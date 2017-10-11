# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 15:51:57 2017

@author: cedric
"""


'''
This script do temporal filtering over sentinel-1 data
'''

'''
IMPORT
'''
import os, shutil
from datetime import datetime
import numpy as np

'''
from S1Lib.S1OwnLib import (ReturnRealCalibrationOTBValue,
                                   GetFileByExtensionFromDirectory,
                                   GetNewDatesFromListFilesInputManifest,
                                   ReprojVector,
                                   CheckAllDifferentRelativeOrbit,
                                   getS1ByTile,
                                   CreateShapeFromPath,
                                   ProcessS1Dataset)
'''

from S1Lib.S1OwnLib import (get_immediate_subdirectories,
                            GetNewDatesComparingOrthoFolderAndTempFiltFolder,
                            ApplyLeePreFiltering,
                            GetInputOutputListFilesForTempFiltering,
                            TileTemporalFiltering,
                            GenerateDualPolColorcompositiondB,
                            GenerateDualPolColorcompositionInt,
                            )
# TO DELETE
'''
NOTES
Things to do

'''

##Sentinel-1 Deforestation Process=group
##2 - Temporal Filtering over area=name
##Input_Data_Folder=folder
##Spatial_Window_Size_for_Temporal_Filter=number 11
##Apply_Lee_Pre_Filtering=boolean False
##Spatial_Window_Size_for_Lee_Filter=number 5
##Looks_Number_for_Lee_Filter=number 5
##Output_in_dB=boolean True
##Output_Data_Folder=folder
##Ram=number 256


'''
Input OF THE PROGRAM
'''

'''
This string have to contain the folder path that contain all the Orthorectified
and calibrate S1 data from previous step
'''
Input_Data_Folder = '/media/cedric/CL/ONFGuyane/Data/Sentinel1/TestScript/Ortho/p120/Full'

'''
This integer contain the window size for the temporal filter
'''
Spatial_Window_Size_for_Temporal_Filter = 11

'''
This boolean enable or not to apply lee pre filtering
'''
Apply_Lee_Pre_Filtering = False

'''
Integer Window size of the Lee pre filtering, used only if enable
'''
Spatial_Window_Size_for_Lee_Filter = 5


'''
Integer containing the equivalent number of looks of the input data
'''
Looks_Number_for_Lee_Filter = 5


'''
Enable to export the data in decibel (if not intensity)
'''
Output_in_dB = True


'''
Output folder that will contain one subfolder for each processed polygon and
each relative orbit
'''
Output_Data_Folder = '/media/cedric/CL/ONFGuyane/Data/Sentinel1/TestScript/Filt'


'''
Amount of allocate ram
In case of use half of the available memory (not the physical memory)
'''
Ram = 2000


'''
THE PROGRAM ITSELF
'''
# Internal variable


'''
Main step
1 - Get new date to process
2 - If enable apply Lee pre filtering
3 - Apply temporel filtering
4 - Generate color composition
'''

#1 - List all Input files
# List of all tif files

NewDates = GetNewDatesComparingOrthoFolderAndTempFiltFolder(Input_Data_Folder,
                                                            Output_Data_Folder)

# Create tmp dir
# Get time
TimeNow = datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
TmpDir = os.path.join(Output_Data_Folder, 'tmp' + TimeNow)
if not os.path.exists(TmpDir):
    os.makedirs(TmpDir)


# Create Temp Filtering Tmp Dir
TmpDirTempFilt = os.path.join(TmpDir, "TempFilt")
if not os.path.exists(TmpDirTempFilt):
    os.makedirs(TmpDirTempFilt)



# List all SENTINEL-1 sub directories
AllSubFolders = get_immediate_subdirectories(Input_Data_Folder)
# Filter to get only new dates
SubFolders = []
for NDate in NewDates:
    for folder in AllSubFolders:
        if NDate in folder:
            SubFolders.append(folder)


#2 - If enable apply Lee pre filtering
if Apply_Lee_Pre_Filtering:
    ApplyLeePreFiltering(SubFolders, TmpDir, Spatial_Window_Size_for_Lee_Filter,
                         Looks_Number_for_Lee_Filter, Ram)
    # Update Input dir to be Tmp folder
    Input_Data_Folder = TmpDir


#3 - Apply temporel filtering
# Get input output list file to filter
# Copol in one side and cross pol
AllInCopolList,AllOutCopolList, AllInCrosspolList, AllOutCrosspolList, \
CopolQueguanFile, CrosspolQueguanFile\
= GetInputOutputListFilesForTempFiltering(Input_Data_Folder,Output_Data_Folder,
                                        NewDates,Output_in_dB, TmpDirTempFilt,
                                        Spatial_Window_Size_for_Temporal_Filter)
                                        
# Apply the temporal filtering
NumDate = len(AllInCopolList)
# Estimate the size of the block to use based on user available ram parameter
BlockSize = int(np.sqrt(float(Ram) * np.power(1024,2) /(4. * 2. *(2. * float(NumDate + 1.)))))

TileTemporalFiltering(Input_Data_Folder, AllInCopolList, AllOutCopolList,
                      CopolQueguanFile, BlockSize,
                      Spatial_Window_Size_for_Temporal_Filter)

# Apply filtering to crosspol
TileTemporalFiltering(Input_Data_Folder, AllInCrosspolList, AllOutCrosspolList,
                      CrosspolQueguanFile, BlockSize,
                      Spatial_Window_Size_for_Temporal_Filter)

#4 - Generate color composition
SubFolders = get_immediate_subdirectories(TmpDirTempFilt)
if Output_in_dB:
	GenerateDualPolColorcompositiondB(SubFolders, Output_Data_Folder, Ram)
else:
	GenerateDualPolColorcompositionInt(SubFolders, Output_Data_Folder, Ram)

# Delete old Quegan File
if os.path.exists(CopolQueguanFile):
	os.remove(CopolQueguanFile)

if os.path.exists(CrosspolQueguanFile):
	os.remove(CrosspolQueguanFile)

# Delete Tmp dir
shutil.rmtree(TmpDir)