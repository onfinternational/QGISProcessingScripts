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
import os, shutil, sys
from datetime import datetime
import numpy as np
from inspect import getsourcefile


'''
This function join path always using unix sep
'''
def modOsJoinPath(alistToJoin):
    joinedPath = os.path.join(*alistToJoin).replace("\\","/")
    
    return joinedPath


'''
This function return the qgis script folder
'''
def getQGISProcessingScriptFolder():
    from qgis.core import QgsApplication
    QGISProcessingScriptFolder = os.path.dirname(QgsApplication.qgisSettingsDirPath())
    QGISProcessingScriptFolder = modOsJoinPath([QGISProcessingScriptFolder,
    'processing', 'scripts'])
    
    return QGISProcessingScriptFolder

'''
This function load the necessary libs in different way
if we are running script in qgis processing or not
'''
def AddS1LibToPath():
    # Boolean to know if in qgis
    import qgis.utils
    inqgis = qgis.utils.iface is not None
    
    if inqgis:
        QGISProcessingScriptFolder = getQGISProcessingScriptFolder()
        
        # Create the S1Lib lib folder path
        ScriptPath = modOsJoinPath([QGISProcessingScriptFolder, 'S1Lib'])
    else:
        LocalScriptFileDir = os.path.dirname(os.path.abspath((getsourcefile(lambda:0)))).replace("\\","/")
        # Create the S1Lib lib folder path
        ScriptPath = modOsJoinPath([LocalScriptFileDir, 'S1Lib'])  
        
    # Add path to sys
    sys.path.insert(0,ScriptPath)
    
    

# Load OTB Libs
AddS1LibToPath()

from S1OwnLib import (get_immediate_subdirectories,
                          GetNewDatesComparingOrthoFolderAndTempFiltFolder,
                          ApplyLeePreFiltering,
                          GetInputOutputListFilesForTempFiltering,
                          GenerateDualPolColorcompositiondB,
                          GenerateDualPolColorcompositionInt,
                          TileTemporalFilteringRIOS,
                          )

'''
NOTES
Things to do

'''

##Sentinel-1 Deforestation Process V2=group
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
Input_Data_Folder = '/media/cedric/CL/ONFGuyane/Data/Sentinel1/TestScript/Ortho/p120'

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
Output_Data_Folder = '/media/cedric/CL/ONFGuyane/Data/Sentinel1/TestScript/Filt_LibRIOS'


'''
Amount of allocate ram
In case of use half of the available memory (not the physical memory)
'''
Ram = 3000


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
TmpDir = modOsJoinPath([Output_Data_Folder, 'tmp' + TimeNow])
if not os.path.exists(TmpDir):
    os.makedirs(TmpDir)


# Create Temp Filtering Tmp Dir
TmpDirTempFilt = modOsJoinPath([TmpDir, "TempFilt"])
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
                                  
# Test if new date if not exit     
if len(AllInCopolList) == 0:
    ExceptionMessage = 'No new dates to process'
    # Delete Tmp dir
    shutil.rmtree(TmpDir)
    
    raise Exception(ExceptionMessage)
                                  
# Apply the temporal filtering
NumDate = len(AllInCopolList)
# Estimate the size of the block to use based on user available ram parameter
BlockSize = int(np.sqrt(float(Ram) * np.power(1024,2) /(4. * 2. *(2. * float(NumDate + 1.)))))


#start = time.time()
# Usion RIOS lib Ongoing DEV                    
TileTemporalFilteringRIOS(Input_Data_Folder, AllInCopolList, AllOutCopolList,
                      CopolQueguanFile, BlockSize,
                      Spatial_Window_Size_for_Temporal_Filter, Output_Data_Folder)

# Usion RIOS lib Ongoing DEV  
TileTemporalFilteringRIOS(Input_Data_Folder, AllInCrosspolList, AllOutCrosspolList,
                      CrosspolQueguanFile, BlockSize,
                      Spatial_Window_Size_for_Temporal_Filter, Output_Data_Folder)
#end = time.time()


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

#print 'Temps execution co et cross pol filtering', end - start