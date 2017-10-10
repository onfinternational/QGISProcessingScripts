# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 18:14:51 2017

@author: cedric
"""

'''
This script apply calibration and orthorectification process of S1 GRD data
'''

'''
IMPORT
'''
import os


from S1Lib.S1OwnLib import (ReturnRealCalibrationOTBValue,
                                   GetFileByExtensionFromDirectory,
                                   GetNewDatesFromListFilesInputManifest,
                                   ReprojVector,
                                   CheckAllDifferentRelativeOrbit,
                                   getS1ByTile,
                                   CreateShapeFromPath,
                                   ProcessS1Dataset)
#from S1Lib.S1OTBLib import 

# TO DELETE
'''
NOTES
Things to do
    - Optimize the ortorectification by using extent of the polygon to replace clip
'''





##Sentinel-1 Deforestation Process=group
##1 - Calibration and Orthorectification over tiles=name
##Input_Data_Folder=folder
##DEM_Folder=folder
##Input_Polygon_File=vector
##Relative_Orbit_Field_Name=field Input_Polygon_File
##Relative_Orbit_To_Process=string 120-47
##Calibration_Type=string Sigma0
##Output_EPSG=crs
##Output_Resolution=number 10
##Output_Data_Folder=folder
##Ram=number 256
# END TO DELETE



'''
Input OF THE PROGRAM
'''

'''
This string have to contain the folder path that contain all the unziped S1
 data to process
'''
Input_Data_Folder = '/media/cedric/CL/ONFGuyane/Data/Sentinel1/TestScript/Input'

'''
This string have to contain the folder path that contain all dem covering the
study area. Tiff format is required
'''
DEM_Folder = '/media/cedric/CL/ONFGuyane/Data/Sentinel1/TestScript/DEM'

'''
Path of the vector file that contain one or more polygon to process
Each polygone will be process independantly like independant tile

This vector have to contain one field with integer type containing
the Relative orbit number to process for the current polygon
'''
Input_Polygon_File = '/media/cedric/CL/ONFGuyane/Data/Sentinel1/TestScript/cliparea.shp'

'''
Name of the field containing the relative orbit
'''
Relative_Orbit_Field_Name = 'id'


'''
Value of all the relative orbit to process. Have to be separated by '-' like
120-47
'''
Relative_Orbit_To_Process = '120'


'''
String containing the calibration type:
Sigma0
Gamma0
Beta0
'''
Calibration_Type = 'Sigma0'

'''
Output EPSG (only in meter)
'''
Output_EPSG = 'EPSG:3857'

'''
Output resolution in meter
'''
Output_Resolution = '10'

'''
Output folder that will contain one subfolder for each processed polygon and
each relative orbit
'''
Output_Data_Folder = '/media/cedric/CL/ONFGuyane/Data/Sentinel1/TestScript/Ortho'


'''
Amount of allocate ram
In case of use half of the available memory (not the physical memory)
'''
Ram = 2000


'''
THE PROGRAM ITSELF
'''
# Internal variable
Noise = 'false'
'''
Main step
1 - Transform Calibration user input into OTB command and get its acronym
2 - List ALL Sentinel-1 Manifest file
3 - Filter the one that not already processed in the output folder
4 - Scan path for all polygon and control that there is different path for all
polygon 
5 - Loop thru relative orbit given by the user and 
    5a - Loop thru polygon
        For each relative orbit loop thru all the polygon of this orbit
        - Create on subfolder per polygon
        - Filter all the data to take the one which intersect the study area
        - Process the data
'''



# 1 - Transform Calibration user input into OTB command and get its acronym
OTBCalibrationType, CalibrationName = ReturnRealCalibrationOTBValue(Calibration_Type)

# 2 - List ALL Sentinel-1 Manifest file
# List all SENTINEL-1 manifest.safe files
ManifestFiles = GetFileByExtensionFromDirectory(Input_Data_Folder, 'manifest.safe')

# 3 - Filter the one that not already processed in the output folder
ManifestFiles = GetNewDatesFromListFilesInputManifest(ManifestFiles, Input_Data_Folder, Output_Data_Folder)

# 4 - Scan path for all polygon and control that there is different path for all
# polygon
# Get Path name list user
PathUserList = Relative_Orbit_To_Process.split('-')
PathUserList = [int(path) for path in PathUserList]
CheckAllDifferentRelativeOrbit(Input_Polygon_File,Relative_Orbit_Field_Name, PathUserList)

# Reproject Shapefile to 4326 to be compliant with S1 raw data
Input_Polygon_FileEPSG4326 = Input_Polygon_File.replace('.shp', 'EPSG4326.shp')
ReprojVector(Input_Polygon_File, Input_Polygon_FileEPSG4326, 4326)


# 5 - Loop thru relative orbit given by the user and 
#    5a - Loop thru polygon
for userPath in PathUserList:
    # Filter files that intersect the required path
    intersectRaster = getS1ByTile(Input_Polygon_FileEPSG4326,
                                  ManifestFiles,
                                  Relative_Orbit_Field_Name,
                                  userPath)
    # If intersect is not empty we create output directory
    if len(intersectRaster) > 0:
        PathDir = os.path.join(Output_Data_Folder, 'p' + str(userPath))
        if not os.path.exists(PathDir):
            os.makedirs(PathDir)
    else: # if no data go to next path
        continue
    
    # Create Shape file of the current path
    PathShape = os.path.join(PathDir, 'p' + str(userPath) + '.shp')
    CreateShapeFromPath(Input_Polygon_FileEPSG4326,
                        Relative_Orbit_Field_Name,
                        str(userPath),
                        PathShape)
                        
    
    
    # Run the process
    ProcessS1Dataset(intersectRaster,
                     PathDir,
                     PathShape,
                     DEM_Folder,
                     Output_Resolution,
                     Calibration_Type,
                     Noise,
                     CalibrationName,
                     OTBCalibrationType,
                     Output_EPSG,
                     Ram)