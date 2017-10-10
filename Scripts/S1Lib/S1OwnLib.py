# -*- coding: utf-8 -*-
"""
S1OwnLibrary.py
Created on Tue Sep 26 18:13:52 2017

@author: cedric
"""

'''
Sentinel-1 own Library
'''
import os, subprocess, shutil

from osgeo import  ogr, osr

from datetime import datetime
'''
Function that take a generic string of available input value and that
return the real OTB argument value (Calibration_Type) and its acronyme Calib_Name
'''
def ReturnRealCalibrationOTBValue(aCalibration_Type='Sigma0'):
    Calibration_Type = ''
    if aCalibration_Type == 'Sigma0':
        Calibration_Type = 'sigma'
        Calib_Name='Sig0'
    elif aCalibration_Type == 'Gamma0':
        Calibration_Type = 'gamma'
        Calib_Name='Gam0'
    elif aCalibration_Type == 'Beta0':
        Calibration_Type = 'beta'
        Calib_Name='Bet0'
    else:
        Calibration_Type='sigma'
        Calib_Name='Sig0'
    
    return Calibration_Type, Calib_Name
    
    
'''
This function list all file with a given file extention in one directory
'''
def GetFileByExtensionFromDirectory(directory = '/Dir/To/Scan', filt_ext = 'abs.safe'):
    list_file = []
    for root, dirnames, filenames in os.walk(directory):
        for filename in filenames:
            if filename.endswith((filt_ext)):
                list_file.append(os.path.join(root, filename))

    list_file.sort()
    return list_file

'''
This function return the relative orbit of a raw S1 manifest file
'''
def getRelativeOrbit(aManifestFile):
    file = open(aManifestFile,"r")
    for line in file:
        if "<safe:relativeOrbitNumber type=\"start\">" in line:
            RelativeOrbit = int(line.replace("            <safe:relativeOrbitNumber type=\"start\">","").replace("</safe:relativeOrbitNumber>",""))
    file.close()
    return RelativeOrbit
   
'''
Extract the acquisition day from raw S1 tif file

'''
def getAcqDayFromRawS1TifFile(aFileName):
    aFileName = aFileName.split("/")[-1]
   
    return aFileName.split("_")[4][0:8]

def getDayFromS1FileOrFolder(aFileName):
    aFileName = aFileName.split("/")[-1]


    pos=max([aFileName.find('_2015'),aFileName.find('_2016'),aFileName.find('_2017'),aFileName.find('_2018'),aFileName.find('_2019'),aFileName.find('_2020')])
    return aFileName[pos+1:pos+9]

'''
Extract the acquisition day on processed data
'''
def getDateFromS1Raster(PathToRaster):
    return PathToRaster.split("/")[-1].split("-")[4][0:8]

'''
This function filter the manifest list file based on date extraction of already
processed data
'''
def GetNewDatesFromListFilesInputManifest(aInputList, aInputDir, aOutputDir):
    # List all unique date already processed (output folder)
    AllOutputFiles = GetFileByExtensionFromDirectory(aOutputDir, 'tiff') +  GetFileByExtensionFromDirectory(aOutputDir, 'tif')
   
    # List all unique date in input file
    AllInputFiles = GetFileByExtensionFromDirectory(aInputDir, 'tiff') +  GetFileByExtensionFromDirectory(aInputDir, 'tif')

    ListDates = [getAcqDayFromRawS1TifFile(RastPath) for RastPath in AllOutputFiles ]
    UniqueOutputDates = list(set(ListDates))
    UniqueOutputDates.sort()
   
    ListDates = [getDateFromS1Raster(RastPath) for RastPath in AllInputFiles ]
    UniqueInputDates = list(set(ListDates))
    UniqueInputDates.sort()

    # Get New dates
    NewDates = [x for x in UniqueInputDates if x not in UniqueOutputDates]

    # Filter the data
    NewList = []
    for NDate in NewDates:
        for item in aInputList:
            if NDate in item:
                NewList.append(item)

    return NewList


'''
This function get new date comparing orthorectification folder and Temporal filtering
folder.
'''
def GetNewDatesComparingOrthoFolderAndTempFiltFolder(aOrthoFolder, aTempFiltFolder):
    AllInputFiles = GetFileByExtensionFromDirectory(aOrthoFolder, 'tif')
    
    # Filter to dont take into account the output
    AllInputFiles = [ afile for afile in AllInputFiles if 'TempProcStack' not in afile]

    # Get all acquisition date
    ListInputDates = [getDayFromS1FileOrFolder(RastPath) for RastPath in AllInputFiles ]
    
    # Get all unique dates
    UniqueInputDates = list(set(ListInputDates))
    UniqueInputDates.sort()
    
    # Get Dates from Output
    AllOutputFiles = GetFileByExtensionFromDirectory(aTempFiltFolder, 'tif')
    ListOutputDates = [getDayFromS1FileOrFolder(RastPath) for RastPath in AllOutputFiles ]
    UniqueOutputDates = list(set(ListOutputDates))
    UniqueOutputDates.sort()
    
    # Get New dates
    NewDates = [x for x in UniqueInputDates if x not in UniqueOutputDates]
    
    return NewDates

def GetInputOutputListFilesForTempFiltering(aInput_Data_Folder,aOutput_Data_Folder,
                                        aNewDates,aOutput_in_dB, aTmpDirTempFilt,
                                        aWindow_Temp_Filtering):
    # Get input files to process
    AllTifFiles = GetFileByExtensionFromDirectory(aInput_Data_Folder, 'tif')
    # Dont take into account already filtered dates
    AllTifFile = []
    for NDate in aNewDates:
        for file in AllTifFiles:
            if NDate in file:
                AllTifFile.append(file)
    
    AllCopolFile = [file for file in AllTifFile if ('VV' in file) or  ('HH' in file)]
    AllCrosspolFile = [file for file in AllTifFile if ('HV' in file) or  ('VH' in file)]
    
    if 'nt' in os.name:
        AllCopolFile = [file.replace("\\","/") for file in AllTifFile if ('VV' in file) or  ('HH' in file)]
        AllCrosspolFile = [file.replace("\\","/") for file in AllTifFile if ('HV' in file) or  ('VH' in file)]
    else:
        AllCopolFile = [file for file in AllTifFile if ('VV' in file) or  ('HH' in file)]
        AllCrosspolFile = [file for file in AllTifFile if ('HV' in file) or  ('VH' in file)]
    
    # Create Output list data
    AllOutCopolFile = []
    AllOutCrosspolFile = []
    
    if not aOutput_in_dB:
        aTmpDirTempFilt = aOutput_Data_Folder
    
    for CopolFile in AllCopolFile:
        # DirName of currentfile
        DirFile =  os.path.dirname(CopolFile)
        DirName = os.path.split(DirFile)[1]
    
        # Create Output subfolder
        OutFolder = os.path.join(aTmpDirTempFilt,DirName)
        if not os.path.exists(OutFolder):
            os.makedirs(OutFolder)
    
        FileName = os.path.basename(os.path.splitext(CopolFile)[0])
        TempFilterFileName = os.path.join(OutFolder,FileName + '_TempFilt_W' + str(aWindow_Temp_Filtering) + '.tif')
        AllOutCopolFile.append(TempFilterFileName.replace("\\","/"))
    
    # Test if previous temporal filtering
    '''
    TO CONTROL
    '''
    AllTifFiles = GetFileByExtensionFromDirectory(aOutput_Data_Folder, 'tif')
    TempProcStackFile = [ Queg for Queg in AllTifFiles if 'TempProcStack' in Queg]
    if len(TempProcStackFile) ==2:
        for QFile in TempProcStackFile:
            if ('VV' or 'HH') in QFile:
                CopolQueguaFile = QFile
            else:
                CrosspolQueguaFile = QFile
    else:
        CopolQueguaFile = ''
        CrosspolQueguaFile = ''
    
    
    for CrosspolFile in AllCrosspolFile:
        # DirName of currentfile
        DirFile =  os.path.dirname(CrosspolFile)
        DirName = os.path.split(DirFile)[1]
    
        # Create Output subfolder
        OutFolder = os.path.join(aTmpDirTempFilt,DirName)
        if not os.path.exists(OutFolder):
            os.makedirs(OutFolder)
    
        FileName = os.path.basename(os.path.splitext(CrosspolFile)[0])
        TempFilterFileName = os.path.join(OutFolder,FileName + '_TempFilt_W' + str(Window_Temp_Filtering) + '.tif')
        AllOutCrosspolFile.append(TempFilterFileName.replace("\\","/"))
    

'''
This function apply lee filtering data to all files in a folder
'''
def ApplyLeePreFiltering(aFolderList, aOutputFolder, aWindowSize, aENL, aRam):
    Radius = int(aWindowSize / 2.)
    # Loop thru different S1 data (folder)
    for Folder in aFolderList:
        # List all tif and tiff files
        AllTifFile = GetFileByExtensionFromDirectory(Folder, 'tif')

        InDirName = os.path.split(Folder)[1]

        # Create Output subfolder
        OutFolder = os.path.join(aOutputFolder,InDirName)
        if not os.path.exists(OutFolder):
            os.makedirs(OutFolder)

        for file in AllTifFile:
            FileName = os.path.basename(os.path.splitext(file)[0])
            OutputFile = os.path.join(OutFolder,  FileName + '_SpkLee_W' + str(aWindowSize) + '_NL' + str(aENL) +'.tif')

            OTBLeeFiltering(file,OutputFile, Radius, aENL, aRam)


'''
This function use OTB to apply Lee Filtering
'''
def OTBLeeFiltering(aInputFile,aOutputFile, aLeeRadius, aENL, aRam):
    cmd = "otbcli_Despeckle"
    cmd += " -in " + aInputFile
    cmd += " -filter lee"
    cmd += " -filter.lee.rad " + str(aLeeRadius)
    cmd += " -filter.lee.nblooks " + str(aENL)
    cmd += " -ram " + str(aRam)
    cmd += " -out " + aOutputFile

    # progress.setInfo(cmd)
    # print cmd

    p1 = subprocess.Popen (cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,  stderr=subprocess.PIPE)
    ret= p1.communicate()[1]
'''
This function control that the input vector contain one different path per polygon
if not  stop and return an alert
In addition it control that the Path user request is compliant with path in vector
'''
def CheckAllDifferentRelativeOrbit(aInputShape,aPathFieldName, aUserPathList):
    ExceptionMessage = ''
    # Initiate Boolean alert
    isAllUniquePath = False
    isUserPathCompliantShape = True
    # Convert to string due to Windows problem
    aPathFieldName = str(aPathFieldName)
    
    # Open the shape
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(aInputShape, 0)
    layer = dataSource.GetLayer()
    
    # Loop thru feature and get path value
    AllPathList = []
    for feature in layer:
        AllPathList.append(feature.GetField(aPathFieldName))
    
    # Close the source
    dataSource = None
    
    # Check if we have all different Path
    UniquePath = list(set(AllPathList))
    if len(AllPathList) == len(UniquePath):
        isAllUniquePath = True
    else:
        ExceptionMessage += 'You not have all different Path in your shape, '
    
    # Check if Path user request correspond to the shape
    for UserPath in aUserPathList:
        if UserPath not in AllPathList:
            isUserPathCompliantShape = False
    
    if not isUserPathCompliantShape:
        ExceptionMessage += 'You ask a path not in your shape'
            
    # Exception
    if not isAllUniquePath or not isUserPathCompliantShape:
        raise Exception(ExceptionMessage)

'''
This function return the list path files that intersect the study area AND
that are in the required orbit
'''
def getS1ByTile(aInputShape,rawRasterList, aPathFieldName, aPath):
    # Initialization of current Path
    PolygonPath = 0
    
    # Convert to string due to Windows problem
    aPathFieldName = str(aPathFieldName)
    
    # We open the shape
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(aInputShape, 0)
    layer = dataSource.GetLayer()
    
    # Loop Thru feature
    for currentPath in layer:
        print 
        if currentPath.GetField(aPathFieldName)==aPath:
            PolygonPath = aPath
            break
        
    # We get the geometry of the user polygon
    pathFootPrint = currentPath.GetGeometryRef()
    
    intersectRaster=[]
    for image in rawRasterList:
        # We get relative orbit
        RasterRelOrbit = getRelativeOrbit(image)
        
        NW,NE,SE,SW = getOrigin(image)

        # We create a new empty geometry to store the S1 raw footprint
        poly = ogr.Geometry(ogr.wkbPolygon)
        
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(NW[1], NW[0],0)
        ring.AddPoint(NE[1], NE[0],0)
        ring.AddPoint(SE[1], SE[0],0)
        ring.AddPoint(SW[1], SW[0],0)
        ring.AddPoint(NW[1], NW[0],0)
        
        poly.AddGeometry(ring)
        
        # Generate intercection between user polygon and curent S1 footprint
        intersection = poly.Intersection(pathFootPrint)

        # test if we have intercection
        if intersection.GetArea()!=0 and RasterRelOrbit == PolygonPath:
            intersectRaster.append(image)
            
        poly = None
        ring = None
        
    dataSource = None
    
    return intersectRaster  


'''
This function call ogr in system command to reproject a vector
'''
def ReprojVector(aInputVector, aOutputVector, aEPSG):
    CmdList = [
        'ogr2ogr -overwrite',
        '-f \"ESRI Shapefile\" %s'     % aOutputVector,
        '%s' % aInputVector,
        '-t_srs EPSG:%s'     % str(aEPSG)
    ]
    p1 = subprocess.Popen (' '.join(CmdList), shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,  stderr=subprocess.PIPE)
    ret= p1.communicate()[1]


'''
This function rasterize a vector using pixel size parameter
'''
def Rasterize(aInputVector, aOutputRaster, aOutPixelSize):
    # Multilook
    cmd = "gdal_rasterize "
    cmd += " -burn 1 "
    cmd += " -of GTiff "
    cmd += " -a_nodata 0 "
    cmd += " -a_srs EPSG:3857 "
    cmd += " -tr " + str(aOutPixelSize) + " "+ str(aOutPixelSize)+ " "
    cmd += aInputVector + " "
    cmd += aOutputRaster

    # print cmd
    p1 = subprocess.Popen (cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,  stderr=subprocess.PIPE)
    ret= p1.communicate()[1]
   
'''
This function extract the user polygon using the path number
'''
def CreateShapeFromPath(aInputShape,aFieldName, aFieldValue, aOutputShape):
    # Convert to string due to Windows problem
    aFieldName = str(aFieldName)
    
    # Open shape file
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(aInputShape, 0)
    layer = dataSource.GetLayer()
    
    # Loop thru polygon
    for currentPath in layer:
        if currentPath.GetField(aFieldName)==aFieldValue:
            break
    pathFootPrint = currentPath.GetGeometryRef()
    
    # Save extent to a new Shapefile
    outShapefile = aOutputShape
    LayerName = os.path.basename(os.path.splitext(aOutputShape)[0])
    outDriver = ogr.GetDriverByName("ESRI Shapefile")
    
    # Remove output shapefile if it already exists
    if os.path.exists(outShapefile):
        outDriver.DeleteDataSource(outShapefile)
        
    # Create the output shapefile
    outDataSource = outDriver.CreateDataSource(outShapefile)
    outLayer = outDataSource.CreateLayer(LayerName, geom_type=ogr.wkbPolygon)
    
    # Create the feature and set values
    featureDefn = outLayer.GetLayerDefn()
    feature = ogr.Feature(featureDefn)
    feature.SetGeometry(pathFootPrint)
    outLayer.CreateFeature(feature)
    feature = None
    
    outDataSource = None
    dataSource = None
    
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    
    # Write the prj file
    OutPrjFile = os.path.splitext(aOutputShape)[0] + '.prj'
    srs.MorphToESRI()
    filePrj = open(OutPrjFile, 'w')
    filePrj.write(srs.ExportToWkt())
    filePrj.close()

'''
This function get the coordinate origin from S1 manifeste file
'''
def getOrigin(aManifestFile):
    with open(aManifestFile,"r") as saveFile:
        for line in saveFile:
            if "<gml:coordinates>" in line:
                coor = line.replace("                <gml:coordinates>","").replace("</gml:coordinates>","").split(" ")
                coord = [(float(val.replace("\n","").split(",")[0]),float(val.replace("\n","").split(",")[1]))for val in coor]

    return coord[0],coord[1],coord[2],coord[3]

'''
This function get the name of the data withour .safe
'''
def getBaseDataNameFromS1Folder(aFolderName):
    return '_'.join(aFolderName.split("/")[-1].split("_")[:5])

'''
This function proces the S1 GRD data
'''
def ProcessS1Dataset(aInputManifestList, aOutputDir, aInputShape, aDemDir, aResolution,
                     aCalibration_Type, aNoise, aCalibName, aCalibType,
                     aEPSG, aRam):
                         
    # Create tmp dir
    # Get time
    TimeNow = datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
    TmpDir = os.path.join(aOutputDir, 'tmp' + TimeNow)
    if not os.path.exists(TmpDir):
        os.makedirs(TmpDir)
        
    # Rasterize the polygone if enable option
    # Repro Vector to EPSG 3857
    ReprojVectorPath = os.path.join(TmpDir,'ReprojVector.shp')
    ReprojVector(aInputShape, ReprojVectorPath, 3857)
    ClipRasterBase = os.path.join(aOutputDir,'RasterArea.tif')
    Rasterize(ReprojVectorPath, ClipRasterBase, 100)


    # List all SENTINEL-1 sub directories
    S1DataFolders = [ os.path.dirname(ManFile) for ManFile in aInputManifestList]
 
    # Loop thru different S1 data (folder)
    for Folder in S1DataFolders:
        TiffFiles = GetFileByExtensionFromDirectory(Folder, 'tiff')
        
        SafeName = os.path.split(Folder)[-1]
        FolderName = SafeName.split('.')[0]
        SatelliteFolder = FolderName.split('_')[0]
        DateFolder = FolderName.split('_')[4]
        DateFolder=DateFolder[0:DateFolder.find('T')]
        
        WorkingFolder = os.path.join(aOutputDir,FolderName)
        if not os.path.exists(WorkingFolder):
            os.makedirs(WorkingFolder)
        
        for TiffFile in TiffFiles:
            if 'grd-vh' in TiffFile:
                SuffixOrtho = '_VH_'

            if 'grd-hv' in TiffFile:
                SuffixOrtho = '_HV_'
			
            if 'grd-vv' in TiffFile:
                SuffixOrtho = '_VV_'

            if 'grd-hh' in TiffFile:
                SuffixOrtho = '_HH_'

		
            SuffixOrtho += aCalibName+'_Ortho.tif'

            OutClipedFile = os.path.join(TmpDir,'S1Clip.tif')
            OTBExtractRoi(TiffFile, ClipRasterBase, OutClipedFile, aRam)
            TiffFile = OutClipedFile

            OutputFile = os.path.join(TmpDir,FolderName,FolderName + SuffixOrtho)
            if not os.path.exists(os.path.dirname(OutputFile)):
                os.makedirs(os.path.dirname(OutputFile))

            OutputFileCliped = os.path.join(aOutputDir,FolderName,
                                            SatelliteFolder + '_' + DateFolder + SuffixOrtho)		

            # Do GRD conversion to Ortho (including calibration)
            GRD2Calib_Ortho(TiffFile,OutputFile,aDemDir, aCalibType, aNoise, aRam, TmpDir)

            GdalClipRasterWithVector(OutputFile,aInputShape ,
                                     OutputFileCliped, aResolution,  0, aEPSG, aRam)
        
            if os.path.exists(OutputFile):
                os.remove(OutputFile)
            
    # Del Temporay Dir
    shutil.rmtree(TmpDir)
    if os.path.exists(ClipRasterBase):
        os.remove(ClipRasterBase)
        
    # Check if we need to mosaic concecutive data
    # List all folder
    SubFolders = get_immediate_subdirectories(aOutputDir)
    concatenateImage(SubFolders, aRam)

'''
This function clip rax sentinel-1 data using reference raster (use it extent)
'''
def OTBExtractRoi(aInputFile, aInputRefFile, aOutputFile, aRam):
    cmd = "otbcli_ExtractROI -in "
    cmd += aInputFile + " "
    cmd += " -ram " + str(aRam)
    cmd += " -out "
    cmd += aOutputFile
    cmd += " -mode fit "
    cmd += " -mode.fit.ref " + aInputRefFile

    # progress.setInfo(cmd)
    # print cmd

    p1 = subprocess.Popen (cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,  stderr=subprocess.PIPE)
    ret= p1.communicate()[1]
    

'''
This function run calibration and next orthorectification
'''
def GRD2Calib_Ortho(aInputFile,aOutputFile, aDemDir, aCalibration_Type, aNoise, aRam, aTmpDir):
    TimeNow = datetime.now().strftime("%Y_%m_%d_%H_%M_%S")

    # Output Intermediate Calibrated file
    OutputCalibratedFile = os.path.join(aTmpDir, TimeNow + "_Calib.tif")

    # Calibration
    OTBSARCalibration(aInputFile, OutputCalibratedFile, aNoise, aCalibration_Type, aRam)

    # Orthorectification
    OTBOrthorectification(OutputCalibratedFile, aOutputFile, aDemDir, aRam)

    if os.path.exists(OutputCalibratedFile):
        os.remove(OutputCalibratedFile)
            

'''
This function calibrate SAR file
aNoise   --> to remove (or not) noise  - True or False
aCalibration_Type    --> to chose aCalibration_Type      - sigma/gamma/beta/dn
'''
def OTBSARCalibration(aInputFile, aOutputFile, aNoise, aCalibration_Type, aRam):
    cmd = "otbcli_SARCalibration -in "
    cmd += aInputFile + " "
    cmd += " -ram " + str(aRam)
    cmd += " -out "
    cmd += aOutputFile
    cmd += " -noise " + aNoise
    cmd += " -lut " + aCalibration_Type
 
    # progress.setInfo(cmd)
    # print cmd
 
    p1 = subprocess.Popen (cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,  stderr=subprocess.PIPE)
    ret= p1.communicate()[1]
    
'''
This function do orthorectification of any optical/radar file
aDEMDir  --> directory that contain DEM

opt.gridspacing is not in parameter function but need to be near 4 to 10 times of output pixel size
For example, with 10m output pixel size, you can chose bettween 40 to 100.
More less is the value more accurate is the results but more long
In addition this value have to be linked to DEM pixel size
'''
def OTBOrthorectification(aInputFile, aOutputFile, aDEMDir , aRam):
    cmd = "otbcli_OrthoRectification -io.in "
    cmd += aInputFile + " "
    cmd += " -opt.ram " + str(aRam)
    cmd += " -io.out "
    cmd += aOutputFile
    cmd += " -elev.dem " + aDEMDir
    cmd += " -elev.geoid " + "./Data/Geoid/egm96.grd"
    cmd += " -opt.gridspacing 40 "

    # progress.setInfo(cmd)
    print cmd

    p1 = subprocess.Popen (cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,  stderr=subprocess.PIPE)
    ret= p1.communicate()[1]

def GdalClipRasterWithVector(aInputRaster, aInputVector, aOutputFile, aOutputResolution, aSRCNoData, aOutEPSG, aRam):
    cmd = "gdalwarp "
    cmd += "--config GDAL_CACHEMAX " + str(aRam) + " -multi -wo NUM_THREADS=val/ALL_CPUS "
    cmd += " -tr " + str(aOutputResolution) + ' ' + str(aOutputResolution) + ' '
    cmd += " -t_srs " + str(aOutEPSG)
    cmd += " -r average -q -multi -crop_to_cutline -cutline "
    cmd += aInputVector
    cmd += ' -co COMPRESS=DEFLATE -co PREDICTOR=2 '
    cmd += " -dstnodata " + str(aSRCNoData)
    cmd += " -of GTiff " + aInputRaster + ' ' + aOutputFile

    # progress.setInfo(cmd)
    # print cmd

    p1 = subprocess.Popen (cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,  stderr=subprocess.PIPE)
    ret= p1.communicate()[1]
    

'''
This function list this immediate subdirectories in one given directory
'''
def get_immediate_subdirectories(a_dir):
    return [os.path.join(a_dir, name) for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]
                
                
def concatenateImage(aInputFolderList, aRam):
    # MainFolder = 
    FirstFile = aInputFolderList[0]

    MainFolder =  list(os.path.split(FirstFile)[:-1])[0]
    aInputFolderList.sort()

    # List all date
    AcqDateList = [getDayFromS1FileOrFolder(FolderName) for FolderName in aInputFolderList]
    UniqueAcqDateList = list(set(AcqDateList))
    for AcDate in UniqueAcqDateList:
        FileToConcatenate = [folder  for folder in aInputFolderList if AcDate in folder]
        if len(FileToConcatenate) >= 2:
            # Create new Directory (same name without hour)
            SFolder = os.path.split(FileToConcatenate[0])
            SFolder = SFolder[-1].split('_')
            SFolder[4] = AcDate
            SFolder = '_'.join(SFolder[0:5])
            NewOutDir = os.path.join(MainFolder,SFolder)
            if not os.path.exists(NewOutDir):
                os.makedirs(NewOutDir)

            # List all tiff file
            TiffFiles = []
            for FileToConcat in FileToConcatenate:
                TiffFiles += GetFileByExtensionFromDirectory(FileToConcat, '.tif')

            # List copol file (VV or HH)
            CopolFile = [ item for item in TiffFiles if ('_VV_' or '_HH_') in item ]
            FileName = os.path.basename(os.path.splitext(CopolFile[0])[0])
            FileName = os.path.split(FileName)
            FileName = FileName[-1].split('_')
            FileName[1] = AcDate
            FileName = '_'.join(FileName)


            OutputConcatCopolFile = os.path.join(NewOutDir,FileName +'.tif')
            # Concatenate Copol
            OTBConcatenate(CopolFile, OutputConcatCopolFile, aRam)

            # List crosspol file (VH or HV)
            CrosspolFile = [ item for item in TiffFiles if ('_VH_' or '_HV_') in item ]
            FileName = os.path.basename(os.path.splitext(CrosspolFile[0])[0])
            FileName = os.path.split(FileName)
            FileName = FileName[-1].split('_')
            FileName[1] = AcDate
            FileName = '_'.join(FileName)


            OutputConcatCrosspolFile = os.path.join(NewOutDir,FileName +'.tif')
            # Concatenate Copol
            OTBConcatenate(CrosspolFile, OutputConcatCrosspolFile, aRam)

            # delete the 2 folders
            for file in CopolFile + CrosspolFile:
                if  os.path.exists(os.path.dirname(file)):
                    shutil.rmtree(os.path.dirname(file))
                    

def OTBConcatenate(aInputRasterList, aOutputRaster, aRam):
    cmd = "otbcli_BandMath -il "
    for rast in aInputRasterList:
        cmd += rast + ' '
    cmd += " -ram " + str(aRam)
    cmd += " -out "
    cmd += aOutputRaster
    cmd += " -exp "
    cmd += "\"max("
    for i in range(len(aInputRasterList)):
        cmd += 'im' + str(i+1) + 'b1,'
    cmd = cmd[:-1]
    cmd += ")\""

    # print cmd
    p1 = subprocess.Popen (cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,  stderr=subprocess.PIPE)
    ret= p1.communicate()[1]