"""
###############################
@author: zhenwei.shi, Maastro##
###############################
"""
from __future__ import print_function

import PyrexReader
import PyrexWithParams
import PyrexOutput
import yaml
from PyrexXNAT import ParseStructure, xnat_collection
#import logging


import os
import pandas
'''
Usage for individual case: python HelloPyrexBatchProcessing.py 

Read parameter file of Pyrex:

# - path:
#    - myWorkingDirectory is the root directory where DICOM files are saved.
#    - exportDir is the directory where results are exported.
# - collectionURL: specify the URL of cloud repository, like 'XNAT'.
# - myProject: specify the name of dataset on cloud reposity, like 'stwstrategyrdr'
# - export_format: specify the format of output, such as rdf or csv.
# - export_name: specify the name of result file.
######################################
'''

#read Params file of Pyradiomics
try:
	paramsFile = os.path.join(os.getcwd(),'ParamsSettings','Pyradiomics_Params.yaml')
except:
	print('Error: Could not find params file of Pyradiomics!')


#reading Params file of Pyrex
try:
    inputFile = os.path.join(os.getcwd(),'ParamsSettings','Pyrex_BatchProcessingParams.yaml')
    with open(inputFile,'r') as params:
        p = yaml.load(params)
        
    #read parameters from configuration file.
    myWorkingDirectory = p['PATH']['myWorkingDirectory']
    collectionURL = p['collectionURL'][0]
    myProjectID = p['myProjectID'][0]
    exportDir = p['PATH']['exportDir']
    export_format = p['export_format'][0]
    export_name = p['export_name'][0]
except:
	print('Error: Could not find params file of Pyrex!')

CaseID = os.listdir(myWorkingDirectory)
#xnat_collection(myWorkingDirectory,collectionURL,myProjectID)
Img_path,RT_path = ParseStructure(myWorkingDirectory) #detect the path of Image and RTstruct

######## log file

######## create a pandas data frame for data
flists = pandas.DataFrame(data= {'CaseID':CaseID}).T
results = pandas.DataFrame()
for entry in flists:
#for entry in range(8,21):
    print(CaseID[entry])
    mask_vol=PyrexReader.Read_RTSTRUCT(RT_path[entry])
    M=mask_vol[0]
    target = []
    for j in range(0,len(M.StructureSetROISequence)):
        target.append(M.StructureSetROISequence[j].ROIName)
    print(target)
    
#    imageFilepath = flists[entry]['Img_path']
#    rtFilepath = flists[entry]['RT_path']
#    label = flists[entry].get('Label', None)

    for k in range(0,len(target)):
        try:
            featureVector = flists[entry]
#            if (imageFilepath is not None) and (rtFilepath is not None):
#                featureVector = flists[entry]
            Image,Mask = PyrexReader.Img_Bimask(Img_path[entry],RT_path[entry],target[k])
            print(k,target[k])
            # sava results in csv
            if export_format == 'csv':
                result = pandas.Series(PyrexWithParams.CalculationRun(Image,Mask,paramsFile))
                featureVector = featureVector.append(result)
                featureVector.name = k
                results = results.join(featureVector, how='outer')
                Image = []
                Mask= []
                result=[]
            # save results in triple stroe
            else:
                featureVector = PyrexWithParams.CalculationRun(Image,Mask,paramsFile) #compute radiomics
                featureVector.update({'patient':CaseID[entry],'contour':target[k]}) #add patient ID and contour                
                PyrexOutput.RadiomicsStore(featureVector,exportDir,CaseID[entry],target[k],export_format,export_name) #store radiomics locally with a specific format 
        except Exception:
            print('FEATURE EXTRACTION FAILED:')

    outputFilepath = os.path.join(exportDir,CaseID[entry]+'.csv')
    results.T.to_csv(outputFilepath, index=False, na_rep='NaN')    