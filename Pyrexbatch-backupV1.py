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
import logging
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
def main():
  outPath = r''
  #read Params file of Pyradiomics
  progress_filename = os.path.join(outPath, 'pyrex_log.txt')
  # Configure logging
  rLogger = logging.getLogger('radiomics')

  # Set logging level
  # rLogger.setLevel(logging.INFO)  # Not needed, default log level of logger is INFO

  # Create handler for writing to log file
  handler = logging.FileHandler(filename=progress_filename, mode='w')
  handler.setFormatter(logging.Formatter('%(levelname)s:%(name)s: %(message)s'))
  rLogger.addHandler(handler)

  # Initialize logging for batch log messages
  logger = rLogger.getChild('batch')

  # Set verbosity level for output to stderr (default level = WARNING)
  #radiomics.setVerbosity(logging.INFO)

  #logger.info('pyradiomics version: %s', radiomics.__version__)
#  logger.info('Reading Params file for pyradiomics')

  # Reading Params file for pyradiomics
  try:
  	paramsFile = os.path.join(os.getcwd(),'ParamsSettings','Pyradiomics_Params.yaml')
  except Exception:
    logger.error('Could not find params file of Pyradiomics!', exc_info=True)
    exit(-1)
#  logger.info('Reading Params file of Pyrex')
  # Reading Params file of Pyrex
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
  except Exception:
    logger.error('Could not find params file of Pyrex!', exc_info=True)
    exit(-1)
  logger.info('Reading Parameters Done')


  CaseID = os.listdir(myWorkingDirectory)
#  logger.info('Parsing DICOM files and RTSTRUCT in working directory')
  #xnat_collection(myWorkingDirectory,collectionURL,myProjectID)
  Img_path,RT_path = ParseStructure(myWorkingDirectory) #detect the path of Image and RTstruct
  logger.info('DICOM and RTSTRUCT Parsing Done')


  ######## create a pandas data frame for data
  flists = pandas.DataFrame(data= {'CaseID':CaseID}).T
  logger.info('Starting Pyrex')
  for entry in flists:
  #for entry in range(8,21):
      results = pandas.DataFrame()
#      logger.info('processing patient: %s', CaseID[entry])
      mask_vol=PyrexReader.Read_RTSTRUCT(RT_path[entry])
      logger.info('Loading RTSTRUCT: %s', RT_path[entry])
      M=mask_vol[0]
      target = []
      for j in range(0,len(M.StructureSetROISequence)):
          target.append(M.StructureSetROISequence[j].ROIName)
      logger.info('VOI: %s', RT_path[entry])
     
      for k in range(0,len(target)):
          try:
              featureVector = flists[entry]
              Image,Mask = PyrexReader.Img_Bimask(Img_path[entry],RT_path[entry],target[k])
              logger.info('%d Processing Radiomics on %s of Patient (%s)',k, target[k],CaseID[entry])
              if export_format == 'csv': # sava results in csv
                  try:
                      result = pandas.Series(PyrexWithParams.CalculationRun(Image,Mask,paramsFile))
                      featureVector = featureVector.append(result)
                      featureVector.name = k
                      results = results.join(featureVector, how='outer')
                      Image = []
                      Mask= []
                      result=[]
                  except Exception:
                    logger.error('FEATURE EXTRACTION FAILED for CSV output:', exc_info=True)
              else:# save results in triple stroe
                  try:
                    featureVector = PyrexWithParams.CalculationRun(Image,Mask,paramsFile) #compute radiomics
                    featureVector.update({'patient':CaseID[entry],'contour':target[k]}) #add patient ID and contour                
                    PyrexOutput.RadiomicsStore(featureVector,exportDir,CaseID[entry],target[k],export_format,export_name) #store radiomics locally with a specific format 
                    logger.info('Extraction complete, writing rdf')
                  except Exception:
                    logger.error('FEATURE EXTRACTION FAILED for RDF output:', exc_info=True)
          except Exception:
              logger.info('FEATURE EXTRACTION FAILED:')
          logger.info('-------------------------------------------')
          print(k,CaseID[entry],target[k])

      logger.info('Extraction complete, writing CSV')          
      outputFilepath = os.path.join(exportDir,CaseID[entry]+'.csv')
      results.T.to_csv(outputFilepath, index=False, na_rep='NaN')
      logger.info('CSV writing complete')
      logger.info('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
      
if __name__ == '__main__':
  main()   