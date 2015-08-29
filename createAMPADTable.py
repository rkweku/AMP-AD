# This program downloads a series of files from Synapse and
# creates one one large table from those files to upload to Synapse
# for querying. Files are merged based on a common identifier
# across the files to merge samples with one another.

import argparse
import pandas as pd
import synapseclient
import sys
import urllib2
import lxml.etree as ET

from Bio import Entrez
from synapseclient import Schema, Column, Table, Row, RowSet, as_table_columns

#########################Execution Variables###################################
# The Synapse ID of the AMPAD project
synID = 'syn2580853'
###############################################################################

# Login to Synapse
syn = synapseclient.Synapse()
syn.login()

def getSynapseFile(synID):
    """
    Download the data from a given synapse ID

    Input:
        synID: The ID of the file to be downloaded
    Returns:
        pandas dataframe of the datafile
    """

    # Download the data from synapse
    entity = syn.get(synID)

    # If the file is a csv file, set delimeter to a comma
    if(entity.path.endswith('.csv')):
        # Save the actual data as a pandas dataframe
        wholeFile = pd.read_csv(entity.path, sep=',')

    # If file is not specified, it is *likely* a tab
    else:
        wholeFile = pd.read_csv(entity.path, sep='\t')

    return(wholeFile)

def createAMPADTable(keyFile, clinicalFile):
    """
    Create the AMP AD table with merged data from keyFile with clinicalFile.
    If any of the supplementary files exist for a particular dataset, change
    the binary classifiers to the synapse ID holding the data and reset 0
    to null for the table.

    Input:
        keyFile: Dataframe with the keys and information regarding what
            exists for each patient
        clinicalFile: Dataframe with clinical data for various patients

    """

    toUpload = []

    clinicalHeader = clinicalFile.columns.values

    #seenList = []
    # Iterate through each project within keyFile
    for i, row in keyFile.iterrows():
        # Create empty list for new row to be added to synapse table
        newRow = []

        # Ignore binary varibles which all end in '_data'
        for item in row.iteritems():
            if(item[0] == 'niagas_data'):
                if(not pd.isnull(row.niagas_data)):
                    newRow.append(arrayExpressionSynID)
                else:
                    newRow.append(float('nan'))

            elif(not item[0].endswith('_data')):
                newRow.append(item[1])
        
        # Check if row has clinical data
        if(row.clinical_data):
            # Create reference to clinicalFile project ID
            clinicalKeyList = clinicalFile['projid']
            
            # get the index of the projID in the clinical file
            index = clinicalKeyList[clinicalKeyList==row.projid].index.tolist()

            if(len(index) == 1):
                index = index[0]
                #seenList.append(row.projid)
                for entry in clinicalFile.iloc[index][1:]:
                    newRow.append(entry)

            # If the length of the idnex is 0, it means the key file thinks
            # there is clinical information for this patient but it does
            # not exist in the clinical file
            elif(len(index) == 0):
                print("Key file indicates that projID %s should have "\
                    "clinical information, but it does not exist in "\
                    "the clinical information file" % row.projid)
                for _ in range(1, len(clinicalHeader)): 
                    newRow.append(float('nan'))

            # If the lengh of index list is greater than 1, that means projID
            # appears more than once in the file. Send warning to user
            else:
                print("projID %s appears more than once in clinical file at "\
                    "positions %s" % (row.projid, index))
                for _ in range(1, len(clinicalHeader)): 
                    newRow.append(float('nan'))

        else:
            for _ in range(1, len(clinicalHeader)): 
                newRow.append(float('nan'))
        
        # Check if row has gwas data
        if(row.gwas_data):
            newRow.append(genotypeSynID)
            newRow.append(imputedGenotypeSynID)
        else:
            newRow.append(float('nan'))
            newRow.append(float('nan'))
        
        if(row.mwas_data):
            newRow.append(methylationSynID)
        else:
            newRow.append(float('nan'))
        
        if(row.mirna_data):
            newRow.append(mirnaSynID)
        else:
            newRow.append(float('nan'))
        
        if(row.mrna_data):
            newRow.append(rnaseqSynID)
        else:
            newRow.append(float('nan'))

        toUpload.append(newRow)

    df = pd.DataFrame(toUpload)
    columns = list(keyFile.columns.values)
    index = columns.index('clinical_data') - 1
    columns.remove('clinical_data')

    idnex = columns.index('gwas_data')
    columns.remove('gwas_data')
    columns.insert(index+1, 'genotype data')
    columns.insert(index+2, 'imputed genotype data')
    
    for i in range(1, len(clinicalHeader)):
        columns.insert(index+i, clinicalHeader[i])

    df.columns = columns

    df.to_csv('mergedTables.csv', encodings='utf-8', index=False)

    print("Uploading to Synapse")
    schema = Schema(name = 'AMP AD Samples Table',
        columns=as_table_columns(df), parent='syn2580853')
    syn.store(Table(schema, df))

def getSynIDs():
    """
    Get the Synapse IDs of the files and folders within the parentID

    Returns:
        Dictionary of all studies with a list of all synIDs and the data
        types associated with each file

    """ 

    # Initialize empty dictionary for the studies
    studies = {}

    # Run a query to find all files that do not have None as a data type or
    # as the study name
    query = "select study, dataType from file where projectId == '%s' and "\
            "dataType != 'None' and study != 'None'" % synID
    results = syn.chunkedQuery(query)

    # Loop through all files to add them to a dictionary to loop through later
    for res in results:
        # Get as strings the study, fileID and data type for each entity
        study = str(res['file.study'][0])
        fileID = str(res['file.id'])
        dataType = str(res['file.dataType'][0])

        # If the study does not exist yet, make an empty dictionary for it
        if(study not in studies.keys()):
            studies[study] = {}

        # If the data type does not exist for the study, make an empty
        # list for it
        if(dataType not in studies[study].keys()):
            studies[study][dataType] = []

        # Append the synapse ID for the data type in this study
        studies[study][dataType].append(fileID)

    return(studies)

def main():

    studies = getSynIDs()

    print studies['ROSMAP'].keys()
    print studies['ROSMAP']
    # Download the key file
    #keyFile = getSynapseFile(keySynID)

    #for study in studies:
        
    # Download clinical file
    clinicalFile = getSynapseFile(clinicalSynID)

    createAMPADTable(keyFile, clinicalFile)
    
#    keyList = keyFile['projid']
#    print(keyList[keyList==3380931].index.tolist())

    

if __name__ == "__main__":
    main()
