#!/usr/bin/python

import csv
import os

def find_common_metadata(mydir):
    # get all filter files
    projectList = os.listdir(mydir)
    metadataDict = {}
    for project in projectList:
        filteredFileName = os.path.join(project, "cor_summary_spearman.txt.filtered")
        metadataList = find_metadata(filteredFileName)    
        # adds to a dictionary with metadata: [list of projects]
        for x in metadataList:
            if x in metadataDict:
                metadataDict[x] = metadataDict[x].append(project)
            else:
                metadatadict[x] = [project]
    shared = open("shared_table.tsv", 'w')
    unique = open("unique_table.tsv", 'w')
    sharedWriter = csv.writer(shared, delimiter='\t')
    uniqueWriter = csv.writer(unique, delimiter='\t')
    # for every key in dictionary: spit out relevent projects
    for key, value in metadataDict.itertems():
        if len(value) == 1:
            uniqueWriter.writerow(value)
        else:
            sharedWriter.writerow(value)

    shared.close()
    unique.close()


def find_metadata(filtered):
    f = open(filtered, 'r')
    fReader = csv.reader(f, delimiter = '\t')
    metadataList = []
    for row in fReader:
        if len(row) >= 2 and row[1] != 'NA':
            metadataList.append(row[0])
    return metadataList
