#!/usr/bin/python

# python /home/ubuntu/git/Ronald-and-Mert/unique_and_sharedv2.py /mnt/single_projects

import csv
import os
import sys

script, directory = sys.argv

def find_common_metadata(mydir):
    #nameFilter = "cor_summary_spearman.txt.filtered.PCO:15.r:0.6.p:0.001"
    nameFilter = "cor_summary_spearman.txt.filtered.PCO:50.r:0.6.p:0.01"
    # get all filter files
    metadataDict = {}
    projDict = {}
    projIndex = 1
    projList = []
    for root, dirnames, filenames in os.walk(mydir):
        for project in dirnames:
            projList.append(project)
            projDict[project] = projIndex
            projIndex += 1
            allfiles = os.listdir(os.path.join(mydir, project))
            for f in allfiles:
                if nameFilter in f:
                    myfile = f
            filteredFileName = os.path.join(mydir, project, myfile)
            metadataList = find_metadata(filteredFileName) 
            # adds to a dictionary with metadata: [[proj1, data1], [proj2, data2]]
            # if list isn't empty
            if metadataList:
                for x in metadataList:
                    # if metadata in list
                    if x[0] in metadataDict:
                        metadataDict[x[0]].append([project, x[1]])
                    else:
                        metadataDict[x[0]] = [[project, x[1]]]
    shared = open(os.path.join(mydir, "shared_table.PCO:50.r:0.6.p:0.01.tsv"), 'w')
    unique = open(os.path.join(mydir, "unique_table.PCO:50.r:0.6.p:0.01.tsv"), 'w')
    sharedWriter = csv.writer(shared, delimiter='\t')
    uniqueWriter = csv.writer(unique, delimiter='\t')
    # write table header
    projList = ['', 'Frequency', 'List of Projects'] + projList
    sharedWriter.writerow(projList)
    uniqueWriter.writerow(projList)
    # for every key in dictionary: spit out relevent projects
    for key, value in metadataDict.iteritems():
        # create list of length proj + 2: col for frequency, col for projects
        row = ['NA'] * (projIndex + 2)
        # frequency
        row[1] = len(value)
        # comma separated value for projects
        # iterate through and append to string
        specificProjList = ''
        for item in value:
            if specificProjList:
                specificProjList += ','
            specificProjList += item[0]
        row[2] = specificProjList
        for pair in value:
            # row[projIndex+2] = data
            row[projDict[pair[0]]+2] = pair[1]
        row[0] = key
        if len(value) == 1:
            uniqueWriter.writerow(row)
        else:
            sharedWriter.writerow(row)

    shared.close()
    unique.close()


def find_metadata(filtered):
    f = open(filtered, 'r')
    fReader = csv.reader(f, delimiter = '\t')
    metadataList = []
    for row in fReader:
        if len(row) >= 2 and row[1] != 'NA':
            # [[metadata1, value1], [metadata2, value2], ...]
            metadataList.append([row[0], row[1]])
    return metadataList

find_common_metadata(directory)
