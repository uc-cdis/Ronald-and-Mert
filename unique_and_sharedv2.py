#!/usr/bin/python

# python /home/ubuntu/git/Ronald-and-Mert/unique_and_sharedv2.py /mnt/single_projects

import csv
import os
import sys

script, directory = sys.argv

def find_common_metadata(mydir = ".", nameFilter = "cor_summary_spearman.txt.filtered.PCO:15.r:0.6.p:0.001"):
	"""
	Iterates through a directory that contains subdirectories with each project
	and produces a file with the significant metadata for each project if
	the filtered file exists (the file that's produced by filter_cor_matrix.r)

	Args:
		mydir: The directory that contains all the project subdirectories
		nameFilter: Part of the name of the filtered file that's spit out
					by filter_cor_matrix.r
					(example: cor_summary_spearman.txt.filtered.PCO:15.r:0.6.p:0.001)

	Returns:
		(None)
		Exports two tables: shared_table.*filter*.tsv and unique_table.*filter*.tsv
		shared_table: table of the pieces of metadata that passed the filter along
					with the number of a projects and which projects were in
					the filtered table
		unique_table: table of the pieces of metadata that passed the filter for
					only a single project
	"""
    #nameFilter = "cor_summary_spearman.txt.filtered.PCO:15.r:0.6.p:0.001"
    #nameFilter = "cor_summary_spearman.txt.filtered.PCO:50.r:0.6.p:0.01"
    # get all filter files
    metadataDict = {}
    projDict = {}
    projIndex = 1
    projList = []
    for root, dirnames, filenames in os.walk(mydir):
        for project in dirnames:
            projList.append(project)
			# projDict used to associate a numerical index with each project
            projDict[project] = projIndex
            projIndex += 1
			# use nameFilter to find the right filtered file
            allfiles = os.listdir(os.path.join(mydir, project))
            for f in allfiles:
                if nameFilter in f:
                    myfile = f
            filteredFileName = os.path.join(mydir, project, myfile)
			# metadataList = [[metadata1, value1], [metadata2, value2], ...]
            metadataList = find_metadata(filteredFileName)
			# adds to a metadataDict with metadata:
			# 	{metadata1: [[proj1, value1], [proj2, value2], ...], metadata2: ...}
            # if list isn't empty
            if metadataList:
                for x in metadataList:
                    # if metadata is already in metadataDict
                    if x[0] in metadataDict:
                        metadataDict[x[0]].append([project, x[1]])
                    else:
                        metadataDict[x[0]] = [[project, x[1]]]

	# need to change the name manually all the time: TODO FIX
    shared = open(os.path.join(mydir, "shared_table.PCO:50.r:0.6.p:0.01.tsv"), 'w')
    unique = open(os.path.join(mydir, "unique_table.PCO:50.r:0.6.p:0.01.tsv"), 'w')
    sharedWriter = csv.writer(shared, delimiter='\t')
    uniqueWriter = csv.writer(unique, delimiter='\t')

    # write table header (all the projects)
    projList = ['', 'Frequency', 'List of Projects'] + projList
    sharedWriter.writerow(projList)
    uniqueWriter.writerow(projList)
    for key, value in metadataDict.iteritems():
		# creates a list "row" with the elements in every row
		#	row[1] is the number of projects where this piece of metadata is signficant
		#	row[2] is a csv of the projects where the metadata is significant
        row = ['NA'] * (projIndex + 2)
        row[1] = len(value)
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
	"""
	Process correlation matrix that has been processed by filtere_cor_matrix.r
	Args:
		filtered: File that shows which pieces of metadata fit filtered requirements
				(file exported by filter_cor_matrix.r)

	Returns:
		List of significant metadata and the PCO:R:P value:
			[[metadata1, value1], [metadata2, value2], ...]
	"""
    f = open(filtered, 'r')
    fReader = csv.reader(f, delimiter = '\t')
    metadataList = []
    for row in fReader:
        if len(row) >= 2 and row[1] != 'NA':
            # [[metadata1, value1], [metadata2, value2], ...]
            metadataList.append([row[0], row[1]])
    return metadataList

find_common_metadata(directory)
