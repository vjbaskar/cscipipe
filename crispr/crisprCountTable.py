#!/usr/bin/env python3
import os
import argparse

def getCounts(fileName):
	i=0
	count_dict = dict()
	with open(fileName, mode = "r") as file:
		for line in file:
			line = line.strip()
			f = line.split(sep="\t")
			lastCol = len(f) - 1
			grna_id = f[0]
			gene = f[1]	
			sequence = f[2]
			libCount = f[3]	
			exptCount = f[lastCol]
			if i == 0:
				exptname = exptCount
				count_dict[exptname] = {}
			else:
				if grna_id != "UN":
					count_dict[exptname][sequence] = exptCount
			i += 1
	return(count_dict)

def getFileList(countlist):
	"""
	Given a file name, returns a 2D list
	"""
	fileList = []
	with open(countlist, mode = "r") as file:
		for line in file:
			fileList.append(line.strip().split())
	return(fileList)

def getCountDictionary(fileList):
	"""
	Given a 2D list
	Loops through the list
	Opens the list[1] and produces a dictionary
	Creates a dictOfdict { expt_name: { grna_seq1: count, grna_seq2: count, ... } }
	"""
	countDict = {}
	for a in fileList:
			name = a[0]
			fileName = a[1]
			print("..." + name)
			temp = getCounts(fileName)
			countDict.update(temp)
	return(countDict)



if __name__ == "__main__":
	# Input Params
	try:
		countlist = os.environ['inputFile']
		metadatafile = os.environ['organism']
		outptFile=os.environ['outpt']
	except:
		print("""
		++ Reqd:
		-f|--inputFile: File containing input files. Typically output from crisprCounter (see below)
		-o|--organism: File containing gRNA sequences [ see ...SHARED/reference/crispr/libs/]
		-O|--outpt): Output count file [outpt.txt]

		++ inputFile format ++
		name1	name1.counts.tsv
		name2	name2.counts.tsv
		""")
		exit(1)	
	#countlist = "samples.counttable.txt"
	#metadatafile = "mouseV2.txt"
	
	
	print("Reading count file list")
	fileList = getFileList(countlist)
	
	print("Reading files ...")
	countDict = getCountDictionary(fileList)

	# Collect sample names
	columns = [ k for k in countDict.keys() ]

	
	# Read in library meta data
	metadata = getFileList(metadatafile)

	with open(outptFile, mode = "w") as of:
		# Print header
		header = metadata[0] + columns
		of.writelines("\t".join(header) + "\n")
		# Print counts for all the gRNAs in the meta file
		for i in metadata[1:]:
			grna_id = i[0]
			gene = i[1]
			sequence = i[2]
			libCount = i[3]
			rowCounts = []
			for expts in columns:
				temp = countDict[expts].get(sequence, "NA")
				rowCounts.append(temp)

			temp = "\t".join(i) + "\t" +  "\t".join(rowCounts)
			of.writelines(temp + "\n")
