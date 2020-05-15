#!/usr/bin/env python3

import sys
import re
import glob
import prettytable
import pandas as pd
import argparse
import os

def readFile(filename):
	fileContents = list()
	with open(filename, "r") as f:
		for line in f:
			line = line.strip()
			fileContents.append(line)
	return fileContents

def getStatusLine(fileContents):
	totallines = len(fileContents)
	statusLine = -1
	for i in range(totallines):
		result = re.match(r"Resource", fileContents[i]) 
		if result:
			statusLine = i - 2
	return(statusLine)

def splitStrip(i):
	temp = i.split(sep=":")[1]
	temp = temp.strip()
	return temp

def getJobDetails(fileContents, lineNumber):
	lines = fileContents[lineNumber:]
	status = lines[0]
	for i in lines:
		if re.match(r"CPU", i):
			cpu_time = splitStrip(i)
		if re.match(r"Max Memory", i):
			max_mem = splitStrip(i)
		if re.match(r"Total Requested Memory", i):
			total_mem = splitStrip(i)
		if re.match(r"Max Processes", i):
			max_proc = splitStrip(i)
		if re.match(r"Max Threads", i):
			max_threads = splitStrip(i)
		if re.match(r"Run time", i):
			run_time = splitStrip(i)
	x = {'cpu_time':cpu_time, 
	'status': status,
	'max_mem': max_mem,
	'total_mem': total_mem,
	'max_proc': max_proc,
	'max_threads': max_threads,
	'run_time': run_time
	}
	return(x)

def getStartEnd(fileContents):
	for i in fileContents:
		if	re.match(r"Started at", i):
			start = i.replace("Started at", "")
		if	re.match(r"Terminated at", i):
			end = i.replace("Terminated at", "")
	x = {"start": start, "end":end}
	return(x)

def pullOutJobData(fileName):
	#print(f"Pulling out job data from ... {fileName}", file = sys.stderr)
	fileContents = readFile(fileName) # read file as a list
	lineNumber = getStatusLine(fileContents) # get status line
	if lineNumber == -1:
		job_details = {"status": "running"}
		return(job_details)
	job_status = fileContents[lineNumber]
	job_start_end = getStartEnd(fileContents)
	job_details = getJobDetails(fileContents, lineNumber)

	if not re.match(r"Successfully completed", job_status):
		jminus1=fileContents[lineNumber -  1]
		job_status = job_status + " - " + jminus1

	job_details.update(job_start_end)
	job_details.update({"status": job_status})
	return(job_details)

class job:
	counter = 0
	def __init__(self, fileName):
		self.fileName = fileName
		temp = pullOutJobData(fileName)
		self.status = temp['status']
		if self.status == "running":
			return
		
		self.cpu_time = temp['cpu_time']
		self.max_mem = temp['max_mem']
		self.total_mem = temp['total_mem']
		self.max_proc = temp['max_proc']
		self.run_time = temp['run_time']
		self.start = temp['start']
		self.end = temp['end']
		job.counter += 1

	def details(self):
		job_details = self.__dict__.items()
		if self.counter == 1:
			for k,v in job_details:
				print("%s" % k, end = "\t")
			print()

		for k,v in job_details:
			print("%s" % v, end = "\t")
		print()
	
	def forTable(self, onlyHeader = False):
		job_details = self.__dict__.items()
		x = list()
		if onlyHeader:
			for k,v in job_details:
				x.append(k)
		if not onlyHeader:
			for k,v in job_details:
				x.append(v)
		return(x)

			
		#print(f"{self.fileName}\t{self.status}\t{self.start}\t{self.end}\t{self.cpu_time}\t{self.max_mem}\t{self.total_mem}")


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("--progArgs", default = "pretty", help="output type: [pretty, csv]")
	parser.add_argument("--comments", default = "", help="filter for names")
	args = parser.parse_args()
	if not os.path.isdir(".bsub/"):
		print("No farm job log found. See if .bsub exists", file = sys.stderr)
		exit(0)
	
	# Search files
	search = ".bsub/*" + args.comments + "*.farm"
	files = glob.glob(search)
	
	c = 0
	lof = list()
	for f in files:
		j = job(f)
		if j.status == "running":
			continue
		if c == 0:
			colnames = j.forTable(onlyHeader = True)
			table = prettytable.PrettyTable(colnames)
			lof.append(colnames)
		l = j.forTable(onlyHeader = False)
		table.add_row(l)
		lof.append(l)
		c += 1
	
	if args.progArgs == "pretty":
		print(table)
	if args.progArgs == "csv":
		df = pd.DataFrame(lof)
		print(df.to_csv(index=False, header = False))



