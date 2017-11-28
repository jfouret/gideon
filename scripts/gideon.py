#!/usr/bin/env python
version='SEDMATCHVERSION'
scripts='SEDMATCHLIB'
author='Julien Fouret'
contact='julien@fouret.me'

import argparse

class _HelpAction(argparse._HelpAction):

	def __call__(self, parser, namespace, values, option_string=None):
		parser.print_help()
		# retrieve subparsers from parser
		subparsers_actions = [
			action for action in parser._actions
			if isinstance(action, argparse._SubParsersAction)]
			# there will probably only be one subparser_action,
			# but better save than sorry
		for subparsers_action in subparsers_actions:
		# get all subparsers and print help
			for choice, subparser in subparsers_action.choices.items():
				print("\n\n\n\n--------"+("-"*len(choice))+"\n"+"### {} ###".format(choice))+"\n"+"--------"+("-"*len(choice)+"\n")
				print(subparser.format_help())
				print("______________________________________________________________________________")
		parser.exit()

parser = argparse.ArgumentParser(description='Wraper for GIDEON analyses',epilog="Version : "+str(version)+"\nAuthor : "+author+" for more informations or enquiries please contact "+contact,add_help=False,formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-h','--help', action=_HelpAction, help='if you need some help')

subparsers=parser.add_subparsers(help='Sub-commands',dest="analysis")

vhelp="Print version and exit"
vparser=subparsers.add_parser('version',description=vhelp,help=vhelp)

checkhelp="Check if requirements are satisfied"
checkparser=subparsers.add_parser('check',description=checkhelp,help=checkhelp)

genePythiaHelp="Knowledge-driven method to identify a set of genes linked with a thematic (defined by terms). PubMed data-mining"
genePythiaParser=subparsers.add_parser('genePythia',description=genePythiaHelp,help=genePythiaHelp)
genePythiaParser.add_argument('-o','--outDir', metavar='Path',required=True, help="outPut directory")
#parse and gene name will be the prefix
genePythiaParser.add_argument('-geneList', metavar='file',required=True, help="2-columns tabular file without header separated by tabulation. 1st col: gene name | 2nd col: comma-separated alias list")
genePythiaParser.add_argument('-term', metavar='syn1,syn2,syn3',required=True, help="commat-separated list of terms to search (synonyms). Please add\" if a term is composed with spaces" )
genePythiaParser.add_argument('-b','--batch_number', metavar='N',default='50',required=False, help="")
genePythiaParser.add_argument('-q','--queue', metavar='name',required=True, help="PBS queue to be used")
genePythiaParser.add_argument('-maxPub', metavar='N',default='80',required=False, help="max number of publications to print")
genePythiaParser.add_argument('-mail', metavar='user@domain.tld',required=True, help="mail for Entrez query")
genePythiaParser.add_argument('-ab',action="store_true",required=False, help="adding the column abstract")


confHelp="Help you to create the conf file in YAML format"
confParse=subparsers.add_parser('configure',description=confHelp,help=confHelp)

analysisHelp="run the GIDEON analysis with the yaml conf file for gene lists and databases"
analysisParse=subparsers.add_parser('analyse',description=analysisHelp,help=analysisHelp)
analysisParse.add_argument('-o','--outDir', metavar='Path',required=True, help="outPut directory")
analysisParse.add_argument('-Tmax', metavar='N',default='90',required=False, help="maximum percentage of target genes in a metagroup to be selected")
analysisParse.add_argument('-Tmin', metavar='N',default='10',required=False, help="minimum percentage of target genes in a metagroup to be selected")
analysisParse.add_argument('-onlyEnrich',action="store_true",required=False, help="Perform only the enrichment test to evaluate enrichment qualities")
analysisParse.add_argument('-multiplex',action="store_true",required=False, help="Integrate several databases as a layer in a multiplex layer. Then the Louvain algorithm is used to build metagroup")

args=parser.parse_args()

def check():
	from subprocess import Popen, PIPE
	import IPython
	res="### Informations related to Python ###\n"
	res+=IPython.sys_info().replace("': '",":\t").replace("': u'",":\t").replace("","").replace("',","").replace("{'","").replace(" '","").replace("'}","")
	res+="\n######\n"
	res+="\n### Informations related to R ###\n"
	process = Popen("Rscript -e 'sessionInfo()'", stdout=PIPE, stderr=PIPE,shell=True)
	stdout, stderr = process.communicate()
	res+=stdout
	res+="\n######\n"
	return(res)

import sys

if args.analysis=="version":
	print(version)
	sys.exit()
elif args.analysis=="check":
	print(check())
	sys.exit()
elif args.analysis=="genePythia":
	import datetime
	nowTime=datetime.datetime.now()
	from jupype import *
	rootedDir=RootDir(args.outDir,pbs=True)
	rootedDir.logs.writeArgs(args)
	with open(rootedDir.logs.path+"/check.txt",'a')as checkFile:
		checkFile.write("system check printed at time: "+str(nowTime)+"\n")
		checkFile.write(check()+"\n\n\n")
	if args.ab:
		ab=" -ab"
	else:
		ab=""
	batch_count=0
	job_num=0
	cmdList=[]
	batch_lim=int(args.batch_number)
	with open(args.geneList) as geneListFile:
		for line in geneListFile.readlines():
			line=line.strip()
			batch_count+=1
			if "\t" in line:
				gene,alias=line.split("\t")
			else:
				gene=line
				alias="None"
			cmdList.append("mkdir -p "+rootedDir.results+"/"+gene+"\ncd "+rootedDir.results+"/"+gene+"\n\n"+scripts+"/genePythia/genePythia.py -gene "+gene+" -alias "+alias+" -term '"+args.term+"' -maxPub "+args.maxPub+" -mail "+args.mail+ab)
			if batch_count>batch_lim:
				batch_count=0
				job_num+=1
				submitQsubWithPBS(createPBS(cmdList,"genePythia_"+str(job_num),queue=args.queue,workdir=rootedDir.results))
				cmdList=[]
	if cmdList!=[]:
		job_num+=1
		submitQsubWithPBS(createPBS(cmdList,"genePythia_"+str(job_num),queue=args.queue,workdir=rootedDir.results))
	saveRoot(rootedDir)
	sys.exit(0)




