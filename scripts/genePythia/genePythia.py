#!/export/prog/python/bin/python

import argparse

version=2.1
year=2017
author='Julien Fouret'
employer='ViroScan3D'
contact='julien@fouret.me'
Licence=" ".join(["LICENCE V3D-intern",
"The author grants you the permission to use the software under the following conditions:",
"(i) Distribution is limited to ViroScan3D or ProfileXpert, usage via network is considered as a distribution",
"(ii) No modification is allowed without the permissions of the author(s) or ViroScan3D",
"(iii) Commercial use is not allowed without the permission of Viroscan3D",
"(iv) Any use must credit the contribution from the author and ViroScan3D",
"(v) NO WARRANTY is granted",
"(vi) ViroScan3D or the author are not liable for any use of this software"])

desc="\n".join(["genePythia is a software allowing automactic search (via https requests) of articles in pubmed database for co-occurence of a theme designed by a list of terms with a gene (or any aliases if precised)",
"Statistics are reported and in addtion all details about results are fetched and formatted"])
dep="Python>=2.7;Biopython>=1.68"

##parse argument
parser = argparse.ArgumentParser(description=desc,epilog="Version : "+str(version)+" | "+str(year)+" | Author : "+author+" | Employer : "+employer+" for more informations or enquiries please contact "+contact+" | DEPENDENCES : "+dep+" | "+Licence,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-gene', metavar='name',required=True, help="name of the gene")
parser.add_argument('-alias', metavar='alias1,alias2,alias',default='None',required=False, help="commat-separated list of")
parser.add_argument('-term', metavar='syn1,syn2,syn3',required=True, help="commat-separated list of terms to search (synonyms). Please add\" if a term is composed with spaces" )
parser.add_argument('-maxPub', metavar='N',default='80',required=False, help="max number of publications to print")
parser.add_argument('-mail', metavar='user@domain.tld',default=contact,required=False, help="mail for Entrez query")
parser.add_argument('-ab',action="store_true",required=False, help="adding the column abstract")
parser.add_argument('-v',action="store_true",required=False, help="verbose")
args=parser.parse_args()
from Bio import Entrez,Medline
import codecs
import sys
import httplib

def search(query,defType="xml"):
	global args
	Entrez.email = args.mail
	handle = Entrez.esearch(db='pubmed',
		sort='relevance',
		retmode=defType,
		retmax='9999999',
		term=query)
	results = Entrez.read(handle)
	return results

def getCitations(pmid):
	try: 
		results = Entrez.read(Entrez.elink(dbfrom="pubmed", db="pmc",LinkName="pubmed_pmc_refs", id=pmid))
		pmc_ids = [link["Id"] for link in results[0]["LinkSetDb"][0]["Link"]]
		return str(len(pmc_ids))
	except:
		return '0'

def getNumberPub(query):
	try:
        	results = search(query)
        	id_list = results['IdList']
		return len(id_list)
	except: return 0

def getFormattedPubList(query):
	global args
	if args.v:
		print("INFO:PUBLIST : starting to extract details about publications")
	try:
		if args.v:
			print("INFO:PUBLIST : searching the query")
		results = search(query)
		if args.v:
                        print("INFO:PUBLIST : extracting the all pmids")
		id_list = results['IdList']
		if args.v:
			print("INFO:PUBLIST : Fetching details of top "+args.maxPub+" from "+str(len(id_list))+" results ranked by relevance")
		ids = ','.join(id_list)
		Entrez.email = args.mail
		handle = Entrez.efetch(db='pubmed',
			retmode='text',
			rettype="medline",
			retmax=args.maxPub,
			id=ids)
		print("INFO:PUBLIST : Parsing fetched details")
		papers = Medline.parse(handle)
		formattedList=list()
		notEmpty=True
		iterTest=1
		print("INFO:PUBLIST : Initiating the iteration to format the details")
		while notEmpty==True:
			errorHTTP=True
			while errorHTTP==True:
				errorHTTP=False
				try:
					paper=next(papers)
					pos=str(iterTest)
					if args.v:
						print('###==> article number: '+pos+"\n")
					iterTest+=1
					title=paper.get("TI", "?")
					pmid=paper.get("PMID", "?")
					authors=';'.join(paper.get("FAU", "?"))
					date=paper.get("DP", "?")
					journal=paper.get("JT", "?")
					cited=getCitations(pmid)
					if args.ab:
						abstract=paper.get("AB", "?")
						formattedList.append([pos,pmid,title,authors,journal,date,cited,abstract])
					else:
						formattedList.append([pos,pmid,title,authors,journal,date,cited])
				except httplib.IncompleteRead:
					errorHTTP=True
					print('ERROR:PUBLIST : httplib.incompletedRead'+"\n")
				except :
					notEmpty=False
					if args.v:
						print('INFO:PUBLIST : Final error because the list of papers is finished'+"\n")
		if args.v:
			print("INFO:PUBLIST : Ending the iteration to format the details")
		return formattedList
	except Exception as e:
		print("ERROR:PUBLIST : <== Unexpected ==> "+"\n"+str(e))
		return None
##############
#		for i, paper in enumerate(papers):
#			#formattedList.append()
#			####Medline.read from here ???
###			title=paper.get("TI", "?")
#			pmid=paper.get("PMID", "?")
#			authors=';'.join(paper.get("FAU", "?"))
#			date=paper.get("DP", "?")
#			journal=paper.get("JT", "?")
#			pos=str(i+1)
#			cited=getCitations(pmid)
#			formattedList.append([pos,pmid,title,authors,journal,date,cited])
#		return formattedList
#	except:
#		return None


#pos pubmedId title journal year cited

if args.alias!='None':
	alias=True
	aliasList=args.alias.split(',')

else:
	alias=False

termList=args.term.split(',')

# start only with gene
if args.ab:
	header=['pos','pmid','title','authors','journal','date','cited','abstract']
else:
	header=['pos','pmid','title','authors','journal','date','cited']
queryGene=args.gene
queryTerm=' OR '.join(termList)
queryBoth='('+queryGene+') AND ('+queryTerm+')'

numGene=getNumberPub(queryGene)
geneNumBoth=getNumberPub(queryBoth)
geneResults=getFormattedPubList(queryBoth)

with open('gene_article.tab','w') as outFile:
	outFile.write("\t".join(header))
	if geneResults!=None:
		for lineList in geneResults:
			outFile.write("\n"+"\t".join(lineList).encode('utf-8','ignore'))

#continue with alias
if alias:
	queryAlias=' OR '.join(aliasList)
	queryBoth='('+queryGene+') AND ('+queryTerm+')'
	aliasNumBoth=getNumberPub(queryBoth)
	numAlias=getNumberPub(queryAlias)
	aliasResults=getFormattedPubList(queryBoth)

	with open('alias_article.tab','w') as outFile:
		outFile.write("\t".join(header))
		if aliasResults!=None:
			for lineList in aliasResults:
				outFile.write("\n"+"\t".join(lineList).encode('utf-8','ignore'))
	with open('biblio.tab','w') as outFile:
		outFile.write("\t".join(['terms','alias','termGene','termAlias','onlyGene','onlyAlias']))
		outFile.write("\n"+"\t".join([args.term,','.join(aliasList),str(geneNumBoth),str(aliasNumBoth),str(numGene),str(numAlias)]))
else : 
	with open('biblio.tab','w') as outFile:
		outFile.write("\t".join(['terms','alias','termGene','termAlias','onlyGene','onlyAlias']))
		outFile.write("\n"+"\t".join([args.term,"NA",str(geneNumBoth),"NA",str(numGene),"NA"]))

sys.exit()
