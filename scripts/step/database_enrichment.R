#!/usr/bin/env Rscript

suppressMessages(library(argparser, quietly=TRUE))

# Create a parser
p <- arg_parser("R program designed to perform an enrichment analysis")

# Add command line arguments
p <- add_argument(p, "i", help="input gene list return-separated with no header.", type="character")
p <- add_argument(p, "u", help="universe gene list return-separated with no header.", type="character")
p <- add_argument(p, "a", help="annotation table. This table is tab delimited with a header. The column genes contains the comma separated list of genes in a term", type="character")
p <- add_argument(p, "o", help="output prefix", type="character")
p <- add_argument(p, "p", help="comma-separated list of regex pattern to remove from gene name",default="",type="character")

# Parse the command line arguments
argv <- parse_args(p)

suppressMessages(library(dplyr, quietly=TRUE))
suppressMessages(library(xlsx, quietly=TRUE))

### GO ENRICHMENT ###
geneInput=argv$i
geneUniverse=argv$u
annSource=argv$a
enrichOutPut=argv$o

filtered$

goVector=filtered$id
      goKeep=logical(length = length(goVector))
      names(goKeep)=filtered$id
      for (goID in goVector){
        offspring=GOffspring[goID]
        if (length(intersect(offspring,goVector))>0){
          goKeep[goID]=F
        }else{
          goKeep[goID]=T
        }
      }


if (argv$=!=""){
  for (patternToRemove in strsplit())
  GI=read.table(geneInput,as.is=T)$V1
  GI=unique(gsub("dup[0-9]+_","",GI,perl=T))
  Universe=unique(gsub("dup[0-9]+_","",as.character(read.table(geneUniverse)$V1),perl=T))
}else{
  GI=as.character(read.table(geneInput)$V1)
  Universe=as.character(read.table(geneUniverse)$V1)
}

if (argv$rmkg){
  print("removing kg prefix")
  GI=unique(gsub("kg_","",GI,perl=T))
  Universe=unique(gsub("kg_","",Universe,perl=T))
}

ann=read.delim(file=annSource,sep="\t",header=T,as.is=T)

ann$pval=NA
ann$GI=''
ann$Ngenes=NA
ann$NGI=NA

N=length(Universe) # Total number of genes considered
k=length(GI)
  
for (i in 1:dim(ann)[1]){
  if (argv$rmdup){
    annGenes=unique(gsub("dup[0-9]+_","",unlist(strsplit(as.character(ann[i,'genes']), split=",")),perl=T))
  }else{
    annGenes=unlist(strsplit(as.character(ann[i,'genes']), split=","))
  }
  if (argv$rmkg){
    annGenes=unique(gsub("kg_","",annGenes,perl=T))
  }
  annGenes=intersect(annGenes,Universe)
  overlap=intersect(annGenes,GI)
  x=length(overlap)
  m=length(annGenes)
  n=N-m
  #1-phyper(q-1,m,n,k)   # one-sided fisher equivalent to hypergeometric
  mat=matrix(c(x,k-x,m-x,N-x-(k-x)-(m-x)),nrow=2)
  ann[i,'pval']=fisher.test(mat,alternative = 'greater')$p.value
  ann[i,'GI']=paste(overlap,collapse = ",")
  ann[i,'genes']=paste(annGenes,collapse = ",")
  ann[i,'NGI']=x
  ann[i,'Ngenes']=m
}

# basic filtering
ann=subset(ann,NGI>0)
ann$qval_BH=p.adjust(ann$pval,method='BH')
ann$qval_bonferroni=p.adjust(ann$pval,method='bonferroni')

#kable(subset(ann,pval<0.05))
ann=arrange(ann,pval)
write.table(ann,file=paste(enrichOutPut,"tab",sep='.'),quote = F,row.names = F,sep="\t")
write.xlsx(ann,file=paste(enrichOutPut,'xlsx',sep='.'))