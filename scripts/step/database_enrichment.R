#!/usr/bin/env Rscript

suppressMessages(library(argparser, quietly=TRUE))

# Create a parser
p <- arg_parser("R program designed to perform an enrichment analysis with topology-based annotation. Over-representation analysis is performed according to the parent-child-intersection method. Multi-testing is done with BH and Bonferoni methods on the subset of term with a minimal p-value below 1e-7 by default. P-values are based on combination analysis as described in the parent-child paper using GNU multiple precision arithmetic libraries and their R implementation")

# Add command line arguments
p <- add_argument(p, "-i", help="input gene list return-separated with no header.", type="character")
p <- add_argument(p, "-u", help="universe gene list return-separated with no header.", type="character")
p <- add_argument(p, "-a", help="annotation table. This table is tab delimited with a header (term,genes,type,name,parents). The column genes contains the comma separated list of genes in a term", type="character")
p <- add_argument(p, "-o", help="output prefix", type="character")
p <- add_argument(p, "-p", help="minimum value for pmin (hihgest possible p-value)", type="character",default = "1e-7")
p <- add_argument(p, "-r", help="regex pattern to remove from gene name",default="-.*",type="character")
p <- add_argument(p, "-t", help="comma-separated list of type to keep (keep all by default)",default="",type="character")

# Parse the command line arguments
argv <- parse_args(p)
suppressMessages(library(tidyr, quietly=TRUE))
suppressMessages(library(stringr, quietly=TRUE))
suppressMessages(library(dplyr, quietly=TRUE))
suppressMessages(library(xlsx, quietly=TRUE))
suppressMessages(library(Rmpfr, quietly=TRUE))

### ENRICHMENT ###
geneInput=argv$i
geneUniverse=argv$u
annSource=argv$a
enrichOutPut=argv$o
lim_pmin=as.numeric(argv$p)
regex=argv$r
type_keep=argv$t

GI=as.character(read.table(geneInput)$V1)
Universe=as.character(read.table(geneUniverse)$V1)

Numbers_GI=c(length(GI))
Numbers_Universe=c(length(Universe))

if (regex!=""){
  GI=unique(gsub(regex,"",GI,perl=T))
  Universe=unique(gsub(regex,"",Universe,perl=T))
}

Numbers_GI=c(Numbers_GI,length(GI))
Numbers_Universe=c(Numbers_Universe,length(Universe))
ann=read.delim(file=annSource,sep="\t",header=T,as.is=T)
if (type_keep!=""){
  list_ann=list()
  type_keeps=str_split(type_keep,",")[[1]]
  
  for (i in seq(1,length(type_keeps))){
    type_keep=type_keeps[i]
    list_ann[[i]]=data.frame(ann %>% group_by(term) %>% filter( type_keep %in% strsplit(type,";")[[1]] ))
  }
  ann = bind_rows(list_ann) %>% distinct()
}

##

allgenes=unique(gsub(regex,"",unique(str_split(paste(ann$genes,collapse = ";"),'[;,]{1}')[[1]]),perl=T))

Universe=intersect(Universe,allgenes)
GI=intersect(GI,allgenes)
##

Numbers_GI=c(Numbers_GI,length(GI))
Numbers_Universe=c(Numbers_Universe,length(Universe))

GI=intersect(GI,Universe)

Numbers_GI=c(Numbers_GI,length(GI))
Numbers_Universe=c(Numbers_Universe,length(Universe))

data_numbers=data.frame(GI=Numbers_GI,Universe=Numbers_Universe)
row.names(data_numbers)=c("Input","Unique Input after gsub","Mapped in annotation space","GI in Universe")

test=knitr::kable(data_numbers,format = "rst")

write(test,file = paste(enrichOutPut,"info",sep='.'))

N=length(Universe) # Total number of genes considered
k=length(GI)
  
filter_glist=function(comma_genes,glist){
  genes=unique(gsub(regex,"",strsplit(comma_genes,'[;,]{1}',perl=T)[[1]],perl=T))
  genes_filtered=intersect(glist,genes)
  return(paste(genes_filtered,collapse = ","))
}


filter_univers=function(comma_genes){
  return(filter_glist(comma_genes,Universe))
}

filter_gi=function(comma_genes){
  return(filter_glist(comma_genes,GI))
}

ann=ann %>% group_by(term) %>% mutate(genes=filter_univers(genes))

ann=ann %>% group_by(term) %>% mutate(GI=filter_gi(genes))

ann=data.frame(ann)
rownames(ann)=ann$term

ann[is.na(ann)]=""

get_count_from_ids=function(terms,name){
  if (terms==""){
    if (name=="genes"){
      return(N)
    }else if(name=="GI"){
      return(k)
    }
  }else{
    vec=strsplit(terms,"[;,]{1}")[[1]]
    elements=strsplit(paste(ann[vec,name],collapse = ","),',')[[1]]
    t=table(elements)
    return(sum(as.integer(t)==length(vec))) # Count intersection parent-child intersection method
  }
}

GI_count_from_ids=function(terms){
  get_count_from_ids(terms,"GI")
}

genes_count_from_ids=function(terms){
  get_count_from_ids(terms,"genes")
}

ann=ann %>% group_by(term) %>% mutate(
  NGI=GI_count_from_ids(term),
  Ngenes=genes_count_from_ids(term),
  parent_NGI=GI_count_from_ids(parents),
  parent_Ngenes=genes_count_from_ids(parents)
)
  
get_pval=function(Ngenes,NGI,parent_Ngenes,parent_NGI){
  numerator_k=function(k){
    return(chooseMpfr(Ngenes,k)*chooseMpfr(parent_Ngenes-Ngenes,parent_NGI-k))
  }
  denominator=chooseMpfr(parent_Ngenes,parent_NGI)
  result = tryCatch({
    as.numeric(signif(sum(numerator_k(seq(NGI,min(Ngenes,parent_NGI))))/denominator,digits=6))
  }, error = function(e) {
    NA
  })
  return(result)
}

get_pmin=function(Ngenes,parent_Ngenes,parent_NGI){
  pval_k=function(k){
    return((chooseMpfr(Ngenes,k)*chooseMpfr(parent_Ngenes-Ngenes,parent_NGI-k))/chooseMpfr(parent_Ngenes,parent_NGI))
  }
  return(as.numeric(signif(pval_k(min(Ngenes,parent_NGI)),digits=6)))
}

ann=ann %>% group_by(term) %>% mutate(pval=get_pval(Ngenes,NGI,parent_Ngenes,parent_NGI),pmin=get_pmin(Ngenes,parent_Ngenes,parent_NGI))

ann=data.frame(ann)
# basic filtering
ann$qval_BH=NA
ann$qval_bonferroni=NA

ann[ann$pmin<(10^(-7)),"qval_BH"]=p.adjust(ann[ann$pmin<(10^(-7)),"pval"],method='BH')

ann[ann$pmin<(10^(-7)),"qval_bonferroni"]=p.adjust(ann[ann$pmin<(10^(-7)),"pval"],method='bonferroni')

#kable(subset(ann,pval<0.05))
ann=arrange(ann,qval_BH)

write.table(ann,file=paste(enrichOutPut,"tab",sep='.'),quote = F,row.names = F,sep="\t")

ann = ann %>% group_by(term) %>% mutate(
  string_db_image=paste('https://string-db.org/api/image/network?identifiers=',gsub(",","%0d",GI),"&species=9606",sep=""),
  string_db_site=paste('https://string-db.org/newstring_cgi/show_network_section.pl?identifiers=',gsub(",","%0d",GI),"&species=9606",sep=""))

write.xlsx(data.frame(ann),file=paste(enrichOutPut,'xlsx',sep='.'),row.names = F)






