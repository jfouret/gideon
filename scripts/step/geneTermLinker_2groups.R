### Function ###

loadData=function(conf,exclude=c(),n=4){
  suppressMessages(require(topGO))
  suppressMessages(require(dplyr))
  GOffspring <- c(as.list(GOBPOFFSPRING),as.list(GOCCOFFSPRING),as.list(GOMFOFFSPRING))
  confdata=read.table(conf,header=T,sep="\t",as.is = T)
  newHeader=c('group','space','id','desc','allGenes','genes','N','pval')
  for (i in 1:dim(confdata)[1]){
    CfileName=confdata[i,'file']
    Cgroup=confdata[i,'group']
    Ctitle=confdata[i,'title']
    Cannspace=confdata[i,'annspace']
    Cid=confdata[i,'id']
    Cdesc=confdata[i,'desc']
    Cgenes=confdata[i,'genes']
    CallGenes=confdata[i,'allGenes']
    Cprob=confdata[i,'prob']
    CprobMax=confdata[i,'probMax']
    CN=confdata[i,'N']
    tmpdata=read.table(CfileName,header=T,sep="\t",as.is=T,quote = '"')
    if (is.na(Cannspace)){
      annspaceName=Ctitle
    }else{
      annspaceName=paste(Ctitle,tmpdata[,Cannspace],sep='_')
    }
    tmpdata$space=annspaceName
    tmpdata$group=Cgroup
    sourceHeader=c('group','space',Cid,Cdesc,CallGenes,Cgenes,CN,Cprob)
    tmpdata=tmpdata[,sourceHeader]
    colnames(tmpdata)=newHeader
    tmpdata=subset(tmpdata,!(space %in% exclude))
    data_threshold=as.data.frame(tmpdata %>% group_by(space) %>% summarise(threshold=mean(N)+n*sd(N)))
    row.names(data_threshold)=data_threshold$space
    filtered=filter(tmpdata,(pval<CprobMax)&(N<data_threshold[space,'threshold'])&(genes!=''))
    if (Ctitle=='go'|Ctitle=='Go'|Ctitle=='GO'){
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
      filtered=filtered[goKeep,]
    }
    if (i==1){
      data=filtered
    }else{
      data=rbind(data,filtered)
    }
  }
  return(data)
}
integrateGroupData=function(data,genesA,genesB,genesG){
  require(dplyr)
  dA=filter(data,group=='A')
  dB=filter(data,group=='B')
  onlyA=anti_join(dA,dB,by='id')
  onlyB=anti_join(dB,dA,by='id')
  both=inner_join(dA,dB,by='id',suffix = c(".A", ".B"))
  ###onlyA$G=0 ###
  for (i in 1:dim(onlyA)[1]){
    genes=unlist(strsplit(as.character(onlyA[i,'allGenes']),','))
    genes.A=intersect(genesA,genes)
    genes.B=intersect(genesB,genes)
    ###onlyA[i,'G']=(length(genes.A)-length(genes.B))/length(genes)
    onlyA[i,'genes']=paste(unique(c(genes.A,genes.B)),collapse = ',')
    
  }
  ###onlyB$G=0
  for (i in 1:dim(onlyB)[1]){
    genes=unlist(strsplit(as.character(onlyB[i,'allGenes']),','))
    genes.A=intersect(genesA,genes)
    genes.B=intersect(genesB,genes)
    ###onlyB[i,'G']=(length(genes.A)-length(genes.B))/length(genes)
    onlyB[i,'genes']=paste(unique(c(genes.A,genes.B)),collapse = ',')
  }
  if (dim(both)[1]>0){
    ##both$G=0
    for (i in 1:dim(both)[1]){
      genes.A=unlist(strsplit(both[i,'genes.A'],','))
      genes.B=unlist(strsplit(both[i,'genes.B'],','))
      genes=unique(c(genes.A,genes.B))
      ##both[i,'G']=(length(genes.A)-length(genes.B))/length(genes)
      if (both[i,"pval.A"]>both[i,"pval.B"]){
        both[i,"group"]="B"
      }else{
        both[i,"group"]="A"
      }
      both[i,'genes']=paste(genes,collapse = ',')
    }
    both=transmute(both,group=group,space=space.A,id=id,desc=desc.A,genes=genes,allGenes=allGenes.A,N=pmin(N.A,N.B),pval=pmin(pval.A,pval.B)) 
    ndata=rbind(both,onlyA,onlyB)
  }else{
    ndata=rbind(onlyA,onlyB)
  }
  ndata$G=0
  for (i in 1:dim(ndata)[1]){
    Gvector=genesG[unlist(strsplit(ndata[i,'genes'],','))]
    ndata[i,'G']=sum(Gvector)/length(Gvector)
  }
  return(ndata)
}
getIncidencePval=function(filtered){
  ### step 2 create an incidence matrix
  genes=unique(strsplit(paste(filtered$genes,collapse=','),',')[[1]])
  Nterm=dim(filtered)[1]
  Ngenes=length(genes)
  geneDict=1:Ngenes
  names(geneDict)=genes
  incid_mat_pval=matrix(0,nrow = Nterm,ncol = Ngenes+1)
  for (i in 1:dim(filtered)[1]){
    termGenes=strsplit(filtered[i,'genes'],',')[[1]]
    for (geneInTerm in termGenes){
      incid_mat_pval[i,geneDict[geneInTerm]]=1
    }
    incid_mat_pval[i,Ngenes+1]=filtered[i,'pval']*filtered[i,'N']
  }
  
  rownames(incid_mat_pval)=filtered$id
  colnames(incid_mat_pval)=c(genes,'pvalxN')
  return(incid_mat_pval)
}
getIncidences=function(data,genes){
  ### step 2 create an incidence matrix
  Nterm=dim(data)[1]
  Ngenes=length(genes)
  geneDict=1:Ngenes
  names(geneDict)=genes
  incid_mat=matrix(0,nrow = Nterm,ncol = Ngenes)
  for (i in 1:dim(data)[1]){
    termGenes=strsplit(data[i,'genes'],',')[[1]]
    for (geneInTerm in termGenes){
      incid_mat[i,geneDict[geneInTerm]]=1
    }
  }
  
  rownames(incid_mat)=data$id
  colnames(incid_mat)=genes
  return(incid_mat)
}
getCluster=function(incid_mat_pval){
  cosineDist <- function(x){
    1 - x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))
  }
  termDist=as.dist(cosineDist(incid_mat_pval))
  cl=hclust(termDist,method='ward.D2')
  return(cl)
}
getMetagroups=function(cl,hlim,filtered,reduced=F){
  metagroups=cutree(cl,h=hlim*max(cl$height))
  filtered$metagroup=metagroups
  if (reduced==F){
    row.names(filtered)=filtered$id
    return(filtered)
  }else{
    filtered_sort=filtered[with(filtered, order(pval)), ]
    init=0
    for (meta in 1:max(metagroups)){
      sub_meta=subset(filtered_sort,metagroup==meta)
      genes_meta=unique(strsplit(paste(sub_meta$genes,collapse=','),',')[[1]])
      Nmeta=length(genes_meta)
      genes_cover=c()
      for (i in 1:dim(sub_meta)[1]){
        genes_cover=unique(c(genes_cover,strsplit(sub_meta[i,'genes'],',')[[1]]))
        if (length(genes_cover)==Nmeta){
          break()
        }
      }
      noredundance_meta=sub_meta[1:i,]
      if (init==0){
        init=1
        data_meta=noredundance_meta
      }else{
        data_meta=rbind(data_meta,noredundance_meta)
      }
    }
    row.names(data_meta)=data_meta$id
    return(data_meta)
  }
}
getMetaIncidTerm=function(data_meta,weight=F){
  require(stringr)
  terms=as.character(data_meta$id)
  metas=unique(data_meta$metagroup)
  incid_meta_terms=matrix(0,nrow=length(metas),ncol=length(terms))
  rownames(incid_meta_terms)=as.character(metas)
  colnames(incid_meta_terms)=terms
  for (term in terms){
    if (weight==T){
      incid_meta_terms[data_meta[term,'metagroup'],term]=str_count(data_meta[term,'genes'],',')+1
    }else{
      incid_meta_terms[data_meta[term,'metagroup'],term]=1
    }
  }
  return(incid_meta_terms)
}
getMetaIncidGenes=function(data_meta){
  terms=as.character(data_meta$id)
  genes=unique(strsplit(paste(data_meta$genes,collapse=','),',')[[1]])
  metas=unique(data_meta$metagroup)
  incid_meta_genes=matrix(0,nrow=length(metas),ncol=length(genes))
  rownames(incid_meta_genes)=as.character(metas)
  colnames(incid_meta_genes)=genes
  for (term in terms){
    metaTerm=data_meta[term,'metagroup']
    for (geneInTerm in strsplit(data_meta[term,'genes'],',')[[1]]){
      incid_meta_genes[metaTerm,geneInTerm]=1
    }
  }
  return(incid_meta_genes)
}
plotMetaBiPart=function(net,plot_title,labels='',labels2='',colorsScale=T,labelSize=0.6){
  require('igraph')
  l <- layout_with_fr(net,dim=2,niter = 1000)
  V(net)$color <- c("steel blue", "orange")[V(net)$type+1]
  if (colorsScale==T){
    V(net)[V(net)$type==F]$color<-rainbow(length(V(net)[V(net)$type==F]$color), alpha = 0.3)
  }
  V(net)$shape <- c("square", "circle")[V(net)$type+1]
  if (labels2==''){
    V(net)$label[V(net)$type==F] <- V(net)$name[V(net)$type==F]
  }else{
    V(net)$label[V(net)$type==F] <- labels2[V(net)$name[V(net)$type==F]]
  }
  if (labels==''){
    V(net)$label[V(net)$type==T] <- V(net)$name[V(net)$type==T]
  }else{
    V(net)$label[V(net)$type==T] <- labels[V(net)$name[V(net)$type==T]]
  }
  V(net)$size=8
  V(net)$label.color='black'
  V(net)$label.cex=labelSize
  V(net)$label.font=2
  plot(net,layout=l,main=plot_title)
}
plotGTBiPart=function(net,plot_title,marks,labels='',labels2='',colorsScale=F){
  require('igraph')
  l <- layout_with_dh(net,maxiter = 100)
  V(net)$color <- c("steel blue", "orange")[V(net)$type+1]
  if (colorsScale==T){
    V(net)[V(net)$type==F]$color<-rainbow(length(V(net)[V(net)$type==F]$color), alpha = 0.3)
  }
  V(net)$shape <- c("square", "circle")[V(net)$type+1]
  if (labels2==''){
    V(net)$label[V(net)$type==F] <- V(net)$name[V(net)$type==F]
  }else{
    V(net)$label[V(net)$type==F] <- labels2[V(net)$name[V(net)$type==F]]
  }
  if (labels==''){
    V(net)$label[V(net)$type==T] <- V(net)$name[V(net)$type==T]
  }else{
    V(net)$label[V(net)$type==T] <- labels[V(net)$name[V(net)$type==T]]
  }
  V(net)$size=8
  V(net)$label.cex=.6
  V(net)$label.font=2
  markGroups=list()
  for (i in 1:length(marks)){
    markGroups[[i]]<-strsplit(marks[i],',')[[1]]
  }
  plot(net,layout=l,main=plot_title,mark.groups=markGroups, mark.col=rainbow(length(marks), alpha = 0.3))
}
getMetaSum=function(data_meta,incid_meta_genes,incid_meta_terms){
  metagroups=unique(data_meta$metagroup)
  matrix_sum=matrix('',nrow=length(metagroups),ncol=3)
  colnames(matrix_sum)=c('genes','allGenes','terms')
  for (meta in metagroups){
    sub_meta=subset(data_meta,metagroup==meta)
    matrix_sum[meta,'terms']=paste(as.character(sub_meta$id),collapse=',')
    matrix_sum[meta,'genes']=paste(unique(strsplit(paste(sub_meta$genes,collapse = ','),',')[[1]]),collapse=',')
    matrix_sum[meta,'allGenes']=paste(unique(strsplit(paste(sub_meta$allGenes,collapse = ','),',')[[1]]),collapse=',')
  }
  summary_meta=as.data.frame(matrix_sum)
  return(summary_meta)
}
plotProjMarked=function(netProj,marks,plot_title){
  require('igraph')
  V(netProj)$color <- 'orange'
  V(netProj)$label.cex=.6
  V(netProj)$label.font=2
  V(netProj)$size=8
  markGroups=list()
  for (i in 1:length(marks)){
    markGroups[[i]]<-strsplit(marks[i],',')[[1]]
  }
  l <- layout_with_fr(netProj,dim=2,niter = 1000)
  plot(netProj,mark.groups=markGroups, mark.col=rainbow(length(marks), alpha = 0.3),mark.border = NA,main=plot_title)
}
analyseNetworks=function(data,pvalMax,n,k,removePrefix=F){
  require('igraph')
  if (removePrefix==T){
    prefix=paste(c('noPrefix_','pvalMax_',pvalMax,'_n_',n,'_k_',k),collapse='')
  }else{
    prefix=paste(c('pvalMax_',pvalMax,'_n_',n,'_k_',k),collapse='')
  }
  filtered=filterTerms(data,n=n,pvalMax=pvalMax,removePrefix=removePrefix)
  incid_mat_pval=getIncidencePval(filtered)
  cl=getCluster(incid_mat_pval)
  data_meta=getMetagroups(cl,k,filtered)
  incid_meta_terms=getMetaIncidTerm(data_meta)
  incid_meta_genes=getMetaIncidGenes(data_meta)
  summary_meta=getMetaSum(data_meta,incid_meta_genes,incid_meta_terms)
  markGenes=as.character(summary_meta$genes)
  markTerms=as.character(summary_meta$terms)
  netMetaTerms=graph_from_incidence_matrix(incid_meta_terms)
  netMetaGenes=graph_from_incidence_matrix(incid_meta_genes)
  labels=as.character(filtered$desc)
  names(labels)=as.character(filtered$id)
  incid_mat=incid_mat_pval[,1:(dim(incid_mat_pval)[2]-1)]
  net_term_genes=graph_from_incidence_matrix(incid_mat)
  net.bg<-bipartite.projection(net_term_genes)
  netGenes <- net.bg$proj2
  netTerms <- net.bg$proj1
  V(netTerms)$label <- labels[V(netTerms)$name]
  pdf(file = paste(prefix,"_metas.pdf",sep=''),width = 20,height = 20)
  par(mfrow=c(2,2))
  plotProjMarked(netTerms,markTerms,'Terms Network with metagroups')
  plotProjMarked(netGenes,markGenes,'Genes Network with metagroups')
  plotMetaBiPart(netMetaTerms,'Terms-metagroups bipartite-network',labels,labelSize=1)
  plotMetaBiPart(netMetaGenes,'Genes-metagroups bipartite-network')
  dev.off()
  pdf(file = paste(prefix,"_GT.pdf",sep=''),width = 15,height = 15)
  par(mfrow=c(1,1))
  plotGTBiPart(net_term_genes,'Genes-terms bipartite-network',markGenes,labels2=labels)
  dev.off()
}
plotCL=function(cl,h){
  source("/export/work/batnipah/juTools/funMining/geneSetEnrichment/A2RPlot.R")
  op = par(bg = "white")
  nb_group=max(cutree(cl,h=h*max(cl$height)))
  cols = rainbow(nb_group)
  return(A2Rplot(cl, k = nb_group, boxes = FALSE, col.up = "black", col.down = cols,show.labels=F))
}

combine.brown=function(data_matrix,pval_list){
  suppressMessages(require("EmpiricalBrownsMethod"))
  if (length(pval_list)==1){
    return(pval_list)
  }else{
    return(empiricalBrownsMethod(data_matrix,pval_list))
  }
}

breaks_meta=function(data_meta,net_terms){
  uniMetas=unique(data_meta$metagroup)
  for (i in uniMetas){
    metaTerms=subset(data_meta,metagroup==i)$id
    subnets=decompose.graph(induced_subgraph(net_terms,metaTerms))
    if (length(subnets)>1){
      for (j in 1:length(subnets)){
        for (term in V(subnets[[j]])$name){
          data_meta[term,"metagroup"]=100000*i+j
        }
      }
    }
  }
  uniMetas=unique(data_meta$metagroup)
  dictMetas=as.character(1:length(uniMetas))
  names(dictMetas)=uniMetas
  data_meta$metagroup=as.integer(dictMetas[as.character(data_meta$metagroup)])
  return(data_meta)
}

analyseNetworks_2groups=function(conf,n,genesAFile,genesBFile,genesAUniversFile,genesBUniversFile,prefixPatternList,printGraph=F,subAnalyse=T,Glim=-0.8,topPer=25,hlim=0.5,pvalLim=0.01){
  suppressMessages(require('igraph'))
  suppressMessages(require('ggplot2'))
  suppressMessages(require('grid'))
  suppressMessages(require('gridExtra'))
  genesA=read.table(genesAFile,header=F,as.is=T)$V1
  genesB=read.table(genesBFile,header=F,as.is=T)$V1
  genesAUnivers=read.table(genesAUniversFile,header=F,as.is=T)$V1
  genesBUnivers=read.table(genesBUniversFile,header=F,as.is=T)$V1
  for (prefixPattern in prefixPatternList){
    genesA=unique(gsub(prefixPattern,"",genesA,perl=T))
    genesB=unique(gsub(prefixPattern,"",genesB,perl=T))
    genesAUnivers=unique(gsub(prefixPattern,"",genesAUnivers,perl=T))
    genesBUnivers=unique(gsub(prefixPattern,"",genesBUnivers,perl=T))
  }
  
  onlyA=setdiff(genesA,genesB)
  onlyB=setdiff(genesB,genesA)
  allGenes=unique(c(genesA,genesB))
  genesG=rep.int(0, length(allGenes))
  names(genesG)=allGenes
  genesG[onlyA]=1
  genesG[onlyB]=-1
  ##genesG[genesG==0]
  filtered=loadData(conf,exclude=c("go_cellular_component","go_molecular_function"))
  print("Data Loaded from config file")
  dataSum=as.data.frame(filtered %>% group_by(group,space) %>% summarise(count=n()))
  print(dataSum)
  filtered=integrateGroupData(filtered,genesA,genesB,genesG)
  term2group=filtered$group
  names(term2group)=filtered$id
  print("Both groups integrated")
  incid_mat_pval=getIncidencePval(filtered)
  print("Incidence term-gene matrix with genes presence absence as 1 and 0 and pval of terms")
  incid_mat=incid_mat_pval[,1:(dim(incid_mat_pval)[2]-1)]
  net_term_genes=graph_from_incidence_matrix(incid_mat)
  print("gene-term network infered")
  tmp_net.bg<-bipartite.projection(net_term_genes)
  tmp_netTerms <- tmp_net.bg$proj1
  independent_subnets=decompose.graph(tmp_netTerms)
  data_meta=NULL
  blank <- grid.rect(gp=gpar(col="white",alpha=0))
  fileName=paste('graphClust_setps.pdf',sep='')
  lm=matrix(2,ncol=8,nrow=2)
  lm[1,3:8]=3
  lm[1,1]=4
  lm[1,7]=5
  lm[1,8]=6
  lm[1,2]=1
  for (subNum in 1:length(independent_subnets)){
    subPath=paste('subnet',subNum,sep='')
    dir.create(subPath,showWarnings = F)
    setwd(subPath)
    if (length(V(independent_subnets[[subNum]]))>1){
      cl=getCluster(incid_mat[V(independent_subnets[[subNum]])$name,])
      print(paste(subPath,"Term clustered using cosine distance and Ward.D2 algorithm"))
      pdf(file=fileName)
      for (i in seq(0.9,0.1,-0.1)){
        if (i>length(V(independent_subnets[[subNum]]))){break}
        plotCL(cl,i)
        cutted=cutree(cl,h=i*max(cl$height))
        sumDecompose=data.frame(mg=1:i)
        sumDecompose$sub=0
        for (j in 1:i){
          mgterms=names(cutted)[cutted==j]
          subnet=induced_subgraph(independent_subnets[[subNum]],mgterms)
          sumDecompose[j,"sub"]=length(decompose.graph(subnet))
        } 
        sumDecompose=subset(sumDecompose,sub!=1)
        if(dim(sumDecompose)[1]==0){
          grob1=blank
        }else{
          grob1=tableGrob(sumDecompose,rows=NULL)
        }
        grid.arrange(grob1,blank,blank,blank,textGrob(paste("h=",i,sep=""), gp=gpar(fontsize=15,font=8)),blank,layout_matrix =lm,newpage = F)
      }
      dev.off()
      print(paste(subPath,"Dendogram graph by number of cluster considered printed"))
      #k <- readline(prompt=paste(subPath,"Enter the number of desired groups: "))
      #k <- as.integer(k)
      prefix=paste(c('h',hlim,'_'),collapse='')
      if (is.null(data_meta)){
        data_meta=getMetagroups(cl,hlim,subset(filtered,id %in% V(independent_subnets[[subNum]])$name))
      }else{
        tmp_data_meta=getMetagroups(cl,hlim,subset(filtered,id %in% V(independent_subnets[[subNum]])$name))
        tmp_data_meta$metagroup=tmp_data_meta$metagroup+max(data_meta$metagroup)
        data_meta=rbind(data_meta,tmp_data_meta)
      }
    }else{
      if (is.null(data_meta)){
        data_meta=subset(filtered,id %in% V(independent_subnets[[subNum]])$name)
        data_meta$metagroup=1
      }else{
        tmp_data_meta=subset(filtered,id %in% V(independent_subnets[[subNum]])$name)
        tmp_data_meta$metagroup=max(data_meta$metagroup)+1
        data_meta=rbind(data_meta,tmp_data_meta)
      }
    }
    setwd("..")
  }
  row.names(data_meta)=data_meta$id
  data_meta=breaks_meta(data_meta,tmp_netTerms)
  k=max(data_meta$metagroup)
  
  info=c(
    paste("A group (G=1)  from:",genesAFile),
    paste("B group (G=-1) from:",genesBFile),
    paste("A group univers (G=1)  from:",genesAUniversFile),
    paste("B group univers (G=-1) from:",genesBUniversFile),
    paste("Enriched term details from the config file:",conf,"\n\n"),
    knitr::kable(dataSum),collapse="\n",
    paste("prefix to remove ? ",paste(prefixPatternList,collapse=" | ")),
    paste("Print network graph ? ",printGraph),
    paste("Height to cut the clustering dendogram : ",hlim),
    paste("Total number of cluster:",k),
    paste("Top ranked (betweenness centrality) percentage : ",topPer),
    paste("n threshold for filtering too big term : ",n),
    paste("Analyse sub-networks ? ",subAnalyse),
    paste("G minimum to analyse metagroups : ",Glim)
  )
  writeLines(info,"info.txt")
  
  
  write.table(data_meta,"term_meta.tab",quote = F,sep = "\t",row.names = F)
  print("Metagroups extracted from cluster with a reduced number of term")
  
  #data_meta_all=getMetagroups(cl,hlim,filtered,reduced=F)
  #write.table(data_meta,"term_meta_all.tab",quote = F,sep = "\t",row.names = F)
  print("Metagroups extracted from cluster ")
  # get G for all metagroups 
  suppressMessages(require(stringr))
  incid_meta_terms=getMetaIncidTerm(data_meta)
  print("Get incidence metagroup-term")
  incid_meta_genes=getMetaIncidGenes(data_meta)
  print("Get incidence metagroup-gene")
  summary_meta=getMetaSum(data_meta,incid_meta_genes,incid_meta_terms)
  summary_meta$G=0
  termMetaAnnot=c()
  print("Summerize metagroup G")
  suppressMessages(require(survcomp))
  data_meta$group=term2group[data_meta$id]
  summary_meta$metagroup=1:k
  incid_matA=getIncidences(subset(data_meta,group=="A"),genesAUnivers)
  incid_matB=getIncidences(subset(data_meta,group=="B"),genesBUnivers)
  print("Gene-terms incidence matrix computed for A and B separately (different universe of genes considered)")
  
  for (i in 1:dim(summary_meta)[1]){
    ### Repeat test
    tmpGeneList=unlist(strsplit(as.character(summary_meta[i,'genes']),','))
    tmpAllGeneList=unlist(strsplit(as.character(summary_meta[i,'allGenes']),','))
    tmpTermList=unlist(strsplit(as.character(summary_meta[i,'terms']),','))
    tmp_annot=rep.int(i,length(tmpTermList))
    names(tmp_annot)=tmpTermList
    termMetaAnnot=c(termMetaAnnot,tmp_annot)
    summary_meta[i,'nGenes']=length(tmpGeneList)
    summary_meta[i,'nTerms']=length(tmpTermList)
    summary_meta[i,'N']=length(tmpAllGeneList)
    summary_meta[i,'G']=sum(genesG[tmpGeneList])/length(tmpGeneList)
    #N=19000
    #K=length(unique(c(genesA,genesB)))
    #m=summary_meta[i,'N']
    #x=summary_meta[i,'nGenes']=length(tmpGeneList)
    #mat=matrix(c(x,K-x,m-x,N-x-(K-x)-(m-x)),nrow=2)
    #summary_meta[i,'pval']=fisher.test(mat,alternative = 'greater')$p.value
    calculAandB=T
    if (length(subset(data_meta,(metagroup==i)&(group=="A"))$pval)==0){
      summary_meta[i,'pval']=combine.brown(incid_matB[subset(data_meta,(metagroup==i)&(group=="B"))$id,],subset(data_meta,(metagroup==i)&(group=="B"))$pval)
      calculAandB=F
    }
    if(length(subset(data_meta,(metagroup==i)&(group=="B"))$pval)==0){
      summary_meta[i,'pval']=combine.brown(incid_matA[subset(data_meta,(metagroup==i)&(group=="A"))$id,],subset(data_meta,(metagroup==i)&(group=="A"))$pval)
      calculAandB=F  
    }
    if (calculAandB) {
      pvalA=combine.brown(incid_matA[subset(data_meta,(metagroup==i)&(group=="A"))$id,],subset(data_meta,(metagroup==i)&(group=="A"))$pval)
      pvalB=combine.brown(incid_matB[subset(data_meta,(metagroup==i)&(group=="B"))$id,],subset(data_meta,(metagroup==i)&(group=="B"))$pval)
      summary_meta[i,'pval']=combine.test(c(pvalA,pvalB), method = "fisher")  #sumlog(c(pvalA,pvalB))
    }
  }
  
  print("p-value and G computed for all metagroups")
  print("p-value was computed for each group A and B per metagroups with an empirical Brown method considering covariance between terms")
  print("Final p-value was computed by combining A and B p-value using Fisher method")
  
  graph_meta=ggplot(summary_meta)+
    geom_bar(aes(metagroup,nGenes,fill=pval),stat="identity")+
    geom_bar(aes(metagroup,-nTerms,fill=pval),stat="identity")+
    geom_hline(yintercept = 0)+
    geom_text(aes(metagroup,nGenes,label=signif(G,digits=3)),hjust=-0.05)+
    geom_label(aes(metagroup,0,label=signif(nGenes/nTerms,digits=3)))+
    labs(title="Metagroups composition",y="number of terms(-) and genes(+)")+
    scale_fill_gradientn(limits=c(0,1),values = scales::rescale(c(0,0.1,1)),colours = c("#50FF50","#FF5050","#5050FF"))+
    scale_y_continuous(breaks=unique(c(seq(-max(summary_meta$nTerms),max(summary_meta$nGenes),by=20),max(summary_meta$nGenes)))) +
    scale_x_continuous(breaks=1:k)+
    coord_flip()+
    theme_bw()
  ggsave(filename = "metagroups_summary.pdf",graph_meta,height = 1+0.5*max(data_meta$metagroup),width = 8)
  write.table(summary_meta,"summary_meta.tab",quote = F,sep = "\t",row.names = F)
  row.names(filtered)=filtered$id
  for (mgi in 1:dim(summary_meta)[1]){
    filename=paste("metagroup_",mgi,'.txt',sep='')
    target_ids=strsplit(as.character(summary_meta$terms[mgi]),',')[[1]]
    write.table(filtered[target_ids,][order(filtered[target_ids,]$pval),],sep ="\t",row.names = F ,file=filename,quote = F)
  }
  print("Graph for metagroups elaborated")
  ### Prepare named dict for annotations ! ! !
  Gmeta=summary_meta$G
  names(Gmeta)=row.names(summary_meta)
  Gterm=data_meta$G
  names(Gterm)=data_meta$id
  DESCterm=data_meta$desc
  names(DESCterm)=data_meta$id
  print("names prepared for annotation")
  #################
  markGenes=as.character(summary_meta$genes)
  markTerms=as.character(summary_meta$terms)
  print("G values in vectors for network annotation")
  weigthed_incid_meta_terms=getMetaIncidTerm(data_meta,T)
  print("number of genes weigthed incidence matrix between term and metagroup created")
  netMetaTerms=graph_from_incidence_matrix(weigthed_incid_meta_terms,weighted=T)
  print("metagroup-term network infered")
  #Annotation of network metagroups + terms
  V(netMetaTerms)$cat='term'
  V(netMetaTerms)[V(netMetaTerms)$type==F]$cat='meta'
  V(netMetaTerms)[V(netMetaTerms)$cat=='meta']$G=Gmeta[V(netMetaTerms)[V(netMetaTerms)$cat=='meta']]
  V(netMetaTerms)[V(netMetaTerms)$cat=='term']$G=Gterm[V(netMetaTerms)[V(netMetaTerms)$cat=='term']]
  V(netMetaTerms)$desc=V(netMetaTerms)$name
  V(netMetaTerms)[V(netMetaTerms)$cat=='term']$desc=DESCterm[V(netMetaTerms)[V(netMetaTerms)$cat=='term']]
  print("metagroup-term network annotated")
  write_graph(netMetaTerms, 'netMetaTerms.graphml', format ="graphml")
  print("metagroup-term network writed")
  netMetaGenes=graph_from_incidence_matrix(incid_meta_genes)
  print("metagroup-genes network infered")
  #Annotation of network metagroups + genes
  V(netMetaGenes)$cat='gene'
  V(netMetaGenes)[V(netMetaGenes)$type==F]$cat='meta'
  V(netMetaGenes)[V(netMetaGenes)$cat=='meta']$G=Gmeta[V(netMetaGenes)[V(netMetaGenes)$cat=='meta']]
  V(netMetaGenes)[V(netMetaGenes)$cat=='gene']$G=genesG[V(netMetaGenes)[V(netMetaGenes)$cat=='gene']]
  V(netMetaGenes)$desc=V(netMetaGenes)$name
  print("metagroup-genes network annotated")
  write_graph(netMetaGenes, 'netMetaGenes.graphml', format ="graphml")
  print("metagroup-genes network writed as netMetaGenes.graphml")
  
  labels=as.character(filtered$desc)
  names(labels)=as.character(filtered$id)
  
  print("gene/term description for vertex annotation prepared")
  
  #Annotation of network metagroups + genes
  V(net_term_genes)$cat='gene'
  V(net_term_genes)[V(net_term_genes)$type==F]$cat='term'
  V(net_term_genes)[V(net_term_genes)$cat=='gene']$G=genesG[V(net_term_genes)[V(net_term_genes)$cat=='gene']$name]
  V(net_term_genes)[V(net_term_genes)$cat=='term']$G=Gterm[V(net_term_genes)[V(net_term_genes)$cat=='term']$name]
  V(net_term_genes)$desc=V(net_term_genes)$name
  V(net_term_genes)[V(net_term_genes)$cat=='term']$desc=DESCterm[V(net_term_genes)[V(net_term_genes)$cat=='term']$name]
  V(net_term_genes)$meta=0
  V(net_term_genes)[V(net_term_genes)$cat=='term']$meta=termMetaAnnot[V(net_term_genes)[V(net_term_genes)$cat=='term']$name]
  print("gene-term network annotated")
  write_graph(net_term_genes, 'net_term_genes.graphml', format ="graphml")
  print("gene-term network writed as net_term_genes.graphml")
  net.bg<-bipartite.projection(net_term_genes)
  print("gene-term network projected")
  netGenes <- net.bg$proj2
  print("Extract gene network from the projection")
  netTerms <- net.bg$proj1
  print("Extract term network from the projection")
  V(netTerms)$label <- labels[V(netTerms)$name]
  print("Term network annotated")
  
  if (subAnalyse){
    print("subanalysis started")
    suppressMessages(require('gplots'))
    suppressMessages(require('VennDiagram'))
    suppressMessages(require('gridExtra'))
    selectGroups=subset(summary_meta,(G>Glim)&(pval<pvalLim))
    print(paste(dim(selectGroups)[1]," metgroups selected for analysis with G > ",Glim,sep=""))
    centralities=list()
    candidates=list()
    topCandidates=list()
    df_candidates=NULL
    for (group in selectGroups$metagroup){
      subnet=induced_subgraph(net_term_genes,V(net_term_genes)[(V(net_term_genes)$meta==group)|(V(net_term_genes)$meta==0)])
      subnet=induced_subgraph(subnet,V(subnet)[igraph::degree(subnet, v = V(subnet))[V(subnet)$name]!=0])
      centrality=betweenness(subnet, v = V(subnet), directed = FALSE, weights = NULL,
                             nobigint = TRUE, normalized = TRUE)
      centralities[[as.character(group)]]=data.frame(name=V(subnet)$name,cat=V(subnet)$cat,G=V(subnet)$G,betweenness=centrality[V(subnet)$name])
      candidates[[as.character(group)]]=subset(centralities[[as.character(group)]],(cat=="gene")&(G!=1))$name
      topCandidates[[as.character(group)]]=subset(subset(centralities[[as.character(group)]][order(centralities[[as.character(group)]]$betweenness,decreasing = T),],(cat=="gene")&(G!=1))[1:as.integer(length(candidates[[as.character(group)]])*topPer/100),],betweenness!=0)$name
      if (is.null(df_candidates)){
        df_candidates=data.frame(name=candidates[[as.character(group)]],mg=rep(group,length(candidates[[as.character(group)]])))
        df_topCandidates=data.frame(name=topCandidates[[as.character(group)]],mg=rep(group,length(topCandidates[[as.character(group)]])))
        df_centrality=centralities[[as.character(group)]]
        df_centrality$mg=group
      }else{
        df_candidates=rbind(df_candidates,data.frame(name=candidates[[as.character(group)]],mg=rep(group,length(candidates[[as.character(group)]]))))
        df_topCandidates=rbind(df_topCandidates,data.frame(name=topCandidates[[as.character(group)]],mg=rep(group,length(topCandidates[[as.character(group)]]))))
        tmp_centrality=centralities[[as.character(group)]]
        tmp_centrality$mg=group
        df_centrality=rbind(df_centrality,tmp_centrality)
      }
    }
    print("Gene centralities calculated for each selected metagroups")
    
    futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
    
    if (length(candidates)>5){
      bestMG=as.character(selectGroups[order(selectGroups$pval),]$metagroup[1:5])
      vcandidates=candidates[bestMG]
      vTopCandidates=topCandidates[bestMG]
    }else{
      vcandidates=candidates
      vTopCandidates=topCandidates
    }
    vp1=venn.diagram(vcandidates,
                     fill=names(vcandidates),
                     filename = NULL,
                     main="Candidates in each metagroups",
                     category.names=paste("mg",names(vcandidates),sep=""))
    vp2=venn.diagram(vTopCandidates,
                     fill=names(vTopCandidates),
                     filename = NULL,
                     main="Top ranked candidates in each metagroups",
                     category.names=paste("mg",names(vTopCandidates),sep=""))
    pdf(file="venn_candidates_metagroups.pdf",width = 8,height = 4)
    grid.arrange(gTree(children=vp1),gTree(children=vp2),ncol=2)
    dev.off()
    print("Venn comparisons of candidates/top candidates amongs metagroups")
    write.table(df_candidates,"candidates.txt",row.names=F,quote=F,sep = "\t")
    write.table(df_topCandidates,"topCandidates.txt",row.names=F,quote=F,sep = "\t")
    write.table(df_centrality,"centralities.txt",row.names=F,quote=F,sep = "\t")
    print("Candidates/top candidates lists printed")
  }
  if (printGraph){
    pdf(file = paste(prefix,"_metas.pdf",sep=''),width = 20,height = 20)
    par(mfrow=c(2,2))
    plotProjMarked(netTerms,markTerms,'Terms Network with metagroups')
    plotProjMarked(netGenes,markGenes,'Genes Network with metagroups')
    plotMetaBiPart(netMetaTerms,'Terms-metagroups bipartite-network',labels,labelSize=1)
    plotMetaBiPart(netMetaGenes,'Genes-metagroups bipartite-network')
    dev.off()
    print("4-plots pdf Terms Network; Genes network ; Terms-metagroups bipartite-network; Genes-metagroups bipartite-network")
    pdf(file = paste(prefix,"_GT.pdf",sep=''),width = 15,height = 15)
    par(mfrow=c(1,1))
    plotGTBiPart(net_term_genes,'Genes-terms bipartite-network',markGenes,labels2=labels)
    dev.off()
    print("Last plot: Genes-terms bipartite-network")
  }
  
}
