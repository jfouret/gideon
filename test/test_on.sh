#!/usr/bin/env bash
mkdir -p on
cd on
for i in IPSc NipahPPI santarius_I_II siRNA_required 
do
	mkdir -p $i
	cd $i
	echo -e "file\tgroup\ttitle\tannspace\tid\tdesc\tallGenes\tgenes\tprob\tprobMax\tN" > config.txt
	echo -e "$(readlink -f ../../de/PS_pteropus/uniprot_keywords/enriched.tab)\tB\tuniprot_keywords\tNA\tterm\tname\tgenes\tGI\tpval\t0.05\tNgenes" >> config.txt
	echo -e "$(readlink -f ../../de/PS_pteropus/reactome/enriched.tab)\tB\treactome\tNA\tterm\tname\tgenes\tGI\tpval\t0.05\tNgenes" >> config.txt
	echo -e "$(readlink -f ../../de/PS_pteropus/interpro/enriched.tab)\tB\tinterpro\tNA\tterm\tname\tgenes\tGI\tpval\t0.001\tNgenes" >> config.txt
	echo -e "$(readlink -f ../../de/${i}/uniprot_keywords/enriched.tab)\tA\tuniprot_keywords\tNA\tterm\tname\tgenes\tGI\tpval\t0.05\tNgenes" >> config.txt
	echo -e "$(readlink -f ../../de/${i}/reactome/enriched.tab)\tA\treactome\tNA\tterm\tname\tgenes\tGI\tpval\t0.05\tNgenes" >> config.txt
	echo -e "$(readlink -f ../../de/${i}/interpro/enriched.tab)\tA\tinterpro\tNA\tterm\tname\tgenes\tGI\tpval\t0.001\tNgenes" >> config.txt
	
	echo "##################################"
	echo "Analyzing occurence network between PS_pteropus and $i"

	echo -e "conf='config.txt'" > run.R
	echo -e "genesAFile='$(readlink -f ../../gi/PS_pteropus.txt)'" >> run.R
	echo -e "genesBFile='$(readlink -f ../../gi/${i}.txt)'" >> run.R
	echo -e "genesAUniversFile='$(readlink -f ../../de/universe.txt)'" >> run.R
	echo -e "genesBUniversFile=genesAUniversFile" >> run.R
	echo -e "prefixPatternList=c('-.*')" >> run.R
	echo -e "genemapping='$(readlink -f ../../databases/spid2gneName_map.txt)'" >> run.R
	echo -e "source('$(readlink -f ../../../scripts/step/geneTermLinker_2groups.R)')" >> run.R
	echo -e "analyseNetworks_2groups(conf,4,genesAFile,genesBFile,genesAUniversFile,genesBUniversFile,prefixPatternList,printGraph=F,geneMap=genemapping,Glim=0.1,hlim=0.4)" >> run.R
	
	echo "### working directory"
	echo "$PWD"
	echo "### config file"
	echo "$(readlink -f config.txt)"
	echo "### R script used"
	echo "$(readlink -f run.R)"

	Rscript run.R

	cd ..
done