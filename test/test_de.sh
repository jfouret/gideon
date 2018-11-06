#!/usr/bin/env bash
cd de
for i in IPSc NipahPPI PS_pteropus santarius_I_II siRNA_required 
do
	mkdir -p $i
	cd $i
	for j in reactome interpro uniprot_keywords
	do
		echo "Analyzing $i list with $j database"
		mkdir -p $j
		cd $j
		if [ $j = "reactome" ]
		then
			type=""
		fi
		if [ $j = "interpro" ]
		then
			type="-t 'Conserved_site,Domain,Family,Homologous_superfamily,Repeat'"
		fi
		if [ $j = "uniprot_keywords" ]
		then
			type="-t 'Biological process'"
		fi

		Rscript ../../../../scripts/step/database_enrichment.R \
			-i ../../../gi/${i}.txt \
			-u ../../universe.txt \
			-a ../../../databases/${j}.tab \
			-o enriched \
			$(echo $type)
	done
	cd ..
done

