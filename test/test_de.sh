#!/usr/bin/env bash
cd de
for i in IPSc NipahPPI PS_pteropus santarius_I_II siRNA_required 
do
	mkdir -p $i
	cd $i
	for j in uniprot_keywords reactome interpro
	do
		echo "##################################"
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
			type="-t \"Biological process\""
		fi
		cmd="Rscript ../../../../scripts/step/database_enrichment.R"
		cmd="$cmd -i ../../../gi/${i}.txt"
		cmd="$cmd -u ../../universe.txt"
		cmd="$cmd -a ../../../databases/${j}.tab"
		cmd="$cmd -o enriched"
		cmd="$cmd $(echo $type)"
		echo "### Working directory ###"
		echo "$PWD"
		echo "### Command ###"
		echo "$cmd"
		echo $cmd | bash
		cd ..
	done
	cd ..
done

