# GIDEON

THIS README IS UNDER CONSTRUCTION

# Driver script
A driver script is being built to simplify the use of the software, meanwhile please find a example below :

# Example
## GI = Gene Identification
It is up to you to perform the gene identification, you might want to use the `gideon/scripts/genePythia/genePythia.py` script ( `-h` for help) to perform pubmed data mining

## DE: Database Enrichment
After gene identification you must have a database to perform enrichment analyses. This step can be done using the script `scripts/step/database_enrichment.R` (`-h` for help)

```
user@machine:path/to/gideon$ Rscript scripts/step/database_enrichment.R -h
usage: scripts/step/database_enrichment.R [--] [--help] [--opts OPTS] [-i -I [-u -U [-a -A [-o -O [-p -P [-r -R [-t -T

R program designed to perform an enrichment analysis with topology-based annotation. Over-representation analysis is performed according to the parent-child-intersection method. Multi-testing is done with BH and Bonferoni methods on the subset of term with a minimal p-value below 1e-7 by default. P-values are based on combination analysis as described in the parent-child paper using GNU multiple precision arithmetic libraries and their R implementation


flags:
  -h, --help                    show this help message and exit

optional arguments:
  -x, --opts OPTS                       RDS file containing argument values
  -i, -i I                      input gene list return-separated with no header.
  -u, -u U                      universe gene list return-separated with no header.
  -a, -a A                      annotation table. This table is tab delimited with a header (term,genes,type,name,parents). The column genes contains the comma separated list of genes in a term
  -o, -o O                      output prefix
  -p, -p P                      minimum value for pmin (hihgest possible p-value) [default: 1e-7]
  -r, -r R                      regex pattern to remove from gene name [default: -.*]
  -t, -t T                      comma-separated list of type to keep (keep all by default) [default: ]
```

## ON: Occurence Network

after enrichment of both A and B gene list you should link results file in a file called `config.txt` to further analysis

### config file format

This need to be a tab-separated file with the following header

* `file` : file with enrichment results
* `group` : A or B
* `title` : title for the annotation space
* `annspace` : leave NA
* `id` : (leave as below) column containing term id in the file
* `desc` : (leave as below) column containing term description in the file
* `allGenes` : (leave as below) column containing all genes of the term in the file
* `genes` : (leave as below)  column containing all genes present in both the term and the gene list analyzed in the file
* `prob` : column containing the p-value in the file, unadjusted, BH-adjusted and Boneferonni adjusted p-values are avalable
* `probMax` : Max p-value to keep the term
* `N` : (leave as below) column containing the amount of genes in the term in the file

```
file    group   title   annspace        id      desc    allGenes        genes   prob    probMax N
B/keywords.tab    B       keywords        NA      term    name    genes   GI      pval    0.05  Ngenes
B/interpro.tab    B       interpro        NA      term    name    genes   GI      pval    0.001 Ngenes
B/reactome.tab    B       reactome        NA      term    name    genes   GI      pval    0.05  Ngenes
A/keywords.tab       A       keywords        NA      term    name    genes   GI      pval    0.1   Ngenes
A/interpro.tab       A       interpro        NA      term    name    genes   GI      pval    0.005 Ngenes
A/reactome.tab       A       reactome        NA      term    name    genes   GI      pval    0.1   Ngenes
```

The final step of gideon is encoded within a R function `analyseNetwork_2groups`

### Example of usage:

```
setwd('.')

conf="config.txt"
genesAFile='A.txt' # list of gene id with no header for group A
genesBFile='B.txt'
genesAUniversFile='universe.txt' # list of gene id with no header that are universe of group A
genesBUniversFile=genesAUniversFile
prefixPatternList=c("-.*") # patterns to remove
genemapping="spid2gneName_map.txt" # tabular file with gene id to gene symbol to have meaningfull network image

source('scripts/step/geneTermLinker_2groups.R')

analyseNetworks_2groups(conf,4,genesAFile,genesBFile,genesAUniversFile,genesBUniversFile,prefixPatternList,printGraph=F,geneMap=genemapping,Glim=0.1,hlim=0.4)

```

# LICENCE
All files in this repository are authored by Julien FOURET and licenced under the GNU Public Licence v3.0

This work have been done during a PhD fellowship co-funded by [ViroScan3D](http://www.viroscan3d.com/) and the DGA (Direction Générale de l'Armement) in the context of a [CIFRE-Défense](https://www.ixarm.com/fr/theses-dga-cifre-defense)

```
/*
 *   GIDEON is a software allowing the mining of hidden functional 
 *   relationship between 2 distinct group of genes 
 *   
 *   Copyright (C) 2018  Julien Fouret
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
```