SEDPREFIX=$(shell echo ${prefix} | sed 's/\//\\\//g')
GITVERSION=0.0.0
GITREPOSED=$(shell pwd | sed 's/\//\\\//g')
#$(shell git describe --tags | sed 's/^v//g')

bin/gideon:
	mkdir -p bin
	sed -e "s/SEDMATCHVERSION/${GITVERSION}/g" scripts/gideon.py | sed -e "s/SEDMATCHLIB/${GITREPOSED}\/scripts/g" > bin/gideon
	chmod 777 bin/gideon

.PHONY : install
install:
	mkdir -p ${prefix}/lib
	-rm ${prefix}/lib/gideon
	cp -r scripts ${prefix}/lib/gideon
	mkdir -p ${prefix}/bin
	sed -e 's/SEDMATCHVERSION//g' scripts/gideon.py | \ 
	sed -e "s/SEDMATCHLIB/${SEDPREFIX}/gideon/g" > ${prefix}/bin/gideon


.PHONY : clean 
clean :
	-rm bin/

.PHONY : remove
remove:
	-rm -r ${prefix}/lib/gideon
	-rm -r ${prefix}/bin

.PHONY : doc
doc : 
	Rscript -e "rmarkdown::render('README.md')"