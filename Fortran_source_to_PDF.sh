#!/usr/bin/env bash

for i in `find ./ -name *.f90` ;
do filename=$(basename $i .f90); enscript --highlight=f90 --color $i -o Source_PDF/$filename.pdf ;
done
#find ./ -name *.f90 | xargs