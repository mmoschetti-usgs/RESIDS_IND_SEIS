#!/bin/bash
# script info
# 
 
# 
for ff in `ls ORIG/ | grep csv`; do 
  cat ORIG/$ff | sed 's/"//g' | sed 's/NA/-999/g' | sed 's/ //g'  > $ff
done
