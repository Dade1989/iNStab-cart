#!/bin/bash

#filename=$1
fileout=$1

rm $fileout

for f in ../q*.f90 ../start_sparse_kit.f90 ../sparse_matrix_operations.f90
do
  echo "Finding the SUBROUTINEs in file $f"
  echo "File $f +++++++++++++++++++++++++++++++++++++++++" >> $fileout
  
  #while read -r word _ 
  while read -r line 
  do
      name=$(echo $line | awk '{print $1;}')
      if [ "$name" == "SUBROUTINE" ]; then
         subname=$(echo $line | awk '{print $2}')
         echo "SUBROUTINE - $subname" >> $fileout
      fi
  done < "$f"
  
done
echo "              ... end of the script"
