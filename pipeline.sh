#!/bin/bash

# Make sure a directory is given as the second parameter
if [[ -z $2 ]]
then
  echo "Second argument missing (please specify folder name)."
  exit 1
fi
mkdir -p "$2/temp_fasta"
mkdir -p "$2/maf"
mkdir -p "$2/filtered_maf"

declare -i id=1;
maflist="";
query="";
queryfile="";
lines=$(cat $1 | wc -l)
echo Number of lines in file $1 : $lines;
#lines=$((lines+1));
while read line
do
  read -a arr <<< $line
  if [ $id == 1 ];
  then
    query=${arr[0]}
    queryfile=${arr[1]}
  else
    echo $queryfile;
    echo ${arr[1]};
    twin=${query}_${arr[0]}
    GSAlign -r ${arr[1]} -q ${queryfile} -o $2/maf/${twin} -no_vcf -t 32 -idy 85 -sen -alen 50 -ind 20 -fmt 1;
    newmaf=${twin}".maf"
    java -jar ./MFbio.jar --task maf2uniquequery --srcdir $2/maf/${newmaf} --destdir $2/temp_fasta/${twin}".fa" --file1 $2/filtered_maf/${newmaf} --p1 50;
    queryfile=$2/temp_fasta/${twin}".fa";
    maflist=${newmaf}","${maflist};
  fi
  id=$((id+1));
  #echo $id;
done <<<$(cat $1)

echo $maflist;
java -jar ~/biotools/MFbio/MFbio.jar --task maf2msa --srcdir $2/filtered_maf --p1 ${maflist} --destdir $2/concatinated_msa.fa --file1 $2/msa.maf --file2 $1
