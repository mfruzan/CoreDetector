#!/bin/bash

if [ "$1" == "-h" ]; then
  echo ""
  echo "CoreDetector pipeline: for further help see https://github.com/mfruzan/MultipleSequenceAlignment/"
  echo ""
  echo -e "Usage:\n      ./pipeline_Minimap.sh <genome_list> <out_dir> <divergence> <ncpus>\n"
  echo -e "Mandatory options:\n\
	genome_list\tText file lists genome names and paths to FASTA files\n\
	out_dir\t\tnamed directory will be created\n\
	divergence\tlevel of genome divergence, int between 1 and 40\n"
  echo -e "Optional:\n\
	cpus\t\tdefault is 4 cpus\n\
	-h\t\tPrint Help (this message) and exit\n"
  exit 0
fi

# Make sure a directory is given as the second parameter
if [[ -z $2 ]]
then
  echo "Second argument missing (please specify folder name)."
  exit 1
fi
mkdir -p "$2/temp_fasta"
mkdir -p "$2/maf"
mkdir -p "$2/filtered_maf"

if [[ -z $3 ]]
then
  echo "Third argument missing (please specify divergance level, between 1 and 40)."
  exit 1
fi

divergance=$(($3));

cores=4
if [[ ! -z $4 ]]
then
  cores=$(($4));
fi



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
  #echo ${arr[0]}, ${arr[1]}, ${arr[2]};

  if [ $id == 1 ];
  then
    query=${arr[0]}
    queryfile=${arr[1]}
  else
    echo $queryfile;
    echo ${arr[1]};
    twin=${query}_${arr[0]}
    if [[ ${arr[2]} ]]
    then
       echo "BWT index specified.";
      minimap2 -k 19 -w 10 -U 50,500 --rmq=yes -r 10k,100k -g 10k -A 1 -B 1 -O 4,10 -E 2,1 -s 400 -z 400 -N 50  -t ${cores}  --cs=long  --secondary=no  ${arr[1]}  ${queryfile} | paftools.js view -f maf - >$2/maf/${twin}.maf;
    else
       #echo "no index file";
       if (( $divergance <= 5 )) 
       then
         minimap2 -x asm5 --cs=long -t ${cores} --secondary=no ${arr[1]}  ${queryfile} | paftools.js view -f maf - >$2/maf/${twin}.maf;
       elif (( $divergance <= 10 )) 
       then
         minimap2 -x asm10 --cs=long -t ${cores} --secondary=no ${arr[1]}  ${queryfile} | paftools.js view -f maf - >$2/maf/${twin}.maf;       
       elif (( $divergance <= 20 )) 
       then
         minimap2 -x asm20 --cs=long -t ${cores} --secondary=no ${arr[1]}  ${queryfile} | paftools.js view -f maf - >$2/maf/${twin}.maf;
       elif (( $divergance <= 30 )) 
       then
         minimap2  minimap2  -k 19 -w 10 -U 50,500 --rmq=yes -r 10k,100k -g 10k -A 1 -B 2 -O 4,10 -E 2,1 -s 200 -z 200 -N 50  -t ${cores}  --cs=long  --secondary=no ${arr[1]}  ${queryfile} | paftools.js view -f maf - >$2/maf/${twin}.maf;
       else
         minimap2  -k 19 -w 10 -U 50,500 --rmq=yes -r 10k,100k -g 10k -A 1 -B 1 -O 4,10 -E 2,1 -s 400 -z 400 -N 50  -t ${cores}  --cs=long  --secondary=no ${arr[1]}  ${queryfile} | paftools.js view -f maf - >$2/maf/${twin}.maf;
       fi
    fi

    newmaf=${twin}".maf"
    java -jar MFbio.jar --task maf2uniquequery --srcdir $2/maf/${newmaf} --destdir $2/temp_fasta/${twin}".fa" --file1 $2/filtered_maf/${newmaf} --p1 50;
    queryfile=$2/temp_fasta/${twin}".fa";
    maflist=${newmaf}","${maflist};
  fi
  id=$((id+1));
  #echo $id;
done <<<$(cat $1)

#echo $maflist;
#get 80% of system available memroy in Gbyte for java
mem=$(awk '/MemAvailable/ { printf "%.3f \n", $2/1024/1024 }' /proc/meminfo)
echo "available memory(Gb)"${mem};
mem=${mem%.*}
mem=$(($mem))
mem=$(($mem*80/100))
if (($mem < 1))
then
 mem=1
fi

java -jar -Xmx${mem}g MFbio.jar --task maf2msa --srcdir $2/filtered_maf --p1 ${maflist} --destdir $2/concatinated_msa.fa --file1 $2/msa.maf --file2 $1
