#!/bin/bash

while getopts g:d:n:o:m:c: option
do 
    case "${option}"
        in
        g)genome=${OPTARG};;
        d)diverg=${OPTARG};;
        n)ncores=${OPTARG};;
        o)outdir=${OPTARG};;
        m)minlen=${OPTARG};;
        c)chrom=${OPTARG};;

    esac
done


if [ "$1" == "-h" ]; then
  echo ""
  echo "CoreDetector pipeline: for further help see https://github.com/mfruzan/CoreDetector.git/"
  echo ""
  echo -e "Usage:\n      ./pipeline_Minimap.sh -g <genome_list> -o <out_dir> -d <divergence> -n <ncpus>  -m <minlength>  -c <chromosome>\n"
  echo -e "Mandatory options:\n\
	genome_list\tText file lists genome names and paths to FASTA files\n\
	out_dir\t\tnamed directory will be created\n\
	divergence\tlevel of genome divergence, int between 1 and 40\n"
  echo -e "Optional:\n\
	ncpus\t\tdefault is 4 cpus\n\
        minlength(Minimum alignment length)\t\tdefault is 200bp\n\
        chromosome(chromosome matching)\t\tdefault is 0 or disabled, to enable set it to 1\n\
	-h\t\tPrint Help (this message) and exit\n"
  exit 0
fi


if [[ -z $genome ]]
then
  echo "genomes file missing (please specify -g)."
  exit 1
fi



if [[ -z $outdir ]]
then
  echo "output directory missing (please specify -o)."
  exit 1
fi
mkdir -p "$outdir/temp_fasta"
mkdir -p "$outdir/maf"
mkdir -p "$outdir/filtered_maf"

if [[ -z $diverg ]]
then
  echo "Divergence level is missing (please specify -d, between 1 and 40)."
  exit 1
fi

divergence=$(($diverg));

cores=4
if [[ ! -z $ncores ]]
then
  cores=$(($ncores));
fi


mlen=200
if [[ ! -z $minlen ]]
then
  mlen=$(($minlen));
fi

chromosome=0
if [[ ! -z $chrom ]]
then
  chromosome=$(($chrom));
fi


declare -i id=1;
declare -i I_param=4;
declare -i K_param=1;
maflist="";
query="";
queryfile="";
lines=$(cat $genome | wc -l)
echo Number of lines in file $genome : $lines;
#lines=$((lines+1));
while read line
do
  read -a arr <<< $line
  #echo ${arr[0]}, ${arr[1]}, ${arr[2]};

  if [ $id == 1 ];
  then
    query=${arr[0]}
    queryfile=${arr[1]}
    #set size of query file
    sz=$(ls -l ${arr[1]} | cut -d ' ' -f 5)
    gb=$((sz/1073741824))
    echo Genome size gb : $gb;
    if (($gb > 3))
    then
      K_param=$((gb+2));
      I_param=$((gb+2));
    fi   
  else
    echo $queryfile;
    echo ${arr[1]};
    twin=${query}_${arr[0]}
    if [[ ${arr[2]} ]]
    then
       echo "BWT index specified.";
      minimap2 -k 19 -w 10 -U 50,500 --rmq=yes -r 10k,100k -g 10k -A 1 -B 1 -O 4,10 -E 2,1 -s 400 -z 400 -N 50  -t ${cores}  --cs=long  --secondary=no  ${arr[1]}  ${queryfile} | paftools.js view -f maf - >$outdir/maf/${twin}.maf;
    else
       #echo "no index file, index is built on the fly.";
       if (( $divergence <= 5 )) 
       then
         minimap2 -x asm5  -I${I_param}g  -K${K_param}g  -H --cs=long -t ${cores} --secondary=no ${arr[1]}  ${queryfile} | paftools.js view -f maf - >$outdir/maf/${twin}.maf;
       elif (( $divergence <= 10 )) 
       then
         minimap2 -x asm10 -I${I_param}g  -K${K_param}g  -H --cs=long -t ${cores} --secondary=no ${arr[1]}  ${queryfile} | paftools.js view -f maf - >$outdir/maf/${twin}.maf;
       elif (( $divergence <= 20 )) 
       then
         minimap2 -x asm20  -I${I_param}g  -K${K_param}g -H  --cs=long -t ${cores} --secondary=no ${arr[1]}  ${queryfile} | paftools.js view -f maf - >$outdir/maf/${twin}.maf;
       elif (( $divergence <= 30 )) 
       then
         minimap2  -I${I_param}g  -K${K_param}g -H -k 19 -w 10 -U 50,500 --rmq=yes -r 10k,100k -g 10k -A 1 -B 2 -O 4,10 -E 2,1 -s 200 -z 200 -N 50  -t ${cores}  --cs=long  --secondary=no ${arr[1]}  ${queryfile} | paftools.js view -f maf - >$outdir/maf/${twin}.maf;
       else
         minimap2  -I${I_param}g  -K${K_param}g  -H -k 19 -w 10 -U 50,500 --rmq=yes -r 10k,100k -g 10k -A 1 -B 1 -O 4,10 -E 2,1 -s 400 -z 400 -N 50  -t ${cores}  --cs=long  --secondary=no ${arr[1]}  ${queryfile} | paftools.js view -f maf - >$outdir/maf/${twin}.maf;
       fi
    fi

    newmaf=${twin}".maf"
    java -jar MFbio.jar --task maf2uniquequery --srcdir $outdir/maf/${newmaf} --destdir $outdir/temp_fasta/${twin}".fa" --file1 $outdir/filtered_maf/${newmaf} --p1 ${mlen}  --p2  ${chromosome};
    queryfile=$outdir/temp_fasta/${twin}".fa";
    maflist=${newmaf}","${maflist};
  fi
  id=$((id+1));
  #echo $id;
done <<<$(cat $genome)

#echo $maflist;
#get 80% of system available memroy in Gbyte for java
mem=$(awk '/MemAvailable/ { printf "%.3f \n", $2/1024/1024 }' /proc/meminfo)
echo available memory ${mem} gb;
mem=${mem%.*}
mem=$(($mem))
mem=$(($mem*80/100))
if (($mem < 1))
then
 mem=1
fi

java -jar -Xmx${mem}g MFbio.jar --task maf2msa --srcdir $outdir/filtered_maf --p1 ${maflist} --destdir $outdir/concatinated_msa.fa --file1 $outdir/msa.maf --file2 $genome ;
