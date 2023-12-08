[![GitHub Downloads](https://img.shields.io/github/downloads/lh3/minimap2/total.svg?style=social&logo=github&label=Download)](https://github.com/lh3/minimap2/releases)

# CoreDetector Multiple Genome Aligner
CoreDetector is a new fast and flexible program that is able to identify the core-genome sequence of larger and more evolutionary diverse genomes. 

## <a name="qstart"></a>Quick start


#### Step 1. Download and install MiniMap2
[![GitHub Downloads](https://img.shields.io/github/downloads/lh3/minimap2/total.svg?style=social&logo=github&label=Download)](https://github.com/lh3/minimap2/releases)

#### Step 2. Download CoreDetector
```bash
git clone https://github.com/mfruzan/CoreDetector.git

# change directory into the CoreDetector directory
cd CoreDetector

# make sure the pipeline is executable
chmod +x pipeline_Minimap.sh
```

#### Step 3. Run pipeline on a list of the genome names and paths 
```bash
# run the example set of genomes, a directory "example_out" is created for the alignment results 
./pipeline_Minimap.sh  -g example/genomes.txt  -o example_out -d 20  -n 16
```
## Table of Contents

- [Quick start](#qstart)
- [User Guide](#userguide)
  - [Dependencies](#depends)
  - [Input formats](#iformat)
  - [Options](#options)

## <a name="userguide"></a>User Guide

The CoreDetector [Manual](https://github.com/mfruzan/CoreDetector/blob/master/Manual.md) explains installation, usage and further analysis examples. 

The CoreDetector pipeline can be run: 

* in the current directory
* in a folder of your choice, just copy the MFbio.jar and pipeline_Minimap.sh files to  the new folder and change directory into the new folder.
* from anywhere by copying MFbio.jar and pipeline_Minimap.sh to an executable bin PATH


**Important**: Ensure you change the path in pipeline_Minimap.sh using a text editor. 

For example: If you copy MFbio.jar and pipeline\_Minimap.sh files into an executable bin PATH "/usr/local/bin/", make sure you change lines 136 and 156 lines in pipeline_Minimap.sh from "MFbio.jar to "/usr/local/bin/MFbio.jar".

```bash
sudo cp MFbio.jar pipeline_Minimap.sh /usr/local/bin/

```

Before running the pipleline_Minimap.sh make sure it has execute permission.

```bash
# make sure the pipeline is executable
chmod +x pipeline_Minimap.sh
```


## <a name="depends"></a>Dependencies

Make sure Java 1.8 or higher is installed. 

If Minimap2 aligner is used make sure minimap2 (https://github.com/lh3/minimap2) and k8 javascript engine (https://github.com/attractivechaos/k8) are installed and their executable files are in system path or PATH variable.
In order to do that, first change directory to where you want minimap2 and k8 being installed and then run following commands:
```bash
# install minimap2
git clone https://github.com/lh3/minimap2
cd minimap2 && make
# install the k8 javascript shell
curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf -
cp k8-0.2.4/k8-`uname -s` k8              # or copy it to a directory on your $PATH
export PATH="$PATH:`pwd`:`pwd`/misc"   
```


If GSAlign pairwise aligner is used make sure GSAlign (https://github.com/hsinnan75/GSAlign) is installed and its executable files are in PATH environment variable. 


## <a name="iformat"></a>Input formats

Is a text file that lists of the name and full path to the FASTA files for each genome. 

The text file has a line for each genome. For example in a text file called 'genomes.txt', each line represents a genome. Each line contains an alias name followed by the full path to its FASTA file separated by a space or Tab. 

```bash
Alg130	example/Alg130.fna
DW5	example/DW5.fna
M4	example/M4.fna
```
Here we can arbitrarily choose genome 1 as the query and the remainder genomes become the subjects. 
## <a name="iformat"></a>Output formats
CoreDetector generates 2 output files in the output folder: <b>msa.maf</b> and <b>concatinated_msa.fa</b> 

msa.maf is standard maf file that each entry of it contains one subject file for each genome. Coordinates and strandness of entries are in respect to the original genome fasta file. This maf file
is approperiate for structural variation detection.

concatinated_msa.fa file is a fasta file that has one entry for each genome file that is constructed by concatinating that genome's subject line from all entries of msa.maf file. The name of each entry is the same as the name of genome introduced in genomes.txt input file. This file is approperiate for phylogenetics tree construction.

## <a name="options"></a>Options

#### CoreDetector pipeline using  Minimap2

```bash
./pipeline_Minimap.sh  -g genomes.txt  -o /output/folder -d 20  -n 16
```
> * -g argument is text file that contains the list of genomes which points to the FASTA sequences. See example/genome.txt
> * -o argument is an output folder where alignment files will be written to. 
>	###### Note: If this folder does not exist it will be created. 
> * -d argument is an integer for the expected genome divergence level and can be any number between 1 and 40. 
> * -n argument (optional) is the number of cores/CPUs (default is 4).
> * -c argument (optional) enables chromosome number matching (1:enable, 0:disable, default is 0) Please note that when enabled then CoreDetector considers a contig name starts with a chromsome number, such as 2B or 14, followed by a white space (or characters '_' , '-' ) 
> * -m argument (optional) sets minimum alignment length in bp (default is 200)
>	###### Note the first three arguments are required. 

#### pipeline using GSAlign 
For the GSAlign pipeline the following GSAlign arguments can be edited and saved in a unix text editor.   

> You can change:
>
> -t (number of threads) 
> -alen (minimum alignment length) 
> -idy (minimum identity between query and subject) 
> -ind (maximum indel length)


