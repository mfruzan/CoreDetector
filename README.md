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

#### Step 3. Set path to CoreDetector Java file lines 73 and 82 in pipeline_Minimap.sh

#### Step 4. Run pipeline on a list of the genome names and paths 
```bash
# run the example set of genomes, a directory "example_out" is created for the alignment results 
./pipeline_Minimap.sh  example/genomes.txt  example_out 20  16
```
## Table of Contents

- [Quick start](#qstart)
- [User Guide](#userguide)
  - [Dependencies](#depends)
  - [Input formats](#iformat)
  - [Options](#options)

## <a name="userguide"></a>User Guide

The coreDetector pipeline can be run in the current directory. 

Or copy MFbio.jar and pipeline_Minimap.sh to a folder of your choice and change directory into that folder.

Before running the pipleline_Minimap.sh make sure it has execute permission.

Important: Ensure you change the path in pipeline_Minimap.sh using a text editor. 

For example: Copy MFbio.jar and pipeline\_Minimap.sh files into an executable bin PATH "/usr/local/bin/". Make sure you change lines 73 and 82 lines in pipeline_Minimap.sh from "~/biotools/MFbio/MFbio.jar to "/usr/local/bin/MFbio.jar"
before you do this.

```bash
sudo cp MFbio.jar pipeline_Minimap.sh /usr/local/bin/

```

## <a name="depends"></a>Dependencies

Make sure Java 1.8 or higher is installed. 
If GSAlign pairwise aligner is used make sure GSAlign (https://github.com/hsinnan75/GSAlign) is installed and its executable files are in PATH environment variable. 

If Minimap2 aligner is used make sure minimap2 (https://lh3.github.io/minimap2/minimap2.html) and K8 javascript engine (https://github.com/attractivechaos/k8) are installed and their executable files are in system path or PATH variable.

## <a name="iformat"></a>Input formats

#### Input 1: Genome sequences in FASTA format


For example 3 genomes in 3 fasta files: genome1.fa, genome2.fa and genome3.fa.


#### Input 2: Is a text file that lists of the name and full path to the FASTA files for each genome. 

The text file has a line for each genome. For example in a text file called 'genomes.txt', each line represents a genome. Each line contains an alias name followed by the full path to its FASTA file separated by a space or Tab. 

```bash
Genome1 /dir/to/fasta/files/genome1.fa
Genome2 /dir/to/fasta/files/genome2.fa
Genome3 /dir/to/fasta/files/genome3.fa
```
Here we can arbitrarily choose genome 1 as the query and the remainder genomes become the subjects. 


## <a name="options"></a>Options

#### CoreDetector pipeline using  Minimap2

```bash
./pipeline_Minimap.sh  genomes.txt  /output/folder 20  16
```
> * First argument is the list of genomes in a text file which points to the FASTA sequences. See example/genome.txt
> * Second argument is a string to create a new output folder or the path to an existing folder. 
>	###### Note: If this folder in the path does not exist it will be created. 
> * Third argument is an integer for the expected genome divergence level and can be any number between 1 and 40. 
> * Fourth argument (optional) is the number of cores/CPUs (default is 4). 
>	###### Note the first three arguments are required. 

#### pipeline using GSAlign 
For the GSAlign pipeline the following GSAlign arguments can be edited and saved in a unix text editor.   

> You can change:
>
> -t (number of threads) 
> -alen (minimum alignment length) 
> -idy (minimum identity between query and subject) 
> -ind (maximum indel length)


