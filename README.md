# CoreDetector Multiple Genome Aligner
<!-- badges -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
![Licence pending](https://img.shields.io/badge/licence-TBA-blue)
![Static Badge](https://img.shields.io/badge/version-pending-80b6ff)
[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Fmfruzan%2FCoreDetector&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com)

CoreDetector is a new fast and flexible program that is able to identify the core-genome sequence of larger and more evolutionary diverse genomes. 

## <a name="qstart"></a>Quick start
Installation and configuration of CoreDetector on Linux-based operating systems proceeds as follows.

#### Step 1. Configure your `$PATH` for CoreDetector binary dependencies
CoreDetector depends on the [Minimap2](https://github.com/lh3/minimap2) versatile pairwise aligner (and its related `paftools.js` utility), as well as the [K8 Javascript shell](https://github.com/attractivechaos/k8). The easiest way is to install these to a prepared folder on the system `$PATH` for them, so that they are always available when CoreDetector runs:
```bash
mkdir -p $HOME/bin
echo "export PATH=$HOME/bin:${PATH}" >> $HOME/.bashrc && source $HOME/.bashrc
```

#### Step 2. Download and install Minimap2 (v2.26)
Grab the v2.26 release of Minimap2 from its GitHub repository [here](https://github.com/lh3/minimap2/releases/tag/v2.26). Alternatively, copy-paste the below commands to automatically download, compile and configure Minimap2. (Note: this compilation requires compiler tools and the zlib development headers to be installed: on Ubuntu 22.04, you can easily install these compilation dependencies with `sudo apt-get -y install build-essential zlib1g-dev`.)
```bash
wget "https://github.com/lh3/minimap2/releases/download/v2.26/minimap2-2.26.tar.bz2"
tar -xjf minimap2-2.26.tar.bz2
cd minimap2-2.26 && make
cp minimap2 misc/paftools.js $HOME/bin/
cd ..
```

#### Step 3. Download and install K8 (v1.0)
Get the v1.0 release of the K8 Javascript shell from its GitHub repository [here](https://github.com/attractivechaos/k8/releases/tag/v1.0). Alternatively, execute the following commands to automatically download and configure the precompiled K8 binary:
```bash
wget "https://github.com/attractivechaos/k8/releases/download/v1.0/k8-1.0.tar.bz2"
tar -xjf k8-1.0.tar.bz2
cp k8-1.0/k8-x86_64-Linux $HOME/bin/k8
```

#### Step 4. Install a Java runtime/development kit
[OpenJDK-11](https://openjdk.org/projects/jdk/11/) (or later versions) have been confirmed to work well with CoreDetector. For most Linux systems, these are easily installed via the package manager. E.g., to install OpenJDK-11 (the default JDK) on Ubuntu 22.04:
```bash
sudo apt-get -y install openjdk-11-jdk  # or default-jdk
```

#### Step 5. Download CoreDetector and run an example pipeline
Finally, pull this GitHub repository to download the CoreDetector tool, and run a test case on the provided example set of genomes to confirm that the tool is working correctly.
```bash
git clone https://github.com/mfruzan/CoreDetector.git
cd CoreDetector
chmod +x pipeline_Minimap.sh

./pipeline_Minimap.sh -g example/quick_genomes.txt -o example_out -d 20 -n 16
```

## <a name="dockerqstart"></a>Quick start (using Docker)
Alternatively, easily set up CoreDetector in a Docker container using the provided Dockerfile, which completely automates the installation. For information about setting up Docker on Windows/Mac/Linux and using containers, see [docs.docker.com](https://docs.docker.com/).
```bash
git clone https://github.com/mfruzan/CoreDetector.git
cd CoreDetector
sudo docker build -t coredetector .
sudo docker run -it -v $(pwd)/example:/example coredetector
```
In the interactive shell for the container, you can immediately run the Multiple Genome Aligner tool:
```bash
./pipeline_Minimap.sh -g example/quick_genomes.txt -o example/output -d 20 -n 16
```

## Table of Contents

- [Quick start](#qstart)
- [Quick start (using Docker)](#dockerqstart)
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
Here we can arbitrarily choose Alg130 as the query and the remainder genomes become the subjects. 
## <a name="iformat"></a>Output formats
CoreDetector generates 2 output files in the output folder: <b>msa.maf</b> and <b>concatinated_msa.fa</b> 

msa.maf is standard maf file that each entry of it contains one subject file for each genome. Coordinates and strandness of entries are in respect to the original genome fasta file. This maf file
is approperiate for structural variation detection.

concatinated_msa.fa file is a fasta file that has one entry for each genome file that is constructed by concatinating that genome's subject line from all entries of msa.maf file. The name of each entry is the same as the name of genome introduced in genomes.txt input file. This file is approperiate for phylogenetics tree construction.

## <a name="options"></a>Options

#### CoreDetector pipeline using  Minimap2

```bash
./pipeline_Minimap.sh  -g example/quick_genomes.txt  -o output_folder -d 20  -n 16
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


