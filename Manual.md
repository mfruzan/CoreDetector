# Computing and analysis of core-genome alignments with CoreDetector

### Mario Fruzangohar & Paula Moolhuijzen
##### mario.fruzangohar@adelaide.edu.au
##### 31-08-2023

The analysis of the core-genome alignments of conserved sequence is important to measure genetic changes between population individuals and show not only the evolutionary relationships within a population but provide further insight into core gene functions and how these may shift over time or geography. This is however complicated by the limitation of current tools containing the functionality to process larger and more diverse speices. CoreDetector is a fast and flexible program that is able to identify the core-genome sequence of larger and more evolutionary diverse genomes. This document explains how to install and use CoreDetector, for core genome alignment, phylogenetic tree reconstruction and gives an example for the comparison of phylogenetic trees.

## Table of Contents

- [1. Installation](#install)
  - [Dependencies](#depends)
- [2. Input formats](#iformat)
- [3. Basic usage](#options)
- [4. Further analyses](#analysis)
  - [Phylogeny](#phylo)
  - [Comparing trees](#comp)


## <a name="install"></a> 1. CoreDetector installation


### <a name="depends"></a>Dependencies

CoreDetector depends on the fast and efficient pairwise alignment tool Minimap2.

**`Step 1.`** Download the precompiled binaries of the latest Minimap2 aligner version. 

[![GitHub Downloads](https://img.shields.io/github/downloads/lh3/minimap2/total.svg?style=social&logo=github&label=Download)](https://github.com/lh3/minimap2/releases)

or

```bash
curl -L https://github.com/lh3/minimap2/releases/download/v2.26/minimap2-2.26_x64-linux.tar.bz2 | tar -jxvf -
./minimap2-2.26_x64-linux/minimap2
```

**`Step 2.`** Make sure the MiniMap2 executable files are in your PATH environment variable. 

**`Step 3.`** Make sure Java 1.8 or higher is installed. 


### <a name="download"></a> Obtaining and setting up CoreDetector

**`Step 4.`** Download the CoreDetector package from GitHub

```bash 
git clone https://github.com/mfruzan/CoreDetector.git
```

**`Step 5.`** Change directory into the CoreDetector directory

```bash
cd CoreDetector
```

**`Step 6.`** Make sure the pipeline is executable by changing the file permissions.

```bash
chmod +x pipeline_Minimap.sh
```

**`Step 7.`** 
The CoreDetector pipeline can be run in the current directory. 

```bash
./pipeline_Minimap.sh
```

Or you can alternatively copy the files MFbio.jar and pipeline_Minimap.sh to a folder of your choice that is in your executable PATH.

For example:
You can sudo copy MFbio.jar and pipeline_Minimap.sh into an executable bin PATH "/usr/local/bin/".

```bash
sudo cp MFbio.jar pipeline_Minimap.sh /usr/local/bin/
```

**`Step 8.`** 
Ensure you change the path to the executables in pipeline_Minimap.sh using a text editor. 

To find the full path of the current directory or where the executables are you can use the pwd command

```bash
pwd .
```

Make sure you change lines 73 and 82 lines in pipeline_Minimap.sh from "~/biotools/MFbio/MFbio.jar to the new path.

For the 'example' the path would be: 

```bash
/usr/local/bin/MFbio.jar
```



## <a name="iformat"></a> 2. Data input formats for CoreDetector

CoreDetector requires two data inputs. The first is fasta formatted genome files (input 1) and the second is a text file that lists of the name and full path to the FASTA files for each genome (input 2). 

Published fasta formatted genomes can be downloaded using NCBI tools datasets and dataformats. Install these tools and we will download a data set of genomes in the next section.

**`Step 1.`** Download and install NCBI tools datasets and dataformat to use in the next section.

```bash
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat'
chmod +x datasets dataformat
sudo cp datasets dataformat /usr/local/bin/

```

**`Step 2.`** The text file list of genomes (input 2).

The genome text file should have a line for each genome. You can create a text file called 'genomes.txt', each line represents a genome. Each line contains an alias name for the genome followed by the full path to its fasta file separated by a either space or Tab. 

```bash
Genome1 /dir/to/fasta/files/g1.fa
Genome2 /dir/to/fasta/files/g2.fa
Genome3 /dir/to/fasta/files/g3.fa
```

We have provided in our GitHub a file that lists the genomes we will use in the next section "example/genomes.txt". In this file genome 1 becomes the query and the remainder genomes become the subjects. 


## <a name="options"></a>3. CoreDetector usage

After CoreDetector is setup, the Minimap2 pipeline can be run.

**`Step 1.`** Let us first view the CoreDetector options (arguments)

The required and mandatoy options for the Minimap2 pipeline (pipeline_Minimap.sh) can be  viewed using the help (-h) argument as follows.

```bash
pipeline_Minimap.sh -h
```
Output:

```
CoreDetector pipeline: for further help see https://github.com/mfruzan/MultipleSequenceAlignment/

Usage:
      ./pipeline_Minimap.sh <genome_list> <out_dir> <divergence> <ncpus>

Mandatory options:
	genome_list	Text file lists genome names and paths to FASTA files
	out_dir		named directory will be created
	divergence	level of genome divergence, int between 1 and 40

Optional:
	cpus		default is 4 cpus
	-h		Print Help (this message) and exit
```

Note the argument order is required and the first three arguments are mandatory:

> * First argument is the list of genomes in a text file which points to the FASTA sequences. See example/genome.txt
> * Second argument is a string to create a new output folder or the path to an existing folder. 
>	 Note: If this folder in the path does not exist it will be created. 
> * Third argument is an integer for the expected genome divergence level and can be any number between 1 and 40. 
> * Fourth argument (optional) is the number of cores/CPUs (default is 4).


**`Step 2.`**  Download the fasta formatted genomes from NCBI 

To see how CoreDetector works we will apply it to a set of 23 fungal genomes using the NCBI datasets tool to download the genomes based on a single BioProject number PRJNA315205

```bash
mkdir genomes
cd genomes
datasets download genome accession PRJNA315205 --filename bioproject_dataset.zip
unzip bioproject_dataset.zip 
ls genomes/ncbi_dataset/data/
```

Output:

```
GCA_003171515.3  GCA_003231365.1  GCA_022544805.1  GCA_022788425.1  GCA_022813025.1
GCA_003171545.1  GCA_003231415.2  GCA_022578365.1  GCA_022788435.1  GCA_022813065.1
GCA_003231325.2  GCA_003231425.2  GCA_022578395.1  GCA_022788445.1  GCA_022837075.1
GCA_003231345.1  GCA_008692205.1  GCA_022788405.1  GCA_022788505.1  assembly_data_report.jsonl
GCA_003231355.1  GCA_022544795.1  GCA_022788415.1  GCA_022788515.1  dataset_catalog.json
```

**`Step 3.`** Run the CoreDetector pipeline using the list of genomes

Note: The genomes.txt is supplied in the example folder. You can move the file to the current path.

```bash
mv example/genomes.txt .
./pipeline_Minimap.sh  genomes.txt  output 20  16
ls output/
```

After CoreDetector has completed (~5mins) you will see the final alignment length and find the results in the output folder. The concatinated_msa.fa file contains the extracted fasta alignment that can then be used for further analysis.


Output:

```
...
start writing to output file...
Length of Alignment 33764393
Total Query Length 32082566
Done!

concatinated_msa.fa  filtered_maf  maf  msa.maf  temp_fasta
```

## <a name="analysis"></a> 4. Futher analysis of the core genome


### <a name="phylo"></a>Phylogeny

Now that we have the core alignment in fasta format we can conduct phylogenetic analysis. You can use a tool of your preference, but here we will use the tool Phylip.

**`Step 1.`** First convert the concatinated_msa.fa to phylip format using fasta\_to\_phylip.py (~1 min).  

```bash
sudo cp scripts/fasta_to_phylip.py /usr/local/bin/
cd output/
fasta_to_phylip.py &>log.convert
ls concatinated_msa.*
```

Output:


```
concatinated_msa.fa  concatinated_msa.phylip
```

**`Step 2.`** The download and installation instuctions for Phylip can be found here https://evolution.genetics.washington.edu/phylip/getme-new1.html

On ubuntu 20.04+ you can install Phylip using apt-get:

```bash
sudo apt-get install phylip
```

**`Step 3.`** Calculate the distances between species. Phylip is interacive on the commandline. Copy the concatinated_msa.phylip to infile (which will be detected by Phylip). You will be prompted to select a distance formula 'D'. Here we selected Jukes-Cantor, then selected 'Y' to accept the selection. The outfile contains the distance matrix for the multiple sequence alignment.

```bash
cp concatinated_msa.phylip infile
phylip dnadist infile
```

```
Nucleic acid sequence Distance Matrix program, version 3.697

Settings for this run:
  D  Distance (F84, Kimura, Jukes-Cantor, LogDet)?  Jukes-Cantor
  G          Gamma distributed rates across sites?  No
  C            One category of substitution rates?  Yes
  W                         Use weights for sites?  No
  L                       Form of distance matrix?  Square
  M                    Analyze multiple data sets?  No
  I                   Input sequences interleaved?  Yes
  0            Terminal type (IBM PC, ANSI, none)?  ANSI
  1             Print out the data at start of run  No
  2           Print indications of progress of run  Yes

  Y to accept these or type the letter for one to change
Y
Distances calculated for species
    Genome1      ......................
    Genome23     .....................
    Genome22     ....................
    Genome21     ...................
    Genome20     ..................
    Genome19     .................
    Genome18     ................
    Genome17     ...............
    Genome16     ..............
    Genome15     .............
    Genome14     ............
    Genome13     ...........
    Genome12     ..........
    Genome11     .........
    Genome10     ........
    Genome9      .......
    Genome8      ......
    Genome7      .....
    Genome6      ....
    Genome5      ...
    Genome4      ..
    Genome3      .
    Genome2   

Distances written to file "outfile"

Done.

```
**`Step 4`** Move the outfile to the new infile to create a phylogenetic tree. We will use the neighbor joing method here.

```bash
mv outfile infile 
phylip neighbor infile
```

```
Neighbor-Joining/UPGMA method version 3.697

Settings for this run:
  N       Neighbor-joining or UPGMA tree?  Neighbor-joining
  O                        Outgroup root?  Yes, at species number 23
  L         Lower-triangular data matrix?  No
  R         Upper-triangular data matrix?  No
  S                        Subreplicates?  No
  J     Randomize input order of species?  No. Use input order
  M           Analyze multiple data sets?  No
  0   Terminal type (IBM PC, ANSI, none)?  ANSI
  1    Print out the data at start of run  No
  2  Print indications of progress of run  Yes
  3                        Print out tree  Yes
  4       Write out trees onto tree file?  Yes


  Y to accept these or type the letter for one to change
Y

Cycle  20: species 3 (   0.00068) joins species 4 (   0.00020)
Cycle  19: species 5 (   0.00063) joins species 7 (   0.00014)
Cycle  18: species 12 (   0.00072) joins species 13 (   0.00045)
Cycle  17: node 5 (   0.00046) joins species 11 (   0.00068)
Cycle  16: node 12 (   0.00045) joins species 15 (   0.00045)
Cycle  15: species 20 (   0.00045) joins species 21 (   0.00035)
Cycle  14: species 19 (   0.00062) joins node 20 (   0.00009)
Cycle  13: node 19 (   0.00011) joins species 23 (   0.00052)
Cycle  12: species 6 (   0.00029) joins species 8 (   0.00092)
Cycle  11: node 3 (   0.00071) joins node 5 (   0.00031)
Cycle  10: node 6 (   0.00007) joins species 10 (   0.00028)
Cycle   9: node 3 (   0.00008) joins node 6 (   0.00002)
Cycle   8: node 3 (   0.00019) joins species 9 (   0.00028)
Cycle   7: node 3 (   0.00056) joins node 12 (   0.00060)
Cycle   6: species 14 (   0.00038) joins node 19 (   0.00047)
Cycle   5: species 17 (   0.00044) joins species 22 (   0.00364)
Cycle   4: species 2 (   0.00289) joins node 3 (   0.00010)
Cycle   3: node 2 (   0.00002) joins node 14 (   0.00011)
Cycle   2: node 2 (   0.00007) joins node 17 (   0.00010)
Cycle   1: node 2 (   0.00002) joins species 18 (   0.00054)
last cycle:
 species 1  (   0.00010) joins node 2  (   0.00001) joins species 16  (   0.00043)

Output written on file "outfile"

Tree written on file "outtree"

Done.
```
**`Step 5.`** The outtree can then be viewed using FigTree. On ubuntu 20.04+ FigTree can be installed using apt-get install.

```bash
sudo apt-get install figtree 
figtree outtree 
```
### <a name="comp"></a>Comparing trees

The CoreDetector phylogenetic tree can then be compared to trees generated using other tools (e.g. Phylonium and parsnp). We have generated trees based on 29 fungal genomes, from our CoreDetector manuscript, available in the cloned folder r_analysis/manuscript

You will be using RStudio for this analysis. If you do not have it already, you can follow the instructions to download, install and launch RStudio 

https://rstudio-education.github.io/hopr/starting.html

**`Step 1`** In the CoreDetector directory open r_analysis/MGA-Analysis.Rmd with RStudio.

**`Step 2`** Run the code chunks **"1. set the path to files"** and **"2. load the libraries required"**.

**`Step 3`** In the next section "Fungal pathogen phylogenic tree comparison"  the code chunk  **"3. Load CoreDetector, Parsnp and Phylonium trees"** loads three trees that are then sorted (using the TreeTools package) on each node into a consistent order, so that node rotation does not obscure similarities between similar trees for comparison.

**`Step 4`** We will now compare the first two trees by running the code chunk for **"4. Compare the generated trees from CoreDetector and Parsnp"** (using ape package) which returns a detatiled report of this comparison. 

![Figure 1](./r\_analysis/compare\_tree1\_tree2.png "Figure 1") 
**Figure 1.** CoreDetector tree (tree 1) with splits incommon to the Parsnp generated tree (left) and Parsnp tree (tree 2) shows splits incommon to the CoreDetector generated tree (right).


**`Step 5`** We next compare the first and third tree set by running the code chunk for **"5. Compare the generated trees from CoreDetector and Phylonium"** (Figure 2).

![Figure 2](./r\_analysis/compare\_tree1\_tree3.png "Figure 2") 
**Figure 2.** CoreDetector tree (tree 1) with splits incommon to the Phylonium generated tree (left) and the Phylonium tree (tree 3) shows splits incommon to CoreDetector generated tree (right).

**`Step 6`** This final step plots the three trees including the geographic origin of the isolates highlighted by running **"6. For each tree highlight the geographic origin of the isolates"**. Figure 3 shows the three trees and the geographic source of the isolates highlighted.

![Figure 3](./r\_analysis/tree\_clades.png "Figure 3")
**Figure 3.** *Pyrenophora tritici-repentis* fungal pathogen phylogenetic tree topology comparisons. Trees generated from CoreDetector (left), Parsnp (centre) and Phylonium (right) show three groups related to geographic locations, Europe (violet), Australia (blue) and North Africa (tan). The Ptr isolate identifiers are shown in all three trees.



