# MultipleSequenceAlignment
Mutliple Sequence Alignment Pipeline for Larger Bacterial and Fungal Genomes. 

Input: Fasta Files of genomes. For example we have 3 fungal genomes in 3 fasta files: genome1.fa, genome2.fa and genome3.fa. Here we can arbitrarily choose genome 1 as query and the rest of genomes as subjects. Then create a text file called 'genomes.txt' that each line contains genome alias name followed by full path to its fasta file separated by a space or Tab. In our example, genomes.txt would be like:
```bash
Genome1 /dir/to/fasta/files/g1.fa
Genome2 /dir/to/fasta/files/g2.fa
Genome3 /dir/to/fasta/files/g3.fa
```
Make sure Java and GSAlign https://github.com/hsinnan75/GSAlign are installed and their executable files are in PATH environment variable.
Then copy MFbio.jar and pipeline.sh to a folder of your choice and CD to this directory. Before running pipleline.sh make sure it has execute permission. Then run

```bash
./pipeline.sh  genomes.txt  /output/folder
```
The first arguments points to the file was created in previous step and second argument is the path to the folder that all output files will be generated inside.
If this folder does not exist it will be created. (both arguments are required).

To change GSAlign arguments just edit pipeline.sh file and save it. You can change -t (number of threads) -alen (minimum alignment length) -idy (minimum identity of query and subject sequences) -ind (maximum indel length) 

To utilize all of CPU/Cores to speed up alignment, set -t parameter of GSAlign to the number of cores in your system.



