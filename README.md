# MultipleSequenceAlignment
Mutliple Sequence Alignment Pipeline for Larger Bacterial and Fungal Genomes. 

Input: Fasta Files of genomes. For example we have 3 fungal genomes in 3 fasta files: genome1.fa, genome2.fa and genome3.fa. Here we can arbitrarily choose genome 1 as query and the rest of genomes as subjects. Then create a text file called 'genomes.txt' that each line contains genome alias name followed by full path to its fasta file separated by a space or Tab. In our example, genomes.txt would be like:
```bash
Genome1 /dir/to/fasta/files/g1.fa
Genome2 /dir/to/fasta/files/g2.fa
Genome3 /dir/to/fasta/files/g3.fa
```
Make sure Java 1.8 or higher is installed. If GSAlign pairwise aligner is used make sure GSAlign https://github.com/hsinnan75/GSAlign is installed and its executable files are in PATH environment variable. If Minimap2 aligner is used make sure minimap2 https://lh3.github.io/minimap2/minimap2.html and K8 javascript engine https://github.com/attractivechaos/k8 are installed and their exutable files are in system path.
Then copy MFbio.jar and pipeline.sh to a folder of your choice and CD to that folder. Before running pipleline.sh make sure it has execute permission. Then run

```bash
./pipeline_Minimap.sh  genomes.txt  /output/folder 20  16
```
The first arguments points to the file was created in previous step and second argument is the path to the folder that all output files will be generated inside.
If this folder does not exist it will be created. Third argument referes to divergance level and can be a numver from 1-40. Fourth argument is number of cores/CPUs(default is 4). The first three arguments are required. 

To change GSAlign arguments just edit pipeline.sh file and save it. You can change -t (number of threads) -alen (minimum alignment length) -idy (minimum identity between query and subject) -ind (maximum indel length).





