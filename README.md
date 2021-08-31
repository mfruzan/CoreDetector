# MultipleSequenceAlignment
Mutliple Sequence Alignment Pipeline for Larger Bacterial and Fungal Genomes. 

Input: Fasta File of whole genomes. For example we have 3 fungal genomes in 3 fasta files: genome1.fa, genome2.fa and genome3.fa. Here we can arbitrarily choose genome 1 as query and the rest of genomes as subjects.

Step 1: Piarwise alignment between genome 1 and 2 using any tools that generates MAF format as output, here we used lastz (http://www.bx.psu.edu/miller_lab/dist/README.lastz-1.02.00/README.lastz-1.02.00a.html), genome 1 as query and genome 2 as subject ,allowing 10% mistmatch and 5% gap:
```bash
./lastz  genome2.fa[multiple]  genome1.fa --ambiguous=n  --ambiguous=iupac  --gfextend --chain --gapped  --identity=90 --continuity=95  --format=maf  --out /dir/to/maf/1_vs_2.maf
```
If we use GSAlign aligner from https://github.com/hsinnan75/GSAlign, the command would be (using 8 threads):
```bash
./GSAlign -r genome2.fa -q genome1 -o /dir/to/maf/1_vs_2 -no_vcf -t 8 -idy 90 -one -alen 50 -ind 10 -fmt 1
```
Step 2: Extract non overlaping query (genome1) from MAF file generated at step 1 and write into new fasta file. Contig names of new fasta file will include start and end offset, separated by ! special character.
```bash
java -jar  MFbio.jar  --task maf2fastaunique  --srcdir /dir/to/maf/1_vs_2.maf  --destdir /dir/to/temp/fasta/1_vs_2.fa --file1 /dir/to/filtered/maf/1_2.maf  --p1 50
```
This job generated 2 outputs, one temporary fasta file that will beused for the next step pairwise alignment, and one filtered maf file that will be  used at last stage backtracking algorithm. --p1 is minimum HSP length that goes into outputs HSP shorter than 50bp will be skipped)

step 3: Pairwise fasta file generated at step 2 into genome 3:
lastz:
```bash
./lastz  genome3.fa[multiple]  /dir/to/temp/fasta/1_vs_2.fa --ambiguous=n  --ambiguous=iupac  --gfextend --chain --gapped  --identity=90 --continuity=95  --format=maf  --out /dir/to/maf/1_vs_2_vs_3.maf
```
GSAlign:
```bash
./GSAlign -r genome3.fa -q /dir/to/temp/fasta/1_vs_2.fa -o /dir/to/maf/1_vs_2_vs_3 -no_vcf -t 8 -idy 90 -one -alen 50 -ind 10 -fmt 1
```

Step 4: Extract non-ovelapping entries from last maf file (in our example 1_vs_2_vs_3.maf) and generate new maf file:
```bash
java -jar  MFbio.jar  showform=no  task=mafuniquequery  --srcdir /dir/to/maf/1_vs_2_vs_3.maf  --destdir /dir/to/temp/fasta/1_vs_2_vs_3.fa --file1 /dir/to/filtered/maf/1_3.maf  --p1 50
```
If we had more genomes for alignments we would just repeated steps 2 and 3 for each genome.

### Now it is time to carry out back tracking algorithm:
step 5: Generate multiple sequence alignment file:
```bash
java -jar  MFbio.jar  --task maf2msa  --srcdir /dir/of/filtered/maf/  --p1 1_3.maf,1_2.maf  --destdir concatinated_msa.fa --file1 msa.maf --file2 genomes.txt
```
where --srcdir is directory where filtered maf files are located, --p1 is comma separated filtered maf file names generated at step 3, ordered from last to the first one.
--file1 is generated multiple sequence alignment file. --destdir is concatinated entries of msa.maf in fasta format. 

--file2 is input tab separated text file in the following format:  Genome_Alias_Name   Genome_fasta_file

one line for each genome, for our example genomes.txt could be like:
```bash
Genome1 /dir/to/fasta/files/g1.fa
Genome2 /dir/to/fasta/files/g2.fa
Genome3 /dir/to/fasta/files/g3.fa
```
