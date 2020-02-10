# MultipleSequenceAlignment
Mutliple Sequence Alignment Pipeline for Larger Bacterial and Fungal Genomes. 

Input: Fasta File of whole genomes. For example we have 3 fungal genomes in 3 fasta files: genome1.fa, genome2.fa and genome3.fa. Here we can arbitrarily choose genome 1 as query and the rest of genomes as subjects.

Step 1: Piarwise alignment between genome 1 and 2 using any tools that generates MAF format as output, here we used lastz (http://www.bx.psu.edu/miller_lab/dist/README.lastz-1.02.00/README.lastz-1.02.00a.html), genome 1 as query and genome 2 as subject ,allowing 10% mistmatch and 5% gap:
```bash
./lastz  genome2.fa[multiple]  genome1.fa --ambiguous=n  --ambiguous=iupac  --gfextend --chain --gapped  --identity=90 --continuity=95  --format=maf  --out 1_vs_2.maf
```
Step 2: Extract non overlaping query (genome1) from MAF file generated at step 1 and write into new fasta file. Contig names of new fasta file will include start and end offset, separated by ! special character.
```bash
java -jar  MFbio.jar  showform=no  task=maf2fastaunique  srcdir=1_vs_2.maf  destdir=1_vs_2.fa
```
step 3: Pairwise fasta file generated at step 2 into genome 3:
```bash
./lastz  genome3.fa[multiple]  1_vs_2.fa --ambiguous=n  --ambiguous=iupac  --gfextend --chain --gapped  --identity=90 --continuity=95  --format=maf  --out 1_vs_2_vs_3.maf
```
If we had more genomes for alignments we would just repeated steps 2 and 3 for each genome.

### Now it is time to carry out back tracking algorithm:
Step 4: Extract non-ovelapping entries from last maf file (in our example 1_vs_2_vs_3.maf) and generate new maf file:
```bash
java -jar  MFbio.jar  showform=no  task=mafuniquequery  srcdir=1_vs_2_vs_3.maf  destdir=1_vs_2_vs_3_unique.maf
```
step 5: Generate multiple sequence alignment file:
```bash
java -jar  MFbio.jar  showform=no  task=maf2msa  srcdir=/dir/of/maf/files  file1=1_vs_2_vs_3.maf,1_vs_2.maf  destdir=output.msa
```
where srcdir is directory where maf files are located, file1 comma separated of maf files generated at step 3, ordered from last to the first one.


