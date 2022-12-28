/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package mfruzan.common;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import mfruzan.algorithm.BlastHSP;
import mfruzan.algorithm.MapUtil;

/**
 *
 * @author a1195806
 */
public class BlastResult {
    BufferedReader reader = null;
    private String fileName = null;
    byte format = 0; 
    public BlastResult(String fileName, byte format){
        this.fileName = fileName;
        this.format = format;
        try{
            if (fileName.endsWith(".gz")){
                InputStream filein = new FileInputStream(fileName);
                InputStream gzipin = new GZIPInputStream(filein);
                Reader ir = new InputStreamReader(gzipin);
                this.reader = new BufferedReader(ir);
            }else
                this.reader = new BufferedReader(new FileReader(this.fileName));
        } catch (Exception ex) {
            System.out.println(ex.getMessage());
        }
    }
    public BlastResult(){
        
    }
    public void mergeWithMutationFile(String mutFile, double min_identity_percent, final boolean exact_match){
        //We blast an unknown sequence (as query) to a mutation reference fasta (as subject), if mutation exists, then we work out position of the mutation in the query.
        // This method, accepts Blast result and MUT file produces a new MUT file (new mutation files are the ones that fall into subject sequence of the Blast)
        // Blast result output format 7 :  qseqid sseqid length pident mismatch gapopen gaps qstart qend sstart send qseq sseq (13 fields)
        // if exact_match is false, as long as contig of mutations matches with contig of subject, it will be written to the file
        // if exact match is true mutation should fall between start and end of subject
        

        Writer writer = null;
        MUTProcessor mut = new MUTProcessor(mutFile, true, true);
        
        try {
            
            File f = new File(this.fileName);
            writer = new FileWriter(new File(f.getParent(), "query_mutations.MUT"));
            String blast_line = reader.readLine();
            while (blast_line != null && !blast_line.isEmpty()){
                if (!blast_line.startsWith("#")){
                    // this line is not a comment; comment lines start with #
                    String[] arr = blast_line.split("\\s");
                    double pident = Double.parseDouble(arr[3]);
                    if (pident >= min_identity_percent){
                        // we process this alignment
                        int qstart = Integer.parseInt(arr[7]);
                        int qend = Integer.parseInt(arr[8]);
                        int sstart = Integer.parseInt(arr[9]);
                        int send = Integer.parseInt(arr[10]);
                        
                        int sleft = 0;
                        int sright = 0;
                        if (sstart <= send){
                            // stand is plus/plus
                            sleft = sstart; 
                            sright = send;
                        }else{
                            // strand is plus/minus
                            sleft = send;
                            sright = sstart;
                        }
                        // in this way sleft always <= sright
                        List<Mutation>  mutations = null;
                        if (exact_match)
                          mutations = mut.getMutations(arr[1], sleft , sright);
                        else 
                          mutations = mut.getMutations(arr[1]);
                        
                        String qseq = arr[11]; // this sequence may contains - for a gap opening in query
                        String sseq = arr[12]; // this sequence may contains - for gap openning in subject
                        
                        // now we have all mutations, 
                        for(Mutation mu : mutations){
                            if(exact_match){
                                // we work out exact mutation position in the query (we know the position in the subject)
                                int counter = sleft-1; // initialize the counter
                                int idx = 0;
                                for (; idx<sseq.length(); idx++ ){

                                    if (sseq.charAt(idx) != '-')  // if character is gap we should skip it without increasing it 
                                        counter++; 
                                    if (counter==mu.position)
                                        break;


                                }
                                // when we get here , idx contains index of position in alignment/sseq where mutation occurs, 
                                //because the length of sseq and qseq are the same(in the alignment) so the same idx position in the qseq is also candidate position for mutation in the query
                                // we only need to skip gaps (-) in the query sequence to work out its pure position in the query
                                counter = qstart -1;
                                for(int i =0; i<=idx; i++){
                                    if (qseq.charAt(i) != '-')  // if character is gap we should skip it without increasing it 
                                        counter++;                                 
                                }
                                // when we get here counter has the exact location of mutation in the query sequence
                                // so now we can write a new MUT file, this time based on query sequence                            
                                writer.write(String.format("%s\t%d\t%d\t%s\n", arr[0] ,counter , counter, mu.print12()));
                            }else // print mutations in the same contig as query alignment                           
                               writer.write(String.format("*%s\t%s\t%d\t%s\n", arr[0] , arr[1], mu.position, mu.print12()));
                        }// for all mutations
                        
                    }
                }                
                blast_line = reader.readLine();
            } //while
            
            
        } catch (Exception ex) {System.out.println("mergeWithMutationFile::" + ex.getMessage());}
        finally{
            try {
                writer.close();
                reader.close();                
            } catch (IOException ex) {}
        }    

    }
   
    
    public void getIndividualSequence(String pileupFile, String pileupIndex, double min_identity_percent){
        // This method reads the blast result format 7 and for each subject sequence builds its equivalent from a pileup file and its index,
        //Output file is blast result file with a new column added to that for individual sequence
        Writer writer = null;
        PileupProcessor pileup = new PileupProcessor();
        
        try {
            
            File f = new File(this.fileName);
            writer = new FileWriter(new File(f.getParent(), "blast_result_individual_seq.txt"));
            String blast_line = reader.readLine();
            while (blast_line != null && !blast_line.isEmpty()){
                if (!blast_line.startsWith("#")){
                    // this line is not a comment; comment lines start with #
                    String[] arr = blast_line.split("\\s");
                    double pident = Double.parseDouble(arr[3]);
                    if (pident >= min_identity_percent){
                        // we process this alignment
                        int sstart = Integer.parseInt(arr[9]);
                        int send = Integer.parseInt(arr[10]);
                        
                        int sleft = 0;
                        int sright = 0;
                        if (sstart <= send){
                            // stand is plus/plus
                            sleft = sstart; 
                            sright = send;
                        }else{
                            // strand is plus/minus
                            sleft = send;
                            sright = sstart;
                        }
                        // we assume pileupfile has only one individual
                        String seq = pileup.getMLESampleSequenceUsingIndex(pileupFile, pileupIndex, 1, 1, 1, arr[1], sleft, sright, 0, null);
                        writer.write(String.format("%s\t%s\n", blast_line, seq));
                        
                    } // if
                }                
                blast_line = reader.readLine();
            } //while
            
            
        } catch (Exception ex) {System.out.println("getIndividualSequence::" + ex.getMessage());}
        finally{
            try {
                writer.close();
                reader.close();
                
            } catch (IOException ex) {}
        }    
        
    }
    public void findUniqueHit(double thereshold){
        //Looking at the hits for each query, a hit can be unique if its score is significantly better that its next hit, and significantly less than its previous hit.
        //Output will go to blast_unique_hits.txt
        //if thershold = 0.2 or 20% then for example score 1000 is significantly bigger than 700 but not from 850
         FileWriter blastWriter = null; // new fasta file with updated contig names is created by this
         
         try{ 
             blastWriter = new FileWriter(new File(new File(fileName).getParent(), "unique_hit_blast.txt"));
             //blastWriter = new FileWriter("c:/data/unique_hit_blast.txt");
             
             String last_line = null;
             String second_last_line = null;
             int last_score = 0;
             int second_last_score = 0;
             String last_qcontig = null;
             String blast_line = reader.readLine();
             while (blast_line != null && !blast_line.isEmpty()){
                 
                 String[] arr = blast_line.split("\\s");
                 String qcontig =  arr[0]; 
                 //double percent = Double.parseDouble(arr[7]);
                 //int alignment_len = Integer.parseInt(arr[6]); // not used at this stage
                 //double evalue = Double.parseDouble(arr[12]);
                if(!qcontig.equals(last_qcontig)){
                     // this means contig has been changed, reset variables
                     last_score = 0;
                     second_last_score = 0;
                     last_line = null;
                     second_last_line = null;
                 }

                 int score  = Integer.parseInt(arr[11]);
                 if((1-thereshold)*last_score > score){
                     if( (second_last_score==0) || ((1-thereshold)*second_last_score > last_score)){
                         //then last_line is good hit 
                         blastWriter.write(last_line+"\n");
                     }
                 }
                 
                 last_qcontig = qcontig;
                 //blastWriter.write(blast_line+"\n");
                 second_last_line = last_line;
                 last_line = blast_line;
                 second_last_score = last_score;
                 last_score = score;
                 blast_line = reader.readLine();
             }         
         } catch (Exception ex) {System.out.println(ex.getMessage());}
        finally{
            //fsf.close();
            try {
                blastWriter.close();
            } catch (IOException ex) {}
        }  
    }
    
    public void addLookupField(String src, String dest, String outFile, int  dest_column_idx, double min_overlap){
        // This method does not need to be here, because it is not specific to blast result file, but because its method is similar to the methods here i write it here
        // It accepts 2 region files both file should have contig name and start and end position, and both files should be ordered based on the contig name (alphabetically) and start position(numerically) (the reason we want everything sorted is that we dont have index for these files and for things to be quick we need them sorted)
        // then for every region in source file if we find ovelapped region (with minumum ovelap percentage) in dest file, then we add column number dest_column_idx (1 based) to the end of source line
        // eventually we write updated src file to output 
        // This method assumes first , second and third columns in both files are contig_name start_index and end_index respectively.
        String srcLine = null;
        String destLine = null;
        FileWriter writer = null;
        BufferedReader src_reader = null;
        BufferedReader dest_reader = null;
        String last_contig = null;
        try{
            //read first destination line
            destLine = dest_reader.readLine();
            String[] arrd = destLine.split("\\s");
            String dcontig = arrd[0];
            int dstart = Integer.parseInt(arrd[1]);
            int dend = Integer.parseInt(arrd[2]);
            srcLine = src_reader.readLine();
            while(srcLine!=null){
                String[] arrs = srcLine.split("\\s");
                String scontig = arrs[0];
                int sstart = Integer.parseInt(arrs[1]);
                int send = Integer.parseInt(arrs[2]);
                // look for this contig in dest file 
                while(scontig.compareTo(dcontig)<0){
                    destLine = dest_reader.readLine();
                    if (destLine == null || destLine.isEmpty())
                        break;
                    arrd = destLine.split("\\s");
                    dcontig = arrd[0];
                    dstart = Integer.parseInt(arrd[1]);
                    dend = Integer.parseInt(arrd[2]);
                }
                // we get here if destLine==null or scontig==dcontig or scontig>dcontig
                if (destLine == null || destLine.isEmpty())
                    break; // we have reached to the end dest file and there is no hope to find further overlap
                if (scontig.equals(dcontig)){
                    // we have to find all the overlaps of current source region with dest region
                }
                
                
                
                srcLine = src_reader.readLine();
            }//while
        } catch (Exception ex) {System.out.println("addLookupField" + ex.getMessage());}
            finally{
            try {
                writer.close();
                src_reader.close();
                dest_reader.close();
            } catch (Exception ex) {}
        }         
        
    }  
    
      
    public void findHomeologAllele(String varFile, String homologFile ,String subjectFastaFile, String subjectFastaIndex, String queryFastaFile, String outFile){
        // This method accepts a var file and a homolog txt file of format (query_contig start_pos_query_contig	end_pos_query_contig	subject_contig	subject_start_pos	subject_end_pos	+/-) that is manually created from the output of findHomologSubject (no header) and shows which regions of query contig is homolog to which region of subject contig..
        // if a polymorphic position in var file falls inside a region in homeologFile, then finds its homolog base and writes it to a new outFile(contig+position+ref-base(e.g. A)+variation(e.g A/C)+homolog base(e.g. T)). 
        // Regions in homolog file are sorted and are not overlapped. positions in varFile are also sorted based on position. It is assumed order of query contigs in varFile is the same queryFastaFile
        // queryFastaFile    must have just one contig and that is the query contig
        String homologLine = null;
        String varLine = null;
        FileWriter writer = null;
        BufferedReader homolog_reader = null;
        BufferedReader var_reader = null;
        int step = 0;
        try{
            //FastaSequenceIndex  qfi = new FastaSequenceIndex(new File(queryFastaIndex));
            //IndexedFastaSequenceFile qf = new IndexedFastaSequenceFile(new File(queryFastaFile), qfi);
            //For query fasta file we dont need index, because there is just one contig inside it
            FastaSequenceFile query_fasta = new FastaSequenceFile(new File(queryFastaFile), true);
            ReferenceSequence first_sequence = query_fasta.nextSequence(); // the only and the first sequence
            String query_sequence = new String(first_sequence.getBases(), "UTF8");
            step = 1;
            FastaSequenceIndex fsi = new FastaSequenceIndex(new File(subjectFastaIndex));
            IndexedFastaSequenceFile ifsi =new IndexedFastaSequenceFile(new File(subjectFastaFile) , fsi);

            writer = new FileWriter(outFile);
            //writer.write("Contig\tPosition\tReference Base\tVariation\tHomolog Base\n");
            homolog_reader = new BufferedReader(new FileReader(homologFile));
            
            homologLine = homolog_reader.readLine(); // read first entry of homolog file (no header)
            
            String[] arrs = homologLine.split("\\s");
            // initialize variable with first entry of homolog file
            String qcontig = arrs[0].trim();// because all query regions refer to just one contig as query, so this variable is not used            
            int qstart = Integer.parseInt(arrs[1]);
            int qend = Integer.parseInt(arrs[2]);
            String scontig = arrs[3].trim();
            int sstart = Integer.parseInt(arrs[4]);
            int send = Integer.parseInt(arrs[5]);
            String orientation = arrs[6].trim();
            // perform global alignment between query and subject regions:
            ReferenceSequence refseq = ifsi.getSubsequenceAt(scontig, sstart, send);
            String subject_region = null;
            if (orientation.equals("+"))
              subject_region =  new String(refseq.getBases(), "UTF-8");
            else
              subject_region =  Helper.getRC(new String(refseq.getBases(), "UTF-8")); // bio java also has api to get RC of a sequence                       

            String query_region = query_sequence.substring(qstart-1, qend);// change 1 based indexes of blast to 0-based indexes 
            //String query_region = new String(refseq.getBases(), "UTF-8");
            // Now do global alignment subject region against query region and show the alignment profile
            DNASequence seq1 = new DNASequence(query_region, AmbiguityDNACompoundSet.getDNACompoundSet());
            DNASequence seq2 = new DNASequence(subject_region, AmbiguityDNACompoundSet.getDNACompoundSet());
            SubstitutionMatrix<NucleotideCompound> matrix = SubstitutionMatrixHelper.getNuc_mf2();  // 5: match score -4: mismatch score -2: there is N
            SimpleGapPenalty gapP = new SimpleGapPenalty();
            gapP.setOpenPenalty((short)4);  
            gapP.setExtensionPenalty((short)4);

            SequencePair<DNASequence, NucleotideCompound> psa = Alignments.getPairwiseAlignment(seq1, seq2, Alignments.PairwiseSequenceAlignerType.GLOBAL, gapP, matrix);

            writer.write(homologLine+"\n");
            writer.write(psa.toString()+"\n");
            String q_aligned = psa.getQuery().getSequenceAsString();// Will contains ACGTN-
            String s_aligned = psa.getTarget().getSequenceAsString();// Will contains ACGTN-
            
            step = 2;
            var_reader = new BufferedReader(new FileReader(varFile)); 
            varLine = var_reader.readLine();
            while (varLine != null && !varLine.isEmpty()){
                step = 3;
                String[] arr = varLine.split("\\s");
                String contig = arr[0].trim();
                int var_pos = Integer.parseInt(arr[1]);
                String refBase = arr[2];
                String variation = arr[3];
                // now we start looking for regions in homologFile until find a region that current polymorphic position is in it or pass by it.               
                while( var_pos >= qstart){
                   if (var_pos <= qend){
                       // we found the right region that contains polymorphic position, we extract query and subject sequence and align them
                       //first we determin where in alignment profile var_pos will be located
                       int cc = qstart-1;
                       int i=0;
                       for(; ;i++){
                           char ch = q_aligned.charAt(i);
                           if (ch!='-')
                               cc++;
                           if(cc== var_pos)
                               break;
                       }
                       // when we get here i+1 will be the index of alignment profile
                       char original = q_aligned.charAt(i);
                       char homolog = s_aligned.charAt(i); // it could be ACGTN-
                       writer.write(String.format("%s\t%d\t%s\t%s\t%c\t%c\n", contig, var_pos, refBase, variation,  original, homolog));
                       //writer.write(String.format("%s[%c]%s\n", q_aligned.substring(0, i), original, q_aligned.substring(i+1)));
                       //writer.write(String.format("%s[%c]%s\n", s_aligned.substring(0, i), homolog, s_aligned.substring(i+1)));
                       break;
                   } 
                   //so we get here if var_pos > qend, so we need to read the next entry/line in homolog file
                   homologLine = homolog_reader.readLine();
                   if (homologLine==null  || homologLine.isEmpty())
                       break;
                   
                   // so we get here we still have entry in the homolog file
                   arrs = homologLine.split("\\s");
                   qcontig = arr[0].trim();
                   qstart = Integer.parseInt(arrs[1]);
                   qend = Integer.parseInt(arrs[2]);
                   scontig = arrs[3].trim();
                   sstart = Integer.parseInt(arrs[4]);
                   send = Integer.parseInt(arrs[5]);
                   orientation = arrs[6].trim();
                   refseq = ifsi.getSubsequenceAt(scontig, sstart, send);
                   subject_region = null;
                   if (orientation.equals("+"))
                      subject_region =  new String(refseq.getBases(), "UTF-8");
                   else
                      subject_region =  Helper.getRC(new String(refseq.getBases(), "UTF-8"));                        

                   
                   query_region = query_sequence.substring(qstart-1, qend);// change 1 based indexes of blast to 0-based indexes 
                    // Now do global alignment subject region against query region and show the alignment profile
                   seq1 = new DNASequence(query_region, AmbiguityDNACompoundSet.getDNACompoundSet());
                   seq2 = new DNASequence(subject_region, AmbiguityDNACompoundSet.getDNACompoundSet());
                   psa = Alignments.getPairwiseAlignment(seq1, seq2, Alignments.PairwiseSequenceAlignerType.GLOBAL, gapP, matrix);
                   writer.write(homologLine+"\n");
                   writer.write(psa.toString()+"\n");
                   q_aligned = psa.getQuery().getSequenceAsString();
                   s_aligned = psa.getTarget().getSequenceAsString();
                   
                }//while homolog file
                step = 4;
                varLine = var_reader.readLine();
               
            }

            step = 5;
 
        } catch (Exception ex) {System.out.println("addHomeologBases" + ex.getMessage() + step+ ":" + varLine );}
            finally{
            try {
                writer.close();
                homolog_reader.close();
                var_reader.close();
            } catch (Exception ex) {}
        }         
    }
    
    
        public void findHomeologAllele2(String varFile, String homologFile , String queryFastaFile, String outFile){
        //findHomeologAllele2 is very similar to findHomeologAllele, but it processes a homolog file that subject alignment block is already there, so we wont need to do pairwise alignment 
        // This method accepts a var file and a homolog txt file of format (q_contig q_start q_end	s_contig	s_start	s_end	+/-     q_alignment_block     s_alignment_block) This file can be produced by MAFProcessor
        // if a polymorphic position in var file falls inside a region in homeologFile, then finds its homolog base and writes it to a new outFile(contig+position+ref-base(e.g. A)+variation(e.g A/C)+homolog base(e.g. T)). 
        // Regions in homolog file are sorted and can be partially overlapped (MAFProcessor sorts but does not guarantee they are distincs). positions in varFile are also sorted . 
        // queryFastaFile    must have just one contig and that is the query contig
        String homologLine = null;
        String varLine = null;
        FileWriter writer = null;
        BufferedReader homolog_reader = null;
        BufferedReader var_reader = null;
        int step = 0;
        try{
            //FastaSequenceIndex  qfi = new FastaSequenceIndex(new File(queryFastaIndex));
            //IndexedFastaSequenceFile qf = new IndexedFastaSequenceFile(new File(queryFastaFile), qfi);
            //For query fasta file we dont need index, because there is just one contig inside it
            FastaSequenceFile query_fasta = new FastaSequenceFile(new File(queryFastaFile), true);
            ReferenceSequence first_sequence = query_fasta.nextSequence(); // the only and the first sequence
            String query_sequence = new String(first_sequence.getBases(), "UTF8");
            step = 1;

            writer = new FileWriter(outFile);
            //writer.write("Contig\tPosition\tReference Base\tVariation\tHomolog Base\n");
            homolog_reader = new BufferedReader(new FileReader(homologFile));
            
            homologLine = homolog_reader.readLine(); // read first entry of homolog file (no header)
            
            String[] arrs = homologLine.split("\\s");
            // initialize variable with first entry of homolog file
            String qcontig = arrs[0].trim();// because all query regions refer to just one contig as query, so this variable is not used            
            int qstart = Integer.parseInt(arrs[1]);
            int qend = Integer.parseInt(arrs[2]);
            String scontig = arrs[3].trim();
            int sstart = Integer.parseInt(arrs[4]);
            int send = Integer.parseInt(arrs[5]);
            String orientation = arrs[6].trim();
            String q_alignment_block = arrs[7].trim(); // Will contains ACGTN-
            String s_alignment_block = arrs[8].trim(); // Will contains ACGTN- 

            String query_region = query_sequence.substring(qstart-1, qend);// change 1 based indexes of blast to 0-based indexes 
            //String query_region = new String(refseq.getBases(), "UTF-8");
            // Now do global alignment subject region against query region and show the alignment profile

            //writer.write(query_region+"\n");
            //writer.write(q_alignment_block+"\n");
            //writer.write(s_alignment_block+"\n");
            
            step = 2;
            var_reader = new BufferedReader(new FileReader(varFile)); 
            varLine = var_reader.readLine();
            while (varLine != null && !varLine.isEmpty()){
                step = 3;
                String[] arr = varLine.split("\\s");
                String contig = arr[0].trim();
                int var_pos = Integer.parseInt(arr[1]);
                String refBase = arr[2];
                String variation = arr[3];
                // now we start looking for regions in homologFile until find a region that current polymorphic position is in it or pass by it.               
                while( var_pos >= qstart){
                   if (var_pos <= qend){
                       // we found the right homolog region that contains polymorphic position,
                       //first we determin where in alignment profile var_pos will be located
                       int cc = qstart-1;
                       int i=0;
                       for(; ;i++){
                           char ch = q_alignment_block.charAt(i);
                           if (ch!='-')
                               cc++;
                           if(cc== var_pos)
                               break;
                       }
                       // when we get here i will be the index of alignment profile
                       char original = q_alignment_block.charAt(i); // it can not be -
                       char homolog_allele = s_alignment_block.charAt(i); // it could be ACGTN-
                       writer.write(String.format("%s\t%d\t%s\t%s\t%c\t%c\n", contig, var_pos, refBase, variation,  original, homolog_allele));
                       //writer.write(String.format("%s[%c]%s\n", q_aligned.substring(0, i), original, q_aligned.substring(i+1)));
                       //writer.write(String.format("%s[%c]%s\n", s_aligned.substring(0, i), homolog, s_aligned.substring(i+1)));
                       break; // break from inner loop
                   } 
                   //so we get here if var_pos > qend, so we need to read the next entry/line in homolog file
                   homologLine = homolog_reader.readLine();
                   if (homologLine==null  || homologLine.isEmpty())
                       break; // break from outer loop
                   
                   // so we get here we still have entry in the homolog file
                   arrs = homologLine.split("\\s");
                   qcontig = arr[0].trim();
                   qstart = Integer.parseInt(arrs[1]);
                   qend = Integer.parseInt(arrs[2]);
                   scontig = arrs[3].trim();
                   sstart = Integer.parseInt(arrs[4]);
                   send = Integer.parseInt(arrs[5]);
                   orientation = arrs[6].trim();
                   q_alignment_block = arrs[7].trim(); // Will contains ACGTN-
                   s_alignment_block = arrs[8].trim(); // Will contains ACGTN- 
                   
                   query_region = query_sequence.substring(qstart-1, qend);// change 1 based indexes of blast to 0-based indexes 
                   //writer.write(query_region+"\n");
                   //writer.write(q_alignment_block+"\n");
                   //writer.write(s_alignment_block+"\n");
                   
                }//while homolog file
                step = 4;
                varLine = var_reader.readLine();
               
            }

            step = 5;
 
        } catch (Exception ex) {System.out.println("addHomeologBases" + ex.getMessage() + step+ ":" + varLine );}
            finally{
            try {
                writer.close();
                homolog_reader.close();
                var_reader.close();
            } catch (Exception ex) {}
        }         
    }

    
    
    
    public void findHomologSubject(String subjectFastaFile, String subjectFastaIndex, String resultFile,   double min_identity, double max_evalue){
        // This method assumes query is just one big contigs that has been aligned to many contigs as subject.In this situation blast sorts the result based on subject contig and for each subject sorts them based on HSP scores from top to bottom
        // The aim of this method is to determine if a subject contig is a true homolg of a subregion on the query contig
        // we set some rules, subject regions must have the same orientation. The following rule in the paranthesis was removed on 15/9/2016 because whole subject can be big enough to accomodate whole query(and from start to the end must account for at least 95% of the length of the contig). 
        // If we measure from the start to the end of the query, it must also be equal to the subject region with 10% variation in length to account for In/Dels.
        // For each subject contig we produce a single line (if it follows above rules) 
        // We also report percentage of whole query contig that is covered and also percentage of subject region (not subject contig, because we dont expect whole query contig will map to whole subject contig)
        
        String blast_line = null;
        String last_subject_contig = null; // blast result is sorted by subject contig
        FileWriter writer = null;
        TreeMap<Integer, Integer> query_posMap_plus = new TreeMap(); // this Map holds start and end positions of each hit in the query, and sorted based on its key (qstart)
        TreeMap<Integer, Integer> query_posMap_minus = new TreeMap(); // this Map holds start and end positions of each hit in the query, and sorted based on its key (qstart)
        
        TreeMap<Integer, Integer> subject_posMap_plus = new TreeMap(); // this Map holds start and end positions of each positive hit in the last_contig, and sorted based on its key (sstart)
        TreeMap<Integer, Integer> subject_posMap_minus = new TreeMap(); // this Map holds start and end positions of each negative hit in the last_contig, and sorted based on its key (sstart)
        int step = 0;
        try{ 
            
            FastaSequenceIndex fsi = new FastaSequenceIndex(new File(subjectFastaIndex));
            //IndexedFastaSequenceFile ifsi =new IndexedFastaSequenceFile(new File(subjectFastaFile) , fsi);

            writer = new FileWriter(resultFile);
            writer.write("Query Start\tQuery End\tQuery Region Coverage\tSubject Contig\tSubject Start\tSubject End\tOrientation\tSubject Region Covered\n");
            blast_line = reader.readLine();
            step = 1;
            while (blast_line != null && !blast_line.isEmpty()){
                String[] arr = blast_line.split("\\s");
                String sseqid = arr[3].trim();
                double percent = Double.parseDouble(arr[7]);
                int alignment_len = Integer.parseInt(arr[6]); // not used at this stage
                double evalue = Double.parseDouble(arr[12]);
                int qstart = Integer.parseInt(arr[1]);
                int qend = Integer.parseInt(arr[2]);
                int sstart = Integer.parseInt(arr[4]);
                int send = Integer.parseInt(arr[5]);
                //if(sseqid.equals("3AS_3321827"))
                 //           alignment_len = 0;
                
                if(!sseqid.equals(last_subject_contig)){
                    // process  last contig posMap
                    if (last_subject_contig != null){
                        // process data related to last_subject_contig
                        // now we investigate 3 posMaps
                        MapUtil util = new MapUtil();
                        int qpluslen = 0;
                        int qminuslen = 0;
                        int biggest_plus_len = 0;
                        int biggest_minus_len = 0;
                        int subject_plus_len = 0;
                        int subject_minus_len = 0;
                        int plus_difference = -1;
                        int minus_difference = -1;
                        int query_plus_covered = 0;
                        int query_minus_covered = 0;
                        int plus_covered = 0;
                        int minus_covered = 0;
                        step = 2;
                        //int subject_conig_size = (int)fsi.getContigSize(last_subject_contig);
                        step = 3;
                       
                        if (query_posMap_plus.size() > 0){
                            qpluslen = util.right_end(query_posMap_plus) - query_posMap_plus.firstKey()+1;
                            query_plus_covered = qpluslen - util.notCovered(query_posMap_plus);                            
                        }
                        if (query_posMap_minus.size() > 0){
                            qminuslen = util.right_end(query_posMap_minus) - query_posMap_minus.firstKey()+1;
                            query_minus_covered = qminuslen - util.notCovered(query_posMap_minus);                            
                        }

                        step = 4;
                        if (subject_posMap_plus.size()>0){
                            subject_plus_len = util.right_end(subject_posMap_plus) - subject_posMap_plus.firstKey() + 1;
                            plus_difference = Math.abs(qpluslen-subject_plus_len);
                            plus_covered = subject_plus_len - util.notCovered(subject_posMap_plus);
                        }
                        step = 5;
                        if (subject_posMap_minus.size()>0){
                            subject_minus_len = util.right_end(subject_posMap_minus) - subject_posMap_minus.firstKey() + 1;
                            minus_difference = Math.abs(qminuslen-subject_minus_len);
                            minus_covered = subject_minus_len - util.notCovered(subject_posMap_minus);
                        }
                        step = 6;
                        //System.out.println(last_subject_contig+":"+subject_conig_size+":"+qlen);
                        biggest_plus_len = Math.max(qpluslen, subject_plus_len);
                        biggest_minus_len = Math.max(qminuslen, subject_minus_len);
                        if(plus_difference>=0  && plus_difference < 0.2*biggest_plus_len ){
                            //plus-plus hits are selected
                            step = 61;
                            writer.write(String.format("%d\t%d\t%.2f\t%s\t%d\t%d\t+\t%.2f\n", query_posMap_plus.firstKey(), util.right_end(query_posMap_plus) , (float)query_plus_covered/qpluslen ,last_subject_contig, subject_posMap_plus.firstKey(), util.right_end(subject_posMap_plus), (float)plus_covered/subject_plus_len));
                            
                        }
                        if (minus_difference>=0   && minus_difference < 0.2*biggest_minus_len  ){
                            
                            //plus-minus hits are selected
                            step = 62;
                            writer.write(String.format("%d\t%d\t%.2f\t%s\t%d\t%d\t-\t%.2f\n", query_posMap_minus.firstKey(), util.right_end(query_posMap_minus) , (float)query_minus_covered/qminuslen , last_subject_contig, subject_posMap_minus.firstKey(), util.right_end(subject_posMap_minus), (float)minus_covered/subject_minus_len));
                        }
                        // else subject contig is not selected as homolog candidate
                        step = 7;
                    }
                    // empty the Map for new contig
                    query_posMap_plus.clear();
                    query_posMap_minus.clear();
                    subject_posMap_plus.clear();
                    subject_posMap_minus.clear();
                    last_subject_contig = sseqid;
                    step = 8;
                }
                 
                if (evalue<max_evalue && percent>min_identity ){
                    
                   if(sstart<send){
                       // plus-plus 
                       subject_posMap_plus.put(sstart, send);
                       query_posMap_plus.put(qstart, qend);
                   }else{
                       //plus-minus
                       subject_posMap_minus.put(send, sstart);
                       query_posMap_minus.put(qstart, qend);
                   }
                }
                 step = 9;
                blast_line = reader.readLine();
                
            }//while loop
             // the last subject contig still needs to be processed
            MapUtil util = new MapUtil();
            int qpluslen = 0;
            int qminuslen = 0;
            int biggest_plus_len = 0;
            int biggest_minus_len = 0;           
            int subject_plus_len = 0;
            int subject_minus_len = 0;
            int plus_difference = -1;
            int minus_difference = -1;
            int query_plus_covered = 0;
            int query_minus_covered = 0;
            int plus_covered = 0;
            int minus_covered = 0;
            
            //int subject_conig_size = (int)fsi.getContigSize(last_subject_contig);

            biggest_plus_len = Math.max(qpluslen, subject_plus_len);
            biggest_minus_len = Math.max(qminuslen, subject_minus_len);
            
            if (query_posMap_plus.size() > 0){
                qpluslen = util.right_end(query_posMap_plus) - query_posMap_plus.firstKey()+1;
                query_plus_covered = qpluslen - util.notCovered(query_posMap_plus);                            
            }
            if (query_posMap_minus.size() > 0){
                qminuslen = util.right_end(query_posMap_minus) - query_posMap_minus.firstKey()+1;
                query_minus_covered = qminuslen - util.notCovered(query_posMap_minus);                            
            }

            if (subject_posMap_plus.size()>0){
                subject_plus_len = util.right_end(subject_posMap_plus) - subject_posMap_plus.firstKey() + 1;
                plus_difference = Math.abs(qpluslen-subject_plus_len);
                plus_covered = subject_plus_len - util.notCovered(subject_posMap_plus);
            }
            if (subject_posMap_minus.size()>0){
                subject_minus_len = util.right_end(subject_posMap_minus) - subject_posMap_minus.firstKey() + 1;
                minus_difference = Math.abs(qminuslen-subject_minus_len);
                minus_covered = subject_minus_len - util.notCovered(subject_posMap_minus);
            }
            //System.out.println(last_subject_contig+":"+subject_conig_size+":"+qlen);
            if(plus_difference>=0  && plus_difference < 0.2*biggest_plus_len ){
                //plus-plus hits are selected
                writer.write(String.format("%d\t%d\t%.2f\t%s\t%d\t%d\t+\t%.2f\n", query_posMap_plus.firstKey(), util.right_end(query_posMap_plus) , (float)query_plus_covered/qpluslen ,last_subject_contig, subject_posMap_plus.firstKey(), util.right_end(subject_posMap_plus), (float)plus_covered/subject_plus_len));

            }
            if (minus_difference>=0  && minus_difference < 0.2*biggest_minus_len ){
                //plus-minus hits are selected
                writer.write(String.format("%d\t%d\t%.2f\t%s\t%d\t%d\t-\t%.2f\n", query_posMap_minus.firstKey(), util.right_end(query_posMap_minus) , (float)query_minus_covered/qminuslen , last_subject_contig, subject_posMap_minus.firstKey(), util.right_end(subject_posMap_minus), (float)minus_covered/subject_minus_len));
            }

             
        } catch (Exception ex) {System.out.println("findHomologSubject" + ex.getMessage() + " at step "+step);}
        finally{
            try {
                writer.close();
            } catch (Exception ex) {}
        } 
        
    }
   

    
    public void findAssemblyHomoeleogRegions(String queryFastaFile, String annotationFile ,String resultFile, double min_identity, double max_evalue){
        // exome capture is blasted to reference 1 assembly, we need to find all 3 copies of each segment in exome capture
        //reference combines with txt annotation file for chromosome 
        //annotationFile has got 2 columns, coulmn one is contig name and column 2 is chromosome name i.e. 1A,1B,...
        
        //Note: i used variable max_evalue as min alignment length, you need to return it to its original state
        BufferedReader annotationReader = null;
        Map<String, String> popseqMap = new HashMap(); // a map of contig names -> Chromosome

        FastaSequenceFile fsf = null;  // we sequentially read query fasta file so we dont need its index
        String blast_line = null;
        FileWriter writer = null;
        FileWriter gff = null;
        String last_query_contig = null; // blast result is sorted by query contig and then subject contig
        int last_query_contig_size = 0;
        long sum_query_contig_size = 0;
        
        int step = 0;
        // Instead of using Treemap we have decided to use another structure. We use LinkedList for each chromosome, that is ordered based on query region
        // There shouldnt be a significant overlap between query regions
        //For each query contig we need at most 3 copies: A,B and D
        LinkedList<BlastHSP> A_genome = new LinkedList();
        LinkedList<BlastHSP> B_genome = new LinkedList();
        LinkedList<BlastHSP> D_genome = new LinkedList();
        int notInsertedHSP_A = 0;
        int notInsertedHSP_B = 0;
        int notInsertedHSP_D = 0;
        
        long total_covered_query_inserted=0;
        long total_covered_query = 0;
        //long total_covered_query_nogap = 0;

        try{ 
            if (annotationFile != null){
                annotationReader = new BufferedReader(new FileReader(annotationFile));
                annotationReader.readLine(); // read header
                String line = annotationReader.readLine();
                while (line != null && !line.isEmpty()){
                    String[] arr = line.split("\t");
                    String contig = arr[0];
                    String chromosome = arr[1].toUpperCase();                    
                    //double position = Double.parseDouble(arr[2]);                   
                    popseqMap.put(contig, chromosome);                    
                    line = annotationReader.readLine();
                } //while             
                annotationReader.close();
            }
            fsf = new FastaSequenceFile(new File(queryFastaFile), true);          

            writer = new FileWriter(resultFile);
            gff = new FileWriter(new File(new File(resultFile).getParent(),"subjects.gff"));
            //writer.write("Query Contig\tQuery Start\tQuery End\tSubject Contig\tSubject Start\tSubject End\tOrientation\n");
            blast_line = reader.readLine();
            step = 1;
            while (blast_line != null && !blast_line.isEmpty()){
                String[] arr = blast_line.split("\\s");
                step = 2;
                String qseqid = arr[0].trim();
                int qstart = Integer.parseInt(arr[1]);
                int qend = Integer.parseInt(arr[2]);
                String sseqid = arr[3].trim();
                int sstart = Integer.parseInt(arr[4]);
                int send = Integer.parseInt(arr[5]);
                int alignment_len = Integer.parseInt(arr[6]); // not used at this stage
                double percent = Double.parseDouble(arr[7]); // normal formulation
                int mismatches = Integer.parseInt(arr[8]);
                int gapopen = Integer.parseInt(arr[9]);
                //double evalue = Double.parseDouble(arr[12]);
                double no_gap_percent = (double)(alignment_len - mismatches - gapopen)/alignment_len; // special formulation
                no_gap_percent *= 100;
                String chromosome = null;
                if (popseqMap.size() > 0)
                    chromosome = popseqMap.get(sseqid);
                else 
                    chromosome = Helper.getWheatChromosome(sseqid); // if unspecified it returns U
                if (chromosome==null)
                    System.out.println(sseqid + " contig was not found in popseq map.");
//if (qseqid.equals("contig00003"))             
step = 3;
                if(!qseqid.equals(last_query_contig)){
                   if (last_query_contig != null){
                       int query_size_A = getTotalQuerySize(A_genome);
                       total_covered_query_inserted += query_size_A;
                       double coverage_A = (double)query_size_A/last_query_contig_size;
                       writer.write(String.format("A:\t%s\t%f\t%s\t%d\t%d\n", last_query_contig, coverage_A, getMajoritySubjectContig(A_genome), notInsertedHSP_A, A_genome.size()));
                       //writer.write (last_query_contig + " A: Covered:" + (double)query_size_A/last_query_contig_size + ", not inserted HSP:" + notInsertedHSP_A + "\n");
                       for(BlastHSP aHSP : A_genome){
                          writer.write(String.format("%s\t%d\t%d\t%s\t%d\t%d\t%f\n", aHSP.q_contig, aHSP.q_start, aHSP.q_end , aHSP.s_contig, aHSP.s_start, aHSP.s_end, aHSP.identity_percent)); 
                          gff.write(aHSP.toGFFLine("")+"\n");
                       }
                       int query_size_B = getTotalQuerySize(B_genome);
                       total_covered_query_inserted += query_size_B;
                       double coverage_B = (double)query_size_B/last_query_contig_size;
                       writer.write(String.format("B:\t%s\t%f\t%s\t%d\t%d\n", last_query_contig, coverage_B, getMajoritySubjectContig(B_genome), notInsertedHSP_B, B_genome.size()));
                       //writer.write (last_query_contig + " B: Covered:" + (double)query_size_B/last_query_contig_size + ", not inserted HSP:" + notInsertedHSP_B + "\n");
                       for(BlastHSP aHSP : B_genome){
                          writer.write(String.format("%s\t%d\t%d\t%s\t%d\t%d\t%f\n", aHSP.q_contig, aHSP.q_start, aHSP.q_end , aHSP.s_contig, aHSP.s_start, aHSP.s_end, aHSP.identity_percent));  
                          gff.write(aHSP.toGFFLine("")+"\n");
                       }
                       int query_size_D = getTotalQuerySize(D_genome);
                       total_covered_query_inserted += query_size_D;
                       double coverage_D= (double)query_size_D/last_query_contig_size;
                       writer.write(String.format("D:\t%s\t%f\t%s\t%d\t%d\n", last_query_contig, coverage_D, getMajoritySubjectContig(D_genome), notInsertedHSP_D, D_genome.size()));                       
                       //writer.write (last_query_contig + " D: Covered:" + (double)query_size_D/last_query_contig_size + ", not inserted HSP:" + notInsertedHSP_D + "\n");
                       for(BlastHSP aHSP : D_genome){
                          writer.write(String.format("%s\t%d\t%d\t%s\t%d\t%d\t%f\n", aHSP.q_contig, aHSP.q_start, aHSP.q_end , aHSP.s_contig, aHSP.s_start, aHSP.s_end, aHSP.identity_percent));   
                          gff.write(aHSP.toGFFLine("")+"\n");
                       }
                       writer.write("#\n"); //this is separator
                       
                       
                   } 

                    //look for qseqid in the query Fasta File contigs (we assume the order of contigs is the same in both files)
                    boolean found = false;                    
                    ReferenceSequence seq = fsf.nextSequence();
                    while(seq != null){
                        // when the name of contig after > is multipart, in blast result only the first part is stored
                        //System.out.println(seq.getName());
                        if (seq.getName().equals(qseqid)){
                            found = true;  
                            last_query_contig_size = seq.length();
                            sum_query_contig_size += last_query_contig_size;
                            break;
                        }      
                        // if qcontig is not found, just output empty lines into out put
                        writer.write(String.format("A:\t%s\t%f\t%s\t%d\t%d\n", seq.getName(), 0.0, null, 0, 0));
                        writer.write(String.format("B:\t%s\t%f\t%s\t%d\t%d\n", seq.getName(), 0.0, null, 0, 0));
                        writer.write(String.format("D:\t%s\t%f\t%s\t%d\t%d\n", seq.getName(), 0.0, null, 0, 0));
                        writer.write("#\n"); //this is separator
                        
                        seq = fsf.nextSequence();
                    }
                    if(!found)
                        throw new Exception("Error: Order of blast result and original query fasta file are not the same. qseqid:" + qseqid);

                   A_genome.clear();
                   B_genome.clear();
                   D_genome.clear();
                   notInsertedHSP_A = 0;
                   notInsertedHSP_B = 0;
                   notInsertedHSP_D = 0;
                   last_query_contig = qseqid;
                   
                }
                step = 4;
                //if (evalue<max_evalue && no_gap_percent>min_identity )
                //    total_covered_query_nogap += qend-qstart+1;
               // if (alignment_len>=max_evalue && no_gap_percent>min_identity ){
                if (alignment_len>=max_evalue && percent>min_identity ){
                    // we put each hit if its query is distinct or just 10bp overlap with existing queries and there is no overlap with their subject
                    // we can not use orientation, because we are not sure the strand of the exome capture contigs
                    total_covered_query += qend-qstart+1;
                    BlastHSP hsp = new BlastHSP(qseqid,qstart,qend,sseqid,sstart,send,alignment_len,percent,0);
                    
                    if (chromosome.matches("\\d+A")){
                        if (!insertIntoBlastHSPList(A_genome, hsp))
                            notInsertedHSP_A++;
                    }else if (chromosome.matches("\\d+B")){
                        if (!insertIntoBlastHSPList(B_genome, hsp))
                            notInsertedHSP_B++;
                    }else if (chromosome.matches("\\d+D")){
                        if (!insertIntoBlastHSPList(D_genome, hsp))
                            notInsertedHSP_D++;
                    }else{
                        //here means chromosome is unspecified- so we try to insert it into A list if unsuccessful, then insert it into B and so forth
                        if (insertIntoBlastHSPList(A_genome, hsp)){
                        }else if (insertIntoBlastHSPList(B_genome, hsp)){
                        }else{
                           insertIntoBlastHSPList(D_genome, hsp); 
                        }
                    }
                    
                }
                step = 5;
                blast_line = reader.readLine();
                
            }//while loop
            
            // The last query contig has not have a chance, so we write it now
           if (last_query_contig != null){
               int query_size_A = getTotalQuerySize(A_genome);
               total_covered_query_inserted += query_size_A;
               double coverage_A = (double)query_size_A/last_query_contig_size;
               writer.write(String.format("A:\t%s\t%f\t%s\t%d\t%d\n", last_query_contig, coverage_A, getMajoritySubjectContig(A_genome), notInsertedHSP_A, A_genome.size()));
               //writer.write (last_query_contig + " A: Covered:" + (double)query_size_A/last_query_contig_size + ", not inserted HSP:" + notInsertedHSP_A + "\n");
               for(BlastHSP aHSP : A_genome){
                  writer.write(String.format("%s\t%d\t%d\t%s\t%d\t%d\t%f\n", aHSP.q_contig, aHSP.q_start, aHSP.q_end , aHSP.s_contig, aHSP.s_start, aHSP.s_end, aHSP.identity_percent)); 
                  gff.write(aHSP.toGFFLine("")+"\n");
               }
               int query_size_B = getTotalQuerySize(B_genome);
               total_covered_query_inserted += query_size_B;
               double coverage_B = (double)query_size_B/last_query_contig_size;
               writer.write(String.format("B:\t%s\t%f\t%s\t%d\t%d\n", last_query_contig, coverage_B, getMajoritySubjectContig(B_genome), notInsertedHSP_B, B_genome.size()));
               //writer.write (last_query_contig + " B: Covered:" + (double)query_size_B/last_query_contig_size + ", not inserted HSP:" + notInsertedHSP_B + "\n");
               for(BlastHSP aHSP : B_genome){
                  writer.write(String.format("%s\t%d\t%d\t%s\t%d\t%d\t%f\n", aHSP.q_contig, aHSP.q_start, aHSP.q_end , aHSP.s_contig, aHSP.s_start, aHSP.s_end, aHSP.identity_percent));  
                  gff.write(aHSP.toGFFLine("")+"\n");
               }
               int query_size_D = getTotalQuerySize(D_genome);
               total_covered_query_inserted += query_size_D;
               double coverage_D= (double)query_size_D/last_query_contig_size;
               writer.write(String.format("D:\t%s\t%f\t%s\t%d\t%d\n", last_query_contig, coverage_D, getMajoritySubjectContig(D_genome), notInsertedHSP_D, D_genome.size()));                       
               //writer.write (last_query_contig + " D: Covered:" + (double)query_size_D/last_query_contig_size + ", not inserted HSP:" + notInsertedHSP_D + "\n");
               for(BlastHSP aHSP : D_genome){
                  writer.write(String.format("%s\t%d\t%d\t%s\t%d\t%d\t%f\n", aHSP.q_contig, aHSP.q_start, aHSP.q_end , aHSP.s_contig, aHSP.s_start, aHSP.s_end, aHSP.identity_percent));   
                  gff.write(aHSP.toGFFLine("")+"\n");
               }
               writer.write("#\n"); //this is separator
           } 
            
            System.out.println("Total length of Query Contigs: " + sum_query_contig_size);
            System.out.println("Total Query Regions qualified for insertion: " + total_covered_query);
            //System.out.println("Total Query Regions qualified for insertion (nogap): " + total_covered_query_nogap);
            System.out.println("Total Query Regions  inserted: " + total_covered_query_inserted);
            
           
            
        } catch (Exception ex) {System.out.println("findAssemblyHomoeleogRegions" + ex.getMessage() + " at step "+step);}
        finally{
            try {
                writer.close();
                gff.close();
            } catch (Exception ex) {}
        } 
        
        
    }

   public void merge2AssemblyHomoeleogRegions(String file1, String file2, String queryFastaFile, String queryFastaIndexFile){
       // It is supposed that file1 and file2 have the same number of records
       //chromosome data is supposed to be extracted from file1
       //Output is called regions.gff
       //it also outputs contigs that were not found in a fasta file.
        

   
       HomeologRecordReader reader1 = null;
       HomeologRecordReader reader2 = null;
       FastaSequenceIndex fsi = null;
       IndexedFastaSequenceFile ifsi = null;
       
        FileWriter gff = null;
        FileWriter notFoundContig = null;
   
        int rec1_count = 0;
        int rec2_count = 0;
        int zero_count = 0;
        long zero_len = 0;
        
        try{ 
 
            fsi = new FastaSequenceIndex(new File(queryFastaIndexFile));
            ifsi =new IndexedFastaSequenceFile(new File(queryFastaFile) , fsi);

            
            gff = new FileWriter(new File(new File(file1).getParent(),"regions.gff"));
            notFoundContig = new FileWriter(new File(new File(file1).getParent(),"notFoundContigs.fa"));
            reader1 = new HomeologRecordReader(file1); //svevo
            reader2 = new HomeologRecordReader(file2); //kronos
            HomeologRecord rec1 = reader1.getNext();
            HomeologRecord rec2 = reader2.getNext();
            
            while (rec1 != null && rec2 != null){
                //rec1 that is svevo in this example is used for chromosome names
                // if one chromosome is U and the other one is non-U then replace U with non-U 
                if (rec1.A_chromosome.equals("U")&& !rec1.B_chromosome.equals("U"))
                    rec1.A_chromosome = rec1.B_chromosome.replace('B', 'A');
                else if (!rec1.A_chromosome.equals("U")&& rec1.B_chromosome.equals("U"))
                    rec1.B_chromosome = rec1.A_chromosome.replace('A', 'B');
                // see total coverage of A and B 
                if((rec1.A_coverage+rec1.B_coverage>0.1) || (rec2.A_coverage+rec2.B_coverage>0.1)){
                    if ( ((rec1.A_coverage+rec1.B_coverage)>=(rec2.A_coverage+rec2.B_coverage))  ){
                        for(BlastHSP aHSP : rec1.A_hsps){
                          gff.write(aHSP.toGFFLine(rec1.A_chromosome)+"\n");
                        }
                        for(BlastHSP aHSP : rec1.B_hsps){
                          gff.write(aHSP.toGFFLine(rec1.B_chromosome)+"\n");
                        }
                        rec1_count++;
                    }else {
                        // rec2 is written into gff file, but chromosome is still taken from rec1 (Svevo has chromosome number in its contig name)
                        for(BlastHSP aHSP : rec2.A_hsps){
                          gff.write(aHSP.toGFFLine(rec1.A_chromosome)+"\n");
                        }
                        for(BlastHSP aHSP : rec2.B_hsps){
                          gff.write(aHSP.toGFFLine(rec1.B_chromosome)+"\n");
                        }                    
                        rec2_count++;
                    }
                }else{
                    zero_count++;
                    int qcontig_size = (int)fsi.getContigSize(rec1.qcontig);
                    zero_len += qcontig_size;
                    notFoundContig.write( String.format(">%s\n%s\n", rec1.qcontig, ifsi.getSequence(rec1.qcontig).getBaseString()) );
                }
                rec1 = reader1.getNext();
                rec2 = reader2.getNext();
            }//while loop
            if((rec1==null) ^ (rec2==null) )
                throw new Exception("Number of records in 2 files are not equal.");
            System.out.println("Total regions selected from File 1 : " + rec1_count);
            System.out.println("Total regions selected from File 2 :" + rec2_count);
            System.out.println("Total regions not covered at all:" + zero_count);
            System.out.println("Total regions not covered size (bp):" + zero_len);
        } catch (Exception ex) {System.out.println("merge2AssemblyHomoeleogRegions:" + ex.getMessage());}
        finally{
            try {
                gff.close();
                notFoundContig.close();
                reader1.close();
                reader2.close();
            } catch (Exception ex) {}
        } 
        
        
    }
    
    
    
    private void insertIntoBlastHSPList2(LinkedList<BlastHSP> list, BlastHSP aHSP){
        // we want to maintain order of HSPs based on the subject start (if HSP is +/- based on subject_end)
        //List contains ordered HSP that their orientation has been normalized.
        
        aHSP.normalizeOrientation(); // if orientation is +/- then swap sstart and send
        
        boolean added = false;
        for(int i=0;i<list.size();i++){
            BlastHSP hsp = list.get(i);
            if(aHSP.s_start < hsp.s_start){
                list.add(i, aHSP);
                added = true;
                break;
            }
        }
        
        if (!added)// then it means that aHSP comese after all the elements of the list
            list.add(aHSP);
    }
    
    public static boolean insertIntoBlastHSPList(LinkedList<BlastHSP> list, BlastHSP aHSP){
        // finds proper position in the list for hsp and insert or replace it in the list
        // Here we dont check for subject overlap. Just query overlaps. For more accurate result we need to check for subject overlap as well.
        //However, because blast automaticcally sort results based on query contig and subject contigs (for each subject contig sort based on the score), there is no need to check for subject contig overlaps
        int i=0;
        for(; i<list.size(); i++){
            BlastHSP hsp = list.get(i);
            if (i==0){
                if (aHSP.comesBefore(hsp, 12)){
                    list.add(0, aHSP);
                    return true;
                }
            }else{
                // i>0, then aHSP should be between element i-1 and i to be able to be inserted into i
                if (aHSP.comesBefore(hsp, 12) && aHSP.comesAfter(list.get(i-1), 12)){
                    list.add(i, aHSP);
                    return true;
                }
            }
        }
        // we get here for 3 reasons
        //1: list is empty
        if (list.size()==0){
            list.add(aHSP);
            return true;
        }
        //2: aHSP comes after last element of the list
        if (aHSP.comesAfter(list.get(list.size()-1), 12)){
            list.add(aHSP);// add it 
            return true;
        }
        //3: aHSP has significant overlap with existing HSPs in the list
        
        return false;
    }
    private int getTotalQuerySize(LinkedList<BlastHSP> list){
        int total = 0;
        for(BlastHSP hsp : list)
            total += hsp.q_end-hsp.q_start + 1 ;
        
        return total;
    }
    private String getMajoritySubjectContig(LinkedList<BlastHSP> list){
        // which subject contig has been repeated the most
        Map<String, Integer> counts = new HashMap();
        if (list.size() == 1)
            return list.get(0).s_contig;
              
        MapUtil util = new MapUtil();
        for(BlastHSP hsp : list)
            util.increment(counts, hsp.s_contig);
            
         String most_frequent_contig = null;
         int max_count = 0;
         for(String contig : counts.keySet()){
             int count = counts.get(contig);
             if(count>max_count || (count==max_count && !Helper.getWheatChromosome(contig).equals("U"))){
                 max_count = count;
                 most_frequent_contig = contig;
             }
         } //for
        
        return most_frequent_contig;
    }
    
    public void findSplittingPoints(String resultFile, double min_identity, int min_alignment_len){
        // This task process a blast result file (blast multi-contig query to multi-contig subject)finds query contigs and their cordinates where 2 separate segments of query contig (with significant space between them) aligned around one point in the reference
        // Such a point on subject contig can potentially be the site of a big insertion point, and the sequence between end of segment1 and start of segment 2 of query contig is inserted sequence
        // Obviously such a point on subject contig can not be exact so we assign 10 bp thereshold for that.
        // Note: these 2 segments in query contig must satisfy min_identity and min_alignment_len and also must align to subject with the same orientation (+/+, +/-)
        String blast_line = null;
        String last_query_contig = null; // blast result is sorted by query contig and then subject contig
        String last_subject_contig = null; // blast result is sorted by subject contig
        
        FileWriter writer = null;
        //TreeMap<Integer, Integer> query_posMap_plus = new TreeMap(); // this Map holds start and end positions of each hit in the query, and sorted based on its key (qstart)
        //TreeMap<Integer, Integer> query_posMap_minus = new TreeMap(); // this Map holds start and end positions of each hit in the query, and sorted based on its key (qstart)
        
        //TreeMap<Integer, Integer> subject_posMap_plus = new TreeMap(); // this Map holds start and end positions of each positive hit in the last_contig, and sorted based on its key (sstart)
        //TreeMap<Integer, Integer> subject_posMap_minus = new TreeMap(); // this Map holds start and end positions of each negative hit in the last_contig, and sorted based on its key (sstart)

        LinkedList<BlastHSP> positive_list = new LinkedList(); // for each query and subject contig stores +/+ HSPs sorted based on subject_start
        LinkedList<BlastHSP> negative_list = new LinkedList();//for each query and subject contigs stores +/- HSPs  sorted based on subject_end
        int step = 0;
        try{
           writer = new FileWriter(resultFile);        
          // writer.write("Query Contig\tQuery Start\tQuery End\tQuery Contig Coverage\tSubject Contig\tSubject Start\tSubject End\tOrientation\tSubject Region Covered\n");
            blast_line = reader.readLine();
            step = 1;
            while (blast_line != null && !blast_line.isEmpty()){
                String[] arr = blast_line.split("\\s");
                String sseqid = arr[3].trim();
                double percent = Double.parseDouble(arr[7]);
                int alignment_len = Integer.parseInt(arr[6]); 
                double evalue = Double.parseDouble(arr[12]);
                String qseqid = arr[0].trim();
                int qstart = Integer.parseInt(arr[1]);
                int qend = Integer.parseInt(arr[2]);
                int sstart = Integer.parseInt(arr[4]);
                int send = Integer.parseInt(arr[5]);
                if(!qseqid.equals(last_query_contig) || !sseqid.equals(last_subject_contig)){
                    // process  last contig posMap
                    if (last_subject_contig != null){//this if makes sure we are not processing the first line of the blast file
                        // process data related to last_query_contig and last_subject_contig
                        // now we investigate 2 LinkedList (sorted HSPs) and we will report any 2 `connected HSPs                       
                        for (int i=0; i<positive_list.size()-1; i++){
                            // see if i and i+1 are connected
                            if (positive_list.get(i).subjetConnected(positive_list.get(i+1), 10)){
                                // report HSP
                                writer.write(positive_list.get(i).toString()+ " -> " + positive_list.get(i+1).toString() + "\n");
                            }
                        }
                        for (int i=0; i<negative_list.size()-1; i++){
                            // see if i and i+1 are connected
                            if (negative_list.get(i).subjetConnected(negative_list.get(i+1), 10)){
                                // report HSP
                                writer.write(negative_list.get(i).toString()+ " -> " + negative_list.get(i+1).toString() + "\n");
                            }
                        }

                    } // if we are not at the first line of blast result
                    
                    // empty the Map for new contig
                    //query_posMap_plus.clear();
                    //query_posMap_minus.clear();
                    //subject_posMap_plus.clear();
                    //subject_posMap_minus.clear();
                    positive_list.clear();
                    negative_list.clear();
                    last_subject_contig = sseqid;
                    last_query_contig = qseqid;
                    step = 8;
                }// if query contig or subject contig has changed
                
                //Note: By sorting HSPs based on subject_start and then check 2 consecutive HSP for connection, we dont guarantee to find all connected HSP but in simple case we find them
                if (alignment_len>=min_alignment_len && percent>=min_identity ){   
                    BlastHSP hsp = new BlastHSP(qseqid,qstart,qend,sseqid,sstart,send,alignment_len,percent,evalue);
                   if(sstart<send){
                       // plus-plus 
                       insertIntoBlastHSPList2(positive_list, hsp);
                   }else{
                       //plus-minus
                       insertIntoBlastHSPList2(negative_list, hsp);
                   }
                }
                 step = 9;
                blast_line = reader.readLine();
                
            }//while there is another line in blast result file

        } catch (Exception ex) {System.out.println("findSplittingPoints" + ex.getMessage() + " at step "+step);}
        finally{
            try {
                writer.close();
            } catch (Exception ex) {}
        } 
    }
    public void findHomologRegions(String queryFastaFile,  String resultFile,   double min_identity, double max_evalue){
        // This method is similar to findHomologSubject, but differes in the query contig, this one accepts a query fasta file with many contigs in it
        // This method assumes order of query contigs in the blast result file is as the same for query fasta file, therefore we wont need the index of query fasta file.
        // The order of blast hits is like this: first based on query contig, for each query contig order by subject contig, for each subject contig hits are ordered based on the descending score 
        // The aim of this method is to determine if a region in a subject contig is a true homolg of a query contig 
        
        // We set some rules,(1) subject HSP regions must have the same orientation (Actually we separate positive from negative strand subject HSPs) and (2) from start to the end of query HSPs must account for at least 95% of the length of the query contig. (3)If we measure from the start to the end of the query contig, it must also be equal to the subject region with 10% (this percentage can be parameterized) variation in length to account for In/Dels.
        // for each query contig and subject contig pair we produce a single line (if it follows above rules) 
        // we also report percentage of whole query contig that is covered and also percentage of subject region (not subject contig, because we dont expect whole query contig will map to whole subject contig)
            
        FastaSequenceFile fsf = null;  // we sequentially read query fasta file so we dont need its index
        String blast_line = null;
        String last_query_contig = null; // blast result is sorted by query contig and then subject contig
        String last_subject_contig = null; // blast result is sorted by subject contig
        int last_query_contig_size = 0;
        FileWriter writer = null;
        TreeMap<Integer, Integer> query_posMap_plus = new TreeMap(); // this Map holds start and end positions of each hit in the query, and sorted based on its key (qstart)
        TreeMap<Integer, Integer> query_posMap_minus = new TreeMap(); // this Map holds start and end positions of each hit in the query, and sorted based on its key (qstart)
        
        TreeMap<Integer, Integer> subject_posMap_plus = new TreeMap(); // this Map holds start and end positions of each positive hit in the last_contig, and sorted based on its key (sstart)
        TreeMap<Integer, Integer> subject_posMap_minus = new TreeMap(); // this Map holds start and end positions of each negative hit in the last_contig, and sorted based on its key (sstart)
        int step = 0;
        try{ 
            
            
            fsf = new FastaSequenceFile(new File(queryFastaFile), true);          

            writer = new FileWriter(resultFile);
            writer.write("Query Contig\tQuery Start\tQuery End\tQuery Contig Coverage\tSubject Contig\tSubject Start\tSubject End\tOrientation\tSubject Region Covered\n");
            blast_line = reader.readLine();
            step = 1;
            while (blast_line != null && !blast_line.isEmpty()){
                String[] arr = blast_line.split("\\s");
                String sseqid = arr[3].trim();
                double percent = Double.parseDouble(arr[7]);
                int alignment_len = Integer.parseInt(arr[6]); // not used at this stage
                double evalue = Double.parseDouble(arr[12]);
                String qseqid = arr[0].trim();
                int qstart = Integer.parseInt(arr[1]);
                int qend = Integer.parseInt(arr[2]);
                int sstart = Integer.parseInt(arr[4]);
                int send = Integer.parseInt(arr[5]);
                //if(sseqid.equals("3AS_3321827"))
                 //           alignment_len = 0;
                
                if(!qseqid.equals(last_query_contig) || !sseqid.equals(last_subject_contig)){
                    // process  last contig posMap
                    if (last_subject_contig != null){//this if makes sure we are not processing the first line of the blast file
                        // process data related to last_query_contig and last_subject_contig
                        // now we investigate 3 posMaps
                        MapUtil util = new MapUtil();
                        int qpluslen = 0;
                        int qminuslen = 0;
                        int biggest_plus_len = 0;
                        int biggest_minus_len = 0;
                        int subject_plus_len = 0;
                        int subject_minus_len = 0;
                        int plus_difference = -1;  // it is very important to be -1, not to be confused by 0
                        int minus_difference = -1;  // it is very important to be -1, not to be confused by 0
                        int query_plus_covered = 0;  // qpluslen is important, not query_plus_covered
                        int query_minus_covered = 0;  // qminuslen is important, not query_minus_covered
                        int subject_plus_covered = 0;
                        int subject_minus_covered = 0;
                       
                        step = 3;
                       
                        if (query_posMap_plus.size() > 0){
                            qpluslen = util.right_end(query_posMap_plus) - query_posMap_plus.firstKey()+1;
                            query_plus_covered = qpluslen - util.notCovered(query_posMap_plus);                            
                        }
                        if (query_posMap_minus.size() > 0){
                            qminuslen = util.right_end(query_posMap_minus) - query_posMap_minus.firstKey()+1;
                            query_minus_covered = qminuslen - util.notCovered(query_posMap_minus);                            
                        }

                        step = 4;
                        if (subject_posMap_plus.size()>0){
                            subject_plus_len = util.right_end(subject_posMap_plus) - subject_posMap_plus.firstKey() + 1;
                            plus_difference = Math.abs(qpluslen-subject_plus_len);
                            subject_plus_covered = subject_plus_len - util.notCovered(subject_posMap_plus);
                        }
                        step = 5;
                        if (subject_posMap_minus.size()>0){
                            subject_minus_len = util.right_end(subject_posMap_minus) - subject_posMap_minus.firstKey() + 1;
                            minus_difference = Math.abs(qminuslen-subject_minus_len);
                            subject_minus_covered = subject_minus_len - util.notCovered(subject_posMap_minus);
                        }
                        step = 6;
                        //System.out.println(last_subject_contig+":"+subject_conig_size+":"+qlen);
                        biggest_plus_len = Math.max(qpluslen, subject_plus_len);
                        biggest_minus_len = Math.max(qminuslen, subject_minus_len);
                        if(plus_difference>=0  && plus_difference < 0.1*biggest_plus_len ){
                            //plus-plus hits are selected
                            step = 61;
                            writer.write(String.format("%s\t%d\t%d\t%.2f\t%s\t%d\t%d\t+\t%.2f\n", last_query_contig, query_posMap_plus.firstKey(), util.right_end(query_posMap_plus) , (float)qpluslen/last_query_contig_size ,last_subject_contig, subject_posMap_plus.firstKey(), util.right_end(subject_posMap_plus), (float)subject_plus_covered/subject_plus_len));
                            
                        }
                        if (minus_difference>=0   && minus_difference < 0.1*biggest_minus_len  ){                            
                            //plus-minus hits are selected
                            step = 62;
                            writer.write(String.format("%s\t%d\t%d\t%.2f\t%s\t%d\t%d\t-\t%.2f\n", last_query_contig, query_posMap_minus.firstKey(), util.right_end(query_posMap_minus) , (float)qminuslen/last_query_contig_size , last_subject_contig, subject_posMap_minus.firstKey(), util.right_end(subject_posMap_minus), (float)subject_minus_covered/subject_minus_len));
                        }
                        // else subject contig is not selected as homolog candidate
                        step = 7;
                    } // if we are not at the first line of blast result
                    
                    if(!qseqid.equals(last_query_contig)){
                        //look for qseqid in the query Fasta File contigs (we assume the order of contigs is the same in both files)
                        boolean found = false;                    
                        ReferenceSequence seq = fsf.nextSequence();
                        while(seq != null){
                            // when the name of contig after > is multipart, in blast result only the first part is stored
                            //System.out.println(seq.getName());
                            if (seq.getName().equals(qseqid)){
                                found = true;  
                                last_query_contig_size = seq.length();
                                break;
                            }                        
                            seq = fsf.nextSequence();
                        }
                        if(!found)
                            throw new Exception("Error: Order of blast result and original query fasta file are not the same. qseqid:" + qseqid);
                    }
                    // empty the Map for new contig
                    query_posMap_plus.clear();
                    query_posMap_minus.clear();
                    subject_posMap_plus.clear();
                    subject_posMap_minus.clear();
                    last_subject_contig = sseqid;
                    last_query_contig = qseqid;
                    step = 8;
                }// if query contig or subject contig has changed
                 
                if (evalue<max_evalue && percent>min_identity ){                    
                   if(sstart<send){
                       // plus-plus 
                       subject_posMap_plus.put(sstart, send);
                       query_posMap_plus.put(qstart, qend);
                   }else{
                       //plus-minus
                       subject_posMap_minus.put(send, sstart);
                       query_posMap_minus.put(qstart, qend);
                   }
                }
                 step = 9;
                blast_line = reader.readLine();
                
            }//while loop
             // the last subject contig still needs to be processed
            MapUtil util = new MapUtil();
            int qpluslen = 0;
            int qminuslen = 0;
            int biggest_plus_len = 0;
            int biggest_minus_len = 0;
            int subject_plus_len = 0;
            int subject_minus_len = 0;
            int plus_difference = -1;  // it is very important to be -1, not to be confused by 0
            int minus_difference = -1;  // it is very important to be -1, not to be confused by 0
            int query_plus_covered = 0; // qpluslen is important, not query_plus_covered
            int query_minus_covered = 0;  // qminuslen is important, not query_minus_covered
            int subject_plus_covered = 0;
            int subject_minus_covered = 0;
            
            if (query_posMap_plus.size() > 0){
                qpluslen = util.right_end(query_posMap_plus) - query_posMap_plus.firstKey()+1;
                query_plus_covered = qpluslen - util.notCovered(query_posMap_plus);                            
            }
            if (query_posMap_minus.size() > 0){
                qminuslen = util.right_end(query_posMap_minus) - query_posMap_minus.firstKey()+1;
                query_minus_covered = qminuslen - util.notCovered(query_posMap_minus);                            
            }

            step = 4;
            if (subject_posMap_plus.size()>0){
                subject_plus_len = util.right_end(subject_posMap_plus) - subject_posMap_plus.firstKey() + 1;
                plus_difference = Math.abs(qpluslen-subject_plus_len);
                subject_plus_covered = subject_plus_len - util.notCovered(subject_posMap_plus);
            }
            step = 5;
            if (subject_posMap_minus.size()>0){
                subject_minus_len = util.right_end(subject_posMap_minus) - subject_posMap_minus.firstKey() + 1;
                minus_difference = Math.abs(qminuslen-subject_minus_len);
                subject_minus_covered = subject_minus_len - util.notCovered(subject_posMap_minus);
            }
            step = 6;
            //System.out.println(last_subject_contig+":"+subject_conig_size+":"+qlen);
            biggest_plus_len = Math.max(qpluslen, subject_plus_len);
            biggest_minus_len = Math.max(qminuslen, subject_minus_len);
            if(plus_difference>=0 && plus_difference < 0.1*biggest_plus_len ){
                //plus-plus hits are selected
                step = 61;
                writer.write(String.format("%s\t%d\t%d\t%.2f\t%s\t%d\t%d\t+\t%.2f\n", last_query_contig, query_posMap_plus.firstKey(), util.right_end(query_posMap_plus) , (float)qpluslen/last_query_contig_size ,last_subject_contig, subject_posMap_plus.firstKey(), util.right_end(subject_posMap_plus), (float)subject_plus_covered/subject_plus_len));

            }
            if (minus_difference>=0   && minus_difference < 0.1*biggest_minus_len  ){                            
                //plus-minus hits are selected
                step = 62;
                writer.write(String.format("%s\t%d\t%d\t%.2f\t%s\t%d\t%d\t-\t%.2f\n", last_query_contig, query_posMap_minus.firstKey(), util.right_end(query_posMap_minus) , (float)qminuslen/last_query_contig_size , last_subject_contig, subject_posMap_minus.firstKey(), util.right_end(subject_posMap_minus), (float)subject_minus_covered/subject_minus_len));
            }

             
        } catch (Exception ex) {System.out.println("findHomologRegions" + ex.getMessage() + " at step "+step);}
        finally{
            try {
                writer.close();
            } catch (Exception ex) {}
        } 
        
    }

    public void getQueryRegionsCovered(String queryFastaFile, double min_identity_percent, int min_alignment_length){
        // This method reads  fasta query file and for each contig in the this file works out how much of it is covered by blast result
        // In other words, this method uses the same query fasta file that was used in blast command, and estimates how much of that is covered by the blast
        // blast result file is of the format 6  qseqid qstart qend sseqid sstart send length pident mismatch gapopen gaps
        // for a blast result entry min_identity and min_alignment_length is criteria to be selected
        // we assume blast result sorted by 'qseqid' but not with 'qstart' and 'qend'
        // The strategy is we first sort all the blast entries in a contig by qstart
        // we assume order of contigs showed up in the query fasta is the same as order of qseqid in blast result file., otherwise this method will raise an Exception
        
            FastaSequenceFile fsf = null;  // we sequentially read query fasta file so we dont need its index
            String last_qseqid = "";
            long last_contig_len = 0;  // length of the contig of last_qseqid
            Map<Integer, Integer> posMap = new TreeMap(); // this Map holds start and end positions of each query hit in the last_contig, and sorted based on its key (qstart)
            //we report following 2 variables
            long total_contig_len = 0; // sum of length of all the contigs in query fasta file
            long total_not_covered = 0; // if a whole contig or some parts of query fasta file that is not mapped 

        try {
            long pbIndex = 0;
            double lastpercent = 0.0;
            long fileSz = new File(this.fileName).length();

            fsf = new FastaSequenceFile(new File(queryFastaFile), true);
            String blast_line = reader.readLine();
            while (blast_line != null && !blast_line.isEmpty()){
                pbIndex += blast_line.getBytes().length; // we add number of bytes read
                double newpercent = ((double)pbIndex/fileSz)*100;
                if (newpercent-lastpercent>=1){
                    lastpercent = newpercent;
                    System.out.println("progress %"+(int)lastpercent);
                }

                
                String[] arr = blast_line.split("\\s");
                String qseqid = arr[0].trim();
                int qstart = Integer.parseInt(arr[1]);
                int qend = Integer.parseInt(arr[2]);
                double percent = Double.parseDouble(arr[7]);
                int alignment_len = Integer.parseInt(arr[6]);
                
                if(!qseqid.equals(last_qseqid)){
                    // process  last contig posMap
                    if (!last_qseqid.isEmpty()){
                        // posMap is ordered based on the starting position
                        Iterator<Map.Entry<Integer, Integer>> it = posMap.entrySet().iterator();
                        int last_end = 0; // end position of the last entry
                        while(it.hasNext()){
                            Map.Entry<Integer,Integer> entry = it.next();
                            int start = entry.getKey();
                            int end = entry.getValue();
                            if (start > last_end)
                                total_not_covered += start-last_end-1;
                            if(end>last_end)
                                last_end=end;
                        }// while in posMap
                        // we add remainder of the contig from last_end until contig_length
                        if (last_end < last_contig_len)
                            total_not_covered += last_contig_len - last_end;

                    }// if last_qseqid is not empty
                    
                    //look for qseqid in the query Fasta File contigs (we assume the order of contigs is the same in both files)
                    boolean found = false;
                    
                    ReferenceSequence seq = fsf.nextSequence();
                    while(seq != null){
                        total_contig_len += seq.length();
                        if (seq.getName().equals(qseqid)){
                            found = true;  
                            last_contig_len = seq.length();
                            break;
                        }else 
                            total_not_covered += seq.length(); // this whole contig is not mapped, so all o its length is not covered
                        
                        seq = fsf.nextSequence();
                    }
                    if(!found)
                        throw new Exception("Error: Order of blast result and original query fasta file are not the same.");
                    // empty the Map for new contig
                    posMap.clear();
                    last_qseqid = qseqid;
                }
                // add the query blast region to the Map if meets the criteria of min_percentage and min_alignment_length
                if (percent >= min_identity_percent  && alignment_len>=min_alignment_length)
                    posMap.put(qstart, qend);
                blast_line = reader.readLine();
            } //while
            // the last contig still needs to be processed
            if(posMap.size() > 0){
                Iterator<Map.Entry<Integer, Integer>> it = posMap.entrySet().iterator();
                int last_end = 0; // end position of the last entry
                while(it.hasNext()){
                    Map.Entry<Integer,Integer> entry = it.next();
                    int start = entry.getKey();
                    int end = entry.getValue();
                    if (start > last_end)
                        total_not_covered += start-last_end-1;
                    if(end>last_end)
                        last_end=end;
                }//while
                // we add remainder of the contig from last_end until contig_length
                if (last_end < last_contig_len)
                    total_not_covered += last_contig_len - last_end;

            }//if last contig has entries
            
            System.out.println("Total contig length: " + total_contig_len);
            System.out.println("Total not covered: " + total_not_covered);
            
            
        } catch (Exception ex) {System.out.println("getQueryRegionsCovered::" + ex.getMessage());}
        finally{
            fsf.close();
        }    
            
        
    }
    
    public void writeQueryRegionsNotCovered(String queryFastaFile, double min_identity_percent, int min_contig_length, String outFasta){
        //Very similar to getQueryRegionsCovered, written on 7/5/2019 for Pragya project in Horsham
        // This method reads  fasta query file and for each contig in the this file writes not mapped sequences as a new fasta file
        // In other words, this method uses the same query fasta file that was used in blast command, and estimates how much of that is covered by the blast
        // blast result file is of the format 6  qseqid qstart qend sseqid sstart send length pident mismatch gapopen gaps
        // for a blast result entry min_identity  is criteria to be selected
        // we assume blast result sorted by 'qseqid' but not with 'qstart' and 'qend'
        // The strategy is we first sort all the blast entries in a contig by qstart
        // we assume order of contigs showed up in the query fasta is the same as order of qseqid in blast result file., otherwise this method will raise an Exception
        
            FastaSequenceFile fsf = null;  // we sequentially read query fasta file so we dont need its index
            String last_qseqid = "";
            long last_contig_len = 0;  // length of the contig of last_qseqid
            String last_qseq = null;
            Map<Integer, Integer> posMap = new TreeMap(); // this Map holds start and end positions of each query hit in the last_contig, and sorted based on its key (qstart)
            //we report following 2 variables
            long total_contig_len = 0; // sum of length of all the contigs in query fasta file
            long total_not_covered = 0; // if a whole contig or some parts of query fasta file that is not mapped 
            FileWriter writer = null;

        try {
            long pbIndex = 0;
            double lastpercent = 0.0;
            long fileSz = new File(this.fileName).length();

            writer = new FileWriter(outFasta);
            
            fsf = new FastaSequenceFile(new File(queryFastaFile), true);
            String blast_line = reader.readLine();
            while (blast_line != null && !blast_line.isEmpty()){
                pbIndex += blast_line.getBytes().length; // we add number of bytes read
                double newpercent = ((double)pbIndex/fileSz)*100;
                if (newpercent-lastpercent>=1){
                    lastpercent = newpercent;
                    System.out.println("progress %"+(int)lastpercent);
                }

                
                String[] arr = blast_line.split("\\s");
                String qseqid = arr[0].trim();
                int qstart = Integer.parseInt(arr[1]);
                int qend = Integer.parseInt(arr[2]);
                double percent = Double.parseDouble(arr[7]);
                int alignment_len = Integer.parseInt(arr[6]);
                
                if(!qseqid.equals(last_qseqid)){
                    // process  last contig posMap
                    int cnt = 1; // this counter is used to be added at the end of query contig name, and reset to 1 for each query contig
                    if (!last_qseqid.isEmpty()){
                        // posMap is ordered based on the starting position
                        Iterator<Map.Entry<Integer, Integer>> it = posMap.entrySet().iterator();
                        int last_end = 0; // end position of the last entry
                        while(it.hasNext()){
                            Map.Entry<Integer,Integer> entry = it.next();
                            int start = entry.getKey();
                            int end = entry.getValue();
                            System.out.println(start + ":" + end);
                            if (start > last_end){
                                if (start - last_end > min_contig_length){
                                    writer.write(String.format(">%s_%d_%d_%d\n%s\n", last_qseqid, cnt++, last_end, start, last_qseq.substring(last_end, start)));
                                }
                                total_not_covered += start-last_end-1;
                            }
                            if(end>last_end)
                                last_end=end;
                        }// while in posMap
                        // we add remainder of the contig from last_end until contig_length
                        if (last_end < last_contig_len){
                            if (last_contig_len - last_end > min_contig_length){
                                writer.write(String.format(">%s_%d_%d_%d\n%s\n", last_qseqid, cnt++, last_end, (int)last_contig_len, last_qseq.substring(last_end, (int)last_contig_len)));
                            }
                            total_not_covered += last_contig_len - last_end;
                        }

                    }// if last_qseqid is not empty
                    
                    //look for qseqid in the query Fasta File contigs (we assume the order of contigs is the same in both files)
                    boolean found = false;
                    
                    ReferenceSequence seq = fsf.nextSequence();
                    while(seq != null){
                        total_contig_len += seq.length();
                        if (seq.getName().equals(qseqid)){
                            found = true;  
                            last_contig_len = seq.length();
                            last_qseq = seq.getBaseString();
                            break;
                        }else{ 
                            total_not_covered += seq.length(); // this whole contig is not mapped, so all of its length is not covered
                            //write whole query contig to output fasta file
                            writer.write(String.format(">%s\n%s\n", seq.getName(),  seq.getBaseString()));
                        }
                        
                        seq = fsf.nextSequence();
                    }
                    if(!found)
                        throw new Exception("Error: Order of blast result and original query fasta file are not the same.");
                    // empty the Map for new contig
                    posMap.clear();
                    last_qseqid = qseqid;
                }
                // add the query blast region to the Map if meets the criteria of min_percentage and min_alignment_length
                if (percent >= min_identity_percent)
                    posMap.put(qstart, qend);
                blast_line = reader.readLine();
            } //while
            // the last contig still needs to be processed
            if(posMap.size() > 0){
                int cnt = 1;
                Iterator<Map.Entry<Integer, Integer>> it = posMap.entrySet().iterator();
                int last_end = 0; // end position of the last entry
                while(it.hasNext()){
                    Map.Entry<Integer,Integer> entry = it.next();
                    int start = entry.getKey();
                    int end = entry.getValue();
                    System.out.println(start + ":" + end);
                    if (start > last_end){
                        if (start - last_end > min_contig_length){
                            writer.write(String.format(">%s_%d\n%s\n", last_qseqid, cnt++, last_qseq.substring(last_end, start)));
                        }                        
                        total_not_covered += start-last_end-1;
                    }
                    if(end>last_end)
                        last_end=end;
                }//while
                // we add remainder of the contig from last_end until contig_length
                if (last_end < last_contig_len){
                    if (last_contig_len - last_end > min_contig_length){
                        writer.write(String.format(">%s_%d\n%s\n", last_qseqid, cnt++, last_qseq.substring(last_end, (int)last_contig_len)));
                    }
                    
                    total_not_covered += last_contig_len - last_end;
                }

            }//if last contig has entries
            
            
            // now see if there is any remaining contigs in query fasta file
            ReferenceSequence seq = fsf.nextSequence();
            while(seq != null){
                total_contig_len += seq.length();
                total_not_covered += seq.length(); // this whole contig is not mapped, so all of its length is not covered
                //write whole query contig to output fasta file
                writer.write(String.format(">%s\n%s\n", seq.getName(),  seq.getBaseString()));
                seq = fsf.nextSequence();
            }            
            System.out.println("Total contig length: " + total_contig_len);
            System.out.println("Total not covered: " + total_not_covered);
            
            
        } catch (Exception ex) {System.out.println("getQueryRegionsCovered::" + ex.getMessage());}
        finally{
            fsf.close();
        }    
            
        
    }
    
    public void processExomeBlastResult(double min_identity_percent, int min_alignment_len){
        //written on 8/11/2017, This is an adhoc method used to find out why some exomes in Nimblegene exome capture reference are not found in the new generated exome reference
        Set<String>  d_genes = new HashSet();// it contains name of genes that are in A and B genome
        Set<String> a_b_genes = new HashSet(); // it contains name of genes that are in A and B genome
        try {
            String blast_line = reader.readLine();
            
            while (blast_line != null && !blast_line.isEmpty()){
                String[] arr = blast_line.split("\\s");
                String sseqid = arr[3].trim();
                double percent = Double.parseDouble(arr[7]); // i.e.  98.5
                int len = Integer.parseInt(arr[6]);

                if(percent>= min_identity_percent && len>= min_alignment_len){
                    //sseqid is something like:chr7D_exon455500_36493287_36497728_TraesCS7D01G065800 or chrUn_exon480299_90992481_90996366_TraesCSU01G103500
                    String[] parts = sseqid.split("_");
                    if(parts[0].indexOf('D')>0){
                        d_genes.add(parts[4]);
                    }else{
                        a_b_genes.add(parts[4]);
                    }
                    
                }
                
                blast_line = reader.readLine();
            }//while
            System.out.println("Number of D genomes " + d_genes.size());
            for(String gene : a_b_genes){
                System.out.println(gene);
            }

        } catch (Exception ex) {System.out.println(ex.getMessage());}
        finally{
           
        }  
        
    }
    
    public void filterResults(double min_identity_percent, double max_evalue, int alignment_len){
        // filters results of blast and creates a new blast result file(This is handy if blast result is too big to open it in excell or other editors for filtering)
        // This task reads a blast result of the format 6 qseqid qstart qend sseqid sstart send length pident mismatch gapopen gaps score evalue qseq sseq
         FileWriter blastWriter = null; // new fasta file with updated contig names is created by this
         FileWriter bedWriter = null; // new fasta file with updated contig names is created by this
         FileWriter gffWriter = null; // new fasta file with updated contig names is created by this
         
         try{ 
             blastWriter = new FileWriter(new File(new File(fileName).getParent(), "filtered_blast.txt"));
             bedWriter = new FileWriter(new File(new File(fileName).getParent(), "subject_bed.txt"));
             gffWriter = new FileWriter(new File(new File(fileName).getParent(), "subject_gff.txt"));
             String blast_line = reader.readLine();
             while (blast_line != null && !blast_line.isEmpty()){

                 String[] arr = blast_line.split("\\s");
                 double percent = Double.parseDouble(arr[7]);
                 int len = Integer.parseInt(arr[6]); // not used at this stage
                 double evalue = Double.parseDouble(arr[12]);
                 String scontig = arr[3];
                 int sstart = Integer.parseInt(arr[4]);
                 int send = Integer.parseInt(arr[5]);
                 String orientation = "+";
                 if(sstart>send){
                     //swap these
                     int tmp = sstart;
                     sstart = send;
                     send = tmp;
                     orientation = "-";
                 }
                 
                 if (evalue<=max_evalue && percent>=min_identity_percent && len>=alignment_len){
                     blastWriter.write(blast_line+"\n");
                     bedWriter.write(String.format("%s\t%d\t%d\n", scontig, sstart-1, send));
                     gffWriter.write(String.format("%s\t.\tHSP\t%d\t%d\t.\t%s\t.\tID=?\n", scontig, sstart, send, orientation));                     
                 }
                 
                 blast_line = reader.readLine();
             }         
         } catch (Exception ex) {System.out.println("Filter blast Result::" + ex.getMessage());}
        finally{
            //fsf.close();
            try {
                blastWriter.close();
                bedWriter.close();
                gffWriter.close();
            } catch (IOException ex) {}
        }  
    }
    public void fixHomoPolymere(String queryFastaFile, String newFastaFile, double min_identity_percent, int min_alignment_length){
        // This task reads a blast result of the format 6 qseqid qstart qend sseqid sstart send length pident mismatch gapopen gaps score evalue qseq sseq
        // If finds a homopolymere in the query that can be replaced with reference using alignment, then replaces it and updates query fasta file, eventually outputs updated query fasta file
        // This task was created to fix the problem of homopolymere in exome capture reference, arising from 454 sequencing technology
        //First we blast exome capture versus survey sequences, then we consider alignments for each query contig if they have at least 97% identity and 100 bp length (alignments should not overlap), then we report how many homopolymere were corrected
            FileWriter fastaWriter = null; // new fasta file with updated contig 
            FastaSequenceFile fsf = null;  // we sequentially read query fasta file so we dont need its index
            String last_qseqid = null;
            int replacement_count = 0 ;
            Map<Integer, Integer> alignmentMap = new HashMap(); // this Map holds start and end positions of each alignment in the current query contig, and sorted based on its key (qstart)

            ReferenceSequence qcontig = null; // always point to the current contig in queryFastaFile
            StringBuffer qbuff = new StringBuffer(); // a buffer to keep sequence of qcontig
            int step = 1;
        try {
            long pbIndex = 0;
            double lastpercent = 0.0;
            long fileSz = new File(this.fileName).length();

            fastaWriter = new FileWriter(newFastaFile);
            fsf = new FastaSequenceFile(new File(queryFastaFile), true);
            String blast_line = reader.readLine();
            step = 2;
            while (blast_line != null && !blast_line.isEmpty()){
                pbIndex += blast_line.getBytes().length; // we add number of bytes read
                double newpercent = ((double)pbIndex/fileSz)*100;
                if (newpercent-lastpercent>=1){
                    lastpercent = newpercent;
                    System.out.println("progress %"+(int)lastpercent);
                }
                
                String[] arr = blast_line.split("\\s");
                String qseqid = arr[0].trim();
                int qstart = Integer.parseInt(arr[1]);
                int qend = Integer.parseInt(arr[2]);
                String qseq = arr[13];
                String sseq = arr[14];               
                double percent = Double.parseDouble(arr[7]);
                int alignment_len = Integer.parseInt(arr[6]);
                if(!qseqid.equals(last_qseqid)){
                    //look for qseqid in the query Fasta File contigs (we assume the order of contigs is the same in both files)
                    boolean found = false;
                    qcontig = fsf.nextSequence();
                    while(qcontig != null){
                        if (qcontig.getName().equals(qseqid)){
                            found = true;  
                            break;
                        }else{
                           fastaWriter.write(String.format(">%s\n", qcontig.getName()));                           
                           fastaWriter.write( new String(qcontig.getBases(), "UTF-8") + "\n"); 
                        } 
                        qcontig= fsf.nextSequence();
                    }
                    if(!found)
                        throw new Exception("Error: Order of blast result and original query fasta file are not the same.");
                    // so qcontig always points to the current contig in the query fasta
                    // process  last query contig 
                    if (last_qseqid != null){
                        // write last query contig in qbuff with updated homopolymere into the output fasta file
                        // First replace any Gap (-) in the buffer due to alignment to reference
                        String new_qcontig = qbuff.toString().replaceAll("-", "");
                        // now write it
                        fastaWriter.write(String.format(">%s\n", last_qseqid));                           
                        fastaWriter.write(new_qcontig + "\n"); 
                    }// if last_qseqid is not null
                    // empty the Map for new contig
                    alignmentMap.clear();
                    // reset query buffer
                    qbuff.delete(0, qbuff.length());
                    // initialize query contig buffer
                    qbuff.append(new String(qcontig.getBases(), "UTF-8") ); 
                    last_qseqid = qseqid;
                }
                // add the query blast region to the Map if meets the criteria of min_percentage and min_alignment_length
                if (percent >= min_identity_percent  && alignment_len>=min_alignment_length){
                    // current alignment has criteria to be processed, check if there is no overlap with previous aligment regions in the current query contig
                    if (!new MapUtil().overlap(alignmentMap, qstart, qend)){
                        // check for N and homopolymere in query seq (to be replaced by reference if reference is any better)
                        Pattern p = Pattern.compile("A{4,}+|C{4,}+|T{4,}+|G{4,}+"); // possessive(+)
                        Matcher m = p.matcher(qseq); // rember that qseq contains Gaps or -, but we only look for continuous base without gap in the middle
                        while(m.find()){
                            // extract the equal part from reference (subject)
                            int sidx = m.start();
                            int eidx = m.end();
                            String ref_part = sseq.substring(sidx, eidx);
                            String query_part = qseq.substring(sidx, eidx);
                            if (Helper.isGood(ref_part) && !ref_part.equals(query_part)){
                                // if reference part is different from 
                                //System.out.println(qseqid + ":" + qcontig.getName() + ":" + qbuff.length() + ":" + qstart + ":" + sidx + ":" + eidx + ":" + ref_part.length() );
                                int gaps = Helper.getGaps(qseq.substring(0, sidx));
                                qbuff.replace(qstart-1 + sidx -gaps, qstart-1 + eidx -gaps, ref_part.toLowerCase());  // remember that ref_part might have Gaps that at the end must be removed
                                replacement_count++;
                            }
                            //System.out.println(txtInput1.getText().substring(m.start(), m.end()));
                        }
                        alignmentMap.put(qstart, qend);
                    }
                }
                blast_line = reader.readLine();
            } //while
            // the last query contig still needs to be processed, because it did not have a chance in while loop
            String new_qcontig_seq = qbuff.toString().replace('-', '\0');
            // now write it
            fastaWriter.write(String.format(">%s\n", qcontig.getName()));                           
            fastaWriter.write(new_qcontig_seq + "\n"); 
           
            //Write rest of queryFastaFile to the new query fasta file
            qcontig = fsf.nextSequence();
            while(qcontig != null){
               fastaWriter.write(String.format(">%s\n", qcontig.getName()));                           
               fastaWriter.write( new String(qcontig.getBases(), "UTF-8") + "\n"); 
               qcontig= fsf.nextSequence();
            }
           System.out.println("Total HomoPolymere Replacements: " + replacement_count);
           // System.out.println("Total not covered: " + total_not_covered);
            
            
        } catch (Exception ex) {System.out.println("fixHomoPolymere::" + ex.getMessage() + " Happened at step:" + step);}
        finally{
            
            try {
                fastaWriter.close();
            } catch (IOException ex) {}
            fsf.close();
        }    
        
    }
    public void getKmereReadsInfo(double min_identity_percent, int min_alignment_length, double max_e_value){
        // This method process a blast result of the format: 6 qseqid qstart qend sseqid sstart send length pident mismatch gapopen gaps score evalue qseq sseq
        // For each query, finds first hit that is normally the best hit (blast sorts output result based on the query), then see if this hit meets min_identity, max_e_value and min_length
        // Then it breaks down query sequence and subject sequence to k-mer sequences (k=100) and in each 100 bp sequence reads number of mismatch, gap openning and total number of the gaps and gaps that are multiple of 3 (3,6,9,....)
        // Then writes output in a file in the same directory of the blast result file.
        
        
        
        
        
    }
    
    public void merge_quey_subject_combined(String queryFile, String queryIndexFile, String subjectFile, String subjectIndexFile, String resultFile , int min_alignment_length){
    // This method merges query and subject fasta file, the resultant file should contain both but not redunduncy
    //This method assumes that blast might report one big perfect match in the form of multiple perfect smaller mathes, thats why uses a TreeMap to store and sort small regions and finally combine them all. 
    // Its application could be where we want to merge 2 assemblies where 2 different tools has been used to assmble the same set of reads
    // For each query contig: 1-Query is not mached to any subject in this case query is outputed as it is 2- Query is contained in a subject, in this case subject is outputed as it is 3- query can extend one subject uniquely 4- query can bridge between 2 subjects uniquely
    // The output is merging result fasta file 
        FastaSequenceFile fsf = null;  // we sequentially read query fasta file with this
        String blast_line = null;
        String last_query_contig = null; // blast result is sorted by query contig and then subject contig
        int last_query_index = 0; //zero-based index of query contig in the fasta file
        int last_query_contig_len = 0;
       
        String last_subject_contig = null; // blast result is sorted by subject contig
        int last_subject_index = 0; //zero-based index of subject contig in the fasta file
        int last_subject_contig_len = 0;
        
        
        FileWriter writer = null;
        TreeMap<Integer, Integer> query_posMap_plus = new TreeMap(); // this Map holds start and end positions of each hit in the query, and sorted based on its key (qstart)
        TreeMap<Integer, Integer> query_posMap_minus = new TreeMap(); // this Map holds start and end positions of each hit in the query, and sorted based on its key (qstart)
        
        TreeMap<Integer, Integer> subject_posMap_plus = new TreeMap(); // this Map holds start and end positions of each positive hit in the last_contig, and sorted based on its key (sstart)
        TreeMap<Integer, Integer> subject_posMap_minus = new TreeMap(); // this Map holds start and end positions of each negative hit in the last_contig, and sorted based on its key (sstart)
        
        List<ContigExtension>  extensions = new ArrayList();
        
        int step = 0;
        FastaSequenceIndex sfsi = new FastaSequenceIndex(new File(subjectIndexFile));
        IndexedFastaSequenceFile sifsi =new IndexedFastaSequenceFile(new File(subjectFile) , sfsi);

        //this is for indexed access to the query file
        FastaSequenceIndex qfsi = new FastaSequenceIndex(new File(queryIndexFile));
        IndexedFastaSequenceFile qifsi =new IndexedFastaSequenceFile(new File(queryFile) , qfsi);
        
        Map<Integer, String> q_id_map = new HashMap(); // a map from positional index of contig (0-based) to the actual contig name
        Map<Integer, String> s_id_map = new HashMap(); // a map from positional index of contig (0-based) to the actual contig name
        
        MapUtil util = new MapUtil();
        try{ 
                     
            writer = new FileWriter(resultFile);
            //read through the blast file
            
            blast_line = reader.readLine();
            step = 1;
            while (blast_line != null && !blast_line.isEmpty()){
                String[] arr = blast_line.split("\\s");
                String sseqid = arr[3].trim();
                double percent = Double.parseDouble(arr[7]);
                int alignment_len = Integer.parseInt(arr[6]); 
                double evalue = Double.parseDouble(arr[12]);
                String qseqid = arr[0].trim();
                int qstart = Integer.parseInt(arr[1]);
                int qend = Integer.parseInt(arr[2]);
                int sstart = Integer.parseInt(arr[4]);
                int send = Integer.parseInt(arr[5]);
              //  if(qseqid.equals("124"))
              //              step = 1;
                
                if(!qseqid.equals(last_query_contig) || !sseqid.equals(last_subject_contig)){
                    // process  last contig posMap
                    if (last_subject_contig != null){//this if makes sure we are not processing the first line of the blast file
                        // process data related to last_query_contig and last_subject_contig
                        // now we investigate 3 posMaps
                        
                        /*
                        int qpluslen = 0;
                        int qminuslen = 0;
                        int biggest_plus_len = 0;
                        int biggest_minus_len = 0;
                        int subject_plus_len = 0;
                        int subject_minus_len = 0;
                        int plus_difference = -1;  // it is very important to be -1, not to be confused by 0
                        int minus_difference = -1;  // it is very important to be -1, not to be confused by 0
                        int query_plus_covered = 0;  // qpluslen is important, not query_plus_covered
                        int query_minus_covered = 0;  // qminuslen is important, not query_minus_covered
                        int subject_plus_covered = 0;
                        int subject_minus_covered = 0;
                       */
                        step = 3;
                        //first process plus alignment
                        if (query_posMap_plus.size() > 0){
                            int query_start = query_posMap_plus.firstKey();
                            int query_end = util.right_end(query_posMap_plus);
                            
                            int subject_start_plus = subject_posMap_plus.firstKey();
                            int subject_end_plus =  util.right_end(subject_posMap_plus);                            
                            if(query_start == 1 && query_end == last_query_contig_len && subject_start_plus>1){
                                //whole query contig is contained inside subject
                                //take action accordingly
                                ContigExtension ext = new ContigExtension();
                                ext.contigA = last_query_index;
                                ext.contigB = last_subject_index;
                                ext.relative_position = Helper.POS_Q_Inside_S;
                                //orientation and index is not important
                                extensions.add(ext);
                                
                            }else if ((subject_start_plus == 1 && subject_end_plus == last_subject_contig_len && query_start>1) ){
                                //whole subject contig is inside query contig
                                ContigExtension ext = new ContigExtension();
                                ext.contigA = last_query_index;
                                ext.contigB = last_subject_index;
                                ext.relative_position = Helper.POS_S_Inside_Q;
                                //orientation and index is not important
                                extensions.add(ext);
                            }else if (query_start == 1 && subject_end_plus == last_subject_contig_len ){
                                //subject is extended by query
                                //take action accordingly
                                ContigExtension ext = new ContigExtension();
                                ext.contigA = last_query_index;
                                ext.contigB = last_subject_index;
                                ext.relative_position = Helper.POS_S_Extended_Q;
                                ext.orientation = '+';
                                ext.index = query_end ; //length to be merged
                                extensions.add(ext);
                                
                            }else if ( subject_start_plus == 1 && query_end == last_query_contig_len ){
                                //query is extended by subject 
                                ContigExtension ext = new ContigExtension();
                                ext.contigA = last_query_index;
                                ext.contigB = last_subject_index;
                                ext.relative_position = Helper.POS_Q_Extended_S;
                                ext.orientation = '+';
                                ext.index = subject_end_plus; //length to be merged
                                extensions.add(ext);
                                
                                //take action accordingly
                            }
                        }
                        
                        //now process minus alignment 
                        if (query_posMap_minus.size() > 0){
                            int query_start = query_posMap_minus.firstKey();
                            int query_end = util.right_end(query_posMap_minus);

                            int subject_start_minus = subject_posMap_minus.firstKey();
                            int subject_end_minus =  util.right_end(subject_posMap_minus);                            
                            
                            if(query_start == 1 && query_end == last_query_contig_len && subject_start_minus>1){
                                //whole query contig is contained inside subject
                                //take action accordingly
                                ContigExtension ext = new ContigExtension();
                                ext.contigA = last_query_index;
                                ext.contigB = last_subject_index;
                                ext.relative_position = Helper.POS_Q_Inside_S;
                                //orientation and index is not important
                                extensions.add(ext);
                                
                            }else if (subject_start_minus == 1 && subject_end_minus == last_subject_contig_len && query_start>1){
                                //whole subject contig is inside query contig
                                ContigExtension ext = new ContigExtension();
                                ext.contigA = last_query_index;
                                ext.contigB = last_subject_index;
                                ext.relative_position = Helper.POS_S_Inside_Q;
                                //orientation and index is not important
                                extensions.add(ext);
                            }else if (query_start == 1 && subject_start_minus == 1){ 
                                //RC of subject is extended by query
                                ContigExtension ext = new ContigExtension();
                                ext.contigA = last_query_index;
                                ext.contigB = last_subject_index;
                                ext.relative_position = Helper.POS_S_Extended_Q;
                                ext.orientation = '-'; //- here implies RC of subject
                                ext.index = query_end; //length to be merged
                                extensions.add(ext);
                                
                            }else if(query_end == last_query_contig_len && subject_end_minus == last_subject_contig_len){
                                //query is extended by RC of subject 
                                ContigExtension ext = new ContigExtension();
                                ext.contigA = last_query_index;
                                ext.contigB = last_subject_index;
                                ext.relative_position = Helper.POS_Q_Extended_S;
                                ext.orientation = '-';//- here implies RC of subject
                                ext.index = query_end - query_start + 1; //length to be merged
                                extensions.add(ext);                                
                            }

                        }
                    } // if we are not at the first line of blast result
                    step = 5;
                    if(!qseqid.equals(last_query_contig)){
                        last_query_contig_len = qifsi.getSequence(qseqid).length();
                        last_query_index = qifsi.getSequence(qseqid).getContigIndex(); //0-based index
                        q_id_map.put(last_query_index, qseqid);
                    }
                    step = 6;
                    if(!sseqid.equals(last_subject_contig)){
                        last_subject_contig_len= sifsi.getSequence(sseqid).length();
                        last_subject_index = sifsi.getSequence(sseqid).getContigIndex(); //0-based index
                        s_id_map.put(last_subject_index, sseqid);
                    }
                    last_subject_contig = sseqid;
                    last_query_contig = qseqid;
                    step = 7;
                    // empty the Map for new contig
                    query_posMap_plus.clear();
                    query_posMap_minus.clear();
                    subject_posMap_plus.clear();
                    subject_posMap_minus.clear();
                    step = 8;
                }// if query contig or subject contig has changed
                 
                if (alignment_len>=min_alignment_length && percent==100){  
                    
                   if(sstart<send){
                       // plus-plus 
                       query_posMap_plus.put(qstart, qend);
                       subject_posMap_plus.put(sstart, send);
                       //query_posMap_plus.put(qstart, qend);
                   }else{
                       //plus-minus
                       query_posMap_minus.put(qstart, qend);
                       subject_posMap_minus.put(send, sstart);                       
                   }
                }
                 step = 9;
                blast_line = reader.readLine();
                
            }//while loop in blast file lines
             // the last subject contig still needs to be processed
            //first process plus alignment
            if (query_posMap_plus.size() > 0){
                int query_start = query_posMap_plus.firstKey();
                int query_end = util.right_end(query_posMap_plus);

                int subject_start_plus = subject_posMap_plus.firstKey();
                int subject_end_plus =  util.right_end(subject_posMap_plus);                            
                if(query_start == 1 && query_end == last_query_contig_len && subject_start_plus>1){
                    //whole query contig is contained inside subject
                    //take action accordingly
                    ContigExtension ext = new ContigExtension();
                    ext.contigA = last_query_index;
                    ext.contigB = last_subject_index;
                    ext.relative_position = Helper.POS_Q_Inside_S;
                    //orientation and index is not important
                    extensions.add(ext);

                }else if ((subject_start_plus == 1 && subject_end_plus == last_subject_contig_len && query_start>1) ){
                    //whole subject contig is inside query contig
                    ContigExtension ext = new ContigExtension();
                    ext.contigA = last_query_index;
                    ext.contigB = last_subject_index;
                    ext.relative_position = Helper.POS_S_Inside_Q;
                    //orientation and index is not important
                    extensions.add(ext);
                }else if (query_start == 1 && subject_end_plus == last_subject_contig_len ){
                    //subject is extended by query
                    //take action accordingly
                    ContigExtension ext = new ContigExtension();
                    ext.contigA = last_query_index;
                    ext.contigB = last_subject_index;
                    ext.relative_position = Helper.POS_S_Extended_Q;
                    ext.orientation = '+';
                    ext.index = query_end ; //length to be merged
                    extensions.add(ext);

                }else if ( subject_start_plus == 1 && query_end == last_query_contig_len ){
                    //query is extended by subject 
                    ContigExtension ext = new ContigExtension();
                    ext.contigA = last_query_index;
                    ext.contigB = last_subject_index;
                    ext.relative_position = Helper.POS_Q_Extended_S;
                    ext.orientation = '+';
                    ext.index = subject_end_plus; //length to be merged
                    extensions.add(ext);

                    //take action accordingly
                }
            }

            //now process minus alignment 
            if (query_posMap_minus.size() > 0){
                int query_start = query_posMap_minus.firstKey();
                int query_end = util.right_end(query_posMap_minus);

                int subject_start_minus = subject_posMap_minus.firstKey();
                int subject_end_minus =  util.right_end(subject_posMap_minus);                            

                if(query_start == 1 && query_end == last_query_contig_len && subject_start_minus>1){
                    //whole query contig is contained inside subject
                    //take action accordingly
                    ContigExtension ext = new ContigExtension();
                    ext.contigA = last_query_index;
                    ext.contigB = last_subject_index;
                    ext.relative_position = Helper.POS_Q_Inside_S;
                    //orientation and index is not important
                    extensions.add(ext);

                }else if (subject_start_minus == 1 && subject_end_minus == last_subject_contig_len && query_start>1){
                    //whole subject contig is inside query contig
                    ContigExtension ext = new ContigExtension();
                    ext.contigA = last_query_index;
                    ext.contigB = last_subject_index;
                    ext.relative_position = Helper.POS_S_Inside_Q;
                    //orientation and index is not important
                    extensions.add(ext);
                }else if (query_start == 1 && subject_start_minus == 1 ){ 
                    //RC of subject is extended by query
                    ContigExtension ext = new ContigExtension();
                    ext.contigA = last_query_index;
                    ext.contigB = last_subject_index;
                    ext.relative_position = Helper.POS_S_Extended_Q;
                    ext.orientation = '-'; //- here implies RC of subject
                    ext.index = query_end; //length to be merged
                    extensions.add(ext);

                }else if(query_end == last_query_contig_len && subject_end_minus == last_subject_contig_len){
                    //query is extended by RC of subject 
                    ContigExtension ext = new ContigExtension();
                    ext.contigA = last_query_index;
                    ext.contigB = last_subject_index;
                    ext.relative_position = Helper.POS_Q_Extended_S;
                    ext.orientation = '-';//- here implies RC of subject
                    ext.index = query_end - query_start + 1; //length to be merged
                    extensions.add(ext);                                
                }

            }
            System.out.println("Number of extensions detected: " + extensions.size());
            step = 10;
            /*
            For example if there is following extension between queries and subjects, then 2 matrixes will be filled like this:
                  ---------qi----                                 -------qm------
            ---sk-----        -------R.C. Sj---------------------------
            
            query_extensions                         subject_Extensions
index       0   |1  |2  |3              index        0  |1  |2  |3  
            --------------                           ------------------
            .............                            ..................
i'th row    j   |20 |0  |k                           ..................
            .............               j'th row     m  |15 |0  |i
m'th row    -1  |-1 |-1 |j                           ..................
            .............               k'th row     i  |30 |1  |-1
                                                     ..................
            
            In the third column 0 means - and 1 means +
            */
            //now we process extensions, There are some rules we apply:
            // (1) every single query or subject contig,  can only participate in just one chain/merging action
            // (2) if one subject or query is extended by multiple query/subject contigs we only consider longest one
            //In order to implement these rules, we prepare 2 arrays one for processed action one for extension 
            int subject_size = sfsi.size();
            int query_size = qfsi.size();
            
            System.out.println("Number of contigs in query fasta file: " + query_size);
            System.out.println("Number of contigs in subject fasta file: " + subject_size);
            int[][] query_extensions = new int[query_size][4];
            int[][] subject_extensions = new int[subject_size][4];
            boolean[] query_processed = new boolean[query_size]; //if a query contig is part of a valid chain then it is set to true
            boolean[] subject_processed = new boolean[subject_size];//if a subject contig is part of a valid chain then it is set to true
            //definition of 4 coulmns [][0]=extended by contig index  [][1]=merging length [][2]= merging orientation 1: forward 0: reverse  [][3]=extends contig index (opposite of [][0])  
            //now initialize arrays:
            for (int i=0;i<subject_size; i++){
                subject_processed[i] = false;
                subject_extensions[i][0]=-1;
                subject_extensions[i][1]=-1;
                subject_extensions[i][2]=-1;
                subject_extensions[i][3]=-1;
            }

            for (int i=0;i<query_size; i++){
                query_processed[i] = false;
                query_extensions[i][0]=-1;
                query_extensions[i][1]=-1;
                query_extensions[i][2]=-1;
                query_extensions[i][3]=-1;
                
            }
            step = 11;
            //now it is time to process extensions list and fill up the extensions arrays
            //For each query and subject there is only one row in the matrix so each contig can only have one action
            //rule: if  contig A is inside contig B, then if A is extended or extends contig C, then we ignore extension action, in other words inside action will override any other extended by action.
            //but if a contig A is inside contig B and contig B is extended by contig C, then both actions (inside and extended by) are stored
            for(ContigExtension ext: extensions){
                if (ext.relative_position == Helper.POS_Q_Inside_S){
                    query_extensions[ext.contigA][0] = ext.contigB;
                    query_extensions[ext.contigA][1] = -1; //index 1 along with index 0 is used to distiguish inside action vs extension action, inside action is valid if index 0 is not -1 and index 1 is -1 
                    query_extensions[ext.contigA][2] = -1;
                    query_extensions[ext.contigA][3] = -1;
                    //coulmns 1,2 and 3 will remain -1, so we will realize that the action is inside
                    //we also stores opposite of this action : contain action
                    //subject_extensions[ext.contigB][3] = ext.contigA;
                }else if (ext.relative_position == Helper.POS_S_Inside_Q){
                    subject_extensions[ext.contigB][0] = ext.contigA;
                    subject_extensions[ext.contigB][1] = -1; //index 1 along with index 0 is used to distiguish inside action vs extension , inside action is valid if index 0 is not -1 and index 1 is -1 
                    subject_extensions[ext.contigB][2] = -1;
                    subject_extensions[ext.contigB][3] = -1;
                    //coulmns 1,2 and 3 will remain -1, so we will realize that the action is inside
                     //we also store opposite of this action : contain action
                    //query_extensions[ext.contigA][3] = ext.contigB;
                }else if (ext.relative_position == Helper.POS_Q_Extended_S){
                    //first make sure this query is not already inside a subject
                    //see if query is extended by something beforehand, if new extension is longer then we replace it
                    if (query_extensions[ext.contigA][0]==-1 ){
                        query_extensions[ext.contigA][0]=ext.contigB;
                        subject_extensions[ext.contigB][3]=ext.contigA; // it is used for back tracking
                        //now store the length of merging
                        query_extensions[ext.contigA][1] = ext.index;
                        //now store orientation of merging
                        query_extensions[ext.contigA][2] = ext.orientation=='+'?1:0;
                    }else if (query_extensions[ext.contigA][1]!=-1 && query_extensions[ext.contigA][1]<ext.index){
                        // we get here when this query has already been extended by another subject contig with smaller ovelap length, so according to the rule we update it
                        // before updating query with new extension, we need to make sure we reset all backtrack flags
                        subject_extensions[query_extensions[ext.contigA][0]][3] = -1;

                        //now update query contig 
                        query_extensions[ext.contigA][0]=ext.contigB;
                        subject_extensions[ext.contigB][3]=ext.contigA; // it is used for back tracking
                        //now store the length of merging
                        query_extensions[ext.contigA][1] = ext.index;
                        //now store orientation of merging
                        query_extensions[ext.contigA][2] = ext.orientation=='+'?1:0;
                        
                    }
                   
                }else if (ext.relative_position == Helper.POS_S_Extended_Q){
                    //first make sure this subject is not already inside a query
                     //see if subject is extended by some query beforehand, if new extension is longer then we replace it
                    if (subject_extensions[ext.contigB][0]==-1 ){
                        subject_extensions[ext.contigB][0]=ext.contigA;
                        query_extensions[ext.contigA][3]=ext.contigB;
                         //now store the length of merging
                        subject_extensions[ext.contigB][1] =  ext.index;
                        //now store orientation of merging
                        subject_extensions[ext.contigB][2] = ext.orientation=='+'?1:0;
                    }else if (subject_extensions[ext.contigB][1] != -1 && subject_extensions[ext.contigB][1]<ext.index){
                        // we get here when this subject has already been extended by another contig with smaller ovelap length, so according to the rule we update it
                        // before updating subject with new extension, we need to make sure we reset all backtrack flags
                        query_extensions[subject_extensions[ext.contigB][0]][3] = -1;
                        //now updating subject_Extension row
                        subject_extensions[ext.contigB][0]=ext.contigA;
                        query_extensions[ext.contigA][3]=ext.contigB;
                         //now store the length of merging
                        subject_extensions[ext.contigB][1] =  ext.index;
                        //now store orientation of merging
                        subject_extensions[ext.contigB][2] = ext.orientation=='+'?1:0;
                        
                    }
                    
                }
            }
            // matrix has been loaded with correct data
            System.out.println("Matrix Extension has been Loaded.");
            //now we go through every contig in query, if it is part of chain we will detect the chain and merge them and report them as one single contig
           step = 12;
           ContigExtender2 ex = new ContigExtender2(query_extensions, subject_extensions, query_processed, subject_processed);
           
           long total_q_bp = 0;
           long merged_bp = 0; //total base pairs have been merged
           long written_bp = 0;
           int qcontig_idx = 0;
           fsf = new FastaSequenceFile(new File(queryFile), true); 
           ReferenceSequence refseq = fsf.nextSequence();
           
           while (refseq != null){
                //String contig_seq = new String(refseq.getBases());
                //contigs.put(idx, new CompactSequence2(contig_seq, idx));
                //contigs.add(contig_seq);
               total_q_bp += refseq.length();
               if (query_extensions[qcontig_idx][0]==-1 && query_extensions[qcontig_idx][3]==-1){
                   //This means this query is not inside any subject and also is not part of any chain, so write it as it is (it is lonely)
                   
                    writer.write(">"+refseq.getName()+"\n");
                    writer.write(refseq.getBaseString() +"\n");  
                    written_bp += refseq.length();
               }else if (query_extensions[qcontig_idx][0]!=-1 && query_extensions[qcontig_idx][1] == -1 ){
                   //this means this is inside action, query is inside a subject, we dont write anything now, we will write subject when we process subject fasta file
                   //this query may extends a subject or be extended by a subject, but we are not interested in them, because nside action overrides other actions for this query
                   merged_bp += refseq.length();
               }else if (query_extensions[qcontig_idx][0]!=-1 || query_extensions[qcontig_idx][3]!=-1) {
                   //Here means query is part of possible extension chain (we may be able go foreward or backward)
                   //now our task is to discover the chain, before this, we need to know if this query contig has been part of a chain already.
                   if (!query_processed[qcontig_idx]){
                       // This means query has not already been part of any chain
                       List<ContigExtension> list = ex.getNextChain(qcontig_idx);
                       StringBuffer chainStr = new StringBuffer();
                       try{
                           
                       if(list.size()==0){
                           // This happens when this contig can be extended by another contig that already been processed, so list will be ampty
                           //output contig as it is
                           step = 120;
                            writer.write(">"+refseq.getName()+"\n");
                            writer.write(refseq.getBaseString()+"\n");
                            merged_bp += refseq.length();
                       }else{
                           // merge elements of the extension chain
                           //because we also go backward in the chain, then a chain may start with a subject
                           for(ContigExtension cext : list)
                               chainStr.append(cext.toString());
                           step = 121;
                           StringBuffer seq_buff = new StringBuffer();
                           StringBuffer name_buff = new StringBuffer();
                           String first_contig = null;
                           if (list.get(0).relative_position == Helper.POS_Q_Extended_S){
                               //first contig is query
                               first_contig = qifsi.getSequence(q_id_map.get(list.get(0).contigA)).getBaseString();
                               name_buff.append("q_" + list.get(0).contigA);
                           }else if (list.get(0).relative_position == Helper.POS_S_Extended_Q){
                               //first contig is subject
                               if (list.get(0).orientation == '+'){
                                   first_contig = sifsi.getSequence(s_id_map.get(list.get(0).contigB)).getBaseString();
                                   name_buff.append("s+_" + list.get(0).contigB);
                               }else if (list.get(0).orientation == '-'){
                                   first_contig = Helper.getRC(sifsi.getSequence(s_id_map.get(list.get(0).contigB)).getBaseString());
                                   name_buff.append("s-_" + list.get(0).contigB);                                   
                               }
                               
                               
                           }
                           seq_buff.append(first_contig);
                        step = 122;
                           for(ContigExtension ext : list){
                               if(ext.relative_position == Helper.POS_Q_Extended_S){
                                   
                                   String ssequence = sifsi.getSequence(s_id_map.get(ext.contigB)).getBaseString();
                                   if (ext.orientation == '+'){
                                       name_buff.append("_s+_" + ext.contigB);
                                       seq_buff.append(ssequence.substring(ext.index));
                                   }else if (ext.orientation == '-'){
                                       name_buff.append("_s-_" + ext.contigB);
                                       seq_buff.append(Helper.getRC(ssequence).substring(ext.index));
                                   }
                               }else if (ext.relative_position == Helper.POS_S_Extended_Q){
                                   name_buff.append("_q_" + ext.contigA);
                                   String qsequence = qifsi.getSequence(q_id_map.get(ext.contigA)).getBaseString();
                                   seq_buff.append(qsequence.substring(ext.index));                                    
                               }
                               
                               merged_bp += ext.index;
                           }
                           step = 123;
                           writer.write(">"+name_buff.toString()+"\n"); 
                           writer.write(seq_buff.toString() +"\n");
                           written_bp += seq_buff.length();

                       }//else chain has any elements
                       }catch(Exception exp){System.out.println(exp.toString() + " STEP:" + step + chainStr.toString());}
                   }//if query contig has not processed before
               }// if query contig can be part of a chain
                refseq = fsf.nextSequence();
                qcontig_idx++;
           }// while  query contig         

           System.out.println("Query contigs has been processed.");
           //now it is time to process subject contigs, for subject contigs we dont need to find chains, becuase all the subjects participated in any chain was already processed and merged and written before
           long total_s_bp = 0;
           int special = 0;
           int scontig_idx = 0;
           fsf.close();
           fsf = new FastaSequenceFile(new File(subjectFile), true); 
           refseq = fsf.nextSequence();           
           while (refseq != null){
                total_s_bp += refseq.length();
                if (subject_extensions[scontig_idx][0]==-1 && subject_extensions[scontig_idx][3]==-1){
                   //This means this subject contig is not inside any query and also is not part of any chain, so write it as it is                   
                    writer.write(">"+refseq.getName()+"\n");
                    writer.write(refseq.getBaseString() +"\n");  
                    written_bp += refseq.length();
               }else if (subject_extensions[scontig_idx][0]!=-1 && subject_extensions[scontig_idx][3] == -1 ){
                   //this means this is inside action, subject is inside a quey, we dont write anything, because it has been written when processed query contigs
                   merged_bp += refseq.length();
               }else if (!subject_processed[scontig_idx]){
                   //for any other reason, that subject has not been processed then we write it
                    writer.write(">"+refseq.getName()+"\n");
                    writer.write(refseq.getBaseString() +"\n");  
                    written_bp += refseq.length();    
                    special++;
               }
                refseq = fsf.nextSequence();
                scontig_idx++;
           }
            
           System.out.println("Total Query contigs bp:" + total_q_bp);
           System.out.println("Total Subject contigs bp:" + total_s_bp);
           System.out.println("Total new contigs written bp:" + written_bp);
           System.out.println("Total merged bp:" + merged_bp);
           System.out.println("Total Subject Contigs with Special cases:" + special);
           
        } catch (Exception ex) {System.out.println("Merging query and subject fasta files " + ex.getMessage() + " at step "+step + blast_line);}
        finally{
            try {
                
                writer.close();
                qifsi.close();
                sifsi.close();
                fsf.close();
            } catch (Exception ex) {}
        } 
        
    }

    
    public void merge_quey_subject_simple(String queryFile, String queryIndexFile, String subjectFile, String subjectIndexFile, String resultFile , int min_alignment_length){
     //Attention:This exactly like merge_quey_subject_combined method with this difference that we dont combine HSP for a specific subject and query , so logic gets much simpler, we assume Blast resul is trustable and HSP s that have overlapped already merged to a bigger HSP by blast itself
    //1) This method merges query and subject fasta file, the resultant file should contain both but not redunduncy
    // 2)Its application could be where we want to merge 2 assemblies where 2 different tools has been used to assmble the same set of reads
    // 3)For each query contig: 1-Query is not matched to any subject: in this case query is outputed as it is 2- Query is contained in a subject: in this case subject is outputed as it is 3- query can extend one subject uniquely 4- query can bridge between 2 subjects uniquely
    // 4)blast result is sorted by query contig and then subject contig and then total alignment score
    // The output is merging result fasta file 
        FastaSequenceFile fsf = null;  // we sequentially read query fasta file with this
        String blast_line = null;
        
        FileWriter writer = null;

        List<ContigExtension>  extensions = new ArrayList();
        
        int step = 0;
        //this is for indexed access to the subject fasta file
        FastaSequenceIndex sfsi = new FastaSequenceIndex(new File(subjectIndexFile));
        IndexedFastaSequenceFile sifsi =new IndexedFastaSequenceFile(new File(subjectFile) , sfsi);

        //this is for indexed access to the query fasta file
        FastaSequenceIndex qfsi = new FastaSequenceIndex(new File(queryIndexFile));
        IndexedFastaSequenceFile qifsi =new IndexedFastaSequenceFile(new File(queryFile) , qfsi);
        
        Map<Integer, String> q_id_map = new HashMap(); // a map from positional index of contig (0-based) to the actual contig name
        Map<Integer, String> s_id_map = new HashMap(); // a map from positional index of contig (0-based) to the actual contig name
        
        MapUtil util = new MapUtil();
        try{ 
                     
            writer = new FileWriter(resultFile);
            //read through the blast file
            
            blast_line = reader.readLine();
            step = 1;
            while (blast_line != null && !blast_line.isEmpty()){
                String[] arr = blast_line.split("\\s");
                String sseqid = arr[3].trim();
                double percent = Double.parseDouble(arr[7]);
                int alignment_len = Integer.parseInt(arr[6]); 
                double evalue = Double.parseDouble(arr[12]);
                String qseqid = arr[0].trim();
                int qstart = Integer.parseInt(arr[1]);
                int qend = Integer.parseInt(arr[2]);
                int sstart = Integer.parseInt(arr[4]);
                int send = Integer.parseInt(arr[5]);
                
                if (alignment_len>=min_alignment_length && percent==100){
                    int last_query_contig_len = qifsi.getSequence(qseqid).length();
                    int last_query_index = qifsi.getSequence(qseqid).getContigIndex(); //0-based index
                    q_id_map.put(last_query_index, qseqid);
                    int last_subject_contig_len= sifsi.getSequence(sseqid).length();
                    int last_subject_index = sifsi.getSequence(sseqid).getContigIndex(); //0-based index
                    s_id_map.put(last_subject_index, sseqid);

                    if (sstart<send){
                      //++
                        int subject_start_plus = sstart;
                        int subject_end_plus  = send;
                        if(qstart == 1 && qend == last_query_contig_len && subject_start_plus>1){
                            //whole query contig is contained inside subject
                            //take action accordingly
                            ContigExtension ext = new ContigExtension();
                            ext.contigA = last_query_index;
                            ext.contigB = last_subject_index;
                            ext.relative_position = Helper.POS_Q_Inside_S;
                            //orientation and index is not important
                            extensions.add(ext);

                        }else if ((subject_start_plus == 1 && subject_end_plus == last_subject_contig_len && qstart>1) ){
                            //whole subject contig is inside query contig
                            ContigExtension ext = new ContigExtension();
                            ext.contigA = last_query_index;
                            ext.contigB = last_subject_index;
                            ext.relative_position = Helper.POS_S_Inside_Q;
                            //orientation and index is not important
                            extensions.add(ext);
                        }else if (qstart == 1 && subject_end_plus == last_subject_contig_len ){
                            //subject is extended by query
                            //take action accordingly
                            ContigExtension ext = new ContigExtension();
                            ext.contigA = last_query_index;
                            ext.contigB = last_subject_index;
                            ext.relative_position = Helper.POS_S_Extended_Q;
                            ext.orientation = '+';
                            ext.index = qend ; //length to be merged
                            extensions.add(ext);

                        }else if ( subject_start_plus == 1 && qend == last_query_contig_len ){
                            //query is extended by subject 
                            ContigExtension ext = new ContigExtension();
                            ext.contigA = last_query_index;
                            ext.contigB = last_subject_index;
                            ext.relative_position = Helper.POS_Q_Extended_S;
                            ext.orientation = '+';
                            ext.index = subject_end_plus; //length to be merged
                            extensions.add(ext);

                            //take action accordingly
                        }
                    }else{
                      //+-
                        int subject_start_minus = send;
                        int subject_end_minus  = sstart;
                        if(qstart == 1 && qend == last_query_contig_len && subject_start_minus>1){
                            //whole query contig is contained inside subject
                            //take action accordingly
                            ContigExtension ext = new ContigExtension();
                            ext.contigA = last_query_index;
                            ext.contigB = last_subject_index;
                            ext.relative_position = Helper.POS_Q_Inside_S;
                            //orientation and index is not important
                            extensions.add(ext);

                        }else if (subject_start_minus == 1 && subject_end_minus == last_subject_contig_len && qstart>1){
                            //whole subject contig is inside query contig
                            ContigExtension ext = new ContigExtension();
                            ext.contigA = last_query_index;
                            ext.contigB = last_subject_index;
                            ext.relative_position = Helper.POS_S_Inside_Q;
                            //orientation and index is not important
                            extensions.add(ext);
                        }else if (qstart == 1 && subject_start_minus == 1){ 
                            //RC of subject is extended by query
                            ContigExtension ext = new ContigExtension();
                            ext.contigA = last_query_index;
                            ext.contigB = last_subject_index;
                            ext.relative_position = Helper.POS_S_Extended_Q;
                            ext.orientation = '-'; //- here implies RC of subject
                            ext.index = qend; //length to be merged
                            extensions.add(ext);

                        }else if(qend == last_query_contig_len && subject_end_minus == last_subject_contig_len){
                            //query is extended by RC of subject 
                            ContigExtension ext = new ContigExtension();
                            ext.contigA = last_query_index;
                            ext.contigB = last_subject_index;
                            ext.relative_position = Helper.POS_Q_Extended_S;
                            ext.orientation = '-';//- here implies RC of subject
                            ext.index = qend - qstart + 1; //length to be merged
                            extensions.add(ext);                                
                        }
                        
                        
                    }
                }
                blast_line = reader.readLine();
                
            }//while loop in blast file lines
            System.out.println("Number of extensions detected: " + extensions.size());
            step = 10;
            /*
            For example if there is following extension between queries and subjects, then 2 matrixes will be filled like this:
                  ---------qi----                                 -------qm------
            ---sk-----        -------R.C. Sj---------------------------
            
            query_extensions                         subject_Extensions
index       0   |1  |2  |3              index        0  |1  |2  |3  
            --------------                           ------------------
            .............                            ..................
i'th row    j   |20 |0  |k                           ..................
            .............               j'th row     m  |15 |0  |i
m'th row    -1  |-1 |-1 |j                           ..................
            .............               k'th row     i  |30 |1  |-1
                                                     ..................
            
            In the third column 0 means - and 1 means +
            */
            //now we process extensions, There are some rules we apply:
            // (1) every single query or subject contig,  can only participate in just one chain/merging action
            // (2) if one subject or query is extended by multiple query/subject contigs we only consider longest one
            //In order to implement these rules, we prepare 2 arrays one for processed action one for extension 
            int subject_size = sfsi.size();
            int query_size = qfsi.size();
            
            System.out.println("Number of contigs in query fasta file: " + query_size);
            System.out.println("Number of contigs in subject fasta file: " + subject_size);
            int[][] query_extensions = new int[query_size][4];
            int[][] subject_extensions = new int[subject_size][4];
            boolean[] query_processed = new boolean[query_size]; //if a query contig is part of a valid chain then it is set to true
            boolean[] subject_processed = new boolean[subject_size];//if a subject contig is part of a valid chain then it is set to true
            //definition of 4 coulmns [][0]=extended by contig index/OR inside contig index  [][1]=merging length [][2]= merging orientation 1: forward 0: reverse  [][3]=extends (opposite of extended by) contig index (opposite of [][0]), used for back tracking   
            //now initialize arrays:
            for (int i=0;i<subject_size; i++){
                subject_processed[i] = false;
                subject_extensions[i][0]=-1;
                subject_extensions[i][1]=-1;
                subject_extensions[i][2]=-1;
                subject_extensions[i][3]=-1;
            }

            for (int i=0;i<query_size; i++){
                query_processed[i] = false;
                query_extensions[i][0]=-1;
                query_extensions[i][1]=-1;
                query_extensions[i][2]=-1;
                query_extensions[i][3]=-1;
                
            }
            step = 11;
            //now it is time to process extensions list and fill up the extensions arrays
            //For each query and subject there is only one row in the matrix so each contig can only have one action
            //rule: if  contig A is inside contig B, then if A is extended or extends contig C, then we ignore extension action, in other words inside action will override any other extended by action (Tehcnically it is not possible, because there is only one place to store 2 facts!!).
            //but if a contig A is inside contig B and contig B is extended by contig C, then both actions (inside and extended by) are stored (Thechnically it is also possible because these facts are stored in 2 different rows)
            for(ContigExtension ext: extensions){
                if (ext.relative_position == Helper.POS_Q_Inside_S){
                     //we also make sure backtrack flag is reset into -1
                    if (query_extensions[ext.contigA][0] > -1)
                       subject_extensions[query_extensions[ext.contigA][0]][3] = -1;

                    query_extensions[ext.contigA][0] = ext.contigB;
                    query_extensions[ext.contigA][1] = -1; //index 1 along with index 0 is used to distiguish inside action vs extension action, inside action is detected if and only if index 0 is not -1 and index 1 is -1 
                    query_extensions[ext.contigA][2] = -1;
                    query_extensions[ext.contigA][3] = -1;
                    //coulmns 1,2 and 3 will remain -1, so we will realize that the action is inside
                }else if (ext.relative_position == Helper.POS_S_Inside_Q){
                     //we also make sure backtrack flag is reset into -1
                    if (query_extensions[ext.contigA][0] > -1)
                       query_extensions[query_extensions[ext.contigA][0]][3] = -1;
                    
                    subject_extensions[ext.contigB][0] = ext.contigA;
                    subject_extensions[ext.contigB][1] = -1; //index 1 along with index 0 is used to distiguish inside action vs extension , inside action is is detected if and only if index 0 is not -1 and index 1 is -1 
                    subject_extensions[ext.contigB][2] = -1;
                    subject_extensions[ext.contigB][3] = -1;
                    //coulmns 1,2 and 3 will remain -1, so we will realize that the action is inside
                }else if (ext.relative_position == Helper.POS_Q_Extended_S){
                    //first make sure this query is not already inside a subject
                    //see if query is extended by something beforehand, if new extension is longer then we replace it
                    if (query_extensions[ext.contigA][0]==-1 ){
                        //it is safe to store the extension
                        query_extensions[ext.contigA][0]=ext.contigB;
                        subject_extensions[ext.contigB][3]=ext.contigA; // it is used for back tracking
                        //now store the length of merging
                        query_extensions[ext.contigA][1] = ext.index;
                        //now store orientation of merging
                        query_extensions[ext.contigA][2] = ext.orientation=='+'?1:0;
                    }else if (query_extensions[ext.contigA][1]!=-1 && query_extensions[ext.contigA][1]<ext.index){
                        // we get here when this query has already been extended by another subject contig with smaller ovelap length, so according to the rule we update it
                        // before updating query with new extension, we need to make sure we reset all backtrack flags
                        subject_extensions[query_extensions[ext.contigA][0]][3] = -1;

                        //now update query contig 
                        query_extensions[ext.contigA][0]=ext.contigB;
                        subject_extensions[ext.contigB][3]=ext.contigA; // it is used for back tracking
                        //now store the length of merging
                        query_extensions[ext.contigA][1] = ext.index;
                        //now store orientation of merging
                        query_extensions[ext.contigA][2] = ext.orientation=='+'?1:0;
                        
                    }
                   
                }else if (ext.relative_position == Helper.POS_S_Extended_Q){
                    //first make sure this subject is not already inside a query
                     //see if subject is extended by some query beforehand, if new extension is longer then we replace it
                    if (subject_extensions[ext.contigB][0]==-1 ){
                        subject_extensions[ext.contigB][0]=ext.contigA;
                        query_extensions[ext.contigA][3]=ext.contigB;
                         //now store the length of merging
                        subject_extensions[ext.contigB][1] =  ext.index;
                        //now store orientation of merging
                        subject_extensions[ext.contigB][2] = ext.orientation=='+'?1:0;
                    }else if (subject_extensions[ext.contigB][1] != -1 && subject_extensions[ext.contigB][1]<ext.index){
                        // we get here when this subject has already been extended by another contig with smaller ovelap length, so according to the rule we update it
                        // before updating subject with new extension, we need to make sure we reset all backtrack flags
                        query_extensions[subject_extensions[ext.contigB][0]][3] = -1;
                        //now updating subject_Extension row
                        subject_extensions[ext.contigB][0]=ext.contigA;
                        query_extensions[ext.contigA][3]=ext.contigB;
                         //now store the length of merging
                        subject_extensions[ext.contigB][1] =  ext.index;
                        //now store orientation of merging
                        subject_extensions[ext.contigB][2] = ext.orientation=='+'?1:0;
                        
                    }
                    
                }
            }
            // matrix has been loaded with correct data
            System.out.println("Matrix Extension has been Loaded.");
            //now we go through every contig in query, if it is part of chain we will detect the chain and merge them and report them as one single contig
           step = 12;
           ContigExtender2 ex = new ContigExtender2(query_extensions, subject_extensions, query_processed, subject_processed);
           
           long total_q_bp = 0;
           long merged_bp = 0; //total base pairs have been merged
           long written_bp = 0;
           int qcontig_idx = 0;
           int total_chains = 0;
           int special = 0;
           fsf = new FastaSequenceFile(new File(queryFile), true); 
           ReferenceSequence refseq = fsf.nextSequence();
           
           while (refseq != null){
                //String contig_seq = new String(refseq.getBases());
                //contigs.put(idx, new CompactSequence2(contig_seq, idx));
                //contigs.add(contig_seq);
               total_q_bp += refseq.length();
               if (query_extensions[qcontig_idx][0]==-1 && query_extensions[qcontig_idx][3]==-1){
                   //This means this query is not inside any subject and also is not part of any chain, so write it as it is (it is lonely contig)
                   
                    writer.write(">q_"+qcontig_idx+"\n");
                    writer.write(refseq.getBaseString() +"\n");  
                    written_bp += refseq.length();
               }else if (query_extensions[qcontig_idx][0]!=-1 && query_extensions[qcontig_idx][1] == -1 ){
                   //this means this is inside action, query is inside a subject, we dont write anything now, we will write subject when we process subject fasta file
                   //this query may extends a subject or be extended by a subject, but we are not interested in them, because inside action overrides other actions for this query
                   merged_bp += refseq.length();
               }else if (query_extensions[qcontig_idx][0]!=-1 || query_extensions[qcontig_idx][3]!=-1) {
                   //Here means query is part of possible extension chain (we may be able go foreward or backward)
                   //now our task is to discover the chain, before this, we need to know if this query contig has been part of a chain already.
                   if (!query_processed[qcontig_idx]){
                       // This means query has not already been part of any chain
                       List<ContigExtension> list = ex.getNextChain(qcontig_idx);
                       StringBuffer chainStr = new StringBuffer(); // this is for debugging purpose
                       try{
                           
                       if(list.size()==0){
                           // This happens when this contig can be extended by another contig that already been processed, 
                           //or query contig is extended by subject contig and subject is inside a bigger query contig, in both cases the list will be empty
                           //output contig as it is
                           step = 120;
                            writer.write(">q_"+qcontig_idx+"\n");
                            writer.write(refseq.getBaseString()+"\n");
                            written_bp += refseq.length();
                       }else{
                           total_chains++;
                           // merge elements of the extension chain
                           
                           //Following 2 lines are for debugging purpose
                           for(ContigExtension cext : list)
                               chainStr.append(cext.toString());
                           step = 121;
                           
                           //because we also go backward in the chain, then a chain may start with a subject
                           StringBuffer seq_buff = new StringBuffer();
                           StringBuffer name_buff = new StringBuffer();
                           String first_contig = null;
                           if (list.get(0).relative_position == Helper.POS_Q_Extended_S){
                               //first contig is query
                               first_contig = qifsi.getSequence(q_id_map.get(list.get(0).contigA)).getBaseString();
                               name_buff.append("q_" + list.get(0).contigA);
                           }else if (list.get(0).relative_position == Helper.POS_S_Extended_Q){
                               //first contig is subject
                               if (list.get(0).orientation == '+'){
                                   first_contig = sifsi.getSequence(s_id_map.get(list.get(0).contigB)).getBaseString();
                                   name_buff.append("s+_" + list.get(0).contigB);
                               }else if (list.get(0).orientation == '-'){
                                   first_contig = Helper.getRC(sifsi.getSequence(s_id_map.get(list.get(0).contigB)).getBaseString());
                                   name_buff.append("s-_" + list.get(0).contigB);                                   
                               }
                               
                               
                           }
                           seq_buff.append(first_contig);
                        step = 122;
                           for(ContigExtension ext : list){
                               if(ext.relative_position == Helper.POS_Q_Extended_S){
                                   
                                   String ssequence = sifsi.getSequence(s_id_map.get(ext.contigB)).getBaseString();
                                   if (ext.orientation == '+'){
                                       name_buff.append("_s+_" + ext.contigB);
                                       seq_buff.append(ssequence.substring(ext.index));
                                   }else if (ext.orientation == '-'){
                                       name_buff.append("_s-_" + ext.contigB);
                                       seq_buff.append(Helper.getRC(ssequence).substring(ext.index));
                                   }
                               }else if (ext.relative_position == Helper.POS_S_Extended_Q){
                                   name_buff.append("_q_" + ext.contigA);
                                   String qsequence = qifsi.getSequence(q_id_map.get(ext.contigA)).getBaseString();
                                   seq_buff.append(qsequence.substring(ext.index));                                    
                               }
                               
                               merged_bp += ext.index;
                           }
                           step = 123;
                           writer.write(">"+name_buff.toString()+"\n"); 
                           writer.write(seq_buff.toString() +"\n");
                           written_bp += seq_buff.length();

                       }//else chain has any elements
                       }catch(Exception exp){System.out.println(exp.toString() + " STEP:" + step + chainStr.toString());}
                   }//if query contig has not processed before
               }// if query contig can be part of a chain
               else if (!query_processed[qcontig_idx]){
                   //for any other reason if query has not been written and is unprocessed then we write it (but statistics shows it is always 0, so we never get here)
                   //we can remove this block safely anyway.
                    special++;
                    writer.write(">q_"+qcontig_idx+"*\n");
                    writer.write(refseq.getBaseString()+"\n");
                    written_bp += refseq.length();  
                    //
               }
                refseq = fsf.nextSequence();
                qcontig_idx++;
           }// while  query contig         

           System.out.println("Query contigs has been processed.");
           //now it is time to process subject contigs, for subject contigs we dont need to find chains, becuase all the subjects participated in any chain was already processed and merged and written before
           long total_s_bp = 0;
           int special2 = 0;
           int scontig_idx = 0;
           fsf.close();
           fsf = new FastaSequenceFile(new File(subjectFile), true); 
           refseq = fsf.nextSequence();           
           while (refseq != null){
                total_s_bp += refseq.length();
                if (subject_extensions[scontig_idx][0]==-1 && subject_extensions[scontig_idx][3]==-1){
                   //This means this subject contig is not inside any query and also is not part of any chain, so write it as it is                   
                    writer.write(">s_"+scontig_idx+"\n");
                    writer.write(refseq.getBaseString() +"\n");  
                    written_bp += refseq.length();
               }else if (subject_extensions[scontig_idx][0]!=-1 && subject_extensions[scontig_idx][1] == -1 ){
                   //this means this is inside action, subject is inside a query, we dont write anything, because it has been written when processed query contigs
                   merged_bp += refseq.length();
               }else if (!subject_processed[scontig_idx]){
                   // for any other reason, that subject has not been in a chain then we write it
                   // The reason we get here is simply because this subject is extended by a query that query is inside another subject, so there would be no chain.
                   // But because we discovered all the chains by starting from a query contig, so these subjects never get a chance being called by getnextChainMethod, so we are left with a bunch of subjects that their index 0 and 3 are not zero, and also their processed flag is still false. therefore we still need to write them.
                   // This situation also occurs to query contigs, but because all query contigs go through getNextChain method, so their processed flag will set to True and they are written as they are to output
                    writer.write(">s_"+scontig_idx+"*\n");
                    writer.write(refseq.getBaseString() +"\n");  
                    written_bp += refseq.length();    
                    special2++;
                    //following 2 lines are for debugging
                    //System.out.print(scontig_idx + ":{" + subject_extensions[scontig_idx][0]+ ", " + subject_extensions[scontig_idx][1] + ", " + subject_extensions[scontig_idx][2] + ", " + subject_extensions[scontig_idx][3]+ "}->");
                    //System.out.println( "{" + query_extensions[subject_extensions[scontig_idx][0]][0]+ ", " + query_extensions[subject_extensions[scontig_idx][0]][1] + ", " + query_extensions[subject_extensions[scontig_idx][0]][2] + ", " + query_extensions[subject_extensions[scontig_idx][0]][3]+ "}");
               }
                refseq = fsf.nextSequence();
                scontig_idx++;
           }
            
           System.out.println("Total Query contigs bp:" + total_q_bp);
           System.out.println("Total Subject contigs bp:" + total_s_bp);
           System.out.println("Total Chains:" + total_chains);
           System.out.println("Total new contigs written bp:" + written_bp);
           System.out.println("Total merged bp:" + merged_bp);
           //System.out.println("Total Query Contigs with Special cases:" + special);
           //System.out.println("Total Subject Contigs with Special cases:" + special2);
           
        } catch (Exception ex) {System.out.println("Merging query and subject fasta files " + ex.getMessage() + " at step "+step + blast_line);}
        finally{
            try {
                
                writer.close();
                qifsi.close();
                sifsi.close();
                fsf.close();
            } catch (Exception ex) {}
        } 
        
    }
    
    
    
     public void intersect(String refFile, String refFileIndex,  String popseqFile, double min_identity_percent, int min_alignment_length){
        // This method merges query fasta file contigs with subject fasta contigs for example if we blast 7A (query) to 7B (subject), and create a new merged fasta file
         // Potenitally, we are looking for query contigs that can uniquely extend a subject contig from either of its ends
         //Query:                |----------------|                      |--------------------------|
         //Subject:     |-----------|   |---|  |---------------------------|
         // refFile file contains query and subject contigs both. It is used to estimate length of each contig
         // blast result file is of the format 6  qseqid qstart qend sseqid sstart send length pident mismatch gapopen gaps
         //Conditions: Every query can extend multiple subjects, but each subject must be extended by maximum 2 subjects (one subject from each side)
        Writer writer = null;
        BufferedReader popseqReader = null;
        Map<String, Double> popseqMap = new HashMap(); // a map of contig names -> cM

        try {
            // first read ref fasta file
            FastaSequenceIndex fsi = new FastaSequenceIndex(new File(refFileIndex));
            //IndexedFastaSequenceFile ifsi =new IndexedFastaSequenceFile(new File(refFile) , fsi);
            // now read popseq file:
            popseqReader = new BufferedReader(new FileReader(popseqFile));
            popseqReader.readLine(); // read header
            String line = popseqReader.readLine();
            while (line != null && !line.isEmpty()){
                String[] arr = line.split("\t");
                String contig = arr[0];
                //String chromosome = arr[1].toLowerCase();                    
                double position = Double.parseDouble(arr[2]);                   
                popseqMap.put(contig, position);                    
                line = popseqReader.readLine();
            } //while                

            long total_alignment_len = 0;
            int extension_num_l = 0;   // number of cases, where a query contig can extend a subject contig from left end
            int extension_num_r = 0;   // number of cases, where a query contig can extend a subject contig from right end
            int extension_num_lr = 0;   // number of cases, where a query contig can extend a subject contig from both ends
            long extension_len = 0; 
            int extension_homeless = 0; // where where query or subject contig position is unknown
            int extension_matched = 0;  // where query and subject positions are the same
            long pbIndex = 0;
            double lastpercent = 0.0;
            long fileSz = new File(this.fileName).length();
            String blast_line = reader.readLine();
            while (blast_line != null && !blast_line.isEmpty()){
                pbIndex += blast_line.getBytes().length; // we add number of bytes read
                double newpercent = ((double)pbIndex/fileSz)*100;
                if (newpercent-lastpercent>=1){
                    lastpercent = newpercent;
                    System.out.println("progress %"+(int)lastpercent);
                }
                
                String[] arr = blast_line.split("\\s");
                String q_contig = arr[0].trim();
                int qstart = Integer.parseInt(arr[1]);
                int qend = Integer.parseInt(arr[2]);
                String s_contig = arr[3].trim();
                int sstart = Integer.parseInt(arr[4]);
                int send = Integer.parseInt(arr[5]);
                int alignment_len = Integer.parseInt(arr[6]);
                double percent = Double.parseDouble(arr[7]);
                
                if (percent >= min_identity_percent  && alignment_len >= min_alignment_length){
                    total_alignment_len += alignment_len;
                    int query_len = (int)fsi.getContigSize(q_contig);
                    int subject_len = (int)fsi.getContigSize(s_contig);
                    if (sstart < send){
                        // this means orientation of alignment from query to the subject is the same
                        if( (query_len > alignment_len) && ((qstart==1 && send==subject_len) || (sstart==1 && qend==query_len) || (sstart==1 && send==subject_len)) ){
                            // This means query contig has got potential to extend subject contig from at least one of its ends
                            if (qstart==1 && send==subject_len)
                                extension_num_l++;
                            else if (sstart==1 && qend==query_len)
                                extension_num_r++;
                            else if (sstart==1 && send==subject_len)
                                extension_num_lr++;
                            // now we need to extract their cM position and see if they are the same
                           
                            extension_len += query_len - alignment_len;  // number of bp of extension
                            Double q_cM = popseqMap.get(q_contig);
                            Double s_cM = popseqMap.get(s_contig);
                            if (q_cM ==null || s_cM==null)
                                extension_homeless++; // position of query or subject is unknown
                            else if (q_cM.doubleValue() == s_cM.doubleValue())
                                extension_matched++;
                            else
                                System.out.println(String.format("%s:%.2f vs. %s:%.2f", q_contig, q_cM, s_contig, s_cM));
                        }
                    }else{
                        // This means query and subject orientations are different
                        if( (query_len > alignment_len) && ((qend==query_len && sstart==subject_len) || (qstart==1 && send==1) || (sstart==subject_len && send==1)) ){
                            // This means query contig has got potential to extend subject contig from at least one of its ends
                            if (qend==query_len && sstart==subject_len)
                                extension_num_l++;
                            else if (qstart==1 && send==1)
                                extension_num_r++;
                            else if (sstart==subject_len && send==1)
                                extension_num_lr++;
                            // now we need to extract their cM position and see if they are the same
                           
                            extension_len += query_len - alignment_len;  // number of bp of extension
                            Double q_cM = popseqMap.get(q_contig);
                            Double s_cM = popseqMap.get(s_contig);
                            if (q_cM ==null || s_cM==null)
                                extension_homeless++; // position of query or subject is unknown
                            else if (q_cM.doubleValue() == s_cM.doubleValue())
                                extension_matched++;
                            else
                                System.out.println(String.format("%s:%.2f vs. %s:%.2f", q_contig, q_cM, s_contig, s_cM));
                        }
                        
                    }
                }
                
                blast_line = reader.readLine();
            }//while

            System.out.println(String.format("Alignment Length:%d Extension Left:%d Extension Right:%d Extension Both:%d Extension Sum:%d Extension Length(bp):%d Extension Matched:%d Extension Homeless:%d", total_alignment_len, extension_num_l, extension_num_r, extension_num_lr, extension_num_l+extension_num_r+extension_num_lr, extension_len, extension_homeless, extension_matched));
            // now we repot all the statistics
        } catch (Exception ex) {System.out.println(ex.getMessage());}
        finally{
            
        }    

    }
    
}
