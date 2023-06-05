/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package mfruzan.common;

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
import java.io.RandomAccessFile;
import java.io.Reader;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.GZIPInputStream;
import mfruzan.algorithm.BlastHSP;
import mfruzan.algorithm.MAFSLine;
import mfruzan.sina.Start;

/**
 *
 * @author mario.fruzangohar
 * This class reads Multiple alignment File format, it also contains reader and writer methods, thats why it looks weird!
 */
public class MAFProcessor {
    BufferedReader buf = null;
    //TreeMap<Integer, BlastHSP> query_positionMap_plus = new TreeMap(); // a map of start position of query -> HSP; if strand is + ; TreeMap naturally ordered by its key
    //TreeMap<Integer, BlastHSP> query_positionMap_minus = new TreeMap();// a map of start position of query -> HSP; if strand is - ; TreeMap naturally ordered by its key
    // MAF results are not sorted based on the score , first we order them
    
    TreeMap<Double, BlastHSP> scoreMap_plus = new TreeMap(); // a map of Score of HSP -> HSP; if strand is + ; TreeMap naturally ordered by its key
    TreeMap<Double, BlastHSP> scoreMap_minus = new TreeMap();// a map of Score of HSP -> HSP; if strand is - ; TreeMap naturally ordered by its key
    
    public List<BlastHSP> all_hsp = new ArrayList();
    public String fName = null;
    public MAFProcessor(String fileName){        
        try{
            if (fileName.endsWith("maf.gz")){
                InputStream filein = new FileInputStream(fileName);
                InputStream gzipin = new GZIPInputStream(filein);
                Reader ir = new InputStreamReader(gzipin);
                this.buf = new BufferedReader(ir);
            }else if (fileName.endsWith(".maf")|| fileName.endsWith(".maff")){
                this.buf = new BufferedReader(new FileReader(fileName));
            }
            fName = new File(fileName).getName();
        }catch(Exception ex){System.out.println(ex.getMessage());}
    } 
    
    public MAFProcessor(){}
    //---------------------------------------------------
    public void load(){
        // this method might need a lot of memory because loads whole file into memory to sort based on HSP score
        System.out.println("Start Loading ...");
        BlastHSP hsp = getNext();
        while(hsp != null){
           all_hsp.add(hsp);
           if (hsp.orientation == 1){
               // positive +/+
               //query_positionMap_plus.put(hsp.q_start, hsp);
               scoreMap_plus.put(hsp.score, hsp);
           }else{
               //Negative +/-
               scoreMap_minus.put(hsp.score, hsp);
           }
           //--------------
           hsp = getNext();
        } 
        
        System.out.println(all_hsp.size() + " was loaded.");
    }
     public void loadBinary(){
         // assumes each MAF entry has only 2 SLines (subject + query) such as lastz output
        // this method uses less memory than load()
        System.out.println("Start Loading ...");
        BlastHSP hsp = getNext();
        while(hsp != null){
           all_hsp.add(hsp);
           hsp = getNext();
        } 
        
        System.out.println(all_hsp.size() + " was loaded.");
    }
     public void loadMulti(){
        // loading multi line MAF file (means each MAF entry might have got more than 2 S line)
        System.out.println("Start Loading ...");
        BlastHSP hsp = getNextMulti();
        while(hsp != null){
           all_hsp.add(hsp);
           hsp = getNextMulti();
        } 
        
        System.out.println(all_hsp.size() + " was loaded.");
    }     
    public BlastHSP searchQuery(String cntgName, int start, int end){
        for(BlastHSP hsp : all_hsp){
            if (hsp.q_contig.equals(cntgName) && hsp.q_start== start && hsp.q_end == end)
                return hsp.clone();
        }
        
        return null;
    }
    public void writeAsFasta(String fastaFile){
        // writes each HSP sequence as fasta file  (writes the squence that is smaller)
        FileWriter writer = null;
        try{
            writer = new FileWriter(fastaFile);
            //write ++
            Map<Double, BlastHSP> map_plus = scoreMap_plus.descendingMap(); // reverse the order from high to low
            int cnt_counter = 1;
            for(Map.Entry<Double, BlastHSP> entry : map_plus.entrySet()) {
                   BlastHSP hsp = entry.getValue();
                   String seq = null;
                   if (hsp.s_len<hsp.q_len)
                       seq = hsp.s_alignment_block;
                   else
                       seq= hsp.q_alignment_block;
                   //now write it as a contig
                   writer.write(String.format(">cntg%d\n", cnt_counter++));
                   writer.write(seq.replaceAll("-", "")+"\n");                   
            } 
            
            // now write +-
            Map<Double, BlastHSP> map_minus = scoreMap_minus.descendingMap();
            for(Map.Entry<Double, BlastHSP> entry : map_minus.entrySet()) {                    
               BlastHSP hsp = entry.getValue();
               String seq = null;
               if (hsp.s_len<hsp.q_len)
                   seq = hsp.s_alignment_block;
               else
                   seq= hsp.q_alignment_block;
               //now write it as a contig
               writer.write(String.format(">cntg%d\n", cnt_counter++));
               writer.write(seq.replaceAll("-", "")+"\n");                   
               
            }                 

            
            
            
        } catch (Exception ex) {System.out.println("writeAsFasta:" + ex.getMessage() );}
        finally{
            try {
                writer.close();
            } catch (Exception ex) {}
        }             
    }
    
    
    public void writeAsFasta2(String fastaFile){
        // retired and replaced by extractNonOverlap
        // lastz output is ordered by query conig ID 
        // writes each HSP sequence as fasta file  (writes the query always)(16/12/2019)
        FileWriter writer = null;
        try{
            writer = new FileWriter(fastaFile);
            String last_contig = null;
            TreeMap<Double, BlastHSP> map = new TreeMap(); 
            for(BlastHSP hsp : all_hsp) {
                   
               //now write it as a contig
               if(last_contig !=null && !last_contig.equals(hsp.q_contig) ){
                   System.out.print(last_contig + ":"  +map.size());
                    //Map<Double, BlastHSP> map_contig = map.descendingMap();
                    // now map_contig contains HSP from high score to low score
                    // we define another List and we put non-overlapping elements of map_contig into this new List
                    List<BlastHSP> list = new ArrayList();
                    for(Map.Entry<Double, BlastHSP> entry : map.descendingMap().entrySet()) {                    
                       BlastHSP aHsp = entry.getValue();
                       aHsp.addToListNoOverlap(list, 5);
                    } 
                    System.out.println(" write:" + list.size());
                    for(BlastHSP bHSP : list){
                       writer.write(String.format(">%s!%d!%d\n", bHSP.q_contig, bHSP.q_start, bHSP.q_end));
                       writer.write(bHSP.q_alignment_block.replaceAll("-", "")+"\n");
                    }

                   map.clear();
               }
               last_contig = hsp.q_contig;
               map.put(hsp.score, hsp);
                               
            } // for all HSPs
            //write the last query contig
            System.out.println(last_contig + ":"  +map.size());
            List<BlastHSP> list = new ArrayList();
            for(Map.Entry<Double, BlastHSP> entry : map.descendingMap().entrySet()) {                    
               BlastHSP aHsp = entry.getValue();
               aHsp.addToListNoOverlap(list, 5);
            } 
            System.out.print(" write:" + list.size());
            for(BlastHSP bHSP : list){
               writer.write(String.format(">%s!%d!%d\n", bHSP.q_contig, bHSP.q_start, bHSP.q_end));
               writer.write(bHSP.q_alignment_block.replaceAll("-", "")+"\n");
            }
            
        } catch (Exception ex) {System.out.println("writeAsFasta2:" + ex.getMessage() );}
        finally{
            try {
                writer.close();
            } catch (Exception ex) {}
        }             
    }
    public void extractNonOverlapQuery(String fastaFile, String outMAFFile, int min_length){
        // This method assumes aligner output MAF file reports MAF entries into order of query contig
        //Updated 26/8/2021, can accept output of lastz ,GSAlign and Minimap2
        // lastz output is ordered by query conig ID , GSAlign output also ordered by query contig ID
        // writes each non-overlapping HSP query sequence as fasta file  (writes the query always)(12/3/2020)
        //Also new filtered MAF file is wrriten
        //It also generated min, max, mean, N50 and Sum of size of query contigs (so save us running fastainfo job)
        int step = 1;
        FileWriter writer = null;
        FileWriter mafw = null;
        List<Integer> q_contig_lengths = new ArrayList();
        int total_len = 0;
        int total_len_ungapped = 0;
        int total_rejected_len = 0;
        try{
            writer = new FileWriter(fastaFile);
            mafw = new FileWriter(outMAFFile);
            String last_contig = null;
            TreeMap<Double, BlastHSP> map = new TreeMap(); // holds entries of a querycontig a treemap can be ordered by its key
            BlastHSP hsp = getNext();
            // in the last iteration of below loop hsp is null, but lastContigWritten will set to true and the loop will exit, because hsp is null and lastContigWritten is set to true
            boolean lastContigWritten = false;
            step = 2;
            while(hsp != null || !lastContigWritten){
                if(hsp == null || (last_contig !=null && !last_contig.equals(hsp.q_contig)) ){
                   //System.out.print(last_contig + ":"  +map.size());
                    //Map<Double, BlastHSP> map_contig = map.descendingMap();
                    // now map_contig contains HSP from high score to low score
                    // we define another List and we put non-overlapping elements of map_contig into this new List
                    List<BlastHSP> list = new ArrayList();
                    for(Map.Entry<Double, BlastHSP> entry : map.descendingMap().entrySet()) {                    
                       BlastHSP aHsp = entry.getValue();
                       //boolean added = aHsp.addToListNoOverlap(list, 5000000);
                       //System.out.println("process:");
                       //aHsp.print4debug();
                       boolean log = false;
                       //if (aHsp.q_contig.equals("JAIRFR010000010.1") && aHsp.q_start==977981 && aHsp.q_end==978756)
                       //    log = true;
                       boolean added = aHsp.addToList(list, log);
                       if (!added)
                           total_rejected_len += aHsp.q_len;
                       else{
                          //System.out.println("added:");
                          //aHsp.printQueryAsMAF();
                          //aHsp.printSubjectAsMAF();
                       }
                       
                    } 
                    step = 3;
                    //System.out.println(" Number of HSP write:" + list.size());
                    for(BlastHSP bHSP : list){
                        if (bHSP.q_alignment_block.length() >= min_length){
                           q_contig_lengths.add(bHSP.q_alignment_block.length());
                           total_len += bHSP.q_alignment_block.length();
                           //because GSAlign will add qry. and ref. to the begining of contig names, we just get rid of them
                           if (bHSP.q_contig.startsWith("qry.")){
                               bHSP.q_contig = bHSP.q_contig.substring(4);
                               bHSP.queryLine = bHSP.queryLine.replaceFirst("qry\\.", "");
                           }
                           if (bHSP.s_contig.startsWith("ref.")){
                               bHSP.s_contig = bHSP.s_contig.substring(4);
                               bHSP.subjectLine = bHSP.subjectLine.replaceFirst("ref\\.", "");
                           }
                           
                           writer.write(String.format(">%s!%d!%d\n", bHSP.q_contig, bHSP.q_start, bHSP.q_end));
                           //always change query string to its original form
                           String qseq = bHSP.q_alignment_block.replaceAll("-", "");
                           total_len_ungapped += qseq.length();
                           if(!bHSP.q_pos_strand)
                               qseq = Helper.getRC(qseq);// remember that getRC change all the characters to uppercase
                           writer.write(qseq+"\n");
                          // now write MAF file
                          writeNext(mafw, bHSP);
                           //mafw.write(bHSP.scoreLine + "\n");
                           //mafw.write(bHSP.subjectLine + "\n");
                           //mafw.write(bHSP.queryLine + "\n");
                           //mafw.write("\n");
                        }
                       
                    }

                    step = 4;
                   map.clear();
                   //---------------
                   if (hsp == null)
                       lastContigWritten = true;
                   //---------------
               }
                step = 5;
               if (hsp != null){
                   last_contig = hsp.q_contig;
                   map.put(hsp.score, hsp);               
                   hsp = getNext();
               }
               step = 6;
            }// while we have more hsp
            
        Collections.sort(q_contig_lengths); // sort ascending order    
        step = 7;
        //System.out.println(contig_lengths.get(0));
        //System.out.println(contig_lengths.get(contig_lengths.size()-1));
        int n50_len = 0;
        long n50 = (int) (0.5*total_len);
        long cumulative_len = 0;
        try{
            for(int i=q_contig_lengths.size()-1; i>=0; i--){
                cumulative_len += q_contig_lengths.get(i);
                if (cumulative_len>=n50){
                    n50_len = q_contig_lengths.get(i);
                    break;
                }
            }           
        }catch(Exception ex){}
        step = 8;
        //System.out.println("Total Rejected(bp) Query HSP, because of overlap:"+ total_rejected_len);
        //System.out.println(String.format("Total HSP:%d, Total Alignment Length(bp) :%d, Max alignment Length(bp):%d, Min alignment Length(bp):%d, N50:%d",  q_contig_lengths.size(), total_len, q_contig_lengths.get(q_contig_lengths.size()-1), q_contig_lengths.get(0),  n50_len));    
        if (q_contig_lengths.size()>0)            
           System.out.println(String.format("Total HSP:%d, Total Rejected(bp):%d, Total Alignment Length(bp) :%d, Total ungapped Alignment Length(bp) :%d ,Max alignment Length(bp):%d, Min alignment Length(bp):%d, N50:%d",  q_contig_lengths.size(), total_rejected_len, total_len, total_len_ungapped,q_contig_lengths.get(q_contig_lengths.size()-1), q_contig_lengths.get(0),  n50_len));    
        } catch (Exception ex) {System.out.println("Error at extractNonOverlapQuery:" + ex.getMessage() + " At step " +step);}
        finally{
            try {
                writer.close();
                mafw.close();
            } catch (Exception ex) {}
        }             
    }
    public void uniqueQuery(String newMAFF){
        //retired, the logic has been merged into extractNonOverlapQuery method
        // lastz output is ordered by query conig ID 
        // writes each HSP sequence as avoid ovelap in query regions(9/1/2020)
        FileWriter writer = null;
        try{
            writer = new FileWriter(newMAFF);
            String last_contig = null;
            TreeMap<Double, BlastHSP> map = new TreeMap(); 
            for(BlastHSP hsp : all_hsp) {
                   
               //now write it as a contig
               if(last_contig !=null && !last_contig.equals(hsp.q_contig) ){
                   System.out.print(last_contig + ":"  +map.size());
                    //Map<Double, BlastHSP> map_contig = map.descendingMap();
                    // now map_contig contains HSP from high score to low score
                    // we define another List and we put non-overlapping elements of map_contig into this new List
                    List<BlastHSP> list = new ArrayList();
                    for(Map.Entry<Double, BlastHSP> entry : map.descendingMap().entrySet()) {                    
                       BlastHSP aHsp = entry.getValue();
                       aHsp.addToListNoOverlap(list, 0);
                    } 
                    System.out.println(" write:" + list.size());
                    for(BlastHSP bHSP : list){
                       writer.write(bHSP.scoreLine + "\n");
                       writer.write(bHSP.subjectLine + "\n");
                       writer.write(bHSP.queryLine + "\n");
                       writer.write("\n");
                    }

                   map.clear();
               }
               last_contig = hsp.q_contig;
               map.put(hsp.score, hsp);
                               
            } // for all HSPs
            //write the last query contig
            System.out.println(last_contig + ":"  +map.size());
            List<BlastHSP> list = new ArrayList();
            for(Map.Entry<Double, BlastHSP> entry : map.descendingMap().entrySet()) {                    
               BlastHSP aHsp = entry.getValue();
               aHsp.addToListNoOverlap(list, 0);
            } 
            System.out.print(" write:" + list.size());
            for(BlastHSP bHSP : list){
               writer.write(bHSP.scoreLine + "\n");
               writer.write(bHSP.subjectLine + "\n");
               writer.write(bHSP.queryLine + "\n");
               writer.write("\n");
            }
            
        } catch (Exception ex) {System.out.println("unique query:" + ex.getMessage() );}
        finally{
            try {
                writer.close();
            } catch (Exception ex) {}
        }             
    }
    public void getQueryAsFasta(String fastaFile){
        // writes each HSP sequence as fasta file  (writes the squence that is smaller)
        FileWriter writer = null;
        try{
            writer = new FileWriter(fastaFile);
            //write ++
            Map<Double, BlastHSP> map_plus = scoreMap_plus.descendingMap(); // reverse the order from high to low
            int cnt_counter = 1;
            for(Map.Entry<Double, BlastHSP> entry : map_plus.entrySet()) {
                   BlastHSP hsp = entry.getValue();
                   String seq = null;
                   if (hsp.s_len<hsp.q_len)
                       seq = hsp.s_alignment_block;
                   else
                       seq= hsp.q_alignment_block;
                   //now write it as a contig
                   writer.write(String.format(">cntg%d\n", cnt_counter++));
                   writer.write(seq.replaceAll("-", "")+"\n");                   
            } 
            
            // now write +-
            Map<Double, BlastHSP> map_minus = scoreMap_minus.descendingMap();
            for(Map.Entry<Double, BlastHSP> entry : map_minus.entrySet()) {                    
               BlastHSP hsp = entry.getValue();
               String seq = null;
               if (hsp.s_len<hsp.q_len)
                   seq = hsp.s_alignment_block;
               else
                   seq= hsp.q_alignment_block;
               //now write it as a contig
               writer.write(String.format(">cntg%d\n", cnt_counter++));
               writer.write(seq.replaceAll("-", "")+"\n");                   
               
            }                 

            
            
            
        } catch (Exception ex) {System.out.println("writeAsFasta:" + ex.getMessage() );}
        finally{
            try {
                writer.close();
            } catch (Exception ex) {}
        }             
    }
    
    public void writeHomologFile(String homologFile, String filter){
        // if filter is null then there is no filter, if it is + then Map_plus is considered if it is - then Map_minus is considered
        //if it is +/- 
        // writes HSP as tab separated file, each line is : q_contig q_start q_end s_contig s_start s_end orientation  query_sequence     subject_sequence(always s_tart<e_nd)
        // output file is sorted based on q start
        LinkedList<BlastHSP> composition = new LinkedList(); //we select some HSP in some specified order and put then into this linked list, to compose target region
        FileWriter writer = null;
        
        try{ 
            writer = new FileWriter(homologFile);
            if (filter!= null && filter.equals("+")){
                Map<Double, BlastHSP> map_plus = scoreMap_plus.descendingMap(); // reverse the order from high to low
                for(Map.Entry<Double, BlastHSP> entry : map_plus.entrySet()) {
                   BlastHSP hsp = entry.getValue();
                   BlastResult.insertIntoBlastHSPList(composition, hsp);
                } 
            }else if (filter!= null && filter.equals("-")){
                Map<Double, BlastHSP> map_minus = scoreMap_minus.descendingMap();
                for(Map.Entry<Double, BlastHSP> entry : map_minus.entrySet()) {                    
                   BlastHSP hsp = entry.getValue();
                   BlastResult.insertIntoBlastHSPList(composition, hsp);
                }                 
            }else if (filter!= null && filter.equals("+-")){
                Map<Double, BlastHSP> map_plus = scoreMap_plus.descendingMap(); // reverse the order from high to low
                Map<Double, BlastHSP> map_minus = scoreMap_minus.descendingMap();
                // Proiroty is with positive HSP, if after adding all of them any space left between query regins, then negative HSP will try to fill them.
                for(Map.Entry<Double, BlastHSP> entry : map_plus.entrySet()) {
                   BlastHSP hsp = entry.getValue();
                   BlastResult.insertIntoBlastHSPList(composition, hsp);
                } 
                for(Map.Entry<Double, BlastHSP> entry : map_minus.entrySet()) {                    
                   BlastHSP hsp = entry.getValue();
                   BlastResult.insertIntoBlastHSPList(composition, hsp);
                }                 
                
            }else if (filter!= null && filter.equals("-+")){
                // Proiroty is with nageative HSP, if after adding all of them any space left between query regins, then positive HSP will try to fill them.
                Map<Double, BlastHSP> map_plus = scoreMap_plus.descendingMap(); // reverse the order from high to low
                Map<Double, BlastHSP> map_minus = scoreMap_minus.descendingMap();

                for(Map.Entry<Double, BlastHSP> entry : map_minus.entrySet()) {                    
                   BlastHSP hsp = entry.getValue();
                   BlastResult.insertIntoBlastHSPList(composition, hsp);
                }                 
                
                for(Map.Entry<Double, BlastHSP> entry : map_plus.entrySet()) {
                   BlastHSP hsp = entry.getValue();
                   BlastResult.insertIntoBlastHSPList(composition, hsp);
                } 
                
            }
            
            // now print composition list
            for(BlastHSP hsp : composition)
                writer.write(String.format("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\n", hsp.q_contig, hsp.q_start, hsp.q_end, hsp.s_contig, hsp.s_start, hsp.s_end, hsp.orientation==1?"+":"-", hsp.q_alignment_block, hsp.s_alignment_block));                           
            
        } catch (Exception ex) {System.out.println("writeHomologFile" + ex.getMessage() );}
        finally{
            try {
                writer.close();
            } catch (Exception ex) {}
        } 
            
        
    }
    
    public void generateMSA(String srcDir, String maffFiles, String fafFile, String mafFile, String nameFile){
       
        // This method generated fa alignment file (this needs enough memory, to attach all the small sequnces relarated to one specie in one big sequnce)
        // the second output is maf file 
        FileWriter fafWriter = null;
        FileWriter mafWriter = null;
        BufferedReader nameReader = null;
        
        try {
            fafWriter = new FileWriter(fafFile);         
            mafWriter = new FileWriter(mafFile);
            System.out.println("Processing Maff Files at : " + srcDir );
            String[] files = maffFiles.split(",");

            List<String> names = new ArrayList();
            List<MAFProcessor> mafList = new ArrayList();
            // the following 2 list have size of mafList + 1, because for n genome we will have n-1 maf files
            List<StringBuffer> buffList = new ArrayList(); // this List needs enough memory 
            List<StringBuffer> partialBuffList = new ArrayList();
            // add first element
            buffList.add(new StringBuffer());
            partialBuffList.add(new StringBuffer());
            // we add maf files from last to first one into mafList
            for(int i=0; i<files.length; i++){
                // this part can be memory intensive, we need to have space for all maf files in the memory
                MAFProcessor maf = new MAFProcessor(srcDir + "/" +files[i]);
                maf.loadBinary();
                mafList.add(maf);
                buffList.add(new StringBuffer());
                partialBuffList.add(new StringBuffer());
            }
            //when we get here buffList and partialBuffList both have entries equal to number of Genomes
            nameReader = new BufferedReader(new FileReader(nameFile));
            String nameLine = nameReader.readLine();
            while(nameLine != null && !nameLine.isEmpty()){
                String name = nameLine.split("\\s")[0];
                names.add(name);
                nameLine = nameReader.readLine();
            }
            nameReader.close();
            if (names.size() != files.length+1)
                throw new Exception("Number of names in does not match number of MAF files");
            //co ordinates of subjects are absulute but query is always relative to previous step for example if query is something like : qcontig!20!15!2!4!1!2
            //the absolute coordinates of query is calculated recursively at each step co-ordinates are added to start co-ordinate of previous step: starting from end 1+2-1=2 and 2+2-1=3, then 2 and 3 are going to be added to the start of step before: 2+20-1=21 3+20-1=22  so absulute values are 21 until 22 of the qcontig

            for(BlastHSP hsp : mafList.get(0).all_hsp){
                //last maf file is the first in the list
                System.out.println("processing Query  :" + hsp.q_contig + ":" + hsp.q_start + "-" + hsp.q_end);
                //partialBuffList are restet for each query region
                //empty all StringBuffers inside partialBuffList  NQWZ01000001.1!6!121918!1!121913!1!121913!36790!121913!1!85124!1!85124!1!85124:22106-85124
                if (hsp.q_contig.equals("ST00Arg_chr_1!2300082!2384753!21455!49311!9796!27857!6720!18062!1!11343!1!11343!1317!11343!1!10023") && hsp.q_start==1 && hsp.q_end==4041)
                    hsp.q_start=1;
                for(StringBuffer bf : partialBuffList)
                    bf.delete(0, bf.length());
                //for the last maf file
                //String first_query_block = hsp.q_alignment_block;
                //first write query
                partialBuffList.get(0).append(hsp.q_alignment_block);
                // a MAF entry will be created
                BlastHSP mafEntry = new BlastHSP();
                mafEntry.mafLines = new ArrayList();
                mafEntry.mafLines.add(new MAFSLine()); // this is only a place holder for query, its values will be calculated and replaced after following loop
                
                //the following 5 variables are responsible to calculate absolute co-ordinates of the query , (calculation is explained above), original query contig name and size
                int absolute_qcontig_start =  hsp.q_start;
                int absolute_qcontig_end = hsp.q_end;
                String absolute_qcontig_name = hsp.q_contig;
                int absolute_qcontig_size = hsp.q_contig_size;
                boolean absolute_qcontig_pos_strand = hsp.q_pos_strand;
                byte absolute_qorientation = 1;
                
                //then the subject
                
                int offset_start = 1;
                int offset_end = hsp.q_end - hsp.q_start + 1;                                 

                String searching_contig = hsp.q_contig;

                int searching_end = hsp.q_end;
                int searching_start = hsp.q_start;

                // we creates a list of StringBuffer and for each extracted subject we insert string into it
                try{
                    for(int i = 0; i<mafList.size(); i++){
                        MAFProcessor maf = mafList.get(i);
                        BlastHSP aHSP = maf.searchQuery(searching_contig, searching_start, searching_end);

                        absolute_qcontig_name = aHSP.q_contig;
                        absolute_qcontig_size = aHSP.q_contig_size;
                        absolute_qcontig_pos_strand = aHSP.q_pos_strand;
                        absolute_qorientation = aHSP.orientation;

                        if (aHSP==null)
                           throw new Exception("HSP not found:" + searching_contig + " From " + searching_start + " To "+ searching_end + " in maf file " + maf.fName);
                        //read from search_start to search_end
                        //String subject_part = null;
                        //if (offset_start == 0)
                        //   subject_part =aHSP.getSubjectPart(offset_start+1, offset_end+1);
                        //else
                          //  subject_part =aHSP.getSubjectPart(offset_start, offset_end, first_query_block);
                      aHSP.alignSubjectPart3(offset_start, offset_end, i+1, partialBuffList, false);
                        //buffList.get(i+1).append(subject_part);

                        if(partialBuffList.get(i+1).length()==0)
                            throw new Exception(""); // this will cause whole alignment for current query will be discarded

                        // now we have information to write MAFSLine for this subject genome, at this stage alignment_block is not valid and it might be updated
                        // This constructor read as BlastHSP format
                        MAFSLine line = new MAFSLine(aHSP.s_contig, aHSP.s_start, aHSP.s_end, aHSP.s_pos_strand, aHSP.s_contig_size, "");
                        mafEntry.mafLines.add(line);
                        //update search_start, search_end, search_contig. offset_start, offset_end
                        if (aHSP.q_contig.indexOf('!')>0){
                            // if there is one level above; we dont enter this 'if' condition for the last MAF file
                            offset_start = aHSP.q_start + offset_start -1;
                            offset_end = aHSP.q_start  + offset_end - 1;
                            String[] arr = aHSP.q_contig.split("!");
                            searching_end = Integer.parseInt(arr[arr.length-1]);
                            searching_start = Integer.parseInt(arr[arr.length-2]);


                            int last_idx = aHSP.q_contig.lastIndexOf('!');
                            int second_last_idx = aHSP.q_contig.substring(0, last_idx).lastIndexOf('!');
                            searching_contig = aHSP.q_contig.substring(0, second_last_idx);

                            absolute_qcontig_start += searching_start -1;
                            absolute_qcontig_end += searching_start -1;


                        }

                    }//for each maf file

                    // now we have correct absolute coordinates of query contig, we add it query to MAFLines
                    MAFSLine line = new MAFSLine(absolute_qcontig_name, absolute_qcontig_start, absolute_qcontig_end, absolute_qcontig_pos_strand,  absolute_qcontig_size, "");
                    mafEntry.mafLines.set(0, line);

                    //Now we update standness of all subject contigs that is depends on standness of query

                    for(int i=1; i<mafEntry.mafLines.size(); i++)
                        mafEntry.mafLines.get(i).pos_strand = mafEntry.mafLines.get(i).orientation==1?absolute_qcontig_pos_strand:!absolute_qcontig_pos_strand;

                    // now we update all the alignment blocks

                    for(int i=0; i<buffList.size(); i++){
                        buffList.get(i).append(partialBuffList.get(i).toString());
                        mafEntry.mafLines.get(i).alignment_block = partialBuffList.get(i).toString();

                    }

                    mafEntry.validateMAFEntry();
                    //now we write whole maf entry
                    writeNextMulti(mafWriter,mafEntry);
                }catch(Exception ex){}
                        
                    
            }// for each query contig in  the last maf file
            System.out.println("start writing to output file...");
            long qlength = buffList.get(0).length();
            for(int i=0; i<buffList.size(); i++){
                String name = null;
                if (i==0)
                    name = names.get(0);
                else 
                    name = names.get(names.size()-i);
                int slength = buffList.get(i).length();
                if(qlength != slength)
                    System.out.println(String.format("subject %s(%d) length was not equal to the query (%d)", name, slength, qlength )) ;
                fafWriter.write(String.format(">%s\n%s\n", name ,buffList.get(i).toString()));                            

            }
            System.out.println("Length of Alignment " + buffList.get(0).length());
            System.out.println("Done!");

        } catch (Exception ex) {
            System.out.println(ex.getMessage()  );
        }finally{
            try {
                fafWriter.close();
                mafWriter.close();
            } catch (IOException ex) {
                //Logger.getLogger(Start.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        
    }

   
    public BlastHSP getNext(){
        //updated on 5/5/22
        // all the cordinates from MAF will be converted to Blast -like
        //In Blast result, query is always + but subject can be -, in MAF it is opposite: Subject (First Line) is always + , but query(Second line) can be -
        // We change orientation of query and subject to be like Blast (like blast start and end are 1-based; unlike blast start is  always less than end even when - orientation )
        BlastHSP hsp = new BlastHSP();
        //System.out.println("Inside getNext()");
        try{

            // read until we get to a line like a score=176450 or end of file
            String line =  buf.readLine();
            while(line!= null){
                
                if(line.matches("a\\s+(score=)?.+")){
                    // extract score from a line
                    hsp.scoreLine = line;
                    String second_part = line.split("\\s+")[1].trim();
                    int equal_idx = second_part.indexOf('=');
                    if(equal_idx>=0)
                        hsp.score = Double.parseDouble(second_part.substring(equal_idx+1));
                    else
                        hsp.score = Double.parseDouble(second_part);
                    // read next line is subject and the line after that is query line
                    line = buf.readLine(); // subject line
                    hsp.subjectLine = line;
                    if (line==null)
                        throw new Exception("Unexpected end of MAF file: no subject line");
                    String[] parts = line.split("\\s+");
                    //String 
                    hsp.s_contig = parts[1];
                    hsp.s_len = Integer.parseInt(parts[3]);                    
                    int zero_start = Integer.parseInt(parts[2]);
                    hsp.s_alignment_block = parts[6];
                    hsp.s_contig_size = Integer.parseInt(parts[5]);
                     
                    
                    if(parts[4].equals("+")){
                        // positive orientation
                        
                        hsp.s_pos_strand = true;
                        hsp.s_start = zero_start+1; // convert zero based into 1 based
                        hsp.s_end = hsp.s_start + hsp.s_len -1;
                    }else{
                        // negative orientation 
                        //In Maf standard, 'start' is 0-based and relative to - strand
                      //In blast standard 'start' and 'end' of hsp is 1-based and relative to positive strand
                        hsp.s_pos_strand = false;
                        
                        hsp.s_end = hsp.s_contig_size - zero_start;
                        hsp.s_start = hsp.s_end - hsp.s_len + 1;
                    }

                    hsp.alignment_length = hsp.s_alignment_block.length();
                    // read next line as query
                    line = buf.readLine(); // query line
                    hsp.queryLine = line;
                    if (line==null)
                        throw new Exception("Unexpected end of MAF file: no query line");
                    parts = line.split("\\s+");
                    //String 
                    hsp.q_contig = parts[1];
                    hsp.q_len = Integer.parseInt(parts[3]);
                    zero_start = Integer.parseInt(parts[2]);
                    hsp.q_contig_size = Integer.parseInt(parts[5]);
                    hsp.q_alignment_block = parts[6];
                    
                    if(parts[4].equals("+")){
                        // positive orientation  
                        
                        hsp.q_pos_strand = true;
                        hsp.q_start = zero_start+1; // convert zero based into 1 based
                        hsp.q_end = hsp.q_start + hsp.q_len -1;
                    }else{
                        // negative orientation  
                        
                        hsp.q_pos_strand = false;
                        hsp.q_end = hsp.q_contig_size - zero_start;
                        hsp.q_start = hsp.q_end - hsp.q_len + 1;
                    }
 
                    return hsp;
                }//if
                
                
                
                line =  buf.readLine();
            }//while

            //firstLine = null;
        }catch(Exception ex){System.out.println(ex.getMessage());}
        
        
        return null;
    }
    private BlastHSP getNextPAF(){
        //updated on 16/6/22
        // all the cordinates from PAF will be converted to Blast -like
        //In Blast result, query is always + but subject can be -, in PAF it is opposite: Subject (First Line) is always + , but query(Second line) can be -
        // We change orientation of query and subject to be like Blast (like blast start and end are 1-based; unlike blast start is  less than end even when - orientation )
        BlastHSP hsp = new BlastHSP();
        //System.out.println("Inside getNext()");
        try{

            // read until we get to a line like a score=176450 or end of file
            String line =  buf.readLine();
            while(line!= null){
                
                if(line.matches("a\\sscore=.+")){
                    // extract score from a line
                    hsp.scoreLine = line;
                    int equal_idx = line.indexOf('=');
                    hsp.score = Double.parseDouble(line.substring(equal_idx+1));
                    // read next line is subject and the line after that is query line
                    line = buf.readLine(); // subject line
                    hsp.subjectLine = line;
                    if (line==null)
                        throw new Exception("Unexpected end of MAF file: no subject line");
                    String[] parts = line.split("\\s+");
                    //String 
                    hsp.s_contig = parts[1];
                    hsp.s_len = Integer.parseInt(parts[3]);                    
                    int zero_start = Integer.parseInt(parts[2]);
                    hsp.s_alignment_block = parts[6];
                    hsp.s_contig_size = Integer.parseInt(parts[5]);
                     
                    
                    if(parts[4].equals("+")){
                        // positive orientation
                        
                        hsp.s_pos_strand = true;
                        hsp.s_start = zero_start+1; // convert zero based into 1 based
                        hsp.s_end = hsp.s_start + hsp.s_len -1;
                    }else{
                        // negative orientation 
                        //In Maf standard, 'start' is 0-based and relative to - strand
                      //In blast standard 'start' and 'end' of hsp is 1-based and relative to positive strand
                        hsp.s_pos_strand = false;
                        
                        hsp.s_end = hsp.s_contig_size - zero_start;
                        hsp.s_start = hsp.s_end - hsp.s_len + 1;
                    }

                    hsp.alignment_length = hsp.s_alignment_block.length();
                    // read next line as query
                    line = buf.readLine(); // query line
                    hsp.queryLine = line;
                    if (line==null)
                        throw new Exception("Unexpected end of MAF file: no query line");
                    parts = line.split("\\s+");
                    //String 
                    hsp.q_contig = parts[1];
                    hsp.q_len = Integer.parseInt(parts[3]);
                    zero_start = Integer.parseInt(parts[2]);
                    hsp.q_contig_size = Integer.parseInt(parts[5]);
                    hsp.q_alignment_block = parts[6];
                    
                    if(parts[4].equals("+")){
                        // positive orientation  
                        
                        hsp.q_pos_strand = true;
                        hsp.q_start = zero_start+1; // convert zero based into 1 based
                        hsp.q_end = hsp.q_start + hsp.q_len -1;
                    }else{
                        // negative orientation  
                        
                        hsp.q_pos_strand = false;
                        hsp.q_end = hsp.q_contig_size - zero_start;
                        hsp.q_start = hsp.q_end - hsp.q_len + 1;
                    }
 
                    return hsp;
                }//if
                
                
                
                line =  buf.readLine();
            }//while

            //firstLine = null;
        }catch(Exception ex){System.out.println(ex.getMessage());}
        
        
        return null;
    }
    
    private BlastHSP getNextMulti(){
        // all the cordinates from MAF will be converted to Blast -like
        //In Blast result, query is always + but subject can be -, in MAF it is opposite: Subject (First Line) is always + , but query(Second line) can be -
        // unlike blast, start_pos always less than end_pos
        BlastHSP hsp = new BlastHSP();
        hsp.mafLines = new ArrayList();
        //System.out.println("Inside getNext()");
        try{

            // read until we get to a line like: a score=176450 or end of file
            String line =  buf.readLine();
            while(line!= null){
    
                if(line.toLowerCase().matches("a\\s*(score=)?.*") ){
                    // extract score from a line
                    hsp.scoreLine = line;
                    //determine the score
                    if (line.split("\\s+").length > 1){
                        String second_part = line.split("\\s+")[1].trim();
                        int equal_idx = second_part.indexOf('=');
                        // find end of equal index
                        if (equal_idx>0)
                              hsp.score = Double.parseDouble(second_part.substring(equal_idx+1));
                        else
                            hsp.score = Double.parseDouble(second_part);
                    }else{
                        hsp.score = 0;
                    }
                            

                    // read next S lines until  
                    line = buf.readLine(); 
                    while(line!= null && line.toLowerCase().startsWith("s")){
                       String[] parts = line.split("\\s+");
                       MAFSLine  mafLine = new MAFSLine(parts[1], Integer.parseInt(parts[2]), Integer.parseInt(parts[3]), parts[4].equals("+")?true:false,Integer.parseInt(parts[5]), parts[6], true);
                       mafLine.sLine = line;
                       hsp.mafLines.add(mafLine);
                       line = buf.readLine(); 
                    }

                    // when we get here , we have reached to a line that does not start with s
                    
 
                    return hsp;
                }//if
                
                line =  buf.readLine();
            }//while

            //firstLine = null;
        }catch(Exception ex){System.out.println(ex.getMessage());}
        
        
        return null;
    }
    public void writeNext(Writer mafFile, BlastHSP hsp ) throws Exception{
        //opposite of getNext
        //writes one hsp into new MAF file
        // This function assumes each hsp is loaded with Blast like HSPs, so coordinates are 1-based and inclusive
        
        
            //First write a Line
            
            mafFile.write(String.format("a  score=%f\n", hsp.score));
           int s_zero_start = 0;
           if(hsp.s_pos_strand){
               s_zero_start = hsp.s_start - 1;
           }else{
               s_zero_start = hsp.s_contig_size - hsp.s_end;
           }
           // now write s line
           mafFile.write(String.format("s\t%s\t%d\t%d\t%s\t%d\t%s\n", hsp.s_contig, s_zero_start, hsp.s_len, hsp.s_pos_strand?"+":"-", hsp.s_contig_size, hsp.s_alignment_block));

           int q_zero_start = 0;
           if(hsp.q_pos_strand){
               q_zero_start = hsp.q_start - 1;
           }else{
               q_zero_start = hsp.q_contig_size - hsp.q_end;
           }
           // now write s line
           mafFile.write(String.format("s\t%s\t%d\t%d\t%s\t%d\t%s\n", hsp.q_contig, q_zero_start, hsp.q_len, hsp.q_pos_strand?"+":"-", hsp.q_contig_size, hsp.q_alignment_block));
           
           //Now write empty line
            mafFile.write("\n");
        
    }
    
    public void writeNextMulti(Writer mafFile, BlastHSP hsp ) throws Exception{
        //opposite of getNextMulti
        //writes one hsp into new MAF file
        // This function assumes each hsp is loaded with Blast like HSPs, so coordinates are 1-based and inclusive
        
        
            //First write a Line
            
            mafFile.write(String.format("a  score=%f\n", hsp.score));
            for(MAFSLine line : hsp.mafLines){
                
               int zero_start = 0;
               int frag_len = 0;
               if(line.pos_strand){
                   zero_start = line.startPOS - 1;
               }else{
                   zero_start = line.contigSize - line.endPOS;
               }
                frag_len = line.endPOS - line.startPOS + 1;
               // now write s line
               mafFile.write(String.format("s  %s %d %d %s %d %s\n", line.contig, zero_start, frag_len, line.pos_strand?"+":"-", line.contigSize, line.alignment_block));
            }
           //Now write empty line
            mafFile.write("\n");
        
    }
    
    public void readWriteMulti(String newMAFFile){
        
        //this is for test only
        FileWriter writer = null;
        try{
            writer = new FileWriter(newMAFFile);
            BlastHSP hsp = getNextMulti();
            int cc = 0;
            while(hsp != null){
               writeNextMulti(writer, hsp);
               cc++;
               hsp = getNextMulti();
               if (hsp.mafLines.size()<2)
                    System.out.println("MAF Entry should have at least 2 Lines for subject and query , but current size is " + hsp.mafLines.size());
            } 
            System.out.println(cc + " records written");
        }catch(Exception ex){System.out.println(ex.getMessage());}
        finally{
            try {
                writer.close();
            } catch (IOException ex) {Logger.getLogger(MAFProcessor.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }
    public void readWriteMulti(String newMAFFile, int min_size){
        // This is adhoc to check cactus output
        //only writes MAF entries with min_size number of  MAFLines
        Set<String> bacteria_names = new HashSet(Arrays.asList("Anc0", "fung1", "fung2", "fung3", "fung4", "fung5", "fung6", "fung7","fung8","fung9"));
// for each bacteria we count amount of overlaps and report it in this map
        Map<String, Integer> overlaps = new HashMap();
        Map<String, List<MAFSLine>> lists = new HashMap();
        lists.put("Anc0", new ArrayList());
        lists.put("fung1", new ArrayList());
        lists.put("fung2", new ArrayList());
        lists.put("fung3", new ArrayList());
        lists.put("fung4", new ArrayList());
        lists.put("fung5", new ArrayList());
        lists.put("fung6", new ArrayList());
        lists.put("fung7", new ArrayList());
        lists.put("fung8", new ArrayList());
        lists.put("fung9", new ArrayList());
        
        
        FileWriter writer = null;
        try{
            writer = new FileWriter(newMAFFile);
            BlastHSP hsp = getNextMulti();
            int cc = 0;
            long len = 0; // total length of alignment blocks that are written to output
            
            while(hsp != null){
                hsp.calculateOverlaps(overlaps, lists, 5);
                if (hsp.getFirstContigPartSet().containsAll(bacteria_names)){
                    
                   cc++; 
                   len += hsp.mafLines.get(0).alignment_block.length();
                   writeNextMulti(writer, hsp);
                }
               
               
               hsp = getNextMulti();
            } //while
            System.out.println(cc + " records written");
            System.out.println(len + " Total Written Alignment Length");
            System.out.println("Calculated overlaps: ");
            for(String key:overlaps.keySet())
                System.out.println(key + " : " + overlaps.get(key));
        }catch(Exception ex){System.out.println(ex.getMessage());}
        finally{
            try {
                writer.close();
            } catch (IOException ex) {Logger.getLogger(MAFProcessor.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }    

    public void analyzeCactus(String fastaFile, String newMAFFile){
        // This is adhoc to check cactus output
        //only writes MAF entries with that contains all species
        //HashSet<String> bacteria_names = new HashSet(Arrays.asList("Drosophila_simulans", "Drosophila_melanogaster", "Drosophila_sechellia", "Drosophila_yakuba", "Drosophila_erecta", "Drosophila_ananassae", "Drosophila_pseudoobscura","Drosophila_persimilis","Drosophila_willistoni","Drosophila_virilis","Drosophila_mojavensis","Drosophila_grimshawi"));
        HashSet<String> bacteria_names = new HashSet(Arrays.asList("mus_musculus",
"Allactaga_bullata",
"Aplodontia_rufa",
"Castor_canadensis",
"Cricetomys_gambianus",
"Cricetulus_griseus",
"Dipodomys_ordii",
"Dipodomys_stephensi",
"Ellobius_lutescens",
"Ellobius_talpinus",
"Glis_Glis",
"Graphiurus_murinus",
"Ictidomys_tridecemlineatus",
"Jaculus_jaculus",
"Marmota_marmota",
"Meriones_unguiculatus",
"Mesocricetus_auratus",
"Microtus_ochrogaster",
"Mus_caroli",
"Mus_pahari",
"Mus_spretus",
"Muscardinus_avellanarius",
"Nannospalax_galili",
"Ondatra_zibethicus",
"Onychomys_torridus",
"Perognathus_longimembris",
"Peromyscus_maniculatus",
"Psammomys_obesus",
"Rattus_norvegicus",
"Sigmodon_hispidus",
"Spermophilus_dauricus",
"Xerus_inauris",
"Zapus_hudsonius",
"Acomys_cahirinus"));

        

        List<Integer> alignment_lengths = new ArrayList();
        Map<String, StringBuffer> seqs = new HashMap();
        //seqs.put("Anc0", new StringBuffer());
        /*seqs.put("Drosophila_simulans", new StringBuffer());
        seqs.put("Drosophila_melanogaster", new StringBuffer());
        seqs.put("Drosophila_sechellia", new StringBuffer());
        seqs.put("Drosophila_yakuba", new StringBuffer());
        seqs.put("Drosophila_erecta", new StringBuffer());
        seqs.put("Drosophila_ananassae", new StringBuffer());
        seqs.put("Drosophila_pseudoobscura", new StringBuffer());
        seqs.put("Drosophila_persimilis", new StringBuffer());
        seqs.put("Drosophila_willistoni", new StringBuffer());
        seqs.put("Drosophila_virilis", new StringBuffer());
        seqs.put("Drosophila_mojavensis", new StringBuffer());
        seqs.put("Drosophila_grimshawi", new StringBuffer());*/
seqs.put("mus_musculus", new StringBuffer());
seqs.put("Allactaga_bullata", new StringBuffer());
seqs.put("Aplodontia_rufa", new StringBuffer());
seqs.put("Castor_canadensis", new StringBuffer());
seqs.put("Cricetomys_gambianus", new StringBuffer());
seqs.put("Cricetulus_griseus", new StringBuffer());
seqs.put("Dipodomys_ordii", new StringBuffer());
seqs.put("Dipodomys_stephensi", new StringBuffer());
seqs.put("Ellobius_lutescens", new StringBuffer());
seqs.put("Ellobius_talpinus", new StringBuffer());
seqs.put("Glis_Glis", new StringBuffer());
seqs.put("Graphiurus_murinus", new StringBuffer());
seqs.put("Ictidomys_tridecemlineatus", new StringBuffer());
seqs.put("Jaculus_jaculus", new StringBuffer());
seqs.put("Marmota_marmota", new StringBuffer());
seqs.put("Meriones_unguiculatus", new StringBuffer());
seqs.put("Mesocricetus_auratus", new StringBuffer());
seqs.put("Microtus_ochrogaster", new StringBuffer());
seqs.put("Mus_caroli", new StringBuffer());
seqs.put("Mus_pahari", new StringBuffer());
seqs.put("Mus_spretus", new StringBuffer());
seqs.put("Muscardinus_avellanarius", new StringBuffer());
seqs.put("Nannospalax_galili", new StringBuffer());
seqs.put("Ondatra_zibethicus", new StringBuffer());
seqs.put("Onychomys_torridus", new StringBuffer());
seqs.put("Perognathus_longimembris", new StringBuffer());
seqs.put("Peromyscus_maniculatus", new StringBuffer());
seqs.put("Psammomys_obesus", new StringBuffer());
seqs.put("Rattus_norvegicus", new StringBuffer());
seqs.put("Sigmodon_hispidus", new StringBuffer());
seqs.put("Spermophilus_dauricus", new StringBuffer());
seqs.put("Xerus_inauris", new StringBuffer());
seqs.put("Zapus_hudsonius", new StringBuffer());
seqs.put("Acomys_cahirinus", new StringBuffer());
        
        
        Map<Integer, Long> lens = new HashMap(); // a map of size of enries (here 2-13) to number of bp 
        
        FileWriter writer = null;
        FileWriter mafwriter = null;
        try{
            mafwriter = new FileWriter(newMAFFile);
            writer = new FileWriter(fastaFile);
            BlastHSP hsp = getNextMulti();
            int cc = 0;
            long len = 0; // total length of alignment blocks that are written to output
            
            while(hsp != null){
                   
                Set ss = hsp.getFirstContigPartSet();
                
                if (ss.containsAll(bacteria_names)&& !hsp.mafLines.get(0).alignment_block.matches("^N+$")){
                   cc++; 
                   alignment_lengths.add(hsp.mafLines.get(0).alignment_block.length());
                   len += hsp.mafLines.get(0).alignment_block.length();
                   // because some hsp entries could have more than one line for one bacteria, we just make sure we add it once
                   HashSet names_cloned = (HashSet)bacteria_names.clone();
                   for(MAFSLine line : hsp.mafLines){
                       String bname = line.contig.split("\\.")[0];
                       if(names_cloned.contains(bname)){
                          seqs.get(bname).append(line.alignment_block);
                          names_cloned.remove(bname);
                       }
                   }
                   writeNextMulti(mafwriter, hsp);
                   if (!hsp.equalLength())
                        System.out.println("not equal length");
                  
                }
               
                int size = hsp.mafLines.size();
                Long lenbp = lens.get(size);
                long lbp = hsp.mafLines.get(0).alignment_block.length();
                if (lenbp != null)
                    lens.put(size, lenbp+lbp);
                else
                    lens.put(size, lbp);
               
               hsp = getNextMulti();
            } //while
            for(String bacteria : seqs.keySet()){
                writer.write(String.format(">%s\n", bacteria));
                //System.out.println("length "+ seqs.get(bacteria).length());
                writer.write(seqs.get(bacteria)+"\n");
            }
            for(int ll : lens.keySet()){
               System.out.println(ll + "=>" + lens.get(ll));                
            }

            
            System.out.println(cc + " records written" + " calculated " + alignment_lengths.size());
            System.out.println(len + " Total Written Alignment Length(bp)");
            // now calculate N50 and 
        Collections.sort(alignment_lengths); // sort ascending order    
        //System.out.println(contig_lengths.get(0));
        //System.out.println(contig_lengths.get(contig_lengths.size()-1));
        int n50_len = 0;
        long n50 = (int) (0.5*len);
        long cumulative_len = 0;
        try{
            for(int i=alignment_lengths.size()-1; i>=0; i--){
                cumulative_len += alignment_lengths.get(i);
                if (cumulative_len>=n50){
                    n50_len = alignment_lengths.get(i);
                    break;
                }
            }           
        }catch(Exception ex){}
          System.out.println("N50:" + n50_len + " Min:" + alignment_lengths.get(0) + " Max:" + alignment_lengths.get(alignment_lengths.size()-1));   
        }catch(Exception ex){System.out.println(ex.getMessage());}
        finally{
            try {
                writer.close();
                mafwriter.close();
            } catch (IOException ex) {Logger.getLogger(MAFProcessor.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }    
    public void analyzeSibeliaz(String fastaFile, String newMAFFile){
        // This is adhoc to check cactus output
        //only writes MAF entries with that contains all species
        //HashSet<String> bacteria_names = new HashSet(Arrays.asList("Drosophila_simulans", "Drosophila_melanogaster", "Drosophila_sechellia", "Drosophila_yakuba", "Drosophila_erecta", "Drosophila_ananassae", "Drosophila_pseudoobscura","Drosophila_persimilis","Drosophila_willistoni","Drosophila_virilis","Drosophila_mojavensis","Drosophila_grimshawi"));
        HashSet<String> bacteria_names = new HashSet(Arrays.asList("sim", "mel", "sec", "yak", "ere", "ana", "pse","per","wil","vir","moj","gri"));
        

        List<Integer> alignment_lengths = new ArrayList();
        Map<String, StringBuffer> seqs = new HashMap();
        //seqs.put("Anc0", new StringBuffer());
        seqs.put("sim", new StringBuffer());
        seqs.put("mel", new StringBuffer());
        seqs.put("sec", new StringBuffer());
        seqs.put("yak", new StringBuffer());
        seqs.put("ere", new StringBuffer());
        seqs.put("ana", new StringBuffer());
        seqs.put("pse", new StringBuffer());
        seqs.put("per", new StringBuffer());
        seqs.put("wil", new StringBuffer());
        seqs.put("vir", new StringBuffer());
        seqs.put("moj", new StringBuffer());
        seqs.put("gri", new StringBuffer());
        
        Map<Integer, Long> lens = new HashMap(); // a map of size of enries (here 2-13) to number of bp 
        
        FileWriter writer = null;
        FileWriter mafwriter = null;
        try{
            mafwriter = new FileWriter(newMAFFile);
            writer = new FileWriter(fastaFile);
            BlastHSP hsp = getNextMulti();
            int cc = 0;
            long len = 0; // total length of alignment blocks that are written to output
            
            while(hsp != null){
                if (hsp.mafLines.get(0).contig.startsWith("Anc0"))
                    hsp.mafLines.remove(0); // remove Anc0 entry                    
                //else
                //    System.out.println("Anc0 not found");
                   
                Set ss = hsp.getFirstContigPartSet2();
                
                if (ss.containsAll(bacteria_names)){
                                
                   cc++; 
                   alignment_lengths.add(hsp.mafLines.get(0).alignment_block.length());
                   len += hsp.mafLines.get(0).alignment_block.length();
                   // because some hsp entries could have more than one line for one bacteria, we just make sure we add it once
                   HashSet names_cloned = (HashSet)bacteria_names.clone();
                   for(MAFSLine line : hsp.mafLines){
                       String bname = line.contig.split("_")[0];
                       if(names_cloned.contains(bname)){
                          seqs.get(bname).append(line.alignment_block);
                          names_cloned.remove(bname);
                       }
                   }
                   writeNextMulti(mafwriter, hsp);
                   if (!hsp.equalLength())
                        System.out.println("not equal length");
                  
                }
               
                int size = hsp.mafLines.size();
                Long lenbp = lens.get(size);
                long lbp = hsp.mafLines.get(0).alignment_block.length();
                if (lenbp != null)
                    lens.put(size, lenbp+lbp);
                else
                    lens.put(size, lbp);
               
               hsp = getNextMulti();
            } //while
            for(String bacteria : seqs.keySet()){
                writer.write(String.format(">%s\n", bacteria));
                //System.out.println("length "+ seqs.get(bacteria).length());
                writer.write(seqs.get(bacteria)+"\n");
            }
            for(int ll : lens.keySet()){
               System.out.println(ll + "=>" + lens.get(ll));                
            }

            
            System.out.println(cc + " records written" + " calculated " + alignment_lengths.size());
            System.out.println(len + " Total Written Alignment Length(bp)");
            // now calculate N50 and 
        Collections.sort(alignment_lengths); // sort ascending order    
        //System.out.println(contig_lengths.get(0));
        //System.out.println(contig_lengths.get(contig_lengths.size()-1));
        int n50_len = 0;
        long n50 = (int) (0.5*len);
        long cumulative_len = 0;
        try{
            for(int i=alignment_lengths.size()-1; i>=0; i--){
                cumulative_len += alignment_lengths.get(i);
                if (cumulative_len>=n50){
                    n50_len = alignment_lengths.get(i);
                    break;
                }
            }           
        }catch(Exception ex){}
          System.out.println("N50:" + n50_len + " Min:" + alignment_lengths.get(0) + " Max:" + alignment_lengths.get(alignment_lengths.size()-1));   
        }catch(Exception ex){System.out.println(ex.getMessage());}
        finally{
            try {
                writer.close();
                mafwriter.close();
            } catch (IOException ex) {Logger.getLogger(MAFProcessor.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }    

    public void checkCordinates(String fastaFileNames){
        // checks if coordinates of contigs match the sequence
        try{
            List<IndexedFastaSequenceFile> fastaList = new ArrayList();
            BufferedReader nameReader = new BufferedReader(new FileReader(fastaFileNames));
            String nameLine = nameReader.readLine();
            while(nameLine != null && !nameLine.isEmpty()){
                String[] names = nameLine.split("\\s+");
                if (names.length<2){
                    System.out.println("Error: Fasta file name and its index must be provided");
                    return;                    
                }
                String filename = names[0];
                String fileindex = names[1];
                FastaSequenceIndex fsi = new FastaSequenceIndex(new File(fileindex));
                IndexedFastaSequenceFile ifsi =new IndexedFastaSequenceFile(new File(filename) , fsi);
                fastaList.add(ifsi);
                
                nameLine = nameReader.readLine();
            }
            nameReader.close();
            if(fastaList.size()<2){
                System.out.println("Error: At least 2 Fasta files must be provided");
                return;
            }

            BlastHSP hsp = getNextMulti();
            if(hsp.mafLines.size()>fastaList.size()){
                System.out.println("Error: Number of MAF entries are more than number of fasta files");
                return;                
            }
            while(hsp != null){ 
                for (int i =0; i<hsp.mafLines.size(); i++){
                    MAFSLine mafline = hsp.mafLines.get(i);
                    ReferenceSequence seq = fastaList.get(i).getSubsequenceAt(mafline.contig, mafline.startPOS, mafline.endPOS);
                    String extracted_seq = seq.getBaseString().toUpperCase();
                    if (!mafline.pos_strand)
                        extracted_seq = Helper.getRC(extracted_seq);
                    
                    String lineseq = mafline.alignment_block.replaceAll("-", "").toUpperCase();
                    if (!extracted_seq.equals(lineseq)){
                        System.out.println("Following S Line not match:");
                        System.out.println(mafline.sLine);
                        System.out.println(extracted_seq);                        
                    }
                    
                }
                
/*               alignment_lengths.add(hsp.q_alignment_block.length());
               total_len += hsp.q_alignment_block.length();
                for(int i=0; i< hsp.q_alignment_block.length(); i++){
                    char q_ch = hsp.q_alignment_block.charAt(i);
                    char s_ch = hsp.s_alignment_block.charAt(i);
                    if (q_ch == 'N' || q_ch=='n' || s_ch == 'N' || s_ch=='n')
                        n_num++;                        
                    else if (q_ch=='-' ){
                        gap_num++;
                        if (i>0 && hsp.q_alignment_block.charAt(i-1)!='-')
                            gap_open++; // this means we have entered a block of gaps, so we increase gap_event
                    }else if (s_ch=='-' ){
                        gap_num++;
                        if (i>0 && hsp.s_alignment_block.charAt(i-1)!='-')
                            gap_open++;
                    }else if (Helper.nt_to_i(q_ch)!=Helper.nt_to_i(s_ch)){
                        mismatch_num++;
                    }
//
                    if (Character.isLowerCase(q_ch) || Character.isLowerCase(s_ch))
                        low_num++;
                }//for
*/
                hsp = getNextMulti();
               
            }// while we have more hsp
            
        } catch (Exception ex) {System.out.println("Error at checkCordinates:" + ex.getMessage() );}
        
    }
    public void info(){
        // prepares like fastainfo, number of entries, total_length, N50, N70, total gap occurence, total gap (bp), gap-compressed identity, gap_uncopressed identity
        List<Integer> alignment_lengths = new ArrayList();
        int total_len = 0;
        int n_num = 0;
        int low_num = 0;
        int gap_num = 0;
        int gap_open = 0; // continuous gaps are couned only 1
        int mismatch_num = 0; //
       
        try{
            BlastHSP hsp = getNext();
            // in the last iteration of below loop hsp is null, but lastContigWritten will set to true and the loop will exit, because hsp is null and lastContigWritten is set to true
            while(hsp != null){ 
               
               alignment_lengths.add(hsp.q_alignment_block.length());
               total_len += hsp.q_alignment_block.length();
                for(int i=0; i< hsp.q_alignment_block.length(); i++){
                    char q_ch = hsp.q_alignment_block.charAt(i);
                    char s_ch = hsp.s_alignment_block.charAt(i);
                    if (q_ch == 'N' || q_ch=='n' || s_ch == 'N' || s_ch=='n')
                        n_num++;                        
                    else if (q_ch=='-' ){
                        gap_num++;
                        if (i>0 && hsp.q_alignment_block.charAt(i-1)!='-')
                            gap_open++; // this means we have entered a block of gaps, so we increase gap_event
                    }else if (s_ch=='-' ){
                        gap_num++;
                        if (i>0 && hsp.s_alignment_block.charAt(i-1)!='-')
                            gap_open++;
                    }else if (Helper.nt_to_i(q_ch)!=Helper.nt_to_i(s_ch)){
                        mismatch_num++;
                    }
//
                    if (Character.isLowerCase(q_ch) || Character.isLowerCase(s_ch))
                        low_num++;
                }//for
               hsp = getNext();
               
            }// while we have more hsp
            
        Collections.sort(alignment_lengths); // sort ascending order    
        //System.out.println(contig_lengths.get(0));
        //System.out.println(contig_lengths.get(contig_lengths.size()-1));
        int n50_len = 0;
        long n50 = (int) (0.5*total_len);
        long cumulative_len = 0;
        try{
            for(int i=alignment_lengths.size()-1; i>=0; i--){
                cumulative_len += alignment_lengths.get(i);
                if (cumulative_len>=n50){
                    n50_len = alignment_lengths.get(i);
                    break;
                }
            }           
        }catch(Exception ex){}
        double gapcompressed_divergance = (double)(mismatch_num+gap_open)/(total_len - n_num);
        double gapuncompressed_divergance = (double)(mismatch_num+gap_num)/(total_len - n_num);
        
        if (alignment_lengths.size()>0)            
           System.out.println(String.format("Total HSP:%d, Total Alignment Length(bp) :%d, Max alignment Length(bp):%d, Min alignment Length(bp):%d, N50:%d, gap Open:%d, gaps:%d, mismatch:%d, N:%d, lower case:%d, divergance %f, gap-compressed divergance %f",  alignment_lengths.size(), total_len, alignment_lengths.get(alignment_lengths.size()-1), alignment_lengths.get(0),  n50_len, gap_open, gap_num, mismatch_num, n_num, low_num, gapuncompressed_divergance, gapcompressed_divergance));    
        } catch (Exception ex) {System.out.println("Error at info:" + ex.getMessage() );}
 
    }
    public void info_multi(){
        // prepares like fastainfo, number of entries, total_length, N50, N70, total gap occurence, total gap (bp), gap-compressed identity, gap_uncopressed identity
        List<Integer> alignment_lengths = new ArrayList();
        int total_len = 0;
        int n_num = 0;
        int low_num = 0;
        int gap_num = 0;
        int gap_open = 0; // continuous gaps are couned only 1
       
        try{
            BlastHSP hsp = getNextMulti();
            // in the last iteration of below loop hsp is null, but lastContigWritten will set to true and the loop will exit, because hsp is null and lastContigWritten is set to true
            while(hsp != null){ 
               MAFSLine line = hsp.mafLines.get(0);
               alignment_lengths.add(line.alignment_block.length());
               total_len += line.alignment_block.length();
                for(int i=0; i< line.alignment_block.length(); i++){
                    char ch = line.alignment_block.charAt(i);
                    //char s_ch = hsp.s_alignment_block.charAt(i);
                    if (ch == 'N' || ch=='n' )
                        n_num++;                        
                    else if (ch=='-' ){
                        gap_num++;
                        if (i>0 && line.alignment_block.charAt(i-1)!='-')
                            gap_open++; // this means we have entered a block of gaps, so we increase gap_event
                    }
//
                    if (Character.isLowerCase(ch) )
                        low_num++;
                }//for
               hsp = getNextMulti();
               
            }// while we have more hsp
            
        Collections.sort(alignment_lengths); // sort ascending order    
        //System.out.println(contig_lengths.get(0));
        //System.out.println(contig_lengths.get(contig_lengths.size()-1));
        int n50_len = 0;
        long n50 = (int) (0.5*total_len);
        long cumulative_len = 0;
        try{
            for(int i=alignment_lengths.size()-1; i>=0; i--){
                cumulative_len += alignment_lengths.get(i);
                if (cumulative_len>=n50){
                    n50_len = alignment_lengths.get(i);
                    break;
                }
            }           
        }catch(Exception ex){}
        if (alignment_lengths.size()>0)            
           System.out.println(String.format("Total HSP:%d, Total Alignment Length(bp) :%d, Max alignment Length(bp):%d, Min alignment Length(bp):%d, N50:%d, gap Open:%d, gaps:%d, N:%d, lower case:%d",  alignment_lengths.size(), total_len, alignment_lengths.get(alignment_lengths.size()-1), alignment_lengths.get(0),  n50_len, gap_open, gap_num,  n_num, low_num));    
        } catch (Exception ex) {System.out.println("Error at info:" + ex.getMessage() );}
 
    }

    public void analyzeMugsy(String fastaFile){
        // This is adhoc to check mugsy output
        //only writes MAF entries with that contains all bacteria
        Set<String> bacteria_names = new HashSet(Arrays.asList("new_race_ARCrossB10_upper","race_1_isolate134_upper","race_1_isolate239_upper","race_1_isolate5213_upper","race_1_isolate11137_upper","race_1_isolateBFP_upper","race_1_isolateM4_upper","race_2_86_124_upper","race_5_isolateDW5_upper"));

        Map<String, StringBuffer> seqs = new HashMap();
        seqs.put("race_5_isolateDW5_upper", new StringBuffer());
        seqs.put("race_2_86_124_upper", new StringBuffer());
        seqs.put("race_1_isolateM4_upper", new StringBuffer());
        seqs.put("race_1_isolateBFP_upper", new StringBuffer());
        seqs.put("race_1_isolate11137_upper", new StringBuffer());
        seqs.put("race_1_isolate5213_upper", new StringBuffer());
        seqs.put("race_1_isolate239_upper", new StringBuffer());
        seqs.put("race_1_isolate134_upper", new StringBuffer());
        seqs.put("new_race_ARCrossB10_upper", new StringBuffer());
        
        Map<Integer, Long> lens = new HashMap(); // a map of size of enries (here 2-9) to number of bp 
        
        FileWriter writer = null;
        try{
            writer = new FileWriter(fastaFile);
            BlastHSP hsp = getNextMulti();
            int cc = 0;
            long len = 0; // total length of alignment blocks that are written to output
            
            while(hsp != null){
                Set ss = hsp.getFirstContigPartSet();
                
                if (ss.containsAll(bacteria_names)){
                    
                   cc++; 
                   len += hsp.mafLines.get(0).alignment_block.length();
                   for(MAFSLine line : hsp.mafLines)
                       seqs.get(line.contig.split("\\.")[0]).append(line.alignment_block);
                  
                }
               
                int size = hsp.mafLines.size();
                Long lenbp = lens.get(size);
                long lbp = hsp.mafLines.get(0).alignment_block.length();
                if (lenbp != null)
                    lens.put(size, lenbp+lbp);
                else
                    lens.put(size, lbp);
                    
               hsp = getNextMulti();
            } //while
            for(String bacteria : seqs.keySet()){
                writer.write(String.format(">%s\n", bacteria));
                writer.write(seqs.get(bacteria)+"\n");
            }
            for(int ll : lens.keySet()){
               System.out.println(ll + "=>" + lens.get(ll));                
            }

            
            System.out.println(cc + " records written");
            System.out.println(len + " Total Written Alignment Length");
        }catch(Exception ex){System.out.println(ex.getMessage());}
        finally{
            try {
                writer.close();
            } catch (IOException ex) {Logger.getLogger(MAFProcessor.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }    
}
