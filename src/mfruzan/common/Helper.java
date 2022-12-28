/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package mfruzan.common;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;
import mfruzan.algorithm.SNP;
import static mfruzan.common.Helper.T;
import htsjdk.variant.variantcontext.VariantContext;
import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import mfruzan.algorithm.MapUtil;
import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.ChromosomeSequence;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.GeneSequence;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.Strand;
import org.biojava3.core.sequence.TranscriptSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;

/**
 *
 * @author a1195806
 */
public class Helper {
        public static final String TruSeqUniversalAdapter =   "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"; // 58 bases
    public static final String RCTruSeqUniversalAdapter = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT";  // 58 bases
    public static final String  TruSeqIndexedAdapter =     "GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG"; // 63 bases
    public static final String  RCTruSeqIndexedAdapter = "CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"; // 64 bases, notice of an extra T at the end
    //These are for truseq HT protocol
    
    public static final String TruSeqUniversalAdapterHT =   "AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT"; // 70 bases
    public static final String RCTruSeqUniversalAdapterHT = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT";  // 70 bases
    public static final String  TruSeqIndexedAdapterHT =     "GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG"; // 65 bases
    public static final String  RCTruSeqIndexedAdapterHT = "CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC"; // 65 bases
    
    
    //These are Nextra protocol   CTGTCTCTTATACACATCT, CTGTCTCTTATACACATCT + AGATGTGTATAAGAGACAG
    public static final String NextraAdapter = "ATGTGTATAAGAGACA";
    public static final String NextraAdapter_V2 = "AGATGTGTATAAGAGACAG";
    public static final String RCNextraAdapter = "CTGTCTCTTATACACATCT";
    public static final String NextraI5Adapter =   "AATGATACGGCGACCACCGAGATCTACACNNNNNNNNTCGTCGGCAGCGTC"; // 51 bases
    public static final String RCNextraI5Adapter = "GACGCTGCCGACGANNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT";  // 51 bases
    public static final String  NextraI7Adapter =     "CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTCTCGTGGGCTCGG"; // 47 bases
    public static final String  RCNextraI7Adapter = "CCGAGCCCACGAGACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG"; // 47 bases
    
    //public static final String  TruSeqIndexedAdapter1 =   "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"; // 33 bases
    // public static final String  TruSeqIndexedAdapter2 =   "ATCTCGTATGCCGTCTTCTGCTTG";    // 24 bases
     
    // public static final String  RCTruSeqIndexedAdapter1 =   "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT";  // 34 bases extra T
    // public static final String  RCTruSeqIndexedAdapter2 =   "CAAGCAGAAGACGGCATACGAGAT";   // 24 bases
     
    public static final String LongReadAdapter = "CCGGTTCTTCCCTGCCGAACCCTATCTTCGTCGGCAGCGTCAGATGTGTATAAGAGACAGTACGCTTGCAT";
    
    // type of mutation
   public static final byte MT_DELETION = 1; // base  deletion
   public static final byte MT_INSERTION = 2; // base  insertion
   public static final byte MT_Substitution = 3; // base  change

    // type of Variation
   public static final byte VR_DELETION = 1; // base  deletion
   public static final byte VR_INSERTION = 2; // base  insertion
   public static final byte VR_SNP = 3; // base  change
   public static final byte VR_Default = 4; // base  change
   
   // type of SNP or polymorphic position
   public static final byte SNP_MultiBase = 1; // different bases, or Deletion
   public static final byte SNP_MultiInsert = 2; // different insertions, or presence or absence of an Insetion
   public static final byte SNP_MultiDelete = 3; // different insertions, or presence or absence of an Insetion
   
   //relative positions of 2 sequences by each other
   public static final byte POS_Q_Extended_S = 1; // Query sequence is extended by the Subject one from right hand side
   public static final byte POS_S_Extended_Q = 2; // one sequence is extended by the second one from right hand side
   public static final byte POS_Q_Inside_S = 3; // one seq is inside the other one
   public static final byte POS_S_Inside_Q = 4; // one seq is inside the other one
   public static final byte POS_Distinct = 5; // 2 sequences have no overlap

   public static final byte MT_Default = 4; // default mutation type
   public static final byte MT_Hidden = 5; // This is subtype of Transition mutation, when ratio test will reveal mutation
  
   public static final byte ZG_UNKNOWN = 0; // Zygosity unknown
   public static final byte ZG_HOMOZYGOUS = 1; // Homozygous
   public static final byte ZG_HETROZYGOUS = 2; // Hetrozygous   
   public static final byte ZG_INVALID = -1; // something is not right!
   //zygosity ratio in tetraploid
   public static final byte ZGR_4_4 = 1; //  out of four alleles all of them are of type mutation 
   public static final byte ZGR_4_2 = 2;  // out of four allels half of them are of type mutation
   public static final byte ZGR_4_1 = 3;  // out of four allels one-forth of them are of type mutation
   
   // number neucleotie letters alphabetically
   public static final byte A = 0;   
   public static final byte C = 1;
   public static final byte G = 2;
   public static final byte T = 3;
   public static final byte GAP = 5; // in pairwise and mutiple alignment this represent a GAP in the alignment
   public static final byte N = 6; // any base
   public static final byte R = 7;
   public static final byte Y = 8;
   public static final byte S = 9;
   public static final byte W = 10;
   public static final byte K = 11;
   public static final byte M = 12;
   public static final byte B = 13;
   public static final byte d = 14;
   public static final byte H = 15;
   public static final byte V = 16;
   public static final byte UNKNOWN = -1;
   public static final byte D = -2; // base is deleted
   
   public static final double PVALUE_NORMAL = 0.05;
   public static final double PVALUE_STRICT = 0.15;
   
   public static final char[] TRIPLET_VAR = {'!','"','#','$','%','&','\'', '(',')','*','+',',','-','.','/','0','1','2','3','4','5','6','7','8','9', ':',';', '<','=','>','?','@','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','[','\\' ,']', '^','_', '`', 'a','b'}; // 64+1 characters   
   public static final char UPPER_BOUND = '~'; // ascii of this character is greator than the ones in TRIPLET_VAR
   
   public static final HashMap<String, String> codonMap = new HashMap();
   static{
       codonMap.put("TTT", "Phe");
   }
   
   public static int getMaximumIndex(int[] arr){
       // returns the index of the element in the array that has maximum number
       int idx = 0;
       int max = Integer.MIN_VALUE;
       for (int i=0; i<arr.length; i++){
           if (arr[i] > max){
               max=arr[i];
               idx = i;
           }
       }
       
       return idx;
   }
   
   public static byte nt_to_i(char nt){
           switch (nt){
                case 'A':
                case 'a':
                    return A;                       
                case 'T':                        
                case 't':
                case 'U':                        
                case 'u':                    
                    return T;                        
                case 'C':
                case 'c':
                    return C;                        
                case 'G':                        
                case 'g':
                    return G;   
                case '-':
                    return GAP;
                case 'n':
                case 'N':
                    return N;
                case 'd':
                case 'D':
                case '*':
                    return D; 
                case 'R':
                case 'r':
                    return A;
                case 'W':
                case 'w':
                    return A;
                case 'M':
                case 'm':
                    return A;
                case 'Y':
                case 'y':
                    return C;
                case 'S':
                case 's':
                    return C;
                case 'K':
                case 'k':
                    return G;
                case 'B':
                case 'b':
                    return C;
                case 'H':
                case 'h':
                    return A;
                case 'V':
                case 'v':
                    return A;
                    
                    
            }
           return UNKNOWN;       
    }  
   public static String relative_pos_string(byte rpos){
       switch(rpos){
           case POS_Q_Extended_S:
               return "Q extended by S";
           case POS_S_Extended_Q:
               return "S extended by Q";
           case  POS_Q_Inside_S:
               return "Q inside S";
           case POS_S_Inside_Q:
               return "S inside Q";
           case POS_Distinct:
               return "Distinct";
       }
       
       return "UNKNOWN";
   }
   public static int getTripletIndex(char triplet){
       for(int i=0; i<TRIPLET_VAR.length; i++)
           if (TRIPLET_VAR[i]== triplet)
               return i;
       
       return -1;
   }
   public static boolean isGood(String seq){
       // a seq is considered good if it is not homopolymere or contains any n or N
       if (seq.indexOf('n')>=0 || seq.indexOf('N')>=0)
           return false;
       //if (seq.matches("[Aa]+|[Tt]+|[Cc]+|[Gg]+"))
              // return false;
       return true;
   }
    public static int seqNNum(String seq){
       // accepts a seq and counts number of N and return 
       int cnt = 0;
       for (int i=0; i<seq.length(); i++){
           if(seq.charAt(i)=='N' ||  seq.charAt(i)=='n')
               cnt++;
       }
       
       return (cnt);
   }

   public static String getReverse(String seq){
       StringBuffer buf = new StringBuffer();
       for(int i=seq.length()-1; i>=0; i--)
           buf.append(seq.charAt(i));
       
       return buf.toString();
   }
   public static String getRC(String seq){
       // return reverse complement of a sequence (this method guarantees that the length of input parameter is the same as output) 
       StringBuffer buff = new StringBuffer();
       for(int i=seq.length()-1; i>=0; i--){
           char ch = seq.charAt(i);
           byte b = nt_to_i(ch);
           if (b==A)
               buff.append('T');
           else if (b==T)
               buff.append('A');
           else if (b==C)
               buff.append('G');
           else if (b==G)
               buff.append('C');
           else if (b==N)
               buff.append('N');
           else if (b==GAP)
               buff.append('-');           
           else if (b==UNKNOWN)
               buff.append('N');
       }
       return buff.toString();
   }
   public static String standardize(String seq){
       // return lower case and IUPAC code with 
       StringBuffer buff = new StringBuffer();
       for(int i=0; i<seq.length(); i++){
           char ch = seq.charAt(i);
           byte b = nt_to_i(ch);
           if (b==A)
               buff.append('A');
           else if (b==T)
               buff.append('T');
           else if (b==C)
               buff.append('C');
           else if (b==G)
               buff.append('G');
           else if (b==N)
               buff.append('N');
           else if (b==GAP)
               buff.append('-');           
           else if (b==UNKNOWN)
               buff.append('N');
       }
       return buff.toString();
   }
   
   public static String codonToAminoAcid(String codon){
       switch(codon){
           case "TTT":
           case "TTC":
           case "UUU":
           case "UUC":
               return "Phe";
           case "TTA":
           case "TTG":
           case "CTT":
           case "CTC":
           case "CTA":
           case "CTG":
           case "UUA":
           case "UUG":
           case "CUU":
           case "CUC":
           case "CUA":
           case "CUG":
               return "Leu";
           case "ATT":
           case "ATC":
           case "ATA":
           case "AUU":
           case "AUC":
           case "AUA":
               return  "Ile" ;
           case "ATG":
           case "AUG":
               return "Met";
           case "GTT":
           case "GTC":
           case "GTA":
           case "GTG":
           case "GUU":
           case "GUC":
           case "GUA":
           case "GUG":
               return "Val";
           case "TCT":
           case "TCC":
           case "TCA":
           case "TCG":
           case "UCU":
           case "UCC":
           case "UCA":
           case "UCG":
           case "AGT":
           case "AGC":
           case "AGU":              
               return "Ser";
           case "CCT":
           case "CCC":
           case "CCA":
           case "CCG":
           case "CCU":               
               return "Pro";
           case "ACT":
           case "ACC":
           case "ACA":
           case "ACG":
           case "ACU":
               return "Thr";
           case "GCT":
           case "GCC":
           case "GCA":
           case "GCG":
           case "GCU":
               return "Ala";
           case "TAT":
           case "TAC":
           case "UAU":
           case "UAC":
               return "Tyr";
           case "CAT":
           case "CAC":
           case "CAU":
               return "His";
           case "CAA":
           case "CAG":
               return "Gln";
           case "AAT":
           case "AAC":
               return "Asn";
           case "AAA":
           case "AAG":
               return "Lys";
           case "GAT":
           case "GAC":
           case "GAU":
               return "Asp";
           case "GAA":
           case "GAG":
               return "Glu";
           case "TGT":
           case "TGC":
           case "UGU":
           case "UGC":
               return "Cys";
           case "TGG":
           case "UGG":
               return "Trp";
           case "CGT":
           case "CGC":
           case "CGA":
           case "CGG":
           case "CGU":
           case "AGA":
           case "AGG":               
               return "Arg";
           case "GGT"  :
           case "GGC":
           case "GGA":
           case "GGG":
           case "GGU":
               return "Gly";
           case "TGA":
           case "TAA":
           case "TAG":
           case "UGA":
           case "UAA":
           case "UAG":
               return "***";
               
       }
       //unrcognized codon
       return "Err";
   }
   public static String translateDNA(String dna,  byte name_type) throws Exception{
       
       //this method converts DNA string to protein string, short or long name format; using BioJava methods
       //name_type 0: short name, name_type 1 : long name
       StringBuffer buff = new StringBuffer();
       ChromosomeSequence seq = new ChromosomeSequence(dna);
       GeneSequence gene = seq.addGene(new AccessionID("gene1"), 1, seq.getLength(), Strand.POSITIVE);
       TranscriptSequence transcript = gene.addTranscript(new AccessionID("transcript1"), 1, gene.getLength());
       transcript.addCDS(new AccessionID("CDS1"), 1, transcript.getLength(), 0);
       
       ProteinSequence prot = transcript.getProteinSequence();
       for (AminoAcidCompound aa : prot){
           if (name_type == 0)
              buff.append(aa.getShortName());
           else if (name_type == 1)
               buff.append(aa.getLongName());
       }
       
       return buff.toString();
   }
   public static String simplify(String seq){
       // remove gaps , replace any non-actgn characters with n
        StringBuffer buff = new StringBuffer();
        for(int i=0; i<seq.length(); i++){
           char ch = seq.charAt(i);
           byte b = nt_to_i(ch);
           if (b==A)
               buff.append('A');
           else if (b==T)
               buff.append('T');
           else if (b==C)
               buff.append('C');
           else if (b==G)
               buff.append('G');
           else if (b==N)
               buff.append('N');
           else if (b==UNKNOWN)
               buff.append('N');
            // if it is * D or - do nothing
        }
       return buff.toString();
   }
   public static String extractGeneID(String column9){
       // from column nine of a gff file
       String[] arr = column9.split(";");
       for (String el : arr){
           //if (el.startsWith("gene_id="))
           if (el.startsWith("ID="))
               return el.substring(3);
       }
       return null;
   }
   // another implementation using biojava
   public static String getReverseComplement(String seq){
        
            DNASequence dna = new DNASequence(seq);
            
            return dna.getReverseComplement().getSequenceAsString(); 
   }
   public static double getSimilarity(String seq1, String seq2){
        
        int cc = 0;
        for(int i=0; i<seq1.length(); i++)
            if (seq1.charAt(i)==seq2.charAt(i) || seq1.charAt(i)=='N' || seq2.charAt(i)=='N')
                cc++;
            
            return (double)cc/seq1.length();
   }
   public static String getN_BlockEquivalant(String seq){
        // this function, accepts something like NNNNNNNN, or NNNNACCTGGG or ACCNNNNNCGGCA
       //policy, if more than 90% of bases are N then returns empty string other wise trims N at one end and returns remaining
       //First Remove starting and ending N
       String left_trim = seq.replaceAll("^N+", "");
       String right_trim = left_trim.replaceAll("N+$", "");
       return right_trim;
      /* int non_N_start_idx = 0; // inclusive
       int non_N_end_idx = seq.length(); // exclusive
       
       if (seq.toUpperCase().startsWith("N"))
           for(int i=0; i<seq.length(); i++)
               if ()
       StringBuffer buff = new StringBuffer();
        int N_count = 0;
        for(char base : seq.toUpperCase().toCharArray()){
            if (base=='N')
                N_count++;
        }
        double ratio = (double)N_count/seq.length();
        if(ratio>=0.9)
            return "";
        return seq.toUpperCase().replaceAll("N+", "");   */
   }
   
   public static <T> T[] subarray(T[] arr, Class<T> cl, int start_idx, int end_idx){
      //start_idx is inclusive and end_idx is exclusive
       //Class<T> clazz;
       T[] out = (T[]) Array.newInstance(cl, end_idx - start_idx);
       int j = 0;
       for(int i=start_idx; i<end_idx; i++)
           out[j++] = arr[i];
       
       return out;
  }
   
   public static String fixLength(String in, char ch, int len) throws Exception{
       // this methods , adds ch to in if in length is less than len
       if (in.length()>len)
           throw new Exception("Exception inside fixLength");
       if (in.length()== len)
           return in;
       StringBuffer out = new StringBuffer(in);
       for(int i=0; i<len-in.length(); i++)
           out.append(ch);
       
       return out.toString();
   }
   public static byte tripletToDecimal(String triplet){
       // accepts a triplet of ACGT and based on chartai (4 base) produces a number from 0-63, if triplet contains N then it returns 64
       byte b = 0;
       // first make sure ther is no N inside 
       //triplet.toLowerCase().replaceAll("n", "");
       if (triplet.toLowerCase().indexOf('n')>=0)
           return 64; // this means triplet contains N, so perhaps not very reliable
       // we get here if triplet has no N at all
       for (int i=0; i<3; i++)
           b += nt_to_i(triplet.charAt(i))* Math.pow(4,2-i);
       
       return b;
   }
   
   public static short quadletToDecimal(String quadlet){
       // accepts a quadlet of ACGT and based on chartai (4 base) produces a number from 0-255, if quadlet contains N then it replaces it with base A
       short sh = 0;
       for (byte i=0; i<4; i++){
           byte b = 0;
           if (quadlet.charAt(i)=='N' || quadlet.charAt(i)=='n')
               b = 0; // treated as base A
           else
               b = nt_to_i(quadlet.charAt(i));
           sh += b* Math.pow(4,3-i);
       }
       
       return sh;
   }
   public static short seqCNum(String seq){
       // accepts a seq and counts number of C and return a number from 0-255
       short cnt = 0;
       for (short i=0; i<seq.length(); i++){
           if(seq.charAt(i)=='C' ||  seq.charAt(i)=='c')
               cnt++;
       }
       
       return (cnt<256?cnt:255);
   }
   
   public static String decimalToQuadlet(char ch){
       // accepts a char or a number and converts it to base 4 quadlet (A:0, C:1, G:2, T:3)
       StringBuffer out = new StringBuffer();
       short sh = (short)ch;
       String quad = Integer.toString(sh, 4);
       // then fill up quad with leading 0 if its length is less than 4  
       for (byte i=0; i<4-quad.length(); i++)
           out.append('A'); // A is 0
       
       for(byte i=0; i<quad.length(); i++){
           switch(quad.charAt(i)){
               case '0':
                   out.append('A');
                   break;
               case '1':
                   out.append('C');
                   break;
               case '2':
                   out.append('G');
                   break;
               case '3':
                   out.append('T');
                   break;
           }
       }
           
       return out.toString();
   } 
   
   public static Short octaToDecimal(String octa){
       // accepts a Octalet of ACGT and based on base 4  produces a number from 0-63
       short b = 0;
       // first make sure ther is no N inside 
       //triplet.toLowerCase().replaceAll("n", "");
       if (octa.toLowerCase().indexOf('n')>=0)
           return null; // this means triplet contains N, so perhaps not very reliable
       // we get here if octa has no N at all
       for (int i=0; i<8; i++)
           b += nt_to_i(octa.charAt(i))* Math.pow(4,i);
       
       return b;
   }
   
   public static List<Short> getReadPrefixOcta(String read){
       
       // work out number of buckets
       int octas = read.length()/8;
       List<Short>  list = new ArrayList();
       
       
       for(int i=0;i<octas; i++){
           // i is index of each octa
           list.add(octaToDecimal(read.substring(i*8, (i+1)*8))) ;
       }
       
       return list;
   }
   
   public static byte[] getReadDigest(String read){
       // accepts a read as string and create a digest (breaks the read into 15 base buckets and encode first triplet of each bucket and writes it into a byte)
       int len = read.length();
       // work out number of buckets
       int buckets = len/15;
       if (len%15 >= 3)
           buckets++;
       
       byte[] output = new byte[buckets];
       
       for(int i=0;i<buckets; i++){
           // i is index of each bucket
           output[i] = tripletToDecimal(read.substring(i*15, i*15+3));
       }
       
       return output;
   }
   
   public static List<Byte> getReadDigest(String read, int bucket_len){
       // accepts a read as string and create a digest (breaks the read into bucket_len base buckets and encode first triplet of each bucket and writes it into a byte)
       int len = read.length();
       // work out number of buckets
       int buckets = len/bucket_len;
       if (len%bucket_len >= 3)
           buckets++;
       
       List<Byte> list = new ArrayList();
       
       for(int i=0;i<buckets; i++){
           // i is index of each bucket
           list.add(tripletToDecimal(read.substring(i*bucket_len, i*bucket_len+3))) ;
       }
       
       return list;
   }
   
   public static List<Short> getReadDigestOcta(String read, int bucket_len){
       // accepts a read as string and create a digest (breaks the read into bucket_len base buckets and encode first octa of each bucket and writes it into a Short)
       int len = read.length();
       // work out number of buckets
       int buckets = len/bucket_len;
       if (len%bucket_len >= 8)
           buckets++;
       
       List<Short> list = new ArrayList();
       
       for(int i=0;i<buckets; i++){
           // i is index of each bucket
           list.add(octaToDecimal(read.substring(i*bucket_len, i*bucket_len+8))) ;
       }
       
       return list;
   }
   
   public static String getReadDigest2(String read, int bucket_len){
       // accepts a read as string and create a digest (breaks the read into 15 base buckets and encode first triplet of each bucket and writes it into a character)
       
       int len = read.length();
       // work out number of buckets
       int buckets = len/bucket_len;
       if (len%bucket_len >= 3)
           buckets++;
       
       char[] chars = new char[buckets];
       
       for(int i=0;i<buckets; i++){
           // i is index of each bucket
           chars[i] = TRIPLET_VAR[tripletToDecimal(read.substring(i*bucket_len, i*bucket_len+3))];
       }
       
       return new String(chars);
   }
   
    public static String getPrefixDigest(String read, int bucket_len){
       // accepts a read as string and create a digest (breaks the read into buckets and encode first quadlet of each bucket and writes it into a character (0-255))
       
       int len = read.length();
       // work out number of buckets
       int buckets = len/bucket_len;
       if (len%bucket_len >= 4)
           buckets++;
       
       char[] chars = new char[buckets];
       
       for(int i=0;i<buckets; i++){
           // i is index of each bucket
           chars[i] = (char)quadletToDecimal(read.substring(i*bucket_len, i*bucket_len + 4));
       }
       
       return new String(chars);
   }
    public static String getPrefixDigest2(String read, int bucket_len){
       // accepts a read as string and create a digest (breaks the read into buckets and encode last quadlet of each bucket and writes it into a character (0-255))
       
       int len = read.length();
       // work out number of buckets
       int buckets = len/bucket_len;
       //if (len%bucket_len >= 4)
       //    buckets++;
       
       char[] chars = new char[buckets];
       
       for(int i=1;i<=buckets; i++){
           // i is index of each bucket
           chars[i-1] = (char)quadletToDecimal(read.substring(i*bucket_len - 4, i*bucket_len));
       }
       //because last bucket may have a lower length of bucket_len, we just grab the last 4 bases of the read
       //chars[buckets] = (char)quadletToDecimal(read.substring(len - 4, len));
       
       return new String(chars);
   }
    public static String getFullDigest(String read, int bucket_len){
       // accepts a read as string and create a digest (breaks the read into buckets and encode some bases from middle  and writes it into a character (0-255))
       
       int len = read.length();
       // work out number of buckets
       int buckets = len/bucket_len;
       //if (len%bucket_len >= 4)
       //    buckets++;
       
       char[] chars = new char[buckets];
       int l = (bucket_len-4)/4;
       for(int i=0;i<buckets; i++){
           // i is index of each bucket
           //start from 5'th base, apart from the first 4 bases of each bucket (first 4 bases were already used for Building Trie), divide the rest into 4 parts and grab the first base in each part and constitute 4 bases
           
           String bucket = read.substring(i*bucket_len , (i+1)*bucket_len);
           
           chars[i] = (char)quadletToDecimal(new StringBuffer().append(bucket.charAt(4)).append(bucket.charAt(4+l)).append(bucket.charAt(4+(2*l))).append(bucket.charAt(4+(3*l))).toString());
       }
       //because last bucket may have a lower length of bucket_len, we just grab the last 4 bases of the read
       //chars[buckets] = (char)quadletToDecimal(read.substring(len - 4, len));
       
       return new String(chars);
   }
    public static String getRatioDigest(String read, int bucket_len, int quadlets){
       // accepts a read as string and create a digest (breaks the read into buckets and encode quadlets from quadlet 1 until quadlets )
       
       int len = read.length();
       // work out number of buckets
       int buckets = len/bucket_len;
       //if (len%bucket_len >= 4)
       //    buckets++;
       
       char[] chars = new char[buckets*quadlets];
       int j = 0; //index of chars
       for(int i=0;i<buckets; i++){
           // i is index of each bucket
           //start from 5'th base, apart from the first 4 bases of each bucket (first 4 bases were already used for Building Trie), divide the rest into 4 parts and grab the first base in each part and constitute 4 bases
           
           String bucket = read.substring(i*bucket_len , (i+1)*bucket_len);
           //for each bucket we read this number of quadlet starting from second quadlet
           for(int k=1; k<=quadlets;k++)
              chars[j++] = (char)quadletToDecimal(bucket.substring(k*4, (k+1)*4));
       }
       //because last bucket may have a lower length of bucket_len, we just grab the last 4 bases of the read
       //chars[buckets] = (char)quadletToDecimal(read.substring(len - 4, len));
       
       return new String(chars);
   }

    public static String getFullEncode(String read){
       // accepts a read as string and create a digest (breaks the read into buckets and encode quadlets from quadlet 1 until quadlets )
        //because length of sequence may not be dividable by 4 , then last remaining bases are stored as they are without encoding
        // so for example if length of sequence is 66 then we need 16 chars to store 64 bases and 2 more chars for last 2 bases (18 in total)
       
       int len = read.length();
       // work out number of buckets
       int quadlets = len/4;
       //int remainder = len%4;
       
       //char[] chars = new char[quadlets+remainder];
       char[] chars = new char[quadlets];
       for(int i=0;i<quadlets; i++){
            chars[i] = (char)quadletToDecimal(read.substring(i*4, (i+1)*4));
       }
       /*for(int i=0;i<remainder; i++){
            chars[quadlets+i] = read.charAt(len-remainder+i);
       }*/
       
       //because last bucket may have a lower length of bucket_len, we just grab the last 4 bases of the read
       return new String(chars);
   }
   public static boolean encodedReadStartsWith(String encoded_read, String read_remainder, String suffix){
       //this method is complement to getFullEncode, it returns true if a read comprises encoded_read+read_remainder starts with suffix
       StringBuffer buff = new StringBuffer();
       //first decode the read
       for(int i=0; i<encoded_read.length(); i++)
           buff.append(decimalToQuadlet(encoded_read.charAt(i)));       
       buff.append(read_remainder);
       return buff.toString().startsWith(suffix);
       
   } 
   public static boolean encodedReadStartsWith(String encoded_read, String read_remainder, String suffix, int mismatch){
       //oveloaded method, allow for mismatches
       //this method is complement to getFullEncode, it returns true if a read comprises encoded_read+read_remainder starts with suffix
       StringBuffer buff = new StringBuffer();
       //first decode the read
       for(int i=0; i<encoded_read.length(); i++)
           buff.append(decimalToQuadlet(encoded_read.charAt(i)));       
       buff.append(read_remainder);
       String read = buff.toString();
       if(read.length()<=suffix.length())
           return false;
       if (hammingDistance(read.substring(0, suffix.length()), suffix) > mismatch)
         return false;
       
       return true;
       
   }
   public static int hammingDistance(String str1, String str2){
       //find the hamming distance and return it, we assume str1 and str2 have the same length
        if (str1.length() != str2.length())
           return -1;
       
        int dist = 0;
        for(int i=0; i<str1.length(); i++){
           char ch1 = str1.charAt(i);
           char ch2 = str2.charAt(i);
           if(ch1!=ch2 && ch1!='N' && ch2!='N')
               dist++; //so for example if at least one of chars is not N then distance will not increase
        }
       
        return dist;
       
   }
   
    public static String decodeRead(String encoded_read, String read_remainder){
       //oveloaded method, allow for mismatches
       //this method is complement to getFullEncode, it returns true if a read comprises encoded_read+read_remainder starts with suffix
       StringBuffer buff = new StringBuffer();
       //first decode the read
       for(int i=0; i<encoded_read.length(); i++)
           buff.append(decimalToQuadlet(encoded_read.charAt(i))); 
       if (read_remainder!=null && !read_remainder.isEmpty())
          buff.append(read_remainder);       
       return buff.toString();
       
    }   
    public static String getPrefixCNum(String read, int bucket_len){
       // accepts a read as string and create a digest (breaks the read into buckets and counts number of CG:)
       
       int len = read.length();
       // work out number of buckets
       int buckets = len/bucket_len;
       //if (len%bucket_len >= 4)
       //    buckets++;
       
       char[] chars = new char[buckets];
       
       for(int i=0;i<buckets; i++){
           // i is index of each bucket
           chars[i] = (char)seqCNum(read.substring(i*bucket_len , (i+1)*bucket_len));
       }
       //because last bucket may have a lower length of bucket_len, we just grab the last 4 bases of the read
       //chars[buckets] = (char)quadletToDecimal(read.substring(len - 4, len));
       
       return new String(chars);
   }
    
   public static String getReadSuffixTriplets(String read){
       // accepts a read as string and breaks down read to triplet bases and encode each triplet, if for the last triplet is one or 2 bases left we just ignore it
       
       int len = read.length()-1; // excluding first base (we just need suffixes)
       // work out number of buckets
       int triplets = len/3;
       char[] chars = new char[triplets];
       
       for(int i=0;i<triplets; i++){
           // i is index of each triplet
           chars[triplets-i-1] = TRIPLET_VAR[tripletToDecimal(read.substring(len-((i+1)*3), len-(i*3)))];
       }
       
       return new String(chars);
   }
   
   public static List<Byte> getReadPrefixTriplets(String read){
       // accepts a read as string and breaks down read to triplet bases and encode each triplet, if for the last triplet is one or 2 bases left we just ignore it
       
       int len = read.length(); // excluding last base (we just need prefixes)
       // work out number of buckets
       int triplets = len/3;
       List<Byte> list = new ArrayList();
       
       for(int i=0;i<triplets; i++){
           // i is index of each triplet
           list.add(tripletToDecimal(read.substring(i*3, (i*3)+3)));
       }
       
       return list;
   }
   
   public static int getGaps(String alignment){
       // it returns number of - exist
       int gaps = 0;
       for (int i=0; i<alignment.length(); i++){
           if (alignment.charAt(i)=='-')
               gaps++;
       }
       return gaps;
   }
   public static TreeMap getGapsMap(String seq){
       // this will return position of the next characer to one or multiple - , for example if sequence is AA--TC-AAT-G then will return 4->2, 7->1, 11->1
       TreeMap<Integer, Integer> out = new TreeMap();
       Pattern p = Pattern.compile("\\-+"); // possessive(+)
        Matcher   m = p.matcher(seq); // remember that qseq contains Gaps or -, but we only look for continuous base without gap in the middle
        while(m.find()){
            // extract the equal part from reference (subject)
            out.put(m.end(), m.end()-m.start());
        } 
        
        return out;
       
   }
   public static TreeMap getGapsMap2(String seq){
       // this will return position of gaps in no gap forma of sequence , for example if sequence is AA--TC-AAT-G then will return 2->2, 4->1, 7->1
       TreeMap<Integer, Integer> map = new TreeMap();
       Pattern p = Pattern.compile("\\-+"); // possessive(+)
        Matcher   m = p.matcher(seq); // remember that qseq contains Gaps or -, but we only look for continuous base without gap in the middle
        while(m.find()){
            // extract the equal part from reference (subject)
            map.put(m.start(), m.end()-m.start());
        } 
        //at this stage out contains 2->2, 6->1, 10->1
        TreeMap<Integer, Integer> out = new TreeMap();        
        int sum = 0;
        for(Entry<Integer,Integer> ent : map.entrySet()){
            out.put(ent.getKey()-sum, ent.getValue());
            sum += ent.getValue();
        }
   //2->2, 4->1,  7->1     
        return out;
       
   }
   
   public static String MT_String(byte mt){
       String str = "";
       switch(mt){
           case MT_DELETION:
               str = "Deletion";
               break;
           case MT_INSERTION:
               str = "Insertion";
               break;
           case MT_Substitution:
               str = "Substitution";
               break;
           case MT_Default:
               str = "Default";
               break;
           case MT_Hidden:
               str = "Hidden";
               break;
               
       }
       return str;
   }
  
   public static String MT_Short_String(byte mt){
       String str = "";
       switch(mt){
           case MT_DELETION:
               str = "DEL";
               break;
           case MT_INSERTION:
               str = "INS";
               break;
           case MT_Substitution:
               str = "SUB";
               break;
           case MT_Default:
               str = "DEF";
               break;
           case MT_Hidden:
               str = "HID";
               break;
               
       }
       return str;
   }


   public static byte MT_String_Byte(String mut){
       if(mut.toLowerCase().startsWith("del"))
           return MT_DELETION;
       if(mut.toLowerCase().startsWith("ins"))
           return MT_INSERTION;
       if(mut.toLowerCase().startsWith("tr"))
           return MT_Substitution;
       if(mut.toLowerCase().startsWith("def"))
           return MT_Default;
       if(mut.toLowerCase().startsWith("hid"))
           return MT_Hidden;
       
       return -1;
       
   }
    public static String VR_Short_String(byte vt){
       String str = "";
       switch(vt){
           case VR_DELETION:
               str = "DEL";
               break;
           case VR_INSERTION:
               str = "INS";
               break;
           case VR_SNP:
               str = "SNP";
               break;
           case VR_Default:
               str = "DEF";
               break;
               
       }
       return str;
   }

   public static String ZG_Short_String(byte zg){
       String str = "";
       switch(zg){
           case ZG_UNKNOWN:
               str = "UNKNOWN";
               break;
           case ZG_HOMOZYGOUS:
               str = "HOM";
               break;
           case ZG_HETROZYGOUS:
               str = "HET";
               break;
           case ZG_INVALID:
               str = "ERR";
               break;           
       }
       return str;
   }
   
   public static String ZG_String(byte zg){
       String str = "";
       switch(zg){
           case ZG_UNKNOWN:
               str = "UNKNOWN";
               break;
           case ZG_HOMOZYGOUS:
               str = "HOMOZYGOUS";
               break;
           case ZG_HETROZYGOUS:
               str = "HETEROZYGOUS";
               break;
           case ZG_INVALID:
               str = "ERROR";
               break;           
       }
       return str;
   }
   
   public static String getVCFAllelesStr(VariantContext vc, int sample_index) throws Exception{
       if (vc.getGenotypes().get(sample_index).getAlleles().size() == 1)// its haploid sucha ax chrX and chrY
           return vc.getGenotypes().get(sample_index).getAllele(0).getBaseString().toUpperCase();
           //throw new Exception("Genotype/Specie is not diploid at contig " + vc.getContig() + ":" + vc.getStart());
       if (vc.getGenotypes().get(sample_index).getAllele(0).getBaseString().equals(vc.getGenotypes().get(sample_index).getAllele(1).getBaseString()))
           return vc.getGenotypes().get(sample_index).getAllele(0).getBaseString().toUpperCase();// homozygous
       else
           return (vc.getGenotypes().get(sample_index).getAllele(0).getBaseString().toUpperCase()+ "/" +vc.getGenotypes().get(sample_index).getAllele(1).getBaseString().toUpperCase());// heterozygous
   }
   public static String getVCFAlternateStr(VariantContext vc) throws Exception{
       StringBuffer buff = new StringBuffer();
       for(htsjdk.variant.variantcontext.Allele al : vc.getAlternateAlleles()){
           buff.append(al.getBaseString().toUpperCase()+"/");
       }
       if (buff.charAt(buff.length()-1)=='/')
           buff.deleteCharAt(buff.length()-1);
       return buff.toString(); 
   }
   
   public static byte ZG_String_Byte(String zg){
       if(zg.toLowerCase().startsWith("hom"))
           return ZG_HOMOZYGOUS;
       if (zg.toLowerCase().startsWith("het"))
           return ZG_HETROZYGOUS;
       return ZG_UNKNOWN;
   }
   public static String Base_String(byte b){
       switch(b){
           case A:
               return "A";
           case T:
               return "T";
           case C:
               return "C";
           case G:
               return "G";
            case N:
               return "N";
            case GAP:
                return "-";
            case D:
                return "D";
            case UNKNOWN:
                return "";
            case R:
                return "A";
            case Y:
                return "C";
            case S:
                return "C";
            case W:
                return "A";
            case K:
                return "G";
            case M:
                return "A";
            case B:
                return "C";
            case H:
                return "AH";
            case V:
                return "A";
                
       }
       
       return " ";
   }
  
   
   public static int scale_to(int scale_to, int cov, int max_cov){
      if (cov==0)
           return 0;

       int new_cov = 0; 
       if(scale_to>1){
           if (cov>=max_cov)
               new_cov = scale_to;
           else{
             int interval = max_cov/scale_to; // it return an integer 
             new_cov = cov/interval + 1;
           }
       }else if (scale_to==1){
           // in this mode we have binary mode, if cov is above max_cov it returns 1 otherwise it is 0
           if(cov>max_cov)
               new_cov = 1;
           else
               new_cov = 0;
       }
       
       
            
       return new_cov;
   }
   public static List<Byte> getProfileColumn(List<String> alignment, int col){
      List<Byte> out = new ArrayList() ;
      for(String alignedseq : alignment)
          out.add(nt_to_i(alignedseq.charAt(col)) );
      
      return out;
   }
   
   public static double phredToProb(int phred){
       return Math.pow(10, (double)phred/-10);
   }
   
   public static double getMAPQ(char ch){
       int ascii = ch; // ch is a character that represent quality of the map
       int phred = ascii - 33;
       return phredToProb(phred);
   }
   public static double getPhred(char ch){
       int ascii = ch; // ch is a character that represent quality of the map
       int phred = ascii - 33;
       return phred;
   }
   public static double getMAPQScale(char ch, int mean_mapq){
       int ascii = ch; // ch is a character that represent quality of the map
       int phred = ascii - 33;
       
       if (phred >= mean_mapq)
           return 1;
       
       double d = (1-phredToProb(phred))/(1-phredToProb(mean_mapq));
       return d;
   }
   
   public static int goToNextInteger(int idx, String[] arr){
       //from index position in the arr sacns the array to get to the first Integer, returns -1 if nothing has been found
       int i = 0;
       for (i=idx; i< arr.length; i++)
           if (arr[i].trim().matches("\\d+"))
               break;
       
       
       if (i==arr.length)
           return -1;
       
       return i;
   }
   public static int getReferencePosition(String polymorphicSequence, int idx, int start_pos){
       //returns reference position of the base at idx (zero based) position in the sequence: AACCDDNNNNNNNNC[+A/AA/AAA]C:ACTT[CD][AC]AADGCNNNNNNNNN
        boolean openBracket = false;
        int pos = start_pos;
        for (int i=0; i<polymorphicSequence.length();i++){
           if (i == idx)
               return pos; 
           char ch = polymorphicSequence.charAt(i);
           switch (ch){
               case '[':
                   openBracket = true;
                   if (polymorphicSequence.charAt(i+1) != '+')
                       pos++;
                   break;
               case ']':
                   openBracket = false;
                   break;
               case 'N':
               case 'A':
               case 'T':
               case 'C':
               case 'G':   
               case 'D':
                   if (!openBracket)
                       pos++;
           }//switch
       } //for
       return -1;
   }
   
   public static int getCharacterCount(char ch, String str){
       // counts number of characters in the string
       int cc = 0;
       for (int i=0; i<str.length(); i++)
           if(str.charAt(i)==ch)
               cc++;
       
       return cc;
   }
   public static int getStartRepeat(String str,int min_len ){
       // if the same character repeats at least min_len time in the begining of the string it will return length of repeat otherwise retunrs 0
       if (str.isEmpty())
           return 0;
       
       int cc = 1;
       char start_ch = str.charAt(0);
       for (int i=1; i<str.length(); i++)
           if(str.charAt(i)==start_ch)
               cc++;
           else 
               break;
       if(cc>= min_len)
         return cc;
       else
           return 0;
   }
   public static int getEndRepeat(String str,int min_len ){
       // if the same character repeats at least min_len time at the end of the string it will return length of repeat otherwise retunrs 0
       if (str.isEmpty())
           return 0;
       
       int cc = 1;
       char end_ch = str.charAt(str.length()-1);
       for (int i=str.length()-2; i>=0; i--)
           if(str.charAt(i)==end_ch)
               cc++;
           else 
               break;
       if(cc>= min_len)
         return cc;
       else
           return 0;
   }
   
   public static String removeLeadingTrailingN(String str){
       // if a sequence is like NNACCCTTNNNNN it will return ACCCTT
       
       int i = 0; //left index 
       for ( ; i<str.length(); i++){
           if (str.charAt(i)!='N')
               break;
       }
       //so when we get here i points to first non-N character
       int j = str.length()-1; //right index 
       for ( ; j>=0; j--){
           if (str.charAt(j)!='N')
               break;
       }
       //j points to first non N character from right hand side
       if (i<=j)
           return str.substring(i, j+1);       
       else
           return "";
   }
   public static Integer sumarr(Collection<Integer> arr){
       int sum = 0;
       for(int n : arr)
           sum += n; 
       return sum;
   }
   
   public static int sequenceDistance(String seq1, String seq2){       
       //this method returns sequence distance between 2 sequences of the same length
       int dist = 0;
       for (int i=0; i<seq1.length(); i++)
           if (nt_to_i(seq1.charAt(i))!= nt_to_i(seq2.charAt(i)) )
               dist++;
       
       return dist; // distance zero means 2 sequences are the same
   }
   
   
   public static List<SNP> getSNPSubList(List<SNP> list, int pos_start, int pos_end){
       // process the list and return part of that (sublist) located in the range from start_pos until pos_end
       //This is very inefficient and takes a long time, use next method
       List<SNP> subList = new ArrayList();
       int sz = list.size();
       if (sz<2)
           return null;
       int pos1 = list.get(0).position;
       int posn = list.get(sz-1).position;
       if (pos_end<pos1)
           return null;
       if (pos_start>posn)
           return null;    
       int start_search = 0;
       //for(int i=start_search)
       for(SNP snp : list){
           if (snp.position>= pos_start && snp.position<=pos_end)                
               subList.add(snp.clone()); // Note: we create a new SNP object for every element of the subList, and we didnt not use the element of main list, why? because our dag changes alleles of SNPs 
       }
       return subList;           
   }
   
      public static List<SNP> getSNPSubList(int start_search_idx, List<SNP> list, int pos_start, int pos_end){
       // process the list and return part of that (sublist) located in the range from start_pos until pos_end
       List<SNP> subList = new ArrayList();
       int sz = list.size();
       if (sz<2)
           return null;
       int pos1 = list.get(0).position;
       int posn = list.get(sz-1).position;
       if (pos_end<pos1)
           return null;
       if (pos_start>posn)
           return null;    
       for(int i=start_search_idx; i<sz; i++){
           SNP snp = list.get(i);
           if (snp.position>= pos_start && snp.position<=pos_end)                
               subList.add(snp.clone()); // Note: we create a new SNP object for every element of the subList, and we didnt not use the element of main list, why? because our dag changes alleles of SNPs 
           else if (snp.position>pos_end)
               break; // we passed by end of range
       }
       return subList;           
   }

      public static Set<Integer> getSNPPositions(List<SNP> list){
       // process the list and return part of that (sublist) located in the range from start_pos until pos_end
       Set<Integer> positions = new HashSet();
       for(SNP snp:list){
           for (int i=snp.position; i<snp.position+snp.ref.length();i++)
               positions.add(i);
       }
       return positions;           
   }

   
   public static void detectBroken(List<SNP> snpList,int max_insert_size){
       int last_snp_pos = 0; // this variable is used to determine breaking point in SNP list
       for(SNP snp: snpList){
            boolean broken = false;
            if (max_insert_size>0 && last_snp_pos>0 && (snp.position - last_snp_pos) > max_insert_size)
                snp.broken = true;
           
           last_snp_pos = snp.position;
       }
   }
   public static void setSNPAlleleFromRead(String read, Cigar ncg, int read_start, SNP snp){
       // Note : This method 2 times was verified and tested, do not change it !
       // for each snp in snpList, see if a snp of type snp.type exist at snp.position then we set snp.allele based on the read and its cigar values.
       // Note: no record will start with I operator in cigar (we tested it by SAMFileInfo.check_insertion_start() method)
       // This method will set snp.allele to A|C|T|G|N|D|+  or +Inserted Bases
       int ref_pos = read_start; // a pointer to reference position
       int read_idx = 0; // a pointer to read position (0 based)
       int cigar_idx = 0; // index for each Cigar Element
       //int list_idx = 0; // always points to the current element in the snpList (we dont loop based on elements of snpList)
       int cigar_size = ncg.numCigarElements();
       try{
       while (cigar_idx<cigar_size){
            CigarElement ce = ncg.getCigarElement(cigar_idx);
            for (int j=0; j<ce.getLength(); j++){
                if (ref_pos == snp.position){
                   // here depending on the type of SNP, we set the value of SNP and return
                   if (snp.type == Helper.SNP_MultiBase){
                       // see what we have in the read (a base or a deletion)
                       if (ce.getOperator().equals(CigarOperator.M)){
                           snp.allele = Helper.Base_String( Helper.nt_to_i(read.charAt(read_idx)) ) ;
                           return;
                       }else if (ce.getOperator().equals(CigarOperator.D)){
                           snp.allele = "D";
                           return;
                       }                       
                   }else if (snp.type == Helper.SNP_MultiInsert){
                       // first we see if we are at the last element of current CigarElement
                       if (j==ce.getLength()-1){
                           if (cigar_idx < cigar_size-1){
                               CigarElement ce_next = ncg.getCigarElement(cigar_idx + 1);
                               if (ce_next.getOperator().equals(CigarOperator.I)){
                                   int insertion_length = ce_next.getLength();
                                   //we read next insertion_length bases from the read
                                   snp.allele = "+" + read.substring(read_idx + 1, read_idx + 1 + insertion_length);// SNP was resolved
                                   return;
                               }else{
                                   // if next cigar element is not an Insertion then allele will be just +
                                    snp.allele = "+" ;  // + allele for MultiInsert SNP means absence of an Insertion in the read (SNP was resolved successfully).
                                    return;                                                      
                               }
                            }else{
                               // we can get here when insertion is at the end of the read then allele will be unresolved, and will be resolved by another read
                                    // do nothing (allele will remain null)                                                   
                           }
                       }else{                           
                            snp.allele = "+" ;  // + allele for MultiInsert SNP means absence of an Insertion in the read (SNP was resolved successfully).
                            return;                                                      
                       }
                        //here we get next cigar element
                   }// if multi insert
               } // if ref_pos == snp.position
                // consume a base from the read if cigar element is not D
                if (!ce.getOperator().equals(CigarOperator.D))
                  read_idx++;
                // consume a base from the reference if cigar element is not I
                if (!ce.getOperator().equals(CigarOperator.I))
                  ref_pos++;
            }// j for each base in a Cigar Element
            cigar_idx++;
       } //while Cigar element
       
   
       }catch(Exception ex){
                   System.out.println("exception inside setSNPAlleleFromRead" + ex.getMessage() + "  snp:" + snp.toString());
       }
   }
   public static Set<Integer> getReadSNPPositions(Cigar ncg, String md, int start_ref_pos){
       //This method returns all the Ref SNP positions (Including InDels) .
       //InDels can be detected using Cigar string and SNPs can be detected from MD string.
       int step = 0;
       Set<Integer> positions = new HashSet();

       try{
           Set<Character> bases = new HashSet(Arrays.asList('A','a','C','c','G','g','T','t'));
           //first use Cigar string to work out InDel positions
           int ref_pos = start_ref_pos;
           for(CigarElement ce : ncg.getCigarElements()){
               if(ce.getOperator().equals(CigarOperator.I))
                   positions.add(ref_pos);
               else if (ce.getOperator().equals(CigarOperator.D))
                   positions.add(ref_pos);

               if (!ce.getOperator().equals(CigarOperator.I))
                  ref_pos += ce.getLength();
           }
step = 1;
           // Now use MD string to workout SNP positions
           if (md == null)
              return positions;
         //  System.out.println(ncg.toString() + "::" + md);

           int md_pos = 0;
           ref_pos = start_ref_pos;
           while(md_pos < md.length()){
               char ch = md.charAt(md_pos);
               if (bases.contains(ch)){
                   //we read from md until we get to a non-ACGT character
                   positions.add(ref_pos);
                   md_pos++;
                   ref_pos++;
               }else if (Character.isDigit(ch)){
                   //we read until we get to a non digit character
                   StringBuffer buff = new StringBuffer();

                   while(Character.isDigit(ch)){
                       buff.append(ch);
                       md_pos++;
                       //exit loop condition can not go into while statement
                       if(md_pos == md.length())
                           break;
                       ch = md.charAt(md_pos);
                   }
                   ref_pos += Integer.parseInt(buff.toString());
                   // we get here when ch is not digit, or end of md string
                   if(md_pos == md.length())
                       break;

               }else if (ch=='^'){
                   // this is a deletion start, we read coming ACGT until we get to a non-ACGT character
                   md_pos++;
                   ch = md.charAt(md_pos);
                   while( bases.contains(ch)){
                       ref_pos++;
                       md_pos++;
                       if(md_pos == md.length())
                           break;

                       ch = md.charAt(md_pos);

                   }
                   if(md_pos == md.length())
                           break;
               }else{
                   //unexpected situation, like N or IUPAC codes
                   ref_pos++;
                   md_pos++;
               }
           }//while(md_pos < md.length())
step = 2;
       }catch(Exception ex){
          System.out.println("Exception raised inside getReadSNPPositions step " + step + " Cigar " + ncg.toString() + " MD " +md);
           
       }
       return positions;
   }
  public static String getSNPAlleleFromRead(String read, Cigar ncg, int read_start, int snp_start, int snp_end){
       //developed 25/6/2019
       // Note : This method 2 times was verified and tested, do not change it !
       //snp_start and end are both 1-based and inclusive
       //On 5/3/2020 it was turned out that it can not handle funny situation like 1D1I   , so i ammended code
       // xIyD never happens in terms of aligner, but DI could happen, for example instead of 1M aligner might reposrt 1D1I, and only because of meeting penalty score threshold.
       //keep in mind that this function may return empty string as allele if snp_start:snp_end falls in a D element 
       StringBuffer buff = new StringBuffer();
       int ref_pos = read_start; // a pointer to reference position
       int read_idx = 0; // a pointer to read position (0 based)
       int cigar_idx = 0; // index for each Cigar Element
       //int list_idx = 0; // always points to the current element in the snpList (we dont loop based on elements of snpList)
       int cigar_size = ncg.numCigarElements();
       try{
       boolean exit_loop = false;
       while (cigar_idx<cigar_size){
            CigarElement ce = ncg.getCigarElement(cigar_idx);
            for (int j=0; j<ce.getLength(); j++){
                if (ref_pos >= snp_start){
                   // here depending on the type of SNP, we set the value of SNP and return
                   if (!ce.getOperator().equals(CigarOperator.D))
                       buff.append(read.charAt(read_idx));
                   // see if we are at the end of current cigar element and also at the end of reference postion
                    if (ref_pos== snp_end && j==ce.getLength()-1)
                       if (cigar_idx < cigar_size-1){
                           CigarElement ce_next = ncg.getCigarElement(cigar_idx + 1);
                           if (ce_next.getOperator().equals(CigarOperator.I)){
                               int insertion_length = ce_next.getLength();
                               //we read next insertion_length bases from the read
                               //added logic on 5/3/2020 
                               if(ce.getOperator().equals(CigarOperator.D))
                                  buff.append(read.substring(read_idx , read_idx  + insertion_length));// then it is like 1D1I , does it work for 2D3I  ? It should
                               else 
                                   buff.append(read.substring(read_idx + 1, read_idx + 1 + insertion_length));// like 50M2I  
                           }
                    }
                   
                } 
                
                if (!ce.getOperator().equals(CigarOperator.D))
                  read_idx++;
                // consume a base from the reference if cigar element is not I
                if (!ce.getOperator().equals(CigarOperator.I))
                  ref_pos++;
                if (ref_pos > snp_end){
                  exit_loop = true;
                  break;
                }

            }// j for each base in a Cigar Element
            if(exit_loop)
                break;
            cigar_idx++;
       } //while Cigar element
       
   
       }catch(Exception ex){
                   System.out.println("exception inside getSNPAlleleFromRead" + ex.getMessage() );
       }
       
       return buff.toString();
   }
 public static long[]  prepareReadCigarArray(String read, Cigar ncg, int read_start) throws Exception{
     //developed 10/2/2022
       long[] ref_pos_arr = new long[read.length()];
       long ref_pos = read_start -1; // a pointer to reference position
       int read_idx = -1; // a pointer to read position (0 based)
       long last_M_ref_pos = -1;
       for (CigarElement ce : ncg.getCigarElements()){
           if (ce.getOperator().equals(CigarOperator.M) || ce.getOperator().equals(CigarOperator.EQ) || ce.getOperator().equals(CigarOperator.X)){
               for (int i=0; i<ce.getLength();i++){
                   ref_pos_arr[++read_idx] = ++ref_pos; // last used  ref_pos
                   last_M_ref_pos = ref_pos;
               }
           }else if (ce.getOperator().equals(CigarOperator.D)){
               for (int i=0; i<ce.getLength();i++)
                   ref_pos++; // last used  ref_pos               
           }else if (ce.getOperator().equals(CigarOperator.I)){
               for (int i=0; i<ce.getLength();i++)
                   ref_pos_arr[++read_idx] = last_M_ref_pos; // last used  ref_pos

           }
       }
       //System.out.println(read_idx);
       if (read_idx <read.length()-1)
           throw new Exception("Exception inside getSNPAlleleFromReadSimple, Cigar string did not match with read length"); 
       return ref_pos_arr;
 }
 public static void setSNPAlleleFromRead(String read, long[]arr, List<SNP> snps) throws Exception{
       //developed 10/2/2022

       // remeber that read.length and arr.length must be equal
       // xIyD never happens in terms of aligner, but DI could happen, for example instead of 1M aligner might reposrt 1D1I, and only because of meeting penalty score threshold.
       //keep in mind that this function may return empty string as allele if snp_start:snp_end falls in a D element 
       
       int idx=0;
       for(SNP snp: snps){
           StringBuffer buff = new StringBuffer();
           int snp_start = snp.position;
           int snp_end = snp.position + snp.ref.length()-1;
           while(idx<arr.length && arr[idx]<snp_start)
               idx++;   
           //we get here when arr[idx]>= snp_start
            while(idx<arr.length && arr[idx]<=snp_end){
                buff.append(read.charAt(idx));
                idx++;
            }   
            
            snp.allele = buff.toString();
           
       }// for each SNP
  
}
  public static Cigar standardizeCigarSimple(Cigar cigar){
    Cigar ncg = new Cigar();
    for (CigarElement ce : cigar.getCigarElements()){
        if (!ce.getOperator().equals(CigarOperator.S)){ // if it is not soft trimmed
            ncg.add(ce);
        }

    }          
    return ncg;
  }
  
  public static Cigar standardizeCigar(Cigar cigar){
    Cigar ncg = new Cigar();
    for (CigarElement ce : cigar.getCigarElements()){
        if (!ce.getOperator().equals(CigarOperator.S)){ // if it is not soft trimmed
            if (ce.getOperator().equals(CigarOperator.EQ) || ce.getOperator().equals(CigarOperator.X))
                ncg.add(new CigarElement(ce.getLength(),CigarOperator.M));
            else 
                ncg.add(ce);
        }

    }          
    //now merge multiple M into one 
    Cigar nncg = new Cigar();
    int M_block_len = 0;
    for (CigarElement ce : ncg.getCigarElements()){
        if (!ce.getOperator().equals(CigarOperator.M)){
            //check if M block is already started
            if (M_block_len > 0){
                // here is end of M block, write it into 
                nncg.add(new CigarElement(M_block_len,CigarOperator.M));
                M_block_len = 0;
            } 
            nncg.add(ce);
        }else
            M_block_len += ce.getLength();

    }
    if (M_block_len > 0)
        nncg.add(new CigarElement(M_block_len,CigarOperator.M));   
    return nncg;
} 
public static Cigar getCigarFromString(String str){
    Cigar cig = new Cigar();
    Pattern p = Pattern.compile("\\d+?[MID]"); // reluctant
    Matcher   m = p.matcher(str); // remember that qseq contains Gaps or -, but we only look for continuous base without gap in the middle
    while(m.find()){
        // extract the equal part from reference (subject)
        String sub = str.substring(m.start(), m.end());
        if (sub.indexOf("I")>0)
            cig.add(new CigarElement(Integer.parseInt(sub.substring(0, sub.indexOf("I"))) , CigarOperator.I));
        else if (sub.indexOf("D")>0)
            cig.add(new CigarElement(Integer.parseInt(sub.substring(0, sub.indexOf("D"))) , CigarOperator.D));
        else if (sub.indexOf("M")>0)
            cig.add(new CigarElement(Integer.parseInt(sub.substring(0, sub.indexOf("M"))) , CigarOperator.M));
        //if (m.end()-m.start()>100)
        //    System.out.println(refseq.getBaseString());
    }       
    
    return cig;
}
   /*
   public static void setSNPAlleleFromReadVCF(String read, Cigar ncg, int read_start, SNP snp){
       //Written on 4/6/2019 for SNP filling from VCF file
       // Note : This method 2 times was verified and tested, do not change it !
       // for each snp in snpList, see if a snp of type snp.type exist at snp.position then we set snp.allele based on the read and its cigar values.
       // Note: no record will start with I operator in cigar (we tested it by SAMFileInfo.check_insertion_start() method)
       // This method will set snp.allele based on VCF file
       int ref_pos = read_start; // a pointer to reference position
       int read_idx = 0; // a pointer to read position (0 based)
       int cigar_idx = 0; // index for each Cigar Element
       //int list_idx = 0; // always points to the current element in the snpList (we dont loop based on elements of snpList)
       int cigar_size = ncg.numCigarElements();
       try{
       while (cigar_idx<cigar_size){
            CigarElement ce = ncg.getCigarElement(cigar_idx);
            for (int j=0; j<ce.getLength(); j++){
                if (ref_pos == snp.position){
                   // here depending on the type of SNP, we set the value of SNP and return
                   if (snp.type == Helper.SNP_MultiBase){
                       // see what we have in the read (a base or a deletion)
                       if (ce.getOperator().equals(CigarOperator.M)){
                           snp.allele = Helper.Base_String( Helper.nt_to_i(read.charAt(read_idx)) ) ;
                           return;
                       }                     
                   }else if (snp.type == Helper.SNP_MultiInsert){
                       // first we see if we are at the last element of current CigarElement
                       if (j==ce.getLength()-1){
                           if (cigar_idx < cigar_size-1){
                               CigarElement ce_next = ncg.getCigarElement(cigar_idx + 1);
                               if (ce_next.getOperator().equals(CigarOperator.I)){
                                   int insertion_length = ce_next.getLength();
                                   //we read current position and next insertion_length bases from the read
                                   snp.allele = read.substring(read_idx , read_idx + insertion_length + 1);// SNP was resolved
                                   return;
                               }else{
                                   // if next cigar element is not an Insertion then allele will be just reference allele
                                    snp.allele = snp.ref ;  
                                    return;                                                      
                               }
                            }else{
                               // we can get here when insertion is at the end of the read then allele will be unresolved, and will be resolved by another read
                                    // do nothing (allele will remain null)                                                   
                           }
                       }else{                           
                            snp.allele = snp.ref ;  // reference allele for MultiInsert SNP means absence of an Insertion in the read (SNP was resolved successfully).
                            return;                                                      
                       }
                        //here we get next cigar element
                   }// if multi insert
                   else if (snp.type == Helper.SNP_MultiDelete){
                       // first we see how far we are from end of current  CigarElement (including current element)
                       // so if ce.length = 4 and j=3 (means we are at the last element of current ce) then distance will be 1
                       //likewide if ce.length = 5 and j=3 (means we are before the last element of current ce) then distance will be 2
                       int distance_to_end = ce.getLength() -j ;
                       //Then we must have a del SNP with allele length equal to distance_to_end and also next CE element must be a D
                       int max_allele_length = snp.maxAlleleLength();
                       
                       if (distance_to_end <= max_allele_length){
                           if (cigar_idx < cigar_size-1){
                               CigarElement ce_next = ncg.getCigarElement(cigar_idx + 1);
                               if (ce_next.getOperator().equals(CigarOperator.D)){
                                  
                                   //From current reference position we read distance_to_end bases
                                   snp.allele = read.substring(read_idx , read_idx + distance_to_end);// SNP was resolved
                                   return;
                               }else{
                                   // if next cigar element is not a deletion then allele will be just reference allele
                                    snp.allele = read.substring(read_idx , read_idx + snp.ref.length()) ;  
                                    return;                                                      
                               }
                            }else{
                               // we can get here when deletion is at the end of the read then allele will be unresolved, and will be resolved by another read
                                    // do nothing (allele will remain null)                                                   
                           }
                       }else{                                     
                            snp.allele = read.substring(read_idx , read_idx + snp.ref.length()); 
                            return;                                                      
                       }                     
                       
                   }//if multi delete
               } // if ref_pos == snp.position
                // consume a base from the read if cigar element is not D
                if (!ce.getOperator().equals(CigarOperator.D))
                  read_idx++;
                // consume a base from the reference if cigar element is not I
                if (!ce.getOperator().equals(CigarOperator.I))
                  ref_pos++;
            }// j for each base in a Cigar Element
            cigar_idx++;
       } //while Cigar element
       
   
       }catch(Exception ex){
                   System.out.println("exception inside setSNPAlleleFromRead" + ex.getMessage() + "  snp:" + snp.toString());
       }
   }*/
   
    public static List<Allele> getDiploidGenotype(String ref, double minor_ratio, int min_number, Map<String, Integer> bases){
        //Revised on 4/5/2021, has to have minimum number of reads, regardless been reference allele or not
        //This method was written in 25/6/2019, it determines zygosity and genotype of a diploid position
        //parameter explain: min_number  we need to have at least min_number of those regardless if it is reference allele or not
        //minor_ratio: is the ratio of minor allele to distinguish zygosity:if minor_ratio is 0.2 then if we have 0.85 to 0.15 ratios then we consider it homo and if it is 0.22 and 0.78 then it is hetero
        //Note: if you change logic here you need change logic in PileupTriplet as well.
        List<Allele> alleles = new ArrayList();
        String top_gt = null;
        String second_top_gt = null;
        double top = 0;
        double second_top = 0;
        int top_count = 0;
        int second_top_count = 0;
        
        int total_count = new MapUtil().totalCount(bases);

        for(String base : bases.keySet()){
            int cnt = bases.get(base);
            double prob = (double)cnt/total_count;
            if (/*!base.equals(ref) &&*/ cnt<min_number)
                continue;
            if ((prob>top) || (prob>0 && base.equals(ref) && prob==top )){
               // shift whatever is in top into second top
                second_top = top;
                second_top_count = top_count;
                second_top_gt = top_gt;

                top = prob;
                top_count = cnt;
                top_gt = base;

            }else if ((prob>second_top)|| (prob>0 && base.equals(ref) && prob==second_top)){
                second_top = prob;
                second_top_count = cnt;
                second_top_gt = base;                    
            }

        }            
        double major_ratio = 1-minor_ratio;

        boolean isHet = false;

        if (top>0 && top <= major_ratio && second_top>=minor_ratio){
            isHet = true; // see if it is Heterozygous maximum difference is 70/30, if for example if it is 80/20 we consider it homozygous
        }else if (top>0 && second_top < minor_ratio){
           isHet = false; 
        }

        if(!isHet){
            if (top_gt != null)
               alleles.add(new Allele(top_gt, top_count));
        }else if (isHet){
            if (top_gt != null)
                alleles.add(new Allele(top_gt, top_count));
            if (second_top_gt != null)
                alleles.add(new Allele(second_top_gt, second_top_count));
        }

        return alleles;   
    }

   public static int getDepth(List<Allele> alleles, String allele){
       for(Allele al : alleles)
           if(al.value.toUpperCase().equals(allele.toUpperCase()))
               return al.depth;
       return 0;
   }
   
   private static Cigar getPartialCigar(SAMRecord rec, int part_start_pos, int part_end_pos){
       // return part of cigar that starts from ref_start_pos and end with ref_end_pos (both inclusive)
       Cigar ncg = new Cigar();
       int ref_pos = rec.getAlignmentStart(); // a pointer to reference position       
       int cigar_idx = 0; // index for each Cigar Element
       int cigar_size = rec.getCigarLength();
       try{
           //System.out.println("Start pos " + ref_start_pos + " end pos " + ref_end_pos + " Cigar is " + cg.toString() );
       while (cigar_idx<cigar_size){
            CigarElement ce = rec.getCigar().getCigarElement(cigar_idx);
            for (int j=0; j<ce.getLength(); j++){
                if (ref_pos >= part_start_pos && ref_pos<=part_end_pos){
                        ncg.add(new CigarElement(1, ce.getOperator()));   
                        //System.out.print(ce.getOperator().toString());
                } // if ref_pos == ref_start_pos
                // consume a base from the reference if cigar element is not I
                if (!ce.getOperator().equals(CigarOperator.I))
                  ref_pos++;
            }// j for each base in a Cigar Element
            cigar_idx++;
       } //while Cigar element
       
       return ncg;
   
       }catch(Exception ex){
                   System.out.println("exception inside getPartialCigar" + ex.getMessage());
       }

       return ncg;       
   }
   private static Cigar trimCigarLeft(Cigar cg, int len){
       // from left side trims SAM read with len and returns new Cigar, is called inside the mergePairedReads
       
       Cigar ncg = new Cigar();
       int cigar_idx = 0; // index for each Cigar Element
       int read_idx = 0; // a pointer to read position (0 based)
       int j=0; // pointer in current cigar element
       int cigar_size = cg.numCigarElements(); 
       boolean reached = false;
       while (cigar_idx<cigar_size){
            CigarElement ce = cg.getCigarElement(cigar_idx);
            
            
            for (j=1; j<=ce.getLength(); j++){
               if (!ce.getOperator().equals(CigarOperator.D))
                  read_idx++;
               if (read_idx == len){
                   reached = true;
                   break; // break from the inner loop
               }
            }//for
            if (reached)
                break; // break from outer loop
            
            cigar_idx++;
       }//while
       // when we get here we have definitely have reached
       if(reached){
           //cigar_idx points to the current cigar element and j points to the current Cigar Element position
           CigarElement ce = cg.getCigarElement(cigar_idx);
           // Generally a Cigar should not start or end with D (it does not make sense), but here we need D because this cigar will attach to anothe cigar
           if (j<ce.getLength()){
                  ncg.add(new CigarElement(ce.getLength()-j, ce.getOperator()));
           }
           // then add remaining cigar elements intio ncg without any change
           if (cigar_idx<cigar_size-1)
                for (int i=cigar_idx+1; i< cigar_size; i++)
                    ncg.add(cg.getCigarElement(i));
       }
       
       return ncg;
   }
   public static SAMRecord mergePairedReads(SAMRecord pair1, SAMRecord pair2) throws Exception{
       //2 paired reads Cigars get joined to make a single Insert fragment, it is assumed that pair1 start alignment is before pair2 start alignment (if bam file is sorted this condition is gauranteed)
       if (pair2.getAlignmentEnd() <= pair1.getAlignmentEnd())
           return pair1;
       if (pair1.getAlignmentStart()==pair2.getAlignmentStart())
           return pair2; // this means pair one is inside pair2 so we just return pair2
       // so we get here when pair1.getAlignmentStart()<pair2.getAlignmentStart() and pair2.getAlignmentEnd() > pair1.getAlignmentEnd()
       SAMRecord frag = pair1.deepCopy(); // you could also use method clone()
       int insert_ftag_start_idx = pair1.getAlignmentStart();
       //int insert_frag_end_idx = pair2.getAlignmentEnd();
       StringBuffer frag_str = new StringBuffer();      
       frag_str.append(pair1.getReadString());
       
       Cigar ncg = new Cigar();
       // we add cigar element of read 1 except last one
       for (int i=0; i<pair1.getCigarLength()-1; i++)
           ncg.add(pair1.getCigar().getCigarElement(i));
                  
       
       //now work out overlap size between 2 reads
       int overlap = 0; 
      
       //System.out.println("pair 1 alinment start:"+ pair1.getAlignmentStart()  + "pair 1 alinment end:"+ pair1.getAlignmentEnd()+ " Pair 2 alignment start:" + pair2.getAlignmentStart() + " pair1 read length "+ pair1.getReadLength() );
       if (pair1.getAlignmentEnd()>= pair2.getAlignmentStart()){
           /*
           reads look like this ------>
                                   <------- 
           */
           //first we grab part of pair1 cigar in the window starting from pair2.getAlignmentStart() upto pair1.getAlignmentEnd()
           Cigar partial_cg = Helper.getPartialCigar(pair1, pair2.getAlignmentStart(), pair1.getAlignmentEnd());
           // then we count total read length associated to the partial_cg, that will be overlap between 2 reads
           
           for (CigarElement ce : partial_cg.getCigarElements()){
               
               if (!ce.getOperator().equals(CigarOperator.D))
                   overlap++;
           }
           //System.out.println("Overlap is " + overlap);
           // we cut overlap length from the begining of pair2 and add it to string buffer
           frag_str.append(pair2.getReadString().substring(overlap));
           // now work out new Cigar
           Cigar right_part = trimCigarLeft(pair2.getCigar(), overlap);
           CigarElement last_ce = pair1.getCigar().getLastCigarElement();
           CigarElement first_ce = right_part.getFirstCigarElement();
           if (last_ce.getOperator().equals(first_ce.getOperator()))
               ncg.add(new CigarElement(last_ce.getLength()+first_ce.getLength(), last_ce.getOperator()));
           else{
               ncg.add(last_ce);
               ncg.add(first_ce);
           }
           // now add remaining of the right part of pair2
           for (int i=1; i<right_part.numCigarElements(); i++)
                   ncg.add(right_part.getCigarElement(i));
           
       }else if (pair1.getAlignmentEnd()< pair2.getAlignmentStart()){
           //reads look like this ------->     <-----
           overlap = pair2.getAlignmentStart() - pair1.getAlignmentEnd() -1 ; 
           if (overlap>0){
               for (int i=0; i<overlap; i++)
                   frag_str.append('N');
               frag_str.append(pair2.getReadString());
               // now work out new Cigar
               CigarElement last_ce = pair1.getCigar().getLastCigarElement();
               CigarElement first_ce = pair2.getCigar().getFirstCigarElement();

               if (last_ce.getOperator().equals(CigarOperator.M) && first_ce.getOperator().equals(CigarOperator.M)){
                   ncg.add(new CigarElement(last_ce.getLength()+overlap+first_ce.getLength() ,CigarOperator.M));
               }else if (last_ce.getOperator().equals(CigarOperator.M)){
                   ncg.add(new CigarElement(last_ce.getLength()+overlap ,CigarOperator.M));
                   ncg.add(first_ce);
               }else if (first_ce.getOperator().equals(CigarOperator.M)){
                   ncg.add(last_ce);
                   ncg.add(new CigarElement(overlap+first_ce.getLength() ,CigarOperator.M));
               }else{
                   // when none of the operators is not M then
                   ncg.add(last_ce);
                   ncg.add(new CigarElement(overlap ,CigarOperator.M));
                   ncg.add(first_ce);

               }
              // now add rest of pair2 cigar elements
              for (int i=1; i<pair2.getCigarLength(); i++)
                       ncg.add(pair2.getCigar().getCigarElement(i));

           }else if (overlap==0){
           //reads look like this -----><----------
           frag_str.append(pair2.getReadString());
          
           CigarElement last_ce = pair1.getCigar().getLastCigarElement();
           CigarElement first_ce = pair2.getCigar().getFirstCigarElement();
           if (!last_ce.getOperator().equals(first_ce.getOperator())){
               ncg.add(last_ce);
               ncg.add(first_ce);
           }else
               ncg.add(new CigarElement(last_ce.getLength()+first_ce.getLength(), last_ce.getOperator()));
               // now add rest of pair2 cigar elements
               for (int i=1; i<pair2.getCigarLength(); i++)
                    ncg.add(pair2.getCigar().getCigarElement(i));               
           }
       }
       
       frag.setCigar(ncg);
       //byte[] arr = frag_str.toString().getBytes("UTF-8");
       frag.setReadBases(frag_str.toString().getBytes("UTF-8"));
      // frag.setReadString(frag_str.toString());
       frag.setAlignmentStart(insert_ftag_start_idx);
       return frag;
   }
   
   public static boolean validateSAMRecord(SAMRecord rec){
       // see if sum of cigar elements length is equal to read size
       int sum_ce = 0;
       for(CigarElement ce : rec.getCigar().getCigarElements())
           if (!ce.getOperator().equals(CigarOperator.D))
              sum_ce += ce.getLength();
       
       if (sum_ce == rec.getReadLength())
           return true;
       
       return false;
   }
   
   public static String getGenotypeVariation(byte[] bases){
       
       Set<Byte> gset = new HashSet();
       for (byte b :bases)
           gset.add(b); // in this way we get rid of duplicate bases
       
       StringBuffer buff = new StringBuffer();
       for (byte b : gset)
           if (b != UNKNOWN)
              buff.append(Base_String(b)+"/");
       
       return buff.substring(0, buff.length()-1);
   }
   
   public static boolean isVariation(String[] genotypes){
       //It checks if there is variation in genotypes of samples
       Set<String> gset = new HashSet();
       for (String g :genotypes)
           if(g!=null)
              gset.add(g); // in this way we get rid of duplicate bases
       
       if (gset.size()>1)
           return true;
       
        return false;  
   }
   
   public static boolean isVariationFromReference(String[] genotypes, String ref){
       //It checks if there is any variation from the reference then return true
       //current logic, if a genotype is empty string it is treated as difference from reference (or deletion); however it may be incporrect;we may change this logic in future.
       for (String g :genotypes)
           if(g!=null && !g.isEmpty() &&  !g.toLowerCase().equals(ref.toLowerCase()))
              return true; 
       
       
        return false;  
   }
   
   public static String isRepeatable(String read){
       //It checks if this read is repeated, like ACCACCACCACCACC..., if it is true, it will return the repeated substring otherwise it returns null
       //Note: If a substring of size n is repeated in the read, then we will have n different repeats , for example in the above example we can say ACC, CCA, and CAC have been repeated.
       //This method returns the first repeat, that is ACC
      
       TreeMap<String, Integer> candidates = new TreeMap(); // a map of candidate repeat to its current index (0-based)
       for(int i=0; i<read.length(); i++){
           //see if character at position i can verify or reject any repeat candidates in the map, if rejects then that entry will be removed from the Map
            Iterator<Entry<String, Integer>> it = candidates.entrySet().iterator();
            while (it.hasNext()){
                Entry<String, Integer> entry = it.next();
                //look for the next character in the candidate
                int idx = entry.getValue();
                String candidate = entry.getKey();
                if (idx == candidate.length()-1)
                    idx = 0; // circular 
                else
                    idx++;
                char next_char = candidate.charAt(idx);
                    
                if (next_char != read.charAt(i))
                    it.remove();
                else //advance the pointer index for this candidate
                    candidates.replace(candidate, idx);
            }
            //at the end we add this prefix-sequence to the Map as candidate repeat, with default pointer index pointing at the end of 
            candidates.put(read.substring(0, i+1), i);
       }//for each character in the read
       
       
       // now see if anything left in the map, we need to report the smallest one, so we report the first element.
       //For example: in the above example ACC, ACCACC, ACCACCACC all are repeated, but we report the first one
       if (!candidates.isEmpty())
          return candidates.firstKey();
       
       //This method never returns null, at worst case it returns the whole sequence as repeat
       return null;
       
   }
   
   public static String getWheatChromosome(String contig_id){
        String chromosome = "U"; // unknown
        contig_id = contig_id.toLowerCase();
        if(contig_id.indexOf("1a")>=0)
            chromosome = "1A";
        else if (contig_id.indexOf("1b")>=0)
            chromosome = "1B";
        else if (contig_id.indexOf("1d")>=0)
            chromosome = "1D";        
        else if (contig_id.indexOf("2a")>=0)
            chromosome = "2A";
        else if (contig_id.indexOf("2b")>=0)
            chromosome = "2B";
        else if (contig_id.indexOf("2d")>=0)
            chromosome = "2D";
        else if (contig_id.indexOf("3a")>=0)
            chromosome = "3A";
        else if (contig_id.indexOf("3b")>=0)
            chromosome = "3B";
        else if (contig_id.indexOf("3d")>=0)
            chromosome = "3D";        
        else if (contig_id.indexOf("4a")>=0)
            chromosome = "4A";
        else if (contig_id.indexOf("4b")>=0)
            chromosome = "4B";
        else if (contig_id.indexOf("4d")>=0)
            chromosome = "4D";
        else if (contig_id.indexOf("5a")>=0)
            chromosome = "5A";
        else if (contig_id.indexOf("5b")>=0)
            chromosome = "5B";
        else if (contig_id.indexOf("5d")>=0)
            chromosome = "5D";
        else if (contig_id.indexOf("6a")>=0)
            chromosome = "6A";
        else if (contig_id.indexOf("6b")>=0)
            chromosome = "6B";
        else if (contig_id.indexOf("6d")>=0)
            chromosome = "6D";
        else if (contig_id.indexOf("7a")>=0)
            chromosome = "7A";
        else if (contig_id.indexOf("7b")>=0)
            chromosome = "7B";
        else if (contig_id.indexOf("7d")>=0)
            chromosome = "7D";       
       return chromosome;
   }
   public static boolean areInTheSameWheatChromosome(String chr1, String chr2){
       if (chr1.startsWith(chr2) || chr2.startsWith(chr1))
           return true;
       // because of translocation some parts of 4B translocated with 7A
       if ((chr1.startsWith("7A")&&chr2.startsWith("4B")) || (chr1.startsWith("4B")&&chr2.startsWith("7A")))
           return true;
       
       return false;
   }
   
   
     public static int readIdxToArrayIdx(int read_idx, int read_num){
        // we need to know if we get read idx (1-based) where in the array it should be stored
        if (read_idx>0)
            return read_idx-1;
        //so we get here when read_idx is negative
        return (read_num + Math.abs(read_idx)-1);
    }
     public static int arrayIdxToReadIdx(int arr_idx, int read_num){
        //we have array_idx, we want to know which read idx (1-based) is in the fastq file
        if (arr_idx<read_num)
            return arr_idx+1;
        //we get here if arr_idx>=read_num
        return (-(arr_idx-read_num+1));
    }
    public static int readIdxToProcessedArrayIdx(int read_idx){
        // we need to know if we get read idx (1-based) whether it has been processed or not
        
        return (Math.abs(read_idx)-1);
    }

}
