/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
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


public class Helper {
  
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