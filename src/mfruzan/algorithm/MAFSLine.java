/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mfruzan.algorithm;

import java.util.List;

/**
 *
 * @author a1195806
 */
public class MAFSLine {
    //encapsulates one line of MAF S Line
    public String contig = null;
    public int startPOS = 0;  // this is 1-based start of fragment in contig (BlastHSP like)
    public int fragment_size = 0;
    public int contigSize = 0;
    public boolean pos_strand = true; // default is positive orientation 
    public String alignment_block = null ;  // this is block of sequence including - 
    //this is caculated field
    public int endPOS = 0;  // this is 1-based end of fragment in contig (BlastHSP like)
 //calculated field (is used in MAF file format)
    public byte orientation = 1; // 1:same (+/+ or -/-) 0:reverse (+/- or -/+)    
    public String sLine = null;
    
    public MAFSLine(){
       // this is default constructor      
    }
    public MAFSLine(String c, int zero_pos,  int fs,boolean ps,  int cs,  String ab, boolean dummy){
        // this constructor reads a MAFSLine in MAF standard, then convert it to BlastHSP standard (update startPOS and endPOS, startPOS always less than endPOS)
        this.contig = c;
        this.fragment_size = fs;
        this.pos_strand = ps;
        this.contigSize = cs;
        

        
        this.alignment_block = ab;
        // work out startPOS and endPOS
       if (this.pos_strand){
          this.startPOS = zero_pos + 1;   // convert zero based into 1 based
          this.endPOS = this.startPOS + this.fragment_size - 1; // 1-based end of fragment
       }else{
          this.endPOS = this.contigSize -  zero_pos; // 1-based end of fragment
          this.startPOS = this.endPOS - this.fragment_size + 1;  //1-based start of fragment
       }    
        
    }
  /* public MAFSLine(String c, int spos, int epos,boolean ps, byte or, int cs,  String ab){
        // this constructor starts a MAFSLine in BLASTHSP standard
        this.contig = c;
        
        this.pos_strand = ps;
        this.orientation = or;
        this.contigSize = cs;
        this.startPOS = spos;
        this.endPOS = epos;
        this.alignment_block = ab;
    }    */
    public MAFSLine(String c, int spos, int epos,boolean ps, int cs,  String ab){
        // this constructor reads a MAFSLine in BLASTHSP standard
        this.contig = c;
        
        this.pos_strand = ps;
        
        this.contigSize = cs;
        this.startPOS = spos; // 1-bases inclusive
        this.endPOS = epos; //1-based inclusive
        //fragment length can be calculated as endPOS-startPOS + 1
        this.alignment_block = ab;
    }    
   
    public int overlap(MAFSLine line){
        // it allows for overlap bp overlap between regions
        if (!line.contig.equals(this.contig))
            return 0;
        if (line.startPOS > this.endPOS  || line.endPOS<this.startPOS)
            return 0;
        //otherwise there is overlap between them 
        //see if one is covered the other one
        if (line.startPOS>= this.startPOS && line.endPOS<=this.endPOS)
            return (line.endPOS - line.startPOS + 1);
        if (this.startPOS>= line.startPOS && this.endPOS<=line.endPOS)
            return (this.endPOS - this.startPOS + 1);
        //we get here when there is partial ovelap
        return Math.min(Math.abs(line.startPOS - this.endPOS), Math.abs(line.endPOS - this.startPOS));
        
    }
    
    public int addToListNoOverlap(List<MAFSLine> list, int max_overlap){
        // tries to add this to the list if there is no major overlapping with elements of the list and return 0
        //otherwise it returns number of overlaps
        for(MAFSLine hsp: list){
            if(overlap(hsp)> max_overlap)
                return overlap(hsp);            
        }
        // we get here when there is no major overlapping 
        list.add(this);
        return 0;
    }    
    
}
