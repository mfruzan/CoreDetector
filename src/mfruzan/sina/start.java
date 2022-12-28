/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package mfruzan.sina;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.security.CodeSource;
import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.ForkJoinPool;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import mfruzan.common.MAFProcessor;

/**
 *
 * @author a1195806
 */
public class start {
       public static void main(String args[]) throws ParseException {
        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
         */
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException ex) {
            java.util.logging.Logger.getLogger(Start.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(Start.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(Start.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(Start.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>
        //We want our code runs by eligible person
        Set<String> allowed_addresses = new HashSet<>(Arrays.asList("00-0E-1E-9C-D1-B0", "00-0E-1E-9D-07-C0", "5C-F3-FC-DB-AB-1C", "90-E2-BA-0A-D0-F4", "00-23-24-3D-7E-FF"));
        //if (!allowed_addresses.contains(getFirstMAC()))
          //  System.exit(0);
        //We want to be sure that our code comes from SINA.jar
        CodeSource src = start.class.getProtectionDomain().getCodeSource();
        String jar = src.getLocation().getFile();
        System.out.print(jar);
        if (!jar.endsWith("/classes/") && !(jar.endsWith("MFbio.jar")|| jar.endsWith("SINA.jar") || jar.endsWith("WHEATbio.jar")))
            System.exit(0);
        DateFormat format = new SimpleDateFormat("yyyy-MM-dd");
       // if (new Date().after(format.parse("2020-09-30")))
        //    System.exit(0); 
        /* Create and display the form */
        
        boolean showForm = false;
        String srcDir = null;
        String destDir = null;
        String task = null;
        String host = null; // db host
        String exporterClass =null;
        String p1 =null;  // this is used for a general purpose parameter
        String p2 =null;  // this is used for a general purpose parameter
        String p3 =null;  // this is used for a general purpose parameter
        String p4 =null;  // this is used for a general purpose parameter
        String p5 =null;  // this is used for a general purpose parameter
        String p6 =null;  // this is used for a general purpose parameter
        String p7 =null;  // this is used for a general purpose parameter
        String p8 =null;  // this is used for a general purpose parameter
 
        
        
        String mate1File = null;
        String mate2File = null;
        String unpairedFile = null;
        
        String file1 = null;
        String file2 = null;
        String file3 = null;
        String file4 = null;
        String file5 = null;

        String bam = null;
        String bamx = null;

        String ref = null;
        String refx = null;
        String vcf = null;
        String out = null;
       String minmapq =null;  
        String maxmapq =null;  
        String seqtype =null;  
        String afl = null; //average fragment length 
        String threads =null;  // this is used for number of threads in thread pool        
        
        
        // developed on 30/9/2020 to be able to read parameters like --showform no rather than showform=no
        Map<String, String> prms = new HashMap(); // will keep entries like showform -> no
        for(int i =0; i<args.length; i++){
            //even args should start with --
            if (i%2==0){
                if (!args[i].startsWith("--")){
                    System.out.println("command line argument should starts with --");  
                    System.exit(1);  
                }
                // so we get here when argument has reached 
                String value = args[i+1];
                if (value.startsWith("--")){
                    System.out.println("argument value is missing "+ args[i]);  
                    System.exit(1);                      
                }
                prms.put(args[i].substring(2), value);
            }
        }
        
        
        /*
        
        Properties prop = new Properties();
        StringBuffer sb = new StringBuffer();
        for(int i =0; i<args.length; i++){
            sb.append(args[i]+"\n");
        }
        try {
            prop.load(new StringReader(sb.toString()));
        } catch (IOException ex) {
             System.out.println(ex.getMessage());
        }
       */
       //--------------------------------------------
        try{//for debugging purpose set below variable to True
            showForm = false;
            if (prms.get("showform")!= null && prms.get("showform").equals("yes"))
                    showForm = true;
           task = prms.get("task");
           host = prms.get("host");
           srcDir = prms.get("srcdir");
           destDir = prms.get("destdir");
           exporterClass = prms.get("class");
           p1 = prms.get("p1");
           p2 = prms.get("p2");
           p3 = prms.get("p3");
           p4 = prms.get("p4");
           p5 = prms.get("p5");
           p6 = prms.get("p6");
           p7 = prms.get("p7");
           p8 = prms.get("p8");
           threads = prms.get("threads");
           
           mate1File = prms.get("1");
           mate2File = prms.get("2");
           unpairedFile = prms.get("U");
           file1 = prms.get("file1");
           file2 = prms.get("file2");
           file3 = prms.get("file3");
           file4 = prms.get("file4");
           file5 = prms.get("file5");
            
           ref = prms.get("ref");
           refx = prms.get("refx");
           bam = prms.get("bam");
           bamx = prms.get("bamx");
           vcf = prms.get("vcf");
           out = prms.get("out");
           minmapq =prms.get("minmapq");  
           maxmapq =prms.get("maxmapq");  
           seqtype =prms.get("seqtype");  
           afl = prms.get("afl"); //averag           
            
            /*
           if(prop.getProperty("showform").equals("no")){
             showForm = false;
           }
           //-------------------------------------------
           task = prop.getProperty("task");
           host = prop.getProperty("host");
           srcDir = prop.getProperty("srcdir");
           destDir = prop.getProperty("destdir");
           exporterClass = prop.getProperty("class");
           p1 = prop.getProperty("p1");
           p2 = prop.getProperty("p2");
           p3 = prop.getProperty("p3");
           p4 = prop.getProperty("p4");
           p5 = prop.getProperty("p5");
           p6 = prop.getProperty("p6");
           p7 = prop.getProperty("p7");
           p8 = prop.getProperty("p8");
           threads = prop.getProperty("threads");
           
           mate1File = prop.getProperty("1");
           mate2File = prop.getProperty("2");
           unpairedFile = prop.getProperty("U");
           file1 = prop.getProperty("file1");
           file2 = prop.getProperty("file2");
           file3 = prop.getProperty("file3");
           file4 = prop.getProperty("file4");
           file5 = prop.getProperty("file5");
           */

        }catch(NullPointerException ex){}

        if(showForm){
            java.awt.EventQueue.invokeLater(new Runnable() {
                public void run() {
                    new start().setVisible(true);
                }
            });
        }else{
            // here we examine parameters
           // if (host != null)
               // DBUtils.host = host.trim();
           // else
              //  DBUtils.host = "localhost";


             if (task.equals("maf2fasta")){
                // This task reads maf file and writes it as fasta file (for each HSP selects smaller sequence)
                //src dir : maf file
                //dest: output fasta file
                if(srcDir == null  || destDir == null  ){
                   System.out.println("'All following parameters are required: 'srcdir', 'destdir''");  
                   System.exit(1);                     
                }  
                System.out.println("Processing File: " + srcDir );
                
                MAFProcessor maf = new MAFProcessor(srcDir); 
                maf.load();
                maf.writeAsFasta(destDir);
                System.out.println("Done.");

            }else if (task.equals("mafinfo")){
                // This task reads binary maf file and prepare stats like fastainfo
                //src dir : maf file
                if(srcDir == null ){
                   System.out.println("'All following parameters are required: 'srcdir''");  
                   System.exit(1);                     
                }  
                System.out.println("Processing File: " + srcDir );
                
                MAFProcessor maf = new MAFProcessor(srcDir); 
                maf.info_multi();
                System.out.println("Done.");

            }else if (task.equals("mafcheck")){
                // This task checks cordinates of maf lines and compares with the sequence
                //src dir : maf file
                //file1: a file containing list of fasta files and their fasta indexes
                if(srcDir == null || file1 == null ){
                   System.out.println("All following parameters are required: 'srcdir', 'file1'");  
                   System.exit(1);                     
                }  
                System.out.println("Processing File: " + srcDir );
                
                MAFProcessor maf = new MAFProcessor(srcDir); 
                maf.checkCordinates(file1);
                System.out.println("Done.");

            }else if (task.equals("maf2uniquequery")){
                // This task reads maf file and writes non-ovelapping query sequences as fasta file also adds start and end positions to the contig names separated by !  
                //src dir : maf file
                //dest: output fasta file
                //file1: output Filtered MAF File
                //p1(optional): min length of query 
                if(srcDir == null  || destDir == null || file1 == null ){
                   System.out.println("'All following parameters are required: 'srcdir', 'destdir', 'file1'");  
                   System.exit(1);                     
                }  
                System.out.println("Processing File: " + srcDir );
                int min_len = 0;
                if (p1!=null)
                    min_len = Integer.parseInt(p1);
                MAFProcessor maf = new MAFProcessor(srcDir); 
                maf.extractNonOverlapQuery(destDir, file1, min_len);
                System.out.println("Done.");

            }else if (task.equals("maf2msa")){
                // wrriten 6/1/2020
                // This task reads maf file and writes it as fasta file (selects query)
                //src dir : directory of binary maf files
                //p1: name of binary maf files , separated by comma (last to first)
                //dest: output msa fasta file
                //file1 : output maf file
                //file2: genome file , alias names : its a tab separated file, Alias_Name  File_Name  (First line is Query the following lines are subjects in order of the alignment)
                
                if(srcDir == null  || destDir == null || file1 == null || file2==null || p1 == null){
                   System.out.println("'All following parameters are required: 'srcdir', 'destdir', 'file1', 'file2', 'p1'");  
                   System.exit(1);                     
                } 
                new MAFProcessor().generateMSA(srcDir, p1, destDir, file1, file2);
                System.out.println("Done.");

            }else if (task.equals("mafuniquequery")){
                // wrriten 6/1/2020
                // This task reads maf file and selects non-ovelap query regions
                //src dir : input maf file
                //dest: output maf file
               
                if(srcDir == null  || destDir == null  ){
                   System.out.println("'All following parameters are required: 'srcdir', 'destdir'");  
                   System.exit(1);                     
                } 
                System.out.println("Processing File: " + srcDir );
                
                MAFProcessor maf = new MAFProcessor(srcDir); 
                maf.load();
                maf.uniqueQuery(destDir);
                System.out.println("Done.");

            }            else{
               System.out.println("Error:Task '" + task + "' Unrecognized."); 
            }
        }

    
}
