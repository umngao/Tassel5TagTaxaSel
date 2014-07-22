/*
* To change this template, choose Tools | Templates
* and open the template in the editor.
*/
package net.maizegenetics.dna.map;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import net.maizegenetics.dna.tag.AbstractTags;
import net.maizegenetics.dna.tag.TagsByTaxa.FilePacking;
import net.maizegenetics.dna.BaseEncoder;

/**
 * Class that hold genetic position of tags, genetic positions are the position of markers associated with tags
 * @author Fei Lu
 */
public class TagsOnGeneticMap extends AbstractTags {
    protected int[] gChr;
    protected int[] gPos;
    protected byte[] ifPAV;
    protected float[] prediction;
    
    /**
     * Construct TOGM from a file
     * @param infileS
     *        File name of input TagsOnGeneticMap file
     * @param format
     *        FilePacking format
     */
    public TagsOnGeneticMap (String infileS, FilePacking format) {
        this.readDistFile(infileS, format);
    }
    
    /**
     * Return chromosome of genetic position
     * @param index
     * @return Chromosome of genetic position
     */
    public int getGChr (int index) {
        return gChr[index];
    }
    
    /**
     * Return site of genetic position
     * @param index
     * @return Genetic position
     */
    public int getGPos (int index) {
        return gPos[index];
    }
    
    /**
     * Return if this tag is PAV
     * @param index
     * @return 
     */
    public byte getIfPAV (int index) {
        return ifPAV[index];
    }
    
    /**
     * Return prediction value from model
     * @param index
     * @return 
     */
    public float getPrediction (int index) {
        return prediction[index];
    }
    
    /**
     * Read tagsOnGeneticMap file
     * @param infileS
     * @param format
     */
    public void readDistFile (String infileS, FilePacking format) {
        System.out.println("Reading TOGM file from " + infileS);
        File infile = new File (infileS);
        switch (format) {
            case Text:
                readTextTOGMFile(infile);
                break;
            default:
                readBinaryTOGMFile(infile);
                break;
        }
        System.out.println("TOGM file read. Tatol: " + this.getTagCount() + " Tags");
    }
    
    /**
     * Read text TOGM file
     * @param infile
     */
    private void readTextTOGMFile (File infile) {
        try {
            BufferedReader br = new BufferedReader (new FileReader(infile), 65536);
            br.readLine();
            tagLengthInLong = br.readLine().split("\t")[0].length()/BaseEncoder.chunkSize;
            int tagNum = 1;
            while ((br.readLine())!=null) {
                tagNum++;
            }
            this.iniMatrix(tagLengthInLong, tagNum);
            br = new BufferedReader (new FileReader(infile), 65536);
            br.readLine();
            for (int i = 0; i < tagNum; i++) {
                String[] temp = br.readLine().split("\t");
                long[] t = BaseEncoder.getLongArrayFromSeq(temp[0]);
                for (int j = 0; j < tagLengthInLong; j++) {
                    tags[j][i] = t[j];
                }
                tagLength[i] = Byte.parseByte(temp[1]);
                gChr[i] = Integer.parseInt(temp[2]);
                gPos[i] = Integer.parseInt(temp[3]);
                ifPAV[i] = Byte.parseByte(temp[4]);
                prediction[i] = Float.parseFloat(temp[5]);
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Read binary TOGM file
     * @param infile
     */
    private void readBinaryTOGMFile (File infile) {
        try {
            DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(infile), 65536));
            tagLengthInLong = dis.readInt();
            int tagNum = (int)(infile.length()/(8*tagLengthInLong)+1+4+4+1+4);
            this.iniMatrix(tagLengthInLong, tagNum);
            for (int i = 0; i < this.getTagCount(); i++) {
                for (int j = 0; j < this.tagLengthInLong; j++) {
                    this.tags[j][i] = dis.readLong();
                }
                this.tagLength[i] = dis.readByte();
                this.gChr[i] = dis.readInt();
                this.gPos[i] = dis.readInt();
                this.ifPAV[i] = dis.readByte();
                this.prediction[i] = dis.readFloat();
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Initialize the matrix of TOGM
     * @param tagLengthInLong
     *        Tag length in Long primitive data type
     * @param tagNum
     *        Total tag number
     */
    protected void iniMatrix (int tagLengthInLong, int tagNum) {
        tags = new long[tagLengthInLong][tagNum];
        tagLength = new byte[tagNum];
        gChr = new int[tagNum];
        gPos = new int[tagNum];
        ifPAV = new byte[tagNum];
        prediction = new float[tagNum];
    }
    
    /**
     * Write TagsOnGeneticMap file
     * @param outfileS
     *        File name of output file
     * @param format
     *        FilePacking format
     */
    public void writeDistFile (String outfileS, FilePacking format) {
        System.out.println("Writing TOGM file to " + outfileS);
        switch (format) {
            case Text:
                writeTextTOGMFile(outfileS);
                break;
            default:
                writeBinaryTOGMFile(outfileS);
                break;
        }
        System.out.println("TOGM file written");
    }
    
    /**
     * Write text TOGM file
     * @param outfileS
     */
    private void writeTextTOGMFile (String outfileS) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
            bw.write("Tag\tTagLength\tGChr\tGPos\tIfPAV\tPredictedDistance");
            bw.newLine();
            long[] temp = new long[this.tagLengthInLong];
            for (int i = 0; i < this.getTagCount(); i++) {
                for (int j = 0; j < temp.length; j++) {
                    temp[j] = tags[j][i];
                }
                bw.write(BaseEncoder.getSequenceFromLong(temp)+"\t"+String.valueOf(this.getTagLength(i))+"\t");
                bw.write(String.valueOf(this.gChr[i])+"\t"+String.valueOf(this.gPos[i])+"\t"+String.valueOf(ifPAV[i])+"\t"+String.valueOf(this.prediction[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Write binary TOGM file
     * @param outfileS
     */
    private void writeBinaryTOGMFile (String outfileS) {
        try {
            DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS), 65536));
            dos.writeInt(tagLengthInLong);
            for (int i = 0; i < this.getTagCount(); i++) {
                for (int j = 0; j < this.tagLengthInLong; j++) {
                    dos.writeLong(this.tags[j][i]);
                }
                dos.writeByte(this.getTagLength(i));
                dos.writeInt(this.getGChr(i));
                dos.writeInt(this.getGPos(i));
                dos.writeByte(this.getIfPAV(i));
                dos.writeFloat(this.getPrediction(i));
            }
            dos.flush();
            dos.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void writeFastQ (String outfileS) {
        String defaultQualityS = null;
        for (int i = 0; i < this.getTagSizeInLong()*32; i++) {
            defaultQualityS+="f";
        }
        long[] t = new long[this.getTagSizeInLong()];
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
            for (int i = 0; i < this.getTagCount(); i++) {
                bw.write("@"+String.valueOf(i));
                bw.newLine();
                for (int j = 0; j < this.getTagSizeInLong(); j++) {
                    t[j] = this.tags[j][i];
                }
                bw.write(BaseEncoder.getSequenceFromLong(t).substring(0, this.getTagLength(i)));
                bw.newLine();
                bw.write("+");
                bw.newLine();
                bw.write(defaultQualityS.substring(0, this.getTagLength(i)));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
