package net.maizegenetics.taxa;


import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.AbstractTableReport;

import java.io.PrintWriter;
import java.io.Serializable;
import java.util.*;

import org.apache.commons.codec.language.DoubleMetaphone;
import org.apache.commons.codec.language.Metaphone;
import org.apache.commons.codec.language.RefinedSoundex;
import org.apache.commons.codec.language.Soundex;


/**
 * User: Ed
 * Date: Mar 30, 2005
 * Time: 1:39:47 PM
 */
public class IdentifierSynonymizer extends AbstractTableReport implements Serializable, TableReport {

    HashMap<String,Integer> idSynonyms = new HashMap<>();    //TODO needs to be entirely updated to new collections
    private TaxaList referenceIDGroup;
    private int unmatchCount = 0;
    private int technique = 0;
    private double globalMin = Double.POSITIVE_INFINITY;
    private double globalMax = Double.NEGATIVE_INFINITY;

    public IdentifierSynonymizer(TaxaList preferredTaxa, TaxaList[] alternateTaxaSets) {
        init(preferredTaxa, alternateTaxaSets);
    }
    public IdentifierSynonymizer(TaxaList preferredTaxa, TaxaList[] alternateTaxaSets,int technique) {
        this.technique = technique;
        init(preferredTaxa, alternateTaxaSets);
    }

    public IdentifierSynonymizer(TaxaList preferredTaxa, TaxaList alternateTaxa) {
        TaxaList[] alternateTaxaSets = new TaxaList[1];
        alternateTaxaSets[0] = alternateTaxa;
        init(preferredTaxa, alternateTaxaSets);
    }

    private void init(TaxaList preferredTaxa, TaxaList[] alternateTaxaSets) {
        //referenceIDGroup=preferredTaxa;
        referenceIDGroup = preferredTaxa;
        Taxon currID;
        //Load up the synonym table with all the known names
        for (int i = 0; i < referenceIDGroup.numberOfTaxa(); i++) {
            idSynonyms.put(referenceIDGroup.taxaName(i), i);
        }
        //Find the unknown names and place them in a list
        for (int a = 0; a < alternateTaxaSets.length; a++) {
            for (int i = 0; i < alternateTaxaSets[a].numberOfTaxa(); i++) {
                currID = alternateTaxaSets[a].get(i);
                if (idSynonyms.containsKey(currID.getName()) == false) {
                    ArrayList<String> theBest = findBestMatch(currID.toString());
                    if (theBest.size() == 1) {
                        String bs = (String) theBest.get(0);
                        int indexOfBest = referenceIDGroup.indexOf(bs);
                        idSynonyms.put(currID.toString(), indexOfBest);
                    } else {
                        idSynonyms.put(currID.toString(), -1);
                        unmatchCount++;
                    }
                }
                else {
                    globalMin = 0;
                }
            }
        }
    }
   

    private ArrayList<String> findBestMatch(String unmatchedString) {
        ArrayList<String> bestMatches = new ArrayList<>();
        double maxScore = -1;
        double minScore = Double.POSITIVE_INFINITY;
        double sm;
        int levelOfRestriction = 0;
        boolean ignoreCase = true, ignoreWhite = false, ignorePunc = false;
        while ((bestMatches.size() != 1) && (levelOfRestriction < 4)) {
            switch (levelOfRestriction) {
                case 1:
                    ignoreCase = true;
                    break;
                case 2:
                    ignoreWhite = true;
                    break;
                case 3:
                    ignorePunc = true;
                    break;
            }
            /*
            for (int i = 0; i < referenceIDGroup.numberOfTaxa(); i++) {
                sm = getScore(referenceIDGroup.taxaName(i), unmatchedString, ignoreCase, ignoreWhite, ignorePunc,technique);
                //sm = scoreMatch(referenceIDGroup.taxaName(i), unmatchedString, ignoreCase, ignoreWhite, ignorePunc);
                if (sm > maxScore) {
                    bestMatches.clear();
                    bestMatches.add(referenceIDGroup.taxaName(i));
                    maxScore = sm;
                } else if (sm == maxScore) {
                    bestMatches.add(referenceIDGroup.taxaName(i));
                }
            }*/
            
            for (int i = 0; i < referenceIDGroup.numberOfTaxa(); i++) {
                sm = getScore(referenceIDGroup.taxaName(i), unmatchedString, ignoreCase, ignoreWhite, ignorePunc,technique);
                
                //sm = scoreMatch(referenceIDGroup.taxaName(i), unmatchedString, ignoreCase, ignoreWhite, ignorePunc);
                if (sm < minScore) {
                    bestMatches.clear();
                    bestMatches.add(referenceIDGroup.taxaName(i));
                    minScore = sm;
                    if(minScore<globalMin) {
                        System.out.println(minScore);
                        globalMin = minScore;
                    } 
                } else if (sm == minScore) {
                    bestMatches.add(referenceIDGroup.taxaName(i));
                }
            }
            if(minScore>globalMax) {
                globalMax = minScore;
            }
            
            levelOfRestriction++;
        }
        return bestMatches;
    }

    public ArrayList<String> findOrderedMatches(String unmatchedString, int levelOfRestriction) {
        SortedMap<Double,String> theSortMap = new TreeMap<>();
        double sm;
        boolean ignoreCase = false, ignoreWhite = false, ignorePunc = false;
        if (levelOfRestriction > 0) {
            ignoreCase = true;
        }
        if (levelOfRestriction > 1) {
            ignoreWhite = true;
        }
        if (levelOfRestriction > 2) {
            ignorePunc = true;
        }
        for (int i = 0; i < referenceIDGroup.numberOfTaxa(); i++) {
            //sm = scoreMatch(referenceIDGroup.taxaName(i), unmatchedString, ignoreCase, ignoreWhite, ignorePunc);
            sm = getScore(referenceIDGroup.taxaName(i), unmatchedString, ignoreCase, ignoreWhite, ignorePunc,technique);
            
            //theSortMap.put(1 - sm - ((double) i / 100000.0), referenceIDGroup.taxaName(i));
            theSortMap.put(sm - ((double) i / 100000.0), referenceIDGroup.taxaName(i));
        }
        return new ArrayList<>(theSortMap.values());
    }

    public static double getScore(String s1, String s2, boolean ignoreCase, boolean ignoreWhite, boolean ignorePunc, int technique) {
        double score = 0.0;
      
        //dice need to do a 1- as high similarity = low distance
        if(technique==0) {
            score = 1.0 - scoreMatch(s1,s2,ignoreCase,ignoreWhite,ignorePunc);
        }
        //String edit
        else if(technique==1) {
            score = editDistanceScoreMatch(s1,s2,ignoreCase,ignoreWhite,ignorePunc);
        }
        //DTW with hamming
        else if(technique==2) {
            score = dtwDist(s1,s2,"hamming",ignoreCase,true,ignorePunc);
        }
        //DTW with keyboard dist
        else if(technique==3) {
            score = dtwDist(s1,s2,"key",ignoreCase,true,ignorePunc);
        }
        //Hamming with soundex
        else if(technique==4) {
            score = hammingDistSoundex(s1,s2,ignoreCase,ignoreWhite,ignorePunc);
        }
        //Dice with metaphone  need to do a 1- as high similarity = low distance
        else if(technique==5) {
            score = 1 - diceWithMetaphone(s1,s2,ignoreCase,ignoreWhite,ignorePunc);
        }
        //Edit Distance with metaphone
        else if(technique==6) {
            score=editWithMetaphone(s1,s2,ignoreCase,ignoreWhite,ignorePunc);
        }
        return score;
    }
    
    public static double hammingDistSoundex(String s1, String s2, boolean ignoreCase, boolean ignoreWhite, boolean ignorePunc) {
        s1 = soundex2(s1, true, true, true);
        s2 = soundex2(s2, true, true, true);
        int sum = 0;
        for(int i = 0; i<s1.length();i++) {
            sum += hammingDist(s1.charAt(i), s2.charAt(i));
        }
        return (double)sum;
    }
    public static double diceWithMetaphone(String s1, String s2, boolean ignoreCase, boolean ignoreWhite, boolean ignorePunc) {
        s1 = metaphone2(s1,true,true,true);
        s2 = metaphone2(s2,true,true,true);
        return scoreMatch(s1,s2,ignoreCase,ignoreWhite,ignorePunc);
    }
    
    public static double editWithMetaphone(String s1, String s2, boolean ignoreCase, boolean ignoreWhite, boolean ignorePunc) {
        s1 = metaphone2(s1,true,true,true);
        s2 = metaphone2(s2,true,true,true);
        return editDistanceScoreMatch(s1,s2,ignoreCase,ignoreWhite,ignorePunc);
    }
    private double scoreMatch2(String s1, String s2, boolean ignoreCase, boolean ignoreWhite, boolean ignorePunc) {
        //idea from http://www.catalysoft.com/articles/StrikeAMatch.html?article=How_to_Strike_a_Match_15
        //this is faster but it can be tricked if there are long runs of characters in s1
        int score = 0;
        double sm;
        s1 = cleanName(s1, ignoreCase, ignoreWhite, ignorePunc);
        s2 = cleanName(s2, ignoreCase, ignoreWhite, ignorePunc);
//        System.out.println("s1="+s1+"  s2="+s2);
        for (int c1 = 0; c1 < (s1.length() - 1); c1++) {
            for (int c2 = 0; c2 < (s2.length() - 1); c2++) {
                if ((s1.charAt(c1) == s2.charAt(c2)) && (s1.charAt(c1 + 1) == s2.charAt(c2 + 1))) {
                    score++;
                    break;
                }
            }
        }
        sm = (2.0 * (double) score) / (s1.length() + s2.length() - 2);
        return sm;
    }

    /** @return lexical similarity value in the range [0,1] */
    public static double scoreMatch(String s1, String s2, boolean ignoreCase, boolean ignoreWhite, boolean ignorePunc) {
        //idea from http://www.catalysoft.com/articles/StrikeAMatch.html?article=How_to_Strike_a_Match_15
        //this is slower but it will not be tricked if there are long runs of characters in s1
        s1 = cleanName(s1, ignoreCase, ignoreWhite, ignorePunc);
        s2 = cleanName(s2, ignoreCase, ignoreWhite, ignorePunc);
        ArrayList<String> pairs1 = letterPairs(s1);
        ArrayList<String> pairs2 = letterPairs(s2);

        int intersection = 0;
        int union = pairs1.size() + pairs2.size();
        for (int i = 0; i < pairs1.size(); i++) {
            Object pair1 = pairs1.get(i);
            for (int j = 0; j < pairs2.size(); j++) {
                Object pair2 = pairs2.get(j);
                if (pair1.equals(pair2)) {
                    intersection++;
                    pairs2.remove(j);
                    break;
                }
            }
        }
        return (2.0 * intersection) / union;
    }
    
    public static double editDistanceScoreMatch(String s1, String s2, boolean ignoreCase, boolean ignoreWhite, boolean ignorePunc) {
        
        s1 = cleanName(s1, ignoreCase, ignoreWhite, ignorePunc);
        s2 = cleanName(s2, ignoreCase, ignoreWhite, ignorePunc);
        
        double[][] editMatrix = new double[s1.length()][s2.length()];
        
        //Init first row and column of editMatrix
        editMatrix[0][0] = (s1.charAt(0) == s2.charAt(0))?0.0:1.0;
        for(int i = 1; i<editMatrix.length; i++) {
            editMatrix[i][0] = (s1.charAt(i)==s2.charAt(0))?editMatrix[i-1][0]: editMatrix[i-1][0]+1;
        }
        for(int i = 1; i<editMatrix[0].length;i++) {
            editMatrix[0][i] = (s1.charAt(0) == s2.charAt(i))? editMatrix[0][i-1]:editMatrix[0][i-1]+1;
        }
        
        //Fill in the rest of the matrix
        for(int i = 1; i < editMatrix.length; i++) {
            for(int j = 1; j< editMatrix[i].length;j++) {
                if(s1.charAt(i) == s2.charAt(j)) {
                    editMatrix[i][j] = editMatrix[i-1][j-1];
                }
                else {
                    double diagCost = editMatrix[i-1][j-1];
                    double upCost = editMatrix[i-1][j];
                    double leftCost = editMatrix[i][j-1];
                    
                    if(diagCost <= upCost && diagCost <= leftCost) {
                        //Substitution
                        editMatrix[i][j] = diagCost + 1;
                    }
                    else if(upCost <= leftCost && upCost <= diagCost) {
                        //Deletion
                        editMatrix[i][j] = upCost + 1;
                    }
                    else {
                        //Insertion
                        editMatrix[i][j] = leftCost + 1;
                    }
                }
            }
        }
        //Return the last position of the matrix as similarity score
        /*
        for(int i = 0; i<editMatrix.length;i++) {
            for(int j = 0; j<editMatrix[i].length;j++) {
                System.out.print(editMatrix[i][j]+",");
            }
            System.out.println();
        }
        System.out.println();
        if(s2.equals("VA22")) {
            System.out.println();
        }
        */
        return editMatrix[editMatrix.length-1][editMatrix[editMatrix.length-1].length-1];
    }

    public static String refinedSoundex(String s1, boolean ignoreCase, boolean ignoreWhite, boolean ignorePunc) {
        s1 = cleanName(s1, true, true, true);
        
        StringBuilder soundexCode = new StringBuilder();
        HashMap<Character, Character> consMap= new HashMap<Character,Character>();
        consMap.put('B', '1');
        consMap.put('C', '3');
        consMap.put('D', '6');
        consMap.put('F', '2');
        consMap.put('G', '4');
        consMap.put('J', '4');
        consMap.put('K', '3');
        consMap.put('L', '7');
        consMap.put('M', '8');
        consMap.put('N', '8');
        consMap.put('P', '1');
        consMap.put('Q', '5');
        consMap.put('R', '9');
        consMap.put('S', '3');
        consMap.put('T', '6');
        consMap.put('V', '2');
        consMap.put('X', '5');
        consMap.put('Z', '5');
        consMap.put('0', 'a');
        consMap.put('1', 'b');
        consMap.put('2', 'c');
        consMap.put('3', 'd');
        consMap.put('4', 'e');
        consMap.put('5', 'f');
        consMap.put('6', 'g');
        consMap.put('7', 'h');
        consMap.put('8', 'i');
        consMap.put('9', 'j');
        soundexCode.append(s1.charAt(0));
        if(s1.charAt(0) != 'H' && s1.charAt(0) != 'W') {
            soundexCode.append(consMap.get(s1.charAt(0)));
        }
        for(int i = 1; i< s1.length();i++) {
            if(s1.charAt(i) != 'W' && s1.charAt(i) != 'H' && s1.charAt(i) != 'A' && s1.charAt(i) != 'E' && s1.charAt(i) != 'I' && s1.charAt(i) != 'O' && s1.charAt(i) != 'U' && s1.charAt(i) != 'Y') {
                char currentCharacter = consMap.get(s1.charAt(i));
                char prevCharacter = soundexCode.charAt(soundexCode.length()-1);
                if(currentCharacter != prevCharacter) {
                    soundexCode.append(currentCharacter);
                }
            }
            else {
                
                char prevCharacter = soundexCode.charAt(soundexCode.length()-1);
                if(prevCharacter != '0') {
                    soundexCode.append('0');
                }
                //soundexCode.append('0');
            }
        }
        return soundexCode.toString();
    }

    public static String metaphone(String s1, boolean ignoreCase, boolean ignoreWhite, boolean ignorePunc) {
        s1 = cleanName(s1, true, true, true);
        
        //Step1 remove all repeated letters except for duplicate C characters
        StringBuilder sb = new StringBuilder();
        char prevChar = s1.charAt(0);
        sb.append(s1.charAt(0));
        for(int i = 1; i<s1.length();i++) {
            if(prevChar != s1.charAt(i) || s1.charAt(i) == 'C') {
                prevChar = s1.charAt(i);
                sb.append(s1.charAt(i));
            }
        }
        //Step2 change first to characters if KN, GN, PN, AE, WR
        s1 = sb.toString();
        sb = new StringBuilder();
        if(s1.startsWith("KN")) {
            s1.replaceFirst("KN","N");
        }
        else if(s1.startsWith("GN")) {
            s1.replaceFirst("GN", "N");
        }
        else if(s1.startsWith("PN")) {
            s1.replaceFirst("PN", "N");
        }
        else if(s1.startsWith("AE")) {
            s1.replaceFirst("AE", "E");
        }
        else if(s1.startsWith("WR")) {
            s1.replaceFirst("WR", "R");
        }
        
        //Step3 remove B character from end if preceded by M
        if(s1.endsWith("MB")) {
            s1 = s1.substring(s1.length()-2, s1.length()-1);
        }
        
        //Step4 Replace occurrences of C depending on what surrounds it
        s1.replaceAll("CIA", "XIA");
        s1.replaceAll("SCH", "SKH");
        s1.replaceAll("CH", "XH");
        s1.replaceAll("CI", "SI");
        s1.replaceAll("CE", "SE");
        s1.replaceAll("CY", "SY");
        s1.replaceAll("C", "K");
        
        //Step 5 Replace D with J or T
        s1.replaceAll("DGE", "JGE");
        s1.replaceAll("DGY", "JGY");
        s1.replaceAll("DGI", "JGY");
        s1.replaceAll("D", "T");
        
        //Step 6 Replace GH with H as long as its not before a vowel or ending the word
        s1.replaceAll("GH[^AEIOU|]","H");
        
        //Step 7 fix GN or GNED ending words
        s1.replaceAll("GN$", "N");
        s1.replaceAll("GNED$", "NED");
        
        //Step 8 replace G characters
        s1.replaceAll("GI", "JI");
        s1.replaceAll("GE", "JE");
        s1.replaceAll("GY", "JY");
        s1.replaceAll("G", "K");
        
        //Step 9 remove H after a vowel but not before
        for(int i = 0; i<s1.length()-2;i++) {
            if(s1.charAt(i) == 'A' || s1.charAt(i) == 'E' || s1.charAt(i) == 'I' || s1.charAt(i) == 'O' || s1.charAt(i) == 'U') {
                if(s1.charAt(i+1) != 'H') {
                    sb.append(s1.charAt(i));
                }
                else {
                    if(sb.charAt(i+2) == 'A' || s1.charAt(i+2) != 'E' || s1.charAt(i+2) == 'I' || s1.charAt(i+2) == 'O' || s1.charAt(i+2) == 'U') {
                       sb.append(s1.charAt(i+1));
                       sb.append(s1.charAt(i+2));
                       i+=2;
                    }
                    else {
                        //skip H character
                        i+=1;
                    }
                }
            }
            else {
                sb.append(s1.charAt(i));
            }
        }
        sb.append(s1.charAt(s1.length()-2));
        if(s1.charAt(s1.length()-1)=='H') {
            if(s1.charAt(s1.length()-2)!='A' && s1.charAt(s1.length()-2)!='E' && s1.charAt(s1.length()-2)!='I' && s1.charAt(s1.length()-2)!='O' && s1.charAt(s1.length()-2)!='U') {
                sb.append('H');
            }
        }
        s1 = sb.toString();
        sb = new StringBuilder();
        
        //Step10 Various transformations
        s1.replaceAll("CK", "K");
        s1.replaceAll("PH", "F");
        s1.replaceAll("Q", "K");
        s1.replaceAll("V", "F");
        s1.replaceAll("Z", "S");
        
        //Step11 S to X
        s1.replaceAll("SH","XH");
        s1.replaceAll("SIO", "XIO");
        s1.replaceAll("SIA", "XIA");
        
        //Step12 replace T
        s1.replaceAll("TIA", "XIA");
        s1.replaceAll("TIO", "XIO");
        s1.replaceAll("TH", "0");
        s1.replaceAll("TCH", "CH");
        
        //Step13 remove WH if it is at the start/ remove all W if no vowel after it
        if(s1.charAt(0)=='W' && s1.charAt(1)=='H') {
            s1="W"+s1.substring(2, s1.length());
        }
        for(int i = 0; i<s1.length()-1;i++) {
            if(s1.charAt(i)=='W') {
                if(s1.charAt(i+1) != 'A' && s1.charAt(i+1) != 'E' && s1.charAt(i+1) != 'I' && s1.charAt(i+1) != 'O' && s1.charAt(i+1) != 'U') {
                    sb.append(s1.charAt(i));
                }
            }
            else {
                sb.append(s1.charAt(i));
            }
        }
        if(s1.charAt(s1.length()-1)!='W') {
            sb.append(s1.charAt(s1.length()-1));
        }
        s1 = sb.toString();
        sb = new StringBuilder();
        
        //Step14 replace X
        if(s1.charAt(0)=='X') {
            s1 = "S"+s1.substring(1, s1.length());
        }
        s1.replaceAll("X", "KS");
        
        //Step15 Remove Y which are not before a vowel
        for(int i = 0; i<s1.length()-1;i++) {
            if(s1.charAt(i)=='Y') {
                if(s1.charAt(i+1) != 'A' && s1.charAt(i+1) != 'E' && s1.charAt(i+1) != 'I' && s1.charAt(i+1) != 'O' && s1.charAt(i+1) != 'U') {
                    sb.append(s1.charAt(i));
                }
            }
            else {
                sb.append(s1.charAt(i));
            }
        }
        if(s1.charAt(s1.length()-1)!='Y') {
            sb.append(s1.charAt(s1.length()-1));
        }
        s1 = sb.toString();
        sb = new StringBuilder();
        
        //Step16 remove all vowels except if word starts with vowel
        sb.append(s1.charAt(0));
        for(int i = 1; i<s1.length();i++) {
            if(s1.charAt(i) != 'A' && s1.charAt(i) != 'E' && s1.charAt(i) != 'I' && s1.charAt(i) != 'O' && s1.charAt(i) != 'U') {
                sb.append(s1.charAt(i));
            }
        }
        s1 = sb.toString();
        return s1;
    }

    public static String metaphone2(String s1, boolean ignoreCase, boolean ignoreWhite, boolean ignorePunc) {
        s1 = cleanName(s1, true, true, true);
        if(s1.equals("")) {
            return "";
        }
        //Parse out numbers i.e. ab12cd3e becomes ["ab","12","cd","3","e"] then metaphone all letter only strings
        ArrayList<String> parsed = new ArrayList<String>();
        String current = "";
        boolean digitMode = false;
        if(Character.isDigit(s1.charAt(0))) {
            current+=s1.charAt(0);
            digitMode = true;
        }
        for(int i = 0; i<s1.length();i++) {
            if(Character.isDigit(s1.charAt(i)) && !digitMode) {
                parsed.add(current);
                current = ""+s1.charAt(i);
                digitMode = true;
            }
            else if(!Character.isDigit(s1.charAt(i)) && digitMode) {
                parsed.add(current);
                current = ""+s1.charAt(i);
                digitMode = false;
            }
            else {
                current += s1.charAt(i);
            }
        }
        parsed.add(current);
        Metaphone metaphone = new Metaphone();
        String encodedString = "";
        for(int i = 0; i<parsed.size();i++) {
            if(!Character.isDigit(parsed.get(i).charAt(0))) {
               encodedString += metaphone.encode(parsed.get(i)); 
            }
        }
        return encodedString;
        //return metaphone.metaphone(s1);
    }
    public static String soundex2(String s1, boolean ignoreCase, boolean ignoreWhite, boolean ignorePunc) {
        s1 = cleanName(s1, true, true, true);
        Soundex soundex = new Soundex();
        return soundex.soundex(s1);
    }
    public static String refinedSoundex2(String s1, boolean ignoreCase, boolean ignoreWhite, boolean ignorePunc) {
        s1 = cleanName(s1, true, true, true);
        RefinedSoundex soundex = new RefinedSoundex();
        return soundex.soundex(s1);
    }
    /*
     * Method to determine keyboard distance.
     * Idea borrowed from http://search.cpan.org/~krburton/String-KeyboardDistance-1.01/KeyboardDistance.pm
     * This needs refractoring it is messy
     */
    private static int keyboardDist(char firstChar, char secondChar) {
        /*    | 0   1   2   3   4   5   6   7   8   9   10  11  12  13
            --+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
            0 | ` | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 0 | - | = |   |
            1 |   | q | w | e | r | t | y | u | i | o | p | [ | ] | \ |
            2 |   | a | s | d | f | g | h | j | k | l | ; | ' |   |   |
            3 |   | z | x | c | v | b | n | m | , | . | / |   |   |   |
            --+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
        */
        HashMap<Character,Integer[]> map = new HashMap<>();
        map.put('`', new Integer[]{0,0});
        map.put('~', new Integer[]{0,0});
        map.put('1', new Integer[]{0,1});
        map.put('!', new Integer[]{0,1});
        map.put('2', new Integer[]{0,2});
        map.put('@', new Integer[]{0,2});
        map.put('3', new Integer[]{0,3});
        map.put('#', new Integer[]{0,3});
        map.put('4', new Integer[]{0,4});
        map.put('$', new Integer[]{0,4});
        map.put('5', new Integer[]{0,5});
        map.put('%', new Integer[]{0,5});
        map.put('6', new Integer[]{0,6});
        map.put('^', new Integer[]{0,6});
        map.put('7', new Integer[]{0,7});
        map.put('&', new Integer[]{0,7});
        map.put('8', new Integer[]{0,8});
        map.put('*', new Integer[]{0,8});
        map.put('9', new Integer[]{0,9});
        map.put('(', new Integer[]{0,9});
        map.put('0', new Integer[]{0,10});
        map.put(')', new Integer[]{0,10});
        map.put('-', new Integer[]{0,11});
        map.put('_', new Integer[]{0,11});
        map.put('=', new Integer[]{0,12});
        map.put('+', new Integer[]{0,12});
        
        map.put('q', new Integer[]{1,1});
        map.put('Q', new Integer[]{1,1});
        map.put('w', new Integer[]{1,2});
        map.put('W', new Integer[]{1,2});
        map.put('e', new Integer[]{1,3});
        map.put('E', new Integer[]{1,3});
        map.put('r', new Integer[]{1,4});
        map.put('R', new Integer[]{1,4});
        map.put('t', new Integer[]{1,5});
        map.put('T', new Integer[]{1,5});
        map.put('y', new Integer[]{1,6});
        map.put('Y', new Integer[]{1,6});
        map.put('u', new Integer[]{1,7});
        map.put('U', new Integer[]{1,7});
        map.put('i', new Integer[]{1,8});
        map.put('I', new Integer[]{1,8});
        map.put('o', new Integer[]{1,9});
        map.put('O', new Integer[]{1,9});
        map.put('p', new Integer[]{1,10});
        map.put('P', new Integer[]{1,10});
        map.put('[', new Integer[]{1,11});
        map.put('{', new Integer[]{1,11});
        map.put(']', new Integer[]{1,12});
        map.put('}', new Integer[]{1,12});
        map.put('\\', new Integer[]{1,13});
        map.put('|', new Integer[]{1,13});
        
        map.put('a', new Integer[]{2,1});
        map.put('A', new Integer[]{2,1});
        map.put('s', new Integer[]{2,2});
        map.put('S', new Integer[]{2,2});
        map.put('d', new Integer[]{2,3});
        map.put('D', new Integer[]{2,3});
        map.put('f', new Integer[]{2,4});
        map.put('F', new Integer[]{2,4});
        map.put('g', new Integer[]{2,5});
        map.put('G', new Integer[]{2,5});
        map.put('h', new Integer[]{2,6});
        map.put('H', new Integer[]{2,6});
        map.put('j', new Integer[]{2,7});
        map.put('J', new Integer[]{2,7});
        map.put('k', new Integer[]{2,8});
        map.put('K', new Integer[]{2,8});
        map.put('l', new Integer[]{2,9});
        map.put('L', new Integer[]{2,9});
        map.put(';', new Integer[]{2,10});
        map.put(':', new Integer[]{2,10});
        map.put('\'', new Integer[]{2,11});
        map.put('"', new Integer[]{2,11});
        
        map.put('z', new Integer[]{3,1});
        map.put('Z', new Integer[]{3,1});
        map.put('x', new Integer[]{3,2});
        map.put('X', new Integer[]{3,2});
        map.put('c', new Integer[]{3,3});
        map.put('C', new Integer[]{3,3});
        map.put('v', new Integer[]{3,4});
        map.put('V', new Integer[]{3,4});
        map.put('b', new Integer[]{3,5});
        map.put('B', new Integer[]{3,5});
        map.put('n', new Integer[]{3,6});
        map.put('N', new Integer[]{3,6});
        map.put('m', new Integer[]{3,7});
        map.put('M', new Integer[]{3,7});
        map.put(',', new Integer[]{3,8});
        map.put('<', new Integer[]{3,8});
        map.put('.', new Integer[]{3,9});
        map.put('>', new Integer[]{3,9});
        map.put('/', new Integer[]{3,10});
        map.put('?', new Integer[]{3,10});
        
       Integer[] coords1 = map.get(firstChar);
       Integer[] coords2 = map.get(secondChar);
       
       //calculate manhattan distance between the characters
       return Math.abs(coords1[0] - coords2[0]) + Math.abs(coords1[1] - coords2[1]);
        
    }
    
    private static int hammingDist(char firstChar, char secondChar) {
        if(firstChar == secondChar) {
            return 0;
        }
        else{
            return 1;
        }   
    }
   
    /*
     * Simple implementation of Dynamic Time Warping distance
     * 
     * Currently uses KeyboardDistance as the distance measurement
     * More will be implemented soon
     */
    
    private static double dtwDist(String str1, String str2, String distMeas,boolean ignoreCase, boolean ignoreWhite, boolean ignorePunc) {
        str1 = cleanName(str1,ignoreCase,ignoreWhite,ignorePunc);
        str2 = cleanName(str2,ignoreCase,ignoreWhite,ignorePunc);
        double[][] costMat = new double[str1.length()+1][str2.length()+1];
        //Initialize arrays
        for(int i = 0; i<costMat.length;i++) {
            costMat[i][0] = Double.POSITIVE_INFINITY;
        }
        for(int i = 0; i<costMat[0].length;i++) {
            costMat[0][i] = Double.POSITIVE_INFINITY;
        }
        
        double currDist = 0.0;
        if(distMeas.equals("key")) {
            costMat[1][1] = (double)keyboardDist(str1.charAt(0),str2.charAt(0));
        }
        else if(distMeas.equals("hamming")) {
            costMat[1][1] = (double)hammingDist(str1.charAt(0), str2.charAt(0));
        }
        for(int i = 2;i<costMat[1].length;i++) {
            if(distMeas.equals("key")) {
                costMat[1][i] = costMat[1][i-1]+(double)keyboardDist(str1.charAt(0),str2.charAt(i-1));
            }
            else if(distMeas.equals("hamming")) {
                costMat[1][i] = costMat[1][i-1]+(double)hammingDist(str1.charAt(0), str2.charAt(i-1));
            }        
        }
        for(int i = 2; i<costMat.length;i++) {
            if(distMeas.equals("key")) {
                costMat[i][1] = costMat[i-1][1]+(double)keyboardDist(str1.charAt(i-1),str2.charAt(0));
            }
            else if(distMeas.equals("hamming")) {
                costMat[i][1] = costMat[i-1][1]+(double)hammingDist(str1.charAt(i-1), str2.charAt(0));
            }
        }
        for(int i = 2; i<costMat.length;i++) {
            for(int j = 2; j<costMat[i].length;j++) {
                currDist = 0.0;
                //get distance
                if(distMeas.equals("key")) {
                    currDist = (double)keyboardDist(str1.charAt(i-1),str2.charAt(j-1));
                    //System.out.println(currDist);
                }
                else if(distMeas.equals("hamming")) {
                    currDist = (double)hammingDist(str1.charAt(i-1), str2.charAt(j-1));
                }
                
                if(costMat[i-1][j-1] < costMat[i-1][j] && costMat[i-1][j-1] < costMat[i][j-1]) {
                    costMat[i][j] = costMat[i-1][j-1] + currDist;
                }
                else if(costMat[i-1][j]<costMat[i-1][j-1] && costMat[i-1][j] < costMat[i][j-1]) {
                    costMat[i][j] = costMat[i-1][j] + currDist;
                }
                else {
                    costMat[i][j] = costMat[i][j-1] + currDist;
                }
            }
        }
        //return end point
        
        /*
        if(distMeas.equals("key")) {
        System.out.println(str2);
        for(int i = 1; i<costMat.length;i++) {
            System.out.print(str1.charAt(i-1)+": ");
            for(int j = 1; j<costMat[i].length;j++) {
                System.out.print(costMat[i][j]+",");
            }
            System.out.println();
        }
        System.out.println("END:"+costMat[costMat.length-1][costMat[costMat.length-1].length-1]);
        }
        */
        return costMat[costMat.length-1][costMat[costMat.length-1].length-1];
    }
    /** @return an array of adjacent letter pairs contained in the input string */
    private static ArrayList<String> letterPairs(String str) {
        ArrayList<String> allPairs = new ArrayList<>();
        //int numPairs = str.length()-1;
        //String[] pairs = new String[numPairs];
        for (int i = 0; i < (str.length() - 1); i++) {
            allPairs.add(str.substring(i, i + 2));
        }
        return allPairs;
    }

    private static String cleanName(String s, boolean ignoreCase, boolean ignoreWhite, boolean ignorePunc) {
        if (ignoreCase) {
            s = s.toUpperCase();
        }
        //StringBuffer sb=new StringBuffer(s);
        //int x;
        if (ignoreWhite) {
            s.replaceAll("\\s", "");
        // while((x=sb.indexOf(" "))>=0) {sb.deleteCharAt(x);}
        }
        if (ignorePunc) {
            //           s=s.replaceAll("\\W","");
            s = s.replaceAll("[^a-zA-Z0-9]", "");
        }
        // sb=new StringBuffer(s);
        return s;
    }

    public void changeAlignmentIdentifiers(TaxaList alternateIdGroups) {
        TaxaList[] aidg = new TaxaList[1];
        aidg[0] = alternateIdGroups;
        changeAlignmentIdentifiers(aidg[0]);
    }

    public void changeAlignmentIdentifiers(TaxaList[] alternateIdGroups) {
        Taxon currID;
        for (int a = 0; a < alternateIdGroups.length; a++) {
            TaxaListBuilder tLB=new TaxaListBuilder();
            for (int i = 0; i < alternateIdGroups[a].numberOfTaxa(); i++) {
                currID = alternateIdGroups[a].get(i);
                if (getPreferredIndex(currID.getName()) > -1) {
                    tLB.add(new Taxon.Builder(getPreferredName(currID.getName())).build());
                } else {
                    tLB.add(new Taxon.Builder(currID).build());
                }
            }
        }
    }

    public String toString() {
        String s = "Synonym Table\n" + idSynonyms.toString() + "\n\n";
        return s;    //To change body of overridden methods use File | Settings | File Templates.
    }

    public String getPreferredName(String theID) {
        int index = getPreferredIndex(theID);
        if (index > -1) {
            return referenceIDGroup.taxaName(index);
        } else {
            return "";
        }
    }

    public int getPreferredIndex(String theID) {
        Object index = idSynonyms.get(theID);
        if (index == null) {
            return -1;
        } else {
            return ((Integer) index).intValue();
        }
    }

    public Taxon getPreferredIdentifier(Taxon theID) {
        int index = getPreferredIndex(theID.getName());
        if (index > -1) {
            return referenceIDGroup.get(index);
        } else {
            return null;
        }
    }

    public void deleteByThreshold(double threshold) {
        Object[] keyArray = idSynonyms.keySet().toArray();
        String synName, realName;
        double score;
        for (int i = 0; i < keyArray.length; i++) {
            synName = "" + (String) keyArray[i];
            if (getPreferredIndex(synName) > -1) {
                realName = "" + getPreferredName(synName);
                score = scoreMatch(synName, realName, true, false, false);
                if (score < threshold) {
                    idSynonyms.put(synName, new Integer(-1));
                }
            }
        }
    }

    public boolean setRealName(String synName, String realName) {
        int synID = getPreferredIndex(synName);
        int realID = referenceIDGroup.indexOf(realName);
        if ((synID > -1) && (realID > -1)) {
            realName = "" + getPreferredName(synName);
            idSynonyms.put(synName, new Integer(realID));
            return true;
        } else {
            return false;
        }
    }

    public boolean setRealID(String synName, int realID) {
        if ((realID <= referenceIDGroup.numberOfTaxa()) && (realID > -2)) {
            idSynonyms.put(synName, new Integer(realID));
            return true;
        } else {
            return false;
        }
    }

    public Object[] getRealNames() {
        Object[] idArray = new Object[referenceIDGroup.numberOfTaxa()];
        for (int i = 0; i < referenceIDGroup.numberOfTaxa(); i++) {
            idArray[i] = referenceIDGroup.get(i).toString();
        }
        return idArray;
    }

    public void report(PrintWriter out) {
        //String s="Synonym Table\n"+idSynonyms.toString()+"\n\n"+"Unmatched\n"+unmatchedIDs.toString();
        out.println("Synonym Table");
        out.println(idSynonyms.size() + " unique matches");
        out.println(unmatchCount + " unmatched:");
    }

    public Object[] getTableColumnNames() {
        String[] cn = new String[4];
        cn[0] = "TaxaSynonym";
        cn[1] = "TaxaRealName";
        cn[2] = "RefIDNum";
        cn[3] = "MatchScore";
        return cn;
    }

    /**
     * Returns specified row.
     *
     * @param row row number
     *
     * @return row
     */
    public Object[] getRow(long rowLong) {

        int row = (int) rowLong;
        Object[] data = new Object[4];
        Object[] keyArray = idSynonyms.keySet().toArray();
        data[0] = (String) keyArray[row];
        data[1] = getPreferredName((String) keyArray[row]);
        data[2] = "" + getPreferredIndex((String) keyArray[row]);
        if(technique==0) {
            data[3] = "" + scoreMatch("" + data[0], "" + data[1], true, false, false);
        }
        else {
            if(Integer.parseInt((String)data[2])==-1) {
                //data[3] = ""+Double.POSITIVE_INFINITY;
                data[3] = ""+0;
            }
            else {
                //System.out.println("GlMin: "+globalMin + " GLMax: "+globalMax+" Val: "+getScore("" + data[0], "" + data[1], true, false, false,technique));
                data[3] = "" + (1.0-((getScore("" + data[0], "" + data[1], true, false, false,technique) - globalMin)/(globalMax-globalMin)));
                //data[3] = "" + getScore("" + data[0], "" + data[1], true, false, false,technique);
            }
        }
        return data;

    }

    public void deleteElements(Object[] key) {
        for (int i = 0; i < key.length; i++) {
            idSynonyms.remove(key[i]);
        }
    }

    public String getTableTitle() {
        return "Taxa Synonym Table";
    }

    // interface ExtendedTableReport
    public int getColumnCount() {
        return 4;
    }

    public long getRowCount() {
        return idSynonyms.size();
    }

    public long getElementCount() {
        return getColumnCount() * getRowCount();
    }
}
