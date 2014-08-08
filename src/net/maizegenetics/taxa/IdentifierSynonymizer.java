package net.maizegenetics.taxa;


import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.Report;
import net.maizegenetics.util.AbstractTableReport;

import java.io.PrintWriter;
import java.io.Serializable;
import java.util.*;


/**
 * User: Ed
 * Date: Mar 30, 2005
 * Time: 1:39:47 PM
 */
public class IdentifierSynonymizer extends AbstractTableReport implements Serializable, Report, TableReport {

    HashMap<String,Integer> idSynonyms = new HashMap<>();    //TODO needs to be entirely updated to new collections
    private TaxaList referenceIDGroup;
    private int unmatchCount = 0;

    public IdentifierSynonymizer(TaxaList preferredTaxa, TaxaList[] alternateTaxaSets) {
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
            }
        }
    }

    private ArrayList<String> findBestMatch(String unmatchedString) {
        ArrayList<String> bestMatches = new ArrayList<>();
        double maxScore = -1;
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
            for (int i = 0; i < referenceIDGroup.numberOfTaxa(); i++) {
                sm = scoreMatch(referenceIDGroup.taxaName(i), unmatchedString, ignoreCase, ignoreWhite, ignorePunc);
                if (sm > maxScore) {
                    bestMatches.clear();
                    bestMatches.add(referenceIDGroup.taxaName(i));
                    maxScore = sm;
                } else if (sm == maxScore) {
                    bestMatches.add(referenceIDGroup.taxaName(i));
                }
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
            sm = scoreMatch(referenceIDGroup.taxaName(i), unmatchedString, ignoreCase, ignoreWhite, ignorePunc);
            theSortMap.put(1 - sm - ((double) i / 100000.0), referenceIDGroup.taxaName(i));
        }
        return new ArrayList<>(theSortMap.values());
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
        data[3] = "" + scoreMatch("" + data[0], "" + data[1], true, false, false);
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
