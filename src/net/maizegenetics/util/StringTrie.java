package net.maizegenetics.util;

/**
 * This is an implementation of a trie (prefix tree) in Java.
 * Supports opeations like searching a string, searching a prefix, searching by prefix etc.
 * @author Janu Verma
 * jv367@cornell.edu
 **/


import net.maizegenetics.analysis.gbs.Barcode;

import java.lang.Character;
import java.lang.Object;
import java.lang.String;
import java.lang.System;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Collection;

public class StringTrie{

    private TrieNode root;



    /** Constructor */

    public StringTrie()
    {
        root = new TrieNode();
    }



    /** Adds a word to the trie
     * @param word
     */
    public void addWord(String word){
        root.addWord(word.toLowerCase());
    }


    /**
     * Add a collection of words to the trie.
     * @param words
     */

    public void addAllWords(Collection<String> words){
        for (String w:words)
            addWord(w);
    }


    /**
     * checks if the String is in the trie.
     * @param s
     * @return true if the string is in the trie
     */
    public boolean contains(String s){
        TrieNode currentNode = root;
        for (int i = 0; i < s.length(); i++){
            char c = s.charAt(i);
            if (currentNode.containsKey(c))
                currentNode = currentNode.getNode(c);
            else
                return false;
        }
        return true;
    }



    /**
     * get the words in the trie with the given
     * prefix
     * @param prefix
     * @return a List contaning String objects containing the words
     * in the Trie with the given prefix.
     */
    public List getWords(String prefix){
        // Find the node which represents the last letter of the prefix.
        TrieNode lastNode = root;
        for (int i = 0; i < prefix.length(); i++){
            lastNode = lastNode.getNode(prefix.charAt(i));

            // If no node matches, then no words exits, return empty list
            if (lastNode == null) return new ArrayList();
        }
        // Return the words which eminate from the last node
        return lastNode.getWords();
    }


    /**
     * Find the longest prefix of a string
     * @param input
     */
    public String longestPrefix(String input){
        String result = "";
        int length = input.length();
        TrieNode crawl = root;
        int level, prevMatch = 0;
        for (level = 0; level < length-1; level++){
            char ch = input.charAt(level);
            TrieNode child = crawl.getNode(ch);  //Get the Node reprsenting the character.
            if (crawl.containsKey(ch)){
                result += ch;
                crawl = child;
                if (crawl.isWord)
                    prevMatch = level + 1;
            }
            else break;
        }
        if (!crawl.isWord)
            return result.substring(0,prevMatch);

        else return result;
    }





    public static void main(String args[]){
        List info;
        String perid;
        TrieNode t;
        List child;
        StringTrie tau;
        char ch = 'a';
        StringTrie trie = new StringTrie();
        trie.addWord("alpha");
        trie.addWord("beta");
        trie.addWord("altare");
        trie.addWord("altarations");
        trie.addWord("bet");
        trie.addWord("al");
        info = trie.getWords("al");
        // System.out.print(info + "\n");
        t = trie.root;
        child = t.getWords();
        Object alpha = child.get(0);
        // System.out.println(alpha);
        //tau = trie.getNode(ch);
        //child = t.children();
        String yaay = "altarationsX";
        String barcode1 ="ACGACAACGACG";
        String barcode2 ="ACCAACGACG";

        String query="ACGACAACGACGACTGATCGATCGATGTACGATCGATCG";

        perid = trie.longestPrefix(yaay);
        System.out.println(perid);
        //boolean gaga = t.isLeaf;
        //System.out.print("gaga: " + gaga);
    }
}

class TrieNode{
    public TrieNode parent;
    public TrieNode[] children;
    public boolean isLeaf; // Quick way to check if any children exist
    public boolean isWord; // does this node represent the last character
    private char character; //character the node represents
    private Barcode theBarcode;


    /**
     * Constructor for top level root node.
     */
    public TrieNode()
    {
        children = new TrieNode[26];
        isLeaf = true;
        isWord = false;
    }


    /**
     * Constructor for the child node.
     */
    public TrieNode(char character){
        this();
        this.character = character;
    }



    /**
     * Adds a word to this node. This method is called recursively and
     * adds child nodes for each successive letter in the word,
     * therefore recursive calls will be made with partial words.
     * @param word - the word to add
     */
    protected void addWord(String word){
        isLeaf = false;
        int charPos = word.charAt(0) - 'a';
        if (children[charPos] == null){
            children[charPos] = new TrieNode(word.charAt(0));
            children[charPos].parent = this;
        }
        if (word.length() > 1){
            children[charPos].addWord(word.substring(1));
        }
        else{
            children[charPos].isWord = true;
        }

    }




    /**
     * Return the child TrieNode representing rthe given char,
     * or null if no node exists.
     * @param c
     * @return TrieNode
     */
    protected TrieNode getNode(char c){
        return children[c-'a'];
    }



    /**
     * checks if the given character is a children.
     * @param c
     * @return true if c is a children
     */
    public boolean containsKey(char c){
        List followers = new ArrayList();
        for (TrieNode x : children) {
            if (x != null) {
                char y = x.character;
                followers.add(y);
            }
        }
        return (followers.contains(c));
    }



    /**
     * Returns a List of String objects which are lower in the
     * hierarchy than this node.
     * @return List of words
     */
    protected List getWords() {
        //Create a list to return.
        List list = new ArrayList();

        // If this node represents a word, add it.
        if (isWord) {
            list.add(toString());
        }
        // if any children
        if (!isLeaf) {
            //Add any words belonging to any children
            for (int i = 0; i < children.length; i++) {
                if (children[i] != null) {
                    list.addAll(children[i].getWords());
                }
            }
        }
        return list;
    }





    /**
     * Gets the string that this node represents.
     * e.g.g if this node represents the charcter t, whose parent
     * represents the character a, whose parent represents the charcter
     * c, then the string would be cat.
     * @return String
     */

    public String toString(){
        if (parent == null){
            return "";
        }
        else{
            return parent.toString() + new String(new char[]{character});
        }
    }
}
