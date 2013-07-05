package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.List;

/**
 * A black/white radix tree implementation.
 * A radix tree is a type of trie where common prefixes are merged together to save space 
 * This implementation employs short arrays rather than maps to exploit the fact that all OPSIN tokens are ASCII.
 * @author dl387
 *
 */
class OpsinRadixTrie {
	final OpsinTrieNode rootNode;

	OpsinRadixTrie() {
		rootNode = new OpsinTrieNode("", false);
	}

	/**
	 * Adds a string to the Trie.
	 * This string should not contain any non ASCII characters
	 * @param token
	 */
	void addToken(String token) {
		int tokenLength =token.length();
		String remaingStr =token;
		OpsinTrieNode currentNode = rootNode;
		for (int i = 0; i < tokenLength;) {
			int charsMatched = currentNode.getNumberOfMatchingCharacters(remaingStr, 0);
			remaingStr = remaingStr.substring(charsMatched);
			i+=charsMatched;
			currentNode = currentNode.add(remaingStr, charsMatched);
		}
		currentNode.setIsEndPoint(true);
	}

	/**
	 * Returns all possible runs of the input string that reached end point nodes in the trie
	 * e.g. ylidene might return 2 ("yl"), 6 ("yliden") and 7 ("ylidene")
	 * Results are given as the index of the end of the match in the chemicalName
	 * Returns null if no runs were possible
	 * @param posInName The point at which to start matching
	 * @param chemicalName
	 * @return
	 */
	List<Integer> findMatches(String chemicalName, int posInName) {
		int untokenisedChemicalNameLength = chemicalName.length();
		List<Integer> lengths = null;
		if (rootNode.isEndPoint()){
			lengths = new ArrayList<Integer>();
			lengths.add(posInName);
		}
		OpsinTrieNode currentNode = rootNode;
		for (int i = posInName; i < untokenisedChemicalNameLength; i++) {
			OpsinTrieNode node = currentNode.getChild(chemicalName.charAt(i));
			if (node != null) {
				currentNode = node;
				int charsMatched = currentNode.getNumberOfMatchingCharacters(chemicalName, i);
				i += (charsMatched - 1);
				if (charsMatched == currentNode.getValue().length()){
					if (currentNode.isEndPoint()) {
						if (lengths == null) {
							lengths = new ArrayList<Integer>();
						}
						lengths.add(i + 1);
					}
				}
				else{
					break;
				}
			} else {
				break;
			}
		}
		return lengths;
	}

	/**
	 * Same as findLengthsOfMatches but the trie has been populated by reversed tokens
	 * @param untokenisedChemicalName
	 * @return
	 */
	List<Integer> findLengthsOfMatchesReadingStringRightToLeft(String untokenisedChemicalName) {
		int untokenisedChemicalNameLength = untokenisedChemicalName.length();
		List<Integer> lengths = null;
		if (rootNode.isEndPoint()){
			lengths = new ArrayList<Integer>();
			lengths.add(0);
		}
		OpsinTrieNode currentNode = rootNode;
		for (int i = untokenisedChemicalNameLength-1; i >=0; i--) {
			OpsinTrieNode node = currentNode.getChild(untokenisedChemicalName.charAt(i));
			if (node != null) {
				currentNode = node;
				int charsMatched =currentNode.getNumberOfMatchingCharactersInReverse(untokenisedChemicalName.substring(0, i + 1));
				i-=(charsMatched-1);
				if (charsMatched == currentNode.getValue().length()){
					if (currentNode.isEndPoint()) {
						if (lengths == null) {
							lengths = new ArrayList<Integer>();
						}
						lengths.add(untokenisedChemicalNameLength -i);
					}
				}
				else{
					break;
				}
			} else {
				break;
			}
		}
		return lengths;
	}
}

class OpsinTrieNode {

	private boolean isEndPoint;
	private String key;
	private OpsinTrieNode[] children = new OpsinTrieNode[128];

	OpsinTrieNode(String key, boolean isEndPoint) {
		this.isEndPoint = isEndPoint;
		this.key = key;
	}

	String getValue() {
		return key;
	}
	
	boolean isEndPoint() {
		return isEndPoint;
	}

	void setIsEndPoint(boolean isEndPoint) {
		this.isEndPoint = isEndPoint;
	}
	
	private void setChildren(OpsinTrieNode[] children) {
		this.children = children;
	}
	
	OpsinTrieNode add(String remaingStr, int charsMatched) {
		if (charsMatched < key.length()){//need to split this Trie node
			OpsinTrieNode newNode = new OpsinTrieNode(key.substring(charsMatched), isEndPoint);
			newNode.setChildren(children);
			children = new OpsinTrieNode[128];
			children[key.charAt(charsMatched)] = newNode;
			key = key.substring(0, charsMatched);
			isEndPoint =false;
		}
		if (remaingStr.length()!=0){
			int charValue = (int) remaingStr.charAt(0);
			if (children[charValue] == null) {
				children[charValue] = new OpsinTrieNode(remaingStr, false);
			}
			return children[charValue];
		}
		return this;
	}

	int getNumberOfMatchingCharacters(String chemicalName, int posInName) {
		int maxLength = Math.min(key.length(), chemicalName.length() - posInName);
		for (int i = 0; i < maxLength; i++) {
			if (key.charAt(i) != chemicalName.charAt(posInName + i)){
				return i;
			}
		}
		return maxLength;
	}
	
	int getNumberOfMatchingCharactersInReverse(String token) {
		int tokenLength = token.length();
		int maxLength = Math.min(key.length(), tokenLength);
		for (int i = 0; i < maxLength; i++) {
			if (key.charAt(i)!=token.charAt(tokenLength-i-1)){
				return i;
			}
		}
		return maxLength;
	}

	OpsinTrieNode getChild(char c) {
		return children[(int) c];
	}
}
