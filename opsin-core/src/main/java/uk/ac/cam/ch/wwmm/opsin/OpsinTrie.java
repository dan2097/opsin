package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.List;

class OpsinTrie {
	private final OpsinTrieNode rootNode;

	OpsinTrie() {
		rootNode = new OpsinTrieNode(' ');
	}

	void addToken(String token) {
		OpsinTrieNode currentNode = rootNode;
		for (char c : token.toCharArray()) {
			currentNode = currentNode.add(c);
		}
		currentNode.add('$');//use a dollar to indicate the end of the token
	}

	List<Integer> findLengthsOfMatches(String untokenisedChemicalName, int untokenisedChemicalNameLength) {
		List<Integer> lengths = null;
		OpsinTrieNode currentNode = rootNode;
		for (int i = 0; i < untokenisedChemicalNameLength; i++) {
			OpsinTrieNode node = currentNode.getChild(untokenisedChemicalName.charAt(i));
			if (node != null) {
				currentNode = node;
				if (currentNode.getChild('$') != null) {
					if (lengths == null) {
						lengths = new ArrayList<Integer>();
					}
					lengths.add(i + 1);
				}
			} else {
				break;
			}
		}
		return lengths;
	}

	List<Integer> findLengthsOfMatchesReadingStringRightToLeft(String untokenisedChemicalName, int untokenisedChemicalNameLength) {
		List<Integer> lengths = null;
		OpsinTrieNode currentNode = rootNode;
		for (int i = untokenisedChemicalNameLength-1; i >=0; i--) {
			OpsinTrieNode node = currentNode.getChild(untokenisedChemicalName.charAt(i));
			if (node != null) {
				currentNode = node;
				if (currentNode.getChild('$') != null) {
					if (lengths == null) {
						lengths = new ArrayList<Integer>();
					}
					lengths.add(untokenisedChemicalNameLength -i);
				}
			} else {
				break;
			}
		}
		return lengths;
	}
}

class OpsinTrieNode {

	private final char character;
	private final OpsinTrieNode[] children = new OpsinTrieNode[128];

	OpsinTrieNode(char c) {
		character = c;
	}

	char getValue() {
		return character;
	}

	OpsinTrieNode add(char c) {
		int charValue = (int) c;
		if (children[charValue] == null) {
			children[charValue] = new OpsinTrieNode(c);
		}
		return children[charValue];
	}

	OpsinTrieNode getChild(char c) {
		return children[(int) c];
	}
}
