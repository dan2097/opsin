package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.List;

/**A "struct" containing data on the parsing of a chemical name.
 *
 * @author ptc24/dl387
 *
 */
class Parse {
	/**The chemical name.*/
	private final String name;
	/**The words within the name, and their parsing data.*/
	private final List<ParseWord> words = new ArrayList<ParseWord>();

	/**
	 * Creates a parse object for a chemicalName
	 * @param chemicalName
	 */
	Parse(String chemicalName) {
		name = chemicalName;
	}

	Parse deepCopy() {
		Parse p = new Parse(name);
		for(ParseWord pw : words) {
			p.words.add(pw.deepCopy());
		}
		return p;
	}

	public String toString() {
		return "[" + name + ", " +  words.toString() + "]";
	}

	List<ParseWord> getWords() {
		return words;
	}

	boolean addWord(ParseWord word) {
		return words.add(word);
	}
	
	boolean removeWord(ParseWord word) {
		return words.remove(word);
	}
	
	ParseWord getWord(int indice) {
		return words.get(indice);
	}

	String getName() {
		return name;
	}

}
