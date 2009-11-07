package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.List;

/**A "struct" containing data on the parsing of a chemical name.
 * 
 * @author ptc24
 *
 */
class Parse {
	/**The chemical name.*/
	String name;
	/**The rule detailing the arrangement of words within the name.*/
	String wordRule;
	/**The words within the name, and their parsing data.*/
	List<ParseWord> words = new ArrayList<ParseWord>();
	
	/**
	 * Creates a parse object for a chemicalName
	 * @param chemicalName
	 */
	public Parse(String chemicalName) {
		this.name = chemicalName;
	}

	Parse deepCopy() {
		Parse p = new Parse(name);
		p.wordRule = wordRule;
		for(ParseWord pw : words) {
			p.words.add(pw.deepCopy());
		}
		return p;
	}
	
	public String toString() {
		return "[" + name + ", " + wordRule + ", " + words.toString() + "]";
	}
}
