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
	List<ParseWord> words;
	
	Parse deepCopy() {
		Parse p = new Parse();
		p.name = name;
		p.wordRule = wordRule;
		p.words = new ArrayList<ParseWord>();
		for(ParseWord pw : words) {
			p.words.add(pw.deepCopy());
		}
		return p;
	}
	
	public String toString() {
		return "[" + name + ", " + wordRule + ", " + words.toString() + "]";
	}
}
