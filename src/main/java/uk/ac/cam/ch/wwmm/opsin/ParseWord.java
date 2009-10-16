package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.List;
/**A "struct" containing data on the parsing of a word in a chemical name.
 * 
 * @author ptc24
 *
 */
class ParseWord {
	/**The word itself.*/
	String word;
	/**The type of word - a literal, a "root" word or a "substituent" word.*/
	String wordType;
	/**All of the possible tokenisations of the word.*/
	List<ParseTokens> parseTokens;
	
	ParseWord deepCopy() {
		ParseWord pw = new ParseWord();
		pw.word = word;
		pw.wordType = wordType;
		if(parseTokens == null) {
			pw.parseTokens = null;
		} else {
			pw.parseTokens = new ArrayList<ParseTokens>();
			for(ParseTokens pt : parseTokens) {
				pw.parseTokens.add(pt);
			}
		}
		return pw;
	}
	
	public String toString() {
		return "[" + word + ", " + wordType + ", " + parseTokens.toString() + "]";
	}

}
