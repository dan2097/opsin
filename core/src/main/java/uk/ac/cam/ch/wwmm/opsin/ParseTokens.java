package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.List;
/**A "struct" containing data a possible tokenisation of a word in a chemical name.
 *
 * @author ptc24
 *
 */
class ParseTokens {
	/**The tokens that the word is made up of.*/
	List<String> tokens;
	/**A list of possible annotations of that token.*/
	List<Character> annotations;

	ParseTokens deepCopy() {
		ParseTokens pt = new ParseTokens();
		pt.tokens = new ArrayList<String>();
		for(String t : tokens) pt.tokens.add(t);
		pt.annotations = new ArrayList<Character>();
		for(Character c : annotations) {
			pt.annotations.add(c);
		}
		return pt;
	}

	public String toString() {
		return "[" + tokens + ", " + annotations + "]";
	}

}
