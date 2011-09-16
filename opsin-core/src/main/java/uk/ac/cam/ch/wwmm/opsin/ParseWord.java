package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.List;
/**A "struct" containing data on the parsing of a word in a chemical name.
 *
 * @author ptc24
 * @author dl387
 *
 */
class ParseWord {
	/**The word itself.*/
	private final String word;
	/**All of the possible tokenisations of the word.*/
	private List<ParseTokens> parseTokens;

	ParseWord deepCopy() {
		return new ParseWord(word, parseTokens);
	}
	
	ParseWord(String word, List<ParseTokens> parseTokens) {
		this.word =word;
		if (parseTokens ==null){
			this.parseTokens =null;
		}
		else{
			this.parseTokens =  new ArrayList<ParseTokens>(parseTokens);
		}
	}

	String getWord() {
		return word;
	}

	List<ParseTokens> getParseTokens() {
		return parseTokens;
	}
	
	void setParseTokens(List<ParseTokens> parseTokens) {
		this.parseTokens = parseTokens;
	}

	public String toString() {
		return "[" + word + ", " + parseTokens.toString() + "]";
	}

}
