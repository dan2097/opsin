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
	private final List<ParseTokens> parseTokens;
	/**Whether an explicit mixture delimiter follows this word.*/
	private boolean endsComponent;

	ParseWord deepCopy() {
		return new ParseWord(word, parseTokens, endsComponent);
	}

	ParseWord(String word, List<ParseTokens> parseTokens) {
		this(word, parseTokens, false);
	}
	
	ParseWord(String word, List<ParseTokens> parseTokens, boolean endsComponent) {
		this.word = word;
		this.parseTokens = new ArrayList<>(parseTokens);
		this.endsComponent = endsComponent;
	}

	String getWord() {
		return word;
	}

	List<ParseTokens> getParseTokens() {
		return parseTokens;
	}
	
	boolean endsComponent() {
		return endsComponent;
	}

	void setEndsComponent(boolean endsComponent) {
		this.endsComponent = endsComponent;
	}

	public String toString() {
		return "[" + word + ", " + parseTokens.toString() + "]";
	}

}
