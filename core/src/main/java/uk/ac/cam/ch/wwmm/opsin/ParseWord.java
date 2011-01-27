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

	static enum WordType{
		full,
		substituent,
		functionalTerm,
	}

	/**The type of word - full, substituent or functionalGroup*/
	private WordType wordType = null;

	ParseWord deepCopy() {
		ParseWord pw = new ParseWord(word, parseTokens);
		pw.wordType = wordType;
		return pw;
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
	
	WordType getWordType() {
		return wordType;
	}

	void setWordType(WordType wordType) {
		this.wordType = wordType;
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
		return "[" + word + ", " + wordType + ", " + parseTokens.toString() + "]";
	}

}
