package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.List;
/**A "struct" containing data a possible tokenisation of a word in a chemical name.
 *
 * @author ptc24/dl387
 *
 */
public class ParseTokens {
	/**The tokens that the word is made up of.*/
	private List<String> tokens = new ArrayList<String>();

	/**A list of possible annotations of that token.*/
	private List<Character> annotations = new ArrayList<Character>();

	ParseTokens deepCopy() {
		ParseTokens pt = new ParseTokens();
		pt.tokens = new ArrayList<String>(tokens);
		pt.annotations = new ArrayList<Character>(annotations);
		return pt;
	}
	
	/**
	 * Creates an empty parseTokens
	 */
	ParseTokens() {
	}
	
	/**
	 * Creates a parseTokens from an existing list of tokens and annotations
	 * The lists should be of identical lengths otherwise an exception is thrown
	 * @throws ParsingException 
	 */
	ParseTokens(List<String> tokens, List<Character> annotations ) throws ParsingException {
		if (tokens.size() != annotations.size()){
			throw new ParsingException("OPSIN bug: mismatch between the sizes of tokens list and annotation list");
		}
		this.tokens = tokens;
		this.annotations = annotations;
	}
	
	void addToken(String tokenString, char annotationSymbol){
		tokens.add(tokenString);
		annotations.add(annotationSymbol);
	}

	public List<String> getTokens() {
		return tokens;
	}

	public List<Character> getAnnotations() {
		return annotations;
	}
	
	public String toString() {
		return "[" + tokens + ", " + annotations + "]";
	}

}
