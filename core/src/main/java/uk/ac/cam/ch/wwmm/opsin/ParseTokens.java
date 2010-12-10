package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
/**A "struct" containing data a possible tokenisation of a word in a chemical name.
 *
 * @author ptc24/dl387
 *
 */
public class ParseTokens {
	/**The tokens that the word is made up of.*/
	private final List<String> tokens;

	/**A list of possible annotations of that token.*/
	private final List<Character> annotations;

	
	/**
	 * Creates a parseTokens from an existing list of tokens and annotations
	 * The lists should be of identical lengths otherwise an exception is thrown
     * @param tokens
     * @param annotations
	 */
	ParseTokens(List<String> tokens, List<Character> annotations ){
		if (tokens.size() != annotations.size()){
			throw new IllegalArgumentException("OPSIN bug: mismatch between the sizes of tokens list and annotation list");
		}
		this.tokens = Collections.unmodifiableList(new ArrayList<String>(tokens));
		this.annotations = Collections.unmodifiableList(new ArrayList<Character>(annotations));
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
	
	@Override
	public boolean equals(Object other) {
		if (this == other) {
			return true;
		}
		if (other instanceof ParseTokens) {
			ParseTokens otherPT = (ParseTokens) other;
			return this.tokens.equals(otherPT.tokens) && this.annotations.equals(otherPT.annotations);
		}
		return false;
	}
	
	@Override
	public int hashCode() {		
		return (3 * this.tokens.hashCode()) * (7 * this.annotations.hashCode());
	}
}
