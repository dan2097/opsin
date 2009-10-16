package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import nu.xom.Elements;
import dk.brics.automaton.RunAutomaton;

/**Performs finite-state allocation of roles ("annotations") to tokens:
 * The chemical name is broken down into tokens e.g. ethyl -->eth yl by applying the chemical grammar in regexes.xml
 * The tokens eth and yl are associated with a letter which is referred to here as an annotation which is the role of the token.
 * These letters are defined in regexes.xml and would in this case have the meaning alkaneStem and inlineSuffix
 *
 * The chemical grammer employs the annoatations associated with the tokens when deciding what may follow what has already been seen
 * e.g. you cannot start a chemical name with yl and an optional e is valid after an arylGroup
 *
 * @author ptc24/dl387
 *
 */
class ParseRules {

	/** A "struct" containing bits of state needed during finite-state parsing. */
	private class AnnotatorState {
		/** The current state of the DFA. */
		int state;
		/** The annotation so far. */
		List<Character> annot;
		/** The strings these annotations correspond to. */
		ArrayList<String> tokens;
	}

	/**
	 * Wrapper class for returning multiple objects
	 */
	final class TwoReturnValues {
		final List<ParseTokens> first;
		final boolean second;

		public TwoReturnValues(List<ParseTokens> first, boolean second) {
			this.first = first;
			this.second = second;
		}

		public List<ParseTokens> getFirst() {
			return first;
		}

		public boolean getSecond() {
			return second;
		}
	}

	/** Maps regular expression names to regular expressions. */
	HashMap<String, String> regexDict;
	/** A DFA for "root" words in chemical names. */
	private RunAutomaton chemAutomaton;
	/** A DFA for "substituent" words in chemical names. */
	private RunAutomaton subAutomaton;

	private final static char endOfSubstituent = '\u00e9';
	private final static char endOfMainGroup = '\u00e2';
	private final static char interSubstituentHyphen = '\u00e4';

	private TokenManager tokenManager;
	private ResourceGetter resourceGetter;

	//private String largestLex="";//holds the smallest amount of unlexed/unlexable text so far

	/** Initialises the finite-state parser, reading in the rules from regexes.xml.
	 * @param tokenManager
	 * @param resourceGetter 
	 *
	 * @throws Exception If the rules file can't be read properly.
	 */
	ParseRules(TokenManager tokenManager, ResourceGetter resourceGetter) throws Exception {
		this.tokenManager = tokenManager;
		this.resourceGetter = resourceGetter;
		regexDict = new HashMap<String, String>();
		Elements regexes = resourceGetter.getXMLDocument("regexes.xml").getRootElement().getChildElements("regex");
		Pattern p = Pattern.compile("%.*?%");
		for(int i=0;i<regexes.size();i++) {
			String name = regexes.get(i).getAttributeValue("name");
			String value = regexes.get(i).getAttributeValue("value");
			Matcher m = p.matcher(value);
			String newValue = "";
			int position = 0;
			while(m.find()) {
				newValue += value.substring(position, m.start());
				if (regexDict.get(m.group())==null){
					throw new ParsingException("Regex entry for: " + m.group() + " missing! Check regexes.xml");
				}
				newValue += regexDict.get(m.group());
				position = m.end();
			}
			newValue += value.substring(position);
			regexDict.put(name, newValue);
		}

		chemAutomaton = getChemicalAutomaton();
		subAutomaton = getSubstituentAutomaton();
	}

	/** Compiles the DFA for a "root" chemical name word.
	 *
	 * @return The DFA for a "root" chemical name word.
	 */
	RunAutomaton getChemicalAutomaton() {
		String re = regexDict.get("%chemical%");
		return AutomatonInitialiser.getAutomaton("chemical", re, resourceGetter);
	}

	/** Compiles the DFA for a "substituent" chemical name word.
	 *
	 * @return The DFA for a "substituent" chemical name word.
	 */
	RunAutomaton getSubstituentAutomaton() {
		String re = regexDict.get("%substituent%");
		return AutomatonInitialiser.getAutomaton("substituent", re+"*", resourceGetter);
	}

	/** Compiles the DFA for a chemical name word.
	 *
	 * @param wr The word rule for the word.
	 * @return The DFA for the chemical name word.
	 */
	RunAutomaton getAutomatonForWordRule(String wr) {
		if(wr.equals("substituent")) return subAutomaton;
		return chemAutomaton;
	}

	/**Determines the possible annotations for a chemical word.
	 *
	 * @param possibleAnnotations A list of (list of possible annotations)s for each token in the word
	 * @param wordRule The word rule for the word.
	 * @return A list of possible annotations for the word.
	 */
	TwoReturnValues getParses(String chemicalWord, String wordRule) {
		//largestLex=chemicalWord;
		RunAutomaton automaton = getAutomatonForWordRule(wordRule);
		AnnotatorState as = new AnnotatorState();
		as.state = automaton.getInitialState();
		as.annot = new ArrayList<Character>();
		as.tokens = new ArrayList<String>();
		List<AnnotatorState> states = new ArrayList<AnnotatorState>();
		moveToNextAnnotation(chemicalWord, chemicalWord.toLowerCase(), as, states, automaton, automaton.getCharIntervals());
		List<ParseTokens> outputList = new ArrayList<ParseTokens>();
		boolean interSubHyphenAllowed = false;//false if name is not lexable or if none of the partial lexes allow such a hyphen as their next token
		for(AnnotatorState aas : states) {
			if(automaton.isAccept(aas.state)) {
				ParseTokens pt =new ParseTokens();
				pt.tokens=aas.tokens;
				pt.annotations=aas.annot;
				outputList.add(pt);
			}
			else{
				if (automaton.step(aas.state, interSubstituentHyphen) != -1){
					interSubHyphenAllowed=true;
				}
			}
		}
		if (outputList.size()==0){
			//System.out.println(largestLex);
		}
		TwoReturnValues output =new TwoReturnValues(outputList, interSubHyphenAllowed);
		return output;
	}

	/**
	 * Recursively attempts to find all valid annotations for the given chemical name
	 * This is done by a depth first search
	 * Successful annotations are stored in successfulAnnotations
	 * @param chemicalWord The chemical name with any part of the name that has been lexed already removed from the front
	 * @param chemicalWordLowerCase The same chemical name but lower case, used for token matching but not regex matching
	 * @param as An AnnotatorState, this contains the current state of the automaton.
	 * @param succesfulAnnotations
	 */
	void moveToNextAnnotation(String chemicalWord, String chemicalWordLowerCase, AnnotatorState as, List<AnnotatorState> succesfulAnnotations, RunAutomaton automaton, char[] stateSymbols ){
		int wordLength = chemicalWordLowerCase.length();
		if (wordLength ==0){
			//attempt to manually step to the endOfMainGroup/endOfSubstituent state if allowed
			if (automaton == chemAutomaton){
				if (automaton.step(as.state, endOfMainGroup)!=-1){
					as.annot.add(endOfMainGroup);
					as.state=automaton.step(as.state, endOfMainGroup);
				}
			}
			else if (automaton == subAutomaton){
				if (automaton.step(as.state, endOfSubstituent)!=-1){
					as.annot.add(endOfSubstituent);
					as.state=automaton.step(as.state, endOfSubstituent);
				}
			}
			succesfulAnnotations.add(as);
		}
		String firstTwoLetters =null;
		if (wordLength >=2){
			firstTwoLetters=chemicalWordLowerCase.substring(0,2);
		}

		for (int j = 0; j < stateSymbols.length; j++) {
			char annotationCharacter =stateSymbols[j];
			int potentialNextState = automaton.step(as.state, annotationCharacter);
			if (potentialNextState != -1) {//-1 means this state is not accessible from the previous state
				HashMap<String, List<String>> possibleTokenisationsMap = tokenManager.symbolTokenNamesDict.get(annotationCharacter);
				if (possibleTokenisationsMap!=null){
					List<String> possibleTokenisations =null;
					if (firstTwoLetters!=null){
						possibleTokenisations= possibleTokenisationsMap.get(firstTwoLetters);
					}
					if (possibleTokenisations!=null){//next could be a token
						for (String possibleTokenisation : possibleTokenisations) {
							if (chemicalWordLowerCase.startsWith(possibleTokenisation)){
								AnnotatorState newAs =new AnnotatorState();
								String newchemicalWord=chemicalWord.substring(possibleTokenisation.length());
								String newchemicalWordLowerCase=chemicalWordLowerCase.substring(possibleTokenisation.length());
								newAs.tokens =new ArrayList<String>(as.tokens);
								newAs.tokens.add(possibleTokenisation);
								newAs.annot = new ArrayList<Character>(as.annot);
								newAs.annot.add(annotationCharacter);
								newAs.state=potentialNextState;
								//System.out.println("tokened " +newchemicalWord);
								moveToNextAnnotation(newchemicalWord, newchemicalWordLowerCase, newAs, succesfulAnnotations, automaton, stateSymbols);
							}
						}
					}
				}
				List<Pattern> possibleRegexes =tokenManager.symbolRegexesDict.get(annotationCharacter);
				if (possibleRegexes!=null){//next could be a regex
					for (Pattern pattern : possibleRegexes) {
						Matcher mat =pattern.matcher(chemicalWord);
						if (mat.lookingAt()){//match at start
							AnnotatorState newAs =new AnnotatorState();
							String newchemicalWord=chemicalWord.substring(mat.group(0).length());
							String newchemicalWordLowerCase=chemicalWordLowerCase.substring(mat.group(0).length());
							newAs.tokens =new ArrayList<String>(as.tokens);
							if (!mat.group(0).equals("")){//special case for endOfSubstituent which we do not want to form a token from
								newAs.tokens.add(mat.group(0));
							}
							newAs.annot = new ArrayList<Character>(as.annot);
							newAs.annot.add(annotationCharacter);
							newAs.state=potentialNextState;
							//System.out.println("neword regex " +newchemicalWord);
							moveToNextAnnotation(newchemicalWord, newchemicalWordLowerCase, newAs, succesfulAnnotations, automaton, stateSymbols);
						}
					}
				}
			}
		}
//		if (chemicalWord.length() < largestLex.length()){
//			largestLex=chemicalWord;
//		}
	}

	/**Groups the token annotations for a given word into substituent/s and/or a maingroup by looking for the endOfSubstituent or endOfMainGroup annotation
	 *
	 * @param annots The annotation for a word.
	 * @return A List of lists of annotations, each list corresponds to a substituent/maingroup
	 */
	List<List<Character>> chunkAnnotations(List<Character> annots) {
		LinkedList<List<Character>> chunkList = new LinkedList<List<Character>>();
		List<Character> currentSubOrRoot = new ArrayList<Character>();
		for (Character annot : annots) {
			if (annot != endOfSubstituent && annot !=endOfMainGroup){
				currentSubOrRoot.add(annot);
			}
			else{
				chunkList.add(currentSubOrRoot);
				currentSubOrRoot = new ArrayList<Character>();
			}
		}
		return chunkList;
	}
}
