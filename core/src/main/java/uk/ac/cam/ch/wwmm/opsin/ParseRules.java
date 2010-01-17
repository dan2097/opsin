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
		/** After an accept state has been reached the remaining chemical name is placed here
		 * Usually only annotator states in which untokenisedName ="" need be considered */
		String untokenisedName = null;
	}
	
	/**
	 * A wrapper for an integer with a getter/setter
	 * @author dl387
	 *
	 */
	private class MutableInteger{
		int value;

		int getValue() {
			return value;
		}

		void setValue(int i) {
			value = i;
		}
		public MutableInteger(int i) {
			value =i;
		}
	}


	/** Maps regular expression names to regular expressions. */
	HashMap<String, String> regexDict;
	/** A DFA encompassing the grammar of a chemical word. */
	private RunAutomaton chemAutomaton;
	/** The allowed symbols in chemAutomaton */
	private char[] stateSymbols;

	private final static char endOfSubstituent = '\u00e9';
	private final static char endOfMainGroup = '\u00e2';
	private final static char endOfFunctionalTerm = '\u00FB';
	//private final static char interSubstituentHyphen = '\u00e4';

	private TokenManager tokenManager;
	private ResourceGetter resourceGetter;

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
		stateSymbols = chemAutomaton.getCharIntervals();
	}

	/** Compiles the DFA for a "root" chemical name word.
	 *
	 * @return The DFA for a "root" chemical name word.
	 */
	RunAutomaton getChemicalAutomaton() {
		String re = regexDict.get("%chemical%");
		return AutomatonInitialiser.getAutomaton("chemical", re, resourceGetter);
	}

	/**Determines the possible annotations for a chemical word
	 * Returns a list of parses and how much of the word could not be interpreted
	 * e.g. usually the list will have only one parse and the string will equal ""
	 * For something like ethyloxime. The list will contain the parse for ethyl and the string will equal "oxime" as it was unparsable
	 * For something like eth no parses would be found and the string will equal "eth"
	 *  
	 * @param chemicalWord
	 * @return
	 * @throws ParsingException
	 */
	TwoReturnValues<List<ParseTokens>, String> getParses(String chemicalWord) throws ParsingException {
		AnnotatorState as = new AnnotatorState();
		as.state = chemAutomaton.getInitialState();
		as.annot = new ArrayList<Character>();
		as.tokens = new ArrayList<String>();
		List<AnnotatorState> states = new ArrayList<AnnotatorState>();
		moveToNextAnnotation(chemicalWord, chemicalWord.toLowerCase(), as, states, new MutableInteger(chemicalWord.length()));
		List<ParseTokens> outputList = new ArrayList<ParseTokens>();
		String unparseableName = chemicalWord;
		if (states.size()>0){//at least some of the name could be interpreted into a substituent/full/functionalTerm
			for(AnnotatorState aas : states) {
				ParseTokens pt =new ParseTokens(aas.tokens, aas.annot);
				outputList.add(pt);
				unparseableName=aas.untokenisedName;//all acceptable annotator states found should have the same untokenisedName
			}
		}
		return new TwoReturnValues<List<ParseTokens>, String>(outputList, unparseableName);
	}

	/**
	 * Recursively attempts to find all valid annotations for the given chemical name
	 * This is done by a depth first search
	 * Successful annotations are stored in successfulAnnotations
	 * @param chemicalWord The chemical name with any part of the name that has been lexed already removed from the front
	 * @param chemicalWordLowerCase The same chemical name but lower case, used for token matching but not regex matching
	 * @param as An AnnotatorState, this contains the current state of the automaton.
	 * @param successfulAnnotations
	 */
	void moveToNextAnnotation(String chemicalWord, String chemicalWordLowerCase, AnnotatorState as, List<AnnotatorState> successfulAnnotations, MutableInteger wordLengthRemainingOnLastSuccessfulAnnotations){
		int wordLength = chemicalWordLowerCase.length();
		String firstTwoLetters =null;
		if (wordLength >=2){
			firstTwoLetters=chemicalWordLowerCase.substring(0,2);
		}
        if (chemAutomaton.isAccept(as.state)){
        	if (wordLength <= wordLengthRemainingOnLastSuccessfulAnnotations.getValue()){//this annotation is worthy of consideration
	        	if (wordLength < wordLengthRemainingOnLastSuccessfulAnnotations.getValue()){//this annotation is longer than any previously found annotation
	        		successfulAnnotations.clear();
	        		wordLengthRemainingOnLastSuccessfulAnnotations.setValue(wordLength);
	        	}
        		as.untokenisedName=chemicalWord;
	        	successfulAnnotations.add(as);
        	}
        }

        for (char annotationCharacter : stateSymbols) {
            int potentialNextState = chemAutomaton.step(as.state, annotationCharacter);
            if (potentialNextState != -1) {//-1 means this state is not accessible from the previous state
                HashMap<String, List<String>> possibleTokenisationsMap = tokenManager.symbolTokenNamesDict.get(annotationCharacter);
                if (possibleTokenisationsMap != null) {
                    List<String> possibleTokenisations = null;
                    if (firstTwoLetters != null) {
                        possibleTokenisations = possibleTokenisationsMap.get(firstTwoLetters);
                    }
                    if (possibleTokenisations != null) {//next could be a token
                        for (String possibleTokenisation : possibleTokenisations) {
                            if (chemicalWordLowerCase.startsWith(possibleTokenisation)) {
                                AnnotatorState newAs = new AnnotatorState();
                                String newchemicalWord = chemicalWord.substring(possibleTokenisation.length());
                                String newchemicalWordLowerCase = chemicalWordLowerCase.substring(possibleTokenisation.length());
                                newAs.tokens = new ArrayList<String>(as.tokens);
                                newAs.tokens.add(possibleTokenisation);
                                newAs.annot = new ArrayList<Character>(as.annot);
                                newAs.annot.add(annotationCharacter);
                                newAs.state = potentialNextState;
                                //System.out.println("tokened " +newchemicalWord);
                                moveToNextAnnotation(newchemicalWord, newchemicalWordLowerCase, newAs, successfulAnnotations, wordLengthRemainingOnLastSuccessfulAnnotations);
                            }
                        }
                    }
                }
                List<Pattern> possibleRegexes = tokenManager.symbolRegexesDict.get(annotationCharacter);
                if (possibleRegexes != null) {//next could be a regex
                    for (Pattern pattern : possibleRegexes) {
                        Matcher mat = pattern.matcher(chemicalWord);
                        if (mat.lookingAt()) {//match at start
                            AnnotatorState newAs = new AnnotatorState();
                            String newchemicalWord = chemicalWord.substring(mat.group(0).length());
                            String newchemicalWordLowerCase = chemicalWordLowerCase.substring(mat.group(0).length());
                            newAs.tokens = new ArrayList<String>(as.tokens);
                            newAs.tokens.add(mat.group(0));
                            newAs.annot = new ArrayList<Character>(as.annot);
                            newAs.annot.add(annotationCharacter);
                            newAs.state = potentialNextState;
                            //System.out.println("neword regex " +newchemicalWord);
                            moveToNextAnnotation(newchemicalWord, newchemicalWordLowerCase, newAs, successfulAnnotations, wordLengthRemainingOnLastSuccessfulAnnotations);
                        }
                    }
                }
            }
        }
	}

	/**Groups the token annotations for a given word into substituent/s and/or a maingroup and/or functionalTerm by
	 * looking for the endOfSubstituent/endOfMainGroup/endOfFunctionalTerm annotations
	 *
	 * @param annots The annotations for a word.
	 * @return A List of lists of annotations, each list corresponds to a substituent/maingroup/functionalTerm
	 */
	List<List<Character>> chunkAnnotations(List<Character> annots) {
		LinkedList<List<Character>> chunkList = new LinkedList<List<Character>>();
		List<Character> currentTerm = new ArrayList<Character>();
		for (Character annot : annots) {
			currentTerm.add(annot);
			if (annot.equals(endOfSubstituent) || annot.equals(endOfMainGroup) || annot.equals(endOfFunctionalTerm)){
				chunkList.add(currentTerm);
				currentTerm = new ArrayList<Character>();
			}
		}
		return chunkList;
	}
}
