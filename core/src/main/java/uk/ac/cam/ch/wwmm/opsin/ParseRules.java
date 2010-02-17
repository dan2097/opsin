package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Stack;
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
		/** Holds the part of chemical name that has yet to be tokenised */
		String untokenisedChemicalName;
		/** As above but all lower case */
		String untokenisedChemicalNameLowerCase;
	}

	/** Maps regular expression names to regular expressions. */
	private final HashMap<String, String> regexDict;
	/** A DFA encompassing the grammar of a chemical word. */
	private RunAutomaton chemAutomaton;
	/** The allowed symbols in chemAutomaton */
	private char[] stateSymbols;

	private final static char endOfSubstituent = '\u00e9';
	private final static char endOfMainGroup = '\u00e2';
	private final static char endOfFunctionalTerm = '\u00FB';

	private final TokenManager tokenManager;

	/** Initialises the finite-state parser, reading in the rules from regexes.xml.
	 * @param tokenManager
	 * @param resourceGetter
	 *
	 * @throws Exception If the rules file can't be read properly.
	 */
	ParseRules(TokenManager tokenManager, ResourceGetter resourceGetter) throws Exception {
		this.tokenManager = tokenManager;
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
		return AutomatonInitialiser.getAutomaton("chemical", re, true);
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
	ThreeReturnValues<List<ParseTokens>, String, String> getParses(String chemicalWord) throws ParsingException {
		AnnotatorState startingAS = new AnnotatorState();
		startingAS.state = chemAutomaton.getInitialState();
		startingAS.annot = new ArrayList<Character>();
		startingAS.tokens = new ArrayList<String>();
		startingAS.untokenisedChemicalName = chemicalWord;
		startingAS.untokenisedChemicalNameLowerCase = chemicalWord.toLowerCase();
		Stack<AnnotatorState> asStack = new Stack<AnnotatorState>();
		asStack.add(startingAS);

		int wordLengthRemainingOnLastSuccessfulAnnotations = chemicalWord.length();
		int wordLengthRemainingOnLongestAnnotation = chemicalWord.length();
		List<AnnotatorState> successfulAnnotations = new ArrayList<AnnotatorState>();
		AnnotatorState longestAnnotation = new AnnotatorState();//this is the longest annotation. It does not necessarily end in an accept state
		longestAnnotation.state = chemAutomaton.getInitialState();
		longestAnnotation.annot = new ArrayList<Character>();
		longestAnnotation.tokens = new ArrayList<String>();
		longestAnnotation.untokenisedChemicalName = chemicalWord;
		longestAnnotation.untokenisedChemicalNameLowerCase = chemicalWord.toLowerCase();

		while (!asStack.isEmpty()) {
			AnnotatorState as = asStack.pop();
			String untokenisedChemicalNameLowerCase = as.untokenisedChemicalNameLowerCase;
			String untokenisedChemicalName = as.untokenisedChemicalName;
			int wordLength = untokenisedChemicalNameLowerCase.length();
			String firstTwoLetters =null;
			if (wordLength >=2){
				firstTwoLetters = untokenisedChemicalNameLowerCase.substring(0,2);
			}
	        if (chemAutomaton.isAccept(as.state)){
	        	if (wordLength <= wordLengthRemainingOnLastSuccessfulAnnotations){//this annotation is worthy of consideration
		        	if (wordLength < wordLengthRemainingOnLastSuccessfulAnnotations){//this annotation is longer than any previously found annotation
		        		successfulAnnotations.clear();
		        		wordLengthRemainingOnLastSuccessfulAnnotations = wordLength;
		        	}
		        	else if (successfulAnnotations.size()>128){
		        		throw new ParsingException("Ambiguity in OPSIN's chemical grammar has produced more than 128 annotations. Parsing has been aborted. Please report this as a bug");
		        	}
		        	successfulAnnotations.add(as);
	        	}
	        }
        	//record the longest annotation found so it can be reported to the user for debugging
          	if (wordLength < wordLengthRemainingOnLongestAnnotation){
          		wordLengthRemainingOnLongestAnnotation = wordLength;
        		longestAnnotation =as;
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
	                            if (untokenisedChemicalNameLowerCase.startsWith(possibleTokenisation)) {
	                                AnnotatorState newAs = new AnnotatorState();
	                                newAs.untokenisedChemicalName = untokenisedChemicalName.substring(possibleTokenisation.length());
	                                newAs.untokenisedChemicalNameLowerCase = untokenisedChemicalNameLowerCase.substring(possibleTokenisation.length());
	                                newAs.tokens = new ArrayList<String>(as.tokens);
	                                newAs.tokens.add(possibleTokenisation);
	                                newAs.annot = new ArrayList<Character>(as.annot);
	                                newAs.annot.add(annotationCharacter);
	                                newAs.state = potentialNextState;
	                                //System.out.println("tokened " + newAs.untokenisedChemicalName);
	                                asStack.add(newAs);
	                            }
	                        }
	                    }
	                }
	                List<RunAutomaton> possibleAutomata = tokenManager.symbolRegexAutomataDict.get(annotationCharacter);
	                if (possibleAutomata != null) {//next could be a regex
	                    for (RunAutomaton automaton : possibleAutomata) {
	                        int matchLength = automaton.run(untokenisedChemicalName, 0);
	                    	if (matchLength != -1){//matchLength = -1 means it did not match at the start of the string.
	                            AnnotatorState newAs = new AnnotatorState();
	                            newAs.untokenisedChemicalName =  untokenisedChemicalName.substring(matchLength);
	                            newAs.untokenisedChemicalNameLowerCase = untokenisedChemicalNameLowerCase.substring(matchLength);
	                            newAs.tokens = new ArrayList<String>(as.tokens);
	                            newAs.tokens.add(untokenisedChemicalName.substring(0, matchLength));
	                            newAs.annot = new ArrayList<Character>(as.annot);
	                            newAs.annot.add(annotationCharacter);
	                            newAs.state = potentialNextState;
	                            //System.out.println("neword automata " + newAs.untokenisedChemicalName);
	                            asStack.add(newAs);
	                        }
	                    }
	                }
	                List<Pattern> possibleRegexes = tokenManager.symbolRegexesDict.get(annotationCharacter);
	                if (possibleRegexes != null) {//next could be a regex
	                    for (Pattern pattern : possibleRegexes) {
	                        Matcher mat = pattern.matcher(untokenisedChemicalName);
	                        if (mat.lookingAt()) {//match at start
	                            AnnotatorState newAs = new AnnotatorState();
	                            newAs.untokenisedChemicalName =  untokenisedChemicalName.substring(mat.group(0).length());
	                            newAs.untokenisedChemicalNameLowerCase = untokenisedChemicalNameLowerCase.substring(mat.group(0).length());
	                            newAs.tokens = new ArrayList<String>(as.tokens);
	                            newAs.tokens.add(mat.group(0));
	                            newAs.annot = new ArrayList<Character>(as.annot);
	                            newAs.annot.add(annotationCharacter);
	                            newAs.state = potentialNextState;
	                            //System.out.println("neword regex " + newAs.untokenisedChemicalName);
	                            asStack.add(newAs);
	                        }
	                    }
	                }
	            }
	        }
		}
		List<ParseTokens> outputList = new ArrayList<ParseTokens>();
		String uninterpretableName = chemicalWord;
		String unparseableName = longestAnnotation.untokenisedChemicalName;
		if (successfulAnnotations.size()>0){//at least some of the name could be interpreted into a substituent/full/functionalTerm
			for(AnnotatorState as : successfulAnnotations) {
				ParseTokens pt =new ParseTokens(as.tokens, as.annot);
				outputList.add(pt);
				uninterpretableName=as.untokenisedChemicalName;//all acceptable annotator states found should have the same untokenisedName
			}
		}
		return new ThreeReturnValues<List<ParseTokens>, String, String>(outputList, uninterpretableName, unparseableName);
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
