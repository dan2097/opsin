package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import dk.brics.automaton.RunAutomaton;

/**
 * Instantiate via NameToStructure.getOpsinParser()
 * 
 * Performs finite-state allocation of roles ("annotations") to tokens:
 * The chemical name is broken down into tokens e.g. ethyl -->eth yl by applying the chemical grammar in regexes.xml
 * The tokens eth and yl are associated with a letter which is referred to here as an annotation which is the role of the token.
 * These letters are defined in regexes.xml and would in this case have the meaning alkaneStem and inlineSuffix
 *
 * The chemical grammar employs the annotations associated with the tokens when deciding what may follow what has already been seen
 * e.g. you cannot start a chemical name with yl and an optional e is valid after an arylGroup
 *
 * @author ptc24
 * @author dl387
 *
 */
public class ParseRules {

	/** A "struct" containing bits of state needed during finite-state parsing. */
	private static class AnnotatorState {
		/** The current state of the DFA. */
		int state;
		/** The annotation so far. */
		List<Character> annot;
		/** The strings these annotations correspond to. */
		ArrayList<String> tokens;
		/** The index of the first char in the chemical name that has yet to be tokenised */
		int posInName;
	}

	/** A DFA encompassing the grammar of a chemical word. */
	private RunAutomaton chemAutomaton;
	/** The allowed symbols in chemAutomaton */
	private char[] stateSymbols;

	private final ResourceManager resourceManager;

	/**
	 * Creates a left to right parser that can parse a substituent/full/functional word
	 * @param resourceManager
	 */
	ParseRules(ResourceManager resourceManager){
		this.resourceManager = resourceManager;
		chemAutomaton = resourceManager.chemicalAutomaton;
		stateSymbols = chemAutomaton.getCharIntervals();
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
	public ParseRulesResults getParses(String chemicalWord) throws ParsingException {
		String chemicalWordLowerCase = chemicalWord.toLowerCase();
		AnnotatorState startingAS = new AnnotatorState();
		startingAS.state = chemAutomaton.getInitialState();
		startingAS.annot = new ArrayList<Character>();
		startingAS.tokens = new ArrayList<String>();
		startingAS.posInName = 0;
		LinkedList<AnnotatorState> asStack = new LinkedList<AnnotatorState>();
		asStack.add(startingAS);

		int posInNameOfLastSuccessfulAnnotations = 0;
		List<AnnotatorState> successfulAnnotations = new ArrayList<AnnotatorState>();
		AnnotatorState longestAnnotation = new AnnotatorState();//this is the longest annotation. It does not necessarily end in an accept state
		longestAnnotation.state = chemAutomaton.getInitialState();
		longestAnnotation.annot = new ArrayList<Character>();
		longestAnnotation.tokens = new ArrayList<String>();
		longestAnnotation.posInName = 0;
		int stateSymbolsSize = stateSymbols.length;
		while (!asStack.isEmpty()) {
			AnnotatorState as = asStack.removeFirst();
			int posInName = as.posInName;
			if (chemAutomaton.isAccept(as.state)){
				if (posInName >= posInNameOfLastSuccessfulAnnotations){//this annotation is worthy of consideration
					if (posInName > posInNameOfLastSuccessfulAnnotations){//this annotation is longer than any previously found annotation
						successfulAnnotations.clear();
						posInNameOfLastSuccessfulAnnotations = posInName;
					}
					else if (successfulAnnotations.size()>128){
						throw new ParsingException("Ambiguity in OPSIN's chemical grammar has produced more than 128 annotations. Parsing has been aborted. Please report this as a bug");
					}
					successfulAnnotations.add(as);
				}
			}
			//record the longest annotation found so it can be reported to the user for debugging
			if (posInName > longestAnnotation.posInName){
				longestAnnotation = as;
			}

			for (int i = 0; i < stateSymbolsSize; i++) {
				char annotationCharacter = stateSymbols[i];
				int potentialNextState = chemAutomaton.step(as.state, annotationCharacter);
				if (potentialNextState != -1) {//-1 means this state is not accessible from the previous state
					OpsinRadixTrie possibleTokenisationsTrie = resourceManager.symbolTokenNamesDict[i];
					if (possibleTokenisationsTrie != null) {
						List<Integer> possibleTokenisations = possibleTokenisationsTrie.findMatches(chemicalWordLowerCase, posInName);
						if (possibleTokenisations != null) {//next could be a token
							for (int j = 0, l = possibleTokenisations.size(); j < l; j++) {//typically list size will be 1 so this is faster than an iterator
								int tokenizationIndex = possibleTokenisations.get(j);
								AnnotatorState newAs = new AnnotatorState();
								newAs.posInName = tokenizationIndex;
								newAs.tokens = new ArrayList<String>(as.tokens);
								newAs.tokens.add(chemicalWordLowerCase.substring(posInName, tokenizationIndex));
								newAs.annot = new ArrayList<Character>(as.annot);
								newAs.annot.add(annotationCharacter);
								newAs.state = potentialNextState;
								//System.out.println("tokened " + chemicalWordLowerCase.substring(posInName, tokenizationIndex));
								asStack.add(newAs);
							}
						}
					}
					List<RunAutomaton> possibleAutomata = resourceManager.symbolRegexAutomataDict[i];
					if (possibleAutomata != null) {//next could be an automaton
						for (int j = 0, l = possibleAutomata.size(); j < l; j++) {
							RunAutomaton automaton = possibleAutomata.get(j);
							int matchLength = automaton.run(chemicalWord, posInName);
							if (matchLength != -1){//matchLength = -1 means it did not match
								AnnotatorState newAs = new AnnotatorState();
								newAs.posInName = posInName + matchLength;
								newAs.tokens = new ArrayList<String>(as.tokens);
								newAs.tokens.add(chemicalWord.substring(posInName, posInName + matchLength));
								newAs.annot = new ArrayList<Character>(as.annot);
								newAs.annot.add(annotationCharacter);
								newAs.state = potentialNextState;
								//System.out.println("neword automata " + chemicalWord.substring(posInName, posInName + matchLength));
								asStack.add(newAs);
							}
						}
					}
					List<Pattern> possibleRegexes = resourceManager.symbolRegexesDict[i];
					if (possibleRegexes != null) {//next could be a regex
						for (int j = 0, l = possibleRegexes.size(); j < l; j++) {
							Pattern pattern = possibleRegexes.get(j);
							Matcher mat = pattern.matcher(chemicalWord).region(posInName, chemicalWord.length());
							if (mat.lookingAt()) {//match at start
								AnnotatorState newAs = new AnnotatorState();
								String matchedString = mat.group(0);
								newAs.posInName = posInName + matchedString.length();
								newAs.tokens = new ArrayList<String>(as.tokens);
								newAs.tokens.add(matchedString);
								newAs.annot = new ArrayList<Character>(as.annot);
								newAs.annot.add(annotationCharacter);
								newAs.state = potentialNextState;
								//System.out.println("neword regex " + matchedString);
								asStack.add(newAs);
							}
						}
					}
				}
			}
		}
		List<ParseTokens> outputList = new ArrayList<ParseTokens>();
		String uninterpretableName = chemicalWord;
		String unparseableName = chemicalWord.substring(longestAnnotation.posInName);
		if (successfulAnnotations.size() > 0){//at least some of the name could be interpreted into a substituent/full/functionalTerm
			int bestAcceptPosInName = -1;
			for(AnnotatorState as : successfulAnnotations) {
				ParseTokens pt = new ParseTokens(as.tokens, as.annot);
				outputList.add(pt);
				bestAcceptPosInName = as.posInName;//all acceptable annotator states found should have the same posInName
			}
			uninterpretableName = chemicalWord.substring(bestAcceptPosInName);
		}
		return new ParseRulesResults(outputList, uninterpretableName, unparseableName);
	}
}
