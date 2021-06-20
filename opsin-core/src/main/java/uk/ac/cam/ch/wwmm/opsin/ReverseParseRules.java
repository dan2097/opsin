package uk.ac.cam.ch.wwmm.opsin;

import java.io.IOException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import dk.brics.automaton.RunAutomaton;

/**
 * The same as ParseRules but works from right to left
 * 
 * Performs finite-state allocation of roles ("annotations") to tokens:
 * The chemical name is broken down into tokens e.g. ethyl -->eth yl by applying the chemical grammar in regexes.xml
 * The tokens eth and yl are associated with a letter which is referred to here as an annotation which is the role of the token.
 * These letters are defined in regexes.xml and would in this case have the meaning alkaneStem and inlineSuffix
 *
 * The chemical grammar employs the annotations associated with the tokens when deciding what may follow what has already been seen
 * e.g. you cannot start a chemical name with yl and an optional e is valid after an arylGroup
 *
 * @author dl387
 *
 */
class ReverseParseRules {

	/** A DFA encompassing the grammar of a chemical word. */
	private final RunAutomaton chemAutomaton;
	/** The allowed symbols in chemAutomaton */
	private final char[] stateSymbols;
	
	private final OpsinRadixTrie[] symbolTokenNamesDictReversed;
	private final RunAutomaton[] symbolRegexAutomataDictReversed;
	private final Pattern[] symbolRegexesDictReversed;

	/** 
	 * Creates a right to left parser that can parse a substituent/full/functional word
	 * @param resourceManager
	 * @throws IOException 
	 */
	ReverseParseRules(ResourceManager resourceManager) throws IOException{
		resourceManager.populatedReverseTokenMappings();
		this.chemAutomaton = resourceManager.getReverseChemicalAutomaton();
		this.symbolTokenNamesDictReversed = resourceManager.getSymbolTokenNamesDictReversed();
		this.symbolRegexAutomataDictReversed = resourceManager.getSymbolRegexAutomataDictReversed();
		this.symbolRegexesDictReversed = resourceManager.getSymbolRegexesDictReversed();
		this.stateSymbols = chemAutomaton.getCharIntervals();
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
		AnnotatorState initialState = new AnnotatorState(chemAutomaton.getInitialState(), '\0', chemicalWord.length(), true, null);
		String chemicalWordLowerCase = StringTools.lowerCaseAsciiString(chemicalWord);
		ArrayDeque<AnnotatorState> asStack = new ArrayDeque<>();
		asStack.add(initialState);

		int posInNameOfLastSuccessfulAnnotations = chemicalWord.length();
		List<AnnotatorState> successfulAnnotations = new ArrayList<>();
		AnnotatorState longestAnnotation = initialState;//this is the longest annotation. It does not necessarily end in an accept state
		int stateSymbolsSize = stateSymbols.length;
		while (!asStack.isEmpty()) {
			AnnotatorState as = asStack.removeFirst();
			int posInName = as.getPosInName();
			if (chemAutomaton.isAccept(as.getState())){
				if (posInName <= posInNameOfLastSuccessfulAnnotations){//this annotation is worthy of consideration
					if (posInName < posInNameOfLastSuccessfulAnnotations){//this annotation is longer than any previously found annotation
						successfulAnnotations.clear();
						posInNameOfLastSuccessfulAnnotations = posInName;
					}
					else if (successfulAnnotations.size() > 128){
						throw new ParsingException("Ambiguity in OPSIN's chemical grammar has produced more than 128 annotations. Parsing has been aborted. Please report this as a bug");
					}
					successfulAnnotations.add(as);
				}
			}
			//record the longest annotation found so it can be reported to the user for debugging
			if (posInName < longestAnnotation.getPosInName()){
				longestAnnotation = as;
			}

			for (int i = 0; i < stateSymbolsSize; i++) {
				char annotationCharacter = stateSymbols[i];
				int potentialNextState = chemAutomaton.step(as.getState(), annotationCharacter);
				if (potentialNextState != -1) {//-1 means this state is not accessible from the previous state
					OpsinRadixTrie possibleTokenisationsTrie = symbolTokenNamesDictReversed[i];
					if (possibleTokenisationsTrie != null) {
						List<Integer> possibleTokenisations = possibleTokenisationsTrie.findMatchesReadingStringRightToLeft(chemicalWordLowerCase, posInName);
						if (possibleTokenisations != null) {//next could be a token
							for (int j = 0, l = possibleTokenisations.size(); j < l; j++) {//typically list size will be 1 so this is faster than an iterator
								int tokenizationIndex = possibleTokenisations.get(j);
								AnnotatorState newAs = new AnnotatorState(potentialNextState, annotationCharacter, tokenizationIndex, false, as);
								//System.out.println("tokened " + chemicalWordLowerCase.substring(tokenizationIndex, posInName));
								asStack.add(newAs);
							}
						}
					}
					RunAutomaton possibleAutomata = symbolRegexAutomataDictReversed[i];
					if (possibleAutomata != null) {//next could be an automaton
						int matchLength = runInReverse(possibleAutomata, chemicalWord, posInName);
						if (matchLength != -1){//matchLength = -1 means it did not match
							int tokenizationIndex = posInName - matchLength;
							AnnotatorState newAs = new AnnotatorState(potentialNextState, annotationCharacter, tokenizationIndex, true, as);
							//System.out.println("neword automata " + chemicalWord.substring(tokenizationIndex, posInName));
							asStack.add(newAs);
						}
					}
					Pattern possibleRegex = symbolRegexesDictReversed[i];
					if (possibleRegex != null) {//next could be a regex
						Matcher mat = possibleRegex.matcher(chemicalWord).region(0, posInName);
						mat.useTransparentBounds(true);
						if (mat.find()) {//match at end (patterns use $ anchor)
							int tokenizationIndex = posInName - mat.group(0).length();
							AnnotatorState newAs = new AnnotatorState(potentialNextState, annotationCharacter, tokenizationIndex, true, as);
							//System.out.println("neword regex " + mat.group(0));
							asStack.add(newAs);
						}
					}
				}
			}
		}
		List<ParseTokens> outputList = new ArrayList<>();
		String uninterpretableName = chemicalWord;
		String unparseableName = chemicalWord.substring(0, longestAnnotation.getPosInName());
		if (successfulAnnotations.size() > 0){//at least some of the name could be interpreted into a substituent/full/functionalTerm
			int bestAcceptPosInName = -1;
			for(AnnotatorState as : successfulAnnotations) {
				outputList.add(convertAnnotationStateToParseTokens(as, chemicalWord, chemicalWordLowerCase));
				bestAcceptPosInName = as.getPosInName();//all acceptable annotator states found should have the same posInName
			}
			uninterpretableName = chemicalWord.substring(0, bestAcceptPosInName);
		}
		return new ParseRulesResults(outputList, uninterpretableName, unparseableName);
	}

	/**
	 * Returns the length of the longest accepted run of the given string
	 * starting at pos in the string and working backwards
	 * @param automaton 
	 * @param s the string
	 * @param indexAfterFirstchar pos in string to start at
	 * @return length of the longest accepted run, -1 if no run is accepted
	 */
	private int runInReverse(RunAutomaton automaton, String s, int indexAfterFirstchar) {
		int state = automaton.getInitialState();
		int max = -1;
		for (int pos = indexAfterFirstchar -1; ; pos--) {
			if (automaton.isAccept(state)){
				max = indexAfterFirstchar - 1 - pos;
			}
			if (pos == -1){
				break;
			}
			state = automaton.step(state, s.charAt(pos));
			if (state == -1){
				break;
			}
		}
		return max;
	}
	
	private ParseTokens convertAnnotationStateToParseTokens(AnnotatorState as, String chemicalWord, String chemicalWordLowerCase) {
		List<String> tokens = new ArrayList<>();
		List<Character> annotations = new ArrayList<>();
		AnnotatorState previousAs;
		while ((previousAs = as.getPreviousAs()) != null) {
			if (as.isCaseSensitive()) {
				tokens.add(chemicalWord.substring(as.getPosInName(), previousAs.getPosInName()));
			}
			else{
				tokens.add(chemicalWordLowerCase.substring(as.getPosInName(), previousAs.getPosInName()));
			}
			annotations.add(as.getAnnot());
			as = previousAs;
		}
		return new ParseTokens(tokens, annotations);
	}
}
