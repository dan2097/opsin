package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
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

	/** A DFA encompassing the grammar of a chemical word. */
	private final RunAutomaton chemAutomaton;
	/** The allowed symbols in chemAutomaton */
	private final char[] stateSymbols;
	
	private final OpsinRadixTrie[] symbolTokenNamesDict;
	private final RunAutomaton[] symbolRegexAutomataDict;
	private final Pattern[] symbolRegexesDict;
	
	private final AnnotatorState initialState;

	/**
	 * Creates a left to right parser that can parse a substituent/full/functional word
	 * @param resourceManager
	 */
	ParseRules(ResourceManager resourceManager){
		this.chemAutomaton = resourceManager.getChemicalAutomaton();
		this.symbolTokenNamesDict = resourceManager.getSymbolTokenNamesDict();
		this.symbolRegexAutomataDict = resourceManager.getSymbolRegexAutomataDict();
		this.symbolRegexesDict = resourceManager.getSymbolRegexesDict();
		this.stateSymbols = chemAutomaton.getCharIntervals();
		this.initialState = new AnnotatorState(chemAutomaton.getInitialState(), '\0', 0, true, null);
	}

	/**Determines the possible annotations for a chemical word
	 * Returns a list of parses and how much of the word could not be interpreted
	 * e.g. usually the list will have only one parse and the string will equal ""
	 * For something like ethyloxime. The list will contain the parse for ethyl and the string will equal "oxime" as it was unparsable
	 * For something like eth no parses would be found and the string will equal "eth"
	 *
	 * @param chemicalWord
	 * @return Results of parsing
	 * @throws ParsingException
	 */
	public ParseRulesResults getParses(String chemicalWord) throws ParsingException {
		String chemicalWordLowerCase = StringTools.lowerCaseAsciiString(chemicalWord);
		ArrayDeque<AnnotatorState> asStack = new ArrayDeque<AnnotatorState>();
		asStack.add(initialState);

		int posInNameOfLastSuccessfulAnnotations = 0;
		List<AnnotatorState> successfulAnnotations = new ArrayList<AnnotatorState>();
		AnnotatorState longestAnnotation = initialState;//this is the longest annotation. It does not necessarily end in an accept state
		int stateSymbolsSize = stateSymbols.length;
		while (!asStack.isEmpty()) {
			AnnotatorState as = asStack.removeFirst();
			int posInName = as.getPosInName();
			if (chemAutomaton.isAccept(as.getState())){
				if (posInName >= posInNameOfLastSuccessfulAnnotations){//this annotation is worthy of consideration
					if (posInName > posInNameOfLastSuccessfulAnnotations){//this annotation is longer than any previously found annotation
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
			if (posInName > longestAnnotation.getPosInName()){
				longestAnnotation = as;
			}

			for (int i = 0; i < stateSymbolsSize; i++) {
				char annotationCharacter = stateSymbols[i];
				int potentialNextState = chemAutomaton.step(as.getState(), annotationCharacter);
				if (potentialNextState != -1) {//-1 means this state is not accessible from the previous state
					OpsinRadixTrie possibleTokenisationsTrie = symbolTokenNamesDict[i];
					if (possibleTokenisationsTrie != null) {
						List<Integer> possibleTokenisations = possibleTokenisationsTrie.findMatches(chemicalWordLowerCase, posInName);
						if (possibleTokenisations != null) {//next could be a token
							for (int j = 0, l = possibleTokenisations.size(); j < l; j++) {//typically list size will be 1 so this is faster than an iterator
								int tokenizationIndex = possibleTokenisations.get(j);
								AnnotatorState newAs = new AnnotatorState(potentialNextState, annotationCharacter, tokenizationIndex, false, as);
								//System.out.println("tokened " + chemicalWordLowerCase.substring(posInName, tokenizationIndex));
								asStack.add(newAs);
							}
						}
					}
					RunAutomaton possibleAutomata = symbolRegexAutomataDict[i];
					if (possibleAutomata != null) {//next could be an automaton
						int matchLength = possibleAutomata.run(chemicalWord, posInName);
						if (matchLength != -1){//matchLength = -1 means it did not match
							int tokenizationIndex = posInName + matchLength;
							AnnotatorState newAs = new AnnotatorState(potentialNextState, annotationCharacter, tokenizationIndex, true, as);
							//System.out.println("neword automata " + chemicalWord.substring(posInName, tokenizationIndex));
							asStack.add(newAs);
						}
					}
					Pattern possibleRegex = symbolRegexesDict[i];
					if (possibleRegex != null) {//next could be a regex
						Matcher mat = possibleRegex.matcher(chemicalWord).region(posInName, chemicalWord.length());
						mat.useTransparentBounds(true);
						if (mat.lookingAt()) {//match at start
							int tokenizationIndex = posInName + mat.group(0).length();
							AnnotatorState newAs = new AnnotatorState(potentialNextState, annotationCharacter, tokenizationIndex, true, as);
							//System.out.println("neword regex " + mat.group(0));
							asStack.add(newAs);
						}
					}
				}
			}
		}
		List<ParseTokens> outputList = new ArrayList<ParseTokens>();
		String uninterpretableName = chemicalWord;
		String unparseableName = chemicalWord.substring(longestAnnotation.getPosInName());
		if (successfulAnnotations.size() > 0){//at least some of the name could be interpreted into a substituent/full/functionalTerm
			int bestAcceptPosInName = -1;
			for(AnnotatorState as : successfulAnnotations) {
				outputList.add(convertAnnotationStateToParseTokens(as, chemicalWord, chemicalWordLowerCase));
				bestAcceptPosInName = as.getPosInName();//all acceptable annotator states found should have the same posInName
			}
			uninterpretableName = chemicalWord.substring(bestAcceptPosInName);
		}
		return new ParseRulesResults(outputList, uninterpretableName, unparseableName);
	}

	private ParseTokens convertAnnotationStateToParseTokens(AnnotatorState as, String chemicalWord, String chemicalWordLowerCase) {
		List<String> tokens = new ArrayList<String>();
		List<Character> annotations = new ArrayList<Character>();
		AnnotatorState previousAs;
		while ((previousAs = as.getPreviousAs()) != null) {
			if (as.isCaseSensitive()) {
				tokens.add(chemicalWord.substring(previousAs.getPosInName(), as.getPosInName()));
			}
			else{
				tokens.add(chemicalWordLowerCase.substring(previousAs.getPosInName(), as.getPosInName()));
			}
			annotations.add(as.getAnnot());
			as = previousAs;
		}
		Collections.reverse(tokens);
		Collections.reverse(annotations);
		return new ParseTokens(tokens, annotations);
	}
}
