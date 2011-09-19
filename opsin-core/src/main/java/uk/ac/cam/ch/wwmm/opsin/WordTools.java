package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import static uk.ac.cam.ch.wwmm.opsin.OpsinTools.*;

/**
 * Tools for dealing uniformly with unusually-formed words.
 */
class WordTools {
	/**
	 * Splits cases where the parseTokensList describes a functionalTerm in addition to another mainGroup/substituent into two parseWords
	 * This occurs if the name is formally missing a space e.g. ethylthiocyanate.
	 * If multiple parses are present then it may be possible to disambiguate between them:
	 * 	parses with omitted spaces are discarded if a parse without omitted space is found
	 * 	parses with shorter functional terms are discarded e.g. ethylthiocyanate is [ethyl] [thiocyanate] not [ethylthio] [cyanate]
	 * @param parseTokensList
	 * @param chemicalName
	 * @return
	 * @throws ParsingException
	 */
	static List<ParseWord> splitIntoParseWords(List<ParseTokens> parseTokensList, String chemicalName) throws ParsingException {
		List<ParseTokens> wellFormedParseTokens = new ArrayList<ParseTokens>();//these are all in the same word as would be expected
		List<List<ParseTokens>> splitParseTokensForEachParseTokens = new ArrayList<List<ParseTokens>>();
		/*
		 * Each ParseTokens is split into the number of words it describes
		 * e.g. ethylchloride has one interpretation so splitParseTokensList will have one entry
		 * This entry will be formed of TWO parseTokens, one for the ethyl and one for the chloride
		 */
		int leastWordsInOmmittedSpaceParse = Integer.MAX_VALUE;//we want the least number of words i.e. less omitted spaces
		int longestFunctionalTermEncountered = 0;//we want the longest functional term
		for (ParseTokens parseTokens : parseTokensList) {
			List<Character> annotations = parseTokens.getAnnotations();
			List<List<Character>> chunkedAnnotations = chunkAnnotations(annotations);//chunked into mainGroup/substituent/functionalTerm
			if (containsOmittedSpace(chunkedAnnotations)){
				List<ParseTokens> omittedWordParseTokens = new ArrayList<ParseTokens>();
				List<String> tokens = parseTokens.getTokens();
				List<Character> newAnnotations = new ArrayList<Character>();
				List<String> newTokens = new ArrayList<String>();
				int currentFunctionalTermLength = 0;
				int annotPos = 0;
				for (List<Character> annotationList : chunkedAnnotations) {
					Character finalAnnotationInList = annotationList.get(annotationList.size() - 1);
					if (finalAnnotationInList.equals(END_OF_FUNCTIONALTERM) && newAnnotations.size() > 0) {
						//create a new parseTokens for the substituent/maingroup preceding the functional term
						//not necessary if the functional term is the first thing to be read e.g. in the case of poly
						omittedWordParseTokens.add(new ParseTokens(newTokens, newAnnotations));
						newAnnotations = new ArrayList<Character>();
						newTokens = new ArrayList<String>();
					}
					for (Character annotation : annotationList) {
						newAnnotations.add(annotation);
						newTokens.add(tokens.get(annotPos++));
					}
					if (finalAnnotationInList.equals(END_OF_FUNCTIONALTERM) || finalAnnotationInList.equals(END_OF_MAINGROUP) || annotPos == tokens.size()) {
						omittedWordParseTokens.add(new ParseTokens(newTokens, newAnnotations));
						if (finalAnnotationInList.equals(END_OF_FUNCTIONALTERM)){
							currentFunctionalTermLength = StringTools.stringListToString(newTokens, "").length();
						}
						newAnnotations = new ArrayList<Character>();
						newTokens = new ArrayList<String>();
					}
				}
				if (omittedWordParseTokens.size() <= leastWordsInOmmittedSpaceParse){
					if (omittedWordParseTokens.size() < leastWordsInOmmittedSpaceParse){
						splitParseTokensForEachParseTokens.clear();
						leastWordsInOmmittedSpaceParse = omittedWordParseTokens.size();
						longestFunctionalTermEncountered = 0;
					}
					if (currentFunctionalTermLength >=longestFunctionalTermEncountered){
						if (currentFunctionalTermLength > longestFunctionalTermEncountered){
							splitParseTokensForEachParseTokens.clear();
							longestFunctionalTermEncountered =currentFunctionalTermLength;
						}
						splitParseTokensForEachParseTokens.add(omittedWordParseTokens);
					}
				}
			} else {
				wellFormedParseTokens.add(parseTokens);
			}
		}
		List<ParseWord> parseWords = new ArrayList<ParseWord>();
		if (!wellFormedParseTokens.isEmpty()) {
			parseWords.add(new ParseWord(chemicalName, wellFormedParseTokens));
		} else {
			for (int i = 0; i < leastWordsInOmmittedSpaceParse; i++) {
				List<ParseTokens> parseTokensForWord = new ArrayList<ParseTokens>();
				for (List<ParseTokens> parseTokens : splitParseTokensForEachParseTokens) {
					if (!parseTokensForWord.contains(parseTokens.get(i))){//if only one word is ambiguous there is no need for the unambiguous word to have multiple identical interpretation
						parseTokensForWord.add(parseTokens.get(i));
					}
				}
				parseWords.add(new ParseWord(StringTools.stringListToString(parseTokensForWord.get(0).getTokens(), ""), parseTokensForWord));
			}
		}
		return parseWords;
	}

	private static boolean containsOmittedSpace(List<List<Character>> chunkedAnnotations) {
		if (chunkedAnnotations.size() > 1){//there are multiple subsitutents/maingroup/functionalterms
			for (List<Character> annotationList : chunkedAnnotations) {
				for (Character annotation : annotationList) {
					if (annotation.equals(END_OF_FUNCTIONALTERM)){
						return true;
					}
				}
			}
		}
		return false;
	}

	/**Groups the token annotations for a given word into substituent/s and/or a maingroup and/or functionalTerm by
	 * looking for the endOfSubstituent/endOfMainGroup/endOfFunctionalTerm annotations
	 *
	 * @param annots The annotations for a word.
	 * @return A List of lists of annotations, each list corresponds to a substituent/maingroup/functionalTerm
	 */
	static List<List<Character>> chunkAnnotations(List<Character> annots) {
		LinkedList<List<Character>> chunkList = new LinkedList<List<Character>>();
		List<Character> currentTerm = new ArrayList<Character>();
		for (Character annot : annots) {
			currentTerm.add(annot);
			if (annot.equals(END_OF_SUBSTITUENT) || annot.equals(END_OF_MAINGROUP) || annot.equals(END_OF_FUNCTIONALTERM)) {
				chunkList.add(currentTerm);
				currentTerm = new ArrayList<Character>();
			}
		}
		return chunkList;
	}
	
	/**
	 * Works left to right removing spaces if there are too many opening brackets
	 * @param name
	 * @return
	 * @throws ParsingException If brackets are unbalanced and cannot be balanced by removing whitespace
	 */
	static String removeWhiteSpaceIfBracketsAreUnbalanced(String name) throws ParsingException {
		int bracketLevel = 0;
		int stringLength = name.length();
		for (int i = 0; i < stringLength; i++) {
			char c = name.charAt(i);
			if (c == '(' || c == '[' || c == '{') {
				bracketLevel++;
			} else if (c == ')' || c == ']' || c == '}') {
				bracketLevel--;
			} else if (c == ' ' && bracketLevel > 0) {//brackets unbalanced and a space has been encountered!
				name = name.substring(0, i) + name.substring(i + 1);
				stringLength = name.length();
				i--;
			}
		}
		if (bracketLevel > 0) {
			throw new ParsingException("Unmatched opening bracket found in :" + name);
		} else if (bracketLevel < 0) {
			throw new ParsingException("Unmatched closing bracket found in :" + name);
		}
		return name;
	}
}
