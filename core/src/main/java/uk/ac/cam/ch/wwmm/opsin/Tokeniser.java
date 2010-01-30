package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.regex.Pattern;

/**
 * Uses OPSIN's DFA based grammar to break a name in to tokens with associated meanings ("annotations").
 * @author dl387
 *
 */
class Tokeniser {


	private ParseRules parseRules;
	private Pattern matchWhitespace =Pattern.compile("\\s+");
	private final static char endOfFunctionalTerm = '\u00FB';
	
	Tokeniser(ParseRules parseRules) {
		this.parseRules = parseRules;
	}

	/**
	 * Master method for tokenizing chemical names into words and within words into tokens
	 * @param name The chemical name.
	 * @return
	 * @throws ParsingException 
	 */
	Parse tokenize(String name) throws ParsingException {
		Parse parse = new Parse(name);
		List<String> words = new LinkedList<String>(Arrays.asList(matchWhitespace.split(name)));//initially assume spaces are hard delimiters
		
		for (int i = 0; i < words.size(); i++) {
			String word = words.get(i);
			
			if (words.size()!=1 && StringTools.isLackingCloseBracket(word)){//erroneous space!!
				//something like [(3,4,5-trihydroxy-6-methyl -phenyl..
				words.set(i, word + words.get(i+1));
				words.remove(i+1);
				i--;
				continue;
			}
			
			/*
			 * Returns
			 * List of parses where at least some of the name was assigned a role
			 * Section of name that was uninterpretable (or "" if none was)
			 * Section of name that was unparsable (or "" if none was). This is always shorter or equal to the above string
			 */
			ThreeReturnValues<List<ParseTokens>, String, String> output= parseRules.getParses(word);
			List<ParseTokens> parseTokens =output.getFirst();
			String unparseableName = output.getSecond();

			if (parseTokens.size()>0 && unparseableName.equals("")){//word was fully interpretable
				//If something like ethylchloride is encountered this should be split back to ethyl chloride and there will be 2 ParseWords returned
				//In cases of properly formed names there will be only one ParseWord
				//If there are two parses one of which assumes a missing space and one of which does not the former is discarded
				List<ParseWord> parseWords = splitIntoParseWords(parseTokens, word);
				words.remove(i);
				for (int j = 0; j < parseWords.size(); j++) {
					words.add(i +j, parseWords.get(j).getWord());
					parse.addWord(parseWords.get(j));
				}
				i=i+parseWords.size()-1;//move the indice onwards if a word has been split
			}
			else{//word is unparsable as is. Try and remove a space and try again
				//TODO add a warning message if this code is invoked. A name invoking this is unambiguously BAD
				if (i +1 < words.size()){//join word with the next word
					words.set(i, word + words.get(i +1));
					words.remove(i+1);
					i--;
				}
				else{
					throw new ParsingException(name + " is unparsable due to the following word being unparsable: " + word+ " (spaces will have been removed as they are assumed to be erroneous if parsing fails with them there)");
				}
			}
		}
		return parse;
	}

	private List<ParseWord> splitIntoParseWords(List<ParseTokens> parseTokensList, String chemicalName) throws ParsingException {
		List<ParseTokens> wellFormedParseTokens = new ArrayList<ParseTokens>();//these are all in the same word as would be expected
		List<List<ParseTokens>> omittedWordParseTokensList = new ArrayList<List<ParseTokens>>();//these are grouped into words e.g. ethylchloride will have a list of parseTokens for the ethyl and chloride
		omittedWordParseTokensList.add(new ArrayList<ParseTokens>());
		omittedWordParseTokensList.add(new ArrayList<ParseTokens>());//only 1 space is allowed to be omitted
		for (ParseTokens parseTokens : parseTokensList) {
			List<Character> annotations = parseTokens.getAnnotations();
			List<List<Character>> chunkedAnnotations = parseRules.chunkAnnotations(annotations);//chunked into mainGroup/substituent/functionalTerm
			if (chunkedAnnotations.size()>1 && annotations.contains(endOfFunctionalTerm)){//must be an omitted space as not allowed to have a functionalTerm and anything else
				List<String> tokens = parseTokens.getTokens();
				List<Character> newAnnotations = new ArrayList<Character>();
				List<String> newTokens = new ArrayList<String>();
				int annotPos =0;
				int wordCounter=0;
				for (List<Character> annotationList : chunkedAnnotations) {
					boolean functionalTermNext =false;
					if (annotationList.get(annotationList.size()-1).equals(endOfFunctionalTerm)){
						functionalTermNext=true;
						if (newAnnotations.size()>0){//create a new parseTokens, unless nothing has been read yet e.g. in the case of poly
							if (wordCounter >=2){
								throw new ParsingException("Name appears to have 2 or more omitted spaces!");
							}
							omittedWordParseTokensList.get(wordCounter++).add(new ParseTokens(newTokens, newAnnotations));
							newAnnotations = new ArrayList<Character>();
							newTokens = new ArrayList<String>();
						}
					}
					for (Character annotation : annotationList) {
						newAnnotations.add(annotation);
						newTokens.add(tokens.get(annotPos));
						annotPos++;
					}
					if (functionalTermNext){
						if (wordCounter >=2){
							throw new ParsingException("Name appears to have 2 or more omitted spaces!");
						}
						if (omittedWordParseTokensList.get(wordCounter).size()==0){//only want one parse for functionTerm as functional terms do not generate multiple interpretations (and more to the point this is the simplest way of achieving that!)
							omittedWordParseTokensList.get(wordCounter++).add(new ParseTokens(newTokens, newAnnotations));
						}
						else{
							wordCounter++;
						}
						newAnnotations = new ArrayList<Character>();
						newTokens = new ArrayList<String>();
					}
				}
				if (newAnnotations.size()>0){
					if (wordCounter >=2){
						throw new ParsingException("Name appears to have 2 or more omitted spaces!");
					}
					omittedWordParseTokensList.get(wordCounter++).add(new ParseTokens(newTokens, newAnnotations));
				}
			}
			else{
				wellFormedParseTokens.add(parseTokens);
			}
		}
		List<ParseWord> parseWords = new ArrayList<ParseWord>();
		if (wellFormedParseTokens.size()>0){
			parseWords.add(new ParseWord(chemicalName, wellFormedParseTokens));
		}
		else{
			for (List<ParseTokens> omittedWordParseTokens : omittedWordParseTokensList) {
				parseWords.add(new ParseWord(StringTools.stringListToString(omittedWordParseTokens.get(0).getTokens(), ""), omittedWordParseTokens));
			}
		}
		return parseWords;
	}
}
