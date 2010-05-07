package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Pattern;

import uk.ac.cam.ch.wwmm.opsin.ParseWord.WordType;

/**
 * Uses OPSIN's DFA based grammar to break a name in to tokens with associated meanings ("annotations").
 * @author dl387
 *
 */
class Tokeniser {


	private final ParseRules parseRules;
	private final static char endOfFunctionalTerm = '\u00FB';
	private final Pattern matchCommaSpace = Pattern.compile(", ");
	private final Pattern matchSpace = Pattern.compile(" ");
	private final Pattern matchAcid = Pattern.compile("acid[\\]\\)\\}]*");
	
	Tokeniser(ParseRules parseRules) {
		this.parseRules = parseRules;
	}

	/**
	 * Master method for tokenizing chemical names into words and within words into tokens
	 * @param name The chemical name.
	 * @param allowRemovalOfWhiteSpace 
	 * @return
	 * @throws ParsingException 
	 */
	Parse tokenize(String name, boolean allowRemovalOfWhiteSpace) throws ParsingException {
		Parse parse = new Parse(name);
		String unparsedName;
		if (allowRemovalOfWhiteSpace){
			unparsedName =  removeWhiteSpaceIfBracketsAreUnbalanced(name);
		}
		else{
			unparsedName = name;
		}
		while (unparsedName.length()>0){
			/*
			 * Returns
			 * List of parses where at least some of the name was assigned a role
			 * Section of name that was uninterpretable (or "" if none was)
			 * Section of name that was unparsable (or "" if none was). This is always shorter or equal to the above string
			 */
			ParseRulesResults results = parseRules.getParses(unparsedName);
			List<ParseTokens> parseTokens =results.getParseTokensList();
			String uninterpretableName = results.getUninterpretableName();
			String parsedName = unparsedName.substring(0, unparsedName.length() - uninterpretableName.length());
			if (parseTokens.size()>0 && (uninterpretableName.equals("") || uninterpretableName.charAt(0) ==' ')){//a word was interpretable
				//If something like ethylchloride is encountered this should be split back to ethyl chloride and there will be 2 ParseWords returned
				//In cases of properly formed names there will be only one ParseWord
				//If there are two parses one of which assumes a missing space and one of which does not the former is discarded
				List<ParseWord> parseWords = splitIntoParseWords(parseTokens, parsedName);
				for (int j = 0; j < parseWords.size(); j++) {
					parse.addWord(parseWords.get(j));
				}
				if (!uninterpretableName.equals("")){
					unparsedName = uninterpretableName.substring(1);//remove white space at start of uninterpretableName
				}
				else{
					unparsedName = uninterpretableName;
				}
			}
			else{//word is unparsable as is.
				if (allowRemovalOfWhiteSpace){
					//Try and remove a space and try again
					//TODO add a warning message if this code is invoked. A name invoking this is unambiguously BAD
					int indexOfSpace = uninterpretableName.indexOf(' ');
					if (indexOfSpace != -1 ){
						unparsedName = parsedName + uninterpretableName.substring(0, indexOfSpace) + uninterpretableName.substring(indexOfSpace +1);
					}
					else{
						if (parsedName.equals("")){
							throw new ParsingException(name + " is unparsable due to the following word being unparsable: " + unparsedName+ " (spaces will have been removed as they are assumed to be erroneous if parsing fails with them there)");
						}
						else {
							throw new ParsingException(name + " is unparsable. This part of the name was interpretable: " + parsedName);
						}
					}
				}
				else{
					if (parsedName.equals("")){
						throw new ParsingException(name + " is unparsable due to the following word being unparsable: " + unparsedName);
					}
					else {
						throw new ParsingException(name + " is unparsable. This part of the name was interpretable: " + parsedName);
					}
				}
			}
		}
		return parse;
	}
	
	/**
	 * Works left to right removing spaces if there are too many opening brackets
	 * @param name
	 * @return
	 * @throws ParsingException If brackets are unbalanced and cannot be balanced by removing whitespace
	 */
	private String removeWhiteSpaceIfBracketsAreUnbalanced(String name) throws ParsingException {
		int bracketLevel = 0;
		int stringLength  = name.length();
		for(int i = 0 ; i < stringLength; i++) {
			char c = name.charAt(i);
			if(c == '(' || c == '[' || c == '{') {
				bracketLevel++;
			}
			else if(c == ')' || c == ']' || c == '}') {
				bracketLevel--;
			}
			else if(c == ' ' && bracketLevel > 0){//brackets unbalanced and a space has been encountered!
				name = name.substring(0, i) +name.substring(i +1);
				stringLength  = name.length();
				i--;
			}
		}
		if (bracketLevel > 0){
			throw new ParsingException("Unmatched opening bracket found in :" + name);
		}
		else if (bracketLevel < 0){
			throw new ParsingException("Unmatched closing bracket found in :" + name);
		}
		return name;
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
	
	/**
	 * Inverts a CAS name.
	 * Throws an exception is OPSIN is unable to determine whether something is a substituent or functional term
	 * or if something unexpected in a CAS name is encountered
	 * @param name
	 * @return
	 * @throws ParsingException
	 */
	String uninvertCASName(String name) throws ParsingException {
		List<String> nameComponents = new ArrayList<String>(Arrays.asList(matchCommaSpace.split(name)));
		List<String> substituents = new ArrayList<String>();
		List<String> seperateWordSubstituents = new ArrayList<String>();
		List<String> functionalTerms = new ArrayList<String>();
		
		String parent = nameComponents.get(0);
		String[] parentNameParts = matchSpace.split(parent);
		if (parentNameParts.length !=1){
			if (parentNameParts.length >2 || !matchAcid.matcher(parentNameParts[1]).matches()){
				throw new ParsingException("Invalid CAS name. Parent compound was followed by an unexpected term");
			}
		}
		boolean addedBracket = false;
		boolean esterEncountered = false;
		for (int i = 1; i < nameComponents.size(); i++) {
			String[] components = matchSpace.split(nameComponents.get(i));
			for (String component : components) {
				if (component.endsWith("-")){
					if (isCloseBracketMissing(component)){
						if (addedBracket){
							throw new ParsingException("Close bracket bracket appears to be missing");
						}
						parent += "]";
						addedBracket =true;
					}
					substituents.add(component);
				}
				else{
					ParseRulesResults results = parseRules.getParses(component);
					List<ParseTokens> parseTokens = results.getParseTokensList();
					if (parseTokens.size() >0){
						if (splitIntoParseWords(parseTokens, component).size()>1){
							throw new ParsingException("Missing space found in name prevents interpetation as CAS index name");
						}
						WordType wordType = OpsinTools.determineWordType(parseTokens.get(0).getAnnotations());
						for (int j = 1; j < parseTokens.size(); j++) {
							if (!wordType.equals(OpsinTools.determineWordType(parseTokens.get(j).getAnnotations()))){
								throw new ParsingException(component + "can be interpeted in multiple ways. For the sake of precision OPSIN has decided not to process this as a CAS name");
							}
						}
						if (wordType.equals(WordType.functionalTerm)){
							if (component.equals("ester")){
								if (esterEncountered){
									throw new ParsingException("ester formation was mentioned more than once in CAS name!");
								}
								parent = uninvertEster(parent);
								esterEncountered=true;
							}
							else{
								functionalTerms.add(component);
							}
						}
						else if (wordType.equals(WordType.substituent)){
							seperateWordSubstituents.add(component);
						}
						else if (wordType.equals(WordType.full)){
							throw new ParsingException("Unable to interpret: " + component +" (as part of a CAS index name)- A full word was encountered where a substituent or functionalTerm was expected");
						}
					}
					else{
						throw new ParsingException("Unable to interpret: " + component +" (as part of a CAS index name)");
					}
				}
			}
		}
		StringBuilder casName = new StringBuilder();
		for (String prefixFunctionalTerm : seperateWordSubstituents) {
			casName.append(prefixFunctionalTerm);
			casName.append(" ");
		}
		for (String substituent : substituents) {
			casName.append(substituent);
		}
		casName.append(parent);
		for (String functionalTerm : functionalTerms) {
			casName.append(" ");
			casName.append(functionalTerm);
		}
		return casName.toString();
	}

	/**
	 * Modifies the name of the parent acid from ic to ate (or ous to ite)
	 * hence allowing the formation of the uninverted ester
	 * @param parent
	 * @return
	 * @throws ParsingException
	 */
	private String uninvertEster(String parent) throws ParsingException {
		int len = parent.length();
		if (len <9){
			throw new ParsingException("Failed to uninvert CAS ester");
		}
		char lastChar = parent.charAt(len-1);
		if (lastChar ==')' || lastChar ==']' || lastChar =='}'){
			if (parent.substring(parent.length()-8).equalsIgnoreCase("ic acid)")){
				parent = parent.substring(0, parent.length()-8) + "ate)";
			}
			else if (parent.substring(parent.length()-9).equalsIgnoreCase("ous acid)")){
				parent = parent.substring(0, parent.length()-9) + "ite)";
			}
			else{
				throw new ParsingException("Failed to uninvert CAS ester");
			}
		}
		else{
			if (parent.substring(parent.length()-7).equalsIgnoreCase("ic acid")){
				parent = parent.substring(0, parent.length()-7) + "ate";
			}
			else if (parent.substring(parent.length()-8).equalsIgnoreCase("ous acid")){
				parent = parent.substring(0, parent.length()-8) + "ite";
			}
			else{
				throw new ParsingException("Failed to uninvert CAS ester");
			}
		}
		return parent;
	}

	private boolean isCloseBracketMissing(String component) {
		char[] characters = component.toCharArray();
		for (int i = characters.length -1; i >=0; i--) {
			char character = characters[i];
			if (character =='(' || character =='[' || character =='{'){
				return true;
			}
			if (character ==')' || character ==']' || character =='}'){
				return false;
			}
		}
		return false;
	}
}
