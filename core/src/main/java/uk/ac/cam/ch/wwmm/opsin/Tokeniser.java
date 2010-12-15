package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import uk.ac.cam.ch.wwmm.opsin.ParseWord.WordType;

/**
 * Uses OPSIN's DFA based grammar to break a name in to tokens with associated meanings ("annotations").
 * @author dl387
 *
 */
class Tokeniser {
	private final ParseRules parseRules;
	private final static char endOfSubstituent = '\u00e9';
	private final static char endOfMainGroup = '\u00e2';
	private final static char endOfFunctionalTerm = '\u00FB';
	private final Pattern matchCommaSpace = Pattern.compile(", ");
	private final Pattern matchSpace = Pattern.compile(" ");
	private final Pattern matchAcid = Pattern.compile("acid[\\]\\)\\}]*");
	private final Pattern matchCasCollectiveIndex = Pattern.compile("([\\[\\(\\{]([1-9][0-9]?[cC][iI][, ]?)+[\\]\\)\\}])+|[1-9][0-9]?[cC][iI]", Pattern.CASE_INSENSITIVE );
	private final Pattern matchCompoundWithPhrase = Pattern.compile("(compd\\. with|compound with|and) ", Pattern.CASE_INSENSITIVE );

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
	TokenizationResult tokenize(String name, boolean allowRemovalOfWhiteSpace) throws ParsingException {
		TokenizationResult result = new TokenizationResult(name, allowRemovalOfWhiteSpace);

		while (!result.isSuccessfullyTokenized()){
			ParseRulesResults results = parseRules.getParses(result.getUnparsedName());
			List<ParseTokens> parseTokens = results.getParseTokensList();
			result.setUninterpretableName(results.getUninterpretableName());
			String parsedName = result.getParsedName();

			if (isWordParsable(parseTokens, result)){
				parseWord(result, parseTokens, parsedName, false);
			}
			else if (!fixWord(result, parsedName, results, allowRemovalOfWhiteSpace)) {
				break;
			}
		}

		return result;
	}
	
	/**
	 * Master method for tokenizing chemical names into words and within words into tokens
	 * This is performed in a right to left manner
	 * @param reverseParseRules
	 * @param name The chemical name.
	 * @param allowRemovalOfWhiteSpace 
	 * @return
	 * @throws ParsingException 
	 */
	TokenizationResult tokenizeRightToLeft(ReverseParseRules reverseParseRules, String name, boolean allowRemovalOfWhiteSpace) throws ParsingException {
		TokenizationResult result = new TokenizationResult(name, allowRemovalOfWhiteSpace);

		//bracket matching is not currently being performed as this the input to this function from the parser will often be what the LR tokenizer couldn't handle, which may not have matching brackets
	
		while (!result.isSuccessfullyTokenized()){
			ParseRulesResults results = reverseParseRules.getParses(result.getUnparsedName());
			List<ParseTokens> parseTokens =results.getParseTokensList();
			result.setUninterpretableName(results.getUninterpretableName());
			String parsedName = result.getUnparsedName().substring(result.getUninterpretableName().length());

			if (isWordParsableInReverse(parseTokens, result)) {
				parseWord(result, parseTokens, parsedName, true);
			}
			else if (!fixWordInReverse(result, parsedName, results, allowRemovalOfWhiteSpace)) {
				break;
			}
		}
		
		Collections.reverse(result.getParse().getWords());

		return result;
	}


	private boolean isWordParsableInReverse(List<ParseTokens> parseTokens, TokenizationResult result) {
		return parseTokens.size()>0 && (result.hasUninterpretableName() || result.getUninterpretableName().charAt(result.getUninterpretableName().length()-1)==' ' || result.getUninterpretableName().charAt(result.getUninterpretableName().length()-1) =='-');
	}

	private boolean isWordParsable(List<ParseTokens> parseTokens, TokenizationResult result) {
		return parseTokens.size()>0 && (result.hasUninterpretableName() || result.getUninterpretableName().charAt(0) ==' ' || result.getUninterpretableName().charAt(0) =='-');
	}
	
	private void parseWord(TokenizationResult result, List<ParseTokens> parseTokens, String parsedName, boolean reverse) throws ParsingException {
		//If something like ethylchloride is encountered this should be split back to ethyl chloride and there will be 2 ParseWords returned
		//In cases of properly formed names there will be only one ParseWord
		//If there are two parses one of which assumes a missing space and one of which does not the former is discarded
		addParseWords(parseTokens, parsedName, result.getParse(), reverse);

		if (result.hasUninterpretableName()) {
			result.setUnparsedName(result.getUninterpretableName());
		} else {
			String name = reverse ? result.getUninterpretableName().substring(0, result.getUninterpretableName().length() - 1) : result.getUninterpretableName().substring(1);
			result.setUnparsedName(name);
		}
	}

	private boolean reverseSpaceRemoval(List<ParseWord> parsedWords, TokenizationResult result) throws ParsingException {
		boolean successful = false;

		if (!parsedWords.isEmpty()) {//first see whether the space before the unparseable word is erroneous
			ParseWord pw = parsedWords.get(parsedWords.size() - 1);
			ParseRulesResults backResults = parseRules.getParses(pw.getWord() + result.getUnparsedName());
			List<ParseTokens> backParseTokens = backResults.getParseTokensList();
			String backUninterpretableName = backResults.getUninterpretableName();
			String backParsedName = pw.getWord() + result.getUnparsedName().substring(0, result.getUnparsedName().length() - backUninterpretableName.length());
			if (backParsedName.length() > pw.getWord().length() && backParseTokens.size() > 0 && (backUninterpretableName.equals("") || backUninterpretableName.charAt(0) == ' ')) {//a word was interpretable
				result.getParse().removeWord(pw);
				List<ParseWord> parseWords = splitIntoParseWords(backParseTokens, backParsedName);
				for (ParseWord parseWord : parseWords) {
					result.getParse().addWord(parseWord);
				}
				if (!backUninterpretableName.equals("")) {
					result.setUnparsedName(backUninterpretableName.substring(1));//remove white space at start of uninterpretableName
				} else {
					result.setUnparsedName(backUninterpretableName);
				}
				successful = true;
			}
		}

		return successful;
	}

	private void addParseWords(List<ParseTokens> parseTokens, String parsedName, Parse parse, boolean reverse) throws ParsingException {
		List<ParseWord> parseWords = splitIntoParseWords(parseTokens, parsedName);

		if (reverse) {
			Collections.reverse(parseWords);//make this set of words back to front as well
		}

		for (ParseWord parseWord : parseWords) {
			parse.addWord(parseWord);
		}
	}

	private boolean fixWord(TokenizationResult result, String parsedName, ParseRulesResults results, boolean allowRemovalOfWhiteSpace) throws ParsingException {
		Matcher m = matchCompoundWithPhrase.matcher(result.getUninterpretableName());
		if (m.lookingAt()) {
			result.setUnparsedName(parsedName + result.getUninterpretableName().substring(m.group().length()));
		} else if (matchCasCollectiveIndex.matcher(result.getUninterpretableName()).matches()) {
			result.setUnparsedName(parsedName);
		} else {
			if (allowRemovalOfWhiteSpace) {
				//TODO add a warning message if this code is invoked. A name invoking this is unambiguously BAD
				List<ParseWord> parsedWords = result.getParse().getWords();
				if (!reverseSpaceRemoval(parsedWords, result)) {
					//Try and remove a space from the right and try again
					int indexOfSpace = result.getUninterpretableName().indexOf(' ');
					if (indexOfSpace != -1) {
						result.setUnparsedName( parsedName + result.getUninterpretableName().substring(0, indexOfSpace) + result.getUninterpretableName().substring(indexOfSpace + 1));
					} else {
						return false;
					}
				}
			} else {
				result.setUnparsableName(results.getUnparseableName());
				return false;
			}
		}
		return true;
	}

	private boolean fixWordInReverse(TokenizationResult result, String parsedName, ParseRulesResults results, boolean allowRemovalOfWhiteSpace) {
		if (allowRemovalOfWhiteSpace) {
			//Try and remove a space and try again
			//TODO add a warning message if this code is invoked. A name invoking this is unambiguously BAD
			int indexOfSpace = result.getUninterpretableName().lastIndexOf(' ');
			if (indexOfSpace != -1) {
				result.setUnparsedName( result.getUninterpretableName().substring(0, indexOfSpace) + result.getUninterpretableName().substring(indexOfSpace + 1) + parsedName);
			} else {
				return false;
			}
		} else {
			result.setUnparsableName(results.getUnparseableName());
			return false;
		}
		return true;
	}

	private List<ParseWord> splitIntoParseWords(List<ParseTokens> parseTokensList, String chemicalName) throws ParsingException {
		List<ParseTokens> wellFormedParseTokens = new ArrayList<ParseTokens>();//these are all in the same word as would be expected
		List<List<ParseTokens>> omittedWordParseTokensList = new ArrayList<List<ParseTokens>>();//these are grouped into words e.g. ethylchloride will have a list of parseTokens for the ethyl and chloride
		omittedWordParseTokensList.add(new ArrayList<ParseTokens>());
		omittedWordParseTokensList.add(new ArrayList<ParseTokens>());//only 1 space is allowed to be omitted
		int longestFunctionalTermEncountered =0;//we want the longest functional term
		int shortestNonFunctionalTermEncountered = Integer.MAX_VALUE;//and the shortest non functional term
		for (ParseTokens parseTokens : parseTokensList) {
			List<Character> annotations = parseTokens.getAnnotations();
			List<List<Character>> chunkedAnnotations = chunkAnnotations(annotations);//chunked into mainGroup/substituent/functionalTerm
			if (chunkedAnnotations.size()>1 && annotations.contains(endOfFunctionalTerm)) {//must be an omitted space as not allowed to have a functionalTerm and anything else
				List<String> tokens = parseTokens.getTokens();
				List<Character> newAnnotations = new ArrayList<Character>();
				List<String> newTokens = new ArrayList<String>();
				int annotPos = 0;
				int wordCounter = 0;
				for (List<Character> annotationList : chunkedAnnotations) {
					boolean functionalTermNext = false;
					if (annotationList.get(annotationList.size()-1).equals(endOfFunctionalTerm)) {
						functionalTermNext = true;
						if (newAnnotations.size()>0){//create a new parseTokens, unless nothing has been read yet e.g. in the case of poly
							ParseTokens newParseTokens = new ParseTokens(newTokens, newAnnotations);
							if (wordCounter >=2){
								throw new ParsingException("Name appears to have 2 or more omitted spaces!");
							}
							int currentNonFunctionalTermLength = StringTools.stringListToString(newTokens, "").length();
							if (currentNonFunctionalTermLength <= shortestNonFunctionalTermEncountered  && !omittedWordParseTokensList.get(wordCounter).contains(newParseTokens)){
								if (currentNonFunctionalTermLength < shortestNonFunctionalTermEncountered) {
									omittedWordParseTokensList.get(wordCounter).clear();
									shortestNonFunctionalTermEncountered =currentNonFunctionalTermLength;
								}
								omittedWordParseTokensList.get(wordCounter).add(newParseTokens);
							}
							wordCounter++;
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
						ParseTokens newParseTokens = new ParseTokens(newTokens, newAnnotations);
						if (wordCounter >=2){
							throw new ParsingException("Name appears to have 2 or more omitted spaces!");
						}
						int currentFunctionalTermLength = StringTools.stringListToString(newTokens, "").length();
						if (currentFunctionalTermLength >= longestFunctionalTermEncountered && !omittedWordParseTokensList.get(wordCounter).contains(newParseTokens)){
							if (currentFunctionalTermLength > longestFunctionalTermEncountered) {
								omittedWordParseTokensList.get(wordCounter).clear();
								longestFunctionalTermEncountered =currentFunctionalTermLength;
							}
							omittedWordParseTokensList.get(wordCounter).add(newParseTokens);
						}
						wordCounter++;
						newAnnotations = new ArrayList<Character>();
						newTokens = new ArrayList<String>();
					}
				}
				if (!newAnnotations.isEmpty()) {
					ParseTokens newParseTokens = new ParseTokens(newTokens, newAnnotations);
					if (wordCounter >= 2){
						throw new ParsingException("Name appears to have 2 or more omitted spaces!");
					}
					int currentNonFunctionalTermLength = StringTools.stringListToString(newTokens, "").length();
					if (currentNonFunctionalTermLength <= shortestNonFunctionalTermEncountered  && !omittedWordParseTokensList.get(wordCounter).contains(newParseTokens)){
						if (currentNonFunctionalTermLength < shortestNonFunctionalTermEncountered) {
							omittedWordParseTokensList.get(wordCounter).clear();
							shortestNonFunctionalTermEncountered =currentNonFunctionalTermLength;
						}
						omittedWordParseTokensList.get(wordCounter).add(newParseTokens);
					}
					wordCounter++;
				}
			}
			else{
				wellFormedParseTokens.add(parseTokens);
			}
		}
		List<ParseWord> parseWords = new ArrayList<ParseWord>();
		if (!wellFormedParseTokens.isEmpty()) {
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
			if (matchCasCollectiveIndex.matcher(parentNameParts[parentNameParts.length-1]).matches()){//CAS collective index description should be ignored
				parent = "";
				for (int i = 0; i < parentNameParts.length-1; i++) {
					parent += parentNameParts[i];
				}
				parentNameParts = matchSpace.split(parent);
			}
			for (int i = 1; i < parentNameParts.length; i++) {
				if (!matchAcid.matcher(parentNameParts[i]).matches()){
					ParseRulesResults results = parseRules.getParses(parentNameParts[i]);
					List<ParseTokens> parseTokens = results.getParseTokensList();
					if (parseTokens.isEmpty()){
						throw new ParsingException("Invalid CAS name. Parent compound was followed by an unexpected term");
					}
				}
			}
		}
		boolean addedBracket = false;
		boolean esterEncountered = false;
		for (int i = 1; i < nameComponents.size(); i++) {
			String nameComponent = nameComponents.get(i);
			Matcher m = matchCompoundWithPhrase.matcher(nameComponent);
			boolean compoundWithcomponent =false;
			if (m.lookingAt()){
				nameComponent = nameComponent.substring(m.group().length());
				compoundWithcomponent =true;
			}
			String[] components = matchSpace.split(nameComponents.get(i));
			for (String component : components) {
				if (compoundWithcomponent){
					functionalTerms.add(component);
					continue;
				}
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
							if (component.equalsIgnoreCase("ester")){
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
							if (StringTools.endsWithCaseInsensitive(component, "ate") || StringTools.endsWithCaseInsensitive(component, "ite")//e.g. Piperazinium, 1,1-dimethyl-, 2,2,2-trifluoroacetate hydrochloride
									|| component.equalsIgnoreCase("hydrofluoride") || component.equalsIgnoreCase("hydrochloride") || component.equalsIgnoreCase("hydrobromide") || component.equalsIgnoreCase("hydroiodide")){
								functionalTerms.add(component);
							}
							else{
								throw new ParsingException("Unable to interpret: " + component +" (as part of a CAS index name)- A full word was encountered where a substituent or functionalTerm was expected");
							}
						}
					}
					else{
						if (!matchCasCollectiveIndex.matcher(component).matches()){//CAS collective index description should be ignored
							throw new ParsingException("Unable to interpret: " + component +" (as part of a CAS index name)");
						}
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
