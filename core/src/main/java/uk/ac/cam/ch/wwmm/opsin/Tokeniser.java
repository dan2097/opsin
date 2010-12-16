package uk.ac.cam.ch.wwmm.opsin;

import java.util.Collections;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Uses OPSIN's DFA based grammar to break a name into tokens with associated meanings ("annotations").
 * @author dl387
 *
 */
class Tokeniser {
	private final ParseRules parseRules;
	private final Pattern matchCasCollectiveIndex = Pattern.compile("([\\[\\(\\{]([1-9][0-9]?[cC][iI][, ]?)+[\\]\\)\\}])+|[1-9][0-9]?[cC][iI]", Pattern.CASE_INSENSITIVE );
	private final Pattern matchCompoundWithPhrase = Pattern.compile("(compd\\. with|compound with|and) ", Pattern.CASE_INSENSITIVE );

	Tokeniser(ParseRules parseRules) {
		this.parseRules = parseRules;
	}

	ParseRules getParseRules() {
		return parseRules;
	}

	/**
	 * Master method for tokenizing chemical names into words and within words into tokens
	 * @param name The chemical name.
	 * @param allowRemovalOfWhiteSpace 
	 * @return
	 * @throws ParsingException 
	 */
	TokenizationResult tokenize(String name, boolean allowRemovalOfWhiteSpace) throws ParsingException {
		TokenizationResult result = allowRemovalOfWhiteSpace ? new TokenizationResult(WordTools.removeWhiteSpaceIfBracketsAreUnbalanced(name)) :new TokenizationResult(name);

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
		TokenizationResult result = new TokenizationResult(name);
		//removeWhiteSpaceIfBracketsAreUnbalanced is not currently employed as this the input to this function from the parser will often be what the LR tokenizer couldn't handle, which may not have matching brackets
	
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
		return parseTokens.size()>0 && (result.isFullyInterpretable() || result.getUninterpretableName().charAt(result.getUninterpretableName().length()-1)==' ' || result.getUninterpretableName().charAt(result.getUninterpretableName().length()-1) =='-');
	}

	private boolean isWordParsable(List<ParseTokens> parseTokens, TokenizationResult result) {
		return parseTokens.size()>0 && (result.isFullyInterpretable() || result.getUninterpretableName().charAt(0) ==' ' || result.getUninterpretableName().charAt(0) =='-');
	}
	
	private void parseWord(TokenizationResult result, List<ParseTokens> parseTokens, String parsedName, boolean reverse) throws ParsingException {
		//If something like ethylchloride is encountered this should be split back to ethyl chloride and there will be 2 ParseWords returned
		//In cases of properly formed names there will be only one ParseWord
		//If there are two parses one of which assumes a missing space and one of which does not the former is discarded
		addParseWords(parseTokens, parsedName, result.getParse(), reverse);

		if (result.isFullyInterpretable()) {
			result.setUnparsedName(result.getUninterpretableName());
		} else {
			String name = reverse ? result.getUninterpretableName().substring(0, result.getUninterpretableName().length() - 1) : result.getUninterpretableName().substring(1);
			result.setUnparsedName(name);
		}
	}

	private void addParseWords(List<ParseTokens> parseTokens, String parsedName, Parse parse, boolean reverse) throws ParsingException {
		List<ParseWord> parseWords = WordTools.splitIntoParseWords(parseTokens, parsedName);

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

	/**
	 * Fixes cases like for example "benzene sulfonamide" -->"benzenesulfonamide"
	 * @param parsedWords
	 * @param result
	 * @return
	 * @throws ParsingException
	 */
	private boolean reverseSpaceRemoval(List<ParseWord> parsedWords, TokenizationResult result) throws ParsingException {
		boolean successful = false;
	
		if (!parsedWords.isEmpty()) {//first see whether the space before the unparseable word is erroneous
			ParseWord pw = parsedWords.get(parsedWords.size() - 1);
			ParseRulesResults backResults = parseRules.getParses(pw.getWord() + result.getUnparsedName());
			List<ParseTokens> backParseTokens = backResults.getParseTokensList();
			String backUninterpretableName = backResults.getUninterpretableName();
			String backParsedName = pw.getWord() + result.getUnparsedName().substring(0, result.getUnparsedName().length() - backUninterpretableName.length());
			if (backParsedName.length() > pw.getWord().length() && backParseTokens.size() > 0 && (backUninterpretableName.equals("") || backUninterpretableName.charAt(0) == ' ' || backUninterpretableName.charAt(0) == '-')) {//a word was interpretable
				result.getParse().removeWord(pw);
				List<ParseWord> parseWords = WordTools.splitIntoParseWords(backParseTokens, backParsedName);
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
}
