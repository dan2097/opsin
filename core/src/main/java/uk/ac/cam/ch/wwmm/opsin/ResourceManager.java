package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import dk.brics.automaton.RunAutomaton;

import nu.xom.Document;
import nu.xom.Element;
import nu.xom.Elements;

/**Holds all of the tokens used in parsing of chemical names.
 * Holds all automata
 * Generates XML Elements for tokens.
 *
 * @author ptc24/dl387
 *
 */
class ResourceManager {

	/**Used to load XML files.*/
	private final ResourceGetter resourceGetter;
	
	/**A mapping between primitive tokens, and annotation->Token object mappings.*/
	final HashMap<String, HashMap<Character, Token>> tokenDict = new HashMap<String, HashMap<Character, Token>>();
	/**A mapping between regex tokens, and annotation->Token object mappings.*/
	final HashMap<Character, Token> reSymbolTokenDict = new HashMap<Character, Token>();


	/**A mapping between annotation symbols and the first two letters of token names applicable to
	 * that annotation symbol which then map to token names (annotation->first two letters of token names ->token names mapping).*/
	final List<HashMap<String, List<String>>> symbolTokenNamesDict = new ArrayList<HashMap<String, List<String>>>(200);
	/**A mapping between annotation symbols and DFAs (annotation->automata mapping).*/
	final List<List<RunAutomaton>> symbolRegexAutomataDict = new ArrayList<List<RunAutomaton>>();
	/**A mapping between annotation symbols and regex patterns (annotation->regex pattern mapping).*/
	final List<List<Pattern>> symbolRegexesDict = new ArrayList<List<Pattern>>();
	
	/**The automaton which describes the grammar of a chemical name from left to right*/
	final RunAutomaton chemicalAutomaton;
	
	
	/**As symbolTokenNamesDict but use the last two letters of token names to map to the token names.*/
	List<HashMap<String, List<String>>> symbolTokenNamesDict_TokensByLastTwoLetters;
	/**As symbolRegexAutomataDict but automata are reversed */
	List<List<RunAutomaton>> symbolRegexAutomataDictReversed;
	/**As symbolRegexesDict but regexes match the end of string */
	List<List<Pattern>> symbolRegexesDictReversed;
	
	/**The automaton which describes the grammar of a chemical name from right to left*/
	RunAutomaton reverseChemicalAutomaton;

	/**Generates the ResourceManager.
	 * @param resourceGetter
	 *
	 * @throws Exception If the XML token and regex files can't be read in properly.
	 */
	ResourceManager(ResourceGetter resourceGetter) throws Exception {
		this.resourceGetter = resourceGetter;
		chemicalAutomaton = processChemicalGrammar(false);
		int grammarSymbolsSize = chemicalAutomaton.getCharIntervals().length;
		for (int i = 0; i < grammarSymbolsSize; i++) {//initialise arrayLists
			symbolTokenNamesDict.add(null);
			symbolRegexAutomataDict.add(null);
			symbolRegexesDict.add(null);
		}
		processTokenFiles(false);
		processRegexTokenFiles(false);
	}

	/**
	 * Processes tokenFiles
	 * @param reversed Should the hashing of
	 * @throws Exception
	 */
	private void processTokenFiles(boolean reversed) throws Exception{
		Document tokenFiles = resourceGetter.getXMLDocument("index.xml");
		Elements files = tokenFiles.getRootElement().getChildElements("tokenFile");
		for(int i=0;i<files.size();i++) {
			Element rootElement = resourceGetter.getXMLDocument(files.get(i).getValue()).getRootElement();
			List<Element> tokenLists;
			if (rootElement.getLocalName().equals("tokenLists")){//support for xml files with one "tokenList" or multiple "tokenList" under a "tokenLists" element
				tokenLists = XOMTools.getChildElementsWithTagName(rootElement, "tokenList");
			}
			else{
				tokenLists =new ArrayList<Element>();
				tokenLists.add(rootElement);
			}
			for (Element tokenList : tokenLists) {
				char symbol = tokenList.getAttributeValue("symbol").charAt(0);
				List<Element> tokenElements = XOMTools.getChildElementsWithTagName(tokenList, "token");
				int index = Arrays.binarySearch(chemicalAutomaton.getCharIntervals(), symbol);
				if (index < 0){
					throw new Exception(symbol +" is associated with a tokenList of tagname " + tokenList.getAttributeValue("tagname") +" however it is not actually used in OPSIN's grammar!!!");
				}
				for (Element tokenElement : tokenElements) {
					String t = tokenElement.getValue();

					if(!tokenDict.containsKey(t)) {
						tokenDict.put(t, new HashMap<Character, Token>());
					}
					tokenDict.get(t).put(symbol, new Token(tokenElement, tokenList));
					if (!reversed){
						if(symbolTokenNamesDict.get(index)==null) {
							symbolTokenNamesDict.set(index, new HashMap<String, List<String>>());
						}
						if(!symbolTokenNamesDict.get(index).containsKey(t.substring(0, 2))) {
							symbolTokenNamesDict.get(index).put(t.substring(0, 2), new ArrayList<String>());
						}
						symbolTokenNamesDict.get(index).get(t.substring(0, 2)).add(t);
					}
					else{
						if(symbolTokenNamesDict_TokensByLastTwoLetters.get(index)==null) {
							symbolTokenNamesDict_TokensByLastTwoLetters.set(index, new HashMap<String, List<String>>());
						}
						if(!symbolTokenNamesDict_TokensByLastTwoLetters.get(index).containsKey(t.substring(t.length()-2))) {
							symbolTokenNamesDict_TokensByLastTwoLetters.get(index).put(t.substring(t.length()-2), new ArrayList<String>());
						}
						symbolTokenNamesDict_TokensByLastTwoLetters.get(index).get(t.substring(t.length()-2)).add(t);
					}
				}
			}
		}
	}

	private void processRegexTokenFiles(boolean reversed) throws Exception {
		Element reTokenList = resourceGetter.getXMLDocument("regexTokens.xml").getRootElement();
		Elements regexEls = reTokenList.getChildElements();
	
		HashMap<String, String> tempRegexes = new HashMap<String, String>();
		Pattern matchRegexReplacement = Pattern.compile("%.*?%");
		for(int i=0;i<regexEls.size();i++) {
			Element regexEl = regexEls.get(i);
			String re = regexEl.getAttributeValue("regex");
			Matcher m = matchRegexReplacement.matcher(re);
			String newValue = "";
			int position = 0;
			while(m.find()) {//replace sections enclosed in %..% with the appropriate regex
				newValue += re.substring(position, m.start());
				if (tempRegexes.get(m.group())==null){
					throw new Exception("Regex entry for: " + m.group() + " missing! Check regexTokens.xml");
				}
				newValue += tempRegexes.get(m.group());
				position = m.end();
			}
			newValue += re.substring(position);
			if (regexEl.getLocalName().equals("regex")){
				if (regexEl.getAttribute("name")==null){
					throw new Exception("Regex entry in regexTokenes.xml with no name. regex: " + newValue);
				}
				tempRegexes.put(regexEl.getAttributeValue("name"), newValue);
				continue;
			}
			//must be a regexToken
	
			Character symbol = regexEl.getAttributeValue("symbol").charAt(0);
			reSymbolTokenDict.put(symbol, new Token(regexEl));
	
			int index = Arrays.binarySearch(chemicalAutomaton.getCharIntervals(), symbol);
			if (index < 0){
				throw new Exception(symbol +" is associated with the regex " + newValue +" however it is not actually used in OPSIN's grammar!!!");
			}
			if (!reversed){
				if (regexEl.getAttribute("determinise")!=null){//should the regex be compiled into a DFA for faster execution?
					if(symbolRegexAutomataDict.get(index)==null) {
						symbolRegexAutomataDict.set(index, new ArrayList<RunAutomaton>());
					}
					symbolRegexAutomataDict.get(index).add(AutomatonInitialiser.getAutomaton(regexEl.getAttributeValue("tagname")+"_"+(int)symbol, newValue, false, false));
				}
				else{
					if(symbolRegexesDict.get(index)==null) {
						symbolRegexesDict.set(index, new ArrayList<Pattern>());
					}
					symbolRegexesDict.get(index).add(Pattern.compile(newValue));
				}
			}
			else{
				if (regexEl.getAttribute("determinise")!=null){//should the regex be compiled into a DFA for faster execution?
					if(symbolRegexAutomataDictReversed.get(index)==null) {
						symbolRegexAutomataDictReversed.set(index, new ArrayList<RunAutomaton>());
					}
					symbolRegexAutomataDictReversed.get(index).add(AutomatonInitialiser.getAutomaton(regexEl.getAttributeValue("tagname")+"_"+(int)symbol, newValue, false, true));
				}
				else{
					if(symbolRegexesDictReversed.get(index)==null) {
						symbolRegexesDictReversed.set(index, new ArrayList<Pattern>());
					}
					symbolRegexesDictReversed.get(index).add(Pattern.compile(newValue +"$"));
				}
			}
		}
	}
	
	private RunAutomaton processChemicalGrammar(boolean reversed) throws Exception {
		Map<String, String> regexDict = new HashMap<String, String>();
		Elements regexes = resourceGetter.getXMLDocument("regexes.xml").getRootElement().getChildElements("regex");
		Pattern matchRegexReplacement = Pattern.compile("%.*?%");
		for(int i=0;i<regexes.size();i++) {
			String name = regexes.get(i).getAttributeValue("name");
			String value = regexes.get(i).getAttributeValue("value");
			Matcher m = matchRegexReplacement.matcher(value);
			String newValue = "";
			int position = 0;
			while(m.find()) {
				newValue += value.substring(position, m.start());
				if (regexDict.get(m.group())==null){
					throw new Exception("Regex entry for: " + m.group() + " missing! Check regexes.xml");
				}
				newValue += regexDict.get(m.group());
				position = m.end();
			}
			newValue += value.substring(position);
			if (regexDict.get(name)!=null){
				throw new Exception("Regex entry: " + name + " has duplicate definitions! Check regexes.xml");
			}
			regexDict.put(name, newValue);
		}
		String re = regexDict.get("%chemical%");
		if (!reversed){
			return AutomatonInitialiser.getAutomaton("chemical", re, true, false);
		}
		else{
			return AutomatonInitialiser.getAutomaton("chemical", re, true, true);
		}
	}

	synchronized void populatedReverseTokenMappings() throws Exception{
		if (reverseChemicalAutomaton ==null){
			reverseChemicalAutomaton = processChemicalGrammar(true);
		}
		int grammarSymbolsSize = reverseChemicalAutomaton.getCharIntervals().length;
		if (symbolTokenNamesDict_TokensByLastTwoLetters ==null){
			symbolTokenNamesDict_TokensByLastTwoLetters = new ArrayList<HashMap<String,List<String>>>();
			for (int i = 0; i < grammarSymbolsSize; i++) {
				symbolTokenNamesDict_TokensByLastTwoLetters.add(null);
			}
			processTokenFiles(true);
		}
		if (symbolRegexAutomataDictReversed ==null && symbolRegexesDictReversed==null){
			symbolRegexAutomataDictReversed = new ArrayList<List<RunAutomaton>>();
			symbolRegexesDictReversed = new ArrayList<List<Pattern>>();
			for (int i = 0; i < grammarSymbolsSize; i++) {
				symbolRegexAutomataDictReversed.add(null);
				symbolRegexesDictReversed.add(null);
			}
			processRegexTokenFiles(true);
		}
	}

	/**Given a token string and an annotation character, makes the XML element for
	 * the token string.
	 * @param token The token string.
	 * @param symbol The annotation character.
	 *
	 * @return The XML element produced.
	 * @throws ParsingException
	 */
	Element makeTokenElement(String token, Character symbol) throws ParsingException {
		if(tokenDict.get(token) != null && tokenDict.get(token).get(symbol) !=null ) {
			return tokenDict.get(token).get(symbol).makeElement(token);
		}
		if (reSymbolTokenDict.get(symbol)!=null){
			return reSymbolTokenDict.get(symbol).makeElement(token);
		}
		throw new ParsingException("Parsing Error: This is a bug in the program. A token element could not be found for token: " + token +" using annotation symbol: " +symbol);
	}
}
