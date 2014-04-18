package uk.ac.cam.ch.wwmm.opsin;

import java.io.IOException;
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
 * @author ptc24
 * @author dl387
 *
 */
class ResourceManager {

	/**Used to load XML files.*/
	private final ResourceGetter resourceGetter;
	
	/**Used to serialise and deserialise automata.*/
	private final AutomatonInitialiser automatonInitialiser;
	
	/**A mapping between primitive tokens, and annotation->Token object mappings.*/
	final HashMap<String, Map<Character, Token>> tokenDict = new HashMap<String, Map<Character, Token>>();
	/**A mapping between regex tokens, and annotation->Token object mappings.*/
	final HashMap<Character, Token> reSymbolTokenDict = new HashMap<Character, Token>();


	/**A mapping between annotation symbols and a trie of tokens.*/
	final OpsinRadixTrie[] symbolTokenNamesDict;
	/**A mapping between annotation symbols and DFAs (annotation->automata mapping).*/
	final RunAutomaton[] symbolRegexAutomataDict;
	/**A mapping between annotation symbols and regex patterns (annotation->regex pattern mapping).*/
	final Pattern[] symbolRegexesDict;
	
	/**The automaton which describes the grammar of a chemical name from left to right*/
	final RunAutomaton chemicalAutomaton;
	
	
	/**As symbolTokenNamesDict but the tokens are reversed*/
	OpsinRadixTrie[] symbolTokenNamesDictReversed;
	/**As symbolRegexAutomataDict but automata are reversed */
	RunAutomaton[] symbolRegexAutomataDictReversed;
	/**As symbolRegexesDict but regexes match the end of string */
	Pattern[] symbolRegexesDictReversed;
	
	/**The automaton which describes the grammar of a chemical name from right to left*/
	RunAutomaton reverseChemicalAutomaton;

	/**Generates the ResourceManager.
	 * This involves reading in the token files, the regexToken file (regexTokens.xml) and the grammar file (regexes.xml).
	 * DFA are built or retrieved for the regexTokens and the chemical grammar.
	 * 
	 * Throws an exception if the XML token and regex files can't be read in properly or the grammar cannot be built.
	 * @param resourceGetter
	 * @throws IOException 
	 */
	ResourceManager(ResourceGetter resourceGetter) throws IOException{
		this.resourceGetter = resourceGetter;
		this.automatonInitialiser = new AutomatonInitialiser(resourceGetter.getResourcePath() + "serialisedAutomata/");
		chemicalAutomaton = processChemicalGrammar(false);
		int grammarSymbolsSize = chemicalAutomaton.getCharIntervals().length;
		symbolTokenNamesDict = new OpsinRadixTrie[grammarSymbolsSize];
		symbolRegexAutomataDict = new RunAutomaton[grammarSymbolsSize];
		symbolRegexesDict = new Pattern[grammarSymbolsSize];
		processTokenFiles(false);
		processRegexTokenFiles(false);
	}

	/**
	 * Processes tokenFiles
	 * @param reversed Should the hashing of
	 * @throws IOException 
	 */
	private void processTokenFiles(boolean reversed) throws IOException{
		Document tokenFiles = resourceGetter.getXMLDocument("index.xml");
		Elements files = tokenFiles.getRootElement().getChildElements("tokenFile");
		for(int i = 0, l = files.size(); i < l; i++) {
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
					throw new RuntimeException(symbol +" is associated with a tokenList of tagname " + tokenList.getAttributeValue("tagname") +" however it is not actually used in OPSIN's grammar!!!");
				}
				for (Element tokenElement : tokenElements) {
					String t = tokenElement.getValue();

					Map<Character, Token> symbolToToken = tokenDict.get(t);
					if(symbolToToken == null) {
						symbolToToken = new HashMap<Character, Token>();
						tokenDict.put(t, symbolToToken);
					}
					symbolToToken.put(symbol, new Token(tokenElement, tokenList));
					if (!reversed){
						if(symbolTokenNamesDict[index]==null) {
							symbolTokenNamesDict[index] = new OpsinRadixTrie();
						}
						symbolTokenNamesDict[index].addToken(t);
					}
					else{
						if(symbolTokenNamesDictReversed[index]==null) {
							symbolTokenNamesDictReversed[index] = new OpsinRadixTrie();
						}
						symbolTokenNamesDictReversed[index].addToken(new StringBuilder(t).reverse().toString());
					}
				}
			}
		}
	}

	private void processRegexTokenFiles(boolean reversed) throws IOException{
		Element reTokenList = resourceGetter.getXMLDocument("regexTokens.xml").getRootElement();
		Elements regexEls = reTokenList.getChildElements();
	
		HashMap<String, String> tempRegexes = new HashMap<String, String>();
		Pattern matchRegexReplacement = Pattern.compile("%.*?%");
		for(int i = 0, l = regexEls.size(); i < l; i++) {
			Element regexEl = regexEls.get(i);
			String re = regexEl.getAttributeValue("regex");
			Matcher m = matchRegexReplacement.matcher(re);
			StringBuilder newValueSB = new StringBuilder();
			int position = 0;
			while(m.find()) {//replace sections enclosed in %..% with the appropriate regex
				newValueSB.append(re.substring(position, m.start()));
				if (tempRegexes.get(m.group())==null){
					throw new RuntimeException("Regex entry for: " + m.group() + " missing! Check regexTokens.xml");
				}
				newValueSB.append(tempRegexes.get(m.group()));
				position = m.end();
			}
			newValueSB.append(re.substring(position));
			if (regexEl.getLocalName().equals("regex")){
				if (regexEl.getAttribute("name")==null){
					throw new RuntimeException("Regex entry in regexTokenes.xml with no name. regex: " + newValueSB.toString());
				}
				tempRegexes.put(regexEl.getAttributeValue("name"), newValueSB.toString());
				continue;
			}
			
			Character symbol = regexEl.getAttributeValue("symbol").charAt(0);
			if (!reversed) {
				//reSymbolTokenDict will be populated when the constructor is called for left-right parsing, hence skip for right-left 
				if (reSymbolTokenDict.get(symbol) != null) {
					throw new RuntimeException(symbol +" is associated with multiple regular expressions. The following expression clashes: " + regexEl.toXML() +" This should be resolved by combining regular expressions that map the same symbol" );
				}
				reSymbolTokenDict.put(symbol, new Token(regexEl));
			}
			
			int index = Arrays.binarySearch(chemicalAutomaton.getCharIntervals(), symbol);
			if (index < 0){
				throw new RuntimeException(symbol +" is associated with the regex " + newValueSB.toString() +" however it is not actually used in OPSIN's grammar!!!");
			}
			if (!reversed){
				if (regexEl.getAttribute("determinise") != null){//should the regex be compiled into a DFA for faster execution?
					symbolRegexAutomataDict[index] = automatonInitialiser.loadAutomaton(regexEl.getAttributeValue("tagname")+"_"+(int)symbol, newValueSB.toString(), false, false);
				}
				else{
					symbolRegexesDict[index] = Pattern.compile(newValueSB.toString());
				}
			}
			else{
				if (regexEl.getAttribute("determinise")!=null){//should the regex be compiled into a DFA for faster execution?
					symbolRegexAutomataDictReversed[index] = automatonInitialiser.loadAutomaton(regexEl.getAttributeValue("tagname")+"_"+(int)symbol, newValueSB.toString(), false, true);
				}
				else{
					symbolRegexesDictReversed[index] = Pattern.compile(newValueSB.toString() +"$");
				}
			}
		}
	}
	
	private RunAutomaton processChemicalGrammar(boolean reversed) throws IOException{
		Map<String, StringBuilder> regexDict = new HashMap<String, StringBuilder>();
		Elements regexes = resourceGetter.getXMLDocument("regexes.xml").getRootElement().getChildElements("regex");
		Pattern matchRegexReplacement = Pattern.compile("%.*?%");
		for(int i = 0, l =regexes.size(); i < l; i++) {
			Element regex =regexes.get(i);
			String name = regex.getAttributeValue("name");
			String value = regex.getAttributeValue("value");
			Matcher m = matchRegexReplacement.matcher(value);
			StringBuilder newValueSB = new StringBuilder();
			int position = 0;
			while(m.find()) {
				newValueSB.append(value.substring(position, m.start()));
				if (regexDict.get(m.group()) == null){
					throw new RuntimeException("Regex entry for: " + m.group() + " missing! Check regexes.xml");
				}
				newValueSB.append(regexDict.get(m.group()));
				position = m.end();
			}
			newValueSB.append(value.substring(position));
			if (regexDict.get(name) != null){
				throw new RuntimeException("Regex entry: " + name + " has duplicate definitions! Check regexes.xml");
			}
			regexDict.put(name, newValueSB);
		}
		String re = regexDict.get("%chemical%").toString();
		if (!reversed){
			return automatonInitialiser.loadAutomaton("chemical", re, true, false);
		}
		else{
			return automatonInitialiser.loadAutomaton("chemical", re, true, true);
		}
	}

	synchronized void populatedReverseTokenMappings() throws IOException{
		if (reverseChemicalAutomaton == null){
			reverseChemicalAutomaton = processChemicalGrammar(true);
		}
		int grammarSymbolsSize = reverseChemicalAutomaton.getCharIntervals().length;
		if (symbolTokenNamesDictReversed == null){
			symbolTokenNamesDictReversed = new OpsinRadixTrie[grammarSymbolsSize];
			processTokenFiles(true);
		}
		if (symbolRegexAutomataDictReversed == null && symbolRegexesDictReversed==null){
			symbolRegexAutomataDictReversed = new RunAutomaton[grammarSymbolsSize];
			symbolRegexesDictReversed = new Pattern[grammarSymbolsSize];
			processRegexTokenFiles(true);
		}
	}

	/**Given a token string and an annotation character, makes the XML element for
	 * the token string.
	 * @param tokenString The token string.
	 * @param symbol The annotation character.
	 *
	 * @return The XML element produced.
	 * @throws ParsingException
	 */
	Element makeTokenElement(String tokenString, Character symbol) throws ParsingException {
		Map<Character, Token> annotationToToken = tokenDict.get(tokenString);
		if(annotationToToken != null){
			Token token = annotationToToken.get(symbol);
			if (token != null ) {
				return token.makeElement(tokenString);
			}
		}
		Token regexToken = reSymbolTokenDict.get(symbol);
		if (regexToken != null){
			return regexToken.makeElement(tokenString);
		}
		throw new ParsingException("Parsing Error: This is a bug in the program. A token element could not be found for token: " + tokenString +" using annotation symbol: " +symbol);
	}
}
