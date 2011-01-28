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
	
	/**A mapping between primitive tokens, and annotation->Token object mappings.*/
	final HashMap<String, HashMap<Character, Token>> tokenDict = new HashMap<String, HashMap<Character, Token>>();
	/**A mapping between regex tokens, and annotation->Token object mappings.*/
	final HashMap<Character, Token> reSymbolTokenDict = new HashMap<Character, Token>();


	/**A mapping between annotation symbols and the first two letters of token names applicable to
	 * that annotation symbol which then map to token names (annotation->first two letters of token names ->token names mapping).*/
	final HashMap<String, List<String>>[] symbolTokenNamesDict;
	/**A mapping between annotation symbols and DFAs (annotation->automata mapping).*/
	final List<RunAutomaton>[] symbolRegexAutomataDict;
	/**A mapping between annotation symbols and regex patterns (annotation->regex pattern mapping).*/
	final List<Pattern>[] symbolRegexesDict;
	
	/**The automaton which describes the grammar of a chemical name from left to right*/
	final RunAutomaton chemicalAutomaton;
	
	
	/**As symbolTokenNamesDict but use the last two letters of token names to map to the token names.*/
	HashMap<String, List<String>>[] symbolTokenNamesDict_TokensByLastTwoLetters;
	/**As symbolRegexAutomataDict but automata are reversed */
	List<RunAutomaton>[] symbolRegexAutomataDictReversed;
	/**As symbolRegexesDict but regexes match the end of string */
	List<Pattern>[] symbolRegexesDictReversed;
	
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
	@SuppressWarnings("unchecked")
	ResourceManager(ResourceGetter resourceGetter) throws IOException{
		this.resourceGetter = resourceGetter;
		chemicalAutomaton = processChemicalGrammar(false);
		int grammarSymbolsSize = chemicalAutomaton.getCharIntervals().length;
		symbolTokenNamesDict = new HashMap[grammarSymbolsSize];
		symbolRegexAutomataDict = new List[grammarSymbolsSize];
		symbolRegexesDict = new List[grammarSymbolsSize];
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
					throw new RuntimeException(symbol +" is associated with a tokenList of tagname " + tokenList.getAttributeValue("tagname") +" however it is not actually used in OPSIN's grammar!!!");
				}
				for (Element tokenElement : tokenElements) {
					String t = tokenElement.getValue();

					if(!tokenDict.containsKey(t)) {
						tokenDict.put(t, new HashMap<Character, Token>());
					}
					tokenDict.get(t).put(symbol, new Token(tokenElement, tokenList));
					if (!reversed){
						if(symbolTokenNamesDict[index]==null) {
							symbolTokenNamesDict[index] = new HashMap<String, List<String>>();
						}
						if(!symbolTokenNamesDict[index].containsKey(t.substring(0, 2))) {
							symbolTokenNamesDict[index].put(t.substring(0, 2), new ArrayList<String>());
						}
						symbolTokenNamesDict[index].get(t.substring(0, 2)).add(t);
					}
					else{
						if(symbolTokenNamesDict_TokensByLastTwoLetters[index]==null) {
							symbolTokenNamesDict_TokensByLastTwoLetters[index] = new HashMap<String, List<String>>();
						}
						if(!symbolTokenNamesDict_TokensByLastTwoLetters[index].containsKey(t.substring(t.length()-2))) {
							symbolTokenNamesDict_TokensByLastTwoLetters[index].put(t.substring(t.length()-2), new ArrayList<String>());
						}
						symbolTokenNamesDict_TokensByLastTwoLetters[index].get(t.substring(t.length()-2)).add(t);
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
		for(int i=0;i<regexEls.size();i++) {
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
			//must be a regexToken
	
			Character symbol = regexEl.getAttributeValue("symbol").charAt(0);
			reSymbolTokenDict.put(symbol, new Token(regexEl));
	
			int index = Arrays.binarySearch(chemicalAutomaton.getCharIntervals(), symbol);
			if (index < 0){
				throw new RuntimeException(symbol +" is associated with the regex " + newValueSB.toString() +" however it is not actually used in OPSIN's grammar!!!");
			}
			if (!reversed){
				if (regexEl.getAttribute("determinise")!=null){//should the regex be compiled into a DFA for faster execution?
					if(symbolRegexAutomataDict[index]==null) {
						symbolRegexAutomataDict[index] = new ArrayList<RunAutomaton>();
					}
					symbolRegexAutomataDict[index].add(AutomatonInitialiser.loadAutomaton(regexEl.getAttributeValue("tagname")+"_"+(int)symbol, newValueSB.toString(), false, false));
				}
				else{
					if(symbolRegexesDict[index]==null) {
						symbolRegexesDict[index] = new ArrayList<Pattern>();
					}
					symbolRegexesDict[index].add(Pattern.compile(newValueSB.toString()));
				}
			}
			else{
				if (regexEl.getAttribute("determinise")!=null){//should the regex be compiled into a DFA for faster execution?
					if(symbolRegexAutomataDictReversed[index]==null) {
						symbolRegexAutomataDictReversed[index] = new ArrayList<RunAutomaton>();
					}
					symbolRegexAutomataDictReversed[index].add(AutomatonInitialiser.loadAutomaton(regexEl.getAttributeValue("tagname")+"_"+(int)symbol, newValueSB.toString(), false, true));
				}
				else{
					if(symbolRegexesDictReversed[index]==null) {
						symbolRegexesDictReversed[index] = new ArrayList<Pattern>();
					}
					symbolRegexesDictReversed[index].add(Pattern.compile(newValueSB.toString() +"$"));
				}
			}
		}
	}
	
	private RunAutomaton processChemicalGrammar(boolean reversed) throws IOException{
		Map<String, String> regexDict = new HashMap<String, String>();
		Elements regexes = resourceGetter.getXMLDocument("regexes.xml").getRootElement().getChildElements("regex");
		Pattern matchRegexReplacement = Pattern.compile("%.*?%");
		for(int i=0;i<regexes.size();i++) {
			String name = regexes.get(i).getAttributeValue("name");
			String value = regexes.get(i).getAttributeValue("value");
			Matcher m = matchRegexReplacement.matcher(value);
			StringBuilder newValueSB = new StringBuilder();
			int position = 0;
			while(m.find()) {
				newValueSB.append(value.substring(position, m.start()));
				if (regexDict.get(m.group())==null){
					throw new RuntimeException("Regex entry for: " + m.group() + " missing! Check regexes.xml");
				}
				newValueSB.append(regexDict.get(m.group()));
				position = m.end();
			}
			newValueSB.append(value.substring(position));
			if (regexDict.get(name)!=null){
				throw new RuntimeException("Regex entry: " + name + " has duplicate definitions! Check regexes.xml");
			}
			regexDict.put(name, newValueSB.toString());
		}
		String re = regexDict.get("%chemical%");
		if (!reversed){
			return AutomatonInitialiser.loadAutomaton("chemical", re, true, false);
		}
		else{
			return AutomatonInitialiser.loadAutomaton("chemical", re, true, true);
		}
	}

	@SuppressWarnings("unchecked")
	synchronized void populatedReverseTokenMappings() throws IOException{
		if (reverseChemicalAutomaton ==null){
			reverseChemicalAutomaton = processChemicalGrammar(true);
		}
		int grammarSymbolsSize = reverseChemicalAutomaton.getCharIntervals().length;
		if (symbolTokenNamesDict_TokensByLastTwoLetters ==null){
			symbolTokenNamesDict_TokensByLastTwoLetters  = new HashMap[grammarSymbolsSize];
			processTokenFiles(true);
		}
		if (symbolRegexAutomataDictReversed ==null && symbolRegexesDictReversed==null){
			symbolRegexAutomataDictReversed = new List[grammarSymbolsSize];
			symbolRegexesDictReversed = new List[grammarSymbolsSize];
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
