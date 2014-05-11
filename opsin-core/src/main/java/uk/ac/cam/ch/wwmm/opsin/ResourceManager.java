package uk.ac.cam.ch.wwmm.opsin;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.xml.stream.XMLStreamConstants;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;

import nu.xom.Elements;
import dk.brics.automaton.RunAutomaton;
import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

/**Holds all of the tokens used in parsing of chemical names.
 * Holds all automata
 * Generates XML Elements for tokens.
 *
 * @author ptc24
 * @author dl387
 *
 */
class ResourceManager {
	private final static Element IGNORE_WHEN_WRITING_PARSE_TREE = new Element("");

	/**Used to load XML files.*/
	private final ResourceGetter resourceGetter;
	
	/**Used to serialise and deserialise automata.*/
	private final AutomatonInitialiser automatonInitialiser;
	
	/**A mapping between primitive tokens, and annotation->Token object mappings.*/
	final HashMap<String, Map<Character, Element>> tokenDict = new HashMap<String, Map<Character, Element>>();
	/**A mapping between regex tokens, and annotation->Token object mappings.*/
	final HashMap<Character, Element> reSymbolTokenDict = new HashMap<Character, Element>();


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
	 * @param reversed Should the tokens be reversed
	 */
	private void processTokenFiles(boolean reversed) {
		try {
			XMLStreamReader filesToProcessReader = resourceGetter.getXMLDocument2("index.xml");
			while (filesToProcessReader.hasNext()) {
				int event = filesToProcessReader.next();
				if (event == XMLStreamConstants.START_ELEMENT && 
						filesToProcessReader.getLocalName().equals("tokenFile")) {
					String fileName = filesToProcessReader.getElementText();
					processTokenFile(fileName, reversed);
				}
			}
			filesToProcessReader.close();
		}
		catch (IOException e){
			throw new NameToStructureException(e.getMessage(), e);
		}
		catch (XMLStreamException e){
			throw new NameToStructureException(e.getMessage(), e);
		}
	}

	private void processTokenFile(String fileName, boolean reversed) throws IOException, XMLStreamException {
		XMLStreamReader reader = resourceGetter.getXMLDocument2(fileName);
		while (reader.hasNext()) {
			switch (reader.next()) {
			case XMLStreamConstants.START_ELEMENT:
				String tagName = reader.getLocalName();
				if (tagName.equals("tokenLists")) {
					while (reader.hasNext()) {
						switch (reader.next()) {
						case XMLStreamConstants.START_ELEMENT:
							if (reader.getLocalName().equals("tokenList")) {
								processTokenList(reader, reversed);
							}
							break;
						}
					}
				}
				else if (tagName.equals("tokenList")) {
					processTokenList(reader, reversed);
				}
				break;
			}
		}
		reader.close();
	}

	private void processTokenList(XMLStreamReader reader, boolean reversed) throws XMLStreamException {
		String tokenTagName = null;
		Character symbol = null;
		String type = null;
		String subType = null;
		boolean ignoreWhenWritingXML = false;
		
		for (int i = 0, l = reader.getAttributeCount(); i < l; i++) {
			String atrName = reader.getAttributeLocalName(i);
			String atrValue = reader.getAttributeValue(i);
			if (atrName.equals("tagname")){
				tokenTagName  = atrValue;
			}
			else if (atrName.equals("symbol")){
				symbol = atrValue.charAt(0);
			}
			else if (atrName.equals(TYPE_ATR)){
				type = atrValue;
			}
			else if (atrName.equals(SUBTYPE_ATR)){
				subType = atrValue;
			}
			else if (atrName.equals("ignoreWhenWritingXML")){
				ignoreWhenWritingXML = atrValue.equals("yes");
			}
			else{
				throw new RuntimeException("Malformed tokenlist");
			}
		}
		if (tokenTagName == null || symbol == null) {
			throw new RuntimeException("Malformed tokenlist");
		}
		
		int index = Arrays.binarySearch(chemicalAutomaton.getCharIntervals(), symbol);
		if (index < 0){
			throw new RuntimeException(symbol +" is associated with a tokenList of tagname " + tokenTagName +" however it is not actually used in OPSIN's grammar!!!");
		}
		
		while (reader.hasNext()) {
			switch (reader.next()) {
			case XMLStreamConstants.START_ELEMENT:
				if (reader.getLocalName().equals("token")) {
					Element el;
					if (ignoreWhenWritingXML){
						el = IGNORE_WHEN_WRITING_PARSE_TREE;
					}
					else{
						el = new Element(tokenTagName);
						if (type != null){
							el.addAttribute(TYPE_ATR, type);
						}
						if (subType != null){
							el.addAttribute(SUBTYPE_ATR, subType);
						}
						for (int i = 0, l = reader.getAttributeCount(); i < l; i++) {
							el.addAttribute(reader.getAttributeLocalName(i), reader.getAttributeValue(i));
						}
					}
					String t = reader.getElementText();
					Map<Character, Element> symbolToToken = tokenDict.get(t);
					if(symbolToToken == null) {
						symbolToToken = new HashMap<Character, Element>();
						tokenDict.put(t, symbolToToken);
					}
					symbolToToken.put(symbol, el);

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
				break;
			case XMLStreamConstants.END_ELEMENT:
				if (reader.getLocalName().equals("tokenList")) {
					return;
				}
				break;
			}
		}
	}

	private void processRegexTokenFiles(boolean reversed) throws IOException{
		XMLStreamReader reader = resourceGetter.getXMLDocument2("regexTokens.xml");
		Map<String, String> tempRegexes = new HashMap<String, String>();
		Pattern matchRegexReplacement = Pattern.compile("%.*?%");
		try {
			while (reader.hasNext()) {
				if (reader.next() == XMLStreamConstants.START_ELEMENT) {
					String localName = reader.getLocalName();
					if (!localName.equals("regex") && !localName.equals("regexToken")){
						continue;
					}
					String re = reader.getAttributeValue(null, "regex");
					Matcher m = matchRegexReplacement.matcher(re);
					StringBuilder newValueSB = new StringBuilder();
					int position = 0;
					while(m.find()) {//replace sections enclosed in %..% with the appropriate regex
						newValueSB.append(re.substring(position, m.start()));
						if (tempRegexes.get(m.group()) == null){
							throw new RuntimeException("Regex entry for: " + m.group() + " missing! Check regexTokens.xml");
						}
						newValueSB.append(tempRegexes.get(m.group()));
						position = m.end();
					}
					newValueSB.append(re.substring(position));
					if (localName.equals("regex")) {
						String regexName = reader.getAttributeValue(null, "name");
						if (regexName == null){
							throw new RuntimeException("Regex entry in regexTokenes.xml with no name. regex: " + newValueSB.toString());
						}
						tempRegexes.put(regexName, newValueSB.toString());
						continue;
					}
					addRegexToken(reader, newValueSB.toString(), reversed);
				}
			}
		}
		catch (XMLStreamException e) {
			throw new IOException("Parsing exception occurred while reading regexTokens.xml", e);
		}
		finally {
			try {
				reader.close();
			} catch (XMLStreamException e) {
				throw new IOException("Parsing exception occurred while reading regexTokens.xml", e);
			}
		}
	}
	
	private void addRegexToken(XMLStreamReader reader, String regex, boolean reversed) {
		String tokenTagName = null;
		Character symbol = null;
		String type = null;
		String subType = null;
		String value = null;
		boolean determinise = false;
		boolean ignoreWhenWritingXML = false;
		
		for (int i = 0, l = reader.getAttributeCount(); i < l; i++) {
			String atrName = reader.getAttributeLocalName(i);
			String atrValue = reader.getAttributeValue(i);
			if (atrName.equals("tagname")){
				tokenTagName  = atrValue;
			}
			else if (atrName.equals("symbol")){
				symbol = atrValue.charAt(0);
			}
			else if (atrName.equals(TYPE_ATR)){
				type = atrValue;
			}
			else if (atrName.equals(SUBTYPE_ATR)){
				subType = atrValue;
			}
			else if (atrName.equals("value")){
				value = atrValue;
			}
			else if (atrName.equals("determinise")){
				determinise = atrValue.equals("yes");
			}
			else if (atrName.equals("ignoreWhenWritingXML")){
				ignoreWhenWritingXML = atrValue.equals("yes");
			}
			else if (!atrName.equals("regex")){
				throw new RuntimeException("Malformed regexToken");
			}
		}
		if (tokenTagName == null || symbol == null) {
			throw new RuntimeException("Malformed regexToken");
		}
		
		if (!reversed) {
			//reSymbolTokenDict will be populated when the constructor is called for left-right parsing, hence skip for right-left 
			if (reSymbolTokenDict.get(symbol) != null) {
				throw new RuntimeException(symbol +" is associated with multiple regular expressions. The following expression clashes: " + regex +" This should be resolved by combining regular expressions that map the same symbol" );
			}

			if (ignoreWhenWritingXML) {
				reSymbolTokenDict.put(symbol, IGNORE_WHEN_WRITING_PARSE_TREE);
			}
			else{
				Element el = new Element(tokenTagName);
				if (type != null){
					el.addAttribute(TYPE_ATR, type);
				}
				if (subType != null){
					el.addAttribute(SUBTYPE_ATR, subType);
				}
				if (value != null){
					el.addAttribute(VALUE_ATR, value);
				}
				reSymbolTokenDict.put(symbol, el);
			}
		}
		
		int index = Arrays.binarySearch(chemicalAutomaton.getCharIntervals(), symbol);
		if (index < 0){
			throw new RuntimeException(symbol +" is associated with the regex " + regex +" however it is not actually used in OPSIN's grammar!!!");
		}
		if (!reversed){
			if (determinise){//should the regex be compiled into a DFA for faster execution?
				symbolRegexAutomataDict[index] = automatonInitialiser.loadAutomaton(tokenTagName + "_" + (int)symbol, regex, false, false);
			}
			else{
				symbolRegexesDict[index] = Pattern.compile(regex);
			}
		}
		else{
			if (determinise){//should the regex be compiled into a DFA for faster execution?
				symbolRegexAutomataDictReversed[index] = automatonInitialiser.loadAutomaton(tokenTagName + "_" + (int)symbol, regex, false, true);
			}
			else{
				symbolRegexesDictReversed[index] = Pattern.compile(regex +"$");
			}
		}
	}

	private RunAutomaton processChemicalGrammar(boolean reversed) throws IOException{
		Map<String, StringBuilder> regexDict = new HashMap<String, StringBuilder>();
		Elements regexes = resourceGetter.getXMLDocument("regexes.xml").getRootElement().getChildElements("regex");
		Pattern matchRegexReplacement = Pattern.compile("%.*?%");
		for(int i = 0, l =regexes.size(); i < l; i++) {
			nu.xom.Element regex =regexes.get(i);
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
		Map<Character, Element> annotationToToken = tokenDict.get(tokenString);
		if(annotationToToken != null){
			Element token = annotationToToken.get(symbol);
			if (token != null) {
				if (token == IGNORE_WHEN_WRITING_PARSE_TREE){
					return null;
				}
				Element tokenInstance = new Element(token);
				tokenInstance.setValue(tokenString);
				return tokenInstance;
			}
		}
		Element regexToken = reSymbolTokenDict.get(symbol);
		if (regexToken != null){
			if (regexToken == IGNORE_WHEN_WRITING_PARSE_TREE){
				return null;
			}
			Element tokenInstance = new Element(regexToken);
			tokenInstance.setValue(tokenString);
			return tokenInstance;
		}
		throw new ParsingException("Parsing Error: This is a bug in the program. A token element could not be found for token: " + tokenString +" using annotation symbol: " +symbol);
	}
	
}
