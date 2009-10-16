package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import nu.xom.Document;
import nu.xom.Element;
import nu.xom.Elements;
import nu.xom.Text;

/**Holds all of the tokens used in parsing of chemical names.
 * Generates XML Elements for tokens.
 *
 * @author ptc24
 *
 */
class TokenManager {

	/**A mapping between primitive tokens, and annotation->Token object mappings.*/
	HashMap<String, HashMap<Character, Token>> tokenDict;
	/**A mapping between regex tokens, and annotation->Token object mappings.*/
	HashMap<Character, Token> reSymbolTokenDict;


	/**A mapping between annotation symbols and the first two letters of token names applicable to
	 * that annotation symbol which then map to token names (annotation->2 letters of token names ->token names mapping).*/
	HashMap<Character, HashMap<String, List<String>>> symbolTokenNamesDict;
	/**A mapping between annotation symbols and regex patterns (annotation->regex pattern mapping).*/
	HashMap<Character, List<Pattern>> symbolRegexesDict;

	/**Generates the TokenManager.
	 * @param resourceGetter 
	 *
	 * @throws Exception If the XML token and regex files can't be read in properly.
	 */
	TokenManager(ResourceGetter resourceGetter) throws Exception {
		tokenDict = new HashMap<String, HashMap<Character, Token>>();
		reSymbolTokenDict = new HashMap<Character, Token>();

		symbolTokenNamesDict = new HashMap<Character, HashMap<String, List<String>>>();
		symbolRegexesDict = new HashMap<Character, List<Pattern>>();

		Document tokenFiles = resourceGetter.getXMLDocument("index.xml");
		Elements files = tokenFiles.getRootElement().getChildElements("tokenFile");
		for(int i=0;i<files.size();i++) {
			Element rootElement = resourceGetter.getXMLDocument(((Text)files.get(i).getChild(0)).getValue()).getRootElement();
			List<Element> tokenLists =new ArrayList<Element>();
			if (rootElement.getLocalName().equals("tokenLists")){//support for xml files with one "tokenList" or multiple "tokenList" under a "tokenLists" element
				Elements children =rootElement.getChildElements();
				for (int j = 0; j <children.size(); j++) {
					tokenLists.add(children.get(j));
				}
			}
			else{
				tokenLists.add(rootElement);
			}
			for (Element tokenList : tokenLists) {
				char symbol = tokenList.getAttributeValue("symbol").charAt(0);
				Elements tokenElements = tokenList.getChildElements("token");
				for(int j=0;j<tokenElements.size();j++) {
					String t = ((Text)tokenElements.get(j).getChild(0)).getValue();

					if(!tokenDict.containsKey(t)) {
						tokenDict.put(t, new HashMap<Character, Token>());
					}
					tokenDict.get(t).put(symbol, new Token(tokenElements.get(j), tokenList));

					if(!symbolTokenNamesDict.containsKey(symbol)) {
						symbolTokenNamesDict.put(symbol, new HashMap<String, List<String>>());
					}
					if(!symbolTokenNamesDict.get(symbol).containsKey(t.substring(0, 2))) {
						symbolTokenNamesDict.get(symbol).put(t.substring(0, 2), new ArrayList<String>());
					}
					symbolTokenNamesDict.get(symbol).get(t.substring(0, 2)).add(t);
				}
			}
		}

		Element reTokenList = resourceGetter.getXMLDocument("regexTokens.xml").getRootElement();
		Elements reTokens = reTokenList.getChildElements("regexToken");

		HashMap<String, String> tempRegexes = new HashMap<String, String>();
		Pattern p = Pattern.compile("%.*?%");
		for(int i=0;i<reTokens.size();i++) {
			Element rt = reTokens.get(i);
			String re = rt.getAttributeValue("regex");
			Matcher m = p.matcher(re);
			String newValue = "";
			int position = 0;
			while(m.find()) {//replace sections enclosed in %..% with the appropriate regex
				newValue += re.substring(position, m.start());
				if (tempRegexes.get(m.group())==null){
					throw new ParsingException("Regex entry for: " + m.group() + " missing! Check regexTokens.xml");
				}
				newValue += tempRegexes.get(m.group());
				position = m.end();
			}
			newValue += re.substring(position);
			if (rt.getAttribute("tagname") ==null){
				if (rt.getAttribute("name")==null){
					throw new Exception("Entry in regexTokenes.xml has neither a tagname or a name. regex: " + newValue);
				}
				tempRegexes.put(rt.getAttributeValue("name"), newValue);
				continue;
			}

			Character symbol = rt.getAttributeValue("symbol").charAt(0);
			reSymbolTokenDict.put(symbol, new Token(rt.getAttributeValue("tagname"), rt.getAttributeValue("type"), rt.getAttributeValue("ignoreWhenWritingXML")));

			if(!symbolRegexesDict.containsKey(symbol)) {
				symbolRegexesDict.put(symbol, new ArrayList<Pattern>());
			}
			symbolRegexesDict.get(symbol).add(Pattern.compile(newValue));
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
