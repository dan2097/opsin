package uk.ac.cam.ch.wwmm.opsin;

import java.io.IOException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Deque;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import nu.xom.Attribute;
import nu.xom.Element;
import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;
import static uk.ac.cam.ch.wwmm.opsin.OpsinTools.*;

/**Conducts finite-state parsing on chemical names.
 * Adds XML annotation to the semantic constituents of the name.
 *
 * @author ptc24/dl387
 *
 */
class Parser {

	/**Uses ParseRules to intelligently parse chemical names into substituents/full terms/functional terms*/
	private final Tokeniser tokeniser;
	/**The rules by which words are grouped together (e.g. in functional class nomenclature)*/
	private final WordRules wordRules;
	/**Holds the various tokens used.*/
	private final ResourceManager resourceManager;
	private final ParseRules parseRules;
	
	private final static Pattern matchSemiColonSpace = Pattern.compile("; ");
	private final static Pattern matchStoichiometryIndication = Pattern.compile("[ ]?[\\{\\[\\(](\\d+|\\?)([:/](\\d+|\\?))+[\\}\\]\\)]$");

	/**
	 * No-argument constructor. Uses ResouceGetter found at
	 * uk/ac/cam/ch/wwmm/opsin/resources/
	 * @throws IOException 
	 */
	Parser() throws IOException{
		ResourceGetter resources = new ResourceGetter("uk/ac/cam/ch/wwmm/opsin/resources/");
		this.wordRules = new WordRules(resources);
		this.resourceManager = new ResourceManager(resources);
		this.parseRules = new ParseRules(this.resourceManager);
		this.tokeniser = new Tokeniser(parseRules);
	}

	/**Initialises the parser.
	 * @param resourceManager
	 * @param tokeniser
	 * @param wordRules
	 */
	Parser(WordRules wordRules, Tokeniser tokeniser, ResourceManager resourceManager) {
		this.wordRules =wordRules;
		this.resourceManager = resourceManager;
		this.tokeniser = tokeniser;
		this.parseRules = tokeniser.getParseRules();
	}

	/**Parses a chemical name to an XML representation of the parse.
	 * @param n2sConfig 
	 *
	 * @param name The name to parse.
	 * @return The parse.
	 * @throws ParsingException If the name is unparsable.
	 */
	List<Element> parse(NameToStructureConfig n2sConfig, String name) throws ParsingException {
		Integer[] componentRatios = null;
		if (name.endsWith(")") || name.endsWith("]") || name.endsWith("}")){
			Matcher m = matchStoichiometryIndication.matcher(name);
			if (m.find()){
				componentRatios = processStoichiometryIndication(m.group());
				name = m.replaceAll("");
			}
		}
		Parse parse = null;
		if (name.contains(", ")){
			try{
				TokenizationResult tokenizationResult = tokeniser.tokenize(CASTools.uninvertCASName(name, parseRules), false);
				if (tokenizationResult.isSuccessfullyTokenized()){
					parse = tokenizationResult.getParse();
				}
			}
			catch (ParsingException ignored) {
			}
		}
		else if (name.contains("; ")){//a mixture, spaces are sufficient for OPSIN to treat as a mixture. These spaces for obvious reasons must not be removed
			TokenizationResult tokenizationResult = tokeniser.tokenize(matchSemiColonSpace.matcher(name).replaceAll(" "), false);
			if (tokenizationResult.isSuccessfullyTokenized()){
				parse = tokenizationResult.getParse();
			}
		}
		boolean allowSpaceRemoval = parse ==null ? true : false;
		if (parse == null){
			TokenizationResult tokenizationResult = tokeniser.tokenize(name , true);
			if (tokenizationResult.isSuccessfullyTokenized()){
				parse = tokenizationResult.getParse();
			}
			else{
				if (n2sConfig.isDetailedFailureAnalysis()){
					generateExactParseFailureReason(tokenizationResult, name);
				}
				else{
					throw new ParsingException(name + " is unparsable due to the following being uninterpretable: " + tokenizationResult.getUninterpretableName()
							+ " The following was not parseable: " +tokenizationResult.getUnparsableName());
				}
			}
		}
		
		List<Parse> parses = generateParseCombinations(parse);
		if (parses.size()==0){
			throw new ParsingException("No parses could be found for " + name);
		}
		
		List<Element> results = new ArrayList<Element>();
		for(Parse pp : parses) {
			Element moleculeEl = new Element(MOLECULE_EL);
			moleculeEl.addAttribute(new Attribute(NAME_ATR, name));
			for(ParseWord pw : pp.getWords()) {
				Element word = new Element(WORD_EL);
				moleculeEl.appendChild(word);
				if (pw.getParseTokens().size() >1){
					throw new ParsingException("OPSIN bug: parseWord had multiple annotations after creating additional parses step");
				}
				
				WordType wordType = OpsinTools.determineWordType(pw.getParseTokens().get(0).getAnnotations());
				word.addAttribute(new Attribute(TYPE_ATR, wordType.toString()));
				if (pw.getWord().startsWith("-")){//we want -acid to be the same as acid
					word.addAttribute(new Attribute(VALUE_ATR, pw.getWord().substring(1)));
				}
				else{
					word.addAttribute(new Attribute(VALUE_ATR, pw.getWord()));
				}
				for(ParseTokens pt : pw.getParseTokens()) {
					writeWordXML(word, pw, pt.getTokens(), WordTools.chunkAnnotations(pt.getAnnotations()));
				}
			}
			/* All words are placed into a wordRule.
			 * Often multiple words in the same wordRule.
			 * WordRules can be nested within each other e.g. in Carbonyl cyanide m-chlorophenyl hydrazone ->
			 * <wr><wr>Carbonyl cyanide</wr> m-chlorophenyl hydrazone </wr>
			 */
			try{
				wordRules.groupWordsIntoWordRules(n2sConfig, moleculeEl, allowSpaceRemoval);
				if (componentRatios!=null){
					applyStoichiometryIndicationToWordRules(moleculeEl, componentRatios);
				}
				results.add(moleculeEl);
			}
			catch (ParsingException e) {
				// Using that parse no word rules matched
			}
		}
		if (results.size()==0){
			throw new ParsingException(name + " could be parsed but OPSIN was unsure of the meaning of the words. This error will occur, by default, if a name is just a substituent");
		}
		
		return results;
	}

	static Integer[] processStoichiometryIndication(String ratioString) throws ParsingException {
		ratioString = ratioString.trim();
		ratioString = ratioString.substring(1, ratioString.length()-1);
		String[] ratioStrings = MATCH_COLON.split(ratioString);
		if (ratioStrings.length ==1){
			ratioStrings = MATCH_SLASH.split(ratioString);
		}
		Integer[] componentRatios = new Integer[ratioStrings.length];
		for (int i = 0; i < ratioStrings.length; i++) {
			String currentRatio = ratioStrings[i];
			if (currentRatio.contains("/")){
				throw new ParsingException("Unexpected / in component ratio declaration");
			}
			if (currentRatio.equals("?")){
				componentRatios[i]=1;
			}
			else{
				componentRatios[i]=Integer.parseInt(currentRatio);
			}
		}
		return componentRatios;
	}

	private void generateExactParseFailureReason(TokenizationResult tokenizationResult, String name) throws ParsingException {
		ReverseParseRules reverseParseRules;
		try {
			reverseParseRules = new ReverseParseRules(resourceManager);
		} catch (IOException e) {
			throw new RuntimeException("Failed to load resources for parsing names from right to left!",e);
		}
		String uninterpretableLR = tokenizationResult.getUninterpretableName();
		String unparseableLR = tokenizationResult.getUnparsableName();
		TokenizationResult reverseTokenizationResult = tokeniser.tokenizeRightToLeft(reverseParseRules, uninterpretableLR, true);
		String uninterpretableRL = reverseTokenizationResult.getUninterpretableName();
		String unparseableRL = reverseTokenizationResult.getUnparsableName();
		int indiceToTruncateUpTo =  uninterpretableLR.length()-unparseableLR.length();
		StringBuilder message = new StringBuilder();
		message.append(name);
		if (!uninterpretableRL.equals("")){
			message.append(" was uninterpretable due to the following section of the name: ");
			message.append(uninterpretableRL);
			if (indiceToTruncateUpTo <= unparseableRL.length()){
				String uninterpretableInContext = unparseableRL.substring(indiceToTruncateUpTo);
				if (!uninterpretableInContext.equals("")){
					message.append("  The following was not understandable in the context it was used: ");
					message.append(uninterpretableInContext);
				}
			}
		}
		else{
			message.append(" has no tokens unknown to OPSIN but does not conform to its grammar. ");
			message.append("From left to right it is unparsable due to the following being uninterpretable:");
			message.append(uninterpretableLR);
			message.append(" The following or which was not parseable: ");
			message.append(unparseableLR);
		}
		throw new ParsingException(message.toString());
	}

	/**
	 * For cases where any of the parse's parseWords contain multiple annotations create a
	 * parse for each possibility. Hence after this process there may be multiple Parse objects but
	 * the parseWords they contain will each only have one parseTokens object.
	 * @param parse
	 * @return
	 * @throws ParsingException 
	 */
	private List<Parse> generateParseCombinations(Parse parse) throws ParsingException {
		int numberOfCombinations = 1;
		List<Integer> parseCounts = new ArrayList<Integer>();
		List<ParseWord> parseWords = parse.getWords();
		for (ParseWord pw : parseWords) {
			int parsesForWord = pw.getParseTokens().size();
			parseCounts.add(parsesForWord);
			numberOfCombinations *= parsesForWord;
			if (numberOfCombinations > 128){//checked here to avoid integer overflow on inappropriate input
				throw new ParsingException("Too many different combinations of word interpretation are possible (>128) i.e. name contains too many terms that OPSIN finds ambiguous to interpret");
			}
		}
		if (numberOfCombinations ==1){
			return Arrays.asList(parse);
		}
		List<Parse> parses = new ArrayList<Parse>();
		
		Deque<Parse> parseQueue = new ArrayDeque<Parse>();
		parseQueue.add(new Parse(parse.getName()));
		while (!parseQueue.isEmpty()){
			Parse currentParse = parseQueue.removeFirst();
			int wordsInCurrentParse = currentParse.getWords().size();
			if(wordsInCurrentParse == parseWords.size()) {
				parses.add(currentParse);
			}
			else {
				ParseWord referenceWord = parseWords.get(wordsInCurrentParse);
				List<ParseTokens> referenceWordParseTokens = referenceWord.getParseTokens();
				for (int i = referenceWordParseTokens.size()-1; i >=0; i--) {
					ParseTokens parseTokens = referenceWordParseTokens.get(i);
					Parse parseWithNextWord = i > 0 ? currentParse.deepCopy() : currentParse;
					ParseWord newParseWord = new ParseWord(referenceWord.getWord(), Arrays.asList(parseTokens));
					parseWithNextWord.addWord(newParseWord);
					parseQueue.add(parseWithNextWord);
				}
			}
		}
		return parses;
	}

	/**Write the XML corresponding to a particular word in a parse.
	 *
	 * @param wordEl The empty XML word element to be written into.
	 * @param pw The ParseWord for the word.
	 * @param tokens The list of tokens.
	 * @param annotations The lists of annotations. This has been divided into a separate list per substituent/root/functionalTerm
	 * @throws ParsingException
	 */
	void writeWordXML(Element wordEl, ParseWord pw, List<String> tokens, List<List<Character>> annotations) throws ParsingException {
		int annotNumber = 0;
		int annotPos = 0;
		Element chunk = new Element(SUBSTITUENT_EL);
		wordEl.appendChild(chunk);
		Element lastTokenElement = null;
        for (String token : tokens) {
            if (annotPos >= annotations.get(annotNumber).size()) {
                annotPos = 0;
                annotNumber++;
                chunk = new Element(SUBSTITUENT_EL);
                wordEl.appendChild(chunk);
                lastTokenElement = null;
            }
            Element tokenElement = resourceManager.makeTokenElement(token,
                    annotations.get(annotNumber).get(annotPos));
            if (tokenElement != null) {//null for tokens that have ignoreWhenWritingXML set
                chunk.appendChild(tokenElement);
                lastTokenElement=tokenElement;
            }
            else if (lastTokenElement!=null && !token.equals("")){
            	if (lastTokenElement.getAttribute(SUBSEQUENTUNSEMANTICTOKEN_ATR)!=null){
            		lastTokenElement.getAttribute(SUBSEQUENTUNSEMANTICTOKEN_ATR).setValue(lastTokenElement.getAttributeValue(SUBSEQUENTUNSEMANTICTOKEN_ATR) + token);
            	}
            	else{
            		lastTokenElement.addAttribute(new Attribute(SUBSEQUENTUNSEMANTICTOKEN_ATR, token));
            	}
            }
            annotPos++;
        }
        WordType wordType = WordType.valueOf(wordEl.getAttributeValue(TYPE_ATR));
		if(wordType == WordType.full) {
			chunk.setLocalName(ROOT_EL);
		}
		else if(wordType == WordType.functionalTerm) {
			chunk.setLocalName(FUNCTIONALTERM_EL);
		}
	}
	
	/**
	 * Assigns an indication of stoichiometry to each child word rule of the moleculeEl.
	 * Throws an exception if there is a mismatch between the number of word rules and ratio.
	 * @param moleculeEl
	 * @param componentRatios
	 * @throws ParsingException
	 */
	private void applyStoichiometryIndicationToWordRules(Element moleculeEl,Integer[] componentRatios) throws ParsingException {
		List<Element> wordRules = XOMTools.getChildElementsWithTagName(moleculeEl, WORDRULE_EL);
		if (wordRules.size()!=componentRatios.length){
			throw new ParsingException("Component and stoichiometry indication indication mismatch. OPSIN believes there to be " +wordRules.size() +" components but " + componentRatios.length +" ratios were given!");
		}
		for (int i = 0; i < componentRatios.length; i++) {
			wordRules.get(i).addAttribute(new Attribute(STOICHIOMETRY_ATR,String.valueOf(componentRatios[i])));
		}
	}

}
