package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import uk.ac.cam.ch.wwmm.opsin.ParseWord.WordType;
import uk.ac.cam.ch.wwmm.opsin.Tokeniser.TokenizationResult;

import nu.xom.Attribute;
import nu.xom.Element;

import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

/** For a list of integers, generate all lists of non-negative integers where all
 * of the list members are strictly lower than their corresponding integer in the
 * input list.
 *
 * @author ptc24
 *
 */
final class Combinations {

	private List<List<Integer>> combList;
	private List<Integer> scheme;

	private Combinations() {

	}

	/**For a list of integers, generate all lists of non-negative integers where
	 * all of the list members are strictly lower than their corresponding
	 * integer in the input list.
	 *
	 * @param scheme The input list of integers.
	 * @return The list of lists.
	 */
	public static List<List<Integer>> makeCombinations(List<Integer> scheme) {
		Combinations c = new Combinations();
		c.combList = new ArrayList<List<Integer>>();
		c.scheme = scheme;
		c.grow(new ArrayList<Integer>());
		return c.combList;
	}

	private void grow(List<Integer> soFar) {
		if(soFar.size() == scheme.size()) {
			combList.add(soFar);
		} else {
			for(int i=0;i<scheme.get(soFar.size());i++) {
				List<Integer> ll = new ArrayList<Integer>();
				ll.addAll(soFar);
				ll.add(i);
				grow(ll);
			}
		}
	}
}

/**Conducts finite-state parsing on chemical names. Preserves the name itself, just adding XML
 * annotation to the various constituents on the name.
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
	
	private final static Pattern matchSemiColonSpace = Pattern.compile("; ");
	private final static Pattern matchStoichiometryIndication = Pattern.compile("[ ]?[\\{\\[\\(](\\d+|\\?)([:/](\\d+|\\?))+[\\}\\]\\)]$");
	private final static Pattern matchColon = Pattern.compile(":");
	private final static Pattern matchForwardSlash = Pattern.compile("/");

	/**
	 * No-argument constructor. Uses ResouceGetter found at
	 * uk/ac/cam/ch/wwmm/opsin/resources/
	 */
	Parser() throws Exception {
		ResourceGetter resources = new ResourceGetter("uk/ac/cam/ch/wwmm/opsin/resources/");
		this.wordRules = new WordRules(resources);
		this.resourceManager = new ResourceManager(resources);
		ParseRules parseRules = new ParseRules(this.resourceManager);
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
				componentRatios = processStoichometryIndication(m.group());
				name = m.replaceAll("");
			}
		}
		Parse parse = null;
		if (name.contains(", ")){
			try{
				TokenizationResult tokenizationResult = tokeniser.tokenize(tokeniser.uninvertCASName(name), false);
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
		if (parses.size()>128){
			throw new ParsingException("Too many parses generated, the current limit is 128: " + parses.size());
		}

		
		List<Element> results = new ArrayList<Element>();
		for(Parse pp : parses) {
			Element moleculeEl = new Element(MOLECULE_EL);
			moleculeEl.addAttribute(new Attribute(NAME_ATR, name));
			for(ParseWord pw : pp.getWords()) {
				Element word = new Element(WORD_EL);
				moleculeEl.appendChild(word);
				if (pw.getParseTokens().size() >1){
					throw new ParsingException("OPSIN bug: parseWord had multiple annotations after creating addition parses step");
				}
				
				pw.setWordType(OpsinTools.determineWordType(pw.getParseTokens().get(0).getAnnotations()));
				word.addAttribute(new Attribute(TYPE_ATR, pw.getWordType().toString()));
				if (pw.getWord().startsWith("-")){//we want -acid to be the same as acid
					word.addAttribute(new Attribute(VALUE_ATR, pw.getWord().substring(1)));
				}
				else{
					word.addAttribute(new Attribute(VALUE_ATR, pw.getWord()));
				}
				for(ParseTokens pt : pw.getParseTokens()) {
					writeWordXML(word, pw, pt.getTokens(), tokeniser.chunkAnnotations(pt.getAnnotations()));
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
					applyStoichometryIndicationToWordRules(moleculeEl, componentRatios);
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

	static Integer[] processStoichometryIndication(String ratioString) throws ParsingException {
		ratioString = ratioString.trim();
		ratioString = ratioString.substring(1, ratioString.length()-1);
		String[] ratioStrings = matchColon.split(ratioString);
		if (ratioStrings.length ==1){
			ratioStrings = matchForwardSlash.split(ratioString);
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
		} catch (Exception e) {
			throw new RuntimeException(e);
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
	 * parse for each possibility. Hence after this process there may be multiple parse objects and
	 * the parseWords they contain will each only have one parseTokens object.
	 * @param parse
	 * @return
	 */
	private List<Parse> generateParseCombinations(Parse parse) {
		List<Integer> parseCounts = new ArrayList<Integer>();
		for (ParseWord pw : parse.getWords()) {
			parseCounts.add(pw.getParseTokens().size());
		}
		List<List<Integer>> combinations = Combinations.makeCombinations(parseCounts);
		List<Parse> parses = new ArrayList<Parse>();
		for(List<Integer> c : combinations) {
			Parse parseCopy = parse.deepCopy();
			for(int i=0; i<c.size(); i++) {
				if(parseCounts.get(i) > 1) {
					ParseWord pw = parseCopy.getWord(i);
					List<ParseTokens> ptl = new ArrayList<ParseTokens>();
					ptl.add(pw.getParseTokens().get(c.get(i)));
					pw.setParseTokens(ptl);
				}
			}
			parses.add(parseCopy);
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
            	if (lastTokenElement.getAttribute(SUBSEQUENTUNSEMANTICTOKEN_EL)!=null){
            		lastTokenElement.getAttribute(SUBSEQUENTUNSEMANTICTOKEN_EL).setValue(lastTokenElement.getAttributeValue(SUBSEQUENTUNSEMANTICTOKEN_EL) + token);
            	}
            	else{
            		lastTokenElement.addAttribute(new Attribute(SUBSEQUENTUNSEMANTICTOKEN_EL, token));
            	}
            }
            annotPos++;
        }
		if(pw.getWordType() == WordType.full) {
			chunk.setLocalName(ROOT_EL);
		}
		else if(pw.getWordType() == WordType.functionalTerm) {
			chunk.setLocalName(FUNCTIONALTERM_EL);
		}
	}
	
	/**
	 * Assigns an indication of stoichometry to each child word rule of the moleculeEl.
	 * Throws an exception if there is a mismatch between the number of word rules and ratio.
	 * @param moleculeEl
	 * @param componentRatios
	 * @throws ParsingException
	 */
	private void applyStoichometryIndicationToWordRules(Element moleculeEl,Integer[] componentRatios) throws ParsingException {
		List<Element> wordRules = XOMTools.getChildElementsWithTagName(moleculeEl, WORDRULE_EL);
		if (wordRules.size()!=componentRatios.length){
			throw new ParsingException("Component and stoichometry indication indication mismatch. OPSIN believes there to be " +wordRules.size() +" components but " + componentRatios.length +" ratios were given!");
		}
		for (int i = 0; i < componentRatios.length; i++) {
			wordRules.get(i).addAttribute(new Attribute(STOICHOMETRY_ATR,String.valueOf(componentRatios[i])));
		}
	}

}
