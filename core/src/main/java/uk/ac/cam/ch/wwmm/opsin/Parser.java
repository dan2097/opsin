package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

import uk.ac.cam.ch.wwmm.opsin.ParseWord.WordType;

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
		Parse parse = null;
		if (name.contains(", ")){
			try{
				parse = tokeniser.tokenize(tokeniser.uninvertCASName(name), false);
			}
			catch (ParsingException ignored) {
			}
		}
		else if (name.contains("; ")){//a mixture, spaces are sufficient for OPSIN to treat as a mixture. These spaces for obvious reasons must not be removed
			parse = tokeniser.tokenize(matchSemiColonSpace.matcher(name).replaceAll(" "), false);
		}
		boolean allowSpaceRemoval = parse ==null ? true : false;
		if (parse == null){
			parse = tokeniser.tokenize(name , true);
		}
		
		/* For cases where any of the parse's parseWords contain multiple annotations create a
		 * parse for each possibility. Hence after this process there may be multiple parse objects and
		 * the parseWords they contain will each only have one parseTokens object.
		 */
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
        for (String token : tokens) {
            if (annotPos >= annotations.get(annotNumber).size()) {
                annotPos = 0;
                annotNumber++;
                chunk = new Element(SUBSTITUENT_EL);
                wordEl.appendChild(chunk);
            }
            Element tokenElement = resourceManager.makeTokenElement(token,
                    annotations.get(annotNumber).get(annotPos));
            if (tokenElement != null) {//null for tokens that have ignoreWhenWritingXML set
                chunk.appendChild(tokenElement);
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

}
