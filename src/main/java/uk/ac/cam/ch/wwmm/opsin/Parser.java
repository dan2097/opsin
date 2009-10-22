package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.List;

import nu.xom.Attribute;
import nu.xom.Element;
import uk.ac.cam.ch.wwmm.opsin.ParseRules.TwoReturnValues;
import uk.ac.cam.ch.wwmm.opsin.PreProcessor.OpsinMode;
import uk.ac.cam.ch.wwmm.ptclib.misc.Combinations;

/**Conducts finite-state parsing on chemical names. Preserves the name itself, just adding XML
 * annotation to the various constituents on the name.
 *
 * @author ptc24/d387
 *
 */
class Parser {

	/**The rules by which names are divided into words.*/
	private WordRules wordRules;
	/**Performs finite-state allocation of roles ("annotations") to tokens.*/
	private ParseRules parseRules;
	/**Holds the various tokens used.*/
	private TokenManager tokenManager;

	/**Initialises the parser.
	 * @param tokenManager 
	 * @param parseRules 
	 * @param wordRules
	 */
	Parser(WordRules wordRules, ParseRules parseRules, TokenManager tokenManager) {
		this.wordRules =wordRules;
		this.parseRules = parseRules;
		this.tokenManager =tokenManager;
	}

	/**Parses a chemical name to an XML representation of the parse.
	 *
	 * @param name The name to parse.
	 * @param mode 
	 * @return The parse.
	 * @throws ParsingException If the name is unparsable.
	 */
	List<Element> parse(String name, OpsinMode mode) throws ParsingException {
		Parse p = new Parse();
		p.name = name;
		p.words = new ArrayList<ParseWord>();
		wordRules.parse(p, mode);

		List<Integer> parseCounts = new ArrayList<Integer>();
		int totalParses = 1;

		for (int j = 0; j < p.words.size(); j++) {
			ParseWord pw = p.words.get(j);
			if(pw.wordType.equals("literal")) {
				parseCounts.add(1);
			} else {
				TwoReturnValues returned = parseRules.getParses(pw.word, pw.wordType);
				pw.parseTokens =returned.getFirst();
				if (pw.parseTokens.size()>128){
					throw new ParsingException("Too many parses generated, the current limit is 128: " + pw.parseTokens.size());
				}
				boolean interSubHyphenAllowed =returned.getSecond();
				if (pw.parseTokens.size()==0 && p.wordRule.equals("binaryOrOther")){
					if (j +1 <p.words.size()){
						ParseWord nextWord =p.words.get(j+1);
						if(interSubHyphenAllowed ==true && !nextWord.word.startsWith("-")){
							pw.word+="-";//only add a - if the word is lexable up to this point and allows a hyphen as it's next character
						}
						pw.word+=nextWord.word;
						p.words.remove(j+1);
						j--;
						continue;
					}
				}

				if(pw.parseTokens.size() == 0) throw new ParsingException("No parses for " + name + " using wordrule: " + p.wordRule);
				parseCounts.add(pw.parseTokens.size());
				totalParses *= pw.parseTokens.size();
			}
		}

		//TODO combinations the Combinations class with this class
		List<List<Integer>> combinations = Combinations.makeCombinations(parseCounts);

		List<Parse> parses = new ArrayList<Parse>();

		for(List<Integer> c : combinations) {
			Parse pp = p.deepCopy();
			for(int i=0;i<c.size();i++) {
				if(parseCounts.get(i) > 1) {
					ParseWord pw = pp.words.get(i);
					List<ParseTokens> ptl = new ArrayList<ParseTokens>();
					ptl.add(pw.parseTokens.get(c.get(i)));
					pw.parseTokens = ptl;
				}
			}
			parses.add(pp);
		}

		List<Element> results = new ArrayList<Element>();
		for(Parse pp : parses) {
			Element elem = new Element("molecule");
			elem.addAttribute(new Attribute("name", name));
			elem.addAttribute(new Attribute("wordRule", pp.wordRule));
			for(ParseWord pw : pp.words) {
				Element word = new Element("word");
				elem.appendChild(word);
				word.addAttribute(new Attribute("type", pw.wordType));
				if(pw.wordType.equals("literal")) {
					word.appendChild(pw.word);
				} else {
					for(ParseTokens pt : pw.parseTokens) {
						writeWordXML(word, pw, pt.tokens, parseRules.chunkAnnotations(pt.annotations));
					}
				}
			}
			results.add(elem);
		}
		return results;
	}

	/**Write the XML corresponding to a particular word in a parse.
	 * Assumes that there is at least one possible annotation, and that only the
	 * first is required.
	 *
	 * @param elem The empty XML word element to be written into.
	 * @param pw The ParseWord for the word.
	 * @param tokens The list of tokens.
	 * @param annotations The list of possible annotations.
	 * @throws ParsingException
	 */
	void writeWordXML(Element elem, ParseWord pw, List<String> tokens, List<List<Character>> annotations) throws ParsingException {
		int annotNumber = 0;
		int annotPos = 0;
		Element chunk = new Element("substituent");
		elem.appendChild(chunk);
		for(int i=0;i<tokens.size();i++) {
			if(annotPos >= annotations.get(annotNumber).size()) {
				annotPos = 0;
				annotNumber++;
				chunk = new Element("substituent");
				elem.appendChild(chunk);
			}
			Element tokenElement = tokenManager.makeTokenElement(tokens.get(i),
					annotations.get(annotNumber).get(annotPos));
			if (tokenElement !=null){//null for tokens that have ignoreWhenWritingXML set
				chunk.appendChild(tokenElement);
			}
			annotPos++;
		}
		if(pw.wordType.equals("full")) {
			chunk.setLocalName("root");
		}
	}

}
