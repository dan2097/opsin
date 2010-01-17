package uk.ac.cam.ch.wwmm.opsin;

import static junit.framework.Assert.assertEquals;

import java.util.List;

import org.junit.Before;
import org.junit.Test;

public class TokenizerTest {

	Tokeniser tokenizer;

	@Before
	public void setUp() throws Exception {
		ResourceGetter rg = new ResourceGetter("uk/ac/cam/ch/wwmm/opsin/resources/");
		tokenizer = new Tokeniser(new ParseRules(new TokenManager(rg),rg));
	}
	
	@Test
	public void ethylChloride() throws Exception {
		Parse parse = tokenizer.tokenize("ethyl chloride");
		assertEquals("Two Words", 2, parse.getWords().size());
		ParseWord w = parse.getWord(0);
		assertEquals("One Parse", 1, w.getParseTokens().size());
		List<String> tokens = w.getParseTokens().get(0).getTokens();
		assertEquals("Three tokens", 3, tokens.size());
		assertEquals("First token: eth", "eth", tokens.get(0));
		assertEquals("Second token: yl", "yl", tokens.get(1));
		assertEquals("Third token: end of substituent", "", tokens.get(2));
		w = parse.getWord(1);
		assertEquals("One Parse", 1, w.getParseTokens().size());
		tokens = w.getParseTokens().get(0).getTokens();
		assertEquals("Two tokens", 2, tokens.size());
		assertEquals("First token: chloride", "chloride", tokens.get(0));
		assertEquals("Second token: end of functionalTerm", "", tokens.get(1));


		parse = tokenizer.tokenize("ethylchloride");//missing space
		assertEquals("Two Words", 2, parse.getWords().size());
		w = parse.getWord(0);
		assertEquals("One Parse", 1, w.getParseTokens().size());
		tokens = w.getParseTokens().get(0).getTokens();
		assertEquals("Three tokens", 3, tokens.size());
		assertEquals("First token: eth", "eth", tokens.get(0));
		assertEquals("Second token: yl", "yl", tokens.get(1));
		assertEquals("Third token: end of substituent", "", tokens.get(2));
		w = parse.getWord(1);
		assertEquals("One Parse", 1, w.getParseTokens().size());
		tokens = w.getParseTokens().get(0).getTokens();
		assertEquals("Two tokens", 2, tokens.size());
		assertEquals("First token: chloride", "chloride", tokens.get(0));
		assertEquals("Second token: end of functionalTerm", "", tokens.get(1));
	}
	
	@Test
	public void hexane() throws Exception {
		Parse parse = tokenizer.tokenize("hexane");
		assertEquals("One Word", 1, parse.getWords().size());
		ParseWord w = parse.getWords().get(0);
		assertEquals("One Parse", 1, w.getParseTokens().size());
		List<String> tokens = w.getParseTokens().get(0).getTokens();
		assertEquals("Three tokens", 3, tokens.size());
		assertEquals("First token: hex", "hex", tokens.get(0));
		assertEquals("Second token: ane", "ane", tokens.get(1));
		assertEquals("Third token: end of main group", "", tokens.get(2));
	}
	
	@Test
	public void hexachlorohexane() throws Exception {
		Parse parse = tokenizer.tokenize("hexachlorohexane");
		assertEquals("One Word", 1, parse.getWords().size());
		ParseWord w = parse.getWords().get(0);
		assertEquals("One Parse", 1, w.getParseTokens().size());
		List<String> tokens = w.getParseTokens().get(0).getTokens();
		assertEquals("Six tokens", 6, tokens.size());
		assertEquals("First token: hexa", "hexa", tokens.get(0));
		assertEquals("Second token: chloro", "chloro", tokens.get(1));
		assertEquals("Third token: end of main substituent", "", tokens.get(2));
		assertEquals("Fourth token: hex", "hex", tokens.get(3));
		assertEquals("Fifth token: ane", "ane", tokens.get(4));
		assertEquals("Sixth token: end of main group", "", tokens.get(5));
	}
	@Test
	public void hexachlorohexaneeeeeee() throws Exception {
		boolean failed;
		try{
			tokenizer.tokenize("hexachlorohexaneeeeeee");
			failed = false;
		}
		catch (ParsingException e){
			failed = true;
		}
		assertEquals("Unparsable", true, failed);
	}

	@Test
	public void bracketedHexachlorohexane() throws Exception {
		Parse parse = tokenizer.tokenize("(hexachloro)hexane");
		assertEquals("One Word", 1, parse.getWords().size());
		ParseWord w = parse.getWords().get(0);
		assertEquals("One Parse", 1, w.getParseTokens().size());
		List<String> tokens = w.getParseTokens().get(0).getTokens();
		assertEquals("Eight tokens", 8, tokens.size());
		assertEquals("First token: {", "(", tokens.get(0));
		assertEquals("Second token: hexa", "hexa", tokens.get(1));
		assertEquals("Third token: chloro", "chloro", tokens.get(2));
		assertEquals("Fourth token: )", ")", tokens.get(3));
		assertEquals("Fifth token: end of main substituent", "", tokens.get(4));
		assertEquals("Sixth token: hex", "hex", tokens.get(5));
		assertEquals("Seventh token: ane", "ane", tokens.get(6));
		assertEquals("Eigth token: end of main group", "", tokens.get(7));
	}
	
	@Test
	public void methyl() throws Exception {
		Parse parse = tokenizer.tokenize("methyl");
		assertEquals("One Word", 1, parse.getWords().size());
		ParseWord w = parse.getWords().get(0);
		assertEquals("One Parse", 1, w.getParseTokens().size());
		List<String> tokens = w.getParseTokens().get(0).getTokens();
		assertEquals("Three tokens", 3, tokens.size());
		assertEquals("First token: meth", "meth", tokens.get(0));
		assertEquals("Second token: yl", "yl", tokens.get(1));
		assertEquals("Third token: end of substituent", "", tokens.get(2));
	}
	
}
