package uk.ac.cam.ch.wwmm.opsin;

import static junit.framework.Assert.assertEquals;

import java.util.List;

import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

public class TokenizerTest {

	static Tokeniser tokenizer;
	static ReverseParseRules reverseParseRules;

	@BeforeClass
	public static void  setUp() throws Exception {
		ResourceGetter rg = new ResourceGetter("uk/ac/cam/ch/wwmm/opsin/resources/");
		ResourceManager rm = new ResourceManager(rg);
		tokenizer = new Tokeniser(new ParseRules(rm));
		reverseParseRules = new ReverseParseRules(rm);
	}
	
	@Test
	public void hexane() throws Exception {
		TokenizationResult result= tokenizer.tokenize("hexane", true);
		assertEquals(true, result.isSuccessfullyTokenized());
		assertEquals(true, result.isFullyInterpretable());
		assertEquals("", result.getUninterpretableName());
		assertEquals("", result.getUnparsableName());
		assertEquals("", result.getUnparsedName());
		Parse parse = result.getParse();
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
		Parse parse = tokenizer.tokenize("hexachlorohexane", true).getParse();
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
	public void ethylChloride() throws Exception {
		Parse parse = tokenizer.tokenize("ethyl chloride", true).getParse();
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


		parse = tokenizer.tokenize("ethylchloride", true).getParse();//missing space
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
	public void hexachlorohexaneeeeeee() throws Exception {
		TokenizationResult result = tokenizer.tokenize("hexachlorohexaneeeeeee", true);
		assertEquals("Unparsable", false, result.isSuccessfullyTokenized());
	}

	@Test
	public void bracketedHexachlorohexane() throws Exception {
		Parse parse = tokenizer.tokenize("(hexachloro)hexane", true).getParse();
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
		Parse parse = tokenizer.tokenize("methyl", true).getParse();
		assertEquals("One Word", 1, parse.getWords().size());
		ParseWord w = parse.getWords().get(0);
		assertEquals("One Parse", 1, w.getParseTokens().size());
		List<String> tokens = w.getParseTokens().get(0).getTokens();
		assertEquals("Three tokens", 3, tokens.size());
		assertEquals("First token: meth", "meth", tokens.get(0));
		assertEquals("Second token: yl", "yl", tokens.get(1));
		assertEquals("Third token: end of substituent", "", tokens.get(2));
	}
	
	@Test
	public void aceticacid() throws Exception {
		Parse parse = tokenizer.tokenize("acetic acid", true).getParse();
		assertEquals("One Word", 1, parse.getWords().size());
		ParseWord w = parse.getWords().get(0);
		assertEquals("One Parse", 1, w.getParseTokens().size());
		List<String> tokens = w.getParseTokens().get(0).getTokens();
		assertEquals("Three tokens", 3, tokens.size());
		assertEquals("First token: acet", "acet", tokens.get(0));
		assertEquals("Second token: ic acid", "ic acid", tokens.get(1));
		assertEquals("Third token: end of main group", "", tokens.get(2));
	}
	
	@Test
	public void CCCP() throws Exception {
		TokenizationResult result = tokenizer.tokenize("Carbonyl cyanide m-chlorophenyl oxime", true);
		assertEquals(true, result.isSuccessfullyTokenized());
		assertEquals(true, result.isFullyInterpretable());
		assertEquals("", result.getUninterpretableName());
		assertEquals("", result.getUnparsableName());
		assertEquals("", result.getUnparsedName());
		Parse parse = result.getParse();
		assertEquals("Four Words", 4, parse.getWords().size());
		ParseWord w = parse.getWords().get(0);
		assertEquals("One Parse", 1, w.getParseTokens().size());
		List<String> tokens = w.getParseTokens().get(0).getTokens();
		assertEquals("Three tokens", 3, tokens.size());
		assertEquals("First token: carbon", "carbon", tokens.get(0));
		assertEquals("Second token: yl", "yl", tokens.get(1));
		assertEquals("Third token: end of  substituent", "", tokens.get(2));
		
		w = parse.getWords().get(1);
		assertEquals("One Parse", 1, w.getParseTokens().size());
		tokens = w.getParseTokens().get(0).getTokens();
		assertEquals("Two tokens", 2, tokens.size());
		assertEquals("First token: cyanide", "cyanide", tokens.get(0));
		assertEquals("Second token: end of functionalTerm", "", tokens.get(1));
		
		w = parse.getWords().get(2);
		assertEquals("One Parse", 1, w.getParseTokens().size());
		tokens = w.getParseTokens().get(0).getTokens();
		assertEquals("Six tokens", 6, tokens.size());
		assertEquals("First token: m-", "m-", tokens.get(0));
		assertEquals("Second token: chloro", "chloro", tokens.get(1));
		assertEquals("Third token: end of  substituent", "", tokens.get(2));
		assertEquals("Fourth token: phen", "phen", tokens.get(3));
		assertEquals("Fifth token: yl", "yl", tokens.get(4));
		assertEquals("Sixth token: end of  substituent", "", tokens.get(5));
		
		w = parse.getWords().get(3);
		assertEquals("One Parse", 1, w.getParseTokens().size());
		tokens = w.getParseTokens().get(0).getTokens();
		assertEquals("Two tokens", 2, tokens.size());
		assertEquals("First token: oxime", "oxime", tokens.get(0));
		assertEquals("Second token: end of functionalTerm", "", tokens.get(1));
	}
	
	@Test
	public void CCCP_RL() throws Exception {
		TokenizationResult result = tokenizer.tokenizeRightToLeft(reverseParseRules, "Carbonyl cyanide m-chlorophenyl oxime", true);
		assertEquals(true, result.isSuccessfullyTokenized());
		assertEquals(true, result.isFullyInterpretable());
		assertEquals("", result.getUninterpretableName());
		assertEquals("", result.getUnparsableName());
		assertEquals("", result.getUnparsedName());
		Parse parse = result.getParse();
		assertEquals("Four Words", 4, parse.getWords().size());
		ParseWord w = parse.getWords().get(0);
		assertEquals("One Parse", 1, w.getParseTokens().size());
		List<String> tokens = w.getParseTokens().get(0).getTokens();
		assertEquals("Three tokens", 3, tokens.size());
		assertEquals("First token: carbon", "carbon", tokens.get(0));
		assertEquals("Second token: yl", "yl", tokens.get(1));
		assertEquals("Third token: end of  substituent", "", tokens.get(2));
		
		w = parse.getWords().get(1);
		assertEquals("One Parse", 1, w.getParseTokens().size());
		tokens = w.getParseTokens().get(0).getTokens();
		assertEquals("Two tokens", 2, tokens.size());
		assertEquals("First token: cyanide", "cyanide", tokens.get(0));
		assertEquals("Second token: end of functionalTerm", "", tokens.get(1));
		
		w = parse.getWords().get(2);
		assertEquals("One Parse", 1, w.getParseTokens().size());
		tokens = w.getParseTokens().get(0).getTokens();
		assertEquals("Six tokens", 6, tokens.size());
		assertEquals("First token: m-", "m-", tokens.get(0));
		assertEquals("Second token: chloro", "chloro", tokens.get(1));
		assertEquals("Third token: end of  substituent", "", tokens.get(2));
		assertEquals("Fourth token: phen", "phen", tokens.get(3));
		assertEquals("Fifth token: yl", "yl", tokens.get(4));
		assertEquals("Sixth token: end of  substituent", "", tokens.get(5));
		
		w = parse.getWords().get(3);
		assertEquals("One Parse", 1, w.getParseTokens().size());
		tokens = w.getParseTokens().get(0).getTokens();
		assertEquals("Two tokens", 2, tokens.size());
		assertEquals("First token: oxime", "oxime", tokens.get(0));
		assertEquals("Second token: end of functionalTerm", "", tokens.get(1));
	}
	
	@Test
	public void partiallyInterpretatableLR() throws Exception {
		TokenizationResult result =tokenizer.tokenize("ethyl-2H-foo|ene", true);
		assertEquals(false, result.isSuccessfullyTokenized());
		assertEquals(false, result.isFullyInterpretable());
		assertEquals("2H-foo|ene", result.getUninterpretableName());
		assertEquals("foo|ene", result.getUnparsableName());
		assertEquals("ethyl-2H-foo|ene", result.getUnparsedName());
	}
	
	@Test
	public void partiallyInterpretatableRL1() throws Exception {
		TokenizationResult result =tokenizer.tokenizeRightToLeft(reverseParseRules, "ethyl-2H-foo|ene", true);
		assertEquals(false, result.isSuccessfullyTokenized());
		assertEquals(false, result.isFullyInterpretable());
		assertEquals("ethyl-2H-foo|ene", result.getUninterpretableName());
		assertEquals("ethyl-2H-foo|", result.getUnparsableName());
		assertEquals("ethyl-2H-foo|ene", result.getUnparsedName());
	}
	
	@Test
	public void partiallyInterpretatableRL2() throws Exception {
		TokenizationResult result =tokenizer.tokenizeRightToLeft(reverseParseRules, "fooylpyridine oxide", true);
		assertEquals(false, result.isSuccessfullyTokenized());
		assertEquals(false, result.isFullyInterpretable());
		assertEquals("fooyl", result.getUninterpretableName());
		assertEquals("f", result.getUnparsableName());//o as in the end of thio then oyl
		assertEquals("fooylpyridine", result.getUnparsedName());
	}

	@Test
	public void tokenizeDoesNotTokenizeUnTokenizableName() throws Exception {
		TokenizationResult result =tokenizer.tokenize("ethyl phen|foo toluene", true);
		assertEquals(false, result.isSuccessfullyTokenized());
	}

	@Test
	public void tokenizePreservesSpacesInUninterpretableName() throws Exception {
		TokenizationResult result =tokenizer.tokenize("ethyl phen|foo toluene", true);
		assertEquals("phen|foo toluene", result.getUninterpretableName());
	}

	@Test
	public void tokenizePreservesSpacesInUnparsableName() throws Exception {
		TokenizationResult result =tokenizer.tokenize("ethyl phen|foo toluene", true);
		assertEquals("|foo toluene", result.getUnparsableName());
	}

	@Test
	public void tokenizePreservesSpacesInUnparsedName() throws Exception {
		TokenizationResult result =tokenizer.tokenize("ethyl phen|foo toluene", true);
		assertEquals("phen|foo toluene", result.getUnparsedName());
	}
	
	@Test
	public void tokenizePreservesSpacesInUninterpretableName2() throws Exception {
		TokenizationResult result =tokenizer.tokenize("eth yl |foo arene", true);
		assertEquals("|foo arene", result.getUninterpretableName());
	}

	@Test
	public void tokenizePreservesSpacesInUnparsableName2() throws Exception {
		TokenizationResult result =tokenizer.tokenize("eth yl |foo arene", true);
		assertEquals("|foo arene", result.getUnparsableName());
	}

	@Test
	public void tokenizePreservesSpacesInUnparsedName2() throws Exception {
		TokenizationResult result =tokenizer.tokenize("eth yl |foo arene", true);
		assertEquals("|foo arene", result.getUnparsedName());
	}
	
	@Test
	@Ignore
	public void partiallyInterpretatableRLWithSpaces() throws Exception {
		TokenizationResult result =tokenizer.tokenizeRightToLeft(reverseParseRules, "ethyl foo|yl toluene", true);
		assertEquals(false, result.isSuccessfullyTokenized());
		assertEquals("ethyl foo|yl", result.getUninterpretableName());
		assertEquals("ethyl foo|", result.getUnparsableName());
		assertEquals("ethyl foo|yl", result.getUnparsedName());
	}
}
