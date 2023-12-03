package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.util.List;

import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

public class TokenizerTest {

	private static Tokeniser tokenizer;
	private static ReverseParseRules reverseParseRules;

	@BeforeAll
	public static void setUp() throws IOException{
		ResourceGetter rg = new ResourceGetter("uk/ac/cam/ch/wwmm/opsin/resources/");
		ResourceManager rm = new ResourceManager(rg);
		tokenizer = new Tokeniser(new ParseRules(rm));
		reverseParseRules = new ReverseParseRules(rm);
	}
	
	@AfterAll
	public static void cleanUp(){
		tokenizer = null;
		reverseParseRules = null;
	}
	
	@Test
	public void hexane() throws ParsingException{
		TokenizationResult result= tokenizer.tokenize("hexane", true);
		assertTrue(result.isSuccessfullyTokenized());
		assertTrue(result.isFullyInterpretable());
		assertEquals("", result.getUninterpretableName());
		assertEquals("", result.getUnparsableName());
		assertEquals("", result.getUnparsedName());
		Parse parse = result.getParse();
		assertEquals(1, parse.getWords().size(), "One Word");
		ParseWord w = parse.getWords().get(0);
		assertEquals(1, w.getParseTokens().size(), "One Parse");
		List<String> tokens = w.getParseTokens().get(0).getTokens();
		assertEquals(3, tokens.size(), "Three tokens");
		assertEquals("hex", tokens.get(0), "First token: hex");
		assertEquals("ane", tokens.get(1), "Second token: ane");
		assertEquals("", tokens.get(2), "Third token: end of main group");
	}
	
	@Test
	public void hexachlorohexane() throws ParsingException{
		Parse parse = tokenizer.tokenize("hexachlorohexane", true).getParse();
		assertEquals(1, parse.getWords().size(), "One Word");
		ParseWord w = parse.getWords().get(0);
		assertEquals(1, w.getParseTokens().size(), "One Parse");
		List<String> tokens = w.getParseTokens().get(0).getTokens();
		assertEquals(7, tokens.size(), "Seven tokens");
		assertEquals("hex", tokens.get(0), "First token: hex");
		assertEquals("a", tokens.get(1), "Second token: a");
		assertEquals("chloro", tokens.get(2), "Third token: chloro");
		assertEquals("", tokens.get(3), "Fourth token: end of main substituent");
		assertEquals("hex", tokens.get(4), "Fifth token: hex");
		assertEquals("ane", tokens.get(5), "Sixth token: ane");
		assertEquals("", tokens.get(6), "Seventh token: end of main group");
	}
	
	@Test
	public void ethylChloride() throws ParsingException {
		Parse parse = tokenizer.tokenize("ethyl chloride", true).getParse();
		assertEquals(2, parse.getWords().size(), "Two Words");
		ParseWord w = parse.getWord(0);
		assertEquals(1, w.getParseTokens().size(), "One Parse");
		List<String> tokens = w.getParseTokens().get(0).getTokens();
		assertEquals(3, tokens.size(), "Three tokens");
		assertEquals("eth", tokens.get(0), "First token: eth");
		assertEquals("yl", tokens.get(1), "Second token: yl");
		assertEquals("", tokens.get(2), "Third token: end of substituent");
		w = parse.getWord(1);
		assertEquals(1, w.getParseTokens().size(), "One Parse");
		tokens = w.getParseTokens().get(0).getTokens();
		assertEquals(2, tokens.size(), "Two tokens");
		assertEquals("chloride", tokens.get(0), "First token: chloride");
		assertEquals("", tokens.get(1), "Second token: end of functionalTerm");


		parse = tokenizer.tokenize("ethylchloride", true).getParse();//missing space
		assertEquals(2, parse.getWords().size(), "Two Words");
		w = parse.getWord(0);
		assertEquals(1, w.getParseTokens().size(), "One Parse");
		tokens = w.getParseTokens().get(0).getTokens();
		assertEquals(3, tokens.size(), "Three tokens");
		assertEquals("eth", tokens.get(0), "First token: eth");
		assertEquals("yl", tokens.get(1), "Second token: yl");
		assertEquals("", tokens.get(2), "Third token: end of substituent");
		w = parse.getWord(1);
		assertEquals(1, w.getParseTokens().size(), "One Parse");
		tokens = w.getParseTokens().get(0).getTokens();
		assertEquals(2, tokens.size(), "Two tokens");
		assertEquals("chloride", tokens.get(0), "First token: chloride");
		assertEquals("", tokens.get(1), "Second token: end of functionalTerm");
	}
	
	@Test
	public void hexachlorohexaneeeeeee() throws ParsingException{
		TokenizationResult result = tokenizer.tokenize("hexachlorohexaneeeeeee", true);
		assertFalse(result.isSuccessfullyTokenized(), "Unparsable");
	}

	@Test
	public void bracketedHexachlorohexane() throws ParsingException{
		Parse parse = tokenizer.tokenize("(hexachloro)hexane", true).getParse();
		assertEquals(1, parse.getWords().size(), "One Word");
		ParseWord w = parse.getWords().get(0);
		assertEquals(1, w.getParseTokens().size(), "One Parse");
		List<String> tokens = w.getParseTokens().get(0).getTokens();
		assertEquals(9, tokens.size(),"Nine tokens");
		assertEquals("(", tokens.get(0), "First token: (");
		assertEquals("hex", tokens.get(1), "Second token: hex");
		assertEquals("a", tokens.get(2), "Third token: a");
		assertEquals("chloro", tokens.get(3), "Fourth token: chloro");
		assertEquals(")", tokens.get(4), "Fifth token: )");
		assertEquals("", tokens.get(5), "Sixth token: end of main substituent");
		assertEquals("hex", tokens.get(6), "Seventh token: hex");
		assertEquals("ane", tokens.get(7), "Eigth token: ane");
		assertEquals("", tokens.get(8), "Ninth token: end of main group");
	}
	
	@Test
	public void methyl() throws ParsingException{
		Parse parse = tokenizer.tokenize("methyl", true).getParse();
		assertEquals(1, parse.getWords().size(), "One Word");
		ParseWord w = parse.getWords().get(0);
		assertEquals(1, w.getParseTokens().size(), "One Parse");
		List<String> tokens = w.getParseTokens().get(0).getTokens();
		assertEquals(3, tokens.size(), "Three tokens");
		assertEquals("meth", tokens.get(0), "First token: meth");
		assertEquals("yl", tokens.get(1), "Second token: yl");
		assertEquals("", tokens.get(2), "Third token: end of substituent");
	}
	
	@Test
	public void aceticacid() throws ParsingException{
		Parse parse = tokenizer.tokenize("acetic acid", true).getParse();
		assertEquals(1, parse.getWords().size(), "One Word");
		ParseWord w = parse.getWords().get(0);
		assertEquals(1, w.getParseTokens().size(), "One Parse");
		List<String> tokens = w.getParseTokens().get(0).getTokens();
		assertEquals(3, tokens.size(), "Three tokens");
		assertEquals("acet", tokens.get(0), "First token: acet");
		assertEquals("ic acid", tokens.get(1), "Second token: ic acid");
		assertEquals("", tokens.get(2), "Third token: end of main group");
	}

	@Test
	public void acceptableInterWordBreaks() throws ParsingException{
		assertTrue(tokenizer.tokenize("methane ethane", false).isSuccessfullyTokenized());
		assertTrue(tokenizer.tokenize("methane-ethane", false).isSuccessfullyTokenized());
		assertTrue(tokenizer.tokenize("methane - ethane", false).isSuccessfullyTokenized());
		assertFalse(tokenizer.tokenize("methane -ethane", false).isSuccessfullyTokenized());
		assertFalse(tokenizer.tokenize("methane - ", false).isSuccessfullyTokenized());
		
		assertTrue(tokenizer.tokenizeRightToLeft(reverseParseRules, "methane ethane", false).isSuccessfullyTokenized());
		assertTrue(tokenizer.tokenizeRightToLeft(reverseParseRules, "methane-ethane", false).isSuccessfullyTokenized());
		assertTrue(tokenizer.tokenizeRightToLeft(reverseParseRules, "methane - ethane", false).isSuccessfullyTokenized());
		assertFalse(tokenizer.tokenizeRightToLeft(reverseParseRules, "methane -ethane", false).isSuccessfullyTokenized());
		assertFalse(tokenizer.tokenizeRightToLeft(reverseParseRules, "methane - ", false).isSuccessfullyTokenized());
	}
	
	@Test
	public void compoundWithValidUse() throws ParsingException{
		TokenizationResult result =tokenizer.tokenize("benzene compound with toluene", true);
		assertTrue(result.isSuccessfullyTokenized());
		TokenizationResult result2 =tokenizer.tokenize("benzene and toluene", true);
		assertTrue(result2.isSuccessfullyTokenized());
	}

	@Test
	public void compoundWithInvalidUse1() throws ParsingException{
		TokenizationResult result =tokenizer.tokenize("ethyl and toluene", true);
		assertFalse(result.isSuccessfullyTokenized());
	}

	@Test
	public void compoundWithInvalidUse2() throws ParsingException{
		TokenizationResult result =tokenizer.tokenize("and benzene", true);
		assertFalse(result.isSuccessfullyTokenized());
	}
	
	@Test
	public void CCCP() throws ParsingException{
		TokenizationResult result = tokenizer.tokenize("Carbonyl cyanide m-chlorophenyl oxime", true);
		assertTrue(result.isSuccessfullyTokenized());
		assertTrue(result.isFullyInterpretable());
		assertEquals("", result.getUninterpretableName());
		assertEquals("", result.getUnparsableName());
		assertEquals("", result.getUnparsedName());
		Parse parse = result.getParse();
		assertEquals(4, parse.getWords().size(), "Four Words");
		ParseWord w = parse.getWords().get(0);
		assertEquals(1, w.getParseTokens().size(), "One Parse");
		List<String> tokens = w.getParseTokens().get(0).getTokens();
		assertEquals(3, tokens.size(), "Three tokens");
		assertEquals("carbon", tokens.get(0), "First token: carbon");
		assertEquals("yl", tokens.get(1), "Second token: yl");
		assertEquals("", tokens.get(2), "Third token: end of  substituent");
		
		w = parse.getWords().get(1);
		assertEquals(1, w.getParseTokens().size(), "One Parse");
		tokens = w.getParseTokens().get(0).getTokens();
		assertEquals(2, tokens.size(), "Two tokens");
		assertEquals("cyanide", tokens.get(0), "First token: cyanide");
		assertEquals("", tokens.get(1), "Second token: end of functionalTerm");
		
		w = parse.getWords().get(2);
		assertEquals(1, w.getParseTokens().size(), "One Parse");
		tokens = w.getParseTokens().get(0).getTokens();
		assertEquals(5, tokens.size(), "Five tokens");
		assertEquals("m-", tokens.get(0), "First token: m-");
		assertEquals("chloro", tokens.get(1), "Second token: chloro");
		assertEquals("", tokens.get(2), "Third token: end of  substituent");
		assertEquals("phenyl", tokens.get(3), "Fourth token: phenyl");
		assertEquals("", tokens.get(4), "Fifth token: end of  substituent");
		
		w = parse.getWords().get(3);
		assertEquals(1, w.getParseTokens().size(), "One Parse");
		tokens = w.getParseTokens().get(0).getTokens();
		assertEquals(2, tokens.size(), "Two tokens");
		assertEquals("oxime", tokens.get(0), "First token: oxime");
		assertEquals("", tokens.get(1), "Second token: end of functionalTerm");
	}
	
	@Test
	public void CCCP_RL() throws ParsingException{
		TokenizationResult result = tokenizer.tokenizeRightToLeft(reverseParseRules, "Carbonyl cyanide m-chlorophenyl oxime", true);
		assertTrue(result.isSuccessfullyTokenized());
		assertTrue(result.isFullyInterpretable());
		assertEquals("", result.getUninterpretableName());
		assertEquals("", result.getUnparsableName());
		assertEquals("", result.getUnparsedName());
		Parse parse = result.getParse();
		assertEquals(4, parse.getWords().size(), "Four Words");
		ParseWord w = parse.getWords().get(0);
		assertEquals(1, w.getParseTokens().size(), "One Parse");
		List<String> tokens = w.getParseTokens().get(0).getTokens();
		assertEquals(3, tokens.size(), "Three tokens");
		assertEquals("carbon", tokens.get(0), "First token: carbon");
		assertEquals("yl", tokens.get(1), "Second token: yl");
		assertEquals("", tokens.get(2), "Third token: end of  substituent");
		
		w = parse.getWords().get(1);
		assertEquals(1, w.getParseTokens().size(), "One Parse");
		tokens = w.getParseTokens().get(0).getTokens();
		assertEquals(2, tokens.size(), "Two tokens");
		assertEquals("cyanide", tokens.get(0), "First token: cyanide");
		assertEquals("", tokens.get(1), "Second token: end of functionalTerm");
		
		w = parse.getWords().get(2);
		assertEquals(1, w.getParseTokens().size(), "One Parse");
		tokens = w.getParseTokens().get(0).getTokens();
		assertEquals(5, tokens.size(), "Five tokens");
		assertEquals("m-", tokens.get(0), "First token: m-");
		assertEquals("chloro", tokens.get(1), "Second token: chloro");
		assertEquals("", tokens.get(2), "Third token: end of  substituent");
		assertEquals("phenyl", tokens.get(3), "Fourth token: phenyl");
		assertEquals("", tokens.get(4), "Fifth token: end of  substituent");
		
		w = parse.getWords().get(3);
		assertEquals(1, w.getParseTokens().size(), "One Parse");
		tokens = w.getParseTokens().get(0).getTokens();
		assertEquals(2, tokens.size(), "Two tokens");
		assertEquals("oxime", tokens.get(0), "First token: oxime");
		assertEquals("", tokens.get(1), "Second token: end of functionalTerm");
	}
	
	@Test
	public void partiallyInterpretatableLR() throws ParsingException{
		TokenizationResult result =tokenizer.tokenize("ethyl-2H-foo|ene", true);
		assertFalse(result.isSuccessfullyTokenized());
		assertFalse(result.isFullyInterpretable());
		assertEquals("2H-foo|ene", result.getUninterpretableName());
		assertEquals("foo|ene", result.getUnparsableName());
		assertEquals("ethyl-2H-foo|ene", result.getUnparsedName());
	}
	
	@Test
	public void partiallyInterpretatableRL1() throws ParsingException{
		TokenizationResult result =tokenizer.tokenizeRightToLeft(reverseParseRules, "ethyl-2H-foo|ene", true);
		assertFalse(result.isSuccessfullyTokenized());
		assertFalse(result.isFullyInterpretable());
		assertEquals("ethyl-2H-foo|ene", result.getUninterpretableName());
		assertEquals("ethyl-2H-foo|", result.getUnparsableName());
		assertEquals("ethyl-2H-foo|ene", result.getUnparsedName());
	}
	
	@Test
	public void partiallyInterpretatableRL2() throws ParsingException{
		TokenizationResult result =tokenizer.tokenizeRightToLeft(reverseParseRules, "fooylpyridine oxide", true);
		assertFalse(result.isSuccessfullyTokenized());
		assertFalse(result.isFullyInterpretable());
		assertEquals("fooyl", result.getUninterpretableName());
		assertEquals("f", result.getUnparsableName());//o as in the end of thio then oyl
		assertEquals("fooylpyridine", result.getUnparsedName());
	}

	@Test
	public void tokenizeDoesNotTokenizeUnTokenizableName() throws ParsingException{
		TokenizationResult result =tokenizer.tokenize("ethyl acet|foo toluene", true);
		assertFalse(result.isSuccessfullyTokenized());
	}

	@Test
	public void tokenizePreservesSpacesInUninterpretableNameLR1() throws ParsingException{
		TokenizationResult result =tokenizer.tokenize("ethyl acet|foo toluene", true);
		assertEquals("acet|foo toluene", result.getUninterpretableName());
	}

	@Test
	public void tokenizePreservesSpacesInUnparsableNameLR1() throws ParsingException{
		TokenizationResult result =tokenizer.tokenize("ethyl acet|foo toluene", true);
		assertEquals("|foo toluene", result.getUnparsableName());
	}

	@Test
	public void tokenizePreservesSpacesInUnparsedNameLR1() throws ParsingException{
		TokenizationResult result =tokenizer.tokenize("ethyl acet|foo toluene", true);
		assertEquals("acet|foo toluene", result.getUnparsedName());
	}
	
	@Test
	public void tokenizePreservesSpacesInUninterpretableNameLR2() throws ParsingException{
		TokenizationResult result =tokenizer.tokenize("eth yl acet|foo toluene", true);
		assertEquals("acet|foo toluene", result.getUninterpretableName());
	}

	@Test
	public void tokenizePreservesSpacesInUnparsableNameLR2() throws ParsingException{
		TokenizationResult result =tokenizer.tokenize("eth yl acet|foo toluene", true);
		assertEquals("|foo toluene", result.getUnparsableName());
	}

	@Test
	public void tokenizePreservesSpacesInUnparsedNameLR2() throws ParsingException{
		TokenizationResult result =tokenizer.tokenize("eth yl acet|foo toluene", true);
		assertEquals("acet|foo toluene", result.getUnparsedName());
	}

	@Test
	public void tokenizePreservesSpacesInUninterpretableNameRL1() throws ParsingException{
		TokenizationResult result =tokenizer.tokenizeRightToLeft(reverseParseRules, "ethyl foo|yl toluene", true);
		assertEquals("ethyl foo|yl", result.getUninterpretableName());
	}

	@Test
	public void tokenizePreservesSpacesInUnparsableNameRL1() throws ParsingException{
		TokenizationResult result =tokenizer.tokenizeRightToLeft(reverseParseRules, "ethyl foo|yl toluene", true);
		assertEquals("ethyl foo|", result.getUnparsableName());
	}

	@Test
	public void tokenizePreservesSpacesInUnparsedNameRL1() throws ParsingException{
		TokenizationResult result =tokenizer.tokenizeRightToLeft(reverseParseRules, "ethyl foo|yl toluene", true);
		assertEquals("ethyl foo|yl", result.getUnparsedName());
	}
	
	@Test
	public void tokenizePreservesSpacesInUninterpretableNameRL2() throws ParsingException{
		TokenizationResult result =tokenizer.tokenizeRightToLeft(reverseParseRules, "ethyl foo|yl tolu ene", true);
		assertEquals("ethyl foo|yl", result.getUninterpretableName());
	}

	@Test
	public void tokenizePreservesSpacesInUnparsableNameRL2() throws ParsingException{
		TokenizationResult result =tokenizer.tokenizeRightToLeft(reverseParseRules, "ethyl foo|yl tolu ene", true);
		assertEquals("ethyl foo|", result.getUnparsableName());
	}

	@Test
	public void tokenizePreservesSpacesInUnparsedNameRL2() throws ParsingException{
		TokenizationResult result =tokenizer.tokenizeRightToLeft(reverseParseRules, "ethyl foo|yl tolu ene", true);
		assertEquals("ethyl foo|yl", result.getUnparsedName());
	}
}
