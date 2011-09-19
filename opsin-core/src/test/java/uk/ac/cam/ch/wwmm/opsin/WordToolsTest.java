package uk.ac.cam.ch.wwmm.opsin;

import static junit.framework.Assert.assertEquals;

import java.util.Arrays;
import java.util.List;

import org.junit.Test;

public class WordToolsTest {

	@Test
	public void testNormalCase() throws ParsingException {
		ParseTokens pTokens = new ParseTokens(Arrays.asList("fooane",""), Arrays.asList('a', OpsinTools.END_OF_MAINGROUP));
		List<ParseWord> parseWords = WordTools.splitIntoParseWords(Arrays.asList(pTokens), "fooane");
		assertEquals(1, parseWords.size());
		assertEquals("fooane", parseWords.get(0).getWord());
		assertEquals(1, parseWords.get(0).getParseTokens().size());
		assertEquals(pTokens, parseWords.get(0).getParseTokens().get(0));
	}
	
	@Test
	public void testNormalCase2() throws ParsingException {
		ParseTokens pTokens = new ParseTokens(Arrays.asList("fooyl","","fooane",""), Arrays.asList('a', OpsinTools.END_OF_SUBSTITUENT,'a', OpsinTools.END_OF_MAINGROUP));
		List<ParseWord> parseWords = WordTools.splitIntoParseWords(Arrays.asList(pTokens), "fooylfooane");
		assertEquals(1, parseWords.size());
		assertEquals("fooylfooane", parseWords.get(0).getWord());
		assertEquals(1, parseWords.get(0).getParseTokens().size());
		assertEquals(pTokens, parseWords.get(0).getParseTokens().get(0));
	}
	
	@Test
	public void testNormalCase3() throws ParsingException {
		ParseTokens pTokens = new ParseTokens(Arrays.asList("functionalfoo",""), Arrays.asList('a', OpsinTools.END_OF_FUNCTIONALTERM));
		List<ParseWord> parseWords = WordTools.splitIntoParseWords(Arrays.asList(pTokens), "functionalfoo");
		assertEquals(1, parseWords.size());
		assertEquals("functionalfoo", parseWords.get(0).getWord());
		assertEquals(1, parseWords.get(0).getParseTokens().size());
		assertEquals(pTokens, parseWords.get(0).getParseTokens().get(0));
	}
	
	@Test
	public void testNormalCase4() throws ParsingException {
		ParseTokens pTokens1 = new ParseTokens(Arrays.asList("fooane",""), Arrays.asList('a', OpsinTools.END_OF_MAINGROUP));
		ParseTokens pTokens2 = new ParseTokens(Arrays.asList("fooane",""), Arrays.asList('b', OpsinTools.END_OF_MAINGROUP));
		List<ParseWord> parseWords = WordTools.splitIntoParseWords(Arrays.asList(pTokens1, pTokens2), "fooane");
		assertEquals(1, parseWords.size());
		assertEquals("fooane", parseWords.get(0).getWord());
		assertEquals(2, parseWords.get(0).getParseTokens().size());
		assertEquals(pTokens1, parseWords.get(0).getParseTokens().get(0));
		assertEquals(pTokens2, parseWords.get(0).getParseTokens().get(1));
	}

	@Test
	public void testStartingFunctionalTerm1() throws ParsingException {
		ParseTokens pTokens = new ParseTokens(Arrays.asList("poly","","foo",""), Arrays.asList('a', OpsinTools.END_OF_FUNCTIONALTERM,'a', OpsinTools.END_OF_SUBSTITUENT));
		List<ParseWord> parseWords = WordTools.splitIntoParseWords(Arrays.asList(pTokens), "polyfoo");
		assertEquals(2, parseWords.size());
		assertEquals("poly", parseWords.get(0).getWord());
		assertEquals(1, parseWords.get(0).getParseTokens().size());
		ParseTokens pTokensFunc = new ParseTokens(Arrays.asList("poly",""), Arrays.asList('a', OpsinTools.END_OF_FUNCTIONALTERM));
		assertEquals(pTokensFunc, parseWords.get(0).getParseTokens().get(0));
		assertEquals("foo", parseWords.get(1).getWord());
		assertEquals(1, parseWords.get(1).getParseTokens().size());
		ParseTokens pTokensGroup = new ParseTokens(Arrays.asList("foo",""), Arrays.asList('a', OpsinTools.END_OF_SUBSTITUENT));
		assertEquals(pTokensGroup, parseWords.get(1).getParseTokens().get(0));
	}
	
	@Test
	public void testStartingFunctionalTerm2() throws ParsingException {
		ParseTokens pTokens = new ParseTokens(Arrays.asList("poly","","foo",""), Arrays.asList('a', OpsinTools.END_OF_FUNCTIONALTERM,'a', OpsinTools.END_OF_MAINGROUP));
		List<ParseWord> parseWords = WordTools.splitIntoParseWords(Arrays.asList(pTokens), "polyfoo");
		assertEquals(2, parseWords.size());
		assertEquals("poly", parseWords.get(0).getWord());
		assertEquals(1, parseWords.get(0).getParseTokens().size());
		ParseTokens pTokensFunc = new ParseTokens(Arrays.asList("poly",""), Arrays.asList('a', OpsinTools.END_OF_FUNCTIONALTERM));
		assertEquals(pTokensFunc, parseWords.get(0).getParseTokens().get(0));
		assertEquals("foo", parseWords.get(1).getWord());
		assertEquals(1, parseWords.get(1).getParseTokens().size());
		ParseTokens pTokensGroup = new ParseTokens(Arrays.asList("foo",""), Arrays.asList('a', OpsinTools.END_OF_MAINGROUP));
		assertEquals(pTokensGroup, parseWords.get(1).getParseTokens().get(0));
	}
	
	@Test
	public void testTerminalFunctionalTerm() throws ParsingException {
		ParseTokens pTokens = new ParseTokens(Arrays.asList("fooyl","","functionalfoo",""), Arrays.asList('a', OpsinTools.END_OF_SUBSTITUENT,'a', OpsinTools.END_OF_FUNCTIONALTERM));
		List<ParseWord> parseWords = WordTools.splitIntoParseWords(Arrays.asList(pTokens), "fooylfunctionalfoo");
		assertEquals(2, parseWords.size());
		assertEquals("fooyl", parseWords.get(0).getWord());
		assertEquals(1, parseWords.get(0).getParseTokens().size());
		ParseTokens pTokensSub = new ParseTokens(Arrays.asList("fooyl",""), Arrays.asList('a', OpsinTools.END_OF_SUBSTITUENT));
		assertEquals(pTokensSub, parseWords.get(0).getParseTokens().get(0));
		assertEquals("functionalfoo", parseWords.get(1).getWord());
		assertEquals(1, parseWords.get(1).getParseTokens().size());
		ParseTokens pTokensFunc = new ParseTokens(Arrays.asList("functionalfoo",""), Arrays.asList('a', OpsinTools.END_OF_FUNCTIONALTERM));
		assertEquals(pTokensFunc, parseWords.get(1).getParseTokens().get(0));
	}
	
	@Test
	public void testMultipleParsesTerminalFunctionalTerm() throws ParsingException {
		ParseTokens pTokens1 = new ParseTokens(Arrays.asList("fooyl","","functionalfoo",""), Arrays.asList('a', OpsinTools.END_OF_SUBSTITUENT,'a', OpsinTools.END_OF_FUNCTIONALTERM));
		ParseTokens pTokens2 = new ParseTokens(Arrays.asList("fooyl","","functionalfoo",""), Arrays.asList('b', OpsinTools.END_OF_SUBSTITUENT,'a', OpsinTools.END_OF_FUNCTIONALTERM));
		List<ParseWord> parseWords = WordTools.splitIntoParseWords(Arrays.asList(pTokens1,pTokens2), "fooylfunctionalfoo");
		assertEquals(2, parseWords.size());
		assertEquals("fooyl", parseWords.get(0).getWord());
		assertEquals(2, parseWords.get(0).getParseTokens().size());
		ParseTokens pTokensSub1 = new ParseTokens(Arrays.asList("fooyl",""), Arrays.asList('a', OpsinTools.END_OF_SUBSTITUENT));
		assertEquals(pTokensSub1, parseWords.get(0).getParseTokens().get(0));
		ParseTokens pTokensSub2 = new ParseTokens(Arrays.asList("fooyl",""), Arrays.asList('b', OpsinTools.END_OF_SUBSTITUENT));
		assertEquals(pTokensSub2, parseWords.get(0).getParseTokens().get(1));
		assertEquals("functionalfoo", parseWords.get(1).getWord());
		assertEquals(1, parseWords.get(1).getParseTokens().size());
		ParseTokens pTokensFunc = new ParseTokens(Arrays.asList("functionalfoo",""), Arrays.asList('a', OpsinTools.END_OF_FUNCTIONALTERM));
		assertEquals(pTokensFunc, parseWords.get(1).getParseTokens().get(0));
	}
	
	@Test
	public void testMultipleParsesAmbiguousWordTokenisationTerminalFunctionalTerm() throws ParsingException {
		ParseTokens pTokens1 = new ParseTokens(Arrays.asList("fooyl","","functionalfoo",""), Arrays.asList('a', OpsinTools.END_OF_SUBSTITUENT,'a', OpsinTools.END_OF_FUNCTIONALTERM));
		ParseTokens pTokens2 = new ParseTokens(Arrays.asList("fooylfunc","","tionalfoo",""), Arrays.asList('b', OpsinTools.END_OF_SUBSTITUENT,'a', OpsinTools.END_OF_FUNCTIONALTERM));
		List<ParseWord> parseWords = WordTools.splitIntoParseWords(Arrays.asList(pTokens1,pTokens2), "fooylfunctionalfoo");
		assertEquals(2, parseWords.size());
		assertEquals("fooyl", parseWords.get(0).getWord());
		assertEquals(1, parseWords.get(0).getParseTokens().size());
		ParseTokens pTokensSub = new ParseTokens(Arrays.asList("fooyl",""), Arrays.asList('a', OpsinTools.END_OF_SUBSTITUENT));
		assertEquals(pTokensSub, parseWords.get(0).getParseTokens().get(0));
		assertEquals("functionalfoo", parseWords.get(1).getWord());
		assertEquals(1, parseWords.get(1).getParseTokens().size());
		ParseTokens pTokensFunc = new ParseTokens(Arrays.asList("functionalfoo",""), Arrays.asList('a', OpsinTools.END_OF_FUNCTIONALTERM));
		assertEquals(pTokensFunc, parseWords.get(1).getParseTokens().get(0));
	}
}
