package uk.ac.cam.ch.wwmm.opsin;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertFalse;
import static junit.framework.Assert.assertNotNull;
import static junit.framework.Assert.assertTrue;

import java.util.List;

import org.junit.Test;



public class ParseRulesTest {

	@Test
	public void testParseRules() throws Exception {
		ResourceGetter rg = new ResourceGetter("uk/ac/cam/ch/wwmm/opsin/resources/");
		ParseRules pr = new ParseRules(new TokenManager(rg),rg);
		assertNotNull("Got ParseRules", pr);
		String regex = pr.regexDict.get("%chemical%");
		assertNotNull("Got regex", regex);
		assertTrue("Regex is a decent size", regex.length() > 10);
		assertFalse("Regex contains no %", regex.contains("%"));

		TwoReturnValues<List<ParseTokens>,Boolean> returned = pr.getParses("hexane", "simple");
		List<ParseTokens> parseTokenList =returned.getFirst();
		boolean interSubHyphenAllowed =returned.getSecond();
		assertEquals("One lex", 1, parseTokenList.size());
		assertEquals("Two tokens", 2, parseTokenList.get(0).tokens.size());
		assertEquals("First token: hex", "hex", parseTokenList.get(0).tokens.get(0));
		assertEquals("Second token: ane", "ane", parseTokenList.get(0).tokens.get(1));
		assertEquals(false, interSubHyphenAllowed);

		returned = pr.getParses("hexachlorohexane", "simple");
		parseTokenList =returned.getFirst();
		interSubHyphenAllowed =returned.getSecond();
		assertEquals("Four tokens", 4, parseTokenList.get(0).tokens.size());
		assertEquals("First token: hexa", "hexa", parseTokenList.get(0).tokens.get(0));
		assertEquals("Second token: chloro", "chloro", parseTokenList.get(0).tokens.get(1));
		assertEquals("Third token: hex", "hex", parseTokenList.get(0).tokens.get(2));
		assertEquals("Fourth token: ane", "ane", parseTokenList.get(0).tokens.get(3));
		assertEquals(false, interSubHyphenAllowed);

		returned = pr.getParses("hexachlorohexaneeeeeee", "simple");
		parseTokenList =returned.getFirst();
		interSubHyphenAllowed =returned.getSecond();
		assertEquals("No Parses", 0, parseTokenList.size());
		assertEquals(false, interSubHyphenAllowed);

		returned = pr.getParses("(hexachloro)hexane", "simple");
		parseTokenList =returned.getFirst();
		interSubHyphenAllowed =returned.getSecond();
		assertEquals("One lex", 1, parseTokenList.size());
		assertEquals("Six tokens", 6, parseTokenList.get(0).tokens.size());
		assertEquals("token", "(", parseTokenList.get(0).tokens.get(0));
		assertEquals("token", "hexa", parseTokenList.get(0).tokens.get(1));
		assertEquals("token", "chloro", parseTokenList.get(0).tokens.get(2));
		assertEquals("token", ")", parseTokenList.get(0).tokens.get(3));
		assertEquals("token", "hex", parseTokenList.get(0).tokens.get(4));
		assertEquals("token", "ane", parseTokenList.get(0).tokens.get(5));
		assertEquals(false, interSubHyphenAllowed);

		returned = pr.getParses("methyl", "substituent");
		parseTokenList =returned.getFirst();
		interSubHyphenAllowed =returned.getSecond();
		assertEquals("One lex", 1, parseTokenList.size());
		assertEquals("Two tokens", 2, parseTokenList.get(0).tokens.size());
		assertEquals("token", "meth", parseTokenList.get(0).tokens.get(0));
		assertEquals("token", "yl", parseTokenList.get(0).tokens.get(1));
		assertEquals(false, interSubHyphenAllowed);

		returned = pr.getParses("methyl", "simple");
		parseTokenList =returned.getFirst();
		interSubHyphenAllowed =returned.getSecond();
		assertEquals("No Lexes", 0, parseTokenList.size());
		assertEquals(true, interSubHyphenAllowed);
	}
}
