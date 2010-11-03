package uk.ac.cam.ch.wwmm.opsin;

import static junit.framework.Assert.assertEquals;

import java.util.List;

import org.junit.BeforeClass;
import org.junit.Test;

public class TokenizerTest {

	static Tokeniser tokenizer;

	@BeforeClass
	public static void  setUp() throws Exception {
		ResourceGetter rg = new ResourceGetter("uk/ac/cam/ch/wwmm/opsin/resources/");
		tokenizer = new Tokeniser(new ParseRules(new ResourceManager(rg)));
	}
	
	@Test
	public void ethylChloride() throws Exception {
		Parse parse = tokenizer.tokenize("ethyl chloride", true);
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


		parse = tokenizer.tokenize("ethylchloride", true);//missing space
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
		Parse parse = tokenizer.tokenize("hexane", true);
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
		Parse parse = tokenizer.tokenize("hexachlorohexane", true);
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
			tokenizer.tokenize("hexachlorohexaneeeeeee", true);
			failed = false;
		}
		catch (ParsingException e){
			failed = true;
		}
		assertEquals("Unparsable", true, failed);
	}

	@Test
	public void bracketedHexachlorohexane() throws Exception {
		Parse parse = tokenizer.tokenize("(hexachloro)hexane", true);
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
		Parse parse = tokenizer.tokenize("methyl", true);
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
		Parse parse = tokenizer.tokenize("acetic acid", true);
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
	public void cas1() throws Exception {
		String name = tokenizer.uninvertCASName("Silane, chloromethyl-");
		assertEquals("chloromethyl-Silane", name);
	}
	
	@Test
	public void cas2() throws Exception {
		String name = tokenizer.uninvertCASName("Acetic acid, 2-ethoxy-2-thioxo-");
		assertEquals("2-ethoxy-2-thioxo-Acetic acid", name);
	}

	@Test
	public void cas3() throws Exception {
		String name = tokenizer.uninvertCASName("Silanol, 1,1'-methylenebis-");
		assertEquals("1,1'-methylenebis-Silanol", name);
	}
	
	@Test
	public void cas4() throws Exception {
		String name = tokenizer.uninvertCASName("Phosphonic acid, P,P'-(8-methylene-3,7,10,14-tetraoxo-4,6,11,13-tetraazahexadecane-1,16-diyl)-bis-, P,P,P',P'-tetramethyl ester");
		assertEquals("P,P,P',P'-tetramethyl P,P'-(8-methylene-3,7,10,14-tetraoxo-4,6,11,13-tetraazahexadecane-1,16-diyl)-bis-Phosphonate", name);
	}

	@Test
	public void cas5() throws Exception {
		String name = tokenizer.uninvertCASName("Benzenamine, 3,3',3''-(1-ethenyl-2-ylidene)tris[6-methyl-");
		assertEquals("3,3',3''-(1-ethenyl-2-ylidene)tris[6-methyl-Benzenamine]", name);
	}
	
	@Test
	public void cas6() throws Exception {
		String name = tokenizer.uninvertCASName("Pyridine, 3,3'-thiobis[6-chloro-");
		assertEquals("3,3'-thiobis[6-chloro-Pyridine]", name);
	}
	
	@Test
	public void cas7() throws Exception {
		String name = tokenizer.uninvertCASName("1-Butanesulfonic acid, 2,4-diamino-3-chloro- 1-ethyl ester");
		assertEquals("1-ethyl 2,4-diamino-3-chloro-1-Butanesulfonate", name);
	}
	
	@Test
	public void cas8() throws Exception {
		String name = tokenizer.uninvertCASName("Benzenecarboximidamide, N'-(1E)-1-propen-1-yl-N-(1Z)-1-propen-1-yl-");
		assertEquals("N'-(1E)-1-propen-1-yl-N-(1Z)-1-propen-1-yl-Benzenecarboximidamide", name);
	}
	
	@Test
	public void cas9() throws Exception {
		String name = tokenizer.uninvertCASName("Phosphoric acid, ethyl dimethyl ester");
		assertEquals("ethyl dimethyl Phosphorate", name);
	}
	
	@Test
	public void cas10() throws Exception {
		String name = tokenizer.uninvertCASName("2-Propanone, oxime");
		assertEquals("2-Propanone oxime", name);
	}
	
	@Test
	public void cas11() throws Exception {
		String name = tokenizer.uninvertCASName("Disulfide, bis(2-chloroethyl)");
		assertEquals("bis(2-chloroethyl) Disulfide", name);
	}

	@Test
	public void cas12() throws Exception {
		String name = tokenizer.uninvertCASName("Ethanimidic acid, N-nitro-, (1Z)-");
		assertEquals("N-nitro-(1Z)-Ethanimidic acid", name);
	}
	
	@Test
	public void cas13() throws Exception {
		String name = tokenizer.uninvertCASName("2(1H)-Pyridinone, hydrazone, (2E)-");
		assertEquals("(2E)-2(1H)-Pyridinone hydrazone", name);
	}
	
	@Test
	public void cas14() throws Exception {
		String name = tokenizer.uninvertCASName("benzoic acid, 4,4'-methylenebis[2-chloro-");
		assertEquals("4,4'-methylenebis[2-chloro-benzoic acid]", name);
	}
	
	@Test
	public void cas15() throws Exception {
		String name = tokenizer.uninvertCASName("peroxide, ethyl methyl");
		assertEquals("ethyl methyl peroxide", name);
	}
	
	@Test
	public void cas16() throws Exception {
		String name = tokenizer.uninvertCASName("Phosphonic diamide, P-phenyl- (8CI9CI)");
		assertEquals("P-phenyl-Phosphonic diamide", name);
	}
	
	@Test
	public void cas17() throws Exception {
		String name = tokenizer.uninvertCASName("piperazinium, 1,1-dimethyl-, 2,2,2-trifluoroacetate hydrochloride");
		assertEquals("1,1-dimethyl-piperazinium 2,2,2-trifluoroacetate hydrochloride", name);
	}
	
	@Test(expected=ParsingException.class)
	public void notCas1() throws Exception {
		tokenizer.uninvertCASName("hexanamine, hexylamine");
	}
	
	@Test(expected=ParsingException.class)
	public void notCas2() throws Exception {
		tokenizer.uninvertCASName("cyclopropane-1,2-diyldicarbonyl diisocyanate, cyclopropane-1,2-diylbis(carbonyl)bisisocyanate");
	}
}
