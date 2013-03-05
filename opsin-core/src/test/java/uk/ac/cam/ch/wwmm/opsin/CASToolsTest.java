package uk.ac.cam.ch.wwmm.opsin;

import static junit.framework.Assert.assertEquals;

import java.io.IOException;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 *
 */
public class CASToolsTest {
	private static ParseRules parseRules;

	@BeforeClass
	public static void setUp() throws IOException{
		ResourceGetter rg = new ResourceGetter("uk/ac/cam/ch/wwmm/opsin/resources/");
		parseRules = new ParseRules(new ResourceManager(rg));
	}
	
	@AfterClass
	public static void cleanUp() {
		parseRules = null;
	}

	@Test
	public void cas1() throws ParsingException{
		String name = CASTools.uninvertCASName("Silane, chloromethyl-", parseRules);
		assertEquals("chloromethyl-Silane", name);
	}

	@Test
	public void cas2() throws ParsingException{
		String name = CASTools.uninvertCASName("Acetic acid, 2-ethoxy-2-thioxo-", parseRules);
		assertEquals("2-ethoxy-2-thioxo-Acetic acid", name);
	}

	@Test
	public void cas3() throws ParsingException{
		String name = CASTools.uninvertCASName("Silanol, 1,1'-methylenebis-", parseRules);
		assertEquals("1,1'-methylenebis-Silanol", name);
	}

	@Test
	public void cas4() throws ParsingException{
		String name = CASTools.uninvertCASName("Phosphonic acid, P,P'-(8-methylene-3,7,10,14-tetraoxo-4,6,11,13-tetraazahexadecane-1,16-diyl)-bis-, P,P,P',P'-tetramethyl ester", parseRules);
		assertEquals("P,P,P',P'-tetramethyl P,P'-(8-methylene-3,7,10,14-tetraoxo-4,6,11,13-tetraazahexadecane-1,16-diyl)-bis-Phosphonate", name);
	}

	@Test
	public void cas5() throws ParsingException{
		String name = CASTools.uninvertCASName("Benzenamine, 3,3',3''-(1-ethenyl-2-ylidene)tris[6-methyl-", parseRules);
		assertEquals("3,3',3''-(1-ethenyl-2-ylidene)tris[6-methyl-Benzenamine]", name);
	}

	@Test
	public void cas6() throws ParsingException{
		String name = CASTools.uninvertCASName("Pyridine, 3,3'-thiobis[6-chloro-", parseRules);
		assertEquals("3,3'-thiobis[6-chloro-Pyridine]", name);
	}

	@Test
	public void cas7() throws ParsingException{
		String name = CASTools.uninvertCASName("1-Butanesulfonic acid, 2,4-diamino-3-chloro- 1-ethyl ester", parseRules);
		assertEquals("1-ethyl 2,4-diamino-3-chloro-1-Butanesulfonate", name);
	}

	@Test
	public void cas8() throws ParsingException{
		String name = CASTools.uninvertCASName("Benzenecarboximidamide, N'-(1E)-1-propen-1-yl-N-(1Z)-1-propen-1-yl-", parseRules);
		assertEquals("N'-(1E)-1-propen-1-yl-N-(1Z)-1-propen-1-yl-Benzenecarboximidamide", name);
	}

	@Test
	public void cas9() throws ParsingException{
		String name = CASTools.uninvertCASName("Phosphoric acid, ethyl dimethyl ester", parseRules);
		assertEquals("ethyl dimethyl Phosphorate", name);
	}

	@Test
	public void cas10() throws ParsingException{
		String name = CASTools.uninvertCASName("2-Propanone, oxime", parseRules);
		assertEquals("2-Propanone oxime", name);
	}

	@Test
	public void cas11() throws ParsingException{
		String name = CASTools.uninvertCASName("Disulfide, bis(2-chloroethyl)", parseRules);
		assertEquals("bis(2-chloroethyl) Disulfide", name);
	}

	@Test
	public void cas12() throws ParsingException{
		String name = CASTools.uninvertCASName("Ethanimidic acid, N-nitro-, (1Z)-", parseRules);
		assertEquals("N-nitro-(1Z)-Ethanimidic acid", name);
	}

	@Test
	public void cas13() throws ParsingException{
		String name = CASTools.uninvertCASName("2(1H)-Pyridinone, hydrazone, (2E)-", parseRules);
		assertEquals("(2E)-2(1H)-Pyridinone hydrazone", name);
	}

	@Test
	public void cas14() throws ParsingException{
		String name = CASTools.uninvertCASName("benzoic acid, 4,4'-methylenebis[2-chloro-", parseRules);
		assertEquals("4,4'-methylenebis[2-chloro-benzoic acid]", name);
	}

	@Test
	public void cas15() throws ParsingException{
		String name = CASTools.uninvertCASName("peroxide, ethyl methyl", parseRules);
		assertEquals("ethyl methyl peroxide", name);
	}

	@Test
	public void cas16() throws ParsingException{
		String name = CASTools.uninvertCASName("Phosphonic diamide, P-phenyl- (8CI9CI)", parseRules);
		assertEquals("P-phenyl-Phosphonic diamide", name);
	}

	@Test
	public void cas17() throws ParsingException{
		String name = CASTools.uninvertCASName("piperazinium, 1,1-dimethyl-, 2,2,2-trifluoroacetate hydrochloride", parseRules);
		assertEquals("1,1-dimethyl-piperazinium 2,2,2-trifluoroacetate hydrochloride", name);
	}
	
	@Test
	public void cas18() throws ParsingException{
		String name = CASTools.uninvertCASName("Acetamide, ethylenebis(((ethyl)amino)-", parseRules);
		assertEquals("ethylenebis(((ethyl)amino)-Acetamide)", name);
	}
	
	@Test
	public void cas19() throws ParsingException{
		String name = CASTools.uninvertCASName("Benzenesulfonic acid, 4-amino-, 1-methylhydrazide", parseRules);
		assertEquals("4-amino-Benzenesulfonic acid 1-methylhydrazide", name);
	}
	
	@Test
	public void cas20() throws ParsingException{
		String name = CASTools.uninvertCASName("Acetaldehyde, O-methyloxime", parseRules);
		assertEquals("Acetaldehyde O-methyloxime", name);
	}
	
	@Test
	public void cas21() throws ParsingException{
		String name = CASTools.uninvertCASName("Acetic acid, 2-amino-2-oxo-, 2-(phenylmethylene)hydrazide", parseRules);
		assertEquals("2-amino-2-oxo-Acetic acid 2-(phenylmethylene)hydrazide", name);
	}
	
	@Test
	public void cas22() throws ParsingException{
		String name = CASTools.uninvertCASName("L-Alanine, N-carboxy-, 1-ethyl ester", parseRules);
		assertEquals("1-ethyl N-carboxy-L-Alaninate", name);
	}

	@Test(expected=ParsingException.class)
	public void notCas1() throws ParsingException{
		CASTools.uninvertCASName("hexanamine, hexylamine", parseRules);
	}

	@Test(expected=ParsingException.class)
	public void notCas2() throws ParsingException{
		CASTools.uninvertCASName("cyclopropane-1,2-diyldicarbonyl diisocyanate, cyclopropane-1,2-diylbis(carbonyl)bisisocyanate", parseRules);
	}
	
	@Test(expected=ParsingException.class)
	public void notCas3() throws ParsingException{
		CASTools.uninvertCASName("benzoic acid, ester", parseRules);
	}
}
