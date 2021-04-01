package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.Assert.*;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import uk.ac.cam.ch.wwmm.opsin.OpsinResult.OPSIN_RESULT_STATUS;

public class InchiOutputTest {

	private static NameToInchi n2i;
	@BeforeClass
	public static void setUp() {
		n2i = new NameToInchi();
	}
	
	@AfterClass
	public static void cleanUp(){
		n2i = null;
	}
	
	@Test
	public void testStaticToInChI() throws StructureBuildingException{
		SMILESFragmentBuilder sBuilder = new SMILESFragmentBuilder(new IDManager());
		Fragment f = sBuilder.build("C([H])([H])([H])C(=O)N([H])[H]");
		OpsinResult result = new OpsinResult(f, OPSIN_RESULT_STATUS.SUCCESS, "", "");
		assertEquals("InChI=1/C2H5NO/c1-2(3)4/h1H3,(H2,3,4)/f/h3H2", NameToInchi.convertResultToInChI(result));
	}
	
	@Test
	public void testStaticToStdInChI() throws StructureBuildingException{
		SMILESFragmentBuilder sBuilder = new SMILESFragmentBuilder(new IDManager());
		Fragment f = sBuilder.build("C([H])([H])([H])C(=O)N([H])[H]");
		OpsinResult result = new OpsinResult(f, OPSIN_RESULT_STATUS.SUCCESS, "", "");
		assertEquals("InChI=1S/C2H5NO/c1-2(3)4/h1H3,(H2,3,4)", NameToInchi.convertResultToStdInChI(result));
	}
	
	@Test
	public void testStaticToStdInChIKey() throws StructureBuildingException{
		SMILESFragmentBuilder sBuilder = new SMILESFragmentBuilder(new IDManager());
		Fragment f = sBuilder.build("C([H])([H])([H])C(=O)N([H])[H]");
		OpsinResult result = new OpsinResult(f, OPSIN_RESULT_STATUS.SUCCESS, "", "");
		assertEquals("DLFVBJFMPXGRIB-UHFFFAOYSA-N", NameToInchi.convertResultToStdInChIKey(result));
	}
	
	@Test
	public void testParseToInChI(){
		assertEquals("InChI=1/C2H5NO/c1-2(3)4/h1H3,(H2,3,4)/f/h3H2", n2i.parseToInchi("acetamide"));
	}
	
	
	@Test
	public void testParseToStdInChI(){
		assertEquals("InChI=1S/C2H5NO/c1-2(3)4/h1H3,(H2,3,4)", n2i.parseToStdInchi("acetamide"));
	}
	
	@Test
	public void testParseToStdInChIKey(){
		assertEquals("DLFVBJFMPXGRIB-UHFFFAOYSA-N", n2i.parseToStdInchiKey("acetamide"));
	}

	@Test
	public void ignoreRacemicStereoInInchi() throws StructureBuildingException{
		assertEquals("InChI=1/C8H10O/c1-7(9)8-5-3-2-4-6-8/h2-7,9H,1H3", n2i.parseToInchi("rac-(R)-1-phenylethan-1-ol"));
	}

	// more than one in same rac group
	@Test
	public void keepRacemicStereoInInchi() throws StructureBuildingException{
		assertEquals("InChI=1/C12H10O2/c13-11-8-5-1-3-7-4-2-6-9(10(7)8)12(11)14/h1-6,11-14H/t11-,12-/m1/s1",
				n2i.parseToInchi("rac-trans-acenaphthene-1,2-diol"));
	}

	@Test
	public void consistency() throws StructureBuildingException{
		assertEquals("InChI=1/C31H29Cl2F2N3O3/c1-30(2,3)14-20-16-38(29(41)37-15-18-4-6-19(7-5-18)28(39)40)27(23-10-8-21(32)12-25(23)34)31(20,17-36)24-11-9-22(33)13-26(24)35/h4-13,20,27H,14-16H2,1-3H3,(H,37,41)(H,39,40)/t20-,27-,31-/m1/s1/f/h37,39H",
				n2i.parseToInchi("4-({[rac-(2S,3S,4S)-2,3-bis-(4-chloro-2-fluoro-phenyl)-3-cyano-4-(2,2-dimethyl-propyl)-pyrrolidine-1-carbonyl]-amino}-methyl)-benzoic acid"));
		assertEquals("InChI=1/C31H29Cl2F2N3O3/c1-30(2,3)14-20-16-38(29(41)37-15-18-4-6-19(7-5-18)28(39)40)27(23-10-8-21(32)12-25(23)34)31(20,17-36)24-11-9-22(33)13-26(24)35/h4-13,20,27H,14-16H2,1-3H3,(H,37,41)(H,39,40)/t20-,27-,31-/m1/s1/f/h37,39H",
				n2i.parseToInchi("rac-4-({[(2S,3S,4S)-2,3-bis-(4-chloro-2-fluoro-phenyl)-3-cyano-4-(2,2-dimethyl-propyl)-pyrrolidine-1-carbonyl]-amino}-methyl)-benzoic acid"));
	}
}
