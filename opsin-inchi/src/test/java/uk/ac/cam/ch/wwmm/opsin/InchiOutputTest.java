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
	public void testParseToInChI(){
		assertEquals("InChI=1/C2H5NO/c1-2(3)4/h1H3,(H2,3,4)/f/h3H2", n2i.parseToInchi("acetamide"));
	}
	
	
	@Test
	public void testParseToStdInChI(){
		assertEquals("InChI=1S/C2H5NO/c1-2(3)4/h1H3,(H2,3,4)", n2i.parseToStdInchi("acetamide"));
	}
}
