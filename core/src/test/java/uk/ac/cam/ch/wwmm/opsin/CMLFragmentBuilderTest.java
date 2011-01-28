package uk.ac.cam.ch.wwmm.opsin;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertNotNull;
import static org.mockito.Mockito.mock;

import java.io.IOException;

import nu.xom.Builder;
import nu.xom.Element;

import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;


public class CMLFragmentBuilderTest {

	private BuildState state;
	private static CMLFragmentBuilder cmlBuilder;

	@BeforeClass
	public static void setUp() throws IOException{
		cmlBuilder = new CMLFragmentBuilder(new ResourceGetter("uk/ac/cam/ch/wwmm/opsin/resources/"));
	}

	@Before
	public void setUpInstance(){
		state = new BuildState(mock(NameToStructureConfig.class), mock(SMILESFragmentBuilder.class), cmlBuilder);
	}
	
	@AfterClass
	public static void cleanUp(){
		cmlBuilder = null;
	}

	@Test
	public void testBuildStringString() throws StructureBuildingException {
		Fragment benzene = state.fragManager.buildCML("benzene", "null", "null");
		StructureBuilder.makeHydrogensExplicit(state);
		assertNotNull("Got benzene", benzene);
		assertEquals("Benzene is correct", "<cml xmlns=\"http://www.xml-cml.org/schema\" " +
				"xmlns:conventions=\"http://www.xml-cml.org/convention/\" " +
				"xmlns:cmlDict=\"http://www.xml-cml.org/dictionary/cml/\" " +
				"xmlns:nameDict=\"http://www.xml-cml.org/dictionary/cml/name/\" " +
				"convention=\"conventions:molecular\"><molecule id=\"m1\">" +
				"<name dictRef=\"nameDict:unknown\">benzene</name>" +
				"<atomArray><atom id=\"a1\" elementType=\"C\"><label value=\"1\" dictRef=\"cmlDict:locant\" /></atom>" +
				"<atom id=\"a2\" elementType=\"C\"><label value=\"2\" dictRef=\"cmlDict:locant\" /></atom>" +
				"<atom id=\"a3\" elementType=\"C\"><label value=\"3\" dictRef=\"cmlDict:locant\" /></atom>" +
				"<atom id=\"a4\" elementType=\"C\"><label value=\"4\" dictRef=\"cmlDict:locant\" /></atom>" +
				"<atom id=\"a5\" elementType=\"C\"><label value=\"5\" dictRef=\"cmlDict:locant\" /></atom>" +
				"<atom id=\"a6\" elementType=\"C\"><label value=\"6\" dictRef=\"cmlDict:locant\" /></atom>" +
	            "<atom id=\"a7\" elementType=\"H\" />" +
	            "<atom id=\"a8\" elementType=\"H\" />" +
	            "<atom id=\"a9\" elementType=\"H\" />" +
	            "<atom id=\"a10\" elementType=\"H\" />" +
	            "<atom id=\"a11\" elementType=\"H\" />" +
	            "<atom id=\"a12\" elementType=\"H\" />" +
				"</atomArray><bondArray>" +
				"<bond id=\"a1_a2\" atomRefs2=\"a1 a2\" order=\"D\" />" +
				"<bond id=\"a2_a3\" atomRefs2=\"a2 a3\" order=\"S\" />" +
				"<bond id=\"a3_a4\" atomRefs2=\"a3 a4\" order=\"D\" />" +
				"<bond id=\"a4_a5\" atomRefs2=\"a4 a5\" order=\"S\" />" +
				"<bond id=\"a5_a6\" atomRefs2=\"a5 a6\" order=\"D\" />" +
				"<bond id=\"a6_a1\" atomRefs2=\"a6 a1\" order=\"S\" />" +
	            "<bond id=\"a1_a7\" atomRefs2=\"a1 a7\" order=\"S\" />" +
	            "<bond id=\"a2_a8\" atomRefs2=\"a2 a8\" order=\"S\" />" +
	            "<bond id=\"a3_a9\" atomRefs2=\"a3 a9\" order=\"S\" />" +
	            "<bond id=\"a4_a10\" atomRefs2=\"a4 a10\" order=\"S\" />" +
	            "<bond id=\"a5_a11\" atomRefs2=\"a5 a11\" order=\"S\" />" +
	            "<bond id=\"a6_a12\" atomRefs2=\"a6 a12\" order=\"S\" />" +
				"</bondArray></molecule></cml>", benzene.toCMLMolecule("benzene").toXML());
	}
	
	@Test
	public void testAtom() throws Exception {
		Element cmlAtom = new Builder().build("<atom id=\"a10\" elementType=\"C\" formalCharge=\"-1\">" +
				"<label value=\"1\" /></atom>", "/localhost").getRootElement();
		Atom atom = cmlBuilder.buildAtomFromCML(new FragmentManager(mock(SMILESFragmentBuilder.class), cmlBuilder, new IDManager()), cmlAtom, mock(Fragment.class));
		assertNotNull("Got atom", atom);
		assertEquals("Id = 1", 1, atom.getID());
		assertEquals("Element = C", "C", atom.getElement());
		assertEquals("Charge = -1", -1, atom.getCharge());
	}

}
