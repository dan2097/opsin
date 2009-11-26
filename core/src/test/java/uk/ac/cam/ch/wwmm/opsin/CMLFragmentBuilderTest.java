package uk.ac.cam.ch.wwmm.opsin;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertNotNull;

import org.junit.Before;
import org.junit.Test;


public class CMLFragmentBuilderTest {

	IDManager idManager;
	CMLFragmentBuilder cmlBuilder;

	@Before
	public void setUp() throws Exception {
		idManager = new IDManager();
		cmlBuilder = new CMLFragmentBuilder(new ResourceGetter("uk/ac/cam/ch/wwmm/opsin/resources/"));
	}

	@Test
	public void testBuildStringString() throws Exception {
		Fragment benzene = cmlBuilder.build("benzene", "null", "null", idManager);//hydrogens are implicit
		assertNotNull("Got benzene", benzene);
		assertEquals("Benzene is correct", "<cml xmlns=\"http://www.xml-cml.org/schema\" " +
				"xmlns:cmlDict=\"http://www.xml-cml.org/dictionary/cml/\" " +
				"xmlns:nameDict=\"http://www.xml-cml.org/dictionary/cml/name/\" " +
				"convention=\"cmlDict:cmllite\"><molecule id=\"m1\">" +
				"<name dictRef=\"nameDict:unknown\">benzene</name>" +
				"<atomArray><atom id=\"a1\" elementType=\"C\"><label value=\"1\" dictRef=\"cmlDict:locant\" /></atom>" +
				"<atom id=\"a2\" elementType=\"C\"><label value=\"2\" dictRef=\"cmlDict:locant\" /></atom>" +
				"<atom id=\"a3\" elementType=\"C\"><label value=\"3\" dictRef=\"cmlDict:locant\" /></atom>" +
				"<atom id=\"a4\" elementType=\"C\"><label value=\"4\" dictRef=\"cmlDict:locant\" /></atom>" +
				"<atom id=\"a5\" elementType=\"C\"><label value=\"5\" dictRef=\"cmlDict:locant\" /></atom>" +
				"<atom id=\"a6\" elementType=\"C\"><label value=\"6\" dictRef=\"cmlDict:locant\" /></atom>" +
				"</atomArray><bondArray>" +
				"<bond id=\"a1_a2\" atomRefs2=\"a1 a2\" order=\"D\" />" +
				"<bond id=\"a2_a3\" atomRefs2=\"a2 a3\" order=\"S\" />" +
				"<bond id=\"a3_a4\" atomRefs2=\"a3 a4\" order=\"D\" />" +
				"<bond id=\"a4_a5\" atomRefs2=\"a4 a5\" order=\"S\" />" +
				"<bond id=\"a5_a6\" atomRefs2=\"a5 a6\" order=\"D\" />" +
				"<bond id=\"a6_a1\" atomRefs2=\"a6 a1\" order=\"S\" />" +
				"</bondArray></molecule></cml>", benzene.toCMLMolecule("benzene").toXML());
	}

}
