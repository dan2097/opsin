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
		assertEquals("Benzene is correct", "<cml xmlns=\"http://www.xml-cml.org/schema\"><molecule id=\"m1\">" +
				"<atomArray><atom id=\"a1\" elementType=\"C\"><label value=\"1\" /></atom>" +
				"<atom id=\"a2\" elementType=\"C\"><label value=\"2\" /></atom>" +
				"<atom id=\"a3\" elementType=\"C\"><label value=\"3\" /></atom>" +
				"<atom id=\"a4\" elementType=\"C\"><label value=\"4\" /></atom>" +
				"<atom id=\"a5\" elementType=\"C\"><label value=\"5\" /></atom>" +
				"<atom id=\"a6\" elementType=\"C\"><label value=\"6\" /></atom>" +
				"</atomArray><bondArray>" +
				"<bond atomRefs2=\"a1 a2\" order=\"2\" />" +
				"<bond atomRefs2=\"a2 a3\" order=\"1\" />" +
				"<bond atomRefs2=\"a3 a4\" order=\"2\" />" +
				"<bond atomRefs2=\"a4 a5\" order=\"1\" />" +
				"<bond atomRefs2=\"a5 a6\" order=\"2\" />" +
				"<bond atomRefs2=\"a6 a1\" order=\"1\" />" +
				"</bondArray></molecule></cml>", benzene.toCMLMolecule().toXML());
	}

}
