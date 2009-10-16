package uk.ac.cam.ch.wwmm.opsin;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertNotNull;
import static junit.framework.Assert.assertNull;
import nu.xom.Element;

import org.junit.Test;


public class NameToStructureTest_FromPapers {

	@Test
	public void testNameToStructure() throws Exception {
		NameToStructure nts = new NameToStructure();
		assertNotNull("Got a name to structure convertor", nts);
	}

	@Test
	public void testParseToCML() throws Exception {
		NameToStructure nts = new NameToStructure();
		Element cml = nts.parseToCML("ethane");
		assertEquals("Parsed ethane OK", "<cml xmlns=\"http://www.xml-cml.org/schema\"><molecule id=\"m1\"><atomArray>" +
				"<atom id=\"a1\" elementType=\"C\"><label value=\"1\" /></atom>" +
				"<atom id=\"a2\" elementType=\"C\"><label value=\"2\" /></atom>" +
				"<atom id=\"a3\" elementType=\"H\" />" +
				"<atom id=\"a4\" elementType=\"H\" />" +
				"<atom id=\"a5\" elementType=\"H\" />" +
				"<atom id=\"a6\" elementType=\"H\" />" +
				"<atom id=\"a7\" elementType=\"H\" />" +
				"<atom id=\"a8\" elementType=\"H\" />" +
				"</atomArray><bondArray>" +
				"<bond atomRefs2=\"a1 a2\" order=\"1\" />" +
	            "<bond atomRefs2=\"a1 a3\" order=\"1\" />" +
	            "<bond atomRefs2=\"a1 a4\" order=\"1\" />" +
	            "<bond atomRefs2=\"a1 a5\" order=\"1\" />" +
	            "<bond atomRefs2=\"a2 a6\" order=\"1\" />" +
	            "<bond atomRefs2=\"a2 a7\" order=\"1\" />" +
	            "<bond atomRefs2=\"a2 a8\" order=\"1\" />" +
				"</bondArray></molecule></cml>", cml.toXML());
		assertNull("Won't parse helloworld", nts.parseToCML("helloworld"));
	}
}
