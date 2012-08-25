package uk.ac.cam.ch.wwmm.opsin;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertNotNull;
import static junit.framework.Assert.assertNull;
import nu.xom.Element;

import org.junit.Test;


public class NameToStructureTest {

	@Test
	public void testNameToStructure() {
		NameToStructure nts = NameToStructure.getInstance();
		assertNotNull("Got a name to structure convertor", nts);
	}

	@Test
	public void testParseToCML() {
		NameToStructure nts = NameToStructure.getInstance();
		Element cml = nts.parseToCML("ethane");
		// output is syntactically valid (schema, dictRefs)
		// labels assigned and is correct.
		// contains a molecule with same connectivity as 'frag of CML'

		assertEquals("Parsing 'ethane'", "<cml xmlns=\"http://www.xml-cml.org/schema\" " +
				"xmlns:conventions=\"http://www.xml-cml.org/convention/\" " +
				"xmlns:cmlDict=\"http://www.xml-cml.org/dictionary/cml/\" " +
				"xmlns:nameDict=\"http://www.xml-cml.org/dictionary/cml/name/\" " +
				"convention=\"conventions:molecular\"><molecule id=\"m1\">" +
				"<name dictRef=\"nameDict:unknown\">ethane</name><atomArray>" +
				"<atom id=\"a1\" elementType=\"C\"><label value=\"1\" dictRef=\"cmlDict:locant\" /><label value=\"alpha\" dictRef=\"cmlDict:locant\" /></atom>" +
				"<atom id=\"a2\" elementType=\"C\"><label value=\"2\" dictRef=\"cmlDict:locant\" /><label value=\"beta\" dictRef=\"cmlDict:locant\" /></atom>" +
				"<atom id=\"a3\" elementType=\"H\" />" +
				"<atom id=\"a4\" elementType=\"H\" />" +
				"<atom id=\"a5\" elementType=\"H\" />" +
				"<atom id=\"a6\" elementType=\"H\" />" +
				"<atom id=\"a7\" elementType=\"H\" />" +
				"<atom id=\"a8\" elementType=\"H\" />" +
				"</atomArray><bondArray>" +
				"<bond id=\"a1_a2\" atomRefs2=\"a1 a2\" order=\"S\" />" +
	            "<bond id=\"a1_a3\" atomRefs2=\"a1 a3\" order=\"S\" />" +
	            "<bond id=\"a1_a4\" atomRefs2=\"a1 a4\" order=\"S\" />" +
	            "<bond id=\"a1_a5\" atomRefs2=\"a1 a5\" order=\"S\" />" +
	            "<bond id=\"a2_a6\" atomRefs2=\"a2 a6\" order=\"S\" />" +
	            "<bond id=\"a2_a7\" atomRefs2=\"a2 a7\" order=\"S\" />" +
	            "<bond id=\"a2_a8\" atomRefs2=\"a2 a8\" order=\"S\" />" +
				"</bondArray></molecule></cml>", cml.toXML());
		assertNull("Won't parse helloworld", nts.parseToCML("helloworld"));
	}
	
	@Test
	public void testParseToSmiles() {
		NameToStructure nts = NameToStructure.getInstance();
		String smiles = nts.parseToSmiles("ethane");
		assertEquals("CC", smiles);
	}
}
