package uk.ac.cam.ch.wwmm.opsin;

import javax.xml.stream.XMLStreamConstants;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import static org.junit.Assert.*;
import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

public class VerifyFragmentsTest {
	
	private static ResourceGetter resourceGetter;
	private static SMILESFragmentBuilder sBuilder;

	@BeforeClass
	public static void setUp() {
		resourceGetter = new ResourceGetter("uk/ac/cam/ch/wwmm/opsin/resources/");
		sBuilder = new SMILESFragmentBuilder(new IDManager());
	}
	
	@AfterClass
	public static void cleanUp(){
		resourceGetter = null;
		sBuilder = null;
	}
	
	@Test
	public void verifySMILES() throws Exception {

		XMLStreamReader indexReader = resourceGetter.getXMLStreamReader("index.xml");
		while (indexReader.hasNext()) {
			if (indexReader.next() == XMLStreamConstants.START_ELEMENT &&
					indexReader.getLocalName().equals("tokenFile")) {
				XMLStreamReader tokenReader = resourceGetter.getXMLStreamReader(indexReader.getElementText());
				while (tokenReader.hasNext()) {
					if (tokenReader.next() == XMLStreamConstants.START_ELEMENT) {
						String tagName = tokenReader.getLocalName();
						if (tagName.equals("tokenLists")) {
							while (tokenReader.hasNext()) {
								switch (tokenReader.next()) {
								case XMLStreamConstants.START_ELEMENT:
									if (tokenReader.getLocalName().equals("tokenList")) {
										verifySmilesInTokenList(tokenReader);
									}
									break;
								}
							}
						}
						else if (tagName.equals("tokenList")) {
							verifySmilesInTokenList(tokenReader);
						}
					}
				}
			}
		}
		indexReader.close();
	}

	private void verifySmilesInTokenList(XMLStreamReader reader) throws XMLStreamException {
		String tagname = reader.getAttributeValue(null, "tagname");
		if (tagname.equals(GROUP_EL) ||
				tagname.equals(FUNCTIONALGROUP_EL) ||
				tagname.equals(HETEROATOM_EL) ||
				tagname.equals(SUFFIXPREFIX_EL)) {
			while (reader.hasNext()) {
				switch (reader.next()) {
				case XMLStreamConstants.START_ELEMENT:
					if (reader.getLocalName().equals("token")) {
						Fragment mol = null;
						String smiles = null;
						try{
							smiles = reader.getAttributeValue(null, VALUE_ATR);
							String type = reader.getAttributeValue(null, TYPE_ATR);
							String subType = reader.getAttributeValue(null, SUBTYPE_ATR);
							String labels =  reader.getAttributeValue(null, LABELS_ATR);
							TokenEl tokenEl = new TokenEl(GROUP_EL);
							if (type != null){
								tokenEl.addAttribute(TYPE_ATR, type);
							}
							if (subType != null){
								tokenEl.addAttribute(SUBTYPE_ATR, subType);
							}

							mol = sBuilder.build(smiles, tokenEl, labels != null ? labels : "");
						}
						catch (Exception e) {
							e.printStackTrace();
						}
						assertNotNull("The following token's SMILES or labels were in error: " + smiles, mol);
						try{
							mol.checkValencies();
						}
						catch (StructureBuildingException e) {
							fail("The following token's SMILES produced a structure with invalid valency: " + smiles);
						}
					}
					break;
				case XMLStreamConstants.END_ELEMENT:
					if (reader.getLocalName().equals("tokenList")) {
						return;
					}
					break;
				}
			}
		}
	}
}
