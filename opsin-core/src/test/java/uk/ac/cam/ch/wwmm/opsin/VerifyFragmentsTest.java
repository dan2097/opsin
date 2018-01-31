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
	private static FragmentManager fm;

	@BeforeClass
	public static void setUp() {
		resourceGetter = new ResourceGetter("uk/ac/cam/ch/wwmm/opsin/resources/");
		IDManager idManager = new IDManager();
		fm = new FragmentManager(new SMILESFragmentBuilder(idManager), idManager);
	}
	
	@AfterClass
	public static void cleanUp(){
		resourceGetter = null;
		fm = null;
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
			String type = reader.getAttributeValue(null, TYPE_ATR);
			String subType = reader.getAttributeValue(null, SUBTYPE_ATR);
			while (reader.hasNext()) {
				switch (reader.next()) {
				case XMLStreamConstants.START_ELEMENT:
					if (reader.getLocalName().equals("token")) {
						String smiles = reader.getAttributeValue(null, VALUE_ATR);
						String labels =  reader.getAttributeValue(null, LABELS_ATR);
						TokenEl tokenEl = new TokenEl(GROUP_EL);
						if (type != null){
							tokenEl.addAttribute(TYPE_ATR, type);
						}
						if (subType != null){
							tokenEl.addAttribute(SUBTYPE_ATR, subType);
						}
						Fragment mol = null;
						try {
							mol = fm.buildSMILES(smiles, tokenEl, labels != null ? labels : "");
							fm.convertSpareValenciesToDoubleBonds();
							fm.makeHydrogensExplicit();
							if (!tagname.equals(HETEROATOM_EL)) {
								//some heteroatom replacements have weird valencues, so only verify valency on more normal fragments
								try{
									mol.checkValencies();
								}
								catch (StructureBuildingException e) {
									fail("The following token's SMILES produced a structure with invalid valency: " + smiles);
								}
							}
						}
						catch (Exception e) {
							e.printStackTrace();
							fail("The following SMILES were in error: " + smiles);
						}
						finally {
							if (mol != null) {
								try {
									fm.removeFragment(mol);
								} catch (StructureBuildingException e) {
									e.printStackTrace();
								}
							}
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
