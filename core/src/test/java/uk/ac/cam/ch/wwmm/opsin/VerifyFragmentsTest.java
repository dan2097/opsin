package uk.ac.cam.ch.wwmm.opsin;

import static junit.framework.Assert.*;

import java.util.ArrayList;
import java.util.List;

import nu.xom.Document;
import nu.xom.Element;
import nu.xom.Elements;

import org.junit.Test;
import static org.mockito.Mockito.mock;
import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

public class VerifyFragmentsTest {
	private static final String RESOURCE_LOCATION = "uk/ac/cam/ch/wwmm/opsin/resources/";
	private static final ResourceGetter resourceGetter = new ResourceGetter(RESOURCE_LOCATION);
	
	@Test
	public void verifySMILES() throws Exception {
		FragmentManager fm = new FragmentManager(new SMILESFragmentBuilder(), mock(CMLFragmentBuilder.class), new IDManager());
		Document tokenFileDoc = resourceGetter.getXMLDocument("index.xml");
		Elements tokenFiles = tokenFileDoc.getRootElement().getChildElements();
		for (int i = 0; i < tokenFiles.size(); i++) {
			Element rootElement = resourceGetter.getXMLDocument(tokenFiles.get(i).getValue()).getRootElement();
			List<Element> tokenLists =new ArrayList<Element>();
			if (rootElement.getLocalName().equals("tokenLists")){//support for xml files with one "tokenList" or multiple "tokenList" under a "tokenLists" element
				Elements children =rootElement.getChildElements();
				for (int j = 0; j <children.size(); j++) {
					tokenLists.add(children.get(j));
				}
			}
			else{
				tokenLists.add(rootElement);
			}
			for (Element tokenList : tokenLists) {
				Elements tokenElements = tokenList.getChildElements("token");
				for(int j=0;j<tokenElements.size();j++) {
					Element token = tokenElements.get(j);
					if (SMILES_VALTYPE_VAL.equals(token.getAttributeValue(VALTYPE_ATR))){
						Fragment mol =null;
						try{
							String smiles = token.getAttributeValue(VALUE_ATR);
							String type = token.getAttribute(TYPE_ATR) !=null ?  token.getAttributeValue(TYPE_ATR) : "";
							String subType = token.getAttribute(SUBTYPE_ATR) !=null ?  token.getAttributeValue(SUBTYPE_ATR) : "";
	
							String labels = token.getAttribute(LABELS_ATR) !=null ?  token.getAttributeValue(LABELS_ATR) : "";
							mol = fm.buildSMILES(smiles, type, subType, labels);
						}
						catch (Exception e) {
							e.printStackTrace();
						}
						assertNotNull("The following token's SMILES or labels were in error: " +token.toXML(), mol);
						try{
							fm.checkValencies();
						}
						catch (StructureBuildingException e) {
							fail("The following token's SMILES produced a structure with invalid valency: " +token.toXML());
						}
						fm.removeFragment(mol);
					}
				}
			}
		}
	}
	
	@Test
	public void verifyCML() throws Exception{
		FragmentManager fm = new FragmentManager(mock(SMILESFragmentBuilder.class), new CMLFragmentBuilder(resourceGetter), new IDManager());
		Document tokenFileDoc = resourceGetter.getXMLDocument("index.xml");
		Elements tokenFiles = tokenFileDoc.getRootElement().getChildElements();
		for (int i = 0; i < tokenFiles.size(); i++) {
			Element rootElement = resourceGetter.getXMLDocument(tokenFiles.get(i).getValue()).getRootElement();
			List<Element> tokenLists =new ArrayList<Element>();
			if (rootElement.getLocalName().equals("tokenLists")){//support for xml files with one "tokenList" or multiple "tokenList" under a "tokenLists" element
				Elements children =rootElement.getChildElements();
				for (int j = 0; j <children.size(); j++) {
					tokenLists.add(children.get(j));
				}
			}
			else{
				tokenLists.add(rootElement);
			}
			for (Element tokenList : tokenLists) {
				Elements tokenElements = tokenList.getChildElements("token");
				for(int j=0;j<tokenElements.size();j++) {
					Element token = tokenElements.get(j);
					if (DBKEY_VALTYPE_VAL.equals(token.getAttributeValue(VALTYPE_ATR))){
						Fragment mol =null;
						try{
							String idStr = token.getAttributeValue(VALUE_ATR);
							String type = token.getAttribute(TYPE_ATR) !=null ?  token.getAttributeValue(TYPE_ATR) : "";
							String subType = token.getAttribute(SUBTYPE_ATR) !=null ?  token.getAttributeValue(SUBTYPE_ATR) : "";
							mol = fm.buildCML(idStr, type, subType);
						}
						catch (Exception e) {
							e.printStackTrace();
						}
						assertNotNull("The following token's CML was in error: " +token.toXML(), mol);
						try{
							fm.checkValencies();
						}
						catch (StructureBuildingException e) {
							fail("The following token's CML produced a structure with invalid valency: " +token.toXML());
						}
						fm.removeFragment(mol);
					}
				}
			}
		}
	}
}
