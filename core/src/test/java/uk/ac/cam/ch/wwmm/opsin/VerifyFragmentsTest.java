package uk.ac.cam.ch.wwmm.opsin;

import static junit.framework.Assert.*;

import java.util.ArrayList;
import java.util.List;

import nu.xom.Document;
import nu.xom.Element;
import nu.xom.Elements;

import org.junit.BeforeClass;
import org.junit.Test;

public class VerifyFragmentsTest {
	private static final String RESOURCE_LOCATION = "uk/ac/cam/ch/wwmm/opsin/resources/";
	private static final ResourceGetter resourceGetter = new ResourceGetter(RESOURCE_LOCATION);
	private static SMILESFragmentBuilder smilesBuilder;
	private static CMLFragmentBuilder cmlBuilder;
	private static FragmentManager fm;
	
	@BeforeClass
	public static void setUp() throws Exception {
		smilesBuilder = new SMILESFragmentBuilder();
		cmlBuilder = new CMLFragmentBuilder(resourceGetter);
		fm = new FragmentManager(smilesBuilder, cmlBuilder, new IDManager());
	}
	
	@Test
	public void verifySMILES() throws Exception {
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
					if (token.getAttribute("valType")!=null && token.getAttributeValue("valType").equals("SMILES")){
						Fragment mol =null;
						try{
							String smiles = token.getAttributeValue("value");
							String type = token.getAttribute("type") !=null ?  token.getAttributeValue("type") : "";
							String subType = token.getAttribute("subType") !=null ?  token.getAttributeValue("subType") : "";
	
							String labels = token.getAttribute("labels") !=null ?  token.getAttributeValue("labels") : "";
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
	public void verifyCML() throws Exception {
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
					if (token.getAttribute("valType")!=null && token.getAttributeValue("valType").equals("dbkey")){
						Fragment mol =null;
						try{
							String idStr = token.getAttributeValue("value");
							String type = token.getAttribute("type") !=null ?  token.getAttributeValue("type") : "";
							String subType = token.getAttribute("subType") !=null ?  token.getAttributeValue("subType") : "";
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
