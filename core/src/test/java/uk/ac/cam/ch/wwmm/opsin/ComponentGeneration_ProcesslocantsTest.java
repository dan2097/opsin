package uk.ac.cam.ch.wwmm.opsin;

import static junit.framework.Assert.*;
import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

import nu.xom.Element;

import org.junit.Before;
import org.junit.Test;

public class ComponentGeneration_ProcesslocantsTest {
	
	private Element locant;
	private Element substituent;
	
	@Before
	public void setUpSubstituent(){
		substituent = new Element(SUBSTITUENT_EL);
		locant = new Element(LOCANT_EL);
		substituent.appendChild(locant);
		substituent.appendChild(new Element(GROUP_EL));//a dummy element to give the locant a potential purpose
	}
	
	@Test
	public void testCardinalNumber() throws ComponentGenerationException {
		locant.appendChild("1");
		ComponentGenerator.processLocants(substituent);
		assertEquals("1", locant.getValue());
	}
	
	@Test
	public void testCardinalNumberWithHyphen() throws ComponentGenerationException {
		locant.appendChild("1-");
		ComponentGenerator.processLocants(substituent);
		assertEquals("1", locant.getValue());
	}
	
	
	@Test
	public void testElementSymbol() throws ComponentGenerationException {
		locant.appendChild("N-");
		ComponentGenerator.processLocants(substituent);
		assertEquals("N", locant.getValue());
	}
	
	@Test
	public void testAminoAcidStyleLocant() throws ComponentGenerationException {
		locant.appendChild("N1-");
		ComponentGenerator.processLocants(substituent);
		assertEquals("N1", locant.getValue());
	}
	
	@Test
	public void testCompoundLocant() throws ComponentGenerationException {
		locant.appendChild("1(10)-");
		ComponentGenerator.processLocants(substituent);
		assertEquals("1(10)", locant.getValue());
	}
	
	@Test
	public void testGreek() throws ComponentGenerationException {
		locant.appendChild("alpha");
		ComponentGenerator.processLocants(substituent);
		assertEquals("alpha", locant.getValue());
	}
	
	@Test
	public void testNotlowercase1() throws ComponentGenerationException {
		locant.appendChild("AlPhA-");
		ComponentGenerator.processLocants(substituent);
		assertEquals("alpha", locant.getValue());
	}
	
	@Test
	public void testNotlowercase2() throws ComponentGenerationException {
		locant.appendChild("NAlPhA-");
		ComponentGenerator.processLocants(substituent);
		assertEquals("Nalpha", locant.getValue());
	}
	
	@Test
	public void testIUPAC2004() throws ComponentGenerationException {
		locant.appendChild("2-N-");
		ComponentGenerator.processLocants(substituent);
		assertEquals("N2", locant.getValue());
	}
	
	@Test
	public void testSuperscript1() throws ComponentGenerationException {
		locant.appendChild("N^(2)");
		ComponentGenerator.processLocants(substituent);
		assertEquals("N2", locant.getValue());
	}
	
	@Test
	public void testSuperscript2() throws ComponentGenerationException {
		locant.appendChild("N^2");
		ComponentGenerator.processLocants(substituent);
		assertEquals("N2", locant.getValue());
	}
	
	@Test
	public void testSuperscript3() throws ComponentGenerationException {
		locant.appendChild("N(2)");
		ComponentGenerator.processLocants(substituent);
		assertEquals("N2", locant.getValue());
	}
	
	@Test
	public void testSuperscript4() throws ComponentGenerationException {
		locant.appendChild("N~12~");
		ComponentGenerator.processLocants(substituent);
		assertEquals("N12", locant.getValue());
	}
	
	@Test
	public void testIndicatedHydrogen() throws ComponentGenerationException {
		locant.appendChild("3(5'H)");
		ComponentGenerator.processLocants(substituent);
		assertEquals("3", locant.getValue());
		assertEquals(ADDEDHYDROGENLOCANT_TYPE_VAL, locant.getAttributeValue(TYPE_ATR));
		Element indicatedHydrogen = (Element) XOMTools.getPreviousSibling(locant);
		assertNotNull(indicatedHydrogen);
		assertEquals(INDICATEDHYDROGEN_EL, indicatedHydrogen.getLocalName());
		assertEquals("5'", indicatedHydrogen.getAttributeValue(LOCANT_ATR));
	}
	
	@Test
	public void testIndicatedHydrogen2() throws ComponentGenerationException {
		locant.appendChild("1,2(2H,7H)");
		ComponentGenerator.processLocants(substituent);
		assertEquals("1,2", locant.getValue());
		assertEquals(ADDEDHYDROGENLOCANT_TYPE_VAL, locant.getAttributeValue(TYPE_ATR));
		Element indicatedHydrogen1 = (Element) XOMTools.getPreviousSibling(locant);
		assertNotNull(indicatedHydrogen1);
		assertEquals(INDICATEDHYDROGEN_EL, indicatedHydrogen1.getLocalName());
		assertEquals("7", indicatedHydrogen1.getAttributeValue(LOCANT_ATR));
		Element indicatedHydrogen2 = (Element) XOMTools.getPreviousSibling(indicatedHydrogen1);
		assertNotNull(indicatedHydrogen2);
		assertEquals(INDICATEDHYDROGEN_EL, indicatedHydrogen2.getLocalName());
		assertEquals("2", indicatedHydrogen2.getAttributeValue(LOCANT_ATR));
	}

	@Test
	public void testStereochemistryInLocant1() throws ComponentGenerationException {
		locant.appendChild("5(R)");
		ComponentGenerator.processLocants(substituent);
		assertEquals("5", locant.getValue());
		Element stereochemistry = (Element) XOMTools.getPreviousSibling(locant);
		assertNotNull(stereochemistry);
		assertEquals(STEREOCHEMISTRY_EL, stereochemistry.getLocalName());
		assertEquals(STEREOCHEMISTRYBRACKET_TYPE_VAL, stereochemistry.getAttributeValue(TYPE_ATR));
		assertEquals("(5R)", stereochemistry.getValue());//will be handled by process stereochemistry function
	}
	
	@Test
	public void testStereochemistryInLocant2() throws ComponentGenerationException {
		locant.appendChild("5-(S)");
		ComponentGenerator.processLocants(substituent);
		assertEquals("5", locant.getValue());
		Element stereochemistry = (Element) XOMTools.getPreviousSibling(locant);
		assertNotNull(stereochemistry);
		assertEquals(STEREOCHEMISTRY_EL, stereochemistry.getLocalName());
		assertEquals(STEREOCHEMISTRYBRACKET_TYPE_VAL, stereochemistry.getAttributeValue(TYPE_ATR));
		assertEquals("(5S)", stereochemistry.getValue());//will be handled by process stereochemistry function
	}
	
	@Test
	public void testStereochemistryInLocant3() throws ComponentGenerationException {
		locant.appendChild("N(3)-(S)");
		ComponentGenerator.processLocants(substituent);
		assertEquals("N3", locant.getValue());
		Element stereochemistry = (Element) XOMTools.getPreviousSibling(locant);
		assertNotNull(stereochemistry);
		assertEquals(STEREOCHEMISTRY_EL, stereochemistry.getLocalName());
		assertEquals(STEREOCHEMISTRYBRACKET_TYPE_VAL, stereochemistry.getAttributeValue(TYPE_ATR));
		assertEquals("(N3S)", stereochemistry.getValue());//will be handled by process stereochemistry function
	}
	
	@Test
	public void testMultipleCardinals() throws ComponentGenerationException {
		locant.appendChild("2,3-");
		ComponentGenerator.processLocants(substituent);
		assertEquals("2,3", locant.getValue());
	}
	
	@Test
	public void testMultipleTypesTogether() throws ComponentGenerationException {
		locant.appendChild("2,N5,GaMMa,3-N,N^3,N(2),N~10~,4(5H),3-N(S),1(6)-");
		ComponentGenerator.processLocants(substituent);
		assertEquals("2,N5,gamma,N3,N3,N2,N10,4,N3,1(6)", locant.getValue());
		assertEquals(ADDEDHYDROGENLOCANT_TYPE_VAL, locant.getAttributeValue(TYPE_ATR));
		Element stereochemistry = (Element) XOMTools.getPreviousSibling(locant);
		assertNotNull(stereochemistry);
		assertEquals(STEREOCHEMISTRY_EL, stereochemistry.getLocalName());
		assertEquals(STEREOCHEMISTRYBRACKET_TYPE_VAL, stereochemistry.getAttributeValue(TYPE_ATR));
		assertEquals("(N3S)", stereochemistry.getValue());
		Element indicatedHydrogen = (Element) XOMTools.getPreviousSibling(stereochemistry);
		assertNotNull(indicatedHydrogen);
		assertEquals(INDICATEDHYDROGEN_EL, indicatedHydrogen.getLocalName());
		assertEquals("5", indicatedHydrogen.getAttributeValue(LOCANT_ATR));
	}
}
