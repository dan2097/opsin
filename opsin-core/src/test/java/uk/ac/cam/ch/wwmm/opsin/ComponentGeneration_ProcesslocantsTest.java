package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.Assert.*;
import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

import org.junit.Before;
import org.junit.Test;

public class ComponentGeneration_ProcesslocantsTest {
	
	private Element locant;
	private Element substituent;
	
	@Before
	public void setUpSubstituent(){
		substituent = new GroupingEl(SUBSTITUENT_EL);
		locant = new TokenEl(LOCANT_EL);
		substituent.addChild(locant);
		substituent.addChild(new TokenEl(GROUP_EL));//a dummy element to give the locant a potential purpose
	}
	
	@Test
	public void testCardinalNumber() throws ComponentGenerationException {
		locant.setValue("1");
		ComponentGenerator.processLocants(substituent);
		assertEquals("1", locant.getValue());
	}
	
	@Test
	public void testCardinalNumberWithHyphen() throws ComponentGenerationException {
		locant.setValue("1-");
		ComponentGenerator.processLocants(substituent);
		assertEquals("1", locant.getValue());
	}
	
	
	@Test
	public void testElementSymbol() throws ComponentGenerationException {
		locant.setValue("N-");
		ComponentGenerator.processLocants(substituent);
		assertEquals("N", locant.getValue());
	}
	
	@Test
	public void testAminoAcidStyleLocant() throws ComponentGenerationException {
		locant.setValue("N1-");
		ComponentGenerator.processLocants(substituent);
		assertEquals("N1", locant.getValue());
	}
	
	@Test
	public void testCompoundLocant() throws ComponentGenerationException {
		locant.setValue("1(10)-");
		ComponentGenerator.processLocants(substituent);
		assertEquals("1(10)", locant.getValue());
	}
	
	@Test
	public void testGreek() throws ComponentGenerationException {
		locant.setValue("alpha");
		ComponentGenerator.processLocants(substituent);
		assertEquals("alpha", locant.getValue());
	}
	
	@Test
	public void testNotlowercase1() throws ComponentGenerationException {
		locant.setValue("AlPhA-");
		ComponentGenerator.processLocants(substituent);
		assertEquals("alpha", locant.getValue());
	}
	
	@Test
	public void testNotlowercase2() throws ComponentGenerationException {
		locant.setValue("NAlPhA-");
		ComponentGenerator.processLocants(substituent);
		assertEquals("Nalpha", locant.getValue());
	}
	
	@Test
	public void testIUPAC2004() throws ComponentGenerationException {
		locant.setValue("2-N-");
		ComponentGenerator.processLocants(substituent);
		assertEquals("N2", locant.getValue());
	}
	
	@Test
	public void testSuperscript1() throws ComponentGenerationException {
		locant.setValue("N^(2)");
		ComponentGenerator.processLocants(substituent);
		assertEquals("N2", locant.getValue());
	}
	
	@Test
	public void testSuperscript2() throws ComponentGenerationException {
		locant.setValue("N^2");
		ComponentGenerator.processLocants(substituent);
		assertEquals("N2", locant.getValue());
	}
	
	@Test
	public void testSuperscript3() throws ComponentGenerationException {
		locant.setValue("N(2)");
		ComponentGenerator.processLocants(substituent);
		assertEquals("N2", locant.getValue());
	}
	
	@Test
	public void testSuperscript4() throws ComponentGenerationException {
		locant.setValue("N~12~");
		ComponentGenerator.processLocants(substituent);
		assertEquals("N12", locant.getValue());
	}
	
	@Test
	public void testSuperscript5() throws ComponentGenerationException {
		locant.setValue("N(alpha)");
		ComponentGenerator.processLocants(substituent);
		assertEquals("Nalpha", locant.getValue());
	}
	
	@Test
	public void testSuperscript6() throws ComponentGenerationException {
		locant.setValue("N^alpha");
		ComponentGenerator.processLocants(substituent);
		assertEquals("Nalpha", locant.getValue());
	}
	
	@Test
	public void testSuperscript7() throws ComponentGenerationException {
		locant.setValue("N*12*");
		ComponentGenerator.processLocants(substituent);
		assertEquals("N12", locant.getValue());
	}
	
	@Test
	public void testAddedHydrogen() throws ComponentGenerationException {
		locant.setValue("3(5'H)");
		ComponentGenerator.processLocants(substituent);
		assertEquals("3", locant.getValue());
		assertEquals(ADDEDHYDROGENLOCANT_TYPE_VAL, locant.getAttributeValue(TYPE_ATR));
		Element addedHydrogen = OpsinTools.getPreviousSibling(locant);
		assertNotNull(addedHydrogen);
		assertEquals(ADDEDHYDROGEN_EL, addedHydrogen.getName());
		assertEquals("5'", addedHydrogen.getAttributeValue(LOCANT_ATR));
	}
	
	@Test
	public void testAddedHydrogen2() throws ComponentGenerationException {
		locant.setValue("1,2(2H,7H)");
		ComponentGenerator.processLocants(substituent);
		assertEquals("1,2", locant.getValue());
		assertEquals(ADDEDHYDROGENLOCANT_TYPE_VAL, locant.getAttributeValue(TYPE_ATR));
		Element addedHydrogen1 = OpsinTools.getPreviousSibling(locant);
		assertNotNull(addedHydrogen1);
		assertEquals(ADDEDHYDROGEN_EL, addedHydrogen1.getName());
		assertEquals("7", addedHydrogen1.getAttributeValue(LOCANT_ATR));
		Element addedHydrogen2 = OpsinTools.getPreviousSibling(addedHydrogen1);
		assertNotNull(addedHydrogen2);
		assertEquals(ADDEDHYDROGEN_EL, addedHydrogen2.getName());
		assertEquals("2", addedHydrogen2.getAttributeValue(LOCANT_ATR));
	}

	@Test
	public void testStereochemistryInLocant1() throws ComponentGenerationException {
		locant.setValue("5(R)");
		ComponentGenerator.processLocants(substituent);
		assertEquals("5", locant.getValue());
		Element stereochemistry = OpsinTools.getPreviousSibling(locant);
		assertNotNull(stereochemistry);
		assertEquals(STEREOCHEMISTRY_EL, stereochemistry.getName());
		assertEquals(STEREOCHEMISTRYBRACKET_TYPE_VAL, stereochemistry.getAttributeValue(TYPE_ATR));
		assertEquals("(5R)", stereochemistry.getValue());//will be handled by process stereochemistry function
	}
	
	@Test
	public void testStereochemistryInLocant2() throws ComponentGenerationException {
		locant.setValue("5-(S)");
		ComponentGenerator.processLocants(substituent);
		assertEquals("5", locant.getValue());
		Element stereochemistry = OpsinTools.getPreviousSibling(locant);
		assertNotNull(stereochemistry);
		assertEquals(STEREOCHEMISTRY_EL, stereochemistry.getName());
		assertEquals(STEREOCHEMISTRYBRACKET_TYPE_VAL, stereochemistry.getAttributeValue(TYPE_ATR));
		assertEquals("(5S)", stereochemistry.getValue());//will be handled by process stereochemistry function
	}
	
	@Test
	public void testStereochemistryInLocant3() throws ComponentGenerationException {
		locant.setValue("N(3)-(S)");
		ComponentGenerator.processLocants(substituent);
		assertEquals("N3", locant.getValue());
		Element stereochemistry = OpsinTools.getPreviousSibling(locant);
		assertNotNull(stereochemistry);
		assertEquals(STEREOCHEMISTRY_EL, stereochemistry.getName());
		assertEquals(STEREOCHEMISTRYBRACKET_TYPE_VAL, stereochemistry.getAttributeValue(TYPE_ATR));
		assertEquals("(N3S)", stereochemistry.getValue());//will be handled by process stereochemistry function
	}
	
	@Test
	public void testStereochemistryInLocant4() throws ComponentGenerationException {
		locant.setValue("5(RS)");
		ComponentGenerator.processLocants(substituent);
		assertEquals("5", locant.getValue());
		Element stereochemistry = OpsinTools.getPreviousSibling(locant);
		assertNotNull(stereochemistry);
		assertEquals(STEREOCHEMISTRY_EL, stereochemistry.getName());
		assertEquals(STEREOCHEMISTRYBRACKET_TYPE_VAL, stereochemistry.getAttributeValue(TYPE_ATR));
		assertEquals("(5RS)", stereochemistry.getValue());//will be handled by process stereochemistry function
	}
	
	@Test
	public void testStereochemistryInLocant5() throws ComponentGenerationException {
		locant.setValue("5(R,S)");
		ComponentGenerator.processLocants(substituent);
		assertEquals("5", locant.getValue());
		Element stereochemistry = OpsinTools.getPreviousSibling(locant);
		assertNotNull(stereochemistry);
		assertEquals(STEREOCHEMISTRY_EL, stereochemistry.getName());
		assertEquals(STEREOCHEMISTRYBRACKET_TYPE_VAL, stereochemistry.getAttributeValue(TYPE_ATR));
		assertEquals("(5RS)", stereochemistry.getValue());//will be handled by process stereochemistry function
	}
	
	@Test
	public void testStereochemistryInLocant6() throws ComponentGenerationException {
		locant.setValue("5(R/S)");
		ComponentGenerator.processLocants(substituent);
		assertEquals("5", locant.getValue());
		Element stereochemistry = OpsinTools.getPreviousSibling(locant);
		assertNotNull(stereochemistry);
		assertEquals(STEREOCHEMISTRY_EL, stereochemistry.getName());
		assertEquals(STEREOCHEMISTRYBRACKET_TYPE_VAL, stereochemistry.getAttributeValue(TYPE_ATR));
		assertEquals("(5RS)", stereochemistry.getValue());//will be handled by process stereochemistry function
	}
	
	@Test
	public void testMultipleCardinals() throws ComponentGenerationException {
		locant.setValue("2,3-");
		ComponentGenerator.processLocants(substituent);
		assertEquals("2,3", locant.getValue());
	}
	
	@Test
	public void testMultipleTypesTogether() throws ComponentGenerationException {
		locant.setValue("2,N5,GaMMa,3-N,N^3,N(2),N~10~,4(5H),3-N(S),1(6)-");
		ComponentGenerator.processLocants(substituent);
		assertEquals("2,N5,gamma,N3,N3,N2,N10,4,N3,1(6)", locant.getValue());
		assertEquals(ADDEDHYDROGENLOCANT_TYPE_VAL, locant.getAttributeValue(TYPE_ATR));
		Element stereochemistry = OpsinTools.getPreviousSibling(locant);
		assertNotNull(stereochemistry);
		assertEquals(STEREOCHEMISTRY_EL, stereochemistry.getName());
		assertEquals(STEREOCHEMISTRYBRACKET_TYPE_VAL, stereochemistry.getAttributeValue(TYPE_ATR));
		assertEquals("(N3S)", stereochemistry.getValue());
		Element addedHydrogen = OpsinTools.getPreviousSibling(stereochemistry);
		assertNotNull(addedHydrogen);
		assertEquals(ADDEDHYDROGEN_EL, addedHydrogen.getName());
		assertEquals("5", addedHydrogen.getAttributeValue(LOCANT_ATR));
	}
	
	@Test
	public void testCarbohydrateStyleLocants() throws ComponentGenerationException {
		//2,4,6-tri-O
		locant.setValue("O");
		Element multiplier = new TokenEl(MULTIPLIER_EL);
		multiplier.addAttribute(new Attribute(VALUE_ATR, "3"));
		OpsinTools.insertBefore(locant, multiplier);
		Element numericLocant = new TokenEl(LOCANT_EL);
		numericLocant.setValue("2,4,6");
		OpsinTools.insertBefore(multiplier, numericLocant);
		ComponentGenerator.processLocants(substituent);
		assertEquals("O2,O4,O6", numericLocant.getValue());
		Element group = OpsinTools.getNextSibling(multiplier);
		assertNotNull(group);
		assertEquals(group.getName(), GROUP_EL);
	}
	
	@Test
	public void testCarbohydrateStyleLocantsNoNumericComponent() throws ComponentGenerationException {
		//tri-O
		locant.setValue("O");
		Element multiplier = new TokenEl(MULTIPLIER_EL);
		multiplier.addAttribute(new Attribute(VALUE_ATR, "3"));
		OpsinTools.insertBefore(locant, multiplier);
		ComponentGenerator.processLocants(substituent);
		Element elBeforeMultiplier = OpsinTools.getPreviousSibling(multiplier);
		assertNotNull("A locant should not be in front of the multiplier", elBeforeMultiplier);
		assertEquals(LOCANT_EL, elBeforeMultiplier.getName());
		assertEquals("O,O',O''", elBeforeMultiplier.getValue());
		Element group = OpsinTools.getNextSibling(multiplier);
		assertNotNull(group);
		assertEquals(group.getName(), GROUP_EL);
	}
	
	@Test
	public void testCarbohydrateStyleLocantsCounterExample() throws ComponentGenerationException {
		//2,4,6-tri-2 (this is not a carbohydrate style locant)
		locant.setValue("2");
		Element multiplier = new TokenEl(MULTIPLIER_EL);
		multiplier.addAttribute(new Attribute(VALUE_ATR, "3"));
		OpsinTools.insertBefore(locant, multiplier);
		Element numericLocant = new TokenEl(LOCANT_EL);
		numericLocant.setValue("2,4,6");
		OpsinTools.insertBefore(multiplier, numericLocant);
		ComponentGenerator.processLocants(substituent);
		assertEquals("2,4,6", numericLocant.getValue());
		Element unmodifiedLocant = OpsinTools.getNextSibling(multiplier);
		assertNotNull(unmodifiedLocant);
		assertEquals(unmodifiedLocant.getName(), LOCANT_EL);
		assertEquals("2", unmodifiedLocant.getValue());
	}
}
