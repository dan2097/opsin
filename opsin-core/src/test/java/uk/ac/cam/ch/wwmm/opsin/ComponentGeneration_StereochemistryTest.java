package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.Assert.*;
import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

import java.util.List;

import org.junit.Test;

public class ComponentGeneration_StereochemistryTest {

	@Test
	public void testUnlocantedS() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "(S)");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(1, children.size());
		Element newStereochemistryEl = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getName());
		assertEquals(null, newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("S", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
	}

	@Test
	public void testMultipleUnLocanted() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "(R,R)");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(2, children.size());
		Element newStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl1.getName());
		assertEquals(null, newStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("R", newStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl1.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl2 = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl2.getName());
		assertEquals(null, newStereochemistryEl2.getAttributeValue(LOCANT_ATR));
		assertEquals("R", newStereochemistryEl2.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl2.getAttributeValue(TYPE_ATR));
	}

	@Test
	public void testLocantedR() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "(1R)");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(1, children.size());
		Element newStereochemistryEl = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getName());
		assertEquals("1", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("R", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testMultipleRorSLocanted() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "(alphaR,3S,7'S)");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(3, children.size());
		Element newStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl1.getName());
		assertEquals("alpha", newStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("R", newStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl1.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl2 = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl2.getName());
		assertEquals("3", newStereochemistryEl2.getAttributeValue(LOCANT_ATR));
		assertEquals("S", newStereochemistryEl2.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl2.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl3 = children.get(2);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl3.getName());
		assertEquals("7'", newStereochemistryEl3.getAttributeValue(LOCANT_ATR));
		assertEquals("S", newStereochemistryEl3.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl3.getAttributeValue(TYPE_ATR));
	}
	
	
	@Test
	public void testUnLocantedE() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "(E)");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(1, children.size());
		Element newStereochemistryEl = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getName());
		assertEquals(null, newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("E", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(E_OR_Z_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testLocantedZ() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "(5Z)");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(1, children.size());
		Element newStereochemistryEl = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getName());
		assertEquals("5", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("Z", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(E_OR_Z_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testMultipleRorSorEorZ() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "(NZ,2E,R)");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(3, children.size());
		Element newStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl1.getName());
		assertEquals("N", newStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("Z", newStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(E_OR_Z_TYPE_VAL, newStereochemistryEl1.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl2 = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl2.getName());
		assertEquals("2", newStereochemistryEl2.getAttributeValue(LOCANT_ATR));
		assertEquals("E", newStereochemistryEl2.getAttributeValue(VALUE_ATR));
		assertEquals(E_OR_Z_TYPE_VAL, newStereochemistryEl2.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl3 = children.get(2);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl3.getName());
		assertEquals(null, newStereochemistryEl3.getAttributeValue(LOCANT_ATR));
		assertEquals("R", newStereochemistryEl3.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl3.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testDashInsteadOfComma() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "(NZ,2E-R)");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(3, children.size());
		Element newStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl1.getName());
		assertEquals("N", newStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("Z", newStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(E_OR_Z_TYPE_VAL, newStereochemistryEl1.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl2 = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl2.getName());
		assertEquals("2", newStereochemistryEl2.getAttributeValue(LOCANT_ATR));
		assertEquals("E", newStereochemistryEl2.getAttributeValue(VALUE_ATR));
		assertEquals(E_OR_Z_TYPE_VAL, newStereochemistryEl2.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl3 = children.get(2);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl3.getName());
		assertEquals(null, newStereochemistryEl3.getAttributeValue(LOCANT_ATR));
		assertEquals("R", newStereochemistryEl3.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl3.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testBracketedLocantedCisTrans() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "(3cis,5trans)");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(2, children.size());
		Element newStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl1.getName());
		assertEquals("3", newStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("cis", newStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(CISORTRANS_TYPE_VAL, newStereochemistryEl1.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl2 = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl2.getName());
		assertEquals("5", newStereochemistryEl2.getAttributeValue(LOCANT_ATR));
		assertEquals("trans", newStereochemistryEl2.getAttributeValue(VALUE_ATR));
		assertEquals(CISORTRANS_TYPE_VAL, newStereochemistryEl2.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testBracketedUnlocantedCisTrans() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "(5S-trans)");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(2, children.size());
		Element newStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl1.getName());
		assertEquals("5", newStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("S", newStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl1.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl2 = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl2.getName());
		assertEquals(null, newStereochemistryEl2.getAttributeValue(LOCANT_ATR));
		assertEquals("trans", newStereochemistryEl2.getAttributeValue(VALUE_ATR));
		assertEquals(CISORTRANS_TYPE_VAL, newStereochemistryEl2.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testBracketedExo() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "(exo)");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);
		
		List<Element> children = substituent.getChildElements();
		assertEquals(1, children.size());
		Element newStereochemistryEl = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getName());
		assertEquals(null, newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("exo", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ENDO_EXO_SYN_ANTI_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testBracketedEndo() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "(3-endo,5S)");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);
		
		List<Element> children = substituent.getChildElements();
		assertEquals(2, children.size());
		Element newStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl1.getName());
		assertEquals("3", newStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("endo", newStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(ENDO_EXO_SYN_ANTI_TYPE_VAL, newStereochemistryEl1.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl2 = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl2.getName());
		assertEquals("5", newStereochemistryEl2.getAttributeValue(LOCANT_ATR));
		assertEquals("S", newStereochemistryEl2.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl2.getAttributeValue(TYPE_ATR));
	}

	@Test
	public void testLocantedCisTrans() throws ComponentGenerationException {
		//XML for 3-cis,5-trans:
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element locant = new TokenEl(LOCANT_EL, "3");
		substituent.addChild(locant);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "cis");
		stereochem.addAttribute(new Attribute(TYPE_ATR, CISORTRANS_TYPE_VAL));
		stereochem.addAttribute(new Attribute(VALUE_ATR, "cis"));
		substituent.addChild(stereochem);
		locant = new TokenEl(LOCANT_EL, "5");
		substituent.addChild(locant);
		stereochem = new TokenEl(STEREOCHEMISTRY_EL, "trans");
		stereochem.addAttribute(new Attribute(TYPE_ATR, CISORTRANS_TYPE_VAL));
		stereochem.addAttribute(new Attribute(VALUE_ATR, "trans"));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(2, children.size());
		Element modifiedStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, modifiedStereochemistryEl1.getName());
		assertEquals("3", modifiedStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("cis", modifiedStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(CISORTRANS_TYPE_VAL, modifiedStereochemistryEl1.getAttributeValue(TYPE_ATR));
		
		Element modifiedStereochemistryEl2 = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, modifiedStereochemistryEl2.getName());
		assertEquals("5", modifiedStereochemistryEl2.getAttributeValue(LOCANT_ATR));
		assertEquals("trans", modifiedStereochemistryEl2.getAttributeValue(VALUE_ATR));
		assertEquals(CISORTRANS_TYPE_VAL, modifiedStereochemistryEl2.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testLocantedExoOn() throws ComponentGenerationException {
		//XML for 3-exobicyclo[2.2.2]oct:
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element locant = new TokenEl(LOCANT_EL, "3");
		substituent.addChild(locant);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "exo");
		stereochem.addAttribute(new Attribute(TYPE_ATR, ENDO_EXO_SYN_ANTI_TYPE_VAL));
		stereochem.addAttribute(new Attribute(VALUE_ATR, "exo"));
		substituent.addChild(stereochem);
		Element multiplier = new TokenEl(MULTIPLIER_EL);
		multiplier.addAttribute(new Attribute(TYPE_ATR, VONBAEYER_TYPE_VAL));
		substituent.addChild(multiplier);
		Element vonBaeyer = new TokenEl(VONBAEYER_EL);
		substituent.addChild(vonBaeyer);
		Element group = new TokenEl(GROUP_EL);
		group.addAttribute(new Attribute(TYPE_ATR, CHAIN_TYPE_VAL));
		group.addAttribute(new Attribute(SUBTYPE_ATR, ALKANESTEM_SUBTYPE_VAL));
		substituent.addChild(group);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(4, children.size());
		Element modifiedStereochemistryEl = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, modifiedStereochemistryEl.getName());
		assertEquals("3", modifiedStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("exo", modifiedStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ENDO_EXO_SYN_ANTI_TYPE_VAL, modifiedStereochemistryEl.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testLocantedExo() throws ComponentGenerationException {
		//XML for 3-exoamino
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element locant = new TokenEl(LOCANT_EL, "3");
		substituent.addChild(locant);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "exo");
		stereochem.addAttribute(new Attribute(TYPE_ATR, ENDO_EXO_SYN_ANTI_TYPE_VAL));
		stereochem.addAttribute(new Attribute(VALUE_ATR, "exo"));
		substituent.addChild(stereochem);
		Element group = new TokenEl(GROUP_EL);
		group.addAttribute(new Attribute(TYPE_ATR, SUBSTITUENT_EL));
		group.addAttribute(new Attribute(SUBTYPE_ATR, SIMPLESUBSTITUENT_SUBTYPE_VAL));
		substituent.addChild(group);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(3, children.size());
		assertEquals(LOCANT_EL, children.get(0).getName());
		Element modifiedStereochemistryEl = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, modifiedStereochemistryEl.getName());
		assertEquals("3", modifiedStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("exo", modifiedStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ENDO_EXO_SYN_ANTI_TYPE_VAL, modifiedStereochemistryEl.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testAnti() throws ComponentGenerationException {
		//XML for anti:
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "anti");
		stereochem.addAttribute(new Attribute(TYPE_ATR, ENDO_EXO_SYN_ANTI_TYPE_VAL));
		stereochem.addAttribute(new Attribute(VALUE_ATR, "anti"));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(1, children.size());
		Element unmodifiedStereochemistryEl = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, unmodifiedStereochemistryEl.getName());
		assertEquals("anti", unmodifiedStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ENDO_EXO_SYN_ANTI_TYPE_VAL, unmodifiedStereochemistryEl.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testCis() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "cis");
		stereochem.addAttribute(new Attribute(TYPE_ATR, CISORTRANS_TYPE_VAL));
		stereochem.addAttribute(new Attribute(VALUE_ATR, "cis"));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);
	
		List<Element> children = substituent.getChildElements();
		assertEquals(1, children.size());
		Element modifiedStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, modifiedStereochemistryEl1.getName());
		assertEquals(null, modifiedStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("cis", modifiedStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(CISORTRANS_TYPE_VAL, modifiedStereochemistryEl1.getAttributeValue(TYPE_ATR));
	}

	@Test
	public void testAxial1() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "(M)");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);
		
		List<Element> children = substituent.getChildElements();
		assertEquals(1, children.size());
		Element newStereochemistryEl = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getName());
		assertEquals(null, newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("M", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(AXIAL_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testAxial2() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "(Ra)");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);
		
		List<Element> children = substituent.getChildElements();
		assertEquals(1, children.size());
		Element newStereochemistryEl = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getName());
		assertEquals(null, newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("Ra", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(AXIAL_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
	}

	@Test
	public void testZUnbracketted() throws ComponentGenerationException {//note that IUPAC mandates brackets
		//XML for Z,Z:
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "Z");
		stereochem.addAttribute(new Attribute(TYPE_ATR, E_OR_Z_TYPE_VAL));
		substituent.addChild(stereochem);
		stereochem = new TokenEl(STEREOCHEMISTRY_EL, "Z");
		stereochem.addAttribute(new Attribute(TYPE_ATR, E_OR_Z_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(2, children.size());
		Element modifiedStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, modifiedStereochemistryEl1.getName());
		assertEquals(null, modifiedStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("Z", modifiedStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(E_OR_Z_TYPE_VAL, modifiedStereochemistryEl1.getAttributeValue(TYPE_ATR));
		
		Element modifiedStereochemistryEl2 = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, modifiedStereochemistryEl2.getName());
		assertEquals(null, modifiedStereochemistryEl2.getAttributeValue(LOCANT_ATR));
		assertEquals("Z", modifiedStereochemistryEl2.getAttributeValue(VALUE_ATR));
		assertEquals(E_OR_Z_TYPE_VAL, modifiedStereochemistryEl2.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testEandZUnbrackettedLocanted() throws ComponentGenerationException {//note that IUPAC mandates brackets
		//XML for 2E,4Z:
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element locant = new TokenEl(LOCANT_EL, "2");
		substituent.addChild(locant);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "E");
		stereochem.addAttribute(new Attribute(TYPE_ATR, E_OR_Z_TYPE_VAL));
		substituent.addChild(stereochem);
		locant = new TokenEl(LOCANT_EL, "4");
		substituent.addChild(locant);
		stereochem = new TokenEl(STEREOCHEMISTRY_EL, "Z");
		stereochem.addAttribute(new Attribute(TYPE_ATR, E_OR_Z_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(2, children.size());
		Element modifiedStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, modifiedStereochemistryEl1.getName());
		assertEquals("2", modifiedStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("E", modifiedStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(E_OR_Z_TYPE_VAL, modifiedStereochemistryEl1.getAttributeValue(TYPE_ATR));
		
		Element modifiedStereochemistryEl2 = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, modifiedStereochemistryEl2.getName());
		assertEquals("4", modifiedStereochemistryEl2.getAttributeValue(LOCANT_ATR));
		assertEquals("Z", modifiedStereochemistryEl2.getAttributeValue(VALUE_ATR));
		assertEquals(E_OR_Z_TYPE_VAL, modifiedStereochemistryEl2.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testEandZUnbrackettedBeforeEne() throws ComponentGenerationException {//not allowed in IUPAC names
		//XML for 2E,4Z-diene:
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element locant = new TokenEl(LOCANT_EL, "2");
		substituent.addChild(locant);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "E");
		stereochem.addAttribute(new Attribute(TYPE_ATR, E_OR_Z_TYPE_VAL));
		substituent.addChild(stereochem);
		locant = new TokenEl(LOCANT_EL, "4");
		substituent.addChild(locant);
		stereochem = new TokenEl(STEREOCHEMISTRY_EL, "Z");
		stereochem.addAttribute(new Attribute(TYPE_ATR, E_OR_Z_TYPE_VAL));
		substituent.addChild(stereochem);
		Element multiplier = new TokenEl(MULTIPLIER_EL, "di");
		multiplier.addAttribute(new Attribute(VALUE_ATR, "2"));
		substituent.addChild(multiplier);
		Element unsaturator = new TokenEl(UNSATURATOR_EL, "ene");
		unsaturator.addAttribute(new Attribute(VALUE_ATR, "2"));
		substituent.addChild(unsaturator);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(5, children.size());
		Element modifiedStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, modifiedStereochemistryEl1.getName());
		assertEquals("2", modifiedStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("E", modifiedStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(E_OR_Z_TYPE_VAL, modifiedStereochemistryEl1.getAttributeValue(TYPE_ATR));
		
		Element modifiedStereochemistryEl2 = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, modifiedStereochemistryEl2.getName());
		assertEquals("4", modifiedStereochemistryEl2.getAttributeValue(LOCANT_ATR));
		assertEquals("Z", modifiedStereochemistryEl2.getAttributeValue(VALUE_ATR));
		assertEquals(E_OR_Z_TYPE_VAL, modifiedStereochemistryEl2.getAttributeValue(TYPE_ATR));
		
		Element newLocant = children.get(2);
		assertEquals(LOCANT_EL, newLocant.getName());
		assertEquals("2,4", newLocant.getValue());
		assertEquals(MULTIPLIER_EL ,children.get(3).getName());
		assertEquals(UNSATURATOR_EL, children.get(4).getName());
	}
	
	@Test
	public void testEandZUnbrackettedBeforeYlidene() throws ComponentGenerationException {//not allowed in IUPAC names
		//XML for 2Z-ylidene:
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element locant = new TokenEl(LOCANT_EL, "2");
		substituent.addChild(locant);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "Z");
		stereochem.addAttribute(new Attribute(TYPE_ATR, E_OR_Z_TYPE_VAL));
		substituent.addChild(stereochem);
		Element suffix = new TokenEl(SUFFIX_EL, "ylidene");
		suffix.addAttribute(new Attribute(VALUE_ATR, "ylidene"));
		substituent.addChild(suffix);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(3, children.size());
		Element modifiedStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, modifiedStereochemistryEl1.getName());
		assertEquals("2", modifiedStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("Z", modifiedStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(E_OR_Z_TYPE_VAL, modifiedStereochemistryEl1.getAttributeValue(TYPE_ATR));

		Element newLocant = children.get(1);
		assertEquals(LOCANT_EL, newLocant.getName());
		assertEquals("2", newLocant.getValue());
		assertEquals(SUFFIX_EL ,children.get(2).getName());
	}

	@Test
	public void testBrackettedAlphaBeta() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "(1a,2b,3bEtA,4alpha,5xi)");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		Element naturalProduct = new TokenEl(GROUP_EL);
		naturalProduct.addAttribute(new Attribute(SUBTYPE_ATR, BIOCHEMICAL_SUBTYPE_VAL));
		naturalProduct.addAttribute(new Attribute(ALPHABETACLOCKWISEATOMORDERING_ATR, ""));
		substituent.addChild(naturalProduct);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(6, children.size());
		Element newStereochemistryEl = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getName());
		assertEquals("1", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("alpha", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ALPHA_OR_BETA_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
		
		newStereochemistryEl = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getName());
		assertEquals("2", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("beta", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ALPHA_OR_BETA_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
		
		newStereochemistryEl = children.get(2);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getName());
		assertEquals("3", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("beta", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ALPHA_OR_BETA_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
		
		newStereochemistryEl = children.get(3);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getName());
		assertEquals("4", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("alpha", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ALPHA_OR_BETA_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
		
		newStereochemistryEl = children.get(4);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getName());
		assertEquals("5", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("xi", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ALPHA_OR_BETA_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
	}
	
	
	@Test
	public void testAlphaBeta() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "3beta,5alpha");
		stereochem.addAttribute(new Attribute(TYPE_ATR, ALPHA_OR_BETA_TYPE_VAL));
		substituent.addChild(stereochem);
		Element naturalProduct = new TokenEl(GROUP_EL);
		naturalProduct.addAttribute(new Attribute(SUBTYPE_ATR, BIOCHEMICAL_SUBTYPE_VAL));
		naturalProduct.addAttribute(new Attribute(ALPHABETACLOCKWISEATOMORDERING_ATR, ""));
		substituent.addChild(naturalProduct);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(3, children.size());
		Element newStereochemistryEl = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getName());
		assertEquals("3", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("beta", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ALPHA_OR_BETA_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
		
		newStereochemistryEl = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getName());
		assertEquals("5", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("alpha", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ALPHA_OR_BETA_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testAlphaBetaNotDirectlyPrecedingANaturalProduct1() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "3beta,5alpha");
		stereochem.addAttribute(new Attribute(TYPE_ATR, ALPHA_OR_BETA_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(3, children.size());
		Element newStereochemistryEl = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getName());
		assertEquals("3", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("beta", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ALPHA_OR_BETA_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
		
		newStereochemistryEl = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getName());
		assertEquals("5", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("alpha", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ALPHA_OR_BETA_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
		
		Element newLocantEl = children.get(2);
		assertEquals(LOCANT_EL, newLocantEl.getName());
		assertEquals("3,5", newLocantEl.getValue());
	}
	
	@Test
	public void testAlphaBetaNotDirectlyPrecedingANaturalProduct2() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "(3beta,5alpha)");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(2, children.size());
		Element newStereochemistryEl = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getName());
		assertEquals("3", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("beta", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ALPHA_OR_BETA_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
		
		newStereochemistryEl = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getName());
		assertEquals("5", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("alpha", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ALPHA_OR_BETA_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testAlphaBetaNotDirectlyPrecedingANaturalProduct3() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element naturalProduct = new TokenEl(GROUP_EL);
		naturalProduct.addAttribute(new Attribute(SUBTYPE_ATR, BIOCHEMICAL_SUBTYPE_VAL));
		naturalProduct.addAttribute(new Attribute(ALPHABETACLOCKWISEATOMORDERING_ATR, ""));
		substituent.addChild(naturalProduct);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "3beta,5alpha");
		stereochem.addAttribute(new Attribute(TYPE_ATR, ALPHA_OR_BETA_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(4, children.size());
		Element newStereochemistryEl = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getName());
		assertEquals("3", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("beta", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ALPHA_OR_BETA_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
		
		newStereochemistryEl = children.get(2);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getName());
		assertEquals("5", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("alpha", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ALPHA_OR_BETA_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
		
		Element newLocantEl = children.get(3);
		assertEquals(LOCANT_EL, newLocantEl.getName());
		assertEquals("3,5", newLocantEl.getValue());
	}
	
	@Test
	public void testAlphaBetaStereoMixedWithNormalLocants() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "3beta,4,10,12alpha");
		stereochem.addAttribute(new Attribute(TYPE_ATR, ALPHA_OR_BETA_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(3, children.size());
		Element newStereochemistryEl = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getName());
		assertEquals("3", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("beta", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ALPHA_OR_BETA_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
		
		newStereochemistryEl = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getName());
		assertEquals("12", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("alpha", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ALPHA_OR_BETA_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
		
		Element newLocantEl = children.get(2);
		assertEquals(LOCANT_EL, newLocantEl.getName());
		assertEquals("3,4,10,12", newLocantEl.getValue());
	}
	
	//relative stereochemistry is currently treated the same as absolute stereochemistry
	@Test
	public void testRelativeStereoChemistry1() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "rel-(1R,3S,4S,7R)");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(4, children.size());
		Element newStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl1.getName());
		assertEquals("1", newStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("R", newStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl1.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl2 = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl2.getName());
		assertEquals("3", newStereochemistryEl2.getAttributeValue(LOCANT_ATR));
		assertEquals("S", newStereochemistryEl2.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl2.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl3 = children.get(2);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl3.getName());
		assertEquals("4", newStereochemistryEl3.getAttributeValue(LOCANT_ATR));
		assertEquals("S", newStereochemistryEl3.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl3.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl4 = children.get(3);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl4.getName());
		assertEquals("7", newStereochemistryEl4.getAttributeValue(LOCANT_ATR));
		assertEquals("R", newStereochemistryEl4.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl4.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testRelativeStereoChemistry2() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "(1R*,3S*,4S*,7R*)");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(4, children.size());
		Element newStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl1.getName());
		assertEquals("1", newStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("R", newStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl1.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl2 = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl2.getName());
		assertEquals("3", newStereochemistryEl2.getAttributeValue(LOCANT_ATR));
		assertEquals("S", newStereochemistryEl2.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl2.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl3 = children.get(2);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl3.getName());
		assertEquals("4", newStereochemistryEl3.getAttributeValue(LOCANT_ATR));
		assertEquals("S", newStereochemistryEl3.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl3.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl4 = children.get(3);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl4.getName());
		assertEquals("7", newStereochemistryEl4.getAttributeValue(LOCANT_ATR));
		assertEquals("R", newStereochemistryEl4.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl4.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testRelativeStereoChemistry3() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "rel-");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(0, children.size());
	}
	
	//relativeCisTrans is only supported sufficiently to get constitutionally correct results i.e. locants extracted from the stereochemistry
	@Test
	public void testRelativeCisTrans() throws ComponentGenerationException {
		//c-4-
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "c-4-");
		stereochem.addAttribute(new Attribute(TYPE_ATR, RELATIVECISTRANS_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);
	
		List<Element> children = substituent.getChildElements();
		assertEquals(2, children.size());
		Element modifiedStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, modifiedStereochemistryEl1.getName());
		assertEquals(null, modifiedStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("c-4-", modifiedStereochemistryEl1.getValue());
		assertEquals(RELATIVECISTRANS_TYPE_VAL, modifiedStereochemistryEl1.getAttributeValue(TYPE_ATR));
		Element locant = children.get(1);
		assertEquals(LOCANT_EL, locant.getName());
		assertEquals("4", locant.getValue());
	}
	
	@Test
	public void testRacemate1() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "rac-(2R)");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(1, children.size());
		Element newStereochemistryEl = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getName());
		assertEquals("2", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("RS", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testRacemate2() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "(RS)");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(1, children.size());
		Element newStereochemistryEl = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getName());
		assertEquals(null, newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("RS", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testRacemate2_ci() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "(rs)");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(1, children.size());
		Element newStereochemistryEl = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getName());
		assertEquals(null, newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("RS", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testRacemate3() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "(SR)");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(1, children.size());
		Element newStereochemistryEl = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getName());
		assertEquals(null, newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("SR", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
	}

	@Test
	public void testRacemate4() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "rac-(2R,4S)");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(2, children.size());
		Element newStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl1.getName());
		assertEquals("2", newStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("RS", newStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl1.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl2 = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl2.getName());
		assertEquals("4", newStereochemistryEl2.getAttributeValue(LOCANT_ATR));
		assertEquals("SR", newStereochemistryEl2.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl2.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testRacemate5() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "(2RS,4SR)");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(2, children.size());
		Element newStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl1.getName());
		assertEquals("2", newStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("RS", newStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl1.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl2 = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl2.getName());
		assertEquals("4", newStereochemistryEl2.getAttributeValue(LOCANT_ATR));
		assertEquals("SR", newStereochemistryEl2.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl2.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testRacemate6() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "rac-");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(0, children.size());
	}
	
	@Test
	public void testRacemate7() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "racem-");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(0, children.size());
	}

	@Test
	public void testRacemate8() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "racemic-");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(0, children.size());
	}
	
	@Test
	public void testRacemate9() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "(R/S)-");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(1, children.size());
		Element newStereochemistryEl = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getName());
		assertEquals(null, newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("RS", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testRacemate10() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "(RAC)");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(0, children.size());
	}

	@Test
	public void testRacemateEz1() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "(EZ)");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		Element modifiedStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, modifiedStereochemistryEl1.getName());
		assertEquals("EZ", modifiedStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(E_OR_Z_TYPE_VAL, modifiedStereochemistryEl1.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testRacemateEz2() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "(2EZ)");
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(1, children.size());
		Element modifiedStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, modifiedStereochemistryEl1.getName());
		assertEquals("2", modifiedStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("EZ", modifiedStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(E_OR_Z_TYPE_VAL, modifiedStereochemistryEl1.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testRacemateEz3_unbracketted() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element locant = new TokenEl(LOCANT_EL, "2");
		substituent.addChild(locant);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "ez");
		stereochem.addAttribute(new Attribute(TYPE_ATR, E_OR_Z_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(1, children.size());
		Element modifiedStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, modifiedStereochemistryEl1.getName());
		assertEquals("2", modifiedStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("EZ", modifiedStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(E_OR_Z_TYPE_VAL, modifiedStereochemistryEl1.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testRacemateEz4_unbracketted() throws ComponentGenerationException {
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element stereochem = new TokenEl(STEREOCHEMISTRY_EL, "EZ");
		stereochem.addAttribute(new Attribute(TYPE_ATR, E_OR_Z_TYPE_VAL));
		substituent.addChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		Element modifiedStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, modifiedStereochemistryEl1.getName());
		assertEquals("EZ", modifiedStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(E_OR_Z_TYPE_VAL, modifiedStereochemistryEl1.getAttributeValue(TYPE_ATR));
	}
	
	private void processStereochemistry(Element subOrRoot) throws ComponentGenerationException {
		new ComponentGenerator(new NameToStructureConfig()).processStereochemistry(subOrRoot);
	}
}
