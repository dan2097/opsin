package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.Assert.*;
import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

import java.util.List;

import org.junit.Test;

public class ComponentGeneration_StereochemistryTest {

	@Test
	public void testUnlocantedS() throws ComponentGenerationException {
		Element substituent = new Element(SUBSTITUENT_EL);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.appendChild(stereochem);
		stereochem.setValue("(S)");
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(1, children.size());
		Element newStereochemistryEl = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getLocalName());
		assertEquals(null, newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("S", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
	}

	@Test
	public void testMultipleUnLocanted() throws ComponentGenerationException {
		Element substituent = new Element(SUBSTITUENT_EL);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.appendChild(stereochem);
		stereochem.setValue("(R,R)");
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(2, children.size());
		Element newStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl1.getLocalName());
		assertEquals(null, newStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("R", newStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl1.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl2 = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl2.getLocalName());
		assertEquals(null, newStereochemistryEl2.getAttributeValue(LOCANT_ATR));
		assertEquals("R", newStereochemistryEl2.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl2.getAttributeValue(TYPE_ATR));
	}

	@Test
	public void testLocantedR() throws ComponentGenerationException {
		Element substituent = new Element(SUBSTITUENT_EL);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.appendChild(stereochem);
		stereochem.setValue("(1R)");
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(1, children.size());
		Element newStereochemistryEl = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getLocalName());
		assertEquals("1", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("R", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testMultipleRorSLocanted() throws ComponentGenerationException {
		Element substituent = new Element(SUBSTITUENT_EL);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.appendChild(stereochem);
		stereochem.setValue("(alphaR,3S,7'S)");
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(3, children.size());
		Element newStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl1.getLocalName());
		assertEquals("alpha", newStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("R", newStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl1.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl2 = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl2.getLocalName());
		assertEquals("3", newStereochemistryEl2.getAttributeValue(LOCANT_ATR));
		assertEquals("S", newStereochemistryEl2.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl2.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl3 = children.get(2);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl3.getLocalName());
		assertEquals("7'", newStereochemistryEl3.getAttributeValue(LOCANT_ATR));
		assertEquals("S", newStereochemistryEl3.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl3.getAttributeValue(TYPE_ATR));
	}
	
	
	@Test
	public void testUnLocantedE() throws ComponentGenerationException {
		Element substituent = new Element(SUBSTITUENT_EL);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.appendChild(stereochem);
		stereochem.setValue("(E)");
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(1, children.size());
		Element newStereochemistryEl = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getLocalName());
		assertEquals(null, newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("E", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(E_OR_Z_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testLocantedZ() throws ComponentGenerationException {
		Element substituent = new Element(SUBSTITUENT_EL);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.appendChild(stereochem);
		stereochem.setValue("(5Z)");
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(1, children.size());
		Element newStereochemistryEl = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getLocalName());
		assertEquals("5", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("Z", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(E_OR_Z_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testMultipleRorSorEorZ() throws ComponentGenerationException {
		Element substituent = new Element(SUBSTITUENT_EL);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.appendChild(stereochem);
		stereochem.setValue("(NZ,2E,R)");
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(3, children.size());
		Element newStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl1.getLocalName());
		assertEquals("N", newStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("Z", newStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(E_OR_Z_TYPE_VAL, newStereochemistryEl1.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl2 = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl2.getLocalName());
		assertEquals("2", newStereochemistryEl2.getAttributeValue(LOCANT_ATR));
		assertEquals("E", newStereochemistryEl2.getAttributeValue(VALUE_ATR));
		assertEquals(E_OR_Z_TYPE_VAL, newStereochemistryEl2.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl3 = children.get(2);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl3.getLocalName());
		assertEquals(null, newStereochemistryEl3.getAttributeValue(LOCANT_ATR));
		assertEquals("R", newStereochemistryEl3.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl3.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testDashInsteadOfComma() throws ComponentGenerationException {
		Element substituent = new Element(SUBSTITUENT_EL);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.appendChild(stereochem);
		stereochem.setValue("(NZ,2E-R)");
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(3, children.size());
		Element newStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl1.getLocalName());
		assertEquals("N", newStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("Z", newStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(E_OR_Z_TYPE_VAL, newStereochemistryEl1.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl2 = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl2.getLocalName());
		assertEquals("2", newStereochemistryEl2.getAttributeValue(LOCANT_ATR));
		assertEquals("E", newStereochemistryEl2.getAttributeValue(VALUE_ATR));
		assertEquals(E_OR_Z_TYPE_VAL, newStereochemistryEl2.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl3 = children.get(2);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl3.getLocalName());
		assertEquals(null, newStereochemistryEl3.getAttributeValue(LOCANT_ATR));
		assertEquals("R", newStereochemistryEl3.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl3.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testBracketedLocantedCisTrans() throws ComponentGenerationException {
		Element substituent = new Element(SUBSTITUENT_EL);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.appendChild(stereochem);
		stereochem.setValue("(3cis,5trans)");
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(2, children.size());
		Element newStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl1.getLocalName());
		assertEquals("3", newStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("cis", newStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(CISORTRANS_TYPE_VAL, newStereochemistryEl1.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl2 = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl2.getLocalName());
		assertEquals("5", newStereochemistryEl2.getAttributeValue(LOCANT_ATR));
		assertEquals("trans", newStereochemistryEl2.getAttributeValue(VALUE_ATR));
		assertEquals(CISORTRANS_TYPE_VAL, newStereochemistryEl2.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testBracketedUnlocantedCisTrans() throws ComponentGenerationException {
		Element substituent = new Element(SUBSTITUENT_EL);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.appendChild(stereochem);
		stereochem.setValue("(5S-trans)");
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(2, children.size());
		Element newStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl1.getLocalName());
		assertEquals("5", newStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("S", newStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl1.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl2 = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl2.getLocalName());
		assertEquals(null, newStereochemistryEl2.getAttributeValue(LOCANT_ATR));
		assertEquals("trans", newStereochemistryEl2.getAttributeValue(VALUE_ATR));
		assertEquals(CISORTRANS_TYPE_VAL, newStereochemistryEl2.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testBracketedExo() throws ComponentGenerationException {
		Element substituent = new Element(SUBSTITUENT_EL);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.appendChild(stereochem);
		stereochem.setValue("(exo)");
		processStereochemistry(substituent);
		
		List<Element> children = substituent.getChildElements();
		assertEquals(1, children.size());
		Element newStereochemistryEl = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getLocalName());
		assertEquals(null, newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("exo", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ENDO_EXO_SYN_ANTI_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testBracketedEndo() throws ComponentGenerationException {
		Element substituent = new Element(SUBSTITUENT_EL);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.appendChild(stereochem);
		stereochem.setValue("(3-endo,5S)");
		processStereochemistry(substituent);
		
		List<Element> children = substituent.getChildElements();
		assertEquals(2, children.size());
		Element newStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl1.getLocalName());
		assertEquals("3", newStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("endo", newStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(ENDO_EXO_SYN_ANTI_TYPE_VAL, newStereochemistryEl1.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl2 = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl2.getLocalName());
		assertEquals("5", newStereochemistryEl2.getAttributeValue(LOCANT_ATR));
		assertEquals("S", newStereochemistryEl2.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl2.getAttributeValue(TYPE_ATR));
	}

	@Test
	public void testLocantedCisTrans() throws ComponentGenerationException {
		//XML for 3-cis,5-trans:
		Element substituent = new Element(SUBSTITUENT_EL);
		Element locant = new Element(LOCANT_EL);
		locant.setValue("3");
		substituent.appendChild(locant);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, CISORTRANS_TYPE_VAL));
		stereochem.addAttribute(new Attribute(VALUE_ATR, "cis"));
		stereochem.setValue("cis");
		substituent.appendChild(stereochem);
		locant = new Element(LOCANT_EL);
		locant.setValue("5");
		substituent.appendChild(locant);
		stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, CISORTRANS_TYPE_VAL));
		stereochem.addAttribute(new Attribute(VALUE_ATR, "trans"));
		stereochem.setValue("trans");
		substituent.appendChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(2, children.size());
		Element modifiedStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, modifiedStereochemistryEl1.getLocalName());
		assertEquals("3", modifiedStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("cis", modifiedStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(CISORTRANS_TYPE_VAL, modifiedStereochemistryEl1.getAttributeValue(TYPE_ATR));
		
		Element modifiedStereochemistryEl2 = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, modifiedStereochemistryEl2.getLocalName());
		assertEquals("5", modifiedStereochemistryEl2.getAttributeValue(LOCANT_ATR));
		assertEquals("trans", modifiedStereochemistryEl2.getAttributeValue(VALUE_ATR));
		assertEquals(CISORTRANS_TYPE_VAL, modifiedStereochemistryEl2.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testLocantedExoOn() throws ComponentGenerationException {
		//XML for 3-exobicyclo[2.2.2]oct:
		Element substituent = new Element(SUBSTITUENT_EL);
		Element locant = new Element(LOCANT_EL);
		locant.setValue("3");
		substituent.appendChild(locant);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, ENDO_EXO_SYN_ANTI_TYPE_VAL));
		stereochem.addAttribute(new Attribute(VALUE_ATR, "exo"));
		stereochem.setValue("exo");
		substituent.appendChild(stereochem);
		Element multiplier = new Element(MULTIPLIER_EL);
		multiplier.addAttribute(new Attribute(TYPE_ATR, VONBAEYER_TYPE_VAL));
		substituent.appendChild(multiplier);
		Element vonBaeyer = new Element(VONBAEYER_EL);
		substituent.appendChild(vonBaeyer);
		Element group = new Element(GROUP_EL);
		group.addAttribute(new Attribute(TYPE_ATR, CHAIN_TYPE_VAL));
		group.addAttribute(new Attribute(SUBTYPE_ATR, ALKANESTEM_SUBTYPE_VAL));
		substituent.appendChild(group);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(4, children.size());
		Element modifiedStereochemistryEl = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, modifiedStereochemistryEl.getLocalName());
		assertEquals("3", modifiedStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("exo", modifiedStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ENDO_EXO_SYN_ANTI_TYPE_VAL, modifiedStereochemistryEl.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testLocantedExo() throws ComponentGenerationException {
		//XML for 3-exoamino
		Element substituent = new Element(SUBSTITUENT_EL);
		Element locant = new Element(LOCANT_EL);
		locant.setValue("3");
		substituent.appendChild(locant);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, ENDO_EXO_SYN_ANTI_TYPE_VAL));
		stereochem.addAttribute(new Attribute(VALUE_ATR, "exo"));
		stereochem.setValue("exo");
		substituent.appendChild(stereochem);
		Element group = new Element(GROUP_EL);
		group.addAttribute(new Attribute(TYPE_ATR, SUBSTITUENT_EL));
		group.addAttribute(new Attribute(SUBTYPE_ATR, SIMPLESUBSTITUENT_SUBTYPE_VAL));
		substituent.appendChild(group);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(3, children.size());
		assertEquals(LOCANT_EL, children.get(0).getLocalName());
		Element modifiedStereochemistryEl = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, modifiedStereochemistryEl.getLocalName());
		assertEquals("3", modifiedStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("exo", modifiedStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ENDO_EXO_SYN_ANTI_TYPE_VAL, modifiedStereochemistryEl.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testAnti() throws ComponentGenerationException {
		//XML for anti:
		Element substituent = new Element(SUBSTITUENT_EL);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, ENDO_EXO_SYN_ANTI_TYPE_VAL));
		stereochem.addAttribute(new Attribute(VALUE_ATR, "anti"));
		stereochem.setValue("anti");
		substituent.appendChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(1, children.size());
		Element unmodifiedStereochemistryEl = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, unmodifiedStereochemistryEl.getLocalName());
		assertEquals("anti", unmodifiedStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ENDO_EXO_SYN_ANTI_TYPE_VAL, unmodifiedStereochemistryEl.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testCis() throws ComponentGenerationException {
		Element substituent = new Element(SUBSTITUENT_EL);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, CISORTRANS_TYPE_VAL));
		stereochem.addAttribute(new Attribute(VALUE_ATR, "cis"));
		stereochem.setValue("cis");
		substituent.appendChild(stereochem);
		processStereochemistry(substituent);
	
		List<Element> children = substituent.getChildElements();
		assertEquals(1, children.size());
		Element modifiedStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, modifiedStereochemistryEl1.getLocalName());
		assertEquals(null, modifiedStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("cis", modifiedStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(CISORTRANS_TYPE_VAL, modifiedStereochemistryEl1.getAttributeValue(TYPE_ATR));
	}

	@Test
	public void testZUnbracketted() throws ComponentGenerationException {//note that IUPAC mandates brackets
		//XML for Z,Z:
		Element substituent = new Element(SUBSTITUENT_EL);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, E_OR_Z_TYPE_VAL));
		stereochem.setValue("Z");
		substituent.appendChild(stereochem);
		stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, E_OR_Z_TYPE_VAL));
		stereochem.setValue("Z");
		substituent.appendChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(2, children.size());
		Element modifiedStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, modifiedStereochemistryEl1.getLocalName());
		assertEquals(null, modifiedStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("Z", modifiedStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(E_OR_Z_TYPE_VAL, modifiedStereochemistryEl1.getAttributeValue(TYPE_ATR));
		
		Element modifiedStereochemistryEl2 = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, modifiedStereochemistryEl2.getLocalName());
		assertEquals(null, modifiedStereochemistryEl2.getAttributeValue(LOCANT_ATR));
		assertEquals("Z", modifiedStereochemistryEl2.getAttributeValue(VALUE_ATR));
		assertEquals(E_OR_Z_TYPE_VAL, modifiedStereochemistryEl2.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testEandZUnbrackettedLocanted() throws ComponentGenerationException {//note that IUPAC mandates brackets
		//XML for 2E,4Z:
		Element substituent = new Element(SUBSTITUENT_EL);
		Element locant = new Element(LOCANT_EL);
		locant.setValue("2");
		substituent.appendChild(locant);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, E_OR_Z_TYPE_VAL));
		stereochem.setValue("E");
		substituent.appendChild(stereochem);
		locant = new Element(LOCANT_EL);
		locant.setValue("4");
		substituent.appendChild(locant);
		stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, E_OR_Z_TYPE_VAL));
		stereochem.setValue("Z");
		substituent.appendChild(stereochem);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(2, children.size());
		Element modifiedStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, modifiedStereochemistryEl1.getLocalName());
		assertEquals("2", modifiedStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("E", modifiedStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(E_OR_Z_TYPE_VAL, modifiedStereochemistryEl1.getAttributeValue(TYPE_ATR));
		
		Element modifiedStereochemistryEl2 = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, modifiedStereochemistryEl2.getLocalName());
		assertEquals("4", modifiedStereochemistryEl2.getAttributeValue(LOCANT_ATR));
		assertEquals("Z", modifiedStereochemistryEl2.getAttributeValue(VALUE_ATR));
		assertEquals(E_OR_Z_TYPE_VAL, modifiedStereochemistryEl2.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testBrackettedAlphaBeta() throws ComponentGenerationException {
		Element substituent = new Element(SUBSTITUENT_EL);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.appendChild(stereochem);
		stereochem.setValue("(1a,2b,3bEtA,4alpha,5xi)");
		Element naturalProduct = new Element(GROUP_EL);
		naturalProduct.addAttribute(new Attribute(SUBTYPE_ATR, BIOCHEMICAL_SUBTYPE_VAL));
		naturalProduct.addAttribute(new Attribute(ALPHABETACLOCKWISEATOMORDERING_ATR, ""));
		substituent.appendChild(naturalProduct);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(6, children.size());
		Element newStereochemistryEl = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getLocalName());
		assertEquals("1", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("alpha", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ALPHA_OR_BETA_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
		
		newStereochemistryEl = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getLocalName());
		assertEquals("2", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("beta", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ALPHA_OR_BETA_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
		
		newStereochemistryEl = children.get(2);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getLocalName());
		assertEquals("3", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("beta", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ALPHA_OR_BETA_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
		
		newStereochemistryEl = children.get(3);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getLocalName());
		assertEquals("4", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("alpha", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ALPHA_OR_BETA_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
		
		newStereochemistryEl = children.get(4);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getLocalName());
		assertEquals("5", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("xi", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ALPHA_OR_BETA_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
	}
	
	
	@Test
	public void testAlphaBeta() throws ComponentGenerationException {
		Element substituent = new Element(SUBSTITUENT_EL);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, ALPHA_OR_BETA_TYPE_VAL));
		substituent.appendChild(stereochem);
		stereochem.setValue("3beta,5alpha");
		Element naturalProduct = new Element(GROUP_EL);
		naturalProduct.addAttribute(new Attribute(SUBTYPE_ATR, BIOCHEMICAL_SUBTYPE_VAL));
		naturalProduct.addAttribute(new Attribute(ALPHABETACLOCKWISEATOMORDERING_ATR, ""));
		substituent.appendChild(naturalProduct);
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(3, children.size());
		Element newStereochemistryEl = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getLocalName());
		assertEquals("3", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("beta", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ALPHA_OR_BETA_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
		
		newStereochemistryEl = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getLocalName());
		assertEquals("5", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("alpha", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ALPHA_OR_BETA_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testAlphaBetaNotDirectlyPrecedingANaturalProduct1() throws ComponentGenerationException {
		Element substituent = new Element(SUBSTITUENT_EL);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, ALPHA_OR_BETA_TYPE_VAL));
		substituent.appendChild(stereochem);
		stereochem.setValue("3beta,5alpha");
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(3, children.size());
		Element newStereochemistryEl = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getLocalName());
		assertEquals("3", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("beta", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ALPHA_OR_BETA_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
		
		newStereochemistryEl = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getLocalName());
		assertEquals("5", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("alpha", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ALPHA_OR_BETA_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
		
		Element newLocantEl = children.get(2);
		assertEquals(LOCANT_EL, newLocantEl.getLocalName());
		assertEquals("3,5", newLocantEl.getValue());
	}
	
	@Test
	public void testAlphaBetaNotDirectlyPrecedingANaturalProduct2() throws ComponentGenerationException {
		Element substituent = new Element(SUBSTITUENT_EL);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.appendChild(stereochem);
		stereochem.setValue("(3beta,5alpha)");
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(2, children.size());
		Element newStereochemistryEl = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getLocalName());
		assertEquals("3", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("beta", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ALPHA_OR_BETA_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
		
		newStereochemistryEl = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getLocalName());
		assertEquals("5", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("alpha", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ALPHA_OR_BETA_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testAlphaBetaNotDirectlyPrecedingANaturalProduct3() throws ComponentGenerationException {
		Element substituent = new Element(SUBSTITUENT_EL);
		Element naturalProduct = new Element(GROUP_EL);
		naturalProduct.addAttribute(new Attribute(SUBTYPE_ATR, BIOCHEMICAL_SUBTYPE_VAL));
		naturalProduct.addAttribute(new Attribute(ALPHABETACLOCKWISEATOMORDERING_ATR, ""));
		substituent.appendChild(naturalProduct);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, ALPHA_OR_BETA_TYPE_VAL));
		substituent.appendChild(stereochem);
		stereochem.setValue("3beta,5alpha");
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(4, children.size());
		Element newStereochemistryEl = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getLocalName());
		assertEquals("3", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("beta", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ALPHA_OR_BETA_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
		
		newStereochemistryEl = children.get(2);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getLocalName());
		assertEquals("5", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("alpha", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ALPHA_OR_BETA_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
		
		Element newLocantEl = children.get(3);
		assertEquals(LOCANT_EL, newLocantEl.getLocalName());
		assertEquals("3,5", newLocantEl.getValue());
	}
	
	@Test
	public void testAlphaBetaStereoMixedWithNormalLocants() throws ComponentGenerationException {
		Element substituent = new Element(SUBSTITUENT_EL);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, ALPHA_OR_BETA_TYPE_VAL));
		substituent.appendChild(stereochem);
		stereochem.setValue("3beta,4,10,12alpha");
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(3, children.size());
		Element newStereochemistryEl = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getLocalName());
		assertEquals("3", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("beta", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ALPHA_OR_BETA_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
		
		newStereochemistryEl = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl.getLocalName());
		assertEquals("12", newStereochemistryEl.getAttributeValue(LOCANT_ATR));
		assertEquals("alpha", newStereochemistryEl.getAttributeValue(VALUE_ATR));
		assertEquals(ALPHA_OR_BETA_TYPE_VAL, newStereochemistryEl.getAttributeValue(TYPE_ATR));
		
		Element newLocantEl = children.get(2);
		assertEquals(LOCANT_EL, newLocantEl.getLocalName());
		assertEquals("3,4,10,12", newLocantEl.getValue());
	}
	
	//relative stereochemistry is currently treated the same as absolute stereochemistry
	@Test
	public void testRelativeStereoChemistry1() throws ComponentGenerationException {
		Element substituent = new Element(SUBSTITUENT_EL);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.appendChild(stereochem);
		stereochem.setValue("rel-(1R,3S,4S,7R)");
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(4, children.size());
		Element newStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl1.getLocalName());
		assertEquals("1", newStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("R", newStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl1.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl2 = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl2.getLocalName());
		assertEquals("3", newStereochemistryEl2.getAttributeValue(LOCANT_ATR));
		assertEquals("S", newStereochemistryEl2.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl2.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl3 = children.get(2);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl3.getLocalName());
		assertEquals("4", newStereochemistryEl3.getAttributeValue(LOCANT_ATR));
		assertEquals("S", newStereochemistryEl3.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl3.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl4 = children.get(3);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl4.getLocalName());
		assertEquals("7", newStereochemistryEl4.getAttributeValue(LOCANT_ATR));
		assertEquals("R", newStereochemistryEl4.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl4.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testRelativeStereoChemistry2() throws ComponentGenerationException {
		Element substituent = new Element(SUBSTITUENT_EL);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.appendChild(stereochem);
		stereochem.setValue("(1R*,3S*,4S*,7R*)");
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(4, children.size());
		Element newStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl1.getLocalName());
		assertEquals("1", newStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("R", newStereochemistryEl1.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl1.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl2 = children.get(1);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl2.getLocalName());
		assertEquals("3", newStereochemistryEl2.getAttributeValue(LOCANT_ATR));
		assertEquals("S", newStereochemistryEl2.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl2.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl3 = children.get(2);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl3.getLocalName());
		assertEquals("4", newStereochemistryEl3.getAttributeValue(LOCANT_ATR));
		assertEquals("S", newStereochemistryEl3.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl3.getAttributeValue(TYPE_ATR));
		
		Element newStereochemistryEl4 = children.get(3);
		assertEquals(STEREOCHEMISTRY_EL, newStereochemistryEl4.getLocalName());
		assertEquals("7", newStereochemistryEl4.getAttributeValue(LOCANT_ATR));
		assertEquals("R", newStereochemistryEl4.getAttributeValue(VALUE_ATR));
		assertEquals(R_OR_S_TYPE_VAL, newStereochemistryEl4.getAttributeValue(TYPE_ATR));
	}
	
	@Test
	public void testRelativeStereoChemistry3() throws ComponentGenerationException {
		Element substituent = new Element(SUBSTITUENT_EL);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.appendChild(stereochem);
		stereochem.setValue("rel-");
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(0, children.size());
	}
	
	//relativeCisTrans is only supported sufficiently to get constitutionally correct results i.e. locants extracted from the stereochemistry
	@Test
	public void testRelativeCisTrans() throws ComponentGenerationException {
		//c-4-
		Element substituent = new Element(SUBSTITUENT_EL);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, RELATIVECISTRANS_TYPE_VAL));
		stereochem.setValue("c-4-");
		substituent.appendChild(stereochem);
		processStereochemistry(substituent);
	
		List<Element> children = substituent.getChildElements();
		assertEquals(2, children.size());
		Element modifiedStereochemistryEl1 = children.get(0);
		assertEquals(STEREOCHEMISTRY_EL, modifiedStereochemistryEl1.getLocalName());
		assertEquals(null, modifiedStereochemistryEl1.getAttributeValue(LOCANT_ATR));
		assertEquals("c-4-", modifiedStereochemistryEl1.getValue());
		assertEquals(RELATIVECISTRANS_TYPE_VAL, modifiedStereochemistryEl1.getAttributeValue(TYPE_ATR));
		Element locant = children.get(1);
		assertEquals(LOCANT_EL, locant.getLocalName());
		assertEquals("4", locant.getValue());
	}
	
	//racemates are currently treated identically to completely undefined
	@Test
	public void testRacemate1() throws ComponentGenerationException {
		Element substituent = new Element(SUBSTITUENT_EL);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.appendChild(stereochem);
		stereochem.setValue("rac-(2R)");
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(0, children.size());
	}
	
	@Test
	public void testRacemate2() throws ComponentGenerationException {
		Element substituent = new Element(SUBSTITUENT_EL);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.appendChild(stereochem);
		stereochem.setValue("(RS)");
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(0, children.size());
	}
	
	@Test
	public void testRacemate3() throws ComponentGenerationException {
		Element substituent = new Element(SUBSTITUENT_EL);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.appendChild(stereochem);
		stereochem.setValue("(RS)");
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(0, children.size());
	}

	@Test
	public void testRacemate4() throws ComponentGenerationException {
		Element substituent = new Element(SUBSTITUENT_EL);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.appendChild(stereochem);
		stereochem.setValue("rac-(2R,4R)");
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(0, children.size());
	}
	
	@Test
	public void testRacemate5() throws ComponentGenerationException {
		Element substituent = new Element(SUBSTITUENT_EL);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.appendChild(stereochem);
		stereochem.setValue("(2RS,4RS)");
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(0, children.size());
	}
	
	@Test
	public void testRacemate6() throws ComponentGenerationException {
		Element substituent = new Element(SUBSTITUENT_EL);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.appendChild(stereochem);
		stereochem.setValue("rac-");
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(0, children.size());
	}
	
	@Test
	public void testRacemate7() throws ComponentGenerationException {
		Element substituent = new Element(SUBSTITUENT_EL);
		Element stereochem = new Element(STEREOCHEMISTRY_EL);
		stereochem.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
		substituent.appendChild(stereochem);
		stereochem.setValue("racemic-");
		processStereochemistry(substituent);

		List<Element> children = substituent.getChildElements();
		assertEquals(0, children.size());
	}
	
	private void processStereochemistry(Element subOrRoot) throws ComponentGenerationException {
		new ComponentGenerator(new NameToStructureConfig()).processStereochemistry(subOrRoot);
	}
}
