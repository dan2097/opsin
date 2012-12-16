package uk.ac.cam.ch.wwmm.opsin;

import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;
import static junit.framework.Assert.*;
import nu.xom.Attribute;
import nu.xom.Element;

import org.junit.Test;

public class ComponentProcessorTest {

	@Test(expected=ComponentGenerationException.class)
	public void testSubtractiveWithNoGroupToAttachTo() throws ComponentGenerationException{
		Element word = new Element(WORD_EL);
		Element substituent = new Element(SUBSTITUENT_EL);
		word.appendChild(substituent);
		Element substractivePrefix = new Element(SUBTRACTIVEPREFIX_EL);
		substituent.appendChild(substractivePrefix);
		ComponentProcessor.removeAndMoveToAppropriateGroupIfSubstractivePrefix(substituent);
	}
	
	@Test
	public void testSubtractiveWithBiochemicalToAttachTo() throws ComponentGenerationException{
		Element word = new Element(WORD_EL);
		Element substituent = new Element(SUBSTITUENT_EL);
		Element substractivePrefix = new Element(SUBTRACTIVEPREFIX_EL);
		substituent.appendChild(substractivePrefix);
		word.appendChild(substituent);
		Element root = new Element(ROOT_EL);
		word.appendChild(root);
		Element group = new Element(GROUP_EL);
		group.addAttribute(new Attribute(SUBTYPE_ATR, BIOCHEMICAL_SUBTYPE_VAL));
		root.appendChild(group);

		ComponentProcessor.removeAndMoveToAppropriateGroupIfSubstractivePrefix(substituent);
		assertEquals("Substractive prefix should of been detached", null, substituent.getParent());
		assertEquals(2, root.getChildElements().size());
		assertEquals(substractivePrefix, root.getChildElements().get(0));
	}
	
	@Test
	public void testSubtractiveRightMostPreferred() throws ComponentGenerationException{
		Element word = new Element(WORD_EL);
		Element substituent = new Element(SUBSTITUENT_EL);
		Element substractivePrefix = new Element(SUBTRACTIVEPREFIX_EL);
		substituent.appendChild(substractivePrefix);
		word.appendChild(substituent);
		Element substituent2 = new Element(SUBSTITUENT_EL);
		Element group1 = new Element(GROUP_EL);
		group1.addAttribute(new Attribute(TYPE_ATR, SIMPLEGROUP_SUBTYPE_VAL));
		group1.addAttribute(new Attribute(SUBTYPE_ATR, SIMPLEGROUP_SUBTYPE_VAL));
		substituent2.appendChild(group1);
		word.appendChild(substituent2);
		Element root = new Element(ROOT_EL);
		word.appendChild(root);
		Element group2 = new Element(GROUP_EL);
		group2.addAttribute(new Attribute(TYPE_ATR, SIMPLEGROUP_SUBTYPE_VAL));
		group2.addAttribute(new Attribute(SUBTYPE_ATR, BIOCHEMICAL_SUBTYPE_VAL));
		root.appendChild(group2);

		ComponentProcessor.removeAndMoveToAppropriateGroupIfSubstractivePrefix(substituent);
		assertEquals("Substractive prefix should of been detached", null, substituent.getParent());
		assertEquals(2, root.getChildElements().size());
		assertEquals(substractivePrefix, root.getChildElements().get(0));
	}
	
	@Test
	public void testSubtractiveBiochemicalPreferredToRightMost() throws ComponentGenerationException{
		Element word = new Element(WORD_EL);
		Element substituent = new Element(SUBSTITUENT_EL);
		Element substractivePrefix = new Element(SUBTRACTIVEPREFIX_EL);
		substituent.appendChild(substractivePrefix);
		word.appendChild(substituent);
		Element substituent2 = new Element(SUBSTITUENT_EL);
		Element group1 = new Element(GROUP_EL);
		group1.addAttribute(new Attribute(SUBTYPE_ATR, BIOCHEMICAL_SUBTYPE_VAL));
		substituent2.appendChild(group1);
		word.appendChild(substituent2);
		Element root = new Element(ROOT_EL);
		word.appendChild(root);
		Element group2 = new Element(GROUP_EL);
		group2.addAttribute(new Attribute(SUBTYPE_ATR, SIMPLEGROUP_SUBTYPE_VAL));
		root.appendChild(group2);

		ComponentProcessor.removeAndMoveToAppropriateGroupIfSubstractivePrefix(substituent);
		assertEquals("Substractive prefix should of been detached", null, substituent.getParent());
		assertEquals(1, root.getChildElements().size());
		assertEquals(2, substituent2.getChildElements().size());
		assertEquals(substractivePrefix, substituent2.getChildElements().get(0));
	}
	
	@Test
	public void testSubtractiveWithMultiplierAndLocants() throws ComponentGenerationException{
		Element word = new Element(WORD_EL);
		Element substituent = new Element(SUBSTITUENT_EL);
		Element locant = new Element(LOCANT_EL);
		substituent.appendChild(locant);
		Element multiplier = new Element(MULTIPLIER_EL);
		substituent.appendChild(multiplier);
		Element substractivePrefix = new Element(SUBTRACTIVEPREFIX_EL);
		substituent.appendChild(substractivePrefix);
		word.appendChild(substituent);
		Element root = new Element(ROOT_EL);
		word.appendChild(root);
		Element group = new Element(GROUP_EL);
		group.addAttribute(new Attribute(SUBTYPE_ATR, BIOCHEMICAL_SUBTYPE_VAL));
		root.appendChild(group);

		ComponentProcessor.removeAndMoveToAppropriateGroupIfSubstractivePrefix(substituent);
		assertEquals("Substractive prefix should of been detached", null, substituent.getParent());
		assertEquals(4, root.getChildElements().size());
		assertEquals(locant, root.getChildElements().get(0));
		assertEquals(multiplier, root.getChildElements().get(1));
		assertEquals(substractivePrefix, root.getChildElements().get(2));
	}

	@Test
	public void testDLStereochemistryLOnAminoAcid() throws ComponentGenerationException, StructureBuildingException{
		FragmentManager fm = new FragmentManager(new SMILESFragmentBuilder(), new IDManager());
		Fragment f = fm.buildSMILES("N[C@@H](C)C");
		int parityBefore = f.getAtomByID(2).getAtomParity().getParity();
		ComponentProcessor.applyDlStereochemistryToAminoAcid(f, "l");
		assertEquals(parityBefore, f.getAtomByID(2).getAtomParity().getParity());
	}
	
	@Test
	public void testDLStereochemistryDOnAminoAcid() throws ComponentGenerationException, StructureBuildingException{
		FragmentManager fm = new FragmentManager(new SMILESFragmentBuilder(), new IDManager());
		Fragment f = fm.buildSMILES("N[C@@H](C)C");
		int parityBefore = f.getAtomByID(2).getAtomParity().getParity();
		ComponentProcessor.applyDlStereochemistryToAminoAcid(f, "d");
		assertEquals(parityBefore, -f.getAtomByID(2).getAtomParity().getParity());
	}
	
	@Test
	public void testDLStereochemistryDLOnAminoAcid() throws ComponentGenerationException, StructureBuildingException{
		FragmentManager fm = new FragmentManager(new SMILESFragmentBuilder(), new IDManager());
		Fragment f = fm.buildSMILES("N[C@@H](C)C");
		ComponentProcessor.applyDlStereochemistryToAminoAcid(f, "dl");
		assertEquals(null, f.getAtomByID(2).getAtomParity());
	}
	
	@Test(expected=ComponentGenerationException.class)
	public void testDLStereochemistryDOnAchiralAminoAcid() throws ComponentGenerationException, StructureBuildingException{
		FragmentManager fm = new FragmentManager(new SMILESFragmentBuilder(), new IDManager());
		Fragment f = fm.buildSMILES("NC(C)C");
		ComponentProcessor.applyDlStereochemistryToAminoAcid(f, "d");
	}
	
	
	@Test
	public void testDLStereochemistryLOnCarbohydrate() throws ComponentGenerationException, StructureBuildingException{
		FragmentManager fm = new FragmentManager(new SMILESFragmentBuilder(), new IDManager());
		Fragment f = fm.buildSMILES("N[C@@H](C)C");
		int parityBefore = f.getAtomByID(2).getAtomParity().getParity();
		ComponentProcessor.applyDlStereochemistryToCarbohydrate(f, "l");
		assertEquals(parityBefore, -f.getAtomByID(2).getAtomParity().getParity());
	}
	
	@Test
	public void testDLStereochemistryDOnCarbohydrate() throws ComponentGenerationException, StructureBuildingException{
		FragmentManager fm = new FragmentManager(new SMILESFragmentBuilder(), new IDManager());
		Fragment f = fm.buildSMILES("N[C@@H](C)C");
		int parityBefore = f.getAtomByID(2).getAtomParity().getParity();
		ComponentProcessor.applyDlStereochemistryToCarbohydrate(f, "d");
		assertEquals(parityBefore, f.getAtomByID(2).getAtomParity().getParity());
	}

	@Test
	public void testDStereochemistryDOnCarbohydratePrefix() throws ComponentGenerationException, StructureBuildingException{
		Element prefix = new Element(STEREOCHEMISTRY_EL);
		prefix.addAttribute(new Attribute(TYPE_ATR, CARBOHYDRATECONFIGURATIONPREFIX_TYPE_VAL));
		prefix.addAttribute(new Attribute(VALUE_ATR, "l/r"));//D-threo
		ComponentProcessor.applyDlStereochemistryToCarbohydrateConfigurationalPrefix(prefix, "d");
		assertEquals("l/r", prefix.getAttributeValue(VALUE_ATR));
	}
	
	@Test
	public void testLStereochemistryDOnCarbohydratePrefix() throws ComponentGenerationException, StructureBuildingException{
		Element prefix = new Element(STEREOCHEMISTRY_EL);
		prefix.addAttribute(new Attribute(TYPE_ATR, CARBOHYDRATECONFIGURATIONPREFIX_TYPE_VAL));
		prefix.addAttribute(new Attribute(VALUE_ATR, "r/l"));
		ComponentProcessor.applyDlStereochemistryToCarbohydrateConfigurationalPrefix(prefix, "l");
		assertEquals("l/r", prefix.getAttributeValue(VALUE_ATR));
	}
	
	@Test
	public void testDLStereochemistryDOnCarbohydratePrefix() throws ComponentGenerationException, StructureBuildingException{
		Element prefix = new Element(STEREOCHEMISTRY_EL);
		prefix.addAttribute(new Attribute(TYPE_ATR, CARBOHYDRATECONFIGURATIONPREFIX_TYPE_VAL));
		prefix.addAttribute(new Attribute(VALUE_ATR, "l/r"));
		ComponentProcessor.applyDlStereochemistryToCarbohydrateConfigurationalPrefix(prefix, "dl");
		assertEquals("?/?", prefix.getAttributeValue(VALUE_ATR));
	}
	
}
