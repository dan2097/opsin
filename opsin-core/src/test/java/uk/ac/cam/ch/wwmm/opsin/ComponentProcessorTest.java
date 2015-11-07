package uk.ac.cam.ch.wwmm.opsin;

import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

import org.junit.Test;

import static org.junit.Assert.*;
import static org.mockito.Mockito.mock;

public class ComponentProcessorTest {

	@Test(expected=ComponentGenerationException.class)
	public void testSubtractiveWithNoGroupToAttachTo() throws ComponentGenerationException{
		Element word = new GroupingEl(WORD_EL);
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		word.addChild(substituent);
		Element substractivePrefix = new TokenEl(SUBTRACTIVEPREFIX_EL);
		substractivePrefix.addAttribute(new Attribute(TYPE_ATR, DEOXY_TYPE_VAL));
		substituent.addChild(substractivePrefix);
		ComponentProcessor.removeAndMoveToAppropriateGroupIfSubtractivePrefix(substituent);
	}
	
	@Test
	public void testSubtractiveWithBiochemicalToAttachTo() throws ComponentGenerationException{
		Element word = new GroupingEl(WORD_EL);
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element substractivePrefix = new TokenEl(SUBTRACTIVEPREFIX_EL);
		substractivePrefix.addAttribute(new Attribute(TYPE_ATR, DEOXY_TYPE_VAL));
		substituent.addChild(substractivePrefix);
		word.addChild(substituent);
		Element root = new GroupingEl(ROOT_EL);
		word.addChild(root);
		Element group = new TokenEl(GROUP_EL);
		group.addAttribute(new Attribute(SUBTYPE_ATR, BIOCHEMICAL_SUBTYPE_VAL));
		root.addChild(group);

		ComponentProcessor.removeAndMoveToAppropriateGroupIfSubtractivePrefix(substituent);
		assertEquals("Substractive prefix should of been detached", null, substituent.getParent());
		assertEquals(2, root.getChildCount());
		assertEquals(substractivePrefix, root.getChildElements().get(0));
	}
	
	@Test
	public void testSubtractiveRightMostPreferred() throws ComponentGenerationException{
		Element word = new GroupingEl(WORD_EL);
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element substractivePrefix = new TokenEl(SUBTRACTIVEPREFIX_EL);
		substractivePrefix.addAttribute(new Attribute(TYPE_ATR, DEOXY_TYPE_VAL));
		substituent.addChild(substractivePrefix);
		word.addChild(substituent);
		Element substituent2 = new GroupingEl(SUBSTITUENT_EL);
		Element group1 = new TokenEl(GROUP_EL);
		group1.addAttribute(new Attribute(TYPE_ATR, SIMPLEGROUP_SUBTYPE_VAL));
		group1.addAttribute(new Attribute(SUBTYPE_ATR, SIMPLEGROUP_SUBTYPE_VAL));
		substituent2.addChild(group1);
		word.addChild(substituent2);
		Element root = new GroupingEl(ROOT_EL);
		word.addChild(root);
		Element group2 = new TokenEl(GROUP_EL);
		group2.addAttribute(new Attribute(TYPE_ATR, SIMPLEGROUP_SUBTYPE_VAL));
		group2.addAttribute(new Attribute(SUBTYPE_ATR, BIOCHEMICAL_SUBTYPE_VAL));
		root.addChild(group2);

		ComponentProcessor.removeAndMoveToAppropriateGroupIfSubtractivePrefix(substituent);
		assertEquals("Substractive prefix should of been detached", null, substituent.getParent());
		assertEquals(2, root.getChildCount());
		assertEquals(substractivePrefix, root.getChildElements().get(0));
	}
	
	@Test
	public void testSubtractiveBiochemicalPreferredToRightMost() throws ComponentGenerationException{
		Element word = new GroupingEl(WORD_EL);
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element substractivePrefix = new TokenEl(SUBTRACTIVEPREFIX_EL);
		substractivePrefix.addAttribute(new Attribute(TYPE_ATR, DEOXY_TYPE_VAL));
		substituent.addChild(substractivePrefix);
		word.addChild(substituent);
		Element substituent2 = new GroupingEl(SUBSTITUENT_EL);
		Element group1 = new TokenEl(GROUP_EL);
		group1.addAttribute(new Attribute(SUBTYPE_ATR, BIOCHEMICAL_SUBTYPE_VAL));
		substituent2.addChild(group1);
		word.addChild(substituent2);
		Element root = new GroupingEl(ROOT_EL);
		word.addChild(root);
		Element group2 = new TokenEl(GROUP_EL);
		group2.addAttribute(new Attribute(SUBTYPE_ATR, SIMPLEGROUP_SUBTYPE_VAL));
		root.addChild(group2);

		ComponentProcessor.removeAndMoveToAppropriateGroupIfSubtractivePrefix(substituent);
		assertEquals("Substractive prefix should of been detached", null, substituent.getParent());
		assertEquals(1, root.getChildCount());
		assertEquals(2, substituent2.getChildCount());
		assertEquals(substractivePrefix, substituent2.getChildElements().get(0));
	}
	
	@Test
	public void testSubtractiveWithMultiplierAndLocants() throws ComponentGenerationException{
		Element word = new GroupingEl(WORD_EL);
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element locant = new TokenEl(LOCANT_EL);
		substituent.addChild(locant);
		Element multiplier = new TokenEl(MULTIPLIER_EL);
		substituent.addChild(multiplier);
		Element substractivePrefix = new TokenEl(SUBTRACTIVEPREFIX_EL);
		substractivePrefix.addAttribute(new Attribute(TYPE_ATR, DEOXY_TYPE_VAL));
		substituent.addChild(substractivePrefix);
		word.addChild(substituent);
		Element root = new GroupingEl(ROOT_EL);
		word.addChild(root);
		Element group = new TokenEl(GROUP_EL);
		group.addAttribute(new Attribute(SUBTYPE_ATR, BIOCHEMICAL_SUBTYPE_VAL));
		root.addChild(group);

		ComponentProcessor.removeAndMoveToAppropriateGroupIfSubtractivePrefix(substituent);
		assertEquals("Substractive prefix should of been detached", null, substituent.getParent());
		assertEquals(4, root.getChildCount());
		assertEquals(locant, root.getChildElements().get(0));
		assertEquals(multiplier, root.getChildElements().get(1));
		assertEquals(substractivePrefix, root.getChildElements().get(2));
	}

	@Test
	public void testDLStereochemistryLOnAminoAcid() throws ComponentGenerationException, StructureBuildingException{
		BuildState state = new BuildState(mock(NameToStructureConfig.class));
		Fragment f = state.fragManager.buildSMILES("N[C@@H](C)C");
		Element aminoAcidEl = new TokenEl(GROUP_EL);
		aminoAcidEl.setFrag(f);
		int parityBefore = f.getAtomByID(2).getAtomParity().getParity();
		ComponentProcessor processor = new ComponentProcessor(mock(SuffixRulesLookup.class), state);
		processor.applyDlStereochemistryToAminoAcid(aminoAcidEl, "l");
		assertEquals(parityBefore, f.getAtomByID(2).getAtomParity().getParity());
	}
	
	@Test
	public void testDLStereochemistryDOnAminoAcid() throws ComponentGenerationException, StructureBuildingException{
		BuildState state = new BuildState(mock(NameToStructureConfig.class));
		Fragment f = state.fragManager.buildSMILES("N[C@@H](C)C");
		Element aminoAcidEl = new TokenEl(GROUP_EL);
		aminoAcidEl.setFrag(f);
		int parityBefore = f.getAtomByID(2).getAtomParity().getParity();
		ComponentProcessor processor = new ComponentProcessor(mock(SuffixRulesLookup.class), state);
		processor.applyDlStereochemistryToAminoAcid(aminoAcidEl, "d");
		assertEquals(parityBefore, -f.getAtomByID(2).getAtomParity().getParity());
	}
	
	@Test
	public void testDLStereochemistryDLOnAminoAcid() throws ComponentGenerationException, StructureBuildingException{
		BuildState state = new BuildState(mock(NameToStructureConfig.class));
		Fragment f = state.fragManager.buildSMILES("N[C@@H](C)C");
		Element aminoAcidEl = new TokenEl(GROUP_EL);
		aminoAcidEl.setFrag(f);
		ComponentProcessor processor = new ComponentProcessor(mock(SuffixRulesLookup.class), state);
		processor.applyDlStereochemistryToAminoAcid(aminoAcidEl, "dl");
		assertEquals(null, f.getAtomByID(2).getAtomParity());
	}
	
	@Test(expected=ComponentGenerationException.class)
	public void testDLStereochemistryDOnAchiralAminoAcid() throws ComponentGenerationException, StructureBuildingException{
		BuildState state = new BuildState(mock(NameToStructureConfig.class));
		Fragment f = state.fragManager.buildSMILES("NC(C)C");
		Element aminoAcidEl = new TokenEl(GROUP_EL);
		aminoAcidEl.setFrag(f);
		ComponentProcessor processor = new ComponentProcessor(mock(SuffixRulesLookup.class), state);
		processor.applyDlStereochemistryToAminoAcid(aminoAcidEl, "d");
	}
	
	@Test
	public void testDLStereochemistryLOnCarbohydrate() throws ComponentGenerationException, StructureBuildingException{
		BuildState state = new BuildState(mock(NameToStructureConfig.class));
		Fragment f = state.fragManager.buildSMILES("N[C@@H](C)C");
		Element carbohydrateEl = new TokenEl(GROUP_EL);
		carbohydrateEl.setFrag(f);
		int parityBefore = f.getAtomByID(2).getAtomParity().getParity();
		ComponentProcessor processor = new ComponentProcessor(mock(SuffixRulesLookup.class), state);
		processor.applyDlStereochemistryToCarbohydrate(carbohydrateEl, "l");
		assertEquals(parityBefore, -f.getAtomByID(2).getAtomParity().getParity());
	}
	
	@Test
	public void testDLStereochemistryDOnCarbohydrate() throws ComponentGenerationException, StructureBuildingException{
		BuildState state = new BuildState(mock(NameToStructureConfig.class));
		Fragment f = state.fragManager.buildSMILES("N[C@@H](C)C");
		Element carbohydrateEl = new TokenEl(GROUP_EL);
		carbohydrateEl.setFrag(f);
		int parityBefore = f.getAtomByID(2).getAtomParity().getParity();
		ComponentProcessor processor = new ComponentProcessor(mock(SuffixRulesLookup.class), state);
		processor.applyDlStereochemistryToCarbohydrate(carbohydrateEl, "d");
		assertEquals(parityBefore, f.getAtomByID(2).getAtomParity().getParity());
	}
	
	@Test
	public void testDLStereochemistryInvertedNaturalOnCarbohydrate1() throws ComponentGenerationException, StructureBuildingException{
		BuildState state = new BuildState(mock(NameToStructureConfig.class));
		Fragment f = state.fragManager.buildSMILES("N[C@@H](C)C");
		Element carbohydrateEl = new TokenEl(GROUP_EL);
		carbohydrateEl.addAttribute(new Attribute(NATURALENTISOPPOSITE_ATR, "yes"));
		carbohydrateEl.setFrag(f);
		int parityBefore = f.getAtomByID(2).getAtomParity().getParity();
		ComponentProcessor processor = new ComponentProcessor(mock(SuffixRulesLookup.class), state);
		processor.applyDlStereochemistryToCarbohydrate(carbohydrateEl, "l");
		assertEquals(parityBefore, f.getAtomByID(2).getAtomParity().getParity());
	}
	
	@Test
	public void testDLStereochemistryInvertedNaturalOnCarbohydrate2() throws ComponentGenerationException, StructureBuildingException{
		BuildState state = new BuildState(mock(NameToStructureConfig.class));
		Fragment f = state.fragManager.buildSMILES("N[C@@H](C)C");
		Element carbohydrateEl = new TokenEl(GROUP_EL);
		carbohydrateEl.addAttribute(new Attribute(NATURALENTISOPPOSITE_ATR, "yes"));
		carbohydrateEl.setFrag(f);
		int parityBefore = f.getAtomByID(2).getAtomParity().getParity();
		ComponentProcessor processor = new ComponentProcessor(mock(SuffixRulesLookup.class), state);
		processor.applyDlStereochemistryToCarbohydrate(carbohydrateEl, "d");
		assertEquals(parityBefore, -f.getAtomByID(2).getAtomParity().getParity());
	}

	@Test
	public void testDStereochemistryDOnCarbohydratePrefix() throws ComponentGenerationException, StructureBuildingException{
		Element prefix = new TokenEl(STEREOCHEMISTRY_EL);
		prefix.addAttribute(new Attribute(TYPE_ATR, CARBOHYDRATECONFIGURATIONPREFIX_TYPE_VAL));
		prefix.addAttribute(new Attribute(VALUE_ATR, "l/r"));//D-threo
		ComponentProcessor.applyDlStereochemistryToCarbohydrateConfigurationalPrefix(prefix, "d");
		assertEquals("l/r", prefix.getAttributeValue(VALUE_ATR));
	}
	
	@Test
	public void testLStereochemistryDOnCarbohydratePrefix() throws ComponentGenerationException, StructureBuildingException{
		Element prefix = new TokenEl(STEREOCHEMISTRY_EL);
		prefix.addAttribute(new Attribute(TYPE_ATR, CARBOHYDRATECONFIGURATIONPREFIX_TYPE_VAL));
		prefix.addAttribute(new Attribute(VALUE_ATR, "r/l"));
		ComponentProcessor.applyDlStereochemistryToCarbohydrateConfigurationalPrefix(prefix, "l");
		assertEquals("l/r", prefix.getAttributeValue(VALUE_ATR));
	}
	
	@Test
	public void testDLStereochemistryDOnCarbohydratePrefix() throws ComponentGenerationException, StructureBuildingException{
		Element prefix = new TokenEl(STEREOCHEMISTRY_EL);
		prefix.addAttribute(new Attribute(TYPE_ATR, CARBOHYDRATECONFIGURATIONPREFIX_TYPE_VAL));
		prefix.addAttribute(new Attribute(VALUE_ATR, "l/r"));
		ComponentProcessor.applyDlStereochemistryToCarbohydrateConfigurationalPrefix(prefix, "dl");
		assertEquals("?/?", prefix.getAttributeValue(VALUE_ATR));
	}
	
}
