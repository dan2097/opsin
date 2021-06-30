package uk.ac.cam.ch.wwmm.opsin;

import java.util.Set;

import org.junit.jupiter.api.Test;

import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.mock;

public class StructureBuildingMethodsTest {
	
	@Test
	public void bracketedPrimeNotSpecialCase() {
		Element word = new GroupingEl(WORD_EL);
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		word.addChild(substituent);
		assertEquals(null, StructureBuildingMethods.checkForBracketedPrimedLocantSpecialCase(substituent, "4"));
		assertEquals(null, StructureBuildingMethods.checkForBracketedPrimedLocantSpecialCase(substituent, "4'"));
		assertEquals(null, StructureBuildingMethods.checkForBracketedPrimedLocantSpecialCase(substituent, "4''"));
	}
	
	@Test
	public void bracketedPrimeSpecialCase1() {
		Element word = new GroupingEl(WORD_EL);
		Element bracket = new GroupingEl(BRACKET_EL);
		word.addChild(bracket);
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		bracket.addChild(substituent);
		assertEquals(null, StructureBuildingMethods.checkForBracketedPrimedLocantSpecialCase(substituent, "4"));
		assertEquals("4", StructureBuildingMethods.checkForBracketedPrimedLocantSpecialCase(substituent, "4'"));
		assertEquals(null, StructureBuildingMethods.checkForBracketedPrimedLocantSpecialCase(substituent, "4''"));
		bracket.addAttribute(new Attribute(TYPE_ATR, IMPLICIT_TYPE_VAL));
		assertEquals(null, StructureBuildingMethods.checkForBracketedPrimedLocantSpecialCase(substituent, "4"));
		assertEquals(null, StructureBuildingMethods.checkForBracketedPrimedLocantSpecialCase(substituent, "4'"));
		assertEquals(null, StructureBuildingMethods.checkForBracketedPrimedLocantSpecialCase(substituent, "4''"));
	}
	
	@Test
	public void bracketedPrimeSpecialCase2() {
		Element word = new GroupingEl(WORD_EL);
		Element bracket = new GroupingEl(BRACKET_EL);
		word.addChild(bracket);
		Element bracket2 = new GroupingEl(BRACKET_EL);
		bracket.addChild(bracket2);
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		bracket2.addChild(substituent);
		assertEquals(null, StructureBuildingMethods.checkForBracketedPrimedLocantSpecialCase(substituent, "4"));
		assertEquals(null, StructureBuildingMethods.checkForBracketedPrimedLocantSpecialCase(substituent, "4'"));
		assertEquals("4", StructureBuildingMethods.checkForBracketedPrimedLocantSpecialCase(substituent, "4''"));
		bracket2.addAttribute(new Attribute(TYPE_ATR, IMPLICIT_TYPE_VAL));
		assertEquals(null, StructureBuildingMethods.checkForBracketedPrimedLocantSpecialCase(substituent, "4"));
		assertEquals("4", StructureBuildingMethods.checkForBracketedPrimedLocantSpecialCase(substituent, "4'"));
		assertEquals(null, StructureBuildingMethods.checkForBracketedPrimedLocantSpecialCase(substituent, "4''"));
	}
	
	@Test
	public void notPhosphoSubstitution() throws StructureBuildingException {
		//standard unlocanted substitution
		BuildState state = new BuildState(mock(NameToStructureConfig.class));
		Element word = new GroupingEl(WORD_EL);
		
		Element amino = new TokenEl(GROUP_EL);
		Fragment aminoFrag = state.fragManager.buildSMILES("-N");
		amino.setFrag(aminoFrag);
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		substituent.addChild(amino);
		
		Element methanol = new TokenEl(GROUP_EL);
		methanol.setFrag(state.fragManager.buildSMILES("CO"));
		Element root = new GroupingEl(ROOT_EL);
		root.addChild(methanol);

		word.addChild(substituent);
		word.addChild(root);
		StructureBuildingMethods.resolveRootOrSubstituentUnLocanted(state, substituent);
		
		Set<Bond> interFragmentBonds =  state.fragManager.getInterFragmentBonds(aminoFrag);
		assertEquals(1, interFragmentBonds.size());
		assertEquals(ChemEl.C, interFragmentBonds.iterator().next().getOtherAtom(aminoFrag.getFirstAtom()).getElement());
	}
	
	@Test
	public void phosphoUnlocantedSubstitution() throws StructureBuildingException {
		BuildState state = new BuildState(mock(NameToStructureConfig.class));
		Element word = new GroupingEl(WORD_EL);
		
		Element phospho = new TokenEl(GROUP_EL);
		phospho.addAttribute(new Attribute(SUBTYPE_ATR, PHOSPHO_SUBTYPE_VAL));
		Fragment phosphoFrag = state.fragManager.buildSMILES("-P(=O)O");
		phospho.setFrag(phosphoFrag);
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		substituent.addChild(phospho);
		
		Element methanol = new TokenEl(GROUP_EL);
		methanol.setFrag(state.fragManager.buildSMILES("CO"));
		Element root = new GroupingEl(ROOT_EL);
		root.addChild(methanol);

		word.addChild(substituent);
		word.addChild(root);
		StructureBuildingMethods.resolveRootOrSubstituentUnLocanted(state, substituent);
		
		Set<Bond> interFragmentBonds =  state.fragManager.getInterFragmentBonds(phosphoFrag);
		assertEquals(1, interFragmentBonds.size());
		assertEquals(ChemEl.O, interFragmentBonds.iterator().next().getOtherAtom(phosphoFrag.getFirstAtom()).getElement());
	}
	
	@Test
	public void phosphoLocantedSubstitution() throws StructureBuildingException {
		BuildState state = new BuildState(mock(NameToStructureConfig.class));
		Element word = new GroupingEl(WORD_EL);
		
		Element phospho = new TokenEl(GROUP_EL);
		phospho.addAttribute(new Attribute(SUBTYPE_ATR, PHOSPHO_SUBTYPE_VAL));
		Fragment phosphoFrag = state.fragManager.buildSMILES("-P(=O)O");
		phospho.setFrag(phosphoFrag);
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		substituent.addAttribute(new Attribute(LOCANT_ATR, "4"));
		substituent.addChild(phospho);
		
		Element methanol = new TokenEl(GROUP_EL);
		methanol.setFrag(state.fragManager.buildSMILES("CCCCO",methanol,"1/2/3/4/"));
		Element root = new GroupingEl(ROOT_EL);
		root.addChild(methanol);

		word.addChild(substituent);
		word.addChild(root);
		StructureBuildingMethods.resolveRootOrSubstituentLocanted(state, substituent);
		
		Set<Bond> interFragmentBonds =  state.fragManager.getInterFragmentBonds(phosphoFrag);
		assertEquals(1, interFragmentBonds.size());
		assertEquals(ChemEl.O, interFragmentBonds.iterator().next().getOtherAtom(phosphoFrag.getFirstAtom()).getElement());
	}
}
