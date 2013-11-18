package uk.ac.cam.ch.wwmm.opsin;

import java.util.Set;

import org.junit.Test;

import nu.xom.Attribute;
import nu.xom.Element;
import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;
import static junit.framework.Assert.assertEquals;
import static org.mockito.Mockito.mock;

public class StructureBuildingMethodsTest {
	
	@Test
	public void bracketedPrimeNotSpecialCase() {
		Element word = new Element(WORD_EL);
		Element substituent = new Element(SUBSTITUENT_EL);
		word.appendChild(substituent);
		assertEquals(null, StructureBuildingMethods.checkForBracketedPrimedLocantSpecialCase(substituent, "4"));
		assertEquals(null, StructureBuildingMethods.checkForBracketedPrimedLocantSpecialCase(substituent, "4'"));
		assertEquals(null, StructureBuildingMethods.checkForBracketedPrimedLocantSpecialCase(substituent, "4''"));
	}
	
	@Test
	public void bracketedPrimeSpecialCase1() {
		Element word = new Element(WORD_EL);
		Element bracket = new Element(BRACKET_EL);
		word.appendChild(bracket);
		Element substituent = new Element(SUBSTITUENT_EL);
		bracket.appendChild(substituent);
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
		Element word = new Element(WORD_EL);
		Element bracket = new Element(BRACKET_EL);
		word.appendChild(bracket);
		Element bracket2 = new Element(BRACKET_EL);
		bracket.appendChild(bracket2);
		Element substituent = new Element(SUBSTITUENT_EL);
		bracket2.appendChild(substituent);
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
		Element word = new Element(WORD_EL);
		
		Element amino = new Element(GROUP_EL);
		Fragment aminoFrag = state.fragManager.buildSMILES("-N");
		state.xmlFragmentMap.put(amino, aminoFrag);
		Element substituent = new Element(SUBSTITUENT_EL);
		substituent.appendChild(amino);
		
		Element methanol = new Element(GROUP_EL);
		state.xmlFragmentMap.put(methanol, state.fragManager.buildSMILES("CO"));
		Element root = new Element(ROOT_EL);
		root.appendChild(methanol);

		word.appendChild(substituent);
		word.appendChild(root);
		StructureBuildingMethods.resolveRootOrSubstituentUnLocanted(state, substituent);
		
		Set<Bond> interFragmentBonds =  state.fragManager.getInterFragmentBonds(aminoFrag);
		assertEquals(1, interFragmentBonds.size());
		assertEquals("C", interFragmentBonds.iterator().next().getOtherAtom(aminoFrag.getFirstAtom()).getElement());
	}
	
	@Test
	public void phosphoUnlocantedSubstitution() throws StructureBuildingException {
		BuildState state = new BuildState(mock(NameToStructureConfig.class));
		Element word = new Element(WORD_EL);
		
		Element phospho = new Element(GROUP_EL);
		phospho.addAttribute(new Attribute(SUBTYPE_ATR, PHOSPHO_SUBTYPE_VAL));
		Fragment phosphoFrag = state.fragManager.buildSMILES("-P(=O)O");
		state.xmlFragmentMap.put(phospho, phosphoFrag);
		Element substituent = new Element(SUBSTITUENT_EL);
		substituent.appendChild(phospho);
		
		Element methanol = new Element(GROUP_EL);
		state.xmlFragmentMap.put(methanol, state.fragManager.buildSMILES("CO"));
		Element root = new Element(ROOT_EL);
		root.appendChild(methanol);

		word.appendChild(substituent);
		word.appendChild(root);
		StructureBuildingMethods.resolveRootOrSubstituentUnLocanted(state, substituent);
		
		Set<Bond> interFragmentBonds =  state.fragManager.getInterFragmentBonds(phosphoFrag);
		assertEquals(1, interFragmentBonds.size());
		assertEquals("O", interFragmentBonds.iterator().next().getOtherAtom(phosphoFrag.getFirstAtom()).getElement());
	}
	
	@Test
	public void phosphoLocantedSubstitution() throws StructureBuildingException {
		BuildState state = new BuildState(mock(NameToStructureConfig.class));
		Element word = new Element(WORD_EL);
		
		Element phospho = new Element(GROUP_EL);
		phospho.addAttribute(new Attribute(SUBTYPE_ATR, PHOSPHO_SUBTYPE_VAL));
		Fragment phosphoFrag = state.fragManager.buildSMILES("-P(=O)O");
		state.xmlFragmentMap.put(phospho, phosphoFrag);
		Element substituent = new Element(SUBSTITUENT_EL);
		substituent.addAttribute(new Attribute(LOCANT_ATR, "4"));
		substituent.appendChild(phospho);
		
		Element methanol = new Element(GROUP_EL);
		state.xmlFragmentMap.put(methanol, state.fragManager.buildSMILES("CCCCO","group","1/2/3/4/"));
		Element root = new Element(ROOT_EL);
		root.appendChild(methanol);

		word.appendChild(substituent);
		word.appendChild(root);
		StructureBuildingMethods.resolveRootOrSubstituentLocanted(state, substituent);
		
		Set<Bond> interFragmentBonds =  state.fragManager.getInterFragmentBonds(phosphoFrag);
		assertEquals(1, interFragmentBonds.size());
		assertEquals("O", interFragmentBonds.iterator().next().getOtherAtom(phosphoFrag.getFirstAtom()).getElement());
	}
}
