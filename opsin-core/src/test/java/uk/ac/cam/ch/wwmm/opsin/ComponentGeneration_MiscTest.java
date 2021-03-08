package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.Assert.assertEquals;
import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

import org.junit.Test;

public class ComponentGeneration_MiscTest {

	@Test(expected=ComponentGenerationException.class)
	public void testRejectSingleComponentSaltComponent() throws ComponentGenerationException {
		//reject "hydrate"
		Element molecule = new GroupingEl(MOLECULE_EL);
		Element wordRule = new GroupingEl(WORDRULE_EL);
		Element word = new GroupingEl(WORD_EL);
		Element root = new GroupingEl(ROOT_EL);
		Element group = new TokenEl(GROUP_EL);
		group.addAttribute(new Attribute(TYPE_ATR, SIMPLEGROUP_TYPE_VAL));
		group.addAttribute(new Attribute(SUBTYPE_ATR, SALTCOMPONENT_SUBTYPE_VAL));
		root.addChild(group);
		word.addChild(root);
		wordRule.addChild(word);
		molecule.addChild(wordRule);
		processComponents(molecule);
	}

	@Test
	public void testNumericallyMultipliedSaltComponent() throws ComponentGenerationException {
		Element molecule = new GroupingEl(MOLECULE_EL);
		molecule.addChild(new GroupingEl(WORDRULE_EL));
		
		Element wordRule = new GroupingEl(WORDRULE_EL);
		Element word = new GroupingEl(WORD_EL);
		Element root = new GroupingEl(ROOT_EL);
		Element group = new TokenEl(GROUP_EL);
		group.addAttribute(new Attribute(TYPE_ATR, SIMPLEGROUP_TYPE_VAL));
		group.addAttribute(new Attribute(SUBTYPE_ATR, SALTCOMPONENT_SUBTYPE_VAL));
		group.setValue("2hcl");
		root.addChild(group);
		word.addChild(root);
		wordRule.addChild(word);
		molecule.addChild(wordRule);
		processComponents(molecule);
		assertEquals(2, root.getChildCount());
		Element multiplier = root.getChild(0);
		assertEquals(MULTIPLIER_EL, multiplier.getName());
		assertEquals("2", multiplier.getAttributeValue(VALUE_ATR));
		assertEquals("2", multiplier.getValue());
		Element updatedGroup = root.getChild(1);
		assertEquals("hcl", updatedGroup.getValue());
	}

	private void processComponents(Element parse) throws ComponentGenerationException {
		new ComponentGenerator(new BuildState(new NameToStructureConfig())).processParse(parse);
	}
	
}
