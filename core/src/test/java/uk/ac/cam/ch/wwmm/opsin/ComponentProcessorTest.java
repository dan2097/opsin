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
		group1.addAttribute(new Attribute(SUBTYPE_ATR, SIMPLEGROUP_SUBTYPE_VAL));
		substituent2.appendChild(group1);
		word.appendChild(substituent2);
		Element root = new Element(ROOT_EL);
		word.appendChild(root);
		Element group2 = new Element(GROUP_EL);
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
}
