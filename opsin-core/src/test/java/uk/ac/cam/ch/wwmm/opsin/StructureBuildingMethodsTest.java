package uk.ac.cam.ch.wwmm.opsin;

import org.junit.Test;

import nu.xom.Attribute;
import nu.xom.Element;
import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;
import static junit.framework.Assert.assertEquals;

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
}
