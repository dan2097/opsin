package uk.ac.cam.ch.wwmm.opsin;

import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;
import static junit.framework.Assert.*;

import org.junit.Test;

import nu.xom.Attribute;
import nu.xom.Element;

public class ComponentGeneration_AmbiguitiesAndIrregularities {
	
	@Test
	public void testCorrectlyTokenisedAlkane(){
		Element substituent = new Element(SUBSTITUENT_EL);
		Element alkaneComponent1 = new Element(ALKANESTEMCOMPONENT);
		alkaneComponent1.addAttribute(new Attribute(VALUE_ATR, "4"));
		Element alkaneComponent2 = new Element(ALKANESTEMCOMPONENT);
		alkaneComponent2.addAttribute(new Attribute(VALUE_ATR, "10"));
		substituent.appendChild(alkaneComponent1);
		substituent.appendChild(alkaneComponent2);
		try{
			ComponentGenerator.resolveAmbiguities(substituent);
		}
		catch (ComponentGenerationException e) {
			fail("alkane was well formed, exception should not be thrown");
		}
	}
	
	@Test
	public void testCorrectlyTokenisedAlkane2(){
		Element substituent = new Element(SUBSTITUENT_EL);
		Element multiplier = new Element(MULTIPLIER_EL);
		multiplier.addAttribute(new Attribute(TYPE_ATR, BASIC_TYPE_VAL));
		multiplier.addAttribute(new Attribute(VALUE_ATR, "2"));
		Element alkaneComponent = new Element(ALKANESTEMCOMPONENT);
		alkaneComponent.addAttribute(new Attribute(VALUE_ATR, "10"));
		substituent.appendChild(multiplier);
		substituent.appendChild(alkaneComponent);
		try{
			ComponentGenerator.resolveAmbiguities(substituent);
		}
		catch (ComponentGenerationException e) {
			fail("alkane was well formed, exception should not be thrown");
		}
	}
	
	
	@Test
	public void testCorrectlyTokenisedAlkane3(){//unambiguously 6 hexanes
		Element substituent = new Element(SUBSTITUENT_EL);
		Element multiplier = new Element(MULTIPLIER_EL);
		multiplier.addAttribute(new Attribute(TYPE_ATR, BASIC_TYPE_VAL));
		multiplier.addAttribute(new Attribute(VALUE_ATR, "6"));
		Element alkaneComponent = new Element(ALKANESTEMCOMPONENT);
		alkaneComponent.addAttribute(new Attribute(VALUE_ATR, "6"));
		substituent.appendChild(multiplier);
		substituent.appendChild(alkaneComponent);
		try{
			ComponentGenerator.resolveAmbiguities(substituent);
		}
		catch (ComponentGenerationException e) {
			fail("alkane was well formed, exception should not be thrown");
		}
	}
	
	@Test(expected=ComponentGenerationException.class)//tetradec is 14 not 4 x 10
	public void testMisTokenisedAlkane() throws ComponentGenerationException{
		Element substituent = new Element(SUBSTITUENT_EL);
		Element erroneousMultiplier = new Element(MULTIPLIER_EL);
		erroneousMultiplier.addAttribute(new Attribute(TYPE_ATR, BASIC_TYPE_VAL));
		erroneousMultiplier.addAttribute(new Attribute(VALUE_ATR, "4"));
		Element alkaneComponent2 = new Element(ALKANESTEMCOMPONENT);
		alkaneComponent2.addAttribute(new Attribute(VALUE_ATR, "10"));
		substituent.appendChild(erroneousMultiplier);
		substituent.appendChild(alkaneComponent2);
		ComponentGenerator.resolveAmbiguities(substituent);
	}

	@Test
	public void testLocantsIndicatingTokenizationIsCorrect(){//should be a group multiplier formally
		Element substituent = new Element(SUBSTITUENT_EL);
		Element locant = new Element(LOCANT_EL);
		locant.appendChild("1,2,3,4");
		substituent.appendChild(locant);
		Element multiplier = new Element(MULTIPLIER_EL);
		multiplier.addAttribute(new Attribute(TYPE_ATR, BASIC_TYPE_VAL));
		multiplier.addAttribute(new Attribute(VALUE_ATR, "4"));
		Element alkaneComponent = new Element(ALKANESTEMCOMPONENT);
		alkaneComponent.addAttribute(new Attribute(VALUE_ATR, "10"));
		substituent.appendChild(multiplier);
		substituent.appendChild(alkaneComponent);
		try{
			ComponentGenerator.resolveAmbiguities(substituent);
		}
		catch (ComponentGenerationException e) {
			fail("alkane was well formed, exception should not be thrown");
		}
	}
	
	@Test(expected=ComponentGenerationException.class)//tetradec is 14 not 4 x 10
	public void testLocantsIndicatingTokenizationIsIncorrect() throws ComponentGenerationException{
		Element substituent = new Element(SUBSTITUENT_EL);
		Element locant = new Element(LOCANT_EL);
		locant.appendChild("1");
		substituent.appendChild(locant);
		Element erroneousMultiplier = new Element(MULTIPLIER_EL);
		erroneousMultiplier.addAttribute(new Attribute(TYPE_ATR, BASIC_TYPE_VAL));
		erroneousMultiplier.addAttribute(new Attribute(VALUE_ATR, "4"));
		Element alkaneComponent = new Element(ALKANESTEMCOMPONENT);
		alkaneComponent.addAttribute(new Attribute(VALUE_ATR, "10"));
		substituent.appendChild(erroneousMultiplier);
		substituent.appendChild(alkaneComponent);
		ComponentGenerator.resolveAmbiguities(substituent);
	}
	
	
	@Test(expected=ComponentGenerationException.class)
	public void testTetraphenShouldBeTetra_Phen1() throws ComponentGenerationException{//tetraphenyl
		Element substituent = new Element(SUBSTITUENT_EL);
		Element multiplier = new Element(MULTIPLIER_EL);
		multiplier.addAttribute(new Attribute(TYPE_ATR, BASIC_TYPE_VAL));
		multiplier.addAttribute(new Attribute(VALUE_ATR, "4"));
		Element phen = new Element(HYDROCARBONFUSEDRINGSYSTEM_EL);
		phen.appendChild("phen");
		Element yl = new Element(SUFFIX_EL);
		yl.appendChild("yl");
		substituent.appendChild(multiplier);
		substituent.appendChild(phen);
		substituent.appendChild(yl);
		ComponentGenerator.resolveAmbiguities(substituent);
	}
	
	@Test(expected=ComponentGenerationException.class)
	public void testTetraphenShouldBeTetra_Phen2() throws ComponentGenerationException{//tetraphenoxy
		Element substituent = new Element(SUBSTITUENT_EL);
		Element multiplier = new Element(MULTIPLIER_EL);
		multiplier.addAttribute(new Attribute(TYPE_ATR, BASIC_TYPE_VAL));
		multiplier.addAttribute(new Attribute(VALUE_ATR, "4"));
		Element phen = new Element(HYDROCARBONFUSEDRINGSYSTEM_EL);
		phen.appendChild("phen");
		Element yl = new Element(SUFFIX_EL);
		yl.appendChild("oxy");
		substituent.appendChild(multiplier);
		substituent.appendChild(phen);
		substituent.appendChild(yl);
		ComponentGenerator.resolveAmbiguities(substituent);
	}
	
	@Test
	public void testTetraphenShouldBeTetraphen1(){//tetrapheneyl
		Element substituent = new Element(SUBSTITUENT_EL);
		Element multiplier = new Element(MULTIPLIER_EL);
		multiplier.addAttribute(new Attribute(TYPE_ATR, BASIC_TYPE_VAL));
		multiplier.addAttribute(new Attribute(VALUE_ATR, "4"));
		Element phen = new Element(HYDROCARBONFUSEDRINGSYSTEM_EL);
		phen.addAttribute(new Attribute(SUBSEQUENTUNSEMANTICTOKEN_EL, "e"));
		phen.appendChild("phen");
		Element yl = new Element(SUFFIX_EL);
		yl.appendChild("yl");
		substituent.appendChild(multiplier);
		substituent.appendChild(phen);
		substituent.appendChild(yl);
		try{
			ComponentGenerator.resolveAmbiguities(substituent);
		}
		catch (ComponentGenerationException e) {
			fail("tetraphene was the intended interpretation");
		}
	}
	
	@Test
	public void testTetraphenShouldBeTetraphen2(){//tetraphen2yl
		Element substituent = new Element(SUBSTITUENT_EL);
		Element multiplier = new Element(MULTIPLIER_EL);
		multiplier.addAttribute(new Attribute(TYPE_ATR, BASIC_TYPE_VAL));
		multiplier.addAttribute(new Attribute(VALUE_ATR, "4"));
		Element phen = new Element(HYDROCARBONFUSEDRINGSYSTEM_EL);
		phen.appendChild("phen");
		Element locant = new Element(LOCANT_EL);
		locant.appendChild("2");
		Element yl = new Element(SUFFIX_EL);
		yl.appendChild("yl");
		substituent.appendChild(multiplier);
		substituent.appendChild(phen);
		substituent.appendChild(locant);
		substituent.appendChild(yl);
		try{
			ComponentGenerator.resolveAmbiguities(substituent);
		}
		catch (ComponentGenerationException e) {
			fail("tetraphen as in tetraphene was the intended interpretation");
		}
	}
	
	@Test
	public void testTetraphenShouldBeTetraphen3(){//2tetraphenyl
		Element substituent = new Element(SUBSTITUENT_EL);
		Element locant = new Element(LOCANT_EL);
		locant.appendChild("2");
		Element multiplier = new Element(MULTIPLIER_EL);
		multiplier.addAttribute(new Attribute(TYPE_ATR, BASIC_TYPE_VAL));
		multiplier.addAttribute(new Attribute(VALUE_ATR, "4"));
		Element phen = new Element(HYDROCARBONFUSEDRINGSYSTEM_EL);
		phen.appendChild("phen");
		Element yl = new Element(SUFFIX_EL);
		yl.appendChild("yl");
		substituent.appendChild(locant);
		substituent.appendChild(multiplier);
		substituent.appendChild(phen);
		substituent.appendChild(yl);
		try{
			ComponentGenerator.resolveAmbiguities(substituent);
		}
		catch (ComponentGenerationException e) {
			fail("tetraphen as in tetraphene was the intended interpretation");
		}
	}
	
//TODO multiplier oxy tests, fusion vs Hw locants,  and handleGroupIrregularities tests
}
