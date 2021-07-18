package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.fail;
import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

import org.junit.jupiter.api.Test;

public class ComponentGeneration_AmbiguitiesAndIrregularitiesTest {
	
	@Test
	public void testCorrectlyTokenisedAlkane(){
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element alkaneComponent1 = new TokenEl(ALKANESTEMCOMPONENT);
		alkaneComponent1.addAttribute(new Attribute(VALUE_ATR, "4"));
		Element alkaneComponent2 = new TokenEl(ALKANESTEMCOMPONENT);
		alkaneComponent2.addAttribute(new Attribute(VALUE_ATR, "10"));
		substituent.addChild(alkaneComponent1);
		substituent.addChild(alkaneComponent2);
		try{
			ComponentGenerator.resolveAmbiguities(substituent);
		}
		catch (ComponentGenerationException e) {
			fail("alkane was well formed, exception should not be thrown");
		}
	}
	
	@Test
	public void testCorrectlyTokenisedAlkane2(){
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element multiplier = new TokenEl(MULTIPLIER_EL);
		multiplier.addAttribute(new Attribute(TYPE_ATR, BASIC_TYPE_VAL));
		multiplier.addAttribute(new Attribute(VALUE_ATR, "2"));
		Element alkaneComponent = new TokenEl(ALKANESTEMCOMPONENT);
		alkaneComponent.addAttribute(new Attribute(VALUE_ATR, "10"));
		substituent.addChild(multiplier);
		substituent.addChild(alkaneComponent);
		try{
			ComponentGenerator.resolveAmbiguities(substituent);
		}
		catch (ComponentGenerationException e) {
			fail("alkane was well formed, exception should not be thrown");
		}
	}
	
	
	@Test
	public void testCorrectlyTokenisedAlkane3(){//unambiguously 6 hexanes
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element multiplier = new TokenEl(MULTIPLIER_EL);
		multiplier.addAttribute(new Attribute(TYPE_ATR, BASIC_TYPE_VAL));
		multiplier.addAttribute(new Attribute(VALUE_ATR, "6"));
		Element alkaneComponent = new TokenEl(ALKANESTEMCOMPONENT);
		alkaneComponent.addAttribute(new Attribute(VALUE_ATR, "6"));
		substituent.addChild(multiplier);
		substituent.addChild(alkaneComponent);
		try{
			ComponentGenerator.resolveAmbiguities(substituent);
		}
		catch (ComponentGenerationException e) {
			fail("alkane was well formed, exception should not be thrown");
		}
	}
	
	@Test() // tetradec is 14 not 4 x 10
	public void testMisTokenisedAlkane() {

		assertThrows(ComponentGenerationException.class, () -> {
			Element substituent = new GroupingEl(SUBSTITUENT_EL);
			Element erroneousMultiplier = new TokenEl(MULTIPLIER_EL);
			erroneousMultiplier.addAttribute(new Attribute(TYPE_ATR, BASIC_TYPE_VAL));
			erroneousMultiplier.addAttribute(new Attribute(VALUE_ATR, "4"));
			Element alkaneComponent2 = new TokenEl(ALKANESTEMCOMPONENT);
			alkaneComponent2.addAttribute(new Attribute(VALUE_ATR, "10"));
			substituent.addChild(erroneousMultiplier);
			substituent.addChild(alkaneComponent2);
			ComponentGenerator.resolveAmbiguities(substituent);
		});
	}

	@Test
	public void testLocantsIndicatingTokenizationIsCorrect(){//should be a group multiplier formally
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element locant = new TokenEl(LOCANT_EL, "1,2,3,4");
		substituent.addChild(locant);
		Element multiplier = new TokenEl(MULTIPLIER_EL);
		multiplier.addAttribute(new Attribute(TYPE_ATR, BASIC_TYPE_VAL));
		multiplier.addAttribute(new Attribute(VALUE_ATR, "4"));
		Element alkaneComponent = new TokenEl(ALKANESTEMCOMPONENT);
		alkaneComponent.addAttribute(new Attribute(VALUE_ATR, "10"));
		substituent.addChild(multiplier);
		substituent.addChild(alkaneComponent);
		try{
			ComponentGenerator.resolveAmbiguities(substituent);
		}
		catch (ComponentGenerationException e) {
			fail("alkane was well formed, exception should not be thrown");
		}
	}
	
	@Test() // tetradec is 14 not 4 x 10
	public void testLocantsIndicatingTokenizationIsIncorrect() {
		assertThrows(ComponentGenerationException.class, () -> {
			Element substituent = new GroupingEl(SUBSTITUENT_EL);
			Element locant = new TokenEl(LOCANT_EL, "1");
			substituent.addChild(locant);
			Element erroneousMultiplier = new TokenEl(MULTIPLIER_EL);
			erroneousMultiplier.addAttribute(new Attribute(TYPE_ATR, BASIC_TYPE_VAL));
			erroneousMultiplier.addAttribute(new Attribute(VALUE_ATR, "4"));
			Element alkaneComponent = new TokenEl(ALKANESTEMCOMPONENT);
			alkaneComponent.addAttribute(new Attribute(VALUE_ATR, "10"));
			substituent.addChild(erroneousMultiplier);
			substituent.addChild(alkaneComponent);
			ComponentGenerator.resolveAmbiguities(substituent);
		});
	}
	
	
	@Test()
	public void testTetraphenShouldBeTetra_Phen1() {// tetraphenyl
		assertThrows(ComponentGenerationException.class, () -> {
			Element substituent = new GroupingEl(SUBSTITUENT_EL);
			Element multiplier = new TokenEl(MULTIPLIER_EL);
			multiplier.addAttribute(new Attribute(TYPE_ATR, BASIC_TYPE_VAL));
			multiplier.addAttribute(new Attribute(VALUE_ATR, "4"));
			Element phen = new TokenEl(HYDROCARBONFUSEDRINGSYSTEM_EL, "phen");
			Element yl = new TokenEl(SUFFIX_EL, "yl");
			substituent.addChild(multiplier);
			substituent.addChild(phen);
			substituent.addChild(yl);
			ComponentGenerator.resolveAmbiguities(substituent);
		});
	}
	
	@Test()
	public void testTetraphenShouldBeTetra_Phen2() {// tetraphenoxy
		assertThrows(ComponentGenerationException.class, () -> {
			Element substituent = new GroupingEl(SUBSTITUENT_EL);
			Element multiplier = new TokenEl(MULTIPLIER_EL);
			multiplier.addAttribute(new Attribute(TYPE_ATR, BASIC_TYPE_VAL));
			multiplier.addAttribute(new Attribute(VALUE_ATR, "4"));
			Element phen = new TokenEl(HYDROCARBONFUSEDRINGSYSTEM_EL, "phen");
			Element yl = new TokenEl(SUFFIX_EL, "oxy");
			substituent.addChild(multiplier);
			substituent.addChild(phen);
			substituent.addChild(yl);
			ComponentGenerator.resolveAmbiguities(substituent);
		});
	}
	
	@Test
	public void testTetraphenShouldBeTetraphen1(){//tetrapheneyl
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element multiplier = new TokenEl(MULTIPLIER_EL);
		multiplier.addAttribute(new Attribute(TYPE_ATR, BASIC_TYPE_VAL));
		multiplier.addAttribute(new Attribute(VALUE_ATR, "4"));
		Element phen = new TokenEl(HYDROCARBONFUSEDRINGSYSTEM_EL, "phen");
		phen.addAttribute(new Attribute(SUBSEQUENTUNSEMANTICTOKEN_ATR, "e"));
		Element yl = new TokenEl(SUFFIX_EL, "yl");
		substituent.addChild(multiplier);
		substituent.addChild(phen);
		substituent.addChild(yl);
		try{
			ComponentGenerator.resolveAmbiguities(substituent);
		}
		catch (ComponentGenerationException e) {
			fail("tetraphene was the intended interpretation");
		}
	}
	
	@Test
	public void testTetraphenShouldBeTetraphen2(){//tetraphen2yl
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element multiplier = new TokenEl(MULTIPLIER_EL);
		multiplier.addAttribute(new Attribute(TYPE_ATR, BASIC_TYPE_VAL));
		multiplier.addAttribute(new Attribute(VALUE_ATR, "4"));
		Element phen = new TokenEl(HYDROCARBONFUSEDRINGSYSTEM_EL, "phen");
		Element locant = new TokenEl(LOCANT_EL, "2");
		Element yl = new TokenEl(SUFFIX_EL, "yl");
		substituent.addChild(multiplier);
		substituent.addChild(phen);
		substituent.addChild(locant);
		substituent.addChild(yl);
		try{
			ComponentGenerator.resolveAmbiguities(substituent);
		}
		catch (ComponentGenerationException e) {
			fail("tetraphen as in tetraphene was the intended interpretation");
		}
	}
	
	@Test
	public void testTetraphenShouldBeTetraphen3(){//2tetraphenyl
		Element substituent = new GroupingEl(SUBSTITUENT_EL);
		Element locant = new TokenEl(LOCANT_EL, "2");
		Element multiplier = new TokenEl(MULTIPLIER_EL);
		multiplier.addAttribute(new Attribute(TYPE_ATR, BASIC_TYPE_VAL));
		multiplier.addAttribute(new Attribute(VALUE_ATR, "4"));
		Element phen = new TokenEl(HYDROCARBONFUSEDRINGSYSTEM_EL, "phen");
		Element yl = new TokenEl(SUFFIX_EL, "yl");
		substituent.addChild(locant);
		substituent.addChild(multiplier);
		substituent.addChild(phen);
		substituent.addChild(yl);
		try{
			ComponentGenerator.resolveAmbiguities(substituent);
		}
		catch (ComponentGenerationException e) {
			fail("tetraphen as in tetraphene was the intended interpretation");
		}
	}
	
//TODO multiplier oxy tests, fusion vs Hw locants,  and handleGroupIrregularities tests
}
