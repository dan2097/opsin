package uk.ac.cam.ch.wwmm.opsin;
import java.util.List;
import nu.xom.Element;
import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertFalse;

import org.junit.AfterClass;
import org.junit.BeforeClass;


import org.junit.Test;

public class ParserTest {
	private static Parser parser;
	private static NameToStructureConfig config;

	@BeforeClass
	public static void setUp(){
		parser = new Parser();
		config = NameToStructureConfig.getDefaultConfigInstance();
	}
	
	@AfterClass
	public static void cleanUp(){
		parser = null;
		config = null;
	}

	@Test(expected=ParsingException.class)
	public void testParseThrowsWhenNameIsUninterpretable() throws ParsingException {
		parser.parse(config, "chunky bacon");
	}

	@Test
	public void testParseUninvertsCASNomenclature() throws ParsingException {
		List<Element> parse = parser.parse(config, "Piperidine, 1-(1-oxopropyl)-");

		assertFalse(parse.isEmpty());
	}

	@Test
	public void testParseReturnsOneWordRuleForEachMixtureComponent() throws ParsingException {
		List<Element> parse = parser.parse(config, "benzene; ethane");

		assertEquals(2, parse.get(0).getChildElements(XmlDeclarations.WORDRULE_EL).size());
	}

	@Test(expected=ParsingException.class)
	public void testParseThrowsWhenNameIsSubstituentOnly() throws ParsingException {
		parser.parse(config, "chloro");
	}

	@Test
	public void testConvertStringToComponentRatios1() throws ParsingException {
		String ratio = "(1:2)";
		Integer[] componentRatios = Parser.processStoichometryIndication(ratio);
		assertEquals(2, componentRatios.length);
		for (int i = 0; i < componentRatios.length; i++) {
			if (i==0){
				assertEquals(1,(int) componentRatios[i]);
			}
			if (i==1){
				assertEquals(2,(int) componentRatios[i]);
			}
		}
	}
	
	@Test
	public void testConvertStringToComponentRatios2() throws ParsingException {
		String ratio = "[1/1/2]";
		Integer[] componentRatios = Parser.processStoichometryIndication(ratio);
		assertEquals(3, componentRatios.length);
		for (int i = 0; i < componentRatios.length; i++) {
			if (i==0){
				assertEquals(1,(int) componentRatios[i]);
			}
			if (i==1){
				assertEquals(1,(int) componentRatios[i]);
			}
			if (i==2){
				assertEquals(2,(int) componentRatios[i]);
			}
		}
	}
	
	@Test
	public void testConvertStringToComponentRatios3() throws ParsingException {
		String ratio = "(1:2:?)";
		Integer[] componentRatios = Parser.processStoichometryIndication(ratio);
		assertEquals(3, componentRatios.length);
		for (int i = 0; i < componentRatios.length; i++) {
			if (i==0){
				assertEquals(1,(int) componentRatios[i]);
			}
			if (i==1){
				assertEquals(2,(int) componentRatios[i]);
			}
			if (i==2){
				assertEquals(1,(int) componentRatios[i]);
			}
		}
	}
}
