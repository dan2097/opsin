package uk.ac.cam.ch.wwmm.opsin;


import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.io.IOException;
import java.util.List;

import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;


public class ParserTest {
	private static Parser parser;
	private static NameToStructureConfig config;

	@BeforeAll
	public static void setUp() throws IOException{
		parser = new Parser();
		config = NameToStructureConfig.getDefaultConfigInstance();
	}
	
	@AfterAll
	public static void cleanUp(){
		parser = null;
		config = null;
	}

	@Test()
	public void testParseThrowsWhenNameIsUninterpretable() throws ParsingException {
		assertThrows(ParsingException.class, () -> {
			parser.parse(config, "chunky bacon");
		});
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

	@Test()
	public void testParseThrowsWhenNameIsSubstituentOnly() {
		assertThrows(ParsingException.class, () -> {
			parser.parse(config, "chloro");
		});
	}
	
	@Test()
	public void testNoParseForOneComponentSalt() {
		assertThrows(ParsingException.class, () -> {
			parser.parse(config, "pyridine salt");
		});
	}

	@Test
	public void testConvertStringToComponentRatios1() throws ParsingException {
		String ratio = "(1:2)";
		Integer[] componentRatios = Parser.processStoichiometryIndication(ratio);
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
		Integer[] componentRatios = Parser.processStoichiometryIndication(ratio);
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
		Integer[] componentRatios = Parser.processStoichiometryIndication(ratio);
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
