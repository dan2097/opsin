package uk.ac.cam.ch.wwmm.opsin;
import java.util.List;
import nu.xom.Element;
import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertFalse;
import static org.junit.Assert.fail;
import org.junit.Before;


import org.junit.Test;

public class ParserTest {
	private Parser parser;
	private NameToStructureConfig config;

	@Before
	public void setUp() throws Exception {
		parser = new Parser();
		config = NameToStructureConfig.getDefaultConfigInstance();
	}

	@Test
	public void testParseThrowsWhenNameIsUninterpretable() {
		try {
			parser.parse(config, "chunky bacon");
			fail("Should throw ParsingException");
		} catch (ParsingException e) {
			// no-op
		}
	}

	@Test
	public void testParseUninvertsCASNomenclature() throws ParsingException {
		List<Element> parse = parser.parse(config, "Piperidine, 1-(1-oxopropyl)-");

		assertFalse(parse.isEmpty());
	}

	@Test
	public void testParseReturnsOneWordRulesForEachMixtureComponent() throws ParsingException {
		List<Element> parse = parser.parse(config, "benzene; ethane");

		assertEquals(2, parse.get(0).getChildElements("wordRule").size());
	}

	@Test
	public void testParseThrowsWhenNameIsSubstituentOnly() {
		try {
			parser.parse(config, "chloro");
		} catch (ParsingException e) {
			// no-op
		}
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
