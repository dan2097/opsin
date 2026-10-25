package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import org.junit.jupiter.api.Test;

public class PreProcessorTest {

	@Test()
	public void testPreProcessBlankThrows() {
		assertThrows(PreProcessingException.class, () -> {
		PreProcessor.preProcess("");
		});
	}

	@Test
	public void testPreProcessConvertsDollarA() throws PreProcessingException {
		assertEquals("alpha-bromo", PreProcessor.preProcess("$a-bromo"), "Convert dollar-a");
	}

	@Test
	public void testPreProcessConvertsDollarB() throws PreProcessingException {
		assertEquals("beta-bromo", PreProcessor.preProcess("$b-bromo"), "Convert dollar-b");
	}

	@Test
	public void testPreProcessConvertsDollarG() throws PreProcessingException {
		assertEquals("gamma-bromo", PreProcessor.preProcess("$g-bromo"), "Convert dollar-g");
	}

	@Test
	public void testPreProcessConvertsDollarD() throws PreProcessingException {
		assertEquals("delta-bromo", PreProcessor.preProcess("$d-bromo"), "Convert dollar-d");
	}

	@Test
	public void testPreProcessConvertsDollarE() throws PreProcessingException {
		assertEquals("epsilon-bromo", PreProcessor.preProcess("$e-bromo"), "Convert dollar-e");
	}

	@Test
	public void testPreProcessConvertsDollarL() throws PreProcessingException {
		assertEquals("lambda-bromo", PreProcessor.preProcess("$l-bromo"), "Convert dollar-l");
	}

	@Test
	public void testPreProcessConvertsGreekLetterToWord() throws PreProcessingException {
		assertEquals("alpha-bromo", PreProcessor.preProcess("\u03b1-bromo"), "Convert greek to word");
	}

	@Test
	public void testPreProcessConvertsSulphToSulf() throws PreProcessingException {
		assertEquals("sulfur dioxide", PreProcessor.preProcess("sulphur dioxide"), "Converts 'sulph' to 'sulph'");
	}
	
	@Test
	public void testRemovalOfDotsFromGreekWords1() throws PreProcessingException {
		assertEquals("alpha-methyl-toluene", PreProcessor.preProcess(".alpha.-methyl-toluene"), "Converts '.alpha.' to 'alpha'");
	}
	
	@Test
	public void testRemovalOfDotsFromGreekWords2() throws PreProcessingException {
		assertEquals("alphabetaeta", PreProcessor.preProcess(".alpha..beta..eta."));
	}

	@Test
	public void testHtmlGreeks() throws PreProcessingException {
		assertEquals("alpha-methyl-toluene", PreProcessor.preProcess("&alpha;-methyl-toluene"));
		assertEquals("beta-methyl-styrene", PreProcessor.preProcess("&BETA;-methyl-styrene"));
	}

	@Test
	public void testSuperscriptDigits() throws PreProcessingException {
		assertEquals("1^4-ol", PreProcessor.preProcess("1\u2074-ol"));
		assertEquals("4^10,5^2-diol", PreProcessor.preProcess("4\u00B9\u2070,5\u00B2-diol"));
		assertEquals("3^1-yl", PreProcessor.preProcess("3\u00B9-yl"));
	}

	@Test
	public void testHtmlSuperscripts() throws PreProcessingException {
		assertEquals("1^4-ol", PreProcessor.preProcess("1<sup>4</sup>-ol"));
		assertEquals("1^4-ol", PreProcessor.preProcess("1<SUP>4</SUP>-ol"));
	}

	@Test
	public void testEmDashSeparatesComponents() throws PreProcessingException {
		for (String separator : new String[] {"\u2014", " \u2014", "\u2014 ", " \u2014 ", " \t\u2014\u00A0 "}) {
			assertEquals("benzene - 2-aminopyridine", PreProcessor.preProcess("benzene" + separator + "2-aminopyridine"));
		}
	}

	@Test
	public void testWhitespaceNormalization() throws PreProcessingException {
		assertEquals("sodium chloride", PreProcessor.preProcess("sodium  \t\n\u000B\f\r \u00A0\u3000chloride"));
		assertEquals("sodium chloride", PreProcessor.preProcess("sodium \u00AD chloride"));
		assertEquals(" ", StringTools.convertNonAsciiAndNormaliseRepresentation("\t \u00A0"));
	}

	@Test
	public void testTerminalEmDashRemainsHyphen() throws PreProcessingException {
		assertEquals("ethyl-", PreProcessor.preProcess("ethyl\u2014"));
		assertEquals("ethyl-", PreProcessor.preProcess("ethyl\u2014 \t"));
	}

	@Test
	public void testOtherDashesRemainHyphens() throws PreProcessingException {
		for (char dash : new char[] {'\u2010', '\u2011', '\u2012', '\u2013', '\u2015', '\u2212'}) {
			assertEquals("2-methylpyridine", PreProcessor.preProcess("2" + dash + "methylpyridine"));
		}
	}
}
