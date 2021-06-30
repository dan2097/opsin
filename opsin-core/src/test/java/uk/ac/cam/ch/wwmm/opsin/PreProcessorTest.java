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
}
