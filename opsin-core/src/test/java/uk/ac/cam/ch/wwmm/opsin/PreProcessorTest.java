package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.Assert.assertEquals;
import org.junit.Test;

public class PreProcessorTest {

	@Test(expected=PreProcessingException.class)
	public void testPreProcessBlankThrows() throws PreProcessingException {
		PreProcessor.preProcess("");
	}

	@Test
	public void testPreProcessConvertsDollarA() throws PreProcessingException {
		assertEquals("Convert dollar-a", "alpha-bromo", PreProcessor.preProcess("$a-bromo"));
	}

	@Test
	public void testPreProcessConvertsDollarB() throws PreProcessingException {
		assertEquals("Convert dollar-b", "beta-bromo", PreProcessor.preProcess("$b-bromo"));
	}

	@Test
	public void testPreProcessConvertsDollarG() throws PreProcessingException {
		assertEquals("Convert dollar-g", "gamma-bromo", PreProcessor.preProcess("$g-bromo"));
	}

	@Test
	public void testPreProcessConvertsDollarD() throws PreProcessingException {
		assertEquals("Convert dollar-d", "delta-bromo", PreProcessor.preProcess("$d-bromo"));
	}

	@Test
	public void testPreProcessConvertsDollarE() throws PreProcessingException {
		assertEquals("Convert dollar-e", "epsilon-bromo", PreProcessor.preProcess("$e-bromo"));
	}

	@Test
	public void testPreProcessConvertsDollarL() throws PreProcessingException {
		assertEquals("Convert dollar-l", "lambda-bromo", PreProcessor.preProcess("$l-bromo"));
	}

	@Test
	public void testPreProcessConvertsGreekLetterToWord() throws PreProcessingException {
		assertEquals("Convert greek to word", "alpha-bromo", PreProcessor.preProcess("\u03b1-bromo"));
	}

	@Test
	public void testPreProcessConvertsSulphToSulf() throws PreProcessingException {
		assertEquals("Converts 'sulph' to 'sulph'", "sulfur dioxide", PreProcessor.preProcess("sulphur dioxide"));
	}
	
	@Test
	public void testRemovalOfDotsFromGreekWords1() throws PreProcessingException {
		assertEquals("Converts '.alpha.' to 'alpha'", "alpha-methyl-toluene", PreProcessor.preProcess(".alpha.-methyl-toluene"));
	}
	
	@Test
	public void testRemovalOfDotsFromGreekWords2() throws PreProcessingException {
		assertEquals("alphabetaeta", PreProcessor.preProcess(".alpha..beta..eta."));
	}
}
