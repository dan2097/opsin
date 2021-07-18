package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvFileSource;

public class AmbiguityDetectionTest {
	
	private static NameToStructure n2s;

	@BeforeAll
	public static void setUp() {
		n2s = NameToStructure.getInstance();
	}
	
	@AfterAll
	public static void cleanUp(){
		n2s = null;
	}

	@ParameterizedTest
	@CsvFileSource(resources ="ambiguous.txt",  delimiter='\t')
	public void testNamesThatShouldBeDetectedAsAmbiguous(String ambiguousName) {
		assertTrue(n2s.parseChemicalName(ambiguousName).nameAppearsToBeAmbiguous(), ambiguousName + " should be considered ambiguous");
	}
	
	@ParameterizedTest
	@CsvFileSource(resources ="unambiguous.txt",  delimiter='\t')
	public void testUnAmbiguousCounterExamples(String unambiguousName) {
		assertFalse(n2s.parseChemicalName(unambiguousName).nameAppearsToBeAmbiguous(), unambiguousName + " should be considered unambiguous");
	}
	
}
