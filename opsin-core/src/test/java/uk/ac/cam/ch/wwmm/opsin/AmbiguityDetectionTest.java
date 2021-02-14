package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.Assert.assertEquals;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

public class AmbiguityDetectionTest {
	
	private static NameToStructure n2s;

	@BeforeClass
	public static void setUp() {
		n2s = NameToStructure.getInstance();
	}
	
	@AfterClass
	public static void cleanUp(){
		n2s = null;
	}

	@Test
	public void testNamesThatShouldBeDetectedAsAmbiguous() throws IOException{
		checkNames("ambiguous.txt", true);
	}
	
	@Test
	public void testUnAmbiguousCounterExamples() throws IOException{
		checkNames("unambiguous.txt", false);
	}
	
	private void checkNames(String file, boolean isAmbiguous) throws IOException{
		try(BufferedReader input = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream(file)))) {
			String line = null;
			while ((line = input.readLine()) != null) {
				if(line.startsWith("//")){
					continue;
				}
				OpsinResult result = n2s.parseChemicalName(line);
				assertEquals(line + " gave unexpected result", isAmbiguous, result.nameAppearsToBeAmbiguous());
			}
		}
	}
	
}
