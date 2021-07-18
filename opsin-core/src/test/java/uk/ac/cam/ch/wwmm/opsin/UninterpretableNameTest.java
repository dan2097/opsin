package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import uk.ac.cam.ch.wwmm.opsin.OpsinResult.OPSIN_RESULT_STATUS;

public class UninterpretableNameTest {
	
	private static NameToStructure n2s;

	@BeforeAll
	public static void setUp() {
		n2s = NameToStructure.getInstance();
	}
	
	@AfterAll
	public static void cleanUp(){
		n2s = null;
	}

	@Test
	public void testNamesThatShoudlBeUninterpretable() throws IOException{
		checkNames("uninterpretable.txt");
	}

	private void checkNames(String file) throws IOException{
		try (BufferedReader input = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream(file)))) {
			String line = null;
			while ((line = input.readLine()) != null) {
				if(line.startsWith("//")){
					continue;
				}
				OpsinResult result = n2s.parseChemicalName(line);
				assertEquals(OPSIN_RESULT_STATUS.FAILURE, result.getStatus(), line + " gave unexpected result");
			}
		}
	}
	
}
