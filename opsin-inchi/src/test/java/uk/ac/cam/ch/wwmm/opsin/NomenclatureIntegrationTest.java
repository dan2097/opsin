package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

import org.apache.commons.io.IOUtils;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

public class NomenclatureIntegrationTest {
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
	public void testRadicals() throws IOException{
		NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
		n2sConfig.setAllowRadicals(true);
		String file = "radicals.txt";
		checkNamesAgainstInChIs(file, n2sConfig);
	}
	
	@Test
	public void testAcetals() throws IOException{
		NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
		n2sConfig.setAllowRadicals(true);
		String file = "acetals.txt";
		checkNamesAgainstInChIs(file, n2sConfig);
	}
	
	@Test
	public void testAlcoholEsters() throws IOException{
		NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
		n2sConfig.setAllowRadicals(true);
		String file = "alcoholEsters.txt";
		checkNamesAgainstInChIs(file, n2sConfig);
	}
	
	@Test
	public void testCarbohydrates() throws IOException{
		NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
		n2sConfig.setAllowRadicals(true);
		String file = "carbohydrates.txt";
		checkNamesAgainstInChIs(file, n2sConfig);
	}

	@Test
	public void testChargeBalancing() throws IOException{
		NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
		n2sConfig.setAllowRadicals(true);
		String file = "chargeBalancing.txt";
		checkNamesAgainstInChIs(file, n2sConfig);
	}
	
	@Test
	public void testConjunctiveNomenclature() throws IOException{
		NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
		n2sConfig.setAllowRadicals(true);
		String file = "conjunctiveNomenclature.txt";
		checkNamesAgainstInChIs(file, n2sConfig);
	}
	
	@Test
	public void testCyclicSuffixes() throws IOException{
		NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
		n2sConfig.setAllowRadicals(true);
		String file = "cyclicSuffixes.txt";
		checkNamesAgainstInChIs(file, n2sConfig);
	}
	
	@Test
	public void testEpoxyLike() throws IOException{
		NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
		n2sConfig.setAllowRadicals(true);
		String file = "epoxyLike.txt";
		checkNamesAgainstInChIs(file, n2sConfig);
	}
	
	@Test
	public void testFunctionalReplacement() throws IOException{
		NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
		n2sConfig.setAllowRadicals(true);
		String file = "functionalReplacement.txt";
		checkNamesAgainstInChIs(file, n2sConfig);
	}
	
	@Test
	public void testIsotopes() throws IOException{
		NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
		n2sConfig.setAllowRadicals(true);
		String file = "isotopes.txt";
		checkNamesAgainstInChIs(file, n2sConfig);
	}
	
	@Test
	public void testAdditiveNomenclature() throws IOException{
		NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
		n2sConfig.setAllowRadicals(true);
		String file = "additiveNomenclature.txt";
		checkNamesAgainstInChIs(file, n2sConfig);
	}
	
	@Test
	public void testMultiplicativeNomenclature() throws IOException{
		NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
		n2sConfig.setAllowRadicals(true);
		String file = "multiplicativeNomenclature.txt";
		checkNamesAgainstInChIs(file, n2sConfig);
	}
	
	@Test
	public void testOmittedSpaces() throws IOException{
		NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
		n2sConfig.setAllowRadicals(true);
		String file = "omittedSpaces.txt";
		checkNamesAgainstInChIs(file, n2sConfig);
	}
	
	@Test
	public void testFunctionalClassNomenclature() throws IOException{
		NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
		n2sConfig.setAllowRadicals(true);
		String file = "functionalClasses.txt";
		checkNamesAgainstInChIs(file, n2sConfig);
	}
	
	@Test
	public void testFusedRingNomenclature() throws IOException{
		NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
		n2sConfig.setAllowRadicals(true);
		String file = "fusedRings.txt";
		checkNamesAgainstInChIs(file, n2sConfig);
	}
	
	@Test
	public void testInorganicNomenclature() throws IOException{
		NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
		n2sConfig.setAllowRadicals(true);
		String file = "inorganics.txt";
		checkNamesAgainstInChIs(file, n2sConfig);
	}
	
	@Test
	public void testIonNomenclature() throws IOException{
		NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
		n2sConfig.setAllowRadicals(true);
		String file = "ions.txt";
		checkNamesAgainstInChIs(file, n2sConfig);
	}
	
	@Test
	public void testSpiroNomenclature() throws IOException{
		NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
		n2sConfig.setAllowRadicals(true);
		String file = "spiro.txt";
		checkNamesAgainstInChIs(file, n2sConfig);
	}
	
	@Test
	public void testOrganoMetallics() throws IOException{
		NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
		n2sConfig.setAllowRadicals(true);
		String file = "organometallics.txt";
		checkNamesAgainstInChIs(file, n2sConfig);
	}

	@Test
	public void testImplicitBracketting() throws IOException{
		NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
		n2sConfig.setAllowRadicals(true);
		String file = "implicitBracketting.txt";
		checkNamesAgainstInChIs(file, n2sConfig);
	}
	
	@Test
	public void testStereochemistry() throws IOException{
		NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
		n2sConfig.setAllowRadicals(true);
		String file = "stereochemistry.txt";
		checkNamesAgainstInChIs(file, n2sConfig);
	}

	@Test
	public void testMiscellany() throws IOException{
		NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
		n2sConfig.setAllowRadicals(true);
		String file = "miscellany.txt";
		checkNamesAgainstInChIs(file, n2sConfig);
	}

	private void checkNamesAgainstInChIs(String file, NameToStructureConfig n2sConfig) throws IOException{
		BufferedReader input = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream(file), "UTF-8"));
		try {
			String line = null;
			while ((line = input.readLine()) != null) {
				if(line.startsWith("//")){
					continue;
				}
				String[] lineArray = line.split("\t");
				String inchi = NameToInchi.convertResultToInChI(n2s.parseChemicalName(lineArray[0], n2sConfig));
				if (inchi!=null){
					String opsinInchi = InchiPruner.mainChargeAndStereochemistryLayers(inchi);
					String referenceInchi = InchiPruner.mainChargeAndStereochemistryLayers(lineArray[1]);
	
					if (!opsinInchi.equals(referenceInchi)){
						fail(lineArray[0] +" was misinterpreted as: " + inchi);
					}
				} else {
					fail(lineArray[0] +" was uninterpretable");
				}
			}
		} finally {
			IOUtils.closeQuietly(input);
		}
	}
}
