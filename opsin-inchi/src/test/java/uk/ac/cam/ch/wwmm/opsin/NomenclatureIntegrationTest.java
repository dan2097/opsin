package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.fail;

import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvFileSource;

public class NomenclatureIntegrationTest {
	
	private static NameToStructure n2s;
	private static NameToStructureConfig n2sConfig;

	@BeforeAll
	public static void setUp() {
		n2s = NameToStructure.getInstance();
		n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
		n2sConfig.setAllowRadicals(true);
	}
	
	@AfterAll
	public static void cleanUp(){
		n2s = null;
		n2sConfig = null;
	}

	@ParameterizedTest
	@CsvFileSource(resources = "radicals.txt", delimiter='\t')
	public void testRadicals(String name, String expectedInchi) {
		checkName(name, expectedInchi);
	}
	
	@ParameterizedTest
	@CsvFileSource(resources = "acetals.txt", delimiter='\t')
	public void testAcetals(String name, String expectedInchi) {
		checkName(name, expectedInchi);
	}
	
	@ParameterizedTest
	@CsvFileSource(resources = "alcoholEsters.txt", delimiter='\t')
	public void testAlcoholEsters(String name, String expectedInchi) {
		checkName(name, expectedInchi);
	}
	
	@ParameterizedTest
	@CsvFileSource(resources = "aminoAcids.txt", delimiter='\t')
	public void testAminoAcids(String name, String expectedInchi) {
		checkName(name, expectedInchi);
	}
	
	@ParameterizedTest
	@CsvFileSource(resources = "carbohydrates.txt", delimiter='\t')
	public void testCarbohydrates(String name, String expectedInchi) {
		checkName(name, expectedInchi);
	}

	@ParameterizedTest
	@CsvFileSource(resources = "chargeBalancing.txt", delimiter='\t')
	public void testChargeBalancing(String name, String expectedInchi) {
		checkName(name, expectedInchi);
	}
	
	@ParameterizedTest
	@CsvFileSource(resources = "conjunctiveNomenclature.txt", delimiter='\t')
	public void testConjunctiveNomenclature(String name, String expectedInchi) {
		checkName(name, expectedInchi);
	}
	
	@ParameterizedTest
	@CsvFileSource(resources = "cyclicSuffixes.txt", delimiter='\t')
	public void testCyclicSuffixes(String name, String expectedInchi) {
		checkName(name, expectedInchi);
	}
	
	@ParameterizedTest
	@CsvFileSource(resources = "epoxyLike.txt", delimiter='\t')
	public void testEpoxyLike(String name, String expectedInchi) {
		checkName(name, expectedInchi);
	}
	
	@ParameterizedTest
	@CsvFileSource(resources = "functionalReplacement.txt", delimiter='\t')
	public void testFunctionalReplacement(String name, String expectedInchi) {
		checkName(name, expectedInchi);
	}
	
	@ParameterizedTest
	@CsvFileSource(resources = "isotopes.txt", delimiter='\t')
	public void testIsotopes(String name, String expectedInchi) {
		checkName(name, expectedInchi);
	}
	
	@ParameterizedTest
	@CsvFileSource(resources = "additiveNomenclature.txt", delimiter='\t')
	public void testAdditiveNomenclature(String name, String expectedInchi) {
		checkName(name, expectedInchi);
	}
	
	@ParameterizedTest
	@CsvFileSource(resources = "multiplicativeNomenclature.txt", delimiter='\t')
	public void testMultiplicativeNomenclature(String name, String expectedInchi) {
		checkName(name, expectedInchi);
	}
	
	@ParameterizedTest
	@CsvFileSource(resources = "omittedSpaces.txt", delimiter='\t')
	public void testOmittedSpaces(String name, String expectedInchi) {
		checkName(name, expectedInchi);
	}
	
	@ParameterizedTest
	@CsvFileSource(resources = "functionalClasses.txt", delimiter='\t')
	public void testFunctionalClassNomenclature(String name, String expectedInchi) {
		checkName(name, expectedInchi);
	}
	
	@ParameterizedTest
	@CsvFileSource(resources = "fusedRings.txt", delimiter='\t')
	public void testFusedRingNomenclature(String name, String expectedInchi) {
		checkName(name, expectedInchi);
	}
	
	@ParameterizedTest
	@CsvFileSource(resources = "inorganics.txt", delimiter='\t')
	public void testInorganicNomenclature(String name, String expectedInchi) {
		checkName(name, expectedInchi);
	}
	
	@ParameterizedTest
	@CsvFileSource(resources = "ions.txt", delimiter='\t')
	public void testIonNomenclature(String name, String expectedInchi) {
		checkName(name, expectedInchi);
	}
	
	@ParameterizedTest
	@CsvFileSource(resources = "spiro.txt", delimiter='\t')
	public void testSpiroNomenclature(String name, String expectedInchi) {
		checkName(name, expectedInchi);
	}
	
	@ParameterizedTest
	@CsvFileSource(resources = "organometallics.txt", delimiter='\t')
	public void testOrganoMetallics(String name, String expectedInchi) {
		checkName(name, expectedInchi);
	}

	@ParameterizedTest
	@CsvFileSource(resources = "implicitBracketting.txt", delimiter='\t')
	public void testImplicitBracketting(String name, String expectedInchi) {
		checkName(name, expectedInchi);
	}
	
	@ParameterizedTest
	@CsvFileSource(resources = "stereochemistry.txt", delimiter='\t')
	public void testStereochemistry(String name, String expectedInchi) {
		checkName(name, expectedInchi);
	}
	
	@ParameterizedTest
	@CsvFileSource(resources = "detachablePrefixes.txt", delimiter='\t')
	public void testDetachablePrefixes(String name, String expectedInchi) {
		checkName(name, expectedInchi);
	}

	@ParameterizedTest
	@CsvFileSource(resources = "lettercasing.txt", delimiter='\t')
	public void testLetterCasing(String name, String expectedInchi) {
		checkName(name, expectedInchi);
	}

	@ParameterizedTest
	@CsvFileSource(resources = "miscellany.txt", delimiter='\t')
	public void testMiscellany(String name, String expectedInchi) {
		checkName(name, expectedInchi);
	}
	
	private void checkName(String name, String expectedInchI) {
		String inchi = NameToInchi.convertResultToInChI(n2s.parseChemicalName(name, n2sConfig));
		if (inchi != null) {
			String opsinInchi = InchiPruner.mainChargeAndStereochemistryLayers(inchi);
			String referenceInchi = InchiPruner.mainChargeAndStereochemistryLayers(expectedInchI);
			assertEquals(referenceInchi, opsinInchi, name + " was misinterpreted as: " + inchi);
		} else {
			fail(name +" was uninterpretable");
		}
	}
}
