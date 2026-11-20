package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;

public class SuffixApplierTest {
	private static SuffixRules suffixRules;
	private BuildState state;
	private SuffixApplier suffixApplier;
	private Element group;

	@BeforeAll
	public static void loadRules() throws IOException {
		suffixRules = new SuffixRules(new ResourceGetter("uk/ac/cam/ch/wwmm/opsin/resources/"));
	}

	@BeforeEach
	public void setup() throws StructureBuildingException {
		state = new BuildState(new NameToStructureConfig());
		suffixApplier = new SuffixApplier(state, suffixRules);
		group = new TokenEl(GROUP_EL, "sil");
		group.setFrag(state.fragManager.buildSMILES("[SiH?]", STANDARDGROUP_TYPE_VAL, NUMERIC_LABELS_VAL));
	}

	@ParameterizedTest
	@CsvSource({"yl, 3", "ylidene, 2", "ylidyne, 1"})
	public void reservesOutwardBondValence(String suffixValue, int expectedCapacity) throws ComponentGenerationException {
		assertEquals(Integer.valueOf(expectedCapacity), capacity(suffix(suffixValue)));
		assertEquals(0, group.getFrag().getOutAtomCount());
	}

	@Test
	public void reservesBothHydroxyAttachmentsWithoutApplyingThem() throws StructureBuildingException, ComponentGenerationException {
		Element ol1 = suffixWithFragment("ol", "[*]O");
		Element ol2 = suffixWithFragment("ol", "[*]O");
		assertEquals(2, suffixApplier.calculateMononuclearParentSubstitutionCapacityAfterSuffixes(group, Arrays.asList(ol1, ol2)));
		assertEquals(0, group.getFrag().getFirstAtom().getBondCount());
	}

	@Test
	public void reservesDoubleBondToSuffix() throws StructureBuildingException, ComponentGenerationException {
		assertEquals(2, capacity(suffixWithFragment("one", "[*]=O")));
	}

	@Test
	public void declinesChargeChanges() throws ComponentGenerationException {
		assertNull(capacity(suffix("ium")));
	}

	private Integer capacity(Element suffix) throws ComponentGenerationException {
		return suffixApplier.calculateMononuclearParentSubstitutionCapacityAfterSuffixes(group, Collections.singletonList(suffix));
	}

	private Element suffix(String value) {
		Element suffix = new TokenEl(SUFFIX_EL, value);
		suffix.addAttribute(new Attribute(VALUE_ATR, value));
		return suffix;
	}

	private Element suffixWithFragment(String value, String smiles) throws StructureBuildingException {
		Element suffix = suffix(value);
		suffix.setFrag(state.fragManager.buildSMILES(smiles));
		return suffix;
	}

	@ParameterizedTest
	@CsvSource(delimiter = ';', value = {
		"4-methylacetophenone; 4'-methylacetophenone",
		"3,4-dimethoxyacetophenone; 3',4'-dimethoxyacetophenone",
		"4-methylpropiophenone; 4'-methylpropiophenone"
	})
	public void unprimedPhenoneRingLocants(String unprimed, String primed) {
		NameToStructure n2s = NameToStructure.getInstance();
		OpsinResult result = n2s.parseChemicalName(unprimed);
		assertEquals(OpsinResult.OPSIN_RESULT_STATUS.SUCCESS, result.getStatus(), result.getMessage());
		assertEquals(n2s.parseChemicalName(primed).getSmiles(), result.getSmiles());
	}

	@ParameterizedTest
	@CsvSource(delimiter = ';', value = {
		"2-bromoacetophenone; BrCC(=O)C1=CC=CC=C1",
		"3-chloropropiophenone; ClCCC(=O)C1=CC=CC=C1",
		"4-chlorobutyrophenone; ClCCCC(=O)C1=CC=CC=C1"
	})
	public void unprimedLocantPresentOnAcylChainStillMeansTheChain(String name, String smiles) {
		assertEquals(smiles, NameToStructure.getInstance().parseChemicalName(name).getSmiles());
	}

	@Test
	public void unprimedTwoMakesRingLocantsAmbiguous() {
		//2 may be the acyl carbon or the ring, so the ring numbering is not guessed
		assertEquals(OpsinResult.OPSIN_RESULT_STATUS.FAILURE, NameToStructure.getInstance().parseChemicalName("2,4-dihydroxyacetophenone").getStatus());
	}
}
