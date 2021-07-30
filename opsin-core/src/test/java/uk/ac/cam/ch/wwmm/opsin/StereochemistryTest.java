package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.mockito.Mockito.mock;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.hamcrest.MatcherAssert;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

import org.hamcrest.CoreMatchers;

import java.util.AbstractMap;
import java.util.Iterator;

import uk.ac.cam.ch.wwmm.opsin.BondStereo.BondStereoValue;
import uk.ac.cam.ch.wwmm.opsin.OpsinResult.OPSIN_RESULT_STATUS;
import uk.ac.cam.ch.wwmm.opsin.OpsinWarning.OpsinWarningType;
import uk.ac.cam.ch.wwmm.opsin.StereoAnalyser.StereoBond;
import uk.ac.cam.ch.wwmm.opsin.StereoAnalyser.StereoCentre;

public class StereochemistryTest {

	private FragmentManager fm;
	
	@BeforeEach
	public void setup() {
		IDManager idManager = new IDManager();
		fm = new FragmentManager(new SMILESFragmentBuilder(idManager), idManager);
	}
	
	private static NameToStructure n2s;

	@BeforeAll
	public static void intialSetup() {
		n2s = NameToStructure.getInstance();
	}
	
	@AfterAll
	public static void cleanUp(){
		n2s = null;
	}
	
	/*
	 * Tests for finding stereo centres
	 */
	@Test
	public void findStereoCentresBromoChloroFluoroMethane() {
		Fragment f = n2s.parseChemicalName("bromochlorofluoromethane").getStructure();
		StereoAnalyser stereoAnalyser = new StereoAnalyser(f);
		assertEquals(1, stereoAnalyser.findStereoCentres().size());
		assertEquals(0, stereoAnalyser.findStereoBonds().size());
		StereoCentre sc = stereoAnalyser.findStereoCentres().get(0);
		assertNotNull(sc.getStereoAtom());
		Atom stereoAtom = sc.getStereoAtom();
		assertEquals(ChemEl.C, stereoAtom.getElement());
		assertEquals(4, stereoAtom.getID());
	}
	
	@Test
	public void findStereoCentresNacetylleucine() throws CipOrderingException {
		Fragment f = n2s.parseChemicalName("N-acetylleucine").getStructure();
		StereoAnalyser stereoAnalyser = new StereoAnalyser(f);
		assertEquals(1, stereoAnalyser.findStereoCentres().size());
		assertEquals(0, stereoAnalyser.findStereoBonds().size());
		StereoCentre sc = stereoAnalyser.findStereoCentres().get(0);
		assertNotNull(sc.getStereoAtom());
		Atom stereoAtom = sc.getStereoAtom();
		assertEquals(ChemEl.C, stereoAtom.getElement());
		List<Atom> neighbours = sc.getCipOrderedAtoms();
		for (int i = 0; i < neighbours.size(); i++) {
			Atom a = neighbours.get(i);
			if (i==0){
				assertEquals(ChemEl.H, a.getElement());
			}
			else if (i==1){
				assertEquals(ChemEl.C, a.getElement());
			}
			else if (i==2){
				assertEquals(ChemEl.C, a.getElement());
			}
			else if (i==3){
				assertEquals(ChemEl.N, a.getElement());
			}
		}
	}
	
	@Test
	public void findStereoCentresBut2ene() {
		Fragment f = n2s.parseChemicalName("but-2-ene").getStructure();
		StereoAnalyser stereoAnalyser = new StereoAnalyser(f);
		assertEquals(0, stereoAnalyser.findStereoCentres().size());
		assertEquals(1, stereoAnalyser.findStereoBonds().size());
		StereoBond sb = stereoAnalyser.findStereoBonds().get(0);
		Bond stereoBond = sb.getBond();
		assertNotNull(stereoBond);
		Atom stereoAtom1 = stereoBond.getFromAtom();
		Atom stereoAtom2 = stereoBond.getToAtom();
		assertNotNull(stereoAtom1);
		assertNotNull(stereoAtom2);
		assertNotSame(stereoAtom1, stereoAtom2);
		if (stereoAtom1.getID() == 2){
			assertEquals(3, stereoAtom2.getID());
		}
		else{
			assertEquals(2, stereoAtom2.getID());
			assertEquals(3, stereoAtom1.getID());
		}
	}
	
	/*
	 * Tests for applying stereochemistry
	 */
	
	@Test
	public void applyStereochemistryLocantedZ() throws StructureBuildingException {
		Fragment f = n2s.parseChemicalName("(2Z)-but-2-ene").getStructure();
		Atom atom2 = f.getAtomByLocant("2");
		Atom atom3 = f.getAtomByLocant("3");
		assertNotNull(atom2);
		assertNotNull(atom3);
		Bond chiralBond = atom2.getBondToAtom(atom3);
		assertNotNull(chiralBond);
		BondStereo bondStereo = chiralBond.getBondStereo();
		assertNotNull(bondStereo);
		assertEquals("1 2 3 4", atomRefsToIdStr(bondStereo.getAtomRefs4()));
		assertEquals(BondStereoValue.CIS, bondStereo.getBondStereoValue());
	}
	
	@Test
	public void applyStereochemistryLocantedE() throws StructureBuildingException {
		Fragment f = n2s.parseChemicalName("(2E)-but-2-ene").getStructure();
		Atom atom2 = f.getAtomByLocant("2");
		Atom atom3 = f.getAtomByLocant("3");
		assertNotNull(atom2);
		assertNotNull(atom3);
		Bond chiralBond = atom2.getBondToAtom(atom3);
		assertNotNull(chiralBond);
		BondStereo bondStereo = chiralBond.getBondStereo();
		assertNotNull(bondStereo);
		assertEquals("1 2 3 4", atomRefsToIdStr(bondStereo.getAtomRefs4()));
		assertEquals(BondStereoValue.TRANS, bondStereo.getBondStereoValue());
	}

	@Test
	public void applyStereochemistryUnlocantedZ() throws StructureBuildingException {
		Fragment f = n2s.parseChemicalName("(Z)-but-2-ene").getStructure();
		Atom atom2 = f.getAtomByLocant("2");
		Atom atom3 = f.getAtomByLocant("3");
		assertNotNull(atom2);
		assertNotNull(atom3);
		Bond chiralBond = atom2.getBondToAtom(atom3);
		assertNotNull(chiralBond);
		BondStereo bondStereo = chiralBond.getBondStereo();
		assertNotNull(bondStereo);
		assertEquals("1 2 3 4", atomRefsToIdStr(bondStereo.getAtomRefs4()));
		assertEquals(BondStereoValue.CIS, bondStereo.getBondStereoValue());
	}
	
	@Test
	public void applyStereochemistryUnlocantedE() throws StructureBuildingException {
		Fragment f = n2s.parseChemicalName("(E)-but-2-ene").getStructure();
		Atom atom2 = f.getAtomByLocant("2");
		Atom atom3 = f.getAtomByLocant("3");
		assertNotNull(atom2);
		assertNotNull(atom3);
		Bond chiralBond = atom2.getBondToAtom(atom3);
		assertNotNull(chiralBond);
		BondStereo bondStereo = chiralBond.getBondStereo();
		assertNotNull(bondStereo);
		assertEquals("1 2 3 4", atomRefsToIdStr(bondStereo.getAtomRefs4()));
		assertEquals(BondStereoValue.TRANS, bondStereo.getBondStereoValue());
	}
	
	@Test
	public void applyStereochemistryCis() throws StructureBuildingException {
		Fragment f = n2s.parseChemicalName("cis-but-2-ene").getStructure();
		Atom atom2 = f.getAtomByLocant("2");
		Atom atom3 = f.getAtomByLocant("3");
		assertNotNull(atom2);
		assertNotNull(atom3);
		Bond chiralBond = atom2.getBondToAtom(atom3);
		assertNotNull(chiralBond);
		BondStereo bondStereo = chiralBond.getBondStereo();
		assertNotNull(bondStereo);
		assertEquals("1 2 3 4", atomRefsToIdStr(bondStereo.getAtomRefs4()));
		assertEquals(BondStereoValue.CIS, bondStereo.getBondStereoValue());
	}
	
	@Test
	public void applyStereochemistryTrans() throws StructureBuildingException {
		Fragment f = n2s.parseChemicalName("trans-but-2-ene").getStructure();
		Atom atom2 = f.getAtomByLocant("2");
		Atom atom3 = f.getAtomByLocant("3");
		assertNotNull(atom2);
		assertNotNull(atom3);
		Bond chiralBond = atom2.getBondToAtom(atom3);
		assertNotNull(chiralBond);
		BondStereo bondStereo = chiralBond.getBondStereo();
		assertNotNull(bondStereo);
		assertEquals("1 2 3 4", atomRefsToIdStr(bondStereo.getAtomRefs4()));
		assertEquals(BondStereoValue.TRANS, bondStereo.getBondStereoValue());
	}
	
	
	@Test
	public void applyStereochemistryLocantedRS() throws StructureBuildingException {
		Fragment f = n2s.parseChemicalName("(1S,2R)-2-(methylamino)-1-phenylpropan-1-ol").getStructure();
		List<Atom> atomList = f.getAtomList();
		List<Atom> stereoAtoms = new ArrayList<Atom>();
		for (Atom atom : atomList) {
			if (atom.getAtomParity() != null){
				stereoAtoms.add(atom);
			}
		}
		assertEquals(2, stereoAtoms.size());
		StereoAnalyser stereoAnalyser = new StereoAnalyser(f);
		List<StereoCentre> stereoCentres = stereoAnalyser.findStereoCentres();
		assertEquals(2, stereoCentres.size());
		if (stereoCentres.get(0).getStereoAtom().equals(stereoAtoms.get(0))){
			assertEquals(stereoCentres.get(1).getStereoAtom(), stereoAtoms.get(1));
		}
		else{
			assertEquals(stereoCentres.get(0).getStereoAtom(), stereoAtoms.get(1));
			assertEquals(stereoCentres.get(1).getStereoAtom(), stereoAtoms.get(0));
		}
	}

	/**
	 * Check the number of stereo atoms in a molecule name are assigned to the
	 * correct groups.
	 * @param name chemical name
	 * @param nRacExp number of stereo centers expected to be racemic
	 * @param nRelExp number of stereo centers expected to be relative
	 * @param nAbsExp number of stereo centers expected to be absolute
	 */
	void assertEnhancedStereo(String name, int nRacExp, int nRelExp, int nAbsExp) {
		Fragment f = n2s.parseChemicalName(name).getStructure();
		int nRacAtoms = 0;
		int nRelAtoms = 0;
		int nAbsAtoms = 0;
		for (Atom atom : f.getAtomList()) {
			if (atom.getAtomParity() != null) {
				if (atom.getStereoGroup() == StereoGroup.Rac)
					nRacAtoms++;
				if (atom.getStereoGroup() == StereoGroup.Rel)
					nRelAtoms++;
				else if (atom.getStereoGroup() == StereoGroup.Abs)
					nAbsAtoms++;
			}
		}
		assertEquals(nRacExp, nRacAtoms, "Incorrect number of racemic stereo centers");
		assertEquals(nRelExp, nRelAtoms, "Incorrect number of relative stereo centers");
		assertEquals(nAbsExp, nAbsAtoms, "Incorrect number of absolute stereo centers");
	}

	@Test
	public void applyStereochemistryLocantedRSracemic() throws StructureBuildingException {
		assertEnhancedStereo("(1RS,2SR)-2-(methylamino)-1-phenylpropan-1-ol", 2, 0, 0);
	}

	@Test
	public void applyStereochemistryLocantedRSrel() throws StructureBuildingException {
		assertEnhancedStereo("(1R*,2S*)-2-(methylamino)-1-phenylpropan-1-ol", 0, 2, 0);
	}

	@Test
	public void applyStereochemistryLocantedPartialRac() throws StructureBuildingException {
		assertEnhancedStereo("(1RS,2R)-2-(methylamino)-1-phenylpropan-1-ol", 1, 0, 1);
	}

	@Test
	public void applyStereochemistryLocantedPartialRel() throws StructureBuildingException {
		assertEnhancedStereo("(1R*,2R)-2-(methylamino)-1-phenylpropan-1-ol", 0, 1, 1);
	}

	@Test
	public void applyStereochemistryRacSlash() throws StructureBuildingException {
		assertEnhancedStereo("(1R/S,2R)-2-(methylamino)-1-phenylpropan-1-ol", 1, 0, 1);
	}

	@Test
	public void applyStereochemistryRelHatStar() throws StructureBuildingException {
		assertEnhancedStereo("(1R^*,2S^*)-2-(methylamino)-1-phenylpropan-1-ol", 0, 2, 0);
	}

	@Test
	public void applyStereochemistryRacemicUnlocanted() throws StructureBuildingException {
		assertEnhancedStereo("rac-1-phenylethan-1-ol", 1, 0, 0);
	}

	@Disabled("not allowed")
	public void applyStereochemistryRacemicMultipleUnlocanted() throws StructureBuildingException {
		assertEnhancedStereo("rac-2-(methylamino)-1-phenylpropan-1-ol", 2, 0, 0);
	}

	@Test
	public void applyStereochemistryRelUnlocanted() throws StructureBuildingException {
		assertEnhancedStereo("rel-1-phenylethan-1-ol", 0, 1, 0);
	}

	@Test
	public void applyStereochemistryRelUnlocanted2() throws StructureBuildingException {
		assertEnhancedStereo("(rac)-1-phenylethan-1-ol", 1, 0, 0);
	}

	@Test
	public void prefixTakesPrecedence() throws StructureBuildingException {
		assertEnhancedStereo("rac-(R*)-1-phenylethan-1-ol", 1, 0, 0);
		assertEnhancedStereo("rel-(RS)-1-phenylethan-1-ol", 0, 1, 0);
	}

	@Test
	public void applyStereochemistryLocantedRorS() throws StructureBuildingException {
		// just for James Davidson, should be rel-(1R)- or 1R*
		assertEnhancedStereo("(1R or S)-1-(1-pentyl-1H-pyrazol-5-yl)ethanol", 0, 1, 0);
	}

	@Test
	public void applyStereochemistryRacCis() throws StructureBuildingException {
		// racemic cis
		assertEnhancedStereo("rac-cis-N4-(2,2-dimethyl-3,4-dihydro-3-oxo-2H-pyrido[3,2-b][1,4]oxazin-6-yl)-N2-[6-[2,6-dimethylmorpholino)pyridin-3-yl]-5-fluoro-2,4-pyrimidinediamine", 2, 0, 0);
	}

	@Test
	public void applyStereochemistryPlusMinus() throws StructureBuildingException {
		assertEnhancedStereo("(+/-)-1-(1-pentyl-1H-pyrazol-5-yl)ethanol", 1, 0, 0);
		assertEnhancedStereo("(±)-1-(1-pentyl-1H-pyrazol-5-yl)ethanol", 1, 0, 0);
	}

	@Test
	public void testBracketNormalisation() throws StereochemistryException {
		MatcherAssert.assertThat(ComponentGenerator.normaliseBinaryBrackets("(R)-and(S)-"),
				CoreMatchers.is("(RS)"));
		MatcherAssert.assertThat(ComponentGenerator.normaliseBinaryBrackets("(R,S)-and(S,R)-"),
				CoreMatchers.is("(RS,SR)"));
		MatcherAssert.assertThat(ComponentGenerator.normaliseBinaryBrackets("(2R,3S)-and(2S,3S)-"),
				CoreMatchers.is("(2RS,3S)"));
		MatcherAssert.assertThat(ComponentGenerator.normaliseBinaryBrackets("(2R,3S)-or(2S,3S)-"),
				CoreMatchers.is("(2R*,3S)"));
	}

	@Test
	public void applyStereochemistryMultipleBrackets() throws StructureBuildingException {
		assertEnhancedStereo("(R)- and (S)-1-(1-pentyl-1H-pyrazol-5-yl)ethanol", 1, 0, 0);
		assertEnhancedStereo("(R)- or (S)-1-(1-pentyl-1H-pyrazol-5-yl)ethanol", 0, 1, 0);
		assertEnhancedStereo("(R,S)- or (S,S)-2-(methylamino)-1-phenylpropan-1-ol", 0, 1, 1);
	}

	@Test
	public void onlyApplyRacToPostfix() throws StructureBuildingException {
		assertEnhancedStereo("alanyl-rac-alanine", 1, 0, 1);
	}

	@Test
	public void remoteRacSpecification() throws StructureBuildingException {
		assertEnhancedStereo("rac-tert-butyl 7-[8-(tert-butoxycarbonylamino)-7-fluoro-3-[[(1S,2S,3R)-3-hydroxy-2,3-dimethyl-cyclobutoxy]carbonylamino]-6-isoquinolyl]-8-methyl-2,3-dihydropyrido[2,3-b][1,4]oxazine-1-carboxylate", 3, 0, 0);
		assertEnhancedStereo("(+-)-tert-butyl 7-[8-(tert-butoxycarbonylamino)-7-fluoro-3-[[(1S,2S,3R)-3-hydroxy-2,3-dimethyl-cyclobutoxy]carbonylamino]-6-isoquinolyl]-8-methyl-2,3-dihydropyrido[2,3-b][1,4]oxazine-1-carboxylate", 3, 0, 0);
	}

	// US20080015199A1_2830
	@Test
	public void applyStereochemistryRelUnlocantedRAndS() throws StructureBuildingException {
		assertEnhancedStereo("(R) and (S)-4-{3-[(4-Carbamimidoylphenylamino)-(3,5-dimethoxyphenyl)methyl]-5-oxo-4,5-dihydro-[1,2,4]triazol-1-yl}thiazole-5-carboxylic acid", 1, 0, 0);
	}

	@Test
	public void racemicPeptides() throws StructureBuildingException {
		Fragment f = n2s.parseChemicalName("DL-alanyl-DL-alanine").getStructure();
		Map<Map.Entry<StereoGroup,Integer>, Integer> counter = new HashMap<>();
		for (Atom atom : f.getAtomList()) {
			if (atom.getAtomParity() != null) {
				Map.Entry<StereoGroup,Integer> key
						= new AbstractMap.SimpleImmutableEntry<>(atom.getAtomParity().getStereoGroup(),
						                                         atom.getAtomParity().getStereoGroupNum());
				Integer count = counter.get(key);
				if (count == null)
					count = 0;
				counter.put(key, count+1);
			}
		}
		assertEquals(2, counter.size());
		Iterator<Integer> iterator = counter.values().iterator();
		assertEquals(1, (int)iterator.next());
		assertEquals(1, (int)iterator.next());
	}

	@Test
	public void racemicCarbohydrates() throws StructureBuildingException {
		Fragment f = n2s.parseChemicalName("4-O-α-DL-Glucopyranosyl-α-DL-glucose").getStructure();
		Map<Map.Entry<StereoGroup,Integer>, Integer> counter = new HashMap<>();
		for (Atom atom : f.getAtomList()) {
			if (atom.getAtomParity() != null &&
				atom.getStereoGroup() == StereoGroup.Rac) {
				Map.Entry<StereoGroup,Integer> key
						= new AbstractMap.SimpleImmutableEntry<>(atom.getAtomParity().getStereoGroup(),
						atom.getAtomParity().getStereoGroupNum());
				Integer count = counter.get(key);
				if (count == null)
					count = 0;
				counter.put(key, count+1);
			}
		}
		assertEquals(2, counter.size());
		Iterator<Integer> iterator = counter.values().iterator();
		assertEquals(4, (int)iterator.next());
		assertEquals(4, (int)iterator.next());
	}

	@Test
	public void avoidCollisionOfRacemicDefinitions() throws StructureBuildingException {
		Fragment f = n2s.parseChemicalName("DL-alanyl-(RS)-butan-2-ol").getStructure();
		Map<Map.Entry<StereoGroup,Integer>, Integer> counter = new HashMap<>();
		for (Atom atom : f.getAtomList()) {
			if (atom.getAtomParity() != null) {
				Map.Entry<StereoGroup,Integer> key
						= new AbstractMap.SimpleImmutableEntry<>(atom.getAtomParity().getStereoGroup(),
						atom.getAtomParity().getStereoGroupNum());
				Integer count = counter.get(key);
				if (count == null)
					count = 0;
				counter.put(key, count+1);
			}
		}
		assertEquals(2, counter.size());
		Iterator<Integer> iterator = counter.values().iterator();
		assertEquals(1, (int)iterator.next());
		assertEquals(1, (int)iterator.next());
	}

	@Test
	public void applyStereochemistryLocantedRandS() throws StructureBuildingException {
		Fragment f         = n2s.parseChemicalName("(1R and S)-1-(1-pentyl-1H-pyrazol-5-yl)ethanol").getStructure();
		int      nRacAtoms = 0;
		for (Atom atom : f.getAtomList()) {
			if (atom.getAtomParity() != null) {
				if (atom.getStereoGroup() == StereoGroup.Rac)
					nRacAtoms++;
			}
		}
		assertEquals(1, nRacAtoms);
	}
	
	@Test
	public void testCIPpriority1() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("C(Br)(F)([H])Cl");
		List<Atom> cipOrdered = new CipSequenceRules(f.getFirstAtom()).getNeighbouringAtomsInCipOrder();
		for (int i = 0; i < cipOrdered.size(); i++) {
			Atom a = cipOrdered.get(i);
			if (i==0){
				assertEquals(ChemEl.H, a.getElement());
			}
			else if (i==1){
				assertEquals(ChemEl.F, a.getElement());
			}
			else if (i==2){
				assertEquals(ChemEl.Cl, a.getElement());
			}
			else if (i==3){
				assertEquals(ChemEl.Br, a.getElement());
			}
		}
	}
	
	@Test
	public void testCIPpriority2() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("C([H])(C1CC1)(C1CCC1)O");
		List<Atom> cipOrdered =  new CipSequenceRules(f.getFirstAtom()).getNeighbouringAtomsInCipOrder();
		for (int i = 0; i < cipOrdered.size(); i++) {
			Atom a = cipOrdered.get(i);
			if (i==0){
				assertEquals(ChemEl.H, a.getElement());
			}
			else if (i==1){
				assertEquals(3, a.getID());
			}
			else if (i==2){
				assertEquals(6, a.getID());
			}
			else if (i==3){
				assertEquals(ChemEl.O, a.getElement());
			}
		}
	}
	
	
	@Test
	public void testCIPpriority3() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("[C](N)(C1=CC(O)=CC=C1)([H])C2=CC=C(O)C=C2");
		List<Atom> cipOrdered =  new CipSequenceRules(f.getFirstAtom()).getNeighbouringAtomsInCipOrder();
		for (int i = 0; i < cipOrdered.size(); i++) {
			Atom a = cipOrdered.get(i);
			if (i==0){
				assertEquals(ChemEl.H, a.getElement());
			}
			else if (i==1){
				assertEquals(11, a.getID());
			}
			else if (i==2){
				assertEquals(3, a.getID());
			}
			else if (i==3){
				assertEquals(ChemEl.N, a.getElement());
			}
		}
	}
	
	@Test
	public void testCIPpriority4() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("[C](N)(C1CC(O)CCC1)([H])C2CCC(O)CC2");
		List<Atom> cipOrdered = new CipSequenceRules(f.getFirstAtom()).getNeighbouringAtomsInCipOrder();
		for (int i = 0; i < cipOrdered.size(); i++) {
			Atom a = cipOrdered.get(i);
			if (i==0){
				assertEquals(ChemEl.H, a.getElement());
			}
			else if (i==1){
				assertEquals(11, a.getID());
			}
			else if (i==2){
				assertEquals(3, a.getID());
			}
			else if (i==3){
				assertEquals(ChemEl.N, a.getElement());
			}
		}
	}

	@Test
	public void testCIPpriority5() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("C1([H])(C(=O)O[H])C([H])([H])SC([H])([H])N([H])1");
		List<Atom> cipOrdered = new CipSequenceRules(f.getFirstAtom()).getNeighbouringAtomsInCipOrder();
		for (int i = 0; i < cipOrdered.size(); i++) {
			Atom a = cipOrdered.get(i);
			if (i==0){
				assertEquals(ChemEl.H, a.getElement());
			}
			else if (i==1){
				assertEquals(3, a.getID());
			}
			else if (i==2){
				assertEquals(7, a.getID());
			}
			else if (i==3){
				assertEquals(ChemEl.N, a.getElement());
			}
		}
	}
	
	@Test
	public void testCIPpriority6() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("C1([H])(O)C([H])(C([H])([H])[H])OC([H])([H])C([H])([H])C1([H])(O[H])");
		List<Atom> cipOrdered =  new CipSequenceRules(f.getFirstAtom()).getNeighbouringAtomsInCipOrder();
		for (int i = 0; i < cipOrdered.size(); i++) {
			Atom a = cipOrdered.get(i);
			if (i==0){
				assertEquals(ChemEl.H, a.getElement());
			}
			else if (i==1){
				assertEquals(17, a.getID());
			}
			else if (i==2){
				assertEquals(4, a.getID());
			}
			else if (i==3){
				assertEquals(ChemEl.O, a.getElement());
			}
		}
	}
	
	@Test
	public void testCIPpriority7() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("[H]OC2([H])(C([H])([H])C([H])([H])C3([H])(C4([H])(C([H])([H])C([H])([H])C1=C([H])C([H])([H])C([H])([H])C([H])([H])C1([H])C4([H])(C([H])([H])C([H])([H])C23(C([H])([H])[H])))))");
		List<Atom> cipOrdered =  new CipSequenceRules(f.getAtomList().get(34)).getNeighbouringAtomsInCipOrder();
		for (int i = 0; i < cipOrdered.size(); i++) {
			Atom a = cipOrdered.get(i);
			if (i==0){
				assertEquals(ChemEl.H, a.getElement());
			}
			else if (i==1){
				assertEquals(37, a.getID());
			}
			else if (i==2){
				assertEquals(13, a.getID());
			}
			else if (i==3){
				assertEquals(33, a.getID());
			}
		}
	}
	
	@Test
	public void testCIPpriority8() throws StructureBuildingException {
		Fragment f = n2s.parseChemicalName("(6aR)-6-phenyl-6,6a-dihydroisoindolo[2,1-a]quinazoline-5,11-dione").getStructure();
		List<Atom> cipOrdered = new CipSequenceRules(f.getAtomByLocant("6a")).getNeighbouringAtomsInCipOrder();
		for (int i = 0; i < cipOrdered.size(); i++) {
			Atom a = cipOrdered.get(i);
			if (i==0){
				assertEquals(ChemEl.H, a.getElement());
			}
			else if (i==1){
				assertEquals(ChemEl.C, a.getElement());
			}
			else if (i==2){
				assertEquals("6", a.getFirstLocant());
			}
			else if (i==3){
				assertEquals("12", a.getFirstLocant());
			}
		}
	}
	
	@Test
	public void testCIPpriority9() throws StructureBuildingException {
	  Fragment f = fm.buildSMILES("C1(C=C)CC1C2=CC=CC=C2");
		fm.makeHydrogensExplicit();
		List<Atom> cipOrdered = new CipSequenceRules(f.getFirstAtom()).getNeighbouringAtomsInCipOrder();
		for (int i = 0; i < cipOrdered.size(); i++) {
			Atom a = cipOrdered.get(i);
			if (i==0){
				assertEquals(ChemEl.H, a.getElement());
			}
			else if (i==1){
				assertEquals(4, a.getID());
			}
			else if (i==2){
				assertEquals(2, a.getID());
			}
			else if (i==3){
				assertEquals(5, a.getID());
			}
		}
	}
	
	@Test
	public void testCIPpriority10() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("C(O[H])([H])(C1([H])C([H])(F)C([H])(Cl)C([H])([H])C([H])(I)C1([H])([H]))C1([H])C([H])(F)C([H])(Br)C([H])([H])C([H])(Cl)C1([H])([H])");
		fm.makeHydrogensExplicit();
		List<Atom> cipOrdered = new CipSequenceRules(f.getFirstAtom()).getNeighbouringAtomsInCipOrder();
		for (int i = 0; i < cipOrdered.size(); i++) {
			Atom a = cipOrdered.get(i);
			if (i==0){
				assertEquals(ChemEl.H, a.getElement());
			}
			else if (i==1){
				assertEquals(5, a.getID());
			}
			else if (i==2){
				assertEquals(22, a.getID());
			}
			else if (i==3){
				assertEquals(ChemEl.O, a.getElement());
			}
		}
	}
	
	@Test
	public void testCIPpriority11() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("C17C=CC23C45OC6C19.O74.O2C3.C5.C6(C)C.C9");
		fm.makeHydrogensExplicit();
		//stereocentres at 1,4,5,7,8
		List<Atom> atomList = f.getAtomList();
		List<Atom> cipOrdered = new CipSequenceRules(atomList.get(0)).getNeighbouringAtomsInCipOrder();
		for (int i = 0; i < cipOrdered.size(); i++) {
			Atom a = cipOrdered.get(i);
			if (i==0){
				assertEquals(ChemEl.H, a.getElement());
			}
			else if (i==1){
				assertEquals(2, a.getID());
			}
			else if (i==2){
				assertEquals(8, a.getID());
			}
			else if (i==3){
				assertEquals(ChemEl.O, a.getElement());
			}
		}
		cipOrdered = new CipSequenceRules(atomList.get(3)).getNeighbouringAtomsInCipOrder();
		for (int i = 0; i < cipOrdered.size(); i++) {
			Atom a = cipOrdered.get(i);
			if (i==0){
				assertEquals(3, a.getID());
			}
			else if (i==1){
				assertEquals(11, a.getID());
			}
			else if (i==2){
				assertEquals(5, a.getID());
			}
			else if (i==3){
				assertEquals(ChemEl.O, a.getElement());
			}
		}
		cipOrdered = new CipSequenceRules(atomList.get(4)).getNeighbouringAtomsInCipOrder();
		for (int i = 0; i < cipOrdered.size(); i++) {
			Atom a = cipOrdered.get(i);
			if (i==0){
				assertEquals(12, a.getID());
			}
			else if (i==1){
				assertEquals(4, a.getID());
			}
			else if (i==2){
				assertEquals(6, a.getID());
			}
			else if (i==3){
				assertEquals(9, a.getID());
			}
		}
		cipOrdered = new CipSequenceRules(atomList.get(6)).getNeighbouringAtomsInCipOrder();
		for (int i = 0; i < cipOrdered.size(); i++) {
			Atom a = cipOrdered.get(i);
			if (i==0){
				assertEquals(ChemEl.H, a.getElement());
			}
			else if (i==1){
				assertEquals(13, a.getID());
			}
			else if (i==2){
				assertEquals(8, a.getID());
			}
			else if (i==3){
				assertEquals(ChemEl.O, a.getElement());
			}
		}
		cipOrdered = new CipSequenceRules(atomList.get(7)).getNeighbouringAtomsInCipOrder();
		for (int i = 0; i < cipOrdered.size(); i++) {
			Atom a = cipOrdered.get(i);
			if (i==0){
				assertEquals(ChemEl.H, a.getElement());
			}
			else if (i==1){
				assertEquals(16, a.getID());
			}
			else if (i==2){
				assertEquals(7, a.getID());
			}
			else if (i==3){
				assertEquals(1,  a.getID());
			}
		}
	}
	
	@Test
	public void testCIPpriority12() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("C1(C)(CCC(=O)N1)CCC(=O)NC(C)C");
		fm.makeHydrogensExplicit();
		List<Atom> cipOrdered = new CipSequenceRules(f.getFirstAtom()).getNeighbouringAtomsInCipOrder();
		for (int i = 0; i < cipOrdered.size(); i++) {
			Atom a = cipOrdered.get(i);
			if (i==0){
				assertEquals(2, a.getID());
			}
			else if (i==1){
				assertEquals(3, a.getID());
			}
			else if (i==2){
				assertEquals(8, a.getID());
			}
			else if (i==3){
				assertEquals(ChemEl.N, a.getElement());
			}
		}
	}
	
	@Test
	public void testCIPpriority13() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("C(O)(C#CC)C1=CC=CC=C1");
		fm.makeHydrogensExplicit();
		List<Atom> cipOrdered = new CipSequenceRules(f.getFirstAtom()).getNeighbouringAtomsInCipOrder();
		for (int i = 0; i < cipOrdered.size(); i++) {
			Atom a = cipOrdered.get(i);
			if (i==0){
				assertEquals(ChemEl.H, a.getElement());
			}
			else if (i==1){
				assertEquals(6, a.getID());
			}
			else if (i==2){
				assertEquals(3, a.getID());
			}
			else if (i==3){
				assertEquals(2, a.getID());
			}
		}
	}
	
	@Test
	public void testCIPpriority14() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("C(Cl)([2H])([3H])[H]");
		List<Atom> cipOrdered = new CipSequenceRules(f.getFirstAtom()).getNeighbouringAtomsInCipOrder();
		for (int i = 0; i < cipOrdered.size(); i++) {
			Atom a = cipOrdered.get(i);
			if (i==0){
				assertEquals(5, a.getID());
			}
			else if (i==1){
				assertEquals(3, a.getID());
			}
			else if (i==2){
				assertEquals(4, a.getID());
			}
			else if (i==3){
				assertEquals(2, a.getID());
			}
		}
	}
	
	@Test
	public void testCIPpriority15() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("C([H])(O)(C(C(F)CCl)CCBr)C(C(F)CF)CCI");
		fm.makeHydrogensExplicit();
		List<Atom> cipOrdered = new CipSequenceRules(f.getFirstAtom()).getNeighbouringAtomsInCipOrder();
		assertEquals(4, cipOrdered.size());
		assertEquals(2, cipOrdered.get(0).getID());
		assertEquals(12, cipOrdered.get(1).getID());
		assertEquals(4, cipOrdered.get(2).getID());
		assertEquals(3, cipOrdered.get(3).getID());
	}
	
	@Test()
	public void testCipUnassignable() {
		assertThrows(CipOrderingException.class, () -> {
			// two sides of ring are identical
			Fragment f = fm.buildSMILES("NC1(O)CCC(CCC2CCCCC2)CC1");
			new CipSequenceRules(f.getAtomList().get(1)).getNeighbouringAtomsInCipOrder();
		});
	}
	
	@Test
	public void testAtomParityEquivalence1() {
		Atom a1= new Atom(1, ChemEl.C, mock(Fragment.class));
		Atom a2= new Atom(2, ChemEl.C, mock(Fragment.class));
		Atom a3= new Atom(3, ChemEl.C, mock(Fragment.class));
		Atom a4= new Atom(4, ChemEl.C, mock(Fragment.class));
		Atom[] atomRefs1 = new Atom[]{a1,a2,a3,a4};
		Atom[] atomRefs2 = new Atom[]{a3,a4,a1,a2};
		//2 swaps (4 by bubble sort)
		assertEquals(true, StereochemistryHandler.checkEquivalencyOfAtomsRefs4AndParity(atomRefs1, 1, atomRefs2, 1));
		assertEquals(false, StereochemistryHandler.checkEquivalencyOfAtomsRefs4AndParity(atomRefs1, 1, atomRefs2, -1));
	}
	
	@Test
	public void testAtomParityEquivalence2() {
		Atom a1= new Atom(1, ChemEl.C, mock(Fragment.class));
		Atom a2= new Atom(2, ChemEl.C, mock(Fragment.class));
		Atom a3= new Atom(3, ChemEl.C, mock(Fragment.class));
		Atom a4= new Atom(4, ChemEl.C, mock(Fragment.class));
		Atom[] atomRefs1 = new Atom[]{a1,a2,a3,a4};
		Atom[] atomRefs2 = new Atom[]{a2,a4,a1,a3};
		//3 swaps
		assertEquals(false, StereochemistryHandler.checkEquivalencyOfAtomsRefs4AndParity(atomRefs1, 1, atomRefs2, 1));
		assertEquals(true, StereochemistryHandler.checkEquivalencyOfAtomsRefs4AndParity(atomRefs1, 1, atomRefs2, -1));
	}
	
	@Test
	public void testCisTransUnambiguous() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("[H]C([H])([H])C([H])=C([H])C([H])([H])[H]");
		assertEquals(true, StereochemistryHandler.cisTransUnambiguousOnBond(f.findBond(5, 7)));
	}
	
	@Test
	public void testCisTransAmbiguous() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("[H]C([H])([H])C(Cl)=C([H])C([H])([H])[H]");
		assertEquals(false, StereochemistryHandler.cisTransUnambiguousOnBond(f.findBond(5, 7)));
	}
	
	@Test
	public void testChiralAtomWhichBecomesAchiral() throws StructureBuildingException {
		Fragment f = n2s.parseChemicalName("alpha-amino-alanine").getStructure();
		StereoAnalyser stereoAnalyser = new StereoAnalyser(f);
		assertEquals(0, stereoAnalyser.findStereoCentres().size());
		assertEquals(0, stereoAnalyser.findStereoBonds().size());
		Atom formerChiralCentre = f.getAtomByLocantOrThrow("alpha");
		assertNull(formerChiralCentre.getAtomParity(), "This atom is no longer a chiral centre and hence should not have an associated atom parity");
	}
	
	@Test
	public void testChiralBondWhichBecomesAchiral() throws StructureBuildingException {
		Fragment f = n2s.parseChemicalName("3-methylcrotonic acid").getStructure();
		StereoAnalyser stereoAnalyser = new StereoAnalyser(f);
		assertEquals(0, stereoAnalyser.findStereoCentres().size());
		assertEquals(0, stereoAnalyser.findStereoBonds().size());
		Bond formerChiralBond = f.getAtomByLocantOrThrow("2").getBondToAtomOrThrow(f.getAtomByLocantOrThrow("3"));
		assertNull(formerChiralBond.getBondStereo(), "This Bond is no longer a chiral centre and hence should not have an associated bond stereo");
	}
	
	@Test
	public void testIsTetrahedral() throws StructureBuildingException {
		assertEquals(true, StereoAnalyser.isKnownPotentiallyStereogenic(fm.buildSMILES("C(N)(O)(Cl)Br").getFirstAtom()));
		assertEquals(true, StereoAnalyser.isKnownPotentiallyStereogenic(fm.buildSMILES("[Si](N)(O)(Cl)Br").getFirstAtom()));
		assertEquals(true, StereoAnalyser.isKnownPotentiallyStereogenic(fm.buildSMILES("[Ge](N)(O)(Cl)Br").getFirstAtom()));
		assertEquals(true, StereoAnalyser.isKnownPotentiallyStereogenic(fm.buildSMILES("[N+](N)(O)(Cl)Br").getFirstAtom()));
		assertEquals(true, StereoAnalyser.isKnownPotentiallyStereogenic(fm.buildSMILES("[P+](N)(O)(Cl)Br").getFirstAtom()));
		assertEquals(true, StereoAnalyser.isKnownPotentiallyStereogenic(fm.buildSMILES("[As+](N)(O)(Cl)Br").getFirstAtom()));
		assertEquals(true, StereoAnalyser.isKnownPotentiallyStereogenic(fm.buildSMILES("[B-](N)(O)(Cl)Br").getFirstAtom()));
		assertEquals(true, StereoAnalyser.isKnownPotentiallyStereogenic(fm.buildSMILES("[Sn](N)(O)(Cl)Br").getFirstAtom()));
		assertEquals(true, StereoAnalyser.isKnownPotentiallyStereogenic(fm.buildSMILES("[N](=N)(O)(Cl)Br").getFirstAtom()));
		assertEquals(true, StereoAnalyser.isKnownPotentiallyStereogenic(fm.buildSMILES("[P](=N)(O)(Cl)Br").getFirstAtom()));
		assertEquals(true, StereoAnalyser.isKnownPotentiallyStereogenic(fm.buildSMILES("[S](=N)(=O)(Cl)Br").getFirstAtom()));
		assertEquals(true, StereoAnalyser.isKnownPotentiallyStereogenic(fm.buildSMILES("[S+](=N)(O)(Cl)Br").getFirstAtom()));
		assertEquals(true, StereoAnalyser.isKnownPotentiallyStereogenic(fm.buildSMILES("[S](=O)(Cl)Br").getFirstAtom()));
		assertEquals(true, StereoAnalyser.isKnownPotentiallyStereogenic(fm.buildSMILES("[S+](O)(Cl)Br").getFirstAtom()));
		assertEquals(true, StereoAnalyser.isKnownPotentiallyStereogenic(fm.buildSMILES("N1(C)(OS1)").getFirstAtom()));
		assertEquals(true, StereoAnalyser.isKnownPotentiallyStereogenic(fm.buildSMILES("[Se](=N)(=O)(Cl)Br").getFirstAtom()));
		assertEquals(true, StereoAnalyser.isKnownPotentiallyStereogenic(fm.buildSMILES("[Se+](=N)(O)(Cl)Br").getFirstAtom()));
		assertEquals(true, StereoAnalyser.isKnownPotentiallyStereogenic(fm.buildSMILES("[Se](=O)(Cl)Br").getFirstAtom()));
		assertEquals(true, StereoAnalyser.isKnownPotentiallyStereogenic(fm.buildSMILES("[Se+](O)(Cl)Br").getFirstAtom()));
	}
	
	@Test
	public void testAchiralDueToResonance() throws StructureBuildingException {
		assertEquals(true, StereoAnalyser.isAchiralDueToResonanceOrTautomerism(fm.buildSMILES("[S](=N)(=O)([O-])Br").getFirstAtom()));
		assertEquals(true, StereoAnalyser.isAchiralDueToResonanceOrTautomerism(fm.buildSMILES("[S](=O)([O-])Br").getFirstAtom()));
		assertEquals(false, StereoAnalyser.isAchiralDueToResonanceOrTautomerism(fm.buildSMILES("[S](=S)([O-])Br").getFirstAtom()));
		assertEquals(false, StereoAnalyser.isAchiralDueToResonanceOrTautomerism(fm.buildSMILES("C(N)([O-])(Cl)Br").getFirstAtom()));
	}
	
	@Test
	public void testAchiralDueToTautomerism() throws StructureBuildingException {
		assertEquals(true, StereoAnalyser.isAchiralDueToResonanceOrTautomerism(fm.buildSMILES("[S](=N)(=O)([OH])Br").getFirstAtom()));
		assertEquals(true, StereoAnalyser.isAchiralDueToResonanceOrTautomerism(fm.buildSMILES("[S](=O)([OH])Br").getFirstAtom()));
		assertEquals(false, StereoAnalyser.isAchiralDueToResonanceOrTautomerism(fm.buildSMILES("[S](=S)([OH])Br").getFirstAtom()));
		assertEquals(false, StereoAnalyser.isAchiralDueToResonanceOrTautomerism(fm.buildSMILES("C(N)([OH])(Cl)Br").getFirstAtom()));
		assertEquals(true, StereoAnalyser.isAchiralDueToResonanceOrTautomerism(fm.buildSMILES("N([H])(CC)(C)").getFirstAtom()));
		assertEquals(false, StereoAnalyser.isAchiralDueToResonanceOrTautomerism(fm.buildSMILES("N1(C)(OS1)").getFirstAtom()));
	}
	
	@Test
	public void testFindPseudoAsymmetricCarbon1() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("OCC(O)C(O)C(O)CO");
		fm.makeHydrogensExplicit();
		StereoAnalyser stereoAnalyser = new StereoAnalyser(f);
		List<StereoCentre> stereoCentres = stereoAnalyser.findStereoCentres();
		assertEquals(3, stereoCentres.size());
		for (int i = 0; i < stereoCentres.size(); i++) {
			StereoCentre stereocentre = stereoCentres.get(i);
			if (i < 2){
				assertEquals(true, stereocentre.isTrueStereoCentre());
			}
			else{
				assertEquals(false, stereocentre.isTrueStereoCentre());
				assertEquals(5, stereocentre.getStereoAtom().getID());
			}
		}
	}
	
	@Test
	public void testFindPseudoAsymmetricCarbon2() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("OCC(O)C(C(Cl)(Br)C)(C(Cl)(Br)C)C(O)CO");
		fm.makeHydrogensExplicit();
		StereoAnalyser stereoAnalyser = new StereoAnalyser(f);
		List<StereoCentre> stereoCentres = stereoAnalyser.findStereoCentres();
		assertEquals(5, stereoCentres.size());
		for (int i = 0; i < stereoCentres.size(); i++) {
			StereoCentre stereocentre = stereoCentres.get(i);
			if (i <4){
				assertEquals(true, stereocentre.isTrueStereoCentre());
			}
			else{
				assertEquals(false, stereocentre.isTrueStereoCentre());
				assertEquals(5, stereocentre.getStereoAtom().getID());
			}
		}
	}
	
	@Test
	public void testAmbiguousStereoTerm() {
		OpsinResult result = n2s.parseChemicalName("trans-N-[2-Chloro-5-(2-methoxyethyl)benzyl]-N-cyclopropyl-4-hydroxy-4-(1-methyl-2-oxo-1,2-dihydro-4-pyridinyl)-3-piperidinecarboxamide");
		assertEquals(OPSIN_RESULT_STATUS.WARNING, result.getStatus());
		assertEquals(1, result.getWarnings().size());
		assertEquals(OpsinWarningType.APPEARS_AMBIGUOUS, result.getWarnings().get(0).getType());
	}

	private String atomRefsToIdStr(Atom[] atomRefs4) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < atomRefs4.length; i++) {
			sb.append(atomRefs4[i].getID());
			if (i + 1 < atomRefs4.length) {
				sb.append(' ');
			}
		}
		return sb.toString();
	}
}
