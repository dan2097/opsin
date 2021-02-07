package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.Assert.*;
import static org.mockito.Mockito.mock;

import java.util.ArrayList;
import java.util.List;

import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

import uk.ac.cam.ch.wwmm.opsin.BondStereo.BondStereoValue;
import uk.ac.cam.ch.wwmm.opsin.StereoAnalyser.StereoBond;
import uk.ac.cam.ch.wwmm.opsin.StereoAnalyser.StereoCentre;

public class StereochemistryTest {

	private FragmentManager fm;
	
	@Before
	public void setup() {
		IDManager idManager = new IDManager();
		fm = new FragmentManager(new SMILESFragmentBuilder(idManager), idManager);
	}
	
	private static NameToStructure n2s;

	@BeforeClass
	public static void intialSetup() {
		n2s = NameToStructure.getInstance();
	}
	
	@AfterClass
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

	@Test
	public void applyStereochemistryLocantedRSracemic() throws StructureBuildingException {
		Fragment f         = n2s.parseChemicalName("(1RS,2SR)-2-(methylamino)-1-phenylpropan-1-ol").getStructure();
		int      nRacAtoms = 0;
		for (Atom atom : f.getAtomList()) {
			if (atom.getAtomParity() != null && atom.getStereoGroup() == StereoGroup.Rac) {
				nRacAtoms++;
			}
		}
		assertEquals(2, nRacAtoms);
	}

	@Test
	public void applyStereochemistryLocantedRSrel() throws StructureBuildingException {
		Fragment f         = n2s.parseChemicalName("(1R*,2S*)-2-(methylamino)-1-phenylpropan-1-ol").getStructure();
		int      nRelAtoms = 0;
		for (Atom atom : f.getAtomList()) {
			if (atom.getAtomParity() != null && atom.getStereoGroup() == StereoGroup.Rel) {
				nRelAtoms++;
			}
		}
		assertEquals(2, nRelAtoms);
	}

	@Test
	public void applyStereochemistryLocantedPartialRac() throws StructureBuildingException {
		Fragment f         = n2s.parseChemicalName("(1RS,2R)-2-(methylamino)-1-phenylpropan-1-ol").getStructure();
		int      nRacAtoms = 0;
		int      nAbsAtoms = 0;
		for (Atom atom : f.getAtomList()) {
			if (atom.getAtomParity() != null) {
				if (atom.getStereoGroup() == StereoGroup.Rac)
					nRacAtoms++;
				else if (atom.getStereoGroup() == StereoGroup.Abs)
					nAbsAtoms++;
			}
		}
		assertEquals(1, nRacAtoms);
		assertEquals(1, nAbsAtoms);
	}

	@Test
	public void applyStereochemistryLocantedPartialRel() throws StructureBuildingException {
		Fragment f         = n2s.parseChemicalName("(1R*,2R)-2-(methylamino)-1-phenylpropan-1-ol").getStructure();
		int      nRelAtoms = 0;
		int      nAbsAtoms = 0;
		for (Atom atom : f.getAtomList()) {
			if (atom.getAtomParity() != null) {
				if (atom.getStereoGroup() == StereoGroup.Rel)
					nRelAtoms++;
				else if (atom.getStereoGroup() == StereoGroup.Abs)
					nAbsAtoms++;
			}
		}
		assertEquals(1, nRelAtoms);
		assertEquals(1, nAbsAtoms);
	}

	@Test
	public void applyStereochemistryRacemicUnlocanted() throws StructureBuildingException {
		Fragment f         = n2s.parseChemicalName("rac-1-phenylethan-1-ol").getStructure();
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
	public void applyStereochemistryRacemicMultipleUnlocanted() throws StructureBuildingException {
		Fragment f = n2s.parseChemicalName("rac-2-(methylamino)-1-phenylpropan-1-ol").getStructure();
		assertNotNull(f);
		int      nRacAtoms = 0;
		for (Atom atom : f.getAtomList()) {
			if (atom.getAtomParity() != null) {
				if (atom.getStereoGroup() == StereoGroup.Rac)
					nRacAtoms++;
			}
		}
		assertEquals(2, nRacAtoms);
	}

	@Test
	public void applyStereochemistryRelUnlocanted() throws StructureBuildingException {
		Fragment f         = n2s.parseChemicalName("rel-1-phenylethan-1-ol").getStructure();
		int      nRelAtoms = 0;
		for (Atom atom : f.getAtomList()) {
			if (atom.getAtomParity() != null) {
				if (atom.getStereoGroup() == StereoGroup.Rel)
					nRelAtoms++;
			}
		}
		assertEquals(1, nRelAtoms);
	}

	@Test
	public void applyStereochemistryRelUnlocanted2() throws StructureBuildingException {
		Fragment f         = n2s.parseChemicalName("(rac)-1-phenylethan-1-ol").getStructure();
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
	public void applyStereochemistryLocantedRorS() throws StructureBuildingException {
		// just for James Davidson, should be 1R*
		Fragment f         = n2s.parseChemicalName("(1R or S)-1-(1-pentyl-1H-pyrazol-5-yl)ethanol").getStructure();
		int      nRacAtoms = 0;
		for (Atom atom : f.getAtomList()) {
			if (atom.getAtomParity() != null) {
				if (atom.getStereoGroup() == StereoGroup.Rel)
					nRacAtoms++;
			}
		}
		assertEquals(1, nRacAtoms);
	}

	// US20080015199A1_2830
	// probably better to support via composite entity like techniques
	@Ignore
	public void applyStereochemistryRelUnlocantedRAndS() throws StructureBuildingException {
		Fragment f         = n2s.parseChemicalName("(R) and (S)-4-{3-[(4-Carbamimidoylphenylamino)-(3,5-dimethoxyphenyl)methyl]-5-oxo-4,5-dihydro-[1,2,4]triazol-1-yl}thiazole-5-carboxylic acid").getStructure();
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
	
	@Test(expected=CipOrderingException.class)
	public void testCipUnassignable() throws StructureBuildingException {
		//two sides of ring are identical
		Fragment f = fm.buildSMILES("NC1(O)CCC(CCC2CCCCC2)CC1");
		new CipSequenceRules(f.getAtomList().get(1)).getNeighbouringAtomsInCipOrder();
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
		assertNull("This atom is no longer a chiral centre and hence should not have an associated atom parity", formerChiralCentre.getAtomParity());
	}
	
	@Test
	public void testChiralBondWhichBecomesAchiral() throws StructureBuildingException {
		Fragment f = n2s.parseChemicalName("3-methylcrotonic acid").getStructure();
		StereoAnalyser stereoAnalyser = new StereoAnalyser(f);
		assertEquals(0, stereoAnalyser.findStereoCentres().size());
		assertEquals(0, stereoAnalyser.findStereoBonds().size());
		Bond formerChiralBond = f.getAtomByLocantOrThrow("2").getBondToAtomOrThrow(f.getAtomByLocantOrThrow("3"));
		assertNull("This Bond is no longer a chiral centre and hence should not have an associated bond stereo", formerChiralBond.getBondStereo());
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
