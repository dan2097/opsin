package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;

import java.util.List;

import static junit.framework.Assert.*;

import nu.xom.Element;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.mockito.Mockito.mock;
import uk.ac.cam.ch.wwmm.opsin.Atom;
import uk.ac.cam.ch.wwmm.opsin.Fragment;
import uk.ac.cam.ch.wwmm.opsin.NameToStructure;
import uk.ac.cam.ch.wwmm.opsin.BondStereo.BondStereoValue;
import uk.ac.cam.ch.wwmm.opsin.StereoAnalyser.StereoBond;
import uk.ac.cam.ch.wwmm.opsin.StereoAnalyser.StereoCentre;

public class StereochemistryTest {

	private static NameToStructure n2s;
	private static SMILESFragmentBuilder sBuilder;

	@BeforeClass
	public static void setup() throws NameToStructureException {
		n2s = NameToStructure.getInstance();
		sBuilder = new SMILESFragmentBuilder();
	}
	
	@AfterClass
	public static void cleanUp(){
		n2s = null;
		sBuilder =null;
	}
	
	/*
	 * Tests for finding stereo centres
	 */
	@Test
	public void findStereoCentresBromoChloroFluoroMethane() throws StructureBuildingException {
		Fragment f = n2s.parseChemicalName("bromochlorofluoromethane").getStructure();
		StereoAnalyser stereoAnalyser = new StereoAnalyser(f);
		assertEquals(1, stereoAnalyser.findStereoCentres().size());
		assertEquals(0, stereoAnalyser.findStereoBonds().size());
		StereoCentre sc = stereoAnalyser.findStereoCentres().get(0);
		assertNotNull(sc.getStereoAtom());
		Atom stereoAtom = sc.getStereoAtom();
		assertEquals("C", stereoAtom.getElement());
		assertEquals(4, stereoAtom.getID());
	}
	
	@Test
	public void findStereoCentresNacetylleucine() throws StructureBuildingException {
		Fragment f = n2s.parseChemicalName("N-acetylleucine").getStructure();
		StereoAnalyser stereoAnalyser = new StereoAnalyser(f);
		assertEquals(1, stereoAnalyser.findStereoCentres().size());
		assertEquals(0, stereoAnalyser.findStereoBonds().size());
		StereoCentre sc = stereoAnalyser.findStereoCentres().get(0);
		assertNotNull(sc.getStereoAtom());
		Atom stereoAtom = sc.getStereoAtom();
		assertEquals("C", stereoAtom.getElement());
		List<Atom> neighbours = sc.getCipOrderedAtoms();
		for (int i = 0; i < neighbours.size(); i++) {
			Atom a = neighbours.get(i);
			if (i==0){
				assertEquals(a.getElement(), "H");
			}
			else if (i==1){
				assertEquals(a.getElement(), "C");
			}
			else if (i==2){
				assertEquals(a.getElement(), "C");
			}
			else if (i==3){
				assertEquals(a.getElement(), "N");
			}
		}
	}
	
	@Test
	public void findStereoCentresBut2ene() throws StructureBuildingException {
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
		Element bondStereo = chiralBond.getBondStereo().toCML();
		assertNotNull(bondStereo);
		assertEquals(XmlDeclarations.CML_BONDSTEREO_EL, bondStereo.getLocalName());
		String atomRefs4 = bondStereo.getAttributeValue(XmlDeclarations.CML_ATOMREFS4_ATR);
		assertEquals("a1 a2 a3 a4", atomRefs4);
		assertEquals("a1 a2 a3 a4", atomRefs4);
		assertEquals(BondStereoValue.CIS.toString(), bondStereo.getValue());
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
		Element bondStereo = chiralBond.getBondStereo().toCML();
		assertNotNull(bondStereo);
		assertEquals(XmlDeclarations.CML_BONDSTEREO_EL, bondStereo.getLocalName());
		String atomRefs4 = bondStereo.getAttributeValue(XmlDeclarations.CML_ATOMREFS4_ATR);
		assertEquals("a1 a2 a3 a4", atomRefs4);
		assertEquals("a1 a2 a3 a4", atomRefs4);
		assertEquals(BondStereoValue.TRANS.toString(), bondStereo.getValue());
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
		Element bondStereo = chiralBond.getBondStereo().toCML();
		assertNotNull(bondStereo);
		assertEquals(XmlDeclarations.CML_BONDSTEREO_EL, bondStereo.getLocalName());
		String atomRefs4 = bondStereo.getAttributeValue(XmlDeclarations.CML_ATOMREFS4_ATR);
		assertEquals("a1 a2 a3 a4", atomRefs4);
		assertEquals("a1 a2 a3 a4", atomRefs4);
		assertEquals(BondStereoValue.CIS.toString(), bondStereo.getValue());
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
		Element bondStereo = chiralBond.getBondStereo().toCML();
		assertNotNull(bondStereo);
		assertEquals(XmlDeclarations.CML_BONDSTEREO_EL, bondStereo.getLocalName());
		String atomRefs4 = bondStereo.getAttributeValue(XmlDeclarations.CML_ATOMREFS4_ATR);
		assertEquals("a1 a2 a3 a4", atomRefs4);
		assertEquals("a1 a2 a3 a4", atomRefs4);
		assertEquals(BondStereoValue.TRANS.toString(), bondStereo.getValue());
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
		Element bondStereo = chiralBond.getBondStereo().toCML();
		assertNotNull(bondStereo);
		assertEquals(XmlDeclarations.CML_BONDSTEREO_EL, bondStereo.getLocalName());
		String atomRefs4 = bondStereo.getAttributeValue(XmlDeclarations.CML_ATOMREFS4_ATR);
		assertEquals("a1 a2 a3 a4", atomRefs4);
		assertEquals("a1 a2 a3 a4", atomRefs4);
		assertEquals(BondStereoValue.CIS.toString(), bondStereo.getValue());
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
		Element bondStereo = chiralBond.getBondStereo().toCML();
		assertNotNull(bondStereo);
		assertEquals(XmlDeclarations.CML_BONDSTEREO_EL, bondStereo.getLocalName());
		String atomRefs4 = bondStereo.getAttributeValue(XmlDeclarations.CML_ATOMREFS4_ATR);
		assertEquals("a1 a2 a3 a4", atomRefs4);
		assertEquals("a1 a2 a3 a4", atomRefs4);
		assertEquals(BondStereoValue.TRANS.toString(), bondStereo.getValue());
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
	public void testCIPpriority1() throws StructureBuildingException {
		FragmentManager fm = new FragmentManager(sBuilder, mock(CMLFragmentBuilder.class), new IDManager());
		Fragment f = sBuilder.build("C(Br)(F)([H])Cl", fm);
		StereoAnalyser stereoAnalyser = new StereoAnalyser(f);
		List<Atom> cipOrdered = stereoAnalyser.getNeighbouringAtomsInCIPOrder(f.getFirstAtom());
		for (int i = 0; i < cipOrdered.size(); i++) {
			Atom a = cipOrdered.get(i);
			if (i==0){
				assertEquals("H", a.getElement());
			}
			else if (i==1){
				assertEquals("F", a.getElement());
			}
			else if (i==2){
				assertEquals("Cl", a.getElement());
			}
			else if (i==3){
				assertEquals("Br", a.getElement());
			}
		}
	}
	
	@Test
	public void testCIPpriority2() throws StructureBuildingException {
		FragmentManager fm = new FragmentManager(sBuilder, mock(CMLFragmentBuilder.class), new IDManager());
		Fragment f = sBuilder.build("C([H])(C1CC1)(C1CCC1)O", fm);
		StereoAnalyser stereoAnalyser = new StereoAnalyser(f);
		List<Atom> cipOrdered = stereoAnalyser.getNeighbouringAtomsInCIPOrder(f.getFirstAtom());
		for (int i = 0; i < cipOrdered.size(); i++) {
			Atom a = cipOrdered.get(i);
			if (i==0){
				assertEquals("H", a.getElement());
			}
			else if (i==1){
				assertEquals(3, a.getID());
			}
			else if (i==2){
				assertEquals(6, a.getID());
			}
			else if (i==3){
				assertEquals("O", a.getElement());
			}
		}
	}
	
	
	@Test
	public void testCIPpriority3() throws StructureBuildingException {
		FragmentManager fm = new FragmentManager(sBuilder, mock(CMLFragmentBuilder.class), new IDManager());
		Fragment f = sBuilder.build("[C](N)(C1=CC(O)=CC=C1)([H])C2=CC=C(O)C=C2", fm);
		StereoAnalyser stereoAnalyser = new StereoAnalyser(f);
		List<Atom> cipOrdered = stereoAnalyser.getNeighbouringAtomsInCIPOrder(f.getFirstAtom());
		for (int i = 0; i < cipOrdered.size(); i++) {
			Atom a = cipOrdered.get(i);
			if (i==0){
				assertEquals("H", a.getElement());
			}
			else if (i==1){
				assertEquals(11, a.getID());
			}
			else if (i==2){
				assertEquals(3, a.getID());
			}
			else if (i==3){
				assertEquals("N", a.getElement());
			}
		}
	}
	
	@Test
	public void testCIPpriority4() throws StructureBuildingException {
		FragmentManager fm = new FragmentManager(sBuilder, mock(CMLFragmentBuilder.class), new IDManager());
		Fragment f = sBuilder.build("[C](N)(C1CC(O)CCC1)([H])C2CCC(O)CC2", fm);
		StereoAnalyser stereoAnalyser = new StereoAnalyser(f);
		List<Atom> cipOrdered = stereoAnalyser.getNeighbouringAtomsInCIPOrder(f.getFirstAtom());
		for (int i = 0; i < cipOrdered.size(); i++) {
			Atom a = cipOrdered.get(i);
			if (i==0){
				assertEquals("H", a.getElement());
			}
			else if (i==1){
				assertEquals(11, a.getID());
			}
			else if (i==2){
				assertEquals(3, a.getID());
			}
			else if (i==3){
				assertEquals("N", a.getElement());
			}
		}
	}

	@Test
	public void testCIPpriority5() throws StructureBuildingException {
		FragmentManager fm = new FragmentManager(sBuilder, mock(CMLFragmentBuilder.class), new IDManager());
		Fragment f = sBuilder.build("C1([H])(C(=O)O[H])C([H])([H])SC([H])([H])N([H])1", fm);
		StereoAnalyser stereoAnalyser = new StereoAnalyser(f);
		List<Atom> cipOrdered = stereoAnalyser.getNeighbouringAtomsInCIPOrder(f.getFirstAtom());
		for (int i = 0; i < cipOrdered.size(); i++) {
			Atom a = cipOrdered.get(i);
			if (i==0){
				assertEquals("H", a.getElement());
			}
			else if (i==1){
				assertEquals(3, a.getID());
			}
			else if (i==2){
				assertEquals(7, a.getID());
			}
			else if (i==3){
				assertEquals("N", a.getElement());
			}
		}
	}
	
	@Test
	public void testCIPpriority6() throws StructureBuildingException {
		FragmentManager fm = new FragmentManager(sBuilder, mock(CMLFragmentBuilder.class), new IDManager());
		Fragment f = sBuilder.build("C1([H])(O)C([H])(C([H])([H])[H])OC([H])([H])C([H])([H])C1([H])(O[H])", fm);
		StereoAnalyser stereoAnalyser = new StereoAnalyser(f);
		List<Atom> cipOrdered = stereoAnalyser.getNeighbouringAtomsInCIPOrder(f.getFirstAtom());
		for (int i = 0; i < cipOrdered.size(); i++) {
			Atom a = cipOrdered.get(i);
			if (i==0){
				assertEquals("H", a.getElement());
			}
			else if (i==1){
				assertEquals(17, a.getID());
			}
			else if (i==2){
				assertEquals(4, a.getID());
			}
			else if (i==3){
				assertEquals("O", a.getElement());
			}
		}
	}
	
	@Test
	public void testCIPpriority7() throws StructureBuildingException {
		FragmentManager fm = new FragmentManager(sBuilder, mock(CMLFragmentBuilder.class), new IDManager());
		Fragment f = sBuilder.build("[H]OC2([H])(C([H])([H])C([H])([H])C3([H])(C4([H])(C([H])([H])C([H])([H])C1=C([H])C([H])([H])C([H])([H])C([H])([H])C1([H])C4([H])(C([H])([H])C([H])([H])C23(C([H])([H])[H])))))", fm);
		StereoAnalyser stereoAnalyser = new StereoAnalyser(f);
		List<Atom> cipOrdered = stereoAnalyser.getNeighbouringAtomsInCIPOrder(f.getAtomList().get(34));
		for (int i = 0; i < cipOrdered.size(); i++) {
			Atom a = cipOrdered.get(i);
			if (i==0){
				assertEquals("H", a.getElement());
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
		StereoAnalyser stereoAnalyser = new StereoAnalyser(f);
		List<Atom> cipOrdered = stereoAnalyser.getNeighbouringAtomsInCIPOrder(f.getAtomByLocant("6a"));
		for (int i = 0; i < cipOrdered.size(); i++) {
			Atom a = cipOrdered.get(i);
			if (i==0){
				assertEquals("H", a.getElement());
			}
			else if (i==1){
				assertEquals("C", a.getElement());
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
		BuildState state  =new BuildState(mock(NameToStructureConfig.class), sBuilder, mock(CMLFragmentBuilder.class));
		Fragment f = state.fragManager.buildSMILES("C1(C=C)CC1C2=CC=CC=C2");
		StructureBuilder.makeHydrogensExplicit(state);
		StereoAnalyser stereoAnalyser = new StereoAnalyser(f);
		List<Atom> cipOrdered = stereoAnalyser.getNeighbouringAtomsInCIPOrder(f.getFirstAtom());
		for (int i = 0; i < cipOrdered.size(); i++) {
			Atom a = cipOrdered.get(i);
			if (i==0){
				assertEquals("H", a.getElement());
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
		BuildState state  =new BuildState(mock(NameToStructureConfig.class), sBuilder, mock(CMLFragmentBuilder.class));
		Fragment f = state.fragManager.buildSMILES("C(O[H])([H])(C1([H])C([H])(F)C([H])(Cl)C([H])([H])C([H])(I)C1([H])([H]))C1([H])C([H])(F)C([H])(Br)C([H])([H])C([H])(Cl)C1([H])([H])");
		StructureBuilder.makeHydrogensExplicit(state);
		StereoAnalyser stereoAnalyser = new StereoAnalyser(f);
		List<Atom> cipOrdered = stereoAnalyser.getNeighbouringAtomsInCIPOrder(f.getFirstAtom());
		for (int i = 0; i < cipOrdered.size(); i++) {
			Atom a = cipOrdered.get(i);
			if (i==0){
				assertEquals("H", a.getElement());
			}
			else if (i==1){
				assertEquals(5, a.getID());
			}
			else if (i==2){
				assertEquals(22, a.getID());
			}
			else if (i==3){
				assertEquals("O", a.getElement());
			}
		}
	}
	
	@Test
	public void testCIPpriority11() throws StructureBuildingException {
		BuildState state  =new BuildState(mock(NameToStructureConfig.class), sBuilder, mock(CMLFragmentBuilder.class));
		Fragment f = state.fragManager.buildSMILES("C17C=CC23C45OC6C19.O74.O2C3.C5.C6(C)C.C9");
		StructureBuilder.makeHydrogensExplicit(state);
		StereoAnalyser stereoAnalyser = new StereoAnalyser(f);
		//stereocentres at 1,4,5,7,8
		List<Atom> atomList = f.getAtomList();
		List<Atom> cipOrdered = stereoAnalyser.getNeighbouringAtomsInCIPOrder(atomList.get(0));
		for (int i = 0; i < cipOrdered.size(); i++) {
			Atom a = cipOrdered.get(i);
			if (i==0){
				assertEquals("H", a.getElement());
			}
			else if (i==1){
				assertEquals(2, a.getID());
			}
			else if (i==2){
				assertEquals(8, a.getID());
			}
			else if (i==3){
				assertEquals("O", a.getElement());
			}
		}
		cipOrdered = stereoAnalyser.getNeighbouringAtomsInCIPOrder(atomList.get(3));
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
				assertEquals("O", a.getElement());
			}
		}
		cipOrdered = stereoAnalyser.getNeighbouringAtomsInCIPOrder(atomList.get(4));
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
		cipOrdered = stereoAnalyser.getNeighbouringAtomsInCIPOrder(atomList.get(6));
		for (int i = 0; i < cipOrdered.size(); i++) {
			Atom a = cipOrdered.get(i);
			if (i==0){
				assertEquals("H", a.getElement());
			}
			else if (i==1){
				assertEquals(13, a.getID());
			}
			else if (i==2){
				assertEquals(8, a.getID());
			}
			else if (i==3){
				assertEquals("O", a.getElement());
			}
		}
		cipOrdered = stereoAnalyser.getNeighbouringAtomsInCIPOrder(atomList.get(7));
		for (int i = 0; i < cipOrdered.size(); i++) {
			Atom a = cipOrdered.get(i);
			if (i==0){
				assertEquals("H", a.getElement());
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
		BuildState state  =new BuildState(mock(NameToStructureConfig.class), sBuilder, mock(CMLFragmentBuilder.class));
		Fragment f = state.fragManager.buildSMILES("C1(C)(CCC(=O)N1)CCC(=O)NC(C)C");
		StructureBuilder.makeHydrogensExplicit(state);
		StereoAnalyser stereoAnalyser = new StereoAnalyser(f);
		List<Atom> cipOrdered = stereoAnalyser.getNeighbouringAtomsInCIPOrder(f.getFirstAtom());
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
				assertEquals("N", a.getElement());
			}
		}
	}
	
	@Test
	public void testAtomParityEquivalence1() throws StructureBuildingException {
		Atom a1= new Atom(1, "C", mock(Fragment.class));
		Atom a2= new Atom(2, "C", mock(Fragment.class));
		Atom a3= new Atom(3, "C", mock(Fragment.class));
		Atom a4= new Atom(4, "C", mock(Fragment.class));
		Atom[] atomRefs1 = new Atom[]{a1,a2,a3,a4};
		Atom[] atomRefs2 = new Atom[]{a3,a4,a1,a2};
		//2 swaps (4 by bubble sort)
		assertEquals(true, StereochemistryHandler.checkEquivalencyOfAtomsRefs4AndParity(atomRefs1, 1, atomRefs2, 1));
		assertEquals(false, StereochemistryHandler.checkEquivalencyOfAtomsRefs4AndParity(atomRefs1, 1, atomRefs2, -1));
	}
	
	@Test
	public void testAtomParityEquivalence2() throws StructureBuildingException {
		Atom a1= new Atom(1, "C", mock(Fragment.class));
		Atom a2= new Atom(2, "C", mock(Fragment.class));
		Atom a3= new Atom(3, "C", mock(Fragment.class));
		Atom a4= new Atom(4, "C", mock(Fragment.class));
		Atom[] atomRefs1 = new Atom[]{a1,a2,a3,a4};
		Atom[] atomRefs2 = new Atom[]{a2,a4,a1,a3};
		//3 swaps
		assertEquals(false, StereochemistryHandler.checkEquivalencyOfAtomsRefs4AndParity(atomRefs1, 1, atomRefs2, 1));
		assertEquals(true, StereochemistryHandler.checkEquivalencyOfAtomsRefs4AndParity(atomRefs1, 1, atomRefs2, -1));
	}
	
	@Test
	public void testCisTransUnambiguous() throws StructureBuildingException {
		FragmentManager fm = new FragmentManager(sBuilder, mock(CMLFragmentBuilder.class), new IDManager());
		Fragment f = sBuilder.build("[H]C([H])([H])C([H])=C([H])C([H])([H])[H]", fm);
		assertEquals(true, StereochemistryHandler.cisTranUnambiguousOnBond(f.findBond(5, 7)));
	}
	
	@Test
	public void testCisTransAmbiguous() throws StructureBuildingException {
		FragmentManager fm = new FragmentManager(sBuilder, mock(CMLFragmentBuilder.class), new IDManager());
		Fragment f = sBuilder.build("[H]C([H])([H])C(Cl)=C([H])C([H])([H])[H]", fm);
		assertEquals(false, StereochemistryHandler.cisTranUnambiguousOnBond(f.findBond(5, 7)));
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
}
