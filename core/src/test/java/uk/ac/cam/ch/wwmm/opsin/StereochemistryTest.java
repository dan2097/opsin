package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.List;

import static junit.framework.Assert.*;

import nu.xom.Element;

import org.junit.BeforeClass;
import org.junit.Test;

import uk.ac.cam.ch.wwmm.opsin.Atom;
import uk.ac.cam.ch.wwmm.opsin.Fragment;
import uk.ac.cam.ch.wwmm.opsin.NameToStructure;
import uk.ac.cam.ch.wwmm.opsin.StereoAnalyser.StereoBond;
import uk.ac.cam.ch.wwmm.opsin.StereoAnalyser.StereoCentre;
import uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.BondStereo;
public class StereochemistryTest {

	private static NameToStructure n2s;
	@BeforeClass
	public static void setup() throws Exception {
		n2s = NameToStructure.getInstance();
	}
	
	/*
	 * Tests for finding stereo centres
	 */
	@Test
	public void findStereoCentresBromoChloroFluoroMethane() throws StructureBuildingException {
		Fragment f = n2s.parseChemicalName("bromochlorofluoromethane", false).getStructure();
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
		Fragment f = n2s.parseChemicalName("N-acetylleucine", false).getStructure();
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
		Fragment f = n2s.parseChemicalName("but-2-ene", false).getStructure();
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
		Fragment f = n2s.parseChemicalName("(2Z)-but-2-ene", false).getStructure();
		Atom atom2 = f.getAtomByLocant("2");
		Atom atom3 = f.getAtomByLocant("3");
		assertNotNull(atom2);
		assertNotNull(atom3);
		Bond chiralBond = f.findBond(atom2, atom3);
		assertNotNull(chiralBond);
		Element bondStereo = chiralBond.getBondStereoElement();
		assertNotNull(bondStereo);
		assertEquals(XmlDeclarations.BONDSTEREO_EL, bondStereo.getLocalName());
		String atomRefs4 = bondStereo.getAttributeValue(XmlDeclarations.ATOMREFS4_ATR);
		assertEquals("a1 a2 a3 a4", atomRefs4);
		assertEquals("a1 a2 a3 a4", atomRefs4);
		assertEquals(BondStereo.CIS.toString(), bondStereo.getValue());
	}
	
	@Test
	public void applyStereochemistryLocantedE() throws StructureBuildingException {
		Fragment f = n2s.parseChemicalName("(2E)-but-2-ene", false).getStructure();
		Atom atom2 = f.getAtomByLocant("2");
		Atom atom3 = f.getAtomByLocant("3");
		assertNotNull(atom2);
		assertNotNull(atom3);
		Bond chiralBond = f.findBond(atom2, atom3);
		assertNotNull(chiralBond);
		Element bondStereo = chiralBond.getBondStereoElement();
		assertNotNull(bondStereo);
		assertEquals(XmlDeclarations.BONDSTEREO_EL, bondStereo.getLocalName());
		String atomRefs4 = bondStereo.getAttributeValue(XmlDeclarations.ATOMREFS4_ATR);
		assertEquals("a1 a2 a3 a4", atomRefs4);
		assertEquals("a1 a2 a3 a4", atomRefs4);
		assertEquals(BondStereo.TRANS.toString(), bondStereo.getValue());
	}

	@Test
	public void applyStereochemistryUnlocantedZ() throws StructureBuildingException {
		Fragment f = n2s.parseChemicalName("(Z)-but-2-ene", false).getStructure();
		Atom atom2 = f.getAtomByLocant("2");
		Atom atom3 = f.getAtomByLocant("3");
		assertNotNull(atom2);
		assertNotNull(atom3);
		Bond chiralBond = f.findBond(atom2, atom3);
		assertNotNull(chiralBond);
		Element bondStereo = chiralBond.getBondStereoElement();
		assertNotNull(bondStereo);
		assertEquals(XmlDeclarations.BONDSTEREO_EL, bondStereo.getLocalName());
		String atomRefs4 = bondStereo.getAttributeValue(XmlDeclarations.ATOMREFS4_ATR);
		assertEquals("a1 a2 a3 a4", atomRefs4);
		assertEquals("a1 a2 a3 a4", atomRefs4);
		assertEquals(BondStereo.CIS.toString(), bondStereo.getValue());
	}
	
	@Test
	public void applyStereochemistryUnlocantedE() throws StructureBuildingException {
		Fragment f = n2s.parseChemicalName("(E)-but-2-ene", false).getStructure();
		Atom atom2 = f.getAtomByLocant("2");
		Atom atom3 = f.getAtomByLocant("3");
		assertNotNull(atom2);
		assertNotNull(atom3);
		Bond chiralBond = f.findBond(atom2, atom3);
		assertNotNull(chiralBond);
		Element bondStereo = chiralBond.getBondStereoElement();
		assertNotNull(bondStereo);
		assertEquals(XmlDeclarations.BONDSTEREO_EL, bondStereo.getLocalName());
		String atomRefs4 = bondStereo.getAttributeValue(XmlDeclarations.ATOMREFS4_ATR);
		assertEquals("a1 a2 a3 a4", atomRefs4);
		assertEquals("a1 a2 a3 a4", atomRefs4);
		assertEquals(BondStereo.TRANS.toString(), bondStereo.getValue());
	}
	
	@Test
	public void applyStereochemistryCis() throws StructureBuildingException {
		Fragment f = n2s.parseChemicalName("cis-but-2-ene", false).getStructure();
		Atom atom2 = f.getAtomByLocant("2");
		Atom atom3 = f.getAtomByLocant("3");
		assertNotNull(atom2);
		assertNotNull(atom3);
		Bond chiralBond = f.findBond(atom2, atom3);
		assertNotNull(chiralBond);
		Element bondStereo = chiralBond.getBondStereoElement();
		assertNotNull(bondStereo);
		assertEquals(XmlDeclarations.BONDSTEREO_EL, bondStereo.getLocalName());
		String atomRefs4 = bondStereo.getAttributeValue(XmlDeclarations.ATOMREFS4_ATR);
		assertEquals("a1 a2 a3 a4", atomRefs4);
		assertEquals("a1 a2 a3 a4", atomRefs4);
		assertEquals(BondStereo.CIS.toString(), bondStereo.getValue());
	}
	
	@Test
	public void applyStereochemistryTrans() throws StructureBuildingException {
		Fragment f = n2s.parseChemicalName("trans-but-2-ene", false).getStructure();
		Atom atom2 = f.getAtomByLocant("2");
		Atom atom3 = f.getAtomByLocant("3");
		assertNotNull(atom2);
		assertNotNull(atom3);
		Bond chiralBond = f.findBond(atom2, atom3);
		assertNotNull(chiralBond);
		Element bondStereo = chiralBond.getBondStereoElement();
		assertNotNull(bondStereo);
		assertEquals(XmlDeclarations.BONDSTEREO_EL, bondStereo.getLocalName());
		String atomRefs4 = bondStereo.getAttributeValue(XmlDeclarations.ATOMREFS4_ATR);
		assertEquals("a1 a2 a3 a4", atomRefs4);
		assertEquals("a1 a2 a3 a4", atomRefs4);
		assertEquals(BondStereo.TRANS.toString(), bondStereo.getValue());
	}
	
	
	@Test
	public void applyStereochemistryLocantedRS() throws StructureBuildingException {
		Fragment f = n2s.parseChemicalName("(1S,2R)-2-(methylamino)-1-phenylpropan-1-ol", false).getStructure();
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
}
