package uk.ac.cam.ch.wwmm.opsin;

import java.util.List;

import junit.framework.Assert;

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
		Assert.assertEquals(1, stereoAnalyser.findStereoCentres().size());
		Assert.assertEquals(0, stereoAnalyser.findStereoBonds().size());
		StereoCentre sc = stereoAnalyser.findStereoCentres().get(0);
		Assert.assertNotNull(sc.getStereoAtom());
		Atom stereoAtom = sc.getStereoAtom();
		Assert.assertEquals("C", stereoAtom.getElement());
		Assert.assertEquals(4, stereoAtom.getID());
	}
	
	@Test
	public void findStereoCentresNacetylleucine() throws StructureBuildingException {
		Fragment f = n2s.parseChemicalName("N-acetylleucine", false).getStructure();
		StereoAnalyser stereoAnalyser = new StereoAnalyser(f);
		Assert.assertEquals(1, stereoAnalyser.findStereoCentres().size());
		Assert.assertEquals(0, stereoAnalyser.findStereoBonds().size());
		StereoCentre sc = stereoAnalyser.findStereoCentres().get(0);
		Assert.assertNotNull(sc.getStereoAtom());
		Atom stereoAtom = sc.getStereoAtom();
		Assert.assertEquals("C", stereoAtom.getElement());
		List<Atom> neighbours = sc.getCipOrderedAtoms();
		for (int i = 0; i < neighbours.size(); i++) {
			Atom a = neighbours.get(i);
			if (i==0){
				Assert.assertEquals(a.getElement(), "H");
			}
			else if (i==1){
				Assert.assertEquals(a.getElement(), "C");
			}
			else if (i==2){
				Assert.assertEquals(a.getElement(), "C");
			}
			else if (i==3){
				Assert.assertEquals(a.getElement(), "N");
			}
		}
	}
	
	@Test
	public void findStereoCentresBut2ene() throws StructureBuildingException {
		Fragment f = n2s.parseChemicalName("but-2-ene", false).getStructure();
		StereoAnalyser stereoAnalyser = new StereoAnalyser(f);
		Assert.assertEquals(0, stereoAnalyser.findStereoCentres().size());
		Assert.assertEquals(1, stereoAnalyser.findStereoBonds().size());
		StereoBond sb = stereoAnalyser.findStereoBonds().get(0);
		Bond stereoBond = sb.getBond();
		Assert.assertNotNull(stereoBond);
		Atom stereoAtom1 = stereoBond.getFromAtom();
		Atom stereoAtom2 = stereoBond.getToAtom();
		Assert.assertNotNull(stereoAtom1);
		Assert.assertNotNull(stereoAtom2);
		Assert.assertNotSame(stereoAtom1, stereoAtom2);
		if (stereoAtom1.getID() == 2){
			Assert.assertEquals(3, stereoAtom2.getID());
		}
		else{
			Assert.assertEquals(2, stereoAtom2.getID());
			Assert.assertEquals(3, stereoAtom1.getID());
		}
	}
	
	/*
	 * Tests for applying stereochemistry
	 */
	
	@Test
	public void applyStereochemistryLocantedZBut2ene() throws StructureBuildingException {
		Fragment f = n2s.parseChemicalName("(2Z)-but-2-ene", false).getStructure();
		Atom atom2 = f.getAtomByLocant("2");
		Atom atom3 = f.getAtomByLocant("3");
		Assert.assertNotNull(atom2);
		Assert.assertNotNull(atom3);
		Bond chiralBond = f.findBond(atom2, atom3);
		Assert.assertNotNull(chiralBond);
		Element bondStereo = chiralBond.getBondStereoElement();
		Assert.assertNotNull(bondStereo);
		Assert.assertEquals(XmlDeclarations.BONDSTEREO_EL, bondStereo.getLocalName());
		String atomRefs4 = bondStereo.getAttributeValue(XmlDeclarations.ATOMREFS4_ATR);
		Assert.assertEquals("a1 a2 a3 a4", atomRefs4);
		Assert.assertEquals("a1 a2 a3 a4", atomRefs4);
		Assert.assertEquals(BondStereo.CIS.toString(), bondStereo.getValue());
	}

	@Test
	public void applyStereochemistryZBut2ene() throws StructureBuildingException {
		Fragment f = n2s.parseChemicalName("(Z)-but-2-ene", false).getStructure();
		Atom atom2 = f.getAtomByLocant("2");
		Atom atom3 = f.getAtomByLocant("3");
		Assert.assertNotNull(atom2);
		Assert.assertNotNull(atom3);
		Bond chiralBond = f.findBond(atom2, atom3);
		Assert.assertNotNull(chiralBond);
		Element bondStereo = chiralBond.getBondStereoElement();
		Assert.assertNotNull(bondStereo);
		Assert.assertEquals(XmlDeclarations.BONDSTEREO_EL, bondStereo.getLocalName());
		String atomRefs4 = bondStereo.getAttributeValue(XmlDeclarations.ATOMREFS4_ATR);
		Assert.assertEquals("a1 a2 a3 a4", atomRefs4);
		Assert.assertEquals("a1 a2 a3 a4", atomRefs4);
		Assert.assertEquals(BondStereo.CIS.toString(), bondStereo.getValue());
	}
	
	@Test
	public void applyStereochemistryTransBut2ene() throws StructureBuildingException {
		Fragment f = n2s.parseChemicalName("trans-but-2-ene", false).getStructure();
		Atom atom2 = f.getAtomByLocant("2");
		Atom atom3 = f.getAtomByLocant("3");
		Assert.assertNotNull(atom2);
		Assert.assertNotNull(atom3);
		Bond chiralBond = f.findBond(atom2, atom3);
		Assert.assertNotNull(chiralBond);
		Element bondStereo = chiralBond.getBondStereoElement();
		Assert.assertNotNull(bondStereo);
		Assert.assertEquals(XmlDeclarations.BONDSTEREO_EL, bondStereo.getLocalName());
		String atomRefs4 = bondStereo.getAttributeValue(XmlDeclarations.ATOMREFS4_ATR);
		Assert.assertEquals("a1 a2 a3 a4", atomRefs4);
		Assert.assertEquals("a1 a2 a3 a4", atomRefs4);
		Assert.assertEquals(BondStereo.TRANS.toString(), bondStereo.getValue());
	}
	
}
