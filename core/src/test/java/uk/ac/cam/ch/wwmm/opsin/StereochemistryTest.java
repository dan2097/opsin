package uk.ac.cam.ch.wwmm.opsin;

import java.util.List;

import junit.framework.Assert;

import org.junit.BeforeClass;
import org.junit.Test;

import sea36.chem.core.CMLAtom;
import sea36.chem.core.CMLBond;
import sea36.chem.stereo.StereoAnalyser;
import sea36.chem.stereo.StereoAnalyser.StereoAnalysis;
import sea36.chem.stereo.StereoAnalyser.StereoBond;
import sea36.chem.stereo.StereoAnalyser.StereoCentre;
import uk.ac.cam.ch.wwmm.opsin.Atom;
import uk.ac.cam.ch.wwmm.opsin.OpsinToChemKitWrapper;
import uk.ac.cam.ch.wwmm.opsin.Fragment;
import uk.ac.cam.ch.wwmm.opsin.NameToStructure;
public class StereochemistryTest {

	private static NameToStructure n2s;
	@BeforeClass
	public static void setup() throws Exception {
		n2s = NameToStructure.getInstance();
	}
	
	@Test
	public void bromoChloroFluoroMethane() {
		Fragment f = n2s.parseChemicalName("bromochlorofluoromethane", false).getStructure();
		OpsinToChemKitWrapper chemKitWrapper  =  new OpsinToChemKitWrapper(f);
		StereoAnalyser stereoAnalyser = new StereoAnalyser();
		stereoAnalyser.setIgnoreRecursiveStereoCentres(true);
		StereoAnalysis analysis = stereoAnalyser.findStereoCentres(chemKitWrapper.getChemKitMolecule());
		Assert.assertEquals(1, analysis.getStereoCentres().size());
		Assert.assertEquals(0, analysis.getStereoBonds().size());
		StereoCentre sc = analysis.getStereoCentres().get(0);
		Assert.assertNotNull(sc.getStereoAtom());
		Assert.assertNotNull(chemKitWrapper.getOpsinAtomFromChemKitAtom(sc.getStereoAtom()));
		Atom opsinAtom = chemKitWrapper.getOpsinAtomFromChemKitAtom(sc.getStereoAtom());
		Assert.assertEquals("C", opsinAtom.getElement());
		Assert.assertEquals(4, opsinAtom.getID());
	}
	
	@Test
	public void Nacetylleucine() throws StructureBuildingException {
		Fragment f = n2s.parseChemicalName("N-acetylleucine", false).getStructure();
		OpsinToChemKitWrapper chemKitWrapper  =  new OpsinToChemKitWrapper(f);
		StereoAnalyser stereoAnalyser = new StereoAnalyser();
		stereoAnalyser.setIgnoreRecursiveStereoCentres(true);
		StereoAnalysis analysis = stereoAnalyser.findStereoCentres(chemKitWrapper.getChemKitMolecule());
		Assert.assertEquals(1, analysis.getStereoCentres().size());
		Assert.assertEquals(0, analysis.getStereoBonds().size());
		StereoCentre sc = analysis.getStereoCentres().get(0);
		Assert.assertNotNull(sc.getStereoAtom());
		Assert.assertNotNull(chemKitWrapper.getOpsinAtomFromChemKitAtom(sc.getStereoAtom()));
		Atom opsinAtom = chemKitWrapper.getOpsinAtomFromChemKitAtom(sc.getStereoAtom());
		Assert.assertEquals("C", opsinAtom.getElement());
		List<Atom> neighbours = opsinAtom.getAtomNeighbours();
		int hydrogens =0;
		int carbons =0;
		int nitrogen =0;
		for (Atom atom : neighbours) {
			if (atom.getElement().equals("C")){
				carbons++;
			}
			else if (atom.getElement().equals("N")){
				nitrogen++;
			}
			else if (atom.getElement().equals("H")){
				hydrogens++;
			}
		}
		Assert.assertEquals(2, carbons);
		Assert.assertEquals(1, nitrogen);
		Assert.assertEquals(1, hydrogens);
	}
	
	@Test
	public void but2ene() {
		Fragment f = n2s.parseChemicalName("but-2-ene", false).getStructure();
		OpsinToChemKitWrapper chemKitWrapper  =  new OpsinToChemKitWrapper(f);
		StereoAnalyser stereoAnalyser = new StereoAnalyser();
		stereoAnalyser.setIgnoreRecursiveStereoCentres(true);
		StereoAnalysis analysis = stereoAnalyser.findStereoCentres(chemKitWrapper.getChemKitMolecule());
		Assert.assertEquals(0, analysis.getStereoCentres().size());
		Assert.assertEquals(1, analysis.getStereoBonds().size());
		StereoBond sb = analysis.getStereoBonds().get(0);
		CMLBond chemKitBond = sb.getBond();
		Assert.assertNotNull(chemKitBond);
		CMLAtom atom1 = chemKitBond.getAtom0();
		CMLAtom atom2 = chemKitBond.getAtom1();
		Assert.assertNotNull(atom1);
		Assert.assertNotNull(atom2);
		Atom opsinAtom1 = chemKitWrapper.getOpsinAtomFromChemKitAtom(atom1);
		Atom opsinAtom2 = chemKitWrapper.getOpsinAtomFromChemKitAtom(atom2);
		Assert.assertNotNull(opsinAtom1);
		Assert.assertNotNull(opsinAtom2);
		Assert.assertNotSame(opsinAtom1, opsinAtom2);
		if (opsinAtom1.getID() == 2){
			Assert.assertEquals(3, opsinAtom2.getID());
		}
		else{
			Assert.assertEquals(2, opsinAtom2.getID());
			Assert.assertEquals(3, opsinAtom1.getID());
		}
	}
}
