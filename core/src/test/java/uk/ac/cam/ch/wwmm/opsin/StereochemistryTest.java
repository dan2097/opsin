package uk.ac.cam.ch.wwmm.opsin;

import java.util.List;

import junit.framework.Assert;

import org.junit.BeforeClass;
import org.junit.Test;

import uk.ac.cam.ch.wwmm.opsin.Atom;
import uk.ac.cam.ch.wwmm.opsin.Fragment;
import uk.ac.cam.ch.wwmm.opsin.NameToStructure;
import uk.ac.cam.ch.wwmm.opsin.StereoAnalyser.StereoBond;
import uk.ac.cam.ch.wwmm.opsin.StereoAnalyser.StereoCentre;
public class StereochemistryTest {

	private static NameToStructure n2s;
	@BeforeClass
	public static void setup() throws Exception {
		n2s = NameToStructure.getInstance();
	}
	
	@Test
	public void bromoChloroFluoroMethane() throws StructureBuildingException {
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
	public void Nacetylleucine() throws StructureBuildingException {
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
	public void but2ene() throws StructureBuildingException {
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
}
