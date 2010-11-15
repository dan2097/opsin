package uk.ac.cam.ch.wwmm.opsin;

import static junit.framework.Assert.assertEquals;
import static org.mockito.Mockito.mock;

import java.util.ArrayList;
import java.util.List;

import org.junit.BeforeClass;
import org.junit.Test;

import uk.ac.cam.ch.wwmm.opsin.Fragment;

//Cycle detection is performed as part of fragment creation so we can just check the output of fragment creation
public class CycleDetectorTest {
	private static FragmentManager fm;
	@BeforeClass
	public static void setup() throws Exception {
		fm = new FragmentManager(new SMILESFragmentBuilder(), mock(CMLFragmentBuilder.class), new IDManager());
	}
	
	@Test
	public void testAssignCyclic1() throws StructureBuildingException {
		Fragment frag = fm.buildSMILES("CCCC");
		for (Atom a : frag.getAtomList()) {
			assertEquals("Should be acylic", false, a.getAtomIsInACycle());
		}
	}
	
	
	@Test
	public void testAssignCyclic2() throws StructureBuildingException {
		Fragment frag = fm.buildSMILES("c1ccccc1");
		for (Atom a : frag.getAtomList()) {
			assertEquals("Should be cylic", true, a.getAtomIsInACycle());
		}
	}
	
	@Test
	public void testAssignCyclic3() throws StructureBuildingException {
		Fragment frag = fm.buildSMILES("c12.c23.c34.c45.c56.c61");
		for (Atom a : frag.getAtomList()) {
			assertEquals("Should be cylic", true, a.getAtomIsInACycle());
		}
	}
	
	@Test
	public void testAssignCyclic4() throws StructureBuildingException {
		Fragment frag = fm.buildSMILES("c1ccccc1CCc1ccccc1");
		List<Atom> atomList = frag.getAtomList();
		for (int i = 0; i < atomList.size(); i++) {
			Atom a = atomList.get(i);
			if (i<=5 || i >=8){
				assertEquals("Should be cylic", true, a.getAtomIsInACycle());
			}
			else{
				assertEquals("Should be acylic", false, a.getAtomIsInACycle());
			}
		}
	}
	
	@Test
	public void testAssignCyclic5() throws StructureBuildingException {
		Fragment frag = fm.buildSMILES("CCc1ccc(O)cc1");
		List<Atom> atomList = frag.getAtomList();
		for (int i = 0; i < atomList.size(); i++) {
			Atom a = atomList.get(i);
			if (i<=1 || i==6){
				assertEquals("Should be acylic", false, a.getAtomIsInACycle());
			}
			else{
				assertEquals("Should be cylic", true, a.getAtomIsInACycle());
			}
		}
	}
	
	@Test
	public void testAssignCyclic6() throws StructureBuildingException {
		Fragment frag = fm.buildSMILES("CC1CC(O1)C");
		List<Atom> atomList = frag.getAtomList();
		for (int i = 0; i < atomList.size(); i++) {
			Atom a = atomList.get(i);
			if (i==0 || i==5){
				assertEquals("Should be acylic", false, a.getAtomIsInACycle());
			}
			else{
				assertEquals("Should be cylic", true, a.getAtomIsInACycle());
			}
		}
	}
	
	@Test
	public void testFindIntraFragmentPaths1() throws StructureBuildingException {
		Fragment frag = fm.buildSMILES("c1ccccc1");
		List<Atom> atomList = frag.getAtomList();
		List<List<Atom>> paths = CycleDetector.getIntraFragmentPathsBetweenAtoms(atomList.get(0), atomList.get(3), frag);
		assertEquals(2, paths.size());
		for (List<Atom> path : paths) {
			assertEquals(2, path.size());
		}
		for (List<Atom> path : paths) {
			if (atomList.indexOf(path.get(0))==1){
				assertEquals(2, atomList.indexOf(path.get(1)));
			}
			else{
				assertEquals(5, atomList.indexOf(path.get(0)));
				assertEquals(4, atomList.indexOf(path.get(1)));
			}
		}
	}
	
	@Test
	public void testFindIntraFragmentPaths2() throws StructureBuildingException {
		Fragment frag = fm.buildSMILES("C1CCCC2CCCCC12");
		List<Atom> atomList = frag.getAtomList();
		List<List<Atom>> paths = CycleDetector.getIntraFragmentPathsBetweenAtoms(atomList.get(4), atomList.get(9), frag);
		assertEquals(3, paths.size());
		
		List<List<Atom>> nonZeroLengthPaths = new ArrayList<List<Atom>>(); 
		for (List<Atom> path : paths) {
			if (path.size()!=0){
				assertEquals(4, path.size());
				nonZeroLengthPaths.add(path);
			}
		}
		assertEquals(2, nonZeroLengthPaths.size());
		for (List<Atom> path : nonZeroLengthPaths) {
			if (atomList.indexOf(path.get(0))==5){
				assertEquals(6, atomList.indexOf(path.get(1)));
				assertEquals(7, atomList.indexOf(path.get(2)));
				assertEquals(8, atomList.indexOf(path.get(3)));
			}
			else{
				assertEquals(3, atomList.indexOf(path.get(0)));
				assertEquals(2, atomList.indexOf(path.get(1)));
				assertEquals(1, atomList.indexOf(path.get(2)));
				assertEquals(0, atomList.indexOf(path.get(3)));
			}
		}
	}
}
	
