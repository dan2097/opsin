package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.Assert.*;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.junit.Test;

import uk.ac.cam.ch.wwmm.opsin.Fragment;

//Cycle detection is performed as part of fragment creation so we can just check the output of fragment creation
public class CycleDetectorTest {
	private SMILESFragmentBuilder sBuilder = new SMILESFragmentBuilder(new IDManager());
	
	@Test
	public void testAssignCyclic1() throws StructureBuildingException {
		Fragment frag = sBuilder.build("CCCC");
		for (Atom a : frag.getAtomList()) {
			assertEquals("Should be acylic", false, a.getAtomIsInACycle());
		}
	}

	@Test
	public void testAssignCyclic2() throws StructureBuildingException {
		Fragment frag = sBuilder.build("c1ccccc1");
		for (Atom a : frag.getAtomList()) {
			assertEquals("Should be cylic", true, a.getAtomIsInACycle());
		}
	}
	
	@Test
	public void testAssignCyclic3() throws StructureBuildingException {
		Fragment frag = sBuilder.build("c12.c23.c34.c45.c56.c61");
		for (Atom a : frag.getAtomList()) {
			assertEquals("Should be cylic", true, a.getAtomIsInACycle());
		}
	}
	
	@Test
	public void testAssignCyclic4() throws StructureBuildingException {
		Fragment frag = sBuilder.build("c1ccccc1CCc1ccccc1");
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
		Fragment frag = sBuilder.build("CCc1ccc(O)cc1");
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
		Fragment frag = sBuilder.build("CC1CC(O1)C");
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
	public void testFindPathBetweenAtoms1() throws StructureBuildingException {
		Fragment frag = sBuilder.build("c1ccccc1");
		List<Atom> atomList = frag.getAtomList();
		List<List<Atom>> paths = CycleDetector.getPathBetweenAtomsUsingBonds(atomList.get(0), atomList.get(3), frag.getBondSet());
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
	public void testFindPathBetweenAtoms2() throws StructureBuildingException {
		Fragment frag = sBuilder.build("C1CCCC2CCCCC12");
		List<Atom> atomList = frag.getAtomList();
		Set<Bond> bonds = new HashSet<Bond>(frag.getBondSet());
		bonds.remove(atomList.get(4).getBondToAtom(atomList.get(9)));
		List<List<Atom>> paths = CycleDetector.getPathBetweenAtomsUsingBonds(atomList.get(4), atomList.get(9), bonds);
		assertEquals(2, paths.size());
		
		List<Atom> pathLeftRing;
		List<Atom> pathRightRing;
		if (atomList.indexOf(paths.get(0).get(0))==3){
			pathLeftRing = paths.get(0);
			pathRightRing = paths.get(1);
		}
		else{
			pathLeftRing = paths.get(1);
			pathRightRing = paths.get(0);
		}
		assertEquals(3, atomList.indexOf(pathLeftRing.get(0)));
		assertEquals(2, atomList.indexOf(pathLeftRing.get(1)));
		assertEquals(1, atomList.indexOf(pathLeftRing.get(2)));
		assertEquals(0, atomList.indexOf(pathLeftRing.get(3)));
		
		assertEquals(5, atomList.indexOf(pathRightRing.get(0)));
		assertEquals(6, atomList.indexOf(pathRightRing.get(1)));
		assertEquals(7, atomList.indexOf(pathRightRing.get(2)));
		assertEquals(8, atomList.indexOf(pathRightRing.get(3)));
	}
	
	@Test
	public void testFindPathBetweenAtoms3() throws StructureBuildingException {
		Fragment frag = sBuilder.build("C1(C)CCCC2C(C)CCCC12");
		List<Atom> atomList = frag.getAtomList();
		Set<Bond> bonds = new HashSet<Bond>(frag.getBondSet());
		bonds.remove(atomList.get(0).getBondToAtom(atomList.get(1)));
		bonds.remove(atomList.get(6).getBondToAtom(atomList.get(7)));
		bonds.remove(atomList.get(5).getBondToAtom(atomList.get(11)));
		List<List<Atom>> paths = CycleDetector.getPathBetweenAtomsUsingBonds(atomList.get(0), atomList.get(6), bonds);
		assertEquals(2, paths.size());

		List<Atom> pathLeftRing;
		List<Atom> pathRightRing;
		if (atomList.indexOf(paths.get(0).get(0))==2){
			pathLeftRing = paths.get(0);
			pathRightRing = paths.get(1);
		}
		else{
			pathLeftRing = paths.get(1);
			pathRightRing = paths.get(0);
		}
		assertEquals(2, atomList.indexOf(pathLeftRing.get(0)));
		assertEquals(3, atomList.indexOf(pathLeftRing.get(1)));
		assertEquals(4, atomList.indexOf(pathLeftRing.get(2)));
		assertEquals(5, atomList.indexOf(pathLeftRing.get(3)));
		
		assertEquals(11, atomList.indexOf(pathRightRing.get(0)));
		assertEquals(10, atomList.indexOf(pathRightRing.get(1)));
		assertEquals(9, atomList.indexOf(pathRightRing.get(2)));
		assertEquals(8, atomList.indexOf(pathRightRing.get(3)));
	}
}
	
