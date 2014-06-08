package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.log4j.Logger;

/**
 * Numbers fusedRings
 * @author aa593
 * @author dl387
 *
 */
class FusedRingNumberer {
	private static final Logger LOG = Logger.getLogger(FusedRingNumberer.class);
	private static class RingConnectivityTable {
		final List<RingShape> ringShapes = new ArrayList<RingShape>();
		final List<Ring> neighbouringRings = new ArrayList<Ring>();
		final List<Integer> directionFromRingToNeighbouringRing = new ArrayList<Integer>();
		final List<Ring> usedRings = new ArrayList<Ring>();

		RingConnectivityTable copy(){
			RingConnectivityTable copy = new RingConnectivityTable();
			copy.ringShapes.addAll(ringShapes);
			copy.neighbouringRings.addAll(neighbouringRings);
			copy.directionFromRingToNeighbouringRing.addAll(directionFromRingToNeighbouringRing);
			copy.usedRings.addAll(usedRings);
			return copy;
		}
	}

	/**
	 * Wrapper for a ring of a fused ring system with the shape that ring is currently being treated as having
	 * @author dl387
	 *
	 */
	private static class RingShape{
		private final Ring ring;
		private final FusionRingShape shape;
		public RingShape(Ring ring, FusionRingShape shape) {
			this.ring = ring;
			this.shape = shape;
		}
		Ring getRing() {
			return ring;
		}
		FusionRingShape getShape() {
			return shape;
		}
	}

	enum FusionRingShape{
		enterFromLeftHouse,//5 membered ring
		enterFromTopLeftHouse,//5 membered ring
		enterFromTopRightHouse,//5 membered ring
		enterFromRightHouse,//5 membered ring
		enterFromLeftSevenMembered,//7 membered ring, modified from 6 membered at bottom
		enterFromRightSevenMembered,//7 membered ring, modified from 6 membered at top
		standard
	}

	private static class Chain {
		private final int length;
		private final int startingX;
		private final int y;

		Chain(int length, int startingX, int y) {
			this.length = length;
			this.startingX = startingX;
			this.y = y;
		}

		int getLength() {
			return length;
		}
		int getStartingX() {
			return startingX;
		}
		int getY() {
			return y;
		}
	}

	/**
	 * Sorts by atomSequences by the IUPAC rules for determining the preferred labelling
	 * The most preferred will be sorted to the back (0th position)
	 * @author dl387
	 *
	 */
	private static class SortAtomSequences implements Comparator<List<Atom>> {

	    public int compare(List<Atom> sequenceA, List<Atom> sequenceB){
	    	if (sequenceA.size() != sequenceB.size()){
	    		//Error in fused ring building. Identified ring sequences not the same lengths!
	    		return 0;
	    	}

	    	int i=0;
	    	int j=0;
	    	//Give low numbers for the heteroatoms as a set.
	    	while(i < sequenceA.size()){
				Atom atomA=sequenceA.get(i);
				boolean isAaHeteroatom = atomA.getElement() != ChemEl.C;


				//bridgehead carbon do not increment numbering
				if (!isAaHeteroatom && atomA.getIncomingValency()>=3){
					i++;
					continue;
				}

				Atom atomB=sequenceB.get(j);
				boolean isBaHeteroatom =atomB.getElement() != ChemEl.C;
				if (!isBaHeteroatom && atomB.getIncomingValency()>=3){
					j++;
					continue;
				}

				if (isAaHeteroatom && !isBaHeteroatom){
					return -1;
				}
				if (isBaHeteroatom && !isAaHeteroatom){
					return 1;
				}
	    		i++;j++;
	    	}

	    	i=0;
	    	j=0;
	    	//Give low numbers for heteroatoms when considered in the order: O, S, Se, Te, N, P, As, Sb, Bi, Si, Ge, Sn, Pb, B, Hg
	    	while(i < sequenceA.size()){
				Atom atomA=sequenceA.get(i);

				//bridgehead carbon do not increment numbering
				if (atomA.getElement() == ChemEl.C && atomA.getIncomingValency()>=3){
					i++;
					continue;
				}

				Atom atomB=sequenceB.get(j);
				if (atomB.getElement() == ChemEl.C && atomB.getIncomingValency()>=3){
					j++;
					continue;
				}

				Integer heteroAtomPriorityA = heteroAtomValues.get(atomA.getElement());
				int atomAElementValue = heteroAtomPriorityA != null ? heteroAtomPriorityA : 0;
				
				Integer heteroAtomPriorityB = heteroAtomValues.get(atomB.getElement());
				int atomBElementValue = heteroAtomPriorityB != null ? heteroAtomPriorityB : 0;

				if (atomAElementValue > atomBElementValue){
					return -1;
				}
				if (atomAElementValue < atomBElementValue){
					return 1;
				}
				i++;j++;
	    	}

	    	//Give low numbers to fusion carbon atoms.
	    	for ( i = 0; i < sequenceA.size(); i++) {
				Atom atomA=sequenceA.get(i);
				Atom atomB=sequenceB.get(i);
				if (atomA.getIncomingValency()>=3 && atomA.getElement() == ChemEl.C){
					if (!(atomB.getIncomingValency()>=3 && atomB.getElement() == ChemEl.C)){
						return -1;
					}
				}
				if (atomB.getIncomingValency()>=3 && atomB.getElement() == ChemEl.C){
					if (!(atomA.getIncomingValency()>=3 && atomA.getElement() == ChemEl.C)){
						return 1;
					}
				}
			}
	    	//Note that any sequences still unsorted at this step will have fusion carbon atoms in the same places
	    	//which means you can go through both sequences without constantly looking for fusion carbons i.e. the variable j is no longer needed

	    	//Give low numbers to fusion rather than non-fusion atoms of the same heteroelement.
	    	for (i = 0; i < sequenceA.size(); i++) {
				Atom atomA=sequenceA.get(i);
				Atom atomB=sequenceB.get(i);
				if (atomA.getIncomingValency()>=3){
					if (!(atomB.getIncomingValency()>=3)){
						return -1;
					}
				}
				if (atomB.getIncomingValency()>=3){
					if (!(atomA.getIncomingValency()>=3)){
						return 1;
					}
				}
			}
	    	//TODO consider heteroatoms FR5.4d
	    	return 0;
	    }
	}

	private static final Map<ChemEl, Integer> heteroAtomValues = new EnumMap<ChemEl, Integer>(ChemEl.class);
	static{
		//unknown heteroatoms or carbon are given a value of 0
		heteroAtomValues.put(ChemEl.Hg, 2);
		heteroAtomValues.put(ChemEl.Tl, 3);
		heteroAtomValues.put(ChemEl.In, 4);
		heteroAtomValues.put(ChemEl.Ga, 5);
		heteroAtomValues.put(ChemEl.Al, 6);
		heteroAtomValues.put(ChemEl.B, 7);
		heteroAtomValues.put(ChemEl.Pb, 8);
		heteroAtomValues.put(ChemEl.Sn, 9);
		heteroAtomValues.put(ChemEl.Ge, 10);
		heteroAtomValues.put(ChemEl.Si, 11);
		heteroAtomValues.put(ChemEl.Bi, 12);
		heteroAtomValues.put(ChemEl.Sb, 13);
		heteroAtomValues.put(ChemEl.As, 14);
		heteroAtomValues.put(ChemEl.P, 15);
		heteroAtomValues.put(ChemEl.N, 16);
		heteroAtomValues.put(ChemEl.Te, 17);
		heteroAtomValues.put(ChemEl.Se, 18);
		heteroAtomValues.put(ChemEl.S, 19);
		heteroAtomValues.put(ChemEl.O, 20);
		heteroAtomValues.put(ChemEl.I, 21);
		heteroAtomValues.put(ChemEl.Br, 22);
		heteroAtomValues.put(ChemEl.Cl, 23);
		heteroAtomValues.put(ChemEl.F, 24);
	}
	/*
	 * The meaning of the integers used is as follows:
	 *        2
	 *    3   ^  1
	 *      \ | /
	 * +-4 <-   -> 0
	 *      / | \
	 *   -3   v  -1
	 *       -2
	 *
	 * They indicate the relative directions between rings
	 * Possibly enums should be used...
	 */

	/**
	 * Numbers the fused ring
	 * Currently only works for a very limited selection of rings
	 * @param fusedRing
	 * @throws StructureBuildingException
	 */
	static void numberFusedRing(Fragment fusedRing) throws StructureBuildingException {
		List<Ring> rings = SSSRFinder.getSetOfSmallestRings(fusedRing);
		if (rings.size() <2){
			throw new StructureBuildingException("Ring perception system found less than 2 rings within input fragment!");
		}
		List<Atom> atomList = fusedRing.getAtomList();
		setupAdjacentFusedRingProperties(rings);
		if (!checkRingApplicability(rings)){
			for (Atom atom : atomList) {
				atom.clearLocants();
			}
			return;
		}
		List<List<Atom>> atomSequences = determinePossiblePeripheryAtomOrders(rings, atomList.size());
		if (atomSequences.size()==0){
			for (Atom atom : atomList) {
				atom.clearLocants();
			}
			return;
		}

		// add missing atoms to each path
		for (List<Atom> path : atomSequences) {//TODO properly support interior atom labelling
			for(Atom atom : atomList) {
				if(!path.contains(atom)) {
					path.add(atom);
				}
			}
		}
		// find the preferred numbering scheme then relabel with this scheme
		Collections.sort( atomSequences, new SortAtomSequences());
		fusedRing.setDefaultInAtom(atomSequences.get(0).get(0));
		FragmentTools.relabelLocantsAsFusedRingSystem(atomSequences.get(0));
		fusedRing.reorderAtomCollection(atomSequences.get(0));
	}

	/**
	 * Populates rings with their neighbouring fused rings and the bonds involved
	 * @param rings
	 */
	static void setupAdjacentFusedRingProperties(List<Ring> rings){
		for (int i = 0, l = rings.size(); i < l; i++) {
			Ring curRing  = rings.get(i);
			bondLoop : for (Bond bond : curRing.getBondList()) {		// go through all the bonds for the current ring
				for (int j = i + 1; j < l; j++) {
					Ring otherRing = rings.get(j);	
					if (otherRing.getBondList().contains(bond)) {		// check if this bond belongs to any other ring
						otherRing.addNeighbour(bond, curRing);
						curRing.addNeighbour(bond, otherRing);			// if so, then associate the bond with the adjacent ring
						continue bondLoop;
					}
				}
			}
		}
	}

	/**
	 * Checks that all the rings are of sizes 3-8 or if larger than 8 are involved in 2 or fewer fused bonds
	 * @param rings
	 * @return
	 */
	private static boolean checkRingApplicability(List<Ring> rings) {
		for (Ring ring : rings) {
			if (ring.size() <=2){
				throw new RuntimeException("Invalid ring size: " +ring.size());
			}
			if (ring.size() >8 && ring.getNumberOfFusedBonds() > 2){
				return false;
			}
		}
		return true;
	}

	/**
	 * Returns possible enumerations of atoms. Currently Interior atoms are not considered.
	 * These enumerations will be compliant with rules FR5.1-FR5.3 of the fused ring nomenclature guidelines
	 * http://www.chem.qmul.ac.uk/iupac/fusedring/FR51.html
	 * @param rings
	 * @param atomCountOfFusedRingSystem 
	 * @return
	 * @throws StructureBuildingException
	 */
	private static List<List<Atom>> determinePossiblePeripheryAtomOrders(List<Ring> rings, int atomCountOfFusedRingSystem) throws StructureBuildingException {
		List<Ring> tRings = findTerminalRings(rings);
		if (tRings.size()<1) {
			throw new RuntimeException("OPSIN bug: Unable to find a terminal ring in fused ring system");
		}
		Ring tRing = tRings.get(0);
		Bond b1 = getStartingNonFusedBond(tRing);
		if(b1 == null) {
			throw new RuntimeException("OPSIN Bug: Non-fused bond from terminal ring not found");
		}

		List<RingConnectivityTable> cts = new ArrayList<RingConnectivityTable>();
		RingConnectivityTable startingCT = new RingConnectivityTable();
		cts.add(startingCT);
		buildRingConnectionTables(tRing, null, 0, b1, b1.getFromAtom(), startingCT, cts);
		//The preference against fusion to elongated edges is built into the construction of the ring table
		
		/* FR 5.1.1/FR 5.1.2 Preferred shapes preferred to distorted shapes */
		removeCTsWithDistortedRingShapes(cts);
		//TODO better implement the corner cases of FR 5.1.3-5.1.5

		/* FR-5.2a. Maximum number of rings in a horizontal row */
		Map<RingConnectivityTable, List<Integer>> horizonalRowDirections = findLongestChainDirections(cts);
		List<Ring[][]> ringMaps = createRingMapsAlignedAlongGivenhorizonalRowDirections(horizonalRowDirections);
		/* FR-5.2b-d */
		return findPossiblePaths(ringMaps, atomCountOfFusedRingSystem);
	}

	/**
	 * Finds the rings with the minimum number of fused bonds
	 * @param rings
	 * @return
	 */
	private static List<Ring> findTerminalRings(List<Ring> rings) {
		List<Ring> tRings = new ArrayList<Ring>();

		int minFusedBonds = Integer.MAX_VALUE;
		for (Ring ring : rings){
			if (ring.getNumberOfFusedBonds() < minFusedBonds) {
				minFusedBonds = ring.getNumberOfFusedBonds();
			}
		}

		for  (Ring ring : rings){
			if (ring.getNumberOfFusedBonds() == minFusedBonds) {
				tRings.add(ring);
			}
		}
		return tRings;
	}

	/**
	 * Recursive function to create the connectivity table of the rings, for each connection includes both directions
	 * @param currentRing
	 * @param previousRing
	 * @param previousDir
	 * @param previousBond
	 * @param atom
	 * @param ct
	 * @param cts
	 * @return
	 */
	private static List<RingConnectivityTable> buildRingConnectionTables(Ring currentRing, Ring previousRing, int previousDir, Bond previousBond, Atom atom, RingConnectivityTable ct, List<RingConnectivityTable> cts) {
		// order atoms and bonds in the ring
		currentRing.makeCyclicLists(previousBond, atom);
		List<RingConnectivityTable> generatedCts = new ArrayList<RingConnectivityTable>();
		List<FusionRingShape> allowedShapes = getAllowedShapesForRing(currentRing, previousBond);
		if (allowedShapes.size()==0){
			throw new RuntimeException("OPSIN limitation, unsupported ring size in fused ring numbering");
		}
		ct.usedRings.add(currentRing);
		for (int i = allowedShapes.size()-1; i >=0; i--) {
			FusionRingShape fusionRingShape = allowedShapes.get(i);
			RingConnectivityTable currentCT;
			if (i==0){
				currentCT = ct;
			}
			else{
				currentCT =ct.copy();
				cts.add(currentCT);
				generatedCts.add(currentCT);
			}
			RingShape ringShape = new RingShape(currentRing, fusionRingShape);
			List<RingConnectivityTable> ctsToExpand = new ArrayList<RingConnectivityTable>();
			ctsToExpand.add(currentCT);//all the cts to consider, the currentCT and generated clones
			for (Ring neighbourRing : currentRing.getNeighbours()){
				//find the directions between the current ring and all neighbouring rings including the previous ring
				// this means that the direction to the previous ring will then be known in both directions

				// find direction
				Bond currentBond = findFusionBond(currentRing, neighbourRing);

				int dir = 0;
				if (neighbourRing == previousRing) {
					dir = getOppositeDirection(previousDir);
				}
				else {
					dir = calculateRingDirection(ringShape, previousBond, currentBond, previousDir);
				}
				//System.out.println(currentRing +"|" +neighbourRing +"|" +dir +"|" +(neighbourRing==previousRing));

				// place into connectivity table, like graph, rings and their connection
				for (RingConnectivityTable ctToExpand : ctsToExpand) {
					ctToExpand.ringShapes.add(ringShape);
					ctToExpand.neighbouringRings.add(neighbourRing);
					ctToExpand.directionFromRingToNeighbouringRing.add(dir);
				}
				if (!currentCT.usedRings.contains(neighbourRing)) {
					List<RingConnectivityTable> newCts = new ArrayList<RingConnectivityTable>();
					for (RingConnectivityTable ctToExpand : ctsToExpand) {
						Atom a = getAtomFromBond(currentRing, currentBond);
						List<RingConnectivityTable> generatedDownStreamCts = buildRingConnectionTables(neighbourRing, currentRing, dir, currentBond, a, ctToExpand, cts);
						newCts.addAll(generatedDownStreamCts);
					}
					ctsToExpand.addAll(newCts);
					generatedCts.addAll(newCts);
				}
			}
		}
		return generatedCts;
	}

	/**
	 * Returns the allowed shapes for the given ring.
	 * The starting bond is required to assured that elongated bonds do not unnecesarily correspond to fusions
	 * Currently only 5 membered rings are considered in multiple orientations but the same
	 * is probably required for 7+ member rings
	 * @param ring
	 * @param startingBond
	 * @return
	 */
	private static List<FusionRingShape> getAllowedShapesForRing(Ring ring, Bond startingBond) {
		List<FusionRingShape> allowedRingShapes = new ArrayList<FusionRingShape>();
		int size = ring.size();
		if (size==5){
			List<Bond> fusedBonds = ring.getFusedBonds();
			int fusedBondCount = fusedBonds.size();
			if (fusedBondCount==1){
				allowedRingShapes.add(FusionRingShape.enterFromLeftHouse);
			}
			else if (fusedBondCount==2 || fusedBondCount==3 || fusedBondCount==4){
				List<Integer> distances = new ArrayList<Integer>();//one distance is likely to be 0
				for (Bond fusedBond : fusedBonds) {
					distances.add(calculateDistanceBetweenBonds(startingBond, fusedBond, ring));
				}
				if (!distances.contains(1)){
					allowedRingShapes.add(FusionRingShape.enterFromLeftHouse);
				}
				if (!distances.contains(4)){
					allowedRingShapes.add(FusionRingShape.enterFromRightHouse);
				}

				if (!distances.contains(2)){
					allowedRingShapes.add(FusionRingShape.enterFromTopLeftHouse);
				}
				else if (!distances.contains(3)){
					allowedRingShapes.add(FusionRingShape.enterFromTopRightHouse);
				}
				allowedRingShapes = removeDegenerateRingShapes(allowedRingShapes, distances);
			}
			else if (fusedBondCount==5){
				allowedRingShapes.add(FusionRingShape.enterFromLeftHouse);
				allowedRingShapes.add(FusionRingShape.enterFromRightHouse);
				//top left and top right are the same other than position of the elongated bond which will invariably be used anyway
				allowedRingShapes.add(FusionRingShape.enterFromTopLeftHouse);
			}
		}
		else if (size==7){
			List<Bond> fusedBonds = ring.getFusedBonds();
			int fusedBondCount = fusedBonds.size();
			if (fusedBondCount==1){
				allowedRingShapes.add(FusionRingShape.enterFromLeftSevenMembered);
			}
			else{
				allowedRingShapes.add(FusionRingShape.enterFromLeftSevenMembered);
				allowedRingShapes.add(FusionRingShape.enterFromRightSevenMembered);
			}
		}
		else{
			allowedRingShapes.add(FusionRingShape.standard);
		}
		return allowedRingShapes;
	}

	/**
	 * Removes the ring shapes that for given distances have identical properties
	 * @param allowedRingShapes
	 * @param distances
	 */
	private static List<FusionRingShape> removeDegenerateRingShapes(List<FusionRingShape> allowedRingShapes, List<Integer> distances) {
		distances = new ArrayList<Integer>(distances);
		distances.remove((Integer)0);//remove distance 0 if present, this invariably comes from the starting bond and is not of interest (and breaks getDirectionFromDist)
		for (int i = allowedRingShapes.size() - 1; i >=0; i--) {
			FusionRingShape shapeToConsiderRemoving = allowedRingShapes.get(i);
			for (int j = i - 1; j >=0; j--) {
				FusionRingShape shapeToCompareWith = allowedRingShapes.get(j);
				boolean foundDifference = false;
				for (Integer distance : distances) {
					if (getDirectionFromDist(shapeToConsiderRemoving, 5, distance) != getDirectionFromDist(shapeToCompareWith, 5, distance)){
						foundDifference = true;
						break;
					}
				}
				if (!foundDifference){
					allowedRingShapes.remove(i);
					break;
				}
			}
		}

		return allowedRingShapes;
	}

	/**
	 * Calculates the direction of the next ring according to the distance between fusion bonds and the previous direction
	 * @param ringShape
	 * @param previousBond
	 * @param currentBond
	 * @param previousDir
	 * @return
	 */
	private static int calculateRingDirection(RingShape ringShape, Bond previousBond, Bond currentBond, int previousDir) {
		// take the ring fused to one from the previous loop step
		Ring ring = ringShape.getRing();
		if (ring.getCyclicBondList() == null ) {
			throw new RuntimeException("OPSIN bug: cyclic bond set should have already been populated");
		}

		int dist = calculateDistanceBetweenBonds(previousBond, currentBond, ring);

		if (dist == 0) {
			throw new RuntimeException("OPSIN bug: Distance between bonds is equal to 0");
		}

		int relativeDir = getDirectionFromDist(ringShape.getShape(), ring.size(), dist);
		return determineAbsoluteDirectionUsingPreviousDirection(ringShape.getShape(), ring.size(), relativeDir, previousDir);
	}

	/**
	 * Given two bonds on a ring returns the distance (in bonds) between them
	 * @param bond1
	 * @param bond2
	 * @param ring
	 * @return
	 */
	private static int calculateDistanceBetweenBonds(Bond bond1, Bond bond2, Ring ring) {
		List<Bond> cyclicBondList =ring.getCyclicBondList();
		int previousBondIndice = cyclicBondList.indexOf(bond1);
		int currentBondIndice = cyclicBondList.indexOf(bond2);
		if (previousBondIndice==-1 || currentBondIndice==-1){
			throw new RuntimeException("OPSIN bug: previous and current bond were not present in the cyclic bond list of the current ring");
		}
		int ringSize =ring.size();
		int dist = (ringSize + currentBondIndice - previousBondIndice) % ringSize;
		return dist;
	}

	/**
	 * Uses the ring shape, the ring size and distance between the incoming and outgoing fused bond to determine
	 * the relative direction between the entry point on the ring and the exit point
	 * @param fusionRingShape
	 * @param ringSize
	 * @param dist
	 * @return
	 */
	private static int getDirectionFromDist(FusionRingShape fusionRingShape, int ringSize, int dist) {
		int dir=0;
		if (ringSize == 3) { // 3 member ring
			if (dist == 1) {
				dir = -1;
			}
			else if (dist == 2) {
				dir = 1;
			}
			else throw new RuntimeException("Impossible distance between bonds for a 3 membered ring");
		}
		else if (ringSize == 4) { // 4 member ring
			if (dist == 2) {
				dir = 0;
			}
			else if (dist ==1) {
				dir = -2;
			}
			else if (dist ==3) {
				dir = 2;
			}
			else throw new RuntimeException("Impossible distance between bonds for a 4 membered ring");
		}
		else if (ringSize == 5) { // 5 member ring
			if (fusionRingShape == FusionRingShape.enterFromLeftHouse){
				if (dist ==1){
					dir = -2;//fusion to an elongated bond
				}
				else if (dist ==2){
					dir = 0;
				}
				else if (dist ==3){
					dir = 1;
				}
				else if (dist ==4){
					dir = 3;
				}
				else throw new RuntimeException("Impossible distance between bonds for a 5 membered ring");
			}
			else if (fusionRingShape == FusionRingShape.enterFromTopLeftHouse){
				if (dist ==1){
					dir = -3;
				}
				else if (dist ==2){
					dir = -1;//fusion to an elongated bond
				}
				else if (dist ==3){
					dir = 1;
				}
				else if (dist ==4){
					dir = 3;
				}
				else throw new RuntimeException("Impossible distance between bonds for a 5 membered ring");
			}
			else if (fusionRingShape == FusionRingShape.enterFromTopRightHouse){
				if (dist ==1){
					dir = -3;
				}
				else if (dist ==2){
					dir = -1;
				}
				else if (dist ==3){
					dir = 1;//fusion to an elongated bond
				}
				else if (dist ==4){
					dir = 3;
				}
				else throw new RuntimeException("Impossible distance between bonds for a 5 membered ring");
			}
			else if (fusionRingShape == FusionRingShape.enterFromRightHouse){
				if (dist ==1){
					dir = -3;
				}
				else if (dist ==2){
					dir = -1;
				}
				else if (dist ==3){
					dir = 0;
				}
				else if (dist ==4){
					dir = 2;//fusion to an elongated bond
				}
				else throw new RuntimeException("Impossible distance between bonds for a 5 membered ring");
			}
			else{
				throw new RuntimeException("OPSIN Bug: Unrecognised fusion ring shape for 5 membered ring");
			}
		}
		else if (ringSize == 7) { // 7 member ring
			if (fusionRingShape == FusionRingShape.enterFromLeftSevenMembered){
				if (dist ==1){
					dir = -3;//fusion to an abnormally angled bond
				}
				else if (dist ==2){
					dir = -2;
				}
				else if (dist ==3){
					dir = -1;//fusion to an abnormally angled bond
				}
				else if (dist ==4){
					dir = 0;
				}
				else if (dist ==5){
					dir = 1;
				}
				else if (dist ==6){
					dir = 3;
				}
				else throw new RuntimeException("Impossible distance between bonds for a 7 membered ring");
			}
			else if (fusionRingShape == FusionRingShape.enterFromRightSevenMembered){
				if (dist ==1){
					dir = -3;
				}
				else if (dist ==2){
					dir = -1;
				}
				else if (dist ==3){
					dir = 0;
				}
				else if (dist ==4){
					dir = 1;//fusion to an abnormally angled bond
				}
				else if (dist ==5){
					dir = 2;
				}
				else if (dist ==6){
					dir = 3;//fusion to an abnormally angled bond
				}
				else throw new RuntimeException("Impossible distance between bonds for a 7 membered ring");
			}
			else{
				throw new RuntimeException("OPSIN Bug: Unrecognised fusion ring shape for 7 membered ring");
			}
		}
		else if (ringSize % 2 == 0) {//general case even number of atoms ring (a 6 membered ring or distortion of)
			if (dist == 1) {
				dir = -3;
			}
			else if (dist == ringSize-1) {
				dir = 3;
			}
			else {
				dir = dist - ringSize/2;
				if (Math.abs(dir) > 2 && ringSize >= 8){// 8 and more neighbours
					dir = -2 * Integer.signum(dir);
				}
			}
		}
		else {// general case odd number of atoms ring (distortion of an even numbered ring by insertion of one atom).
			if (dist == 1) {
				dir = -3;
			}
			else if (dist == ringSize/2 || dist == ringSize/2 + 1) {//0 in both cases as effectively we are using a different depiction of the ring system. See FR-5.1.1 (this is done to give the longest horizontal row)
				dir = 0;
			}
			else if (dist == ringSize-1) {
				dir = 3;
			}
			else if(dist < ringSize/2) {
				dir = -2;
			}
			else if(dist > ringSize/2+1) {
				dir = 2;
			}
			else{
				throw new RuntimeException("OPSIN Bug: Unable to determine direction between odd number of atoms ring and next ring");
			}
		}
		return dir;
	}

	private static void removeCTsWithDistortedRingShapes(List<RingConnectivityTable> cts) {
		Map<RingConnectivityTable, List<Integer>> ctToDistortedRings = new HashMap<RingConnectivityTable, List<Integer>>();
		for (RingConnectivityTable ct : cts) {
			List<Integer> distortedRingSizes = new ArrayList<Integer>();
			ctToDistortedRings.put(ct, distortedRingSizes);
			List<RingShape> ringShapes = ct.ringShapes;
			for (int i = 0; i < ringShapes.size(); i++) {
				Ring r1 = ringShapes.get(i).getRing();
				Ring r2 = ct.neighbouringRings.get(i);
				for (int j = i +1; j < ringShapes.size(); j++) {
					if (ringShapes.get(j).getRing().equals(r2) && ct.neighbouringRings.get(j).equals(r1)){//look for the reverse entry in the ring connection table
						int expectedDir = getOppositeDirection(ct.directionFromRingToNeighbouringRing.get(i));
						if (expectedDir != ct.directionFromRingToNeighbouringRing.get(j)){
							distortedRingSizes.add(r2.size());
						}
					}
				}
			}
		}
		int minDistortedRings = Integer.MAX_VALUE;//find the minimum number of distorted rings
		for (List<Integer> distortedRingSizes : ctToDistortedRings.values()) {
			if (distortedRingSizes.size() < minDistortedRings){
				minDistortedRings = distortedRingSizes.size();
			}
		}
		for (int i = cts.size()-1; i>=0; i--) {
			if (ctToDistortedRings.get(cts.get(i)).size()>minDistortedRings){
				cts.remove(i);
			}
		}
	}

	/**
	 * Given a list of cts find the longest chain of rings in a line. This can be used a possible horizontal row
	 * The output is a map between the connection tables and the directions which give the longest chains
	 * Some cts may have no directions that give a chain of rings of this length
	 *
	 * @param cts
	 * @return
	 */
	private static Map<RingConnectivityTable, List<Integer>> findLongestChainDirections(List<RingConnectivityTable> cts){
		Map<RingConnectivityTable, List<Integer>> horizonalRowDirections = new LinkedHashMap<RingConnectivityTable, List<Integer>>();
		int maxChain = 0;
		for (RingConnectivityTable ct : cts) {
			if (ct.ringShapes.size() != ct.neighbouringRings.size() || ct.neighbouringRings.size() != ct.directionFromRingToNeighbouringRing.size()) {
				throw new RuntimeException("OPSIN Bug: Sizes of arrays in fused ring numbering connection table are not equal");
			}
			int ctEntriesSize =ct.ringShapes.size();
			List<Integer> directions = new  ArrayList<Integer>();
			horizonalRowDirections.put(ct, directions);
			for (int i=0; i< ctEntriesSize; i++){
				Ring neighbour = ct.neighbouringRings.get(i);
				int curChain = 1;
				int curDir = ct.directionFromRingToNeighbouringRing.get(i);

				nextRingInChainLoop: for (int k = 0; k <= ct.usedRings.size(); k++) {//<= rather than < so buggy behaviour can be caught
					int indexOfNeighbour = indexOfCorrespondingRingshape(ct.ringShapes, neighbour);

					if (indexOfNeighbour >= 0) {
						for (int j=indexOfNeighbour; j < ctEntriesSize; j++)	{
							if (ct.ringShapes.get(j).getRing() == neighbour && ct.directionFromRingToNeighbouringRing.get(j) == curDir) {
								curChain++;
								neighbour = ct.neighbouringRings.get(j);
								continue nextRingInChainLoop;
							}
						}
					}
					else{
						throw new RuntimeException("OPSIN bug: fused ring numbering: Ring missing from connection table");
					}
					if (curChain >= maxChain ) {
						int oDir = getOppositeDirection(curDir);
						if(curChain > maxChain){//new longest chain found
							for (List<Integer> previousDirections: horizonalRowDirections.values()) {
								previousDirections.clear();
							}
						}
						// if we has this direction before or its opposite, it is the same orientation
						if(curChain > maxChain || (!directions.contains(curDir) && !directions.contains(oDir))) {
							directions.add(curDir);
						}
						maxChain = curChain;
					}
					break;
				}
				if (maxChain > ct.usedRings.size()){
					throw new RuntimeException("OPSIN bug: fused ring layout contained a loop: more rings in a chain than there were rings!");
				}
			}
		}
		return horizonalRowDirections;
	}

	/**
	 * Given a list of ringShapes finds the indice of the ringShape corresponding to the given ring
	 * returns -1 if this is not possible
	 * @param ringShapes
	 * @param ring
	 * @return
	 */
	private static int indexOfCorrespondingRingshape(List<RingShape> ringShapes, Ring ring) {
		for (int i = 0; i < ringShapes.size(); i++) {
			if (ringShapes.get(i).getRing().equals(ring)){
				return i;
			}
		}
		return -1;
	}
	
	
	/**
	 * For each RingConnectivityTable and for each horizontal row direction creates a ringMap aligned along the given horizontal row direction
	 * @param horizonalRowDirectionsMap
	 * @return
	 * @throws StructureBuildingException
	 */
	private static List<Ring[][]> createRingMapsAlignedAlongGivenhorizonalRowDirections(Map<RingConnectivityTable, List<Integer>> horizonalRowDirectionsMap) throws StructureBuildingException {
		List<Ring[][]> ringMaps = new ArrayList<Ring[][]>();
		for (Entry<RingConnectivityTable, List<Integer>> entry : horizonalRowDirectionsMap.entrySet()) {
			RingConnectivityTable ct = entry.getKey();
			if ( ct.ringShapes.size() != ct.neighbouringRings.size() || ct.neighbouringRings.size() != ct.directionFromRingToNeighbouringRing.size() || ct.ringShapes.size() <= 0) {
				throw new RuntimeException("OPSIN Bug: Sizes of arrays in fused ring numbering connection table are not equal");
			}
			int ctEntriesSize = ct.ringShapes.size();
			for (Integer horizonalRowDirection : entry.getValue()) {
				int[] directionFromRingToNeighbouringRing = new int[ctEntriesSize];
				// turn the ring system such as to be aligned along the horizonalRowDirection
				for(int i=0; i<ctEntriesSize; i++){
					RingShape ringShape = ct.ringShapes.get(i);
					directionFromRingToNeighbouringRing[i] = determineAbsoluteDirectionUsingPreviousDirection(ringShape.getShape(), ringShape.getRing().size(), ct.directionFromRingToNeighbouringRing.get(i), -horizonalRowDirection);
				}
				Ring[][] ringMap = generateRingMap(ct, directionFromRingToNeighbouringRing);
				if (ringMap !=null){//null if overlapping bonds rings present
					ringMaps.add(ringMap);
				}
			}
		}
		if (ringMaps.size()==0){
			throw new StructureBuildingException("Fused ring systems with overlapping rings such as in helices cannot currently be numbered");
		}
		return ringMaps;
	}

	/**
	 * Applies FR5.2 B, C and D to determine the preferred orientation and returns lists of potential peripheral atom orderings
	 * @param ringMaps
	 * @param atomCountOfFusedRingSystem 
	 * @return
	 */
	private static List<List<Atom>> findPossiblePaths(List<Ring[][]> ringMaps, int atomCountOfFusedRingSystem){
		List<Double[]> chainQs = new ArrayList<Double[]>();
		List<Ring[][]> correspondingRingMap = new ArrayList<Ring[][]>();
		for (Ring[][] ringMap : ringMaps) {
			List<Chain> chains = findChainsOfMaximumLengthInHorizontalDir(ringMap);
			// For each chain count the number of rings in each quadrant
			for (Chain chain : chains) {
				int midChainXcoord = chain.getLength() + chain.getStartingX() - 1;//Remember the X axis is measured in 1/2s so don't need to 1/2 length

				Double[] qs = countQuadrants(ringMap, midChainXcoord, chain.getY());
				chainQs.add(qs);
				correspondingRingMap.add(ringMap);
			}
		}

		/*
		 * The quadrant numbers are as follows:
		 *
		 *  1  |  0
		 * ----+----
		 *  2  |  3
		 *
		 *  But at this stage it is not known what the mapping between these numbers and the/a preferred orientation of the structure is
		 */
		//  order for each right corner candidates for each chain
		List<List<Integer>> allowedUpperRightQuadrantsForEachChain =rulesBCD(chainQs);

		List<List<Atom>> paths = new ArrayList<List<Atom>> ();
		for (int c=0; c < chainQs.size(); c++) {
			Ring[][] ringMap = correspondingRingMap.get(c);
			List<Integer> allowedUpperRightQuadrants = allowedUpperRightQuadrantsForEachChain.get(c);

			for (Integer upperRightQuadrant : allowedUpperRightQuadrants) {
				Ring[][] qRingMap = transformQuadrantToUpperRightOfRingMap(ringMap, upperRightQuadrant);
				if (LOG.isTraceEnabled()){
					debugRingMap(qRingMap);
				}
				boolean inverseAtoms = (upperRightQuadrant == 2 || upperRightQuadrant == 0);
				List<Atom> peripheralAtomPath = orderAtoms(qRingMap, inverseAtoms, atomCountOfFusedRingSystem);
				paths.add(peripheralAtomPath);
			}
		}

		return paths;
	}

	private static Ring[][] generateRingMap(RingConnectivityTable ct, int[] directionFromRingToNeighbouringRing) {
		int ctEntriesSize = ct.ringShapes.size();
		// Find max and min coordinates for ringMap
		// we put the first ring into takenRings to start with it in the connection table
		int nRings = ct.usedRings.size();
		int[][] coordinates = new int[nRings][]; // correspondent to usedRings
		Ring[] takenRings = new Ring[nRings];
		int takenRingsCnt = 0;
		int maxX = 0;
		int minX = 0;
		int maxY = 0;
		int minY = 0;

		takenRings[takenRingsCnt++] = ct.ringShapes.get(0).getRing();
		coordinates[0] = new int[]{0,0};

		// Go through the rings in a system
		// Find the rings connected to them and assign coordinates according to the direction
		// Each time we go to the ring, whose coordinates were already identified.
		for(int tr=0; tr<nRings-1; tr++) {
			Ring currentRing = takenRings[tr];
			if (currentRing == null){
				throw new RuntimeException("OPSIN bug: Unexpected null ring in fused ring numbering");
			}

			int indexOfCurrentRing = indexOfCorrespondingRingshape(ct.ringShapes, currentRing);

			int xy[] = coordinates[tr]; // find the correspondent coordinates for the ring

			if (indexOfCurrentRing >= 0) {
				for (int j=indexOfCurrentRing; j< ctEntriesSize; j++) {
					if (ct.ringShapes.get(j).getRing() == currentRing) {
						Ring neighbour = ct.neighbouringRings.get(j);
						if (arrayContains(takenRings, neighbour)) {
							continue;
						}

						int[] newXY = new int[2];
						newXY[0] = xy[0] + Math.round(2 * countDX(directionFromRingToNeighbouringRing[j]));
						newXY[1] = xy[1] + countDY(directionFromRingToNeighbouringRing[j]);

						if(takenRingsCnt > takenRings.length) {
							throw new RuntimeException("OPSIN Bug: Fused ring numbering bug");
						}
						takenRings[takenRingsCnt] = neighbour;
						coordinates[takenRingsCnt] = newXY;
						takenRingsCnt++;

						if (newXY[0] > maxX){
							maxX = newXY[0];
						}
						else if (newXY[0] < minX) {
							minX = newXY[0];
						}

						if (newXY[1] > maxY){
							maxY = newXY[1];
						}
						else if (newXY[1] < minY) {
							minY = newXY[1];
						}
					}
				}
			}
			else{
				throw new RuntimeException("OPSIN bug: fused ring numbering: Ring missing from connection table");
			}
		}
		// the height and the width of the map
		int h = maxY - minY + 1;
		int w = maxX - minX + 1;

		Ring[][] ringMap = new Ring[w][h];

		// Map rings using coordinates calculated in the previous step, and transform them according to found minX and minY

		int ix = -minX;
		int iy = -minY;
		if (ix >= w || iy >= h) {
			throw new RuntimeException("OPSIN Bug: Fused ring numbering bug, Coordinates have been calculated wrongly");
		}

		int curX = 0;
		int curY = 0;
		for (int ti = 0; ti < takenRings.length; ti++){
			int[] xy = coordinates[ti];
			curX = xy[0] - minX;
			curY = xy[1] - minY;
			if(curX <0 || curX > w || curY < 0 || curY > h) {
				throw new RuntimeException("OPSIN Bug: Fused ring numbering bug, Coordinates have been calculated wrongly");
			}
			if (ringMap[curX][curY] != null){
				return null;
			}
			ringMap[curX][curY] = takenRings[ti];
		}
		return ringMap;
	}

	/**
	 * Finds all the chains of maximum length for the current direction
	 * @param ringMap
	 * @return
	 */
	private static List<Chain> findChainsOfMaximumLengthInHorizontalDir(Ring[][] ringMap){
		int w = ringMap.length;
		int h = ringMap[0].length;

		List<Chain> chains = new ArrayList<Chain>();

		int maxChain = 0;
		int chain = 0;

		// Find the longest chain
		for (int j=0; j<h; j++)	{
			for (int i=0; i<w; i++)	 {
				if(ringMap[i][j] != null) {
					chain = 1;
					while(i + 2*chain < w && ringMap[i + 2*chain][j] != null ) {
						chain++; // *2 because along the x axis the step is 2
					}
					if (chain > maxChain){
						chains.clear();
						maxChain = chain;
					}
					if(chain >= maxChain) {
						chains.add(new Chain(chain, i, j));
					}
					i += 2*chain;
				}
			}
		}
		return chains;
	}

	/**
	 * Counts number of rings in each quadrant
	 * @param ringMap
	 * @param midChainXcoord
	 * @param yChain
	 * @return
	 */
	private static Double[] countQuadrants(Ring[][] ringMap, int midChainXcoord, int yChain){
		Double[] qs = new Double[4];
		qs[0] = 0d;
		qs[1] = 0d;
		qs[2] = 0d;
		qs[3] = 0d;
		int w = ringMap.length;
		int h = ringMap[0].length;

		// Count rings in each quadrants
		for (int x=0; x<w; x++)	 {
			for (int y=0; y<h; y++)	{
				if (ringMap[x][y] == null) {
					continue;
				}

				if (x == midChainXcoord || y == yChain ) {// if the ring is on the axis
					if( x == midChainXcoord && y > yChain ) {
						qs[0]+=0.5;
						qs[1]+=0.5;
					}
					else if( x == midChainXcoord && y < yChain ) {
						qs[2]+=0.5;
						qs[3]+=0.5;
					}
					else if( x < midChainXcoord && y == yChain ) {
						qs[1]+=0.5;
						qs[2]+=0.5;
					}
					else if( x > midChainXcoord && y == yChain ) {
						qs[0]+=0.5;
						qs[3]+=0.5;
					}
					if (x==midChainXcoord && y==yChain ){
						qs[0]+=0.25;
						qs[1]+=0.25;
						qs[2]+=0.25;
						qs[3]+=0.25;
					}
				}
				else if(x > midChainXcoord && y > yChain) {
					qs[0]++;
				}
				else if(x < midChainXcoord && y > yChain) {
					qs[1]++;
				}
				else if(x < midChainXcoord && y < yChain) {
					qs[2]++;
				}
				else if(x > midChainXcoord && y < yChain) {
					qs[3]++;
				}
			}
		}

		return qs;
	}

	/**
	 * Applying rules FR5.2 B, C and D to the ring system.
	 * Return a list of possible upper right quadrants for each chain given. A chain may have multiple possible upper right quadrants (due to symmetry)
	 * or none if other chains can be shown to be preferable by application of the rules
	 * @param chainQs - array with number of ring in each quadrant for each chain.
	 */
	private static List<List<Integer>> rulesBCD(List<Double[]> chainQs) {
		List<List<Integer>> possibleUpperRightQuadrantsForEachChain = new ArrayList<List<Integer>>();
		int nChains = chainQs.size();
		if (nChains==0){
			throw new RuntimeException("OPSIN Bug: Fused ring numbering, no chains found?");
		}

		// Rule B: Maximum number of rings in upper right quadrant. Upper right corner candidates (it is not at this stage known which quadrant is the upper right one)
		double qmax = 0;

		for (Double[] chainQ : chainQs) {
			for (int j = 0; j < 4; j++)	{
				Double q = chainQ[j];
				if(q > qmax) {
					qmax = q;
				}
			}
		}

		for (Double[] chainQ : chainQs) {
			List<Integer> allowedUpperRightQuadrants = new ArrayList<Integer>();
			for (int j = 0; j < 4; j++){
				if (chainQ[j] == qmax) {
					allowedUpperRightQuadrants.add(j);
				}
			}
			possibleUpperRightQuadrantsForEachChain.add(allowedUpperRightQuadrants);
		}

		// Rule C: Minimum number of rings in lower left quadrant
		double qmin = Double.MAX_VALUE;

		for (int c = 0; c < nChains; c++) {
			List<Integer> possibleUpperRightQuadrant = possibleUpperRightQuadrantsForEachChain.get(c);
			for (Integer upperRightQuad : possibleUpperRightQuadrant) {
				int qdiagonal = (upperRightQuad + 2) % 4;
				if (chainQs.get(c)[qdiagonal] < qmin){
					qmin = chainQs.get(c)[qdiagonal];
				}
			}
		}
		for (int c = 0; c < nChains; c++) {
			List<Integer> possibleUpperRightQuadrant = possibleUpperRightQuadrantsForEachChain.get(c);
			List<Integer> allowedUpperRightQuadrants = new ArrayList<Integer>();
			for (Integer upperRightQuad : possibleUpperRightQuadrant) {
				int qdiagonal = (upperRightQuad + 2) % 4;
				if (chainQs.get(c)[qdiagonal]==qmin) {
					allowedUpperRightQuadrants.add(upperRightQuad);
				}
			}
			possibleUpperRightQuadrantsForEachChain.set(c, allowedUpperRightQuadrants);
		}

		// Rule D: Maximum number of rings above the horizontal row
		double rMax = 0;
		for (int c = 0; c < nChains; c++) {
			List<Integer> possibleUpperRightQuadrant = possibleUpperRightQuadrantsForEachChain.get(c);
			for (Integer upperRightQuad : possibleUpperRightQuadrant) {
				int upperLeftQuad;
				if (upperRightQuad % 2 == 0) {
					upperLeftQuad = upperRightQuad + 1;
				}
				else {
					upperLeftQuad = upperRightQuad - 1;
				}

				if (chainQs.get(c)[upperLeftQuad] + chainQs.get(c)[upperRightQuad] > rMax) {
					rMax = chainQs.get(c)[upperLeftQuad] + chainQs.get(c)[upperRightQuad];
				}
			}
		}
		for (int c = 0; c < nChains; c++) {
			List<Integer> possibleUpperRightQuadrant = possibleUpperRightQuadrantsForEachChain.get(c);
			List<Integer> allowedUpperRightQuadrants = new ArrayList<Integer>();
			for (Integer upperRightQuad : possibleUpperRightQuadrant) {
				int upperLeftQuad;
				if (upperRightQuad % 2 == 0) {
					upperLeftQuad = upperRightQuad + 1;
				}
				else {
					upperLeftQuad = upperRightQuad - 1;
				}

				if (chainQs.get(c)[upperLeftQuad] + chainQs.get(c)[upperRightQuad] == rMax) {
					allowedUpperRightQuadrants.add(upperRightQuad);
				}
			}
			possibleUpperRightQuadrantsForEachChain.set(c, allowedUpperRightQuadrants);
		}
		return possibleUpperRightQuadrantsForEachChain;
	}

	/**
	 * Enumerates the peripheral atoms in a system in accordance with FR-5.3:
	 * First finds the uppermost right ring, takes the next neighbour in the clockwise direction, and so on until the starting atom is reached
	 * @param ringMap
	 * @param inverseAtoms The direction in which the periphery atoms should be enumerated. Anticlockwise by default
	 * @param atomCountOfFusedRingSystem 
	 * @return
	 */
	private static List<Atom> orderAtoms(Ring[][] ringMap, boolean inverseAtoms, int atomCountOfFusedRingSystem){
		int w = ringMap.length;
		int h = ringMap[0].length;

		// find upper right ring
		Ring upperRightRing = null;
		for (int i=w-1; i>=0; i--) {
			if (ringMap[i][h-1] != null) {
				upperRightRing = ringMap[i][h-1];
				break;
			}
		}
		if (upperRightRing == null) {
			throw new RuntimeException("OPSIN Bug: Upper right ring not found when performing fused ring numbering");
		}
		List<Ring> visitedRings = new ArrayList<Ring>();
		visitedRings.add(upperRightRing);
		while (isEntirelyFusionAtoms(upperRightRing)){//c.f cyclopropa[de]anthracene
			upperRightRing = findClockwiseRingFromUpperRightRing(ringMap, upperRightRing, visitedRings);
			if (upperRightRing==null){
				throw new RuntimeException("OPSIN Bug: Unabled to find clockwise ring without fusion atoms");
			}
			visitedRings.add(upperRightRing);
		}

		Ring prevRing = findUpperLeftNeighbourOfUpperRightRing(ringMap, upperRightRing);
		Bond prevBond = findFusionBond(upperRightRing, prevRing);
		Bond nextBond = null;

		Ring currentRing = upperRightRing;
		Ring nextRing = null;
		List<Atom> atomPath = new ArrayList<Atom>();
		int count = 0;
		mainLoop: for (; count <= atomCountOfFusedRingSystem; count++) {
			int ringSize = currentRing.size();

			int startingBondIndex = currentRing.getBondIndex(prevBond) ;

			List<Bond> cyclicBonds = currentRing.getCyclicBondList();
			List<Bond> fusedBonds = currentRing.getFusedBonds();
			if (!inverseAtoms) {
				for(int bondIndex = 0; bondIndex < ringSize; bondIndex++) {
					int i = (startingBondIndex + bondIndex + 1) % ringSize; // +1 because we start from the bond next to stBond and end with it
					// if this bond is fused then it indicates the next ring to move to
					Bond bond = cyclicBonds.get(i);
					if(fusedBonds.contains(bond)) {
						nextBond = bond;
						break;
					}
				}
			}
			else {
				for(int bondIndex = 0; bondIndex < ringSize; bondIndex++) {
					int i = (startingBondIndex - bondIndex -1 + ringSize) % ringSize; // -1 because we start from the bond next to stBond and end with it
					// if this bond is fused then it indicates the next ring to move to
					Bond bond = cyclicBonds.get(i);
					if(fusedBonds.contains(bond)) {
						nextBond = bond;
						break;
					}
				}
			}
			if (nextBond == null) {
				throw new RuntimeException("OPSIN Bug: None of the bonds from this ring were fused, but this is not possible ");
			}

			// next ring
			nextRing = currentRing.getNeighbourOfFusedBond(nextBond);

			int endNumber = currentRing.getBondIndex(nextBond) ;

			// Add atoms in order, considering inverse or not inverse
			if (!inverseAtoms) {
				// if distance between prev bond and cur bond = 1 (it means that fused bonds are next to each other) i.e. come under interior atom numbering
				// we don't add that atom, cause it was added already
				if ( (endNumber - startingBondIndex + ringSize) % ringSize != 1) {
					startingBondIndex = (startingBondIndex + 1) % ringSize;
					endNumber = (endNumber - 1 + ringSize ) % ringSize;
					if (startingBondIndex > endNumber) {
						endNumber += ringSize;
					}

					// start from the atom next to fusion
					for (int j = startingBondIndex; j <= endNumber; j++) {
						Atom atom = currentRing.getCyclicAtomList().get(j % ringSize);
						if (atomPath.contains(atom)) {
							break mainLoop;
						}
						atomPath.add(atom);
					}
				}
			}
			else {
				if ( ( startingBondIndex - endNumber + ringSize) % ringSize != 1) {
					startingBondIndex = (startingBondIndex - 2 + ringSize ) % ringSize;
					endNumber = endNumber % ringSize;
					if (startingBondIndex < endNumber) {
						startingBondIndex += ringSize;
					}

					for (int j = startingBondIndex; j >= endNumber; j-- ) {
						Atom atom = currentRing.getCyclicAtomList().get(j % ringSize);
						if (atomPath.contains(atom)) {
							break mainLoop;
						}
						atomPath.add(atom);
					}
				}
			}
			prevBond = nextBond;
			prevRing = currentRing;
			currentRing = nextRing;
		}
		if (count ==atomCountOfFusedRingSystem){
			throw new RuntimeException("OPSIN Bug: Fused ring numbering may have been stuck in an infinite loop while enumerating peripheral numbering");
		}
		return atomPath;
	}

	private static boolean isEntirelyFusionAtoms(Ring upperRightRing) {
		List<Atom> atomList = upperRightRing.getAtomList();
		for (Atom atom : atomList) {
			if (atom.getBonds().size() < 3){
				return false;
			}
		}
		return true;
	}

	/**
	 * Finds the neighbour ring, which is the clockwise of the given ring.
	 * @param ringMap
	 * @param upperRightRing
	 * @param visitedRings 
	 * @return
	 */
	private static Ring findClockwiseRingFromUpperRightRing (Ring[][] ringMap, Ring upperRightRing, List<Ring> visitedRings){
		Ring clockwiseRing = null;
		int maxX = 0;
		int maxY = 0;

		for (Ring ring : upperRightRing.getNeighbours()) {
			if (visitedRings.contains(ring)){
				continue;
			}
			int xy[] = findRingPosition(ringMap, ring);
			if (xy==null) {
				throw new RuntimeException("OPSIN Bug: Ring not found in ringMap when performing fused ring numbering");
			}

			if (xy[0] > maxX  ||  xy[0] == maxX && xy[1] > maxY ) {
				maxX = xy[0];
				maxY = xy[1];
				clockwiseRing = ring;
			}
		}
		return clockwiseRing;
	}

	/**
	 * Finds the neighbour ring, which is the uppermost and on the left side from the given ring. Used to find previous bond for the uppermost right ring, from which we start to enumerate
	 * @param ringMap
	 * @param upperRightRing
	 * @return
	 */
	private static Ring findUpperLeftNeighbourOfUpperRightRing (Ring[][] ringMap, Ring upperRightRing){
		Ring nRing = null;
		int minX = Integer.MAX_VALUE;
		int maxY = 0;

		for (Ring ring : upperRightRing.getNeighbours()) {
			// upper left would be previous ring
			int xy[] = findRingPosition(ringMap, ring);
			if (xy==null) {
				throw new RuntimeException("OPSIN Bug: Ring not found in ringMap when performing fused ring numbering");
			}

			if (xy[1] > maxY  ||  xy[1] == maxY && xy[0] < minX ) {
				minX = xy[0];
				maxY = xy[1];
				nRing = ring;
			}
		}
		return nRing;
	}

	/**
	 * Finds the position(i,j) of the ring in the map
	 * @param ringMap
	 * @param ring
	 * @return
	 */
	private static int[] findRingPosition(Ring[][] ringMap, Ring ring) {
		int w = ringMap.length;
		int h = ringMap[0].length;

		for(int i=0; i<w; i++) {
			for(int j=0; j<h; j++) {
				if (ringMap[i][j] == ring) {
					return new int[]{i,j};
				}
			}
		}
		return null;
	}

	/**
	 * Transform the map such that the candidate upper right quadrant actually is in the upper right corner
	 * @param ringMap
	 * @param upperRightQuadrant
	 * @return
	 */
	private static Ring[][] transformQuadrantToUpperRightOfRingMap(Ring[][] ringMap, int upperRightQuadrant){
		int w = ringMap.length;
		int h = ringMap[0].length;

		Ring[][] rearrangedMap = new Ring[w][h];
		for (int i=0; i < w; i++) {
			for (int j=0; j < h; j++) {
				if (upperRightQuadrant == 0) {//already is in the upper right
					rearrangedMap[i][j] = ringMap[i][j];
				}
				if(upperRightQuadrant == 1) {//flip in y axis
					rearrangedMap[w-i-1][j] = ringMap[i][j];
				}
				else if(upperRightQuadrant == 2) {//flip in x and y axes
					rearrangedMap[w-i-1][h-j-1] = ringMap[i][j];
				}
				else if(upperRightQuadrant == 3) {//flip in x axis
					rearrangedMap[i][h-j-1] = ringMap[i][j];
				}
			}
		}

		return rearrangedMap;
	}

	/**
	 * Checks if array contains an object
	 * @param array
	 * @param obj
	 * @return
	 */
	private static boolean arrayContains(Object[] array, Object obj) {
		for (Object arrObj : array) {
			if (arrObj == obj) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Returns a bond which is not a bond that is in two rings
	 * Preference is given to a bond that is at least a bond away from a fused bond to avoid problems with 5 member rings starting in bad orientations
	 * @param tRing
	 * @return
	 */
	private static Bond getStartingNonFusedBond(Ring tRing){
		List<Bond> allBonds = new ArrayList<Bond>(tRing.getBondList());
		for (Bond fusedBond : tRing.getFusedBonds()) {
			List<Bond> neighbouringBonds = fusedBond.getFromAtom().getBonds();
			for (Bond bond : neighbouringBonds) {
				allBonds.remove(bond);
			}
			neighbouringBonds = fusedBond.getToAtom().getBonds();
			for (Bond bond : neighbouringBonds) {
				allBonds.remove(bond);
			}
		}
		if (allBonds.size() > 0){
			return allBonds.get(0);
		}
		for (Bond bond : tRing.getBondList()) {
			if(tRing.getNeighbourOfFusedBond(bond) == null){
				// return a non-fused bond
				return bond;
			}
		}
		return null;
	}

	/**
	 * Given the direction of the bond from ring1 to ring2, returns the opposite direction: from ring2 to ring1
	 * @param prevDir
	 * @return
	 */
	static int getOppositeDirection(int prevDir) {
		int dir;
		if (prevDir == 0) {
			dir = 4;
		}
		else if (Math.abs(prevDir) == 4){
			dir =0;
		}
		else if (Math.abs(prevDir) == 2){
			dir = 2 * -1 * Integer.signum(prevDir);
		}
		else if (Math.abs(prevDir) == 1){
			dir = 3 * -1 * Integer.signum(prevDir);
		}
		else {//prevDir will be +-3
			dir = 1 * -1 * Integer.signum(prevDir);
		}
		return dir;
	}

	/**
	 * Finds the atom connected to the bond, takes into account the order of the bonds and atoms in the ring
	 * @param ring
	 * @param curBond
	 * @return
	 */
	private static Atom getAtomFromBond(Ring ring, Bond curBond) {
		if (ring.getCyclicBondList() == null) {
			throw new RuntimeException("The cyclic bond list should already have been generated");
		}
		int bondIndice= ring.getCyclicBondList().indexOf(curBond);
		int atomIndice = ( bondIndice - 1 + ring.size() ) % ring.size();
		return ring.getCyclicAtomList().get(atomIndice);
	}

	/**
	 * Finds the fusion bond between 2 rings
	 * @param r1
	 * @param r2
	 * @return
	 */
	private static Bond findFusionBond (Ring r1, Ring r2) {
		List<Bond> b2 = r2.getBondList();
		for(Bond bond : r1.getBondList()){
			if (b2.contains(bond)) {
				return bond;
			}
		}
		return null;
	}

	/**
	 * Counts delta x distance between previous and next rings
	 * @param val
	 * @return
	 */
	private static float countDX (int val) {
		float dX = 0;
		if (Math.abs(val) == 1) {
			dX += 0.5f;
		}
		else if (Math.abs(val) == 3) {
			dX -= 0.5f;
		}
		else if (Math.abs(val) == 0) {
			dX += 1f;
		}
		else if (Math.abs(val) == 4) {
			dX -= 1f;
		}
		return dX;
	}

	/**
	 * Counts delta y distance (height) between previous and next rings
	 * @param val
	 * @return
	 */

	private static int countDY (int val) {
		int dY = 0;
		if (Math.abs(val) != 4) {
			if (val > 0) {
				dY = 1;
			}
			if (val < 0) {
				dY = -1;
			}
		}
		return dY;
	}

	/**
	 * Take into account the previous direction to convert the given relative direction into a direction that is absolute for the fused ring system
	 * @param fusionRingShape
	 * @param ringSize
	 * @param relativeDirection
	 * @param previousDir
	 * @return
	 */
	static int determineAbsoluteDirectionUsingPreviousDirection(FusionRingShape fusionRingShape, int ringSize, int relativeDirection, int previousDir){
		int interimDirection;
		if (Math.abs(previousDir) == 4) {
			if (relativeDirection == 0) {
				interimDirection = 4;
			}
			else {
				interimDirection = relativeDirection + 4 * -1 * Integer.signum(relativeDirection); // if dir<0 we add 4, if dir>0 we add -4
			}
		}
		else {
			interimDirection = relativeDirection + previousDir;
		}

		if (Math.abs(interimDirection)>4) {// Added
			interimDirection = (8 - Math.abs(interimDirection)) *  Integer.signum(interimDirection) * -1;
		}
		//TODO investigate this function and unit test
		 /* Even numbered rings when angled do not have direction 2.
		 * Almost true for 5 member except for corner case where fusion to elongated bond occurs
		 */
		if (Math.abs(interimDirection) == 2 && ((ringSize % 2 ==0) || ringSize==5)) {
			// if (one of them equal to 1 and another is equal to 3, we decrease absolute value and conserve the sign)
			if (Math.abs(relativeDirection)==1 && Math.abs(previousDir)==3  ||  Math.abs(relativeDirection)==3 && Math.abs(previousDir)==1) {
				interimDirection = 1 * Integer.signum(interimDirection);
			}
			// if both are equal to 1
			else if(Math.abs(relativeDirection)==1 && Math.abs(previousDir)==1 ) {
				interimDirection = 3 * Integer.signum(interimDirection);
			}
			// if both are equal to 3
			else if(Math.abs(relativeDirection)==3 && Math.abs(previousDir)==3 ) {
				interimDirection = 3 * Integer.signum(interimDirection);
			}
			// else it is correctly 2
		}

		if (interimDirection == -4) {
			interimDirection = 4;
		}

		return interimDirection;
	}

	private static void debugRingMap(Ring[][] ringMap) {
		Ring[][] yxOrdered = new Ring[ringMap[0].length][ringMap.length];
		for (int x = 0; x < ringMap.length; x++) {
			Ring[] yRings = ringMap[x];
			for (int y = 0; y < yRings.length; y++) {
				yxOrdered[y][x] =yRings[y];
			}
		}
		for (int y = yxOrdered.length-1; y >=0 ; y--) {
			Ring[] xRings = yxOrdered[y];
			StringBuilder sb = new StringBuilder();
			for (Ring ring : xRings) {
				if (ring!=null){
					int size = ring.size();
					if (size>9){
						if (size==10){
							sb.append("0");
						}
						else if (size % 2 ==0){
							sb.append("2");
						}
						else{
							sb.append("1");
						}
					}
					else{
						sb.append(size);
					}
				}
				else{
					sb.append(" ");
				}
			}
			LOG.trace(sb.toString());
		}
		LOG.trace("#########");

	}
}
