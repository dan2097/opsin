package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Deque;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

import uk.ac.cam.ch.wwmm.opsin.Bond.SMILES_BOND_DIRECTION;
import uk.ac.cam.ch.wwmm.opsin.BondStereo.BondStereoValue;

/**
 * Writes an isomeric SMILES serialisation of an OPSIN fragment
 * @author dl387
 *
 */
class SMILESWriter {

	/**The organic atoms and their allowed implicit valences in SMILES */
	private static final Map<ChemEl,Integer[]> organicAtomsToStandardValencies = new EnumMap<ChemEl, Integer[]>(ChemEl.class);

	/**Closures 1-9, %10-99, 0 */
	private static final  List<String> closureSymbols = new ArrayList<String>();


	/**The available ring closure symbols, ordered from start to end in the preferred order for use.*/
	private final Deque<String> availableClosureSymbols = new ArrayDeque<String>(closureSymbols);

	/**Maps between bonds and the ring closure to use when the atom that ends the bond is encountered.*/
	private final HashMap<Bond, String> bondToClosureSymbolMap = new HashMap<Bond, String>();

	/**Maps between bonds and the atom that this bond will go to in the SMILES. Populated in the order the bonds are to be made */
	private final HashMap<Bond, Atom> bondToNextAtomMap = new LinkedHashMap<Bond, Atom>();

	/**The structure to be converted to SMILES*/
	private final Fragment structure;

	/**Holds the SMILES string which is under construction*/
	private final StringBuilder smilesBuilder = new StringBuilder();
	
	/**Should extended SMILES be output*/
	private final boolean outputExtendedSmiles;
	
	/**The order atoms were traversed when creating the SMILES*/
	private List<Atom> smilesOutputOrder;

	static {
		organicAtomsToStandardValencies.put(ChemEl.B, new Integer[]{3});
		organicAtomsToStandardValencies.put(ChemEl.C, new Integer[]{4});
		organicAtomsToStandardValencies.put(ChemEl.N, new Integer[]{3,5});//note that OPSIN doesn't accept valency 5 nitrogen without the lambda convention
		organicAtomsToStandardValencies.put(ChemEl.O, new Integer[]{2});
		organicAtomsToStandardValencies.put(ChemEl.P, new Integer[]{3,5});
		organicAtomsToStandardValencies.put(ChemEl.S, new Integer[]{2,4,6});
		organicAtomsToStandardValencies.put(ChemEl.F, new Integer[]{1});
		organicAtomsToStandardValencies.put(ChemEl.Cl, new Integer[]{1});
		organicAtomsToStandardValencies.put(ChemEl.Br, new Integer[]{1});
		organicAtomsToStandardValencies.put(ChemEl.I, new Integer[]{1});
		
		organicAtomsToStandardValencies.put(ChemEl.R, new Integer[]{1,2,3,4,5,6,7,8,9});

		for (int i = 1; i <=9; i++) {
			closureSymbols.add(String.valueOf(i));
		}
		for (int i = 10; i <=99; i++) {
			closureSymbols.add("%"+i);
		}
		closureSymbols.add("0");
	}

	/**
	 * Creates a SMILES writer for the given fragment
	 * @param structure
	 * @param outputExtendedSmiles 
	 */
	private SMILESWriter(Fragment structure, boolean outputExtendedSmiles) {
		this.structure = structure;
		this.outputExtendedSmiles = outputExtendedSmiles;
	}

	/**
	 * Generates SMILES for the given fragment
	 * The following assumptions are currently made:
	 * 	The fragment contains no bonds to atoms outside the fragment
	 * 	Hydrogens are all explicit
	 * 	Spare valency has been converted to double bonds
	 * @return SMILES String
	 */
	static String generateSmiles(Fragment structure) {
		return new SMILESWriter(structure, false).writeSmiles();
	}

	/**
	 * Generates extended SMILES for the given fragment
	 * The following assumptions are currently made:
	 * 	The fragment contains no bonds to atoms outside the fragment
	 * 	Hydrogens are all explicit
	 * 	Spare valency has been converted to double bonds
	 * @return Extended SMILES String
	 */
	static String generateExtendedSmiles(Fragment structure) {
		return new SMILESWriter(structure, true).writeSmiles();
	}

	String writeSmiles() {
		assignSmilesOrder();
		assignDoubleBondStereochemistrySlashes();

		List<Atom> atomList = structure.getAtomList();
		smilesOutputOrder = new ArrayList<Atom>(atomList.size());

		boolean isEmpty = true;
		for (Atom currentAtom : atomList) {
			Integer visitedDepth = currentAtom.getProperty(Atom.VISITED);
			if (visitedDepth != null && visitedDepth ==0) {//new component
				if (!isEmpty){
					smilesBuilder.append('.');
				}
				traverseSmiles(currentAtom);
				isEmpty = false;
			}
		}
		
		if (outputExtendedSmiles) {
			writeExtendedSmilesLayer();
		}

		return smilesBuilder.toString();
	}

	private void writeExtendedSmilesLayer() {
		List<String> atomLabels = new ArrayList<String>();
		List<String> positionVariationBonds = new ArrayList<String>();
		Integer lastLabel = null;
		int attachmentPointCounter = 1;
		List<Atom> polymerAttachPoints = structure.getPolymerAttachmentPoints();
		boolean isPolymer = polymerAttachPoints != null && polymerAttachPoints.size() > 0;
		for (int i = 0, l = smilesOutputOrder.size(); i < l; i++) {
			Atom a = smilesOutputOrder.get(i);
			String homologyGroup = a.getProperty(Atom.HOMOLOGY_GROUP);
			if (homologyGroup != null) {
				homologyGroup = escapeExtendedSmilesLabel(homologyGroup);
				if (homologyGroup.startsWith("_")) {
					atomLabels.add(homologyGroup);
				}
				else {
					atomLabels.add(homologyGroup + "_p");
				}
				lastLabel = i;
			}
			else if (a.getElement() == ChemEl.R){
				if (isPolymer) {
					atomLabels.add("star_e");
				}
				else {
					atomLabels.add("_AP" + String.valueOf(attachmentPointCounter++));
				}
				lastLabel = i;
			}
			else {
				atomLabels.add("");
			}
			List<Atom> atomsInPositionVariationBond = a.getProperty(Atom.POSITION_VARIATION_BOND);
			if (atomsInPositionVariationBond != null) {
				StringBuilder sb = new StringBuilder();
				sb.append(i);
				for (int j = 0; j < atomsInPositionVariationBond.size(); j++) {
					sb.append(j==0 ? ':' : '.');
					Atom referencedAtom = atomsInPositionVariationBond.get(j);
					int referencedAtomIndex = smilesOutputOrder.indexOf(referencedAtom);
					if (referencedAtomIndex == -1){
						throw new RuntimeException("OPSIN Bug: Failed to resolve position variation bond atom");
					}
					sb.append(referencedAtomIndex);
				}
				positionVariationBonds.add(sb.toString());
			}
		}
		List<String> extendedSmiles = new ArrayList<String>(2);
		if (lastLabel != null) {
			extendedSmiles.add("$" + StringTools.stringListToString(atomLabels.subList(0, lastLabel + 1), ";") + "$" );
		}
		if (positionVariationBonds.size() > 0) {
			extendedSmiles.add("m:" + StringTools.stringListToString(positionVariationBonds, ","));
		}
		if (isPolymer) {
			StringBuilder sruContents = new StringBuilder();
			sruContents.append("Sg:n:");
			boolean appendDelimiter = false;
			for (int i = 0, l = smilesOutputOrder.size(); i < l; i++) {
				if (smilesOutputOrder.get(i).getElement() != ChemEl.R) {
					if (appendDelimiter) {
						sruContents.append(',');
					}
					sruContents.append(i);
					appendDelimiter = true;
				}
			}
			sruContents.append("::ht");
			extendedSmiles.add(sruContents.toString());
		}
		if (extendedSmiles.size() > 0) {
			smilesBuilder.append(" |");
			smilesBuilder.append(StringTools.stringListToString(extendedSmiles, ","));
			smilesBuilder.append('|');
		}
	}
	
	private String escapeExtendedSmilesLabel(String str) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0, len = str.length(); i < len; i++) {
			char ch = str.charAt(i);
			if ((ch >= 'a' && ch <= 'z') || 
			   (ch >= 'A' && ch <= 'Z')  ||
			   (ch >= '0' && ch <= '9') ) {
				sb.append(ch);
			}
			else {
				sb.append("&#");
				sb.append(String.valueOf((int)ch));
				sb.append(';');
			}
		}
		return sb.toString();
	}

	/**
	 * Walks through the fragment populating the Atom.VISITED property indicating how many bonds
	 * an atom is from the start of the fragment walk. A new walk will be started for each disconnected component of the fragment
	 */
	private void assignSmilesOrder() {
		List<Atom> atomList = structure.getAtomList();
		for (Atom atom : atomList) {
			atom.setProperty(Atom.VISITED, null);
		}
		for (Atom a : atomList) {
			if(a.getProperty(Atom.VISITED) == null && !isSmilesImplicitProton(a)){//true for only the first atom in a fully connected molecule
				traverseMolecule(a);
			}
		}
	}
	
	private static class TraversalState {
		private final Atom atom;
		private final Bond bondTaken;
		private final int depth;

		private TraversalState(Atom atom, Bond bondTaken, int depth) {
			this.atom = atom;
			this.bondTaken = bondTaken;
			this.depth = depth;
		}
	}
	
	/**
	 * Iterative function for populating the Atom.VISITED property 
	 * Also populates the bondToNextAtom Map
	 * @param startingAtom
	 * @return
	 */
	private void traverseMolecule(Atom startingAtom){
		Deque<TraversalState> stack = new ArrayDeque<TraversalState>();
		stack.add(new TraversalState(startingAtom, null, 0));
		while (!stack.isEmpty()){
			TraversalState currentstate = stack.removeLast();
			Atom currentAtom = currentstate.atom;
			Bond bondtaken = currentstate.bondTaken;
			if (bondtaken != null) {
				bondToNextAtomMap.put(bondtaken, currentAtom);
			}	
			if(currentAtom.getProperty(Atom.VISITED) != null){
				continue;
			}
			int depth = currentstate.depth;
			currentAtom.setProperty(Atom.VISITED, depth);
			List<Bond> bonds = currentAtom.getBonds();
			for (int i = bonds.size() - 1; i >=0; i--) {
				Bond bond = bonds.get(i);
				if (bond.equals(bondtaken)){
					continue;
				}
				Atom neighbour = bond.getOtherAtom(currentAtom);
				if (isSmilesImplicitProton(neighbour)){
					continue;
				}
				stack.add(new TraversalState(neighbour, bond, depth + 1));
			}
		}
	}

	private boolean isSmilesImplicitProton(Atom atom) {
		if (atom.getElement() != ChemEl.H){
			//not hydrogen
			return false;
		}
		if (atom.getIsotope() != null && atom.getIsotope() != 1){
			//deuterium/tritium
			return false;
		}
		List<Atom> neighbours = atom.getAtomNeighbours();
		int neighbourCount = neighbours.size();
		if (neighbourCount > 1){
			//bridging hydrogen
			return false;
		}
		if (neighbourCount == 0){
			//just a hydrogen atom
			return false;
		}
		
		Atom neighbour = neighbours.get(0);
		ChemEl chemEl = neighbour.getElement();
		if (chemEl == ChemEl.H || chemEl == ChemEl.R) {
			//only connects to hydrogen or an R-group
			return false;
		}
		if (chemEl == ChemEl.N){
			List<Bond> bondsFromNitrogen = neighbour.getBonds();
			if (bondsFromNitrogen.size() == 2){
				for (Bond bond : bondsFromNitrogen) {
					if (bond.getBondStereo() != null){
						//special case where hydrogen is connected to a nitrogen with imine double bond stereochemistry
						return false;
					}
				}
			}
		}
		return true;
	}

	/**
	 * Goes through the bonds with BondStereo in the order the are to be created in the SMILES
	 * The bondStereo is used to set whether the bonds to non-implicit hydrogens that are adjacent to this bond
	 * should be be represented by / or \ in the SMILES. If this method has already set the slash on some bonds
	 * e.g. in a conjugated system this must be taken into account when setting the next slashes so as to not
	 * create a contradictory double bond stereochemistry definition.
	 */
	private void assignDoubleBondStereochemistrySlashes() {
		Set<Bond> bonds = bondToNextAtomMap.keySet();
		for (Bond bond : bonds) {
			bond.setSmilesStereochemistry(null);
		}
		for (Bond bond : bonds) {
			BondStereo bondStereo =bond.getBondStereo();
			if (bondStereo!=null){
				Atom[] atomRefs4 = bondStereo.getAtomRefs4();
				Bond bond1 = atomRefs4[0].getBondToAtom(atomRefs4[1]);
				Bond bond2 = atomRefs4[2].getBondToAtom(atomRefs4[3]);
				if (bond1 ==null || bond2==null){
					throw new RuntimeException("OPSIN Bug: Bondstereo described atoms that are not bonded");
				}
				Atom bond1ToAtom = bondToNextAtomMap.get(bond1);
				Atom bond2ToAtom = bondToNextAtomMap.get(bond2);
				SMILES_BOND_DIRECTION bond1Slash = bond1.getSmilesStereochemistry();//null except in conjugated systems
				SMILES_BOND_DIRECTION bond2Slash = bond2.getSmilesStereochemistry();
				SMILES_BOND_DIRECTION bond1Direction = SMILES_BOND_DIRECTION.LSLASH;
				SMILES_BOND_DIRECTION bond2Direction = SMILES_BOND_DIRECTION.LSLASH;
				if (bondStereo.getBondStereoValue().equals(BondStereoValue.CIS)){
					bond2Direction = bond2Direction.equals(SMILES_BOND_DIRECTION.LSLASH) ? SMILES_BOND_DIRECTION.RSLASH : SMILES_BOND_DIRECTION.LSLASH;//flip the slash type to be used from \ to /
				}
				if (!bond1ToAtom.equals(atomRefs4[1])){
					bond1Direction = bond1Direction.equals(SMILES_BOND_DIRECTION.LSLASH) ? SMILES_BOND_DIRECTION.RSLASH : SMILES_BOND_DIRECTION.LSLASH;
				}
				if (!bond2ToAtom.equals(atomRefs4[3])){
					bond2Direction = bond2Direction.equals(SMILES_BOND_DIRECTION.LSLASH) ? SMILES_BOND_DIRECTION.RSLASH : SMILES_BOND_DIRECTION.LSLASH;
				}
				
				//One of the bonds may have already have a defined slash from a previous bond stereo. If so make sure that we don't change it.
				if (bond1Slash !=null && !bond1Slash.equals(bond1Direction)){
					bond1Direction = bond1Direction.equals(SMILES_BOND_DIRECTION.LSLASH) ? SMILES_BOND_DIRECTION.RSLASH : SMILES_BOND_DIRECTION.LSLASH;
					bond2Direction = bond2Direction.equals(SMILES_BOND_DIRECTION.LSLASH) ? SMILES_BOND_DIRECTION.RSLASH : SMILES_BOND_DIRECTION.LSLASH;
				}
				else if (bond2Slash !=null && !bond2Slash.equals(bond2Direction)){
					bond1Direction = bond1Direction.equals(SMILES_BOND_DIRECTION.LSLASH) ? SMILES_BOND_DIRECTION.RSLASH : SMILES_BOND_DIRECTION.LSLASH;
					bond2Direction = bond2Direction.equals(SMILES_BOND_DIRECTION.LSLASH) ? SMILES_BOND_DIRECTION.RSLASH : SMILES_BOND_DIRECTION.LSLASH;
				}
	
				//Also need to investigate the bonds which are implicitly set by the bondStereo
				//F   Cl
				// C=C
				//N   O
				//e.g. the bonds from the C-N and C-O (the higher priority atoms will always be used for bond1/2)
				Bond bond1Other =null;
				Bond bond2Other =null;
				SMILES_BOND_DIRECTION bond1OtherDirection =null;
				SMILES_BOND_DIRECTION bond2OtherDirection =null;
				
				List<Bond> bondsFrom2ndAtom = new ArrayList<Bond>(atomRefs4[1].getBonds());
				bondsFrom2ndAtom.remove(bond1);
				bondsFrom2ndAtom.remove(bond);
				if (bondsFrom2ndAtom.size()==1){//can be 0 for imines
					if (bondToNextAtomMap.containsKey(bondsFrom2ndAtom.get(0))){//ignore bonds to implicit hydrogen
						bond1Other = bondsFrom2ndAtom.get(0);
						bond1OtherDirection =  bond1Direction.equals(SMILES_BOND_DIRECTION.LSLASH) ? SMILES_BOND_DIRECTION.RSLASH : SMILES_BOND_DIRECTION.LSLASH;
						if (!bond1ToAtom.equals(atomRefs4[1])){
							bond1OtherDirection = bond1OtherDirection.equals(SMILES_BOND_DIRECTION.LSLASH) ? SMILES_BOND_DIRECTION.RSLASH : SMILES_BOND_DIRECTION.LSLASH;
						}
						if (!bondToNextAtomMap.get(bond1Other).equals(atomRefs4[1])){
							bond1OtherDirection = bond1OtherDirection.equals(SMILES_BOND_DIRECTION.LSLASH) ? SMILES_BOND_DIRECTION.RSLASH : SMILES_BOND_DIRECTION.LSLASH;
						}
					}
				}
				
				List<Bond> bondsFrom3rdAtom= new ArrayList<Bond>(atomRefs4[2].getBonds());
				bondsFrom3rdAtom.remove(bond2);
				bondsFrom3rdAtom.remove(bond);
				if (bondsFrom3rdAtom.size()==1){
					if (bondToNextAtomMap.containsKey(bondsFrom3rdAtom.get(0))){
						bond2Other = bondsFrom3rdAtom.get(0);
						bond2OtherDirection =  bond2Direction.equals(SMILES_BOND_DIRECTION.LSLASH) ? SMILES_BOND_DIRECTION.RSLASH : SMILES_BOND_DIRECTION.LSLASH;
						if (!bond2ToAtom.equals(atomRefs4[3])){
							bond2OtherDirection = bond2OtherDirection.equals(SMILES_BOND_DIRECTION.LSLASH) ? SMILES_BOND_DIRECTION.RSLASH : SMILES_BOND_DIRECTION.LSLASH;
						}
						if (!bondToNextAtomMap.get(bond2Other).equals(bond2Other.getOtherAtom(atomRefs4[2]))){
							bond2OtherDirection = bond2OtherDirection.equals(SMILES_BOND_DIRECTION.LSLASH) ? SMILES_BOND_DIRECTION.RSLASH : SMILES_BOND_DIRECTION.LSLASH;
						}
					}
				}
				
				//One of the bonds may have already have a defined slash from a previous bond stereo. If so make sure that we don't change it.
				if (bond1Other !=null && bond1Other.getSmilesStereochemistry() !=null && !bond1Other.getSmilesStereochemistry().equals(bond1OtherDirection)){
					bond1Direction = bond1Direction.equals(SMILES_BOND_DIRECTION.LSLASH) ? SMILES_BOND_DIRECTION.RSLASH : SMILES_BOND_DIRECTION.LSLASH;
					bond2Direction = bond2Direction.equals(SMILES_BOND_DIRECTION.LSLASH) ? SMILES_BOND_DIRECTION.RSLASH : SMILES_BOND_DIRECTION.LSLASH;
					bond1OtherDirection = bond1OtherDirection.equals(SMILES_BOND_DIRECTION.LSLASH) ? SMILES_BOND_DIRECTION.RSLASH : SMILES_BOND_DIRECTION.LSLASH;
					if (bond2Other!=null){
						bond2OtherDirection = bond2OtherDirection.equals(SMILES_BOND_DIRECTION.LSLASH) ? SMILES_BOND_DIRECTION.RSLASH : SMILES_BOND_DIRECTION.LSLASH;
					}
				}
				else if (bond2Other !=null && bond2Other.getSmilesStereochemistry() !=null && !bond2Other.getSmilesStereochemistry().equals(bond2OtherDirection)){
					bond1Direction = bond1Direction.equals(SMILES_BOND_DIRECTION.LSLASH) ? SMILES_BOND_DIRECTION.RSLASH : SMILES_BOND_DIRECTION.LSLASH;
					bond2Direction = bond2Direction.equals(SMILES_BOND_DIRECTION.LSLASH) ? SMILES_BOND_DIRECTION.RSLASH : SMILES_BOND_DIRECTION.LSLASH;
					bond2OtherDirection = bond2OtherDirection.equals(SMILES_BOND_DIRECTION.LSLASH) ? SMILES_BOND_DIRECTION.RSLASH : SMILES_BOND_DIRECTION.LSLASH;
					if (bond1Other!=null){
						bond1OtherDirection = bond1OtherDirection.equals(SMILES_BOND_DIRECTION.LSLASH) ? SMILES_BOND_DIRECTION.RSLASH : SMILES_BOND_DIRECTION.LSLASH;
					}
				}
				
				//Set slashes for all bonds that are not to implicit hydrogen
				//In non conjugated systems this will yield redundant, but consistent, information
				bond1.setSmilesStereochemistry(bond1Direction);
				bond2.setSmilesStereochemistry(bond2Direction);
	
				if (bond1Other!=null){
					bond1Other.setSmilesStereochemistry(bond1OtherDirection);
				}
				if (bond2Other!=null){
					bond2Other.setSmilesStereochemistry(bond2OtherDirection);
				}
			}
		}
	}

	
	private static final TraversalState startBranch = new TraversalState(null, null, -1);
	private static final TraversalState endBranch = new TraversalState(null, null, -1);
	
	/**
	 * Generates the SMILES starting from the currentAtom, iteratively exploring
	 * in the same order as {@link SMILESWriter#traverseMolecule(Atom)}
	 * @param startingAtom
	 */
	private void traverseSmiles(Atom startingAtom){
		Deque<TraversalState> stack = new ArrayDeque<TraversalState>();
		stack.add(new TraversalState(startingAtom, null, 0));
		while (!stack.isEmpty()){
			TraversalState currentstate = stack.removeLast();
			if (currentstate == startBranch){
				smilesBuilder.append('(');
				continue;
			}
			if (currentstate == endBranch){
				smilesBuilder.append(')');
				continue;
			}
			Atom currentAtom = currentstate.atom;
			Bond bondtaken = currentstate.bondTaken;
			if (bondtaken != null){
				smilesBuilder.append(bondToSmiles(bondtaken));
			}
			int depth = currentstate.depth;

			smilesBuilder.append(atomToSmiles(currentAtom, depth, bondtaken));
			smilesOutputOrder.add(currentAtom);
			List<Bond> bonds = currentAtom.getBonds();
			List<String> newlyAvailableClosureSymbols = null;
			for (Bond bond : bonds) {//ring closures
				if (bond.equals(bondtaken)) {
					continue;
				}
				Atom neighbour = bond.getOtherAtom(currentAtom);
				Integer nDepth = neighbour.getProperty(Atom.VISITED);
				if (nDepth != null && nDepth <= depth){
					String closure = bondToClosureSymbolMap.get(bond);
					smilesBuilder.append(closure);
					if (newlyAvailableClosureSymbols == null){
						newlyAvailableClosureSymbols = new ArrayList<String>();
					}
					newlyAvailableClosureSymbols.add(closure);
				}
			}
			for (Bond bond : bonds) {//ring openings
				Atom neighbour = bond.getOtherAtom(currentAtom);
				Integer nDepth = neighbour.getProperty(Atom.VISITED);
				if (nDepth != null && nDepth > (depth +1)){
					String closure = availableClosureSymbols.removeFirst();
					bondToClosureSymbolMap.put(bond, closure);
					smilesBuilder.append(bondToSmiles(bond));
					smilesBuilder.append(closure);
				}
			}

			if (newlyAvailableClosureSymbols != null) {
				//By not immediately adding to availableClosureSymbols we avoid using the same digit 
				//to both close and open on the same atom
				for (int i = newlyAvailableClosureSymbols.size() -1; i >=0; i--) {
					availableClosureSymbols.addFirst(newlyAvailableClosureSymbols.get(i));
				}
			}
	
			boolean seenFirstBranch = false;
			for (int i = bonds.size() - 1; i >=0; i--) {
				//adjacent atoms which have not been previously written
				Bond bond = bonds.get(i);
				Atom neighbour = bond.getOtherAtom(currentAtom);
				Integer nDepth = neighbour.getProperty(Atom.VISITED);
				if (nDepth != null && nDepth == depth + 1){
					if (!seenFirstBranch){
						stack.add(new TraversalState(neighbour, bond, depth + 1));
						seenFirstBranch = true;
					}
					else {
						stack.add(endBranch);
						stack.add(new TraversalState(neighbour, bond, depth + 1));
						stack.add(startBranch);
					}
				}
			}
		}
	}

	/**
	 * Returns the SMILES describing the given atom.
	 * Where possible square brackets are not included to give more readable SMILES
	 * @param atom
	 * @param depth
	 * @param bondtaken
	 * @return
	 */
	private String atomToSmiles(Atom atom, int depth, Bond bondtaken) {
		StringBuilder atomSmiles = new StringBuilder();
		int hydrogenCount = calculateNumberOfBondedExplicitHydrogen(atom);
		boolean needsSquareBrackets = determineWhetherAtomNeedsSquareBrackets(atom, hydrogenCount);
		if (needsSquareBrackets) {
			atomSmiles.append('[');
		}
		if (atom.getIsotope() != null) {
			atomSmiles.append(atom.getIsotope());
		}
		ChemEl chemEl = atom.getElement();
		if (chemEl == ChemEl.R) {//used for polymers
			atomSmiles.append('*');
		}
		else{
			if (atom.hasSpareValency()) {//spare valency corresponds directly to lower case SMILES in OPSIN's SMILES reader
				atomSmiles.append(chemEl.toString().toLowerCase(Locale.ROOT));
			}
			else{
				atomSmiles.append(chemEl.toString());
			}
		}
		if (atom.getAtomParity() != null){
			atomSmiles.append(atomParityToSmiles(atom, depth, bondtaken));
		}
		if (hydrogenCount != 0 && needsSquareBrackets && chemEl != ChemEl.H){
			atomSmiles.append('H');
			if (hydrogenCount != 1){
				atomSmiles.append(String.valueOf(hydrogenCount));
			}
		}
		int charge = atom.getCharge();
	    if (charge != 0){
	    	if (charge == 1){
	    		atomSmiles.append('+');
	    	}
	    	else if (charge == -1){
	    		atomSmiles.append('-');
	    	}
	    	else{
	    		if (charge > 0){
	    			atomSmiles.append('+');
	    		}
	    		atomSmiles.append(charge);
	    	}
	    }
	    if (needsSquareBrackets) {
	    	Integer atomClass = atom.getProperty(Atom.ATOM_CLASS);
			if (atomClass != null) {
				atomSmiles.append(':');
				atomSmiles.append(String.valueOf(atomClass));
			}
	    	atomSmiles.append(']');
	    }
		return atomSmiles.toString();
	}

	private int calculateNumberOfBondedExplicitHydrogen(Atom atom) {
		List<Atom> neighbours = atom.getAtomNeighbours();
		int count = 0;
		for (Atom neighbour : neighbours) {
			if (neighbour.getProperty(Atom.VISITED) == null){
				count++;
			}
		}
		return count;
	}

	private boolean determineWhetherAtomNeedsSquareBrackets(Atom atom, int hydrogenCount) {
		Integer[] expectedValencies = organicAtomsToStandardValencies.get(atom.getElement());
		if (expectedValencies == null){
			return true;
		}
		if (atom.getCharge() != 0){
			return true;
		}
		if (atom.getIsotope() != null){
			return true;
		}
		if (atom.getAtomParity() != null){
			return true;
		}

		int valency = atom.getIncomingValency();
		boolean valencyCanBeDescribedImplicitly = Arrays.binarySearch(expectedValencies, valency) >= 0;
		int targetImplicitValency =valency;
		if (valency > expectedValencies[expectedValencies.length-1]){
			valencyCanBeDescribedImplicitly = true;
		}
		if (!valencyCanBeDescribedImplicitly){
			return true;
		}
	
		int nonHydrogenValency = valency - hydrogenCount;
		int implicitValencyThatWouldBeGenerated = nonHydrogenValency;
		for (int i = expectedValencies.length - 1; i >= 0; i--) {
			if (expectedValencies[i] >= nonHydrogenValency){
				implicitValencyThatWouldBeGenerated =expectedValencies[i];
			}
		}
		if (targetImplicitValency != implicitValencyThatWouldBeGenerated){
			return true;
		}
		if (atom.getProperty(Atom.ATOM_CLASS) != null) {
			return true;
		}
		return false;
	}

	private String atomParityToSmiles(Atom currentAtom, int depth, Bond bondtaken) {
		AtomParity atomParity = currentAtom.getAtomParity();
		Atom[] atomRefs4 = atomParity.getAtomRefs4().clone();

		List<Atom> atomrefs4Current = new ArrayList<Atom>();

		if (bondtaken != null) {//previous atom
			Atom neighbour = bondtaken.getOtherAtom(currentAtom);
			atomrefs4Current.add(neighbour);
		}

		for (Atom atom : atomRefs4) {//lone pair as in tetrahedral sulfones
			if (atom.equals(currentAtom)){
				atomrefs4Current.add(currentAtom);
			}
		}

		List<Bond> bonds = currentAtom.getBonds();
		for (Bond bond : bonds) {//implicit hydrogen
			Atom neighbour = bond.getOtherAtom(currentAtom);
			if (neighbour.getProperty(Atom.VISITED) == null){
				atomrefs4Current.add(currentAtom);
			}
		}
		for (Bond bond : bonds) {//ring closures
			if (bond.equals(bondtaken)){
				continue;
			}
			Atom neighbour = bond.getOtherAtom(currentAtom);
			if (neighbour.getProperty(Atom.VISITED) == null){
				continue;
			}
			if (neighbour.getProperty(Atom.VISITED) <= depth){
				atomrefs4Current.add(neighbour);
			}
		}
		for (Bond bond : bonds) {//ring openings
			Atom neighbour = bond.getOtherAtom(currentAtom);
			if (neighbour.getProperty(Atom.VISITED) == null){
				continue;
			}
			if (neighbour.getProperty(Atom.VISITED) > (depth +1)){
				atomrefs4Current.add(neighbour);
			}

		}
		for (Bond bond : bonds) {//next atom/s
			Atom neighbour = bond.getOtherAtom(currentAtom);
			if (neighbour.getProperty(Atom.VISITED) == null){
				continue;
			}
			if (neighbour.getProperty(Atom.VISITED) == depth + 1){
				atomrefs4Current.add(neighbour);
			}
		}
		Atom[] atomrefs4CurrentArr = new Atom[4];
		for (int i = 0; i < atomrefs4Current.size(); i++) {
			atomrefs4CurrentArr[i] = atomrefs4Current.get(i);
		}
		for (int i = 0; i < atomRefs4.length; i++) {//replace mentions of explicit hydrogen with the central atom the hydrogens are attached to, to be consistent with the SMILES representation
			if (atomRefs4[i].getProperty(Atom.VISITED) == null){
				atomRefs4[i] = currentAtom;
			}
		}

		boolean equivalent = StereochemistryHandler.checkEquivalencyOfAtomsRefs4AndParity(atomRefs4, atomParity.getParity(), atomrefs4CurrentArr, 1);
		if (equivalent){
			return "@@";
		}
		else{
			return "@";
		}
	}

	/**
	 * Generates the SMILES description of the bond
	 * In the case of cis/trans stereochemistry this relies on the {@link SMILESWriter#assignDoubleBondStereochemistrySlashes}
	 * having been run to setup the smilesBondDirection attribute
	 * @param bond
	 * @return
	 */
	private String bondToSmiles(Bond bond){
		String bondSmiles = "";
		int bondOrder = bond.getOrder();
		if (bondOrder == 2){
			bondSmiles = "=";
		}
		else if (bondOrder == 3){
			bondSmiles = "#";
		}
		else if (bond.getSmilesStereochemistry() != null){
			if (bond.getSmilesStereochemistry() == SMILES_BOND_DIRECTION.RSLASH){
				bondSmiles ="/";
			}
			else{
				bondSmiles ="\\";
			}
		}
		return bondSmiles;
	}
}
