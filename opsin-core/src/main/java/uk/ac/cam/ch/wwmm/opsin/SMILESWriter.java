package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
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
	private static final Map<String,Integer[]> organicAtomsToStandardValencies = new HashMap<String, Integer[]>();

	/**Closures 1-9, %10-99, 0 */
	private static final  LinkedList<String> closureSymbols = new LinkedList<String>();


	/**The available ring closure symbols, ordered from start to end in the preferred order for use.*/
	private final LinkedList<String> availableClosureSymbols = new LinkedList<String>(closureSymbols);

	/**Maps between bonds and the ring closure to use when the atom that ends the bond is encountered.*/
	private final HashMap<Bond, String> bondToClosureSymbolMap = new HashMap<Bond, String>();

	/**Maps between bonds and the atom that this bond will go to in the SMILES. Populated in the order the bonds are to be made */
	private final HashMap<Bond, Atom> bondToNextAtomMap= new LinkedHashMap<Bond, Atom>();

	/**The structure to be converted to SMILES*/
	private final Fragment structure;

	/**Holds the SMILES string which is under construction*/
	private final StringBuilder smilesBuilder = new StringBuilder();

	static {
		organicAtomsToStandardValencies.put("B", new Integer[]{3});
		organicAtomsToStandardValencies.put("C", new Integer[]{4});
		organicAtomsToStandardValencies.put("N", new Integer[]{3,5});//note that OPSIN doesn't accept valency 5 nitrogen without the lambda convention
		organicAtomsToStandardValencies.put("O", new Integer[]{2});
		organicAtomsToStandardValencies.put("P", new Integer[]{3,5});
		organicAtomsToStandardValencies.put("S", new Integer[]{2,4,6});
		organicAtomsToStandardValencies.put("F", new Integer[]{1});
		organicAtomsToStandardValencies.put("Cl", new Integer[]{1});
		organicAtomsToStandardValencies.put("Br", new Integer[]{1});
		organicAtomsToStandardValencies.put("I", new Integer[]{1});

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
	 */
	SMILESWriter(Fragment structure) {
		this.structure =structure;
	}

	/**
	 * Generates SMILES from the fragment the SMILESWriter was created with
	 * The following assumptions are currently made:
	 * 	The fragment contains no bonds to atoms outside the fragment
	 * 	Hydrogens are all explicit
	 * 	Spare valency has been converted to double bonds
	 * @return
	 */
	String generateSmiles() {
		assignSmilesOrder();
		assignDoubleBondStereochemistrySlashes();

		List<Atom> atomList = structure.getAtomList();
		List<Atom> nonProtonAtomList = createNonProtonAtomList(atomList);
		int nonProtonCount = nonProtonAtomList.size();
		boolean isEmpty =true;
		for (int i = 0; i < nonProtonCount; i++) {
			Atom currentAtom = nonProtonAtomList.get(i);
			if(currentAtom.getProperty(Atom.VISITED)==0){//new component
				if (!isEmpty){
					smilesBuilder.append('.');
				}
				traverseSmiles(currentAtom, null, 0);
				isEmpty =false;
			}
		}

		return smilesBuilder.toString();
	}

	/**
	 * Walks through the fragment populating the Atom.VISITED property indicating how many bonds
	 * an atom is from the start of the fragment walk. A new walk will be started for each disconnected component of the fragment
	 */
	private void assignSmilesOrder() {
		List<Atom> atomList =structure.getAtomList();
		for (Atom atom : atomList) {
			atom.setProperty(Atom.VISITED, null);
		}
		for (Atom a : atomList) {
			if(a.getProperty(Atom.VISITED)==null && !isSmilesImplicitProton(a)){//typically for all but the first atom this will be true
				traverseMolecule(a, null, 0);
			}
		}
	}

	/**
	 * Recursive function for populating the Atom.VISITED property 
	 * Also populates the bondToNextAtom Map
	 * @param currentAtom
	 * @param previousAtom
	 * @param depth
	 * @return
	 */
	private void traverseMolecule(Atom currentAtom, Atom previousAtom, int depth){
		if(currentAtom.getProperty(Atom.VISITED)!=null){
			return;
		}
		currentAtom.setProperty(Atom.VISITED, depth);
		Set<Bond> bonds = currentAtom.getBonds();
		for (Bond bond : bonds) {
			Atom neighbour = bond.getOtherAtom(currentAtom);
			if (isSmilesImplicitProton(neighbour)){
				continue;
			}
			if (neighbour.equals(previousAtom)){
				continue;
			}
			bondToNextAtomMap.put(bond, neighbour);
			traverseMolecule(neighbour, currentAtom, depth+1);
		}
	}

	private boolean isSmilesImplicitProton(Atom atom) {
		if (!atom.getElement().equals("H") || (atom.getIsotope()!=null && atom.getIsotope()!=1) ){
			return false;
		}
		else{
			List<Atom> neighbours = atom.getAtomNeighbours();
			//special case where hydrogen is bridging
			if (neighbours.size() > 1){
				return false;
			}
			//special case where hydrogen is a counter ion or only connects to other hydrogen
			boolean foundNonHydrogenNeighbour =false;
			for (Atom neighbour : neighbours) {
				if (!neighbour.getElement().equals("H")){
					foundNonHydrogenNeighbour =true;
				}
			}
			if (!foundNonHydrogenNeighbour){
				return false;
			}
			
			//special case where hydrogen is connected to a nitrogen with imine double bond stereochemistry
			if (neighbours.get(0).getElement().equals("N")){
				Set<Bond> bondsFromNitrogen = neighbours.get(0).getBonds();
				if (bondsFromNitrogen.size()==2){
					for (Bond bond : bondsFromNitrogen) {
						if (bond.getBondStereo()!=null){
							return false;
						}
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
					throw new RuntimeException("Bondstereo described atoms that are not bonded:" +bondStereo.toCML().toXML());
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

	private List<Atom> createNonProtonAtomList(List<Atom> atomList) {
		List<Atom> nonProtonAtomList = new ArrayList<Atom>();
		for (Atom atom : atomList) {
			if (atom.getProperty(Atom.VISITED)!=null){
				nonProtonAtomList.add(atom);
			}
		}
		return nonProtonAtomList;
	}

	/**
	 * Generates the SMILES for the currentAtom and its bonds and then is called recursively to explore the atom's neighbours
	 * @param currentAtom
	 * @param previousAtom
	 * @param depth
	 */
	private void traverseSmiles(Atom currentAtom, Atom previousAtom, int depth){
		smilesBuilder.append(atomToSmiles(currentAtom, depth, previousAtom));
		Set<Bond> bonds = currentAtom.getBonds();
		LinkedList<String> newlyAvailableClosureSymbols = null;
		for (Bond bond : bonds) {//ring closures
			Atom neighbour = bond.getOtherAtom(currentAtom);
			Integer nDepth = neighbour.getProperty(Atom.VISITED);
			if (nDepth!=null && nDepth<=depth && !neighbour.equals(previousAtom)){
				String closure = bondToClosureSymbolMap.get(bond);
				smilesBuilder.append(closure);
				if (newlyAvailableClosureSymbols == null){
					newlyAvailableClosureSymbols = new LinkedList<String>();
				}
				newlyAvailableClosureSymbols.addFirst(closure);
			}
		}
		for (Bond bond : bonds) {//ring openings
			Atom neighbour = bond.getOtherAtom(currentAtom);
			Integer nDepth = neighbour.getProperty(Atom.VISITED);
			if (nDepth!=null && nDepth > (depth +1)){
				String closure = availableClosureSymbols.removeFirst();
				bondToClosureSymbolMap.put(bond, closure);
				smilesBuilder.append(bondToSmiles(bond));
				smilesBuilder.append(closure);
			}
		}
		if (newlyAvailableClosureSymbols != null){
			for (String closure : newlyAvailableClosureSymbols) {
				availableClosureSymbols.addFirst(closure);
			}
		}

		// count outgoing edges
		int count = 0;
		for (Bond bond : bonds) {
			Atom neighbour = bond.getOtherAtom(currentAtom);
			Integer nDepth = neighbour.getProperty(Atom.VISITED);
			if (nDepth!=null && nDepth==depth+1){
				count++;
			}
		}
	
		for (Bond bond : bonds) {//adjacent atoms which have not been previously written
			Atom neighbour = bond.getOtherAtom(currentAtom);
			Integer nDepth = neighbour.getProperty(Atom.VISITED);
			if (nDepth!=null && nDepth==depth+1){
				if (count > 1){
				  smilesBuilder.append('(');
				}
				smilesBuilder.append(bondToSmiles(bond));
				traverseSmiles(neighbour,currentAtom,depth+1);
				if (count > 1){
					smilesBuilder.append(')');
					count--;
				}
			}
		}
	}

	/**
	 * Returns the SMILES describing the given atom.
	 * Where possible square brackets are not included to give more readable SMILES
	 * @param atom
	 * @param depth
	 * @param previousAtom
	 * @return
	 */
	private String atomToSmiles(Atom atom, int depth, Atom previousAtom) {
		StringBuilder atomSmiles = new StringBuilder();
		int hydrogen =calculateNumberOfBondedExplicitHydrogen(atom);
		boolean needsSquareBrackets = determineWhetherAtomNeedsSquareBrackets(atom, hydrogen);
		if (needsSquareBrackets){
			atomSmiles.append('[');
		}
		if (atom.getIsotope()!=null){
			atomSmiles.append(atom.getIsotope());
		}
		String elementSymbol =atom.getElement();
		if (atom.hasSpareValency()){//spare valency corresponds directly to lower case SMILES in OPSIN's SMILES reader
			atomSmiles.append(elementSymbol.toLowerCase());
		}
		else{
			if (elementSymbol.equals("R")){//used for polymers
				atomSmiles.append('*');
			}
			else{
				atomSmiles.append(elementSymbol);
			}
		}
		if (atom.getAtomParity()!=null){
			atomSmiles.append(atomParityToSmiles(atom, depth, previousAtom));
		}
		if (hydrogen !=0 && needsSquareBrackets && !elementSymbol.equals("H")){
			atomSmiles.append('H');
			if (hydrogen !=1){
				atomSmiles.append(String.valueOf(hydrogen));
			}
		}
		int charge = atom.getCharge();
	    if (charge != 0){
	    	if (charge==1){
	    		atomSmiles.append('+');
	    	}
	    	else if (charge==-1){
	    		atomSmiles.append('-');
	    	}
	    	else{
	    		if (charge>0){
	    			atomSmiles.append('+');
	    		}
	    		atomSmiles.append(charge);
	    	}
	    }
	    //atomSmiles.append("[id:"+atom.getID()+"]");
	    if (needsSquareBrackets){
	    	atomSmiles.append(']');
	    }
		return atomSmiles.toString();
	}

	private int calculateNumberOfBondedExplicitHydrogen(Atom atom) {
		List<Atom> neighbours = atom.getAtomNeighbours();
		int count =0;
		for (Atom neighbour : neighbours) {
			if (neighbour.getProperty(Atom.VISITED)==null){
				count++;
			}
		}
		return count;
	}

	private boolean determineWhetherAtomNeedsSquareBrackets(Atom atom, int hydrogenCount) {
		if (!organicAtomsToStandardValencies.containsKey(atom.getElement())){
			return true;
		}
		if (atom.getCharge()!=0){
			return true;
		}
		if (atom.getIsotope()!=null){
			return true;
		}
		if (atom.getAtomParity()!=null){
			return true;
		}
	
		List<Integer> expectedValencies = Arrays.asList(organicAtomsToStandardValencies.get(atom.getElement()));
		int valency = atom.getIncomingValency();
		boolean valencyCanBeDescribedImplicitly = expectedValencies.contains(valency);
		int targetImplicitValency =valency;
		if (valency > expectedValencies.get(expectedValencies.size()-1)){
			valencyCanBeDescribedImplicitly =true;
		}
		if (!valencyCanBeDescribedImplicitly){
			return true;
		}
	
		int nonHydrogenValency = valency -hydrogenCount;
		int implicitValencyThatWouldBeGenerated = nonHydrogenValency;
		for (int i = expectedValencies.size()-1; i>=0; i--) {
			if (expectedValencies.get(i) >= nonHydrogenValency){
				implicitValencyThatWouldBeGenerated =expectedValencies.get(i);
			}
		}
		if (targetImplicitValency != implicitValencyThatWouldBeGenerated){
			return true;
		}
		return false;
	}

	private String atomParityToSmiles(Atom currentAtom, int depth, Atom previousAtom) {
		AtomParity atomParity = currentAtom.getAtomParity();
		StringBuilder tetrahedralStereoChem = new StringBuilder();
		Atom[] atomRefs4 = atomParity.getAtomRefs4().clone();

		List<Atom> atomrefs4Current = new ArrayList<Atom>();

		Set<Bond> bonds = currentAtom.getBonds();
		for (Bond bond : bonds) {//previous atom
			Atom neighbour = bond.getOtherAtom(currentAtom);
			if (neighbour.getProperty(Atom.VISITED)!=null && neighbour.equals(previousAtom) ){
				atomrefs4Current.add(neighbour);
			}
		}
		for (Atom atom : atomRefs4) {//lone pair as in tetrahedral sulfones
			if (atom.equals(currentAtom)){
				atomrefs4Current.add(currentAtom);
			}
		}
		for (Bond bond : bonds) {//implicit hydrogen
			Atom neighbour = bond.getOtherAtom(currentAtom);
			if (neighbour.getProperty(Atom.VISITED)==null){
				atomrefs4Current.add(currentAtom);
			}
		}
		for (Bond bond : bonds) {//ring closures
			Atom neighbour = bond.getOtherAtom(currentAtom);
			if (neighbour.getProperty(Atom.VISITED)==null){
				continue;
			}
			if (neighbour.getProperty(Atom.VISITED)<=depth && !neighbour.equals(previousAtom) ){
				atomrefs4Current.add(neighbour);
			}
		}
		for (Bond bond : bonds) {//ring openings
			Atom neighbour = bond.getOtherAtom(currentAtom);
			if (neighbour.getProperty(Atom.VISITED)==null){
				continue;
			}
			if (neighbour.getProperty(Atom.VISITED)> (depth +1)){
				atomrefs4Current.add(neighbour);
			}

		}
		for (Bond bond : bonds) {//next atom/s
			Atom neighbour = bond.getOtherAtom(currentAtom);
			if (neighbour.getProperty(Atom.VISITED)==null){
				continue;
			}
			if (neighbour.getProperty(Atom.VISITED)==depth+1){
				atomrefs4Current.add(neighbour);
			}
		}
		Atom[] atomrefs4CurrentArr = new Atom[4];
		for (int i = 0; i < atomrefs4Current.size(); i++) {
			atomrefs4CurrentArr[i] = atomrefs4Current.get(i);
		}
		for (int i = 0; i < atomRefs4.length; i++) {//replace mentions of explicit hydrogen with the central atom the hydrogens are attached to, to be consistent with the SMILES representation
			if (atomRefs4[i].getProperty(Atom.VISITED)==null){
				atomRefs4[i] = currentAtom;
			}
		}

		boolean equivalent = StereochemistryHandler.checkEquivalencyOfAtomsRefs4AndParity(atomRefs4, atomParity.getParity(), atomrefs4CurrentArr, 1);
		if (equivalent){
			tetrahedralStereoChem.append("@@");
		}
		else{
			tetrahedralStereoChem.append("@");
		}
		return tetrahedralStereoChem.toString();
	}

	/**
	 * Generates the SMILES description of the bond
	 * In the case of cis/trans stereochemistry this relies on the assignDoubleBondStereochemistrySlashes
	 * having been run to setup the smilesBondDirection attribute
	 * @param bond
	 * @return
	 */
	private String bondToSmiles(Bond bond){
		String bondSmiles ="";
		int bondOrder = bond.getOrder();
		if (bondOrder==2){
			bondSmiles ="=";
		}
		else if (bondOrder==3){
			bondSmiles ="#";
		}
		else if (bond.getSmilesStereochemistry()!=null){
			if (bond.getSmilesStereochemistry()==SMILES_BOND_DIRECTION.RSLASH){
				bondSmiles ="/";
			}
			else{
				bondSmiles ="\\";
			}
		}
		return bondSmiles;
	}
}
