package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

class SMILESWriter {
	
	/**The organic atoms and their allowed implicit valencies in SMILES */
	private static final Map<String,Integer[]> organicAtomsToStandardValencies = new HashMap<String, Integer[]>();
	/**Closures 1-9, %10-99, 0 */
	private static final  LinkedList<String> closureSymbols = new LinkedList<String>();
	
	/**The available ring closure symbols, ordered from start to end in the preferred order for use.*/
	private final LinkedList<String> availableClosureSymbols = new LinkedList<String>(closureSymbols);
	/**Maps between bonds and the ring closure to use when the atom that ends the bond is encountered.*/
	private final HashMap<Bond, String> bondToClosureSymbolMap = new HashMap<Bond, String>();
	private final Fragment structure;

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
	 * Generates SMILES from the fragment the SMILESWriter was created with (currently achiral)
	 * The following assumptions are currently made:
	 * 	The fragment contains no bonds to atoms outside the fragment
	 * 	Hydrogens are all explicit
	 * 	Spare valency has been converted to double bonds.
	 * 
	 * The following is guaranteed:
	 * 	The order of output in the SMILES will match the order of atoms in the fragment
	 * 	Hydrogens are all treated as being implicit or as properties of a non hydrogen atom
	 * @return
	 */
	String generateSmiles() {
		StringBuilder sb = new StringBuilder();
		List<Atom> atomList = structure.getAtomList();
		List<Atom> nonProtonAtomList = createNonProtonAtomList(atomList);
		int nonProtonCount = nonProtonAtomList.size();
		for (int i = 0; i < nonProtonCount; i++) {
			Atom currentAtom = nonProtonAtomList.get(i);
			sb.append(atomToSmiles(currentAtom));
			Atom nextAtom = (i+1) < nonProtonCount ? nonProtonAtomList.get(i+1) : null;
			Atom previousAtom = (i-1) >=0 ? nonProtonAtomList.get(i-1) : null;
			sb.append(bondsFromCurrentAtomToOtherAtoms(currentAtom, nextAtom, previousAtom, nonProtonAtomList));
		}
		return sb.toString();
	}

	private List<Atom> createNonProtonAtomList(List<Atom> atomList) {
		List<Atom> nonProtonAtomList = new ArrayList<Atom>();
		for (Atom atom : atomList) {
			if (!atom.getElement().equals("H") || (atom.getIsotope()!=null && atom.getIsotope()!=1) ){
				nonProtonAtomList.add(atom);
			}
			else{
				//special case where hydrogen is a counter ion or only connects to other hydrogen
				List<Atom> neighbours = atom.getAtomNeighbours();
				boolean foundNonHydrogenNeighbour =false;
				for (Atom neighbour : neighbours) {
					if (!neighbour.getElement().equals("H")){
						foundNonHydrogenNeighbour =true;
					}
				}
				if (!foundNonHydrogenNeighbour || neighbours.size()>1){
					nonProtonAtomList.add(atom);
				}
			}
		}
		return nonProtonAtomList;
	}

	private String atomToSmiles(Atom atom) {
		StringBuilder atomSmiles = new StringBuilder();
		int hydrogen =calculateNumberOfBondedExplicitHydrogen(atom);
		boolean needsSquareBrackets = determineWhetherAtomNeedsSquareBrackets(atom, hydrogen);
		if (needsSquareBrackets){
			atomSmiles.append('[');
		}
		if (atom.getIsotope()!=null){
			atomSmiles.append(atom.getIsotope());
		}
		if (atom.hasSpareValency()){//spare valency corresponds directly to lower case SMILES in OPSIN's SMILES reader
			atomSmiles.append(atom.getElement().toLowerCase());
		}
		else{
			atomSmiles.append(atom.getElement());
		}
		if (hydrogen !=0 && needsSquareBrackets && !atom.getElement().equals("H")){
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
	    if (needsSquareBrackets){
	    	atomSmiles.append(']');
	    }
		return atomSmiles.toString();
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

	private int calculateNumberOfBondedExplicitHydrogen(Atom atom) {
		List<Atom> neighbours = atom.getAtomNeighbours();
		int count =0;
		for (Atom neighbour : neighbours) {
			if (neighbour.getElement().equals("H")){
				count++;
			}
		}
		return count;
	}
	
	/**
	 * Returns a string describing all the bonds from the current atom to its neighbours.
	 * This will take the form or a series of ring opening/closures followed by either nothing or a dot or the bond order to the next adjacent atom
	 * @param currentAtom
	 * @param nextAtom
	 * @param previousAtom
	 * @param nonProtonAtomList 
	 * @return
	 */
	private String bondsFromCurrentAtomToOtherAtoms(Atom currentAtom, Atom nextAtom, Atom previousAtom, List<Atom> nonProtonAtomList) {
		Set<Bond> bonds = currentAtom.getBonds();
		StringBuilder ringClosures = new StringBuilder();
		String bondToNextAtom = nextAtom==null ? "" : ".";
		for (Bond bond : bonds) {
			Atom toAtom =null;
			if (bond.getFromAtom() ==currentAtom){
				toAtom = bond.getToAtom();
			}
			else{
				toAtom = bond.getFromAtom();
			}
			if (nonProtonAtomList.contains(toAtom)){
				String bondSmiles = bondToSmiles(bond);
				if (nextAtom!=null && toAtom.equals(nextAtom)){
					bondToNextAtom = bondSmiles;
					continue;
				}
				if (bondToClosureSymbolMap.containsKey(bond)){
					String closure = bondToClosureSymbolMap.get(bond);
					availableClosureSymbols.addFirst(closure);
					ringClosures.append(closure);
				}
				else{
					if (toAtom.equals(previousAtom)){
						continue;//already bonded
					}
					String closure = availableClosureSymbols.removeFirst();
					bondToClosureSymbolMap.put(bond, closure);
					ringClosures.append(bondSmiles);
					ringClosures.append(closure);
				}
			}
		}
		return ringClosures.toString() + bondToNextAtom;
	}
	
	private String bondToSmiles(Bond bond){
		String bondSmiles ="";
		if (bond.getOrder()==2){
			bondSmiles ="=";
		}
		else if (bond.getOrder()==3){
			bondSmiles ="#";
		}
//TODO implement double bond stereochemistry
//		else if (bond.getSmilesStereochemistry()!=null){
//			if (bond.getSmilesStereochemistry()==SMILES_BOND_DIRECTION.RSLASH){
//				bondSmiles ="/";
//			}
//			else{
//				bondSmiles ="\\";
//			}
//		}
		return bondSmiles;
		
	}
}
