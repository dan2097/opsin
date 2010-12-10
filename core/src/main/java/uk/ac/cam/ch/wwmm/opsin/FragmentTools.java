package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

/**
 * Sorts a list of atoms such that their order agrees with the order symbolic locants are typically assigned
 * 
 * Preferred atoms are sorted to the START of the list
 * @author dl387
 *
 */
class SortAtomsForElementSymbols implements Comparator<Atom> {

	final Map<Atom, Bond> atomToPreviousBondMap;
    public SortAtomsForElementSymbols(Map<Atom, Bond> atomToPreviousBondMap) {
    	this.atomToPreviousBondMap = atomToPreviousBondMap;
	}

	public int compare(Atom a, Atom b){
    	Bond bondA = atomToPreviousBondMap.get(a);
    	Bond bondB = atomToPreviousBondMap.get(b);
    	if (bondA.getOrder() > bondB.getOrder()){//lower order bond is preferred
    		return 1;
    	}
    	if (bondA.getOrder() < bondB.getOrder()){
    		return -1;
    	}
    	
    	if (a.getOutValency() > b.getOutValency()){//prefer atoms with outValency
    		return -1;
    	}
    	if (a.getOutValency() < b.getOutValency()){
    		return 1;
    	}

    	int expectedHydrogenA = StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(a);
    	int expectedHydrogenB = StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(b);
    	
    	
    	if (expectedHydrogenA > expectedHydrogenB){//prefer atoms with more hydrogen
    		return -1;
    	}
    	if (expectedHydrogenA < expectedHydrogenA){
    		return 1;
    	}
    	return 0;
    }
}

/**
 * Performs a very crude sort of atoms such that those that are more likely to be substitued are preferred for low locants
 * Preferred atoms are sorted to the START of the list
 * @author dl387
 *
 */
class SortAtomsForMainGroupElementSymbols implements Comparator<Atom> {

    public int compare(Atom a, Atom b){
    	int compare =a.getElement().compareTo(b.getElement());
    	if (compare !=0){//only bother comparing properly if elements are the same
    		return compare;
    	}

    	int aExpectedHydrogen = StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(a);
    	int bExpectedHydrogen =StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(b);
    	if (aExpectedHydrogen >0 &&  bExpectedHydrogen ==0){//having substitutable hydrogen preferred
    		return -1;
    	}
    	if (aExpectedHydrogen ==0 && bExpectedHydrogen >0){
    		return 1;
    	}
    	List<String> locantsA = a.getLocants();
    	List<String> locantsB = b.getLocants();
    	if (locantsA.size() ==0 &&  locantsB.size() >0){//having no locants preferred
    		return -1;
    	}
    	if (locantsA.size() >0 &&  locantsB.size() ==0){
    		return 1;
    	}
    	return 0;
    }
}

class FragmentTools {
	private static final Pattern matchElementSymbolLocant =Pattern.compile("[A-Z][a-z]?'*");
	private static final Pattern matchNumericLocant =Pattern.compile("\\d+[a-z]?'*");
	
	/**
	 * Sorts by number, then by letter e.g. 4,3,3b,5,3a,2 -->2,3,3a,3b,4,5
	 * @author dl387
	 *
	 */
	static class SortByLocants implements Comparator<Atom> {
      	static final Pattern locantSegmenter =Pattern.compile("(\\d+)([a-z]?)('*)");

	    public int compare(Atom atoma, Atom atomb){
	    	if (atoma.getType().equals(SUFFIX_TYPE_VAL) && !atomb.getType().equals(SUFFIX_TYPE_VAL)){//suffix atoms go to the back
	    		return 1;
	    	}
	    	if (atomb.getType().equals(SUFFIX_TYPE_VAL) && !atoma.getType().equals(SUFFIX_TYPE_VAL)){
	    		return -1;
	    	}

	    	String locanta =atoma.getFirstLocant();
	    	String locantb =atomb.getFirstLocant();
	    	if (locanta==null|| locantb==null){
	    		return 0;
	    	}

	    	Matcher m1  =locantSegmenter.matcher(locanta);
	    	Matcher m2  =locantSegmenter.matcher(locantb);
	    	if (!m1.matches()|| !m2.matches()){//inappropriate locant
	    		return 0;
	    	}
        	String locantaPrimes = m1.group(3);
        	String locantbPrimes = m2.group(3);
            if (locantaPrimes.compareTo(locantbPrimes)>=1) {
                return 1;//e.g. 1'' vs 1'
            } else if (locantbPrimes.compareTo(locantaPrimes)>=1) {
                return -1;//e.g. 1' vs 1''
            }
            else{
		    	int locantaNumber = Integer.parseInt(m1.group(1));
		    	int locantbNumber = Integer.parseInt(m2.group(1));
	
		        if (locantaNumber >locantbNumber) {
		            return 1;//e.g. 3 vs 2 or 3a vs 2
		        } else if (locantbNumber >locantaNumber) {
		            return -1;//e.g. 2 vs 3 or 2 vs 3a
		        }
		        else{
		        	String locantaLetter = m1.group(2);
		        	String locantbLetter = m2.group(2);
		            if (locantaLetter.compareTo(locantbLetter)>=1) {
		                return 1;//e.g. 1b vs 1a
		            } else if (locantbLetter.compareTo(locantaLetter)>=1) {
		                return -1;//e.g. 1a vs 1b
		            }
			        return 0;
		        }
            }
	    }
	}
	
	/**
	 * Assign element locants to groups/suffixes. These are in addition to any numerical locants that are present.
	 * Adds primes to make each locant unique.
	 * For groups a locant is not given to carbon atoms
	 * If an element appears in a suffix then element locants are not assigned to occurrences of that element in the parent group
	 * HeteroAtoms in acidStems connected to the first Atom of the fragment are treated as if they were suffix atoms
	 * @param suffixableFragment
	 * @param suffixFragments
	 * @throws StructureBuildingException 
	 */
	static void assignElementLocants(Fragment suffixableFragment, List<Fragment> suffixFragments) throws StructureBuildingException {
		
		HashMap<String,Integer> elementCount =new HashMap<String,Integer>();//keeps track of how many times each element has been seen
		HashSet<Atom> atomsToIgnore = new HashSet<Atom>();//atoms which already have a symbolic locant
		
		ArrayList<Fragment> allFragments =new ArrayList<Fragment>(suffixFragments);
		allFragments.add(suffixableFragment);
		/*
		 * First check whether any element locants have already been assigned, these will take precedence
		 */
		for (Fragment fragment : allFragments) {
			List<Atom> atomList =fragment.getAtomList();
			for (Atom atom : atomList) {
				List<String> elementSymbolLocants =atom.getElementSymbolLocants();
				for (String locant : elementSymbolLocants) {
					int primeCount =0;
					for(int i=0;i<locant.length();i++) {
						if(locant.charAt(i) == '\'') primeCount++;
					}
					String element =locant.substring(0, locant.length()-primeCount);
					if (elementCount.get(element)==null || (elementCount.get(element) < primeCount +1)){
						elementCount.put(element, primeCount +1);
					}
					atomsToIgnore.add(atom);
				}
			}
		}
		
		HashSet<String> elementToIgnore = new HashSet<String>();
		for (String element : elementCount.keySet()) {
			elementToIgnore.add(element);
		}
		for (Fragment fragment : allFragments) {
			List<Atom> atomList =fragment.getAtomList();
			for (Atom atom : atomList) {
				if (elementToIgnore.contains(atom.getElement())){
					atomsToIgnore.add(atom);
				}
			}
		}
		
		
		String fragType = suffixableFragment.getType();
		if (fragType.equals(NONCARBOXYLICACID_TYPE_VAL) || fragType.equals(CHALCOGENACIDSTEM_TYPE_VAL)){
			if (suffixFragments.size()!=0){
				throw new StructureBuildingException("No suffix fragments were expected to be present on non carboxylic acid");
			}
			processNonCarboxylicAcidLabelling(suffixableFragment, elementCount, atomsToIgnore);
		}
		else{
			if (suffixFragments.size()>0){
				processSuffixLabelling(suffixFragments, elementCount, atomsToIgnore);
				if (elementCount.get("N")!=null &&  elementCount.get("N")>1){//look for special case violation of IUPAC rule, =(N)=(NN) is N//N' in practice rather than N/N'/N''
					//this method will put both locants on the N with substituable hydrogen
					detectAndCorrectHydrazoneDerivativeViolation(suffixFragments);
				}
			}
			processMainGroupLabelling(suffixableFragment, elementCount, atomsToIgnore);
		}
	}


	private static void detectAndCorrectHydrazoneDerivativeViolation(List<Fragment> suffixFragments) {
		fragmentLoop: for (Fragment suffixFrag : suffixFragments) {
			List<Atom> atomList = suffixFrag.getAtomList();
			for (Atom atom : atomList) {
				if (atom.getElement().equals("N") && atom.getIncomingValency() ==3 ){
					List<String> locants =atom.getLocants();
					if (locants.size()==1 && matchElementSymbolLocant.matcher(locants.get(0)).matches()){
						List<Atom> neighbours = atom.getAtomNeighbours();
						for (Atom neighbour : neighbours) {
							if (neighbour.getElement().equals("N") && neighbour.getIncomingValency()==1){
								String locantToAdd = locants.get(0);
								atom.clearLocants();
								neighbour.addLocant(locantToAdd);
								continue fragmentLoop;
							}
						}
					}
				}
			}
		}
	}


	private static void processMainGroupLabelling(Fragment suffixableFragment, HashMap<String, Integer> elementCount, HashSet<Atom> atomsToIgnore) {
		HashSet<String> elementToIgnore = new HashSet<String>();
		for (String element : elementCount.keySet()) {
			elementToIgnore.add(element);
		}
		List<Atom> atomList =suffixableFragment.getAtomList();
		Collections.sort(atomList, new SortAtomsForMainGroupElementSymbols());
		Atom atomToAddCLabelTo = null;//only add a C label if there is only one C in the main group
		boolean seenMoreThanOneC = false;
		for (Atom atom : atomList) {
			if (atomsToIgnore.contains(atom)){
				continue;
			}
			String element =atom.getElement();
			if (elementToIgnore.contains(element)){
				continue;
			}
			if (element.equals("C")){
				if (seenMoreThanOneC){
					continue;
				}
				if (atomToAddCLabelTo !=null){
					atomToAddCLabelTo = null;
					seenMoreThanOneC =true;
				}
				else{
					atomToAddCLabelTo =atom;
				}
				continue;
			}
			if (elementCount.get(element)==null){
				atom.addLocant(element);
				elementCount.put(element,1);
			}
			else{
				int count =elementCount.get(element);
				atom.addLocant(element + StringTools.multiplyString("'", count));
				elementCount.put(element, count +1);
			}
		}
		if (atomToAddCLabelTo !=null){
			atomToAddCLabelTo.addLocant("C");
		}
	}


	private static void processSuffixLabelling(List<Fragment> suffixFragments, HashMap<String, Integer> elementCount, HashSet<Atom> atomsToIgnore) throws StructureBuildingException {
		LinkedList<Atom> startingAtoms = new LinkedList<Atom>();
		Map<Atom, Bond> atomPreviousBondMap = new HashMap<Atom, Bond>();
		Set<Atom> atomsVisited = new HashSet<Atom>();
		for (Fragment fragment : suffixFragments) {
			List<Atom> suffixAtomList =fragment.getAtomList();
			Atom rAtom = suffixAtomList.get(0);
			LinkedList<Atom> nextAtoms = new LinkedList<Atom>(rAtom.getAtomNeighbours());
			for (Atom nextAtom : nextAtoms) {
				atomsVisited.add(nextAtom);
				atomPreviousBondMap.put(nextAtom, fragment.findBondOrThrow(rAtom, nextAtom));
			}
			startingAtoms.addAll(nextAtoms);
		}
		Collections.sort(startingAtoms, new SortAtomsForElementSymbols(atomPreviousBondMap));

		while (startingAtoms.size() > 0){
			assignLocantsAndExploreNeighbours(elementCount, atomsToIgnore, atomsVisited, startingAtoms);
		}
	}


	private static void processNonCarboxylicAcidLabelling(Fragment suffixableFragment, HashMap<String, Integer> elementCount,HashSet<Atom> atomsToIgnore) throws StructureBuildingException {
		Set<Atom> atomsVisited = new HashSet<Atom>();
		List<Atom> atomList =suffixableFragment.getAtomList();
		Atom firstAtom = atomList.get(0);
		LinkedList<Atom> nextAtoms = new LinkedList<Atom>(firstAtom.getAtomNeighbours());
		Map<Atom, Bond> atomPreviousBondMap = new HashMap<Atom, Bond>();
		for (Atom nextAtom : nextAtoms) {
			atomPreviousBondMap.put(nextAtom, suffixableFragment.findBondOrThrow(firstAtom, nextAtom));
		}
		Collections.sort(nextAtoms, new SortAtomsForElementSymbols(atomPreviousBondMap));
		atomsVisited.add(firstAtom);
		while (nextAtoms.size() > 0){
			assignLocantsAndExploreNeighbours(elementCount, atomsToIgnore, atomsVisited, nextAtoms);
		}
		if (!atomsToIgnore.contains(firstAtom) && firstAtom.determineValency(true) > firstAtom.getIncomingValency()){
			//e.g. carbonimidoyl the carbon has locant C
			assignLocant(firstAtom, elementCount);
		}
	}


	private static void assignLocantsAndExploreNeighbours(HashMap<String, Integer> elementCount, HashSet<Atom> atomsToIgnore, Set<Atom> atomsVisited, LinkedList<Atom> nextAtoms) throws StructureBuildingException {
		Atom atom = nextAtoms.removeFirst();
		atomsVisited.add(atom);
		if (!atomsToIgnore.contains(atom)){//assign locant
			assignLocant(atom, elementCount);
		}
		List<Atom> atomNeighbours = atom.getAtomNeighbours();
		for (int i = atomNeighbours.size() -1; i >=0; i--) {
			if (atomsVisited.contains(atomNeighbours.get(i))){
				atomNeighbours.remove(i);
			}
		}
		Map<Atom, Bond> atomPreviousBondMap = new HashMap<Atom, Bond>();
		for (Atom atomNeighbour : atomNeighbours) {
			atomPreviousBondMap.put(atomNeighbour, atom.getFrag().findBondOrThrow(atom, atomNeighbour));
		}
		Collections.sort(atomNeighbours, new SortAtomsForElementSymbols(atomPreviousBondMap));
		nextAtoms.addAll(0, atomNeighbours);
	}


	private static void assignLocant(Atom atom, HashMap<String, Integer> elementCount) {
		String element =atom.getElement();
		if (elementCount.get(element)==null){
			atom.addLocant(element);
			elementCount.put(element,1);
		}
		else{
			int count =elementCount.get(element);
			atom.addLocant(element + StringTools.multiplyString("'", count));
			elementCount.put(element, count +1);
		}
	}


	/** Adjusts the order of a bond in a fragment.
	 *
	 * @param fromAtom The lower-numbered atom in the bond
	 * @param bondOrder The new bond order
	 * @param fragment The fragment
     * @throws StructureBuildingException
	 */
	static void unsaturate(Atom fromAtom, int bondOrder, Fragment fragment) throws StructureBuildingException {
		List<Atom> atomList = fragment.getAtomList();
		Atom toAtom = null;
		int initialIndice = atomList.indexOf(fromAtom);
		if (initialIndice +1 < atomList.size()){
			toAtom = atomList.get(initialIndice +1);
		}
		if (toAtom==null || toAtom.getType().equals(SUFFIX_TYPE_VAL)){//allows something like cyclohexan-6-ene, something like butan-4-ene will still fail
			List<Atom> neighbours = fromAtom.getAtomNeighbours();
			if (neighbours.size() >=2){
				Atom firstAtomInFrag = atomList.get(0);
				for (Atom a : neighbours) {
					if (a ==firstAtomInFrag){
						toAtom = a;
						break;
					}
				}
			}
		}
		Bond b = fragment.findBondOrThrow(fromAtom, toAtom);
		b.setOrder(bondOrder);
	}

	/** Adjusts the order of a bond in a fragment.
	 *
	 * @param fromAtom The first atom in the bond
	 * @param locantTo The locant of the other atom in the bond
	 * @param bondOrder The new bond order
	 * @param fragment The fragment
     * @throws StructureBuildingException
	 */
	static void unsaturate(Atom fromAtom, String locantTo, int bondOrder, Fragment fragment) throws StructureBuildingException {
		Atom toAtom = fragment.getAtomByLocantOrThrow(locantTo);
		Bond b = fragment.findBondOrThrow(fromAtom, toAtom);
		b.setOrder(bondOrder);
	}
	
	/** Works out where to put an "one", if this is unspecified. position 2 for propanone
	 * and higher, else 1. Position 2 is assumed to be 1 higher than the atomIndice given.
	 *
	 * @param fragment The fragment
	 * @param atomIndice
     * @return the appropriate atom indice
	 * @throws StructureBuildingException 
	 */
	static int findKetoneAtomIndice(Fragment fragment, int atomIndice) throws StructureBuildingException {
		if(fragment.getChainLength() < 3){
			return atomIndice;
		}
		else {
			if (atomIndice +1>=fragment.getAtomList().size()){
				return 1;//this probably indicates a problem with the input name but nonetheless 1 is a better answer than an indice which isn't even in the range of the fragment
			}
			else{
				return atomIndice +1;
			}
		}
	}
	
	/**Adjusts the labeling on a fused ring system, such that bridgehead atoms
	 * have locants endings in 'a' or 'b' etc. Example: naphthalene
	 * 1,2,3,4,5,6,7,8,9,10->1,2,3,4,4a,5,6,7,8,8a
     * @param fusedring
	 */
	static void relabelFusedRingSystem(Fragment fusedring){
		relabelFusedRingSystem(fusedring.getAtomList());
	}

	/**Adjusts the labeling on a fused ring system, such that bridgehead atoms
	 * have locants endings in 'a' or 'b' etc. Example: naphthalene
	 * 1,2,3,4,5,6,7,8,9,10->1,2,3,4,4a,5,6,7,8,8a
     * @param atomList
	 */
	static void relabelFusedRingSystem(List<Atom> atomList) {
		int locantVal = 0;
		char locantLetter = 'a';
		for (Atom atom : atomList) {
			atom.clearLocants();
		}
		for (Atom atom : atomList) {
			if(!atom.getElement().equals("C") || atom.getBonds().size() < 3) {
				locantVal++;
				locantLetter = 'a';
				atom.addLocant(Integer.toString(locantVal));
			} else {
				atom.addLocant(Integer.toString(locantVal) + locantLetter);
				locantLetter++;
			}
		}
	}
	
	/**
	 * Adds the given string to all the locants of the atoms.
	 * @param atomList
	 * @param stringToAdd
	 */
	static void relabelLocants(List<Atom> atomList, String stringToAdd) {
		for (Atom atom : atomList) {
			List<String> locants = new ArrayList<String>(atom.getLocants());
			atom.clearLocants();
			for (String locant : locants) {
				atom.addLocant(locant + stringToAdd);
			}
		}
	}
	
	/**
	 * Adds the given string to all the numeric locants of the atoms.
	 * @param atomList
	 * @param stringToAdd
	 */
	static void relabelNumericLocants(List<Atom> atomList, String stringToAdd) {
		for (Atom atom : atomList) {
			List<String> locants = new ArrayList<String>(atom.getLocants());
			for (String locant : locants) {
				if (matchNumericLocant.matcher(locant).matches()){
					atom.removeLocant(locant);
					atom.addLocant(locant + stringToAdd);
				}
			}
		}
	}


	static void splitOutAtomIntoValency1OutAtoms(OutAtom outAtom) {
		Fragment frag =outAtom.getAtom().getFrag();
		for (int i = 1; i < outAtom.getValency(); i++) {
			frag.addOutAtom(outAtom.getAtom(), 1, outAtom.isSetExplicitly());
		}
		outAtom.setValency(1);
	}

	/**
	 * Checks if the specified Nitrogen is potentially involved in [NH]C=N <-> N=C[NH] tautomerism
	 * Given the starting nitrogen returns the other nitrogen or null if that nitrogen does not appear to be involved in such tautomerism
	 * @param nitrogen
	 * @return null or the other nitrogen
	 */
	static Atom detectSimpleNitrogenTautomer(Atom nitrogen) {
		if (nitrogen.getElement().equals("N") && nitrogen.getAtomIsInACycle()){
			List<Atom> neighbours = nitrogen.getAtomNeighbours();
			for (Atom neighbour : neighbours) {
				if (neighbour.hasSpareValency()){
					List<Atom> distance2Neighbours = neighbour.getAtomNeighbours();
					for (Atom distance2Neighbour : distance2Neighbours) {
						if (distance2Neighbour.hasSpareValency() && distance2Neighbour.getElement().equals("N") && distance2Neighbour.getCharge()==0){
							return distance2Neighbour;
						}
					}
				}
			}
		}
		return null;
	}
	

	/** Looks for the atom that would have had a hydrogen indicated,
	 * adds a spareValency to that atom, and sets indicatedHydrogen.
     * @param frag
	 */
	static void pickUpIndicatedHydrogens(Fragment frag) {
		if (frag.getIndicatedHydrogen().size()>0){
			return;
		}
		List<Atom> atomList = frag.getAtomList();
		int svCount = 0;
		int dvCount = 0; /* Divalent, like O */
		for(Atom a : atomList) {
			svCount += a.hasSpareValency() ? 1 : 0;
			if(a.getElement().equals("O")) dvCount++;
			if(a.getElement().equals("S")) dvCount++;
			if(a.getElement().equals("Se")) dvCount++;
			if(a.getElement().equals("Te")) dvCount++;
		}
		if(svCount == atomList.size() - (1 + dvCount) && svCount > 0) {
			// Now we're looking for a trivalent C, or failing that a divalent N/P/As/Sb
			Atom nCandidate = null;
			for(Atom a : atomList) {
				if(a.getAtomIsInACycle() && !a.hasSpareValency()) {
					String element =a.getElement();
					if (element.equals("C")){
						frag.addIndicatedHydrogen(a);
						a.setSpareValency(true);
						return;
					} else if(element.equals("N")) {
						nCandidate = a;
					}
				}
			}
			if(nCandidate != null) {
				frag.addIndicatedHydrogen(nCandidate);
				nCandidate.setSpareValency(true);
			}
		}
	}

	/**Looks for double and higher bonds, converts them to single bonds
	 * and adds corresponding spareValencies to the atoms they join.
	 * @param frag
	 */
	static void convertHighOrderBondsToSpareValencies(Fragment frag)  {
		Set<Bond> bondSet = frag.getBondSet();
		for(Bond b : bondSet) {
			if(b.getOrder() == 2) {
				Atom firstAtom =b.getFromAtom();
				Atom secondAtom =b.getToAtom();
				if (firstAtom.getAtomIsInACycle() && secondAtom.getAtomIsInACycle()){
					b.setOrder(1);
					firstAtom.setSpareValency(true);
					secondAtom.setSpareValency(true);
				}
			}
		}
		FragmentTools.pickUpIndicatedHydrogens(frag);
	}

	/**Increases the order of bonds joining atoms with spareValencies,
	 * and uses up said spareValencies.
	 * @param frag
     * @throws StructureBuildingException If the algorithm can't work out where to put the bonds
	 */
	static void convertSpareValenciesToDoubleBonds(Fragment frag) throws StructureBuildingException {
		List<Atom> atomCollection =frag.getAtomList();
		/* pick atom, getAtomNeighbours, decideIfTerminal, resolve */

		/*
		 * Correct spare valency by looking at valencyState of atom
		 *
		 */
		for(Atom a : atomCollection) {
			a.ensureSVIsConsistantWithValency(true);
		}

		/*
		 * Remove spare valency on atoms which may not form higher order bonds
		 */
		atomLoop: for(Atom a : atomCollection) {
			if(a.hasSpareValency()) {
				for(Atom aa : frag.getIntraFragmentAtomNeighbours(a)) {
					if(aa.hasSpareValency()){
						continue atomLoop;
					}
				}
				a.setSpareValency(false);
			}
		}

		/*
		 Reduce valency of atoms which cannot possibly have any of their bonds converted to double bonds
		 pick an atom which definitely does have spare valency to be the indicated hydrogen.
		*/
		Atom atomToReduceValencyAt =null;
		List<Atom> indicatedHydrogen = frag.getIndicatedHydrogen();
		for (int i = indicatedHydrogen.size() -1; i >=0; i--) {
			if (!indicatedHydrogen.get(i).hasSpareValency()){
				indicatedHydrogen.remove(i);
			}
		}
		if (indicatedHydrogen.size()>0){
			if (indicatedHydrogen.size()>1){
				if (indicatedHydrogen.size()==2){//fix for guanine like purine derivatives
					for (Atom indicatedAtom : indicatedHydrogen) {
						if (!indicatedAtom.hasSpareValency()){
							continue;
						}
						List<Atom> neighbours = indicatedAtom.getAtomNeighbours();
						for (Atom neighbour : neighbours) {
							if (neighbour.hasSpareValency()){
								List<Atom> neighbours2 = neighbour.getAtomNeighbours();
								for (Atom secondNeighbour : neighbours2) {
									if (secondNeighbour==indicatedAtom){
										continue;
									}
									if (secondNeighbour.hasSpareValency() && secondNeighbour.getElement().equals("N")){
										indicatedAtom.setSpareValency(false);
									}
								}
							}
						}
					}
				}
				else{
					for (Atom indicatedAtom : indicatedHydrogen) {
						indicatedAtom.setSpareValency(false);
					}
				}
			}
			else{
				atomToReduceValencyAt = indicatedHydrogen.get(0);
			}
		}
		
		int svCount = 0;
		for(Atom a : atomCollection) {
			svCount += a.hasSpareValency() ? 1 :0;
		}

		if((svCount % 2) == 1) {
			if (atomToReduceValencyAt ==null){
				for(Atom a : atomCollection) {//try and find an atom with SV that neighbours only one atom with SV
					if(a.hasSpareValency()) {
						int atomsWithSV =0;
						for(Atom aa : frag.getIntraFragmentAtomNeighbours(a)) {
							if(aa.hasSpareValency()) {
								atomsWithSV++;
							}
						}
						if (atomsWithSV==1){
							atomToReduceValencyAt=a;
							break;
						}
					}
				}
				if (atomToReduceValencyAt==null){
					atomLoop: for(Atom a : atomCollection) {//try and find an atom with bridgehead atoms with SV on both sides c.f. phenoxastibinine ==10H-phenoxastibinine
						if(a.hasSpareValency()) {
							List<Atom> neighbours =frag.getIntraFragmentAtomNeighbours(a);
							if (neighbours.size()==2){
								for(Atom aa : neighbours) {
									if(frag.getIntraFragmentAtomNeighbours(aa).size() < 3){
										continue atomLoop;
									}
								}
								atomToReduceValencyAt=a;
								break;
							}
						}
					}
					if (atomToReduceValencyAt==null){//Prefer nitrogen to carbon e.g. get NHC=C rather than N=CCH
						for(Atom a : atomCollection) {
							if(a.hasSpareValency()) {
								if (atomToReduceValencyAt==null){
									atomToReduceValencyAt=a;//else just go with the first atom with SV encountered
								}
								if (!a.getElement().equals("C")){
									atomToReduceValencyAt=a;
									break;
								}
							}
						}
					}
				}
			}
			atomToReduceValencyAt.setSpareValency(false);
			svCount--;
		}

		while(svCount > 0) {
			boolean foundTerminalFlag = false;
			boolean foundNonBridgeHeadFlag = false;
			boolean foundBridgeHeadFlag = false;
			for(Atom a : atomCollection) {
				if(a.hasSpareValency()) {
					int count = 0;
					for(Atom aa : frag.getIntraFragmentAtomNeighbours(a)) {
						if(aa.hasSpareValency()) {
							count++;
						}
					}
					if(count == 1) {
						for(Atom aa : frag.getIntraFragmentAtomNeighbours(a)) {
							if(aa.hasSpareValency()) {
								foundTerminalFlag = true;
								a.setSpareValency(false);
								aa.setSpareValency(false);
								frag.findBondOrThrow(a, aa).addOrder(1);
								svCount -= 2;//Two atoms where for one of them this bond is the only double bond it can possible form
								break;
							}
						}
					}
				}
			}
			if(!foundTerminalFlag) {
				for(Atom a : atomCollection) {
					List<Atom> neighbours =frag.getIntraFragmentAtomNeighbours(a);
					if(a.hasSpareValency() && neighbours.size() < 3) {
						for(Atom aa : neighbours) {
							if(aa.hasSpareValency()) {
								foundNonBridgeHeadFlag = true;
								a.setSpareValency(false);
								aa.setSpareValency(false);
								frag.findBondOrThrow(a, aa).addOrder(1);
								svCount -= 2;//Two atoms where one of them is not a bridge head
								break;
							}
						}
					}
					if(foundNonBridgeHeadFlag) break;
				}
				if(!foundNonBridgeHeadFlag){
					for(Atom a : atomCollection) {
						List<Atom> neighbours =frag.getIntraFragmentAtomNeighbours(a);
						if(a.hasSpareValency()) {
							for(Atom aa : neighbours) {
								if(aa.hasSpareValency()) {
									foundBridgeHeadFlag = true;
									a.setSpareValency(false);
									aa.setSpareValency(false);
									frag.findBondOrThrow(a, aa).addOrder(1);
									svCount -= 2;//Two atoms where both of them are a bridge head e.g. necessary for something like coronene
									break;
								}
							}
						}
						if(foundBridgeHeadFlag) break;
					}
					if(!foundBridgeHeadFlag){
						throw new StructureBuildingException("Could not assign all higher order bonds.");
					}
				}
			}
		}
	}


	static Atom getAtomByAminoAcidStyleLocant(Atom backboneAtom, String elementSymbol, String primes) throws StructureBuildingException {
		//Search for appropriate atom by using the same algorithm as is used to assign locants initially
		Fragment frag = backboneAtom.getFrag();

		LinkedList<Atom> nextAtoms = new LinkedList<Atom>();
		Map<Atom, Bond> atomPreviousBondMap = new HashMap<Atom, Bond>();
		Set<Atom> atomsVisited = new HashSet<Atom>();
		List<Atom> neighbours = backboneAtom.getAtomNeighbours();
		mainLoop: for (Atom neighbour : neighbours) {
			atomsVisited.add(neighbour);
			if (!neighbour.getType().equals(SUFFIX_TYPE_VAL)){
				for (String neighbourLocant : neighbour.getLocants()) {
					if (matchNumericLocant.matcher(neighbourLocant).matches()){//gone to an inappropriate atom
						continue mainLoop;
					}
				}
			}
			nextAtoms.add(neighbour);
			atomPreviousBondMap.put(neighbour, frag.findBondOrThrow(backboneAtom, neighbour));
		}

		Collections.sort(nextAtoms, new SortAtomsForElementSymbols(atomPreviousBondMap));
		HashMap<String,Integer> elementCount =new HashMap<String,Integer>();//keeps track of how many times each element has been seen
	
		boolean hydrazoneSpecialCase =false;//look for special case violation of IUPAC rule where the locant of the =N- atom is skipped. This flag is set when =N- is encountered
		while (nextAtoms.size() > 0){
			Atom atom = nextAtoms.removeFirst();
			atomsVisited.add(atom);
			int primesOnPossibleAtom =0;
			String element =atom.getElement();
			if (elementCount.get(element)==null){
				elementCount.put(element,1);
			}
			else{
				int count =elementCount.get(element);
				primesOnPossibleAtom =count;
				elementCount.put(element, count +1);
			}
			if (hydrazoneSpecialCase){
				if (element.equals(elementSymbol) && primes.length() == primesOnPossibleAtom -1){
					return atom;
				}
				hydrazoneSpecialCase =false;
			}

			List<Atom> atomNeighbours = atom.getAtomNeighbours();
			for (int i = atomNeighbours.size() -1; i >=0; i--) {
				Atom neighbour = atomNeighbours.get(i);
				if (atomsVisited.contains(neighbour)){
					atomNeighbours.remove(i);
				}
				if (!neighbour.getType().equals(SUFFIX_TYPE_VAL)){
					for (String neighbourLocant : neighbour.getLocants()) {
						if (matchNumericLocant.matcher(neighbourLocant).matches()){//gone to an inappropriate atom
							atomNeighbours.remove(i);
							break;
						}
					}
				}
			}
			if (atom.getElement().equals("N") && atom.getIncomingValency() ==3 && atom.getCharge()==0 
					&& atomNeighbours.size()==1 && atomNeighbours.get(0).getElement().equals("N")){
				hydrazoneSpecialCase =true;
			}
			else{
				if (element.equals(elementSymbol)){
					if (primes.length() == primesOnPossibleAtom){
						return atom;
					}
				}
			}
			atomPreviousBondMap = new HashMap<Atom, Bond>();
			for (Atom atomNeighbour : atomNeighbours) {
				atomPreviousBondMap.put(atomNeighbour, frag.findBondOrThrow(atom, atomNeighbour));
			}
			Collections.sort(atomNeighbours, new SortAtomsForElementSymbols(atomPreviousBondMap));
			nextAtoms.addAll(0, atomNeighbours);
		}

		if (primes.equals("") && backboneAtom.getElement().equals(elementSymbol)){//maybe it meant the starting atom
			return backboneAtom;
		}
		return null;
	}
	
	
	/**
	 * Determines whether the bond between two elements is likely to be covalent
	 * This is crudely determined based on whether the combination of elements fall outside the ionic and
	 * metallic sections of a van Arkel diagram
	 * @param element1
	 * @param element2
	 * @return
	 */
	static boolean isCovalent(String element1, String element2) {
		Double atom1Electrongegativity = AtomProperties.elementToPaulingElectronegativity.get(element1);
		Double atom2Electrongegativity = AtomProperties.elementToPaulingElectronegativity.get(element2);
		if (atom1Electrongegativity!=null && atom2Electrongegativity !=null){
			double halfSum = (atom1Electrongegativity + atom2Electrongegativity)/2;
			double difference = Math.abs(atom1Electrongegativity - atom2Electrongegativity);
			if (halfSum < 1.6){
				return false;//probably metallic
			}
			if (difference < 1.39* halfSum -2.2){
				return true;
			}			
		}
		return false;
	}
	
	/**
	 * Is the atom a suffix atom or an aldehyde atom or a chalcogen functional atom
	 * @param atom
	 * @return
	 */
	static boolean isCharacteristicAtom(Atom atom) {
		if (atom.getType().equals(SUFFIX_TYPE_VAL)){
			return true;
		}
		if (atom.getProperty(Atom.ISALDEHYDE)!=null && atom.getProperty(Atom.ISALDEHYDE)){//substituting an aldehyde would make it no longer an aldehyde
			return true;
		}
		
		String element =atom.getElement();
		if (element.equals("O")|| element.equals("S") || element.equals("Se") || element.equals("Te")){//potential chalcogen functional atom
			boolean isFunctionalAtom =false;
			for (FunctionalAtom funcAtom : atom.getFrag().getFunctionalAtoms()) {
				if (atom.equals(funcAtom.getAtom())){
					isFunctionalAtom =true;
					break;
				}
			}
			if (isFunctionalAtom){
				return true;
			}
		}
		return false;
	}
}
