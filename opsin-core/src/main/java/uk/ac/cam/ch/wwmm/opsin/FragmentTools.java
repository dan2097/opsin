package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static uk.ac.cam.ch.wwmm.opsin.OpsinTools.*;
import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

/**
 * Sorts a list of atoms such that their order agrees with the order symbolic locants are typically assigned
 * 
 * Preferred atoms are sorted to the START of the list
 * @author dl387
 *
 */
class SortAtomsForElementSymbols implements Comparator<Atom> {

	public int compare(Atom a, Atom b){
		int bondOrderA = a.getProperty(Atom.VISITED);
		int bondOrderB = b.getProperty(Atom.VISITED);
    	if (bondOrderA > bondOrderB) {//lower order bond is preferred
    		return 1;
    	}
    	if (bondOrderA < bondOrderB) {
    		return -1;
    	}
    	
    	if (a.getOutValency() > b.getOutValency()) {//prefer atoms with outValency
    		return -1;
    	}
    	if (a.getOutValency() < b.getOutValency()) {
    		return 1;
    	}

    	int expectedHydrogenA = StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(a);
    	int expectedHydrogenB = StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(b);

    	if (expectedHydrogenA > expectedHydrogenB) {//prefer atoms with more hydrogen
    		return -1;
    	}
    	if (expectedHydrogenA < expectedHydrogenB) {
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
    	int compare = a.getElement().compareTo(b.getElement());
    	if (compare != 0) {//only bother comparing properly if elements are the same
    		return compare;
    	}

    	int aExpectedHydrogen = StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(a);
    	int bExpectedHydrogen = StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(b);
    	if (aExpectedHydrogen > 0 &&  bExpectedHydrogen == 0) {//having substitutable hydrogen preferred
    		return -1;
    	}
    	if (aExpectedHydrogen == 0 && bExpectedHydrogen > 0) {
    		return 1;
    	}
    	List<String> locantsA = a.getLocants();
    	List<String> locantsB = b.getLocants();
    	if (locantsA.size() == 0 &&  locantsB.size() > 0) {//having no locants preferred
    		return -1;
    	}
    	if (locantsA.size() > 0 &&  locantsB.size() == 0) {
    		return 1;
    	}
    	return 0;
    }
}

class FragmentTools {
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
		
		Map<String,Integer> elementCount = new HashMap<String,Integer>();//keeps track of how many times each element has been seen
		Set<Atom> atomsToIgnore = new HashSet<Atom>();//atoms which already have a symbolic locant
		
		List<Fragment> allFragments = new ArrayList<Fragment>(suffixFragments);
		allFragments.add(suffixableFragment);
		/*
		 * First check whether any element locants have already been assigned, these will take precedence
		 */
		for (Fragment fragment : allFragments) {
			List<Atom> atomList = fragment.getAtomList();
			for (Atom atom : atomList) {
				List<String> elementSymbolLocants = atom.getElementSymbolLocants();
				for (String locant : elementSymbolLocants) {
					int primeCount = StringTools.countTerminalPrimes(locant);
					String element = locant.substring(0, locant.length() - primeCount);
					Integer seenCount = elementCount.get(element);
					if (seenCount == null || (seenCount < primeCount + 1)){
						elementCount.put(element, primeCount + 1);
					}
					atomsToIgnore.add(atom);
				}
			}
		}
		
		{
			Set<String> elementsToIgnore = elementCount.keySet();
	
			for (Fragment fragment : allFragments) {
				List<Atom> atomList = fragment.getAtomList();
				for (Atom atom : atomList) {
					if (elementsToIgnore.contains(atom.getElement().toString())){
						atomsToIgnore.add(atom);
					}
				}
			}
		}

		String fragType = suffixableFragment.getType();
		if (fragType.equals(NONCARBOXYLICACID_TYPE_VAL) || fragType.equals(CHALCOGENACIDSTEM_TYPE_VAL)){
			if (suffixFragments.size() != 0){
				throw new StructureBuildingException("No suffix fragments were expected to be present on non carboxylic acid");
			}
			processNonCarboxylicAcidLabelling(suffixableFragment, elementCount, atomsToIgnore);
		}
		else{
			if (suffixFragments.size() > 0){
				processSuffixLabelling(suffixFragments, elementCount, atomsToIgnore);
				Integer seenCount = elementCount.get("N");
				if (seenCount != null && seenCount > 1){//look for special case violation of IUPAC rule, =(N)=(NN) is N//N' in practice rather than N/N'/N''
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
				if (atom.getElement() == ChemEl.N && atom.getIncomingValency() ==3 ){
					List<String> locants =atom.getLocants();
					if (locants.size()==1 && MATCH_ELEMENT_SYMBOL_LOCANT.matcher(locants.get(0)).matches()){
						List<Atom> neighbours = atom.getAtomNeighbours();
						for (Atom neighbour : neighbours) {
							if (neighbour.getElement() == ChemEl.N && neighbour.getIncomingValency()==1){
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

	private static void processMainGroupLabelling(Fragment suffixableFragment, Map<String, Integer> elementCount, Set<Atom> atomsToIgnore) {
		Set<String> elementToIgnore = new HashSet<String>(elementCount.keySet());
		List<Atom> atomList = suffixableFragment.getAtomList();
		Collections.sort(atomList, new SortAtomsForMainGroupElementSymbols());
		Atom atomToAddCLabelTo = null;//only add a C label if there is only one C in the main group
		boolean seenMoreThanOneC = false;
		for (Atom atom : atomList) {
			if (atomsToIgnore.contains(atom)){
				continue;
			}
			ChemEl chemEl = atom.getElement();
			if (elementToIgnore.contains(chemEl.toString())){
				continue;
			}
			if (chemEl == ChemEl.C) {
				if (seenMoreThanOneC) {
					continue;
				}
				if (atomToAddCLabelTo != null){
					atomToAddCLabelTo = null;
					seenMoreThanOneC = true;
				}
				else{
					atomToAddCLabelTo = atom;
				}
			}
			else{
				assignLocant(atom, elementCount);
			}
		}
		if (atomToAddCLabelTo != null){
			atomToAddCLabelTo.addLocant("C");
		}
	}

	private static void processSuffixLabelling(List<Fragment> suffixFragments, Map<String, Integer> elementCount, Set<Atom> atomsToIgnore) {
		List<Atom> startingAtoms = new ArrayList<Atom>();
		Set<Atom> atomsVisited = new HashSet<Atom>();
		for (Fragment fragment : suffixFragments) {
			Atom rAtom = fragment.getFirstAtom();
			List<Atom> nextAtoms = getIntraFragmentNeighboursAndSetVisitedBondOrder(rAtom);
			atomsVisited.addAll(nextAtoms);
			startingAtoms.addAll(nextAtoms);
		}
		Collections.sort(startingAtoms, new SortAtomsForElementSymbols());

		Deque<Atom> atomsToConsider = new ArrayDeque<Atom>(startingAtoms);
		while (atomsToConsider.size() > 0){
			assignLocantsAndExploreNeighbours(elementCount, atomsToIgnore, atomsVisited, atomsToConsider);
		}
	}

	private static void processNonCarboxylicAcidLabelling(Fragment suffixableFragment, Map<String, Integer> elementCount, Set<Atom> atomsToIgnore) {
		Set<Atom> atomsVisited = new HashSet<Atom>();
		Atom firstAtom = suffixableFragment.getFirstAtom();
		List<Atom> startingAtoms = getIntraFragmentNeighboursAndSetVisitedBondOrder(firstAtom);
		
		Collections.sort(startingAtoms, new SortAtomsForElementSymbols());
		atomsVisited.add(firstAtom);
		Deque<Atom> atomsToConsider = new ArrayDeque<Atom>(startingAtoms);
		while (atomsToConsider.size() > 0){
			assignLocantsAndExploreNeighbours(elementCount, atomsToIgnore, atomsVisited, atomsToConsider);
		}
		if (!atomsToIgnore.contains(firstAtom) && firstAtom.determineValency(true) > firstAtom.getIncomingValency()) {
			//e.g. carbonimidoyl the carbon has locant C
			assignLocant(firstAtom, elementCount);
		}
	}

	private static void assignLocantsAndExploreNeighbours(Map<String, Integer> elementCount, Set<Atom> atomsToIgnore, Set<Atom> atomsVisited, Deque<Atom> atomsToConsider) {
		Atom atom = atomsToConsider.removeFirst();
		atomsVisited.add(atom);
		if (!atomsToIgnore.contains(atom)) {//assign locant
			assignLocant(atom, elementCount);
		}
		List<Atom> atomsToExplore = getIntraFragmentNeighboursAndSetVisitedBondOrder(atom);
		atomsToExplore.removeAll(atomsVisited);
		Collections.sort(atomsToExplore, new SortAtomsForElementSymbols());
		for (int i = atomsToExplore.size() - 1; i >= 0; i--) {
			atomsToConsider.addFirst(atomsToExplore.get(i));
		}
	}

	/**
	 * Gets the neighbours of an atom that claim to be within the same frag
	 * The order of bond taken to get to the neighbour is set on the neighbours Atom.VISITED property
	 * @param atom
	 * @return
	 */
	private static List<Atom> getIntraFragmentNeighboursAndSetVisitedBondOrder(Atom atom) {
		List<Atom> atomsToExplore = new ArrayList<Atom>();
		List<Bond> bonds = atom.getBonds();
		for (Bond bond : bonds) {
			Atom neighbour = bond.getOtherAtom(atom);
			if (neighbour.getFrag().equals(atom.getFrag())) {
				atomsToExplore.add(neighbour);
				neighbour.setProperty(Atom.VISITED, bond.getOrder());
			}
		}
		return atomsToExplore;
	}

	private static void assignLocant(Atom atom, Map<String, Integer> elementCount) {
		String element = atom.getElement().toString();
		Integer count = elementCount.get(element);
		if (count == null){
			atom.addLocant(element);
			elementCount.put(element, 1);
		}
		else{
			atom.addLocant(element + StringTools.multiplyString("'", count));
			elementCount.put(element, count + 1);
		}
	}

	/** Adjusts the order of a bond in a fragment.
	 *
	 * @param fromAtom The lower-numbered atom in the bond
	 * @param bondOrder The new bond order
	 * @param fragment The fragment
	 * @return The bond that was unsaturated
     * @throws StructureBuildingException
	 */
	static Bond unsaturate(Atom fromAtom, int bondOrder, Fragment fragment) throws StructureBuildingException {
		Atom toAtom = null;
		Integer locant = null;
		try{
			String primes ="";
			String locantStr = fromAtom.getFirstLocant();
			int numberOfPrimes = StringTools.countTerminalPrimes(locantStr);
			locant = Integer.parseInt(locantStr.substring(0, locantStr.length()-numberOfPrimes));
			primes = StringTools.multiplyString("'", numberOfPrimes);
			Atom possibleToAtom = fragment.getAtomByLocant(String.valueOf(locant +1)+primes);
			if (possibleToAtom !=null && fromAtom.getBondToAtom(possibleToAtom)!=null){
				toAtom = possibleToAtom;
			}
			else if (possibleToAtom ==null && fromAtom.getAtomIsInACycle()){//allow something like cyclohexan-6-ene, something like butan-4-ene will still fail
				possibleToAtom = fragment.getAtomByLocant("1" + primes);
				if (possibleToAtom !=null && fromAtom.getBondToAtom(possibleToAtom)!=null){
					toAtom =possibleToAtom;
				}
			}
		}
		catch (Exception e) {
			List<Atom> atomList = fragment.getAtomList();
			int initialIndice = atomList.indexOf(fromAtom);
			if (initialIndice +1 < atomList.size() && fromAtom.getBondToAtom(atomList.get(initialIndice +1))!=null){
				toAtom = atomList.get(initialIndice +1);
			}
		}
		if (toAtom==null){
			if (locant!=null){
				throw new StructureBuildingException("Could not find bond to unsaturate starting from the atom with locant: " +locant);
			}
			else{
				throw new StructureBuildingException("Could not find bond to unsaturate");
			}
		}
		Bond b = fromAtom.getBondToAtomOrThrow(toAtom);
		b.setOrder(bondOrder);
		return b;
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
		Bond b = fromAtom.getBondToAtomOrThrow(toAtom);
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
			if (atomIndice +1>=fragment.getAtomCount()){
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
     * @param atomList
	 */
	static void relabelLocantsAsFusedRingSystem(List<Atom> atomList) {
		int locantVal = 0;
		char locantLetter = 'a';
		for (Atom atom : atomList) {
			atom.clearLocants();
		}
		for (Atom atom : atomList) {
			if(atom.getElement() != ChemEl.C || atom.getBondCount() < 3) {
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
				if (MATCH_NUMERIC_LOCANT.matcher(locant).matches()){
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
		if (nitrogen.getElement() == ChemEl.N && nitrogen.getAtomIsInACycle()){
			for (Atom neighbour : nitrogen.getAtomNeighbours()) {
				if (neighbour.hasSpareValency() && neighbour.getElement() == ChemEl.C && neighbour.getAtomIsInACycle()){
					List<Atom> distance2Neighbours = neighbour.getAtomNeighbours();
					distance2Neighbours.remove(nitrogen);
					for (Atom distance2Neighbour : distance2Neighbours) {
						if (distance2Neighbour.hasSpareValency() && distance2Neighbour.getElement() == ChemEl.N && distance2Neighbour.getAtomIsInACycle() && distance2Neighbour.getCharge()==0){
							return distance2Neighbour;
						}
					}
				}
			}
		}
		return null;
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
		List<Atom> originalIndicatedHydrogen = frag.getIndicatedHydrogen();
		List<Atom> indicatedHydrogen = new ArrayList<Atom>(originalIndicatedHydrogen);
		for (int i = indicatedHydrogen.size() -1; i >=0; i--) {
			if (!indicatedHydrogen.get(i).hasSpareValency()){
				indicatedHydrogen.remove(i);
			}
		}
		if (indicatedHydrogen.size()>0){
			if (indicatedHydrogen.size()>1){
				for (Atom indicatedAtom : indicatedHydrogen) {
					boolean couldBeInvolvedInSimpleNitrogenTautomerism = false;//fix for guanine like purine derivatives
					if (indicatedAtom.getElement() == ChemEl.N && indicatedAtom.getAtomIsInACycle()){
						atomloop : for (Atom neighbour : indicatedAtom.getAtomNeighbours()) {
							if (neighbour.getElement() == ChemEl.C && neighbour.getAtomIsInACycle()){
								List<Atom> distance2Neighbours = neighbour.getAtomNeighbours();
								distance2Neighbours.remove(indicatedAtom);
								for (Atom distance2Neighbour : distance2Neighbours) {
									if (distance2Neighbour.getElement() == ChemEl.N && distance2Neighbour.getAtomIsInACycle() && !originalIndicatedHydrogen.contains(distance2Neighbour)){
										couldBeInvolvedInSimpleNitrogenTautomerism =true;
										break atomloop;
									}
								}
							}
						}
					}
					//retain spare valency if has the cyclic [NH]C=N moiety but substitution has meant that this tautomerism doesn't actually occur cf. 8-oxoguanine
					if (!couldBeInvolvedInSimpleNitrogenTautomerism || detectSimpleNitrogenTautomer(indicatedAtom) != null){
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
								if (a.getElement() != ChemEl.C) {
									atomToReduceValencyAt = a;
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
								a.getBondToAtomOrThrow(aa).addOrder(1);
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
								a.getBondToAtomOrThrow(aa).addOrder(1);
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
									a.getBondToAtomOrThrow(aa).addOrder(1);
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


	static Atom getAtomByAminoAcidStyleLocant(Atom backboneAtom, String elementSymbol, String primes) {
		//Search for appropriate atom by using the same algorithm as is used to assign locants initially

		List<Atom> startingAtoms = new ArrayList<Atom>();
		Set<Atom> atomsVisited = new HashSet<Atom>();
		List<Atom> neighbours = getIntraFragmentNeighboursAndSetVisitedBondOrder(backboneAtom);
		mainLoop: for (Atom neighbour : neighbours) {
			atomsVisited.add(neighbour);
			if (!neighbour.getType().equals(SUFFIX_TYPE_VAL)){
				for (String neighbourLocant : neighbour.getLocants()) {
					if (MATCH_NUMERIC_LOCANT.matcher(neighbourLocant).matches()){//gone to an inappropriate atom
						continue mainLoop;
					}
				}
			}
			startingAtoms.add(neighbour);
		}

		Collections.sort(startingAtoms, new SortAtomsForElementSymbols());
		Map<String,Integer> elementCount = new HashMap<String,Integer>();//keeps track of how many times each element has been seen
	
		Deque<Atom> atomsToConsider = new ArrayDeque<Atom>(startingAtoms);
		boolean hydrazoneSpecialCase =false;//look for special case violation of IUPAC rule where the locant of the =N- atom is skipped. This flag is set when =N- is encountered
		while (atomsToConsider.size() > 0){
			Atom atom = atomsToConsider.removeFirst();
			atomsVisited.add(atom);
			int primesOnPossibleAtom =0;
			String element =atom.getElement().toString();
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

			List<Atom> atomNeighbours = getIntraFragmentNeighboursAndSetVisitedBondOrder(atom);
			atomNeighbours.removeAll(atomsVisited);
			for (int i = atomNeighbours.size() -1; i >=0; i--) {
				Atom neighbour = atomNeighbours.get(i);
				if (!neighbour.getType().equals(SUFFIX_TYPE_VAL)){
					for (String neighbourLocant : neighbour.getLocants()) {
						if (MATCH_NUMERIC_LOCANT.matcher(neighbourLocant).matches()){//gone to an inappropriate atom
							atomNeighbours.remove(i);
							break;
						}
					}
				}
			}
			if (atom.getElement() == ChemEl.N && atom.getIncomingValency() ==3 && atom.getCharge()==0 
					&& atomNeighbours.size()==1 && atomNeighbours.get(0).getElement() == ChemEl.N){
				hydrazoneSpecialCase =true;
			}
			else{
				if (element.equals(elementSymbol)){
					if (primes.length() == primesOnPossibleAtom){
						return atom;
					}
				}
			}

			Collections.sort(atomNeighbours, new SortAtomsForElementSymbols());
			for (int i = atomNeighbours.size() - 1; i >= 0; i--) {
				atomsToConsider.addFirst(atomNeighbours.get(i));
			}
		}

		if (primes.equals("") && backboneAtom.getElement().toString().equals(elementSymbol)){//maybe it meant the starting atom
			return backboneAtom;
		}
		return null;
	}
	
	
	/**
	 * Determines whether the bond between two elements is likely to be covalent
	 * This is crudely determined based on whether the combination of elements fall outside the ionic and
	 * metallic sections of a van Arkel diagram
	 * @param chemEl1
	 * @param chemEl2
	 * @return
	 */
	static boolean isCovalent(ChemEl chemEl1, ChemEl chemEl2) {
		Double atom1Electrongegativity = AtomProperties.getPaulingElectronegativity(chemEl1);
		Double atom2Electrongegativity = AtomProperties.getPaulingElectronegativity(chemEl2);
		if (atom1Electrongegativity!=null && atom2Electrongegativity !=null){
			double halfSum = (atom1Electrongegativity + atom2Electrongegativity)/2;
			double difference = Math.abs(atom1Electrongegativity - atom2Electrongegativity);
			if (halfSum < 1.6){
				return false;//probably metallic
			}
			if (difference < 1.76 * halfSum - 3.03){
				return true;
			}			
		}
		return false;
	}
	
	/**
	 * Is the atom a suffix atom/carbon of an aldehyde atom/chalcogen functional atom/hydroxy (or chalcogen equivalent)
	 * (by special step heterostems are not considered hydroxy e.g. disulfane)
	 * @param atom
	 * @return
	 */
	static boolean isCharacteristicAtom(Atom atom) {
		if (atom.getType().equals(SUFFIX_TYPE_VAL) || 
			(atom.getElement().isChalcogen() && !HETEROSTEM_SUBTYPE_VAL.equals(atom.getFrag().getSubType()) &&
			atom.getIncomingValency() == 1 &&
			atom.getOutValency() == 0 && atom.getCharge() == 0)) {
			return true;
		}
		return isFunctionalAtomOrAldehyde(atom);
	}
	
	/**
	 * Is the atom an aldehyde atom or a chalcogen functional atom
	 * @param atom
	 * @return
	 */
	static boolean isFunctionalAtomOrAldehyde(Atom atom) {
		if (Boolean.TRUE.equals(atom.getProperty(Atom.ISALDEHYDE))){//substituting an aldehyde would make it no longer an aldehyde
			return true;
		}
		return isFunctionalAtom(atom);
	}
	
	/**
	 * Is the atom a chalcogen functional atom
	 * @param atom
	 * @return
	 */
	static boolean isFunctionalAtom(Atom atom) {
		ChemEl chemEl = atom.getElement();
		if (chemEl.isChalcogen()) {//potential chalcogen functional atom
			Fragment frag = atom.getFrag();
			for (int i = 0, l = frag.getFunctionalAtomCount(); i < l; i++) {
				if (atom.equals(frag.getFunctionalAtom(i).getAtom())){
					return true;
				}
			}
		}
		return false;
	}


	/**
	 * Checks that all atoms in a ring appear to be equivalent
	 * @param ring
	 * @return true if all equivalent, else false
	 */
	static boolean  allAtomsInRingAreIdentical(Fragment ring){
		List<Atom> atomList = ring.getAtomList();
		Atom firstAtom = atomList.get(0);
		ChemEl chemEl = firstAtom.getElement();
		int valency = firstAtom.getIncomingValency();
		boolean spareValency = firstAtom.hasSpareValency();
		for (Atom atom : atomList) {
			if (atom.getElement() != chemEl){
				return false;
			}
			if (atom.getIncomingValency() != valency){
				return false;
			}
			if (atom.hasSpareValency() != spareValency){
				return false;
			}
		}
		return true;
	}


	/**
	 * Removes a terminal atom of a particular element e.g. oxygen
	 * A locant may be specified to indicate what atom is adjacent to the atom to be removed
	 * Formally the atom is replaced by hydrogen, hence stereochemistry is intentionally preserved
	 * @param state 
	 * @param fragment
	 * @param chemEl
	 * @param locant A locant or null
	 * @throws StructureBuildingException 
	 */
	static void removeHydroxyLikeTerminalAtom(BuildState state, Fragment fragment, ChemEl chemEl, String locant) throws StructureBuildingException {
		List<Atom> applicableTerminalAtoms;
		if (locant!=null){
			Atom adjacentAtom = fragment.getAtomByLocantOrThrow(locant);
			applicableTerminalAtoms = findHydroxyLikeTerminalAtoms(adjacentAtom.getAtomNeighbours(), chemEl);
			if (applicableTerminalAtoms.isEmpty()){
				throw new StructureBuildingException("Unable to find terminal atom of type: " + chemEl + " at locant "+ locant +" for subtractive nomenclature");
			}
		}
		else{
			applicableTerminalAtoms = findHydroxyLikeTerminalAtoms(fragment.getAtomList(), chemEl);
			if (applicableTerminalAtoms.isEmpty()){
				throw new StructureBuildingException("Unable to find terminal atom of type: " + chemEl + " for subtractive nomenclature");
			}
		}
		Atom atomToRemove = applicableTerminalAtoms.get(0);
		if (isFunctionalAtom(atomToRemove)){//This can occur with aminoglycosides where the anomeric OH is removed by deoxy
			for (int i = 0, l = fragment.getFunctionalAtomCount(); i < l; i++) {
				if (atomToRemove.equals(fragment.getFunctionalAtom(i).getAtom())){
					fragment.removeFunctionalAtom(i);
					break;
				}
			}
			fragment.addFunctionalAtom(atomToRemove.getFirstBond().getOtherAtom(atomToRemove));
		}
		removeTerminalAtom(state, atomToRemove);
	}
	
	static void removeTerminalAtom(BuildState state, Atom atomToRemove) {
		AtomParity atomParity =  atomToRemove.getAtomNeighbours().get(0).getAtomParity();
		if (atomParity!=null){//replace reference to atom with reference to implicit hydrogen
			Atom[] atomRefs4= atomParity.getAtomRefs4();
			for (int i = 0; i < atomRefs4.length; i++) {
				if (atomRefs4[i]==atomToRemove){
					atomRefs4[i] = AtomParity.deoxyHydrogen;
					break;
				}
			}
		}
		state.fragManager.removeAtomAndAssociatedBonds(atomToRemove);
	}


	/**
	 * Finds terminal atoms of the given element type from the list given
	 * The terminal atoms be single bonded, not radicals and uncharged
	 * @param atoms
	 * @param chemEl
	 * @return 
	 */
	static List<Atom> findHydroxyLikeTerminalAtoms(List<Atom> atoms, ChemEl chemEl) {
		List<Atom> matches =new ArrayList<Atom>();
		for (Atom atom : atoms) {
			if (atom.getElement() == chemEl && atom.getIncomingValency() == 1 &&
				atom.getOutValency() == 0 && atom.getCharge() == 0){
				matches.add(atom);
			}
		}
		return matches;
	}

	/**
	 * Checks whether a bond is part of a 6 member or smaller ring.
	 * This is necessary as such double bonds are assumed to not be capable of having E/Z stereochemistry
	 * @param bond
	 * @return true unless in a 6 member or smaller rings
	 */
	static boolean notIn6MemberOrSmallerRing(Bond bond) {
		Atom fromAtom =bond.getFromAtom();
		Atom toAtom = bond.getToAtom();
		if (fromAtom.getAtomIsInACycle() && toAtom.getAtomIsInACycle()){//obviously both must be in rings
			//attempt to get from the fromAtom to the toAtom in 6 or fewer steps.
			List<Atom> visitedAtoms = new ArrayList<Atom>();
			Deque<Atom> atomsToInvestigate = new ArrayDeque<Atom>();//A queue is not used as I need to make sure that only up to depth 6 is investigated
			List<Atom> neighbours =fromAtom.getAtomNeighbours();
			neighbours.remove(toAtom);
			for (Atom neighbour : neighbours) {
				atomsToInvestigate.add(neighbour);
			}
			visitedAtoms.add(fromAtom);
			for (int i = 0; i < 5; i++) {//up to 5 bonds from the neighbours of the fromAtom i.e. up to ring size 6
				if (atomsToInvestigate.isEmpty()){
					break;
				}
				Deque<Atom> atomsToInvestigateNext = new ArrayDeque<Atom>();
				while (!atomsToInvestigate.isEmpty()) {
					Atom currentAtom =atomsToInvestigate.removeFirst();
					if (currentAtom == toAtom){
						return false;
					}
					visitedAtoms.add(currentAtom);
					neighbours =currentAtom.getAtomNeighbours();
					for (Atom neighbour : neighbours) {
						if (!visitedAtoms.contains(neighbour) && neighbour.getAtomIsInACycle()){
							atomsToInvestigateNext.add(neighbour);
						}
					}
				}
				atomsToInvestigate = atomsToInvestigateNext;
			}
		}
		return true;
	}
	
	/**
	 * Finds the hydroxy atom of all hydroxy functional groups in a fragment
	 * i.e. not in carboxylic acid or oxime
	 * @param frag
	 * @return
	 * @throws StructureBuildingException 
	 */
	static List<Atom> findHydroxyGroups(Fragment frag) throws StructureBuildingException {
		List<Atom> hydroxyAtoms = new ArrayList<Atom>();
		List<Atom> atoms = frag.getAtomList();
		for (Atom atom : atoms) {
			if (atom.getElement() == ChemEl.O && atom.getIncomingValency() == 1 && atom.getOutValency() == 0 && atom.getCharge() == 0){
				Atom adjacentAtom = atom.getAtomNeighbours().get(0);
				List<Atom> neighbours = adjacentAtom.getAtomNeighbours();
				if (adjacentAtom.getElement() == ChemEl.C){
					neighbours.remove(atom);
					if (neighbours.size() >= 1 && neighbours.get(0).getElement() == ChemEl.O && adjacentAtom.getBondToAtomOrThrow(neighbours.get(0)).getOrder()==2){
						continue;
					}
					if (neighbours.size() >= 2 && neighbours.get(1).getElement() == ChemEl.O && adjacentAtom.getBondToAtomOrThrow(neighbours.get(1)).getOrder()==2){
						continue;
					}
					hydroxyAtoms.add(atom);
				}
			}
		}
		return hydroxyAtoms;
	}
	
	static List<Atom> findnAtomsForSubstitution(List<Atom> atomList, Atom preferredAtom, int numberOfSubstitutionsRequired, int bondOrder, boolean takeIntoAccountOutValency) {
		int atomCount = atomList.size();
		int startingIndex = preferredAtom != null ? atomList.indexOf(preferredAtom) : 0;
		if (startingIndex < 0){
			throw new IllegalArgumentException("OPSIN Bug: preferredAtom should be part of the list of atoms to search through");
		}
		CyclicAtomList atoms = new CyclicAtomList(atomList, startingIndex - 1);//next() will retrieve the atom at the startingIndex
		List<Atom> substitutableAtoms = new ArrayList<Atom>();
		for (int i = 0; i < atomCount; i++) {//aromaticity preserved, standard valency assumed, characteristic atoms ignored
			Atom atom = atoms.next();
			if (!FragmentTools.isCharacteristicAtom(atom) || (numberOfSubstitutionsRequired == 1 && atom == preferredAtom)) {
				int currentExpectedValency = atom.determineValency(takeIntoAccountOutValency);
				int usedValency = atom.getIncomingValency() + (atom.hasSpareValency() ? 1 : 0) + (takeIntoAccountOutValency ? atom.getOutValency() : 0);
				int timesAtomCanBeSubstitued = ((currentExpectedValency - usedValency)/ bondOrder);
				for (int j = 1; j <= timesAtomCanBeSubstitued; j++) {
					substitutableAtoms.add(atom);
				}
			}
		}
		if (substitutableAtoms.size() >= numberOfSubstitutionsRequired){
			return substitutableAtoms;
		}
		substitutableAtoms.clear();
		for (int i = 0; i < atomCount; i++) {//aromaticity preserved, standard valency assumed, functional suffixes ignored
			Atom atom = atoms.next();
			if (!FragmentTools.isFunctionalAtomOrAldehyde(atom) || (numberOfSubstitutionsRequired == 1 && atom == preferredAtom)) {
				int currentExpectedValency = atom.determineValency(takeIntoAccountOutValency);
				int usedValency = atom.getIncomingValency() + (atom.hasSpareValency() ? 1 : 0) + (takeIntoAccountOutValency ? atom.getOutValency() : 0);
				int timesAtomCanBeSubstitued = ((currentExpectedValency - usedValency)/ bondOrder);
				for (int j = 1; j <= timesAtomCanBeSubstitued; j++) {
					substitutableAtoms.add(atom);
				}
			}
		}
		if (substitutableAtoms.size() >= numberOfSubstitutionsRequired){
			return substitutableAtoms;
		}
		substitutableAtoms.clear();
		
		for (int i = 0; i < atomCount; i++) {//aromaticity preserved, any sensible valency allowed, anything substitutable
			Atom atom = atoms.next();
			Integer maximumValency = ValencyChecker.getMaximumValency(atom);
			if (maximumValency != null) {
				int usedValency = atom.getIncomingValency() + (atom.hasSpareValency() ? 1 : 0) + (takeIntoAccountOutValency ? atom.getOutValency() : 0);
				int timesAtomCanBeSubstitued = ((maximumValency - usedValency)/ bondOrder);
				for (int j = 1; j <= timesAtomCanBeSubstitued; j++) {
					substitutableAtoms.add(atom);
				}
			}
			else{
				for (int j = 0; j < numberOfSubstitutionsRequired; j++) {
					substitutableAtoms.add(atom);
				}
			}
		}
		if (substitutableAtoms.size() >= numberOfSubstitutionsRequired){
			return substitutableAtoms;
		}
		substitutableAtoms.clear();

		for (int i = 0; i < atomCount; i++) {//aromaticity dropped, any sensible valency allowed, anything substitutable
			Atom atom = atoms.next();
			Integer maximumValency = ValencyChecker.getMaximumValency(atom);
			if (maximumValency != null) {
				int usedValency = atom.getIncomingValency() + (takeIntoAccountOutValency ? atom.getOutValency() : 0);
				int timesAtomCanBeSubstitued = ((maximumValency - usedValency)/ bondOrder);
				for (int j = 1; j <= timesAtomCanBeSubstitued; j++) {
					substitutableAtoms.add(atom);
				}
			}
			else {
				for (int j = 0; j < numberOfSubstitutionsRequired; j++) {
					substitutableAtoms.add(atom);
				}
			}
		}
		if (substitutableAtoms.size() >= numberOfSubstitutionsRequired){
			return substitutableAtoms;
		}
		return null;
	}
	
	static List<Atom> findnAtomsForSubstitution(Fragment frag, Atom preferredAtom, int numberOfSubstitutionsRequired, int bondOrder, boolean takeIntoAccountOutValency) {
		return findnAtomsForSubstitution(frag.getAtomList(), preferredAtom, numberOfSubstitutionsRequired, bondOrder, takeIntoAccountOutValency);
	}
	
	/**
	 * Returns a list of atoms of size >= numberOfSubstitutionsDesired (or null if this not possible)
	 * An atom must have have sufficient valency to support a substituent requiring a bond of order bondOrder
	 * If an atom can support multiple substituents it will appear in the list multiple times
	 * This method iterates over the the fragment atoms attempting to fulfil these requirements with incrementally more lenient constraints:
	 * aromaticity preserved, standard valency assumed, characteristic atoms ignored
	 * aromaticity preserved, standard valency assumed, functional suffixes ignored
	 * aromaticity preserved, any sensible valency allowed, anything substitutable
	 * aromaticity dropped, any sensible valency allowed, anything substitutable
	 * 
	 * Iteration starts from the defaultInAtom (if applicable, else the first atom) i.e. the defaultInAtom if substitutable will be the first atom in the list
	 * @param frag
	 * @param numberOfSubstitutionsRequired
	 * @param bondOrder
	 * @return
	 */
	static List<Atom> findnAtomsForSubstitution(Fragment frag, int numberOfSubstitutionsRequired, int bondOrder) {
		return findnAtomsForSubstitution(frag.getAtomList(), frag.getDefaultInAtom(), numberOfSubstitutionsRequired, bondOrder, true);
	}
	
	/**
	 * Returns a list of the most preferable atoms for substitution (empty list if none are)
	 * An atom must have have sufficient valency to support a substituent requiring a bond of order bondOrder
	 * If an atom can support multiple substituents it will appear in the list multiple times
	 * This method iterates over the the fragment atoms attempting to fulfil these requirements with incrementally more lenient constraints:
	 * aromaticity preserved, standard valency assumed, characteristic atoms ignored
	 * aromaticity preserved, standard valency assumed, functional suffixes ignored
	 * aromaticity preserved, any sensible valency allowed, anything substitutable
	 * aromaticity dropped, any sensible valency allowed, anything substitutable
	 * 
	 * Iteration starts from the defaultInAtom (if applicable, else the first atom) i.e. the defaultInAtom if substitutable will be the first atom in the list
	 * @param frag
	 * @param bondOrder
	 * @return
	 */
	static List<Atom> findSubstituableAtoms(Fragment frag, int bondOrder) {
		List<Atom> potentialAtoms = findnAtomsForSubstitution(frag, 1, bondOrder);
		if (potentialAtoms == null) {
			return Collections.emptyList();
		}
		return potentialAtoms;
	}
	
	/**
	 * Returns the first substitutable atom found using the same criteria as {@link FragmentTools#findnAtomsForSubstitution(Fragment, int, int)}
	 * The defaultInAtom in this case is the provided preferredAtom
	 * takeIntoAccountOutValency may be set to false for cases where the position of radicals is being determined
	 * Return null if no suitable atom can be found
	 * @param frag
	 * @param preferredAtom
	 * @param bondOrder
	 * @param takeIntoAccountOutValency
	 * @return
	 */
	static Atom findAtomForSubstitution(Fragment frag, Atom preferredAtom, int bondOrder, boolean takeIntoAccountOutValency) {
		List<Atom> atoms = findnAtomsForSubstitution(frag.getAtomList(), preferredAtom, 1, bondOrder, takeIntoAccountOutValency);
		if (atoms != null) {
			return atoms.get(0);
		}
		return null;
	}
	
	/**
	 * Returns the first substitutable atom found using the same criteria as {@link FragmentTools#findnAtomsForSubstitution(Fragment, int, int)}
	 * The defaultInAtom in this case is the provided preferredAtom
	 * takeIntoAccountOutValency may be set to false for cases where the position of radicals is being determined
	 * Throws an exception if no suitable atom can be found
	 * @param frag
	 * @param preferredAtom
	 * @param bondOrder
	 * @param takeIntoAccountOutValency
	 * @return
	 * @throws StructureBuildingException 
	 */
	static Atom findAtomForSubstitutionOrThrow(Fragment frag, Atom preferredAtom, int bondOrder, boolean takeIntoAccountOutValency) throws StructureBuildingException {
		Atom a = findAtomForSubstitution(frag, preferredAtom, bondOrder, takeIntoAccountOutValency);
		if (a == null) {
			throw new StructureBuildingException("No suitable atom found");
		}
		return a;
	}

}
