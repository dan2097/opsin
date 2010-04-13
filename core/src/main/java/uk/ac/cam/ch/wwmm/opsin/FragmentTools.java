package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

class FragmentTools {
	/**
	 * Sorts a list of atoms such that their order agrees with the order symbolic locants are typically assigned
	 * @author dl387
	 *
	 */
	static class SortAtomsForElementSymbols implements Comparator<Atom> {

	    public int compare(Atom a, Atom b){
	    	int compare =a.getElement().compareTo(b.getElement());
	    	if (compare !=0){//only bother comparing properly if elements are the same
	    		return compare;
	    	}

	    	if (a.getBonds().size() >b.getBonds().size()){//less bonds is preferred
	    		return 1;
	    	}
	    	if (a.getBonds().size() <b.getBonds().size()){
	    		return -1;
	    	}

	    	int aTotalBondOrder =0;
	    	for (Bond bonda : a.getBonds()) {
				aTotalBondOrder +=bonda.getOrder();
			}

	    	int bTotalBondOrder =0;
	    	for (Bond bondb : b.getBonds()) {
				bTotalBondOrder +=bondb.getOrder();
			}

	    	aTotalBondOrder +=a.getOutValency();//take into account the bond/s which hasn't been added yet which will go to a main fragment
	    	bTotalBondOrder +=b.getOutValency();//take into account the bond/s which hasn't been added yet which will go to a main fragment

	    	if (aTotalBondOrder >bTotalBondOrder){//lower order bonds are preferred
	    		return 1;
	    	}
	    	if (aTotalBondOrder <bTotalBondOrder){
	    		return -1;
	    	}

	    	return 0;
	    }
	}
	/**
	 * Sorts by number, then by letter e.g. 4,3,3b,5,3a,2 -->2,3,3a,3b,4,5
	 * @author dl387
	 *
	 */
	static class SortByLocants implements Comparator<Atom> {
		static final Pattern matchdigits = Pattern.compile("(\\d+).*");
    	static final Pattern matchletters = Pattern.compile(".*([a-z]+)");

	    public int compare(Atom atoma, Atom atomb){
	    	String locanta =atoma.getFirstLocant();
	    	String locantb =atomb.getFirstLocant();

	    	Matcher m1  =matchdigits.matcher(locanta);
	    	int locantaNumber=0;
	    	if (m1.matches()){
	    		locantaNumber=Integer.parseInt(m1.group(1));
	    	}
	    	else{
	    		return 0;//invalid locant (could be intentionally invalid)
	    	}

	    	Matcher m2  =matchdigits.matcher(locantb);
	    	int locantbNumber=0;
	    	if (m2.matches()){
	    		locantbNumber=Integer.parseInt(m2.group(1));
	    	}
	    	else{
	    		return 0;//invalid locant (could be intentionally invalid)
	    	}

	        if (locantaNumber >locantbNumber) {
	            return 1;//e.g. 3 vs 2 or 3a vs 2
	        } else if (locantbNumber >locantaNumber) {
	            return -1;//e.g. 2 vs 3 or 2 vs 3a
	        }
	        else{
	        	m1  =matchletters.matcher(locanta);
	        	String locantaLetter="";
	        	if (m1.matches()){
	        		locantaLetter=m1.group(1);
	        	}
	        	else{
	        		return -1;// e.g. 1 vs 1a
	        	}

	        	m2  =matchletters.matcher(locantb);
	        	String locantbLetter="";
	        	if (m2.matches()){
	        		locantbLetter=m2.group(1);
	        	}
	        	else{
	        		return 1;//e.g. 1a vs 1
	        	}

	            if (locantaLetter.compareTo(locantbLetter)>=1) {
	                return 1;//e.g. 1b vs 1a
	            } else if (locantbLetter.compareTo(locantaLetter)>=1) {
	                return -1;//e.g. 1a vs 1b
	            }
	            return 0;
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
	 */
	static void assignElementLocants(Fragment suffixableFragment, ArrayList<Fragment> suffixFragments) {
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

		for (Fragment fragment : suffixFragments) {
			List<Atom> atomList =fragment.getAtomList();

			/*
			 * Sort them by empirical rules to agree with IUPAC numbering e.g. single bonded atoms are preffered to doubled bonded atoms
			 */
			Collections.sort(atomList, new SortAtomsForElementSymbols());
			for (Atom atom : atomList) {//add the locants
				if (atomsToIgnore.contains(atom)){continue;}
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
		}
		HashSet<String> elementToIgnore = new HashSet<String>(elementCount.keySet());
		elementCount =new HashMap<String,Integer>();
		List<Atom> atomList =suffixableFragment.getAtomList();
		Atom atomToAddCLabelTo=null;//only add a C label if there is only one C in the main group
		for (Atom atom : atomList) {
			if (atomsToIgnore.contains(atom)){continue;}
			String element =atom.getElement();
			if (elementToIgnore.contains(element)){
				continue;
			}
			if (element.equals("C")){
				if (atomToAddCLabelTo !=null){
					elementToIgnore.add("C");
					atomToAddCLabelTo=null;
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
	

	/** Adjusts the order of a bond in a fragment.
	 *
	 * @param fromAtomID The id of the lower-numbered atom in the bond
	 * @param bondOrder The new bond order
	 * @param fragment The fragment
	 */
	static void unsaturate(int fromAtomID, int bondOrder, Fragment fragment) throws StructureBuildingException {
		int toAtomID = fromAtomID + 1;
		if (fragment.getAtomByID(toAtomID)==null || fragment.getAtomByID(toAtomID).getType().equals("suffix")){//allows something like cyclohexan-6-ene, something like butan-4-ene will still fail
			List<Atom> neighbours =fragment.getAtomByIDOrThrow(fromAtomID).getAtomNeighbours();
			if (neighbours.size() >=2){
				int firstID =fragment.getIdOfFirstAtom();
				for (Atom a : neighbours) {
					if (a.getID() ==firstID){
						toAtomID=firstID;
						break;
					}
				}
			}
		}
		Bond b = fragment.findBondOrThrow(fromAtomID, toAtomID);
		b.setOrder(bondOrder);
	}

	/** Adjusts the order of a bond in a fragment.
	 *
	 * @param fromAtomID The id of the first atom in the bond
	 * @param locantTo The locant of the other atom in the bond
	 * @param bondOrder The new bond order
	 * @param fragment The fragment
	 */
	static void unsaturate(int fromAtomID, String locantTo, int bondOrder, Fragment fragment) throws StructureBuildingException {
		int toAtomID = fragment.getIDFromLocantOrThrow(locantTo);
		Bond b = fragment.findBondOrThrow(fromAtomID, toAtomID);
		b.setOrder(bondOrder);
	}
	
	/** Works out where to put an "one", if this is unspecified. position 2 for propanone
	 * and higher, else 1. Position 2 is assumed to be 1 higher than the atomIndice given.
	 *
	 * @param fragment The fragment
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
	 */
	static void relabelFusedRingSystem(Fragment fusedring){
		relabelFusedRingSystem(fusedring.getAtomList());
	}

	/**Adjusts the labeling on a fused ring system, such that bridgehead atoms
	 * have locants endings in 'a' or 'b' etc. Example: naphthalene
	 * 1,2,3,4,5,6,7,8,9,10->1,2,3,4,4a,5,6,7,8,8a
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
	 * Adds the given string to all the first locants of the atoms.
	 * Other locants are removed
	 * @param atomList
	 * @param stringToAdd
	 */
	static void relabelLocants(List<Atom> atomList, String stringToAdd) {
		for (Atom atom : atomList) {
			atom.replaceLocant(atom.getFirstLocant() + stringToAdd);
		}
	}
}
