package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import nu.xom.Element;
import nu.xom.Elements;

/**
 *
 * @author aa593/dl387
 *
 */
class FusedRingBuilder {
	/**
	 * Sorts by number, then by letter e.g. 4,3,3b,5,3a,2 -->2,3,3a,3b,4,5
	 * @author dl387
	 *
	 */
	class SortByLocants implements Comparator<Atom> {
		Pattern matchdigits = Pattern.compile("(\\d+).*");
    	Pattern matchletters = Pattern.compile(".*([a-z]+)");

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
	private Pattern matchSlash = Pattern.compile("/");
	private Pattern matchC = Pattern.compile("C");

	FusedRingBuilder() {
	}

	/**
	 * Master method for processing fused rings. If 2 groups are present will attempt to fuse them
	 * Returns the substituent/root with the 2 groups fused together into 1 group
	 * @param state: contains the current id and fragment manager
	 * @param subOrRoot Element (substituent or root)
	 * @throws PostProcessingException
	 * @throws StructureBuildingException
	 */
	Element processFusedRings(BuildState state, Element subOrRoot) throws PostProcessingException, StructureBuildingException {
		Elements groups =subOrRoot.getChildElements("group");
		if (groups.size() < 2){return subOrRoot;}//nothing to fuse
		if (groups.size() >2){
			throw new PostProcessingException("Unsupported fused ring system: more than 2 groups");
		}
		Elements fusions =subOrRoot.getChildElements("fusion");
		if (fusions.size() >1){
			throw new PostProcessingException("Unsupported fused ring system: more than 1 fusion");
		}
		List<Fragment> rings = new ArrayList<Fragment>();

		/*
		 * Apply any nonstandard ring numbering
		 */
		for(int i=0;i<groups.size();i++) {
			Element group=groups.get(i);
			Fragment ring =state.xmlFragmentMap.get(group);
			if (i==groups.size()-1){
				//perform a quick check that every atom in this group is infact cyclic. Fusion components are enumerated and hence all guarenteed to be purely cyclic
				List<Atom> atomList = ring.getAtomList();
				for (Atom atom : atomList) {
					if (!atom.getAtomIsInACycle()){
						throw new PostProcessingException("Inappropriate group used in fusion nomenclature. Only groups composed entirely of atoms in cycles may be used. i.e. not: " + group.getValue());
					}
				}
				if (group.getAttribute("fusedRingNumbering")!=null){
					String[] standardNumbering = matchSlash.split(group.getAttributeValue("fusedRingNumbering"),-1);
					for (int j = 0; j < standardNumbering.length; j++) {
						atomList.get(j).replaceLocant(standardNumbering[j]);
					}
				}
			}
			rings.add(ring);
		}

		List<String> numericalLocantsOfChild = null;
		List<String> letterLocantsOfParent = null;
		if (fusions.size() > 0){
			Element fusion = fusions.get(0);
			String[] fusionArray=fusion.getValue().split("-");
			if (fusionArray.length ==2){
				String[] locantsOfChildTemp = fusionArray[0].replaceFirst("\\[", "").split(",");
				numericalLocantsOfChild = Arrays.asList(locantsOfChildTemp);
				char[] tempLetterLocantsOfParent = fusionArray[1].replaceFirst("\\]", "").toCharArray();
				letterLocantsOfParent = new ArrayList<String>();
				for (int i = 0; i < tempLetterLocantsOfParent.length; i++) {
					letterLocantsOfParent.add(String.valueOf(tempLetterLocantsOfParent[i]));
				}
			}
			else{
				String tempContents = fusionArray[0].replaceFirst("\\[", "").replaceFirst("\\]", "");
				if (tempContents.contains(",")){//only has digits
					String[] numericalLocantsOfChildTemp =tempContents.split(",");
					numericalLocantsOfChild = Arrays.asList(numericalLocantsOfChildTemp);
				}
				else{//only has letters
					char[] tempLetterLocantsOfParentCharArray = tempContents.toCharArray();
					letterLocantsOfParent = new ArrayList<String>();
					for (int i = 0; i < tempLetterLocantsOfParentCharArray.length; i++) {
						letterLocantsOfParent.add(String.valueOf(tempLetterLocantsOfParentCharArray[i]));
					}
				}
			}
		}

		int edgeLength =1;
		if (numericalLocantsOfChild != null){
			if (numericalLocantsOfChild.size() <=1){
				throw new StructureBuildingException("At least two numerical locants must be provided to perform fusion!");
			}
			edgeLength = numericalLocantsOfChild.size()-1;
		}
		else if (letterLocantsOfParent != null){
			edgeLength = letterLocantsOfParent.size();
		}

		if (numericalLocantsOfChild == null){
			numericalLocantsOfChild = findPossibleNumericalLocants(rings.get(0), edgeLength);
		}

		if (letterLocantsOfParent == null){
			letterLocantsOfParent = findPossibleLetterLocants(rings.get(1), edgeLength);
		}
		if (numericalLocantsOfChild == null || letterLocantsOfParent ==null){
			throw new StructureBuildingException("Unable to find bond to form fused ring system. Some information for forming fused ring system was only supplyed implicitly");
		}

		Fragment fusedRing =fuseRings(state, rings.get(0), rings.get(1), numericalLocantsOfChild, letterLocantsOfParent);//fuse the rings
		String fusedRingName=groups.get(0).getValue();
		if (fusedRingName.equals("benz") || fusedRingName.equals("benzo")){
			benzoSpecificAssignHeteroAtomsUsingLocants(state, groups.get(0), fusedRing);
		}
		for(int i=0;i<fusions.size();i++) {
			fusedRingName+=fusions.get(i).getValue();
		}
		fusedRingName+=groups.get(1).getValue();

		Element fusedRingEl =groups.get(groups.size()-1);//reuse this element to save having to remap suffixes...
		fusedRingEl.getAttribute("value").setValue(fusedRingName);
		fusedRingEl.getAttribute("valType").setValue("generatedFragment");
		fusedRingEl.getAttribute("type").setValue("ring");
		fusedRingEl.getAttribute("subType").setValue("fusedRing");
		fusedRingEl.removeChildren();
		fusedRingEl.appendChild(fusedRingName);

		state.xmlFragmentMap.put(fusedRingEl, fusedRing);

		for(int i=0;i<groups.size() -1;i++) {
			groups.get(i).detach();
		}
		for(int i=0;i<fusions.size();i++) {
			groups.get(i).detach();
		}
		return fusedRingEl;
	}

	/**
	 * Takes a ring an returns and array with one letter corresponding to a side/s
	 * that contains two adjacent non bridgehead carbons
	 * The number of sides is specified by edgeLength
	 * @param ring
	 * @param edgeLength The number of bonds to be fused along
	 * @return
	 */
	private List<String> findPossibleLetterLocants(Fragment ring, int edgeLength) {
		List<Atom> atomlist = sortFragmentAtomListByLocant(ring);
		List<String> letterLocantsOfParent = null;
		List<Atom> carbonAtoms = new ArrayList<Atom>();
		atomlist.add(0, atomlist.get(atomlist.size()-1));//this atomList is a copy so we can safely do this
		for (int i =atomlist.size() -1; i >=0; i--) {//iterate backwards in list to use highest locanted edge in preference.
			//this retains what is currently locant 1 on the parent ring as locant 1 if the first two atoms found match
			Atom atom = atomlist.get(i);
			if (matchC.matcher(atom.getElement()).matches()){
				if (atom.getIncomingValency()>=3){
					carbonAtoms.clear();
					continue;//don't want bridgehead carbons
				}
				carbonAtoms.add(atom);
				if (carbonAtoms.size() ==edgeLength +1 ){//as many in a row as edgelength ->use this side
					letterLocantsOfParent = new ArrayList<String>();
					Collections.reverse(carbonAtoms);
					atomlist.remove(0);
					for (int j = 0; j < edgeLength; j++) {
						letterLocantsOfParent.add(String.valueOf((char)(97 +atomlist.indexOf(carbonAtoms.get(j)))));//97 is ascii for a	
					}
					break;
				}
			}
			else{
				carbonAtoms.clear();
			}
		}
		return letterLocantsOfParent;
	}

	/**
	 * Takes a ring and returns an array of numbers corresponding to a side/s
	 * that contains two adjacent non bridgehead carbons
	 * The number of sides is specified by edgeLength
	 * @param ring
	 * @param edgeLength The number of bonds to be fused along
	 * @return
	 */
	private List<String> findPossibleNumericalLocants(Fragment ring, int edgeLength) {
		List<Atom> atomlist = sortFragmentAtomListByLocant(ring);
		List<String> numericalLocantsOfChild = null;
		List<String> carbonLocants = new ArrayList<String>();
		atomlist.add(atomlist.get(0));//this atomList is a copy so we can safely do this
		for (Atom atom : atomlist) {
			if (matchC.matcher(atom.getElement()).matches()){
				if (atom.getIncomingValency()>=3){
					carbonLocants.clear();
					continue;//don't want bridgehead carbons
				}
				carbonLocants.add(atom.getFirstLocant());
				if (carbonLocants.size()==edgeLength +1){//as many in a row as edgelength ->use this side
					numericalLocantsOfChild = new ArrayList<String>();
					for (String locant : carbonLocants) {
						numericalLocantsOfChild.add(locant);
					}
					break;
				}
			}
			else{
				carbonLocants.clear();
			}
		}
		return numericalLocantsOfChild;
	}

	/**
	 * Performs a single ring fusion using the values in numericalLocantsOfChild/letterLocantsOfParent
	 * @param state
	 * @param childRing
	 * @param parentRing
	 * @param numericalLocantsOfChild
	 * @param letterLocantsOfParent
	 * @return The fused ring fragment. Locants are not representative of final fused ring numbering
	 * @throws StructureBuildingException
	 */
	private Fragment fuseRings(BuildState state, Fragment childRing, Fragment parentRing, List<String> numericalLocantsOfChild, List<String> letterLocantsOfParent) throws StructureBuildingException {
		List<Atom> sortedAtomsInChild =sortFragmentAtomListByLocant(childRing);
		int indexfirst = sortedAtomsInChild.indexOf(childRing.getAtomByLocantOrThrow(numericalLocantsOfChild.get(0)));
		int indexfinal = sortedAtomsInChild.indexOf(childRing.getAtomByLocantOrThrow(numericalLocantsOfChild.get(numericalLocantsOfChild.size()-1)));
		CyclicAtomList cyclicListAtomsInChild = new CyclicAtomList(sortedAtomsInChild, indexfirst);
		List<Atom> childAtoms = null;
		
		List<Atom> possibleChildAtoms = new ArrayList<Atom>();
		possibleChildAtoms.add(cyclicListAtomsInChild.getCurrent());
		while (cyclicListAtomsInChild.getIndice() != indexfinal){//assume numbers are ascending
			possibleChildAtoms.add(cyclicListAtomsInChild.getNext());
		}
		if (letterLocantsOfParent.size() +1 == possibleChildAtoms.size()){
			boolean notInPossibleChildAtoms =false;
			for (int i =1; i < numericalLocantsOfChild.size()-1 ; i ++){
				if (!possibleChildAtoms.contains(childRing.getAtomByLocantOrThrow(numericalLocantsOfChild.get(i)))){
					notInPossibleChildAtoms =true;
				}
			}
			if (!notInPossibleChildAtoms){
				childAtoms = possibleChildAtoms;
			}
		}
		
		if (childAtoms ==null){//that didn't work, so try assuming the numbers are descending
			cyclicListAtomsInChild.setIndice(indexfirst);
			possibleChildAtoms.clear();
			possibleChildAtoms.add(cyclicListAtomsInChild.getCurrent());
			while (cyclicListAtomsInChild.getIndice() != indexfinal){//assume numbers are ascending
				possibleChildAtoms.add(cyclicListAtomsInChild.getPrevious());
			}
			if (letterLocantsOfParent.size() +1 == possibleChildAtoms.size()){
				boolean notInPossibleChildAtoms =false;
				for (int i =1; i < numericalLocantsOfChild.size()-1 ; i ++){
					if (!possibleChildAtoms.contains(childRing.getAtomByLocantOrThrow(numericalLocantsOfChild.get(i)))){
						notInPossibleChildAtoms =true;
					}
				}
				if (!notInPossibleChildAtoms){
					childAtoms = possibleChildAtoms;
				}
			}
		}
		if (childAtoms ==null){
			throw new StructureBuildingException("Malformed fusion bracket!");
		}

		List<List<Atom>> neighboursOfToBeReplacedChildAtoms= new ArrayList<List<Atom>>();//this list is in the same order as childAtoms
		for (Atom atom : childAtoms) {
			List<Atom> neighboursToBeAttachedToParentRing = new ArrayList<Atom>();
			List<Atom> neighbours = childRing.getAtomNeighbours(atom);
			for (Atom neighbour : neighbours) {
				if (!childAtoms.contains(neighbour)){
					neighboursToBeAttachedToParentRing.add(neighbour);
				}
			}
			neighboursOfToBeReplacedChildAtoms.add(neighboursToBeAttachedToParentRing);
		}
		//remove the childAtoms
		for (Atom atom : childAtoms) {
			state.fragManager.removeAtomAndAssociatedBonds(atom);
		}

		List<Atom> parentAtoms = new ArrayList<Atom>();
		List<Atom> sortedAtomsInParent =sortFragmentAtomListByLocant(parentRing);
		CyclicAtomList cyclicListAtomsInParent = new CyclicAtomList(sortedAtomsInParent, (int)letterLocantsOfParent.get(0).charAt(0) -97);//convert from lower case character through ascii to 0-23
		parentAtoms.add(cyclicListAtomsInParent.getCurrent());
		for (int i = 0; i < letterLocantsOfParent.size(); i++) {
			parentAtoms.add(cyclicListAtomsInParent.getNext());
		}
		if (parentAtoms.size()!=childAtoms.size()){
			throw new StructureBuildingException("Problem with fusion descriptors: Parent atoms specified: " + parentAtoms.size() +" Child atoms specified: " + childAtoms.size() + " These should have been identical!");
		}

		parentRing.importFrag(childRing);
		state.fragManager.removeFragment(childRing);
		for (int i = 0; i < parentAtoms.size(); i++) {
			Atom parentAtom = parentAtoms.get(i);
			if (!parentAtom.getElement().equals(childAtoms.get(i).getElement())){
				throw new StructureBuildingException("Invalid fusion descriptor: Heteroatom placement is ambigous as it is not present in both components of the fusion");
			}
			for (Atom atom : neighboursOfToBeReplacedChildAtoms.get(i)) {
				//System.out.println("Atom ID " + atom.getID() +" bonded to " +  parentAtom.get(i));
				state.fragManager.createBond(parentAtom, atom, 1);
			}
		}
		FusedRingNumberer.numberFusedRing(parentRing);//numbers the fused ring;
		Fragment fusedRing =state.fragManager.copyAndRelabel(parentRing);//makes sure the IDs are continuous
		state.fragManager.removeFragment(parentRing);
		return fusedRing;
	}


	/**
	 * Given a fragment returns it's atom list sorted by locant. e.g. 1,2,3,3a,3b,4
	 * @param frag
	 * @return
	 */
	private List<Atom> sortFragmentAtomListByLocant(Fragment frag) {
		List<Atom> atomsInFragment =frag.getAtomList();
		Collections.sort(atomsInFragment, new SortByLocants());
		return atomsInFragment;
	}

	/**
	 * Uses locants in front of the benz/benzo group to assign heteroatoms on the now numbered and used fused ring system
	 * @param state
	 * @param benzoEl
	 * @param fusedRing
	 * @throws StructureBuildingException
	 * @throws PostProcessingException
	 */
	private void benzoSpecificAssignHeteroAtomsUsingLocants(BuildState state, Element benzoEl, Fragment fusedRing) throws StructureBuildingException, PostProcessingException {
		Element previous = (Element) XOMTools.getPreviousSibling(benzoEl);
		if (previous!=null && previous.getLocalName().equals("multiplier")){//e.g. dibenzothiophene
			throw new StructureBuildingException("multiple benzo groups cannot currently be fused");
		}
		LinkedList<Element> locants =new LinkedList<Element>();
		while(previous != null && previous.getLocalName().equals("locant")) {
			locants.add(previous);
			previous=(Element) XOMTools.getPreviousSibling(previous);
		}
		if (locants.size() >0){
			Elements suffixes=((Element)benzoEl.getParent()).getChildElements("suffix");
			int suffixesWithoutLocants =0;
			for (int i = 0; i < suffixes.size(); i++) {
				if (suffixes.get(i).getAttribute("locant")==null){
					suffixesWithoutLocants++;
				}
			}
			if (locants.size() != suffixesWithoutLocants){//In preference locants will be assigned to suffixes rather than to this nomenclature
				List<Atom> atomList =fusedRing.getAtomList();
				LinkedList<Atom> heteroatoms =new LinkedList<Atom>();
				LinkedList<String> elementOfHeteroAtom =new LinkedList<String>();
				for (Atom atom : atomList) {//this iterates in the same order as the numbering system
					if (!atom.getElement().equals("C")){
						heteroatoms.add(atom);
						elementOfHeteroAtom.add(atom.getElement());
					}
				}
				if (locants.size() >=heteroatoms.size()){//atleast as many locants as there are heteroatoms to assign
					for (Atom atom : heteroatoms) {
						atom.setElement("C");
					}
					fusedRing.pickUpIndicatedHydrogen();
					for (int i=0; i< heteroatoms.size(); i ++){
						String elementSymbol =elementOfHeteroAtom.removeLast();
						Element locant =locants.removeFirst();
						fusedRing.getAtomByLocantOrThrow(locant.getAttributeValue("value")).setElement(elementSymbol);
						locant.detach();
					}
				}
				else if (locants.size()>1){
					throw new PostProcessingException("Unable to assign all locants to benzo-fused ring or multiplier was mising");
				}
			}
		}
	}
}
