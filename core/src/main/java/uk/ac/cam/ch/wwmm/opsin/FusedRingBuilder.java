package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.regex.Pattern;

import nu.xom.Element;
import nu.xom.Elements;

/**
 *
 * @author dl387
 *
 */
class FusedRingBuilder {
	private Pattern matchColon = Pattern.compile(":");
	private Pattern matchComma = Pattern.compile(",");
	private Pattern matchDash = Pattern.compile("-");
	private Pattern matchSlash = Pattern.compile("/");
	private Pattern matchC = Pattern.compile("C");
	Pattern matchPrime = Pattern.compile("'");

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
		Element lastGroup = groups.get(groups.size()-1);
		/*
		 * Apply any nonstandard ring numbering and sort atomOrder by locant
		 */
		for(int i=0;i<groups.size();i++) {
			Element group=groups.get(i);
			Fragment ring =state.xmlFragmentMap.get(group);
			if (group == lastGroup){
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
				else{
					ring.sortAtomListByLocant();//for those where the order the locants are in is sensible					}
				}
				for (Atom atom : atomList) {
					atom.clearLocants();//the parentRing does not have locants, letters are used to indicate the edges
				}
			}
			else if (group.getAttribute("fusedRingNumbering")==null){
				ring.sortAtomListByLocant();//for those where the order the locants are in is sensible
			}
		}
		Element previous = (Element) XOMTools.getPreviousSibling(groups.get(0));//TODO support this
		if (previous!=null && previous.getLocalName().equals("multiplier")){//e.g. dibenzothiophene
			throw new StructureBuildingException("multiplied components cannot currently be fused");
		}
		for(int i= groups.size() -2;i >=0; i--) {
			if (groups.get(i).getValue().equals("benz") || groups.get(i).getValue().equals("benzo")){
				Element possibleFusionbracket = (Element) XOMTools.getNextSibling(groups.get(i));
				if (!possibleFusionbracket.getLocalName().equals("fusion")){
					//e.g. 2-benzofuran. Fused rings of this type are a special case treated as being a single component
					//and have a special convention for indicating the position of heteroatoms 
					benzoSpecificFusion(state, groups.get(i), groups.get(i+1));
					groups.get(i).detach();
					groups =subOrRoot.getChildElements("group");
				}
			}
		}
		List<Element> nameComponents  = XOMTools.getChildElementsWithTagNames(subOrRoot, new String[]{"fusion","group"});
		nameComponents.remove(lastGroup);
		
		/*
		 * The number of primes on the component to be connected. 
		 * This is initially 0 indicating fusion of unprimed locants with the letter locants of the parentRing
		 * Subsequently it will switch to 1 indicating fusion of a second order component (primed locants) with a 
		 * first order component (unprimed locants)
		 * Next would be double primed fusing to single primed locants etc.
		 * 
		 */
		Fragment parentRing = state.xmlFragmentMap.get(lastGroup);
		int fusionLevel = 0;
		List<Fragment> fragmentInScopeForEachFusionLevel = new ArrayList<Fragment>();
		fragmentInScopeForEachFusionLevel.add(0, parentRing);//theLast ring fused or parentRing
		for (int i = nameComponents.size()-1; i>=0; i--) {
			Element fusion = null;
			if (nameComponents.get(i).getLocalName().equals("fusion")){
				fusion = nameComponents.get(i--);
			}
			if (i <0 || !nameComponents.get(i).getLocalName().equals("group")){
				throw new StructureBuildingException("Group not found where group expected. This is probably a bug");
			}
			Fragment nextComponent = state.xmlFragmentMap.get(nameComponents.get(i));
			if (fusion==null || matchColon.split(fusion.getValue()).length==1){
				if (fusion!=null){//A fusion bracket without a colon is always used when applying to the parent component
					//check for case of ommitted locant from a higher order fusion bracket e.g. cyclopenta[4,5]pyrrolo[2,3-c]pyridine
					if (matchDash.split(fusion.getValue()).length==1 && 
							matchComma.split(fusion.getValue()).length >1 &&
							allAtomsAreIdentical(nextComponent)){
						List<String> numericalLocantsOfParent = Arrays.asList(matchComma.split(fusion.getValue().substring(1, fusion.getValue().length()-1)));
						List<String> numericalLocantsOfChild = findPossibleNumericalLocants(nextComponent, numericalLocantsOfParent.size()-1);
						processHigherOrderFusionDescriptors(state, nextComponent, fragmentInScopeForEachFusionLevel.get(fusionLevel), numericalLocantsOfChild, numericalLocantsOfParent);
					}
					else{
						fusionLevel = 0;
						performSimpleFusion(state, fusion, nextComponent, fragmentInScopeForEachFusionLevel.get(fusionLevel));//e.g. pyrano[3,2-b]imidazo[4,5-e]pyridine where both are level 0 fusions
					}
				}
				else{
					if (fusionLevel > 0){
						FragmentTools.relabelLocants(nextComponent.getAtomList(), StringTools.multiplyString("'", fusionLevel));
					}
					performSimpleFusion(state, fusion, nextComponent, fragmentInScopeForEachFusionLevel.get(fusionLevel));
				}
			}
			else{
				String firstLocant = matchComma.split(fusion.getValue())[0];
				int numberOfPrimes =0;//determine number of primes in fusor and hence determine fusion level
				for(int j = firstLocant.length() -1; j>0; j--){
					if (firstLocant.charAt(j)=='\''){
						numberOfPrimes++;
					}
				}
				if (numberOfPrimes != fusionLevel){
					if (fusionLevel == numberOfPrimes +1){
						fusionLevel = numberOfPrimes;
					}
					else{
						throw new StructureBuildingException("Incorrect number of primes in fusion bracket: " +fusion.getValue());
					}
				}
				FragmentTools.relabelLocants(nextComponent.getAtomList(), StringTools.multiplyString("'", fusionLevel));
				performHigherOrderFusion(state, fusion, nextComponent, fragmentInScopeForEachFusionLevel.get(fusionLevel));
			}
			fusionLevel++;
			if (fragmentInScopeForEachFusionLevel.size() <= fusionLevel){
				fragmentInScopeForEachFusionLevel.add(fusionLevel, nextComponent);
			}
			else{
				fragmentInScopeForEachFusionLevel.set(fusionLevel, nextComponent);
			}
		}

		for(int i=0; i< groups.size() - 1 ;i++) {
			Fragment ring =state.xmlFragmentMap.get(groups.get(i));
			state.fragManager.incorporateFragment(ring, parentRing);
		}
		FusedRingNumberer.numberFusedRing(parentRing);//numbers the fused ring;
		state.fragManager.removeFragment(parentRing);
		Fragment fusedRing =state.fragManager.copyAndRelabel(parentRing);//makes sure the IDs are continuous

		StringBuilder fusedRingName = new StringBuilder();
		for (Element element : nameComponents) {
			fusedRingName.append(element.getValue());
		}
		fusedRingName.append(groups.get(groups.size()-1).getValue());

		Element fusedRingEl =groups.get(groups.size()-1);//reuse this element to save having to remap suffixes...
		fusedRingEl.getAttribute("value").setValue(fusedRingName.toString());
		fusedRingEl.getAttribute("valType").setValue("generatedFragment");
		fusedRingEl.getAttribute("type").setValue("ring");
		fusedRingEl.getAttribute("subType").setValue("fusedRing");
		fusedRingEl.removeChildren();
		fusedRingEl.appendChild(fusedRingName.toString());

		state.xmlFragmentMap.put(fusedRingEl, fusedRing);

		for (Element element : nameComponents) {
			element.detach();
		}
		return fusedRingEl;
	}

	private void performSimpleFusion(BuildState state, Element fusion, Fragment childRing, Fragment parentRing) throws StructureBuildingException {
		List<String> numericalLocantsOfChild = null;
		List<String> letterLocantsOfParent = null;
		if (fusion != null){
			String[] fusionArray = matchDash.split(fusion.getValue().substring(1, fusion.getValue().length()-1));
			if (fusionArray.length ==2){
				numericalLocantsOfChild = Arrays.asList(matchComma.split(fusionArray[0]));
				char[] tempLetterLocantsOfParent = fusionArray[1].toCharArray();
				letterLocantsOfParent = new ArrayList<String>();
				for (int i = 0; i < tempLetterLocantsOfParent.length; i++) {
					letterLocantsOfParent.add(String.valueOf(tempLetterLocantsOfParent[i]));
				}
			}
			else{
				if (fusionArray[0].contains(",")){//only has digits
					String[] numericalLocantsOfChildTemp = matchComma.split(fusionArray[0]);
					numericalLocantsOfChild = Arrays.asList(numericalLocantsOfChildTemp);
				}
				else{//only has letters
					char[] tempLetterLocantsOfParentCharArray = fusionArray[0].toCharArray();
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
			numericalLocantsOfChild = findPossibleNumericalLocants(childRing, edgeLength);
		}

		if (letterLocantsOfParent == null){
			letterLocantsOfParent = findPossibleLetterLocants(parentRing, edgeLength);
		}
		if (numericalLocantsOfChild == null || letterLocantsOfParent ==null){
			throw new StructureBuildingException("Unable to find bond to form fused ring system. Some information for forming fused ring system was only supplyed implicitly");
		}

		processFirstOrderFusionDescriptors(state, childRing, parentRing, numericalLocantsOfChild, letterLocantsOfParent);//fuse the rings
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
		List<Atom> atomlist = ring.getAtomList();
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
		List<Atom> atomlist = ring.getAtomList();
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
	 * @throws StructureBuildingException
	 */
	private void processFirstOrderFusionDescriptors(BuildState state, Fragment childRing, Fragment parentRing, List<String> numericalLocantsOfChild, List<String> letterLocantsOfParent) throws StructureBuildingException {
		List<Atom> childAtomList = childRing.getAtomList();
		int indexfirst = childAtomList.indexOf(childRing.getAtomByLocantOrThrow(numericalLocantsOfChild.get(0)));
		int indexfinal = childAtomList.indexOf(childRing.getAtomByLocantOrThrow(numericalLocantsOfChild.get(numericalLocantsOfChild.size()-1)));
		CyclicAtomList cyclicListAtomsInChild = new CyclicAtomList(childAtomList, indexfirst);
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
			while (cyclicListAtomsInChild.getIndice() != indexfinal){//assume numbers are descending
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

		List<Atom> parentAtoms = new ArrayList<Atom>();
		List<Atom> parentAtomList = parentRing.getAtomList();
		CyclicAtomList cyclicListAtomsInParent = new CyclicAtomList(parentAtomList, (int)letterLocantsOfParent.get(0).charAt(0) -97);//convert from lower case character through ascii to 0-23
		parentAtoms.add(cyclicListAtomsInParent.getCurrent());
		for (int i = 0; i < letterLocantsOfParent.size(); i++) {
			parentAtoms.add(cyclicListAtomsInParent.getNext());
		}
		fuseRings(state, childRing, parentRing, childAtoms, parentAtoms);
	}
	
	/**
	 * pyrido[1'',2'':1',2']imidazo[4',5':5,6]pyrazino[2,3-b]phenazine
	 * @param state
	 * @param higherOrderFusion
	 * @param nextComponent
	 * @param fusedRing
	 * @throws StructureBuildingException 
	 */
	private void performHigherOrderFusion(BuildState state, Element higherOrderFusion, Fragment nextComponent, Fragment fusedRing) throws StructureBuildingException {
		List<String> numericalLocantsOfChild = null;
		List<String> numericalLocantsOfParent = null;
		String[] fusionArray = matchColon.split(higherOrderFusion.getValue().substring(1, higherOrderFusion.getValue().length()-1));
		if (fusionArray.length ==2){
			numericalLocantsOfChild = Arrays.asList(matchComma.split(fusionArray[0]));
			numericalLocantsOfParent = Arrays.asList(matchComma.split(fusionArray[1]));
		}
		else{
			throw new StructureBuildingException("Malformed fusion bracket: This is an OPSIN bug, check regexTokens.xml");
		}
		processHigherOrderFusionDescriptors(state, nextComponent, fusedRing, numericalLocantsOfChild, numericalLocantsOfParent);//fuse the rings
	}

	/**
	 * Performs a single ring fusion using the values in numericalLocantsOfChild/numericalLocantsOfParent
	 * @param state
	 * @param childRing
	 * @param parentRing
	 * @param numericalLocantsOfChild
	 * @param letterLocantsOfParent
	 * @throws StructureBuildingException
	 */
	private void processHigherOrderFusionDescriptors(BuildState state, Fragment childRing, Fragment parentRing, List<String> numericalLocantsOfChild, List<String> numericalLocantsOfParent) throws StructureBuildingException {
		List<Atom> childAtomList = childRing.getAtomList();
		int indexfirst = childAtomList.indexOf(childRing.getAtomByLocantOrThrow(numericalLocantsOfChild.get(0)));
		int indexfinal = childAtomList.indexOf(childRing.getAtomByLocantOrThrow(numericalLocantsOfChild.get(numericalLocantsOfChild.size()-1)));
		CyclicAtomList cyclicListAtomsInChild = new CyclicAtomList(childAtomList, indexfirst);
		List<Atom> childAtoms = null;
		
		//need to determine whether locants are ascending or descending...locants may be omitted!
		List<Atom> potentialChildAtomsAscending = new ArrayList<Atom>();
		potentialChildAtomsAscending.add(cyclicListAtomsInChild.getCurrent());
		while (cyclicListAtomsInChild.getIndice() != indexfinal){//assume numbers are ascending
			potentialChildAtomsAscending.add(cyclicListAtomsInChild.getNext());
		}
		boolean notInPotentialChildAtoms =false;
		for (int i =1; i < numericalLocantsOfChild.size()-1 ; i ++){
			if (!potentialChildAtomsAscending.contains(childRing.getAtomByLocantOrThrow(numericalLocantsOfChild.get(i)))){
				notInPotentialChildAtoms =true;
			}
		}
		if (notInPotentialChildAtoms){
			potentialChildAtomsAscending  = null;
		}
		
		//and again for descending
		cyclicListAtomsInChild.setIndice(indexfirst);
		List<Atom> potentialChildAtomsDescending = new ArrayList<Atom>();
		potentialChildAtomsDescending.add(cyclicListAtomsInChild.getCurrent());
		while (cyclicListAtomsInChild.getIndice() != indexfinal){//assume numbers are descending
			potentialChildAtomsDescending.add(cyclicListAtomsInChild.getPrevious());
		}
		notInPotentialChildAtoms =false;
		for (int i =1; i < numericalLocantsOfChild.size()-1 ; i ++){
			if (!potentialChildAtomsDescending.contains(childRing.getAtomByLocantOrThrow(numericalLocantsOfChild.get(i)))){
				notInPotentialChildAtoms =true;
			}
		}
		if (notInPotentialChildAtoms){
			potentialChildAtomsDescending = null;
		}
		if (potentialChildAtomsAscending ==null){
			childAtoms = potentialChildAtomsDescending;
		}
		else if (potentialChildAtomsDescending ==null){
			childAtoms = potentialChildAtomsAscending;
		}
		else{
			if (potentialChildAtomsDescending.size() < potentialChildAtomsAscending.size() ){
				childAtoms = potentialChildAtomsDescending;
			}
			else{
				childAtoms = potentialChildAtomsAscending;
			}
		}
		
		if (childAtoms ==null){
			throw new StructureBuildingException("Malformed fusion bracket!");
		}

		List<Atom> parentAtomList = parentRing.getAtomList();
		indexfirst = parentAtomList.indexOf(parentRing.getAtomByLocantOrThrow(numericalLocantsOfParent.get(0)));
		indexfinal = parentAtomList.indexOf(parentRing.getAtomByLocantOrThrow(numericalLocantsOfParent.get(numericalLocantsOfParent.size()-1)));
		CyclicAtomList cyclicListAtomsInParent = new CyclicAtomList(parentAtomList, indexfirst);
		List<Atom> parentAtoms = null;
		
		List<Atom> potentialParentAtoms = new ArrayList<Atom>();
		potentialParentAtoms.add(cyclicListAtomsInParent.getCurrent());
		while (cyclicListAtomsInParent.getIndice() != indexfinal){//assume numbers are ascending
			potentialParentAtoms.add(cyclicListAtomsInParent.getNext());
		}
		if (childAtoms.size() == potentialParentAtoms.size()){
			boolean notInPotentialParentAtoms =false;
			for (int i =1; i < numericalLocantsOfParent.size()-1 ; i ++){
				if (!potentialParentAtoms.contains(parentRing.getAtomByLocantOrThrow(numericalLocantsOfParent.get(i)))){
					notInPotentialParentAtoms =true;
				}
			}
			if (!notInPotentialParentAtoms){
				parentAtoms = potentialParentAtoms;
			}
		}
		
		if (parentAtoms ==null){//that didn't work, so try assuming the numbers are descending
			cyclicListAtomsInParent.setIndice(indexfirst);
			potentialParentAtoms.clear();
			potentialParentAtoms.add(cyclicListAtomsInParent.getCurrent());
			while (cyclicListAtomsInParent.getIndice() != indexfinal){//assume numbers are ascending
				potentialParentAtoms.add(cyclicListAtomsInParent.getPrevious());
			}
			if (childAtoms.size() == potentialParentAtoms.size()){
				boolean notInPotentialParentAtoms =false;
				for (int i =1; i < numericalLocantsOfParent.size()-1 ; i ++){
					if (!potentialParentAtoms.contains(parentRing.getAtomByLocantOrThrow(numericalLocantsOfParent.get(i)))){
						notInPotentialParentAtoms =true;
					}
				}
				if (!notInPotentialParentAtoms){
					parentAtoms = potentialParentAtoms;
				}
			}
		}
		if (parentAtoms ==null){
			throw new StructureBuildingException("Malformed fusion bracket!");
		}
		fuseRings(state, childRing, parentRing, childAtoms, parentAtoms);
	}

	/**
	 * Fuses two rings together returning a fragment containing the fusedRing
	 * @param state
	 * @param childRing
	 * @param parentRing
	 * @param childAtoms
	 * @param parentAtoms
	 * @throws StructureBuildingException
	 */
	private void fuseRings(BuildState state, Fragment childRing, Fragment parentRing, List<Atom> childAtoms, List<Atom> parentAtoms) throws StructureBuildingException {
		if (parentAtoms.size()!=childAtoms.size()){
			throw new StructureBuildingException("Problem with fusion descriptors: Parent atoms specified: " + parentAtoms.size() +" Child atoms specified: " + childAtoms.size() + " These should have been identical!");
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

		for (int i = 0; i < parentAtoms.size(); i++) {
			Atom parentAtom = parentAtoms.get(i);
			if (!parentAtom.getElement().equals(childAtoms.get(i).getElement())){
				throw new StructureBuildingException("Invalid fusion descriptor: Heteroatom placement is ambigous as it is not present in both components of the fusion");
			}
			for (Atom atom : neighboursOfToBeReplacedChildAtoms.get(i)) {
				//System.out.println("Atom ID " + atom.getID() +" bonded to " +  parentAtom.getID());
				state.fragManager.createBond(parentAtom, atom, 1);
			}
		}
	}

	/**
	 * Fuse the benzo with the subsequent ring
	 * Uses locants in front of the benz/benzo group to assign heteroatoms on the now numbered fused ring system
	 * @param state
	 * @param benzoEl
	 * @param parentEl
	 * @throws StructureBuildingException
	 * @throws PostProcessingException
	 */
	private void benzoSpecificFusion(BuildState state, Element benzoEl, Element parentEl) throws StructureBuildingException, PostProcessingException {
		
		/*
		 * Perform the fusion, number it and associate it with the parentEl
		 */
		Fragment benzoRing = state.xmlFragmentMap.get(benzoEl);
		Fragment parentRing = state.xmlFragmentMap.get(parentEl);
		performSimpleFusion(state, null, benzoRing , parentRing);
		state.fragManager.incorporateFragment(benzoRing, parentRing);
		FusedRingNumberer.numberFusedRing(parentRing);//numbers the fused ring;
		state.fragManager.removeFragment(parentRing);
		Fragment fusedRing =state.fragManager.copyAndRelabel(parentRing);//makes sure the IDs are continuous
		state.xmlFragmentMap.put(parentEl, fusedRing);

		/*
		 * Check for locants and use these to set the heteroatom positions
		 */
		Element previous = (Element) XOMTools.getPreviousSibling(benzoEl);
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
	
	/**
	 * Checks that all atoms in a ring appear to be equivalent
	 * @param ring
	 * @return true if all equivalent, else false
	 */
	private boolean  allAtomsAreIdentical(Fragment ring){
		List<Atom> atomList = ring.getAtomList();
		Atom firstAtom = atomList.get(0);
		String element = firstAtom.getElement();
		int valency = firstAtom.getIncomingValency();
		boolean spareValency = firstAtom.hasSpareValency();
		for (Atom atom : atomList) {
			if (!atom.getElement().equals(element)){
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
}
