package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static uk.ac.cam.ch.wwmm.opsin.OpsinTools.*;
import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

import nu.xom.Element;
import nu.xom.Elements;

/**
 * Assembles fused rings named using fusion nomenclature
 * @author dl387
 *
 */
class FusedRingBuilder {
	private final BuildState state;
	private final List<Element> groupsInFusedRing;
	private final Element lastGroup;
	private final Fragment parentRing;
	private final Map<Integer,Fragment> fragmentInScopeForEachFusionLevel = new HashMap<Integer,Fragment>();
	private Map<Atom, Atom> atomsToRemoveToReplacementAtom = new HashMap<Atom, Atom>();

	private FusedRingBuilder(BuildState state, List<Element> groupsInFusedRing) {
		this.state = state;
		this.groupsInFusedRing = groupsInFusedRing;
		lastGroup = groupsInFusedRing.get(groupsInFusedRing.size()-1);
		parentRing = state.xmlFragmentMap.get(lastGroup);
		fragmentInScopeForEachFusionLevel.put(0, parentRing);
	}

	/**
	 * Master method for processing fused rings. Fuses groups together
	 * @param state: contains the current id and fragment manager
	 * @param subOrRoot Element (substituent or root)
	 * @throws StructureBuildingException
	 */
	static void processFusedRings(BuildState state, Element subOrRoot) throws  StructureBuildingException {
		List<Element> groups = XOMTools.getChildElementsWithTagName(subOrRoot, GROUP_EL);
		if (groups.size() < 2){
			return;//nothing to fuse
		}
		List<Element> groupsInFusedRing =new ArrayList<Element>();
		for (int i = groups.size()-1; i >=0; i--) {//group groups into fused rings
			Element group =groups.get(i);
			groupsInFusedRing.add(0, group);
			if (i!=0){
				Element startingEl = group;
				if ((group.getValue().equals("benz") || group.getValue().equals("benzo")) && FUSIONRING_SUBTYPE_VAL.equals(group.getAttributeValue(SUBTYPE_ATR))){
					Element beforeBenzo = (Element) XOMTools.getPreviousSibling(group);
					if (beforeBenzo !=null && beforeBenzo.getLocalName().equals(LOCANT_EL)){
						startingEl = beforeBenzo;
					}
				}
				Element possibleGroup = XOMTools.getPreviousSiblingIgnoringCertainElements(startingEl, new String[]{MULTIPLIER_EL, FUSION_EL});
				if (!groups.get(i-1).equals(possibleGroup)){//end of fused ring system
					if (groupsInFusedRing.size()>=2){
						//This will be invoked in cases where there are multiple fused ring systems in the same subOrRoot such as some spiro systems
						new FusedRingBuilder(state, groupsInFusedRing).buildFusedRing();
					}
					groupsInFusedRing.clear();
				}
			}
		}
		if (groupsInFusedRing.size()>=2){
			new FusedRingBuilder(state, groupsInFusedRing).buildFusedRing();
		}
	}

	/**
	 * Combines the groups given in the {@link FusedRingBuilder} constructor to destructively create the fused ring system
	 * This fused ring is then numbered
	 * @throws StructureBuildingException
	 */
	void buildFusedRing() throws StructureBuildingException{
		/*
		 * Apply any nonstandard ring numbering, sorts atomOrder by locant
		 * Aromatises appropriate cycloalkane rings, Rejects groups with acyclic atoms
		 */
        processRingNumberingAndIrregularities();
		processBenzoFusions();//FR-2.2.8  e.g. in 2H-[1,3]benzodioxino[6',5',4':10,5,6]anthra[2,3-b]azepine  benzodioxino is one component
		List<Element> nameComponents = formNameComponentList();
		nameComponents.remove(lastGroup);

		List<Fragment> componentFragments = new ArrayList<Fragment>();//all the ring fragments (other than the parentRing). These will later be merged into the parentRing
		List<Fragment> parentFragments = new ArrayList<Fragment>();
		parentFragments.add(parentRing);
		
		int numberOfParents = 1;
		Element possibleMultiplier = (Element) XOMTools.getPreviousSibling(lastGroup);
		if (nameComponents.size()>0 && possibleMultiplier !=null && possibleMultiplier.getLocalName().equals(MULTIPLIER_EL)){
			numberOfParents = Integer.parseInt(possibleMultiplier.getAttributeValue(VALUE_ATR));
			possibleMultiplier.detach();
			for (int j = 1; j < numberOfParents; j++) {
				Fragment copyOfParentRing =state.fragManager.copyFragment(parentRing);
				parentFragments.add(copyOfParentRing);
				componentFragments.add(copyOfParentRing);
			}
		}

		/*The indice from nameComponents to use next. Work from right to left i.e. starts at nameComponents.size()-1*/
		int ncIndice = processMultiParentSystem(parentFragments, nameComponents, componentFragments);//handle multiparent systems
		/*
		 * The number of primes on the component to be connected. 
		 * This is initially 0 indicating fusion of unprimed locants with the letter locants of the parentRing
		 * Subsequently it will switch to 1 indicating fusion of a second order component (primed locants) with a 
		 * first order component (unprimed locants)
		 * Next would be double primed fusing to single primed locants etc.
		 * 
		 */
		int fusionLevel = (nameComponents.size()-1 -ncIndice)/2;
		for (; ncIndice>=0; ncIndice--) {
			Element fusion = null;
			if (nameComponents.get(ncIndice).getLocalName().equals(FUSION_EL)){
				fusion = nameComponents.get(ncIndice--);
			}
			if (ncIndice <0 || !nameComponents.get(ncIndice).getLocalName().equals(GROUP_EL)){
				throw new StructureBuildingException("Group not found where group expected. This is probably a bug");
			}
			Fragment nextComponent = state.xmlFragmentMap.get(nameComponents.get(ncIndice));
			int multiplier = 1;
			Element possibleMultiplierEl = (Element) XOMTools.getPreviousSibling(nameComponents.get(ncIndice));//e.g. the di of difuro
			if (possibleMultiplierEl != null && possibleMultiplierEl.getLocalName().equals(MULTIPLIER_EL)){
				multiplier = Integer.parseInt(possibleMultiplierEl.getAttributeValue(VALUE_ATR));
			}
			String[] fusionDescriptors =null;
			if (fusion !=null){
				String fusionDescriptorString = fusion.getValue().toLowerCase().substring(1, fusion.getValue().length()-1);
				if (multiplier ==1){
					fusionDescriptors = new String[]{fusionDescriptorString};
				}
				else{
					if (MATCH_SEMICOLON.split(fusionDescriptorString).length >1){
						fusionDescriptors = MATCH_SEMICOLON.split(fusionDescriptorString);
					}
					else if (MATCH_COLON.split(fusionDescriptorString).length >1){
						fusionDescriptors = MATCH_COLON.split(fusionDescriptorString);
					}
					else if (MATCH_COMMA.split(fusionDescriptorString).length >1){
						fusionDescriptors = MATCH_COMMA.split(fusionDescriptorString);
					}
					else{//multiplier does not appear to mean multiplied component. Could be indicating multiplication of the whole fused ring system
						if (ncIndice!=0){
							throw new StructureBuildingException("Unexpected multiplier: " + possibleMultiplierEl.getValue() +" or incorrect fusion descriptor: " + fusionDescriptorString);
						}
						multiplier =1;
						fusionDescriptors = new String[]{fusionDescriptorString};
					}
				}
			}
			if (multiplier >1){
				possibleMultiplierEl.detach();
			}
			Fragment[] fusionComponents = new Fragment[multiplier];
			for (int j = 0; j < multiplier; j++) {
				if (j>0){
					fusionComponents[j] = state.fragManager.copyAndRelabelFragment(nextComponent,  j);
				}
				else{
					fusionComponents[j] = nextComponent;
				}
			}
			
			for (int j = 0; j < multiplier; j++) {
				Fragment component = fusionComponents[j];
				componentFragments.add(component);
				if (fusion !=null){
					if (MATCH_COLON.split(fusionDescriptors[j]).length==1){//A fusion bracket without a colon is used when applying to the parent component (except in a special case where locants are ommitted)
						//check for case of omitted locant from a higher order fusion bracket e.g. cyclopenta[4,5]pyrrolo[2,3-c]pyridine
						if (MATCH_DASH.split(fusionDescriptors[j]).length==1 && 
								MATCH_COMMA.split(fusionDescriptors[j]).length >1 &&
								FragmentTools.allAtomsInRingAreIdentical(component)){
							relabelAccordingToFusionLevel(component, fusionLevel);
							List<String> numericalLocantsOfParent = Arrays.asList(MATCH_COMMA.split(fusionDescriptors[j]));
							List<String> numericalLocantsOfChild = findPossibleNumericalLocants(component, determineAtomsToFuse(fragmentInScopeForEachFusionLevel.get(fusionLevel), numericalLocantsOfParent, null).size()-1);
							processHigherOrderFusionDescriptors(component, fragmentInScopeForEachFusionLevel.get(fusionLevel), numericalLocantsOfChild, numericalLocantsOfParent);
						}
						else{
							fusionLevel = 0;
							relabelAccordingToFusionLevel(component, fusionLevel);
							String fusionDescriptor = fusionDescriptors[j];
							String[] fusionArray = determineNumericalAndLetterComponents(fusionDescriptor);
							int numberOfPrimes =0;
							if (!fusionArray[1].equals("")){
								numberOfPrimes =StringTools.countTerminalPrimes(fusionArray[1]);
								if (fusionArray[0].equals("")){
									fusionDescriptor = fusionArray[1].replaceAll("'", "");
								}
								else{
									fusionDescriptor = fusionArray[0]+ "-" +fusionArray[1].replaceAll("'", "");
								}
								if (numberOfPrimes >= parentFragments.size()){
									throw new StructureBuildingException("Unexpected prime in fusion descriptor");
								}
							}
							performSimpleFusion(fusionDescriptor, component, parentFragments.get(numberOfPrimes));//e.g. pyrano[3,2-b]imidazo[4,5-e]pyridine where both are level 0 fusions
						}
					}
					else{
						String firstLocant = MATCH_COMMA.split(fusionDescriptors[j])[0];
						int numberOfPrimes = -j;//determine number of primes in fusor and hence determine fusion level
						for(int k = firstLocant.length() -1; k>0; k--){
							if (firstLocant.charAt(k)=='\''){
								numberOfPrimes++;
							}
						}
						if (numberOfPrimes != fusionLevel){
							if (fusionLevel == numberOfPrimes +1){
								fusionLevel = numberOfPrimes;
							}
							else{
								throw new StructureBuildingException("Incorrect number of primes in fusion bracket: " +fusionDescriptors[j]);
							}
						}
						relabelAccordingToFusionLevel(component, fusionLevel);
						performHigherOrderFusion(fusionDescriptors[j], component, fragmentInScopeForEachFusionLevel.get(fusionLevel));
					}
				}
				else{
					relabelAccordingToFusionLevel(component, fusionLevel);
					performSimpleFusion(null, component, fragmentInScopeForEachFusionLevel.get(fusionLevel));
				}
			}
			fusionLevel++;
			if (multiplier ==1){//multiplied components may not be substituted onto
				fragmentInScopeForEachFusionLevel.put(fusionLevel, fusionComponents[0]);
			}
		}
		for (Fragment ring : componentFragments) {
			state.fragManager.incorporateFragment(ring, parentRing);
		}
		removeMergedAtoms();

		FusedRingNumberer.numberFusedRing(parentRing);//numbers the fused ring;

		StringBuilder fusedRingName = new StringBuilder();
		for (Element element : nameComponents) {
			fusedRingName.append(element.getValue());
		}
		fusedRingName.append(lastGroup.getValue());

		Element fusedRingEl =lastGroup;//reuse this element to save having to remap suffixes...
		fusedRingEl.getAttribute(VALUE_ATR).setValue(fusedRingName.toString());
		fusedRingEl.removeAttribute(fusedRingEl.getAttribute(VALTYPE_ATR));
		fusedRingEl.getAttribute(TYPE_ATR).setValue(RING_TYPE_VAL);
		fusedRingEl.getAttribute(SUBTYPE_ATR).setValue(FUSEDRING_SUBTYPE_VAL);
		XOMTools.setTextChild(fusedRingEl, fusedRingName.toString());

		for (Element element : nameComponents) {
			element.detach();
		}
	}

	private void removeMergedAtoms() {
		for (Atom a : atomsToRemoveToReplacementAtom.keySet()) {
			state.fragManager.removeAtomAndAssociatedBonds(a);
		}
		atomsToRemoveToReplacementAtom.clear();
	}

	/**
	 * Forms a list a list of all group and fusion elements between the first and last group in the fused ring
	 * @return
	 */
	private List<Element> formNameComponentList() {
		List<Element> nameComponents  = new ArrayList<Element>();
		Element currentEl = groupsInFusedRing.get(0);
		while(currentEl != lastGroup){
			if (currentEl.getLocalName().equals(GROUP_EL) || currentEl.getLocalName().equals(FUSION_EL)){
				nameComponents.add(currentEl);
			}
			currentEl = (Element) XOMTools.getNextSibling(currentEl);
		}
		return nameComponents;
	}

	private void processRingNumberingAndIrregularities() throws StructureBuildingException {
		for (Element group : groupsInFusedRing) {
            Fragment ring = state.xmlFragmentMap.get(group);
            if (ALKANESTEM_SUBTYPE_VAL.equals(group.getAttributeValue(SUBTYPE_ATR))){
            	aromatiseCyclicAlkane(group);
            }
            if (group == lastGroup) {
                //perform a quick check that every atom in this group is infact cyclic. Fusion components are enumerated and hence all guaranteed to be purely cyclic
                List<Atom> atomList = ring.getAtomList();
                for (Atom atom : atomList) {
                    if (!atom.getAtomIsInACycle()) {
                        throw new StructureBuildingException("Inappropriate group used in fusion nomenclature. Only groups composed entirely of atoms in cycles may be used. i.e. not: " + group.getValue());
                    }
                }
                if (group.getAttribute(FUSEDRINGNUMBERING_ATR) != null) {
                    String[] standardNumbering = MATCH_SLASH.split(group.getAttributeValue(FUSEDRINGNUMBERING_ATR), -1);
                    for (int j = 0; j < standardNumbering.length; j++) {
                        atomList.get(j).replaceLocants(standardNumbering[j]);
                    }
                } else {
                    ring.sortAtomListByLocant();//for those where the order the locants are in is sensible					}
                }
                for (Atom atom : atomList) {
                    atom.clearLocants();//the parentRing does not have locants, letters are used to indicate the edges
                }
            } else if (group.getAttribute(FUSEDRINGNUMBERING_ATR) == null) {
                ring.sortAtomListByLocant();//for those where the order the locants are in is sensible
            }
        }
	}

	/**
	 * Given a cyclicAlkaneGroup determines whether or not it should be aromatised. Unlocanted ene will be detached if it is an aromatisation hint
	 * No unsaturators -->aromatise
	 * Just ane -->don't
	 * More than 1 ene or locants on ene -->don't
	 * yne --> don't
	 * @param cyclicAlkaneGroup
	 */
	private void aromatiseCyclicAlkane(Element cyclicAlkaneGroup) {
		Element next = (Element) XOMTools.getNextSibling(cyclicAlkaneGroup);
		List<Element> unsaturators = new ArrayList<Element>();
		while (next!=null && next.getLocalName().equals(UNSATURATOR_EL)){
			unsaturators.add(next);
			next = (Element) XOMTools.getNextSibling(next);
		}
		boolean conjugate =true;
		if (unsaturators.size()==1){
			int value = Integer.parseInt(unsaturators.get(0).getAttributeValue(VALUE_ATR));
			if (value !=2){
				conjugate =false;
			}
			else if (unsaturators.get(0).getAttribute(LOCANT_ATR)!=null){
				conjugate =false;
			}
		}
		else if (unsaturators.size()==2){
			int value1 = Integer.parseInt(unsaturators.get(0).getAttributeValue(VALUE_ATR));
			if (value1 !=1){
				conjugate =false;
			}
			else{
				int value2 = Integer.parseInt(unsaturators.get(1).getAttributeValue(VALUE_ATR));
				if (value2 !=2 || unsaturators.get(1).getAttribute(LOCANT_ATR)!=null){
					conjugate =false;
				}
			}
		}
		else if (unsaturators.size() >2){
			conjugate =false;
		}
		if (conjugate){
			for (Element unsaturator : unsaturators) {
				unsaturator.detach();
			}
			List<Atom> atomList =state.xmlFragmentMap.get(cyclicAlkaneGroup).getAtomList();
		    for (Atom atom : atomList) {
		        atom.setSpareValency(true);
		    }
		}
	}

	private int processMultiParentSystem(List<Fragment> parentFragments, List<Element> nameComponents, List<Fragment> componentFragments) throws StructureBuildingException {
		int i = nameComponents.size()-1;
		int fusionLevel =0;
		if (i>=0 && parentFragments.size()>1){
			List<Fragment> previousFusionLevelFragments = parentFragments;
			for (; i>=0; i--) {
				if (previousFusionLevelFragments.size()==1){//completed multi parent system
					fragmentInScopeForEachFusionLevel.put(fusionLevel, previousFusionLevelFragments.get(0));
					break;
				}
				Element fusion = null;
				if (nameComponents.get(i).getLocalName().equals(FUSION_EL)){
					fusion = nameComponents.get(i--);
				}
				else{
					throw new StructureBuildingException("Fusion bracket not found where fusion bracket expected");
				}
				if (i <0 || !nameComponents.get(i).getLocalName().equals(GROUP_EL)){
					throw new StructureBuildingException("Group not found where group expected. This is probably a bug");
				}
				Fragment nextComponent = state.xmlFragmentMap.get(nameComponents.get(i));
				relabelAccordingToFusionLevel(nextComponent, fusionLevel);
				int multiplier = 1;
				Element possibleMultiplierEl = (Element) XOMTools.getPreviousSibling(nameComponents.get(i));
				if (possibleMultiplierEl != null && possibleMultiplierEl.getLocalName().equals(MULTIPLIER_EL)){
					multiplier = Integer.parseInt(possibleMultiplierEl.getAttributeValue(VALUE_ATR));
					possibleMultiplierEl.detach();
				}
				List<Fragment> fusionComponents = new ArrayList<Fragment>();
				for (int j = 0; j < multiplier; j++) {
					if (j>0){
						fusionComponents.add(state.fragManager.copyAndRelabelFragment(nextComponent,  j));
					}
					else{
						fusionComponents.add(nextComponent);
					}
				}
				fusionLevel+=multiplier;
				if (multiplier>1 && multiplier != previousFusionLevelFragments.size()){
					throw new StructureBuildingException("Mismatch between number of components and number of parents in fused ring system");
				}
				String fusionDescriptorString = fusion.getValue().toLowerCase().substring(1, fusion.getValue().length()-1);
				String[] fusionDescriptors =null;
				if (MATCH_SEMICOLON.split(fusionDescriptorString).length >1){
					fusionDescriptors = MATCH_SEMICOLON.split(fusionDescriptorString);
				}
				else if (MATCH_COLON.split(fusionDescriptorString).length >1){
					fusionDescriptors = MATCH_COLON.split(fusionDescriptorString);
				}
				else if (MATCH_COMMA.split(fusionDescriptorString).length >1){
					fusionDescriptors = MATCH_COMMA.split(fusionDescriptorString);
				}
				else{
					throw new StructureBuildingException("Invalid fusion descriptor: " + fusionDescriptorString);
				}
				if (fusionDescriptors.length != previousFusionLevelFragments.size()){
					throw new StructureBuildingException("Invalid fusion descriptor: "+fusionDescriptorString +"(Number of locants disagrees with number of parents)");
				}
				for (int j = 0; j < fusionDescriptors.length; j++) {
					String fusionDescriptor = fusionDescriptors[j];
					Fragment component = multiplier>1 ? fusionComponents.get(j) : nextComponent;
					Fragment parentToUse = previousFusionLevelFragments.get(j);
					boolean simpleFusion = MATCH_COLON.split(fusionDescriptor).length <= 1;
					if (simpleFusion){
						String[] fusionArray = determineNumericalAndLetterComponents(fusionDescriptor);
						if (fusionArray[1].length() != 0){
							int numberOfPrimes =StringTools.countTerminalPrimes(fusionArray[1]);
							if (fusionArray[0].length() == 0){
								fusionDescriptor = fusionArray[1].replaceAll("'", "");
							}
							else{
								fusionDescriptor = fusionArray[0]+ "-" +fusionArray[1].replaceAll("'", "");
							}
							if (numberOfPrimes !=j){//check the number of primes on the letter part agree with the parent to use e.g.[4,5-bcd:1,2-c']difuran
								throw new StructureBuildingException("Incorrect number of primes in fusion descriptor: " + fusionDescriptor);
							}
						}
						performSimpleFusion(fusionDescriptor, component, parentToUse);
					}
					else{
						performHigherOrderFusion(fusionDescriptor, component, parentToUse);
					}
				}
				previousFusionLevelFragments = fusionComponents;
				componentFragments.addAll(fusionComponents);
			}
			if (previousFusionLevelFragments.size()!=1){
				throw new StructureBuildingException("Invalid fused ring system. Incomplete multiparent system");
			}
		}
		return i;
	}

	/**
	 * Splits a first order fusion component into it's numerical and letter parts
	 * Either one of these can be the blank string as they may have been omitted
	 * The first entry in the array is the numbers and the second the letters
	 * @param fusionDescriptor
	 * @return
	 */
	private static String[] determineNumericalAndLetterComponents(String fusionDescriptor) {
		String[] fusionArray = MATCH_DASH.split(fusionDescriptor);
		if (fusionArray.length ==2){
			return fusionArray;
		}
		else{
			String[] components = new String[2];
			if (fusionArray[0].contains(",")){//the digit section
				components[0]=fusionArray[0];
				components[1]="";
			}
			else{
				components[0]="";
				components[1]=fusionArray[0];
			}
			return components;
		}
	}

	/**
	 * Searches groups for benz(o) components and fuses them in accordance with
	 * FR-2.2.8 Heterobicyclic components with a benzene ring
	 * @throws StructureBuildingException
	 */
	private void processBenzoFusions() throws StructureBuildingException {
		for(int i= groupsInFusedRing.size() -2;i >=0; i--) {
			if (groupsInFusedRing.get(i).getValue().equals("benz") || groupsInFusedRing.get(i).getValue().equals("benzo")){
				Element possibleFusionbracket = (Element) XOMTools.getNextSibling(groupsInFusedRing.get(i));
				if (!possibleFusionbracket.getLocalName().equals(FUSION_EL)){
					Element possibleMultiplier = (Element) XOMTools.getPreviousSibling(groupsInFusedRing.get(i));
					if (possibleMultiplier==null || !possibleMultiplier.getLocalName().equals(MULTIPLIER_EL)|| possibleMultiplier.getAttributeValue(TYPE_ATR).equals(GROUP_TYPE_VAL)){
						//e.g. 2-benzofuran. Fused rings of this type are a special case treated as being a single component
						//and have a special convention for indicating the position of heteroatoms 
						benzoSpecificFusion(groupsInFusedRing.get(i), groupsInFusedRing.get(i+1));
						groupsInFusedRing.get(i).detach();
						groupsInFusedRing.remove(i);
					}
				}
			}
		}
	}

	/**
	 * Modifies nextComponent's locants according to the fusionLevel.
	 * @param component
	 * @param fusionLevel
	 */
	private void relabelAccordingToFusionLevel(Fragment component, int fusionLevel)  {
		if (fusionLevel > 0){
			FragmentTools.relabelNumericLocants(component.getAtomList(), StringTools.multiplyString("'", fusionLevel));
		}
	}

	/**
	 * Handles fusion between components where the fusion descriptor is of the form:
	 * comma separated locants dash letters
	 * e.g imidazo[4,5-d]pyridine
	 * The fusionDescriptor may be given as null or the letter/numerical part omitted.
	 * Sensible defaults will be found instead
	 * @param fusionDescriptor
	 * @param childRing
	 * @param parentRing
	 * @throws StructureBuildingException
	 */
	private void performSimpleFusion(String fusionDescriptor, Fragment childRing, Fragment parentRing) throws StructureBuildingException {
		List<String> numericalLocantsOfChild = null;
		List<String> letterLocantsOfParent = null;
		if (fusionDescriptor != null){
			String[] fusionArray = MATCH_DASH.split(fusionDescriptor);
			if (fusionArray.length ==2){
				numericalLocantsOfChild = Arrays.asList(MATCH_COMMA.split(fusionArray[0]));
				char[] tempLetterLocantsOfParent = fusionArray[1].toCharArray();
				letterLocantsOfParent = new ArrayList<String>();
                for (char letterLocantOfParent : tempLetterLocantsOfParent) {
                    letterLocantsOfParent.add(String.valueOf(letterLocantOfParent));
                }
			}
			else{
				if (fusionArray[0].contains(",")){//only has digits
					String[] numericalLocantsOfChildTemp = MATCH_COMMA.split(fusionArray[0]);
					numericalLocantsOfChild = Arrays.asList(numericalLocantsOfChildTemp);
				}
				else{//only has letters
					char[] tempLetterLocantsOfParentCharArray = fusionArray[0].toCharArray();
					letterLocantsOfParent = new ArrayList<String>();
                    for (char letterLocantOfParentCharArray : tempLetterLocantsOfParentCharArray) {
                        letterLocantsOfParent.add(String.valueOf(letterLocantOfParentCharArray));
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

		processFirstOrderFusionDescriptors(childRing, parentRing, numericalLocantsOfChild, letterLocantsOfParent);//fuse the rings
	}

	/**
	 * Takes a ring an returns and array with one letter corresponding to a side/s
	 * that contains two adjacent non bridgehead carbons
	 * The number of sides is specified by edgeLength
	 * @param ring
	 * @param edgeLength The number of bonds to be fused along
	 * @return
	 */
	private static List<String> findPossibleLetterLocants(Fragment ring, int edgeLength) {
		List<Atom> atomlist = ring.getAtomList();
		List<String> letterLocantsOfParent = null;
		List<Atom> carbonAtoms = new ArrayList<Atom>();
		atomlist.add(0, atomlist.get(atomlist.size()-1));//this atomList is a copy so we can safely do this
		for (int i =atomlist.size() -1; i >=0; i--) {//iterate backwards in list to use highest locanted edge in preference.
			//this retains what is currently locant 1 on the parent ring as locant 1 if the first two atoms found match
			Atom atom = atomlist.get(i);
			if (atom.getElement().equals("C")){
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
	private static List<String> findPossibleNumericalLocants(Fragment ring, int edgeLength) {
		List<Atom> atomlist = ring.getAtomList();
		List<String> numericalLocantsOfChild = null;
		List<String> carbonLocants = new ArrayList<String>();
		atomlist.add(atomlist.get(0));//this atomList is a copy so we can safely do this
		for (Atom atom : atomlist) {
			if (atom.getElement().equals("C")){
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
	 * @param childRing
	 * @param parentRing
	 * @param numericalLocantsOfChild
	 * @param letterLocantsOfParent
	 * @throws StructureBuildingException
	 */
	private void processFirstOrderFusionDescriptors(Fragment childRing, Fragment parentRing, List<String> numericalLocantsOfChild, List<String> letterLocantsOfParent) throws StructureBuildingException {
		List<Atom> childAtoms = determineAtomsToFuse(childRing, numericalLocantsOfChild, letterLocantsOfParent.size() +1);
		if (childAtoms ==null){
			throw new StructureBuildingException("Malformed fusion bracket!");
		}

		List<Atom> parentAtoms = new ArrayList<Atom>();
		List<Atom> parentAtomList = parentRing.getAtomList();
		//find the indice of the last atom on the surface of the ring. This obviously connects to the first atom. The objective is to exclude any interior atoms.
		List<Atom> neighbours = parentAtomList.get(0).getAtomNeighbours();
		int indice = Integer.MAX_VALUE;
		for (Atom atom : neighbours) {
			int indexOfAtom =parentAtomList.indexOf(atom);
			if (indexOfAtom ==1){//not the next atom
				continue;
			}
			else if (indexOfAtom ==-1){//not in parentRing
				continue;
			}
			if (parentAtomList.indexOf(atom)< indice){
				indice = indexOfAtom;
			}
		}
		CyclicAtomList cyclicListAtomsOnSurfaceOfParent = new CyclicAtomList(parentAtomList.subList(0, indice +1), (int)letterLocantsOfParent.get(0).charAt(0) -97);//convert from lower case character through ascii to 0-23
		parentAtoms.add(cyclicListAtomsOnSurfaceOfParent.getCurrent());
		for (int i = 0; i < letterLocantsOfParent.size(); i++) {
			parentAtoms.add(cyclicListAtomsOnSurfaceOfParent.getNext());
		}
		fuseRings(childAtoms, parentAtoms);
	}
	
	/**
	 * Handles fusion between components where the fusion descriptor is of the form:
	 * comma separated locants colon comma separated locants
	 * e.g pyrido[1'',2'':1',2']imidazo
	 * @param fusionDescriptor
	 * @param nextComponent
	 * @param fusedRing
	 * @throws StructureBuildingException 
	 */
	private void performHigherOrderFusion(String fusionDescriptor, Fragment nextComponent, Fragment fusedRing) throws StructureBuildingException {
		List<String> numericalLocantsOfChild = null;
		List<String> numericalLocantsOfParent = null;
		String[] fusionArray = MATCH_COLON.split(fusionDescriptor);
		if (fusionArray.length ==2){
			numericalLocantsOfChild = Arrays.asList(MATCH_COMMA.split(fusionArray[0]));
			numericalLocantsOfParent = Arrays.asList(MATCH_COMMA.split(fusionArray[1]));
		}
		else{
			throw new StructureBuildingException("Malformed fusion bracket: This is an OPSIN bug, check regexTokens.xml");
		}
		processHigherOrderFusionDescriptors(nextComponent, fusedRing, numericalLocantsOfChild, numericalLocantsOfParent);//fuse the rings
	}

	/**
	 * Performs a single ring fusion using the values in numericalLocantsOfChild/numericalLocantsOfParent
	 * @param childRing
	 * @param parentRing
	 * @param numericalLocantsOfChild
	 * @param numericalLocantsOfParent
	 * @throws StructureBuildingException
	 */
	private void processHigherOrderFusionDescriptors(Fragment childRing, Fragment parentRing, List<String> numericalLocantsOfChild, List<String> numericalLocantsOfParent) throws StructureBuildingException {
		List<Atom> childAtoms =determineAtomsToFuse(childRing, numericalLocantsOfChild, null);
		if (childAtoms ==null){
			throw new StructureBuildingException("Malformed fusion bracket!");
		}

		List<Atom> parentAtoms = determineAtomsToFuse(parentRing, numericalLocantsOfParent, childAtoms.size());
		if (parentAtoms ==null){
			throw new StructureBuildingException("Malformed fusion bracket!");
		}
		fuseRings(childAtoms, parentAtoms);
	}

	/**
	 * Determines which atoms on a ring should be used for fusion given a set of numerical locants.
	 * If from the other ring involved in the fusion it is known how many atoms are expected to be found this should be provided
	 * If this is not known it should be set to null and the smallest number of fusion atoms will be returned.
	 * @param ring
	 * @param numericalLocantsOnRing
	 * @param expectedNumberOfAtomsToBeUsedForFusion
	 * @return
	 * @throws StructureBuildingException
	 */
	private static List<Atom> determineAtomsToFuse(Fragment ring, List<String> numericalLocantsOnRing, Integer expectedNumberOfAtomsToBeUsedForFusion) throws StructureBuildingException {
		List<Atom> parentAtomList = ring.getAtomList();
		int indexfirst = parentAtomList.indexOf(ring.getAtomByLocantOrThrow(numericalLocantsOnRing.get(0)));
		int indexfinal = parentAtomList.indexOf(ring.getAtomByLocantOrThrow(numericalLocantsOnRing.get(numericalLocantsOnRing.size()-1)));
		CyclicAtomList cyclicRingAtomList = new CyclicAtomList(parentAtomList, indexfirst);
		List<Atom> fusionAtoms = null;
		
		List<Atom> potentialFusionAtomsAscending = new ArrayList<Atom>();
		potentialFusionAtomsAscending.add(cyclicRingAtomList.getCurrent());
		while (cyclicRingAtomList.getIndice() != indexfinal){//assume numbers are ascending
			potentialFusionAtomsAscending.add(cyclicRingAtomList.getNext());
		}
		if (expectedNumberOfAtomsToBeUsedForFusion ==null ||expectedNumberOfAtomsToBeUsedForFusion == potentialFusionAtomsAscending.size()){
			boolean notInPotentialParentAtoms =false;
			for (int i =1; i < numericalLocantsOnRing.size()-1 ; i ++){
				if (!potentialFusionAtomsAscending.contains(ring.getAtomByLocantOrThrow(numericalLocantsOnRing.get(i)))){
					notInPotentialParentAtoms =true;
				}
			}
			if (!notInPotentialParentAtoms){
				fusionAtoms = potentialFusionAtomsAscending;
			}
		}
		
		if (fusionAtoms ==null || expectedNumberOfAtomsToBeUsedForFusion ==null){//that didn't work, so try assuming the numbers are descending
			cyclicRingAtomList.setIndice(indexfirst);
			List<Atom> potentialFusionAtomsDescending = new ArrayList<Atom>();
			potentialFusionAtomsDescending.add(cyclicRingAtomList.getCurrent());
			while (cyclicRingAtomList.getIndice() != indexfinal){//assume numbers are descending
				potentialFusionAtomsDescending.add(cyclicRingAtomList.getPrevious());
			}
			if (expectedNumberOfAtomsToBeUsedForFusion ==null || expectedNumberOfAtomsToBeUsedForFusion == potentialFusionAtomsDescending.size()){
				boolean notInPotentialParentAtoms =false;
				for (int i =1; i < numericalLocantsOnRing.size()-1 ; i ++){
					if (!potentialFusionAtomsDescending.contains(ring.getAtomByLocantOrThrow(numericalLocantsOnRing.get(i)))){
						notInPotentialParentAtoms =true;
					}
				}
				if (!notInPotentialParentAtoms){
					if (fusionAtoms!=null && expectedNumberOfAtomsToBeUsedForFusion ==null){
						//prefer less fusion atoms
						if (potentialFusionAtomsDescending.size()< fusionAtoms.size()){
							fusionAtoms = potentialFusionAtomsDescending;
						}
					}
					else{
						fusionAtoms = potentialFusionAtomsDescending;
					}
				}
			}
		}
		return fusionAtoms;
	}

	/**
	 * Creates the bonds required to fuse two rings together.
	 * The child atoms are recorded as atoms that should be removed later
	 * @param childAtoms
	 * @param parentAtoms
	 * @throws StructureBuildingException
	 */
	private void fuseRings(List<Atom> childAtoms, List<Atom> parentAtoms) throws StructureBuildingException {
		if (parentAtoms.size()!=childAtoms.size()){
			throw new StructureBuildingException("Problem with fusion descriptors: Parent atoms specified: " + parentAtoms.size() +" Child atoms specified: " + childAtoms.size() + " These should have been identical!");
		}
		//replace parent atoms if the atom has already been used in fusion with the original atom
		//This will occur if fusion has resulted in something resembling a spiro centre e.g. cyclopenta[1,2-b:5,1-b']bis[1,4]oxathiine
		for (int i = parentAtoms.size() -1; i >=0; i--) {
			if (atomsToRemoveToReplacementAtom.get(parentAtoms.get(i))!=null){
				parentAtoms.set(i, atomsToRemoveToReplacementAtom.get(parentAtoms.get(i)));
			}
			if (atomsToRemoveToReplacementAtom.get(childAtoms.get(i))!=null){
				childAtoms.set(i, atomsToRemoveToReplacementAtom.get(childAtoms.get(i)));
			}
		}
		
		//sync spareValency and check that element type matches
		for (int i = 0; i < childAtoms.size(); i++) {
			Atom parentAtom = parentAtoms.get(i);
			Atom childAtom = childAtoms.get(i);
			if (childAtom.hasSpareValency()){
				parentAtom.setSpareValency(true);
			}
			if (!parentAtom.getElement().equals(childAtom.getElement())){
				throw new StructureBuildingException("Invalid fusion descriptor: Heteroatom placement is ambigous as it is not present in both components of the fusion");
			}
			atomsToRemoveToReplacementAtom.put(childAtom, parentAtom);
		}
		
		Set<Bond> fusionEdgeBonds  = new HashSet<Bond>();//these bonds already exist in both the child and parent atoms
		for (int i = 0; i < childAtoms.size() -1; i++) {
			fusionEdgeBonds.add(childAtoms.get(i).getBondToAtomOrThrow(childAtoms.get(i+1)));
			fusionEdgeBonds.add(parentAtoms.get(i).getBondToAtomOrThrow(parentAtoms.get(i+1)));
		}

		Set<Bond> bondsToAddToParentAtoms = new HashSet<Bond>();
		for (Atom childAtom : childAtoms) {
			for (Bond b : childAtom.getBonds()) {
				if (!fusionEdgeBonds.contains(b)){
					bondsToAddToParentAtoms.add(b);
				}
			}
		}
		
		Set<Bond> bondsToAddToChildAtoms = new HashSet<Bond>();
		for (Atom parentAtom : parentAtoms) {
			for (Bond b : parentAtom.getBonds()) {
				if (!fusionEdgeBonds.contains(b)){
					bondsToAddToChildAtoms.add(b);
				}
			}
		}
		
		for (Bond bond : bondsToAddToParentAtoms) {
			Atom from = bond.getFromAtom();
			int indiceInChildAtoms = childAtoms.indexOf(from);
			if (indiceInChildAtoms !=-1){
				from = parentAtoms.get(indiceInChildAtoms);
			}
			Atom to = bond.getToAtom();
			indiceInChildAtoms = childAtoms.indexOf(to);
			if (indiceInChildAtoms !=-1){
				to = parentAtoms.get(indiceInChildAtoms);
			}
			state.fragManager.createBond(from, to, 1);
		}

		for (Bond bond : bondsToAddToChildAtoms) {
			Atom from = bond.getFromAtom();
			int indiceInParentAtoms = parentAtoms.indexOf(from);
			if (indiceInParentAtoms !=-1){
				from = childAtoms.get(indiceInParentAtoms);
			}
			Atom to = bond.getToAtom();
			indiceInParentAtoms = parentAtoms.indexOf(to);
			if (indiceInParentAtoms !=-1){
				to = childAtoms.get(indiceInParentAtoms);
			}
			Bond newBond = new Bond(from, to, 1);
			if (childAtoms.contains(from)){
				from.addBond(newBond);
			}
			else{
				to.addBond(newBond);
			}
		}
	}

	/**
	 * Fuse the benzo with the subsequent ring
	 * Uses locants in front of the benz/benzo group to assign heteroatoms on the now numbered fused ring system
	 * @param benzoEl
	 * @param parentEl
	 * @throws StructureBuildingException
	 */
	private void benzoSpecificFusion(Element benzoEl, Element parentEl) throws StructureBuildingException {
		
		/*
		 * Perform the fusion, number it and associate it with the parentEl
		 */
		Fragment benzoRing = state.xmlFragmentMap.get(benzoEl);
		Fragment parentRing = state.xmlFragmentMap.get(parentEl);
		performSimpleFusion(null, benzoRing , parentRing);
		state.fragManager.incorporateFragment(benzoRing, parentRing);
		removeMergedAtoms();
		FusedRingNumberer.numberFusedRing(parentRing);//numbers the fused ring;
		Fragment fusedRing =parentRing;

		/*
		 * Check for locants and use these to set the heteroatom positions
		 */
		Element locantEl = (Element) XOMTools.getPreviousSibling(benzoEl);
		if (locantEl != null && locantEl.getLocalName().equals(LOCANT_EL)) {
			String[] locants = MATCH_COMMA.split(locantEl.getValue());
			Elements suffixes=((Element)benzoEl.getParent()).getChildElements(SUFFIX_EL);
			int suffixesWithoutLocants =0;
			for (int i = 0; i < suffixes.size(); i++) {
				if (suffixes.get(i).getAttribute(LOCANT_ATR)==null){
					suffixesWithoutLocants++;
				}
			}
			if (locants.length != suffixesWithoutLocants){//In preference locants will be assigned to suffixes rather than to this nomenclature
				List<Atom> atomList =fusedRing.getAtomList();
				LinkedList<Atom> heteroatoms =new LinkedList<Atom>();
				LinkedList<String> elementOfHeteroAtom =new LinkedList<String>();
				for (Atom atom : atomList) {//this iterates in the same order as the numbering system
					if (!atom.getElement().equals("C")){
						heteroatoms.add(atom);
						elementOfHeteroAtom.add(atom.getElement());
					}
				}
				if (locants.length == heteroatoms.size()){//as many locants as there are heteroatoms to assign
					for (Atom atom : heteroatoms) {
						atom.setElement("C");
					}
					FragmentTools.pickUpIndicatedHydrogens(fusedRing);
					for (int i=0; i< heteroatoms.size(); i ++){
						String elementSymbol =elementOfHeteroAtom.get(i);
						fusedRing.getAtomByLocantOrThrow(locants[i]).setElement(elementSymbol);
					}
					locantEl.detach();
				}
				else if (locants.length > 1){
					throw new StructureBuildingException("Unable to assign all locants to benzo-fused ring or multiplier was mising");
				}
			}
		}
	}
}
