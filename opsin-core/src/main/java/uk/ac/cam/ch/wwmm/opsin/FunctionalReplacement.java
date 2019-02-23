package uk.ac.cam.ch.wwmm.opsin;

import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;

/**
 * Methods for performing functional replacement
 * @author dl387
 *
 */
class FunctionalReplacement {

	/**
	 * Sorts infix transformations by the number of acceptable inputs for the transformation.
	 * e.g. thio ends up towards the end of the list as it accepts both -O or =O whilst say imido only accepts =O
	 * @author dl387
	 *
	 */
	private static class SortInfixTransformations implements Comparator<String> {
		public int compare(String infixTransformation1, String infixTransformation2) {
			int allowedInputs1 = infixTransformation1.split(",").length;
			int allowedInputs2 = infixTransformation2.split(",").length;
			if (allowedInputs1 < allowedInputs2){//infixTransformation1 preferred
				return -1;
			}
			if (allowedInputs1 > allowedInputs2){//infixTransformation2 preferred
				return 1;
			}
			else{
				return 0;
			}
		}
	}
	private static enum PREFIX_REPLACEMENT_TYPE{
		chalcogen,//ambiguous
		halideOrPseudoHalide,//only mean functional replacement when applied to non carboxylic acids
		dedicatedFunctionalReplacementPrefix,//no ambiguity exists
		hydrazono,//ambiguous, only applies to non carboxylic acid
		peroxy//ambiguous, also applies to etheric oxygen
	}	

	private static final Pattern matchChalcogenReplacement= Pattern.compile("thio|seleno|telluro");
	
	private final BuildState state;

	FunctionalReplacement(BuildState state) {
		this.state = state;
	}

	/**
	 * Applies the effects of acid replacing functional class nomenclature
	 * This must be performed early so that prefix/infix functional replacement is performed correctly
	 * and so that element symbol locants are assigned appropriately
	 * @param finalSubOrRootInWord
	 * @param word
	 * @throws ComponentGenerationException 
	 * @throws StructureBuildingException 
	 */
	void processAcidReplacingFunctionalClassNomenclature(Element finalSubOrRootInWord, Element word) throws ComponentGenerationException, StructureBuildingException {
		Element wordRule = OpsinTools.getParentWordRule(word);
		if (WordRule.valueOf(wordRule.getAttributeValue(WORDRULE_ATR)) == WordRule.acidReplacingFunctionalGroup){
			Element parentWordRule = word.getParent();
			if (parentWordRule.indexOf(word)==0){
				for (int i = 1, l = parentWordRule.getChildCount(); i < l ; i++) {
					Element acidReplacingWord = parentWordRule.getChild(i);
					if (!acidReplacingWord.getName().equals(WORD_EL)) {
						throw new RuntimeException("OPSIN bug: problem with acidReplacingFunctionalGroup word rule");
					}
					String type = acidReplacingWord.getAttributeValue(TYPE_ATR);
					if (type.equals(WordType.full.toString())) {
						//case where functionalTerm is substituted
						//as words are processed from right to left in cases like phosphoric acid tri(ethylamide) this will be phosphoric acid ethylamide ethylamide ethylamide
						processAcidReplacingFunctionalClassNomenclatureFullWord(finalSubOrRootInWord, acidReplacingWord);
					}
					else if (type.equals(WordType.functionalTerm.toString())) {
						processAcidReplacingFunctionalClassNomenclatureFunctionalWord(finalSubOrRootInWord, acidReplacingWord);
					}
					else {
						throw new RuntimeException("OPSIN bug: problem with acidReplacingFunctionalGroup word rule");
					}
				}
			}
		}
	}

	/**
	 * Performs prefix functional replacement e.g. thio in thioacetic acid replaces an O with S
	 * Prefixes will present themselves as substituents. There is potential ambiguity between usage as a substituent
	 * and as a functional replacement term in some cases. If the substituent is deemed to indicate functional replacement
	 * it will be detached and its effects applied to the subsequent group
	 * 
	 * The list of groups and substituents given to this method will be mutated in the process.
	 * 
	 * For heterocyclic rings functional replacement should technically be limited to :
	 * pyran, morpholine, chromene, isochromene and xanthene, chromane and isochromane.
	 * but this is not currently enforced
	 * @param groups
	 * @param substituents
	 * @return boolean: has any functional replacement occurred
	 * @throws StructureBuildingException
	 * @throws ComponentGenerationException 
	 */
	boolean processPrefixFunctionalReplacementNomenclature(List<Element> groups, List<Element> substituents) throws StructureBuildingException, ComponentGenerationException {
		int originalNumberOfGroups = groups.size();
		for (int i = originalNumberOfGroups-1; i >=0; i--) {
			Element group =groups.get(i);
			String groupValue = group.getValue();
			PREFIX_REPLACEMENT_TYPE replacementType = null;
			if (matchChalcogenReplacement.matcher(groupValue).matches() && !isChalcogenSubstituent(group) || groupValue.equals("thiono")){
				replacementType =PREFIX_REPLACEMENT_TYPE.chalcogen;
			}
			else if (HALIDEORPSEUDOHALIDE_SUBTYPE_VAL.equals(group.getAttributeValue(SUBTYPE_ATR))){
				replacementType =PREFIX_REPLACEMENT_TYPE.halideOrPseudoHalide;
			}
			else if (DEDICATEDFUNCTIONALREPLACEMENTPREFIX_SUBTYPE_VAL.equals(group.getAttributeValue(SUBTYPE_ATR))){
				replacementType =PREFIX_REPLACEMENT_TYPE.dedicatedFunctionalReplacementPrefix;
			}
			else if (groupValue.equals("hydrazono")){
				replacementType =PREFIX_REPLACEMENT_TYPE.hydrazono;
			}
			else if (groupValue.equals("peroxy")){
				replacementType =PREFIX_REPLACEMENT_TYPE.peroxy;
			}
			if (replacementType != null) {
				//need to check whether this is an instance of functional replacement by checking the substituent/root it is applying to
				Element substituent = group.getParent();
				Element nextSubOrBracket = OpsinTools.getNextSibling(substituent);
				if (nextSubOrBracket!=null && (nextSubOrBracket.getName().equals(ROOT_EL) || nextSubOrBracket.getName().equals(SUBSTITUENT_EL))){
					Element groupToBeModified = nextSubOrBracket.getFirstChildElement(GROUP_EL);
					if (groupPrecededByElementThatBlocksPrefixReplacementInterpetation(groupToBeModified)) {
						if (replacementType == PREFIX_REPLACEMENT_TYPE.dedicatedFunctionalReplacementPrefix){
							throw new ComponentGenerationException("dedicated Functional Replacement Prefix used in an inappropriate position :" + groupValue);
						}
						continue;//not 2,2'-thiodipyran
					}
					Element locantEl = null;//null unless a locant that agrees with the multiplier is present
					Element multiplierEl = null;
					int numberOfAtomsToReplace = 1;//the number of atoms to be functionally replaced, modified by a multiplier e.g. dithio
					Element possibleMultiplier = OpsinTools.getPreviousSibling(group);
					if (possibleMultiplier != null) {
						Element possibleLocant;
						if (possibleMultiplier.getName().equals(MULTIPLIER_EL)) {
							numberOfAtomsToReplace = Integer.valueOf(possibleMultiplier.getAttributeValue(VALUE_ATR));
							possibleLocant = OpsinTools.getPreviousSibling(possibleMultiplier);
							multiplierEl = possibleMultiplier;
						}
						else{
							possibleLocant = possibleMultiplier;
						}
						if (possibleLocant !=null && possibleLocant.getName().equals(LOCANT_EL) && possibleLocant.getAttribute(TYPE_ATR) == null) {
							int numberOfLocants = possibleLocant.getValue().split(",").length;
							if (numberOfLocants == numberOfAtomsToReplace){//locants and number of replacements agree
								locantEl = possibleLocant;
							}
							else if (numberOfAtomsToReplace > 1) {//doesn't look like prefix functional replacement
								if (replacementType  == PREFIX_REPLACEMENT_TYPE.dedicatedFunctionalReplacementPrefix){
									throw new ComponentGenerationException("dedicated Functional Replacement Prefix used in an inappropriate position :" + groupValue);
								}
								continue;
							}
						}
					}

					int oxygenReplaced;
					if (replacementType == PREFIX_REPLACEMENT_TYPE.chalcogen) {
						oxygenReplaced = performChalcogenFunctionalReplacement(groupToBeModified, locantEl, numberOfAtomsToReplace, group.getAttributeValue(VALUE_ATR));
					}
					else if (replacementType == PREFIX_REPLACEMENT_TYPE.peroxy) {
						if (nextSubOrBracket.getName().equals(SUBSTITUENT_EL)) {
							continue;
						}
						oxygenReplaced = performPeroxyFunctionalReplacement(groupToBeModified, locantEl, numberOfAtomsToReplace);
					}
					else if (replacementType == PREFIX_REPLACEMENT_TYPE.dedicatedFunctionalReplacementPrefix){
						if (!groupToBeModified.getAttributeValue(TYPE_ATR).equals(NONCARBOXYLICACID_TYPE_VAL)
								&& !(groupToBeModified.getValue().equals("form") && groupValue.equals("imido"))){
							throw new ComponentGenerationException("dedicated Functional Replacement Prefix used in an inappropriate position :" + groupValue);
						}
						oxygenReplaced = performFunctionalReplacementOnAcid(groupToBeModified, locantEl, numberOfAtomsToReplace, group.getAttributeValue(VALUE_ATR));
						if (oxygenReplaced==0){
							throw new ComponentGenerationException("dedicated Functional Replacement Prefix used in an inappropriate position :" + groupValue);
						}
					}
					else if (replacementType == PREFIX_REPLACEMENT_TYPE.hydrazono || replacementType == PREFIX_REPLACEMENT_TYPE.halideOrPseudoHalide){
						Fragment acidFrag = groupToBeModified.getFrag();
						if (!groupToBeModified.getAttributeValue(TYPE_ATR).equals(NONCARBOXYLICACID_TYPE_VAL) ||
								acidHasSufficientHydrogenForSubstitutionInterpretation(acidFrag, group.getFrag().getOutAtom(0).getValency(), locantEl)){
							//hydrazono replacement only applies to non carboxylic acids e.g. hydrazonooxalic acid
							//need to be careful to note that something like chlorophosphonic acid isn't functional replacement
							continue;
						}
						oxygenReplaced = performFunctionalReplacementOnAcid(groupToBeModified, locantEl, numberOfAtomsToReplace, group.getAttributeValue(VALUE_ATR));
					}
					else{
						throw new StructureBuildingException("OPSIN bug: Unexpected prefix replacement type");
					}
					if (oxygenReplaced>0){
						state.fragManager.removeFragment(group.getFrag());
						substituent.removeChild(group);
						groups.remove(group);
						List<Element> remainingChildren =substituent.getChildElements();//there may be a locant that should be moved
						for (int j = remainingChildren.size()-1; j>=0; j--){
							Element child =substituent.getChild(j);
							child.detach();
							nextSubOrBracket.insertChild(child, 0);
						}
						substituents.remove(substituent);
						substituent.detach();
						if (oxygenReplaced>1){
							multiplierEl.detach();
						}
					}
				}
				else if (replacementType  == PREFIX_REPLACEMENT_TYPE.dedicatedFunctionalReplacementPrefix){
					throw new ComponentGenerationException("dedicated Functional Replacement Prefix used in an inappropriate position :" + groupValue);
				}
			}
		}
		return groups.size() != originalNumberOfGroups;
	}
	
	private boolean isChalcogenSubstituent(Element group) {
		//Is this group followed by a hyphen and directly preceded by a substituent i.e. no multiplier/locant
		//e.g. methylthio-
		Element next = OpsinTools.getNextSibling(group);
		if (next != null && next.getName().equals(HYPHEN_EL) &&
				OpsinTools.getPreviousSibling(group) == null) {
			Element previousGroup = OpsinTools.getPreviousGroup(group);
			if (previousGroup != null) {
				//TODO We actually want to know if a carbon atom is the attachment point... but we don't know the attachment point locations at this point
				Element suffix = OpsinTools.getNextSibling(previousGroup, SUFFIX_EL);
				if (suffix == null || suffix.getFrag() == null) {
					for (Atom a : previousGroup.getFrag()) {
						if (a.getElement() == ChemEl.C) {
							return true;
						}
					}
				}
			}
		}
		return false;
	}

	/**
	 * Currently prefix replacement terms must be directly adjacent to the groupToBeModified with an exception made
	 * for carbohydrate stereochemistry prefixes e.g. 'gluco' and for substractive prefixes e.g. 'deoxy'
	 * @param groupToBeModified
	 * @return
	 */
	private boolean groupPrecededByElementThatBlocksPrefixReplacementInterpetation(Element groupToBeModified) {
		Element previous = OpsinTools.getPreviousSibling(groupToBeModified);
		while (previous !=null && (previous.getName().equals(SUBTRACTIVEPREFIX_EL)
				|| (previous.getName().equals(STEREOCHEMISTRY_EL) && previous.getAttributeValue(TYPE_ATR).equals(CARBOHYDRATECONFIGURATIONPREFIX_TYPE_VAL)))){
			previous = OpsinTools.getPreviousSibling(previous);
		}
		return previous != null;
	}
	
	
	/*
	 * 
	 */
	
	/**
	 * Performs functional replacement using infixes e.g. thio in ethanthioic acid replaces an O with S
	 * @param suffixFragments May be modified if a multiplier is determined to mean multiplication of a suffix, usually untouched
	 * @param suffixes The suffix elements  May be modified if a multiplier is determined to mean multiplication of a suffix, usually untouched
	 * @throws StructureBuildingException
	 * @throws ComponentGenerationException
	 */
	void processInfixFunctionalReplacementNomenclature(List<Element> suffixes, List<Fragment> suffixFragments) throws StructureBuildingException, ComponentGenerationException {
		for (int i = 0; i < suffixes.size(); i++) {
			Element suffix = suffixes.get(i);
			if (suffix.getAttribute(INFIX_ATR) != null){
				Fragment fragToApplyInfixTo = suffix.getFrag();
				Element possibleAcidGroup = OpsinTools.getPreviousSiblingIgnoringCertainElements(suffix, new String[]{MULTIPLIER_EL, INFIX_EL, SUFFIX_EL});
				if (possibleAcidGroup !=null && possibleAcidGroup.getName().equals(GROUP_EL) && 
						(possibleAcidGroup.getAttributeValue(TYPE_ATR).equals(NONCARBOXYLICACID_TYPE_VAL)|| possibleAcidGroup.getAttributeValue(TYPE_ATR).equals(CHALCOGENACIDSTEM_TYPE_VAL))){
					fragToApplyInfixTo = possibleAcidGroup.getFrag();
				}
				if (fragToApplyInfixTo ==null){
					throw new ComponentGenerationException("infix has erroneously been assigned to a suffix which does not correspond to a suffix fragment. suffix: " + suffix.getValue());
				}
				//e.g. =O:S,-O:S (which indicates replacing either a double or single bonded oxygen with S)
				//This is semicolon delimited for each infix
				List<String> infixTransformations = StringTools.arrayToList(suffix.getAttributeValue(INFIX_ATR).split(";"));

				List<Atom> atomList =fragToApplyInfixTo.getAtomList();
				LinkedList<Atom> singleBondedOxygen = new LinkedList<Atom>();
				LinkedList<Atom> doubleBondedOxygen = new LinkedList<Atom>();
				populateTerminalSingleAndDoubleBondedOxygen(atomList, singleBondedOxygen, doubleBondedOxygen);
				int oxygenAvailable = singleBondedOxygen.size() +doubleBondedOxygen.size();

				/*
				 * Modifies suffixes, suffixFragments, suffix and infixTransformations as appropriate
				 */
				disambiguateMultipliedInfixMeaning(suffixes, suffixFragments, suffix, infixTransformations, oxygenAvailable);

				/*
				 * Sort infixTransformations so more specific transformations are performed first
				 * e.g. ethanthioimidic acid-->ethanimidthioic acid as imid can only apply to the double bonded oxygen
				 */
				Collections.sort(infixTransformations, new SortInfixTransformations());

				for (String infixTransformation : infixTransformations) {
					String[] transformationArray = infixTransformation.split(":");
					if (transformationArray.length !=2){
						throw new StructureBuildingException("Atom to be replaced and replacement not specified correctly in infix: " + infixTransformation);
					}
					String[] transformations = transformationArray[0].split(",");
					String replacementSMILES = transformationArray[1];
					boolean acceptDoubleBondedOxygen = false;
					boolean acceptSingleBondedOxygen = false;
					boolean nitrido =false;
					for (String transformation : transformations) {
						if (transformation.startsWith("=")){
							acceptDoubleBondedOxygen = true;
						}
						else if (transformation.startsWith("-")){
							acceptSingleBondedOxygen = true;
						}
						else if (transformation.startsWith("#")){
							nitrido =true;
						}
						else{
							throw new StructureBuildingException("Malformed infix transformation. Expected to start with either - or =. Transformation was: " +transformation);
						}
						if (transformation.length()<2 || transformation.charAt(1)!='O'){
							throw new StructureBuildingException("Only replacement by oxygen is supported. Check infix defintions");
						}
					}
					boolean infixAssignmentAmbiguous =false;
					if ((acceptSingleBondedOxygen ||nitrido)  && !acceptDoubleBondedOxygen){
						if (singleBondedOxygen.size() ==0){
							throw new StructureBuildingException("Cannot find single bonded oxygen for infix with SMILES: "+ replacementSMILES+ " to modify!");
						}
						if (singleBondedOxygen.size() !=1){
							infixAssignmentAmbiguous=true;
						}
					}
					if (!acceptSingleBondedOxygen && (acceptDoubleBondedOxygen | nitrido)){
						if (doubleBondedOxygen.size()==0){
							throw new StructureBuildingException("Cannot find double bonded oxygen for infix with SMILES: "+ replacementSMILES+ " to modify!");
						}
						if (doubleBondedOxygen.size() != 1){
							infixAssignmentAmbiguous=true;
						}
					}
					if (acceptSingleBondedOxygen && acceptDoubleBondedOxygen){
						if (oxygenAvailable ==0){
							throw new StructureBuildingException("Cannot find oxygen for infix with SMILES: "+ replacementSMILES+ " to modify!");
						}
						if (oxygenAvailable !=1){
							infixAssignmentAmbiguous=true;
						}
					}

					Set<Atom> ambiguousElementAtoms = new LinkedHashSet<Atom>();
					Atom atomToUse = null;
					if ((acceptDoubleBondedOxygen || nitrido) && doubleBondedOxygen.size()>0 ){
						atomToUse = doubleBondedOxygen.removeFirst();
					}
					else if (acceptSingleBondedOxygen && singleBondedOxygen.size()>0 ){
						atomToUse = singleBondedOxygen.removeFirst();
					}
					else{
						throw new StructureBuildingException("Cannot find oxygen for infix with SMILES: "+ replacementSMILES+ " to modify!");//this would be a bug
					}
					Fragment replacementFrag = state.fragManager.buildSMILES(replacementSMILES, SUFFIX_TYPE_VAL, NONE_LABELS_VAL);
					if (replacementFrag.getOutAtomCount()>0){//SMILES include an indication of the bond order the replacement fragment will have, this is not intended to be an outatom
						replacementFrag.removeOutAtom(0);
					}
					Atom atomThatWillReplaceOxygen =replacementFrag.getFirstAtom();
					if (replacementFrag.getAtomCount()==1 && atomThatWillReplaceOxygen.getElement().isChalcogen()){
						atomThatWillReplaceOxygen.setCharge(atomToUse.getCharge());
						atomThatWillReplaceOxygen.setProtonsExplicitlyAddedOrRemoved(atomToUse.getProtonsExplicitlyAddedOrRemoved());
					}
					removeOrMoveObsoleteFunctionalAtoms(atomToUse, replacementFrag);//also will move charge if necessary
					moveObsoleteOutAtoms(atomToUse, replacementFrag);//if the replaced atom was an outatom the fragments outatom list need to be corrected
					if (nitrido){
						atomToUse.getFirstBond().setOrder(3);
						Atom removedHydroxy = singleBondedOxygen.removeFirst();
						state.fragManager.removeAtomAndAssociatedBonds(removedHydroxy);
						removeAssociatedFunctionalAtom(removedHydroxy);
					}
					state.fragManager.incorporateFragment(replacementFrag, atomToUse.getFrag());
					state.fragManager.replaceAtomWithAnotherAtomPreservingConnectivity(atomToUse, atomThatWillReplaceOxygen);
					if (infixAssignmentAmbiguous){
						ambiguousElementAtoms.add(atomThatWillReplaceOxygen);
						if (atomThatWillReplaceOxygen.getProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT)!=null){
							ambiguousElementAtoms.addAll(atomThatWillReplaceOxygen.getProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT));
						}
					}
					if (infixAssignmentAmbiguous){//record what atoms could have been replaced. Often this ambiguity is resolved later e.g. S-methyl ethanthioate
						for (Atom a : doubleBondedOxygen) {
							ambiguousElementAtoms.add(a);
							if (a.getProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT)!=null){
								ambiguousElementAtoms.addAll(a.getProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT));
							}
						}
						for (Atom a : singleBondedOxygen) {
							ambiguousElementAtoms.add(a);
							if (a.getProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT)!=null){
								ambiguousElementAtoms.addAll(a.getProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT));
							}
						}
						for (Atom atom : ambiguousElementAtoms) {
							atom.setProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT, ambiguousElementAtoms);
						}
					}
				}
			}
		}
	}

	/*
	 * Functional class nomenclature
	 */
	
	/**
	 * Replaces the appropriate number of functional oxygen atoms with the corresponding fragment
	 * @param acidContainingRoot
	 * @param acidReplacingWord
	 * @throws ComponentGenerationException
	 * @throws StructureBuildingException
	 */
	private void processAcidReplacingFunctionalClassNomenclatureFullWord(Element acidContainingRoot, Element acidReplacingWord) throws ComponentGenerationException, StructureBuildingException {
		String locant = acidReplacingWord.getAttributeValue(LOCANT_ATR);
		Element acidReplacingGroup = StructureBuildingMethods.findRightMostGroupInBracket(acidReplacingWord);
		if (acidReplacingGroup ==null){
			throw new ComponentGenerationException("OPSIN bug: acid replacing group not found where one was expected for acidReplacingFunctionalGroup wordRule");
		}
		String functionalGroupName = acidReplacingGroup.getValue();
		Fragment acidReplacingFrag = acidReplacingGroup.getFrag();
		if (acidReplacingGroup.getParent().getChildCount() != 1){
			throw new ComponentGenerationException("Unexpected qualifier to: " + functionalGroupName);
		}
		
		Element groupToBeModified = acidContainingRoot.getFirstChildElement(GROUP_EL);
		List<Atom> oxygenAtoms = findFunctionalOxygenAtomsInApplicableSuffixes(groupToBeModified);
		if (oxygenAtoms.size() == 0){
			oxygenAtoms = findFunctionalOxygenAtomsInGroup(groupToBeModified);
		}
		if (oxygenAtoms.size() == 0){
			List<Element> conjunctiveSuffixElements =OpsinTools.getNextSiblingsOfType(groupToBeModified, CONJUNCTIVESUFFIXGROUP_EL);
			for (Element conjunctiveSuffixElement : conjunctiveSuffixElements) {
				oxygenAtoms.addAll(findFunctionalOxygenAtomsInGroup(conjunctiveSuffixElement));
			}
		}
		if (oxygenAtoms.size() < 1){
			throw new ComponentGenerationException("Insufficient oxygen to replace with " + functionalGroupName +"s in " + acidContainingRoot.getFirstChildElement(GROUP_EL).getValue());
		}
		
		boolean isAmide = functionalGroupName.equals("amide") || functionalGroupName.equals("amid");
		if (isAmide) {
			if (acidReplacingFrag.getAtomCount()!=1){
				throw new ComponentGenerationException("OPSIN bug: " + functionalGroupName + " not found where expected");
			}
			Atom amideNitrogen = acidReplacingFrag.getFirstAtom();
			amideNitrogen.neutraliseCharge();
			amideNitrogen.clearLocants();
			acidReplacingFrag.addMappingToAtomLocantMap("N", amideNitrogen);
		}
		Atom chosenOxygen = locant != null ? removeOxygenWithAppropriateLocant(oxygenAtoms, locant) : oxygenAtoms.get(0);
		state.fragManager.replaceAtomWithAnotherAtomPreservingConnectivity(chosenOxygen, acidReplacingFrag.getFirstAtom());
		removeAssociatedFunctionalAtom(chosenOxygen);
	}
	

	/**
	 * Replaces the appropriate number of functional oxygen atoms with the corresponding fragment
	 * @param acidContainingRoot
	 * @param functionalWord
	 * @throws ComponentGenerationException
	 * @throws StructureBuildingException
	 */
	private void processAcidReplacingFunctionalClassNomenclatureFunctionalWord(Element acidContainingRoot, Element functionalWord) throws ComponentGenerationException, StructureBuildingException {
		if (functionalWord !=null && functionalWord.getAttributeValue(TYPE_ATR).equals(WordType.functionalTerm.toString())){
			Element functionalTerm = functionalWord.getFirstChildElement(FUNCTIONALTERM_EL);
			if (functionalTerm ==null){
				throw new ComponentGenerationException("OPSIN bug: functionalTerm word not found where one was expected for acidReplacingFunctionalGroup wordRule");
			}
			Element acidReplacingGroup = functionalTerm.getFirstChildElement(FUNCTIONALGROUP_EL);
			String functionalGroupName = acidReplacingGroup.getValue();
			Element possibleLocantOrMultiplier = OpsinTools.getPreviousSibling(acidReplacingGroup);
			int numberOfAcidicHydroxysToReplace = 1;
			String[] locants = null;
			if (possibleLocantOrMultiplier != null){
				if (possibleLocantOrMultiplier.getName().equals(MULTIPLIER_EL)){
					numberOfAcidicHydroxysToReplace = Integer.parseInt(possibleLocantOrMultiplier.getAttributeValue(VALUE_ATR));
					possibleLocantOrMultiplier.detach();
					possibleLocantOrMultiplier = OpsinTools.getPreviousSibling(acidReplacingGroup);
				}
				if (possibleLocantOrMultiplier != null){
					if (possibleLocantOrMultiplier.getName().equals(LOCANT_EL)){
						locants = StringTools.removeDashIfPresent(possibleLocantOrMultiplier.getValue()).split(",");
						possibleLocantOrMultiplier.detach();
					}
					else {
						throw new ComponentGenerationException("Unexpected qualifier to acidReplacingFunctionalGroup functionalTerm");
					}
				}
			}
			if (functionalTerm.getChildCount() != 1){
				throw new ComponentGenerationException("Unexpected qualifier to acidReplacingFunctionalGroup functionalTerm");
			}
			
			Element groupToBeModified = acidContainingRoot.getFirstChildElement(GROUP_EL);
			List<Atom> oxygenAtoms = findFunctionalOxygenAtomsInApplicableSuffixes(groupToBeModified);
			if (oxygenAtoms.size()==0) {
				oxygenAtoms = findFunctionalOxygenAtomsInGroup(groupToBeModified);
			}
			if (oxygenAtoms.size()==0) {
				List<Element> conjunctiveSuffixElements =OpsinTools.getNextSiblingsOfType(groupToBeModified, CONJUNCTIVESUFFIXGROUP_EL);
				for (Element conjunctiveSuffixElement : conjunctiveSuffixElements) {
					oxygenAtoms.addAll(findFunctionalOxygenAtomsInGroup(conjunctiveSuffixElement));
				}
			}
			if (numberOfAcidicHydroxysToReplace > oxygenAtoms.size()){
				throw new ComponentGenerationException("Insufficient oxygen to replace with nitrogen in " + acidContainingRoot.getFirstChildElement(GROUP_EL).getValue());
			}
			boolean isAmide = functionalGroupName.equals("amide") || functionalGroupName.equals("amid");
			if (isAmide) {
				for (int i = 0; i < numberOfAcidicHydroxysToReplace; i++) {
					Atom functionalOxygenToReplace = locants != null ? removeOxygenWithAppropriateLocant(oxygenAtoms, locants[i]) : oxygenAtoms.get(i);
					removeAssociatedFunctionalAtom(functionalOxygenToReplace);
					functionalOxygenToReplace.setElement(ChemEl.N);
				}
			}
			else{
				String groupValue = acidReplacingGroup.getAttributeValue(VALUE_ATR);
				String labelsValue = acidReplacingGroup.getAttributeValue(LABELS_ATR);
				Fragment acidReplacingFrag = state.fragManager.buildSMILES(groupValue, SUFFIX_TYPE_VAL, labelsValue != null ? labelsValue : NONE_LABELS_VAL);
				Fragment acidFragment = groupToBeModified.getFrag();
				if (acidFragment.hasLocant("2")){//prefer numeric locants on group to those of replacing group
					for (Atom atom : acidReplacingFrag.getAtomList()) {
						atom.clearLocants();
					}
				}
				Atom firstFunctionalOxygenToReplace = locants != null ? removeOxygenWithAppropriateLocant(oxygenAtoms, locants[0]) : oxygenAtoms.get(0);
				state.fragManager.replaceAtomWithAnotherAtomPreservingConnectivity(firstFunctionalOxygenToReplace, acidReplacingFrag.getFirstAtom());
				removeAssociatedFunctionalAtom(firstFunctionalOxygenToReplace);
				for (int i = 1; i < numberOfAcidicHydroxysToReplace; i++) {
					Fragment clonedHydrazide = state.fragManager.copyAndRelabelFragment(acidReplacingFrag, i);
					Atom functionalOxygenToReplace = locants != null ? removeOxygenWithAppropriateLocant(oxygenAtoms, locants[i]) : oxygenAtoms.get(i);
					state.fragManager.replaceAtomWithAnotherAtomPreservingConnectivity(functionalOxygenToReplace, clonedHydrazide.getFirstAtom());
					state.fragManager.incorporateFragment(clonedHydrazide, functionalOxygenToReplace.getFrag());
					removeAssociatedFunctionalAtom(functionalOxygenToReplace);
				}
				state.fragManager.incorporateFragment(acidReplacingFrag, firstFunctionalOxygenToReplace.getFrag());
			}
		}
		else{
			throw new ComponentGenerationException("amide word not found where expected, bug?");
		}
	}
	
	private Atom removeOxygenWithAppropriateLocant(List<Atom> oxygenAtoms, String locant) throws ComponentGenerationException {
		for (Iterator<Atom> iterator = oxygenAtoms.iterator(); iterator.hasNext();) {
			Atom atom = iterator.next();
			if (atom.hasLocant(locant)) {
				iterator.remove();
				return atom;
			}	
		}
		//Look for the case whether the locant refers to the backbone
		for (Iterator<Atom> iterator = oxygenAtoms.iterator(); iterator.hasNext();) {
			Atom atom = iterator.next();
			if (OpsinTools.depthFirstSearchForNonSuffixAtomWithLocant(atom, locant) != null){
				iterator.remove();
				return atom;
			}	
		}
		throw new ComponentGenerationException("Failed to find acid group at locant: " + locant);
	}
	

	/*
	 * Prefix functional replacement nomenclature
	 */

	
	private boolean acidHasSufficientHydrogenForSubstitutionInterpretation(Fragment acidFrag, int hydrogenRequiredForSubstitutionInterpretation, Element locantEl) {
		List<Atom> atomsThatWouldBeSubstituted = new ArrayList<Atom>();
		if (locantEl !=null){
			String[] possibleLocants = locantEl.getValue().split(",");
			for (String locant : possibleLocants) {
				Atom atomToBeSubstituted = acidFrag.getAtomByLocant(locant);
				if (atomToBeSubstituted !=null){
					atomsThatWouldBeSubstituted.add(atomToBeSubstituted);
				}
				else{
					atomsThatWouldBeSubstituted.clear();
					atomsThatWouldBeSubstituted.add(acidFrag.getDefaultInAtomOrFirstAtom());
					break;
				}
			}
		}
		else{
			atomsThatWouldBeSubstituted.add(acidFrag.getDefaultInAtomOrFirstAtom());
		}
		for (Atom atom : atomsThatWouldBeSubstituted) {
			if (StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(atom) < hydrogenRequiredForSubstitutionInterpretation){
				return false;//insufficient hydrogens for substitution interpretation
			}
		}
		return true;
	}

	/**
	 * Performs replacement of oxygen atoms by chalogen atoms
	 * If this is ambiguous e.g. thioacetate then Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT is populated
	 * @param groupToBeModified
	 * @param locantEl
	 * @param numberOfAtomsToReplace
	 * @param replacementSmiles
	 * @return
	 * @throws StructureBuildingException
	 */
	private int performChalcogenFunctionalReplacement(Element groupToBeModified, Element locantEl, int numberOfAtomsToReplace, String replacementSmiles) throws StructureBuildingException {
		List<Atom> oxygenAtoms = findOxygenAtomsInApplicableSuffixes(groupToBeModified);
		if (oxygenAtoms.size() == 0) {
			oxygenAtoms = findOxygenAtomsInGroup(groupToBeModified);
		}
		if (locantEl != null) {//locants are used to indicate replacement on trivial groups
			List<Atom> oxygenWithAppropriateLocants = pickOxygensWithAppropriateLocants(locantEl, oxygenAtoms);
			if(oxygenWithAppropriateLocants.size() < numberOfAtomsToReplace) {
				numberOfAtomsToReplace = 1;
				//e.g. -1-thioureidomethyl
			}
			else{
				locantEl.detach();
				oxygenAtoms = oxygenWithAppropriateLocants;
			}
		}
		List<Atom> replaceableAtoms = new ArrayList<Atom>();
		if (replacementSmiles.startsWith("=")) {
			//e.g. thiono
			replacementSmiles = replacementSmiles.substring(1);
			for (Atom oxygen : oxygenAtoms) {
				int incomingValency = oxygen.getIncomingValency();
				int bondCount = oxygen.getBondCount();
				if (bondCount == 1 && incomingValency == 2) {
					replaceableAtoms.add(oxygen);
				}
			}
		}
		else {
			List<Atom> doubleBondedOxygen = new ArrayList<Atom>();
			List<Atom> singleBondedOxygen = new ArrayList<Atom>();
			List<Atom> ethericOxygen = new ArrayList<Atom>();
			for (Atom oxygen : oxygenAtoms) {
				int incomingValency = oxygen.getIncomingValency();
				int bondCount = oxygen.getBondCount();
				if (bondCount == 1 && incomingValency ==2 ) {
					doubleBondedOxygen.add(oxygen);
				}
				else if (bondCount == 1 && incomingValency == 1) {
					singleBondedOxygen.add(oxygen);
				}
				else if (bondCount == 2 && incomingValency == 2) {
					ethericOxygen.add(oxygen);
				}
			}
			replaceableAtoms.addAll(doubleBondedOxygen);
			replaceableAtoms.addAll(singleBondedOxygen);
			replaceableAtoms.addAll(ethericOxygen);
		}

		int totalOxygen = replaceableAtoms.size();
		if (numberOfAtomsToReplace >1){
			if (totalOxygen < numberOfAtomsToReplace){
				numberOfAtomsToReplace=1;
			}
		}

		int atomsReplaced =0;
		if (totalOxygen >=numberOfAtomsToReplace){//check that there atleast as many oxygens as requested replacements
			boolean prefixAssignmentAmbiguous =false;
			Set<Atom> ambiguousElementAtoms = new LinkedHashSet<Atom>();
			if (totalOxygen != numberOfAtomsToReplace){
				prefixAssignmentAmbiguous=true;
			}

			for (Atom atomToReplace : replaceableAtoms) {
				if (atomsReplaced == numberOfAtomsToReplace){
					ambiguousElementAtoms.add(atomToReplace);
					continue;
				}
				else{
					state.fragManager.replaceAtomWithSmiles(atomToReplace, replacementSmiles);
					if (prefixAssignmentAmbiguous){
						ambiguousElementAtoms.add(atomToReplace);
					}
				}
				atomsReplaced++;
			}

			if (prefixAssignmentAmbiguous){//record what atoms could have been replaced. Often this ambiguity is resolved later e.g. S-methyl thioacetate
				for (Atom atom : ambiguousElementAtoms) {
					atom.setProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT, ambiguousElementAtoms);
				}
			}
		}
		return atomsReplaced;
	}


	/**
	 * Converts functional oxygen to peroxy e.g. peroxybenzoic acid
	 * Returns the number of oxygen replaced
	 * @param groupToBeModified
	 * @param locantEl
	 * @param numberOfAtomsToReplace
	 * @return
	 * @throws StructureBuildingException
	 */
	private int performPeroxyFunctionalReplacement(Element groupToBeModified, Element locantEl, int numberOfAtomsToReplace) throws StructureBuildingException {
		List<Atom> oxygenAtoms = findFunctionalOxygenAtomsInApplicableSuffixes(groupToBeModified);
		if (oxygenAtoms.size()==0){
			oxygenAtoms = findEthericOxygenAtomsInGroup(groupToBeModified);
			oxygenAtoms.addAll(findFunctionalOxygenAtomsInGroup(groupToBeModified));
		}
		if (locantEl !=null){
			List<Atom> oxygenWithAppropriateLocants = pickOxygensWithAppropriateLocants(locantEl, oxygenAtoms);
			if(oxygenWithAppropriateLocants.size() < numberOfAtomsToReplace){
				numberOfAtomsToReplace =1;
			}
			else{
				locantEl.detach();
				oxygenAtoms = oxygenWithAppropriateLocants;
			}
		}
		if (numberOfAtomsToReplace >1 && oxygenAtoms.size() < numberOfAtomsToReplace){
			numberOfAtomsToReplace=1;
		}
		int atomsReplaced = 0;
		if (oxygenAtoms.size() >=numberOfAtomsToReplace){//check that there atleast as many oxygens as requested replacements
			atomsReplaced = numberOfAtomsToReplace;
			for (int j = 0; j < numberOfAtomsToReplace; j++) {
				Atom oxygenToReplace = oxygenAtoms.get(j);
				if (oxygenToReplace.getBondCount()==2){//etheric oxygen
					Fragment newOxygen = state.fragManager.buildSMILES("O", SUFFIX_TYPE_VAL, NONE_LABELS_VAL);
					Bond bondToRemove = oxygenToReplace.getFirstBond();
					Atom atomToAttachTo = bondToRemove.getFromAtom() == oxygenToReplace ?  bondToRemove.getToAtom() :  bondToRemove.getFromAtom();
					state.fragManager.createBond(atomToAttachTo, newOxygen.getFirstAtom(), 1);
					state.fragManager.createBond(newOxygen.getFirstAtom(), oxygenToReplace, 1);
					state.fragManager.removeBond(bondToRemove);
					state.fragManager.incorporateFragment(newOxygen, groupToBeModified.getFrag());
				}
				else{
					Fragment replacementFrag = state.fragManager.buildSMILES("OO", SUFFIX_TYPE_VAL, NONE_LABELS_VAL);
					removeOrMoveObsoleteFunctionalAtoms(oxygenToReplace, replacementFrag);
					state.fragManager.replaceAtomWithAnotherAtomPreservingConnectivity(oxygenToReplace, replacementFrag.getFirstAtom());
					state.fragManager.incorporateFragment(replacementFrag, groupToBeModified.getFrag());
				}
			}
		}
		return atomsReplaced;
	}
	
	/**
	 * Replaces double bonded oxygen and/or single bonded oxygen depending on the input SMILES
	 * SMILES with a valency 1 outAtom replace -O, SMILES with a valency 2 outAtom replace =O
	 * SMILES with a valency 3 outAtom replace -O and =O (nitrido)
	 * Returns the number of oxygen replaced
	 * @param groupToBeModified
	 * @param locantEl
	 * @param numberOfAtomsToReplace
	 * @param replacementSmiles
     * @return
	 * @throws StructureBuildingException
	 */
	private int performFunctionalReplacementOnAcid(Element groupToBeModified, Element locantEl, int numberOfAtomsToReplace, String replacementSmiles) throws StructureBuildingException {
		int outValency;
		if (replacementSmiles.startsWith("-")){
			outValency =1;
		}
		else if (replacementSmiles.startsWith("=")){
			outValency =2;
		}
		else if (replacementSmiles.startsWith("#")){
			outValency =3;
		}
		else{
			throw new StructureBuildingException("OPSIN bug: Unexpected valency on fragment for prefix functional replacement");
		}
		replacementSmiles = replacementSmiles.substring(1);
		List<Atom> oxygenAtoms = findOxygenAtomsInApplicableSuffixes(groupToBeModified);
		if (oxygenAtoms.size()==0){
			oxygenAtoms = findOxygenAtomsInGroup(groupToBeModified);
		}
		if (locantEl !=null){//locants are used to indicate replacement on trivial groups
			List<Atom> oxygenWithAppropriateLocants = pickOxygensWithAppropriateLocants(locantEl, oxygenAtoms);
			List<Atom> singleBondedOxygen = new ArrayList<Atom>();
			List<Atom> terminalDoubleBondedOxygen = new ArrayList<Atom>();
			populateTerminalSingleAndDoubleBondedOxygen(oxygenWithAppropriateLocants, singleBondedOxygen, terminalDoubleBondedOxygen);
			if (outValency ==1){
				oxygenWithAppropriateLocants.removeAll(terminalDoubleBondedOxygen);
			}
			else if (outValency ==2){
				oxygenWithAppropriateLocants.removeAll(singleBondedOxygen);
			}
			if(oxygenWithAppropriateLocants.size() < numberOfAtomsToReplace){
				numberOfAtomsToReplace =1;
				//e.g. -1-thioureidomethyl
			}
			else{
				locantEl.detach();
				oxygenAtoms = oxygenWithAppropriateLocants;
			}
		}
		List<Atom> singleBondedOxygen = new ArrayList<Atom>();
		List<Atom> terminalDoubleBondedOxygen = new ArrayList<Atom>();
		populateTerminalSingleAndDoubleBondedOxygen(oxygenAtoms, singleBondedOxygen, terminalDoubleBondedOxygen);
		if (outValency ==1){
			oxygenAtoms.removeAll(terminalDoubleBondedOxygen);
		}
		else if (outValency ==2){
			oxygenAtoms.removeAll(singleBondedOxygen);
			//favour bridging oxygen over double bonded oxygen c.f. imidodicarbonate
			oxygenAtoms.removeAll(terminalDoubleBondedOxygen);
			oxygenAtoms.addAll(terminalDoubleBondedOxygen);
		}
		else {
			if (singleBondedOxygen.size()==0 || terminalDoubleBondedOxygen.size()==0){
				throw new StructureBuildingException("Both a -OH and =O are required for nitrido prefix functional replacement");
			}
			oxygenAtoms.removeAll(singleBondedOxygen);
		}
		if (numberOfAtomsToReplace >1 && oxygenAtoms.size() < numberOfAtomsToReplace){
			numberOfAtomsToReplace=1;
		}

		int atomsReplaced =0;
		if (oxygenAtoms.size() >=numberOfAtomsToReplace){//check that there atleast as many oxygens as requested replacements
			for (Atom atomToReplace : oxygenAtoms) {
				if (atomsReplaced == numberOfAtomsToReplace){
					continue;
				}
				else{
					Fragment replacementFrag = state.fragManager.buildSMILES(replacementSmiles, atomToReplace.getFrag().getTokenEl(), NONE_LABELS_VAL);
					if (outValency ==3){//special case for nitrido
						atomToReplace.getFirstBond().setOrder(3);
						Atom removedHydroxy = singleBondedOxygen.remove(0);
						state.fragManager.removeAtomAndAssociatedBonds(removedHydroxy);
						removeAssociatedFunctionalAtom(removedHydroxy);
					}
					state.fragManager.replaceAtomWithAnotherAtomPreservingConnectivity(atomToReplace, replacementFrag.getFirstAtom());
					if (outValency ==1){
						removeOrMoveObsoleteFunctionalAtoms(atomToReplace, replacementFrag);
					}
					moveObsoleteOutAtoms(atomToReplace, replacementFrag);
					state.fragManager.incorporateFragment(replacementFrag, atomToReplace.getFrag());
				}
				atomsReplaced++;
			}
		}
		return atomsReplaced;
	}
	
	/*
	 * Infix functional replacement nomenclature
	 */

	/**
	 * This block handles infix multiplication. Unless brackets are provided this is ambiguous without knowledge of the suffix that is being modified
	 * For example butandithione could be intepreted as butandi(thione) or butan(dithi)one.
	 * Obviously the latter is wrong in this case but it is the correct interpretation for butandithiate
	 * @param suffixes
	 * @param suffixFragments
	 * @param suffix
	 * @param infixTransformations
	 * @param oxygenAvailable
	 * @throws ComponentGenerationException
	 * @throws StructureBuildingException
	 */
	private void disambiguateMultipliedInfixMeaning(List<Element> suffixes,
			List<Fragment> suffixFragments,Element suffix, List<String> infixTransformations, int oxygenAvailable)
			throws ComponentGenerationException, StructureBuildingException {
		Element possibleInfix = OpsinTools.getPreviousSibling(suffix);
		if (possibleInfix.getName().equals(INFIX_EL)){//the infix is only left when there was ambiguity
			Element possibleMultiplier = OpsinTools.getPreviousSibling(possibleInfix);
			if (possibleMultiplier.getName().equals(MULTIPLIER_EL)){
				int multiplierValue =Integer.parseInt(possibleMultiplier.getAttributeValue(VALUE_ATR));
				if (infixTransformations.size() + multiplierValue-1 <=oxygenAvailable){//multiplier means multiply the infix e.g. butandithiate
					for (int j = 1; j < multiplierValue; j++) {
						infixTransformations.add(0, infixTransformations.get(0));
					}
				}
				else{
					Element possibleLocant = OpsinTools.getPreviousSibling(possibleMultiplier);
					String[] locants = null;
					if (possibleLocant.getName().equals(LOCANT_EL)) {
						locants = possibleLocant.getValue().split(",");
					}
					if (locants !=null){
						if (locants.length!=multiplierValue){
							throw new ComponentGenerationException("Multiplier/locant disagreement when multiplying infixed suffix");
						}
					    suffix.addAttribute(new Attribute(LOCANT_ATR, locants[0]));
					}
					suffix.addAttribute(new Attribute(MULTIPLIED_ATR, "multiplied"));
					for (int j = 1; j < multiplierValue; j++) {//multiplier means multiply the infixed suffix e.g. butandithione
						Element newSuffix = suffix.copy();
						Fragment newSuffixFrag = state.fragManager.copyFragment(suffix.getFrag());
						newSuffix.setFrag(newSuffixFrag);
						suffixFragments.add(newSuffixFrag);
						OpsinTools.insertAfter(suffix, newSuffix);
						suffixes.add(newSuffix);
						if (locants !=null){//assign locants if available
							newSuffix.getAttribute(LOCANT_ATR).setValue(locants[j]);
						}
					}
					if (locants!=null){
						possibleLocant.detach();
					}
				}
				possibleMultiplier.detach();
				possibleInfix.detach();
			}
			else{
				throw new ComponentGenerationException("Multiplier expected in front of ambiguous infix");
			}
		}
	}
	
	/*
	 * Convenience Methods
	 */
	
	/**
	 * Given an atom that is to be replaced by a functional replacement fragment
	 * determines whether this atom is a functional atom and, if it is, performs the following processes:
	 * The functionalAtom is removed. If the the replacement fragment is an atom of O/S/Se/Te or the
	 * the terminal atom of the fragment is a single bonded O/S/Se/Te a functionAom is added to this atom.
	 * @param atomToBeReplaced
	 * @param replacementFrag
	 */
	private void removeOrMoveObsoleteFunctionalAtoms(Atom atomToBeReplaced, Fragment replacementFrag){
		List<Atom> replacementAtomList = replacementFrag.getAtomList();
		Fragment origFrag = atomToBeReplaced.getFrag();
		for (int i = origFrag.getFunctionalAtomCount() - 1; i >=0; i--) {
			FunctionalAtom functionalAtom = origFrag.getFunctionalAtom(i);
			if (atomToBeReplaced.equals(functionalAtom.getAtom())){
				atomToBeReplaced.getFrag().removeFunctionalAtom(i);
				Atom terminalAtomOfReplacementFrag = replacementAtomList.get(replacementAtomList.size()-1);
				if ((terminalAtomOfReplacementFrag.getIncomingValency() ==1 || replacementAtomList.size()==1)&& terminalAtomOfReplacementFrag.getElement().isChalcogen()){
					replacementFrag.addFunctionalAtom(terminalAtomOfReplacementFrag);
					terminalAtomOfReplacementFrag.setCharge(atomToBeReplaced.getCharge());
					terminalAtomOfReplacementFrag.setProtonsExplicitlyAddedOrRemoved(atomToBeReplaced.getProtonsExplicitlyAddedOrRemoved());
				}
				atomToBeReplaced.neutraliseCharge();
			}
		}
	}
	
	/**
	 * Given an atom that is to be replaced by a functional replacement fragment
	 * determines whether this atom has outvalency and if it does removes the outatom from the atom's fragment
	 * and adds an outatom to the replacementFrag
	 * @param atomToBeReplaced
	 * @param replacementFrag
	 */
	private void moveObsoleteOutAtoms(Atom atomToBeReplaced, Fragment replacementFrag){
		if (atomToBeReplaced.getOutValency() >0){//this is not known to occur in well formed IUPAC names but would occur in thioxy (as a suffix)
			List<Atom> replacementAtomList = replacementFrag.getAtomList();
			Fragment origFrag = atomToBeReplaced.getFrag();
			for (int i = origFrag.getOutAtomCount() - 1; i >=0; i--) {
				OutAtom outAtom = origFrag.getOutAtom(i);
				if (atomToBeReplaced.equals(outAtom.getAtom())){
					atomToBeReplaced.getFrag().removeOutAtom(i);
					Atom terminalAtomOfReplacementFrag = replacementAtomList.get(replacementAtomList.size()-1);
					replacementFrag.addOutAtom(terminalAtomOfReplacementFrag, outAtom.getValency(), outAtom.isSetExplicitly());
				}
			}
		}
	}
	
	private void removeAssociatedFunctionalAtom(Atom atomWithFunctionalAtom) throws StructureBuildingException {
		Fragment frag = atomWithFunctionalAtom.getFrag();
		for (int i = frag.getFunctionalAtomCount() - 1; i >=0; i--) {
			FunctionalAtom functionalAtom = frag.getFunctionalAtom(i);
			if (atomWithFunctionalAtom.equals(functionalAtom.getAtom())){
				atomWithFunctionalAtom.getFrag().removeFunctionalAtom(i);
				return;
			}
		}
		throw new StructureBuildingException("OPSIN bug: Unable to find associated functionalAtom");
	}
	

	/**
	 * Returns the subset of oxygenAtoms that possess one of the locants in locantEl
	 * Searches for locant on nearest non suffix atom in case of suffixes
	 * @param locantEl
	 * @param oxygenAtoms
	 * @return
	 */
	private List<Atom> pickOxygensWithAppropriateLocants(Element locantEl, List<Atom> oxygenAtoms) {
		String[] possibleLocants = locantEl.getValue().split(",");
		boolean pLocantSpecialCase = allLocantsP(possibleLocants);
		List<Atom> oxygenWithAppropriateLocants = new ArrayList<Atom>();
		for (Atom atom : oxygenAtoms) {
			List<String> atomlocants = atom.getLocants();
			if (atomlocants.size() > 0) {
				for (String locantVal : possibleLocants) {
					if (atomlocants.contains(locantVal)) {
						 oxygenWithAppropriateLocants.add(atom);
						 break;
					}
				}
			}
			else if (pLocantSpecialCase) {
				for (Atom neighbour : atom.getAtomNeighbours()) {
					if (neighbour.getElement() == ChemEl.P) {
						 oxygenWithAppropriateLocants.add(atom);
						 break;
					}
				}
			}
			else {
				Atom atomWithNumericLocant = OpsinTools.depthFirstSearchForAtomWithNumericLocant(atom);
				if (atomWithNumericLocant != null) {
					List<String> atomWithNumericLocantLocants = atomWithNumericLocant.getLocants();
					for (String locantVal : possibleLocants) {
						if (atomWithNumericLocantLocants.contains(locantVal)) {
							 oxygenWithAppropriateLocants.add(atom);
							 break;
						}
					}
				}
			}
		}
		return oxygenWithAppropriateLocants;
	}

	private boolean allLocantsP(String[] locants) {
		if (locants.length == 0) {
			return false;
		}
		for (String locant : locants) {
			if (!locant.equals("P")) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Returns oxygen atoms in suffixes with functionalAtoms
	 * @param groupToBeModified
	 * @return
	 */
	private List<Atom> findFunctionalOxygenAtomsInApplicableSuffixes(Element groupToBeModified) {
		List<Element> suffixElements =OpsinTools.getNextSiblingsOfType(groupToBeModified, SUFFIX_EL);
		List<Atom> oxygenAtoms = new ArrayList<Atom>();
		for (Element suffix : suffixElements) {
			Fragment suffixFrag = suffix.getFrag();
			if (suffixFrag != null) {//null for non carboxylic acids
				for (int i = 0, l = suffixFrag.getFunctionalAtomCount(); i < l; i++) {
					Atom a = suffixFrag.getFunctionalAtom(i).getAtom();
					if (a.getElement() == ChemEl.O) {
						oxygenAtoms.add(a);
					}
				}
			}
		}
		return oxygenAtoms;
	}

	/**
	 * Returns functional oxygen atoms in groupToBeModified
	 * @param groupToBeModified
	 * @return
	 */
	private List<Atom> findFunctionalOxygenAtomsInGroup(Element groupToBeModified) {
		List<Atom> oxygenAtoms = new ArrayList<Atom>();
		Fragment frag = groupToBeModified.getFrag();
		for (int i = 0, l = frag.getFunctionalAtomCount(); i < l; i++) {
			Atom a = frag.getFunctionalAtom(i).getAtom();
			if (a.getElement() == ChemEl.O){
				oxygenAtoms.add(a);
			}
		}
		return oxygenAtoms;
	}
	
	
	/**
	 * Returns etheric oxygen atoms in groupToBeModified
	 * @param groupToBeModified
	 * @return
	 */
	private List<Atom> findEthericOxygenAtomsInGroup(Element groupToBeModified) {
		List<Atom> oxygenAtoms = new ArrayList<Atom>();
		List<Atom> atomList = groupToBeModified.getFrag().getAtomList();
		for (Atom a: atomList) {
			if (a.getElement() == ChemEl.O && a.getBondCount()==2 && a.getCharge()==0 && a.getIncomingValency()==2){
				oxygenAtoms.add(a);
			}
		}
		return oxygenAtoms;
	}
	
	
	/**
	 * Returns oxygen atoms in suffixes with functionalAtoms or acidStem suffixes or aldehyde suffixes (1979 C-531)
	 * @param groupToBeModified
	 * @return
	 */
	private List<Atom> findOxygenAtomsInApplicableSuffixes(Element groupToBeModified) {
		List<Element> suffixElements =OpsinTools.getNextSiblingsOfType(groupToBeModified, SUFFIX_EL);
		List<Atom> oxygenAtoms = new ArrayList<Atom>();
		for (Element suffix : suffixElements) {
			Fragment suffixFrag = suffix.getFrag();
			if (suffixFrag != null) {//null for non carboxylic acids
				if (suffixFrag.getFunctionalAtomCount() > 0 || groupToBeModified.getAttributeValue(TYPE_ATR).equals(ACIDSTEM_TYPE_VAL) || suffix.getAttributeValue(VALUE_ATR).equals("aldehyde")) {
					List<Atom> atomList = suffixFrag.getAtomList();
					for (Atom a : atomList) {
						if (a.getElement() == ChemEl.O) {
							oxygenAtoms.add(a);
						}
					}
				}
			}
		}
		return oxygenAtoms;
	}
	
	/**
	 * Returns oxygen atoms in groupToBeModified
	 * @param groupToBeModified
	 * @return
	 */
	private List<Atom> findOxygenAtomsInGroup(Element groupToBeModified) {
		List<Atom> oxygenAtoms = new ArrayList<Atom>();
		List<Atom> atomList = groupToBeModified.getFrag().getAtomList();
		for (Atom a : atomList) {
			if (a.getElement() == ChemEl.O){
				oxygenAtoms.add(a);
			}
		}
		return oxygenAtoms;
	}
	

	private void populateTerminalSingleAndDoubleBondedOxygen(List<Atom> atomList, List<Atom> singleBondedOxygen, List<Atom> doubleBondedOxygen) throws StructureBuildingException {
		for (Atom a : atomList) {
			if (a.getElement() == ChemEl.O){//find terminal oxygens
				if (a.getBondCount()==1){
					int incomingValency = a.getIncomingValency();
					if (incomingValency ==2){
						doubleBondedOxygen.add(a);
					}
					else if (incomingValency ==1){
						singleBondedOxygen.add(a);
					}
					else{
						throw new StructureBuildingException("Unexpected bond order to oxygen; excepted 1 or 2 found: " +incomingValency);
					}

				}
			}
		}
	}
}
