package uk.ac.cam.ch.wwmm.opsin;

import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;
import static uk.ac.cam.ch.wwmm.opsin.OpsinTools.*;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;


import nu.xom.Attribute;
import nu.xom.Element;
import nu.xom.Elements;
import nu.xom.Node;

/**
 * Master methods and convenience methods for performing functional replacement
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
			int allowedInputs1 = MATCH_COMMA.split(infixTransformation1).length;
			int allowedInputs2 = MATCH_COMMA.split(infixTransformation2).length;
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
	private final static Pattern matchChalcogen = Pattern.compile("O|S|Se|Te");
	private static final Pattern matchChalcogenReplacement= Pattern.compile("thio|seleno|telluro");
	
	enum PREFIX_REPLACEMENT_TYPE{
		chalcogen,//ambiguous
		halideOrPseudoHalide,//only mean functional replacement when applied to non carboxylic acids
		dedicatedFunctionalReplacementPrefix,//no ambiguity exists
		hydrazono,//ambiguous, only applies to non carboxylic acid
		peroxy//ambiguous, also applies to etheric oxygen
	}

	/**
	 * Applies the effects of amide or hydrazide functional class nomenclature
	 * This must be performed early so that prefix/infix funcional replacement is performed correctly
	 * and so that element symbol locants are assigned appropriately
	 * @param state 
	 * @param finalSubOrRootInWord
	 * @param word
	 * @throws ComponentGenerationException 
	 * @throws StructureBuildingException 
	 */
	static void processAmideOrHydrazideFunctionalClassNomenclature(BuildState state, Element finalSubOrRootInWord, Element word) throws ComponentGenerationException, StructureBuildingException {
		Element wordRule = OpsinTools.getParentWordRule(word);
		WordRule wr = WordRule.valueOf(wordRule.getAttributeValue(WORDRULE_ATR));
		if (wr == WordRule.hydrazide){
			processHydrazideFunctionalClassNomenclature(state, finalSubOrRootInWord, ((Element) XOMTools.getNextSibling(word)));
		}
		else if (wr == WordRule.amide){
			Element parentWordRule = (Element) word.getParent();
			if (parentWordRule.indexOf(word)==0){
				List<Element> amideFullWords = XOMTools.getChildElementsWithTagNameAndAttribute(parentWordRule, WORD_EL, TYPE_ATR, WordType.full.toString());
				amideFullWords.remove(word);
				if (amideFullWords.size()>0){
					//as words are processed from right to left in cases like phosphoric acid tri(ethylamide) this will be phosphoric acid ethylamide ethylamide ethylamide
					for (Element amideWord : amideFullWords) {
						processAmideFunctionalClassNomenclatureFullWord(state, finalSubOrRootInWord, amideWord);
					}
				}
				else if (parentWordRule.getChildElements().size()==2) {
					processAmideFunctionalClassNomenclatureFunctionalWord(state, finalSubOrRootInWord, ((Element) XOMTools.getNextSibling(word)));
				}
				else{
					throw new ComponentGenerationException("OPSIN bug: problem with amide word rule");
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
	 * @param state
	 * @param groups
	 * @param substituents
	 * @return boolean: has any functional replacement occurred
	 * @throws StructureBuildingException
	 * @throws ComponentGenerationException 
	 */
	static boolean processPrefixFunctionalReplacementNomenclature(BuildState state, List<Element> groups, List<Element> substituents) throws StructureBuildingException, ComponentGenerationException {
		int originalNumberOfGroups = groups.size();
		for (int i = originalNumberOfGroups-1; i >=0; i--) {
			Element group =groups.get(i);
			String groupValue = group.getValue();
			PREFIX_REPLACEMENT_TYPE replacementType = null;
			if (matchChalcogenReplacement.matcher(groupValue).matches()){
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
			if (replacementType!=null){
				//need to check whether this is an instance of functional replacement by checking the substituent/root it is applying to
				Element substituent =(Element) group.getParent();
				Element nextSubOrBracket = (Element) XOMTools.getNextSibling(substituent);
				if (nextSubOrBracket!=null && (nextSubOrBracket.getLocalName().equals(ROOT_EL) || nextSubOrBracket.getLocalName().equals(SUBSTITUENT_EL))){
					Element groupToBeModified = nextSubOrBracket.getFirstChildElement(GROUP_EL);
					if (XOMTools.getPreviousSibling(groupToBeModified)!=null){
						if (replacementType  == PREFIX_REPLACEMENT_TYPE.dedicatedFunctionalReplacementPrefix){
							throw new ComponentGenerationException("dedicated Functional Replacement Prefix used in an inappropriate position :" + groupValue);
						}
						continue;//not 2,2'-thiodipyran
					}
					Element locantEl =null;//null unless a locant that agrees with the multiplier is present
					Element multiplierEl =null;
					int numberOfAtomsToReplace =1;//the number of atoms to be functionally replaced, modified by a multiplier e.g. dithio
					Element possibleMultiplier = (Element) XOMTools.getPreviousSibling(group);
					if (possibleMultiplier !=null){
						Element possibleLocant;
						if (possibleMultiplier.getLocalName().equals(MULTIPLIER_EL)){
							numberOfAtomsToReplace =Integer.valueOf(possibleMultiplier.getAttributeValue(VALUE_ATR));
							possibleLocant = (Element) XOMTools.getPreviousSibling(possibleMultiplier);
							multiplierEl = possibleMultiplier;
						}
						else{
							possibleLocant = possibleMultiplier;
						}
						if (possibleLocant !=null && possibleLocant.getLocalName().equals(LOCANT_EL) && possibleLocant.getAttribute(TYPE_ATR)==null) {
							int numberOfLocants = MATCH_COMMA.split(possibleLocant.getValue()).length;
							if (numberOfLocants == numberOfAtomsToReplace){//locants and number of replacements agree
								locantEl = possibleLocant;
							}
							else if (numberOfAtomsToReplace >1){//doesn't look like prefix functional replacement
								if (replacementType  == PREFIX_REPLACEMENT_TYPE.dedicatedFunctionalReplacementPrefix){
									throw new ComponentGenerationException("dedicated Functional Replacement Prefix used in an inappropriate position :" + groupValue);
								}
								continue;
							}
						}
					}

					int oxygenReplaced;
					if (replacementType == PREFIX_REPLACEMENT_TYPE.chalcogen){
						oxygenReplaced = performChalcogenFunctionalReplacement(state, groupToBeModified, locantEl, numberOfAtomsToReplace, group.getAttributeValue(VALUE_ATR));
					}
					else if (replacementType == PREFIX_REPLACEMENT_TYPE.peroxy){
						if (nextSubOrBracket.getLocalName().equals(SUBSTITUENT_EL)){
							continue;
						}
						oxygenReplaced = performPeroxyFunctionalReplacement(state, groupToBeModified, locantEl, numberOfAtomsToReplace);
					}
					else if (replacementType == PREFIX_REPLACEMENT_TYPE.dedicatedFunctionalReplacementPrefix){
						if (!groupToBeModified.getAttributeValue(TYPE_ATR).equals(NONCARBOXYLICACID_TYPE_VAL)){
							throw new ComponentGenerationException("dedicated Functional Replacement Prefix used in an inappropriate position :" + groupValue);
						}
						oxygenReplaced = performFunctionalReplacementOnAcid(state, groupToBeModified, locantEl, numberOfAtomsToReplace, group.getAttributeValue(VALUE_ATR));
						if (oxygenReplaced==0){
							throw new ComponentGenerationException("dedicated Functional Replacement Prefix used in an inappropriate position :" + groupValue);
						}
					}
					else if (replacementType == PREFIX_REPLACEMENT_TYPE.hydrazono || replacementType == PREFIX_REPLACEMENT_TYPE.halideOrPseudoHalide){
						boolean appropriate = checkGroupIsAnAppropriateNonCarboxylicAcid(state, groupToBeModified, state.xmlFragmentMap.get(group).getOutAtom(0).getValency());
						if (!appropriate){
							continue;
						}
						oxygenReplaced = performFunctionalReplacementOnAcid(state, groupToBeModified, locantEl, numberOfAtomsToReplace, group.getAttributeValue(VALUE_ATR));
					}
					else{
						throw new StructureBuildingException("OPSIN bug: Unexpected prefix replacement type");
					}
					if (oxygenReplaced>0){
						state.fragManager.removeFragment(state.xmlFragmentMap.get(group));
						substituent.removeChild(group);
						groups.remove(group);
						Elements remainingChildren =substituent.getChildElements();//there may be a locant that should be moved
						for (int j = remainingChildren.size()-1; j>=0; j--){
							Node child =substituent.getChild(j);
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
	
	
	/*
	 * 
	 */
	
	/**
	 * Performs functional replacement using infixes e.g. thio in ethanthioic acid replaces an O with S
	 * @param state
	 * @param suffixFragments May be modified if a multiplier is determined to mean multiplication of a suffix, usually untouched
	 * @param suffixes The suffix elements  May be modified if a multiplier is determined to mean multiplication of a suffix, usually untouched
	 * @throws StructureBuildingException
	 * @throws ComponentGenerationException
	 */
	static void processInfixFunctionalReplacementNomenclature(BuildState state, List<Element> suffixes, List<Fragment> suffixFragments) throws StructureBuildingException, ComponentGenerationException {
		for (int i = 0; i < suffixes.size(); i++) {
			Element suffix = suffixes.get(i);
			if (suffix.getAttribute(INFIX_ATR)!=null){
				Fragment fragToApplyInfixTo = state.xmlFragmentMap.get(suffix);
				Element possibleAcidGroup = XOMTools.getPreviousSiblingIgnoringCertainElements(suffix, new String[]{MULTIPLIER_EL, INFIX_EL, SUFFIX_EL});
				if (possibleAcidGroup !=null && possibleAcidGroup.getLocalName().equals(GROUP_EL) && 
						(possibleAcidGroup.getAttributeValue(TYPE_ATR).equals(NONCARBOXYLICACID_TYPE_VAL)|| possibleAcidGroup.getAttributeValue(TYPE_ATR).equals(CHALCOGENACIDSTEM_TYPE_VAL))){
					fragToApplyInfixTo = state.xmlFragmentMap.get(possibleAcidGroup);
				}
				if (fragToApplyInfixTo ==null){
					throw new ComponentGenerationException("infix has erroneously been assigned to a suffix which does not correspond to a suffix fragment. suffix: " + suffix.getValue());
				}
				//e.g. =O:S,-O:S (which indicates replacing either a double or single bonded oxygen with S)
				//This is semicolon delimited for each infix
				List<String> infixTransformations = StringTools.arrayToList(MATCH_SEMICOLON.split(suffix.getAttributeValue(INFIX_ATR)));

				List<Atom> atomList =fragToApplyInfixTo.getAtomList();
				LinkedList<Atom> singleBondedOxygen =new LinkedList<Atom>();
				LinkedList<Atom> doubleBondedOxygen =new LinkedList<Atom>();
				populateTerminalSingleAndDoubleBondedOxygen(atomList, singleBondedOxygen, doubleBondedOxygen);
				int oxygenAvailable = singleBondedOxygen.size() +doubleBondedOxygen.size();

				/*
				 * Modifies suffixes, suffixFragments, suffix and infixTransformations as appropriate
				 */
				disambiguateMultipliedInfixMeaning(state, suffixes, suffixFragments, suffix, fragToApplyInfixTo, infixTransformations, oxygenAvailable);

				/*
				 * Sort infixTransformations so more specific transformations are performed first
				 * e.g. ethanthioimidic acid-->ethanimidthioic acid as imid can only apply to the double bonded oxygen
				 */
				Collections.sort(infixTransformations, new SortInfixTransformations());

				for (String infixTransformation : infixTransformations) {
					String[] transformationArray = MATCH_COLON.split(infixTransformation);
					if (transformationArray.length !=2){
						throw new StructureBuildingException("Atom to be replaced and replacement not specified correctly in infix: " + infixTransformation);
					}
					String[] transformations = MATCH_COMMA.split(transformationArray[0]);
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

					Set<Atom> ambiguousElementAtoms = new HashSet<Atom>();
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
					Fragment replacementFrag =state.fragManager.buildSMILES(replacementSMILES, SUFFIX_TYPE_VAL, NONE_LABELS_VAL);
					if (replacementFrag.getOutAtoms().size()>0){
						replacementFrag.removeOutAtom(0);
					}
					Atom atomThatWillReplaceOxyen =replacementFrag.getFirstAtom();
					if (replacementFrag.getAtomList().size()==1 && matchChalcogen.matcher(atomThatWillReplaceOxyen.getElement()).matches()){
						atomThatWillReplaceOxyen.setCharge(atomToUse.getCharge());
						atomThatWillReplaceOxyen.setProtonsExplicitlyAddedOrRemoved(atomToUse.getProtonsExplicitlyAddedOrRemoved());
					}
					removeOrMoveObsoleteFunctionalAtoms(atomToUse, replacementFrag);//also will move charge if necessary
					if (nitrido){
						atomToUse.getFirstBond().setOrder(3);
						state.fragManager.removeAtomAndAssociatedBonds(singleBondedOxygen.removeFirst());
					}
					state.fragManager.incorporateFragment(replacementFrag, atomToUse.getFrag());
					state.fragManager.replaceAtomWithAnotherAtomPreservingConnectivity(atomToUse, atomThatWillReplaceOxyen);
					if (infixAssignmentAmbiguous){
						ambiguousElementAtoms.add(atomThatWillReplaceOxyen);
						if (atomThatWillReplaceOxyen.getProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT)!=null){
							ambiguousElementAtoms.addAll(atomThatWillReplaceOxyen.getProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT));
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
	 * Replaces the appropriate number of functional oxygen atoms with hydrazide fragments
	 * @param state
	 * @param acidContainingRoot
	 * @param functionalWord
	 * @throws ComponentGenerationException
	 * @throws StructureBuildingException
	 */
	private static void processHydrazideFunctionalClassNomenclature(BuildState state, Element acidContainingRoot, Element functionalWord) throws ComponentGenerationException, StructureBuildingException {
		if (functionalWord !=null && functionalWord.getAttributeValue(TYPE_ATR).equals(WordType.functionalTerm.toString())){
			Element functionalTerm = functionalWord.getFirstChildElement(FUNCTIONALTERM_EL);
			if (functionalTerm ==null){
				throw new ComponentGenerationException("OPSIN bug: functionalTerm word not found where one was expected for hydrazide wordRule");
			}
			Element hydrazideGroup = functionalTerm.getFirstChildElement(FUNCTIONALGROUP_EL);
			Fragment hydrazide = ComponentProcessor.resolveGroup(state, hydrazideGroup);
			Element possibleMultiplier = (Element) XOMTools.getPreviousSibling(hydrazideGroup);
			int hydrazides =1;
			if (possibleMultiplier!=null){
				if (!possibleMultiplier.getLocalName().equals(MULTIPLIER_EL)){
					throw new ComponentGenerationException("OPSIN bug: non multiplier found where only a multiplier was expected in hydrazide wordRule");
				}
				hydrazides = Integer.parseInt(possibleMultiplier.getAttributeValue(VALUE_ATR));
				possibleMultiplier.detach();
			}
			if (functionalTerm.getChildElements().size()!=1){
				throw new ComponentGenerationException("Unexpected qualifier to hydrazide functionalTerm");
			}
			
			Element groupToBeModified = acidContainingRoot.getFirstChildElement(GROUP_EL);
			List<Atom> oxygenAtoms = findFunctionalOxygenAtomsInApplicableSuffixes(state, groupToBeModified);
			if (oxygenAtoms.size()==0){
				oxygenAtoms = findFunctionalOxygenAtomsInGroup(state, groupToBeModified);
			}
			if (oxygenAtoms.size()==0){
				List<Element> conjunctiveSuffixElements =XOMTools.getNextSiblingsOfType(groupToBeModified, CONJUNCTIVESUFFIXGROUP_EL);
				for (Element conjunctiveSuffixElement : conjunctiveSuffixElements) {
					oxygenAtoms.addAll(findFunctionalOxygenAtomsInGroup(state, conjunctiveSuffixElement));
				}
			}
			if (hydrazides > oxygenAtoms.size()){
				throw new ComponentGenerationException("Insufficient oxygen to replace with hydrazides in " + acidContainingRoot.getFirstChildElement(GROUP_EL).getValue());
			}
			
			Fragment acidFragment = state.xmlFragmentMap.get(groupToBeModified);
			if (acidFragment.hasLocant("2")){//prefer numeric locants on group to those of hydrazide
				for (Atom atom : hydrazide.getAtomList()) {
					atom.clearLocants();
				}
			}
			state.fragManager.replaceAtomWithAnotherAtomPreservingConnectivity(oxygenAtoms.get(0), hydrazide.getFirstAtom());
			removeAssociatedFunctionalAtom(oxygenAtoms.get(0));
			for (int i = 1; i < hydrazides; i++) {
				Fragment clonedHydrazide = state.fragManager.copyAndRelabelFragment(hydrazide, i);
				state.fragManager.replaceAtomWithAnotherAtomPreservingConnectivity(oxygenAtoms.get(i), clonedHydrazide.getFirstAtom());
				state.fragManager.incorporateFragment(clonedHydrazide, oxygenAtoms.get(i).getFrag());
				removeAssociatedFunctionalAtom(oxygenAtoms.get(i));
			}
			state.fragManager.incorporateFragment(hydrazide, oxygenAtoms.get(0).getFrag());
		}
		else{
			throw new ComponentGenerationException("hydrazide word not found where expected, bug?");
		}
	}
	
	private static void processAmideFunctionalClassNomenclatureFullWord(BuildState state, Element acidContainingRoot, Element amideWord) throws ComponentGenerationException, StructureBuildingException {
		Element amideGroup = StructureBuildingMethods.findRightMostGroupInBracket(amideWord);
		if (amideGroup ==null){
			throw new ComponentGenerationException("OPSIN bug: amide group not found where one was expected for amide wordRule");
		}
		Fragment amide = state.xmlFragmentMap.get(amideGroup);
		if (((Element)amideGroup.getParent()).getChildElements().size()!=1){
			throw new ComponentGenerationException("Unexpected qualifier to amide");
		}
		
		Element groupToBeModified = acidContainingRoot.getFirstChildElement(GROUP_EL);
		List<Atom> oxygenAtoms = findFunctionalOxygenAtomsInApplicableSuffixes(state, groupToBeModified);
		if (oxygenAtoms.size()==0){
			oxygenAtoms = findFunctionalOxygenAtomsInGroup(state, groupToBeModified);
		}
		if (oxygenAtoms.size()==0){
			List<Element> conjunctiveSuffixElements =XOMTools.getNextSiblingsOfType(groupToBeModified, CONJUNCTIVESUFFIXGROUP_EL);
			for (Element conjunctiveSuffixElement : conjunctiveSuffixElements) {
				oxygenAtoms.addAll(findFunctionalOxygenAtomsInGroup(state, conjunctiveSuffixElement));
			}
		}
		if (oxygenAtoms.size()<1){
			throw new ComponentGenerationException("Insufficient oxygen to replace with amides in " + acidContainingRoot.getFirstChildElement(GROUP_EL).getValue());
		}
		if (amide.getAtomList().size()!=1){
			throw new ComponentGenerationException("OPSIN bug: amide not found where expected");
		}
		Atom amideNitrogen = amide.getFirstAtom();
		amideNitrogen.neutraliseCharge();
		amideNitrogen.clearLocants();
		amide.addMappingToAtomLocantMap("N", amideNitrogen);
		state.fragManager.replaceAtomWithAnotherAtomPreservingConnectivity(oxygenAtoms.get(0), amide.getFirstAtom());
		state.fragManager.incorporateFragment(amide, oxygenAtoms.get(0).getFrag());
		removeAssociatedFunctionalAtom(oxygenAtoms.get(0));
	}

	/**
	 * Replaces the appropriate number of functional oxygen atoms with nitrogen atoms
	 * @param state
	 * @param acidContainingRoot
	 * @param functionalWord
	 * @throws ComponentGenerationException
	 * @throws StructureBuildingException
	 */
	private static void processAmideFunctionalClassNomenclatureFunctionalWord(BuildState state, Element acidContainingRoot, Element functionalWord) throws ComponentGenerationException, StructureBuildingException {
		if (functionalWord !=null && functionalWord.getAttributeValue(TYPE_ATR).equals(WordType.functionalTerm.toString())){
			Element functionalTerm = functionalWord.getFirstChildElement(FUNCTIONALTERM_EL);
			if (functionalTerm ==null){
				throw new ComponentGenerationException("OPSIN bug: functionalTerm word not found where one was expected for amide wordRule");
			}
			Element amideGroup = functionalTerm.getFirstChildElement(FUNCTIONALGROUP_EL);
			Element possibleMultiplier = (Element) XOMTools.getPreviousSibling(amideGroup);
			int numbrOfAmidesToForm =1;
			if (possibleMultiplier!=null){
				if (!possibleMultiplier.getLocalName().equals(MULTIPLIER_EL)){
					throw new ComponentGenerationException("OPSIN bug: non multiplier found where only a multiplier was expected in amide wordRule");
				}
				numbrOfAmidesToForm = Integer.parseInt(possibleMultiplier.getAttributeValue(VALUE_ATR));
				possibleMultiplier.detach();
			}
			if (functionalTerm.getChildElements().size()!=1){
				throw new ComponentGenerationException("Unexpected qualifier to amide functionalTerm");
			}
			
			Element groupToBeModified = acidContainingRoot.getFirstChildElement(GROUP_EL);
			List<Atom> oxygenAtoms = findFunctionalOxygenAtomsInApplicableSuffixes(state, groupToBeModified);
			if (oxygenAtoms.size()==0){
				oxygenAtoms = findFunctionalOxygenAtomsInGroup(state, groupToBeModified);
			}
			if (oxygenAtoms.size()==0){
				List<Element> conjunctiveSuffixElements =XOMTools.getNextSiblingsOfType(groupToBeModified, CONJUNCTIVESUFFIXGROUP_EL);
				for (Element conjunctiveSuffixElement : conjunctiveSuffixElements) {
					oxygenAtoms.addAll(findFunctionalOxygenAtomsInGroup(state, conjunctiveSuffixElement));
				}
			}
			if (numbrOfAmidesToForm > oxygenAtoms.size()){
				throw new ComponentGenerationException("Insufficient oxygen to replace with nitrogen in " + acidContainingRoot.getFirstChildElement(GROUP_EL).getValue());
			}
			for (int i = 0; i < numbrOfAmidesToForm; i++) {
				oxygenAtoms.get(i).setElement("N");
				removeAssociatedFunctionalAtom(oxygenAtoms.get(i));
			}
		}
		else{
			throw new ComponentGenerationException("amide word not found where expected, bug?");
		}
	}
	

	/*
	 * Prefix functional replacement nomenclature
	 */

	
	private static boolean checkGroupIsAnAppropriateNonCarboxylicAcid(BuildState state, Element groupToBeModified, int maxSubstitutableHydrogenOnOneAtom) throws StructureBuildingException {
		if (!groupToBeModified.getAttributeValue(TYPE_ATR).equals(NONCARBOXYLICACID_TYPE_VAL)){
			return false;//hydrazono replacement only applies to non carboxylic acids e.g. hydrazonooxalic acid
		}
		return calculateMaxSubstitutableHydrogenAttachedToAnAcidCentre(state, groupToBeModified) < maxSubstitutableHydrogenOnOneAtom;
	}
	
	private static int calculateMaxSubstitutableHydrogenAttachedToAnAcidCentre(BuildState state, Element group) {
		Fragment fragToBeModified = state.xmlFragmentMap.get(group);
		if (group.getAttribute(SUFFIXAPPLIESTO_ATR)!=null){
			int maxSubstitutableHydrogen =0;
			String[] atomIndices = MATCH_COMMA.split(group.getAttributeValue(SUFFIXAPPLIESTO_ATR));
			List<Atom> atomList = fragToBeModified.getAtomList();
			for (String atomIndice : atomIndices) {
				int substitutableHydrogen = StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(atomList.get(Integer.parseInt(atomIndice)-1));
				if (substitutableHydrogen > maxSubstitutableHydrogen){
					maxSubstitutableHydrogen = substitutableHydrogen;
				}
			}
			return maxSubstitutableHydrogen;
		}
		else{
			return StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(fragToBeModified.getFirstAtom());
		}
	}


	/**
	 * Performs replacement of oxygen atoms by chalogen atoms
	 * If this is ambiguous e.g. thioacetate then Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT is populated
	 * @param state
	 * @param groupToBeModified
	 * @param locantEl
	 * @param numberOfAtomsToReplace
	 * @param replacementSmiles
	 * @return
	 * @throws StructureBuildingException
	 */
	private static int performChalcogenFunctionalReplacement(BuildState state, Element groupToBeModified, Element locantEl, int numberOfAtomsToReplace, String replacementSmiles) throws StructureBuildingException {
		List<Atom> oxygenAtoms = findOxygenAtomsInApplicableSuffixes(state, groupToBeModified);
		if (oxygenAtoms.size()==0){
			oxygenAtoms = findOxygenAtomsInGroup(state, groupToBeModified);
		}
		if (locantEl !=null){//locants are used to indicate replacement on trivial groups
			List<Atom> oxygenWithAppropriateLocants = pickOxygensWithAppropriateLocants(locantEl, oxygenAtoms);
			if(oxygenWithAppropriateLocants.size() < numberOfAtomsToReplace){
				numberOfAtomsToReplace =1;
				//e.g. -1-thioureidomethyl
			}
			else{
				locantEl.detach();
				oxygenAtoms = oxygenWithAppropriateLocants;
			}
		}

		List<Atom> doubleBondedOxygen = new ArrayList<Atom>();
		List<Atom> singleBondedOxygen = new ArrayList<Atom>();
		List<Atom> ethericOxygen = new ArrayList<Atom>();
		for (Atom oxygen : oxygenAtoms) {
			int incomingValency = oxygen.getIncomingValency();
			int bondCount = oxygen.getBonds().size();
			if (bondCount==1 && incomingValency==2){
				doubleBondedOxygen.add(oxygen);
			}
			else if (bondCount==1 && incomingValency==1){
				singleBondedOxygen.add(oxygen);
			}
			else if (bondCount==2 && incomingValency==2){
				ethericOxygen.add(oxygen);
			}
		}
		List<Atom> replaceableAtoms = new LinkedList<Atom>();
		replaceableAtoms.addAll(doubleBondedOxygen);
		replaceableAtoms.addAll(singleBondedOxygen);
		replaceableAtoms.addAll(ethericOxygen);
		int totalOxygen = replaceableAtoms.size();
		if (numberOfAtomsToReplace >1){
			if (totalOxygen < numberOfAtomsToReplace){
				numberOfAtomsToReplace=1;
			}
		}

		int atomsReplaced =0;
		if (totalOxygen >=numberOfAtomsToReplace){//check that there atleast as many oxygens as requested replacements
			boolean prefixAssignmentAmbiguous =false;
			Set<Atom> ambiguousElementAtoms = new HashSet<Atom>();
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
	 * @param state
	 * @param groupToBeModified
	 * @param locantEl
	 * @param numberOfAtomsToReplace
	 * @return
	 * @throws StructureBuildingException
	 */
	private static int performPeroxyFunctionalReplacement(BuildState state, Element groupToBeModified, Element locantEl, int numberOfAtomsToReplace) throws StructureBuildingException {
		List<Atom> oxygenAtoms = findFunctionalOxygenAtomsInApplicableSuffixes(state, groupToBeModified);
		if (oxygenAtoms.size()==0){
			oxygenAtoms = findEthericOxygenAtomsInGroup(state, groupToBeModified);
			oxygenAtoms.addAll(findFunctionalOxygenAtomsInGroup(state, groupToBeModified));
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
				if (oxygenToReplace.getBonds().size()==2){//etheric oxygen
					Fragment newOxygen = state.fragManager.buildSMILES("O", SUFFIX_TYPE_VAL, NONE_LABELS_VAL);
					Bond bondToRemove = oxygenToReplace.getFirstBond();
					Atom atomToAttachTo = bondToRemove.getFromAtom() == oxygenToReplace ?  bondToRemove.getToAtom() :  bondToRemove.getFromAtom();
					state.fragManager.createBond(atomToAttachTo, newOxygen.getFirstAtom(), 1);
					state.fragManager.createBond(newOxygen.getFirstAtom(), oxygenToReplace, 1);
					state.fragManager.removeBond(bondToRemove);
					state.fragManager.incorporateFragment(newOxygen, state.xmlFragmentMap.get(groupToBeModified));
				}
				else{
					Fragment replacementFrag = state.fragManager.buildSMILES("OO", SUFFIX_TYPE_VAL, NONE_LABELS_VAL);
					removeOrMoveObsoleteFunctionalAtoms(oxygenToReplace, replacementFrag);
					state.fragManager.replaceAtomWithAnotherAtomPreservingConnectivity(oxygenToReplace, replacementFrag.getFirstAtom());
					state.fragManager.incorporateFragment(replacementFrag, state.xmlFragmentMap.get(groupToBeModified));
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
	 * @param state
	 * @param groupToBeModified
	 * @param locantEl
	 * @param numberOfAtomsToReplace
	 * @param replacementSmiles
     * @return
	 * @throws StructureBuildingException
	 */
	private static int performFunctionalReplacementOnAcid(BuildState state, Element groupToBeModified, Element locantEl, int numberOfAtomsToReplace, String replacementSmiles) throws StructureBuildingException {
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
		List<Atom> oxygenAtoms = findOxygenAtomsInApplicableSuffixes(state, groupToBeModified);
		if (oxygenAtoms.size()==0){
			oxygenAtoms = findOxygenAtomsInGroup(state, groupToBeModified);
		}
		if (locantEl !=null){//locants are used to indicate replacement on trivial groups
			List<Atom> oxygenWithAppropriateLocants = pickOxygensWithAppropriateLocants(locantEl, oxygenAtoms);
			LinkedList<Atom> singleBondedOxygen =new LinkedList<Atom>();
			LinkedList<Atom> terminalDoubleBondedOxygen =new LinkedList<Atom>();
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
		LinkedList<Atom> singleBondedOxygen =new LinkedList<Atom>();
		LinkedList<Atom> terminalDoubleBondedOxygen =new LinkedList<Atom>();
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
					Fragment replacementFrag = state.fragManager.buildSMILES(replacementSmiles, atomToReplace.getFrag().getType(), NONE_LABELS_VAL);
					if (outValency ==3){//special case for nitrido
						atomToReplace.getFirstBond().setOrder(3);
						state.fragManager.removeAtomAndAssociatedBonds(singleBondedOxygen.removeFirst());
					}
					state.fragManager.replaceAtomWithAnotherAtomPreservingConnectivity(atomToReplace, replacementFrag.getFirstAtom());
					if (outValency ==1){
						removeOrMoveObsoleteFunctionalAtoms(atomToReplace, replacementFrag);
					}
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
	 * @param state
	 * @param suffixes
	 * @param suffixFragments
	 * @param suffix
	 * @param suffixFrag
	 * @param infixTransformations
	 * @param oxygenAvailable
	 * @throws ComponentGenerationException
	 * @throws StructureBuildingException
	 */
	private static void disambiguateMultipliedInfixMeaning(BuildState state,List<Element> suffixes,
			List<Fragment> suffixFragments,Element suffix, Fragment suffixFrag, List<String> infixTransformations, int oxygenAvailable)
			throws ComponentGenerationException, StructureBuildingException {
		Element possibleInfix =(Element) XOMTools.getPreviousSibling(suffix);
		if (possibleInfix.getLocalName().equals(INFIX_EL)){//the infix is only left when there was ambiguity
			Element possibleMultiplier =(Element) XOMTools.getPreviousSibling(possibleInfix);
			if (possibleMultiplier.getLocalName().equals(MULTIPLIER_EL)){
				int multiplierValue =Integer.parseInt(possibleMultiplier.getAttributeValue(VALUE_ATR));
				if (infixTransformations.size() + multiplierValue-1 <=oxygenAvailable){//multiplier means multiply the infix e.g. butandithiate
					for (int j = 1; j < multiplierValue; j++) {
						infixTransformations.add(0, infixTransformations.get(0));
					}
				}
				else{
					Element possibleLocant =(Element)XOMTools.getPreviousSibling(possibleMultiplier);
					String[] locants = null;
					if (possibleLocant.getLocalName().equals(LOCANT_EL)) {
						locants = MATCH_COMMA.split(possibleLocant.getValue());
					}
					if (locants !=null){
						if (locants.length!=multiplierValue){
							throw new ComponentGenerationException("Multiplier/locant disagreement when multiplying infixed suffix");
						}
					    suffix.addAttribute(new Attribute(LOCANT_ATR, locants[0]));
					}
					suffix.addAttribute(new Attribute(MULTIPLIED_ATR, "multiplied"));
					for (int j = 1; j < multiplierValue; j++) {//multiplier means multiply the infixed suffix e.g. butandithione
						Element newSuffix =new Element(suffix);
						Fragment newSuffixFrag =state.fragManager.copyFragment(suffixFrag);
						state.xmlFragmentMap.put(newSuffix, newSuffixFrag);
						suffixFragments.add(newSuffixFrag);
						XOMTools.insertAfter(suffix, newSuffix);
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
	private static void removeOrMoveObsoleteFunctionalAtoms(Atom atomToBeReplaced, Fragment replacementFrag){
		List<Atom> replacementAtomList = replacementFrag.getAtomList();
		List<FunctionalAtom> functionalAtoms = atomToBeReplaced.getFrag().getFunctionalAtoms();
		for (int j = functionalAtoms.size()-1; j >=0; j--) {
			FunctionalAtom functionalAtom = functionalAtoms.get(j);
			if (atomToBeReplaced.equals(functionalAtom.getAtom())){
				atomToBeReplaced.getFrag().removeFunctionalAtom(j);
				Atom terminalAtomOfReplacementFrag = replacementAtomList.get(replacementAtomList.size()-1);
				if ((terminalAtomOfReplacementFrag.getIncomingValency() ==1 || replacementAtomList.size()==1)&& matchChalcogen.matcher(terminalAtomOfReplacementFrag.getElement()).matches()){
					replacementFrag.addFunctionalAtom(terminalAtomOfReplacementFrag);
					terminalAtomOfReplacementFrag.setCharge(atomToBeReplaced.getCharge());
					terminalAtomOfReplacementFrag.setProtonsExplicitlyAddedOrRemoved(atomToBeReplaced.getProtonsExplicitlyAddedOrRemoved());
				}
				atomToBeReplaced.neutraliseCharge();
			}
		}
	}
	
	private static void removeAssociatedFunctionalAtom(Atom atomWithFunctionalAtom) throws ComponentGenerationException {
		List<FunctionalAtom> functionalAtoms = atomWithFunctionalAtom.getFrag().getFunctionalAtoms();
		for (int j = functionalAtoms.size()-1; j >=0; j--) {
			FunctionalAtom functionalAtom = functionalAtoms.get(j);
			if (atomWithFunctionalAtom.equals(functionalAtom.getAtom())){
				atomWithFunctionalAtom.getFrag().removeFunctionalAtom(j);
				return;
			}
		}
		throw new ComponentGenerationException("OPSIN bug: Unable to find associated functionalAtom");
	}
	

	/**
	 * Returns the subset of oxygenAtoms that possess one of the locants in locantEl
	 * Searches for locant on nearest non suffix atom in case of suffixes
	 * @param locantEl
	 * @param oxygenAtoms
	 * @return
	 */
	private static List<Atom> pickOxygensWithAppropriateLocants(Element locantEl, List<Atom> oxygenAtoms) {
		String[] possibleLocants = MATCH_COMMA.split(locantEl.getValue());
		List<Atom> oxygenWithAppropriateLocants = new ArrayList<Atom>();
		for (Atom atom : oxygenAtoms) {
			List<String> atomlocants = atom.getLocants();
			if (atomlocants.size()>0){
				for (String locantVal : possibleLocants) {
					if (atomlocants.contains(locantVal)){
						 oxygenWithAppropriateLocants.add(atom);
						 break;
					}
				}
			}
			else{
				Atom atomWithNumericLocant = OpsinTools.depthFirstSearchForAtomWithNumericLocant(atom);
				if (atomWithNumericLocant!=null){
					List<String> atomWithNumericLocantLocants = atomWithNumericLocant.getLocants();
					for (String locantVal : possibleLocants) {
						if (atomWithNumericLocantLocants.contains(locantVal)){
							 oxygenWithAppropriateLocants.add(atom);
							 break;
						}
					}
				}
			}
		}
		return oxygenWithAppropriateLocants;
	}

	/**
	 * Returns oxygen atoms in suffixes with functionalAtoms
	 * @param state
	 * @param groupToBeModified
	 * @return
	 */
	private static List<Atom> findFunctionalOxygenAtomsInApplicableSuffixes(BuildState state, Element groupToBeModified) {
		List<Element> suffixElements =XOMTools.getNextSiblingsOfType(groupToBeModified, SUFFIX_EL);
		List<Atom> oxygenAtoms = new ArrayList<Atom>();
        for (Element suffix : suffixElements) {
            Fragment suffixFrag = state.xmlFragmentMap.get(suffix);
            if (suffixFrag != null) {//null for non carboxylic acids
                List<FunctionalAtom> functionalAtoms = suffixFrag.getFunctionalAtoms();
                for (FunctionalAtom funcA : functionalAtoms) {
                    Atom a = funcA.getAtom();
                    if (a.getElement().equals("O")) {
                        oxygenAtoms.add(funcA.getAtom());
                    }
                }
            }
        }
		return oxygenAtoms;
	}
	
	/**
	 * Returns functional oxygen atoms in groupToBeModified
	 * @param state
	 * @param groupToBeModified
	 * @return
	 */
	private static List<Atom> findFunctionalOxygenAtomsInGroup(BuildState state, Element groupToBeModified) {
		List<Atom> oxygenAtoms = new ArrayList<Atom>();
		List<FunctionalAtom> functionalAtoms = state.xmlFragmentMap.get(groupToBeModified).getFunctionalAtoms();
		for (FunctionalAtom funcA: functionalAtoms) {
			Atom a = funcA.getAtom();
			if (a.getElement().equals("O")){
				oxygenAtoms.add(funcA.getAtom());
			}
		}
		return oxygenAtoms;
	}
	
	
	/**
	 * Returns etheric oxygen atoms in groupToBeModified
	 * @param state
	 * @param groupToBeModified
	 * @return
	 */
	private static List<Atom> findEthericOxygenAtomsInGroup(BuildState state, Element groupToBeModified) {
		List<Atom> oxygenAtoms = new ArrayList<Atom>();
		List<Atom> atomList = state.xmlFragmentMap.get(groupToBeModified).getAtomList();
		for (Atom a: atomList) {
			if (a.getElement().equals("O") && a.getBonds().size()==2 && a.getCharge()==0 && a.getIncomingValency()==2){
				oxygenAtoms.add(a);
			}
		}
		return oxygenAtoms;
	}
	
	
	/**
	 * Returns oxygen atoms in suffixes with functionalAtoms or acidStem suffixes or aldehyde suffixes (1979 C-531)
	 * @param state
	 * @param groupToBeModified
	 * @return
	 */
	private static List<Atom> findOxygenAtomsInApplicableSuffixes(BuildState state, Element groupToBeModified) {
		List<Element> suffixElements =XOMTools.getNextSiblingsOfType(groupToBeModified, SUFFIX_EL);
		List<Atom> oxygenAtoms = new ArrayList<Atom>();
        for (Element suffix : suffixElements) {
            Fragment suffixFrag = state.xmlFragmentMap.get(suffix);
            if (suffixFrag != null) {//null for non carboxylic acids
                if (suffixFrag.getFunctionalAtoms().size() > 0 || groupToBeModified.getAttributeValue(TYPE_ATR).equals(ACIDSTEM_TYPE_VAL) || suffix.getAttributeValue(VALUE_ATR).equals("aldehyde")) {
                    List<Atom> atomList = suffixFrag.getAtomList();
                    for (Atom a : atomList) {
                        if (a.getElement().equals("O")) {
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
	 * @param state
	 * @param groupToBeModified
	 * @return
	 */
	private static List<Atom> findOxygenAtomsInGroup(BuildState state, Element groupToBeModified) {
		List<Atom> oxygenAtoms = new ArrayList<Atom>();
		List<Atom> atomList = state.xmlFragmentMap.get(groupToBeModified).getAtomList();
		for (Atom a : atomList) {
			if (a.getElement().equals("O")){
				oxygenAtoms.add(a);
			}
		}
		return oxygenAtoms;
	}
	

	private static void populateTerminalSingleAndDoubleBondedOxygen(List<Atom> atomList, LinkedList<Atom> singleBondedOxygen,LinkedList<Atom> doubleBondedOxygen) throws StructureBuildingException {
		for (Atom a : atomList) {
			if (a.getElement().equals("O")){//find terminal oxygens
				if (a.getBonds().size()==1){
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
