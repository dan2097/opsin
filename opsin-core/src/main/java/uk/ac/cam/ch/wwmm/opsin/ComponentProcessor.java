package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;
import static uk.ac.cam.ch.wwmm.opsin.OpsinTools.*;

/**Performs structure-aware destructive procedural parsing on parser results.
*
* @author dl387
*
*/

class ComponentProcessor {
	private final static Pattern matchAddedHydrogenBracket =Pattern.compile("[\\[\\(\\{]([^\\[\\(\\{]*)H[\\]\\)\\}]");
	private final static Pattern matchElementSymbolOrAminoAcidLocant = Pattern.compile("[A-Z][a-z]?'*(\\d+[a-z]?'*)?");
	private final static Pattern matchChalcogenReplacement= Pattern.compile("thio|seleno|telluro");
	private final static Pattern matchGroupsThatAreAlsoInlineSuffixes = Pattern.compile("carbon|oxy|sulfen|sulfin|sulfon|selenen|selenin|selenon|telluren|tellurin|telluron");
	private final static String[] traditionalAlkanePositionNames =new String[]{"alpha", "beta", "gamma", "delta", "epsilon", "zeta"};
	
	private final FunctionalReplacement functionalReplacement;
	private final SuffixRulesLookup suffixRulesLookup;
	private final BuildState state;
	
	//rings that look like HW rings but have other meanings. For the HW like inorganics the true meaning is given
	private static final Map<String, String[]> specialHWRings = new HashMap<String, String[]>();
	static{
		//The first entry of the array is a special instruction e.g. blocked or saturated. The correct order of the heteroatoms follows
		//terminal e is ignored from all of the keys as it is optional in the input name
		specialHWRings.put("oxin", new String[]{"blocked"});
		specialHWRings.put("azin", new String[]{"blocked"});
		
		specialHWRings.put("selenin", new String[]{"not_icacid", "Se","C","C","C","C","C"});
		specialHWRings.put("tellurin", new String[]{"not_icacid", "Te","C","C","C","C","C"});

		specialHWRings.put("thiol", new String[]{"not_nothingOrOlate", "S","C","C","C","C"});
		specialHWRings.put("selenol", new String[]{"not_nothingOrOlate", "Se","C","C","C","C"});
		specialHWRings.put("tellurol", new String[]{"not_nothingOrOlate", "Te","C","C","C","C"});
		
		specialHWRings.put("oxazol", new String[]{"","O","C","N","C","C"});
		specialHWRings.put("thiazol", new String[]{"","S","C","N","C","C"});
		specialHWRings.put("selenazol", new String[]{"","Se","C","N","C","C"});
		specialHWRings.put("tellurazol", new String[]{"","Te","C","N","C","C"});
		specialHWRings.put("oxazolidin", new String[]{"","O","C","N","C","C"});
		specialHWRings.put("thiazolidin", new String[]{"","S","C","N","C","C"});
		specialHWRings.put("selenazolidin", new String[]{"","Se","C","N","C","C"});
		specialHWRings.put("tellurazolidin", new String[]{"","Te","C","N","C","C"});
		specialHWRings.put("oxazolid", new String[]{"","O","C","N","C","C"});
		specialHWRings.put("thiazolid", new String[]{"","S","C","N","C","C"});
		specialHWRings.put("selenazolid", new String[]{"","Se","C","N","C","C"});
		specialHWRings.put("tellurazolid", new String[]{"","Te","C","N","C","C"});
		specialHWRings.put("oxazolin", new String[]{"","O","C","N","C","C"});
		specialHWRings.put("thiazolin", new String[]{"","S","C","N","C","C"});
		specialHWRings.put("selenazolin", new String[]{"","Se","C","N","C","C"});
		specialHWRings.put("tellurazolin", new String[]{"","Te","C","N","C","C"});

		specialHWRings.put("oxoxolan", new String[]{"","O","C","O","C","C"});
		specialHWRings.put("oxoxan", new String[]{"","O","C","C","O","C","C"});
		specialHWRings.put("oxoxin", new String[]{"","O","C","C","O","C","C"});

		specialHWRings.put("boroxin", new String[]{"saturated","O","B","O","B","O","B"});
		specialHWRings.put("borazin", new String[]{"saturated","N","B","N","B","N","B"});
		specialHWRings.put("borthiin", new String[]{"saturated","S","B","S","B","S","B"});
	}

	ComponentProcessor(SuffixRulesLookup suffixRules, BuildState state) {
		this.suffixRulesLookup = suffixRules;
		this.state = state;
		this.functionalReplacement = new FunctionalReplacement(state);
	}

	/**
	* Processes a parse result that has already gone through the ComponentGenerator.
	 * At this stage one can expect all substituents/roots to have at least 1 group.
	 * Multiple groups are present in, for example, fusion nomenclature. By the end of this function there will be exactly 1 group
	 * associated with each substituent/root. Multiplicative nomenclature can result in there being multiple roots
	 * @param parse 
	 * @throws ComponentGenerationException
	 * @throws StructureBuildingException
	 */
	void processParse(Element parse) throws ComponentGenerationException, StructureBuildingException {
		List<Element> words =OpsinTools.getDescendantElementsWithTagName(parse, WORD_EL);
		int wordCount =words.size();
		for (int i = wordCount -1; i>=0; i--) {
			Element word = words.get(i);
			String wordRule = OpsinTools.getParentWordRule(word).getAttributeValue(WORDRULE_EL);
			state.currentWordRule = WordRule.valueOf(wordRule);
			if (word.getAttributeValue(TYPE_ATR).equals(WordType.functionalTerm.toString())){
				continue;//functionalTerms are handled on a case by case basis by wordRules
			}

			List<Element> roots = OpsinTools.getDescendantElementsWithTagName(word, ROOT_EL);
			if (roots.size() > 1){
				throw new ComponentGenerationException("Multiple roots, but only 0 or 1 were expected. Found: " + roots.size());
			}
			List<Element> substituents = OpsinTools.getDescendantElementsWithTagName(word, SUBSTITUENT_EL);
			List<Element> substituentsAndRoot = OpsinTools.combineElementLists(substituents, roots);
			List<Element> brackets = OpsinTools.getDescendantElementsWithTagName(word, BRACKET_EL);
			List<Element> substituentsAndRootAndBrackets = OpsinTools.combineElementLists(substituentsAndRoot, brackets);
			List<Element> groups = OpsinTools.getDescendantElementsWithTagName(word, GROUP_EL);

			for (Element group : groups) {
				Fragment thisFrag = resolveGroup(state, group);
				processChargeAndOxidationNumberSpecification(group, thisFrag);//e.g. mercury(2+) or mercury(II)
			}
			
			for (Element subOrRoot : substituentsAndRoot) {
				applyDLPrefixes(subOrRoot);
				processCarbohydrates(subOrRoot);//e.g. glucopyranose (needs to be done before determineLocantMeaning to cope with alpha,beta for undefined anomer stereochemistry)
			}
			
			for (int j = substituents.size() -1; j >=0; j--) {
				Element substituent = substituents.get(j);
				boolean removed = removeAndMoveToAppropriateGroupIfHydroSubstituent(substituent);//this REMOVES a substituent just containing hydro/perhydro elements and moves these elements in front of an appropriate ring
				if (!removed){
					removed = removeAndMoveToAppropriateGroupIfSubtractivePrefix(substituent);
				}
				if (removed){
					substituents.remove(j);
					substituentsAndRoot.remove(substituent);
					substituentsAndRootAndBrackets.remove(substituent);
				}
			}
			
			Element finalSubOrRootInWord = word.getChild(word.getChildCount() - 1);
			while (!finalSubOrRootInWord.getName().equals(ROOT_EL) && !finalSubOrRootInWord.getName().equals(SUBSTITUENT_EL)){
				List<Element> children = OpsinTools.getChildElementsWithTagNames(finalSubOrRootInWord, new String[]{ROOT_EL, SUBSTITUENT_EL, BRACKET_EL});
				if (children.size() == 0){
					throw new ComponentGenerationException("Unable to find finalSubOrRootInWord");
				}
				finalSubOrRootInWord = children.get(children.size() - 1);
			}

			for (Element subOrRootOrBracket : substituentsAndRootAndBrackets) {
				determineLocantMeaning(subOrRootOrBracket, finalSubOrRootInWord);
			}

			for (Element subOrRoot : substituentsAndRoot) {
				processMultipliers(subOrRoot);
				detectConjunctiveSuffixGroups(subOrRoot, groups);
				matchLocantsToDirectFeatures(subOrRoot);

				List<Element> groupsOfSubOrRoot = subOrRoot.getChildElements(GROUP_EL);
				Element lastGroupInSubOrRoot =groupsOfSubOrRoot.get(groupsOfSubOrRoot.size() - 1);
				preliminaryProcessSuffixes(lastGroupInSubOrRoot, subOrRoot.getChildElements(SUFFIX_EL));
			}
			functionalReplacement.processAcidReplacingFunctionalClassNomenclature(finalSubOrRootInWord, word);

			if (functionalReplacement.processPrefixFunctionalReplacementNomenclature(groups, substituents)){//true if functional replacement performed, 1 or more substituents will have been removed
				substituentsAndRoot = OpsinTools.combineElementLists(substituents, roots);
				substituentsAndRootAndBrackets =OpsinTools.combineElementLists(substituentsAndRoot, brackets);
			}
			
			handleGroupIrregularities(groups);

			for (Element subOrRoot : substituentsAndRoot) {
				processHW(subOrRoot);//hantzch-widman rings
				FusedRingBuilder.processFusedRings(state, subOrRoot);
				processFusedRingBridges(subOrRoot);
				assignElementSymbolLocants(subOrRoot);
				processRingAssemblies(subOrRoot);
				processPolyCyclicSpiroNomenclature(subOrRoot);
			}

			for (Element subOrRoot : substituentsAndRoot) {
				applyLambdaConvention(subOrRoot);
				handleMultiRadicals(subOrRoot);
			}

			addImplicitBracketsToAminoAcids(groups, brackets);
			for (Element substituent : substituents) {
				matchLocantsToIndirectFeatures(substituent);
				addImplicitBracketsInCaseWhereSubstituentHasTwoLocants(substituent, brackets);
				implicitlyBracketToPreviousSubstituentIfAppropriate(substituent, brackets);
			}
			for (Element root : roots) {
				matchLocantsToIndirectFeatures(root);
			}

			for (Element subOrRoot : substituentsAndRoot) {
				assignImplicitLocantsToDiTerminalSuffixes(subOrRoot);
				processConjunctiveNomenclature(subOrRoot);
				resolveSuffixes(subOrRoot.getFirstChildElement(GROUP_EL), subOrRoot.getChildElements(SUFFIX_EL));
			}

			moveErroneouslyPositionedLocantsAndMultipliers(brackets);//e.g. (tetramethyl)azanium == tetra(methyl)azanium
			List<Element> children = OpsinTools.getChildElementsWithTagNames(word, new String[]{ROOT_EL, SUBSTITUENT_EL, BRACKET_EL});
			while (children.size() == 1){
				children = OpsinTools.getChildElementsWithTagNames(children.get(0), new String[]{ROOT_EL, SUBSTITUENT_EL, BRACKET_EL});
			}
			if (children.size() > 0){
				assignLocantsToMultipliedRootIfPresent(children.get(children.size() - 1));//multiplicative nomenclature e.g. methylenedibenzene or 3,4'-oxydipyridine
			}
			substituentsAndRootAndBrackets =OpsinTools.combineElementLists(substituentsAndRoot, brackets);//implicit brackets may have been created
			for (Element subBracketOrRoot : substituentsAndRootAndBrackets) {
				assignLocantsAndMultipliers(subBracketOrRoot);
			}
			processBiochemicalLinkageDescriptors(substituents, brackets);
			processWordLevelMultiplierIfApplicable(word, wordCount);
		}
		new WordRulesOmittedSpaceCorrector(state, parse).correctOmittedSpaces();//TODO where should this go?
	}

	/**Resolves the contents of a group element
	 *
	 * @param group The group element
	 * @return The fragment specified by the group element.
	 * @throws StructureBuildingException If the group can't be built.
	 * @throws ComponentGenerationException
	 */
	static Fragment resolveGroup(BuildState state, Element group) throws StructureBuildingException, ComponentGenerationException {
		String groupValue = group.getAttributeValue(VALUE_ATR);

		Fragment thisFrag;
		if (group.getAttribute(LABELS_ATR) != null){
			thisFrag = state.fragManager.buildSMILES(groupValue, group, group.getAttributeValue(LABELS_ATR));
		}
		else{
			thisFrag = state.fragManager.buildSMILES(groupValue, group, "");
		}
		group.setFrag(thisFrag);

		//processes groups like cymene and xylene whose structure is determined by the presence of a locant in front e.g. p-xylene
		processXyleneLikeNomenclature(state, group, thisFrag);

		setFragmentDefaultInAtomIfSpecified(thisFrag, group);
		setFragmentFunctionalAtomsIfSpecified(group, thisFrag);
		applyTraditionalAlkaneNumberingIfAppropriate(group, thisFrag); 
		return thisFrag;
	}

	/**
	 * Checks for groups with the addGroup/addBond/addHeteroAtom attributes. For the addGroup attribute adds the group defined by the SMILES described within
	 * e.g. for xylene  this function would add two methyls. Xylene is initially generated using the structure of benzene!
	 * See tokenList dtd for more information on the syntax of these attributes if it is not clear from the code
	 * @param state
	 * @param group: The group element
	 * @param parentFrag: The fragment that has been generated from the group element
	 * @throws StructureBuildingException
	 * @throws ComponentGenerationException
	 */
	private static void processXyleneLikeNomenclature(BuildState state, Element group, Fragment parentFrag) throws StructureBuildingException, ComponentGenerationException {
		if(group.getAttribute(ADDGROUP_ATR) != null) {
			String addGroupInformation = group.getAttributeValue(ADDGROUP_ATR);
			String[] groupsToBeAdded = MATCH_SEMICOLON.split(addGroupInformation);//typically only one, but 2 in the case of xylene and quinones
			List<Map<String, String>> allGroupInformation = new ArrayList<Map<String, String>>();
			for (String groupToBeAdded : groupsToBeAdded) {//populate allGroupInformation list
				String[] tempArray = MATCH_SPACE.split(groupToBeAdded);
				Map<String, String> groupInformation = new HashMap<String, String>();
				if (tempArray.length != 2 && tempArray.length != 3) {
					throw new ComponentGenerationException("malformed addGroup tag");
				}
				groupInformation.put("SMILES", tempArray[0]);
				if (tempArray[1].startsWith("id")) {
					groupInformation.put("atomReferenceType", "id");
					groupInformation.put("atomReference", tempArray[1].substring(2));
				} else if (tempArray[1].startsWith("locant")) {
					groupInformation.put("atomReferenceType", "locant");
					groupInformation.put("atomReference", tempArray[1].substring(6));
				} else {
					throw new ComponentGenerationException("malformed addGroup tag");
				}
				if (tempArray.length == 3) {//labels may optionally be specified for the group to be added
					groupInformation.put("labels", tempArray[2]);
				}
				allGroupInformation.add(groupInformation);
			}
			Element previousEl = OpsinTools.getPreviousSibling(group);
			if (previousEl !=null && previousEl.getName().equals(LOCANT_EL)){//has the name got specified locants to override the default ones
				List<String> locantValues = StringTools.arrayToList(MATCH_COMMA.split(previousEl.getValue()));
				if ((locantValues.size() == groupsToBeAdded.length || locantValues.size() +1 == groupsToBeAdded.length) && locantAreAcceptableForXyleneLikeNomenclatures(locantValues, group)){//one locant can be implicit in some cases
					boolean assignlocants =true;
					if (locantValues.size()!=groupsToBeAdded.length){
						//check that the firstGroup by default will be added to the atom with locant 1. If this is not the case then as many locants as there were groups should of been specified
						//or no locants should have been specified, which is what will be assumed (i.e. the locants will be left unassigned)
						Map<String, String> groupInformation = allGroupInformation.get(0);
						String locant;
						if (groupInformation.get("atomReferenceType").equals("locant")){
							locant =parentFrag.getAtomByLocantOrThrow(groupInformation.get("atomReference")).getFirstLocant();
						}
						else if (groupInformation.get("atomReferenceType").equals("id") ){
							locant =parentFrag.getAtomByIDOrThrow(parentFrag.getIdOfFirstAtom() + Integer.parseInt(groupInformation.get("atomReference")) -1 ).getFirstLocant();
						}
						else{
							throw new ComponentGenerationException("malformed addGroup tag");
						}
						if (locant ==null || !locant.equals("1")){
							assignlocants=false;
						}
					}
					if (assignlocants){
						for (int i = groupsToBeAdded.length -1; i >=0 ; i--) {
							//if less locants than expected are specified the locants of only the later groups will be changed
							//e.g. 4-xylene will transform 1,2-xylene to 1,4-xylene
							Map<String, String> groupInformation =allGroupInformation.get(i);
							if (locantValues.size() >0){
								groupInformation.put("atomReferenceType", "locant");
								groupInformation.put("atomReference", locantValues.get(locantValues.size()-1));
								locantValues.remove(locantValues.size()-1);
							}
							else{
								break;
							}
						}
						group.removeAttribute(group.getAttribute(FRONTLOCANTSEXPECTED_ATR));
						previousEl.detach();
					}
				}
			}

			for (int i = 0; i < groupsToBeAdded.length; i++) {
				Map<String, String> groupInformation = allGroupInformation.get(i);
				String smilesOfGroupToBeAdded = groupInformation.get("SMILES");
				Fragment newFrag;
				if (groupInformation.get("labels")!=null){
					newFrag = state.fragManager.buildSMILES(smilesOfGroupToBeAdded, group, groupInformation.get("labels"));
				}
				else{
					newFrag = state.fragManager.buildSMILES(smilesOfGroupToBeAdded, group, NONE_LABELS_VAL);
				}

				Atom atomOnParentFrag =null;
				if (groupInformation.get("atomReferenceType").equals("locant")){
					atomOnParentFrag=parentFrag.getAtomByLocantOrThrow(groupInformation.get("atomReference"));
				}
				else if (groupInformation.get("atomReferenceType").equals("id") ){
					atomOnParentFrag= parentFrag.getAtomByIDOrThrow(parentFrag.getIdOfFirstAtom() + Integer.parseInt(groupInformation.get("atomReference")) -1);
				}
				else{
					throw new ComponentGenerationException("malformed addGroup tag");
				}
				if (newFrag.getOutAtomCount() >1){
					throw new ComponentGenerationException("too many outAtoms on group to be added");
				}
				if (newFrag.getOutAtomCount() ==1) {
					OutAtom newFragOutAtom = newFrag.getOutAtom(0);
					newFrag.removeOutAtom(newFragOutAtom);
					state.fragManager.incorporateFragment(newFrag, newFragOutAtom.getAtom(), parentFrag, atomOnParentFrag, newFragOutAtom.getValency());
				}
				else{
					Atom atomOnNewFrag = newFrag.getDefaultInAtom();
					state.fragManager.incorporateFragment(newFrag, atomOnNewFrag, parentFrag, atomOnParentFrag, 1);
				}
			}
		}

		if(group.getAttributeValue(ADDHETEROATOM_ATR)!=null) {
			String addHeteroAtomInformation=group.getAttributeValue(ADDHETEROATOM_ATR);
			String[] heteroAtomsToBeAdded = MATCH_SEMICOLON.split(addHeteroAtomInformation);
			ArrayList<HashMap<String, String>> allHeteroAtomInformation = new ArrayList<HashMap<String, String>>();
			for (String heteroAtomToBeAdded : heteroAtomsToBeAdded) {//populate allHeteroAtomInformation list
				String[] tempArray = MATCH_SPACE.split(heteroAtomToBeAdded);
				HashMap<String, String> heteroAtomInformation = new HashMap<String, String>();
				if (tempArray.length != 2) {
					throw new ComponentGenerationException("malformed addHeteroAtom tag");
				}
				heteroAtomInformation.put("SMILES", tempArray[0]);
				if (tempArray[1].startsWith("id")) {
					heteroAtomInformation.put("atomReferenceType", "id");
					heteroAtomInformation.put("atomReference", tempArray[1].substring(2));
				} else if (tempArray[1].startsWith("locant")) {
					heteroAtomInformation.put("atomReferenceType", "locant");
					heteroAtomInformation.put("atomReference", tempArray[1].substring(6));
				} else {
					throw new ComponentGenerationException("malformed addHeteroAtom tag");
				}
				allHeteroAtomInformation.add(heteroAtomInformation);
			}
			Element previousEl = OpsinTools.getPreviousSibling(group);
			if (previousEl !=null && previousEl.getName().equals(LOCANT_EL)){//has the name got specified locants to override the default ones
				List<String> locantValues =StringTools.arrayToList(MATCH_COMMA.split(previousEl.getValue()));
				if (locantValues.size() ==heteroAtomsToBeAdded.length && locantAreAcceptableForXyleneLikeNomenclatures(locantValues, group)){
					for (int i = heteroAtomsToBeAdded.length -1; i >=0 ; i--) {//all heteroatoms must have a locant or default locants will be used
						HashMap<String, String> groupInformation =allHeteroAtomInformation.get(i);
						groupInformation.put("atomReferenceType", "locant");
						groupInformation.put("atomReference", locantValues.get(locantValues.size()-1));
						locantValues.remove(locantValues.size()-1);
					}
					group.removeAttribute(group.getAttribute(FRONTLOCANTSEXPECTED_ATR));
					previousEl.detach();
				}
			}

			for (int i = 0; i < heteroAtomsToBeAdded.length; i++) {
				HashMap<String, String> heteroAtomInformation =allHeteroAtomInformation.get(i);
				Atom atomOnParentFrag =null;
				if (heteroAtomInformation.get("atomReferenceType").equals("locant")){
					atomOnParentFrag=parentFrag.getAtomByLocantOrThrow(heteroAtomInformation.get("atomReference"));
				}
				else if (heteroAtomInformation.get("atomReferenceType").equals("id") ){
					atomOnParentFrag= parentFrag.getAtomByIDOrThrow(parentFrag.getIdOfFirstAtom() + Integer.parseInt(heteroAtomInformation.get("atomReference")) -1);
				}
				else{
					throw new ComponentGenerationException("malformed addHeteroAtom tag");
				}
				state.fragManager.replaceAtomWithSmiles(atomOnParentFrag, heteroAtomInformation.get("SMILES"));
			}
		}

		if(group.getAttributeValue(ADDBOND_ATR)!=null && !HANTZSCHWIDMAN_SUBTYPE_VAL.equals(group.getAttributeValue(SUBTYPE_ATR))) {//HW add bond is handled later
			String addBondInformation=group.getAttributeValue(ADDBOND_ATR);
			String[] bondsToBeAdded = MATCH_SEMICOLON.split(addBondInformation);
			ArrayList<HashMap<String, String>> allBondInformation = new ArrayList<HashMap<String, String>>();
			for (String bondToBeAdded : bondsToBeAdded) {//populate allBondInformation list
				String[] tempArray = MATCH_SPACE.split(bondToBeAdded);
				HashMap<String, String> bondInformation = new HashMap<String, String>();
				if (tempArray.length != 2) {
					throw new ComponentGenerationException("malformed addBond tag");
				}
				bondInformation.put("bondOrder", tempArray[0]);
				if (tempArray[1].startsWith("id")) {
					bondInformation.put("atomReferenceType", "id");
					bondInformation.put("atomReference", tempArray[1].substring(2));
				} else if (tempArray[1].startsWith("locant")) {
					bondInformation.put("atomReferenceType", "locant");
					bondInformation.put("atomReference", tempArray[1].substring(6));
				} else {
					throw new ComponentGenerationException("malformed addBond tag");
				}
				allBondInformation.add(bondInformation);
			}
			boolean locanted = false;
			Element previousEl = OpsinTools.getPreviousSibling(group);
			if (previousEl !=null && previousEl.getName().equals(LOCANT_EL)){//has the name got specified locants to override the default ones
				List<String> locantValues =StringTools.arrayToList(MATCH_COMMA.split(previousEl.getValue()));
				if (locantValues.size() ==bondsToBeAdded.length && locantAreAcceptableForXyleneLikeNomenclatures(locantValues, group)){
					for (int i = bondsToBeAdded.length -1; i >=0 ; i--) {//all bond order changes must have a locant or default locants will be used
						HashMap<String, String> bondInformation =allBondInformation.get(i);
						bondInformation.put("atomReferenceType", "locant");
						bondInformation.put("atomReference", locantValues.get(locantValues.size()-1));
						locantValues.remove(locantValues.size()-1);
					}
					group.removeAttribute(group.getAttribute(FRONTLOCANTSEXPECTED_ATR));
					previousEl.detach();
					locanted = true;
				}
			}

			for (int i = 0; i < bondsToBeAdded.length; i++) {
				HashMap<String, String> bondInformation =allBondInformation.get(i);
				Atom atomOnParentFrag =null;
				if (bondInformation.get("atomReferenceType").equals("locant")){
					atomOnParentFrag=parentFrag.getAtomByLocantOrThrow(bondInformation.get("atomReference"));
				}
				else if (bondInformation.get("atomReferenceType").equals("id") ){
					atomOnParentFrag= parentFrag.getAtomByIDOrThrow(parentFrag.getIdOfFirstAtom() + Integer.parseInt(bondInformation.get("atomReference")) -1);
				}
				else{
					throw new ComponentGenerationException("malformed addBond tag");
				}
				
				Bond b = FragmentTools.unsaturate(atomOnParentFrag, Integer.parseInt(bondInformation.get("bondOrder")) , parentFrag);
				if (!locanted && b.getOrder() ==2 && 
						parentFrag.getAtomCount()==5 &&
						b.getFromAtom().getAtomIsInACycle() &&
						b.getToAtom().getAtomIsInACycle()){
					//special case just that substitution of groups like imidazoline may actually remove the double bond...
					b.setOrder(1);
					b.getFromAtom().setSpareValency(true);
					b.getToAtom().setSpareValency(true);
				}
			}
		}
	}

	/**
	 * Checks that all locants are present within the front locants expected attribute of the group
	 * @param locantValues
	 * @param group
	 * @return
	 */
	private static boolean locantAreAcceptableForXyleneLikeNomenclatures(List<String> locantValues, Element group) {
		if (group.getAttribute(FRONTLOCANTSEXPECTED_ATR) == null){
			throw new IllegalArgumentException("Group must have frontLocantsExpected to implement xylene-like nomenclature");
		}
		List<String> allowedLocants = Arrays.asList(MATCH_COMMA.split(group.getAttributeValue(FRONTLOCANTSEXPECTED_ATR)));
		for (String locant : locantValues) {
			if (!allowedLocants.contains(locant)){
				return false;
			}
		}
		return true;
	}


	/**
	 * Looks for the presence of {@link XmlDeclarations#DEFAULTINLOCANT_ATR} and {@link XmlDeclarations#DEFAULTINID_ATR} on the group and applies them to the fragment
	 * @param thisFrag
	 * @param group
	 * @throws StructureBuildingException
	 */
	private static void setFragmentDefaultInAtomIfSpecified(Fragment thisFrag, Element group) throws StructureBuildingException {
		String defaultInLocant = group.getAttributeValue(DEFAULTINLOCANT_ATR);
		String defaultInId = group.getAttributeValue(DEFAULTINID_ATR);
		if (defaultInLocant != null){//sets the atom at which substitution will occur to by default
			thisFrag.setDefaultInAtom(thisFrag.getAtomByLocantOrThrow(defaultInLocant));
		}
		else if (defaultInId != null){
			thisFrag.setDefaultInAtom(thisFrag.getAtomByIDOrThrow(thisFrag.getIdOfFirstAtom() + Integer.parseInt(defaultInId) - 1));
		}
	}


	/**
	 * Looks for the presence of FUNCTIONALIDS_ATR on the group and applies them to the fragment
	 * @param group
	 * @param thisFrag
	 * @throws StructureBuildingException
	 */
	private static void setFragmentFunctionalAtomsIfSpecified(Element group, Fragment thisFrag) throws StructureBuildingException {
		if (group.getAttribute(FUNCTIONALIDS_ATR)!=null){
			String[] functionalIDs = MATCH_COMMA.split(group.getAttributeValue(FUNCTIONALIDS_ATR));
			for (String functionalID : functionalIDs) {
				thisFrag.addFunctionalAtom(thisFrag.getAtomByIDOrThrow(thisFrag.getIdOfFirstAtom() + Integer.parseInt(functionalID) - 1));
			}
		}
	}


	private static void applyTraditionalAlkaneNumberingIfAppropriate(Element group, Fragment thisFrag)  {
		String groupType  = group.getAttributeValue(TYPE_ATR);
		if (groupType.equals(ACIDSTEM_TYPE_VAL)){
			List<Atom> atomList = thisFrag.getAtomList();
			Atom startingAtom = thisFrag.getFirstAtom();
			if (group.getAttribute(SUFFIXAPPLIESTO_ATR) != null){
				String suffixAppliesTo = group.getAttributeValue(SUFFIXAPPLIESTO_ATR);
				String suffixAppliesToArr[] = MATCH_COMMA.split(suffixAppliesTo);
				if (suffixAppliesToArr.length != 1){
					return;
				}
				startingAtom = atomList.get(Integer.parseInt(suffixAppliesToArr[0]) - 1);
			}
			List<Atom> neighbours = startingAtom.getAtomNeighbours();
			int counter = -1;
			Atom previousAtom = startingAtom;
			for (int i = neighbours.size() - 1; i >=0; i--) {//only consider carbon atoms
				if (neighbours.get(i).getElement() != ChemEl.C){
					neighbours.remove(i);
				}
			}
			while (neighbours.size() == 1){
				counter++;
				if (counter > 5){
					break;
				}
				Atom nextAtom = neighbours.get(0);
				if (nextAtom.getAtomIsInACycle()){
					break;
				}
				String traditionalLocant = traditionalAlkanePositionNames[counter];
				if (!nextAtom.hasLocant(traditionalLocant)){
					nextAtom.addLocant(traditionalLocant);
				}
				neighbours = nextAtom.getAtomNeighbours();
				neighbours.remove(previousAtom);
				for (int i = neighbours.size()-1; i >=0; i--) {//only consider carbon atoms
					if (neighbours.get(i).getElement() != ChemEl.C){
						neighbours.remove(i);
					}
				}
				previousAtom = nextAtom;
			}
		}
		else if (groupType.equals(CHAIN_TYPE_VAL) && ALKANESTEM_SUBTYPE_VAL.equals(group.getAttributeValue(SUBTYPE_ATR))){
			List<Atom> atomList = thisFrag.getAtomList();
			if (atomList.size() == 1){
				return;
			}
			Element possibleSuffix = OpsinTools.getNextSibling(group, SUFFIX_EL);
			Boolean terminalSuffixWithNoSuffixPrefixPresent = false;
			if (possibleSuffix!=null && TERMINAL_SUBTYPE_VAL.equals(possibleSuffix.getAttributeValue(SUBTYPE_ATR)) && possibleSuffix.getAttribute(SUFFIXPREFIX_ATR) == null){
				terminalSuffixWithNoSuffixPrefixPresent = true;
			}
			for (Atom atom : atomList) {
				String firstLocant = atom.getFirstLocant();
				if (!atom.getAtomIsInACycle() && firstLocant != null && firstLocant.length() == 1 && Character.isDigit(firstLocant.charAt(0))){
					int locantNumber = Integer.parseInt(firstLocant);
					if (terminalSuffixWithNoSuffixPrefixPresent){
						if (locantNumber > 1 && locantNumber <= 7){
							atom.addLocant(traditionalAlkanePositionNames[locantNumber - 2]);
						}
					}
					else{
						if (locantNumber > 0 && locantNumber <= 6){
							atom.addLocant(traditionalAlkanePositionNames[locantNumber - 1]);
						}
					}
				}
			}
		}
	}
	
	private void processChargeAndOxidationNumberSpecification(Element group, Fragment frag)  {
		Element nextEl = OpsinTools.getNextSibling(group);
		if (nextEl != null){
			if(nextEl.getName().equals(CHARGESPECIFIER_EL)) {
				frag.getFirstAtom().setCharge(Integer.parseInt(nextEl.getAttributeValue(VALUE_ATR)));
				nextEl.detach();
			}
			if(nextEl.getName().equals(OXIDATIONNUMBERSPECIFIER_EL)) {
				frag.getFirstAtom().setProperty(Atom.OXIDATION_NUMBER, Integer.parseInt(nextEl.getAttributeValue(VALUE_ATR)));
				nextEl.detach();
			}
		}
	}

	/**
	 * Removes substituents which are just a hydro/perhydro element and moves their contents to be in front of the next in scope ring
	 * @param substituent
	 * @return true is the substituent was a hydro substituent and hence was removed
	 * @throws ComponentGenerationException
	 */
	private boolean removeAndMoveToAppropriateGroupIfHydroSubstituent(Element substituent) throws ComponentGenerationException {
		List<Element> hydroElements = substituent.getChildElements(HYDRO_EL);
		if (hydroElements.size() > 0 && substituent.getChildElements(GROUP_EL).size() == 0){
			Element hydroSubstituent = substituent;
			if (hydroElements.size() != 1){
				throw new ComponentGenerationException("Unexpected number of hydro elements found in substituent");
			}
			Element hydroElement = hydroElements.get(0);
			String hydroValue = hydroElement.getValue();
			if (hydroValue.equals("hydro")){
				Element multiplier = OpsinTools.getPreviousSibling(hydroElement);
				if (multiplier == null || !multiplier.getName().equals(MULTIPLIER_EL) ){
					throw new ComponentGenerationException("Multiplier expected but not found before hydro subsituent");
				}
				if (Integer.parseInt(multiplier.getAttributeValue(VALUE_ATR)) %2 !=0){
					throw new ComponentGenerationException("Hydro can only be added in pairs but multiplier was odd: " + multiplier.getAttributeValue(VALUE_ATR));
				}
			}
			Element targetRing = null;
			Element nextSubOrRootOrBracket = OpsinTools.getNextSibling(hydroSubstituent);
			if (nextSubOrRootOrBracket == null){
				throw new ComponentGenerationException("Cannot find ring for hydro substituent to apply to");
			}
			//first check adjacent substituent/root. If the hydroelement has one locant or the ring is locantless then we can assume the hydro is acting as a nondetachable prefix
			Element potentialRing = nextSubOrRootOrBracket.getFirstChildElement(GROUP_EL);
			if (potentialRing != null && containsCyclicAtoms(potentialRing)){
				Element possibleLocantInFrontOfHydro = OpsinTools.getPreviousSiblingIgnoringCertainElements(hydroElement, new String[]{MULTIPLIER_EL});
				if (possibleLocantInFrontOfHydro != null && possibleLocantInFrontOfHydro.getName().equals(LOCANT_EL) && MATCH_COMMA.split(possibleLocantInFrontOfHydro.getValue()).length == 1){
					//e.g.4-decahydro-1-naphthalenyl
					targetRing =potentialRing;
				}
				else{
					Element possibleLocantInFrontOfRing = OpsinTools.getPreviousSibling(potentialRing, LOCANT_EL);
					if (possibleLocantInFrontOfRing != null){
						if (potentialRing.getAttribute(FRONTLOCANTSEXPECTED_ATR) != null){//check whether the group was expecting a locant e.g. 2-furyl
							String locantValue = possibleLocantInFrontOfRing.getValue();
							String[] expectedLocants = MATCH_COMMA.split(potentialRing.getAttributeValue(FRONTLOCANTSEXPECTED_ATR));
							for (String expectedLocant : expectedLocants) {
								if (locantValue.equals(expectedLocant)){
									targetRing =potentialRing;
									break;
								}
							}
						}
						//check whether the group is a HW system e.g. 1,3-thiazole
						if (potentialRing.getAttributeValue(SUBTYPE_ATR).equals(HANTZSCHWIDMAN_SUBTYPE_VAL)){
							String locantValue = possibleLocantInFrontOfRing.getValue();
							int locants = MATCH_COMMA.split(locantValue).length;
							int heteroCount = 0;
							Element currentElem = OpsinTools.getNextSibling(possibleLocantInFrontOfRing);
							while(!currentElem.equals(potentialRing)){
								if(currentElem.getName().equals(HETEROATOM_EL)) {
									heteroCount++;
								} else if (currentElem.getName().equals(MULTIPLIER_EL)){
									heteroCount += Integer.parseInt(currentElem.getAttributeValue(VALUE_ATR)) -1;
								}
								currentElem = OpsinTools.getNextSibling(currentElem);
							}
							if (heteroCount==locants){//number of locants must match number
								targetRing =potentialRing;
							}
						}
						//check whether the group is a benzofused ring e.g. 1,4-benzodioxin
						if (FUSIONRING_SUBTYPE_VAL.equals(potentialRing.getAttributeValue(SUBTYPE_ATR)) && 
								(potentialRing.getValue().equals("benzo")|| potentialRing.getValue().equals("benz")) &&
								!OpsinTools.getNextSibling(potentialRing).getName().equals(FUSION_EL)){
							targetRing =potentialRing;
						}
					}
					else{
						targetRing =potentialRing;
					}
				}
			}
	
			//that didn't match so the hydro appears to be a detachable prefix. detachable prefixes attach in preference to the rightmost applicable group so search any remaining substituents/roots from right to left
			if (targetRing == null){
				Element nextSubOrRootOrBracketFromLast = hydroSubstituent.getParent().getChild(hydroSubstituent.getParent().getChildCount()-1);//the last sibling
				while (!nextSubOrRootOrBracketFromLast.equals(hydroSubstituent)){
					potentialRing = nextSubOrRootOrBracketFromLast.getFirstChildElement(GROUP_EL);
					if (potentialRing!=null && containsCyclicAtoms(potentialRing)){
						targetRing =potentialRing;
						break;
					}
					else{
						nextSubOrRootOrBracketFromLast = OpsinTools.getPreviousSibling(nextSubOrRootOrBracketFromLast);
					}
				}
			}
			if (targetRing == null){
				throw new ComponentGenerationException("Cannot find ring for hydro substituent to apply to");
			}
			//move the children of the hydro substituent
			List<Element> children = hydroSubstituent.getChildElements();
			for (int i = children.size()-1; i >=0 ; i--) {
				Element child = children.get(i);
				if (!child.getName().equals(HYPHEN_EL)){
					child.detach();
					targetRing.getParent().insertChild(child, 0);
				}
			}
			hydroSubstituent.detach();
			return true;
		}
		return false;
	}
	

	/**
	 * Removes substituents which are just a subtractivePrefix element e.g. deoxy and moves their contents to be in front of the next in scope biochemical fragment (or failing that group)
	 * @param substituent
	 * @return true is the substituent was a subtractivePrefix substituent and hence was removed
	 * @throws ComponentGenerationException
	 */
	static boolean removeAndMoveToAppropriateGroupIfSubtractivePrefix(Element substituent) throws ComponentGenerationException {
		List<Element> subtractivePrefixes = substituent.getChildElements(SUBTRACTIVEPREFIX_EL);
		if (subtractivePrefixes.size() > 0){
			if (subtractivePrefixes.size() != 1){
				throw new RuntimeException("Unexpected number of subtractive prefixes found in substituent");
			}
			Element subtractivePrefix = subtractivePrefixes.get(0);
			Element biochemicalGroup = null;//preferred
			Element standardGroup = null;
			Element nextSubOrRootOrBracket = OpsinTools.getNextSibling(substituent);
			if (nextSubOrRootOrBracket == null){
				throw new ComponentGenerationException("Unable to find group for: " + subtractivePrefix.getValue() +" to apply to!");
			}
			//prefer the nearest (unlocanted) biochemical group or the rightmost standard group
			while (nextSubOrRootOrBracket != null){
				Element groupToConsider = nextSubOrRootOrBracket.getFirstChildElement(GROUP_EL);
				if (groupToConsider!=null){
					if (BIOCHEMICAL_SUBTYPE_VAL.equals(groupToConsider.getAttributeValue(SUBTYPE_ATR)) || groupToConsider.getAttributeValue(TYPE_ATR).equals(CARBOHYDRATE_TYPE_VAL)){
						biochemicalGroup = groupToConsider;
						if (OpsinTools.getPreviousSiblingsOfType(biochemicalGroup, LOCANT_EL).size() == 0){
							break;
						}
					}
					else {
						standardGroup = groupToConsider;
					}
				}
				nextSubOrRootOrBracket = OpsinTools.getNextSibling(nextSubOrRootOrBracket);
			}
			
			Element targetGroup = biochemicalGroup!=null ? biochemicalGroup : standardGroup;
			if (targetGroup == null){
				throw new ComponentGenerationException("Unable to find group for: " + subtractivePrefix.getValue() +" to apply to!");
			}
			if (subtractivePrefix.getAttributeValue(TYPE_ATR).equals(ANHYDRO_TYPE_VAL)){
				Element locant = OpsinTools.getPreviousSibling(subtractivePrefix);
				if (locant == null || !locant.getName().equals(LOCANT_EL)){
					throw new ComponentGenerationException("Two locants are required before an anhydro prefix");
				}
				String locantStr = locant.getValue();
				if (MATCH_COMMA.split(locantStr).length != 2){
					throw new ComponentGenerationException("Two locants are required before an anhydro prefix, but found: "+ locantStr);
				}
				subtractivePrefix.addAttribute(new Attribute(LOCANT_ATR, locantStr));
				locant.detach();
			}
			
			//move the children of the subtractivePrefix substituent
			List<Element> children =substituent.getChildElements();
			for (int i = children.size()-1; i >=0 ; i--) {
				Element child =children.get(i);
				if (!child.getName().equals(HYPHEN_EL)){
					child.detach();
					targetGroup.getParent().insertChild(child, 0);
				}
			}
			substituent.detach();
			return true;
		}
		return false;
	}


	private boolean containsCyclicAtoms(Element potentialRing) {
		Fragment potentialRingFrag = potentialRing.getFrag();
		List<Atom> atomList = potentialRingFrag.getAtomList();
		for (Atom atom : atomList) {
			if (atom.getAtomIsInACycle()){
				return true;
			}
		}
		return false;
	}

	/**
	 * Checks for agreement between the number of locants and multipliers.
	 * If a locant element contains multiple elements and is not next to a multiplier the various cases where this is the case will be checked for
	 * This may result in a locant being moved if it is more convenient for subsequent processing
	 * @param subOrBracketOrRoot The substituent/root/bracket to looks for locants in.
	 * @param finalSubOrRootInWord : used to check if a locant is referring to the root as in multiplicative nomenclature
	 * @throws ComponentGenerationException
	 * @throws StructureBuildingException
	 */
	private void determineLocantMeaning(Element subOrBracketOrRoot, Element finalSubOrRootInWord) throws StructureBuildingException, ComponentGenerationException {
		List<Element> locants = subOrBracketOrRoot.getChildElements(LOCANT_EL);
		Element group = subOrBracketOrRoot.getFirstChildElement(GROUP_EL);//will be null if element is a bracket
		for (Element locant : locants) {
			String[] locantValues = MATCH_COMMA.split(locant.getValue());
			if(locantValues.length > 1) {
				Element afterLocant = OpsinTools.getNextSibling(locant);
				int structuralBracketDepth = 0;
				Element multiplierEl = null;
				while (afterLocant != null){
					String elName = afterLocant.getName();
					if (elName.equals(STRUCTURALOPENBRACKET_EL)){
						structuralBracketDepth++;
					}
					else if (elName.equals(STRUCTURALCLOSEBRACKET_EL)){
						structuralBracketDepth--;
					}
					if (structuralBracketDepth != 0){
						afterLocant = OpsinTools.getNextSibling(afterLocant);
						continue;
					}
					if(elName.equals(LOCANT_EL)) {
						break;
					}
					else if (elName.equals(MULTIPLIER_EL)){
						if (locantValues.length == Integer.parseInt(afterLocant.getAttributeValue(VALUE_ATR))){
							if (afterLocant.equals(OpsinTools.getNextSiblingIgnoringCertainElements(locant, new String[]{INDICATEDHYDROGEN_EL}))){
								//direct locant, typical case. An exception is made for indicated hydrogen e.g. 1,2,4-1H-triazole
								multiplierEl = afterLocant;
								break;
							}
							else{
								Element afterMultiplier = OpsinTools.getNextSibling(afterLocant);
								if (afterMultiplier!=null && (afterMultiplier.getName().equals(SUFFIX_EL) || afterMultiplier.getName().equals(INFIX_EL)
										|| afterMultiplier.getName().equals(UNSATURATOR_EL) || afterMultiplier.getName().equals(GROUP_EL))){
									multiplierEl = afterLocant; //indirect locant
									break;
								}
							}
						}
						if (afterLocant.equals(OpsinTools.getNextSibling(locant))){//if nothing better can be found report this as a locant/multiplier mismatch
							multiplierEl = afterLocant;
						}
					}
					else if (elName.equals(RINGASSEMBLYMULTIPLIER_EL) && afterLocant.equals(OpsinTools.getNextSibling(locant))){//e.g. 1,1'-biphenyl
						multiplierEl = afterLocant;
						if (!FragmentTools.allAtomsInRingAreIdentical(group.getFrag())){//if all atoms are identical then the locant may refer to suffixes
							break;
						}
					}
					else if (elName.equals(FUSEDRINGBRIDGE_EL)&& locantValues.length ==2 && afterLocant.equals(OpsinTools.getNextSibling(locant))){//e.g. 1,8-methano
						break;
					}
					afterLocant = OpsinTools.getNextSibling(afterLocant);
				}
				if(multiplierEl != null) {
					if(Integer.parseInt(multiplierEl.getAttributeValue(VALUE_ATR)) == locantValues.length ) {
						// number of locants and multiplier agree
						boolean locantModified = false;//did determineLocantMeaning do something?
						if (locantValues[locantValues.length-1].endsWith("'") && group!=null && subOrBracketOrRoot.indexOf(group) > subOrBracketOrRoot.indexOf(locant)){//quite possible that this is referring to a multiplied root
							if (group.getAttribute(OUTIDS_ATR)!=null && MATCH_COMMA.split(group.getAttributeValue(OUTIDS_ATR)).length>1){
								locantModified = checkSpecialLocantUses(locant, locantValues, finalSubOrRootInWord);
							}
							else{
								Element afterGroup = OpsinTools.getNextSibling(group);
								int inlineSuffixCount =0;
								int multiplier = 1;
								while (afterGroup != null){
									if(afterGroup.getName().equals(MULTIPLIER_EL)){
										multiplier =Integer.parseInt(afterGroup.getAttributeValue(VALUE_ATR));
									}
									else if(afterGroup.getName().equals(SUFFIX_EL) && afterGroup.getAttributeValue(TYPE_ATR).equals(INLINE_TYPE_VAL)){
										inlineSuffixCount +=(multiplier);
										multiplier=1;
									}
									afterGroup = OpsinTools.getNextSibling(afterGroup);
								}
								if (inlineSuffixCount >=2){
									locantModified = checkSpecialLocantUses(locant, locantValues, finalSubOrRootInWord);
								}
							}
						}
						if (!locantModified && !OpsinTools.getNextSibling(locant).equals(multiplierEl)){//the locants apply indirectly the multiplier e.g. 2,3-butandiol
							//move the locant to be next to the multiplier.
							locant.detach();
							OpsinTools.insertBefore(multiplierEl, locant);
						}
					} else {
						if(!checkSpecialLocantUses(locant, locantValues, finalSubOrRootInWord)) {
							throw new ComponentGenerationException("Mismatch between locant and multiplier counts (" + Integer.toString(locantValues.length) + " and " + multiplierEl.getAttributeValue(VALUE_ATR) + "):" + locant.getValue());
						}
					}
				} else {
					/* Multiple locants without a multiplier */
					if(!checkSpecialLocantUses(locant, locantValues, finalSubOrRootInWord)) {
						throw new ComponentGenerationException("Multiple locants without a multiplier: " + locant.toXML());
					}
				}
			}
		}
	}


	/**Looks for Hantzch-Widman systems, and sees if the number of locants
	 * agrees with the number of heteroatoms.
	 * If this is not the case alternative possibilities are tested:
	 * 	The locants could be intended to indicate the position of outAtoms e.g. 1,4-phenylene
	 * 	The locants could be intended to indicate the attachement points of the root groups in multiplicative nomenclature e.g. 4,4'-methylenedibenzoic acid
	 * @param locant The element corresponding to the locant group to be tested
	 * @param locantValues The locant values;
	 * @param finalSubOrRootInWord : used to check if a locant is referring to the root as in multiplicative nomenclatures)
	 * @return true if there's a HW system, and agreement; or if the locants conform to one of the alternative possibilities, otherwise false.
	 */
	private boolean checkSpecialLocantUses(Element locant, String[] locantValues, Element finalSubOrRootInWord) {
		int count = locantValues.length;
		Element currentElem = OpsinTools.getNextSibling(locant);
		int heteroCount = 0;
		int multiplierValue = 1;
		while(currentElem != null && !currentElem.getName().equals(GROUP_EL)){
			if(currentElem.getName().equals(HETEROATOM_EL)) {
				heteroCount += multiplierValue;
				multiplierValue = 1;
			} else if (currentElem.getName().equals(MULTIPLIER_EL)){
				multiplierValue = Integer.parseInt(currentElem.getAttributeValue(VALUE_ATR));
			}
			else{
				break;
			}
			currentElem = OpsinTools.getNextSibling(currentElem);
		}
		if(currentElem != null && currentElem.getName().equals(GROUP_EL)){
			if (currentElem.getAttributeValue(SUBTYPE_ATR).equals(HANTZSCHWIDMAN_SUBTYPE_VAL)) {
				if(heteroCount == count) {
					return true;
				} else if (heteroCount > 1){
					return false;//there is a case where locants don't apply to heteroatoms in a HW system, but in that case only one locant is expected so this function would not be called
				}
			}
			if (heteroCount == 0 && currentElem.getAttribute(OUTIDS_ATR) != null ) {//e.g. 1,4-phenylene
				String[] outIDs = MATCH_COMMA.split(currentElem.getAttributeValue(OUTIDS_ATR), -1);
				Fragment groupFragment = currentElem.getFrag();
				if (count ==outIDs.length && groupFragment.getAtomCount() > 1){//things like oxy do not need to have their outIDs specified
					int idOfFirstAtomInFrag =groupFragment.getIdOfFirstAtom();
					boolean foundLocantNotPresentOnFragment = false;
					for (int i = outIDs.length - 1; i >=0; i--) {
						Atom a =groupFragment.getAtomByLocant(locantValues[i]);
						if (a == null){
							foundLocantNotPresentOnFragment = true;
							break;
						}
						outIDs[i] = Integer.toString(a.getID() - idOfFirstAtomInFrag + 1);//convert to relative id
					}
					if (!foundLocantNotPresentOnFragment){
						currentElem.getAttribute(OUTIDS_ATR).setValue(StringTools.arrayToString(outIDs, ","));
						locant.detach();
						return true;
					}
				}
			}
			else if(currentElem.getValue().equals("benz") || currentElem.getValue().equals("benzo")){
				Element potentialGroupAfterBenzo = OpsinTools.getNextSibling(currentElem, GROUP_EL);//need to make sure this isn't benzyl
				if (potentialGroupAfterBenzo!=null){
					return true;//e.g. 1,2-benzothiazole
				}
			}
		}
		if(currentElem != null) {
			String name = currentElem.getName();
			if (name.equals(POLYCYCLICSPIRO_EL)){
				return true;
			}
			else if (name.equals(FUSEDRINGBRIDGE_EL) && count == 2){
				return true;
			}
			else if (name.equals(SUFFIX_EL) && CYCLEFORMER_SUBTYPE_VAL.equals(currentElem.getAttributeValue(SUBTYPE_ATR)) && count == 2){
				currentElem.addAttribute(new Attribute(LOCANT_ATR, locant.getValue()));
				locant.detach();
				return true;
			}
		}
		boolean detectedMultiplicativeNomenclature = detectMultiplicativeNomenclature(locant, locantValues, finalSubOrRootInWord);
		if (detectedMultiplicativeNomenclature){
			return true;
		}
		if (currentElem != null && count ==2 && currentElem.getName().equals(GROUP_EL) && EPOXYLIKE_SUBTYPE_VAL.equals(currentElem.getAttributeValue(SUBTYPE_ATR))){
			return true;
		}
		Element parentElem = locant.getParent();
		if (count == 2 && parentElem.getName().equals(BRACKET_EL)){//e.g. 3,4-(dichloromethylenedioxy) this is changed to (dichloro3,4-methylenedioxy)
			List<Element> substituents = parentElem.getChildElements(SUBSTITUENT_EL);
			if (substituents.size() > 0){
				Element finalSub  = substituents.get(substituents.size() - 1);
				Element group = finalSub.getFirstChildElement(GROUP_EL);
				if (EPOXYLIKE_SUBTYPE_VAL.equals(group.getAttributeValue(SUBTYPE_ATR))){
					locant.detach();
					OpsinTools.insertBefore(group, locant);
					return true;
				}
			}
		}
		
		return false;
	}


	/**
	 * Detects multiplicative nomenclature. If it does then the locant will be moved, changed to a multiplicative locant and true will be returned
	 * @param locant
	 * @param locantValues
	 * @param finalSubOrRootInWord
	 * @return
	 */
	private boolean detectMultiplicativeNomenclature(Element locant, String[] locantValues, Element finalSubOrRootInWord) {
		int count =locantValues.length;
		Element multiplier = finalSubOrRootInWord.getChild(0);
		if (finalSubOrRootInWord.getParent().getName().equals(BRACKET_EL)){//e.g. 1,1'-ethynediylbis(1-cyclopentanol)
			if (!multiplier.getName().equals(MULTIPLIER_EL)){
				multiplier = finalSubOrRootInWord.getParent().getChild(0);
			}
			else{
				Element elAfterMultiplier = OpsinTools.getNextSibling(multiplier);
				String elName = elAfterMultiplier.getName();
				if (elName.equals(HETEROATOM_EL) || elName.equals(SUBTRACTIVEPREFIX_EL)|| (elName.equals(HYDRO_EL) && !elAfterMultiplier.getValue().startsWith("per"))|| elName.equals(FUSEDRINGBRIDGE_EL)) {
					multiplier = finalSubOrRootInWord.getParent().getChild(0);
				}
			}
		}
		Element commonParent = locant.getParent().getParent();//this should be a common parent of the multiplier in front of the root. If it is not, then this locant is in a different scope
		Element parentOfMultiplier =multiplier.getParent();
		while (parentOfMultiplier!=null){
			if (commonParent.equals(parentOfMultiplier)){
				if (locantValues[count-1].endsWith("'")  &&
						multiplier.getName().equals(MULTIPLIER_EL) && !OpsinTools.getNextSibling(multiplier).getName().equals(MULTIPLICATIVELOCANT_EL) &&
						Integer.parseInt(multiplier.getAttributeValue(VALUE_ATR)) == count ){//multiplicative nomenclature
					locant.setName(MULTIPLICATIVELOCANT_EL);
					locant.detach();
					OpsinTools.insertAfter(multiplier, locant);
					return true;
				}
			}
			parentOfMultiplier = parentOfMultiplier.getParent();
		}
		return false;
	}

	private void applyDLPrefixes(Element subOrRoot) throws ComponentGenerationException {
		List<Element> dlStereochemistryEls = subOrRoot.getChildElements(DLSTEREOCHEMISTRY_EL);
		for (Element dlStereochemistry : dlStereochemistryEls) {
			String dlStereochemistryValue = dlStereochemistry.getAttributeValue(VALUE_ATR);
			Element elementToApplyTo = OpsinTools.getNextSibling(dlStereochemistry);
			if (elementToApplyTo == null){
				throw new RuntimeException("OPSIN bug: DL stereochemistry found in inappropriate position");
			}
			if (AMINOACID_TYPE_VAL.equals(elementToApplyTo.getAttributeValue(TYPE_ATR))){
				applyDlStereochemistryToAminoAcid(elementToApplyTo, dlStereochemistryValue);
			}
			else if (elementToApplyTo.getAttributeValue(TYPE_ATR).equals(CARBOHYDRATE_TYPE_VAL)){
				applyDlStereochemistryToCarbohydrate(elementToApplyTo, dlStereochemistryValue);
			}
			else if (CARBOHYDRATECONFIGURATIONPREFIX_TYPE_VAL.equals(elementToApplyTo.getAttributeValue(TYPE_ATR))){
				applyDlStereochemistryToCarbohydrateConfigurationalPrefix(elementToApplyTo, dlStereochemistryValue);
			}
			else{
				throw new RuntimeException("OPSIN bug: Unrecognised element after DL stereochemistry: " +elementToApplyTo.toXML());
			}
			dlStereochemistry.detach();
		}
	}

	void applyDlStereochemistryToAminoAcid(Element aminoAcidEl, String dlStereochemistryValue) throws ComponentGenerationException {
		Fragment aminoAcid = aminoAcidEl.getFrag();
		List<Atom> atomList = aminoAcid.getAtomList();
		List<Atom> atomsWithParities = new ArrayList<Atom>();
		for (Atom atom : atomList) {
			if (atom.getAtomParity()!=null){
				atomsWithParities.add(atom);
			}
		}
		if (atomsWithParities.isEmpty()){
			throw new ComponentGenerationException("D/L stereochemistry :" +dlStereochemistryValue + " found before achiral amino acid");
		}
		if (dlStereochemistryValue.equals("dl")){
			for (Atom atom : atomsWithParities) {
				atom.setAtomParity(null);
			}
		}
		else{
			boolean invert;
			if (dlStereochemistryValue.equals("l") || dlStereochemistryValue.equals("ls")){
				invert = false;
			} else if (dlStereochemistryValue.equals("d") || dlStereochemistryValue.equals("ds")){
				invert = true;
			} else{
				throw new ComponentGenerationException("Unexpected value for D/L stereochemistry found before amino acid: " + dlStereochemistryValue );
			}
			if ("yes".equals(aminoAcidEl.getAttributeValue(NATURALENTISOPPOSITE_ATR))){
				invert = !invert;
			}
			
			if (invert) {
				for (Atom atom : atomsWithParities) {
					atom.getAtomParity().setParity(-atom.getAtomParity().getParity());
				}
			}
		}
	}

	void applyDlStereochemistryToCarbohydrate(Element carbohydrateEl, String dlStereochemistryValue) throws ComponentGenerationException {
		Fragment carbohydrate = carbohydrateEl.getFrag();
		List<Atom> atomList = carbohydrate.getAtomList();
		List<Atom> atomsWithParities = new ArrayList<Atom>();
		for (Atom atom : atomList) {
			if (atom.getAtomParity()!=null){
				atomsWithParities.add(atom);
			}
		}
		if (atomsWithParities.isEmpty()){
			throw new ComponentGenerationException("D/L stereochemistry :" + dlStereochemistryValue + " found before achiral carbohydrate");//sounds like a vocab bug...
		}
		
		if (dlStereochemistryValue.equals("dl")){
			for (Atom atom : atomsWithParities) {
				atom.setAtomParity(null);
			}
		}
		else{
			boolean invert;
			if (dlStereochemistryValue.equals("d") || dlStereochemistryValue.equals("dg")){
				invert = false;
			} else if (dlStereochemistryValue.equals("l") || dlStereochemistryValue.equals("lg")){
				invert = true;
			} else{
				throw new ComponentGenerationException("Unexpected value for D/L stereochemistry found before carbohydrate: " + dlStereochemistryValue );
			}
			if ("yes".equals(carbohydrateEl.getAttributeValue(NATURALENTISOPPOSITE_ATR))){
				invert = !invert;
			}

			if (invert) {
				for (Atom atom : atomsWithParities) {
					atom.getAtomParity().setParity(-atom.getAtomParity().getParity());
				}
			}
		}
	}

	static void applyDlStereochemistryToCarbohydrateConfigurationalPrefix(Element elementToApplyTo, String dlStereochemistryValue) throws ComponentGenerationException {
		if (dlStereochemistryValue.equals("d") || dlStereochemistryValue.equals("dg")){
			//do nothing, D- is implicit
		}
		else if (dlStereochemistryValue.equals("l") || dlStereochemistryValue.equals("lg")){
			String[] values = MATCH_SLASH.split(elementToApplyTo.getAttributeValue(VALUE_ATR), -1);
			StringBuilder sb = new StringBuilder();
			for (String value : values) {
				if (value.equals("r")){
					sb.append("l");
				}
				else if (value.equals("l")){
					sb.append("r");
				}
				else{
					throw new RuntimeException("OPSIN Bug: Invalid carbohydrate prefix value: " + elementToApplyTo.getAttributeValue(VALUE_ATR));
				}
				sb.append("/");
			}
			String newVal = sb.toString().substring(0, sb.length()-1);
			elementToApplyTo.getAttribute(VALUE_ATR).setValue(newVal);
		}
		else  if (dlStereochemistryValue.equals("dl")){
			String[] values = MATCH_SLASH.split(elementToApplyTo.getAttributeValue(VALUE_ATR));
			String newVal = "?" + StringTools.multiplyString("/?", values.length-1);
			elementToApplyTo.getAttribute(VALUE_ATR).setValue(newVal);
		}
		else{
			throw new ComponentGenerationException("Unexpected value for D/L stereochemistry found before carbohydrate prefix: " + dlStereochemistryValue );
		}
	}

	/**
	 * Cyclises carbohydrates and regularises their suffixes
	 * @param subOrRoot
	 * @throws StructureBuildingException
	 */
	private void processCarbohydrates(Element subOrRoot) throws StructureBuildingException {
		List<Element> carbohydrates = OpsinTools.getChildElementsWithTagNameAndAttribute(subOrRoot, GROUP_EL, TYPE_ATR, CARBOHYDRATE_TYPE_VAL);
		for (Element carbohydrate : carbohydrates) {
			String subtype = carbohydrate.getAttributeValue(SUBTYPE_ATR);
			boolean isAldose;
			if (CARBOHYDRATESTEMKETOSE_SUBTYPE_VAL.equals(subtype)){
				isAldose = false;
			}
			else if (CARBOHYDRATESTEMALDOSE_SUBTYPE_VAL.equals(subtype) || SYSTEMATICCARBOHYDRATESTEMALDOSE_SUBTYPE_VAL.equals(subtype)){
				isAldose = true;
			}
			else{
				//trivial carbohydrates don't have suffixes
				continue;
			}
			boolean cyclisationPerformed = false;
			Fragment carbohydrateFrag = carbohydrate.getFrag();
			Attribute anomericId = carbohydrate.getAttribute(SUFFIXAPPLIESTO_ATR);
			if (anomericId == null){
				throw new StructureBuildingException("OPSIN bug: Missing suffixAppliesTo on: " + carbohydrate.getValue());
			}
			Atom potentialCarbonyl = carbohydrateFrag.getAtomByID(carbohydrateFrag.getIdOfFirstAtom() + Integer.parseInt(anomericId.getValue()) -1);
			if (potentialCarbonyl == null){
				throw new StructureBuildingException("OPSIN bug: " +  anomericId.getValue() + " did not point to an atom on: " + carbohydrate.getValue());
			}
			carbohydrate.removeAttribute(anomericId);
			
			Element nextSibling = OpsinTools.getNextSibling(carbohydrate);
			while (nextSibling !=null){
				Element nextNextSibling = OpsinTools.getNextSibling(nextSibling);
				String elName = nextSibling.getName();
				if (elName.equals(SUFFIX_EL)){
					Element suffix = nextSibling;
					String value = suffix.getAttributeValue(VALUE_ATR);
					if (value.equals("dialdose") || value.equals("aric acid") || value.equals("arate")){
						if (!isAldose){
							throw new StructureBuildingException(value + " may only be used with aldoses");
						}
						if (cyclisationPerformed){
							throw new StructureBuildingException("OPSIN bug: " +  value + " not expected after carbohydrate cycliser");
						}
						processAldoseDiSuffix(value, carbohydrate, potentialCarbonyl);
						suffix.detach();
					}
					else if (value.startsWith("uron")){
						//strictly these are also aldose di suffixes but in practice they are also used on ketoses
						suffix.addAttribute(new Attribute(LOCANT_ATR, String.valueOf(carbohydrateFrag.getChainLength())));			
					}
					else if (!cyclisationPerformed && (value.equals("ulose") || value.equals("osulose"))){
						if (value.equals("ulose")){
							isAldose = false;
							if (SYSTEMATICCARBOHYDRATESTEMALDOSE_SUBTYPE_VAL.equals(subtype)){
								carbohydrate.getAttribute(SUBTYPE_ATR).setValue(SYSTEMATICCARBOHYDRATESTEMKETOSE_SUBTYPE_VAL);
							}
						}
						potentialCarbonyl = processUloseSuffix(carbohydrate, suffix, potentialCarbonyl);
						suffix.detach();
					}
					else if (value.equals("itol") || value.equals("yl") || value.equals("glycoside")){
						suffix.addAttribute(new Attribute(LOCANT_ATR, potentialCarbonyl.getFirstLocant()));
						if (value.equals("glycoside") && OpsinTools.getParentWordRule(subOrRoot).getAttributeValue(WORDRULE_ATR).equals(WordRule.simple.toString())){
							throw new StructureBuildingException("A glycoside requires a space seperated substituent e.g. methyl alpha-D-glucopyranoside");
						}
					}
				}
				else if (elName.equals(CARBOHYDRATERINGSIZE_EL)){
					if (cyclisationPerformed){
						throw new StructureBuildingException("OPSIN bug: Carbohydate cyclised twice!");
					}
					Element ringSize = nextSibling;
					cycliseCarbohydrate(carbohydrate, ringSize, potentialCarbonyl);
					ringSize.detach();
					cyclisationPerformed = true;
				}
				else if (!elName.equals(LOCANT_EL) &&
						!elName.equals(MULTIPLIER_ATR) &&
						!elName.equals(UNSATURATOR_EL) &&
						!elName.equals(COLONORSEMICOLONDELIMITEDLOCANT_EL)){
					break;
				}
				nextSibling = nextNextSibling;
			}
			if (!cyclisationPerformed){
				applyUnspecifiedRingSizeCyclisationIfPresent(carbohydrate, potentialCarbonyl);
			}
		}
	}

	private void applyUnspecifiedRingSizeCyclisationIfPresent(Element group, Atom potentialCarbonyl) throws StructureBuildingException {
		boolean cyclise = false;
		Element possibleYl = OpsinTools.getNextSibling(group);
		if (possibleYl != null && possibleYl.getName().equals(SUFFIX_EL) && possibleYl.getValue().equals("yl")){
			cyclise = true;
		}
		else{
			Element alphaOrBetaLocantEl = OpsinTools.getPreviousSiblingIgnoringCertainElements(group, new String[]{STEREOCHEMISTRY_EL});
			if (alphaOrBetaLocantEl != null && alphaOrBetaLocantEl.getName().equals(LOCANT_EL) ){
				String value = alphaOrBetaLocantEl.getValue();
				if (value.equals("alpha") || value.equals("beta") || value.equals("alpha,beta") || value.equals("beta,alpha")){
					cyclise = true;
				}
			}
		}
		if (cyclise) {
			Element ringSize = new TokenEl(CARBOHYDRATERINGSIZE_EL);
			String sugarStem = group.getValue();
			if (group.getFrag().hasLocant("5") && !sugarStem.equals("rib") && !sugarStem.equals("fruct")){
				ringSize.addAttribute(new Attribute(VALUE_ATR, "6"));
			}
			else{
				ringSize.addAttribute(new Attribute(VALUE_ATR, "5"));
			}
			OpsinTools.insertAfter(group, ringSize);
			cycliseCarbohydrate(group, ringSize, potentialCarbonyl);
			ringSize.detach();
		}
	}

	/**
	 * Indicates that the compound is a ketose.
	 * This may take the form of replacement of the aldose functionality with ketose functionality or the additon of ketose functionality
	 * The carbonyl may be subsequently used in cyclisation e.g. non-2-ulopyranose
	 * A potentialcarbonyl is returned
	 * @param group
	 * @param suffix
	 * @param potentialCarbonyl 
	 * @return 
	 * @throws StructureBuildingException 
	 */
	private Atom processUloseSuffix(Element group, Element suffix, Atom potentialCarbonyl) throws StructureBuildingException {
		List<String> locantsToConvertToKetones = new ArrayList<String>();
		Element potentialLocantOrMultiplier = OpsinTools.getPreviousSibling(suffix);
		if (potentialLocantOrMultiplier.getName().equals(MULTIPLIER_ATR)){
			int multVal = Integer.parseInt(potentialLocantOrMultiplier.getAttributeValue(VALUE_ATR));
			Element locant = OpsinTools.getPreviousSibling(potentialLocantOrMultiplier);
			if (locant != null && locant.getName().equals(LOCANT_EL)){
				String[] locantStrs = MATCH_COMMA.split(locant.getValue());
				if (locantStrs.length != multVal) {
					throw new StructureBuildingException("Mismatch between locant and multiplier counts (" + locantStrs.length + " and " + multVal + "):" + locant.getValue());
				}
				Collections.addAll(locantsToConvertToKetones, locantStrs);
				locant.detach();
			}
			else{
				for (int i = 0; i < multVal; i++) {
					locantsToConvertToKetones.add(String.valueOf(i + 2));
				}
			}
			potentialLocantOrMultiplier.detach();
		}
		else {
			Element locant = potentialLocantOrMultiplier;
			if (!locant.getName().equals(LOCANT_EL)){
				locant = OpsinTools.getPreviousSibling(group);
			}
			if (locant !=null && locant.getName().equals(LOCANT_EL)){
				String locantStr = locant.getValue();
				if (MATCH_COMMA.split(locantStr).length==1){
					locantsToConvertToKetones.add(locantStr);
				}
				else{
					throw new StructureBuildingException("Incorrect number of locants for ul suffix: " + locantStr);
				}
				locant.detach();
			}
			else{
				locantsToConvertToKetones.add("2");
			}
		}
		Fragment frag = group.getFrag();
		if (suffix.getAttributeValue(VALUE_ATR).equals("ulose")) {//convert aldose to ketose
			Atom aldehydeAtom = potentialCarbonyl;
			boolean foundBond = false;
			for (Bond bond : aldehydeAtom.getBonds()) {
				if (bond.getOrder() ==2){
					Atom otherAtom = bond.getOtherAtom(aldehydeAtom);
					if (otherAtom.getElement() == ChemEl.O && otherAtom.getCharge()==0 && otherAtom.getBondCount()==1){
						bond.setOrder(1);
						foundBond = true;
						break;
					}
				}
			}
			if (!foundBond){
				throw new StructureBuildingException("OPSIN bug: Unable to convert aldose to ketose");
			}
			Atom backboneAtom = frag.getAtomByLocantOrThrow(locantsToConvertToKetones.get(0));
			potentialCarbonyl = backboneAtom;
		}
		for (String locantStr : locantsToConvertToKetones) {
			Atom backboneAtom = frag.getAtomByLocantOrThrow(locantStr);
			boolean foundBond = false;
			for (Bond bond : backboneAtom.getBonds()) {
				if (bond.getOrder() ==1){
					Atom otherAtom = bond.getOtherAtom(backboneAtom);
					if (otherAtom.getElement() == ChemEl.O && otherAtom.getCharge()==0 && otherAtom.getBondCount()==1){
						bond.setOrder(2);
						foundBond = true;
						break;
					}
				}
			}
			if (!foundBond){
				throw new StructureBuildingException("Failed to find hydroxy group at position:" + locantStr);
			}
			backboneAtom.setAtomParity(null);
		}
		return potentialCarbonyl;
	}

	/**
	 * Cyclises carbohydrate configuration prefixes according to the ring size indicator
	 * Alpha/beta stereochemistry is then applied if present
	 * @param carbohydrateGroup
	 * @param ringSize
	 * @param potentialCarbonyl 
	 * @throws StructureBuildingException
	 */
	private void cycliseCarbohydrate(Element carbohydrateGroup, Element ringSize, Atom potentialCarbonyl) throws StructureBuildingException {
		Fragment frag = carbohydrateGroup.getFrag();
		String ringSizeVal = ringSize.getAttributeValue(VALUE_ATR);
		Element potentialLocant = OpsinTools.getPreviousSibling(ringSize);
		Atom carbonylCarbon = null;
		Atom atomToJoinWith = null;
		if (potentialLocant.getName().equals(LOCANT_EL)){
			String[] locants = MATCH_COMMA.split(potentialLocant.getValue());
			if (locants.length != 2){
				throw new StructureBuildingException("Expected 2 locants in front of sugar ring size specifier but found: " + potentialLocant.getValue());
			}
			try{
				int firstLocant = Integer.parseInt(locants[0]);
				int secondLocant = Integer.parseInt(locants[1]);
				if (Math.abs(secondLocant - firstLocant) != (Integer.parseInt(ringSizeVal) -2)){
					throw new StructureBuildingException("Mismatch between ring size: " +  ringSizeVal + " and ring size specified by locants: " + (Math.abs(secondLocant - firstLocant) + 2) );
				}
			}
			catch (NumberFormatException e){
				throw new StructureBuildingException("Locants for ring should be numeric but were: " + potentialLocant.getValue());
			}
			carbonylCarbon = frag.getAtomByLocantOrThrow(locants[0]);
			atomToJoinWith = frag.getAtomByLocantOrThrow("O" + locants[1]);
			potentialLocant.detach();
		}
		
		if (carbonylCarbon == null){
			carbonylCarbon = potentialCarbonyl;
			if (carbonylCarbon ==null){
				throw new RuntimeException("OPSIN bug: Could not find carbonyl carbon in carbohydrate");
			}
		}

		for (Bond b: carbonylCarbon.getBonds()) {
			if (b.getOrder()==2){
				b.setOrder(1);
				break;
			}
		}
		int locantOfCarbonyl;
		try{
			locantOfCarbonyl = Integer.parseInt(carbonylCarbon.getFirstLocant());
		}
		catch (Exception e) {
			throw new RuntimeException("OPSIN bug: Could not determine locant of carbonyl carbon in carbohydrate", e);
		}
		if (atomToJoinWith ==null){
			String locantToJoinWith = String.valueOf(locantOfCarbonyl + Integer.parseInt(ringSizeVal) -2);
			atomToJoinWith =frag.getAtomByLocant("O" +locantToJoinWith);
			if (atomToJoinWith ==null){
				throw new StructureBuildingException("Carbohydrate was not an inappropriate length to form a ring of size: " + ringSizeVal);
			}
		}
		state.fragManager.createBond(carbonylCarbon, atomToJoinWith, 1);
		CycleDetector.assignWhetherAtomsAreInCycles(frag);
		Element alphaOrBetaLocantEl = OpsinTools.getPreviousSiblingIgnoringCertainElements(carbohydrateGroup, new String[]{STEREOCHEMISTRY_EL});
		if (alphaOrBetaLocantEl !=null && alphaOrBetaLocantEl.getName().equals(LOCANT_EL)){
			Element stereoPrefixAfterAlphaBeta = OpsinTools.getNextSibling(alphaOrBetaLocantEl);
			Atom anomericReferenceAtom = getAnomericReferenceAtom(frag);
			if (anomericReferenceAtom ==null){
				throw new RuntimeException("OPSIN bug: Unable to determine anomeric reference atom in: " +carbohydrateGroup.getValue());
			}
			applyAnomerStereochemistryIfPresent(alphaOrBetaLocantEl, carbonylCarbon, anomericReferenceAtom);
			if (carbonylCarbon.getAtomParity() !=null && (SYSTEMATICCARBOHYDRATESTEMALDOSE_SUBTYPE_VAL.equals(carbohydrateGroup.getAttributeValue(SUBTYPE_ATR)) ||
					SYSTEMATICCARBOHYDRATESTEMKETOSE_SUBTYPE_VAL.equals(carbohydrateGroup.getAttributeValue(SUBTYPE_ATR)))){
				//systematic chains only have their stereochemistry defined after structure building to account for the fact that some stereocentres may be removed
				//hence inspect the stereoprefix to see if it is L and flip if so
				String val = stereoPrefixAfterAlphaBeta.getAttributeValue(VALUE_ATR);
				if (val.substring(val.length() -1 , val.length()).equals("l")){//"r" if D, "l" if L
					//flip if L
					AtomParity atomParity = carbonylCarbon.getAtomParity();
					atomParity.setParity(-atomParity.getParity());
				}
			}
		}
		carbonylCarbon.setProperty(Atom.ISANOMERIC, true);
	}

	private void processAldoseDiSuffix(String suffixValue, Element group, Atom aldehydeAtom) throws StructureBuildingException {
		Fragment frag = group.getFrag();
		Atom alcoholAtom = frag.getAtomByLocantOrThrow(String.valueOf(frag.getChainLength()));
		
		if (suffixValue.equals("aric acid") || suffixValue.equals("arate")){
			removeTerminalOxygen(alcoholAtom, 1);
			Fragment f = state.fragManager.buildSMILES("O", group, NONE_LABELS_VAL);
			state.fragManager.incorporateFragment(f, f.getFirstAtom(), frag, alcoholAtom, 2);
			
			f = state.fragManager.buildSMILES("O", group, NONE_LABELS_VAL);
			Atom hydroxyAtom = f.getFirstAtom();
			if (suffixValue.equals("arate")){
				hydroxyAtom.addChargeAndProtons(-1, -1);
			}
			state.fragManager.incorporateFragment(f, f.getFirstAtom(), frag, alcoholAtom, 1);
			frag.addFunctionalAtom(hydroxyAtom);
			
			f = state.fragManager.buildSMILES("O", group, NONE_LABELS_VAL);
			hydroxyAtom = f.getFirstAtom();
			if (suffixValue.equals("arate")){
				hydroxyAtom.addChargeAndProtons(-1, -1);
			}
			state.fragManager.incorporateFragment(f, f.getFirstAtom(), frag, aldehydeAtom, 1);
			frag.addFunctionalAtom(hydroxyAtom);
		}

		else if (suffixValue.equals("dialdose")){
			removeTerminalOxygen(alcoholAtom, 1);
			Fragment f = state.fragManager.buildSMILES("O", group, NONE_LABELS_VAL);
			state.fragManager.incorporateFragment(f, f.getFirstAtom(), frag, alcoholAtom, 2);
		}
		else{
			throw new IllegalArgumentException("OPSIN Bug: Unexpected suffix value: " + suffixValue);
		}		
	}

	/**
	 * Gets the configurationalAtom currently i.e. the defined stereocentre with the highest locant
	 * @param frag
	 * @return
	 */
	private Atom getAnomericReferenceAtom(Fragment frag){
		List<Atom> atomList = frag.getAtomList();
		int highestLocantfound = Integer.MIN_VALUE;
		Atom configurationalAtom = null;
		for (Atom a : atomList) {
			if (a.getAtomParity()==null){
				continue;
			}
			try{
				String locant = a.getFirstLocant();
				int intVal = Integer.parseInt(locant);
				if (intVal > highestLocantfound){
					highestLocantfound = intVal;
					configurationalAtom = a;
				}
			}
			catch (Exception e) {
				//may throw null pointer exceptions or number format exceptions
			}
		}
		return configurationalAtom;
	}

	private void applyAnomerStereochemistryIfPresent(Element alphaOrBetaLocantEl, Atom anomericAtom, Atom anomericReferenceAtom) {
		String value = alphaOrBetaLocantEl.getValue();
		if (value.equals("alpha") || value.equals("beta")){
			Atom[] referenceAtomRefs4 = getDeterministicAtomRefs4ForReferenceAtom(anomericReferenceAtom);
			boolean flip = StereochemistryHandler.checkEquivalencyOfAtomsRefs4AndParity(referenceAtomRefs4, 1, anomericReferenceAtom.getAtomParity().getAtomRefs4(), anomericReferenceAtom.getAtomParity().getParity());
			Atom[] atomRefs4 = getDeterministicAtomRefs4ForAnomericAtom(anomericAtom);
			if (flip){
				if (value.equals("alpha")){
					anomericAtom.setAtomParity(atomRefs4, 1);
				}
				else{
					anomericAtom.setAtomParity(atomRefs4, -1);
				}
			}
			else{
				if (value.equals("alpha")){
					anomericAtom.setAtomParity(atomRefs4, -1);
				}
				else{
					anomericAtom.setAtomParity(atomRefs4, 1);
				}
			}
			alphaOrBetaLocantEl.detach();
		}
		else if (value.equals("alpha,beta") || value.equals("beta,alpha")){
			//unspecified stereochemistry
			alphaOrBetaLocantEl.detach();
		}
	}

	private Atom[] getDeterministicAtomRefs4ForReferenceAtom(Atom referenceAtom) {
		List<Atom> neighbours = referenceAtom.getAtomNeighbours();
		if (neighbours.size()!=3){
			throw new RuntimeException("OPSIN bug: Unexpected number of atoms connected to anomeric reference atom of carbohydrate");
		}
		String nextLowestLocant = String.valueOf(Integer.parseInt(referenceAtom.getFirstLocant()) -1);
		Atom[] atomRefs4 = new Atom[4];
		for (Atom neighbour : neighbours) {
			if (neighbour.getElement() == ChemEl.O) {
				atomRefs4[0] = neighbour;
			}
			else if (neighbour.getElement() == ChemEl.C) {
				if (neighbour.getFirstLocant().equals(nextLowestLocant)){
					atomRefs4[1] = neighbour;
				}
				else {
					atomRefs4[2] = neighbour;
				}
			}
			else{
				throw new RuntimeException("OPSIN bug: Unexpected atom element type connected to for anomeric reference atom");
			}
		}
		atomRefs4[3] = AtomParity.hydrogen;
		for (Atom atom : atomRefs4) {
			if (atom ==null){
				throw new RuntimeException("OPSIN bug: Unable to determine atomRefs4 for anomeric reference atom");
			}
		}
		return atomRefs4;
	}

	private Atom[] getDeterministicAtomRefs4ForAnomericAtom(Atom anomericAtom) {
		List<Atom> neighbours = anomericAtom.getAtomNeighbours();
		if (neighbours.size()!=3 && neighbours.size()!=4){
			throw new RuntimeException("OPSIN bug: Unexpected number of atoms connected to anomeric atom of carbohydrate");
		}
		Atom[] atomRefs4 = new Atom[4];
		for (Atom neighbour : neighbours) {
			if (neighbour.getElement() == ChemEl.C){
				if (neighbour.getAtomIsInACycle()){
					atomRefs4[0] = neighbour;
				}
				else{
					atomRefs4[3] = neighbour;
				}
			}
			else if (neighbour.getElement() == ChemEl.O){
				int incomingVal =neighbour.getIncomingValency();
				if (incomingVal ==1){
					atomRefs4[1] = neighbour;
				}
				else if (incomingVal ==2){
					atomRefs4[2] = neighbour;
				}
				else{
					throw new RuntimeException("OPSIN bug: Unexpected valency on oxygen in carbohydrate");
				}
			}
			else{
				throw new RuntimeException("OPSIN bug: Unexpected atom element type connected to anomeric atom of carbohydrate");
			}
		}
		if (atomRefs4[3]==null){
			atomRefs4[3] = AtomParity.hydrogen;
		}
		for (Atom atom : atomRefs4) {
			if (atom ==null){
				throw new RuntimeException("OPSIN bug: Unable to assign anomeric carbon stereochemistry on carbohydrate");
			}
		}
		return atomRefs4;
	}

	/** Look for multipliers, and multiply out suffixes/unsaturators/heteroatoms/hydros.
	 * Locants are assigned if the number of locants matches the multiplier
	 * associated with them. Eg. triol - > ololol.
	 * Note that infix multiplication is handled seperately as resolution of suffixes is required to perform this unambiguously
	 * @param subOrRoot The substituent/root to looks for multipliers in.
	 */
	private void processMultipliers(Element subOrRoot) {
		List<Element> multipliers = subOrRoot.getChildElements(MULTIPLIER_EL);
		for (Element multiplier : multipliers) {
			Element possibleLocant = OpsinTools.getPreviousSibling(multiplier);
			String[] locants = null;
			if (possibleLocant != null){
				String possibleLocantElName = possibleLocant.getName();
				if (possibleLocantElName.equals(LOCANT_EL)){
					locants = MATCH_COMMA.split(possibleLocant.getValue());
				}
				else if (possibleLocantElName.equals(COLONORSEMICOLONDELIMITEDLOCANT_EL)){
					locants = MATCH_COLON.split(StringTools.removeDashIfPresent(possibleLocant.getValue()));
				}
			}
			Element featureToMultiply = OpsinTools.getNextSibling(multiplier);
			String nextName = featureToMultiply.getName();
			if(nextName.equals(UNSATURATOR_EL) ||
					nextName.equals(SUFFIX_EL) ||
					nextName.equals(SUBTRACTIVEPREFIX_EL) ||
					(nextName.equals(HETEROATOM_EL) && !GROUP_TYPE_VAL.equals(multiplier.getAttributeValue(TYPE_ATR))) ||
					nextName.equals(HYDRO_EL)) {
				int mvalue = Integer.parseInt(multiplier.getAttributeValue(VALUE_ATR));
				if (mvalue>1){
					featureToMultiply.addAttribute(new Attribute(MULTIPLIED_ATR, "multiplied"));
				}
				for(int i= mvalue -1; i >=1; i--) {
					Element newElement = featureToMultiply.copy();
					if (locants !=null && locants.length==mvalue){
						newElement.addAttribute(new Attribute(LOCANT_ATR, locants[i]));
					}
					OpsinTools.insertAfter(featureToMultiply, newElement);
				}
				multiplier.detach();
				if (locants !=null && locants.length==mvalue){
					featureToMultiply.addAttribute(new Attribute(LOCANT_ATR, locants[0]));
					possibleLocant.detach();
				}
			}
		}
	}


	/**
	 * Converts group elements that are identified as being conjunctive suffixes to CONJUNCTIVESUFFIXGROUP_EL
	 * and labels them appropriately. Any suffixes that the conjunctive suffix may have are resolved onto it
	 * @param subOrRoot
	 * @param allGroups
	 * @throws ComponentGenerationException
	 * @throws StructureBuildingException
	 */
	private void detectConjunctiveSuffixGroups(Element subOrRoot, List<Element> allGroups) throws ComponentGenerationException, StructureBuildingException {
		List<Element> groups = subOrRoot.getChildElements(GROUP_EL);
		if (groups.size()>1){
			List<Element> conjunctiveGroups = new ArrayList<Element>();
			Element ringGroup =null;
			for (int i = groups.size() -1 ; i >=0; i--) {
				Element group =groups.get(i);
				if (!group.getAttributeValue(TYPE_ATR).equals(RING_TYPE_VAL)){//e.g. the methanol in benzenemethanol.
					conjunctiveGroups.add(group);
				}
				else{
					ringGroup =group;
					break;
				}
			}
			if (conjunctiveGroups.size() ==0){
				return;
			}
			if (ringGroup ==null){
				throw new ComponentGenerationException("OPSIN bug: unable to find ring associated with conjunctive suffix group");
			}
			if (conjunctiveGroups.size()!=1){
				throw new ComponentGenerationException("OPSIN Bug: Two groups exactly should be present at this point when processing conjunctive nomenclature");
			}
			Element primaryConjunctiveGroup =conjunctiveGroups.get(0);
			Fragment primaryConjunctiveFrag = primaryConjunctiveGroup.getFrag();
			//remove all locants
			List<Atom> atomList = primaryConjunctiveFrag.getAtomList();
			for (Atom atom : atomList) {
				atom.clearLocants();
			}
			List<Element> suffixes = new ArrayList<Element>();
			Element possibleSuffix = OpsinTools.getNextSibling(primaryConjunctiveGroup);
			while (possibleSuffix !=null){
				if (possibleSuffix.getName().equals(SUFFIX_EL)){
					suffixes.add(possibleSuffix);
				}
				possibleSuffix = OpsinTools.getNextSibling(possibleSuffix);
			}
			preliminaryProcessSuffixes(primaryConjunctiveGroup, suffixes);
			resolveSuffixes(primaryConjunctiveGroup, suffixes);
			for (Element suffix : suffixes) {
				suffix.detach();
			}
			primaryConjunctiveGroup.setName(CONJUNCTIVESUFFIXGROUP_EL);
			allGroups.remove(primaryConjunctiveGroup);
			
			Element possibleMultiplier = OpsinTools.getPreviousSibling(primaryConjunctiveGroup);
			//label atoms appropriately
			boolean alphaIsPosition1 = atomList.get(0).getIncomingValency() < 3;
			int counter =0;
			for (int i = (alphaIsPosition1 ? 0 : 1); i < atomList.size(); i++) {
				Atom a = atomList.get(i);
				if (counter==0){
					a.addLocant("alpha");
				}
				else if (counter==1){
					a.addLocant("beta");
				}
				else if (counter==2){
					a.addLocant("gamma");
				}
				else if (counter==3){
					a.addLocant("delta");
				}
				else if (counter==4){
					a.addLocant("epsilon");
				}
				else if (counter==5){
					a.addLocant("zeta");
				}
				else if (counter==6){
					a.addLocant("eta");
				}
				counter++;
			}
			if (MULTIPLIER_EL.equals(possibleMultiplier.getName())){
				int multiplier = Integer.parseInt(possibleMultiplier.getAttributeValue(VALUE_ATR));
				for (int i = 1; i < multiplier; i++) {
					Element conjunctiveSuffixGroup = primaryConjunctiveGroup.copy();
					Fragment newFragment = state.fragManager.copyAndRelabelFragment(primaryConjunctiveFrag, i);
					newFragment.setTokenEl(conjunctiveSuffixGroup);
					conjunctiveSuffixGroup.setFrag(newFragment);
					conjunctiveGroups.add(conjunctiveSuffixGroup);
					OpsinTools.insertAfter(primaryConjunctiveGroup, conjunctiveSuffixGroup);
				}
				Element possibleLocant = OpsinTools.getPreviousSibling(possibleMultiplier);
				possibleMultiplier.detach();
				if (possibleLocant.getName().equals(LOCANT_EL)){
					String[] locants = MATCH_COMMA.split(possibleLocant.getValue());
					if (locants.length!=multiplier){
						throw new ComponentGenerationException("mismatch between number of locants and multiplier in conjunctive nomenclature routine");
					}
					for (int i = 0; i < locants.length; i++) {
						conjunctiveGroups.get(i).addAttribute(new Attribute(LOCANT_ATR, locants[i]));
					}
					possibleLocant.detach();
				}
			}
		}
	}


	/** Match each locant to the next applicable "feature". Assumes that processLocants
	 * has done a good job and rejected cases where no match can be made.
	 * Handles cases where the locant is next to the feature it refers to
	 *
	 * @param subOrRoot The substituent/root to look for locants in.
	 * @throws ComponentGenerationException
	 */
	private void matchLocantsToDirectFeatures(Element subOrRoot) throws ComponentGenerationException {
		List<Element> locants =  subOrRoot.getChildElements(LOCANT_EL);
		List<Element> groups = subOrRoot.getChildElements(GROUP_EL);
		for (Element group : groups) {
			if (group.getAttributeValue(SUBTYPE_ATR).equals(HANTZSCHWIDMAN_SUBTYPE_VAL)){//handle Hantzch-widman systems
				if (group.getAttribute(ADDBOND_ATR)!=null){//special case for partunsatring
					//exception for where a locant is supposed to indicate the location of a double bond...
					List<Element> deltas = subOrRoot.getChildElements(DELTA_EL);
					if (deltas.size()==0){
						Element delta =new TokenEl(DELTA_EL);
						Element appropriateLocant = OpsinTools.getPreviousSiblingIgnoringCertainElements(group, new String[]{HETEROATOM_EL, MULTIPLIER_EL});
						if (appropriateLocant !=null && appropriateLocant.getName().equals(LOCANT_EL) && MATCH_COMMA.split(appropriateLocant.getValue()).length == 1){
							delta.setValue(appropriateLocant.getValue());
							OpsinTools.insertBefore(appropriateLocant, delta);
							appropriateLocant.detach();
							locants.remove(appropriateLocant);
						}
						else{
							delta.setValue("");
							subOrRoot.insertChild(delta, 0);//no obvious attempt to set double bond position, potentially ambiguous, valency will be used to choose later
						}
					}
				}
				if (locants.size()>0 ){
					Element locantBeforeHWSystem = null;
					List<Element> heteroAtoms = new ArrayList<Element>();
					int indexOfGroup = subOrRoot.indexOf(group);
					for (int j = indexOfGroup -1; j >= 0; j--) {
						String elName = subOrRoot.getChild(j).getName();
						if (elName.equals(LOCANT_EL)){
							locantBeforeHWSystem = subOrRoot.getChild(j);
							break;
						}
						else if(elName.equals(HETEROATOM_EL)){
							Element heteroAtom = subOrRoot.getChild(j);
							heteroAtoms.add(heteroAtom);
							if (heteroAtom.getAttribute(LOCANT_ATR)!=null){//locants already assigned, assumedly by process multipliers
								break;
							}
						}
						else{
							break;
						}
					}
					Collections.reverse(heteroAtoms);
					if (locantBeforeHWSystem !=null){
						String[] locantValues = MATCH_COMMA.split(locantBeforeHWSystem.getValue());
						//detect a solitary locant in front of a HW system and prevent it being assigned.
						//something like 1-aziridin-1-yl never means the N is at position 1 as it is at position 1 by convention
						//this special case is not applied to pseudo HW like systems e.g. [1]oxacyclotetradecine
						if (locantValues.length ==1 && group.getFrag().getAtomCount() <=10){
							locants.remove(locantBeforeHWSystem);//don't assign this locant
						}
						else {
							if (locantValues.length == heteroAtoms.size()){
								for (int j = 0; j < locantValues.length; j++) {
									String locantValue = locantValues[j];
									heteroAtoms.get(j).addAttribute(new Attribute(LOCANT_ATR, locantValue));
								}
								locantBeforeHWSystem.detach();
								locants.remove(locantBeforeHWSystem);
							}
							else if (heteroAtoms.size()>1){
								throw new ComponentGenerationException("Mismatch between number of locants and HW heteroatoms");
							}
						}
					}
				}
			}
		}
		assignSingleLocantsToAdjacentFeatures(locants);
	}


	/**
	 * Looks for a suffix/suffix/heteroatom/hydro element adjacent to the given locant
	 * and if the locant element describes just 1 locant asssigns it
	 * @param locants
	 */
	private void assignSingleLocantsToAdjacentFeatures(List<Element> locants) {
		for (Element locant : locants) {
			String[] locantValues = MATCH_COMMA.split(locant.getValue());
			Element referent = OpsinTools.getNextSibling(locant);
			if (referent != null && locantValues.length == 1){
				String refName = referent.getName();
				//Only assigning locants to elements that were not created by a multiplier
				if(referent.getAttribute(LOCANT_ATR) == null && referent.getAttribute(MULTIPLIED_ATR) == null && (refName.equals(UNSATURATOR_EL) ||
						refName.equals(SUFFIX_EL) ||
						refName.equals(HETEROATOM_EL) ||
						refName.equals(CONJUNCTIVESUFFIXGROUP_EL) ||
						refName.equals(SUBTRACTIVEPREFIX_EL) ||
						(refName.equals(HYDRO_EL) && !referent.getValue().startsWith("per") ))) {//not perhydro
					referent.addAttribute(new Attribute(LOCANT_ATR, locantValues[0]));
					locant.detach();
				}
			}
		}
	}


	/**
	 * Handles suffixes, passes them to resolveGroupAddingSuffixes.
	 * Processes the suffixAppliesTo command which multiplies a suffix and attaches the suffixes to the atoms described by the given IDs
	 * @param group
	 * @param suffixes 
	 * @throws ComponentGenerationException
	 * @throws StructureBuildingException
	 */
	private void preliminaryProcessSuffixes(Element group, List<Element> suffixes) throws ComponentGenerationException, StructureBuildingException{
		Fragment suffixableFragment = group.getFrag();

		if (group.getAttribute(SUFFIXAPPLIESTO_ATR)!=null){//typically a trivial polyAcid or aminoAcid
			processSuffixAppliesTo(group, suffixes,suffixableFragment);
		}
		else{
			for (Element suffix : suffixes) {
				if (suffix.getAttribute(ADDITIONALVALUE_ATR)!=null){
					throw new ComponentGenerationException("suffix: " + suffix.getValue() + " used on an inappropriate group");
				}
			}
		}
		applyDefaultLocantsToSuffixesIfApplicable(group, suffixableFragment);

		List<Fragment> suffixFragments =resolveGroupAddingSuffixes(suffixes, suffixableFragment);
		state.xmlSuffixMap.put(group, suffixFragments);
		boolean suffixesResolved =false;
		if (group.getAttributeValue(TYPE_ATR).equals(CHALCOGENACIDSTEM_TYPE_VAL)){//merge the suffix into the chalcogen acid stem e.g sulfonoate needs to be one fragment for infix replacement
			resolveSuffixes(group, suffixes);
			suffixesResolved =true;
		}
		processSuffixPrefixes(suffixes);//e.g. carbox amide
		functionalReplacement.processInfixFunctionalReplacementNomenclature(suffixes, suffixFragments);
		processRemovalOfHydroxyGroupsRules(suffixes, suffixableFragment);

		if (group.getValue().equals("oxal")){//oxalic acid is treated as a non carboxylic acid for the purposes of functional replacment. See P-65.2.3
			resolveSuffixes(group, suffixes);
			group.getAttribute(TYPE_ATR).setValue(NONCARBOXYLICACID_TYPE_VAL);
			suffixesResolved =true;
		}
		if (suffixesResolved){
			//suffixes have already been resolved so need to be detached to avoid being passed to resolveSuffixes later
			for (int i = suffixes.size() -1; i>=0; i--) {
				Element suffix =suffixes.remove(i);
				suffix.detach();
			}
		}
		if (group.getAttribute(NUMBEROFFUNCTIONALATOMSTOREMOVE_ATR)!=null){
			int numberToRemove = Integer.parseInt(group.getAttributeValue(NUMBEROFFUNCTIONALATOMSTOREMOVE_ATR));
			if (numberToRemove > suffixableFragment.getFunctionalAtomCount()){
				throw new ComponentGenerationException("Too many hydrogen for the number of positions on non carboxylic acid");
			}
			for (int i = 0; i< numberToRemove; i++) {
				Atom functionalAtom = suffixableFragment.removeFunctionalAtom(0).getAtom();
				functionalAtom.neutraliseCharge();
			}		
		}
	}

	private void applyDefaultLocantsToSuffixesIfApplicable(Element group, Fragment suffixableFragment)  {
		String defaultLocantsAtrValue = group.getAttributeValue(SUFFIXAPPLIESTOBYDEFAULT_ATR);
		if (defaultLocantsAtrValue != null){
			String[] suffixInstructions = MATCH_COMMA.split(defaultLocantsAtrValue);
			int firstIdInFragment = suffixableFragment.getIdOfFirstAtom();
			Element suffix = OpsinTools.getNextNonChargeSuffix(group);
			for (String suffixInstruction : suffixInstructions) {
				if (suffix !=null){
					suffix.addAttribute(new Attribute(DEFAULTLOCANTID_ATR, Integer.toString(firstIdInFragment + Integer.parseInt(suffixInstruction) -1)));
					suffix = OpsinTools.getNextNonChargeSuffix(suffix);
				}
			}
		}
	}

	/**
	 * Processes the effects of the suffixAppliesTo attribute
	 * @param group
	 * @param suffixes
	 * @param suffixableFragment
	 * @throws ComponentGenerationException
	 */
	private void processSuffixAppliesTo(Element group, List<Element> suffixes, Fragment suffixableFragment) throws ComponentGenerationException {
		//suffixAppliesTo attribute contains instructions for number/positions of suffix
		//this is of the form comma sepeated ids with the number of ids corresponding to the number of instances of the suffix
		Element suffix =OpsinTools.getNextNonChargeSuffix(group);
		if (suffix ==null){
			if (group.getAttributeValue(TYPE_ATR).equals(ACIDSTEM_TYPE_VAL)){
				throw new ComponentGenerationException("No suffix where suffix was expected");
			}
		}
		else{
			if (suffixes.size()>1 && group.getAttributeValue(TYPE_ATR).equals(ACIDSTEM_TYPE_VAL)){
				throw new ComponentGenerationException("More than one suffix detected on trivial polyAcid. Not believed to be allowed");
			}

			String suffixInstruction =group.getAttributeValue(SUFFIXAPPLIESTO_ATR);
			String[] suffixInstructions = MATCH_COMMA.split(suffixInstruction);
			if (CYCLEFORMER_SUBTYPE_VAL.equals(suffix.getAttributeValue(SUBTYPE_ATR))){
				if (suffixInstructions.length !=2){
					throw new ComponentGenerationException("suffix: " + suffix.getValue() + " used on an inappropriate group");
				}
				suffix.addAttribute(new Attribute(LOCANTID_ATR, suffixInstruction));
				return;
			}
			boolean symmetricSuffixes =true;
			if (suffix.getAttribute(ADDITIONALVALUE_ATR)!=null){//handles amic, aldehydic, anilic and amoyl suffixes properly
				if (suffixInstructions.length < 2){
					throw new ComponentGenerationException("suffix: " + suffix.getValue() + " used on an inappropriate group");
				}
				symmetricSuffixes = false;
			}

			int firstIdInFragment=suffixableFragment.getIdOfFirstAtom();
			if (suffix.getAttribute(LOCANT_ATR)==null){
				suffix.addAttribute(new Attribute(LOCANTID_ATR, Integer.toString(firstIdInFragment + Integer.parseInt(suffixInstructions[0]) -1)));
			}
			for (int i = 1; i < suffixInstructions.length; i++) {
				Element newSuffix = new TokenEl(SUFFIX_EL);
				if (symmetricSuffixes){
					newSuffix.addAttribute(new Attribute(VALUE_ATR, suffix.getAttributeValue(VALUE_ATR)));
					newSuffix.addAttribute(new Attribute(TYPE_ATR,  suffix.getAttributeValue(TYPE_ATR)));
					if (suffix.getAttribute(SUBTYPE_ATR)!=null){
						newSuffix.addAttribute(new Attribute(SUBTYPE_ATR,  suffix.getAttributeValue(SUBTYPE_ATR)));
					}
					if (suffix.getAttribute(INFIX_ATR)!=null && suffix.getAttributeValue(INFIX_ATR).startsWith("=")){//clone infixes that effect double bonds but not single bonds e.g. maleamidate still should have one functional atom
						newSuffix.addAttribute(new Attribute(INFIX_ATR,  suffix.getAttributeValue(INFIX_ATR)));
					}
				}
				else{
					newSuffix.addAttribute(new Attribute(VALUE_ATR, suffix.getAttributeValue(ADDITIONALVALUE_ATR)));
					newSuffix.addAttribute(new Attribute(TYPE_ATR, ROOT_EL));
				}
				newSuffix.addAttribute(new Attribute(LOCANTID_ATR, Integer.toString(firstIdInFragment + Integer.parseInt(suffixInstructions[i]) -1)));
				OpsinTools.insertAfter(suffix, newSuffix);
				suffixes.add(newSuffix);
			}
		}
	}


	/**Processes a suffix and returns any fragment the suffix intends to add to the molecule
	 * @param suffixes The suffix elements for a fragment.
	 * @param frag The fragment to which the suffix will be applied
	 * @return An arrayList containing the generated fragments
	 * @throws StructureBuildingException If the suffixes can't be resolved properly.
	 * @throws ComponentGenerationException
	 */
	private List<Fragment> resolveGroupAddingSuffixes(List<Element> suffixes, Fragment frag) throws StructureBuildingException, ComponentGenerationException {
		List<Fragment> suffixFragments =new ArrayList<Fragment>();
		String groupType = frag.getType();
		String subgroupType = frag.getSubType();

		String suffixTypeToUse =null;
		if (suffixRulesLookup.isGroupTypeWithSpecificSuffixRules(groupType)){
			suffixTypeToUse =groupType;
		}
		else{
			suffixTypeToUse = STANDARDGROUP_TYPE_VAL;
		}

		for (Element suffix : suffixes) {
			String suffixValue = suffix.getAttributeValue(VALUE_ATR);

			boolean cyclic;//needed for addSuffixPrefixIfNonePresentAndCyclic rule
			Atom atomLikelyToBeUsedBySuffix = null;
			
			String locant = suffix.getAttributeValue(LOCANT_ATR);
			String locantId = suffix.getAttributeValue(LOCANTID_ATR);

			if (locant != null && locant.indexOf(',') == -1) {
				atomLikelyToBeUsedBySuffix = frag.getAtomByLocant(locant);
			}
			else if (locantId != null && locantId.indexOf(',') == -1) {
				atomLikelyToBeUsedBySuffix = frag.getAtomByIDOrThrow(Integer.parseInt(locantId));
			}
			if (atomLikelyToBeUsedBySuffix==null){
				//a locant has not been specified
				//also can happen in the cases of things like fused rings where the final numbering is not available so lookup by locant fails (in which case all the atoms will be cyclic anyway)
				atomLikelyToBeUsedBySuffix = frag.getFirstAtom();
			}
			cyclic = atomLikelyToBeUsedBySuffix.getAtomIsInACycle();

			List<SuffixRule> suffixRules = suffixRulesLookup.getSuffixRuleTags(suffixTypeToUse, suffixValue, subgroupType);
			Fragment suffixFrag = null;
			/*
			 * Temp fragments are build for each addGroup rule and then merged into suffixFrag
			 */
			for (SuffixRule suffixRule : suffixRules) {
				switch (suffixRule.getType()) {
				case addgroup:
					String labels = suffixRule.getAttributeValue(SUFFIXRULES_LABELS_ATR);
					if (labels == null) {
						labels = NONE_LABELS_VAL;
					}
					suffixFrag = state.fragManager.buildSMILES(suffixRule.getAttributeValue(SUFFIXRULES_SMILES_ATR), SUFFIX_TYPE_VAL, labels);
					List<Atom> atomList = suffixFrag.getAtomList();
					String functionalIdsAtr = suffixRule.getAttributeValue(SUFFIXRULES_FUNCTIONALIDS_ATR);
					if (functionalIdsAtr != null) {
						String[] relativeIdsOfFunctionalAtoms = MATCH_COMMA.split(functionalIdsAtr);
						for (String relativeId : relativeIdsOfFunctionalAtoms) {
							int atomIndice = Integer.parseInt(relativeId) -1;
							if (atomIndice >=atomList.size()){
								throw new StructureBuildingException("Check suffixRules.xml: Atom requested to have a functionalAtom was not within the suffix fragment");
							}
							suffixFrag.addFunctionalAtom(atomList.get(atomIndice));
						}
					}
					String outIdsAtr = suffixRule.getAttributeValue(SUFFIXRULES_OUTIDS_ATR);
					if (outIdsAtr != null) {
						String[] relativeIdsOfOutAtoms = MATCH_COMMA.split(outIdsAtr);
						for (String relativeId : relativeIdsOfOutAtoms) {
							int atomIndice = Integer.parseInt(relativeId) -1;
							if (atomIndice >=atomList.size()){
								throw new StructureBuildingException("Check suffixRules.xml: Atom requested to have a outAtom was not within the suffix fragment");
							}
							suffixFrag.addOutAtom(atomList.get(atomIndice), 1 , true);
						}
					}
					break;
				case addSuffixPrefixIfNonePresentAndCyclic:
					if (cyclic && suffix.getAttribute(SUFFIXPREFIX_ATR)==null){
						suffix.addAttribute(new Attribute(SUFFIXPREFIX_ATR, suffixRule.getAttributeValue(SUFFIXRULES_SMILES_ATR)));
					}
					break;
				case addFunctionalAtomsToHydroxyGroups:
					if (suffixFrag != null){
						throw new ComponentGenerationException("addFunctionalAtomsToHydroxyGroups is not currently compatable with the addGroup suffix rule");
					}
					addFunctionalAtomsToHydroxyGroups(atomLikelyToBeUsedBySuffix);
					break;
				case chargeHydroxyGroups:
					if (suffixFrag != null){
						throw new ComponentGenerationException("chargeHydroxyGroups is not currently compatable with the addGroup suffix rule");
					}
					chargeHydroxyGroups(atomLikelyToBeUsedBySuffix);
					break;
				case removeTerminalOxygen:
					if (suffixFrag != null){
						throw new ComponentGenerationException("removeTerminalOxygen is not currently compatible with the addGroup suffix rule");
					}
					int bondOrder = Integer.parseInt(suffixRule.getAttributeValue(SUFFIXRULES_ORDER_ATR));
					removeTerminalOxygen(atomLikelyToBeUsedBySuffix, bondOrder);
					break;
				default:
					break;
				}
			}
			if (suffixFrag != null) {
				suffixFragments.add(suffixFrag);
				suffix.setFrag(suffixFrag);
			}
		}
		return suffixFragments;
	}
	
	/**Processes any convertHydroxyGroupsToOutAtoms and convertHydroxyGroupsToPositiveCharge instructions
	 * This is not handled as part of resolveGroupAddingSuffixes as something like carbonochloridoyl involves infix replacement
	 * on a hydroxy that would otherwise actually be removed by this rule!
	 * @param suffixes The suffix elements for a fragment.
	 * @param frag The fragment to which the suffix will be applied
	 * @throws ComponentGenerationException
	 */
	private void processRemovalOfHydroxyGroupsRules(List<Element> suffixes, Fragment frag) throws ComponentGenerationException {
		String groupType = frag.getType();
		String subgroupType = frag.getSubType();
		String suffixTypeToUse =null;
		if (suffixRulesLookup.isGroupTypeWithSpecificSuffixRules(groupType)) {
			suffixTypeToUse =groupType;
		}
		else{
			suffixTypeToUse = STANDARDGROUP_TYPE_VAL;
		}
		for (Element suffix : suffixes) {
			String suffixValue = suffix.getAttributeValue(VALUE_ATR);
			List<SuffixRule> suffixRules = suffixRulesLookup.getSuffixRuleTags(suffixTypeToUse, suffixValue, subgroupType);
			for (SuffixRule suffixRule : suffixRules) {
				SuffixRuleType type =suffixRule.getType();
				if (type == SuffixRuleType.convertHydroxyGroupsToOutAtoms) {
					convertHydroxyGroupsToOutAtoms(frag);
				}
				else if (type == SuffixRuleType.convertHydroxyGroupsToPositiveCharge) {
					convertHydroxyGroupsToPositiveCharge(frag);
				}
			}
		}
	}

	/**
	 * Finds all hydroxy groups connected to a given atom and adds a functionalAtom to each of them
	 * @param atom
	 * @throws StructureBuildingException
	 */
	private void addFunctionalAtomsToHydroxyGroups(Atom atom) throws StructureBuildingException {
		List<Atom> neighbours = atom.getAtomNeighbours();
		for (Atom neighbour : neighbours) {
			if (neighbour.getElement() == ChemEl.O && neighbour.getCharge() == 0 && neighbour.getAtomNeighbours().size() == 1 && atom.getBondToAtomOrThrow(neighbour).getOrder() == 1){
				neighbour.getFrag().addFunctionalAtom(neighbour);
			}
		}
	}

	/**
	 * Finds all hydroxy groups connected to a given atom and makes them negatively charged
	 * @param atom
	 * @throws StructureBuildingException
	 */
	private void chargeHydroxyGroups(Atom atom) throws StructureBuildingException {
		List<Atom> neighbours = atom.getAtomNeighbours();
		for (Atom neighbour : neighbours) {
			if (neighbour.getElement() == ChemEl.O && neighbour.getCharge()==0 && neighbour.getAtomNeighbours().size()==1 && atom.getBondToAtomOrThrow(neighbour).getOrder()==1){
				neighbour.addChargeAndProtons(-1, -1);
			}
		}
	}

	/**
	 * Removes a terminal oxygen from the atom 
	 * An exception is thrown if no suitable oxygen could be found connected to the atom
	 * Note that [N+][O-] is treated as N=O
	 * @param atom
	 * @throws StructureBuildingException
	 */
	private void removeTerminalOxygen(Atom atom, int desiredBondOrder) throws StructureBuildingException {
		//TODO prioritise [N+][O-]
		List<Atom> neighbours = atom.getAtomNeighbours();
		for (Atom neighbour : neighbours) {
			if (neighbour.getElement() == ChemEl.O && neighbour.getAtomNeighbours().size()==1){
				Bond b = atom.getBondToAtomOrThrow(neighbour);
				if (b.getOrder()==desiredBondOrder && neighbour.getCharge()==0){
					FragmentTools.removeTerminalAtom(state, neighbour);
					if (atom.getLambdaConventionValency()!=null){//corrects valency for phosphin/arsin/stibin
						atom.setLambdaConventionValency(atom.getLambdaConventionValency()-desiredBondOrder);
					}
					if (atom.getMinimumValency()!=null){//corrects valency for phosphin/arsin/stibin
						atom.setMinimumValency(atom.getMinimumValency()-desiredBondOrder);
					}
					return;
				}
				else if (neighbour.getCharge() ==-1 && b.getOrder()==1 && desiredBondOrder == 2){
					if (atom.getCharge() ==1 && atom.getElement() == ChemEl.N){
						FragmentTools.removeTerminalAtom(state, neighbour);
						atom.neutraliseCharge();
						return;
					}
				}
			}
		}
		if (desiredBondOrder ==2){
			throw new StructureBuildingException("Double bonded oxygen not found at suffix attachment position. Perhaps a suffix has been used inappropriately");
		}
		else if (desiredBondOrder ==1){
			throw new StructureBuildingException("Hydroxy oxygen not found at suffix attachment position. Perhaps a suffix has been used inappropriately");
		}
		else {
			throw new StructureBuildingException("Suitable oxygen not found at suffix attachment position Perhaps a suffix has been used inappropriately");
		}

	}
	
	/**
	 * Given a fragment removes all hydroxy groups and adds a valency 1 outAtom to the adjacent atom for each hydroxy group
	 * Note that O[OH] is not considered a hydroxy c.f. carbonoperoxoyl
	 * @param frag
	 */
	private void convertHydroxyGroupsToOutAtoms(Fragment frag) {
		List<Atom> atomList = frag.getAtomList();
		for (Atom atom : atomList) {
			if (atom.getElement() == ChemEl.O && atom.getCharge()==0 && atom.getBondCount()==1  &&
					atom.getFirstBond().getOrder()==1 && atom.getOutValency() == 0){
				Atom adjacentAtom = atom.getAtomNeighbours().get(0);
				if (adjacentAtom.getElement() != ChemEl.O){
					state.fragManager.removeAtomAndAssociatedBonds(atom);
					frag.addOutAtom(adjacentAtom, 1, true);
				}
			}
		}
	}

	/**
	 * Given a fragment removes all hydroxy groups and applies ylium to the adjacent atom (+1 charge -1 proton)
	 * Note that O[OH] is not considered a hydroxy
	 * @param frag
	 */
	private void convertHydroxyGroupsToPositiveCharge(Fragment frag) {
		List<Atom> atomList = frag.getAtomList();
		for (Atom atom : atomList) {
			if (atom.getElement() == ChemEl.O && atom.getCharge()==0 && atom.getBondCount()==1  &&
					atom.getFirstBond().getOrder()==1 && atom.getOutValency() == 0){
				Atom adjacentAtom = atom.getAtomNeighbours().get(0);
				if (adjacentAtom.getElement() != ChemEl.O){
					state.fragManager.removeAtomAndAssociatedBonds(atom);
					adjacentAtom.addChargeAndProtons(1, -1);
				}
			}
		}
	}

	/**
	 * Searches for suffix elements with the suffixPrefix attribute set
	 * A suffixPrefix is something like sulfon in sulfonamide. It would in this case take the value S(=O)
	 * @param suffixes
	 * @throws StructureBuildingException
	 */
	private void processSuffixPrefixes(List<Element> suffixes) throws StructureBuildingException{
		for (Element suffix : suffixes) {
			if (suffix.getAttribute(SUFFIXPREFIX_ATR)!=null){
				Fragment suffixPrefixFrag = state.fragManager.buildSMILES(suffix.getAttributeValue(SUFFIXPREFIX_ATR), SUFFIX_TYPE_VAL, NONE_LABELS_VAL);
				addFunctionalAtomsToHydroxyGroups(suffixPrefixFrag.getFirstAtom());
				if (suffix.getValue().endsWith("ate") || suffix.getValue().endsWith("at")){
					chargeHydroxyGroups(suffixPrefixFrag.getFirstAtom());
				}
				Atom firstAtomOfPrefix = suffixPrefixFrag.getFirstAtom();
				firstAtomOfPrefix.addLocant("X");//make sure this atom is not given a locant
				Fragment suffixFrag = suffix.getFrag();
				state.fragManager.incorporateFragment(suffixPrefixFrag, suffixFrag);
				
				//manipulate suffixFrag such that all the bonds to the first atom (the R)  go instead to the first atom of suffixPrefixFrag.
				//Then reconnect the R to that atom
				Atom theR = suffixFrag.getFirstAtom();
				List<Atom> neighbours = theR.getAtomNeighbours();
				for (Atom neighbour : neighbours) {
					Bond b = theR.getBondToAtomOrThrow(neighbour);
					state.fragManager.removeBond(b);
					state.fragManager.createBond(neighbour, firstAtomOfPrefix, b.getOrder());
				}
				state.fragManager.createBond(firstAtomOfPrefix, theR, 1);
			}
		}
	}

	/**
	 * Checks through the groups accesible from the startingElement taking into account brackets (i.e. those that it is feasible that the group of the startingElement could substitute onto).
	 * It is assumed that one does not intentionally locant onto something in a deeper level of bracketting (not implicit bracketing). e.g. 2-propyl(ethyl)ammonia will give prop-2-yl
	 * @param startingElement
	 * @param locant: the locant string to check for the presence of
	 * @return whether the locant was found
	 * @throws StructureBuildingException
	 */
	private boolean checkLocantPresentOnPotentialRoot(Element startingElement, String locant) throws StructureBuildingException {
		boolean foundSibling =false;
		Deque<Element> s = new ArrayDeque<Element>();
		s.add(startingElement);
		boolean doneFirstIteration =false;//check on index only done on first iteration to only get elements with an index greater than the starting element
		while (s.size()>0){
			Element currentElement =s.removeLast();
			Element parent = currentElement.getParent();
			List<Element> siblings = OpsinTools.getChildElementsWithTagNames(parent, new String[]{BRACKET_EL, SUBSTITUENT_EL, ROOT_EL});
			int indexOfCurrentElement =parent.indexOf(currentElement);

			for (Element bracketOrSub : siblings) {
				if (!doneFirstIteration && parent.indexOf(bracketOrSub) <= indexOfCurrentElement){
					continue;
				}
				if (bracketOrSub.getName().equals(BRACKET_EL)){//only want to consider implicit brackets, not proper brackets
					if (bracketOrSub.getAttribute(TYPE_ATR)==null){
						continue;
					}
					s.add(bracketOrSub.getChild(0));
				}
				else{
					Element group = bracketOrSub.getFirstChildElement(GROUP_EL);
					Fragment groupFrag = group.getFrag();
					if (groupFrag.hasLocant(locant)){
						return true;
					}
					List<Fragment> suffixes =state.xmlSuffixMap.get(group);
					if (suffixes!=null){
						for (Fragment suffix : suffixes) {
							if (suffix.hasLocant(locant)){
								return true;
							}
						}
					}
					List<Element> conjunctiveGroups = OpsinTools.getNextSiblingsOfType(group, CONJUNCTIVESUFFIXGROUP_EL);
					for (Element conjunctiveGroup : conjunctiveGroups) {
						if (conjunctiveGroup.getFrag().hasLocant(locant)){
							return true;
						}
					}
				}
				foundSibling =true;
			}
			doneFirstIteration =true;
		}

		if (!foundSibling){//Special case: anything the group could potentially substitute onto is in a bracket. The bracket is checked recursively
			s = new ArrayDeque<Element>();
			s.add(startingElement);
			doneFirstIteration =false;//check on index only done on first iteration to only get elements with an index greater than the starting element
			while (s.size()>0){
				Element currentElement =s.removeLast();
				Element parent = currentElement.getParent();
				List<Element> siblings = OpsinTools.getChildElementsWithTagNames(parent, new String[]{BRACKET_EL, SUBSTITUENT_EL, ROOT_EL});
				int indexOfCurrentElement =parent.indexOf(currentElement);

				for (Element bracketOrSub : siblings) {
					if (!doneFirstIteration && parent.indexOf(bracketOrSub) <= indexOfCurrentElement){
						continue;
					}
					if (bracketOrSub.getName().equals(BRACKET_EL)){
						s.add(bracketOrSub.getChild(0));
					}
					else{
						Element group = bracketOrSub.getFirstChildElement(GROUP_EL);
						Fragment groupFrag = group.getFrag();
						if (groupFrag.hasLocant(locant)){
							return true;
						}
						List<Fragment> suffixes =state.xmlSuffixMap.get(group);
						if (suffixes!=null){
							for (Fragment suffix : suffixes) {
								if (suffix.hasLocant(locant)){
									return true;
								}
							}
						}
						List<Element> conjunctiveGroups = OpsinTools.getNextSiblingsOfType(group, CONJUNCTIVESUFFIXGROUP_EL);
						for (Element conjunctiveGroup : conjunctiveGroups) {
							if (conjunctiveGroup.getFrag().hasLocant(locant)){
								return true;
							}
						}
					}
				}
				doneFirstIteration =true;
			}
		}

		return false;
	}

	/** Handles special cases in IUPAC nomenclature that are most elegantly solved by modification of the fragment
	 *  Also sets the default in atom for alkanes so that say methylethyl is prop-2-yl rather than propyl
	 * @param groups
	 * @throws StructureBuildingException
	 * @throws ComponentGenerationException 
	 */
	private void handleGroupIrregularities(List<Element> groups) throws StructureBuildingException, ComponentGenerationException{
		for (Element group : groups) {
			String groupValue =group.getValue();
			if (groupValue.equals("porphyrin")|| groupValue.equals("porphin")){
				List<Element> hydrogenAddingEls = group.getParent().getChildElements(INDICATEDHYDROGEN_EL);
				boolean implicitHydrogenExplicitlySet =false;
				for (Element hydrogenAddingEl : hydrogenAddingEls) {
					String locant = hydrogenAddingEl.getAttributeValue(LOCANT_ATR);
					if (locant !=null && (locant.equals("21") || locant.equals("22") || locant.equals("23") || locant.equals("24"))){
						implicitHydrogenExplicitlySet =true;
					}	
				}
				if (!implicitHydrogenExplicitlySet){
					//porphyrins implicitly have indicated hydrogen at the 21/23 positions
					//directly modify the fragment to avoid problems with locants in for example ring assemblies
					Fragment frag = group.getFrag();
					frag.getAtomByLocantOrThrow("21").setSpareValency(false);
					frag.getAtomByLocantOrThrow("23").setSpareValency(false);
				}
			}
			else if (groupValue.equals("xanthate") || groupValue.equals("xanthic acid") || groupValue.equals("xanthicacid")){
				//This test needs to be here rather in the ComponentGenerator to correctly reject non substituted thioxanthates
				Element wordRule = OpsinTools.getParentWordRule(group);
				if (wordRule.getAttributeValue(WORDRULE_ATR).equals(WordRule.simple.toString())){
					if (OpsinTools.getDescendantElementsWithTagName(wordRule, SUBSTITUENT_EL).size()==0){
						throw new ComponentGenerationException(groupValue +" describes a class of compounds rather than a particular compound");
					}
				}
			}

			if ("yes".equals(group.getAttributeValue(USABLEASJOINER_ATR)) 
					&& group.getAttribute(DEFAULTINID_ATR) == null
					&& group.getAttribute(DEFAULTINLOCANT_ATR) == null) {
				//makes linkers by default attach end to end
				Fragment frag = group.getFrag();
				int chainLength = frag.getChainLength();
				if (chainLength > 1){
					boolean connectEndToEndWithPreviousSub = true;
					if (group.getAttributeValue(TYPE_ATR).equals(CHAIN_TYPE_VAL) && ALKANESTEM_SUBTYPE_VAL.equals(group.getAttributeValue(SUBTYPE_ATR))){//don't do this if the group is preceded by another alkaneStem e.g. methylethyl makes more sense as prop-2-yl rather than propyl
						Element previousSubstituent = OpsinTools.getPreviousSibling(group.getParent());
						if (previousSubstituent != null){
							List<Element> previousSubstGroups = previousSubstituent.getChildElements(GROUP_EL);
							if (previousSubstGroups.size() == 1){
								Element previousGroup = previousSubstGroups.get(0);
								if (previousGroup.getAttributeValue(TYPE_ATR).equals(CHAIN_TYPE_VAL) && ALKANESTEM_SUBTYPE_VAL.equals(previousGroup.getAttributeValue(SUBTYPE_ATR))){
									Element suffixAfterGroup = OpsinTools.getNextSibling(previousGroup, SUFFIX_EL);
									if (suffixAfterGroup == null || suffixAfterGroup.getFrag() == null || suffixAfterGroup.getFrag().getOutAtomCount() == 0){
										connectEndToEndWithPreviousSub = false;
									}
								}
							}
						}
					}
					if (connectEndToEndWithPreviousSub){
						Element parent = group.getParent();
						while (parent.getName().equals(BRACKET_EL)){
							parent = parent.getParent();
						}
						if (parent.getName().equals(ROOT_EL)){
							Element previous = OpsinTools.getPrevious(group);
							if (previous == null || !previous.getName().equals(MULTIPLIER_EL)){
								connectEndToEndWithPreviousSub = false;
							}
						}
					}
					if (connectEndToEndWithPreviousSub){
						group.addAttribute(new Attribute(DEFAULTINID_ATR, Integer.toString(chainLength)));
						frag.setDefaultInAtom(frag.getAtomByLocantOrThrow(Integer.toString(chainLength)));
					}
				}
			}
		}
	}


	/**
	 * Handles Hantzsch-Widman rings. Adds SMILES to the group corresponding to the ring's structure
	 * @param subOrRoot
	 * @throws StructureBuildingException
	 * @throws ComponentGenerationException
	 */
	private void processHW(Element subOrRoot) throws StructureBuildingException, ComponentGenerationException{
		List<Element> hwGroups = OpsinTools.getChildElementsWithTagNameAndAttribute(subOrRoot, GROUP_EL, SUBTYPE_ATR, HANTZSCHWIDMAN_SUBTYPE_VAL);
		for (Element group : hwGroups) {
			Fragment hwRing = group.getFrag();
			List<Atom> atomList =hwRing.getAtomList();
			Element prev = OpsinTools.getPreviousSibling(group);
			List<Element> prevs = new ArrayList<Element>();
			boolean noLocants = true;
			while(prev != null && prev.getName().equals(HETEROATOM_EL)) {
				prevs.add(prev);
				if(prev.getAttribute(LOCANT_ATR) != null) {
					noLocants = false;
				}
				prev = OpsinTools.getPreviousSibling(prev);
			}
			if (atomList.size() == 6 && group.getValue().equals("an")){
				boolean hasNitrogen = false;
				boolean hasSiorGeorSnorPb = false;
				boolean saturatedRing = true;
				for(Element heteroatom : prevs){
					String heteroAtomElement =heteroatom.getAttributeValue(VALUE_ATR);
					Matcher m = MATCH_ELEMENT_SYMBOL.matcher(heteroAtomElement);
					if (!m.find()){
						throw new ComponentGenerationException("Failed to extract element from HW heteroatom");
					}
					heteroAtomElement = m.group();
					if (heteroAtomElement.equals("N")){
						hasNitrogen = true;
					}
					if (heteroAtomElement.equals("Si") ||
						heteroAtomElement.equals("Ge") ||
						heteroAtomElement.equals("Sn") ||
						heteroAtomElement.equals("Pb") ){
						hasSiorGeorSnorPb = true;
					}
				}
				for (Atom a: atomList) {
					if (a.hasSpareValency()){
						saturatedRing = false;
					}
				}
				if (saturatedRing && !hasNitrogen && hasSiorGeorSnorPb){
					throw new ComponentGenerationException("Blocked HW system (6 member saturated ring with no nitrogen but has Si/Ge/Sn/Pb)");
				}
			}
			StringBuilder nameSB = new StringBuilder();
			Collections.reverse(prevs);
			for(Element heteroatom : prevs){
				nameSB.append(heteroatom.getValue());
			}
			nameSB.append(group.getValue());
			String name = nameSB.toString().toLowerCase(Locale.ROOT);
			if(noLocants && prevs.size() > 0) {
				if(specialHWRings.containsKey(name)) {
					String[] specialRingInformation = specialHWRings.get(name);
					String specialInstruction =specialRingInformation[0];
					if (!specialInstruction.equals("")){
						if (specialInstruction.equals("blocked")){
							throw new ComponentGenerationException("Blocked HW system");
						}
						else if (specialInstruction.equals("saturated")){
							for (Atom a: hwRing.getAtomList()) {
								a.setSpareValency(false);
							}
						}
						else if (specialInstruction.equals("not_icacid")){
							if (group.getAttribute(SUBSEQUENTUNSEMANTICTOKEN_ATR) == null){
								Element nextEl = OpsinTools.getNextSibling(group);
								if (nextEl != null && nextEl.getName().equals(SUFFIX_EL) && nextEl.getAttribute(LOCANT_ATR) == null && nextEl.getAttributeValue(VALUE_ATR).equals("ic")){
									throw new ComponentGenerationException(name + nextEl.getValue() +" appears to be a generic class name, not a HW ring");
								}
							}
						}
						else if (specialInstruction.equals("not_nothingOrOlate")){
							if (group.getAttribute(SUBSEQUENTUNSEMANTICTOKEN_ATR) == null){
								Element nextEl = OpsinTools.getNextSibling(group);
								if (nextEl==null || (nextEl!=null && nextEl.getName().equals(SUFFIX_EL) && nextEl.getAttribute(LOCANT_ATR)==null && nextEl.getAttributeValue(VALUE_ATR).equals("ate"))){
									throw new ComponentGenerationException(name +" has the syntax for a HW ring but probably does not mean that in this context");
								}
							}
						}
						else{
							throw new ComponentGenerationException("OPSIN Bug: Unrecognised special HW ring instruction");
						}
					}
					//something like oxazole where by convention locants go 1,3 or a inorganic HW-like system
					for (int j = 1; j < specialRingInformation.length; j++) {
						Atom a = hwRing.getAtomByLocantOrThrow(Integer.toString(j));
						a.setElement(ChemEl.valueOf(specialRingInformation[j]));
					}
					for(Element p : prevs){
						p.detach();
					}
					prevs.clear();
				}
			}
			Set<Element> elementsToRemove = new HashSet<Element>();
			for(Element heteroatom : prevs){//add locanted heteroatoms
				if (heteroatom.getAttribute(LOCANT_ATR) !=null){
					String locant = heteroatom.getAttributeValue(LOCANT_ATR);
					String elementReplacement = heteroatom.getAttributeValue(VALUE_ATR);
					Matcher m = MATCH_ELEMENT_SYMBOL.matcher(elementReplacement);
					if (!m.find()){
						throw new ComponentGenerationException("Failed to extract element from HW heteroatom");
					}
					elementReplacement = m.group();
					Atom a = hwRing.getAtomByLocantOrThrow(locant);
					a.setElement(ChemEl.valueOf(elementReplacement));
					if (heteroatom.getAttribute(LAMBDA_ATR) != null){
						a.setLambdaConventionValency(Integer.parseInt(heteroatom.getAttributeValue(LAMBDA_ATR)));
					}
					heteroatom.detach();
					elementsToRemove.add(heteroatom);
				}
			}
			for(Element p : elementsToRemove){
				prevs.remove(p);
			}

			//add unlocanted heteroatoms
			int defaultLocant = 1;
			for(Element heteroatom : prevs){
				String elementReplacement =heteroatom.getAttributeValue(VALUE_ATR);
				Matcher m = MATCH_ELEMENT_SYMBOL.matcher(elementReplacement);
				if (!m.find()){
					throw new ComponentGenerationException("Failed to extract element from HW heteroatom");
				}
				elementReplacement = m.group();

				while (hwRing.getAtomByLocantOrThrow(Integer.toString(defaultLocant)).getElement() != ChemEl.C){
					defaultLocant++;
				}
				Atom a = hwRing.getAtomByLocantOrThrow(Integer.toString(defaultLocant));
				a.setElement(ChemEl.valueOf(elementReplacement));
				if (heteroatom.getAttribute(LAMBDA_ATR)!=null){
					a.setLambdaConventionValency(Integer.parseInt(heteroatom.getAttributeValue(LAMBDA_ATR)));
				}
				heteroatom.detach();
			}

			List<Element> deltaEls = subOrRoot.getChildElements(DELTA_EL);//add specified double bonds
			for (Element deltaEl : deltaEls) {
				String locantOfDoubleBond = deltaEl.getValue();
				if (locantOfDoubleBond.equals("")){
					Element newUnsaturator = new TokenEl(UNSATURATOR_EL);
					newUnsaturator.addAttribute(new Attribute(VALUE_ATR, "2"));
					OpsinTools.insertAfter(group, newUnsaturator);
				}
				else{
					Atom firstInDoubleBond = hwRing.getAtomByLocantOrThrow(locantOfDoubleBond);
					Atom secondInDoubleBond = hwRing.getAtomByIDOrThrow(firstInDoubleBond.getID() + 1);
					Bond b = firstInDoubleBond.getBondToAtomOrThrow(secondInDoubleBond);
					b.setOrder(2);
				}
				deltaEl.detach();
			}
			group.setValue(name);
		}
	}


	/**
	 * Assigns Element symbols to groups, suffixes and conjunctive suffixes.
	 * Suffixes have preference.
	 * @param subOrRoot
	 * @throws StructureBuildingException 
	 */
	private void assignElementSymbolLocants(Element subOrRoot) throws StructureBuildingException {
		List<Element> groups = subOrRoot.getChildElements(GROUP_EL);
		Element lastGroupElementInSubOrRoot =groups.get(groups.size()-1);
		List<Fragment> suffixFragments = new ArrayList<Fragment>(state.xmlSuffixMap.get(lastGroupElementInSubOrRoot));
		Fragment suffixableFragment = lastGroupElementInSubOrRoot.getFrag();
		//treat conjunctive suffixesas if they were suffixes
		List<Element> conjunctiveGroups = subOrRoot.getChildElements(CONJUNCTIVESUFFIXGROUP_EL);
		for (Element group : conjunctiveGroups) {
			suffixFragments.add(group.getFrag());
		}
		FragmentTools.assignElementLocants(suffixableFragment, suffixFragments);
		for (int i = groups.size()-2; i>=0; i--) {
			FragmentTools.assignElementLocants(groups.get(i).getFrag(), new ArrayList<Fragment>());
		}
	}


	/**
	 * Processes constructs such as biphenyl, 1,1':4',1''-Terphenyl, 2,2'-Bipyridylium, m-Quaterphenyl
	 * @param subOrRoot
	 * @throws ComponentGenerationException
	 * @throws StructureBuildingException
	 */
	private void processRingAssemblies(Element subOrRoot) throws ComponentGenerationException, StructureBuildingException {
		List<Element> ringAssemblyMultipliers = subOrRoot.getChildElements(RINGASSEMBLYMULTIPLIER_EL);
		for (Element multiplier : ringAssemblyMultipliers) {
			int mvalue = Integer.parseInt(multiplier.getAttributeValue(VALUE_ATR));

			/*
			 * Populate locants with locants. Two locants are required for every pair of rings to be joined.
			 * e.g. bi requires 2, ter requires 4 etc.
			 */
			List<List<String>> ringJoiningLocants =new ArrayList<List<String>>();
			Element potentialLocant = OpsinTools.getPreviousSibling(multiplier);
			Element group = OpsinTools.getNextSibling(multiplier, GROUP_EL);
			if (potentialLocant!=null && (potentialLocant.getName().equals(COLONORSEMICOLONDELIMITEDLOCANT_EL)||potentialLocant.getName().equals(LOCANT_EL)) ){//a locant appears to have been provided to indicate how to connect the rings of the ringAssembly
				if (ORTHOMETAPARA_TYPE_VAL.equals(potentialLocant.getAttributeValue(TYPE_ATR))){//an OMP locant has been provided to indicate how to connect the rings of the ringAssembly
					String locant2 =potentialLocant.getValue();
					String locant1 ="1";
					ArrayList<String> locantArrayList =new ArrayList<String>();
					locantArrayList.add("1");
					locantArrayList.add("1'");
					ringJoiningLocants.add(locantArrayList);
					for (int j = 1; j < mvalue -1; j++) {
						locantArrayList =new ArrayList<String>();
						locantArrayList.add(locant2 + StringTools.multiplyString("'", j));
						locantArrayList.add(locant1 + StringTools.multiplyString("'", j+1));
						ringJoiningLocants.add(locantArrayList);
					}
					potentialLocant.detach();
				}
				else{
					String locantText =StringTools.removeDashIfPresent(potentialLocant.getValue());
					//locantText might be something like 1,1':3',1''
					String[] perRingLocantArray = MATCH_COLONORSEMICOLON.split(locantText);
					if (perRingLocantArray.length !=(mvalue -1)){
						throw new ComponentGenerationException("Disagreement between number of locants(" + locantText +") and ring assembly multiplier: " + mvalue);
					}
					if (perRingLocantArray.length!=1 || MATCH_COMMA.split(perRingLocantArray[0]).length!=1){//not for the case of a single locant
						for (String ringLocantArray : perRingLocantArray) {
							String[] locantArray = MATCH_COMMA.split(ringLocantArray);
							if (locantArray.length !=2){
								throw new ComponentGenerationException("missing locant, expected 2 locants: " + ringLocantArray);
							}
							ringJoiningLocants.add(Arrays.asList(locantArray));
						}
						potentialLocant.detach();
					}
				}
			}

			Fragment fragmentToResolveAndDuplicate = group.getFrag();
			Element elementToResolve;//temporary element containing elements that should be resolved before the ring is duplicated
			Element nextEl = OpsinTools.getNextSibling(multiplier);
			if (nextEl.getName().equals(STRUCTURALOPENBRACKET_EL)){//brackets have been provided to aid disambiguation. These brackets are detached e.g. bi(cyclohexyl)
				elementToResolve = new GroupingEl(SUBSTITUENT_EL);
				Element currentEl =nextEl;
				nextEl = OpsinTools.getNextSibling(currentEl);
				currentEl.detach();
				while (nextEl !=null && !nextEl.getName().equals(STRUCTURALCLOSEBRACKET_EL)){
					currentEl =nextEl;
					nextEl = OpsinTools.getNextSibling(currentEl);
					currentEl.detach();
					elementToResolve.addChild(currentEl);
				}
				if (nextEl!=null){
					nextEl.detach();
				}
			}
			else{
				elementToResolve = determineElementsToResolveIntoRingAssembly(multiplier, ringJoiningLocants.size(), fragmentToResolveAndDuplicate.getOutAtomCount());
			}

			List<Element> suffixes = elementToResolve.getChildElements(SUFFIX_EL);
			resolveSuffixes(group, suffixes);
			StructureBuildingMethods.resolveLocantedFeatures(state, elementToResolve);
			StructureBuildingMethods.resolveUnLocantedFeatures(state, elementToResolve);
			group.detach();
			OpsinTools.insertAfter(multiplier, group);

			int bondOrder = 1;
			if (fragmentToResolveAndDuplicate.getOutAtomCount()>0){//e.g. bicyclohexanylidene
				bondOrder =fragmentToResolveAndDuplicate.getOutAtom(0).getValency();
			}
			if (fragmentToResolveAndDuplicate.getOutAtomCount()>1){
				throw new StructureBuildingException("Ring assembly fragment should have one or no OutAtoms; not more than one!");
			}

			List<Fragment> clonedFragments = new ArrayList<Fragment>();
			for (int j = 1; j < mvalue; j++) {
				clonedFragments.add(state.fragManager.copyAndRelabelFragment(fragmentToResolveAndDuplicate, j));
			}
			for (int j = 0; j < mvalue-1; j++) {
				Fragment clone =clonedFragments.get(j);
				Atom atomOnParent;
				Atom atomOnLatestClone;
				if (ringJoiningLocants.size()>0){//locants defined
					atomOnParent = fragmentToResolveAndDuplicate.getAtomByLocantOrThrow(ringJoiningLocants.get(j).get(0));
					String secondLocant = ringJoiningLocants.get(j).get(1);
					if (mvalue ==2 && !secondLocant.endsWith("'")){
						//Allow prime to be (incorrectly) omitted on second locant in bi ring assemblies e.g. 2,2-bipyridine
						try {
							atomOnLatestClone = clone.getAtomByLocantOrThrow(secondLocant);
						}
						catch (StructureBuildingException e){
							atomOnLatestClone = clone.getAtomByLocant(secondLocant + "'");
							if (atomOnLatestClone == null){
								throw e;
							}
						}
					}
					else{
						atomOnLatestClone = clone.getAtomByLocantOrThrow(secondLocant);
					}
				}
				else{
					if (fragmentToResolveAndDuplicate.getOutAtomCount()>0 && mvalue==2){
						atomOnParent = fragmentToResolveAndDuplicate.getOutAtom(0).getAtom();
						atomOnLatestClone = clone.getOutAtom(0).getAtom();
					}
					else{
						atomOnParent =fragmentToResolveAndDuplicate.getAtomOrNextSuitableAtomOrThrow(fragmentToResolveAndDuplicate.getDefaultInAtom(), bondOrder, true);
						atomOnLatestClone = clone.getAtomOrNextSuitableAtomOrThrow(clone.getDefaultInAtom(), bondOrder, true);
					}
				}
				if (fragmentToResolveAndDuplicate.getOutAtomCount()>0){
					fragmentToResolveAndDuplicate.removeOutAtom(0);
				}
				if (clone.getOutAtomCount()>0){
					clone.removeOutAtom(0);
				}
				state.fragManager.incorporateFragment(clone, atomOnLatestClone, fragmentToResolveAndDuplicate, atomOnParent, bondOrder);
				fragmentToResolveAndDuplicate.setDefaultInAtom(clone.getDefaultInAtom());
			}
			group.setValue(multiplier.getValue() +group.getValue());
			Element possibleOpenStructuralBracket = OpsinTools.getPreviousSibling(multiplier);
			if (possibleOpenStructuralBracket!=null && possibleOpenStructuralBracket.getName().equals(STRUCTURALOPENBRACKET_EL)){//e.g. [2,2'-bipyridin].
				//To emphasise there can actually be two sets of structural brackets e.g. [1,1'-bi(cyclohexyl)]
				OpsinTools.getNextSibling(possibleOpenStructuralBracket, STRUCTURALCLOSEBRACKET_EL).detach();
				possibleOpenStructuralBracket.detach();
			}
			multiplier.detach();
		}
	}

	/**
	 * Given the element after the ring assembly multiplier determines which siblings should be resolved by adding them to elementToResolve
	 * @param multiplier
	 * @param ringJoiningLocants
	 * @param outAtomCount
	 * @return
	 * @throws ComponentGenerationException
	 */
	private Element determineElementsToResolveIntoRingAssembly(Element multiplier, int ringJoiningLocants, int outAtomCount) throws ComponentGenerationException {
		Element elementToResolve = new GroupingEl(SUBSTITUENT_EL);
		boolean groupFound = false;
		boolean inlineSuffixSeen = outAtomCount > 0;
		Element currentEl = OpsinTools.getNextSibling(multiplier);
		while (currentEl !=null){
			Element nextEl = OpsinTools.getNextSibling(currentEl);
			if (!groupFound || currentEl.getName().equals(SUFFIX_EL) && currentEl.getAttributeValue(TYPE_ATR).equals(CHARGE_TYPE_VAL)|| currentEl.getName().equals(UNSATURATOR_EL)){
				currentEl.detach();
				elementToResolve.addChild(currentEl);
			}
			else if (currentEl.getName().equals(SUFFIX_EL)){
				if (!inlineSuffixSeen && currentEl.getAttributeValue(TYPE_ATR).equals(INLINE_TYPE_VAL) && currentEl.getAttributeValue(MULTIPLIED_ATR) ==null
						&& (currentEl.getAttribute(LOCANT_ATR)==null || ("2".equals(multiplier.getAttributeValue(VALUE_ATR)) && ringJoiningLocants==0)) && currentEl.getFrag()==null){
					inlineSuffixSeen = true;
					currentEl.detach();
					elementToResolve.addChild(currentEl);
				}
				else{
					break;
				}
			}
			else{
				break;
			}
			if (currentEl.getName().equals(GROUP_EL)){
				groupFound = true;
			}
			currentEl = nextEl;
		}
		Element parent = multiplier.getParent();
		if (!parent.getName().equals(SUBSTITUENT_EL) && OpsinTools.getChildElementsWithTagNameAndAttribute(parent, SUFFIX_EL, TYPE_ATR, INLINE_TYPE_VAL).size()!=0){
			throw new ComponentGenerationException("Unexpected radical adding suffix on ring assembly");			
		}
		return elementToResolve;
	}


	/**
	 * Proccess any polycyclic spiro systems present in subOrRoot
	 * It is assumed that at this stage all hantzch widman rings/fused rings have been resolved to single groups allowing them to be simply spiro fused
	 * 
	 * http://www.chem.qmul.ac.uk/iupac/spiro/ (SP-2 through SP-6)
	 * @param subOrRoot
	 * @throws ComponentGenerationException
	 * @throws StructureBuildingException
	 */
	private void processPolyCyclicSpiroNomenclature(Element subOrRoot) throws ComponentGenerationException, StructureBuildingException {
		List<Element> polyCyclicSpiros = subOrRoot.getChildElements(POLYCYCLICSPIRO_EL);
		if (polyCyclicSpiros.size()>0){
			Element polyCyclicSpiroDescriptor = polyCyclicSpiros.get(0);
			String value = polyCyclicSpiroDescriptor.getAttributeValue(VALUE_ATR);
			if (value.equals("spiro")){
				if (polyCyclicSpiros.size()!=1){
					throw new ComponentGenerationException("Nested polyspiro systems are not supported");
				}
				processNonIdenticalPolyCyclicSpiro(polyCyclicSpiroDescriptor);
			}
			else if (value.equals("spiroOldMethod")){
				processOldMethodPolyCyclicSpiro(polyCyclicSpiros);
			}
			else if (value.equals("spirobi")){
				if (polyCyclicSpiros.size()!=1){
					throw new ComponentGenerationException("Nested polyspiro systems are not supported");
				}
				processSpiroBiOrTer(polyCyclicSpiroDescriptor, 2);
			}
			else if (value.equals("spiroter")){
				if (polyCyclicSpiros.size()!=1){
					throw new ComponentGenerationException("Nested polyspiro systems are not supported");
				}
				processSpiroBiOrTer(polyCyclicSpiroDescriptor, 3);
			}
			else if (value.equals("dispiroter")){
				if (polyCyclicSpiros.size()!=1){
					throw new ComponentGenerationException("Nested polyspiro systems are not supported");
				}
				processDispiroter(polyCyclicSpiroDescriptor);
			}
			else{
				throw new ComponentGenerationException("Unsupported spiro system encountered");
			}
			polyCyclicSpiroDescriptor.detach();
		}
	}

	private void processNonIdenticalPolyCyclicSpiro(Element polyCyclicSpiroDescriptor) throws ComponentGenerationException, StructureBuildingException {
		Element subOrRoot = polyCyclicSpiroDescriptor.getParent();
		List<Element> groups = subOrRoot.getChildElements(GROUP_EL);
		if (groups.size()<2){
			throw new ComponentGenerationException("OPSIN Bug: Atleast two groups were expected in polycyclic spiro system");
		}
		Element openBracket = OpsinTools.getNextSibling(polyCyclicSpiroDescriptor);
		if (!openBracket.getName().equals(STRUCTURALOPENBRACKET_EL)){
			throw new ComponentGenerationException("OPSIN Bug: Open bracket not found where open bracket expeced");
		}
		List<Element> spiroBracketElements = OpsinTools.getSiblingsUpToElementWithTagName(openBracket, STRUCTURALCLOSEBRACKET_EL);
		Element closeBracket = OpsinTools.getNextSibling(spiroBracketElements.get(spiroBracketElements.size()-1));
		if (closeBracket == null || !closeBracket.getName().equals(STRUCTURALCLOSEBRACKET_EL)){
			throw new ComponentGenerationException("OPSIN Bug: Open bracket not found where open bracket expeced");
		}
		
		Element firstGroup = groups.get(0);
		List<Element> firstGroupEls = new ArrayList<Element>();
		int indexOfOpenBracket = subOrRoot.indexOf(openBracket);
		int indexOfFirstGroup = subOrRoot.indexOf(firstGroup);
		for (int i =indexOfOpenBracket +1; i < indexOfFirstGroup; i++) {
			firstGroupEls.add(subOrRoot.getChild(i));
		}
		firstGroupEls.add(firstGroup);
		firstGroupEls.addAll(OpsinTools.getNextAdjacentSiblingsOfType(firstGroup, UNSATURATOR_EL));
		resolveFeaturesOntoGroup(firstGroupEls);
		Set<Atom> spiroAtoms = new HashSet<Atom>();
		for (int i = 1; i < groups.size(); i++) {
			Element nextGroup =groups.get(i);
			Element locant = OpsinTools.getNextSibling(groups.get(i-1), SPIROLOCANT_EL);
			if (locant ==null){
				throw new ComponentGenerationException("Unable to find locantEl for polycyclic spiro system");
			}
			
			List<Element> nextGroupEls = new ArrayList<Element>();
			int indexOfLocant = subOrRoot.indexOf(locant);
			int indexOfNextGroup = subOrRoot.indexOf(nextGroup);
			for (int j =indexOfLocant +1; j < indexOfNextGroup; j++) {
				nextGroupEls.add(subOrRoot.getChild(j));
			}
			nextGroupEls.add(nextGroup);
			nextGroupEls.addAll(OpsinTools.getNextAdjacentSiblingsOfType(nextGroup, UNSATURATOR_EL));
			resolveFeaturesOntoGroup(nextGroupEls);
			
			String[] locants = MATCH_COMMA.split(StringTools.removeDashIfPresent(locant.getValue()));
			if (locants.length!=2){
				throw new ComponentGenerationException("Incorrect number of locants found before component of polycyclic spiro system");
			}
			for (int j = 0; j < locants.length; j++) {
				String locantText= locants[j];
				Matcher m = matchAddedHydrogenBracket.matcher(locantText);
				if (m.find()){
					Element addedHydrogenElement=new TokenEl(ADDEDHYDROGEN_EL);
					addedHydrogenElement.addAttribute(new Attribute(LOCANT_ATR, m.group(1)));
					OpsinTools.insertBefore(locant, addedHydrogenElement);
					locant.addAttribute(new Attribute(TYPE_ATR, ADDEDHYDROGENLOCANT_TYPE_VAL));
					locants[j] = m.replaceAll("");
				}
			}
			locant.detach();
			Fragment nextFragment = nextGroup.getFrag();
			FragmentTools.relabelNumericLocants(nextFragment.getAtomList(), StringTools.multiplyString("'", i));
			String secondLocant = locants[1];
			Atom atomToBeReplaced;
			if (secondLocant.endsWith("'")){
				atomToBeReplaced = nextFragment.getAtomByLocantOrThrow(locants[1]);
			}
			else{
				//for simple spiro fusions the prime is often forgotten
				atomToBeReplaced = nextFragment.getAtomByLocantOrThrow(locants[1] + "'");
			}
			Atom atomOnParentFrag = null;
			for (int j = 0; j < i; j++) {
				atomOnParentFrag = groups.get(j).getFrag().getAtomByLocant(locants[0]);
				if (atomOnParentFrag!=null){
					break;
				}
			}
			if (atomOnParentFrag==null){
				throw new ComponentGenerationException("Could not find the atom with locant " + locants[0] +" for use in polycyclic spiro system");
			}
			spiroAtoms.add(atomOnParentFrag);
			state.fragManager.replaceAtomWithAnotherAtomPreservingConnectivity(atomToBeReplaced, atomOnParentFrag);
			if (atomToBeReplaced.hasSpareValency()){
				atomOnParentFrag.setSpareValency(true);
			}
		}
		if (spiroAtoms.size()>1){
			Element expectedMultiplier = OpsinTools.getPreviousSibling(polyCyclicSpiroDescriptor);
			if (expectedMultiplier!=null && expectedMultiplier.getName().equals(MULTIPLIER_EL) && Integer.parseInt(expectedMultiplier.getAttributeValue(VALUE_ATR))==spiroAtoms.size()){
				expectedMultiplier.detach();
			}
		}
		Element rootGroup = groups.get(groups.size()-1);
		Fragment rootFrag = rootGroup.getFrag();
		String name = rootGroup.getValue();
		for (int i = 0; i < groups.size() -1; i++) {
			Element group =groups.get(i);
			state.fragManager.incorporateFragment(group.getFrag(), rootFrag);
			name = group.getValue() + name;
			group.detach();
		}
		rootGroup.setValue(polyCyclicSpiroDescriptor.getValue() + name);
		openBracket.detach();
		closeBracket.detach();
	}
	

	/**
	 * Processes spiro systems described using the now deprectated method described in the 1979 guidelines Rule A-42
	 * @param spiroElements
	 * @throws ComponentGenerationException
	 * @throws StructureBuildingException
	 */
	private void processOldMethodPolyCyclicSpiro(List<Element> spiroElements) throws ComponentGenerationException, StructureBuildingException {
		Element firstSpiro =spiroElements.get(0);
		Element subOrRoot = firstSpiro.getParent();
		Element firstEl = subOrRoot.getChild(0);
		List<Element> elementsToResolve = OpsinTools.getSiblingsUpToElementWithTagName(firstEl, POLYCYCLICSPIRO_EL);
		elementsToResolve.add(0, firstEl);
		resolveFeaturesOntoGroup(elementsToResolve);
		
		for (int i = 0; i < spiroElements.size(); i++) {
			Element currentSpiro = spiroElements.get(i);
			Element previousGroup = OpsinTools.getPreviousSibling(currentSpiro, GROUP_EL);
			if (previousGroup==null){
				throw new ComponentGenerationException("OPSIN bug: unable to locate group before polycylic spiro descriptor");
			}
			Element nextGroup = OpsinTools.getNextSibling(currentSpiro, GROUP_EL);
			if (nextGroup==null){
				throw new ComponentGenerationException("OPSIN bug: unable to locate group after polycylic spiro descriptor");
			}
			Fragment parentFrag = nextGroup.getFrag();
			Fragment previousFrag = previousGroup.getFrag();
			FragmentTools.relabelNumericLocants(parentFrag.getAtomList(), StringTools.multiplyString("'",i+1));
			elementsToResolve = OpsinTools.getSiblingsUpToElementWithTagName(currentSpiro, POLYCYCLICSPIRO_EL);
			resolveFeaturesOntoGroup(elementsToResolve);
			
			String locant1 =null;
			Element possibleFirstLocant = OpsinTools.getPreviousSibling(currentSpiro);
			if (possibleFirstLocant !=null && possibleFirstLocant.getName().equals(LOCANT_EL)){
				if (MATCH_COMMA.split(possibleFirstLocant.getValue()).length==1){
					locant1 = possibleFirstLocant.getValue();
					possibleFirstLocant.detach();
				}
				else{
					throw new ComponentGenerationException("Malformed locant before polycyclic spiro descriptor");
				}
			}
			Atom atomToBeReplaced;
			if (locant1 != null){
				atomToBeReplaced = previousFrag.getAtomByLocantOrThrow(locant1);
			}
			else{
				atomToBeReplaced = previousFrag.getAtomOrNextSuitableAtomOrThrow(previousFrag.getFirstAtom(), 2, true);
			}
			Atom atomOnParentFrag;
			String locant2 =null;
			Element possibleSecondLocant = OpsinTools.getNextSibling(currentSpiro);
			if (possibleSecondLocant !=null && possibleSecondLocant.getName().equals(LOCANT_EL)){
				if (MATCH_COMMA.split(possibleSecondLocant.getValue()).length==1){
					locant2 = possibleSecondLocant.getValue();
					possibleSecondLocant.detach();
				}
				else{
					throw new ComponentGenerationException("Malformed locant after polycyclic spiro descriptor");
				}
			}
			if (locant2!=null){
				atomOnParentFrag = parentFrag.getAtomByLocantOrThrow(locant2);
			}
			else{
				atomOnParentFrag = parentFrag.getAtomOrNextSuitableAtomOrThrow(parentFrag.getFirstAtom(), 2, true);
			}
			state.fragManager.replaceAtomWithAnotherAtomPreservingConnectivity(atomToBeReplaced, atomOnParentFrag);
			if (atomToBeReplaced.hasSpareValency()){
				atomOnParentFrag.setSpareValency(true);
			}
			if (atomToBeReplaced.getCharge()!=0 && atomOnParentFrag.getCharge()==0){
				atomOnParentFrag.setCharge(atomToBeReplaced.getCharge());
				atomOnParentFrag.setProtonsExplicitlyAddedOrRemoved(atomToBeReplaced.getProtonsExplicitlyAddedOrRemoved());
			}
			state.fragManager.incorporateFragment(previousFrag, parentFrag);
			nextGroup.setValue(previousGroup.getValue() + currentSpiro.getValue() + nextGroup.getValue());
			previousGroup.detach();
		}
	}


	/**
	 * Two or three copies of the fragment after polyCyclicSpiroDescriptor are spiro fused at one centre
	 * @param polyCyclicSpiroDescriptor
	 * @param components
	 * @throws ComponentGenerationException
	 * @throws StructureBuildingException
	 */
	private void processSpiroBiOrTer(Element polyCyclicSpiroDescriptor, int components) throws ComponentGenerationException, StructureBuildingException {
		Element locant = OpsinTools.getPreviousSibling(polyCyclicSpiroDescriptor);
		String locantText;
		if (locant ==null || !locant.getName().equals(LOCANT_EL)){
			if (components==2){
				locantText ="1,1'";
			}
			else{
				throw new ComponentGenerationException("Unable to find locant indicating atoms to form polycyclic spiro system!");
			}
		}
		else{
			locantText = locant.getValue();
			locant.detach();
		}
		String[] locants = MATCH_COMMA.split(locantText);
		if (locants.length!=components){
			throw new ComponentGenerationException("Mismatch between spiro descriptor and number of locants provided");
		}
		Element group = OpsinTools.getNextSibling(polyCyclicSpiroDescriptor, GROUP_EL);
		if (group==null){
			throw new ComponentGenerationException("Cannot find group to which spirobi/ter descriptor applies");
		}

		determineFeaturesToResolveInSingleComponentSpiro(polyCyclicSpiroDescriptor);
		Fragment fragment = group.getFrag();
		List<Fragment> clones = new ArrayList<Fragment>();
		for (int i = 1; i < components ; i++) {
			clones.add(state.fragManager.copyAndRelabelFragment(fragment, i));
		}
		
		Atom atomOnOriginalFragment = fragment.getAtomByLocantOrThrow(locants[0]);
		for (int i = 1; i < components ; i++) {
			Fragment clone = clones.get(i - 1);
			Atom atomToBeReplaced;
			if (components ==2 && !locants[i].endsWith("'")){
				//Allow prime to be (incorrectly) omitted on second locant in spirobi
				try {
					atomToBeReplaced = clone.getAtomByLocantOrThrow(locants[i]);
				}
				catch (StructureBuildingException e){
					atomToBeReplaced = clone.getAtomByLocant(locants[i] + "'");
					if (atomToBeReplaced == null){
						throw e;
					}
				}
			}
			else{
				atomToBeReplaced = clone.getAtomByLocantOrThrow(locants[i]);
			}

			state.fragManager.replaceAtomWithAnotherAtomPreservingConnectivity(atomToBeReplaced, atomOnOriginalFragment);
			if (atomToBeReplaced.hasSpareValency()){
				atomOnOriginalFragment.setSpareValency(true);
			}
		}
		for (Fragment clone : clones) {
			state.fragManager.incorporateFragment(clone, fragment);
		}
		
		group.setValue(polyCyclicSpiroDescriptor.getValue() + group.getValue());
	}

	/**
	 * Three copies of the fragment after polyCyclicSpiroDescriptor are spiro fused at two centres
	 * @param polyCyclicSpiroDescriptor
	 * @throws StructureBuildingException
	 * @throws ComponentGenerationException
	 */
	private void processDispiroter(Element polyCyclicSpiroDescriptor) throws StructureBuildingException, ComponentGenerationException {
		String value = polyCyclicSpiroDescriptor.getValue();
		value = value.substring(0, value.length()-10);//remove dispiroter
		value = StringTools.removeDashIfPresent(value);
		String[] locants = MATCH_COLON.split(value);
		Element group = OpsinTools.getNextSibling(polyCyclicSpiroDescriptor, GROUP_EL);
		if (group==null){
			throw new ComponentGenerationException("Cannot find group to which dispiroter descriptor applies");
		}
		determineFeaturesToResolveInSingleComponentSpiro(polyCyclicSpiroDescriptor);
		Fragment fragment = group.getFrag();
		List<Fragment> clones = new ArrayList<Fragment>();
		for (int i = 1; i < 3 ; i++) {
			clones.add(state.fragManager.copyAndRelabelFragment(fragment, i));
		}
		for (Fragment clone : clones) {
			state.fragManager.incorporateFragment(clone, fragment);
		}
		
		Atom atomOnLessPrimedFragment = fragment.getAtomByLocantOrThrow(MATCH_COMMA.split(locants[0])[0]);
		Atom atomToBeReplaced = fragment.getAtomByLocantOrThrow(MATCH_COMMA.split(locants[0])[1]);
		state.fragManager.replaceAtomWithAnotherAtomPreservingConnectivity(atomToBeReplaced, atomOnLessPrimedFragment);
		if (atomToBeReplaced.hasSpareValency()){
			atomOnLessPrimedFragment.setSpareValency(true);
		}
		
		atomOnLessPrimedFragment = fragment.getAtomByLocantOrThrow(MATCH_COMMA.split(locants[1])[0]);
		atomToBeReplaced = fragment.getAtomByLocantOrThrow(MATCH_COMMA.split(locants[1])[1]);
		state.fragManager.replaceAtomWithAnotherAtomPreservingConnectivity(atomToBeReplaced, atomOnLessPrimedFragment);
		if (atomToBeReplaced.hasSpareValency()){
			atomOnLessPrimedFragment.setSpareValency(true);
		}

		group.setValue("dispiroter" + group.getValue());
	}

	/**
	 * The features between the polyCyclicSpiroDescriptor and the first group element, or beween the STRUCTURALOPENBRACKET_EL and STRUCTURALCLOSEBRACKET_EL
	 * are found and then passed to resolveFeaturesOntoGroup
	 * @param polyCyclicSpiroDescriptor
	 * @throws StructureBuildingException
	 * @throws ComponentGenerationException 
	 */
	private void determineFeaturesToResolveInSingleComponentSpiro(Element polyCyclicSpiroDescriptor) throws StructureBuildingException, ComponentGenerationException {
		Element possibleOpenBracket = OpsinTools.getNextSibling(polyCyclicSpiroDescriptor);
		List<Element> elementsToResolve;
		if (possibleOpenBracket.getName().equals(STRUCTURALOPENBRACKET_EL)){
			possibleOpenBracket.detach();
			elementsToResolve = OpsinTools.getSiblingsUpToElementWithTagName(polyCyclicSpiroDescriptor, STRUCTURALCLOSEBRACKET_EL);
			OpsinTools.getNextSibling(elementsToResolve.get(elementsToResolve.size()-1)).detach();//detach close bracket
		}
		else{
			elementsToResolve = OpsinTools.getSiblingsUpToElementWithTagName(polyCyclicSpiroDescriptor, GROUP_EL);
		}
		resolveFeaturesOntoGroup(elementsToResolve);
	}
	
	/**
	 * Given some elements including a group element resolves all locanted and unlocanted features.
	 * If suffixes are present these are resolved and detached
	 * @param elementsToResolve
	 * @throws StructureBuildingException 
	 * @throws ComponentGenerationException 
	 */
	private void resolveFeaturesOntoGroup(List<Element> elementsToResolve) throws StructureBuildingException, ComponentGenerationException{
		if (elementsToResolve.size()==0){
			return;
		}
		Element substituentToResolve = new GroupingEl(SUBSTITUENT_EL);//temporary element containing elements that should be resolved before the ring is cloned
		Element parent = elementsToResolve.get(0).getParent();
		int index = parent.indexOf(elementsToResolve.get(0));
		Element group =null;
		List<Element> suffixes = new ArrayList<Element>();
		Element locant =null;
		for (Element element : elementsToResolve) {
			String elName =element.getName();
			if (elName.equals(GROUP_EL)){
				group = element;
			}
			else if (elName.equals(SUFFIX_EL)){
				suffixes.add(element);
			}
			else if (elName.equals(LOCANT_EL) && group==null){
				locant = element;
			}
			element.detach();
			substituentToResolve.addChild(element);
		}
		if (group ==null){
			throw new ComponentGenerationException("OPSIN bug: group element should of been given to method");
		}
		if (locant !=null){//locant is probably an indirect locant, try and assign it
			List<Element> locantAble = findElementsMissingIndirectLocants(substituentToResolve, locant);
			String[] locantValues = MATCH_COMMA.split(locant.getValue());
			if (locantAble.size() >= locantValues.length){
				for (int i = 0; i < locantValues.length; i++) {
					String locantValue = locantValues[i];
					locantAble.get(i).addAttribute(new Attribute(LOCANT_ATR, locantValue));
				}
				locant.detach();
			}
		}
		if (!suffixes.isEmpty()){
			resolveSuffixes(group, suffixes);
			for (Element suffix : suffixes) {
				suffix.detach();
			}
		}
		if (substituentToResolve.getChildCount() != 0){
			StructureBuildingMethods.resolveLocantedFeatures(state, substituentToResolve);
			StructureBuildingMethods.resolveUnLocantedFeatures(state, substituentToResolve);
			List<Element> children = substituentToResolve.getChildElements();
			for (int i = children.size() -1; i>=0; i--) {
				Element child = children.get(i);
				child.detach();
				parent.insertChild(child, index);
			}
		}
	}
	
	private static class SortBridgesByHighestLocantedBridgehead implements Comparator<Fragment>{

		private final Map<Fragment, Atom[]> bridgeToRingAtoms;
		SortBridgesByHighestLocantedBridgehead(Map<Fragment, Atom[]> bridgeToRingAtoms) {
			this.bridgeToRingAtoms = bridgeToRingAtoms;
		}

		public int compare(Fragment bridge1, Fragment bridge2) {
			Atom[] ringAtoms1 = bridgeToRingAtoms.get(bridge1);
			int bridge1HighestRingLocant = Math.max(getLocantNumber(ringAtoms1[0]),getLocantNumber(ringAtoms1[1]));
			
			Atom[] ringAtoms2 = bridgeToRingAtoms.get(bridge2);
			int bridge2HighestRingLocant = Math.max(getLocantNumber(ringAtoms2[0]),getLocantNumber(ringAtoms2[1]));
			
			if (bridge1HighestRingLocant > bridge2HighestRingLocant){
				return -1;
			}
			if (bridge1HighestRingLocant < bridge2HighestRingLocant){
				return 1;
			}
			return 0;
		}
		
	}

	/**
	 * Processes bridges e.g. 4,7-methanoindene
	 * Resolves and attaches said bridges to the adjacent ring fragment
	 * Numbers the bridges in accordance with FR-8.6/FR-8.7
	 * @param subOrRoot
	 * @throws StructureBuildingException 
	 */
	private void processFusedRingBridges(Element subOrRoot) throws StructureBuildingException {
		List<Element> bridges = subOrRoot.getChildElements(FUSEDRINGBRIDGE_EL);
		int bridgeCount = bridges.size();
		if (bridgeCount == 0) {
			return;
		}
		Element groupEl = OpsinTools.getNextSibling(bridges.get(bridgeCount - 1), GROUP_EL);
		Fragment ringFrag = groupEl.getFrag();
		Map<Fragment, Atom[]> bridgeToRingAtoms = new LinkedHashMap<Fragment, Atom[]>();
		for (Element bridge : bridges) {
			Fragment bridgeFrag = state.fragManager.buildSMILES(bridge.getAttributeValue(VALUE_ATR), groupEl, NONE_LABELS_VAL);
			List<Atom> bridgeAtomList = bridgeFrag.getAtomList();
			bridgeFrag.addOutAtom(bridgeAtomList.get(0), 1, true);
			bridgeFrag.addOutAtom(bridgeAtomList.get(bridgeAtomList.size() - 1), 1, true);
			Element possibleLocant = OpsinTools.getPreviousSibling(bridge);
			Atom[] ringAtoms;
			if (possibleLocant != null && possibleLocant.getName().equals(LOCANT_EL)){
				String[] locantArray = MATCH_COMMA.split(possibleLocant.getValue());
				if (locantArray.length == 2){
					bridgeFrag.getOutAtom(0).setLocant(locantArray[0]);
					bridgeFrag.getOutAtom(1).setLocant(locantArray[1]);
					possibleLocant.detach();
				}
				ringAtoms = StructureBuildingMethods.formEpoxide(state, bridgeFrag, ringFrag.getDefaultInAtom());
			}
			else{
				ringAtoms = StructureBuildingMethods.formEpoxide(state, bridgeFrag, ringFrag.getAtomOrNextSuitableAtomOrThrow(ringFrag.getDefaultInAtom(), 1, true));
			}
			bridgeToRingAtoms.put(bridgeFrag, ringAtoms);
			state.fragManager.incorporateFragment(bridgeFrag, ringFrag);
			Element unsaturator = OpsinTools.getNextSibling(bridge);
			if (unsaturator.getName().equals(UNSATURATOR_EL)) {
				unsaturator.detach();
			}
			bridge.detach();
		}
		int highestLocant = getHighestNumericLocant(ringFrag);
		List<Fragment> bridgeFragments = new ArrayList<Fragment>(bridgeToRingAtoms.keySet());
		Collections.sort(bridgeFragments, new SortBridgesByHighestLocantedBridgehead(bridgeToRingAtoms));
		for (Fragment bridgeFragment: bridgeFragments) {
			List<Atom> bridgeFragmentAtoms = bridgeFragment.getAtomList();
			Atom[] ringAtoms = bridgeToRingAtoms.get(bridgeFragment);
			if (getLocantNumber(ringAtoms[0]) <= getLocantNumber(ringAtoms[1])){
				for (int i = bridgeFragmentAtoms.size() - 1; i >=0; i--) {
					bridgeFragmentAtoms.get(i).addLocant(String.valueOf(++highestLocant));
				}
			}
			else{
				for (Atom atom : bridgeFragmentAtoms) {
					atom.addLocant(String.valueOf(++highestLocant));
				}
			}
		}
	}

	private static int getLocantNumber(Atom atom) {
		String locant = atom.getFirstLocant();
		if (locant != null) {
			Matcher m = MATCH_NUMERIC_LOCANT.matcher(locant);
			if (m.matches()){
				return Integer.parseInt(m.group(1));
			}
		}
		return 0;
	}

	private int getHighestNumericLocant(Fragment ringFrag) {
		for (int i = 1; ; i++) {
			if (ringFrag.getAtomByLocant(String.valueOf(i)) == null){
				return i - 1;
			}
		}
	}

	/**
	 * Searches for lambdaConvention elements and applies the valency they specify to the atom
	 * they specify on the substituent/root's fragment
	 * @param subOrRoot
	 * @throws StructureBuildingException
	 */
	private void applyLambdaConvention(Element subOrRoot) throws StructureBuildingException {
		List<Element> lambdaConventionEls = subOrRoot.getChildElements(LAMBDACONVENTION_EL);
		for (Element lambdaConventionEl : lambdaConventionEls) {
			Fragment frag = subOrRoot.getFirstChildElement(GROUP_EL).getFrag();
			if (lambdaConventionEl.getAttribute(LOCANT_ATR)!=null){
				frag.getAtomByLocantOrThrow(lambdaConventionEl.getAttributeValue(LOCANT_ATR)).setLambdaConventionValency(Integer.parseInt(lambdaConventionEl.getAttributeValue(LAMBDA_ATR)));
			}
			else{
				if (frag.getAtomCount()!=1){
					throw new StructureBuildingException("Ambiguous use of lambda convention. Fragment has more than 1 atom but no locant was specified for the lambda");
				}
				frag.getFirstAtom().setLambdaConventionValency(Integer.parseInt(lambdaConventionEl.getAttributeValue(LAMBDA_ATR)));
			}
			lambdaConventionEl.detach();
		}
	}


	/**
	 * Uses the number of outAtoms that are present to assign the number of outAtoms on substituents that can have a variable number of outAtoms
	 * Hence at this point it can be determined if a multi radical susbtituent is present in the name
	 * This would be expected in multiplicative nomenclature and is noted in the state so that the StructureBuilder knows to resolve the
	 * section of the name from that point onwards in a left to right manner rather than right to left
	 * @param subOrRoot: The sub/root to look in
	 * @throws ComponentGenerationException
	 * @throws StructureBuildingException
	 */
	private void handleMultiRadicals(Element subOrRoot) throws ComponentGenerationException, StructureBuildingException{
		Element group =subOrRoot.getFirstChildElement(GROUP_EL);
		String groupValue =group.getValue();
		Fragment thisFrag = group.getFrag();
		if (groupValue.equals("methylene") || groupValue.equals("oxy")|| matchChalcogenReplacement.matcher(groupValue).matches()){//resolves for example trimethylene to propan-1,3-diyl or dithio to disulfan-1,2-diyl. Locants may not be specified before the multiplier
			Element beforeGroup = OpsinTools.getPreviousSibling(group);
			if (beforeGroup!=null && beforeGroup.getName().equals(MULTIPLIER_ATR) && beforeGroup.getAttributeValue(TYPE_ATR).equals(BASIC_TYPE_VAL) && OpsinTools.getPreviousSibling(beforeGroup)==null){
				int multiplierVal = Integer.parseInt(beforeGroup.getAttributeValue(VALUE_ATR));
				if (!unsuitableForFormingChainMultiradical(group, beforeGroup)){
					if (groupValue.equals("methylene")){
						group.getAttribute(VALUE_ATR).setValue(StringTools.multiplyString("C", multiplierVal));
					}
					else if (groupValue.equals("oxy")){
						group.getAttribute(VALUE_ATR).setValue(StringTools.multiplyString("O", multiplierVal));
					}
					else if (groupValue.equals("thio")){
						group.getAttribute(VALUE_ATR).setValue(StringTools.multiplyString("S", multiplierVal));
					}
					else if (groupValue.equals("seleno")){
						group.getAttribute(VALUE_ATR).setValue(StringTools.multiplyString("[SeH?]", multiplierVal));
					}
					else if (groupValue.equals("telluro")){
						group.getAttribute(VALUE_ATR).setValue(StringTools.multiplyString("[TeH?]", multiplierVal));
					}
					else{
						throw new ComponentGenerationException("unexpected group value");
					}
					group.getAttribute(OUTIDS_ATR).setValue("1,"+Integer.parseInt(beforeGroup.getAttributeValue(VALUE_ATR)));
					group.setValue(beforeGroup.getValue() + groupValue);
					beforeGroup.detach();
					if (group.getAttribute(LABELS_ATR)!=null){//use numeric numbering
						group.getAttribute(LABELS_ATR).setValue(NUMERIC_LABELS_VAL);
					}
					else{
						group.addAttribute(new Attribute(LABELS_ATR, NUMERIC_LABELS_VAL));
					}
					state.fragManager.removeFragment(thisFrag);
					thisFrag = resolveGroup(state, group);
					group.removeAttribute(group.getAttribute(USABLEASJOINER_ATR));
				}
			}
		}

		if (group.getAttribute(OUTIDS_ATR)!=null){//adds outIDs at the specified atoms
			String[] radicalPositions = MATCH_COMMA.split(group.getAttributeValue(OUTIDS_ATR));
			int firstIdInFrag =thisFrag.getIdOfFirstAtom();
			for (String radicalID : radicalPositions) {
				thisFrag.addOutAtom(firstIdInFrag + Integer.parseInt(radicalID) - 1, 1, true);
			}
		}
		int outAtomCount = thisFrag.getOutAtomCount();
		if (outAtomCount >=2){
			if (groupValue.equals("amine")){//amine is a special case as it shouldn't technically be allowed but is allowed due to it's common usage in EDTA
				Element previousGroup = OpsinTools.getPreviousGroup(group);
				Element nextGroup = OpsinTools.getNextGroup(group);
				if (previousGroup==null || previousGroup.getFrag().getOutAtomCount() < 2 || nextGroup==null){//must be preceded by a multi radical
					throw new ComponentGenerationException("Invalid use of amine as a substituent!");
				}
			}
			if (state.currentWordRule == WordRule.polymer){
				if (outAtomCount >=3){//In poly mode nothing may have more than 2 outAtoms e.g. nitrilo is -N= or =N-
					int valency =0;
					for (int i = 2; i < outAtomCount; i++) {
						OutAtom nextOutAtom = thisFrag.getOutAtom(i);
						valency += nextOutAtom.getValency();
						thisFrag.removeOutAtom(nextOutAtom);
					}
					thisFrag.getOutAtom(1).setValency(thisFrag.getOutAtom(1).getValency() + valency);
				}
			}
		}
		
		if (outAtomCount ==2 && EPOXYLIKE_SUBTYPE_VAL.equals(group.getAttributeValue(SUBTYPE_ATR))){
			Element possibleLocant = OpsinTools.getPreviousSibling(group);
			if (possibleLocant !=null){
				String[] locantValues = MATCH_COMMA.split(possibleLocant.getValue());
				if (locantValues.length==2){
					thisFrag.getOutAtom(0).setLocant(locantValues[0]);
					thisFrag.getOutAtom(1).setLocant(locantValues[1]);
					possibleLocant.detach();
					subOrRoot.addAttribute(new Attribute(LOCANT_ATR, locantValues[0]));
				}
			}
		}

		int totalOutAtoms = outAtomCount + calculateOutAtomsToBeAddedFromInlineSuffixes(group, subOrRoot.getChildElements(SUFFIX_EL));
		if (totalOutAtoms >= 2){
			group.addAttribute(new Attribute (ISAMULTIRADICAL_ATR, Integer.toString(totalOutAtoms)));
		}
	}

	/**
	 * Checks for cases where multiplier(methylene) or multiplier(thio) and the like should not be interpreted as one fragment
	 * Something like nitrilotrithiotriacetic acid or oxetane-3,3-diyldimethylene
	 * @param group
	 * @param multiplierBeforeGroup 
	 * @return
	 */
	private boolean unsuitableForFormingChainMultiradical(Element group, Element multiplierBeforeGroup) {
		Element previousGroup = OpsinTools.getPreviousGroup(group);
		if (previousGroup!=null){
			if (previousGroup.getAttribute(ISAMULTIRADICAL_ATR)!=null){
				if (previousGroup.getAttributeValue(ACCEPTSADDITIVEBONDS_ATR)!=null && OpsinTools.getPreviousSibling(previousGroup.getParent())!=null){
					return false;
				}
				//the initial multiplier is proceded by another multiplier e.g. bis(dithio)
				if (OpsinTools.getPrevious(multiplierBeforeGroup).getName().equals(MULTIPLIER_EL)){
					return false;
				}
				if (previousGroup.getAttributeValue(ISAMULTIRADICAL_ATR).equals(multiplierBeforeGroup.getAttributeValue(VALUE_ATR))){
					return true;//probably multiplicative
				}
				else{
					return false;
				}
			}
			else if (OpsinTools.getPreviousSibling(previousGroup, MULTIPLIER_EL)==null){
				//This is a 99% solution to the detection of cases such as ethylidenedioxy == ethan-1,1-diyldioxy
				Fragment previousGroupFrag = previousGroup.getFrag();
				int outAtomValency =0;
				if (previousGroupFrag.getOutAtomCount()==1){
					outAtomValency = previousGroupFrag.getOutAtom(0).getValency();
				}
				else{
					Element suffix = OpsinTools.getNextSibling(previousGroup, SUFFIX_EL);
					if (suffix!=null && suffix.getAttributeValue(VALUE_ATR).equals("ylidene")){
						outAtomValency =2;
					}
					if (suffix!=null && suffix.getAttributeValue(VALUE_ATR).equals("ylidyne")){
						outAtomValency =3;
					}
				}
				if (outAtomValency==Integer.parseInt(multiplierBeforeGroup.getAttributeValue(VALUE_ATR))){
					return true;
				}
			}
		}
		return false;
	}


	/**
	 * Calculates number of OutAtoms that the resolveSuffixes method will add.
	 * @param group
	 * @param suffixes
	 * @return numberOfOutAtoms that will be added by resolveSuffixes
	 * @throws ComponentGenerationException
	 */
	private int calculateOutAtomsToBeAddedFromInlineSuffixes(Element group, List<Element> suffixes) throws  ComponentGenerationException {
		int outAtomsThatWillBeAdded = 0;
		Fragment frag = group.getFrag();
		String groupType = frag.getType();
		String subgroupType = frag.getSubType();
		String suffixTypeToUse =null;
		if (suffixRulesLookup.isGroupTypeWithSpecificSuffixRules(groupType)){
			suffixTypeToUse =groupType;
		}
		else{
			suffixTypeToUse = STANDARDGROUP_TYPE_VAL;
		}

		List<Fragment> suffixList =state.xmlSuffixMap.get(group);

		for (Fragment suffix : suffixList) {
			outAtomsThatWillBeAdded += suffix.getOutAtomCount();
		}
		for (Element suffix : suffixes) {
			String suffixValue = suffix.getAttributeValue(VALUE_ATR);
			List<SuffixRule> suffixRules = suffixRulesLookup.getSuffixRuleTags(suffixTypeToUse, suffixValue, subgroupType);
			for (SuffixRule suffixRule : suffixRules) {
				if(suffixRule.getType() == SuffixRuleType.setOutAtom) {
					outAtomsThatWillBeAdded += 1;
				}
			}
		}
		return outAtomsThatWillBeAdded;
	}

	/**
	 * Corrects something like L-alanyl-L-glutaminyl-L-arginyl-O-phosphono-L-seryl-L-alanyl-L-proline to:
	 * ((((L-alanyl-L-glutaminyl)-L-arginyl)-O-phosphono-L-seryl)-L-alanyl)-L-proline
	 * i.e. substituents go onto the last mentioned amino acid; amino acids chain together to form peptides
	 * @param groups
	 * @param brackets
	 */
	private void addImplicitBracketsToAminoAcids(List<Element> groups, List<Element> brackets) {
		for (int i = groups.size() -1; i >=0; i--) {
			Element group = groups.get(i);
			if (group.getAttributeValue(TYPE_ATR).equals(AMINOACID_TYPE_VAL) && OpsinTools.getNextGroup(group)!=null){
				Element possibleLocant = OpsinTools.getPreviousSiblingIgnoringCertainElements(group, new String[]{MULTIPLIER_EL});
				if (possibleLocant != null && possibleLocant.getName().equals(LOCANT_EL)){
					continue;
				}
				Element subOrRoot = group.getParent();
				
				//now find the brackets/substituents before this element
				Element previous = OpsinTools.getPreviousSibling(subOrRoot);
				List<Element> previousElements = new ArrayList<Element>();
				while( previous !=null){
					if (!previous.getName().equals(SUBSTITUENT_EL) && !previous.getName().equals(BRACKET_EL)){
						break;
					}
					previousElements.add(previous);
					previous = OpsinTools.getPreviousSibling(previous);
				}
				if (previousElements.size()>0){//an implicit bracket is needed
					Collections.reverse(previousElements);
					Element bracket = new GroupingEl(BRACKET_EL);
					bracket.addAttribute(new Attribute(TYPE_ATR, IMPLICIT_TYPE_VAL));
					Element parent = subOrRoot.getParent();
					int indexToInsertAt = parent.indexOf(previousElements.get(0));
					for (Element element : previousElements) {
						element.detach();
						bracket.addChild(element);
					}

					subOrRoot.detach();
					bracket.addChild(subOrRoot);
					parent.insertChild(bracket, indexToInsertAt);
					brackets.add(bracket);
				}
			}
		}
	}

	/**
	 * Looks for whether this substituent should be bracketed to the substituent before it
	 * E.g. dimethylaminobenzene -> (dimethylamino)benzene when the substituent is the amino
	 * The list of brackets is modified if the method does something
	 * @param substituent
	 * @param brackets
	 * @throws StructureBuildingException
	 * @throws ComponentGenerationException
	 */
	private void implicitlyBracketToPreviousSubstituentIfAppropriate(Element substituent, List<Element> brackets) throws ComponentGenerationException, StructureBuildingException {
		String firstElInSubName = substituent.getChild(0).getName();
		if (firstElInSubName.equals(LOCANT_EL) ||firstElInSubName.equals(MULTIPLIER_EL)){
			return;
		}

		Element substituentGroup = substituent.getFirstChildElement(GROUP_EL);
		//Only some substituents are valid joiners (e.g. no rings are valid joiners). Need to be atleast bivalent.
		if (substituentGroup.getAttribute(USABLEASJOINER_ATR)==null){
			return;
		}

		//checks that the element before is a substituent or a bracket which will obviously include substituent/s
		//this makes sure there's more than just a substituent in the bracket
		Element elementBeforeSubstituent = OpsinTools.getPreviousSibling(substituent);
		if (elementBeforeSubstituent ==null||
				!elementBeforeSubstituent.getName().equals(SUBSTITUENT_EL) &&
				!elementBeforeSubstituent.getName().equals(BRACKET_EL)){
			return;
		}
		
		Element elementAftersubstituent = OpsinTools.getNextSibling(substituent);
		if (elementAftersubstituent != null){
			//Not preceded and succeded by a bracket e.g. Not (benzyl)methyl(phenyl)amine	c.f. P-16.4.1.3 (draft 2004)
			if (elementBeforeSubstituent.getName().equals(BRACKET_EL) && !IMPLICIT_TYPE_VAL.equals(elementBeforeSubstituent.getAttributeValue(TYPE_ATR)) && elementAftersubstituent.getName().equals(BRACKET_EL)){
				Element firstChildElementOfElementAfterSubstituent = elementAftersubstituent.getChild(0);
				if ((firstChildElementOfElementAfterSubstituent.getName().equals(SUBSTITUENT_EL) || firstChildElementOfElementAfterSubstituent.getName().equals(BRACKET_EL))
					&& !OpsinTools.getPrevious(firstChildElementOfElementAfterSubstituent).getName().equals(HYPHEN_EL)){
					return;
				}
			}
		}
		//there must be an element after the substituent (or the substituent is being used for locanted ester formation) for the implicit bracket to be required
		if (elementAftersubstituent ==null ||
				!elementAftersubstituent.getName().equals(SUBSTITUENT_EL) &&
				!elementAftersubstituent.getName().equals(BRACKET_EL) &&
				!elementAftersubstituent.getName().equals(ROOT_EL)){
			if (!(elementAftersubstituent == null && locantedEsterImplicitBracketSpecialCase(substituent, elementBeforeSubstituent))) {
				return;
			}
		}

		//look for hyphen between substituents, this seems to indicate implicit bracketing was not desired e.g. dimethylaminomethane vs dimethyl-aminomethane
		Element elementDirectlyBeforeSubstituent = OpsinTools.getPrevious(substituent.getChild(0));//can't return null as we know elementBeforeSubstituent is not null
		if (elementDirectlyBeforeSubstituent.getName().equals(HYPHEN_EL)){
			return;
		}
		
		Fragment frag = substituentGroup.getFrag();
		String theSubstituentSubType = substituentGroup.getAttributeValue(SUBTYPE_ATR);
		String theSubstituentType = substituentGroup.getAttributeValue(TYPE_ATR);

		//prevents alkyl chains being bracketed together e.g. ethylmethylamine
		//...unless it's something like 2-methylethyl where the first appears to be locanted onto the second
		List<Element> groupElements = OpsinTools.getDescendantElementsWithTagName(elementBeforeSubstituent, GROUP_EL);//one for a substituent, possibly more for a bracket
		Element lastGroupOfElementBeforeSub =groupElements.get(groupElements.size() - 1);
		if (lastGroupOfElementBeforeSub == null) {
			throw new ComponentGenerationException("No group where group was expected");
		}
		if (theSubstituentType.equals(CHAIN_TYPE_VAL) && theSubstituentSubType.equals(ALKANESTEM_SUBTYPE_VAL) &&
				lastGroupOfElementBeforeSub.getAttributeValue(TYPE_ATR).equals(CHAIN_TYPE_VAL) && lastGroupOfElementBeforeSub.getAttributeValue(SUBTYPE_ATR).equals(ALKANESTEM_SUBTYPE_VAL)){
			boolean placeInImplicitBracket = false;

			Element suffixAfterGroup = OpsinTools.getNextSibling(lastGroupOfElementBeforeSub, SUFFIX_EL);
			//if the alkane ends in oxy, sulfinyl, sulfonyl etc. it's not a pure alkane
			//the outatom check rules out things like "oyl" which don't extend the chain
			if (suffixAfterGroup !=null && suffixAfterGroup.getFrag() != null && suffixAfterGroup.getFrag().getOutAtomCount() > 0){
				placeInImplicitBracket = true;
			}
			//look for locants and check whether they appear to be referring to the other chain
			if (!placeInImplicitBracket){
				List<Element> childrenOfElementBeforeSubstituent = elementBeforeSubstituent.getChildElements();
				Boolean foundLocantNotReferringToChain = null;
				for (Element childOfElBeforeSub : childrenOfElementBeforeSubstituent) {
					String currentElementName = childOfElBeforeSub.getName();
					if (currentElementName.equals(LOCANT_EL)){
						String locantText = childOfElBeforeSub.getValue();
						if(!frag.hasLocant(locantText)){
							foundLocantNotReferringToChain = true;
							break;
						}
						else{
							foundLocantNotReferringToChain = false;
						}
					}
					else if (currentElementName.equals(STEREOCHEMISTRY_EL)){
					}
					else{
						break;
					}
				}
				if (foundLocantNotReferringToChain !=null && !foundLocantNotReferringToChain){//a locant was found and it appeared to refer to the other chain
					placeInImplicitBracket = true;
				}
			}
			if (!placeInImplicitBracket){
				return;
			}
		}

		//prevent bracketing to multi radicals unless through substitution they are likely to cease being multiradicals
		if (lastGroupOfElementBeforeSub.getAttribute(ISAMULTIRADICAL_ATR)!=null && lastGroupOfElementBeforeSub.getAttribute(ACCEPTSADDITIVEBONDS_ATR)==null && lastGroupOfElementBeforeSub.getAttribute(IMINOLIKE_ATR)==null){
			return;
		}
		if (substituentGroup.getAttribute(ISAMULTIRADICAL_ATR)!=null && substituentGroup.getAttribute(ACCEPTSADDITIVEBONDS_ATR)==null && substituentGroup.getAttribute(IMINOLIKE_ATR)==null){
			return;
		}
		if (lastGroupOfElementBeforeSub.getAttribute(IMINOLIKE_ATR)!=null && substituentGroup.getAttribute(IMINOLIKE_ATR)!=null){
			return;//possibly a multiplicative additive operation
		}
		if (substituentGroup.getValue().equals("sulf") && frag.getAtomCount() == 1){
			Element suffix = OpsinTools.getNextSiblingIgnoringCertainElements(substituentGroup, new String[]{UNSATURATOR_EL});
			if (suffix != null && suffix.getAttributeValue(VALUE_ATR).equals("ylidene")) {
				substituentGroup.removeAttribute(substituentGroup.getAttribute(USABLEASJOINER_ATR));
				//TODO resolve suffixes as early as can be done unambiguously
				//e.g. it should be possible to know that sulfanylidene has 0 hydrogen by azaniumylidyne has 1
				return;
			}
		}
		
		//prevent bracketting perhalogeno terms 
		if (PERHALOGENO_SUBTYPE_VAL.equals(lastGroupOfElementBeforeSub.getAttributeValue(SUBTYPE_ATR))){
			return;
		}

		/*
		 * locant may need to be moved. This occurs when the group in elementBeforeSubstituent is not supposed to be locanted onto
		 *  theSubstituentGroup
		 *  e.g. 2-aminomethyl-1-chlorobenzene where the 2 refers to the benzene NOT the methyl
		 */
		List<Element> locantRelatedElements = new ArrayList<Element>();//sometimes moved
		String[] locantValues = null;
		ArrayList<Element> stereoChemistryElements =new ArrayList<Element>();//always moved if bracketing occurs
		List<Element> childrenOfElementBeforeSubstituent = elementBeforeSubstituent.getChildElements();
		for (Element childOfElBeforeSub : childrenOfElementBeforeSubstituent) {
			String currentElementName = childOfElBeforeSub.getName();
			if (currentElementName.equals(STEREOCHEMISTRY_EL)){
				stereoChemistryElements.add(childOfElBeforeSub);
			}
			else if (currentElementName.equals(LOCANT_EL)){
				if (locantValues !=null){
					break;
				}
				locantRelatedElements.add(childOfElBeforeSub);
				locantValues = MATCH_COMMA.split(childOfElBeforeSub.getValue());
			}
			else{
				break;
			}
		}

		//either all locants will be moved, or none
		boolean moveLocants = false;
		if (locantValues != null){
			Element elAfterLocant = OpsinTools.getNextSibling(locantRelatedElements.get(0));
			for (String locantText : locantValues) {
				
				//Check the right fragment in the bracket:
				//if it only has 1 then assume locanted substitution onto it not intended. Or if doesn't have the required locant
				if (frag.getAtomCount()==1 || !frag.hasLocant(locantText) || matchElementSymbolOrAminoAcidLocant.matcher(locantText).find()
						|| (locantValues.length ==1 && elAfterLocant.getName().equals(MULTIPLIER_EL))){
					if (checkLocantPresentOnPotentialRoot(substituent, locantText)){
						moveLocants = true;//locant location is present elsewhere
						break;
					}
					else {
						if( frag.getAtomCount()==1 && frag.hasLocant(locantText)){
							//1 locant was intended to locant onto fragment with 1 atom
						}
						else{
							moveLocants = true;//the fragment adjacent to the locant doesn't have this locant or doesn't need any indirect locants. Assume it will appear elsewhere later
							break;
						}
					}
				}
			}

			if (moveLocants && locantValues.length > 1){
				if (elAfterLocant != null && elAfterLocant.getName().equals(MULTIPLIER_EL)){
					Element shouldBeAGroupOrSubOrBracket = OpsinTools.getNextSiblingIgnoringCertainElements(elAfterLocant, new String[]{MULTIPLIER_EL});
					if (shouldBeAGroupOrSubOrBracket != null){
						if ((shouldBeAGroupOrSubOrBracket.getName().equals(GROUP_EL) && elAfterLocant.getAttributeValue(TYPE_ATR).equals(GROUP_TYPE_VAL))//e.g. 2,5-bisaminothiobenzene --> 2,5-bis(aminothio)benzene
								|| (matchGroupsThatAreAlsoInlineSuffixes.matcher(substituentGroup.getValue()).matches())){//e.g. 4,4'-dimethoxycarbonyl-2,2'-bioxazole --> 4,4'-di(methoxycarbonyl)-2,2'-bioxazole
							locantRelatedElements.add(elAfterLocant);//e.g. 1,5-bis-(4-methylphenyl)sulfonyl --> 1,5-bis-((4-methylphenyl)sulfonyl)
						}
						else if (ORTHOMETAPARA_TYPE_VAL.equals(locantRelatedElements.get(0).getAttributeValue(TYPE_ATR))){//e.g. p-dimethylamino[ring]
							locantRelatedElements.get(0).setValue(locantValues[1]);
						}
						else if (frag.getAtomCount() == 1) {//e.g. 1,3,4-trimethylthiobenzene --> 1,3,4-tri(methylthio)benzene
							locantRelatedElements.add(elAfterLocant);
						}
						else{//don't bracket other complex multiplied substituents (name hasn't given enough hints if indeed bracketing was expected)
							return;
						}
					}
					else{
						moveLocants = false;
					}
				}
				else{
					moveLocants = false;
				}
			}
		}

		Element bracket = new GroupingEl(BRACKET_EL);
		bracket.addAttribute(new Attribute(TYPE_ATR, IMPLICIT_TYPE_VAL));

		for (Element stereoChemistryElement : stereoChemistryElements) {
			stereoChemistryElement.detach();
			bracket.addChild(stereoChemistryElement);
		}
		if (moveLocants){
			for (Element locantElement : locantRelatedElements) {
				locantElement.detach();
				bracket.addChild(locantElement);
			}
		}

		/*
		 * Case when a multiplier should be moved
		 * e.g. tripropan-2-yloxyphosphane -->tri(propan-2-yloxy)phosphane or trispropan-2-ylaminophosphane --> tris(propan-2-ylamino)phosphane
		 */
		if (locantRelatedElements.size() == 0){
			Element possibleMultiplier =childrenOfElementBeforeSubstituent.get(0);
			if (possibleMultiplier.getName().equals(MULTIPLIER_EL) && (
					matchGroupsThatAreAlsoInlineSuffixes.matcher(substituentGroup.getValue()).matches() || possibleMultiplier.getAttributeValue(TYPE_ATR).equals(GROUP_TYPE_VAL))){
				Element desiredGroup = OpsinTools.getNextSiblingIgnoringCertainElements(possibleMultiplier, new String[]{MULTIPLIER_EL});
				if (desiredGroup !=null && desiredGroup.getName().equals(GROUP_EL)){
					childrenOfElementBeforeSubstituent.get(0).detach();
					bracket.addChild(childrenOfElementBeforeSubstituent.get(0));
				}
			}
		}

		Element parent = substituent.getParent();
		int startIndex=parent.indexOf(elementBeforeSubstituent);
		int endIndex=parent.indexOf(substituent);
		for(int i = 0 ; i <= (endIndex-startIndex);i++) {
			Element n = parent.getChild(startIndex);
			n.detach();
			bracket.addChild(n);
		}
		parent.insertChild(bracket, startIndex);
		brackets.add(bracket);
	}


	/**
	 * Retrusn true in the case that:
	 * the given substituent is a direct child of a word element
	 * The preceding substituent/bracket is the first element in the word element
	 * The current word rule invovles locanted ester like linkages
	 * @param substituent
	 * @param elementBeforeSubstituent
	 * @return
	 */
	private boolean locantedEsterImplicitBracketSpecialCase(Element substituent, Element elementBeforeSubstituent) {
		if (substituent.getParent().getName().equals(WORD_EL) &&
				OpsinTools.getPreviousSibling(elementBeforeSubstituent) == null &&
				(state.currentWordRule == WordRule.ester || state.currentWordRule == WordRule.functionalClassEster || state.currentWordRule == WordRule.multiEster || state.currentWordRule == WordRule.acetal)){
			return true;
		}
		return false;
	}

	/** 
	 * Attempts to match locants to non adjacent suffixes/unsatuators
	 * e.g.  2-propanol, 3-furyl, 2'-Butyronaphthone
	 * @param subOrRoot The substituent/root to look for locants in.
	 * @throws StructureBuildingException
	 */
	private void matchLocantsToIndirectFeatures(Element subOrRoot) throws  StructureBuildingException {
		/* Root fragments (or the root in a bracket) can have prefix-locants
		 * that work on suffixes - (2-furyl), 2-propanol, (2-propylmethyl), (2-propyloxy), 2'-Butyronaphthone.
		 */
		List<Element> locantEls = findLocantsThatCouldBeIndirectLocants(subOrRoot);

		if (locantEls.size()>0){
			Element group =subOrRoot.getFirstChildElement(GROUP_EL);
			Element lastLocant = locantEls.get(locantEls.size()-1);//the locant that may apply to an unsaturator/suffix
			String[] locantValues = MATCH_COMMA.split(lastLocant.getValue());
			if (locantValues.length==1 && group.getAttribute(FRONTLOCANTSEXPECTED_ATR)!=null){//some trivial retained names like 2-furyl expect locants to be in front of them. For these the indirect intepretation will always be used rather than checking whether 2-(furyl) even makes sense
				String[] allowedLocants = MATCH_COMMA.split(group.getAttributeValue(FRONTLOCANTSEXPECTED_ATR));
				for (String allowedLocant : allowedLocants) {
					if (locantValues[0].equals(allowedLocant)){
						Element expectedSuffix = OpsinTools.getNextSibling(group);
						if (expectedSuffix!=null && expectedSuffix.getName().equals(SUFFIX_EL) && expectedSuffix.getAttribute(LOCANT_ATR)==null){
							expectedSuffix.addAttribute(new Attribute(LOCANT_ATR, locantValues[0]));
							lastLocant.detach();
							return;
						}
						break;
					}
				}
			}
			boolean allowIndirectLocants =true;
			if(state.currentWordRule == WordRule.multiEster && !ADDEDHYDROGENLOCANT_TYPE_VAL.equals(lastLocant.getAttributeValue(TYPE_ATR))){//special case e.g. 1-benzyl 4-butyl terephthalate (locants do not apply to yls)
				Element parentEl = subOrRoot.getParent();
				if (parentEl.getName().equals(WORD_EL) && parentEl.getAttributeValue(TYPE_ATR).equals(SUBSTITUENT_EL) && parentEl.getChildCount()==1 &&
						locantValues.length==1 && !ORTHOMETAPARA_TYPE_VAL.equals(lastLocant.getAttributeValue(TYPE_ATR))){
					allowIndirectLocants =false;
				}
			}
			Fragment fragmentAfterLocant = group.getFrag();
			if (fragmentAfterLocant.getAtomCount()<=1){
				allowIndirectLocants =false;//e.g. prevent 1-methyl as meth-1-yl is extremely unlikely to be the intended result
			}

			if (allowIndirectLocants){
				/* The first locant is most likely a locant indicating where this subsituent should be attached.
				 * If the locant cannot be found on a potential root this cannot be the case though (assuming the name is valid of course)
				 */
				if (!ADDEDHYDROGENLOCANT_TYPE_VAL.equals(lastLocant.getAttributeValue(TYPE_ATR)) && locantEls.size() ==1 && group.getAttribute(ISAMULTIRADICAL_ATR)==null &&
						locantValues.length == 1 && checkLocantPresentOnPotentialRoot(subOrRoot, locantValues[0]) && OpsinTools.getPreviousSibling(lastLocant, LOCANT_EL)==null){
					return;
				}
				boolean assignableToIndirectFeatures =true;
				List<Element> locantAble  =findElementsMissingIndirectLocants(subOrRoot, lastLocant);
				if (locantAble.size() < locantValues.length){
					assignableToIndirectFeatures =false;
				}
				else{
					for (String locantValue : locantValues) {
						if (!fragmentAfterLocant.hasLocant(locantValue)){//locant is not available on the group adjacent to the locant!
							assignableToIndirectFeatures =false;
						}
					}
				}
				
				if (!assignableToIndirectFeatures){//usually indicates the name will fail unless the suffix has the locant or heteroatom replacement will create the locant
					if (locantValues.length==1){
						List<Fragment> suffixes =state.xmlSuffixMap.get(group);
						//I do not want to assign element locants as in locants on the suffix as I currently know of no examples where this actually occurs
						if (matchElementSymbolOrAminoAcidLocant.matcher(locantValues[0]).matches()){
							return;
						}
						for (Fragment suffix : suffixes) {
							if (suffix.hasLocant(locantValues[0])){//e.g. 2'-Butyronaphthone
								Atom dummyRAtom =suffix.getFirstAtom();
								List<Atom> neighbours =dummyRAtom.getAtomNeighbours();
								Bond b =null;
								atomLoop: for (Atom atom : neighbours) {
									List<String> neighbourLocants = atom.getLocants();
									for (String neighbourLocant : neighbourLocants) {
										if (MATCH_NUMERIC_LOCANT.matcher(neighbourLocant).matches()){
											b = dummyRAtom.getBondToAtomOrThrow(atom);
											break atomLoop;
										}
									}
								}
								if (b!=null){
									state.fragManager.removeBond(b);//the current bond between the dummy R and the suffix
									state.fragManager.createBond(dummyRAtom, suffix.getAtomByLocantOrThrow(locantValues[0]), b.getOrder());
									lastLocant.detach();
								}
							}
						}
					}
				}
				else{
					for (int i = 0; i < locantValues.length; i++) {
						String locantValue = locantValues[i];
						locantAble.get(i).addAttribute(new Attribute(LOCANT_ATR, locantValue));
					}
					lastLocant.detach();
				}
			}
		}
	}


	/**
	 * Finds locants that are before a group element and not immediately followed by a multiplier
	 * @param subOrRoot
	 * @return
	 */
	private List<Element> findLocantsThatCouldBeIndirectLocants(Element subOrRoot) {
		List<Element> children = subOrRoot.getChildElements();
		List<Element> locantEls = new ArrayList<Element>();
		for (Element el : children) {
			if (el.getName().equals(LOCANT_EL)){
				Element afterLocant = OpsinTools.getNextSibling(el);
				if (afterLocant!=null && afterLocant.getName().equals(MULTIPLIER_EL)){//locant should not be followed by a multiplier. c.f. 1,2,3-tributyl 2-acetyloxypropane-1,2,3-tricarboxylate
					continue;
				}
				locantEls.add(el);
			}
			else if (el.getName().equals(GROUP_EL)){
				break;
			}
		}
		return locantEls;
	}


	/**
	 * Find elements that can have indirect locants but don't currently
	 * This requirement excludes hydro and heteroatoms as it is assumed that locants for these are always adjacent (or handled by the special HW code in the case of heteroatoms)
	 * @param subOrRoot The subOrRoot of interest
	 * @param locantEl the locant, only elements after it will be considered
	 * @return An arrayList of locantable elements
	 */
	private List<Element> findElementsMissingIndirectLocants(Element subOrRoot,Element locantEl) {
		List<Element> locantAble = new ArrayList<Element>();
		List<Element> childrenOfSubOrBracketOrRoot=subOrRoot.getChildElements();
		for (Element el : childrenOfSubOrBracketOrRoot) {
			String name =el.getName();
			if (name.equals(SUFFIX_EL) || name.equals(UNSATURATOR_EL) || name.equals(CONJUNCTIVESUFFIXGROUP_EL)){
				if (el.getAttribute(LOCANT_ATR) ==null && el.getAttribute(LOCANTID_ATR) ==null && el.getAttribute(MULTIPLIED_ATR)==null){// shouldn't already have a locant or be multiplied (should of already had locants assignd to it if that were the case)
					if (subOrRoot.indexOf(el)>subOrRoot.indexOf(locantEl)){
						if (name.equals(SUFFIX_EL)){//check a few special cases that must not be locanted
							Element group = OpsinTools.getPreviousSibling(el, GROUP_EL);
							String type = group.getAttributeValue(TYPE_ATR);
							if ((type.equals(ACIDSTEM_TYPE_VAL) && !CYCLEFORMER_SUBTYPE_VAL.equals(el.getAttributeValue(SUBTYPE_ATR)))||
									type.equals(NONCARBOXYLICACID_TYPE_VAL) || type.equals(CHALCOGENACIDSTEM_TYPE_VAL)){
								continue;
							}
						}
						locantAble.add(el);
					}
				}
			}
		}
		return locantAble;
	}


	/**
	 * Put di-carbon modifying suffixes e.g. oic acids, aldehydes on opposite ends of chain
	 * @param subOrRoot
	 * @throws StructureBuildingException
	 */
	private void assignImplicitLocantsToDiTerminalSuffixes(Element subOrRoot) throws StructureBuildingException {
		Element terminalSuffix1 = subOrRoot.getFirstChildElement(SUFFIX_EL);
		if (terminalSuffix1!=null){
			if (isATerminalSuffix(terminalSuffix1) && OpsinTools.getNextSibling(terminalSuffix1) != null){
				Element terminalSuffix2 =OpsinTools.getNextSibling(terminalSuffix1);
				if (isATerminalSuffix(terminalSuffix2)){
					Element hopefullyAChain = OpsinTools.getPreviousSibling(terminalSuffix1, GROUP_EL);
					if (hopefullyAChain != null && hopefullyAChain.getAttributeValue(TYPE_ATR).equals(CHAIN_TYPE_VAL)){
						int chainLength = hopefullyAChain.getFrag().getChainLength();
						if (chainLength >=2){
							terminalSuffix1.addAttribute(new Attribute(LOCANT_ATR, "1"));
							terminalSuffix2.addAttribute(new Attribute(LOCANT_ATR, Integer.toString(chainLength)));
						}
					}
				}
			}
		}
	}


	/**
	 * Checks whether a suffix element is:
	 * a suffix, an inline suffix OR terminal root suffix, has no current locant
	 * @param suffix
	 * @return
	 */
	private boolean isATerminalSuffix(Element suffix){
		return suffix.getName().equals(SUFFIX_EL) &&
				suffix.getAttribute(LOCANT_ATR) == null &&
				(suffix.getAttributeValue(TYPE_ATR).equals(INLINE_TYPE_VAL) || TERMINAL_SUBTYPE_VAL.equals(suffix.getAttributeValue(SUBTYPE_ATR)));
		}

	private void processConjunctiveNomenclature(Element subOrRoot) throws ComponentGenerationException, StructureBuildingException {
		List<Element> conjunctiveGroups = subOrRoot.getChildElements(CONJUNCTIVESUFFIXGROUP_EL);
		if (conjunctiveGroups.size()>0){
			Element ringGroup = subOrRoot.getFirstChildElement(GROUP_EL);
			Fragment ringFrag = ringGroup.getFrag();
			if (ringFrag.getOutAtomCount()!=0 ){
				throw new ComponentGenerationException("OPSIN Bug: Ring fragment should have no radicals");
			}
			List<Fragment> conjunctiveFragments = new ArrayList<Fragment>();
			for (Element group : conjunctiveGroups) {
				Fragment frag = group.getFrag();
				conjunctiveFragments.add(frag);
			}
			for (int i = 0; i < conjunctiveFragments.size(); i++) {
				Fragment conjunctiveFragment = conjunctiveFragments.get(i);
				if (conjunctiveGroups.get(i).getAttribute(LOCANT_ATR)!=null){
					state.fragManager.createBond(lastNonSuffixCarbonWithSufficientValency(conjunctiveFragment), ringFrag.getAtomByLocantOrThrow(conjunctiveGroups.get(i).getAttributeValue(LOCANT_ATR)) , 1);
				}
				else{
					state.fragManager.createBond(lastNonSuffixCarbonWithSufficientValency(conjunctiveFragment), ringFrag.getAtomOrNextSuitableAtomOrThrow(ringFrag.getFirstAtom(), 1, true) , 1);
				}
				state.fragManager.incorporateFragment(conjunctiveFragment, ringFrag);
			}
		}
	}


	private Atom lastNonSuffixCarbonWithSufficientValency(Fragment conjunctiveFragment) throws ComponentGenerationException {
		List<Atom> atomList = conjunctiveFragment.getAtomList();
		for (int i = atomList.size()-1; i >=0; i--) {
			Atom a = atomList.get(i);
			if (a.getType().equals(SUFFIX_TYPE_VAL)){
				continue;
			}
			if (a.getElement() != ChemEl.C){
				continue;
			}
			if (ValencyChecker.checkValencyAvailableForBond(a, 1)){
				return a;
			}
		}
		throw new ComponentGenerationException("OPSIN Bug: Unable to find non suffix carbon with sufficient valency");
	}


	/**Process the effects of suffixes upon a fragment. 
	 * Unlocanted non-terminal suffixes are not attached yet. All other suffix effects are performed
	 * @param group The group element for the fragment to which the suffixes will be added
	 * @param suffixes The suffix elements for a fragment.
	 * @throws StructureBuildingException If the suffixes can't be resolved properly.
	 * @throws ComponentGenerationException
	 */
	private void resolveSuffixes(Element group, List<Element> suffixes) throws StructureBuildingException, ComponentGenerationException {
		Fragment frag = group.getFrag();
		List<Atom> atomList = frag.getAtomList();//this instance of atomList will not change even once suffixes are merged into the fragment
		int defaultAtom = 0;//index in atomList
		String groupType = frag.getType();
		String subgroupType = frag.getSubType();
		String suffixTypeToUse = suffixRulesLookup.isGroupTypeWithSpecificSuffixRules(groupType) ? groupType : STANDARDGROUP_TYPE_VAL;
		
		List<Fragment> associatedSuffixFrags = state.xmlSuffixMap.get(group);
		if (associatedSuffixFrags != null) {//null for non-final group in polycyclic spiro systems
			associatedSuffixFrags.clear();
		}
		for (Element suffix : suffixes) {
			String suffixValue = suffix.getAttributeValue(VALUE_ATR);
			
			Fragment suffixFrag = null;
			List<SuffixRule> suffixRules = suffixRulesLookup.getSuffixRuleTags(suffixTypeToUse, suffixValue, subgroupType);
			for (SuffixRule suffixRule : suffixRules) {
				if (defaultAtom >= atomList.size()) {
					defaultAtom = 0;
				}
				switch (suffixRule.getType()) {
				case addgroup:
					if (suffixFrag == null) {
						suffixFrag = suffix.getFrag();
						if (suffixFrag == null) {
							throw new RuntimeException("OPSIN Bug: Suffix was expected to have an associated fragment but it wasn't found");
						}
						Atom firstAtomInSuffix = suffixFrag.getFirstAtom();
						if (firstAtomInSuffix.getBondCount() <= 0) {
							throw new ComponentGenerationException("OPSIN Bug: Dummy atom in suffix should have at least one bond to it");
						}
						if (CYCLEFORMER_SUBTYPE_VAL.equals(suffix.getAttributeValue(SUBTYPE_ATR))){
							processCycleFormingSuffix(suffixFrag, frag, suffix);
						}
						else{
							int bondOrderRequired = firstAtomInSuffix.getIncomingValency();
							Atom fragAtomToUse = getFragAtomToUse(frag, suffix, suffixTypeToUse);
							if (fragAtomToUse == null) {
								if ("yes".equals(suffixRule.getAttributeValue(SUFFIXRULES_KETONELOCANT_ATR)) && !atomList.get(defaultAtom).getAtomIsInACycle()) {
									if (defaultAtom == 0) {
										defaultAtom = FragmentTools.findKetoneAtomIndice(frag, defaultAtom);
									}
									fragAtomToUse = atomList.get(defaultAtom);
									defaultAtom++;
								} else {
									fragAtomToUse = atomList.get(defaultAtom);
								}
								fragAtomToUse = frag.getAtomOrNextSuitableAtomOrThrow(fragAtomToUse, bondOrderRequired, true);
								if (FragmentTools.isCharacteristicAtom(fragAtomToUse)){
									throw new StructureBuildingException("No suitable atom found to attach suffix");
								}
							}
	
							//create a new bond and associate it with the suffixfrag and both atoms. Remember the suffixFrag has not been imported into the frag yet
							List<Bond> bonds = new ArrayList<Bond>(firstAtomInSuffix.getBonds());
							for (Bond bondToSuffix : bonds) {
								Atom suffixAtom = bondToSuffix.getOtherAtom(firstAtomInSuffix);
								state.fragManager.createBond(fragAtomToUse, suffixAtom, bondToSuffix.getOrder());
								state.fragManager.removeBond(bondToSuffix);
								if (fragAtomToUse.getIncomingValency() > 2 && (suffixValue.equals("aldehyde") || suffixValue.equals("al")|| suffixValue.equals("aldoxime"))){//formaldehyde/methanal are excluded as they are substitutable
									if("X".equals(suffixAtom.getFirstLocant())){//carbaldehyde
										suffixAtom.setProperty(Atom.ISALDEHYDE, true);
									}
									else{
										fragAtomToUse.setProperty(Atom.ISALDEHYDE, true);
									}
								}
							}
						}
					}
					else{
						throw new ComponentGenerationException("OPSIN bug: Suffix may only have one addgroup rule: " + suffix.getValue());
					}
					break;
				case changecharge:
					int chargeChange = Integer.parseInt(suffixRule.getAttributeValue(SUFFIXRULES_CHARGE_ATR));
					int protonChange = Integer.parseInt(suffixRule.getAttributeValue(SUFFIXRULES_PROTONS_ATR));
					if (suffix.getAttribute(SUFFIXPREFIX_ATR) == null) {
						Atom fragAtomToUse = getFragAtomToUse(frag, suffix, suffixTypeToUse);
						if (fragAtomToUse != null) {
							fragAtomToUse.addChargeAndProtons(chargeChange, protonChange);
						}
						else{
							applyUnlocantedChargeModification(atomList, chargeChange, protonChange);
						}
					}
					else {//a suffix prefixed acylium suffix
						if (suffixFrag == null) {
							throw new StructureBuildingException("OPSIN bug: ordering of elements in suffixRules.xml wrong; changeCharge found before addGroup");
						}
						Set<Bond> bonds = state.fragManager.getInterFragmentBonds(suffixFrag);
						if (bonds.size() != 1) {
							throw new StructureBuildingException("OPSIN bug: Wrong number of bonds between suffix and group");
						}
						for (Bond bond : bonds) {
							if (bond.getFromAtom().getFrag() == suffixFrag) {
								bond.getFromAtom().addChargeAndProtons(chargeChange, protonChange);
							} else {
								bond.getToAtom().addChargeAndProtons(chargeChange, protonChange);
							}
						}
					}
					break;
				case setOutAtom:
					String outValencyAtr = suffixRule.getAttributeValue(SUFFIXRULES_OUTVALENCY_ATR);
					int outValency = outValencyAtr != null ? Integer.parseInt(outValencyAtr) : 1;
					if (suffix.getAttribute(SUFFIXPREFIX_ATR) == null) {
						Atom fragAtomToUse = getFragAtomToUse(frag, suffix, suffixTypeToUse);
						if (fragAtomToUse != null) {
							frag.addOutAtom(fragAtomToUse, outValency, true);
						} else {
							frag.addOutAtom(frag.getFirstAtom(), outValency, false);
						}
					} else {//something like oyl on a ring, which means it is now carbonyl and the outAtom is on the suffix and not frag
						if (suffixFrag == null) {
							throw new StructureBuildingException("OPSIN bug: ordering of elements in suffixRules.xml wrong; setOutAtom found before addGroup");
						}
						Set<Bond> bonds = state.fragManager.getInterFragmentBonds(suffixFrag);
						if (bonds.size() != 1) {
							throw new StructureBuildingException("OPSIN bug: Wrong number of bonds between suffix and group");
						}
						for (Bond bond : bonds) {
							if (bond.getFromAtom().getFrag() == suffixFrag) {
								suffixFrag.addOutAtom(bond.getFromAtom(), outValency, true);
							} else {
								suffixFrag.addOutAtom(bond.getToAtom(), outValency, true);
							}
						}
					}
					break;
				case setAcidicElement:
					ChemEl chemEl = ChemEl.valueOf(suffixRule.getAttributeValue(SUFFIXRULES_ELEMENT_ATR));
					swapElementsSuchThatThisElementIsAcidic(suffixFrag, chemEl);
					break;
				case addSuffixPrefixIfNonePresentAndCyclic:
				case addFunctionalAtomsToHydroxyGroups:
				case chargeHydroxyGroups:
				case removeTerminalOxygen:
				case convertHydroxyGroupsToOutAtoms:
				case convertHydroxyGroupsToPositiveCharge:
					//already processed
					break;
				}
			}

			if (suffixFrag != null) {//merge suffix frag and parent fragment
				state.fragManager.removeAtomAndAssociatedBonds(suffixFrag.getFirstAtom());//the dummy R atom
				Set<String> suffixLocants = new HashSet<String>(suffixFrag.getLocants());
				for (String suffixLocant : suffixLocants) {
					if (Character.isDigit(suffixLocant.charAt(0))){//check that numeric locants do not conflict with the parent fragment e.g. hydrazide 2' with biphenyl 2'
						if (frag.hasLocant(suffixLocant)){
							suffixFrag.getAtomByLocant(suffixLocant).removeLocant(suffixLocant);
						}
					}
				}
				state.fragManager.incorporateFragment(suffixFrag, frag);
				if (CYCLEFORMER_SUBTYPE_VAL.equals(suffix.getAttributeValue(SUBTYPE_ATR))){
					CycleDetector.assignWhetherAtomsAreInCycles(frag);
				}
			}
		}
	}

	private Atom getFragAtomToUse(Fragment frag, Element suffix, String suffixTypeToUse) throws StructureBuildingException {
		Atom fragAtomToUse = null;
		String locant = suffix.getAttributeValue(LOCANT_ATR);
		if (locant != null) {
			fragAtomToUse = frag.getAtomByLocantOrThrow(locant);
		}
		else{
			String locantId = suffix.getAttributeValue(LOCANTID_ATR);
			if (locantId != null) {
				fragAtomToUse = frag.getAtomByIDOrThrow(Integer.parseInt(locantId));
			}
			else{
				String defaultLocantId = suffix.getAttributeValue(DEFAULTLOCANTID_ATR);
				if (defaultLocantId != null) {
					fragAtomToUse = frag.getAtomByIDOrThrow(Integer.parseInt(defaultLocantId));
				}
				else if (suffixTypeToUse.equals(ACIDSTEM_TYPE_VAL) || suffixTypeToUse.equals(NONCARBOXYLICACID_TYPE_VAL) || suffixTypeToUse.equals(CHALCOGENACIDSTEM_TYPE_VAL)) {//means that e.g. sulfonyl, has an explicit outAtom
					fragAtomToUse = frag.getFirstAtom();
				}
			}
		}
		return fragAtomToUse;
	}


	private void processCycleFormingSuffix(Fragment suffixFrag, Fragment suffixableFragment, Element suffix) throws StructureBuildingException, ComponentGenerationException {
		List<Atom> rAtoms = new ArrayList<Atom>();
		for (Atom a : suffixFrag.getAtomList()) {
			if (a.getElement() == ChemEl.R){
				rAtoms.add(a);
			}
		}
		if (rAtoms.size() != 2){
			throw new ComponentGenerationException("OPSIN bug: Incorrect number of R atoms associated with cyclic suffix");
		}
		if (rAtoms.get(0).getBondCount() <= 0 || rAtoms.get(1).getBondCount() <= 0) {
			throw new ComponentGenerationException("OPSIN Bug: Dummy atoms in suffix should have at least one bond to them");
		}
		
		Atom parentAtom1;
		Atom parentAtom2;

		String locant = suffix.getAttributeValue(LOCANT_ATR);
		String locantId = suffix.getAttributeValue(LOCANTID_ATR);
		if (locant != null){
			String[] locants = MATCH_COMMA.split(locant);
			if (locants.length ==2){
				parentAtom1 = suffixableFragment.getAtomByLocantOrThrow(locants[0]);
				parentAtom2 = suffixableFragment.getAtomByLocantOrThrow(locants[1]);
			}
			else if (locants.length ==1){
				parentAtom1 = suffixableFragment.getAtomByLocantOrThrow("1");
				parentAtom2 = suffixableFragment.getAtomByLocantOrThrow(locants[0]);
			}
			else{
				throw new ComponentGenerationException("Incorrect number of locants associated with cycle forming suffix, expected 2 found: " + locants.length);
			}
		}
		else if (locantId !=null) {
			String[] locantIds = MATCH_COMMA.split(locantId);
			if (locantIds.length !=2){
				throw new ComponentGenerationException("OPSIN bug: Should be exactly 2 locants associated with a cyclic suffix");
			}
			int firstIdInFragment = suffixableFragment.getIdOfFirstAtom(); 
			parentAtom1 = suffixableFragment.getAtomByIDOrThrow(firstIdInFragment + Integer.parseInt(locantIds[0]) -1);
			parentAtom2 = suffixableFragment.getAtomByIDOrThrow(firstIdInFragment + Integer.parseInt(locantIds[1]) -1);
		}
		else{
			int chainLength = suffixableFragment.getChainLength();
			if (chainLength > 1 && chainLength == suffixableFragment.getAtomCount()){
				parentAtom1 = suffixableFragment.getAtomByLocantOrThrow("1");
				parentAtom2 = suffixableFragment.getAtomByLocantOrThrow(String.valueOf(chainLength));
			}
			else{
				throw new ComponentGenerationException("cycle forming suffix: " + suffix.getValue() +" should be locanted!");
			}
		}
		if (parentAtom1.equals(parentAtom2)){
			throw new ComponentGenerationException("cycle forming suffix: " + suffix.getValue() +" attempted to form a cycle involving the same atom twice!");
		}
		
		if (suffixableFragment.getType().equals(CARBOHYDRATE_TYPE_VAL)){
			removeTerminalOxygen(parentAtom1, 2);
			removeTerminalOxygen(parentAtom1, 1);
			List<Atom> chainHydroxy = FragmentTools.findHydroxyLikeTerminalAtoms(parentAtom2.getAtomNeighbours(), ChemEl.O);
			if (chainHydroxy.size() == 1){
				FragmentTools.removeTerminalAtom(state, chainHydroxy.get(0));//make sure to retain stereochemistry
			}
			else{
				throw new ComponentGenerationException("The second locant of a carbohydrate lactone should point to a carbon in the chain with a hydroxyl group");
			}
		}
		else{
			if (parentAtom2.getElement() == ChemEl.O){//cyclic suffixes like lactone formally indicate the removal of hydroxy cf. 1979 rule 472.1
				//...although in most cases they are used on structures that don't actually have a hydroxy group
				List<Atom> neighbours = parentAtom2.getAtomNeighbours();
				if (neighbours.size()==1){
					List<Atom> suffixNeighbours = rAtoms.get(1).getAtomNeighbours();
					if (suffixNeighbours.size()==1 && suffixNeighbours.get(0).getElement() == ChemEl.O){
						state.fragManager.removeAtomAndAssociatedBonds(parentAtom2);
						parentAtom2 = neighbours.get(0);
					}
				}
			}
		}
		makeBondsToSuffix(parentAtom1, rAtoms.get(0));
		makeBondsToSuffix(parentAtom2, rAtoms.get(1));
		state.fragManager.removeAtomAndAssociatedBonds(rAtoms.get(1));
	}
	
	/**
	 * Creates bonds between the parentAtom and the atoms connected to the R atoms.
	 * Removes bonds to the R atom
	 * @param parentAtom
	 * @param suffixRAtom
	 */
	private void makeBondsToSuffix(Atom parentAtom, Atom suffixRAtom) {
		List<Bond> bonds = new ArrayList<Bond>(suffixRAtom.getBonds());
		for (Bond bondToSuffix : bonds) {
			Atom suffixAtom = bondToSuffix.getOtherAtom(suffixRAtom);
			state.fragManager.createBond(parentAtom, suffixAtom, bondToSuffix.getOrder());
			state.fragManager.removeBond(bondToSuffix);
		}
	}

	/**
	 * Preference is given to mono cation/anions as they are expected to be more likely
	 * Additionally, Typically if a locant has not been specified then it was intended to refer to a nitrogen even if the nitrogen is not at locant 1 e.g. isoquinolinium
	 * Hence preference is given to nitrogen atoms and then to non carbon atoms
	 * @param atomList
	 * @param chargeChange
	 * @param protonChange
	 */
	private void applyUnlocantedChargeModification(List<Atom> atomList, int chargeChange, int protonChange) {
		Atom likelyAtom = null;
		Atom possibleHeteroatom = null;
		Atom possibleCarbonAtom = null;
		Atom possibleDiOrHigherIon = null;
		for (Atom a : atomList) {
			ChemEl chemEl = a.getElement();
			Integer[] stableValencies = ValencyChecker.getPossibleValencies(chemEl, a.getCharge() + chargeChange);
			if (stableValencies == null) {//unstable valency so seems unlikely
				continue;
			}
			int resultantExpectedValency = (a.getLambdaConventionValency() ==null ? ValencyChecker.getDefaultValency(chemEl) : a.getLambdaConventionValency()) + a.getProtonsExplicitlyAddedOrRemoved() + protonChange;
			boolean matched = false;
			for (Integer stableValency : stableValencies) {
				if (stableValency ==resultantExpectedValency){
					matched =true;
					break;
				}
			}
			if (!matched){//unstable valency so seems unlikely
				continue;
			}
			if (protonChange <0 && StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(a)<=0){
				continue;
			}
			if (Math.abs(a.getCharge())==0){
				if (chemEl == ChemEl.N){
					likelyAtom = a;
					break;
				}
				else if (possibleHeteroatom ==null && chemEl != ChemEl.C){
					possibleHeteroatom= a;
				}
				else if (possibleCarbonAtom ==null){
					possibleCarbonAtom = a;
				}
			}
			else if (possibleDiOrHigherIon ==null){
				possibleDiOrHigherIon = a;
			}
		}
		if (likelyAtom == null) {
			if (possibleHeteroatom !=null){
				likelyAtom = possibleHeteroatom;
			}
			else if (possibleCarbonAtom !=null){
				likelyAtom = possibleCarbonAtom;
			}
			else if (possibleDiOrHigherIon !=null){
				likelyAtom = possibleDiOrHigherIon;
			}
			else{
				likelyAtom = atomList.get(0);
			}
		}
		likelyAtom.addChargeAndProtons(chargeChange, protonChange);
	}
	
	/**
	 * e.g. if element is "S" changes C(=S)O -->C(=O)S
	 * @param frag
	 * @param chemEl
	 * @throws StructureBuildingException 
	 */
	private void swapElementsSuchThatThisElementIsAcidic(Fragment frag, ChemEl chemEl) throws StructureBuildingException {
		for (int i = 0, l =frag.getFunctionalAtomCount(); i < l; i++) {
			Atom atom = frag.getFunctionalAtom(i).getAtom();
			Set<Atom> ambiguouslyElementedAtoms = atom.getProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT);
			if (ambiguouslyElementedAtoms != null) {
				Atom atomToSwapWith = null;
				for (Atom ambiguouslyElementedAtom : ambiguouslyElementedAtoms) {
					if (ambiguouslyElementedAtom.getElement() == chemEl){
						atomToSwapWith = ambiguouslyElementedAtom;
						break;	
					}
				}
				if (atomToSwapWith != null) {
					if (atomToSwapWith != atom) {
						//swap locants and element type
						List<String> tempLocants1 = new ArrayList<String>(atom.getLocants());
						List<String> tempLocants2 = new ArrayList<String>(atomToSwapWith.getLocants());
						atom.clearLocants();
						atomToSwapWith.clearLocants();
						for (String locant : tempLocants1) {
							atomToSwapWith.addLocant(locant);
						}
						for (String locant : tempLocants2) {
							atom.addLocant(locant);
						}
						ChemEl a2ChemEl = atomToSwapWith.getElement();
						atomToSwapWith.setElement(atom.getElement());
						atom.setElement(a2ChemEl);
						ambiguouslyElementedAtoms.remove(atomToSwapWith);
					}
					ambiguouslyElementedAtoms.remove(atom);
					return;
				}
			}
		}
		throw new StructureBuildingException("Unable to find potential acidic atom with element: " + chemEl);
	}
	
	/**
	 * Converts a biochemical linkage description e.g. (1->4) into an O[1-9] locant
	 * If the carbohydrate is preceded by substituents these are placed into a bracket and the bracket locanted
	 * @param substituents
	 * @param brackets
	 * @throws StructureBuildingException
	 */
	private void processBiochemicalLinkageDescriptors(List<Element> substituents, List<Element> brackets) throws StructureBuildingException {
		for (Element substituent : substituents) {
			List<Element> bioLinkLocants = substituent.getChildElements(BIOCHEMICALLINKAGE_EL);
			if (bioLinkLocants.size() > 0){
				if (bioLinkLocants.size() > 1){
					throw new RuntimeException("OPSIN Bug: More than 1 biochemical linkage locant associated with subsituent");
				}
				Element bioLinkLocant = bioLinkLocants.get(0);
				String bioLinkLocantStr = bioLinkLocant.getValue();
				bioLinkLocantStr = bioLinkLocantStr.substring(1, bioLinkLocantStr.length() -1);//strip brackets
				checkAndApplyFirstLocantOfBiochemicalLinkage(substituent, bioLinkLocantStr);
				int secondLocantStartPos = Math.max(bioLinkLocantStr.lastIndexOf('>'), bioLinkLocantStr.lastIndexOf('-')) + 1;
				String locantToConnectTo = bioLinkLocantStr.substring(secondLocantStartPos);
				Element parent = substituent.getParent();
				Attribute locantAtr = new Attribute(LOCANT_ATR, "O" + locantToConnectTo);

				Element elementAfterSubstituent = OpsinTools.getNextSibling(substituent);				
				boolean hasAdjacentGroupToSubstitute = (elementAfterSubstituent !=null &&
						(elementAfterSubstituent.getName().equals(SUBSTITUENT_EL) ||
						elementAfterSubstituent.getName().equals(BRACKET_EL) ||
						elementAfterSubstituent.getName().equals(ROOT_EL)));
				

				/* If a biochemical is not at the end of a scope but is preceded by substituents/brackets
				 * these are bracketted and the locant assigned to the bracket.
				 * Else If the group is the only thing in a bracket the locant is assigned to the bracket (this is used to describe branches)
				 * Else the locant is assigned to the substituent
				 */
				boolean bracketAdded =false;
				if (hasAdjacentGroupToSubstitute){
					//now find the brackets/substituents before this element
					Element previous = OpsinTools.getPreviousSibling(substituent);
					List<Element> previousElements = new ArrayList<Element>();
					while( previous !=null){
						if (!previous.getName().equals(SUBSTITUENT_EL) && !previous.getName().equals(BRACKET_EL)){
							break;
						}
						previousElements.add(previous);
						previous = OpsinTools.getPreviousSibling(previous);
					}
					if (previousElements.size() > 0 ){//an explicit bracket is needed
						Collections.reverse(previousElements);
						Element bracket = new GroupingEl(BRACKET_EL);
						bracket.addAttribute(locantAtr);
						int indexToInsertAt = parent.indexOf(previousElements.get(0));
						for (Element element : previousElements) {
							element.detach();
							bracket.addChild(element);
						}

						substituent.detach();
						bracket.addChild(substituent);
						parent.insertChild(bracket, indexToInsertAt);
						brackets.add(bracket);
						bracketAdded = true;
					}
				}
				
				if (!bracketAdded) {
					Element elToAddAtrTo;
					if (parent.getName().equals(BRACKET_EL) && !hasAdjacentGroupToSubstitute){
						elToAddAtrTo = parent;
					}
					else{
						elToAddAtrTo = substituent;
					}
					if (elToAddAtrTo.getAttribute(LOCANT_ATR) !=null){
						throw new StructureBuildingException("Substituent with biochemical linkage descriptor should not also have a locant: " + elToAddAtrTo.getAttributeValue(LOCANT_ATR));
					}
					elToAddAtrTo.addAttribute(locantAtr);
				}
				bioLinkLocant.detach();
			}
		}
		
		for (Element bracket : brackets) {
			List<Element> bioLinkLocants = bracket.getChildElements(BIOCHEMICALLINKAGE_EL);
			if (bioLinkLocants.size() > 0){
				if (bioLinkLocants.size() > 1){
					throw new RuntimeException("OPSIN Bug: More than 1 biochemical linkage locant associated with bracket");
				}
				Element bioLinkLocant = bioLinkLocants.get(0);
				Element substituent = OpsinTools.getPreviousSibling(bioLinkLocant);
				if (substituent == null || !substituent.getName().equals(SUBSTITUENT_EL)){
					throw new RuntimeException("OPSIN Bug: Substituent expected before biochemical linkage locant");
				}
				String bioLinkLocantStr = bioLinkLocant.getValue();
				bioLinkLocantStr = bioLinkLocantStr.substring(1, bioLinkLocantStr.length() -1);
				checkAndApplyFirstLocantOfBiochemicalLinkage(substituent, bioLinkLocantStr);
				int secondLocantStartPos = Math.max(bioLinkLocantStr.lastIndexOf('>'), bioLinkLocantStr.lastIndexOf('-')) + 1;
				String locantToConnectTo = bioLinkLocantStr.substring(secondLocantStartPos);
				if (bracket.getAttribute(LOCANT_ATR) !=null){
					throw new StructureBuildingException("Substituent with biochemical linkage descriptor should not also have a locant: " + bracket.getAttributeValue(LOCANT_ATR));
				}
				bracket.addAttribute(new Attribute(LOCANT_ATR, "O" + locantToConnectTo));
				bioLinkLocant.detach();
			}
		}
	}

	private void checkAndApplyFirstLocantOfBiochemicalLinkage(Element substituent, String biochemicalLinkage) throws StructureBuildingException {
		Element group = substituent.getFirstChildElement(GROUP_EL);
		Fragment frag = group.getFrag();
		String firstLocant = biochemicalLinkage.substring(0, biochemicalLinkage.indexOf('-'));
		if (group.getAttributeValue(TYPE_ATR).equals(CARBOHYDRATE_TYPE_VAL)) {
			Atom anomericAtom = frag.getAtomByLocantOrThrow(firstLocant);
			boolean anomericIsOutAtom = false;
			for (int i = 0; i < frag.getOutAtomCount(); i++) {
				if (frag.getOutAtom(i).getAtom().equals(anomericAtom)){
					anomericIsOutAtom = true;
				}
			}
			if (!anomericIsOutAtom){
				throw new StructureBuildingException("Invalid glycoside linkage descriptor. Locant: " + firstLocant + " should point to the anomeric carbon");
			}
		}
		else{
			Atom positionOfPhospho = frag.getAtomByLocantOrThrow("O" + firstLocant);
			if (positionOfPhospho.getBondCount() !=1){
				throw new StructureBuildingException(firstLocant + " should be the carbon to which a hydroxy group is attached!");
			}
			if (frag.getOutAtomCount()==1){
				Atom atomToConnect = frag.getOutAtom(0).getAtom();
				state.fragManager.createBond(positionOfPhospho, atomToConnect, 1);
			}
			else{
				throw new RuntimeException("OPSIN Bug: Biochemical linkage only expected on groups with 1 OutAtom");
			}
		}
		
		if (OpsinTools.getNextGroup(group)==null){
			throw new StructureBuildingException("Biochemical linkage descriptor should be followed by another biochemical: " + biochemicalLinkage);
		}
	}

	/**
	 * Moves a multiplier out of a bracket if the bracket contains only one substituent
	 * e.g. (trimethyl) --> tri(methyl).
	 * The multiplier may have locants e.g. [N,N-bis(2-hydroxyethyl)]
	 * This is done because OPSIN has no idea what to do with (trimethyl) as there is nothing within the scope to substitute onto!
	 * @param brackets
	 */
	private void moveErroneouslyPositionedLocantsAndMultipliers(List<Element> brackets) {
		for (int i = brackets.size()-1; i >=0; i--) {
			Element bracket =brackets.get(i);
			List<Element> childElements = bracket.getChildElements();
			boolean hyphenPresent = false;
			int childCount = childElements.size();
			if (childCount==2){
				for (int j = childCount -1; j >=0; j--) {
					if (childElements.get(j).getName().equals(HYPHEN_EL)){
						hyphenPresent=true;
					}
				}
			}
			if (childCount==1 || hyphenPresent && childCount==2){
				List<Element> substituentContent = childElements.get(0).getChildElements();
				if (substituentContent.size()>=2){
					Element locant =null;
					Element multiplier =null;
					Element possibleMultiplier = substituentContent.get(0);
					if (substituentContent.get(0).getName().equals(LOCANT_EL)){//probably erroneous locant
						locant = substituentContent.get(0);
						possibleMultiplier = substituentContent.get(1);
					}
					if (possibleMultiplier.getName().equals(MULTIPLIER_EL)){//erroneously placed multiplier present
						multiplier = possibleMultiplier;
					}
					if (locant!=null){
						if (multiplier==null || MATCH_COMMA.split(locant.getValue()).length == Integer.parseInt(multiplier.getAttributeValue(VALUE_ATR))){
							locant.detach();
							OpsinTools.insertBefore(childElements.get(0), locant);
						}
						else{
							continue;
						}
					}
					if (multiplier !=null){
						multiplier.detach();
						OpsinTools.insertBefore(childElements.get(0), multiplier);
					}
				}
			}
		}
	}


	/**
	 * Given the right most child of a word:
	 * Checks whether this is multiplied e.g. methylenedibenzene
	 * If it is then it checks for the presence of locants e.g. 4,4'-oxydibenzene which has been changed to oxy-4,4'-dibenzene
	 * An attribute called inLocants is then added that is either INLOCANTS_DEFAULT or a comma seperated list of locants
	 * @param rightMostElement
	 * @throws ComponentGenerationException
	 * @throws StructureBuildingException
	 */
	private void assignLocantsToMultipliedRootIfPresent(Element rightMostElement) throws ComponentGenerationException, StructureBuildingException {
		List<Element> multipliers = rightMostElement.getChildElements(MULTIPLIER_EL);
		if(multipliers.size() == 1) {
			Element multiplier =multipliers.get(0);
			if (OpsinTools.getPrevious(multiplier)==null){
				throw new StructureBuildingException("OPSIN bug: Unacceptable input to function");
			}
			List<Element> locants = rightMostElement.getChildElements(MULTIPLICATIVELOCANT_EL);
			if (locants.size()>1){
				throw new ComponentGenerationException("OPSIN bug: Only none or one multiplicative locant expected");
			}
			int multiVal = Integer.parseInt(multiplier.getAttributeValue(VALUE_ATR));
			if (locants.size()==0){
				rightMostElement.addAttribute(new Attribute(INLOCANTS_ATR, INLOCANTS_DEFAULT));
			}
			else{
				Element locantEl = locants.get(0);
				String[] locantValues = MATCH_COMMA.split(locantEl.getValue());
				if (locantValues.length == multiVal){
					rightMostElement.addAttribute(new Attribute(INLOCANTS_ATR, locantEl.getValue()));
					locantEl.detach();
				}
				else{
					throw new ComponentGenerationException("Mismatch between number of locants and number of roots");
				}
			}
		}
		else if (rightMostElement.getName().equals(BRACKET_EL)){
			assignLocantsToMultipliedRootIfPresent(rightMostElement.getChild(rightMostElement.getChildCount()-1));
		}
	}


	/**
	 * Adds an implicit bracket in the case where two locants have been given.
	 * One for the locanting of substituent on to the next substituent and one
	 * for the locanting of this combined substituent onto a parent group
	 * e.g. 5-p-hydroxyphenyl-1,2-dithiole-3-thione --> e.g. 5-(p-hydroxyphenyl)-1,2-dithiole-3-thione
	 * @param substituent
	 * @param brackets
	 */
	private void addImplicitBracketsInCaseWhereSubstituentHasTwoLocants(Element substituent, List<Element> brackets) {
		Element siblingSubstituent = OpsinTools.getNextSibling(substituent);
		if (siblingSubstituent != null && siblingSubstituent.getName().equals(SUBSTITUENT_EL)){
			List<Element> locants = getLocantsAtStartOfSubstituent(substituent);
			if (locants.size() == 2 && locantsAreSingular(locants)
					&& getLocantsAtStartOfSubstituent(siblingSubstituent).size() == 0){
				Element bracket = new GroupingEl(BRACKET_EL);
				bracket.addAttribute(new Attribute(TYPE_ATR, IMPLICIT_TYPE_VAL));
				Element parent = substituent.getParent();
				int indexToInsertAt = parent.indexOf(substituent);
				int elementsToMove = substituent.indexOf(locants.get(0)) + 1;
				for (int i = 0; i < elementsToMove; i++) {
					Element locantOrStereoToMove = substituent.getChild(0);
					locantOrStereoToMove.detach();
					bracket.addChild(locantOrStereoToMove);
				}
				substituent.detach();
				siblingSubstituent.detach();
				bracket.addChild(substituent);
				bracket.addChild(siblingSubstituent);
				parent.insertChild(bracket, indexToInsertAt);
				brackets.add(bracket);
			}
		}
	}

	/**
	 * Retrieves the first elements of a substituent which are locants skipping over stereochemistry elements
	 * @param substituent
	 * @return
	 */
	private List<Element> getLocantsAtStartOfSubstituent(Element substituent) {
		List<Element> locants = new ArrayList<Element>();
		List<Element> children = substituent.getChildElements();
		for (Element child : children) {
			String currentElementName = child.getName();
			if (currentElementName.equals(LOCANT_EL)){
				locants.add(child);
			}
			else if (currentElementName.equals(STEREOCHEMISTRY_EL)){
				//ignore
			}
			else{
				break;
			}
		}
		return locants;
	}

	/**
	 * Checks that none of the locants contain commas
	 * @param locants
	 * @return
	 */
	private boolean locantsAreSingular(List<Element> locants) {
		for (Element locant : locants) {
			if (MATCH_COMMA.split(locant.getValue()).length > 1){
				return false;
			}
		}
		return true;
	}

	/**
	 * Assigns locants and multipliers to substituents/brackets
	 * If both locants and multipliers are present a final check is done that the number of them agree.
	 * WordLevel multipliers are processed e.g. diethyl ethanoate
	 * Adding a locant to a root or any other group that cannot engage in substitive nomenclature will result in an exception being thrown
	 * An exception is made for cases where the locant could be referring to a position on another word
	 * @param subOrBracket
	 * @throws ComponentGenerationException
	 * @throws StructureBuildingException 
	 */
	private void assignLocantsAndMultipliers(Element subOrBracket) throws ComponentGenerationException, StructureBuildingException {
		List<Element> locants = subOrBracket.getChildElements(LOCANT_EL);
		int multiplier =1;
		List<Element> multipliers =  subOrBracket.getChildElements(MULTIPLIER_EL);
		Element parentElem = subOrBracket.getParent();
		boolean oneBelowWordLevel = parentElem.getName().equals(WORD_EL);
		Element groupIfPresent = subOrBracket.getFirstChildElement(GROUP_EL);
		if (multipliers.size()>0){
			if (multipliers.size()>1){
				throw new ComponentGenerationException(subOrBracket.getName() +" has multiple multipliers, unable to determine meaning!");
			}
			if (oneBelowWordLevel &&
					OpsinTools.getNextSibling(subOrBracket) == null &&
					OpsinTools.getPreviousSibling(subOrBracket) == null) {
				return;//word level multiplier
			}
			multiplier = Integer.parseInt(multipliers.get(0).getAttributeValue(VALUE_ATR));
			subOrBracket.addAttribute(new Attribute(MULTIPLIER_ATR, multipliers.get(0).getAttributeValue(VALUE_ATR)));
			//multiplier is INTENTIONALLY not detached. As brackets/subs are only multiplied later on it is neccesary at that stage to determine what elements (if any) are in front of the multiplier
			if (groupIfPresent !=null && PERHALOGENO_SUBTYPE_VAL.equals(groupIfPresent.getAttributeValue(SUBTYPE_ATR))){
				throw new StructureBuildingException(groupIfPresent.getValue() +" cannot be multiplied");
			}
		}
		if(locants.size() > 0) {
			if (multiplier==1 && oneBelowWordLevel && OpsinTools.getPreviousSibling(subOrBracket)==null){//locant might be word Level locant
				if (wordLevelLocantsAllowed(subOrBracket, locants.size())){//something like S-ethyl or S-(2-ethylphenyl) or S-4-tert-butylphenyl
					Element locant = locants.remove(0);
					if (MATCH_COMMA.split(locant.getValue()).length!=1){
						throw new ComponentGenerationException("Multiplier and locant count failed to agree; All locants could not be assigned!");
					}
					parentElem.addAttribute(new Attribute(LOCANT_ATR, locant.getValue()));
					locant.detach();
					if (locants.size()==0){
						return;
					}
				}
			}
			if (subOrBracket.getName().equals(ROOT_EL)){
				locantsToDebugString(locants);
				throw new ComponentGenerationException(locantsToDebugString(locants));
			}
			if (locants.size()!=1){
				throw new ComponentGenerationException(locantsToDebugString(locants));
			}
			Element locantEl = locants.get(0);
			String[] locantValues = MATCH_COMMA.split(locantEl.getValue());
			if (multiplier != locantValues.length){
				throw new ComponentGenerationException("Multiplier and locant count failed to agree; All locants could not be assigned!");
			}

			Element parent = subOrBracket.getParent();
			//attempt to find cases where locant will not be utilised. A special case is made for carbonyl derivatives //e.g. 1H-2-benzopyran-1,3,4-trione 4-[N-(4-chlorophenyl)hydrazone]
			if (!parent.getName().equals(WORD_EL) || !parent.getAttributeValue(TYPE_ATR).equals(WordType.full.toString()) || !state.currentWordRule.equals(WordRule.carbonylDerivative)){
				List<Element> children =parent.getChildElements();
				boolean foundSomethingToSubstitute =false;
				for (int i = parent.indexOf(subOrBracket) +1 ; i < children.size(); i++) {
					if (!children.get(i).getName().equals(HYPHEN_EL)){
						foundSomethingToSubstitute = true;
					}
				}
				if (!foundSomethingToSubstitute){
					throw new ComponentGenerationException(locantsToDebugString(locants));
				}
			}
			if (groupIfPresent !=null && PERHALOGENO_SUBTYPE_VAL.equals(groupIfPresent.getAttributeValue(SUBTYPE_ATR))){
				throw new StructureBuildingException(groupIfPresent.getValue() +" cannot be locanted");
			}
			subOrBracket.addAttribute(new Attribute(LOCANT_ATR, locantEl.getValue()));
			locantEl.detach();
		}
	}
	
	private String locantsToDebugString(List<Element> locants) {
		StringBuilder message = new StringBuilder("Unable to assign all locants. ");
		message.append((locants.size() > 1) ? "These locants " : "This locant ");
		message.append((locants.size() > 1) ? "were " : "was ");
		message.append("not assigned: ");
		for(Element locant : locants) {
			message.append(locant.getValue());
			message.append(" ");
		}
		return message.toString();
	}



	private boolean wordLevelLocantsAllowed(Element subOrBracket, int numberOflocants) {
		Element parentElem = subOrBracket.getParent();
		if (WordType.valueOf(parentElem.getAttributeValue(TYPE_ATR))==WordType.substituent
				&& (OpsinTools.getNextSibling(subOrBracket)==null || numberOflocants>=2)){
			if (state.currentWordRule == WordRule.ester || state.currentWordRule == WordRule.functionalClassEster || state.currentWordRule == WordRule.multiEster || state.currentWordRule == WordRule.acetal){
				return true;
			}
		}
		if ((state.currentWordRule == WordRule.potentialBiochemicalEster || 
				(state.currentWordRule == WordRule.ester &&  (OpsinTools.getNextSibling(subOrBracket)==null || numberOflocants>=2)))
				&& parentElem.getName().equals(WORD_EL)){
			Element wordRule = parentElem.getParent();
			List<Element> words = wordRule.getChildElements(WORD_EL);
			Element ateWord = words.get(words.size()-1);
			if (parentElem==ateWord){
				return true;
			}
		}
			
		return false;
	}

	/**
	 * If a word level multiplier is present e.g. diethyl butandioate then this is processed to ethyl ethyl butandioate
	 * If wordCount is 1 then an exception is thrown if a multiplier is encountered e.g. triphosgene parsed as tri phosgene
	 * @param word
	 * @param wordCount 
	 * @throws StructureBuildingException
	 * @throws ComponentGenerationException
	 */
	private void processWordLevelMultiplierIfApplicable(Element word, int wordCount) throws StructureBuildingException, ComponentGenerationException {
		if (word.getChildCount()==1){
			Element subOrBracket = word.getChild(0);
			Element multiplier = subOrBracket.getFirstChildElement(MULTIPLIER_EL);
			if (multiplier !=null){
				int multiVal =Integer.parseInt(multiplier.getAttributeValue(VALUE_ATR));
				List<Element> locants =subOrBracket.getChildElements(LOCANT_EL);
				boolean assignLocants =false;
				boolean wordLevelLocants = wordLevelLocantsAllowed(subOrBracket, locants.size());//something like O,S-dimethyl phosphorothioate
				if (locants.size()>1){
					throw new ComponentGenerationException("Unable to assign all locants");
				}
				String[] locantValues = null;
				if (locants.size()==1){
					locantValues = MATCH_COMMA.split(locants.get(0).getValue());
					if (locantValues.length == multiVal){
						assignLocants=true;
						locants.get(0).detach();
						if (wordLevelLocants){
							word.addAttribute(new Attribute(LOCANT_ATR, locantValues[0]));
						}
						else{
							throw new ComponentGenerationException(locantsToDebugString(locants));
						}
					}
					else{
						throw new ComponentGenerationException("Unable to assign all locants");
					}
				}
				checkForNonConfusedWithNona(multiplier);
				if (wordCount ==1){
					if (!isMonoFollowedByElement(multiplier, multiVal)){
						throw new StructureBuildingException("Unexpected multiplier found at start of word. Perhaps the name is trivial e.g. triphosgene");
					}
				}
				if (multiVal ==1){//mono
					return;
				}
				List<Element> elementsNotToBeMultiplied = new ArrayList<Element>();//anything before the multiplier
				for (int i = subOrBracket.indexOf(multiplier) -1 ; i >=0 ; i--) {
					Element el = subOrBracket.getChild(i);
					el.detach();
					elementsNotToBeMultiplied.add(el);
				}
				multiplier.detach();
				for(int i=multiVal -1; i>=1; i--) {
					Element clone = state.fragManager.cloneElement(state, word);
					if (assignLocants){
						clone.getAttribute(LOCANT_ATR).setValue(locantValues[i]);
					}
					OpsinTools.insertAfter(word, clone);
				}
				for (Element el : elementsNotToBeMultiplied) {//re-add anything before multiplier to original word
					subOrBracket.insertChild(el, 0);
				}
			}
		}
	}

	private void checkForNonConfusedWithNona(Element multiplier) throws StructureBuildingException {
		if (multiplier.getValue().equals("non")){
			String subsequentUnsemanticToken = multiplier.getAttributeValue(SUBSEQUENTUNSEMANTICTOKEN_ATR);
			if (subsequentUnsemanticToken !=null && subsequentUnsemanticToken.toLowerCase(Locale.ROOT).startsWith("a")){
				return;
			}
			throw new StructureBuildingException("\"non\" probably means \"not\". If a multiplier of value 9 was intended \"nona\" should be used");
		}
	}

	/**
	 * Names like monooxygen may be used to emphasise that a molecule is not being described
	 * @param multiplier
	 * @param multiVal
	 * @return
	 */
	private boolean isMonoFollowedByElement(Element multiplier, int multiVal) {
		if (multiVal ==1){
			Element possibleElement = OpsinTools.getNextSibling(multiplier);
			if (possibleElement != null && possibleElement.getName().equals(GROUP_EL) &&
					(ELEMENTARYATOM_SUBTYPE_VAL.equals(possibleElement.getAttributeValue(SUBTYPE_ATR)) || possibleElement.getValue().equals("hydrogen"))){
				return true;
			}
		}
		return false;
	}
}
