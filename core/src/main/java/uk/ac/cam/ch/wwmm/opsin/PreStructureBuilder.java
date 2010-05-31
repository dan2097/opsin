package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.Stack;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import uk.ac.cam.ch.wwmm.opsin.ParseWord.WordType;
import uk.ac.cam.ch.wwmm.opsin.WordRules.WordRule;
import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;


import nu.xom.Attribute;
import nu.xom.Document;
import nu.xom.Element;
import nu.xom.Elements;
import nu.xom.Node;

/**Performs structure-aware destructive procedural parsing on parser results.
*
* @author dl387
*
*/

class PreStructureBuilder {

	/**
	 * Sorts infix transformations by the number of acceptable inputs for the transformation.
	 * e.g. thio ends up towards the end of the list as it accepts both -O or =O whilst say imido only accepts =O
	 * @author dl387
	 *
	 */
	private class SortInfixTransformations implements Comparator<String> {
		public int compare(String infixTransformation1, String infixTransformation2) {
			int allowedInputs1 = matchComma.split(infixTransformation1).length;
			int allowedInputs2 = matchComma.split(infixTransformation2).length;
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
	
	private final FusedRingBuilder fusedRingBuilder;

	private final Pattern matchCompoundLocant =Pattern.compile("[\\[\\(\\{](\\d+[a-z]?'*)[\\]\\)\\}]");
	private final Pattern matchIndicatedHydrogen =Pattern.compile("(\\d+[a-z]?'*)H");
	private final Pattern matchBracketedEntryInLocant =Pattern.compile("[\\[\\(\\{].*[\\]\\)\\}]");
	private final Pattern matchCisTransInLocants =Pattern.compile("[rct]-");
	private final Pattern matchColon =Pattern.compile(":");
	private final Pattern matchSemiColon =Pattern.compile(";");
	private final Pattern matchComma =Pattern.compile(",");
	private final Pattern matchSpace =Pattern.compile(" ");
	private final Pattern matchElementSymbolOrAminoAcidLocant = Pattern.compile("[A-Z][a-z]?'*(\\d+[a-z]?'*)?");
	private final Pattern matchElementSymbol = Pattern.compile("[A-Z][a-z]?");
	private final Pattern matchOrtho =Pattern.compile("[oO]");
	private final Pattern matchMeta =Pattern.compile("[mM]");
	private final Pattern matchPara =Pattern.compile("[pP]");
	private final Pattern matchNumericLocant =Pattern.compile("\\d+[a-z]?'*");
	private final Pattern matchChalcogen = Pattern.compile("O|S|Se|Te");
	private final Pattern matchChalogenReplacment= Pattern.compile("thio|seleno|telluro");
	private final Pattern matchInlineSuffixesThatAreAlsoGroups = Pattern.compile("carbon|oxy|sulfen|sulfin|sulfon|selenen|selenin|selenon|telluren|tellurin|telluron");

	/*Holds the rules on how suffixes are interpreted. Convenience methods are available to use them*/
	private HashMap<String, HashMap<String, List<Element>>> suffixApplicability;
	private HashMap<String, Element> suffixRules;
	
	//rings that look like HW rings but have other meanings. For the HW like inorganics the true meaning is given
	private static HashMap<String, String[]> specialHWRings = new HashMap<String, String[]>();
	static{
		//The first entry of the array is a special instruction e.g. blocked or saturated. The correct order of the heteroatoms follows
		//terminal e is ignored from all of the keys as it is optional in the input name
		specialHWRings.put("oxin", new String[]{"blocked"});
		specialHWRings.put("azin", new String[]{"blocked"});

		specialHWRings.put("oxazol", new String[]{"","O","C","N","C","C"});
		specialHWRings.put("thiazol", new String[]{"","S","C","N","C","C"});
		specialHWRings.put("selenazol", new String[]{"","Se","C","N","C","C"});
		specialHWRings.put("tellurazol", new String[]{"","Te","C","N","C","C"});
		specialHWRings.put("oxazolidin", new String[]{"","O","C","N","C","C"});
		specialHWRings.put("thiazolidin", new String[]{"","S","C","N","C","C"});
		specialHWRings.put("selenazolidin", new String[]{"","Se","C","N","C","C"});
		specialHWRings.put("tellurazolidin", new String[]{"","Te","C","N","C","C"});
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

	PreStructureBuilder(FusedRingBuilder fusedRingBuilder, ResourceGetter resourceGetter) throws Exception {
		this.fusedRingBuilder = fusedRingBuilder;
		//Populate suffix rules/applicability hashes
		Document suffixApplicabilityDoc = resourceGetter.getXMLDocument("suffixApplicability.xml");
		Document suffixRulesDoc = resourceGetter.getXMLDocument("suffixRules.xml");
		suffixApplicability = new HashMap<String, HashMap<String,List<Element>>>();
		suffixRules = new HashMap<String, Element>();
		Elements groupTypes = suffixApplicabilityDoc.getRootElement().getChildElements("groupType");
		for (int i = 0; i < groupTypes.size(); i++) {
			Element groupType =groupTypes.get(i);
			Elements suffixes = groupType.getChildElements("suffix");
			HashMap<String, List<Element>> suffixToRuleMap= new HashMap<String, List<Element>>();
			for (int j = 0; j < suffixes.size(); j++) {
				Element suffix =suffixes.get(j);
				String suffixValue= suffix.getAttributeValue(VALUE_ATR);
				if (suffixToRuleMap.get(suffixValue)!=null){//can have multiple entries if subType attribute is set
					suffixToRuleMap.get(suffixValue).add(suffix);
				}
				else{
					ArrayList<Element> suffixList =new ArrayList<Element>();
					suffixList.add(suffix);
					suffixToRuleMap.put(suffixValue, suffixList);
				}
			}
			suffixApplicability.put(groupType.getAttributeValue(TYPE_ATR), suffixToRuleMap);
		}

		Elements rules = suffixRulesDoc.getRootElement().getChildElements("rule");
		for (int i = 0; i < rules.size(); i++) {
			Element rule =rules.get(i);
			String ruleValue=rule.getAttributeValue(VALUE_ATR);
			if (suffixRules.get(ruleValue)!=null){
				throw new Exception("Suffix: " +ruleValue +" appears multiple times in suffixRules.xml");
			}
			suffixRules.put(ruleValue, rule);
		}
	}


	/** The master method, postprocesses a parse result. At this stage one can except all substituents/roots to have at least 1 group.
	 * Multiple groups are present in, for example, fusion nomenclature. By the end of this function there will be exactly 1 group
	 * associated with each substituent/root. Multiplicative nomenclature can result in there being multiple roots
	 *
	 * @param elem The element to postprocess.
	 * @return
	 * @throws PostProcessingException 
	 * @throws StructureBuildingException 
	 */
	void postProcess(BuildState state, Element elem) throws PostProcessingException, StructureBuildingException {
		List<Element> words =XOMTools.getDescendantElementsWithTagName(elem, WORD_EL);
		int wordCount =words.size();
		for (Element word : words) {
			String wordRule = OpsinTools.getParentWordRule(word).getAttributeValue(WORDRULE_EL);
			state.currentWordRule = WordRule.valueOf(wordRule);
			if (word.getAttributeValue(TYPE_ATR).equals(WordType.functionalTerm.toString())){
				continue;//functionalTerms are handled on a case by case basis by wordRules
			}

			List<Element> roots = XOMTools.getDescendantElementsWithTagName(word, ROOT_EL);
			if (roots.size() >1){
				throw new PostProcessingException("Multiple roots, but only 0 or 1 were expected. Found: " +roots.size());
			}
			List<Element> substituents = XOMTools.getDescendantElementsWithTagName(word, SUBSTITUENT_EL);
			List<Element> substituentsAndRoot = OpsinTools.combineElementLists(substituents, roots);
			List<Element> brackets =  XOMTools.getDescendantElementsWithTagName(word, BRACKET_EL);
			List<Element> substituentsAndRootAndBrackets =OpsinTools.combineElementLists(substituentsAndRoot, brackets);
			List<Element> groups =  XOMTools.getDescendantElementsWithTagName(word, GROUP_EL);

			for (Element subOrBracketOrRoot : substituentsAndRootAndBrackets) {
				processLocantFeatures(subOrBracketOrRoot);
			}

			for (Element group : groups) {
				Fragment thisFrag = resolveGroup(state, group);
				state.xmlFragmentMap.put(group, thisFrag);
			}

			Element finalSubOrRootInWord =(Element) word.getChild(word.getChildElements().size()-1);
			while (!finalSubOrRootInWord.getLocalName().equals(ROOT_EL) && !finalSubOrRootInWord.getLocalName().equals(SUBSTITUENT_EL)){
				List<Element> children = XOMTools.getChildElementsWithTagNames(finalSubOrRootInWord, new String[]{ROOT_EL, SUBSTITUENT_EL, BRACKET_EL});
				if (children.size()==0){
					throw new PostProcessingException("Unable to find finalSubOrRootInWord");
				}
				finalSubOrRootInWord = children.get(children.size()-1);
			}
		
			for (Element subOrRootOrBracket : substituentsAndRootAndBrackets) {
				checkAndConvertToSingleLocants(state, subOrRootOrBracket, finalSubOrRootInWord);
			}

			for (Element subOrRoot : substituentsAndRoot) {
				processMultipliers(subOrRoot);
				detectConjunctiveSuffixGroups(state, subOrRoot, groups);
				matchLocantsToDirectFeatures(subOrRoot);

				Elements groupsOfSubOrRoot = subOrRoot.getChildElements(GROUP_EL);
				Element lastGroupInSubOrRoot =groupsOfSubOrRoot.get(groupsOfSubOrRoot.size()-1);
				preliminaryProcessSuffixes(state, lastGroupInSubOrRoot, XOMTools.getChildElementsWithTagName(subOrRoot, SUFFIX_EL));
			}

			if (processPrefixFunctionalReplacementNomenclature(state, groups, substituents)){//true if functional replacement performed, 1 or more substituents will have been removed
				substituentsAndRoot = OpsinTools.combineElementLists(substituents, roots);
				substituentsAndRootAndBrackets =OpsinTools.combineElementLists(substituentsAndRoot, brackets);
			}

			for (Element subOrRoot : substituentsAndRoot) {
				processHW(state, subOrRoot);//hantzch-widman rings
				fusedRingBuilder.processFusedRings(state, subOrRoot);
				assignElementSymbolLocants(state, subOrRoot);
				processRingAssemblies(state, subOrRoot);
				processPolyCyclicSpiroNomenclature(state, subOrRoot);
			}

			for (Element subOrRoot : substituentsAndRoot) {
				applyLambdaConvention(state, subOrRoot);
				handleMultiRadicals(state, subOrRoot);
			}

			//System.out.println(new XOMFormatter().elemToString(elem));
			addImplicitBracketsToAminoAcids(state, groups, brackets);
			findAndStructureImplictBrackets(state, substituents, brackets);

			substituentsAndRootAndBrackets =OpsinTools.combineElementLists(substituentsAndRoot, brackets);//findAndStructureImplictBrackets may have created new brackets

			for (Element subOrRoot : substituentsAndRoot) {
				matchLocantsToIndirectFeatures(state, subOrRoot);
				processConjunctiveNomenclature(state, subOrRoot);
				resolveSuffixes(state, subOrRoot.getFirstChildElement(GROUP_EL), XOMTools.getChildElementsWithTagName(subOrRoot, SUFFIX_EL));
			}

			removeClarifyingBrackets(brackets, substituentsAndRootAndBrackets);//e.g. (tetramethyl)azanium == tetramethylazanium

			if (word.getChildCount()>1){
				assignLocantsToMultipliedRootIfPresent(state, (Element) word.getChild(word.getChildCount()-1));//multiplicative nomenclature e.g. methylenedibenzene or 3,4'-oxydipyridine
			}

			for (Element subBracketOrRoot : substituentsAndRootAndBrackets) {
				assignLocantsAndMultipliers(state, subBracketOrRoot);
			}
			processWordLevelMultiplierIfApplicable(state, word, wordCount);

		}
	}

	/**Handles special features of locants e.g. ortho/meta/para, indicated hydrogen, cis/trans in locant
	 *
	 * @param elem The substituent/root/bracket to looks for locants in.
	 * @throws PostProcessingException
	 */
	private void processLocantFeatures(Element elem) throws PostProcessingException{
		Elements ompLocants = elem.getChildElements("orthoMetaPara");
		for(int i=0;i<ompLocants.size();i++) {
			Element locant = ompLocants.get(i);
			String locantText = locant.getValue();
			locantText = locantText.substring(0, 1);
			Element afterOmpLocant = (Element)XOMTools.getNextSibling(locant);
			locant.setLocalName("locant");
			locant.removeChildren();
			locant.addAttribute(new Attribute("type", "orthoMetaPara"));
			if(afterOmpLocant.getLocalName().equals("multiplier") || (afterOmpLocant.getAttribute("outIDs")!=null && matchComma.split(afterOmpLocant.getAttributeValue("outIDs")).length>1) ) {
				if (matchOrtho.matcher(locantText).matches()){
					locant.appendChild("1,2");
				}
				else if (matchMeta.matcher(locantText).matches()){
					locant.appendChild("1,3");
				}
				else if (matchPara.matcher(locantText).matches()){
					locant.appendChild("1,4");
				}
				else{
					throw new PostProcessingException(locantText + " was not identified as being either ortho, meta or para but according to the chemical grammar it should of been");
				}
			}
			else{
				if (matchOrtho.matcher(locantText).matches()){
					locant.appendChild("2");
				}
				else if (matchMeta.matcher(locantText).matches()){
					locant.appendChild("3");
				}
				else if (matchPara.matcher(locantText).matches()){
					locant.appendChild("4");
				}
				else{
					throw new PostProcessingException(locantText + " was not identified as being either ortho, meta or para but according to the chemical grammar it should of been");
				}
			}
		}

		Elements locants = elem.getChildElements("locant");
		for(int i=0;i<locants.size();i++) {
			Element locant = locants.get(i);
			String locantText = locant.getValue();

			//If the indicatedHydrogen has been specified create a tag for it and remove it from the list of locants
			//e.g. 1(9H),5,7 -->indicatedHydrogen tag value (9H) and 1,5,7
			//can get as complicated as 1,2(2H,7H)
			Matcher matches =matchIndicatedHydrogen.matcher(locantText);
			if (matches.find()){
				do {
					Element indicatedHydrogenElement=new Element("indicatedHydrogen");
					indicatedHydrogenElement.addAttribute(new Attribute(LOCANT_ATR, matches.group(1)));
					XOMTools.insertBefore(locant, indicatedHydrogenElement);
				}
				while (matches.find());
				locantText =matchBracketedEntryInLocant.matcher(locantText).replaceAll("");
			}

			/*
			 * Strip out cis/trans information built into locant - currently unhandled
			 */
			matches =matchCisTransInLocants.matcher(locantText);
			if (matches.find()){
				do {
					//currently do nothing
				}
				while (matches.find());
				locantText =matches.replaceAll("");
			}
			XOMTools.setTextChild(locant, locantText);

			Element afterLocants = (Element)XOMTools.getNextSibling(locant);
			if(afterLocants == null) throw new PostProcessingException("Nothing after locant tag: " + locant.toXML());
		}
	}


	/**Resolves the contents of a &lt;group&gt; tag.
	 *
	 * @param group The &lt;group&gt; tag.
	 * @return The fragment specified by the tag.
	 * @throws StructureBuildingException If the group can't be built.
	 * @throws PostProcessingException
	 */
	private Fragment resolveGroup(BuildState state, Element group) throws StructureBuildingException, PostProcessingException {
		String groupType = group.getAttributeValue(TYPE_ATR);
		String groupSubType = group.getAttributeValue("subType");
		String groupValue = group.getAttributeValue(VALUE_ATR);
		String groupValType = group.getAttributeValue("valType");
		Fragment thisFrag =null;
		if(groupValType.equals("chain")) {
			int alkaneLength = new Integer(groupValue);
			String smiles = StringTools.multiplyString("C", alkaneLength);
			thisFrag = state.fragManager.buildSMILES(smiles, groupType, groupSubType, "");
		} else if(groupValType.equals("ring") || groupValType.equals("partunsatring")) {
			int alkaneLength = new Integer(groupValue);
			String smiles = "C1";
			smiles += StringTools.multiplyString("C", alkaneLength-1);
			smiles += "1";
			thisFrag = state.fragManager.buildSMILES(smiles, groupType, groupSubType, "");
		} else if(groupValType.equals("unsatring")) {
			int alkaneLength = new Integer(groupValue);
			String smiles = "c1";
			smiles += StringTools.multiplyString("c", alkaneLength-1);
			smiles += "1";
			thisFrag = state.fragManager.buildSMILES(smiles, groupType, groupSubType, "");
		} else if(groupValType.equals("SMILES")) {
			if (group.getAttribute("labels")!=null){
				thisFrag = state.fragManager.buildSMILES(groupValue, groupType, groupSubType, group.getAttributeValue("labels"));
			}
			else{
				thisFrag = state.fragManager.buildSMILES(groupValue, groupType, groupSubType, "");
			}
		} else if(groupValType.equals("dbkey")) {
			thisFrag = state.fragManager.buildCML(groupValue, groupType, groupSubType);
		} else if(groupValType.equals("atom")) {
			thisFrag = state.fragManager.buildSMILES("[" + groupValue + "]", groupType, groupSubType, "none");
		}
		else{
			throw new StructureBuildingException("Group tag has bad or missing valType: " + group.toXML());
		}
		if (thisFrag ==null){
			throw new StructureBuildingException("null fragment returned from the following xml: " + group.toXML());
		}

		//processes groups like cymene and xylene whose structure is determined by the presence of a locant in front e.g. p-xylene
		processXyleneLikeNomenclature(state, group, thisFrag);

		FragmentTools.convertHighOrderBondsToSpareValencies(thisFrag);//only applied to cyclic bonds

		if (group.getAttribute("defaultInLocant")!=null){//sets the atom at which substitution will occur to by default
			thisFrag.setDefaultInAtom(thisFrag.getAtomByLocantOrThrow(group.getAttributeValue("defaultInLocant")));
		}
		else if (group.getAttribute("defaultInID")!=null){
			thisFrag.setDefaultInAtom(thisFrag.getAtomByIDOrThrow(thisFrag.getIdOfFirstAtom() + Integer.parseInt(group.getAttributeValue("defaultInID")) -1));
		}
		else if (group.getAttribute("usableAsAJoiner") != null && group.getAttributeValue("usableAsAJoiner").equals("yes") && group.getAttribute("suffixAppliesTo")==null){//makes linkers by default attach end to end
			int chainLength =thisFrag.getChainLength();
			if (chainLength >1){
				boolean connectEndToEndWithPreviousSub =true;
				if (groupSubType.equals("alkaneStem")){//don't do this if you the group is preceded by another alkaneStem e.g. methylethyl makes more sense as prop-2-yl rather than propyl
					Element previousSubstituent =(Element) XOMTools.getPreviousSibling(group.getParent());
					if (previousSubstituent!=null){
						Elements groups = previousSubstituent.getChildElements("group");
						if (groups.size()==1 && groups.get(0).getAttributeValue("subType").equals("alkaneStem") && !groups.get(0).getAttributeValue(TYPE_ATR).equals("ring")){
							connectEndToEndWithPreviousSub = false;
						}
					}
				}
				if (connectEndToEndWithPreviousSub){
					group.addAttribute(new Attribute("defaultInID",Integer.toString(chainLength)));
					thisFrag.setDefaultInAtom(thisFrag.getAtomByLocantOrThrow(Integer.toString(chainLength)));
				}
			}
		}

		if (group.getAttribute("functionalIDs")!=null){
			String[] functionalIDs = matchComma.split(group.getAttributeValue("functionalIDs"));
            for (String functionalID : functionalIDs) {
                thisFrag.addFunctionalAtom(thisFrag.getAtomByIDOrThrow(thisFrag.getIdOfFirstAtom() + Integer.parseInt(functionalID) - 1));
            }
		}
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
	 * @throws PostProcessingException
	 */
	private void processXyleneLikeNomenclature(BuildState state, Element group, Fragment parentFrag) throws StructureBuildingException, PostProcessingException {
		if(group.getAttributeValue("addGroup")!=null) {
			String addGroupInformation=group.getAttributeValue("addGroup");
			String[] groupsToBeAdded = matchSemiColon.split(addGroupInformation);//typically only one, but 2 in the case of xylene and quinones
			ArrayList<HashMap<String, String>> allGroupInformation = new ArrayList<HashMap<String, String>>();
            for (String groupToBeAdded : groupsToBeAdded) {//populate allGroupInformation list
                String[] tempArray = matchSpace.split(groupToBeAdded);
                HashMap<String, String> groupInformation = new HashMap<String, String>();
                if (tempArray.length != 2 && tempArray.length != 3) {
                    throw new PostProcessingException("malformed addGroup tag");
                }
                groupInformation.put("SMILES", tempArray[0]);
                if (tempArray[1].startsWith("id")) {
                    groupInformation.put("atomReferenceType", "id");
                    groupInformation.put("atomReference", tempArray[1].substring(2));
                } else if (tempArray[1].startsWith("locant")) {
                    groupInformation.put("atomReferenceType", "locant");
                    groupInformation.put("atomReference", tempArray[1].substring(6));
                } else {
                    throw new PostProcessingException("malformed addGroup tag");
                }
                if (tempArray.length == 3) {//labels may optionally be specified for the group to be added
                    groupInformation.put("labels", tempArray[2]);
                }
                allGroupInformation.add(groupInformation);
            }
			Element previousEl =(Element) XOMTools.getPreviousSibling(group);
			if (previousEl !=null && previousEl.getLocalName().equals("locant")){//has the name got specified locants to override the default ones
				List<String> locantValues =StringTools.arrayToList(matchComma.split(previousEl.getValue()));
				boolean assignlocants =true;
				if (locantValues.size()<groupsToBeAdded.length){
					if (locantValues.size() +1 <groupsToBeAdded.length ){//only one locant can be implicit
						assignlocants=false;
					}
					else {//check that the firstGroup by default will be added to the atom with locant 1. If this is not the case then as many locants as there were groups should of been specified
						//or no locants should have been specified, which is what will be assumed (i.e. the locants will be left unassigned)
						HashMap<String, String> groupInformation =allGroupInformation.get(0);
						String locant;
						if (groupInformation.get("atomReferenceType").equals("locant")){
							locant =parentFrag.getAtomByLocantOrThrow(groupInformation.get("atomReference")).getFirstLocant();
						}
						else if (groupInformation.get("atomReferenceType").equals("id") ){
							locant =parentFrag.getAtomByIDOrThrow(parentFrag.getIdOfFirstAtom() + Integer.parseInt(groupInformation.get("atomReference")) -1 ).getFirstLocant();
						}
						else{
							throw new PostProcessingException("malformed addGroup tag");
						}
						if (locant ==null || !locant.equals("1")){
							assignlocants=false;
						}
					}

				}
				if (assignlocants){
					for (int i = groupsToBeAdded.length -1; i >=0 ; i--) {
						//if less locants than expected are specified the locants of only the later groups will be changed
						//e.g. 4-xylene will transform 1,2-xylene to 1,4-xylene
						HashMap<String, String> groupInformation =allGroupInformation.get(i);
						if (locantValues.size() >0){
							groupInformation.put("atomReferenceType", "locant");
							groupInformation.put("atomReference", locantValues.get(locantValues.size()-1));
							locantValues.remove(locantValues.size()-1);
						}
						else{
							break;
						}
					}
					if (locantValues.size() ==0){
						previousEl.detach();
					}
					else{
						XOMTools.setTextChild(previousEl, StringTools.stringListToString(locantValues, ","));
					}
				}
			}

			for (int i = 0; i < groupsToBeAdded.length; i++) {
				HashMap<String, String> groupInformation =allGroupInformation.get(i);
				String smilesOfGroupToBeAdded = groupInformation.get("SMILES");
				Fragment newFrag;
				if (groupInformation.get("labels")!=null){
					newFrag = state.fragManager.buildSMILES(smilesOfGroupToBeAdded, parentFrag.getType(), parentFrag.getSubType(), groupInformation.get("labels"));
				}
				else{
					newFrag = state.fragManager.buildSMILES(smilesOfGroupToBeAdded, parentFrag.getType(), parentFrag.getSubType(), "none");
				}

				Atom atomOnParentFrag =null;
				if (groupInformation.get("atomReferenceType").equals("locant")){
					atomOnParentFrag=parentFrag.getAtomByLocantOrThrow(groupInformation.get("atomReference"));
				}
				else if (groupInformation.get("atomReferenceType").equals("id") ){
					atomOnParentFrag= parentFrag.getAtomByIDOrThrow(parentFrag.getIdOfFirstAtom() + Integer.parseInt(groupInformation.get("atomReference")) -1);
				}
				else{
					throw new PostProcessingException("malformed addGroup tag");
				}
				if (newFrag.getOutAtoms().size() >1){
					throw new PostProcessingException("too many outAtoms on group to be added");
				}
				if (newFrag.getOutAtoms().size() ==1) {
					OutAtom newFragOutAtom = newFrag.getOutAtom(0);
					newFrag.removeOutAtom(newFragOutAtom);
					state.fragManager.incorporateFragment(newFrag, newFragOutAtom.getAtom().getID(), parentFrag, atomOnParentFrag.getID(), newFragOutAtom.getValency());
				}
				else{
					Atom atomOnNewFrag = newFrag.getDefaultInAtom();
					state.fragManager.incorporateFragment(newFrag, atomOnNewFrag.getID(), parentFrag, atomOnParentFrag.getID(), 1);
				}
			}
		}

		if(group.getAttributeValue("addHeteroAtom")!=null) {
			String addHeteroAtomInformation=group.getAttributeValue("addHeteroAtom");
			String[] heteroAtomsToBeAdded = matchSemiColon.split(addHeteroAtomInformation);
			ArrayList<HashMap<String, String>> allHeteroAtomInformation = new ArrayList<HashMap<String, String>>();
            for (String heteroAtomToBeAdded : heteroAtomsToBeAdded) {//populate allHeteroAtomInformation list
                String[] tempArray = matchSpace.split(heteroAtomToBeAdded);
                HashMap<String, String> heteroAtomInformation = new HashMap<String, String>();
                if (tempArray.length != 2) {
                    throw new PostProcessingException("malformed addHeteroAtom tag");
                }
                heteroAtomInformation.put("SMILES", tempArray[0]);
                if (tempArray[1].startsWith("id")) {
                    heteroAtomInformation.put("atomReferenceType", "id");
                    heteroAtomInformation.put("atomReference", tempArray[1].substring(2));
                } else if (tempArray[1].startsWith("locant")) {
                    heteroAtomInformation.put("atomReferenceType", "locant");
                    heteroAtomInformation.put("atomReference", tempArray[1].substring(6));
                } else {
                    throw new PostProcessingException("malformed addHeteroAtom tag");
                }
                allHeteroAtomInformation.add(heteroAtomInformation);
            }
			Element previousEl =(Element) XOMTools.getPreviousSibling(group);
			if (previousEl !=null && previousEl.getLocalName().equals("locant")){//has the name got specified locants to override the default ones
				List<String> locantValues =StringTools.arrayToList(matchComma.split(previousEl.getValue()));
				if (locantValues.size() >=heteroAtomsToBeAdded.length){
					for (int i = heteroAtomsToBeAdded.length -1; i >=0 ; i--) {//all heteroatoms must have a locant or default locants will be used
						HashMap<String, String> groupInformation =allHeteroAtomInformation.get(i);
						groupInformation.put("atomReferenceType", "locant");
						groupInformation.put("atomReference", locantValues.get(locantValues.size()-1));
						locantValues.remove(locantValues.size()-1);
					}
					if (locantValues.size() ==0){
						previousEl.detach();
					}
					else{
						XOMTools.setTextChild(previousEl, StringTools.stringListToString(locantValues, ","));
					}
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
					throw new PostProcessingException("malformed addHeteroAtom tag");
				}
				String atomSymbol = heteroAtomInformation.get("SMILES");
				if(atomSymbol.startsWith("[")) {
					Fragment f = state.fragManager.buildSMILES(atomSymbol);
					Atom referenceAtom = f.getAtomList().get(0);
					atomSymbol =referenceAtom.getElement();
					atomOnParentFrag.setCharge(referenceAtom.getCharge());
					state.fragManager.removeFragment(f);
				}
				atomOnParentFrag.setElement(atomSymbol);
			}
		}

		if(group.getAttributeValue("addBond")!=null) {
			String addBondInformation=group.getAttributeValue("addBond");
			String[] bondsToBeAdded = matchSemiColon.split(addBondInformation);
			ArrayList<HashMap<String, String>> allBondInformation = new ArrayList<HashMap<String, String>>();
            for (String bondToBeAdded : bondsToBeAdded) {//populate allBondInformation list
                String[] tempArray = matchSpace.split(bondToBeAdded);
                HashMap<String, String> bondInformation = new HashMap<String, String>();
                if (tempArray.length != 2) {
                    throw new PostProcessingException("malformed addBond tag");
                }
                bondInformation.put("bondOrder", tempArray[0]);
                if (tempArray[1].startsWith("id")) {
                    bondInformation.put("atomReferenceType", "id");
                    bondInformation.put("atomReference", tempArray[1].substring(2));
                } else if (tempArray[1].startsWith("locant")) {
                    bondInformation.put("atomReferenceType", "locant");
                    bondInformation.put("atomReference", tempArray[1].substring(6));
                } else {
                    throw new PostProcessingException("malformed addBond tag");
                }
                allBondInformation.add(bondInformation);
            }
			Element previousEl =(Element) XOMTools.getPreviousSibling(group);
			if (previousEl !=null && previousEl.getLocalName().equals("locant")){//has the name got specified locants to override the default ones
				List<String> locantValues =StringTools.arrayToList(matchComma.split(previousEl.getValue()));
				if (locantValues.size() >=bondsToBeAdded.length){
					for (int i = bondsToBeAdded.length -1; i >=0 ; i--) {//all heteroatoms must have a locant or default locants will be used
						HashMap<String, String> bondInformation =allBondInformation.get(i);
						bondInformation.put("atomReferenceType", "locant");
						bondInformation.put("atomReference", locantValues.get(locantValues.size()-1));
						locantValues.remove(locantValues.size()-1);
					}
					if (locantValues.size() ==0){
						previousEl.detach();
					}
					else{
						XOMTools.setTextChild(previousEl, StringTools.stringListToString(locantValues, ","));
					}
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
					throw new PostProcessingException("malformed addBond tag");
				}
				FragmentTools.unsaturate(atomOnParentFrag.getID(), Integer.parseInt(bondInformation.get("bondOrder")) , parentFrag);
			}
		}
	}


	/**Converts locantgroups to individual locant tags, and checks for agreement
	 * between the number of locants, and multipliers. If there is disagreement the locant is checked against some special cases.
	 * If this fails an exception is thrown.
	 * @param state
	 * @param subOrBracketOrRoot The substituent/root/bracket to looks for locants in.
	 * @param finalSubOrRootInWord : used to check if a locant is referring to the root as in multiplicative nomenclature
	 * @throws PostProcessingException
	 * @throws StructureBuildingException
	 * @throws PostProcessingException If there is a disagreement.
	 */
	private void checkAndConvertToSingleLocants(BuildState state, Element subOrBracketOrRoot, Element finalSubOrRootInWord) throws StructureBuildingException, PostProcessingException {
		Elements locants = subOrBracketOrRoot.getChildElements("locant");
		Element group =subOrBracketOrRoot.getFirstChildElement(GROUP_EL);//will be null if element is a bracket
		for(int i=0;i<locants.size();i++) {
			Element locant = locants.get(i);
			String [] locantValues = matchComma.split(locant.getValue());

			Element afterLocants = (Element)XOMTools.getNextSibling(locant);
			if(locantValues.length > 1) {
				while (afterLocants !=null){
					if(afterLocants.getLocalName().equals("multiplier") || afterLocants.getLocalName().equals("locant")) {
						break;
					}
					afterLocants = (Element)XOMTools.getNextSibling(afterLocants);
				}
				if(afterLocants != null && afterLocants.getLocalName().equals("multiplier")) {
					if(Integer.parseInt(afterLocants.getAttributeValue(VALUE_ATR)) == locantValues.length ) {
						// number of locants and multiplier agree
						boolean locantModified =false;//did determineLocantMeaning do something?
						if (locantValues[locantValues.length-1].endsWith("'") && group!=null && subOrBracketOrRoot.indexOf(group) > subOrBracketOrRoot.indexOf(locant)){//quite possible that this is referring to a multiplied root

							if (group.getAttribute("outIDs")!=null && matchComma.split(group.getAttributeValue("outIDs")).length>1){
								locantModified=determineLocantMeaning(state, locant, locantValues, finalSubOrRootInWord);
							}
							else{
								Element afterGroup = (Element)XOMTools.getNextSibling(group);
								int inlineSuffixCount =0;
								int multiplier=1;
								while (afterGroup !=null){
									if(afterGroup.getLocalName().equals("multiplier")){
										multiplier =Integer.parseInt(afterGroup.getAttributeValue(VALUE_ATR));
									}
									else if(afterGroup.getLocalName().equals("suffix") && afterGroup.getAttributeValue(TYPE_ATR).equals("inline")){
										inlineSuffixCount +=(multiplier);
										multiplier=1;
									}
									afterGroup = (Element)XOMTools.getNextSibling(afterGroup);
								}
								if (inlineSuffixCount >=2){
									locantModified=determineLocantMeaning(state, locant, locantValues, finalSubOrRootInWord);
								}
							}
						}
						if (!locantModified && !XOMTools.getNextSibling(locant).equals(afterLocants)){//the locants apply indirectly the multiplier e.g. 2,3-butandiol
							//move the locant to be next to the multiplier.
							locant.detach();
							XOMTools.insertBefore(afterLocants, locant);
						}
					} else {
						if(!determineLocantMeaning(state, locant, locantValues, finalSubOrRootInWord)) throw new PostProcessingException("Mismatch between locant and multiplier counts (" +
								Integer.toString(locantValues.length) + " and " + afterLocants.getAttributeValue(VALUE_ATR) + "):" + locant.toXML());
					}
				} else {
					/* Multiple locants without a multiplier */
					if(!determineLocantMeaning(state, locant, locantValues, finalSubOrRootInWord)) throw new PostProcessingException("Multiple locants without a multiplier: " + locant.toXML());
				}
			}

			//checks that determineLocantMeaning hasn't moved the locant or already assigned the locant a special meaning and detached it
			if (locant.getParent()!=null && locant.getParent().equals(subOrBracketOrRoot)){

                for (String locantValue : locantValues) {
                    Element singleLocant = new Element("locant");
                    String locantType = null;
                    if (locant.getAttribute("type") != null) {
                        locantType = locant.getAttributeValue(TYPE_ATR);
                        singleLocant.addAttribute(new Attribute("type", locantType));
                    }
                    Matcher matches = matchCompoundLocant.matcher(locantValue);
                    if (matches.find()) {
                        singleLocant.addAttribute(new Attribute("compoundLocant", matches.group(1)));
                        locantValue = matches.replaceAll("");
                    }
                    singleLocant.addAttribute(new Attribute("value", locantValue));
                    XOMTools.insertBefore(locant, singleLocant);
                }
				locant.detach();
			}
		}
	}


	/**Looks for Hantzch-Widman systems, and sees if the number of locants
	 * agrees with the number of heteroatoms.
	 * If this is not the case alternative possibilities are tested:
	 * 	The locants could be intended to indicate the position of outAtoms e.g. 1,4-phenylene
	 * 	The locants could be intended to indicate the attachement points of the root groups in multiplicative nomenclature e.g. 4,4'-methylenedibenzoic acid
	 * @param state
	 * @param locant The element corresponding to the locant group to be tested
	 * @param locantValues The locant values;
	 * @param finalSubOrRootInWord : used to check if a locant is referring to the root as in multiplicative nomenclatures)
	 * @return true if there's a HW system, and agreement; or if the locants conform to one of the alternative possibilities, otherwise false.
	 * @throws StructureBuildingException
	 */
	private boolean determineLocantMeaning(BuildState state, Element locant, String[] locantValues, Element finalSubOrRootInWord) throws StructureBuildingException {
		int count =locantValues.length;
		Element currentElem = (Element)XOMTools.getNextSibling(locant);
		int heteroCount = 0;
		int multiplierValue = 1;
		while(currentElem != null && !currentElem.getLocalName().equals(GROUP_EL)){
			if(currentElem.getLocalName().equals(HETEROATOM_EL)) {
				heteroCount+=multiplierValue;
				multiplierValue =1;
			} else if (currentElem.getLocalName().equals(MULTIPLIER_EL)){
				multiplierValue = Integer.parseInt(currentElem.getAttributeValue(VALUE_ATR));
			}
			else{
				break;
			}
			currentElem = (Element)XOMTools.getNextSibling(currentElem);
		}
		if(currentElem != null && currentElem.getLocalName().equals(GROUP_EL)){
			if (currentElem.getAttributeValue(SUBTYPE_ATR).equals(HANTZSCHWIDMAN_SUBTYPE_VAL)) {
				if(heteroCount == count) {
					return true;
				} else {
					return false;//there is a case where locants don't apply to heteroatoms in a HW system, but in that case only one locant is expected so this function would not be called
				}
			}
			else if (heteroCount==0 && currentElem.getAttribute(OUTIDS_ATR)!=null ) {//e.g. 1,4-phenylene
				String[] outIDs = matchComma.split(currentElem.getAttributeValue(OUTIDS_ATR), -1);
				Fragment groupFragment =state.xmlFragmentMap.get(currentElem);
				if (count ==outIDs.length && groupFragment.getAtomList().size()>1){//things like oxy do not need to have their outIDs specified
					int idOfFirstAtomInFrag =groupFragment.getIdOfFirstAtom();
					boolean foundLocantNotPresentOnFragment = false;
					for (int i = outIDs.length-1; i >=0; i--) {
						Atom a =groupFragment.getAtomByLocant(locantValues[i]);
						if (a==null){
							foundLocantNotPresentOnFragment = true;
							break;
						}
						outIDs[i]=Integer.toString(a.getID() -idOfFirstAtomInFrag +1);//convert to relative id
					}
					if (!foundLocantNotPresentOnFragment){
						currentElem.getAttribute(OUTIDS_ATR).setValue(StringTools.arrayToString(outIDs, ","));
						locant.detach();
						return true;
					}
				}
			}
			else if(currentElem.getValue().equals("benz") || currentElem.getValue().equals("benzo")){
				Node potentialGroupAfterBenzo = XOMTools.getNextSibling(currentElem, GROUP_EL);//need to make sure this isn't benzyl
				if (potentialGroupAfterBenzo!=null){
					return true;//e.g. 1,2-benzothiazole
				}
			}
		}
		if(currentElem != null && currentElem.getLocalName().equals(POLYCYCLICSPIRO_EL)){
			return true;
		}
		Element multiplier =(Element) finalSubOrRootInWord.getChild(0);
		if (!multiplier.getLocalName().equals(MULTIPLIER_EL) && ((Element)finalSubOrRootInWord.getParent()).getLocalName().equals(BRACKET_EL)){//e.g. 1,1'-ethynediylbis(1-cyclopentanol)
			multiplier =(Element) finalSubOrRootInWord.getParent().getChild(0);
		}
		Node commonParent =locant.getParent().getParent();//this should be a common parent of the multiplier in front of the root. If it is not, then this locant is in a different scope
		Node parentOfMultiplier =multiplier.getParent();
		while (parentOfMultiplier!=null){
			if (commonParent.equals(parentOfMultiplier)){
				if (locantValues[count-1].endsWith("'")  &&
						multiplier.getLocalName().equals(MULTIPLIER_EL) && multiplier.getAttribute("locantsAssigned")==null &&
						Integer.parseInt(multiplier.getAttributeValue(VALUE_ATR)) == count ){//multiplicative nomenclature
					multiplier.addAttribute(new Attribute ("locantsAssigned",""));
					locant.detach();
					for(int i=locantValues.length-1; i>=0; i--) {
						Element singleLocant = new Element("multiplicativeLocant");
						singleLocant.addAttribute(new Attribute(VALUE_ATR, locantValues[i]));
						XOMTools.insertAfter(multiplier, singleLocant);
					}
					return true;
				}
			}
			parentOfMultiplier=parentOfMultiplier.getParent();
		}
		return false;
	}


	/** Look for multipliers, and multiply out suffixes/unsaturators/heteroatoms/hydros.
	 * Locants are assigned if the number of locants matches the multiplier
	 * associated with them. Eg. triol - > ololol.
	 * Note that infix multiplication is handled seperately as resolution of suffixes is required to perform this unambiguously
	 * @param elem The substituent/root to looks for multipliers in.
	 */
	private void processMultipliers(Element elem) {
		Elements multipliers = elem.getChildElements(MULTIPLIER_EL);
		for(int i=0;i<multipliers.size();i++) {
			Element m = multipliers.get(i);

			ArrayList<Element> locants = new ArrayList<Element>();
			Element possibleLocant =(Element)XOMTools.getPreviousSibling(m);
			while(possibleLocant!=null && possibleLocant.getLocalName().equals(LOCANT_EL)) {
				locants.add(possibleLocant);
				possibleLocant = (Element)XOMTools.getPreviousSibling(possibleLocant);
			}
			Element nextElem = (Element)XOMTools.getNextSibling(m);
			String nextName = nextElem.getLocalName();
			if(nextName.equals(UNSATURATOR_EL) ||
					nextName.equals(SUFFIX_EL) ||
					nextName.equals(HETEROATOM_EL) ||
					nextName.equals(HYDRO_EL)) {
				int mvalue = Integer.parseInt(m.getAttributeValue(VALUE_ATR));
				for(int j=0;j<mvalue;j++) {
					Element newElement =null;
					if (j>0){
						newElement =new Element(nextElem);
					}
					Element referent;
					if (newElement!=null){
						referent=newElement;
					}
					else{
						referent=nextElem;
					}
					if (locants.size()==mvalue){
						Element locant =locants.get(j);
						String locantValue =locant.getAttributeValue(VALUE_ATR);
						//If a compound Locant e.g. 1(6) is detected add a compound locant attribute
						if (locant.getAttribute(COMPOUNDLOCANT_ATR)!=null){
							referent.addAttribute(new Attribute(COMPOUNDLOCANT_ATR, locant.getAttributeValue(COMPOUNDLOCANT_ATR)));
						}
						referent.addAttribute(new Attribute(LOCANT_ATR, locantValue));
						locant.detach();
					}
					else{
						referent.addAttribute(new Attribute("multiplied", "multiplied"));
					}
					if (newElement!=null){
						XOMTools.insertAfter(m, newElement);
					}
				}
				m.detach();
			}
		}
	}


	private void detectConjunctiveSuffixGroups(BuildState state, Element subOrRoot, List<Element> allGroups) throws PostProcessingException, StructureBuildingException {
		List<Element> groups = XOMTools.getChildElementsWithTagName(subOrRoot, GROUP_EL);
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
				throw new PostProcessingException("OPSIN bug: unable to find ring associated with conjunctive suffix group");
			}
			if (conjunctiveGroups.size()!=1){
				throw new PostProcessingException("OPSIN Bug: Two groups exactly should be present at this point when processing conjunctive nomenclature");
			}
			Element primaryConjunctiveGroup =conjunctiveGroups.get(0);
			Fragment primaryConjunctiveFrag = state.xmlFragmentMap.get(primaryConjunctiveGroup);
			//remove all locants
			List<Atom> atomList = primaryConjunctiveFrag.getAtomList();
			for (Atom atom : atomList) {
				atom.clearLocants();
			}
			List<Element> suffixes = new ArrayList<Element>();
			Element possibleSuffix = (Element) XOMTools.getNextSibling(primaryConjunctiveGroup);
			while (possibleSuffix !=null){
				if (possibleSuffix.getLocalName().equals(SUFFIX_EL)){
					suffixes.add(possibleSuffix);
				}
				possibleSuffix = (Element) XOMTools.getNextSibling(possibleSuffix);
			}
			preliminaryProcessSuffixes(state, primaryConjunctiveGroup, suffixes);
			resolveSuffixes(state, primaryConjunctiveGroup, suffixes);
			for (int i = 0; i < suffixes.size(); i++) {
				suffixes.get(i).detach();
			}
			primaryConjunctiveGroup.setLocalName(CONJUNCTIVESUFFIXGROUP_EL);
			allGroups.remove(primaryConjunctiveGroup);
			
			Element possibleMultiplier = (Element) XOMTools.getPreviousSibling(primaryConjunctiveGroup);
			ArrayList<Element> locants = new ArrayList<Element>();
			if (MULTIPLIER_EL.equals(possibleMultiplier.getLocalName())){
				int multiplier = Integer.parseInt(possibleMultiplier.getAttributeValue(VALUE_ATR));
				for (int i = 1; i < multiplier; i++) {
					Element conjunctiveSuffixGroup = new Element(primaryConjunctiveGroup);
					Fragment newFragment = state.fragManager.copyFragment(primaryConjunctiveFrag);
					state.xmlFragmentMap.put(conjunctiveSuffixGroup, newFragment);
					conjunctiveGroups.add(conjunctiveSuffixGroup);
					XOMTools.insertAfter(primaryConjunctiveGroup, conjunctiveSuffixGroup);
				}
				Element possibleLocant =(Element)XOMTools.getPreviousSibling(possibleMultiplier);
				possibleMultiplier.detach();
				while(LOCANT_EL.equals(possibleLocant.getLocalName())) {
					locants.add(possibleLocant);
					possibleLocant = (Element)XOMTools.getPreviousSibling(possibleLocant);
				}
				if (locants.size()>0){
					if (locants.size()!=multiplier){
						throw new PostProcessingException("mismatch between number of locants and multiplier in conjunctive nomenclature routine");
					}
					for (int i = 0; i < locants.size(); i++) {
						conjunctiveGroups.get(i).addAttribute(new Attribute(LOCANT_ATR, locants.get(i).getAttributeValue(VALUE_ATR)));
						locants.get(i).detach();
					}
				}
			}
		}
	}


	/** Match each locant to the next applicable "feature". Assumes that processLocants
	 * has done a good job and rejected cases where no match can be made.
	 * Handles cases where the locant is next to the feature it refers to
	 *
	 * @param subOrRoot The substituent/root to look for locants in.
	 * @throws PostProcessingException
	 */
	private void matchLocantsToDirectFeatures(Element subOrRoot) throws PostProcessingException {
		ArrayList<Element> locants = OpsinTools.elementsToElementArrayList(subOrRoot.getChildElements("locant"));
		Elements groups = subOrRoot.getChildElements("group");
		for (int i = 0; i < groups.size(); i++) {
			Element group = groups.get(i);
			if (group.getAttributeValue("subType").equals("hantzschWidman")){//handle Hantzch-widman systems
				if (group.getAttributeValue("valType")!=null && group.getAttributeValue("valType").equals("partunsatring")){//special case for partunsatring
					//exception for where a locant is supposed to indicate the location of a double bond...
					Elements deltas = subOrRoot.getChildElements("delta");
					if (deltas.size()==0){
						Element delta =new Element("delta");
						if (locants.size()>0 && subOrRoot.indexOf(locants.get(0))< subOrRoot.indexOf(group)){//locant is in front of group
							Element locant=locants.get(0);
							delta.appendChild(locant.getAttributeValue(VALUE_ATR));
							XOMTools.insertBefore(locant, delta);
							locant.detach();
							locants.remove(locant);
						}
						else{
							delta.appendChild("");
							subOrRoot.insertChild(delta, 0);//no obvious attempt to set double bond position, potentially ambiguous, valency will be used to choose later
						}
					}
					group.getAttribute("valType").setValue("ring");
				}
				if (locants.size()>0 ){
					ArrayList<Element> locantsBeforeHWSystem = new ArrayList<Element>();
					ArrayList<Element> heteroAtoms = new ArrayList<Element>();
					int indexOfGroup =subOrRoot.indexOf(group);
					for (int j = indexOfGroup -1; j >= 0; j--) {
						String elName=((Element)subOrRoot.getChild(j)).getLocalName();
						if (elName.equals("locant")){
							locantsBeforeHWSystem.add((Element)subOrRoot.getChild(j));
						}
						else if(elName.equals("heteroatom")){
							heteroAtoms.add((Element)subOrRoot.getChild(j));
						}
						else{
							break;
						}
					}
					if (locantsBeforeHWSystem.size()!=0){
						//detect a solitary locant in front of a HW system and prevent it being assigned.
						//something like 1-aziridin-1-yl never means the N is at position 1 as it is at position 1 by convention
						//this special case is not applied to pseudo HW like systems e.g. [1]oxacyclotetradecine
						if (locantsBeforeHWSystem.size() ==1 && Integer.parseInt(group.getAttributeValue(VALUE_ATR)) <=10){
							locants.remove(locantsBeforeHWSystem.get(0));//don't assign this locant
						}
						else {
							if (locantsBeforeHWSystem.size() == heteroAtoms.size() +1){// e.g. 6-[1,2,3]triazol
								Element locantForSubstituent = locantsBeforeHWSystem.remove(locantsBeforeHWSystem.size()-1);
								locants.remove(locantForSubstituent);
							}
							if (locantsBeforeHWSystem.size() == heteroAtoms.size()){//general case
								for (int j = 0; j < locantsBeforeHWSystem.size(); j++) {
									Element locant =locantsBeforeHWSystem.get(j);
									heteroAtoms.get(j).addAttribute(new Attribute(LOCANT_ATR, locant.getAttributeValue(VALUE_ATR)));
									locant.detach();
									locants.remove(locant);
								}
							}
							else {
								throw new PostProcessingException("Mismatch between number of locants and HW heteroatoms");
							}
						}
					}
				}
			}
		}
	
		for (Element locant : locants) {
			String locantValue =locant.getAttributeValue(VALUE_ATR);
			Element referent = (Element)XOMTools.getNextSibling(locant);
			while(referent.getLocalName().equals("locant") ||
					referent.getAttribute("locant") != null ) {
				referent = (Element)XOMTools.getNextSibling(referent);
			}
			String refName = referent.getLocalName();
			//Only assigning locants to elements that are unsaturator, suffix, heteroatom, hydro and were not created by a multiplier
			if(referent.getAttribute("multiplied") == null && (refName.equals("unsaturator") ||
					refName.equals("suffix") ||
					refName.equals("heteroatom") ||
					refName.equals(CONJUNCTIVESUFFIXGROUP_EL) ||
					(refName.equals("hydro") && !referent.getValue().startsWith("per") ))) {
	
				//If a compound Locant e.g. 1(6) is detected add a compound locant attribute
				if (locant.getAttribute("compoundLocant")!=null){
					referent.addAttribute(new Attribute("compoundLocant", locant.getAttributeValue("compoundLocant")));
				}
				referent.addAttribute(new Attribute(LOCANT_ATR, locantValue));
				locant.detach();
			}
		}
	}


	/**
	 * Handles suffixes, passes them to resolveGroupAddingSuffixes.
	 * Processes the suffixAppliesTo command which multiplies a suffix and attaches the suffixes to the atoms described by the given IDs
	 * @param state
	 * @param group
	 * @param suffixes 
	 * @throws PostProcessingException
	 * @throws StructureBuildingException
	 */
	private void preliminaryProcessSuffixes(BuildState state, Element group, List<Element> suffixes) throws PostProcessingException, StructureBuildingException{
		Fragment suffixableFragment =state.xmlFragmentMap.get(group);

		boolean imideSpecialCase =false;
		if (group.getAttribute(SUFFIXAPPLIESTO_ATR)!=null){//typically a trivial polyAcid or aminoAcid
			//attribute contains instructions for number/positions of suffix
			//this is of the form comma sepeated ids with the number of ids corresponding to the number of instances of the suffix
			Element suffix =OpsinTools.getNextNonChargeSuffix(group);
			if (suffix ==null){
				throw new PostProcessingException("No suffix where suffix was expected");
			}
			if (suffixes.size()>1 && group.getAttributeValue(TYPE_ATR).equals(ACIDSTEM_TYPE_VAL)){
				throw new PostProcessingException("More than one suffix detected on trivial polyAcid. Not believed to be allowed");
			}
			String suffixInstruction =group.getAttributeValue(SUFFIXAPPLIESTO_ATR);
			String[] suffixInstructions = matchComma.split(suffixInstruction);
			boolean symmetricSuffixes =true;
			if (suffix.getAttribute(ADDITIONALVALUE_ATR)!=null){//handles amic, aldehydic, anilic and amoyl suffixes properly
				if (suffixInstructions.length != 2){
					throw new PostProcessingException("suffix: " + suffix.getValue() + " used on an inappropriate group");
				}
				symmetricSuffixes = false;
				if (suffix.getValue().equals("imide")|| suffix.getValue().equals("imido") || suffix.getValue().equals("imidium")  || suffix.getValue().equals("imidylium")){
					imideSpecialCase =true;//prematurely resolve the two suffixes and explicitly join them to form a cyclic imide
				}
			}

			int firstIdInFragment=suffixableFragment.getIdOfFirstAtom();
			if (suffix.getAttribute(LOCANT_ATR)==null){
				suffix.addAttribute(new Attribute("locantID", Integer.toString(firstIdInFragment + Integer.parseInt(suffixInstructions[0]) -1)));
			}
			for (int i = 1; i < suffixInstructions.length; i++) {
				Element newSuffix = new Element(SUFFIX_EL);
				if (symmetricSuffixes){
					newSuffix.addAttribute(new Attribute(VALUE_ATR, suffix.getAttributeValue(VALUE_ATR)));
					newSuffix.addAttribute(new Attribute(TYPE_ATR,  suffix.getAttributeValue(TYPE_ATR)));
					if (suffix.getAttribute(SUBTYPE_ATR)!=null){
						newSuffix.addAttribute(new Attribute(SUBTYPE_ATR,  suffix.getAttributeValue(SUBTYPE_ATR)));
					}
				}
				else{
					newSuffix.addAttribute(new Attribute(VALUE_ATR, suffix.getAttributeValue(ADDITIONALVALUE_ATR)));
					newSuffix.addAttribute(new Attribute(TYPE_ATR, ROOT_EL));
				}
				newSuffix.addAttribute(new Attribute("locantID", Integer.toString(firstIdInFragment + Integer.parseInt(suffixInstructions[i]) -1)));
				XOMTools.insertAfter(suffix, newSuffix);
				suffixes.add(newSuffix);
			}
		}
		else{
			for (Element suffix : suffixes) {
				if (suffix.getAttribute(ADDITIONALVALUE_ATR)!=null){
					throw new PostProcessingException("suffix: " + suffix.getValue() + " used on an inappropriate group");
				}
			}
		}

		ArrayList<Fragment> suffixFragments =resolveGroupAddingSuffixes(state, suffixes, suffixableFragment);
		state.xmlSuffixMap.put(group, suffixFragments);
		boolean suffixesResolved =false;
		if (group.getAttributeValue(TYPE_ATR).equals(CHALCOGENACIDSTEM_TYPE_VAL)){//merge the suffix into the chalcogen acid stem e.g sulfonoate needs to be one fragment for infix replacment
	    	resolveSuffixes(state, group, suffixes);
	    	suffixesResolved =true;
	    }
		processSuffixPrefixes(state, suffixes);//e.g. carbox amide
		processInfixFunctionalReplacementNomenclature(state, suffixes, suffixFragments);
		processConvertHydroxyGroupsToOutAtomsRule(state, suffixes, suffixableFragment);

		if (imideSpecialCase){//Pretty horrible hack to allow cyclic imides
			if (suffixes.size() !=2){
				throw new PostProcessingException("Expected two suffixes fragments for cyclic imide");
			}
			Atom nitrogen =null;
			for (Atom a : suffixFragments.get(0).getAtomList()) {
				if (a.getElement().equals("N")){//amide
					nitrogen =a;
				}
			}
			if (nitrogen ==null){
				throw new PostProcessingException("Nitrogen not found where nitrogen expected");
			}
			Atom carbon = suffixableFragment.getAtomByIDOrThrow(Integer.parseInt(suffixes.get(1).getAttributeValue("locantID")));
			if (!carbon.getElement().equals("C")){
				throw new PostProcessingException("Carbon not found where carbon expected");
			}
			resolveSuffixes(state, group, suffixes);
			suffixesResolved = true;
			state.fragManager.createBond(nitrogen, carbon, 1);//join the N of the amide to the carbon of the acid to form the cyclic imide
		}
		if (suffixesResolved){
			//suffixes have already been resolved so need to be detached to avoid being passed to resolveSuffixes later
			for (int i = suffixes.size() -1; i>=0; i--) {
				Element suffix =suffixes.remove(i);
				suffix.detach();
			}
		}
	}


	/**Processes a suffix and returns any fragment the suffix intends to add to the molecule
	 * @param state
	 * @param suffixes The suffix elements for a fragment.
	 * @param frag The fragment to which the suffix will be applied
	 * @return An arrayList containing the generated fragments
	 * @throws StructureBuildingException If the suffixes can't be resolved properly.
	 * @throws PostProcessingException
	 */
	private ArrayList<Fragment> resolveGroupAddingSuffixes(BuildState state, List<Element> suffixes, Fragment frag) throws StructureBuildingException, PostProcessingException {
		ArrayList<Fragment> suffixFragments =new ArrayList<Fragment>();
		String groupType = frag.getType();
		String subgroupType = frag.getSubType();

		String suffixTypeToUse =null;
		if (suffixApplicability.containsKey(groupType)){
			suffixTypeToUse =groupType;
		}
		else{
			suffixTypeToUse = "standardGroup";
		}

        for (Element suffix : suffixes) {
            String suffixValue = suffix.getAttributeValue(VALUE_ATR);

            boolean cyclic;//needed for addSuffixPrefixIfNonePresentAndCyclic rule
            Atom atomLikelyToBeUsedBySuffix = null;
            if (suffix.getAttribute(LOCANT_ATR) != null) {
            	atomLikelyToBeUsedBySuffix = frag.getAtomByLocant(suffix.getAttributeValue(LOCANT_ATR));
            }
            else if (suffix.getAttribute("locantID") != null) {
            	atomLikelyToBeUsedBySuffix = frag.getAtomByIDOrThrow(Integer.parseInt(suffix.getAttributeValue("locantID")));
            }
            if (atomLikelyToBeUsedBySuffix==null){
            	//a locant has not been specified
            	//also can happen in the cases of things like fused rings where the final numbering is not available so lookup by locant fails (in which case all the atoms will be cyclic anyway)
            	atomLikelyToBeUsedBySuffix = frag.getFirstAtom();
            }
            cyclic = atomLikelyToBeUsedBySuffix.getAtomIsInACycle();

            Elements suffixRuleTags = getSuffixRuleTags(suffixTypeToUse, suffixValue, subgroupType);
            Fragment suffixFrag = null;
            /*
             * Temp fragments are build for each addGroup rule and then merged into suffixFrag
             */
            for (int j = 0; j < suffixRuleTags.size(); j++) {
                Element suffixRuleTag = suffixRuleTags.get(j);
                String suffixRuleTagName = suffixRuleTag.getLocalName();
                if (suffixRuleTagName.equals("addgroup")) {
                    String labels = NONE_LABELS_VAL;
                    if (suffixRuleTag.getAttribute("labels") != null) {
                        labels = suffixRuleTag.getAttributeValue("labels");
                    }
                    suffixFrag = state.fragManager.buildSMILES(suffixRuleTag.getAttributeValue("SMILES"), "suffix", "suffix", labels);
                    List<Atom> atomList = suffixFrag.getAtomList();
                    if (suffixRuleTag.getAttribute("functionalIDs") != null) {
                        String[] relativeIdsOfFunctionalAtoms = matchComma.split(suffixRuleTag.getAttributeValue("functionalIDs"));
                        for (String relativeId : relativeIdsOfFunctionalAtoms) {
                        	int atomIndice = Integer.parseInt(relativeId) -1;
                        	if (atomIndice >=atomList.size()){
                        		throw new StructureBuildingException("Check suffixRules.xml: Atom requested to have a functionalAtom was not within the suffix fragment");
                        	}
                        	suffixFrag.addFunctionalAtom(atomList.get(atomIndice));
						}
                    }
                    if (suffixRuleTag.getAttribute("outIDs") != null) {
                        String[] relativeIdsOfOutAtoms = matchComma.split(suffixRuleTag.getAttributeValue("outIDs"));
                        for (String relativeId : relativeIdsOfOutAtoms) {
                        	int atomIndice = Integer.parseInt(relativeId) -1;
                        	if (atomIndice >=atomList.size()){
                        		throw new StructureBuildingException("Check suffixRules.xml: Atom requested to have a outAtom was not within the suffix fragment");
                        	}
                        	suffixFrag.addOutAtom(atomList.get(atomIndice), 1 , true);
						}
                    }
                }
                else if (suffixRuleTagName.equals("addSuffixPrefixIfNonePresentAndCyclic")){
                	if (cyclic && suffix.getAttribute("suffixPrefix")==null){
                		suffix.addAttribute(new Attribute("suffixPrefix", suffixRuleTag.getAttributeValue("SMILES")));
                	}
                }
				else if (suffixRuleTagName.equals("addFunctionalAtomsToHydroxyGroups")){
					if (suffixFrag != null){
						throw new PostProcessingException("addFunctionalAtomsToHydroxyGroups is not currently compatable with the addGroup suffix rule");
					}
					addFunctionalAtomsToHydroxyGroups(atomLikelyToBeUsedBySuffix);
				}
				else if (suffixRuleTagName.equals("chargeHydroxyGroups")){
					if (suffixFrag != null){
						throw new PostProcessingException("chargeHydroxyGroups is not currently compatable with the addGroup suffix rule");
					}
					chargeHydroxyGroups(atomLikelyToBeUsedBySuffix);
					
				}
				else if (suffixRuleTagName.equals("removeOneDoubleBondedOxygen")){
					if (suffixFrag != null){
						throw new PostProcessingException("removeOneDoubleBondedOxygen is not currently compatable with the addGroup suffix rule");
					}
					removeOneDoubleBondedOxygen(state, atomLikelyToBeUsedBySuffix);
				}
            }
            if (suffixFrag != null) {
				suffixFragments.add(suffixFrag);
				state.xmlFragmentMap.put(suffix, suffixFrag);
            }
        }
		return suffixFragments;
	}
	
	/**Processes any convertHydroxyGroupsToOutAtoms instructions
	 * This is not handled as part of resolveGroupAddingSuffixes as something like carbonochloridoyl involves infix replacement
	 * on a hydroxy that would otherwise actually be removed by this rule!
	 * @param state
	 * @param suffixes The suffix elements for a fragment.
	 * @param frag The fragment to which the suffix will be applied
	 * @return An arrayList containing the generated fragments
	 * @throws PostProcessingException
	 * @throws StructureBuildingException 
	 */
	private void processConvertHydroxyGroupsToOutAtomsRule(BuildState state, List<Element> suffixes, Fragment frag) throws PostProcessingException, StructureBuildingException{
		String groupType = frag.getType();
		String subgroupType = frag.getSubType();
		String suffixTypeToUse =null;
		if (suffixApplicability.containsKey(groupType)){
			suffixTypeToUse =groupType;
		}
		else{
			suffixTypeToUse = "standardGroup";
		}
        for (Element suffix : suffixes) {
            String suffixValue = suffix.getAttributeValue(VALUE_ATR);
            Elements suffixRuleTags = getSuffixRuleTags(suffixTypeToUse, suffixValue, subgroupType);
            for (int j = 0; j < suffixRuleTags.size(); j++) {
                Element suffixRuleTag = suffixRuleTags.get(j);
                String suffixRuleTagName = suffixRuleTag.getLocalName();
                if (suffixRuleTagName.equals("convertHydroxyGroupsToOutAtoms")){
					convertHydroxyGroupsToOutAtoms(state, frag);
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
			if (neighbour.getElement().equals("O") && neighbour.getCharge()==0 && neighbour.getAtomNeighbours().size()==1 && atom.getFrag().findBondOrThrow(atom, neighbour).getOrder()==1){
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
			if (neighbour.getElement().equals("O") && neighbour.getCharge()==0 && neighbour.getAtomNeighbours().size()==1 && atom.getFrag().findBondOrThrow(atom, neighbour).getOrder()==1){
				neighbour.setCharge(-1);
			}
		}
	}

	/**
	 * Removes a double bonded Oxygen from the atom (an [N+][O-] is treated as N=O)
	 * An exception is thrown if no double bonded oxygen could be found connected to the atom
	 * @param state
	 * @param atom
	 * @throws StructureBuildingException
	 */
	private void removeOneDoubleBondedOxygen(BuildState state, Atom atom) throws StructureBuildingException {
		List<Atom> neighbours = atom.getAtomNeighbours();
		for (Atom neighbour : neighbours) {
			if (neighbour.getElement().equals("O") && neighbour.getAtomNeighbours().size()==1){
				Bond b = atom.getFrag().findBondOrThrow(atom, neighbour);
				if (b.getOrder()==2 && neighbour.getCharge()==0){
					state.fragManager.removeAtomAndAssociatedBonds(neighbour);
					if (atom.getLambdaConventionValency()!=null){//corrects valency for phosphin/arsin/stibin
						atom.setLambdaConventionValency(atom.getLambdaConventionValency()-2);
					}
					if (atom.getMinimumValency()!=null){//corrects valency for phosphin/arsin/stibin
						atom.setMinimumValency(atom.getMinimumValency()-2);
					}
					return;
				}
				else if (neighbour.getCharge() ==-1 && b.getOrder()==1){
					if (atom.getCharge() ==1 && atom.getElement().equals("N")){
						state.fragManager.removeAtomAndAssociatedBonds(neighbour);
						atom.setCharge(0);
						return;
					}
				}
			}
		}
		throw new StructureBuildingException("Double bonded oxygen not found in fragment. Perhaps a suffix has been used inappropriately");
	}
	
	/**
	 * 
	 * @param state
	 * @param frag
	 * @throws StructureBuildingException
	 */
	private void convertHydroxyGroupsToOutAtoms(BuildState state, Fragment frag) throws StructureBuildingException {
		List<Atom> atomList = frag.getAtomList();
		for (Atom atom : atomList) {
			List<Atom> neighbours = atom.getAtomNeighbours();
			if (atom.getElement().equals("O") && atom.getCharge()==0 && neighbours.size()==1 && frag.findBondOrThrow(atom, neighbours.get(0)).getOrder()==1){
				state.fragManager.removeAtomAndAssociatedBonds(atom);
				frag.addOutAtom(neighbours.get(0), 1, true);
			}
		}
	}

	/**
	 * Returns the appropriate suffixRule tags for the given arguements.
	 * The suffix rule tags are the children of the appropriate rule in suffixRules.xml
	 * @param suffixTypeToUse
	 * @param suffixValue
	 * @param subgroupType
	 * @return
	 * @throws PostProcessingException
	 */
	private Elements getSuffixRuleTags(String suffixTypeToUse,String suffixValue, String subgroupType) throws PostProcessingException {
		HashMap<String, List<Element>> groupToSuffixMap = suffixApplicability.get(suffixTypeToUse);
		if (groupToSuffixMap==null){
			throw new PostProcessingException("Suffix Type: "+ suffixTypeToUse + " does not have a corresponding groupType entry in suffixApplicability.xml");
		}
		List<Element> potentiallyApplicableSuffixes =groupToSuffixMap.get(suffixValue);
		if(potentiallyApplicableSuffixes==null || potentiallyApplicableSuffixes.size()==0 ) {
			throw new PostProcessingException("Suffix: " +suffixValue +" does not apply to the group it was associated with (type: "+  suffixTypeToUse + ")according to suffixApplicability.xml");
		}
		Element chosenSuffix=null;
        for (Element suffix : potentiallyApplicableSuffixes) {
            if (suffix.getAttribute("subType") != null) {
                if (!suffix.getAttributeValue("subType").equals(subgroupType)) {
                    continue;
                }
            }
            if (chosenSuffix != null) {
                throw new PostProcessingException("Suffix: " + suffixValue + " appears multiple times in suffixApplicability.xml");
            }
            chosenSuffix = suffix;
        }
		if (chosenSuffix==null){
			throw new PostProcessingException("Suffix: " +suffixValue +" does not apply to the group it was associated with (type: "+  suffixTypeToUse + ")due to the group's subType: "+ subgroupType +" according to suffixApplicability.xml");
		}
		Element rule =suffixRules.get(chosenSuffix.getValue());
		if(rule ==null) {
			throw new PostProcessingException("Suffix: " +chosenSuffix.getValue() +" does not have a rule associated with it in suffixRules.xml");
		}
		return rule.getChildElements();
	}


	/**
	 * Searches for suffix elements with the suffixPrefix attribute set
	 * A suffixPrefix is something like sulfon in sulfonamide. It would in this case take the value S(=O)
	 * @param state
	 * @param suffixes
	 * @throws StructureBuildingException
	 */
	private void processSuffixPrefixes(BuildState state, List<Element> suffixes) throws StructureBuildingException{
		for (Element suffix : suffixes) {
			if (suffix.getAttribute(SUFFIXPREFIX_ATR)!=null){
				Fragment suffixPrefixFrag = state.fragManager.buildSMILES(suffix.getAttributeValue(SUFFIXPREFIX_ATR), SUFFIX_TYPE_VAL, NONE_LABELS_VAL);
				addFunctionalAtomsToHydroxyGroups(suffixPrefixFrag.getFirstAtom());
				if (suffix.getValue().endsWith("ate")){
					chargeHydroxyGroups(suffixPrefixFrag.getFirstAtom());
				}
				Atom firstAtomOfPrefix = suffixPrefixFrag.getFirstAtom();
				firstAtomOfPrefix.addLocant("X");
				Fragment suffixFrag = state.xmlFragmentMap.get(suffix);
				state.fragManager.incorporateFragment(suffixPrefixFrag, suffixFrag);
				
				//manipulate suffixFrag such that all the bonds to the first atom (the R)  go instead to the first atom of suffixPrefixFrag.
				//Then reconnect the R to that atom
				Atom theR = suffixFrag.getFirstAtom();
				List<Atom> neighbours = theR.getAtomNeighbours();
				for (Atom neighbour : neighbours) {
					Bond b = suffixFrag.findBondOrThrow(theR, neighbour);
					state.fragManager.removeBond(b);
					state.fragManager.createBond(neighbour, firstAtomOfPrefix, b.getOrder());
				}
				state.fragManager.createBond(firstAtomOfPrefix, theR, 1);
			}
		}
	}
	

	/**
	 * Performs funcional replacement e.g. thio in ethanthioic acid replaces an O with S
	 * @param state
	 * @param suffixFragments Modifed when a special case is encountered, usually untouched
	 * @param suffixes The suffix elements
	 * @throws StructureBuildingException
	 * @throws PostProcessingException
	 */
	private void processInfixFunctionalReplacementNomenclature(BuildState state, List<Element> suffixes, ArrayList<Fragment> suffixFragments) throws StructureBuildingException, PostProcessingException {
		for (int i = 0; i < suffixes.size(); i++) {
			Element suffix = suffixes.get(i);
			if (suffix.getAttribute("infix")!=null){
				Fragment fragToApplyInfixTo = state.xmlFragmentMap.get(suffix);
				Element possibleAcidGroup = XOMTools.getPreviousSiblingIgnoringCertainElements(suffix, new String[]{"multiplier", "infix"});
				if (possibleAcidGroup !=null && possibleAcidGroup.getLocalName().equals(GROUP_EL) && 
						(possibleAcidGroup.getAttributeValue(TYPE_ATR).equals(NONCARBOXYLICACID_TYPE_VAL)|| possibleAcidGroup.getAttributeValue(TYPE_ATR).equals(CHALCOGENACIDSTEM_TYPE_VAL))){
					fragToApplyInfixTo = state.xmlFragmentMap.get(possibleAcidGroup);
				}
				if (fragToApplyInfixTo ==null){
					throw new PostProcessingException("infix has erroneously been assigned to a suffix which does not correspond to a suffix fragment. suffix: " + suffix.getValue());
				}
				//e.g. =O:S,-O:S (which indicates replacing either a double or single bonded oxygen with S)
				//This is semicolon delimited for each infix
				List<String> infixTransformations = StringTools.arrayToList(matchSemiColon.split(suffix.getAttributeValue("infix")));

				List<Atom> atomList =fragToApplyInfixTo.getAtomList();
				LinkedList<Atom> singleBondedOxygen =new LinkedList<Atom>();
				LinkedList<Atom> doubleBondedOxygen =new LinkedList<Atom>();
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
					String[] transformationArray = matchColon.split(infixTransformation);
					if (transformationArray.length !=2){
						throw new StructureBuildingException("Atom to be replaced and replacement not specified correctly in infix: " + infixTransformation);
					}
					String[] transformations = matchComma.split(transformationArray[0]);
					String replacementSMILES = transformationArray[1];
					boolean acceptDoubleBondedOxygen = false;
					boolean acceptSingleBondedOxygen = false;
					for (String transformation : transformations) {
						if (transformation.startsWith("=")){
							acceptDoubleBondedOxygen = true;
						}
						else if (transformation.startsWith("-")){
							acceptSingleBondedOxygen = true;
						}
						else{
							throw new StructureBuildingException("Malformed infix transformation. Expected to start with either - or =. Transformation was: " +transformation);
						}
						if (transformation.length()<2 || transformation.charAt(1)!='O'){
							throw new StructureBuildingException("Only replacement by oxygen is supported. Check infix defintions");
						}
					}
					boolean infixAssignmentAmbiguous =false;
					if (acceptSingleBondedOxygen && !acceptDoubleBondedOxygen){
						if (singleBondedOxygen.size() ==0){
							throw new StructureBuildingException("Cannot find single bonded oxygen for infix with SMILES: "+ replacementSMILES+ " to modify!");
						}
						if (singleBondedOxygen.size() !=1){
							infixAssignmentAmbiguous=true;
						}
					}
					if (!acceptSingleBondedOxygen && acceptDoubleBondedOxygen){
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
					if (acceptDoubleBondedOxygen && doubleBondedOxygen.size()>0 ){
						atomToUse = doubleBondedOxygen.removeFirst();
					}
					else if (acceptSingleBondedOxygen && singleBondedOxygen.size()>0 ){
						atomToUse = singleBondedOxygen.removeFirst();
					}
					else{
						throw new StructureBuildingException("Cannot find oxygen for infix with SMILES: "+ replacementSMILES+ " to modify!");//this would be a bug
					}
					Fragment replacementFrag =state.fragManager.buildSMILES(replacementSMILES, "suffix", "none");
					int bondOrder =0;
					if (replacementFrag.getOutAtoms().size()>0){
						bondOrder = replacementFrag.getOutAtom(0).getValency();
						replacementFrag.removeOutAtom(0);
						//e.g. in nitrido the replacement is #N so need to set bond order to Oxygen appropriately
						if (atomToUse.getAtomNeighbours().size() >0){
							atomToUse.getFrag().findBondOrThrow(atomToUse.getAtomNeighbours().get(0),atomToUse).setOrder(bondOrder);
						}
						else if (atomToUse.getFrag().getInAtoms().size()>0){
							boolean flag = false;
							for (InAtom inAtom : atomToUse.getFrag().getInAtoms()){
								if (inAtom.getAtom().equals(atomToUse)){
									inAtom.setValency(bondOrder);
									flag =true;
									break;
								}
							}
							if (!flag){
								throw new StructureBuildingException("Cannot find inAtom associated with atom in suffix");
							}
						}
						else{
							throw new StructureBuildingException("OPSIN bug: Could not find inAtom for unconnected suffix and atom has no neighbours!");
						}
					}
					removeOrMoveObsoleteFunctionalAtoms(atomToUse, replacementFrag);
					int charge = atomToUse.getCharge();
					state.fragManager.replaceTerminalAtomWithFragment(atomToUse, replacementFrag.getFirstAtom());
					atomToUse.setCharge(charge);
					
					if (infixAssignmentAmbiguous){
						ambiguousElementAtoms.add(atomToUse);
						if (atomToUse.getProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT)!=null){
							ambiguousElementAtoms.addAll(atomToUse.getProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT));
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

	/**
	 * Given an atom that is to be replaced by a functional replacement fragment
	 * determines whether this atom is a functional atom and if it is whether it is:
	 * if the replacement fragment is only 1 atom and is not O/S/Se/Te the functionalAtom is removed
	 * otherwise a check is done on the last atom in the replacement fragment. e.g. peroxy
	 * If this is O/S/Se/Te the functionalAtom is moved to that atom otherwise it is still removed
	 * @param atomToBeReplaced
	 * @param replacementFrag
	 * @throws StructureBuildingException
	 */
	private void removeOrMoveObsoleteFunctionalAtoms(Atom atomToBeReplaced, Fragment replacementFrag)throws StructureBuildingException {
		List<Atom> replacementAtomList = replacementFrag.getAtomList();
		List<FunctionalAtom> functionalAtoms = atomToBeReplaced.getFrag().getFunctionalAtoms();
		for (int j = functionalAtoms.size()-1; j >=0; j--) {
			FunctionalAtom functionalAtom = functionalAtoms.get(j);
			if (atomToBeReplaced.equals(functionalAtom.getAtom())){
				if (replacementFrag.getAtomList().size()>1){
					Atom terminalAtomOfReplacementFrag = replacementAtomList.get(replacementAtomList.size()-1);
					if (matchChalcogen.matcher(terminalAtomOfReplacementFrag.getElement()).matches()){
						functionalAtom.setAtom(terminalAtomOfReplacementFrag);
						terminalAtomOfReplacementFrag.setCharge(atomToBeReplaced.getCharge());
					}
					else{
						atomToBeReplaced.getFrag().removeFunctionalAtom(j);
					}
					atomToBeReplaced.setCharge(0);
				}
				else{
					if (!matchChalcogen.matcher(replacementAtomList.get(0).getElement()).matches()){
						atomToBeReplaced.getFrag().removeFunctionalAtom(j);
						atomToBeReplaced.setCharge(0);
					}
				}
			}
		}
	}


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
	 * @throws PostProcessingException
	 * @throws StructureBuildingException
	 */
	private void disambiguateMultipliedInfixMeaning(BuildState state,
			List<Element> suffixes, ArrayList<Fragment> suffixFragments,
			Element suffix, Fragment suffixFrag,
			List<String> infixTransformations, int oxygenAvailable)
			throws PostProcessingException, StructureBuildingException {
		Element possibleInfix =(Element) XOMTools.getPreviousSibling(suffix);
		if (possibleInfix.getLocalName().equals("infix")){//the infix is only left when there was ambiguity
			Element possibleMultiplier =(Element) XOMTools.getPreviousSibling(possibleInfix);
			if (possibleMultiplier.getLocalName().equals("multiplier")){
				int multiplierValue =Integer.parseInt(possibleMultiplier.getAttributeValue(VALUE_ATR));
				if (infixTransformations.size() + multiplierValue-1 <=oxygenAvailable){//multiplier means multiply the infix e.g. butandithiate
					for (int j = 1; j < multiplierValue; j++) {
						infixTransformations.add(0, infixTransformations.get(0));
					}
				}
				else{
					LinkedList<Element> locants = new LinkedList<Element>();
					Element possibleLocant =(Element)XOMTools.getPreviousSibling(possibleMultiplier);
					while(possibleLocant!=null && possibleLocant.getLocalName().equals("locant")) {
						locants.add(possibleLocant);
						possibleLocant = (Element)XOMTools.getPreviousSibling(possibleLocant);
					}
					if (locants.size()>0){
						if (locants.size()!=multiplierValue){
							throw new PostProcessingException("Multiplier/locant disagreement when multiplying infixed suffix");
						}
						Element locant =locants.removeLast();
					    suffix.addAttribute(new Attribute(LOCANT_ATR, locant.getAttributeValue(VALUE_ATR)));
						suffix.addAttribute(new Attribute("multiplied", "multiplied"));
						locant.detach();
					}
					for (int j = 1; j < multiplierValue; j++) {//multiplier means multiply the infixed suffix e.g. butandithione
						Element newSuffix =new Element(suffix);
						Fragment newSuffixFrag =state.fragManager.copyFragment(suffixFrag);
						state.xmlFragmentMap.put(newSuffix, newSuffixFrag);
						suffixFragments.add(newSuffixFrag);
						XOMTools.insertAfter(suffix, newSuffix);
						suffixes.add(newSuffix);
						if (locants.size()>0){//assign locants if available
							Element locant =locants.removeFirst();
							newSuffix.getAttribute("locant").setValue(locant.getAttributeValue(VALUE_ATR));
							locant.detach();
						}
					}
				}
				possibleMultiplier.detach();
				possibleInfix.detach();
			}
			else{
				throw new PostProcessingException("Multiplier expected in front of ambiguous infix");
			}
		}
	}


	/**
	 * Performs funcional replacement e.g. thio in thioacetic acid replaces an O with S
	 * For heterocyclic rings this should realy be limited to :
	 * pyran, morpholine, chromene, isochromene and xanthene, chromane and isochromane.
	 * @param state
	 * @param groups
	 * @param substituents
	 * @return boolean: has any functional replacement occured
	 * @throws StructureBuildingException
	 */
	private boolean processPrefixFunctionalReplacementNomenclature(BuildState state, List<Element> groups, List<Element> substituents) throws StructureBuildingException {
		boolean doneSomething =false;
		for (int i = groups.size()-1; i >=0; i--) {
			Element group =groups.get(i);
			if (matchChalogenReplacment.matcher(group.getValue()).matches()|| group.getValue().equals("peroxy")){
				//need to check whether this is an instance of functional replacement by checking the substituent/root it is applying to
				Element substituent =(Element) group.getParent();
				Element nextSubOrBracket = (Element) XOMTools.getNextSibling(substituent);
				if (nextSubOrBracket!=null && (nextSubOrBracket.getLocalName().equals("substituent") ||nextSubOrBracket.getLocalName().equals("root"))){
					Element groupToBeModified = nextSubOrBracket.getFirstChildElement(GROUP_EL);
					if (XOMTools.getPreviousSibling(groupToBeModified)!=null){
						continue;//not 2,2'-thiodipyran
					}
					Element possibleMultiplier = (Element) XOMTools.getPreviousSibling(group);
					List<String> possibleLocants = new ArrayList<String>();//usually empty, contains any locant values
					List<Element> locantsToRemove =new ArrayList<Element>();//usually empty, contains any locant elements
					int numberOfAtomsToReplace =1;//the number of atoms to be functionally replaced, modified by a multiplier e.g. dithio
					if (possibleMultiplier !=null){
						Element possibleLocant ;
						if (possibleMultiplier.getLocalName().equals("multiplier")){
							numberOfAtomsToReplace =Integer.valueOf(possibleMultiplier.getAttributeValue(VALUE_ATR));
							possibleLocant = (Element) XOMTools.getPreviousSibling(possibleMultiplier);
						}
						else{
							possibleLocant = possibleMultiplier;
						}
						while (possibleLocant !=null && possibleLocant.getLocalName().equals("locant") && possibleLocant.getAttribute("type")==null) {
							possibleLocants.add(possibleLocant.getAttributeValue(VALUE_ATR));
							locantsToRemove.add(possibleLocant);
							possibleLocant = (Element) XOMTools.getPreviousSibling(possibleLocant);
						}
						if (possibleLocants.size() >0 && possibleLocants.size() !=numberOfAtomsToReplace){//locants and number of replacements disagree
							if (possibleLocants.size()>1){
								continue;
							}
							possibleLocants.clear();
							locantsToRemove.clear();
						}
					}
					int atomCount = state.xmlFragmentMap.get(group).getAtomList().size();
					String replacementSMILES = group.getAttributeValue(VALUE_ATR);
					Fragment frag = state.xmlFragmentMap.get(groupToBeModified);
					ArrayList<Fragment> suffixes = state.xmlSuffixMap.get(groupToBeModified);
					Set<Atom> replaceableAtoms =new LinkedHashSet<Atom>();

					if (atomCount==1){
						for (Atom atom : frag.getAtomList()) {
							if (atom.getElement().equals("O")){
								replaceableAtoms.add(atom);
							}
						}
						if (possibleLocants.size() >0){//locants are used to indicate replacement on trivial groups
							Set<Atom> atomsToUse =new LinkedHashSet<Atom>();
							for (Atom atom : replaceableAtoms) {
								boolean mainGroupAtom =false;//is the atom part of a chain/ring or like an oxo group of a chain/ring
								for (String locant  : atom.getLocants()) {
									if (matchNumericLocant.matcher(locant).matches()){
										mainGroupAtom =true;
										for (String locantVal : possibleLocants) {
											if (locant.equals(locantVal)){
												 atomsToUse.add(atom);
											}
										}
										break;
									}
								}
								if (!mainGroupAtom){
									for (String locantVal : possibleLocants) {
										 if (OpsinTools.depthFirstSearchForNonSuffixAtomWithLocant(atom, locantVal) != null){
											 atomsToUse.add(atom);
										 }
									}
								}
							}
							if(atomsToUse.size() != numberOfAtomsToReplace){
								if (possibleLocants.size()>1){
									throw new StructureBuildingException("Failed to find the correct number of oxygen using locants:" + possibleLocants);
								}
								possibleLocants.clear();
								locantsToRemove.clear();
								//e.g. -1-thioureidomethyl
							}
							else{
								replaceableAtoms =atomsToUse;
							}
						}
						if (possibleLocants.size() ==0){
							//suffixes are ignored in the case of fused ring systems due to suffixes being null and for all suffixes not on acid stems (except phen-->thiophenol)
							if (suffixes!=null){
								Set<Atom> suffixAtomsToUse =new LinkedHashSet<Atom>();
								ArrayList<Fragment> applicableSuffixes = new ArrayList<Fragment>(suffixes);
								if (!groupToBeModified.getAttributeValue(TYPE_ATR).equals("acidStem") && !groupToBeModified.getValue().equals("phen")){
									//remove all non acid suffixes
									for (Fragment fragment : suffixes) {
										Element suffix = state.xmlFragmentMap.getElement(fragment);
										if (!suffix.getAttributeValue(VALUE_ATR).equals("ic") && !suffix.getAttributeValue(VALUE_ATR).equals("ous")){
											applicableSuffixes.remove(fragment);
										}
									}
								}
								for (Fragment suffixFrag : applicableSuffixes) {
									for (Atom atom : suffixFrag.getAtomList()) {
										if (atom.getElement().equals("O") && atom.getBonds().size()==1 && atom.getIncomingValency()==2){
											suffixAtomsToUse.add(atom);
										}
									}
								}
								for (Fragment suffixFrag : applicableSuffixes) {
									for (Atom atom : suffixFrag.getAtomList()) {
										if (atom.getElement().equals("O") && atom.getBonds().size()==1 && atom.getIncomingValency()==1){
											suffixAtomsToUse.add(atom);
										}
									}
								}
								suffixAtomsToUse.addAll(replaceableAtoms);//suffix atoms are preferable
								replaceableAtoms = suffixAtomsToUse;
							}
						}
					}
					else{
						if (possibleLocants.size()>0){
							continue;//you do not locant peroxy
						}
						//only consider ic suffixes when it's not simple atom replacement e.g. peroxy.
						ArrayList<Element> suffixElements =XOMTools.getNextSiblingsOfType(groupToBeModified, "suffix");
						for (int j = 0; j < suffixElements.size(); j++) {
							if (suffixElements.get(j).getAttributeValue(VALUE_ATR).equals("ic")||suffixElements.get(j).getAttributeValue(VALUE_ATR).equals("ous")){
								Fragment suffixFrag =suffixes.get(j);
								List<Atom> atomList =suffixFrag.getAtomList();
								for (Atom a : atomList) {
									if (a.getElement().equals("O")){
										if (a.getIncomingValency()==1){
											replaceableAtoms.add(a);
										}
									}
								}
							}
						}
					}
					boolean multiplierUsed = false;
					if (numberOfAtomsToReplace >1){
						if (replaceableAtoms.size() >= numberOfAtomsToReplace){
							multiplierUsed =true;
						}
						else{
							numberOfAtomsToReplace=1;
						}
					}
					if (replaceableAtoms.size() >=numberOfAtomsToReplace){//check that there atleast as many oxygens as requested replacements
						boolean prefixAssignmentAmbiguous =false;
						Set<Atom> ambiguousElementAtoms = new HashSet<Atom>();
						if (replaceableAtoms.size() != numberOfAtomsToReplace){
							prefixAssignmentAmbiguous=true;
						}

						int replacementsDone =0;
						for (Atom atomToReplace : replaceableAtoms) {
							if (replacementsDone == numberOfAtomsToReplace){
								ambiguousElementAtoms.add(atomToReplace);
								continue;
							}
							if (atomCount>1){//something like peroxy
								Fragment replacementFrag = state.fragManager.buildSMILES(replacementSMILES, "suffix", "none");
								removeOrMoveObsoleteFunctionalAtoms(atomToReplace, replacementFrag);
								int charge = atomToReplace.getCharge();
								state.fragManager.replaceTerminalAtomWithFragment(atomToReplace, replacementFrag.getFirstAtom());
								atomToReplace.setCharge(charge);
							}
							else{
								state.fragManager.makeHeteroatom(atomToReplace, replacementSMILES, false);
								if (prefixAssignmentAmbiguous){
									ambiguousElementAtoms.add(atomToReplace);
								}
							}
							replacementsDone++;
						}

						if (prefixAssignmentAmbiguous){//record what atoms could have been replaced. Often this ambiguity is resolved later e.g. S-methyl thioacetate
							for (Atom atom : ambiguousElementAtoms) {
								atom.setProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT, ambiguousElementAtoms);
							}
						}

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
						if (multiplierUsed){
							possibleMultiplier.detach();
						}
						for (Element locant : locantsToRemove) {
							locant.detach();
						}
						doneSomething=true;
					}
				}
			}
		}
		return doneSomething;
	}

	/**
	 * Checks through the groups accesible from the startingElement taking into account brackets (i.e. those that it is feasible that the group of the startingElement could substitute onto).
	 * It is assumed that one does not intentionally locant onto something in a deeper level of bracketting (not implicit bracketing). e.g. 2-propyl(ethyl)ammonia will give prop-2-yl
	 * @param state
	 * @param startingElement
	 * @param locant: the locant string to check for the presence of
	 * @return whether the locant was found
	 * @throws StructureBuildingException
	 */
	private boolean checkLocantPresentOnPotentialRoot(BuildState state, Element startingElement, String locant) throws StructureBuildingException {
		boolean foundSibling =false;
		Stack<Element> s = new Stack<Element>();
		s.add(startingElement);
		boolean doneFirstIteration =false;//check on index only done on first iteration to only get elements with an index greater than the starting element
		while (s.size()>0){
			Element currentElement =s.pop();
			Element parent = (Element)currentElement.getParent();
			List<Element> siblings = XOMTools.getChildElementsWithTagNames(parent, new String[]{"bracket", "substituent", "root"});
			int indexOfCurrentElement =parent.indexOf(currentElement);

			for (Element bracketOrSub : siblings) {
				if (!doneFirstIteration && parent.indexOf(bracketOrSub) <= indexOfCurrentElement){
					continue;
				}
				if (bracketOrSub.getLocalName().equals("bracket")){//only want to consider implicit brackets, not proper brackets
					if (bracketOrSub.getAttribute("type")==null){
						continue;
					}
					s.push((Element)bracketOrSub.getChild(0));
				}
				else{
					Element group = bracketOrSub.getFirstChildElement(GROUP_EL);
					Fragment groupFrag =state.xmlFragmentMap.get(group);
					if (groupFrag.hasLocant(locant)){
						return true;
					}
					ArrayList<Fragment> suffixes =state.xmlSuffixMap.get(group);
					if (suffixes!=null){
						for (Fragment suffix : suffixes) {
							if (suffix.hasLocant(locant)){
								return true;
							}
						}
					}
				}
				foundSibling =true;
			}
			doneFirstIteration =true;
		}

		if (!foundSibling){//Special case: anything the group could potentially substitute onto is in a bracket. The bracket is checked recursively
			s = new Stack<Element>();
			s.add(startingElement);
			doneFirstIteration =false;//check on index only done on first iteration to only get elements with an index greater than the starting element
			while (s.size()>0){
				Element currentElement =s.pop();
				Element parent = (Element)currentElement.getParent();
				List<Element> siblings = XOMTools.getChildElementsWithTagNames(parent, new String[]{"bracket", "substituent", "root"});
				int indexOfCurrentElement =parent.indexOf(currentElement);

				for (Element bracketOrSub : siblings) {
					if (!doneFirstIteration && parent.indexOf(bracketOrSub) <= indexOfCurrentElement){
						continue;
					}
					if (bracketOrSub.getLocalName().equals("bracket")){
						s.push((Element)bracketOrSub.getChild(0));
					}
					else{
						Element group = bracketOrSub.getFirstChildElement(GROUP_EL);
						Fragment groupFrag =state.xmlFragmentMap.get(group);
						if (groupFrag.hasLocant(locant)){
							return true;
						}
						ArrayList<Fragment> suffixes =state.xmlSuffixMap.get(group);
						if (suffixes!=null){
							for (Fragment suffix : suffixes) {
								if (suffix.hasLocant(locant)){
									return true;
								}
							}
						}
					}
				}
				doneFirstIteration =true;
			}
		}

		return false;
	}

	/**
	 * Handles Hantzsch-Widman rings. Adds SMILES to the group corresponding to the ring's structure
	 * @param state
	 * @param subOrRoot
	 * @throws StructureBuildingException
	 * @throws PostProcessingException
	 */
	private void processHW(BuildState state, Element subOrRoot) throws StructureBuildingException, PostProcessingException{
		Elements groups = subOrRoot.getChildElements("group");
		for(int i=0; i<groups.size(); i++) {
			Element group = groups.get(i);
			if (!group.getAttributeValue("subType").equals("hantzschWidman")){
				continue;
			}
			String ringType =group.getAttributeValue("valType");
			int ringSize = Integer.parseInt(group.getAttributeValue(VALUE_ATR));
			Element prev = (Element) XOMTools.getPreviousSibling(group);
			ArrayList<Element> prevs = new ArrayList<Element>();
			boolean noLocants = true;
			while(prev != null && prev.getLocalName().equals(HETEROATOM_EL)) {
				prevs.add(prev);
				if(prev.getAttribute("locant") != null) {
					noLocants = false;
				}
				prev = (Element) XOMTools.getPreviousSibling(prev);
			}
			boolean hasNitrogen = false;
			boolean hasSiorGeorSborPb=false;
			for(Element heteroatom : prevs){
				String heteroAtomElement =heteroatom.getAttributeValue(VALUE_ATR);
				if (heteroAtomElement.startsWith("[") && heteroAtomElement.endsWith("]")){
					heteroAtomElement=heteroAtomElement.substring(1, heteroAtomElement.length()-1);
				}
				if (heteroAtomElement.equals("N")){
					hasNitrogen=true;
				}
				if (heteroAtomElement.equals("Si") ||
					heteroAtomElement.equals("Ge") ||
					heteroAtomElement.equals("Sb") ||
					heteroAtomElement.equals("Pb") ){
					hasSiorGeorSborPb =true;
				}
			}
			if (ringSize == 6 && ringType.equals("ring") && !hasNitrogen && hasSiorGeorSborPb && (group.getValue().equals("in") ||group.getValue().equals("an"))){
				throw new PostProcessingException("Blocked HW system (6 member saturated ring with no nitrogen but has Si/Ge/Sb/Pb)");
			}
			String name = "";
			Collections.reverse(prevs);
			for(Element heteroatom : prevs) name += heteroatom.getValue();
			name += group.getValue();
			name = name.toLowerCase();
			Fragment hwRing =state.xmlFragmentMap.get(group);
			if(noLocants && prevs.size() > 0) {
				if(specialHWRings.containsKey(name)) {
					String[] specialRingInformation =specialHWRings.get(name);
					if (specialRingInformation[0].equals("blocked")){
						throw new PostProcessingException("Blocked HW system");
					}
					else if (specialRingInformation[0].equals("saturated")){
						for (Atom a: hwRing.getAtomList()) {
							a.setSpareValency(false);
						}
					}//something like oxazole where by convention locants go 1,3 or a inorganic HW-like system
					for (int j = 1; j < specialRingInformation.length; j++) {
						Atom a =hwRing.getAtomByLocantOrThrow(Integer.toString(j));
						a.setElement(specialRingInformation[j]);
					}
					for(Element p : prevs){
						p.detach();
					}
					prevs.clear();
				}
			}
			HashSet<Element> elementsToRemove =new HashSet<Element>();
			for(Element heteroatom : prevs){//add locanted heteroatoms
				if (heteroatom.getAttribute("locant") !=null){
					String locant =heteroatom.getAttributeValue("locant");
					String elementReplacement =heteroatom.getAttributeValue(VALUE_ATR);
					Matcher m = matchElementSymbol.matcher(elementReplacement);
					if (!m.find()){
						throw new PostProcessingException("Failed to extract element from HW heteroatom");
					}
					elementReplacement = m.group();
					Atom a =hwRing.getAtomByLocantOrThrow(locant);
					a.setElement(elementReplacement);
					if (heteroatom.getAttribute("lambda")!=null){
						a.setLambdaConventionValency(Integer.parseInt(heteroatom.getAttributeValue("lambda")));
					}
					heteroatom.detach();
					elementsToRemove.add(heteroatom);
				}
			}
			for(Element p : elementsToRemove){
				prevs.remove(p);
			}

			//add unlocanted heteroatoms
			int defaultLocant=1;
			for(Element heteroatom : prevs){
				String elementReplacement =heteroatom.getAttributeValue(VALUE_ATR);
				Matcher m = matchElementSymbol.matcher(elementReplacement);
				if (!m.find()){
					throw new PostProcessingException("Failed to extract element from HW heteroatom");
				}
				elementReplacement = m.group();

				while (!hwRing.getAtomByLocantOrThrow(Integer.toString(defaultLocant)).getElement().equals("C")){
					defaultLocant++;
				}
				Atom a =hwRing.getAtomByLocantOrThrow(Integer.toString(defaultLocant));
				a.setElement(elementReplacement);
				if (heteroatom.getAttribute("lambda")!=null){
					a.setLambdaConventionValency(Integer.parseInt(heteroatom.getAttributeValue("lambda")));
				}
				heteroatom.detach();
			}

			Elements deltas = subOrRoot.getChildElements("delta");//add specified double bonds
			for (int j = 0; j < deltas.size(); j++) {
				String locantOfDoubleBond = deltas.get(j).getValue();
				Atom firstInDoubleBond;
				Atom secondInDoubleBond;
				if (locantOfDoubleBond.equals("")){
					int defaultId=hwRing.getIdOfFirstAtom();
					firstInDoubleBond =hwRing.getAtomByIDOrThrow(defaultId);
					secondInDoubleBond =hwRing.getAtomByIDOrThrow(defaultId +1);
					while (firstInDoubleBond.hasSpareValency() || !ValencyChecker.checkValencyAvailableForBond(firstInDoubleBond, 1) ||
							secondInDoubleBond.hasSpareValency() || !ValencyChecker.checkValencyAvailableForBond(secondInDoubleBond, 1)){
						defaultId++;
						firstInDoubleBond =hwRing.getAtomByIDOrThrow(defaultId);
						secondInDoubleBond =hwRing.getAtomByIDOrThrow(defaultId +1);
						if (firstInDoubleBond.getType().equals("suffix") || secondInDoubleBond.getType().equals("suffix")){
							throw new StructureBuildingException("No suitable atom found");
						}
					}
				}
				else{
					firstInDoubleBond = hwRing.getAtomByLocantOrThrow(locantOfDoubleBond);
					secondInDoubleBond = hwRing.getAtomByIDOrThrow(firstInDoubleBond.getID() +1);
				}
				Bond b =hwRing.findBond(firstInDoubleBond, secondInDoubleBond);
				b.setOrder(2);
				deltas.get(j).detach();
			}
			XOMTools.setTextChild(group, name);
		}
	}


	/**
	 * Assigns Element symbols to groups and suffixes.
	 * Suffixes have preference.
	 * @param state
	 * @param subOrRoot
	 * @throws StructureBuildingException 
	 */
	private void assignElementSymbolLocants(BuildState state, Element subOrRoot) throws StructureBuildingException {
		Elements groupsOfSubOrRoot = subOrRoot.getChildElements(GROUP_EL);
		Element lastGroupElementInSubOrRoot =groupsOfSubOrRoot.get(groupsOfSubOrRoot.size()-1);
		List<Fragment> suffixFragments =state.xmlSuffixMap.get(lastGroupElementInSubOrRoot);
		Fragment suffixableFragment =state.xmlFragmentMap.get(lastGroupElementInSubOrRoot);
		FragmentTools.assignElementLocants(suffixableFragment, suffixFragments);
	}


	/**
	 * Processes constructs such as biphenyl, 1,1':4',1''-Terphenyl, 2,2'-Bipyridylium, m-Quaterphenyl
	 * @param state
	 * @param subOrRoot
	 * @throws PostProcessingException
	 * @throws StructureBuildingException
	 */
	private void processRingAssemblies(BuildState state, Element subOrRoot) throws PostProcessingException, StructureBuildingException {
		Elements ringAssemblyMultipliers = subOrRoot.getChildElements("ringAssemblyMultiplier");
		for (int i = 0; i < ringAssemblyMultipliers.size(); i++) {
			Element multiplier = ringAssemblyMultipliers.get(i);
			int mvalue = Integer.parseInt(multiplier.getAttributeValue(VALUE_ATR));

			/*
			 * Populate locants with locants. Two locants are required for every pair of rings to be joined.
			 * e.g. bi requires 2, ter requires 4 etc.
			 */
			ArrayList<List<String>> ringJoiningLocants =new ArrayList<List<String>>();
			Element previousEl =(Element)XOMTools.getPreviousSibling(multiplier);
			Element group =(Element)XOMTools.getNextSibling(multiplier, "group");
			if (previousEl!=null && previousEl.getLocalName().equals("ringAssemblyLocant")){//a locant appears to have provided to indicate how to connect the rings of the ringAssembly
				String locantText =StringTools.removeDashIfPresent(previousEl.getValue());
				//special cases where often locants are meant to apply to suffixes rather than being a description of where the rings connect to each other
				if (group.getValue().equals("phen") || group.getValue().equals("hex") || group.getValue().equals("benz")){
					//Find elements that can have locants but don't currently
					List<Element> locantAble = XOMTools.getChildElementsWithTagNames(subOrRoot, new String[]{"suffix", "unsaturator", "heteroatom", "hydro"});
					int locantAbleElements=locantAble.size();
					for(int j=locantAbleElements -1;j >= 0;j--) {
						if ((locantAble.get(j)).getAttribute("locant") !=null){
							locantAble.remove(j);
						}
					}
					if(2 <= locantAble.size()) {
						throw new PostProcessingException("Most likely the ringAssemblyLocant: " + previousEl.getValue() + " is actually a normal locant that is supposed to apply to elements after the ring assembly");
					}
				}
				//locantText might be something like 1,1':3',1''
				String[] perRingLocantArray =matchColon.split(locantText);
				if (perRingLocantArray.length !=(mvalue -1)){
					throw new PostProcessingException("Disagreement between number of locants(" + locantText +") and ring assembly multiplier: " + mvalue);
				}
				for (int j = 0; j < perRingLocantArray.length; j++) {
					String[] locantArray = matchComma.split(perRingLocantArray[j]);
					if (locantArray.length !=2){
						throw new PostProcessingException("missing locant, expected 2 locants: " + perRingLocantArray[j]);
					}
					ringJoiningLocants.add(Arrays.asList(locantArray));
				}
				previousEl.detach();
			}
			else if (previousEl!=null && previousEl.getLocalName().equals("locant")){
				if (previousEl.getAttribute("type")!=null && previousEl.getAttributeValue(TYPE_ATR).equals("orthoMetaPara")){//an OMP locant appears to have provided to indicate how to connect the rings of the ringAssembly
					String locant2 =previousEl.getAttributeValue(VALUE_ATR);
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
					previousEl.detach();
				}
			}

			Element elementToResolve = new Element("substituent");//temporary element containing elements that should be resolved before the ring is duplicated
			Element nextEl =(Element) XOMTools.getNextSibling(multiplier);
			if (nextEl.getLocalName().equals("structuralOpenBracket")){//brackets have been provided to aid disambiguation. These brackets are detached
				Element currentEl =nextEl;
				nextEl = (Element) XOMTools.getNextSibling(currentEl);
				currentEl.detach();
				while (nextEl !=null && !nextEl.getLocalName().equals("structuralCloseBracket")){
					currentEl =nextEl;
					nextEl = (Element) XOMTools.getNextSibling(currentEl);
					currentEl.detach();
					elementToResolve.appendChild(currentEl);
				}
				if (nextEl!=null){
					nextEl.detach();
				}
			}
			else{
				int groupFound = 0;
				int inlineSuffixSeen = 0;
				while (nextEl !=null){
					Element currentEl =nextEl;
					nextEl = (Element) XOMTools.getNextSibling(currentEl);
					if (groupFound==0 ||
							(inlineSuffixSeen == 0 && currentEl.getLocalName().equals("suffix") && currentEl.getAttributeValue(TYPE_ATR).equals("inline") && currentEl.getAttribute("locant")==null)||
							(currentEl.getLocalName().equals("suffix") && currentEl.getAttributeValue(TYPE_ATR).equals("charge"))){
						currentEl.detach();
						elementToResolve.appendChild(currentEl);
					}
					else{
						break;
					}
					if (currentEl.getLocalName().equals(GROUP_EL)){
						groupFound = 1;
					}
					if ((currentEl.getLocalName().equals("suffix") && currentEl.getAttributeValue(TYPE_ATR).equals("inline"))){
						inlineSuffixSeen = 1;
					}
				}
			}

			List<Element> suffixes = XOMTools.getChildElementsWithTagName(elementToResolve, SUFFIX_EL);
			Fragment fragmentToResolveAndDuplicate =state.xmlFragmentMap.get(group);
			resolveSuffixes(state, group, suffixes);
			StructureBuildingMethods.resolveLocantedFeatures(state, elementToResolve);
			StructureBuildingMethods.resolveUnLocantedFeatures(state, elementToResolve);
			group.detach();
			XOMTools.insertAfter(multiplier, group);

			int bondOrder = 1;
			if (fragmentToResolveAndDuplicate.getOutAtoms().size()>0){//e.g. bicyclohexanylidene
				bondOrder =fragmentToResolveAndDuplicate.getOutAtom(0).getValency();
				fragmentToResolveAndDuplicate.removeOutAtom(0);
			}
			if (fragmentToResolveAndDuplicate.getOutAtoms().size()>0){
				throw new StructureBuildingException("Ring assembly fragment should have one or no OutAtoms; not more than one!");
			}

			ArrayList<Fragment> clonedFragments = new ArrayList<Fragment>();
			for (int j = 1; j < mvalue; j++) {
				clonedFragments.add(state.fragManager.copyAndRelabelFragment(fragmentToResolveAndDuplicate, StringTools.multiplyString("'", j)));
			}
			for (int j = 0; j < mvalue-1; j++) {
				Fragment clone =clonedFragments.get(j);
				Atom atomOnParent;
				Atom atomOnLatestClone;
				if (ringJoiningLocants.size()>0){//locants defined
					atomOnParent = fragmentToResolveAndDuplicate.getAtomByLocantOrThrow(ringJoiningLocants.get(j).get(0));
					atomOnLatestClone = clone.getAtomByLocantOrThrow(ringJoiningLocants.get(j).get(1));
				}
				else{
					atomOnParent =fragmentToResolveAndDuplicate.getAtomOrNextSuitableAtomOrThrow(fragmentToResolveAndDuplicate.getDefaultInAtom(), bondOrder);
					atomOnLatestClone = clone.getAtomOrNextSuitableAtomOrThrow(clone.getDefaultInAtom(), bondOrder);
				}
				state.fragManager.incorporateFragment(clone, atomOnLatestClone.getID(), fragmentToResolveAndDuplicate, atomOnParent.getID(), bondOrder);
			}
			XOMTools.setTextChild(group, multiplier.getValue() +group.getValue());
			multiplier.detach();
		}
	}


	private void processPolyCyclicSpiroNomenclature(BuildState state, Element subOrRoot) throws PostProcessingException, StructureBuildingException {
		List<Element> polyCyclicSpiros = XOMTools.getChildElementsWithTagName(subOrRoot, POLYCYCLICSPIRO_EL);
		if (polyCyclicSpiros.size()>0){
			if (polyCyclicSpiros.size()!=1){
				throw new PostProcessingException("Nested polyspiro systems are not supported");
			}
			Element polyCyclicSpiroDescriptor = polyCyclicSpiros.get(0);
			String value = polyCyclicSpiroDescriptor.getAttributeValue(VALUE_ATR);
			String[] tempArray = matchComma.split(value);
			int spiros = Integer.parseInt(tempArray[0]);
			boolean identicalComponents = tempArray[1].equals("identicalComponents");
			if (!identicalComponents){
				throw new PostProcessingException("Only identical components supported so far");
			}
			List<Element> locants = new ArrayList<Element>();
			Element expectedLocant = (Element) XOMTools.getPreviousSibling(polyCyclicSpiroDescriptor);
			while (expectedLocant!=null && expectedLocant.getLocalName().equals(LOCANT_EL)){
				locants.add(expectedLocant);
				expectedLocant = (Element) XOMTools.getPreviousSibling(expectedLocant);
			}
			Collections.reverse(locants);
			if (locants.size()!=spiros){
				throw new PostProcessingException("Mismatch between spiro descriptor and number of locants provided");
			}
			Element group = (Element) XOMTools.getNextSibling(polyCyclicSpiroDescriptor, GROUP_EL);
			if (group==null){
				throw new PostProcessingException("Cannot find group to which spiro descriptor applies");
			}
			Element substituentToResolve = new Element(SUBSTITUENT_EL);//temporary element containing elements that should be resolved before the ring is cloned

			Element possibleOpenBracket = (Element) XOMTools.getNextSibling(polyCyclicSpiroDescriptor);
			List<Element> elementsToResolve;
			if (possibleOpenBracket.getLocalName().equals(STRUCTURALOPENBRACKET_EL)){
				possibleOpenBracket.detach();
				elementsToResolve = XOMTools.getSiblingsUpToElementWithTagName(polyCyclicSpiroDescriptor, STRUCTURALCLOSEBRACKET_EL);
				XOMTools.getNextSibling(elementsToResolve.get(elementsToResolve.size()-1)).detach();//detach close bracket
			}
			else{
				elementsToResolve = XOMTools.getSiblingsUpToElementWithTagName(polyCyclicSpiroDescriptor, GROUP_EL);
			}
			for (Element element : elementsToResolve) {
				element.detach();
				substituentToResolve.appendChild(element);
			}
			if (substituentToResolve.getChildElements().size()!=0){
				group.detach();
				substituentToResolve.appendChild(group);
				StructureBuildingMethods.resolveLocantedFeatures(state, substituentToResolve);
				StructureBuildingMethods.resolveUnLocantedFeatures(state, substituentToResolve);
				group.detach();
				XOMTools.insertAfter(polyCyclicSpiroDescriptor, group);
			}
			Fragment fragment = state.xmlFragmentMap.get(group);
			List<Fragment> clones = new ArrayList<Fragment>();
			for (int i = 1; i < spiros ; i++) {
				clones.add(state.fragManager.copyAndRelabelFragment(fragment, StringTools.multiplyString("'", i)));
			}
			for (Fragment clone : clones) {
				state.fragManager.incorporateFragment(clone, fragment);
			}
			
			Atom atomOnOriginalFragment = fragment.getAtomByLocantOrThrow(locants.get(0).getAttributeValue(VALUE_ATR));
			for (int i = 1; i < spiros ; i++) {
				Atom atomToBeReplaced = fragment.getAtomByLocantOrThrow(locants.get(i).getAttributeValue(VALUE_ATR));
				state.fragManager.replaceAtomWithAnotherAtomPreservingConnectivity(atomToBeReplaced, atomOnOriginalFragment);
			}
			for (Element locant : locants) {
				locant.detach();
			}
			XOMTools.setTextChild(group, polyCyclicSpiroDescriptor.getValue() + group.getValue());
			polyCyclicSpiroDescriptor.detach();
		}
	}

	/**
	 * Searches for lambdaConvention elements and applies the valency they specify to the atom
	 * they specify on the substituent/root's fragment
	 * @param state
	 * @param subOrRoot
	 * @throws StructureBuildingException
	 */
	private void applyLambdaConvention(BuildState state, Element subOrRoot) throws StructureBuildingException {
		List<Element> lambdaConventionEls = XOMTools.getChildElementsWithTagName(subOrRoot, LAMBDACONVENTION_EL);
		for (Element lambdaConventionEl : lambdaConventionEls) {
			Fragment frag = state.xmlFragmentMap.get(subOrRoot.getFirstChildElement(GROUP_EL));
			if (lambdaConventionEl.getAttribute(LOCANT_ATR)!=null){
				frag.getAtomByLocantOrThrow(lambdaConventionEl.getAttributeValue(LOCANT_ATR)).setLambdaConventionValency(Integer.parseInt(lambdaConventionEl.getAttributeValue(LAMBDA_ATR)));
			}
			else{
				if (frag.getAtomList().size()!=1){
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
	 * @param state
	 * @param subOrRoot: The sub/root to look in
	 * @throws PostProcessingException
	 * @throws StructureBuildingException
	 */
	private void handleMultiRadicals(BuildState state, Element subOrRoot) throws PostProcessingException, StructureBuildingException{
		Element group =subOrRoot.getFirstChildElement(GROUP_EL);
		String groupValue =group.getValue();
		Fragment thisFrag = state.xmlFragmentMap.get(group);
		if (groupValue.equals("methylene") || matchChalogenReplacment.matcher(groupValue).matches()){//resolves for example trimethylene to propan-1,3-diyl or dithio to disulfan-1,2-diyl. Locants may not be specified before the multiplier
			Element beforeGroup =(Element) XOMTools.getPreviousSibling(group);
			if (beforeGroup!=null && beforeGroup.getLocalName().equals("multiplier") && beforeGroup.getAttributeValue(TYPE_ATR).equals("basic") && XOMTools.getPreviousSibling(beforeGroup)==null){
				int multiplierVal = Integer.parseInt(beforeGroup.getAttributeValue(VALUE_ATR));
				Element afterGroup = (Element) XOMTools.getNext(group);
				if (((afterGroup !=null && afterGroup.getLocalName().equals("multiplier")&& Integer.parseInt(afterGroup.getAttributeValue(VALUE_ATR)) == multiplierVal) || afterGroup==null) &&
						OpsinTools.getPreviousGroup(group)!=null && ((Element)OpsinTools.getPreviousGroup(group)).getAttribute("isAMultiRadical")!=null &&
						!((Element)XOMTools.getPrevious(beforeGroup)).getLocalName().equals("multiplier")){
					//Something like nitrilotrithiotriacetic acid or oxetane-3,3-diyldimethylene:
					//preceeded by a multiplier that is equal to the multiplier that follows it (or nothing follows it) and the initial multiplier is not proceded by another multiplier e.g. bis(dithio)
				}
				else{
					if (groupValue.equals("methylene")){
						group.getAttribute("value").setValue(StringTools.multiplyString("C", multiplierVal));
					}
					else if (groupValue.equals("thio")){
						group.getAttribute("value").setValue(StringTools.multiplyString("S", multiplierVal));
					}
					else if (groupValue.equals("seleno")){
						group.getAttribute("value").setValue(StringTools.multiplyString("[Se]", multiplierVal));
					}
					else if (groupValue.equals("telluro")){
						group.getAttribute("value").setValue(StringTools.multiplyString("[Te]", multiplierVal));
					}
					else{
						throw new PostProcessingException("unexpected group value");
					}
					group.getAttribute("outIDs").setValue("1,"+Integer.parseInt(beforeGroup.getAttributeValue(VALUE_ATR)));
					XOMTools.setTextChild(group, beforeGroup.getValue() + groupValue);
					beforeGroup.detach();
					if (group.getAttribute("labels")!=null){//use standard numbering
						group.getAttribute("labels").detach();
					}
					state.fragManager.removeFragment(thisFrag);
					thisFrag =resolveGroup(state, group);
					state.xmlFragmentMap.put(group, thisFrag);
					group.removeAttribute(group.getAttribute("usableAsAJoiner"));
				}
			}
		}

		if (group.getAttribute("outIDs")!=null){//adds outIDs at the specified atoms
			String[] radicalPositions = matchComma.split(group.getAttributeValue("outIDs"));
			int firstIdInFrag =thisFrag.getIdOfFirstAtom();
            for (String radicalID : radicalPositions) {
                thisFrag.addOutAtom(firstIdInFrag + Integer.parseInt(radicalID) - 1, 1, true);
            }
		}
		int outAtomCount = thisFrag.getOutAtoms().size();
		if (outAtomCount >=2){
			if (groupValue.equals("amine")){//amine is a special case as it shouldn't technically be allowed but is allowed due to it's common usage in EDTA
				Element previousGroup =(Element) OpsinTools.getPreviousGroup(group);
				if (previousGroup==null || state.xmlFragmentMap.get(previousGroup).getOutAtoms().size() < 2){//must be preceded by a multi radical
					throw new PostProcessingException("Invalid use of amine as a substituent!");
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

		int totalOutAtoms = outAtomCount + calculateOutAtomsToBeAddedFromInlineSuffixes(state, group, subOrRoot.getChildElements("suffix"));
		if (totalOutAtoms >= 2){
			group.addAttribute(new Attribute ("isAMultiRadical", Integer.toString(totalOutAtoms)));
		}
	}

	/**
	 * Calculates number of OutAtoms that the resolveSuffixes method will add.
	 * @param state
	 * @param group
	 * @param suffixes
	 * @return numberOfOutAtoms that will be added by resolveSuffixes
	 * @throws PostProcessingException
	 */
	private int calculateOutAtomsToBeAddedFromInlineSuffixes(BuildState state, Element group, Elements suffixes) throws  PostProcessingException {
		int outAtomsThatWillBeAdded = 0;
		Fragment frag = state.xmlFragmentMap.get(group);
		String groupType = frag.getType();
		String subgroupType = frag.getSubType();
		String suffixTypeToUse =null;
		if (suffixApplicability.containsKey(groupType)){
			suffixTypeToUse =groupType;
		}
		else{
			suffixTypeToUse = "standardGroup";
		}

		List<Fragment> suffixList =state.xmlSuffixMap.get(group);

		for (Fragment suffix : suffixList) {
			outAtomsThatWillBeAdded += suffix.getOutAtoms().size();
		}
		for(int i=0;i<suffixes.size();i++) {
			Element suffix = suffixes.get(i);
			String suffixValue = suffix.getAttributeValue(VALUE_ATR);
			Elements suffixRuleTags =getSuffixRuleTags(suffixTypeToUse, suffixValue, subgroupType);
			for(int j=0;j<suffixRuleTags.size();j++) {
				Element suffixRuleTag = suffixRuleTags.get(j);
				String suffixRuleTagName =suffixRuleTag.getLocalName();
				if(suffixRuleTagName.equals("setOutAtom")) {
					outAtomsThatWillBeAdded +=1;
				}
			}
		}
		return outAtomsThatWillBeAdded;
	}

	/**
	 * Corrects something like L-alanyl-L-glutaminyl-L-arginyl-O-phosphono-L-seryl-L-alanyl-L-proline to:
	 * ((((L-alanyl-L-glutaminyl)-L-arginyl)-O-phosphono-L-seryl)-L-alanyl)-L-proline
	 * i.e. substituents go onto the last mentioned amino acid; amino acids chain together to form peptides
	 * @param state
	 * @param groups
	 * @param brackets
	 */
	private void addImplicitBracketsToAminoAcids(BuildState state, List<Element> groups, List<Element> brackets) {
		for (int i = groups.size() -1; i >=0; i--) {
			Element group = groups.get(i);
			if (group.getAttributeValue(TYPE_ATR).equals(AMINOACID_TYPE_VAL)){
				Element subOrRoot = (Element) group.getParent();
				
				//now find the brackets/substituents before this element
				Element previous = (Element) XOMTools.getPreviousSibling(subOrRoot);
				List<Element> previousElements = new ArrayList<Element>();
				while( previous !=null){
					if (!previous.getLocalName().equals(SUBSTITUENT_EL) && !previous.getLocalName().equals(BRACKET_EL)){
						break;
					}
					previousElements.add(previous);
					previous = (Element) XOMTools.getPreviousSibling(previous);
				}
				if (previousElements.size()>0){//an implicit bracket is needed
					Collections.reverse(previousElements);
					Element bracket = new Element(BRACKET_EL);
					bracket.addAttribute(new Attribute(TYPE_ATR, "implicit"));
					Element parent = (Element) subOrRoot.getParent();
					int indexToInsertAt = parent.indexOf(previousElements.get(0));
					for (Element element : previousElements) {
						element.detach();
						bracket.appendChild(element);
					}

					subOrRoot.detach();
					bracket.appendChild(subOrRoot);
					parent.insertChild(bracket, indexToInsertAt);
					brackets.add(bracket);
				}
			}
		}
	}


	/**Looks for places where brackets should have been, and does the same
	 * as findAndStructureBrackets. E.g. dimethylaminobenzene -> (dimethylamino)benzene.
	 * The bracketting in the above case occurs when the substituent that is being procesed is the amino group
	 * @param state
	 * @param brackets
	 * @param substituents: An arraylist of substituent elements
	 * @return Whether the method did something, and so needs to be called again.
	 * @throws StructureBuildingException
	 */
	private void findAndStructureImplictBrackets(BuildState state, List<Element> substituents, List<Element> brackets) throws PostProcessingException, StructureBuildingException {

		for (Element substituent : substituents) {//will attempt to bracket this substituent with the substituent before it
			String firstElInSubName =((Element)substituent.getChild(0)).getLocalName();
			if (firstElInSubName.equals("locant") ||firstElInSubName.equals("multiplier")){
				continue;
			}

			Element substituentGroup = substituent.getFirstChildElement(GROUP_EL);
			String theSubstituentSubType = substituentGroup.getAttributeValue("subType");
			String theSubstituentType = substituentGroup.getAttributeValue(TYPE_ATR);
			//Only some substituents are valid joiners (e.g. no rings are valid joiners). Need to be atleast bivalent.
			if (substituentGroup.getAttribute("usableAsAJoiner")==null){
				continue;
			}
			Fragment frag =state.xmlFragmentMap.get(substituentGroup);

			//there must be an element after the substituent for the implicit bracket to be required
			Element elementAftersubstituent =(Element)XOMTools.getNextSibling(substituent);
			if (elementAftersubstituent ==null ||
					!elementAftersubstituent.getLocalName().equals("substituent") &&
					!elementAftersubstituent.getLocalName().equals("bracket") &&
					!elementAftersubstituent.getLocalName().equals("root")){
				continue;
			}

			//checks that the element before is a substituent or a bracket which will obviously include substituent/s
			//this makes sure there's more than just a substituent in the bracket
			Element elementBeforeSubstituent =(Element)XOMTools.getPreviousSibling(substituent);
			if (elementBeforeSubstituent ==null||
					!elementBeforeSubstituent.getLocalName().equals("substituent") &&
					!elementBeforeSubstituent.getLocalName().equals("bracket")){
				continue;
			}

			//look for hyphen between substituents, this seems to indicate implicit bracketing was not desired e.g. dimethylaminomethane vs dimethyl-aminomethane
			Element elementDirectlyBeforeSubstituent = (Element) XOMTools.getPrevious(substituent.getChild(0));//can't return null as we know elementBeforeSubstituent is not null
			if (elementDirectlyBeforeSubstituent.getLocalName().equals("hyphen")){
				continue;
			}

			//prevents alkyl chains being bracketed together e.g. ethylmethylamine
			//...unless it's something like 2-methylethyl where the first appears to be locanted onto the second
			List<Element> groupElements  = XOMTools.getDescendantElementsWithTagName(elementBeforeSubstituent, "group");//one for a substituent, possibly more for a bracket
			Element lastGroupOfElementBeforeSub =groupElements.get(groupElements.size()-1);
			if (lastGroupOfElementBeforeSub==null){throw new PostProcessingException("No group where group was expected");}
			if (theSubstituentType.equals("chain") && theSubstituentSubType.equals("alkaneStem") &&
					lastGroupOfElementBeforeSub.getAttributeValue(TYPE_ATR).equals("chain") && lastGroupOfElementBeforeSub.getAttributeValue("subType").equals("alkaneStem")){
				boolean placeInImplicitBracket =false;

				Element suffixAfterGroup=(Element)XOMTools.getNextSibling(lastGroupOfElementBeforeSub, "suffix");
				//if the alkane ends in oxy, sulfinyl, sulfonyl etc. it's not a pure alkane (other suffixes don't need to be considered as they would produce silly structures)
				if (suffixAfterGroup !=null && matchInlineSuffixesThatAreAlsoGroups.matcher(suffixAfterGroup.getValue()).matches()){
					placeInImplicitBracket =true;
				}
				//look for locants and check whether they appear to be referring to the other chain
				if (!placeInImplicitBracket){
					Elements childrenOfElementBeforeSubstituent  =elementBeforeSubstituent.getChildElements();
					Boolean foundLocantNotReferringToChain =null;
					for (int i = 0; i < childrenOfElementBeforeSubstituent.size(); i++) {
						String currentElementName = childrenOfElementBeforeSubstituent.get(i).getLocalName();
						if (currentElementName.equals("locant")){
							String locantText =childrenOfElementBeforeSubstituent.get(i).getAttributeValue(VALUE_ATR);
							if(!frag.hasLocant(locantText)){
								foundLocantNotReferringToChain=true;
								break;
							}
							else{
								foundLocantNotReferringToChain=false;
							}
						}
						else if (currentElementName.equals("stereoChemistry")){
						}
						else{
							break;
						}
					}
					if (foundLocantNotReferringToChain !=null && !foundLocantNotReferringToChain){//a locant was found and it appeared to refer to the other chain
						placeInImplicitBracket=true;
					}
				}
				if (!placeInImplicitBracket){
					continue;
				}
			}

			//prevent bracketing to multi radicals unless through substitution they are likely to cease being multiradicals
			if (lastGroupOfElementBeforeSub.getAttribute("isAMultiRadical")!=null && lastGroupOfElementBeforeSub.getAttribute("acceptsAdditiveBonds")==null){
				continue;
			}

			/*
			 * locant may need to be moved. This occurs when the group in elementBeforeSubstituent is not supposed to be locanted onto
			 *  theSubstituentGroup
			 *  e.g. 2-aminomethyl-1-chlorobenzene where the 2 refers to the benzene NOT the methyl
			 */
			ArrayList<Element> locantElements =new ArrayList<Element>();//sometimes moved
			ArrayList<Element> stereoChemistryElements =new ArrayList<Element>();//always moved if bracketing occurs
			Elements childrenOfElementBeforeSubstituent = elementBeforeSubstituent.getChildElements();
			for (int i = 0; i < childrenOfElementBeforeSubstituent.size(); i++) {
				String currentElementName = childrenOfElementBeforeSubstituent.get(i).getLocalName();
				if (currentElementName.equals("stereoChemistry")){
					stereoChemistryElements.add(childrenOfElementBeforeSubstituent.get(i));
				}
				else if (currentElementName.equals("locant")){
					locantElements.add(childrenOfElementBeforeSubstituent.get(i));
				}
				else{
					break;
				}
			}

			//either all locants will be moved, or none
			Boolean moveLocants = false;
			for (Element locant : locantElements) {
				String locantText =locant.getAttributeValue(VALUE_ATR);
				if (lastGroupOfElementBeforeSub.getAttribute("frontLocantsExpected")!=null){
					StringTools.arrayToList(matchComma.split(lastGroupOfElementBeforeSub.getAttributeValue("frontLocantsExpected"))).contains(locantText);
					continue;
				}
				
				//Check the right fragment in the bracket:
				//if it only has 1 then assume locanted substitution onto it not intended. Or if doesn't have the required locant
				if (frag.getAtomList().size()==1 ||	!frag.hasLocant(locantText) || matchElementSymbolOrAminoAcidLocant.matcher(locantText).find()){
					if (checkLocantPresentOnPotentialRoot(state, substituent, locantText)){
						moveLocants =true;//locant location is present elsewhere
					}
					else if (findElementsMissingIndirectLocants(elementBeforeSubstituent, locant).size()==0 || !state.xmlFragmentMap.get(lastGroupOfElementBeforeSub).hasLocant(locantText)){
						moveLocants =true;//the fragment adjacent to the locant doesn't have this locant or doesn't need any indirect locants. Assume it will appear elsewhere later
					}
				}
			}

			if (moveLocants && locantElements.size() > 1){
				Element shouldBeAMultiplierNode = (Element)XOMTools.getNextSibling(locantElements.get(locantElements.size()-1));
				if (shouldBeAMultiplierNode !=null && shouldBeAMultiplierNode.getLocalName().equals("multiplier")){
					Element shouldBeAGroupOrSubOrBracket = (Element)XOMTools.getNextSibling(shouldBeAMultiplierNode);
					if (shouldBeAGroupOrSubOrBracket !=null){
						if (shouldBeAGroupOrSubOrBracket.getLocalName().equals(GROUP_EL) && 
								(matchInlineSuffixesThatAreAlsoGroups.matcher(substituentGroup.getValue()).matches()//e.g. 4,4'-dimethoxycarbonyl-2,2'-bioxazole --> 4,4'-di(methoxycarbonyl)-2,2'-bioxazole
								|| shouldBeAMultiplierNode.getAttributeValue(TYPE_ATR).equals("group"))){//e.g. 2,5-bisaminothiobenzene --> 2,5-bis(aminothio)benzene
							locantElements.add(shouldBeAMultiplierNode);
						}
						else{//don't bracket other complex multiplied substituents (name hasn't given enough hints if indeed bracketing was expected)
							continue;
						}
					}
					else{
						moveLocants =false;
					}
				}
				else{
					moveLocants =false;
				}
			}

			Element bracket = new Element("bracket");
			bracket.addAttribute(new Attribute("type", "implicit"));

            for (Element stereoChemistryElement : stereoChemistryElements) {
            	stereoChemistryElement.detach();
                bracket.appendChild(stereoChemistryElement);
            }
			if (moveLocants){
                for (Element locantElement : locantElements) {
                    locantElement.detach();
                    bracket.appendChild(locantElement);
                }
			}

			/*
			 * Case when a multiplier should be moved
			 * e.g. tripropan-2-yloxyphosphane -->tri(propan-2-yloxy)phosphane or trispropan-2-ylaminophosphane --> tris(propan-2-ylamino)phosphane
			 */
			if (locantElements.size()==0){
				Element possibleMultiplier =childrenOfElementBeforeSubstituent.get(0);
				if (possibleMultiplier.getLocalName().equals("multiplier") && (
						matchInlineSuffixesThatAreAlsoGroups.matcher(substituentGroup.getValue()).matches() || possibleMultiplier.getAttributeValue(TYPE_ATR).equals("group"))){
					if (childrenOfElementBeforeSubstituent.get(1).getLocalName().equals(GROUP_EL)){
						childrenOfElementBeforeSubstituent.get(0).detach();
						bracket.appendChild(childrenOfElementBeforeSubstituent.get(0));
					}
				}
			}

			Element parent = (Element)substituent.getParent();
			int startIndex=parent.indexOf(elementBeforeSubstituent);
			int endIndex=parent.indexOf(substituent);
			for(int i = 0 ; i <= (endIndex-startIndex);i++) {
				Node n = parent.getChild(startIndex);
				n.detach();
				bracket.appendChild(n);
			}
			parent.insertChild(bracket, startIndex);
			brackets.add(bracket);
		}
	}


	/** Match each locant to the next applicable "feature". Assumes that processLocants
	 * has done a good job and rejected cases where no match can be made.
	 * Handles cases where the locant is in front of the group but the feature is after the group
	 * @param state
	 * @param subOrRoot The substituent/root to look for locants in.
	 * @throws StructureBuildingException
	 */
	private void matchLocantsToIndirectFeatures(BuildState state, Element subOrRoot) throws  StructureBuildingException {
		/* Root fragments (or the root in a bracket) can have prefix-locants
		 * that work on suffixes - (2-furyl), 2-propanol, (2-propylmethyl), (2-propyloxy), 2'-Butyronaphthone.
		 */
		Elements locants = subOrRoot.getChildElements("locant");
		List<Element> locantList = new ArrayList<Element>();//needed because some locants can readily be determined to be needed for other purposes
		for (int i = 0; i < locants.size(); i++) {
			Element locant =locants.get(i);
			Element afterLocants =XOMTools.getNextSiblingIgnoringCertainElements(locant, new String[]{"locant"});
			if (afterLocants!=null && afterLocants.getLocalName().equals("multiplier")){//locant should not be followed by a multiplier. c.f. 1,2,3-tributyl 2-acetyloxypropane-1,2,3-tricarboxylate
				continue;
			}
			locantList.add(locant);
		}
		Element group =subOrRoot.getFirstChildElement(GROUP_EL);
		if (locantList.size()==1 && group!=null && group.getAttribute("frontLocantsExpected")!=null){//some trivial retained names like 2-furyl expect locants to be in front of them. For these the indirect intepretation will always be used rather than checking whether 2-(furyl) even makes sense
			Element locant =locantList.get(0);
			String locantValue =locant.getAttributeValue(VALUE_ATR);
			String[] allowedLocants=matchComma.split(group.getAttributeValue("frontLocantsExpected"));
			for (String allowedLocant : allowedLocants) {
				if (locantValue.equals(allowedLocant)){
					Element expectedSuffix =(Element) XOMTools.getNextSibling(group);
					if (expectedSuffix!=null && expectedSuffix.getLocalName().equals("suffix") && expectedSuffix.getAttribute("locant")==null){
						locantList.clear();
						expectedSuffix.addAttribute(new Attribute(LOCANT_ATR, locantValue));
						locant.detach();
					}
					break;
				}
			}
		}
		if (locantList.size()>0 && group!=null){
			boolean allowIndirectLocants =true;
			if(state.currentWordRule == WordRule.multiEster){//special case e.g. 1-benzyl 4-butyl terephthalate (locants do not apply to yls)
				Element parentEl=(Element) subOrRoot.getParent();
				if (parentEl.getLocalName().equals("word") && parentEl.getAttributeValue(TYPE_ATR).equals("substituent") && parentEl.getChildCount()==1 &&
						locants.size()==1 && (locants.get(0).getAttribute("type")==null || !locants.get(0).getAttributeValue(TYPE_ATR).equals("orthoMetaPara"))){
					allowIndirectLocants =false;
				}
			}
			Fragment fragmentAfterLocant =state.xmlFragmentMap.get(group);
			if (fragmentAfterLocant.getAtomList().size()<=1){
				allowIndirectLocants =false;//e.g. prevent 1-methyl as meth-1-yl is extremely unlikely to be the intended result
			}

			ArrayList<Element> locantsToAssignToIndirectFeatures = new ArrayList<Element>();
			ArrayList<Element> locantAble =null;
			if (allowIndirectLocants){
				for (int i = locants.size()-1; i >=0 ; i--) {
					Element locant =locants.get(i);
					String locantValue =locant.getAttributeValue(VALUE_ATR);
					if (i >0 || !checkLocantPresentOnPotentialRoot(state, subOrRoot, locantValue)){//the first locant is most likely a locant indicating where this subsituent should be attached. If the locant cannot be found on a potential root this cannot be the case though (assuming the name is valid of course)
						if (fragmentAfterLocant.hasLocant(locantValue)){//locant not available elsewhere and is available on the group associated with this element
							if (locantAble ==null){
								locantAble = findElementsMissingIndirectLocants(subOrRoot, locant);
							}
							if (locantsToAssignToIndirectFeatures.size() < locantAble.size()){
								locantsToAssignToIndirectFeatures.add(0, locant);//last in first out
							}
						}
						else{//usually indicates the name will fail unless the suffix has the locant or heteroatom replacement will create the locant
							List<Fragment> suffixes =state.xmlSuffixMap.get(group);
							//I do not want to assign element locants as in locants on the suffix as I currently know of no examples where this actually occurs
							if (matchElementSymbolOrAminoAcidLocant.matcher(locantValue).matches() || locant.getAttribute("compoundLocant")!=null){
								continue;
							}
							for (Fragment suffix : suffixes) {
								if (suffix.hasLocant(locantValue)){//e.g. 2'-Butyronaphthone
									Atom dummyRAtom =suffix.getFirstAtom();
									List<Atom> neighbours =dummyRAtom.getAtomNeighbours();
									Bond b =null;
									atomLoop: for (Atom atom : neighbours) {
										List<String> neighbourLocants = atom.getLocants();
										for (String neighbourLocant : neighbourLocants) {
											if (matchNumericLocant.matcher(neighbourLocant).matches()){
												b=suffix.findBondOrThrow(dummyRAtom, atom);
												break atomLoop;
											}
										}
									}
									if (b!=null){
										state.fragManager.removeBond(b);//the current bond between the dummy R and the suffix
										state.fragManager.createBond(dummyRAtom, suffix.getAtomByLocantOrThrow(locantValue), b.getOrder());
										locant.detach();
									}
								}
							}
						}
					}
				}
			}
			for (Element locant : locantsToAssignToIndirectFeatures) {
				Element locantAbleElement =locantAble.get(0);
				locantAble.remove(0);

				String locantValue =locant.getAttributeValue(VALUE_ATR);
				//If a compound Locant e.g. 1(6) is detected add a compound locant attribute
				if (locant.getAttribute("compoundLocant")!=null){
					locantAbleElement.addAttribute(new Attribute("compoundLocant", locant.getAttributeValue("compoundLocant")));
				}
				locantAbleElement.addAttribute(new Attribute(LOCANT_ATR, locantValue));
				locant.detach();
			}
		}

		//put di-carbon modifying suffixes e.g. oic acids, aldehydes on opposite ends of chain
		Elements suffixEls = subOrRoot.getChildElements("suffix");
		for (int i = 0; i < suffixEls.size()-1; i++) {
			Element terminalSuffix1 = suffixEls.get(i);
			if (isATerminalSuffix(terminalSuffix1) && XOMTools.getNextSibling(terminalSuffix1) != null){
				Element terminalSuffix2 =(Element)XOMTools.getNextSibling(terminalSuffix1);
				if (isATerminalSuffix(terminalSuffix2)){
					Element hopefullyAChain = (Element) XOMTools.getPreviousSibling((Element)terminalSuffix1, "group");
					if (hopefullyAChain != null && hopefullyAChain.getAttributeValue(TYPE_ATR).equals("chain")){
						terminalSuffix1.addAttribute(new Attribute(LOCANT_ATR, "1"));
						terminalSuffix2.addAttribute(new Attribute(LOCANT_ATR, Integer.toString(state.xmlFragmentMap.get(hopefullyAChain).getChainLength())));
						break;
					}
				}
			}
		}
	}


	/**
	 * Find elements that can have indirect locants but don't currently
	 * This requirement excludes hydro and heteroatoms as it is assumed that locants for these are always adjacent (or handled by the special HW code in the case of heteroatoms)
	 * @param subOrRoot The subOrRoot of interest
	 * @param locantEl the locant, only elements after it will be considered
	 * @return An arrayList of locantable elements
	 */
	private ArrayList<Element> findElementsMissingIndirectLocants(Element subOrRoot,Element locantEl) {
		ArrayList<Element> locantAble = new ArrayList<Element>();
		Elements childrenOfSubOrBracketOrRoot=subOrRoot.getChildElements();
		for (int j = 0; j < childrenOfSubOrBracketOrRoot.size(); j++) {
			Element el =childrenOfSubOrBracketOrRoot.get(j);
			String name =el.getLocalName();
			if (name.equals("suffix") || name.equals("unsaturator") || name.equals(CONJUNCTIVESUFFIXGROUP_EL)){
				if (el.getAttribute("locant") ==null && el.getAttribute("multiplied")==null){// shouldn't already have a locant or be multiplied (should of already had locants assignd to it if that were the case)
					if (subOrRoot.indexOf(el)>subOrRoot.indexOf(locantEl)){
						locantAble.add(el);
					}
				}
			}
		}
		return locantAble;
	}
	
	/**
	 * Checks whether a suffix element is:
	 * a suffix, an inline suffix OR terminal root suffix, has no current locant
	 * @param suffix
	 * @return
	 */
	private boolean isATerminalSuffix(Element suffix){
        return suffix.getLocalName().equals(SUFFIX_EL) &&
                suffix.getAttribute(LOCANT_ATR) == null &&
                (suffix.getAttributeValue(TYPE_ATR).equals(INLINE_TYPE_VAL) || TERMINAL_SUBTYPE_VAL.equals(suffix.getAttributeValue(SUBTYPE_ATR)));
		}

	private void processConjunctiveNomenclature(BuildState state, Element subOrRoot) throws PostProcessingException, StructureBuildingException {
		List<Element> conjunctiveGroups = XOMTools.getChildElementsWithTagName(subOrRoot, CONJUNCTIVESUFFIXGROUP_EL);
		if (conjunctiveGroups.size()>0){
			Element ringGroup = subOrRoot.getFirstChildElement(GROUP_EL);
			Fragment ringFrag = state.xmlFragmentMap.get(ringGroup);
			if (ringFrag.getOutAtoms().size()!=0 ){
				throw new PostProcessingException("OPSIN Bug: Ring fragment should have no radicals");
			}
			List<Fragment> conjunctiveFragments = new ArrayList<Fragment>();
			for (Element group : conjunctiveGroups) {
				Fragment frag = state.xmlFragmentMap.get(group);
				conjunctiveFragments.add(frag);
			}
			if (conjunctiveGroups.size()==1){
				//label atoms appropriately
				List<Atom> atomList = conjunctiveFragments.get(0).getAtomList();
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
			}
			for (int i = 0; i < conjunctiveFragments.size(); i++) {
				Fragment conjunctiveFragment = conjunctiveFragments.get(i);
				if (conjunctiveGroups.get(i).getAttribute(LOCANT_ATR)!=null){
					state.fragManager.createBond(lastNonSuffixCarbonWithSufficientValency(conjunctiveFragment), ringFrag.getAtomByLocantOrThrow(conjunctiveGroups.get(i).getAttributeValue(LOCANT_ATR)) , 1);
				}
				else{
					state.fragManager.createBond(lastNonSuffixCarbonWithSufficientValency(conjunctiveFragment), ringFrag.getAtomByIdOrNextSuitableAtomOrThrow(ringFrag.getIdOfFirstAtom(), 1) , 1);
				}
				state.fragManager.incorporateFragment(conjunctiveFragment, ringFrag);
			}
		}
	}


	private Atom lastNonSuffixCarbonWithSufficientValency(Fragment conjunctiveFragment) throws PostProcessingException {
		List<Atom> atomList = conjunctiveFragment.getAtomList();
		for (int i = atomList.size()-1; i >=0; i--) {
			Atom a = atomList.get(i);
			if (a.getType().equals(SUFFIX_TYPE_VAL)){
				continue;
			}
			if (!a.getElement().equals("C")){
				continue;
			}
			if (ValencyChecker.checkValencyAvailableForBond(a, 1)){
				return a;
			}
		}
		throw new PostProcessingException("OPSIN Bug: Unable to find non suffix carbon with sufficient valency");
	}


	/**Process the effects of suffixes upon a fragment. 
	 * Unlocanted non-terminal suffixes are not attached yet. All other suffix effects are performed
	 * @param state
	 * @param group The group element for the fragment to which the suffixes will be added
	 * @param suffixes The suffix elements for a fragment.
	 * @throws StructureBuildingException If the suffixes can't be resolved properly.
	 * @throws PostProcessingException
	 */
	private void resolveSuffixes(BuildState state, Element group, List<Element> suffixes) throws StructureBuildingException, PostProcessingException {
		Fragment frag = state.xmlFragmentMap.get(group);
		int firstAtomID = frag.getIdOfFirstAtom();//typically equivalent to locant 1
		List<Atom> atomList =frag.getAtomList();//this instance of atomList will not change even once suffixes are merged into the fragment
		int defaultAtom =0;//indice in atomList
		String groupType = frag.getType();
		String subgroupType = frag.getSubType();
		String suffixTypeToUse =null;
		if (suffixApplicability.containsKey(groupType)){
			suffixTypeToUse =groupType;
		}
		else{
			suffixTypeToUse ="standardGroup";
		}

		List<Fragment> suffixList = state.xmlSuffixMap.get(group);
		for(int i=0;i<suffixes.size();i++) {
			Element suffix = suffixes.get(i);
			String suffixValue = suffix.getAttributeValue(VALUE_ATR);

			String locant = StructureBuildingMethods.getLocant(suffix);
			int idOnParentFragToUse=0;
			if (!locant.equals("0")){
				idOnParentFragToUse =frag.getIDFromLocantOrThrow(locant);
			}
			if (idOnParentFragToUse==0 && suffix.getAttribute("locantID")!=null){
				idOnParentFragToUse = Integer.parseInt(suffix.getAttributeValue("locantID"));
			}
			if (idOnParentFragToUse==0 && (suffixTypeToUse.equals(ACIDSTEM_TYPE_VAL) || suffixTypeToUse.equals(NONCARBOXYLICACID_TYPE_VAL)|| suffixTypeToUse.equals(CHALCOGENACIDSTEM_TYPE_VAL))){//means that e.g. sulfonyl has an explicit outAtom
				idOnParentFragToUse = firstAtomID;
			}

			if (suffix.getAttribute("compoundLocant")!=null){
				throw new StructureBuildingException("Error: Compound locant assigned to suffix!");
			}

			Fragment suffixFrag =null;
			Elements suffixRuleTags =getSuffixRuleTags(suffixTypeToUse, suffixValue, subgroupType);
			for(int j=0;j<suffixRuleTags.size();j++) {
				Element suffixRuleTag = suffixRuleTags.get(j);
				String suffixRuleTagName =suffixRuleTag.getLocalName();
				if (defaultAtom >= atomList.size()){
					defaultAtom=0;
				}
				if(suffixRuleTagName.equals("addgroup")) {
					if (suffixFrag==null){
						if (suffixList.size() <=0){
							throw new PostProcessingException("OPSIN Bug: Suffixlist should not be empty");
						}
						suffixFrag = suffixList.remove(0);//take the first suffix out of the list, it should of been added in the same order that it is now being read.
						
						if (suffixFrag.getFirstAtom().getBonds().size() <=0){
							throw new PostProcessingException("OPSIN Bug: Dummy atom in suffix should have at least one bond to it");
						}
						int bondOrderRequired = suffixFrag.getFirstAtom().getIncomingValency();
						if(idOnParentFragToUse==0) {
							if(suffixRuleTag.getAttribute("ketoneLocant") != null && !atomList.get(defaultAtom).getAtomIsInACycle()) {
								if(defaultAtom == 0) defaultAtom = FragmentTools.findKetoneAtomIndice(frag, defaultAtom);
								idOnParentFragToUse = atomList.get(defaultAtom).getID();
								defaultAtom++;
							}
							else{
								idOnParentFragToUse =atomList.get(defaultAtom).getID();
							}
							idOnParentFragToUse =frag.getAtomByIdOrNextSuitableAtomOrThrow(idOnParentFragToUse, bondOrderRequired).getID();
						}
						
						//create a new bond and associate it with the suffixfrag and both atoms. Remember the suffixFrag has not been imported into the frag yet
						List<Bond> bonds = new ArrayList<Bond>(suffixFrag.getFirstAtom().getBonds());
						for (Bond bondToSuffix : bonds) {
							Atom suffixAtom;
							if (bondToSuffix.getToAtom().getElement().equals("R")){
								suffixAtom = bondToSuffix.getFromAtom();
							}
							else{
								suffixAtom = bondToSuffix.getToAtom();
							}
							Atom parentfragAtom = frag.getAtomByIDOrThrow(idOnParentFragToUse);
							state.fragManager.createBond(parentfragAtom, suffixAtom, bondToSuffix.getOrder());
							if (suffixValue.equals("one") && groupType.equals("ring")){//special case: one acts in a similar way to the hydro tag c.f. tetrahydrobenzen-1,4-dione
								parentfragAtom.setProperty(Atom.KETONE_SUFFIX_ATTACHED, true);
							}
							state.fragManager.removeBond(bondToSuffix);
						}
					}
				} else if(suffixRuleTagName.equals("changecharge")) {
					int chargeChange =Integer.parseInt(suffixRuleTag.getAttributeValue("charge"));
					int protonChange =Integer.parseInt(suffixRuleTag.getAttributeValue("protons"));
					if(idOnParentFragToUse==0){
						//Typically if a locant has not been specified then it was intended to refer to a nitrogen even if the nitrogen is not at locant 1 e.g. isoquinolinium
						//It's extremely rare to want a carbocation so any heteroatom is preferred with preference given to N
						Atom possibleAtom =null;
						for (Atom a : atomList) {
							if (ValencyChecker.getPossibleValencies(a.getElement(), a.getCharge() + chargeChange)==null){//unstable valency so seems unlikely
								continue;
							}
							String element =a.getElement();
							if (element.equals("N")){
								possibleAtom =a;
								break;
							}
							else if (!element.equals("C")){
								if (possibleAtom == null){
									possibleAtom =a;
								}
							}
						}
						if (possibleAtom==null){
							idOnParentFragToUse =atomList.get(defaultAtom).getID();
							defaultAtom++;
						}
						else{
							idOnParentFragToUse =possibleAtom.getID();
						}
					}
					frag.getAtomByIDOrThrow(idOnParentFragToUse).addChargeAndProtons(chargeChange, protonChange);
				}else if(suffixRuleTagName.equals("setOutAtom")) {
					int outValency = suffixRuleTag.getAttribute("outValency") != null ? Integer.parseInt(suffixRuleTag.getAttributeValue("outValency")) : 1;
					if (suffix.getAttribute("suffixPrefix")==null){
						if(idOnParentFragToUse!=0){
							frag.addOutAtom(idOnParentFragToUse, outValency, true);
						}
						else{
							frag.addOutAtom(firstAtomID, outValency, false);
						}
					}
					else{//something like oyl on a ring, which means it is now carbonyl and the outAtom is on the suffix and not frag
						if (suffixFrag ==null){
							throw new StructureBuildingException("OPSIN bug: ordering of elements in suffixRules.xml wrong; setOutAtom found before addGroup");
						}
						Set<Bond> bonds = state.fragManager.getInterFragmentBonds(suffixFrag);
						if (bonds.size()!=1){
							throw new StructureBuildingException("OPSIN bug: Wrong number of bonds between suffix and group");
						}
						for (Bond bond : bonds) {
							if (bond.getFromAtom().getFrag()==suffixFrag){
								suffixFrag.addOutAtom(bond.getFromAtom(), outValency, true);
							}
							else{
								suffixFrag.addOutAtom(bond.getToAtom(), outValency, true);
							}
						}
					}
				}
				else if (suffixRuleTagName.equals("addSuffixPrefixIfNonePresentAndCyclic")){
					//already processed
				}
				else if (suffixRuleTagName.equals("addFunctionalAtomsToHydroxyGroups")){
					//already processed
				}
				else if (suffixRuleTagName.equals("chargeHydroxyGroups")){
					//already processed
				}
				else if (suffixRuleTagName.equals("removeOneDoubleBondedOxygen")){
					//already processed
				}
				else if (suffixRuleTagName.equals("convertHydroxyGroupsToOutAtoms")){
					//already processed
				}
				else{
					throw new StructureBuildingException("Unknown suffix rule:" + suffixRuleTagName);
				}
			}

			if (suffixFrag!=null){//merge suffix frag and parent fragment
				state.fragManager.removeAtomAndAssociatedBonds(suffixFrag.getFirstAtom());//the dummy R atom
				state.fragManager.incorporateFragment(suffixFrag, frag);
			}
		}
	}

	/**
	 * Removes brackets that only contain one element.
	 * Removed brackets are reflected in brackets and substituentsAndRootAndBrackets
	 * @param brackets
	 * @param substituentsAndRootAndBrackets
	 */
	private void removeClarifyingBrackets(List<Element> brackets, List<Element> substituentsAndRootAndBrackets) {
		for (int i = brackets.size()-1; i >=0; i--) {
			Element bracket =brackets.get(i);
			Elements childElements = bracket.getChildElements();
			boolean hyphenPresent = false;
			if (childElements.size()==2){
				for (int j = childElements.size() -1; j >=0; j--) {
					if (childElements.get(j).getLocalName().equals(HYPHEN_EL)){
						hyphenPresent=true;
					}
				}
			}
			if (childElements.size()==1 || hyphenPresent && childElements.size()==2){
				//this bracket is now unnecesary as implicit brackets have already been added and OPSIN by default substitutes onto the rightmost element
				for (int j = childElements.size() -1; j >=0; j--) {
					Element elToBeMoved = childElements.get(j);
					elToBeMoved.detach();
					XOMTools.insertAfter(bracket, elToBeMoved);
				}
				bracket.detach();
				brackets.remove(i);
				substituentsAndRootAndBrackets.remove(bracket);
			}
		}
	}


	/**
	 * Given the right most child of a word:
	 * Checks whether this is multiplied e.g. methylenedibenzene
	 * If it is then it checks for the presence of locants e.g. 4,4'-oxydibenzene which has been changed to oxy-4,4'-dibenzene
	 * An attribute called inLocants is then added that is either "default" or a comma seperated list of locants
	 * @param state
	 * @param rightMostElement
	 * @throws PostProcessingException
	 * @throws StructureBuildingException
	 */
	private void assignLocantsToMultipliedRootIfPresent(BuildState state, Element rightMostElement) throws PostProcessingException, StructureBuildingException {
		Elements multipliers = rightMostElement.getChildElements("multiplier");
		if(multipliers.size() == 1) {
			Element multiplier =multipliers.get(0);
			if (XOMTools.getPrevious(multiplier)==null){
				throw new StructureBuildingException("OPSIN bug: Unacceptable input to function");
			}
			List<Element> locants = XOMTools.getChildElementsWithTagName(rightMostElement, "multiplicativeLocant");
			int multiVal = Integer.parseInt(multiplier.getAttributeValue(VALUE_ATR));
			if (locants.size()==0){
				rightMostElement.addAttribute(new Attribute("inLocants","default"));
			}
			else if (locants.size() == multiVal){
				String locantString="";
				for (int i = 0; i < locants.size(); i++) {
					if (i!=0){
						locantString+=",";
					}
					Element locant = locants.get(i);
					locantString += locant.getAttributeValue(VALUE_ATR);
					locant.detach();
				}
				rightMostElement.addAttribute(new Attribute("inLocants",locantString));
			}
			else{
				throw new PostProcessingException("Mismatch between number of locants and number of roots");
			}
			Element word =(Element) rightMostElement.getParent();
			if(!word.getLocalName().equals("word")){
				throw new StructureBuildingException("Excpected input to function was the child of a word (OPSIN bug)");
			}
		}
	}


	/**
	 * Assigns locants and multipliers to substituents/brackets
	 * If both locants and multipliers are present a final check is done that the number of them agree.
	 * WordLevel multipliers are processed e.g. diethyl ethanoate
	 * Adding a locant to a root or any other group that cannot engage in substitive nomenclature will result in an exception being thrown
	 * An exception is made for cases where the locant could be referring to a position on another word
	 * @param state 
	 * @param subOrBracket
	 * @throws PostProcessingException
	 */
	private void assignLocantsAndMultipliers(BuildState state, Element subOrBracket) throws PostProcessingException {
		List<Element> locants = XOMTools.getChildElementsWithTagName(subOrBracket,LOCANT_EL);
		int multiplier =1;
		Element possibleMultiplier = subOrBracket.getFirstChildElement("multiplier");
		Element parentElem =(Element)subOrBracket.getParent();
		boolean oneBelowWordLevel = parentElem.getLocalName().equals("word");
		if (possibleMultiplier!=null){
				if (oneBelowWordLevel &&
					XOMTools.getNextSibling(subOrBracket) == null &&
					XOMTools.getPreviousSibling(subOrBracket) == null) {
				return;//word level multiplier
			}
			multiplier = Integer.parseInt(possibleMultiplier.getAttributeValue(VALUE_ATR));
			subOrBracket.addAttribute(new Attribute("multiplier", possibleMultiplier.getAttributeValue(VALUE_ATR)));
			//multiplier is INTENTIONALLY not detached. As brackets/subs are only multiplied later on it is neccesary at that stage to determine what elements (if any) are in front of the multiplier
		}
		if(locants.size() > 0) {
			if (subOrBracket.getLocalName().equals("root")){
				throw new PostProcessingException("Unable to assign all locants");
			}
			if (multiplier==1 && oneBelowWordLevel){//locant might be word Level locant
				if (WordType.valueOf(parentElem.getAttributeValue(TYPE_ATR))==WordType.substituent && (XOMTools.getNextSibling(subOrBracket)==null || locants.size()==2)){//something like S-ethyl or S-(2-ethylphenyl) or S-4-tert-butylphenyl
					if (state.currentWordRule == WordRule.ester || state.currentWordRule == WordRule.functionalClassEster || state.currentWordRule == WordRule.multiEster){
						Element locant = locants.remove(0);
						parentElem.addAttribute(new Attribute(LOCANT_ATR, locant.getAttributeValue(VALUE_ATR)));
						locant.detach();
						if (locants.size()==0){
							return;
						}
					}
				}
			}
			if (multiplier !=locants.size()){
				throw new PostProcessingException("Multiplier and locant count failed to agree; All locants could not be assigned!");
			}

			Element parent =(Element) subOrBracket.getParent();
			if (!parent.getLocalName().equals("word")){//attempt to find cases where locant will not be utilised. This if statement allows the use of locants for ester formation
				Elements children =parent.getChildElements();
				boolean foundSomethingToSubstitute =false;
				for (int i = parent.indexOf(subOrBracket) +1 ; i < children.size(); i++) {
					if (!children.get(i).getLocalName().equals("hyphen")){
						foundSomethingToSubstitute = true;
					}
				}
				if (!foundSomethingToSubstitute){
					throw new PostProcessingException("Unable to assign all locants");
				}
			}
			String locantString="";
			for (int i = 0; i < locants.size(); i++) {
				Element locant = locants.get(i);
				if(locant.getAttribute("compoundLocant")!=null){
					throw new PostProcessingException("A compound locant cannot be used to locant a sub/bracket!");
				}
				if (i!=0){
					locantString += ",";
				}
				locantString += locant.getAttributeValue(VALUE_ATR);
				locant.detach();
			}
			subOrBracket.addAttribute(new Attribute(LOCANT_ATR, locantString));
		}
	}

	/**
	 * If a word level multiplier is present e.g. diethyl butandioate then this is processed to ethyl ethyl butandioate
	 * If wordCount is 1 then an exception is thrown if a multiplier is encountered e.g. triphosgene parsed as tri phosgene
	 * @param state
	 * @param word
	 * @param wordCount 
	 * @throws StructureBuildingException
	 * @throws PostProcessingException
	 */
	private void processWordLevelMultiplierIfApplicable(BuildState state, Element word, int wordCount) throws StructureBuildingException, PostProcessingException {
		if (word.getChildCount()==1){
			Element subOrBracket = (Element) word.getChild(0);
			Element multiplier = subOrBracket.getFirstChildElement("multiplier");
			if (multiplier !=null){
				int multiVal =Integer.parseInt(multiplier.getAttributeValue(VALUE_ATR));
				if (multiVal ==1){return;}
				Elements locants =subOrBracket.getChildElements("locant");
				boolean assignLocants =false;
				boolean wordLevelLocants =false;
				if (XOMTools.getNextSibling(subOrBracket)==null && WordType.valueOf(word.getAttributeValue(TYPE_ATR))==WordType.substituent){//something like O,S-dimethyl phosphorothioate
					if (state.currentWordRule == WordRule.ester || state.currentWordRule == WordRule.functionalClassEster || state.currentWordRule == WordRule.multiEster){
						wordLevelLocants =true;
					}
				}
				if (locants.size()==multiVal){
					assignLocants=true;
					for (int i = 0; i < locants.size(); i++) {
						if(locants.get(i).getAttribute("compoundLocant")!=null){
							throw new PostProcessingException("A compound locant cannot be used to locant a sub/bracket!");
						}
						locants.get(i).detach();
					}
					if (wordLevelLocants){
						word.addAttribute(new Attribute(LOCANT_ATR, locants.get(0).getAttributeValue(VALUE_ATR)));
					}
					else{
						subOrBracket.addAttribute(new Attribute(LOCANT_ATR, locants.get(0).getAttributeValue(VALUE_ATR)));
					}
				}
				else if (locants.size()!=0){
					throw new PostProcessingException("Unable to assign all locants");
				}
				List<Element> elementsNotToBeMultiplied = new ArrayList<Element>();//anything before the multiplier
				for (int i = subOrBracket.indexOf(multiplier) -1 ; i >=0 ; i--) {
					Element el = (Element) subOrBracket.getChild(i);
					el.detach();
					elementsNotToBeMultiplied.add(el);
				}
				multiplier.detach();
				if (wordCount ==1){
					throw new StructureBuildingException("Unexpected multiplier found at start of word. Perhaps the name is trivial e.g. triphosgene");
				}
				for(int i=multiVal -1; i>=1; i--) {
					Element clone = state.fragManager.cloneElement(state, word);
					if (assignLocants){
						if (wordLevelLocants){
							clone.getAttribute("locant").setValue(locants.get(i).getAttributeValue(VALUE_ATR));
						}
						else{
							((Element) clone.getChild(0)).getAttribute("locant").setValue(locants.get(i).getAttributeValue(VALUE_ATR));
						}
					}
					XOMTools.insertAfter(word, clone);
				}
				for (Element el : elementsNotToBeMultiplied) {//re-add anything before multiplier to original word
					subOrBracket.insertChild(el, 0);
				}
			}
		}
	}
}
