package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Stack;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import uk.ac.cam.ch.wwmm.ptclib.string.StringTools;
import uk.ac.cam.ch.wwmm.ptclib.xml.XOMTools;

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

public class PreStructureBuilder {

	private FusedRingBuilder fusedRingBuilder;
	
	private Pattern matchCompoundLocant =Pattern.compile("[\\[\\(\\{](\\d+[a-z]?'*)[\\]\\)\\}]");
	private Pattern matchIndicatedHydrogen =Pattern.compile("(\\d+[a-z]?'*)H");
	private Pattern matchBracketedEntryInLocant =Pattern.compile("[\\[\\(\\{].*[\\]\\)\\}]");
	private Pattern matchCisTransInLocants =Pattern.compile("[rct]-");
	private Pattern matchColon =Pattern.compile(":");
	private Pattern matchSemiColon =Pattern.compile(";");
	private Pattern matchComma =Pattern.compile(",");
	private Pattern matchSpace =Pattern.compile(" ");
	private Pattern matchElementSymbol = Pattern.compile("[A-Z].?");
	private Pattern matchUpperCase = Pattern.compile("[A-Z]");
	private Pattern matchOrtho =Pattern.compile("[oO]");
	private Pattern matchMeta =Pattern.compile("[mM]");
	private Pattern matchPara =Pattern.compile("[pP]");
	private Pattern matchChalogenReplacment= Pattern.compile("thio|seleno|telluro|peroxy");
	private Pattern matchInlineSuffixesThatAreAlsoGroups = Pattern.compile("carbon|oxy|sulfin|sulfon");
	private Pattern matchSuffixesThatGoAtEndOfChainsByDefault = Pattern.compile("al|amide|ate|hydrazide|hydrazonic|hydroxamic|hydroximic|ic|imidic|nitrile|oyl|sulfonate");

	//rings that look like HW rings but have other meanings. For the HW like inorganics the true meaning is given
	private HashMap<String, String[]> specialHWRings;

	/*Holds the rules on how suffixes are interpreted. Convenience methods are available to use them*/
	private HashMap<String, HashMap<String, List<Element>>> suffixApplicability;
	private HashMap<String, Element> suffixRules;

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
				String suffixValue= suffix.getAttributeValue("value");
				if (suffixToRuleMap.get(suffixValue)!=null){//can have multiple entries if subType attribute is set
					suffixToRuleMap.get(suffixValue).add(suffix);
				}
				else{
					ArrayList<Element> suffixList =new ArrayList<Element>();
					suffixList.add(suffix);
					suffixToRuleMap.put(suffixValue, suffixList);
				}
			}
			suffixApplicability.put(groupType.getAttributeValue("type"), suffixToRuleMap);
		}

		Elements rules = suffixRulesDoc.getRootElement().getChildElements("rule");
		for (int i = 0; i < rules.size(); i++) {
			Element rule =rules.get(i);
			String ruleValue=rule.getAttributeValue("value");
			if (suffixRules.get(ruleValue)!=null){
				throw new Exception("Suffix: " +ruleValue +" appears multiple times in suffixRules.xml");
			}
			suffixRules.put(ruleValue, rule);
		}


		//The first entry of the array is a special instruction e.g. blocked or saturated. The correct order of the heteroatoms follows
		specialHWRings = new HashMap<String, String[]>();//terminal e is ignored
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

		specialHWRings.put("boroxin", new String[]{"saturated","O","B","O","B","O","B"});
		specialHWRings.put("borazin", new String[]{"saturated","N","B","N","B","N","B"});
		specialHWRings.put("borthiin", new String[]{"saturated","S","B","S","B","S","B"});
	}


	/** The master method, postprocesses a parse result. At this stage one can except all substituents/roots to have at least 1 group.
	 * Multiple groups are present in, for example, fusion nomenclature. By the end of this function there will be exactly 1 group
	 * associated with each substituent/root. Multiplicative nomenclature can result in there being multiple roots
	 *
	 * @param elem The element to postprocess.
	 * @return The postprocessed element. The same as elem.
	 * @throws Exception
	 */
	Element postProcess(Element elem, BuildState state) throws Exception {
		Elements words =elem.getChildElements("word");
		for (int i = 0; i < words.size(); i++) {
			Element word =words.get(i);
			if (word.getAttributeValue("type").equals("literal")){
				continue;
			}

			List<Element> roots = OpsinTools.findDescendantElementsWithTagName(word, "root");
			if (roots.size() >1){
				throw new PostProcessingException("Multiple roots, but only 0 or 1 were expected. Found: " +roots.size());
			}
			List<Element> substituents = OpsinTools.findDescendantElementsWithTagName(word, "substituent");
			List<Element> substituentsAndRoot = OpsinTools.combineElementLists(substituents, roots);
			List<Element> brackets =  OpsinTools.findDescendantElementsWithTagName(word, "bracket");
			List<Element> substituentsAndRootAndBrackets =OpsinTools.combineElementLists(substituentsAndRoot, brackets);
			List<Element> groups =  OpsinTools.findDescendantElementsWithTagName(word, "group");

			Element root =null;
			if (roots.size() ==1) root=roots.get(0);

			for (Element subOrBracketOrRoot : substituentsAndRootAndBrackets) {
				processLocants(subOrBracketOrRoot);
			}

			for (Element group : groups) {
				Fragment thisFrag = resolveGroup(state, group);
				state.xmlFragmentMap.put(group, thisFrag);
			}

			for (Element e : substituentsAndRootAndBrackets) {
				checkAndConvertToSingleLocants(state, e, root);
			}

			for (Element subOrRoot : substituentsAndRoot) {
				processMultipliers(subOrRoot);
				matchLocantsToDirectFeatures(subOrRoot);
				preliminaryProcessSuffixes(state, subOrRoot);
			}

			if (processPrefixFunctionalReplacementNomenclature(state, groups, substituents)){//true if functional replacement performed, 1 or more substituents will have been removed
				substituentsAndRoot = OpsinTools.combineElementLists(substituents, roots);
				substituentsAndRootAndBrackets =OpsinTools.combineElementLists(substituentsAndRoot, brackets);
			}
			groups=null;//groups is potentially out of sync with xml at this point

			for (Element subOrRoot : substituentsAndRoot) {
				processHW(state, subOrRoot);//hantzch-widman rings
				fusedRingBuilder.processFusedRings(state, subOrRoot);
				assignElementSymbolLocants(state, subOrRoot);
				processRingAssemblies(state, subOrRoot);
				processComplicatedSpiroNomenclature(state, subOrRoot);
			}

			//System.out.println(new XOMFormatter().elemToString(elem));
			findAndStructureImplictBrackets(state, substituents, brackets);

			substituentsAndRootAndBrackets =OpsinTools.combineElementLists(substituentsAndRoot, brackets);//findAndStructureImplictBrackets may have created new brackets

			for (Element subOrRoot : substituentsAndRoot) {
				matchLocantsToIndirectFeatures(state, subOrRoot);
				resolveRemainingSuffixes(state, subOrRoot);
			}

			removePointlessBrackets(brackets, substituentsAndRootAndBrackets);//e.g. (tetramethyl)azanium == tetramethylazanium

			for (Element subOrRoot : substituentsAndRoot) {
				handleMultiRadicals(state, subOrRoot);
			}
			
			if (word.getChildCount()>1){
				assignLocantsToMultipliedRootIfPresent(state, (Element) word.getChild(word.getChildCount()-1));//multiplicative nomenclature e.g. methylenedibenzene or 3,4'-oxydipyridine
			}
			
			for (Element subBracketOrRoot : substituentsAndRootAndBrackets) {
				assignLocantsAndMultipliers(subBracketOrRoot);
			}
			processWordLevelMultiplierIfApplicable(state, word);
			
		}

		return elem;
	}

	/**Handles special features of locants e.g. ortho/meta/para, indicated hydrogen, cis/trans in locant
	 *
	 * @param elem The substituent/root/bracket to looks for locants in.
	 * @throws PostProcessingException
	 */
	private void processLocants(Element elem) throws PostProcessingException{
		Elements ompLocants = elem.getChildElements("orthoMetaPara");
		for(int i=0;i<ompLocants.size();i++) {
			Element locant = ompLocants.get(i);
			String locantText = locant.getValue();
			locantText = locantText.substring(0, 1);
			Element afterOmpLocant = (Element)XOMTools.getNextSibling(locant);
			locant.setLocalName("locant");
			locant.removeChildren();
			locant.addAttribute(new Attribute("type", "orthoMetaPara"));
			if(afterOmpLocant.getLocalName().equals("multiplier") || afterOmpLocant.getAttribute("outIDs")!=null) {
				if (matchOrtho.matcher(locantText).matches()){
					locant.appendChild("1,2-");
				}
				else if (matchMeta.matcher(locantText).matches()){
					locant.appendChild("1,3-");
				}
				else if (matchPara.matcher(locantText).matches()){
					locant.appendChild("1,4-");
				}
				else{
					throw new PostProcessingException(locantText + " was not identified as being either ortho, meta or para but according to the chemical grammar it should of been");
				}
			}
			else{
				if (matchOrtho.matcher(locantText).matches()){
					locant.appendChild("2-");
				}
				else if (matchMeta.matcher(locantText).matches()){
					locant.appendChild("3-");
				}
				else if (matchPara.matcher(locantText).matches()){
					locant.appendChild("4-");
				}
				else{
					throw new PostProcessingException(locantText + " was not identified as being either ortho, meta or para but according to the chemical grammar it should of been");
				}
			}
		}

		Elements locants = elem.getChildElements("locant");
		for(int i=0;i<locants.size();i++) {
			Element locant = locants.get(i);
			String locantText = OpsinTools.removeDashIfPresent(locant.getValue());

			//If the indicatedHydrogen has been specified create a tag for it and remove it from the list of locants
			//e.g. 1(9H),5,7 -->indicatedHydrogen tag value (9H) and 1,5,7
			//can get as complicated as 1,2(2H,7H)
			Matcher matches =matchIndicatedHydrogen.matcher(locantText);
			if (matches.find()==true){
				do {
					Element indicatedHydrogenElement=new Element("indicatedHydrogen");
					indicatedHydrogenElement.addAttribute(new Attribute("locant", matches.group(1)));
					XOMTools.insertBefore(locant, indicatedHydrogenElement);
				}
				while (matches.find());
				locantText =matchBracketedEntryInLocant.matcher(locantText).replaceAll("");
			}

			/*
			 * Strip out cis/trans information built into locant - currently unhandled
			 */
			matches =matchCisTransInLocants.matcher(locantText);
			if (matches.find()==true){
				do {
					//currently do nothing
				}
				while (matches.find());
				locantText =matches.replaceAll("");
			}
			OpsinTools.setTextChild(locant, locantText);

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
		String groupType = group.getAttributeValue("type");
		String groupSubType = group.getAttributeValue("subType");
		String groupValue = group.getAttributeValue("value");
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
		CycleDetector.assignWhetherAtomsAreInCycles(thisFrag);
		thisFrag.convertHighOrderBondsToSpareValencies();//only applied to cyclic bonds

		if (group.getAttribute("defaultInLocant")!=null){//sets the atom at which substitution will occur to by default
			thisFrag.setDefaultInID(thisFrag.getAtomByLocantOrThrow(group.getAttributeValue("defaultInLocant")).getID());
		}
		else if (group.getAttribute("defaultInID")!=null){
			thisFrag.setDefaultInID(thisFrag.getIdOfFirstAtom() + Integer.parseInt(group.getAttributeValue("defaultInID")) -1);
		}
		else if (group.getAttribute("usableAsAJoiner") != null && group.getAttributeValue("usableAsAJoiner").equals("yes")){//makes linkers by default attach end to end
			int chainLength =thisFrag.getChainLength();
			if (chainLength >1){
				boolean connectEndToEndWithPreviousSub =true;
				if (groupSubType.equals("alkaneStem")){//don't do this if you the group is preceded by another alkanStem e.g. methylethyl makes more sense as prop-2-yl rather than propyl
					Element previousSubstituent =(Element) XOMTools.getPreviousSibling(group.getParent());
					if (previousSubstituent!=null){
						Elements groups = previousSubstituent.getChildElements("group");
						if (groups.size()==1 && groups.get(0).getAttributeValue("subType").equals("alkaneStem") && !groups.get(0).getAttributeValue("type").equals("ring")){
							connectEndToEndWithPreviousSub = false;
						}
					}
				}
				if (connectEndToEndWithPreviousSub){
					group.addAttribute(new Attribute("defaultInID",Integer.toString(chainLength)));
					thisFrag.setDefaultInID(thisFrag.getIDFromLocantOrThrow(Integer.toString(chainLength)));
				}
			}
		}

		if (group.getAttribute("functionalIDs")!=null){
			String[] functionalIDs = matchComma.split(group.getAttributeValue("functionalIDs"));
			for (int i = 0; i < functionalIDs.length; i++) {
				thisFrag.addFunctionalID(thisFrag.getIdOfFirstAtom() +Integer.parseInt(functionalIDs[i]) -1);
			}
		}

		if (thisFrag.getOutIDs().size()==0 && groupType.equals("substituent") && groupSubType.equals("simpleSubstituent")){
			//simple substituents implicitly will be given an outID assuming the SMILESbuilder hasn't already given them one
			thisFrag.addOutID(thisFrag.getIdOfFirstAtom(), 1, true);
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
			for (int i = 0; i < groupsToBeAdded.length; i++) {//populate allGroupInformation list
				String groupToBeAdded = groupsToBeAdded[i];
				String[] tempArray =matchSpace.split(groupToBeAdded);
				HashMap<String, String> groupInformation = new HashMap<String, String>();
				if (tempArray.length!=2 && tempArray.length!=3){
					throw new PostProcessingException("malformed addGroup tag");
				}
				groupInformation.put("SMILES", tempArray[0]);
				if (tempArray[1].startsWith("id")){
					groupInformation.put("atomReferenceType", "id");
					groupInformation.put("atomReference", tempArray[1].substring(2));
				}
				else if (tempArray[1].startsWith("locant")){
					groupInformation.put("atomReferenceType", "locant");
					groupInformation.put("atomReference", tempArray[1].substring(6));
				}
				else{
					throw new PostProcessingException("malformed addGroup tag");
				}
				if (tempArray.length==3){//labels may optionally be specified for the group to be added
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
				if (assignlocants==true){
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
						OpsinTools.setTextChild(previousEl, StringTools.objectListToString(locantValues, ","));
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
				if (newFrag.getOutIDs().size() >1){
					throw new PostProcessingException("too many outIDs on group to be added");
				}
				if (newFrag.getOutIDs().size() ==1) {
					OutID newFragOutID = newFrag.getOutID(0);
					state.fragManager.incorporateFragment(newFrag, newFragOutID.id, parentFrag, atomOnParentFrag.getID(), newFragOutID.valency);
				}
				else{
					Atom atomOnNewFrag =newFrag.getAtomByIDOrThrow(newFrag.getDefaultInID());
					state.fragManager.incorporateFragment(newFrag, atomOnNewFrag.getID(), parentFrag, atomOnParentFrag.getID(), 1);
				}
			}
		}

		if(group.getAttributeValue("addHeteroAtom")!=null) {
			String addHeteroAtomInformation=group.getAttributeValue("addHeteroAtom");
			String[] heteroAtomsToBeAdded = matchSemiColon.split(addHeteroAtomInformation);
			ArrayList<HashMap<String, String>> allHeteroAtomInformation = new ArrayList<HashMap<String, String>>();
			for (int i = 0; i < heteroAtomsToBeAdded.length; i++) {//populate allHeteroAtomInformation list
				String heteroAtomToBeAdded = heteroAtomsToBeAdded[i];
				String[] tempArray =matchSpace.split(heteroAtomToBeAdded);
				HashMap<String, String> heteroAtomInformation = new HashMap<String, String>();
				if (tempArray.length!=2){
					throw new PostProcessingException("malformed addHeteroAtom tag");
				}
				heteroAtomInformation.put("SMILES", tempArray[0]);
				if (tempArray[1].startsWith("id")){
					heteroAtomInformation.put("atomReferenceType", "id");
					heteroAtomInformation.put("atomReference", tempArray[1].substring(2));
				}
				else if (tempArray[1].startsWith("locant")){
					heteroAtomInformation.put("atomReferenceType", "locant");
					heteroAtomInformation.put("atomReference", tempArray[1].substring(6));
				}
				else{
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
						OpsinTools.setTextChild(previousEl, StringTools.objectListToString(locantValues, ","));
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
			for (int i = 0; i < bondsToBeAdded.length; i++) {//populate allBondInformation list
				String bondToBeAdded = bondsToBeAdded[i];
				String[] tempArray =matchSpace.split(bondToBeAdded);
				HashMap<String, String> bondInformation = new HashMap<String, String>();
				if (tempArray.length!=2){
					throw new PostProcessingException("malformed addBond tag");
				}
				bondInformation.put("bondOrder", tempArray[0]);
				if (tempArray[1].startsWith("id")){
					bondInformation.put("atomReferenceType", "id");
					bondInformation.put("atomReference", tempArray[1].substring(2));
				}
				else if (tempArray[1].startsWith("locant")){
					bondInformation.put("atomReferenceType", "locant");
					bondInformation.put("atomReference", tempArray[1].substring(6));
				}
				else{
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
						OpsinTools.setTextChild(previousEl, StringTools.objectListToString(locantValues, ","));
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
				state.fragManager.unsaturate(atomOnParentFrag.getID(), Integer.parseInt(bondInformation.get("bondOrder")) , parentFrag);
			}
		}
	}


	/**Converts locantgroups to individual locant tags, and checks for agreement
	 * between the number of locants, and multipliers. If there is disagreement the locant is checked against some special cases.
	 * If this fails an exception is thrown.
	 * @param state
	 * @param subOrBracketOrRoot The substituent/root/bracket to looks for locants in.
	 * @param root : used to check if a locant is referring to the root as in multiplicative nomenclature (root can be null for substituents)
	 * @throws PostProcessingException
	 * @throws StructureBuildingException
	 * @throws PostProcessingException If there is a disagreement.
	 */
	private void checkAndConvertToSingleLocants(BuildState state, Element subOrBracketOrRoot, Element root) throws StructureBuildingException, PostProcessingException {
		Elements locants = subOrBracketOrRoot.getChildElements("locant");
		Element group =subOrBracketOrRoot.getFirstChildElement("group");//will be null if element is a bracket
		for(int i=0;i<locants.size();i++) {
			Element locant = locants.get(i);
			String locantText = OpsinTools.removeDashIfPresent(locant.getValue());

			String [] locantValues = locantText.split(",");

			Element afterLocants = (Element)XOMTools.getNextSibling(locant);
			if(locantValues.length > 1) {
				while (afterLocants !=null){
					if(afterLocants.getLocalName().equals("multiplier") || afterLocants.getLocalName().equals("locant")) {
						break;
					}
					afterLocants = (Element)XOMTools.getNextSibling(afterLocants);
				}
				if(afterLocants != null && afterLocants.getLocalName().equals("multiplier")) {
					if(Integer.parseInt(afterLocants.getAttributeValue("value")) == locantValues.length ) {
						// number of locants and multiplier agree
						boolean locantModified =false;//did determineLocantMeaning do something?
						if (locantValues[locantValues.length-1].endsWith("'") && group!=null && subOrBracketOrRoot.indexOf(group) > subOrBracketOrRoot.indexOf(locant)){//quite possible that this is referring to a multiplied root

							if (group.getAttribute("outIDs")!=null){
								locantModified=determineLocantMeaning(state, locant, locantValues, root);
							}
							else{
								Element afterGroup = (Element)XOMTools.getNextSibling(group);
								int inlineSuffixCount =0;
								int multiplier=1;
								while (afterGroup !=null){
									if(afterGroup.getLocalName().equals("multiplier")){
										multiplier =Integer.parseInt(afterGroup.getAttributeValue("value"));
									}
									else if(afterGroup.getLocalName().equals("suffix") && afterGroup.getAttributeValue("type").equals("inline")){
										inlineSuffixCount +=(multiplier);
										multiplier=1;
									}
									afterGroup = (Element)XOMTools.getNextSibling(afterGroup);
								}
								if (inlineSuffixCount >=2){
									locantModified=determineLocantMeaning(state, locant, locantValues, root);
								}
							}
						}
						if (!locantModified && !XOMTools.getNextSibling(locant).equals(afterLocants)){//the locants apply indirectly the multiplier e.g. 2,3-butandiol
							//move the locant to be next to the multiplier.
							locant.detach();
							XOMTools.insertBefore(afterLocants, locant);
						}
					} else {
						if(!determineLocantMeaning(state, locant, locantValues, root)) throw new PostProcessingException("Mismatch between locant and multiplier counts (" +
								Integer.toString(locantValues.length) + " and " + afterLocants.getAttributeValue("value") + "):" + locant.toXML());
					}
				} else {
					/* Multiple locants without a multiplier */
					if(!determineLocantMeaning(state, locant, locantValues, root)) throw new PostProcessingException("Multiple locants without a multiplier: " + locant.toXML());
				}
			}

			//checks that determineLocantMeaning hasn't moved the locant or already assigned the locant a special meaning and detached it
			if (locant.getParent()!=null && locant.getParent().equals(subOrBracketOrRoot)){

				for(int j=0;j<locantValues.length;j++) {
					Element singleLocant = new Element("locant");
					String locantType = null;
					if (locant.getAttribute("type")!=null){
						locantType =locant.getAttributeValue("type");
						singleLocant.addAttribute(new Attribute("type", locantType));
					}
					String locantValue =locantValues[j];
					Matcher matches =matchCompoundLocant.matcher(locantValue);
					if (matches.find()==true){
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
	 * 	The locants could be intended to indicate the position of outIDs e.g. 1,4-phenylene
	 * 	The locants could be intended to indicate the attachement points of the root groups in multiplicative nomenclature e.g. 4,4'-methylenedibenzoic acid
	 * @param state
	 * @param locant The element corresponding to the locant group before the HW system.
	 * @param locantValues The locant values;
	 * @param root : used to check if a locant is referring to the root as in multiplicative nomenclature (root can be null for substituents)
	 * @return true if there's a HW system, and agreement; or if the locants conform to one of the alternative possibilities, otherwise false.
	 * @throws StructureBuildingException
	 */
	private boolean determineLocantMeaning(BuildState state, Element locant, String[] locantValues, Element root) throws StructureBuildingException {
		if (locant.getAttribute("type")!=null && locant.getAttributeValue("type").equals("multiplicativeNomenclature")) return true;//already known function (the locant must have been been previously moved by this method to an element that checkAndConvertToSingleLocants had not yet encountered)
		int count =locantValues.length;
		Element currentElem = (Element)XOMTools.getNextSibling(locant);
		int heteroCount = 0;
		while(currentElem != null && !currentElem.getLocalName().equals("group")){
			if(currentElem.getLocalName().equals("heteroatom")) {
				heteroCount++;
			} else if (currentElem.getLocalName().equals("multiplier")){
				heteroCount += Integer.parseInt(currentElem.getAttributeValue("value")) -1;
			}
			currentElem = (Element)XOMTools.getNextSibling(currentElem);
		}
		if(currentElem != null && currentElem.getLocalName().equals("group")){
			if (currentElem.getAttributeValue("subType").equals("hantzschWidman")) {
				if(heteroCount == count) {
					return true;
				} else {
					return false;//there is a case where locants don't apply to heteroatoms in a HW system, but in that case only one locant is expected so this function would not be called
				}
			}
			else if (heteroCount==0 && currentElem.getAttribute("outIDs")!=null) {//e.g. 1,4-phenylene
				String[] outIDs = matchComma.split(currentElem.getAttributeValue("outIDs"), -1);
				if (count ==outIDs.length && currentElem.getAttribute("substitutionRemovesOutIDs")==null){//things like oxy do not need to have their outIDs specified
					Fragment groupFragment =state.xmlFragmentMap.get(currentElem);
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
						currentElem.getAttribute("outIDs").setValue(StringTools.arrayToString(outIDs, ","));
						locant.detach();
						return true;
					}
				}
			}
			else if(currentElem.getValue().equals("benz") || currentElem.getValue().equals("benzo")){
				Node potentialGroupAfterBenzo = XOMTools.getNextSibling(currentElem, "group");//need to make sure this isn't benzyl
				if (potentialGroupAfterBenzo!=null){
					return true;//e.g. 1,2-benzothiazole
				}
			}
		}
		if (root!=null){
			Element multiplier =(Element) root.getChild(0);
			if (!multiplier.getLocalName().equals("multiplier") && ((Element)root.getParent()).getLocalName().equals("bracket")){//e.g. 1,1'-ethynediylbis(1-cyclopentanol)
				multiplier =(Element) root.getParent().getChild(0);
			}
			Node commonParent =locant.getParent().getParent();//this should be a common parent of the multiplier in front of the root. If it is not, then this locant is in a different scope
			Node parentOfMultiplier =multiplier.getParent();
			while (parentOfMultiplier!=null){
				if (commonParent.equals(parentOfMultiplier)){
					if (locantValues[count-1].endsWith("'")  &&
							multiplier.getLocalName().equals("multiplier") && multiplier.getAttribute("locantsAssigned")==null &&
							Integer.parseInt(multiplier.getAttributeValue("value")) == count ){//multiplicative nomenclature
						multiplier.addAttribute(new Attribute ("locantsAssigned",""));
						locant.detach();
						for(int i=locantValues.length-1; i>=0; i--) {
							Element singleLocant = new Element("multiplicativeLocant");
							singleLocant.addAttribute(new Attribute("value", locantValues[i]));
							XOMTools.insertAfter(multiplier, singleLocant);
						}
						return true;
					}
				}
				parentOfMultiplier=parentOfMultiplier.getParent();
			}
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
		Elements multipliers = elem.getChildElements("multiplier");
		for(int i=0;i<multipliers.size();i++) {
			Element m = multipliers.get(i);

			ArrayList<Element> locants = new ArrayList<Element>();
			Element possibleLocant =(Element)XOMTools.getPreviousSibling(m);
			while(possibleLocant!=null && possibleLocant.getLocalName().equals("locant")) {
				locants.add(possibleLocant);
				possibleLocant = (Element)XOMTools.getPreviousSibling(possibleLocant);
			}
			Element nextElem = (Element)XOMTools.getNextSibling(m);
			String nextName = nextElem.getLocalName();
			if(nextName.equals("unsaturator") ||
					nextName.equals("suffix") ||
					nextName.equals("heteroatom") ||
					nextName.equals("hydro")) {
				int mvalue = Integer.parseInt(m.getAttributeValue("value"));
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
						String locantValue =locant.getAttributeValue("value");
						//If a compound Locant e.g. 1(6) is detected add a compound locant attribute
						if (locant.getAttribute("compoundLocant")!=null){
							referent.addAttribute(new Attribute("compoundLocant", locant.getAttributeValue("compoundLocant")));
						}
						referent.addAttribute(new Attribute("locant", locantValue));
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


	/**
	 * Handles suffixes, passes them to resolveGroupAddingSuffixes.
	 * Processes the suffixAppliesTo command which multiplies a suffix and attaches the suffixes to the atoms described by the given IDs
	 * @param state
	 * @param subOrRoot
	 * @throws PostProcessingException
	 * @throws StructureBuildingException
	 */
	private void preliminaryProcessSuffixes(BuildState state, Element subOrRoot) throws PostProcessingException, StructureBuildingException{
		Elements groupsOfSubOrRoot = subOrRoot.getChildElements("group");
		Element lastGroupElementInSubOrRoot =groupsOfSubOrRoot.get(groupsOfSubOrRoot.size()-1);
		Fragment suffixableFragment =state.xmlFragmentMap.get(lastGroupElementInSubOrRoot);
		ArrayList<Element> suffixes =OpsinTools.elementsToElementArrayList(subOrRoot.getChildElements("suffix"));

		boolean imideSpecialCase =false;
		if (lastGroupElementInSubOrRoot.getAttribute("suffixAppliesTo")!=null){//trivial polyAcid
			//attribute contains instructions for number/positions of suffix
			//this is of the form comma sepeated ids with the number of ids corresponding to the number of instances of the suffix
			Element suffix =OpsinTools.getNextNonChargeSuffix(lastGroupElementInSubOrRoot);
			if (suffix ==null){
				throw new PostProcessingException("No suffix where suffix was expected");
			}
			if (suffixes.size()>1){
				throw new PostProcessingException("More than one suffix detected on trivial polyAcid. Not believed to be allowed");
			}
			String suffixInstruction =lastGroupElementInSubOrRoot.getAttributeValue("suffixAppliesTo");
			String[] suffixInstructions = matchComma.split(suffixInstruction);
			boolean symmetricSuffixes =true;
			if (suffix.getAttribute("additionalValue")!=null){//handles amic, aldehydic, anilic and amoyl suffixes properly
				if (suffixInstructions.length != 2){
					throw new PostProcessingException("suffix: " + suffix.getValue() + " used on an inappropriate group");
				}
				symmetricSuffixes = false;
				if (suffix.getValue().equals("imide")|| suffix.getValue().equals("imido") || suffix.getValue().equals("imidium")){
					imideSpecialCase =true;//prematurely resolve the two suffixes and explicitly join them to form a cyclic imide
				}
			}

			int firstIdInFragment=suffixableFragment.getIdOfFirstAtom();
			suffix.addAttribute(new Attribute("locantID", Integer.toString(firstIdInFragment + Integer.parseInt(suffixInstructions[0]) -1)));
			for (int i = 1; i < suffixInstructions.length; i++) {
				Element newSuffix = new Element("suffix");
				if (symmetricSuffixes){
					newSuffix.addAttribute(new Attribute("value", suffix.getAttributeValue("value")));
					newSuffix.addAttribute(new Attribute("type",  suffix.getAttributeValue("type")));
				}
				else{
					newSuffix.addAttribute(new Attribute("value", suffix.getAttributeValue("additionalValue")));
					newSuffix.addAttribute(new Attribute("type", "root"));
				}
				newSuffix.addAttribute(new Attribute("locantID", Integer.toString(firstIdInFragment + Integer.parseInt(suffixInstructions[i]) -1)));
				XOMTools.insertAfter(suffix, newSuffix);
				suffixes.add(newSuffix);
			}
		}
		else{
			for (Element suffix : suffixes) {
				if (suffix.getAttribute("additionalValue")!=null){
					throw new PostProcessingException("suffix: " + suffix.getValue() + " used on an inappropriate group");
				}
			}
		}

		ArrayList<Fragment> suffixFragments =resolveGroupAddingSuffixes(state, suffixes, suffixableFragment, lastGroupElementInSubOrRoot);
		processInfixFunctionalReplacementNomenclature(state, suffixes, suffixFragments);
		state.xmlSuffixMap.put(lastGroupElementInSubOrRoot, suffixFragments);
		
		if (imideSpecialCase){//Pretty horrible hack to allow cyclic imides
			if (suffixes.size() !=2){
				throw new PostProcessingException("Expected two suffixes fragments for cyclic imide");
			}
			Atom nitrogen =suffixFragments.get(0).getAtomList().get(1);//amide
			if (!nitrogen.getElement().equals("N")){
				throw new PostProcessingException("Nitrogen not found where nitrogen expected");
			}
			Atom carbon = suffixableFragment.getAtomByIDOrThrow(Integer.parseInt(suffixes.get(1).getAttributeValue("locantID")));
			if (!carbon.getElement().equals("C")){
				throw new PostProcessingException("Carbon not found where carbon expected");
			}
			resolveSuffixes(state, suffixableFragment, subOrRoot.getChildElements("suffix"), lastGroupElementInSubOrRoot);
			for (Element suffix : suffixes) {//suffixes have already been resolved so need to be detached to avoid being passed to resolveSuffixes later
				suffix.detach();
			}
			suffixableFragment.addBond(new Bond(nitrogen, carbon, 1));//join the N of the amide to the carbon of the acid to form the cyclic imide
		}
	}


	/**Processes a suffix and returns any fragment the suffix intends to add to the molecule
	 * @param state
	 * @param suffixes The suffix elements for a fragment.
	 * @param frag The fragment to which the suffix will be applied
	 * @param group The group in the XML
	 * @return An arrayList containing the generated fragments
	 * @throws StructureBuildingException If the suffixes can't be resolved properly.
	 * @throws PostProcessingException
	 */
	private ArrayList<Fragment> resolveGroupAddingSuffixes(BuildState state, ArrayList<Element> suffixes, Fragment frag, Element group) throws StructureBuildingException, PostProcessingException {
		ArrayList<Fragment> suffixFragments =new ArrayList<Fragment>();
		String groupType = frag.getType();
		String subgroupType = frag.getSubType();

		String suffixTypeToUse =null;
		boolean suffixTypeDetermined =false;
		if (suffixApplicability.containsKey(groupType)){
			suffixTypeToUse =groupType;
			suffixTypeDetermined=true;
		}
		else if (groupType.equals("simpleGroup")){
			suffixTypeToUse="substituent";
			suffixTypeDetermined=true;
		}

		//if suffixTypeToUse is still null then it is type cyclic or acyclic as determined by the atom property atomIsInACycle
		for(int i=0;i<suffixes.size();i++) {
			Element suffix = suffixes.get(i);
			String suffixValue = suffix.getAttributeValue("value");

			if (!suffixTypeDetermined){
				boolean cyclic;
				if (suffix.getAttribute("locant")!=null){
					Atom a =frag.getAtomByLocant(suffix.getAttributeValue("locant"));
					if (a!=null){
						cyclic=a.getAtomIsInACycle();
					}
					else{//can happen in the cases of things like fused rings where the final numbering is not available (in which case all the atoms will be cyclic anyway)
						cyclic = frag.getAtomByIDOrThrow(frag.getIdOfFirstAtom()).getAtomIsInACycle();
					}
				}
				else{
					cyclic = frag.getAtomByIDOrThrow(frag.getIdOfFirstAtom()).getAtomIsInACycle();
				}
				if (cyclic){
					suffixTypeToUse="cyclic";
				}
				else{
					suffixTypeToUse="acyclic";
				}
			}
			Elements suffixRuleTags =getSuffixRuleTags(suffixTypeToUse, suffixValue, subgroupType);

			if (suffixTypeToUse.equals("acidStem")){//special IUPAC rules for acids
				if (suffix.getAttributeValue("type").equals("inline")){//handles cases such as carbonyl which are diradicals and have lost a hydroxy from their acid description
					removeHydroxyGroupsAndAddOutIDs(state, group, frag.getAtomByIDOrThrow(frag.getIdOfFirstAtom()));
				}
				else if (suffix.getAttributeValue("value").equals("ate")|| suffix.getAttributeValue("value").equals("ite")){//handles cases such as phosphonate and carbonate that have multiple functional IDs and O- atoms
					chargeHydroxyGroupsAndAddOutIDs(frag, frag.getAtomByIDOrThrow(frag.getIdOfFirstAtom()));
				}
			}
			Fragment suffixFrag =null;
			/*
			 * Temp fragments are build for each addGroup rule and then merged into suffixFrag
			 */
			for(int j=0;j<suffixRuleTags.size();j++) {
				Element suffixRuleTag = suffixRuleTags.get(j);
				if(suffixRuleTag.getLocalName().equals("addgroup")) {
					Fragment tempSuffixFrag =null;
					String bondOrderStr = suffixRuleTag.getAttributeValue("bondOrder");
					int bondOrder = 1;
					if(bondOrderStr != null) bondOrder = Integer.parseInt(bondOrderStr);
					String labels="none";
					if (suffixRuleTag.getAttribute("labels")!=null){
						labels =suffixRuleTag.getAttributeValue("labels");
					}

					if(suffixRuleTag.getAttribute("setsOutID") != null) {
						tempSuffixFrag= state.fragManager.buildSMILES(suffixRuleTag.getAttributeValue("SMILES"), "suffix", "outSuffix", labels);
						if(tempSuffixFrag.getOutIDs().size() == 0) {
							if(suffixRuleTag.getAttribute("outValency") != null) {
								tempSuffixFrag.addOutID(tempSuffixFrag.getIdOfFirstAtom(), Integer.parseInt(suffixRuleTag.getAttributeValue("outValency")), true);
							}
							else{
								tempSuffixFrag.addOutID(tempSuffixFrag.getIdOfFirstAtom(), 1, true);
							}
						}
					}
					else if(suffixRuleTag.getAttribute("setsDefaultInID") != null) {
						tempSuffixFrag= state.fragManager.buildSMILES(suffixRuleTag.getAttributeValue("SMILES"), "suffix", "inSuffix", labels);
					}
					else if	(suffixRuleTag.getAttribute("setsFunctionalID") != null) {
						tempSuffixFrag= state.fragManager.buildSMILES(suffixRuleTag.getAttributeValue("SMILES"), "suffix", "functionalSuffix", labels);
						if(tempSuffixFrag.getFunctionalIDs().size() == 0) {
							tempSuffixFrag.addFunctionalID(tempSuffixFrag.getIdOfFirstAtom());
						}
					}
					else{
						tempSuffixFrag= state.fragManager.buildSMILES(suffixRuleTag.getAttributeValue("SMILES"), "suffix", "suffix", labels);
					}


					if (suffixFrag==null){
						suffixFrag=tempSuffixFrag;
					}
					else{
						suffixFrag.importFrag(tempSuffixFrag);
						if(suffixRuleTag.getAttribute("setsOutID") != null) {
							suffixFrag.addOutIDs(tempSuffixFrag.getOutIDs());
						}
						else if(suffixRuleTag.getAttribute("setsDefaultInID") != null) {
							suffixFrag.setDefaultInID(tempSuffixFrag.getDefaultInID());
						}
						else if	(suffixRuleTag.getAttribute("setsFunctionalID") != null) {
							suffixFrag.addFunctionalIDs(tempSuffixFrag.getFunctionalIDs());
						}
						state.fragManager.removeFragment(tempSuffixFrag);
					}
					suffixFrag.addInID(tempSuffixFrag.getDefaultInID(), bondOrder);
				}
			}
			if (suffixFrag!=null){
				suffixFragments.add(suffixFrag);
				state.xmlFragmentMap.put(suffix,suffixFrag);
			}
		}
		return suffixFragments;
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
		if(potentiallyApplicableSuffixes==null || potentiallyApplicableSuffixes.size()==0 ) throw new PostProcessingException("Suffix: " +suffixValue +" does not apply to the group it was associated with (type: "+  suffixTypeToUse + ")according to suffixApplicability.xml");
		Element chosenSuffix=null;
		for (int i = 0; i < potentiallyApplicableSuffixes.size(); i++) {
			Element suffix =potentiallyApplicableSuffixes.get(i);
			if (suffix.getAttribute("subType")!=null){
				if (!suffix.getAttributeValue("subType").equals(subgroupType)){
					continue;
				}
			}
			if (chosenSuffix!=null){
				throw new PostProcessingException("Suffix: " +suffixValue +" appears multiple times in suffixApplicability.xml");
			}
			chosenSuffix=suffix;
		}
		if (chosenSuffix==null){
			throw new PostProcessingException("Suffix: " +suffixValue +" does not apply to the group it was associated with (type: "+  suffixTypeToUse + ")due to the group's subType: "+ subgroupType +" according to suffixApplicability.xml");
		}
		Element rule =suffixRules.get(chosenSuffix.getValue());
		if(rule ==null) throw new PostProcessingException("Suffix: " +chosenSuffix.getValue() +" does not have a rule associated with it in suffixRules.xml");
		return rule.getChildElements();
	}


	/**
	 * Removes surrounding hydroxy groups from given atom and adds that many optional outIDs to the group
	 * @param group
	 * @param atom
	 * @throws StructureBuildingException
	 */
	private void removeHydroxyGroupsAndAddOutIDs(BuildState state, Element group, Atom atom) throws StructureBuildingException {
		List<Atom> neighbours = atom.getAtomNeighbours();
		Fragment frag =atom.getFrag();
		int atomsRemoved =0;
		for (Atom neighbour : neighbours) {
			if (neighbour.getElement().equals("O") && neighbour.getCharge()==0 && neighbour.getAtomNeighbours().size()==1 && neighbour.getFrag().equals(frag)){
				if (frag.findBond(atom, neighbour).getOrder()==1){
					frag.removeAtom(neighbour, state.fragManager);
					atomsRemoved++;
				}
			}
		}
		if (atomsRemoved >0){
			if (group.getAttribute("outIDs")!=null){
				throw new StructureBuildingException("Group should not already have outIDs");
			}
			String outIDs =StringTools.multiplyString(atom.getID() -frag.getIdOfFirstAtom() +1 +",", atomsRemoved);//add an id relative to this frag
			outIDs =outIDs.substring(0, outIDs.length()-1);//remove extra comma
			group.addAttribute(new Attribute("outIDs", outIDs));
			if (group.getAttribute("substitutionRemovesOutIDs")==null){
				group.addAttribute(new Attribute("substitutionRemovesOutIDs", "yes"));//c.f. carbonyl
			}
		}
	}


	/**
	 * Adds a negative to all hydroxy groups surrounding the given atom and adds that many functional outIDs to the fragment
	 * @param frag
	 * @param atom
	 * @throws StructureBuildingException
	 */
	private void chargeHydroxyGroupsAndAddOutIDs(Fragment frag, Atom atom) throws StructureBuildingException {
		List<Atom> neighbours = atom.getAtomNeighbours();
		for (Atom neighbour : neighbours) {
			if (neighbour.getElement().equals("O") && neighbour.getCharge()==0 && neighbour.getAtomNeighbours().size()==1 && neighbour.getFrag().equals(frag)){
				if (frag.findBond(atom, neighbour).getOrder()==1){
					neighbour.setCharge(-1);
					frag.addFunctionalID(neighbour.getID());
				}
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
	private void processInfixFunctionalReplacementNomenclature(BuildState state, ArrayList<Element> suffixes, ArrayList<Fragment> suffixFragments) throws StructureBuildingException, PostProcessingException {
		for (int i = 0; i < suffixes.size(); i++) {
			Element suffix = suffixes.get(i);
			if (suffix.getAttribute("infix")!=null){
				Fragment suffixFrag = state.xmlFragmentMap.get(suffix);
				if (suffixFrag==null){
					throw new PostProcessingException("infix has erroneously been assigned to a suffix which does not correspond to a suffix fragment. suffix: " + suffix.getValue());
				}
				String replacementSMILES = suffix.getAttributeValue("infix");
				int infixCount = Integer.parseInt(suffix.getAttributeValue("infixCount"));
				List<Atom> atomList =suffixFrag.getAtomList();
				LinkedList<Atom> singleBondedOxygen =new LinkedList<Atom>();
				LinkedList<Atom> doubleBondedOxygen =new LinkedList<Atom>();
				for (Atom a : atomList) {
					if (a.getElement().equals("O")){//find terminal oxygens
						if ((a.getBonds().size()==1 && a.getOutValency()==0)|| (a.getBonds().size()==0)){
							int incomingValency =a.getOutValency() +a.getIncomingValency();
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
				Matcher m =matchUpperCase.matcher(replacementSMILES);
				int atomCount=0;
				while(m.find()) {
					atomCount++;//assumption made that number of upper case letter = number of atoms the SMILES is describing
				}
				int oxygenAvailable = singleBondedOxygen.size() +doubleBondedOxygen.size();

				/* This block handles infix multiplication. Unless brackets are provided this is ambiguous without knowledge of the suffix that is being modified
				 * For example butandithione could be intepreted as butandi(thione) or butan(dithi)one. Obviously the latter is wrong in this case
				 * but it is the correct interpretation for butandithiate
				 */
				Element possibleInfix =(Element) XOMTools.getPreviousSibling(suffix);
				if (possibleInfix.getLocalName().equals("infix")){//the infix is only left when there was ambiguity
					Element possibleMultiplier =(Element) XOMTools.getPreviousSibling(possibleInfix);
					if (possibleMultiplier.getLocalName().equals("multiplier")){
						int multiplierValue =Integer.parseInt(possibleMultiplier.getAttributeValue("value"));
						if (infixCount*multiplierValue <=oxygenAvailable){//multiplier means multiply the infix e.g. butandithiate
							infixCount =infixCount*multiplierValue;
							suffix.getAttribute("infixCount").setValue(Integer.toString(infixCount));
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
							    suffix.addAttribute(new Attribute("locant", locant.getAttributeValue("value")));
								suffix.addAttribute(new Attribute("multiplied", "multiplied"));
								locant.detach();
							}
							for (int j = 1; j < multiplierValue; j++) {//multiplier means multiply the infixed suffix e.g. butandithione
								Element newSuffix =new Element(suffix);
								Fragment newSuffixFrag =state.fragManager.copyAndRelabel(suffixFrag);
								state.xmlFragmentMap.put(newSuffix, newSuffixFrag);
								suffixFragments.add(newSuffixFrag);
								XOMTools.insertAfter(suffix, newSuffix);
								suffixes.add(newSuffix);
								if (locants.size()>0){
									Element locant =locants.removeFirst();
									newSuffix.getAttribute("locant").setValue(locant.getAttributeValue("value"));
									locant.detach();
								}
							}
						}
						possibleMultiplier.detach();
						possibleInfix.detach();
					}
					else{
						throw new PostProcessingException("Multiplier expected in front of infix");
					}
				}

				if (atomCount>1){//something like peroxy
					if (infixCount > singleBondedOxygen.size()){
						throw new StructureBuildingException("Cannot find single bonded oxygen for infix with SMILES: "+ replacementSMILES+ " to modify!");
					}
				}
				else{
					if (infixCount > oxygenAvailable){
						throw new StructureBuildingException("Cannot find oxygen for infix with SMILES: "+ replacementSMILES+ " to modify!");
					}
				}

				boolean infixAssignmentAmbiguous =false;
				ArrayList<Atom> ambiguousElementAtoms = new ArrayList<Atom>();
				if (infixCount != oxygenAvailable){
					infixAssignmentAmbiguous=true;
				}
				for (int j = 0; j < infixCount; j++) {
					if (atomCount>1){//something like peroxy
						state.fragManager.replaceTerminalAtomWithFragment(singleBondedOxygen.get(0), state.fragManager.buildSMILES(replacementSMILES, "suffix", "none").getFirstAtom());
					}
					else{
						Atom a;
						if (doubleBondedOxygen.size()>0){
							a=doubleBondedOxygen.removeFirst();
						}
						else{
							a=singleBondedOxygen.removeFirst();
						}
						state.fragManager.makeHeteroatom(a, replacementSMILES, false);
						if (infixAssignmentAmbiguous){
							ambiguousElementAtoms.add(a);
						}
					}
				}
				if (infixAssignmentAmbiguous){//record what atoms could have been replaced. Often this ambiguity is resolved later e.g. S-methyl ethanthioate
					for (Atom a : doubleBondedOxygen) {
						ambiguousElementAtoms.add(a);
					}
					for (Atom a : singleBondedOxygen) {
						ambiguousElementAtoms.add(a);
					}
					String atomIDsString="";
					for (Atom atom : ambiguousElementAtoms) {
						atomIDsString+=atom.getID();
						atomIDsString+=",";
					}
					atomIDsString= atomIDsString.substring(0, atomIDsString.length()-1);
					for (Atom atom : ambiguousElementAtoms) {
						atom.setNote("ambiguousElementAssignment", atomIDsString);
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
	 * @throws StructureBuildingException
	 */
	private void matchLocantsToDirectFeatures(Element subOrRoot) throws PostProcessingException, StructureBuildingException {
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
							delta.appendChild(locant.getAttributeValue("value"));
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
						if (locantsBeforeHWSystem.size() ==1 && Integer.parseInt(group.getAttributeValue("value")) <=10){
							locants.remove(locantsBeforeHWSystem.get(0));//don't assign this locant
						}
						else if (locantsBeforeHWSystem.size() == heteroAtoms.size()){//general case
							for (int j = 0; j < locantsBeforeHWSystem.size(); j++) {
								Element locant =locantsBeforeHWSystem.get(j);
								heteroAtoms.get(j).addAttribute(new Attribute("locant", locant.getAttributeValue("value")));
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

		for (Element locant : locants) {
			String locantValue =locant.getAttributeValue("value");
			if (matchElementSymbol.matcher(locantValue).matches()){//element symbol locant
				continue;
			}
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
					refName.equals("hydro"))) {

				//If a compound Locant e.g. 1(6) is detected add a compound locant attribute
				if (locant.getAttribute("compoundLocant")!=null){
					referent.addAttribute(new Attribute("compoundLocant", locant.getAttributeValue("compoundLocant")));
				}
				referent.addAttribute(new Attribute("locant", locantValue));
				locant.detach();
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
	 * @throws PostProcessingException 
	 */
	private boolean processPrefixFunctionalReplacementNomenclature(BuildState state, List<Element> groups, List<Element> substituents) throws StructureBuildingException, PostProcessingException {
		boolean doneSomething =false;
		for (Element group : groups) {
			if (matchChalogenReplacment.matcher(group.getValue()).matches()){
				//need to check whether this is an instance of functional replacement by checking the substituent/root it is applying to
				Element substituent =(Element) group.getParent();
				Element firstElInSubstituent=(Element)substituent.getChild(0);
				Element nextSubOrBracket = (Element) XOMTools.getNextSibling(substituent);
				if (nextSubOrBracket!=null && (nextSubOrBracket.getLocalName().equals("substituent") ||nextSubOrBracket.getLocalName().equals("root"))){
					if (nextSubOrBracket.getLocalName().equals("root") && firstElInSubstituent.getLocalName().equals("locant")){
						continue;//the substituent appears to be a prefix and not a case of functional replacement as it appears to be locanted on this root
					}
					Element groupToBeModified =(Element) nextSubOrBracket.getChild(0);
					String nameOfGroupToBeModified =groupToBeModified.getValue();
					if (XOMTools.getNextSibling(group)!=null || !groupToBeModified.getLocalName().equals("group")){
						continue;//there should be no gap between the functional replacement and group c.f. not 2,2'-thiodipyran
					}
					int numberOfAtomsToReplace =1;
					if (firstElInSubstituent.getLocalName().equals("multiplier")){
						numberOfAtomsToReplace =Integer.valueOf(firstElInSubstituent.getAttributeValue("value"));
					}

					String replacementSMILES = group.getAttributeValue("value");
					Matcher m =matchUpperCase.matcher(replacementSMILES);
					int atomCount=0;
					while(m.find()) {
						atomCount++;//assumption made that number of upper case letter = number of atoms the SMILES is describing
					}
					Fragment frag = state.xmlFragmentMap.get(groupToBeModified);
					ArrayList<Fragment> suffixes = state.xmlSuffixMap.get(groupToBeModified);
					ArrayList<Atom> replaceableAtoms =new ArrayList<Atom>();

					if (atomCount==1){
						for (Atom atom : frag.getAtomList()) {
							if (atom.getElement().equals("O")){
								replaceableAtoms.add(atom);
							}
						}
						//suffixes are ignored in the case of fused ring systems due to suffixes being null and for all suffixes not on acid stems (except phen-->thiophenol)
						if (suffixes!=null){
							ArrayList<Fragment> applicableSuffixes = new ArrayList<Fragment>(suffixes);
							if (!groupToBeModified.getAttributeValue("type").equals("acidStem") && !nameOfGroupToBeModified.equals("phen")){
								//remove all non acid suffixes
								for (Fragment fragment : suffixes) {
									Element suffix = state.xmlFragmentMap.getElement(fragment);
									if (!suffix.getAttributeValue("value").equals("ic") && !suffix.getAttributeValue("value").equals("ous")){
										applicableSuffixes.remove(fragment);
									}
								}
							}
							boolean endOfSuffixesFound=false;
							int i=0;
							while (!endOfSuffixesFound){//makes oxygen atoms alternate between suffixes c.f. dithioterephthalic acid
								endOfSuffixesFound =true;
								for (Fragment fragment : applicableSuffixes) {
									if (i < fragment.getAtomList().size()){
										Atom a =fragment.getAtomList().get(i);
										if (a.getElement().equals("O")){
											replaceableAtoms.add(a);
										}
										endOfSuffixesFound=false;
									}
								}
								i++;
							}
						}
					}
					else{
						//only consider ic suffixes when it's not simple atom replacement e.g. peroxy.
						ArrayList<Element> suffixElements =OpsinTools.getLaterSiblingElementsOfType(groupToBeModified, "suffix");
						for (int i = 0; i < suffixElements.size(); i++) {
							if (suffixElements.get(i).getAttributeValue("value").equals("ic")||suffixElements.get(i).getAttributeValue("value").equals("ous")){
								Fragment suffixFrag =suffixes.get(i);
								List<Atom> atomList =suffixFrag.getAtomList();
								for (Atom a : atomList) {
									if (a.getElement().equals("O")){
								    	int totalBondOrder =a.getIncomingValency();
								    	totalBondOrder += a.getOutValency();//take into account the bond which hasn't been added yet to a main fragment
										if (totalBondOrder==1 || (atomCount==1 && totalBondOrder==2)){
											replaceableAtoms.add(a);
										}
									}
								}
							}
						}
					}
					if (replaceableAtoms.size() >=numberOfAtomsToReplace){//check that there atleast as many oxygens as requested replacements
						boolean prefixAssignmentAmbiguous =false;
						ArrayList<Atom> ambiguousElementAtoms = new ArrayList<Atom>();
						if (replaceableAtoms.size() != numberOfAtomsToReplace){
							prefixAssignmentAmbiguous=true;
						}

						for (int i = 0; i < numberOfAtomsToReplace; i++) {
							Atom a =replaceableAtoms.get(i);
							if (atomCount>1){//something like peroxy
								state.fragManager.replaceTerminalAtomWithFragment(a, state.fragManager.buildSMILES(replacementSMILES, "suffix", "none").getFirstAtom());
							}
							else{
								state.fragManager.makeHeteroatom(a, replacementSMILES, false);
								if (prefixAssignmentAmbiguous){
									ambiguousElementAtoms.add(a);
								}
							}
						}
						
						if (prefixAssignmentAmbiguous){//record what atoms could have been replaced. Often this ambiguity is resolved later e.g. S-methyl thioacetate
							for (Atom a : replaceableAtoms) {
								ambiguousElementAtoms.add(a);
							}
							String atomIDsString="";
							for (Atom atom : ambiguousElementAtoms) {
								atomIDsString+=atom.getID();
								atomIDsString+=",";
							}
							atomIDsString= atomIDsString.substring(0, atomIDsString.length()-1);
							for (Atom atom : ambiguousElementAtoms) {
								atom.setNote("ambiguousElementAssignment", atomIDsString);
							}
						}
						
						state.fragManager.removeFragment(state.xmlFragmentMap.get(group));
						substituent.removeChild(group);
						Elements remainingChildren =substituent.getChildElements();//there may be a locant that should be moved
						for (int i = remainingChildren.size()-1; i>=0; i--){
							Node child =substituent.getChild(i);
							child.detach();
							nextSubOrBracket.appendChild(child);
						}
						substituents.remove(substituent);
						substituent.detach();
						if (firstElInSubstituent.getLocalName().equals("multiplier")){
							firstElInSubstituent.detach();
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
			List<Element> siblings = OpsinTools.findChildElementsWithTagNames(parent, new String[]{"bracket", "substituent", "root"});
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
					Element group = bracketOrSub.getFirstChildElement("group");
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

		if (foundSibling ==false){//Special case: anything the group could potentially substitute onto is in a bracket. The bracket is checked recursively
			s = new Stack<Element>();
			s.add(startingElement);
			doneFirstIteration =false;//check on index only done on first iteration to only get elements with an index greater than the starting element
			while (s.size()>0){
				Element currentElement =s.pop();
				Element parent = (Element)currentElement.getParent();
				List<Element> siblings = OpsinTools.findChildElementsWithTagNames(parent, new String[]{"bracket", "substituent", "root"});
				int indexOfCurrentElement =parent.indexOf(currentElement);

				for (Element bracketOrSub : siblings) {
					if (!doneFirstIteration && parent.indexOf(bracketOrSub) <= indexOfCurrentElement){
						continue;
					}
					if (bracketOrSub.getLocalName().equals("bracket")){
						s.push((Element)bracketOrSub.getChild(0));
					}
					else{
						Element group = bracketOrSub.getFirstChildElement("group");
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
	 * Assigns Element symbols to groups and suffixes.
	 * Suffixes have preference.
	 * @param state
	 * @param subOrRoot
	 */
	private void assignElementSymbolLocants(BuildState state, Element subOrRoot) {
		Elements groupsOfSubOrRoot = subOrRoot.getChildElements("group");
		Element lastGroupElementInSubOrRoot =groupsOfSubOrRoot.get(groupsOfSubOrRoot.size()-1);
		ArrayList<Fragment> suffixFragments =state.xmlSuffixMap.get(lastGroupElementInSubOrRoot);
		Fragment suffixableFragment =state.xmlFragmentMap.get(lastGroupElementInSubOrRoot);
		FragmentManager.assignElementLocants(suffixableFragment, suffixFragments);
	}


	/**
	 * Handles Hantzsch-Widman rings. Adds SMILES to the group corresponding to the ring's structure
	 * @param state
	 * @param elem
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
			int ringSize = Integer.parseInt(group.getAttributeValue("value"));
			Element prev = (Element) XOMTools.getPreviousSibling(group);
			ArrayList<Element> prevs = new ArrayList<Element>();
			boolean noLocants = true;
			while(prev != null && prev.getLocalName().equals("heteroatom")) {
				prevs.add(prev);
				if(prev.getAttribute("locant") != null) {
					noLocants = false;
				}
				prev = (Element) XOMTools.getPreviousSibling(prev);
			}
			boolean hasNitrogen = false;
			boolean hasSiorGeorSborPb=false;
			for(Element heteroatom : prevs){
				String heteroAtomElement =heteroatom.getAttributeValue("value");
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
			if (ringSize == 6 && ringType.equals("ring") && hasNitrogen ==false && hasSiorGeorSborPb ==true && (group.getValue().equals("in") ||group.getValue().equals("an"))){
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
							a.setSpareValency(0);
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
					String elementReplacement =heteroatom.getAttributeValue("value");
					if (elementReplacement.startsWith("[") && elementReplacement.endsWith("]")){
						elementReplacement=elementReplacement.substring(1, elementReplacement.length()-1);
					}
					Atom a =hwRing.getAtomByLocantOrThrow(locant);
					a.setElement(elementReplacement);
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
				String elementReplacement =heteroatom.getAttributeValue("value");
				if (elementReplacement.startsWith("[") && elementReplacement.endsWith("]")){
					elementReplacement=elementReplacement.substring(1, elementReplacement.length()-1);
				}

				while (!hwRing.getAtomByLocantOrThrow(Integer.toString(defaultLocant)).getElement().equals("C")){
					defaultLocant++;
				}
				Atom a =hwRing.getAtomByLocantOrThrow(Integer.toString(defaultLocant));
				a.setElement(elementReplacement);
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
					while (firstInDoubleBond.getSpareValency() != 0 || ValencyChecker.checkValencyAvailableForBond(firstInDoubleBond, 1) != true ||
							secondInDoubleBond.getSpareValency() != 0 || ValencyChecker.checkValencyAvailableForBond(secondInDoubleBond, 1) != true){
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
			OpsinTools.setTextChild(group, name);
		}
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
			int mvalue = Integer.parseInt(multiplier.getAttributeValue("value"));

			/*
			 * Populate locants with locants. Two locants are required for every pair of rings to be joined.
			 * e.g. bi requires 2, ter requires 4 etc.
			 */
			ArrayList<List<String>> ringJoiningLocants =new ArrayList<List<String>>();
			Element previousEl =(Element)XOMTools.getPreviousSibling(multiplier);
			Element group =(Element)XOMTools.getNextSibling(multiplier, "group");
			if (previousEl!=null && previousEl.getLocalName().equals("ringAssemblyLocant")){//a locant appears to have provided to indicate how to connect the rings of the ringAssembly
				String locantText =OpsinTools.removeDashIfPresent(previousEl.getValue());
				//special cases where often locants are meant to apply to suffixes rather than being a description of where the rings connect to each other
				if (group.getValue().equals("phen") || group.getValue().equals("hex") || group.getValue().equals("benz")){
					//Find elements that can have locants but don't currently
					List<Element> locantAble = OpsinTools.findChildElementsWithTagNames(subOrRoot, new String[]{"suffix", "unsaturator", "heteroatom", "hydro"});
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
				if (previousEl.getAttribute("type")!=null && previousEl.getAttributeValue("type").equals("orthoMetaPara")){//an OMP locant appears to have provided to indicate how to connect the rings of the ringAssembly
					String locant2 =previousEl.getAttributeValue("value");
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
							(inlineSuffixSeen == 0 && currentEl.getLocalName().equals("suffix") && currentEl.getAttributeValue("type").equals("inline") && currentEl.getAttribute("locant")==null)||
							(currentEl.getLocalName().equals("suffix") && currentEl.getAttributeValue("type").equals("charge"))){
						currentEl.detach();
						elementToResolve.appendChild(currentEl);
					}
					else{
						break;
					}
					if (currentEl.getLocalName().equals("group")){
						groupFound = 1;
					}
					if ((currentEl.getLocalName().equals("suffix") && currentEl.getAttributeValue("type").equals("inline"))){
						inlineSuffixSeen = 1;
					}
				}
			}

			Elements suffixes =elementToResolve.getChildElements("suffix");
			Fragment fragmentToResolveAndDuplicate =state.xmlFragmentMap.get(group);
			resolveSuffixes(state, fragmentToResolveAndDuplicate, suffixes, group);
			StructureBuildingMethods.resolveLocantedFeatures(state, elementToResolve);
			StructureBuildingMethods.resolveUnLocantedFeatures(state, elementToResolve);
			group.detach();
			XOMTools.insertAfter(multiplier, group);

			int bondOrder = 1;
			if (fragmentToResolveAndDuplicate.getOutIDs().size()>0){//e.g. bicyclohexanylidene
				bondOrder =fragmentToResolveAndDuplicate.getOutID(0).valency;
				fragmentToResolveAndDuplicate.removeOutID(0);
			}
			if (fragmentToResolveAndDuplicate.getOutIDs().size()>0){
				throw new StructureBuildingException("Ring assembly fragment should have one or no OutIDs; not more than one!");
			}

			ArrayList<Fragment> clonedFragments = new ArrayList<Fragment>();
			for (int j = 1; j < mvalue; j++) {
				clonedFragments.add(state.fragManager.copyAndRelabel(fragmentToResolveAndDuplicate, StringTools.multiplyString("'", j)));
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
					atomOnParent =fragmentToResolveAndDuplicate.getAtomByIdOrNextSuitableAtomOrThrow(fragmentToResolveAndDuplicate.getDefaultInID(), bondOrder);
					atomOnLatestClone = clone.getAtomByIdOrNextSuitableAtomOrThrow(clone.getDefaultInID(), bondOrder);
				}
				state.fragManager.incorporateFragment(clone, atomOnLatestClone.getID(), fragmentToResolveAndDuplicate, atomOnParent.getID(), bondOrder);
			}
			OpsinTools.setTextChild(group, multiplier.getValue() +group.getValue());
			multiplier.detach();
		}
	}


	private void processComplicatedSpiroNomenclature(BuildState state, Element subOrRoot) {
		// TODO Auto-generated method stub

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

		for (Element theSubstituent : substituents) {
			String firstElInSubName =((Element)theSubstituent.getChild(0)).getLocalName();
			if (firstElInSubName.equals("locant") ||firstElInSubName.equals("multiplier")){
				continue;
			}

			Element theSubstituentGroup = theSubstituent.getFirstChildElement("group");
			String theSubstituentSubType = theSubstituentGroup.getAttributeValue("subType");
			String theSubstituentType = theSubstituentGroup.getAttributeValue("type");

			//Only some substituents are valid joiners (e.g. no rings are valid joiners). Need to be atleast bivalent
			if (theSubstituentGroup.getAttribute("usableAsAJoiner")==null){
				continue;
			}
			Fragment theSubstituentGroupFragment =state.xmlFragmentMap.get(theSubstituentGroup);

			//there must be an element after the substituent for the implicit bracket to be required
			Element elementAftersubstituent =(Element)XOMTools.getNextSibling(theSubstituent);
			if (elementAftersubstituent ==null ||
					!elementAftersubstituent.getLocalName().equals("substituent") &&
					!elementAftersubstituent.getLocalName().equals("bracket") &&
					!elementAftersubstituent.getLocalName().equals("root")){
				continue;
			}

			//checks that the element before is a substituent or a bracket which will obviously include substituent/s
			//this makes sure there's more than just a substituent in the bracket
			Element elementBeforeSubstituent =(Element)XOMTools.getPreviousSibling(theSubstituent);
			if (elementBeforeSubstituent ==null||
					!elementBeforeSubstituent.getLocalName().equals("substituent") &&
					!elementBeforeSubstituent.getLocalName().equals("bracket")){
				continue;
			}

			//look for hyphen between substituents, this seems to indicate implicit bracketing was not desired e.g. dimethylaminomethane vs dimethyl-aminomethane
			Element elementDirectlyBeforeSubstituent = (Element) OpsinTools.getPrevious(theSubstituent.getChild(0));//can't return null as we know elementBeforeSubstituent is not null
			if (elementDirectlyBeforeSubstituent.getLocalName().equals("hyphen")){
				continue;
			}

			//prevents alkyl chains being bracketed together e.g. ethylmethylamine
			//...unless it's something like 2-methylethyl where the first appears to be locanted onto the second
			List<Element> groupElements  = OpsinTools.findDescendantElementsWithTagName(elementBeforeSubstituent, "group");//one for a substituent, possibly more for a bracket
			Element group =groupElements.get(groupElements.size()-1);
			if (group==null){throw new PostProcessingException("No group where group was expected");}
			if (theSubstituentType.equals("chain") && theSubstituentSubType.equals("alkaneStem") &&
					group.getAttributeValue("type").equals("chain") && (group.getAttributeValue("subType").equals("alkaneStem") || group.getAttributeValue("subType").equals("alkaneStem-irregular"))){
				boolean placeInImplicitBracket =false;

				Element suffixAfterGroup=(Element)XOMTools.getNextSibling(group, "suffix");
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
						if (currentElementName.equals("stereoChemistry")){
						}
						else if (currentElementName.equals("locant")){
							String locantText =childrenOfElementBeforeSubstituent.get(i).getAttributeValue("value");
							if(!theSubstituentGroupFragment.hasLocant(locantText)){
								foundLocantNotReferringToChain=true;
								break;
							}
							else{
								foundLocantNotReferringToChain=false;
							}
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
			
			

			Element twoElementsBeforeSubstituent =(Element)XOMTools.getPreviousSibling(group.getParent());
			//prevent bracketing to multi radicals
			//This is really tacky...
			//TODO do this properly...
			if (group.getAttribute("outIDs")!=null && (group.getAttribute("usableAsAJoiner")==null || twoElementsBeforeSubstituent==null )){
				continue;
			}
			Element afterGroup = (Element)XOMTools.getNextSibling(group);
			int inlineSuffixCount =0;
			int multiplier=1;
			while (afterGroup !=null){
				if(afterGroup.getLocalName().equals("multiplier")){
					multiplier =Integer.parseInt(afterGroup.getAttributeValue("value"));
				}
				else if(afterGroup.getLocalName().equals("suffix") && afterGroup.getAttributeValue("type").equals("inline")){
					inlineSuffixCount +=(multiplier);
					multiplier=1;
				}
				afterGroup = (Element)XOMTools.getNextSibling(afterGroup);
			}
			if (inlineSuffixCount >=2){
				continue;
			}

			Element bracket = new Element("bracket");
			bracket.addAttribute(new Attribute("type", "implicit"));

			/*
			 * locant may need to be moved. This occurs when the group in elementBeforeSubstituent is not supposed to be locanted onto
			 *  theSubstituentGroup
			 *  e.g. 2-aminomethyl-1-chlorobenzene where the 2 refers to the benzene NOT the methyl
			 */
			ArrayList<Element> locantElements =new ArrayList<Element>();
			Elements childrenOfElementBeforeSubstituent  =elementBeforeSubstituent.getChildElements();
			int nonStereoChemistryLocants =0;
			for (int i = 0; i < childrenOfElementBeforeSubstituent.size(); i++) {
				String currentElementName = childrenOfElementBeforeSubstituent.get(i).getLocalName();
				if (currentElementName.equals("stereoChemistry")){
					locantElements.add(childrenOfElementBeforeSubstituent.get(i));
				}
				else if (currentElementName.equals("locant")){
					locantElements.add(childrenOfElementBeforeSubstituent.get(i));
					nonStereoChemistryLocants++;
				}
				else{
					break;
				}
			}

			//either all locants will be moved, or none
			Boolean moveLocants=null;
			int flag=0;
			for (Element locant : locantElements) {
				if (!locant.getLocalName().equals("locant")){
					continue;//ignore stereochemistry elements, at least for the moment
				}
				String locantText = locant.getAttributeValue("value");

				if (theSubstituentGroupFragment.hasLocant("2")){//if only has locant 1 then assume substitution onto it not intended
					if(!theSubstituentGroupFragment.hasLocant(locantText)){
						flag=1;
					}
				}
				else{
					flag=1;
				}
			}
			if (flag==0){
				moveLocants =false;//if the locant applies to the theSubstituentGroup then don't move
			}

			flag =0;
			if (moveLocants ==null){
				for (Element locant : locantElements) {
					if (!locant.getLocalName().equals("locant")){
						continue;//ignore stereochemistry elements, atleast for the moment
					}
					String locantText = locant.getAttributeValue("value");
					if (!checkLocantPresentOnPotentialRoot(state, theSubstituent, locantText)){
						flag =1;
					}
				}
				if (flag==0){
					moveLocants = true;//if the locant applies to a group which is not theSubstituentGroup then move
				}
				else{
					moveLocants =false;
				}
			}
			if (moveLocants && nonStereoChemistryLocants>1){
				Element shouldBeAMultiplierNode = (Element)XOMTools.getNextSibling(locantElements.get(locantElements.size()-1));
				if (shouldBeAMultiplierNode instanceof Element){
					Element shouldBeAGroupOrSubOrBracket = (Element)XOMTools.getNextSibling(shouldBeAMultiplierNode);
					if (shouldBeAGroupOrSubOrBracket instanceof Element && shouldBeAMultiplierNode.getLocalName().equals("multiplier")
							&& shouldBeAGroupOrSubOrBracket.getLocalName().equals("group") ||
							shouldBeAGroupOrSubOrBracket.getLocalName().equals("substituent") ||
							shouldBeAGroupOrSubOrBracket.getLocalName().equals("bracket")){
						if (matchInlineSuffixesThatAreAlsoGroups.matcher(theSubstituentGroup.getValue()).matches()){//e.g. 4, 4'-dimethoxycarbonyl-2, 2'-bioxazole
							locantElements.add(shouldBeAMultiplierNode);
						}
						else{//don't bracket complex multiplied substituents
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
			if (moveLocants){
				for (int i = 0; i < locantElements.size(); i++) {
					locantElements.get(i).detach();
					bracket.appendChild(locantElements.get(i));
				}
			}

			/*
			 * A special case when a multiplier should be moved
			 * e.g. tripropan-2-yloxyphosphane -->tri(propan-2-yloxy)phosphane
			 */
			if (locantElements.size()==0 && matchInlineSuffixesThatAreAlsoGroups.matcher(theSubstituentGroup.getValue()).matches()){
				Element possibleMultiplier =childrenOfElementBeforeSubstituent.get(0);
				if (possibleMultiplier.getLocalName().equals("multiplier")){
					if (childrenOfElementBeforeSubstituent.get(1).getLocalName().equals("group")){
						childrenOfElementBeforeSubstituent.get(0).detach();
						bracket.appendChild(childrenOfElementBeforeSubstituent.get(0));
					}
				}
			}

			Element parent = (Element)theSubstituent.getParent();
			int startIndex=parent.indexOf(elementBeforeSubstituent);
			int endIndex=parent.indexOf(theSubstituent);
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
	 * @throws PostProcessingException
	 * @throws StructureBuildingException
	 */
	private void matchLocantsToIndirectFeatures(BuildState state, Element subOrRoot) throws PostProcessingException, StructureBuildingException {
		/* Root fragments (or the root in a bracket) can have prefix-locants
		 * that work on suffixes - (2-furyl), 2-propanol, (2-propylmethyl), (2-propyloxy), 2'-Butyronaphthone.
		 */
		Elements locants = subOrRoot.getChildElements("locant");
		List<Element> locantList = new ArrayList<Element>();//needed because some locants can readily be determined to be needed for other purposes
		for (int i = 0; i < locants.size(); i++) {
			Element locant =locants.get(i);
			Element afterLocants =OpsinTools.getNextIgnoringCertainElements(locant, new String[]{"locant"});
			if (afterLocants!=null && afterLocants.getLocalName().equals("multiplier")){//locant should not be followed by a multiplier. c.f. 1,2,3-tributyl 2-acetyloxypropane-1,2,3-tricarboxylate
				continue;
			}
			locantList.add(locant);
		}
		Element group =subOrRoot.getFirstChildElement("group");
		if (locantList.size()==1 && group!=null && group.getAttribute("frontLocantsExpected")!=null){//some trivial retained names like 2-furyl expect locants to be in front of them. For these the indirect intepretation will always be used rather than checking whether 2-(furyl) even makes sense
			Element locant =locantList.get(0);
			String locantValue =locant.getAttributeValue("value");
			String[] allowedLocants=matchComma.split(group.getAttributeValue("frontLocantsExpected"));
			for (String allowedLocant : allowedLocants) {
				if (locantValue.equals(allowedLocant)){
					Element expectedSuffix =(Element) XOMTools.getNextSibling(group);
					if (expectedSuffix!=null && expectedSuffix.getLocalName().equals("suffix") && expectedSuffix.getAttribute("locant")==null){
						locantList.clear();
						expectedSuffix.addAttribute(new Attribute("locant", locantValue));
						locant.detach();
					}
					break;
				}
			}
		}
		if (locantList.size()>0 && group!=null){
			boolean allowIndirectLocants =true;
			if(state.wordRule.equals("diester")){//special case e.g. 1-benzyl 4-butyl terephthalate (locants do not apply to yls)
				Element parentEl=(Element) subOrRoot.getParent();
				if (parentEl.getLocalName().equals("word") && parentEl.getAttributeValue("type").equals("substituent") && parentEl.getChildCount()==1 &&
						locants.size()==1 && (locants.get(0).getAttribute("type")==null || !locants.get(0).getAttributeValue("type").equals("orthoMetaPara"))){
					allowIndirectLocants =false;
				}
			}

			ArrayList<Element> locantsToAssignToIndirectFeatures = new ArrayList<Element>();
			ArrayList<Element> locantAble =null;
			if (allowIndirectLocants){
				Fragment thisFrag =state.xmlFragmentMap.get(group);
				for (int i = locants.size()-1; i >=0 ; i--) {
					Element locant =locants.get(i);
					String locantValue =locant.getAttributeValue("value");
					if (!checkLocantPresentOnPotentialRoot(state, subOrRoot, locantValue)){
						if (thisFrag.hasLocant(locantValue)){//locant not available elsewhere and is available on the group associated with this element
							if (locantAble ==null){
								//Find elements that can have locants but don't currently
								Elements childrenOfSubOrBracketOrRoot=subOrRoot.getChildElements();
								locantAble=new ArrayList<Element>();
								for (int j = 0; j < childrenOfSubOrBracketOrRoot.size(); j++) {
									Element el =childrenOfSubOrBracketOrRoot.get(j);
									String name =el.getLocalName();
									if (name.equals("suffix") || name.equals("unsaturator") || name.equals("heteroatom") || name.equals("hydro")){
										if (el.getAttribute("locant") ==null && (el.getAttribute("multiplied")==null)){// shouldn't already have a locant or be multiplied (should of already had locants assignd to it if that were the case)
											if (subOrRoot.indexOf(el)>subOrRoot.indexOf(locant)){
												locantAble.add(el);
											}
										}
									}
								}
							}
							if (locantsToAssignToIndirectFeatures.size() < locantAble.size()){
								locantsToAssignToIndirectFeatures.add(0, locant);//last in first out
							}
						}
						else{//usually indicates the name will fail unless the suffix has the locant or heteroatom replacement will create the locant
							List<Fragment> suffixes =state.xmlSuffixMap.get(group);
							//I do not want to assign element locants as in locants on the suffix as I currently know of no examples where this actually occurs
							if (matchElementSymbol.matcher(locantValue).matches()==true || locant.getAttribute("compoundLocant")!=null){
								continue;
							}
							for (Fragment suffix : suffixes) {
								if (suffix.hasLocant(locantValue)){
									suffix.setDefaultInID(suffix.getAtomByLocantOrThrow(locantValue).getID());
									locant.detach();
								}
							}
						}
					}
				}
			}
			for (Element locant : locantsToAssignToIndirectFeatures) {
				Element locantAbleElement =locantAble.get(0);
				locantAble.remove(0);

				String locantValue =locant.getAttributeValue("value");
				//If a compound Locant e.g. 1(6) is detected add a compound locant attribute
				if (locant.getAttribute("compoundLocant")!=null){
					locantAbleElement.addAttribute(new Attribute("compoundLocant", locant.getAttributeValue("compoundLocant")));
				}
				locantAbleElement.addAttribute(new Attribute("locant", locantValue));
				locant.detach();
			}
		}

		//put di-carbon modifying suffixes e.g. oic acids, aldehydes on opposite ends of chain
		Elements suffixEls = subOrRoot.getChildElements("suffix");
		for (int i = 0; i < suffixEls.size()-1; i++) {
			Element diCarbonModifyingSuffix1 = suffixEls.get(i);
			if (matchSuffixesThatGoAtEndOfChainsByDefault.matcher(diCarbonModifyingSuffix1.getAttributeValue("value")).matches()){
				if ((diCarbonModifyingSuffix1).getAttribute("locant")==null && XOMTools.getNextSibling(diCarbonModifyingSuffix1) != null){
					Element diCarbonModifyingSuffix2 =(Element)XOMTools.getNextSibling(diCarbonModifyingSuffix1);
					if (diCarbonModifyingSuffix2.getLocalName().equals("suffix") &&
							matchSuffixesThatGoAtEndOfChainsByDefault.matcher(diCarbonModifyingSuffix2.getAttributeValue("value")).matches() &&
							diCarbonModifyingSuffix2.getAttribute("locant")==null){
						Element hopefullyAChain = (Element) XOMTools.getPreviousSibling((Element)diCarbonModifyingSuffix1, "group");
						if (hopefullyAChain != null && hopefullyAChain.getAttributeValue("type").equals("chain")){
							((Element)diCarbonModifyingSuffix1).addAttribute(new Attribute("locant", "1"));
							((Element)diCarbonModifyingSuffix2).addAttribute(new Attribute("locant", Integer.toString(state.xmlFragmentMap.get(hopefullyAChain).getChainLength())));
							break;
						}
					}
				}
			}
		}
	}


	/**
	 * Resolves the effects of remaining suffixes and attaches resolved suffixes.
	 * @param state
	 * @param subOrRoot: The sub/root to look in
	 * @throws StructureBuildingException
	 * @throws PostProcessingException
	 */
	private void resolveRemainingSuffixes(BuildState state, Element subOrRoot) throws StructureBuildingException, PostProcessingException{
		Element group =subOrRoot.getFirstChildElement("group");
		Fragment thisFrag = state.xmlFragmentMap.get(group);
		Elements suffixes =subOrRoot.getChildElements("suffix");
		resolveSuffixes(state, thisFrag, suffixes, group);
	}

	/**Process the effects of suffixes upon a fragment.
	 * @param state
	 * @param frag The fragment to which to add the suffixes
	 * @param suffixes The suffix elements for a fragment.
	 * @param group The element for the group that the fragment will attach to
	 * @throws StructureBuildingException If the suffixes can't be resolved properly.
	 * @throws PostProcessingException
	 */
	private void resolveSuffixes(BuildState state, Fragment frag, Elements suffixes,  Element group) throws StructureBuildingException, PostProcessingException {
		int firstAtomID = frag.getIdOfFirstAtom();//typically equivalent to locant 1
		List<Atom> atomList =frag.getAtomList();//this instance of atomList will not change even once suffixes are merged into the fragment
		int defaultAtom =0;//indice in atomList
		String groupType = frag.getType();
		String subgroupType = frag.getSubType();
		String suffixTypeToUse =null;
		boolean suffixTypeDetermined =false;
		if (suffixApplicability.containsKey(groupType)){
			suffixTypeToUse =groupType;
			suffixTypeDetermined=true;
		}
		else if (groupType.equals("simpleGroup")){
			suffixTypeToUse="substituent";
			suffixTypeDetermined=true;
		}
		List<Fragment> suffixList =state.xmlSuffixMap.get(group);
		for(int i=0;i<suffixes.size();i++) {
			Element suffix = suffixes.get(i);
			String suffixValue = suffix.getAttributeValue("value");

			if (!suffixTypeDetermined){
				boolean cyclic;
				if (suffix.getAttribute("locant")!=null){
					Atom a =frag.getAtomByLocant(suffix.getAttributeValue("locant"));
					if (a!=null){
						cyclic=a.getAtomIsInACycle();
					}
					else{//can happen in the cases of things like fused rings where the final numbering is not available (in which case all the atoms will be cyclic anyway)
						cyclic = frag.getFirstAtom().getAtomIsInACycle();
					}
				}
				else{
					cyclic = frag.getFirstAtom().getAtomIsInACycle();
				}
				if (cyclic){
					suffixTypeToUse="cyclic";
				}
				else{
					suffixTypeToUse="acyclic";
				}
			}

			String locant = StructureBuildingMethods.getLocant(suffix);
			int idOnParentFragToUse=0;
			if (!locant.equals("0")){
				idOnParentFragToUse =frag.getIDFromLocantOrThrow(locant);
			}
			if (idOnParentFragToUse==0 && suffix.getAttribute("locantID")!=null){
				idOnParentFragToUse = Integer.parseInt(suffix.getAttributeValue("locantID"));
			}

			String compoundLocant = null;
			if (suffix.getAttribute("compoundLocant")!=null){
				compoundLocant=suffix.getAttributeValue("compoundLocant");
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
						suffixFrag = suffixList.get(0);//take the first suffix out of the list, it should of been added in the same order that it is now being read.
						suffixList.remove(0);
					}
					InID connectingInfo = suffixFrag.getInID(0);
					int bondOrder =connectingInfo.valency;
					suffixFrag.removeInID(0);

					if(idOnParentFragToUse==0) {
						if(suffixRuleTag.getAttribute("ketoneLocant") != null && suffixRuleTag.getAttributeValue("ketoneLocant").equals("yes")) {
							if(defaultAtom == 0) defaultAtom = state.fragManager.findKetoneAtomIndice(frag, defaultAtom);
							idOnParentFragToUse = atomList.get(defaultAtom).getID();
							defaultAtom++;
						}
						else{
							idOnParentFragToUse =atomList.get(defaultAtom).getID();
						}
						idOnParentFragToUse =frag.getAtomByIdOrNextSuitableAtomOrThrow(idOnParentFragToUse, bondOrder).getID();
					}
					//create a new bond and associate it with the suffixfrag and both atoms. Remember the suffixFrag has not been imported into the frag yet
					Atom suffixAtom = suffixFrag.getAtomByIDOrThrow(connectingInfo.id);
					Atom parentfragAtom = frag.getAtomByIDOrThrow(idOnParentFragToUse);
					if (parentfragAtom.getCharge()==1 && parentfragAtom.getElement().equals("N") && suffixAtom.getElement().equals("O") && suffixAtom.getCharge()==0 && bondOrder ==2){//special case to cope with azinic acid and the like
						suffixAtom.setCharge(-1);
						bondOrder =1;
					}
					Bond newBond =new Bond(suffixAtom,parentfragAtom, bondOrder);
					suffixFrag.getAtomByIDOrThrow(connectingInfo.id).addBond(newBond);
					frag.getAtomByIDOrThrow(idOnParentFragToUse).addBond(newBond);
					suffixFrag.addBond(newBond,false);
					if(suffixRuleTag.getAttribute("setsDefaultInID") != null) {
						frag.setDefaultInID(suffixFrag.getDefaultInID());
					}
					if (suffixValue.equals("one") && groupType.equals("ring")){//special case: one acts in a similar way to the hydro tag c.f. tetrahydrobenzen-1,4-dione
						frag.getAtomByIDOrThrow(idOnParentFragToUse).setNote("OneSuffixAttached", "1");
					}
				} else if(suffixRuleTagName.equals("doublebond")) {
					if(idOnParentFragToUse==0){
						idOnParentFragToUse = atomList.get(defaultAtom).getID();
						defaultAtom +=2;
					}
					if (compoundLocant!=null){
						state.fragManager.unsaturate(idOnParentFragToUse, compoundLocant, 2, frag);
					}
					else{
						state.fragManager.unsaturate(idOnParentFragToUse, 2, frag);
					}
				} else if(suffixRuleTagName.equals("triplebond")) {
					if(idOnParentFragToUse==0){
						idOnParentFragToUse = atomList.get(defaultAtom).getID();
						defaultAtom +=2;
					}
					if (compoundLocant!=null){
						state.fragManager.unsaturate(idOnParentFragToUse, compoundLocant, 3, frag);
					}
					else{
						state.fragManager.unsaturate(idOnParentFragToUse, 3, frag);
					}
				} else if(suffixRuleTagName.equals("changecharge")) {
					int chargeChange =Integer.parseInt(suffixRuleTag.getAttributeValue("charge"));
					if(idOnParentFragToUse==0){
						//Typically if a locant has not been specified and it's an ium then it was intended to refer to a nitrogen even if the nitrogen is not at locant 1 e.g. isoquinolinium
						//It's extremely rare to want a carbocation
						if (chargeChange == 1){
							Atom possibleAtom =null;
							for (Atom a : atomList) {
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
						else{
							idOnParentFragToUse =atomList.get(defaultAtom).getID();
							defaultAtom++;
						}
					}
					state.fragManager.changeCharge(idOnParentFragToUse, chargeChange, frag);
				}else if(suffixRuleTagName.equals("setOutID")) {
					if(suffixRuleTag.getAttribute("outValency") != null) {
						if(idOnParentFragToUse!=0){
							frag.addOutID(idOnParentFragToUse, Integer.parseInt(suffixRuleTag.getAttributeValue("outValency")), true);
						}
						else{
							frag.addOutID(firstAtomID, Integer.parseInt(suffixRuleTag.getAttributeValue("outValency")), false);
						}
					}
					else{
						if(idOnParentFragToUse!=0){
							frag.addOutID(idOnParentFragToUse, 1, true);
						}
						else{
							frag.addOutID(firstAtomID, 1, false);
						}
					}
				}else if(suffixRuleTagName.equals("setOptionalOutID")) {
					if(suffixRuleTag.getAttribute("outValency") != null) {
						if(idOnParentFragToUse!=0){
							frag.addOutID(idOnParentFragToUse, Integer.parseInt(suffixRuleTag.getAttributeValue("outValency")), true);
						}
						else{
							frag.addOutID(firstAtomID, Integer.parseInt(suffixRuleTag.getAttributeValue("outValency")), false);
						}
					}
					else{
						if(idOnParentFragToUse!=0){
							frag.addOutID(idOnParentFragToUse, 1, true);
						}
						else{
							frag.addOutID(firstAtomID, 1, false);
						}
					}
					if (group.getAttribute("substitutionRemovesOutIDs")==null){
						group.addAttribute(new Attribute("substitutionRemovesOutIDs", "yes"));//c.f. sulfonyl
					}
				}
				else{
					throw new StructureBuildingException("Unknown suffix rule:" + suffixRuleTagName);
				}
			}

			if (suffixFrag!=null){//merge suffix frag and parent fragment
				frag.importFrag(suffixFrag);
				state.fragManager.removeFragment(suffixFrag);
			}
		}
	}

	/**
	 * Removes brackets that only contain one element.
	 * Removed brackets are reflected in brackets and substituentsAndRootAndBrackets
	 * @param brackets
	 * @param substituentsAndRootAndBrackets
	 */
	private void removePointlessBrackets(List<Element> brackets, List<Element> substituentsAndRootAndBrackets) {
		for (int i = brackets.size()-1; i >=0; i--) {
			Element bracket =brackets.get(i);
			if (bracket.getChildCount()==1){//pointless bracket
				Element elToBeMoved = (Element) bracket.getChild(0);
				elToBeMoved.detach();
				XOMTools.insertAfter(bracket, elToBeMoved);
				bracket.detach();
				brackets.remove(i);
				substituentsAndRootAndBrackets.remove(bracket);
			}
		}
	}


	/**
	 * Uses the number of outIDs that are present to assign the number of outIDs on substituents that can have a variable number of outIDs
	 * Hence at this point it can be determined if a multi radical susbtituent is present in the name
	 * This would be expected in multiplicative nomenclature and is noted in the state so that the StructureBuilder knows to resolve the
	 * section of the name from that point onwards in a left to right manner rather than right to left
	 * @param state
	 * @param subOrRoot: The sub/root to look in
	 * @throws PostProcessingException
	 * @throws StructureBuildingException
	 */
	private void handleMultiRadicals(BuildState state, Element subOrRoot) throws PostProcessingException, StructureBuildingException{
		Element group =subOrRoot.getFirstChildElement("group");
		Fragment thisFrag = state.xmlFragmentMap.get(group);
		if (group.getAttribute("outIDs")!=null){//adds outIDs at the specified atoms
			String[] radicalPositions = matchComma.split(group.getAttributeValue("outIDs"));
			int firstIdInFrag =thisFrag.getIdOfFirstAtom();
			for (int i = 0; i < radicalPositions.length; i++) {
				String radicalID =radicalPositions[i];
				thisFrag.addOutID(firstIdInFrag + Integer.parseInt(radicalID) -1, 1, true);
			}
		}
		int outIDsSize = thisFrag.getOutIDs().size();
		if (outIDsSize >=2){
			group.addAttribute(new Attribute ("isAMultiRadical", "yes"));
			Element previousGroup =(Element) OpsinTools.getPreviousGroup(group);
			if (group.getValue().equals("amine")){//amine is a special case as it shouldn't technically be allowed but is allowed due to it's common usage in EDTA
				if (previousGroup==null || state.xmlFragmentMap.get(previousGroup).getOutIDs().size() < 2){//must be preceded by a multi radical
					throw new PostProcessingException("Invalid use of amine as a substituent!");
				}
			}
		}

		//Element possibleMultiplier= (Element)OpsinTools.getNext(subOrRoot.getChild(subOrRoot.getChildCount()-1));

//		System.out.println("OUDIDS" +thisFrag.getOutIDs().size());
//		System.out.println(subOrRoot.toXML());
//		if (possibleMultiplier !=null){//i.e. this is not the root -there is something after the group
//			if (outIDsSize>=2 && !possibleMultiplier.getLocalName().equals("multiplier")){//look for cases like methylenecyclohexane. need to convert 2 outIDs into 1 outID valency 2
//				List<OutID> outIDs =thisFrag.getOutIDs();
//				int id = outIDs.get(0).id;
//				boolean allOnSameAtom =true;
//				for (OutID outID : outIDs) {//all outIDs must be on the same atom
//					if (id!=outID.id){
//						allOnSameAtom=false;
//					}
//				}
//				if (allOnSameAtom){
//					if (state.mode.equals(OpsinMode.poly)){
//						if (outIDsSize>2){//e.g. nitrilo
//							for (int i = outIDs.size() -1; i >=1 ; i--) {
//								thisFrag.removeOutID(outIDs.get(i));
//							}
//							//set the outIds so the one with valency greater than 2 is the first and the other is second
//							OutID singleBondoutId = thisFrag.getOutID(0);
//							thisFrag.removeOutID(0);
//							thisFrag.addOutID(id, outIDsSize -1, true);
//							thisFrag.addOutID(singleBondoutId);
//						}
//					}else{
//						for (int i = outIDs.size() -1; i >=0 ; i--) {
//							thisFrag.removeOutID(outIDs.get(i));
//						}
//						thisFrag.addOutID(id, outIDsSize, true);
//					}
//				}
//			}
//			else if (possibleMultiplier.getLocalName().equals("multiplier") && outIDsSize==1){//special case where something like benzylidene is being used as if it meant benzdiyl for multiplicative nomenclature
//				//this is allowed in the IUPAC 79 recommendations but not recommended in the current recommendations
//				OutID outID =thisFrag.getOutID(0);
//				if (outID.valency == Integer.parseInt(possibleMultiplier.getAttributeValue("value"))){
//					Element parentWord =OpsinTools.getParentWord(group);
//					Element root =parentWord.getFirstChildElement("root");
//					if (!parentWord.getAttributeValue("type").equals("full")|| ((Element)root.getChild(0)).getLocalName().equals("multiplier")){//checks that the name appears to be multiplicative
//						int value =outID.valency;
//						if (outID.setExplicitly){
//							thisFrag.getAtomByIDOrThrow(outID.id).addOutValency(-(outID.valency-1));//correct outValency
//						}
//						outID.valency=1;
//						for (int i = 1; i < value; i++) {
//							thisFrag.addOutID(outID.id, 1, outID.setExplicitly);
//						}
//					}
//				}
//			}
//		}
//		System.out.println("OUDIDS" +thisFrag.getOutIDs().size());
//		System.out.println(subOrRoot.toXML());
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
			if (OpsinTools.getPrevious(multiplier)==null){
				throw new StructureBuildingException("OPSIN bug: Unacceptable input to function");
			}
			List<Element> locants = OpsinTools.findChildElementsWithTagNames(rightMostElement, new String[] {"multiplicativeLocant"});
			int multiVal = Integer.parseInt(multiplier.getAttributeValue("value"));
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
					locantString += locant.getAttributeValue("value");
					locant.detach();
				}
				rightMostElement.addAttribute(new Attribute("inLocants",locantString));
			}
			else{
				throw new PostProcessingException("Mismatch between number of locants and number of roots");
			}
		}
	}


	/**
	 * Assigns locants and multipliers to substituents/brackets
	 * If both locants and multipliers are present a final check is done that the number of them agree.
	 * WordLevel multipliers are processed e.g. diethyl ethanoate
	 * Adding a locant to a root or any other group that cannot engage in substitive nomenclature will result in an exception being thrown
	 * An exception is made for cases where the locant could be referring to a position on another word
	 * @param subOrBracket
	 * @throws PostProcessingException
	 */
	private void assignLocantsAndMultipliers(Element subOrBracket) throws PostProcessingException {
		Elements locants =subOrBracket.getChildElements("locant");
		int multiplier =1;
		Element possibleMultiplier = subOrBracket.getFirstChildElement("multiplier");
		if (possibleMultiplier!=null){
			Element parentElem =(Element)subOrBracket.getParent();
			if(parentElem.getLocalName().equals("word") &&
					XOMTools.getNextSibling(subOrBracket) == null &&
					XOMTools.getPreviousSibling(subOrBracket) == null) {
				return;//word level multiplier
			}
			multiplier = Integer.parseInt(possibleMultiplier.getAttributeValue("value"));
			subOrBracket.addAttribute(new Attribute("multiplier", possibleMultiplier.getAttributeValue("value")));
			possibleMultiplier.detach();
		}
		if(locants.size() > 0) {
			if (subOrBracket.getLocalName().equals("root")){
				throw new PostProcessingException("Unable to assign all locants");
			}
			if (multiplier !=locants.size()){
				throw new PostProcessingException("Multiplier and locant count failed to agree. This is an OPSIN bug as such a check should of been performed earlier!");
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
				locantString += locant.getAttributeValue("value");
				locant.detach();
			}
			subOrBracket.addAttribute(new Attribute("locant", locantString));
		}
	}
	
	/**
	 * If a word level multiplier is present e.g. diethyl butandioate then this is processed to ethyl ethyl butandioate
	 * @param state
	 * @param word
	 * @throws StructureBuildingException 
	 * @throws PostProcessingException 
	 */
	private void processWordLevelMultiplierIfApplicable(BuildState state, Element word) throws StructureBuildingException, PostProcessingException {
		if (word.getChildCount()==1){
			Element subOrBracket = (Element) word.getChild(0);
			Element multiplier = subOrBracket.getFirstChildElement("multiplier");
			if (multiplier !=null){
				int multiVal =Integer.parseInt(multiplier.getAttributeValue("value"));
				multiplier.detach();
				if (multiVal ==1){return;}
				Elements locants =subOrBracket.getChildElements("locant");
				boolean assignLocants =false;
				if (locants.size()==multiVal){
					assignLocants=true;
					for (int i = 0; i < locants.size(); i++) {
						if(locants.get(i).getAttribute("compoundLocant")!=null){
							throw new PostProcessingException("A compound locant cannot be used to locant a sub/bracket!");
						}
					}
					subOrBracket.addAttribute(new Attribute("locant", locants.get(0).getAttributeValue("value")));
				}
				else if (locants.size()!=0){
					throw new PostProcessingException("Unable to assign all locants");
				}
				for(int i=multiVal -1; i>=1; i--) {
					Element clone = state.fragManager.cloneElement(word, state);
					if (assignLocants){
						subOrBracket.addAttribute(new Attribute("locant", locants.get(i).getAttributeValue("value")));
					}
					XOMTools.insertAfter(word, clone);
				}
			}
		}
	}
}
