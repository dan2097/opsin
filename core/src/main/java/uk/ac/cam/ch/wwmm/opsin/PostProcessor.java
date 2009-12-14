package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


import nu.xom.Attribute;
import nu.xom.Element;
import nu.xom.Elements;
import nu.xom.Node;

/**Does destructive procedural parsing on parser results.
 *
 * @author ptc24/dl387
 *
 */
class PostProcessor {

	/**
	 * Sort bridges such as the highest priority secondary bridges come first
	 * e.g. 1^(6,3).1^(15,13)
	 * rearranged to 1^(15,13).1^(6,3)
	 * @author dl387
	 *
	 */
	private class VonBaeyerSecondaryBridgeSort implements Comparator<HashMap<String, Integer>> {

	    public int compare(HashMap<String, Integer> bridge1, HashMap<String, Integer> bridge2){
	    	//first we compare the larger coordinate, due to an earlier potential swapping of coordinates this is always in  "AtomId_Larger"
	    	int largerCoordinate1 = bridge1.get("AtomId_Larger");
	    	int largerCoordinate2 = bridge2.get("AtomId_Larger");
			if (largerCoordinate1 >largerCoordinate2) {
				return -1;
			}
			else if (largerCoordinate2 >largerCoordinate1) {
				return 1;
			}
			//tie
	    	int smallerCoordinate1 = bridge1.get("AtomId_Smaller");
	    	int smallerCoordinate2 = bridge2.get("AtomId_Smaller");
			if (smallerCoordinate1 >smallerCoordinate2) {
				return -1;
			}
			else if (smallerCoordinate2 >smallerCoordinate1) {
				return 1;
			}
			//tie
	    	int bridgelength1 = bridge1.get("Bridge Length");
	    	int bridgelength2 = bridge2.get("Bridge Length");
			if (bridgelength1 >bridgelength2) {
				return -1;
			}
			else if (bridgelength2 >bridgelength1) {
				return 1;
			}
			else{
				return 0;
			}
	    }
	}

	//match a fusion bracket with only numerical locants. If this is followed by a HW group it probably wasn't a fusion bracket
	private Pattern matchNumberLocantsOnlyFusionBracket = Pattern.compile("\\[\\d+[a-z]?(,\\d+[a-z]?)*\\]");
	private Pattern matchVonBaeyer = Pattern.compile("(\\d+\\^?[\\({]?\\d*,?\\d*[\\)}]?\\^?\\^?)");
	private Pattern matchAnnulene = Pattern.compile("\\[([1-9]\\d*)\\]annulen");
	private Pattern matchStereochemistry = Pattern.compile("((?:\\d+[a-z]?|[A-Z][a-z]?)'*)([RSEZrsez])");
	private Pattern matchRS = Pattern.compile("[RSrs]");
	private Pattern matchEZ = Pattern.compile("[EZez]");
	private Pattern matchLambdaConvention = Pattern.compile("(\\S+)?lambda\\D*(\\d+)\\D*");
	private Pattern matchComma =Pattern.compile(",");
	private Pattern matchDot =Pattern.compile("\\.");
	private Pattern matchNonDigit =Pattern.compile("\\D+");

	private TokenManager tokenManager;

	PostProcessor(TokenManager tokenManager) {
		this.tokenManager =tokenManager;
	}

	/** The master method, postprocesses a parse result.
	 *
	 * @param elem The element to postprocess.
	 * @param state
	 * @return
	 * @throws Exception
	 */
	void postProcess(Element elem, BuildState state) throws Exception {
		state.wordRule=elem.getAttributeValue("wordRule");
		/* Throws exceptions for occurrences that are ambiguous and this parse has picked the incorrect interpretation */
		resolveAmbiguities(elem);

		List<Element> substituentsAndRoot = XOMTools.getDescendantElementsWithTagNames(elem, new String[]{"substituent", "root"});

		for (Element subOrRoot: substituentsAndRoot) {
			processHeterogenousHydrides(subOrRoot);
			processIndicatedHydrogens(subOrRoot);
			processStereochemistry(subOrRoot);
			processInfixes(subOrRoot);
			processLambdaConvention(subOrRoot);
		}
		List<Element> groups =  XOMTools.getDescendantElementsWithTagName(elem, "group");

		processHydroCarbonRings(elem);
		for (Element group : groups) {
			processRings(group);//processes cyclo, von baeyer and spiro tokens
			handleIrregularities(group, state);//handles benzyl, diethylene glycol, phenanthrone and other awkward bits of nomenclature
		}

		/* Converts open/close bracket elements to bracket elements and
		 *  places the elements inbetween within the newly created bracket */
		while(findAndStructureBrackets(substituentsAndRoot));

		for (Element group : groups) {
			processHydroSubstituents(group);//this REMOVES hydro substituents and adds hydro elements in front of an appropriate ring
		}

		addOmittedSpaces(state, elem);//e.g. change ethylmethyl ether to ethyl methyl ether
	}

	/**
	 * Resolves common ambiguities e.g. tetradeca being 4x10carbon chain rather than 14carbon chain
	 * @param elem
	 * @throws PostProcessingException
	 */
	private void resolveAmbiguities(Element elem) throws PostProcessingException {
		List<Element> multipliers = XOMTools.getDescendantElementsWithTagName(elem, "multiplier");
		for (Element apparentMultiplier : multipliers) {
			Element nextEl = (Element)XOMTools.getNextSibling(apparentMultiplier);
			if(nextEl !=null && nextEl.getLocalName().equals("group")){//detects ambiguous use of things like tetradeca
				String multiplierAndGroup =apparentMultiplier.getValue() + nextEl.getValue();
				HashMap<String, HashMap<Character, Token>> tokenDict =tokenManager.tokenDict;
				HashMap<Character,Token> tokenMap = tokenDict.get(multiplierAndGroup);
				if (tokenMap !=null){
					Element isThisALocant =(Element)XOMTools.getPreviousSibling(apparentMultiplier);
					if (isThisALocant == null ||
							!isThisALocant.getLocalName().equals("locant") ||
							isThisALocant.getValue().split(",").length != Integer.parseInt(apparentMultiplier.getAttributeValue("value"))){
						throw new PostProcessingException(multiplierAndGroup +" should not have been lexed as two tokens!");
					}
				}
			}

			if (nextEl !=null && nextEl.getLocalName().equals("hydrocarbonFusedRingSystem")&& nextEl.getValue().equals("phen")){
				Element possibleSuffix = (Element) XOMTools.getNextSibling(nextEl);
				if (possibleSuffix!=null){//null if not used as substituent
					String multiplierAndGroup =apparentMultiplier.getValue() + nextEl.getValue();
					if (possibleSuffix.getValue().equals("yl")){
						throw new PostProcessingException(multiplierAndGroup +" should not have been lexed as one token!");
					}
					Element isThisALocant =(Element)XOMTools.getPreviousSibling(apparentMultiplier);
					if (isThisALocant != null && isThisALocant.getLocalName().equals("locant") && isThisALocant.getValue().split(",").length == Integer.parseInt(apparentMultiplier.getAttributeValue("value"))){
						throw new PostProcessingException(multiplierAndGroup +" should not have been lexed as one token!");
					}
				}
			}
		}

		List<Element> fusions = XOMTools.getDescendantElementsWithTagName(elem, "fusion");
		for (Element fusion : fusions) {
			if (matchNumberLocantsOnlyFusionBracket.matcher(fusion.getValue()).matches()){
				Element nextGroup =(Element) XOMTools.getNextSibling(fusion, "group");
				if (nextGroup !=null && nextGroup.getAttributeValue("subType").equals("hantzschWidman")){
					throw new PostProcessingException("This fusion bracket is in fact more likely to be a description of the locants of a HW ring");
				}
			}
		}
	}

	/**Form heterogeneous hydrides/substituents
	 * These are chains of one heteroatom or alternating heteroatoms and are expressed using SMILES
	 * They are typically treated in an analagous way to alkanes
	 * @param elem The root/substituents
	 */
	private void processHeterogenousHydrides(Element elem)  {
		Elements multipliers = elem.getChildElements("multiplier");
		for(int i=0;i<multipliers.size();i++) {
			Element m = multipliers.get(i);
			int mvalue = Integer.parseInt(m.getAttributeValue("value"));
			Element multipliedElem = (Element)XOMTools.getNextSibling(m);

			Element possiblyALocant = (Element)XOMTools.getPreviousSibling(m);
			if(possiblyALocant !=null && possiblyALocant.getLocalName().equals("locant")&& mvalue==matchComma.split(possiblyALocant.getValue()).length){
				continue;//something like 1,2-disulfanylpropane
			}

			if(multipliedElem.getLocalName().equals("group") &&
					multipliedElem.getAttribute("subType")!=null &&
					multipliedElem.getAttributeValue("subType").equals("heteroStem")) {
				//chain of heteroatoms
				String smiles=multipliedElem.getAttributeValue("value");
				multipliedElem.getAttribute("value").setValue(StringTools.multiplyString(smiles, mvalue));
				m.detach();
			}
		}
		Elements groups = elem.getChildElements("group");

		if (groups.size()==0){
			for(int i=0;i<multipliers.size();i++) {
				Element m = multipliers.get(i);
				Element multipliedElem = (Element)XOMTools.getNextSibling(m);
				if(multipliedElem.getLocalName().equals("heteroatom")){
					Element possiblyAnotherHeteroAtom = (Element)XOMTools.getNextSibling(multipliedElem);
					if (possiblyAnotherHeteroAtom !=null && possiblyAnotherHeteroAtom.getLocalName().equals("heteroatom")){
						Element possiblyAnUnsaturator = XOMTools.getNextSiblingIgnoringCertainElements(possiblyAnotherHeteroAtom, new String[]{"locant", "multiplier"});//typically ane but can be ene or yne e.g. triphosphaza-1,3-diene
						if (possiblyAnUnsaturator !=null && possiblyAnUnsaturator.getLocalName().equals("unsaturator")){
							//chain of alternating heteroatoms
							int mvalue = Integer.parseInt(m.getAttributeValue("value"));
							String smiles="";
							Element possiblyARingFormingEl = (Element)XOMTools.getPreviousSibling(m);
							boolean heteroatomChainWillFormARing =false;
							if (possiblyARingFormingEl!=null && (possiblyARingFormingEl.getLocalName().equals("cyclo") || possiblyARingFormingEl.getLocalName().equals("vonBaeyer") || possiblyARingFormingEl.getLocalName().equals("spiro"))){
								heteroatomChainWillFormARing=true;
								//will be cyclised later.
								//FIXME sort based on order in HW system (also add check to HW stuff that the heteroatoms in that are in the correct order so that incorrectly ordered systems may be rejected.
								for (int j = 0; j < mvalue; j++) {
									smiles+=possiblyAnotherHeteroAtom.getAttributeValue("value");
									smiles+=multipliedElem.getAttributeValue("value");
								}
							}
							else{
								for (int j = 0; j < mvalue -1; j++) {
									smiles+=multipliedElem.getAttributeValue("value");
									smiles+=possiblyAnotherHeteroAtom.getAttributeValue("value");
								}
								smiles+=multipliedElem.getAttributeValue("value");
							}
							multipliedElem.detach();

							Element addedGroup=new Element("group");
							addedGroup.addAttribute(new Attribute("value", smiles));
							addedGroup.addAttribute(new Attribute("valType", "SMILES"));
							addedGroup.addAttribute(new Attribute("type", "chain"));
							addedGroup.addAttribute(new Attribute("subType", "heteroStem"));
							if (!heteroatomChainWillFormARing){
								addedGroup.addAttribute(new Attribute("usableAsAJoiner", "yes"));
							}
							addedGroup.appendChild(smiles);
							XOMTools.insertAfter(possiblyAnotherHeteroAtom, addedGroup);

							possiblyAnotherHeteroAtom.detach();
							m.detach();
						}
					}
				}
			}
		}
	}

	/** Handle 1H- in 1H-pyrrole etc.
	 *
	 * @param elem The substituent/root to looks for indicated hydrogens in.
	 */
	private void processIndicatedHydrogens(Element elem) {
		Elements hydrogens = elem.getChildElements("hydrogen");
		for(int i=0;i<hydrogens.size();i++) {
			Element hydrogen = hydrogens.get(i);
			String txt = hydrogen.getChild(0).getValue();
			String[] hydrogenLocants =txt.split(",");
            for (String hydrogenLocant : hydrogenLocants) {
                if (hydrogenLocant.endsWith("H-")) {
                    Element newHydrogenElement = new Element("hydrogen");
                    newHydrogenElement.addAttribute(new Attribute("locant", hydrogenLocant.substring(0, hydrogenLocant.length() - 2)));
                    XOMTools.insertAfter(hydrogen, newHydrogenElement);
                } else if (hydrogenLocant.endsWith("H")) {
                    Element newHydrogenElement = new Element("hydrogen");
                    newHydrogenElement.addAttribute(new Attribute("locant", hydrogenLocant.substring(0, hydrogenLocant.length() - 1)));
                    XOMTools.insertAfter(hydrogen, newHydrogenElement);
                }
            }
			hydrogen.detach();
		}
	}

	/** Handles stereoChemistry (R/Z/E/Z/cis/trans)
	 *  Will assign a locant to a stereoChemistry element if one was specified/available
	 *
	 * @param elem The substituent/root to looks for stereoChemistry in.
	 * @throws PostProcessingException
	 */
	private void processStereochemistry(Element elem) throws PostProcessingException {
		Elements stereoChemistryElements = elem.getChildElements("stereoChemistry");
		for(int i=0;i<stereoChemistryElements.size();i++) {
			Element stereoChemistryElement = stereoChemistryElements.get(i);
			if (stereoChemistryElement.getAttributeValue("type").equals("stereochemistryBracket")){
				String txt = stereoChemistryElement.getValue();
				if (txt.startsWith("rel-") || txt.contains("*")){
					//currently unsupported
				}
				else{
					txt =txt.substring(1, txt.length()-1);//remove opening and closing bracket.
					String[] stereoChemistryDescriptors = matchComma.split(txt);
                    for (String stereoChemistryDescriptor : stereoChemistryDescriptors) {
                        if (stereoChemistryDescriptor.length() > 1) {
                            Matcher m = matchStereochemistry.matcher(stereoChemistryDescriptor);
                            if (m.matches()) {
                                Element stereoChemEl = new Element("stereoChemistry");
                                stereoChemEl.addAttribute(new Attribute("locant", m.group(1)));
                                stereoChemEl.addAttribute(new Attribute("value", m.group(2).toUpperCase()));
                                stereoChemEl.appendChild(stereoChemistryDescriptor);
                                XOMTools.insertAfter(stereoChemistryElement, stereoChemEl);
                                if (matchRS.matcher(m.group(2)).matches()) {
                                    stereoChemEl.addAttribute(new Attribute("type", "RorS"));
                                } else {
                                    stereoChemEl.addAttribute(new Attribute("type", "EorZ"));
                                }
                            } else {
                                throw new PostProcessingException("Malformed stereochemistry element: " + stereoChemistryElement.getValue());
                            }
                        } else {
                            Element stereoChemEl = new Element("stereoChemistry");
                            stereoChemEl.addAttribute(new Attribute("value", stereoChemistryDescriptor.toUpperCase()));
                            stereoChemEl.appendChild(stereoChemistryDescriptor);
                            XOMTools.insertAfter(stereoChemistryElement, stereoChemEl);
                            if (matchRS.matcher(stereoChemistryDescriptor).matches()) {
                                stereoChemEl.addAttribute(new Attribute("type", "RorS"));
                            } else if (matchEZ.matcher(stereoChemistryDescriptor).matches()) {
                                stereoChemEl.addAttribute(new Attribute("type", "EorZ"));
                            } else {
                                throw new PostProcessingException("Malformed stereochemistry element: " + stereoChemistryElement.getValue());
                            }
                        }
                    }
				}
				stereoChemistryElement.detach();
			}
			else if (stereoChemistryElement.getAttributeValue("type").equals("cisOrTrans")){//assign a locant if one is directly before the cis/trans
				Element possibleLocant = (Element) XOMTools.getPrevious(stereoChemistryElement);
				if (possibleLocant !=null && possibleLocant.getLocalName().equals("locant") && matchComma.split(possibleLocant.getValue()).length==1){
					stereoChemistryElement.addAttribute(new Attribute("locant", StringTools.removeDashIfPresent(possibleLocant.getValue())));
					possibleLocant.detach();
				}
			}
		}
	}

	/**
	 * Looks for infixes and assigns them to the next suffix
	 * If the infix/suffix block has been bracketed e.g (dithioate) then the infixCount is set appropriately
	 * If this is not the case then it is ambiguous as to whether the multiplier is referring to the infix or the infixed suffix
	 * This ambiguity is resolved in processInfixFunctionalReplacementNomenclature by looking at the structure of the suffix to be modified
	 * @param elem
	 * @throws PostProcessingException
	 */
	private void processInfixes(Element elem) throws PostProcessingException {
		List<Element> infixes = XOMTools.getDescendantElementsWithTagName(elem, "infix");
		for (Element infix : infixes) {
			Element suffix = (Element) XOMTools.getNextSibling(infix);
			if (suffix ==null || !suffix.getLocalName().equals("suffix")){
				throw new PostProcessingException("No suffix found next next to infix: "+ infix.getValue());
			}
			suffix.addAttribute(new Attribute("infix", infix.getAttributeValue("value")));
			suffix.addAttribute(new Attribute("infixCount", "1"));
			Element possibleMultiplier = (Element) XOMTools.getPreviousSibling(infix);
			Element possibleBracket;
			if (possibleMultiplier.getLocalName().equals("multiplier")){
				possibleBracket  = (Element) XOMTools.getPreviousSibling(possibleMultiplier);
			}
			else{
				possibleBracket=possibleMultiplier;
				possibleMultiplier=null;
				infix.detach();
			}
			if (possibleBracket.getLocalName().equals("structuralOpenBracket")){
				Element bracket = (Element) XOMTools.getNextSibling(suffix);
				if (!bracket.getLocalName().equals("structuralCloseBracket")){
					throw new PostProcessingException("Matching closing bracket not found around infix/suffix block");
				}
				if (possibleMultiplier!=null){
					suffix.getAttribute("infixCount").setValue(possibleMultiplier.getAttributeValue("value"));
					possibleMultiplier.detach();
					infix.detach();
				}
				possibleBracket.detach();
				bracket.detach();
			}
		}
	}

	/**
	 * Identifies lambdaConvention elements.
	 * The elementsValue is expected to be a comma seperated lambda values and 0 or more locants. Where a lambda value has the following form:
	 * optional locant, the word lambda and then a number which is the valency specified (with possibly some attempt to indicate this number is superscripted)
	 * If the element is followed by heteroatoms (possibly multiplied) they are multiplied and the locant/lambda assigned to them
	 * Otherwise a new lambdaConvention element is created with the valency specified by the lambda convention taking the attribute "lambda"
	 * In the case where heteroatoms belong to a fused ring system a new lambdaConvention element is also created. The original locants are retained in the benzo specific fused ring nomenclature:
	 * 2H-5lambda^5-phosphinino[3,2-b]pyran --> 2H 5lambda^5 phosphinino[3,2-b]pyran  BUT
	 * 1lambda^4,5-Benzodithiepin  --> 1lambda^4 1,5-Benzodithiepin
	 * @param subOrRoot
	 * @throws PostProcessingException
	 */
	private void processLambdaConvention(Element subOrRoot) throws PostProcessingException {
		List<Element> lambdaConventionEls = XOMTools.getChildElementsWithTagNames(subOrRoot, new String[]{"lambdaConvention"});
		boolean fusedRingPresent = false;
		if (lambdaConventionEls.size()>0){
			if (subOrRoot.getChildElements("group").size()>1){
				fusedRingPresent = true;
			}
		}
		for (Element lambdaConventionEl : lambdaConventionEls) {
			boolean benzoFusedRing =false;//Is the lambdaConvention el followed by benz/benzo of a fused ring system (these have front locants which correspond to the final fused rings numbering)
			String[] lambdaValues = matchComma.split(StringTools.removeDashIfPresent(lambdaConventionEl.getValue()));
			Element possibleHeteroatomOrMultiplier = (Element) XOMTools.getNextSibling(lambdaConventionEl);
			int heteroCount = 0;
			int multiplierValue = 1;
			while(possibleHeteroatomOrMultiplier != null){
				if(possibleHeteroatomOrMultiplier.getLocalName().equals("heteroatom")) {
					heteroCount+=multiplierValue;
					multiplierValue =1;
				} else if (possibleHeteroatomOrMultiplier.getLocalName().equals("multiplier")){
					multiplierValue = Integer.parseInt(possibleHeteroatomOrMultiplier.getAttributeValue("value"));
				}
				else{
					break;
				}
				possibleHeteroatomOrMultiplier = (Element)XOMTools.getNextSibling(possibleHeteroatomOrMultiplier);
			}
			boolean assignLambdasToHeteroAtoms =false;
			if (lambdaValues.length==heteroCount){//heteroatom and number of locants +lambdas must match
				if (fusedRingPresent && possibleHeteroatomOrMultiplier!=null && possibleHeteroatomOrMultiplier.getLocalName().equals("group") && possibleHeteroatomOrMultiplier.getAttributeValue("subType").equals("hantzschWidman")){
					//You must not set the locants of a HW system which forms a component of a fused ring system. The locant specified corresponds to the complete fused ring system.
				}
				else{
					assignLambdasToHeteroAtoms =true;
				}
			}
			else if(heteroCount==0 && fusedRingPresent &&
					possibleHeteroatomOrMultiplier!=null && possibleHeteroatomOrMultiplier.getLocalName().equals("group") &&
					possibleHeteroatomOrMultiplier.getValue().equals("benzo")||possibleHeteroatomOrMultiplier.getValue().equals("benz")){
				benzoFusedRing = true;
			}
			List<Element> heteroAtoms = new ArrayList<Element>();//contains the heteroatoms to apply the lambda values too. Can be empty if the values are applied to a group directly rather than to a heteroatom
			if (assignLambdasToHeteroAtoms){//populate heteroAtoms, multiplied heteroatoms are multiplied out
				Element multiplier = null;
				Element heteroatomOrMultiplier = (Element) XOMTools.getNextSibling(lambdaConventionEl);
				while(heteroatomOrMultiplier != null){
					if(heteroatomOrMultiplier.getLocalName().equals("heteroatom")) {
						heteroAtoms.add(heteroatomOrMultiplier);
						if (multiplier!=null){
							for (int i = 1; i < Integer.parseInt(multiplier.getAttributeValue("value")); i++) {
								Element newHeteroAtom = new Element(heteroatomOrMultiplier);
								XOMTools.insertAfter(heteroatomOrMultiplier, newHeteroAtom);
								heteroAtoms.add(newHeteroAtom);
							}
							multiplier.detach();
							multiplier=null;
						}
					} else if (heteroatomOrMultiplier.getLocalName().equals("multiplier")){
						if (multiplier !=null){
							break;
						}
						else{
							multiplier = heteroatomOrMultiplier;
						}
					}
					else{
						break;
					}
					heteroatomOrMultiplier = (Element)XOMTools.getNextSibling(heteroatomOrMultiplier);
				}
			}

			for (int i = 0; i < lambdaValues.length; i++) {//assign all the lambdas to heteroatoms or to newly created lambdaConvention elements
				String lambdaValue = lambdaValues[i];
				Matcher m = matchLambdaConvention.matcher(lambdaValue);
				if (m.matches()){//a lambda
					Attribute valencyChange = new Attribute("lambda", m.group(2));
					Attribute locantAtr = null;
					if (m.group(1)!=null){
						locantAtr = new Attribute("locant", m.group(1));
					}
					if (benzoFusedRing){
						if (m.group(1)==null){
							throw new PostProcessingException("Locant not found for lambda convention before a benzo fused ring system");
						}
						lambdaValues[i] = m.group(1);
					}
					if (assignLambdasToHeteroAtoms){
						Element heteroAtom = heteroAtoms.get(i);
						heteroAtom.addAttribute(valencyChange);
						if (locantAtr!=null){
							heteroAtom.addAttribute(locantAtr);
						}
					}
					else{
						Element newLambda = new Element("lambdaConvention");
						newLambda.addAttribute(valencyChange);
						if (locantAtr!=null){
							newLambda.addAttribute(locantAtr);
						}
						XOMTools.insertBefore(lambdaConventionEl, newLambda);
					}
				}
				else{//just a locant e.g 1,3lambda5
					if (!assignLambdasToHeteroAtoms){
						if (!benzoFusedRing){
							throw new PostProcessingException("Lambda convention not specified for locant: " + lambdaValue);
						}
					}
					else{
						Element heteroAtom = heteroAtoms.get(i);
						heteroAtom.addAttribute(new Attribute("locant", lambdaValue));
					}
				}
			}
			if (!benzoFusedRing){
				lambdaConventionEl.detach();
			}
			else{
				lambdaConventionEl.setLocalName("locant");
				XOMTools.setTextChild(lambdaConventionEl, StringTools.arrayToString(lambdaValues, ","));
			}
		}
	}

	/**Looks for annulen/polyacene/polyaphene/polyalene/polyphenylene/polynaphthylene/polyhelicene tags and replaces them with a group with appropriate SMILES.
	 * @param elem The element to look for tags in
	 * @throws PostProcessingException
	 */
	private void processHydroCarbonRings(Element elem) throws PostProcessingException {
		List<Element> annulens = XOMTools.getDescendantElementsWithTagName(elem, "annulen");
		for (Element annulen : annulens) {
			String annulenValue =annulen.getValue();
	        Matcher match = matchAnnulene.matcher(annulenValue);
	        match.matches();
	        if (match.groupCount() !=1){
	        	throw new PostProcessingException("Invalid annulen tag");
	        }

	        int annulenSize=Integer.valueOf(match.group(1));
	        if (annulenSize <3){
	        	throw new PostProcessingException("Invalid annulen tag");
	        }

			//build [annulenSize]annulene ring as SMILES
			String SMILES = "c1" +StringTools.multiplyString("c", annulenSize -1);
			SMILES += "1";

			Element group =new Element("group");
			group.addAttribute(new Attribute("value", SMILES));
			group.addAttribute(new Attribute("valType", "SMILES"));
			group.addAttribute(new Attribute("type", "ring"));
			group.addAttribute(new Attribute("subType", "annulene"));
			group.appendChild(annulenValue);
			annulen.getParent().replaceChild(annulen, group);
		}

		List<Element> hydrocarbonFRSystems = XOMTools.getDescendantElementsWithTagName(elem, "hydrocarbonFusedRingSystem");
		for (Element hydrocarbonFRSystem : hydrocarbonFRSystems) {
			Element multiplier = (Element)XOMTools.getPreviousSibling(hydrocarbonFRSystem);
			if(multiplier != null && multiplier.getLocalName().equals("multiplier")) {
				int multiplierValue =Integer.parseInt(multiplier.getAttributeValue("value"));
				String classOfHydrocarbonFRSystem =hydrocarbonFRSystem.getAttributeValue("value");
				Element newGroup =new Element("group");
				String SMILES="";
				if (classOfHydrocarbonFRSystem.equals("polyacene")){
					if (multiplierValue <=3){
						throw new PostProcessingException("Invalid polyacene");
					}
					SMILES= "c1ccc";
					for (int j = 2; j <= multiplierValue; j++) {
						SMILES+="c"+ringClosure(j);
						SMILES+="c";
					}
					SMILES+= "ccc";
					for (int j = multiplierValue; j >2; j--) {
						SMILES+="c"+ringClosure(j);
						SMILES+="c";
					}
					SMILES+="c12";
				}else if (classOfHydrocarbonFRSystem.equals("polyaphene")){
					if (multiplierValue <=3){
						throw new PostProcessingException("Invalid polyaphene");
					}
					SMILES= "c1ccc";

					int ringsAbovePlane;
					int ringsOnPlane;
					int ringOpeningCounter=2;
					if (multiplierValue %2==0){
						ringsAbovePlane =(multiplierValue-2)/2;
						ringsOnPlane =ringsAbovePlane +1;
					}
					else{
						ringsAbovePlane =(multiplierValue-1)/2;
						ringsOnPlane =ringsAbovePlane;
					}

					for (int j = 1; j <= ringsAbovePlane; j++) {
						SMILES+="c"+ringClosure(ringOpeningCounter++);
						SMILES+="c";
					}

					for (int j = 1; j <= ringsOnPlane; j++) {
						SMILES+="c";
						SMILES+="c"+ringClosure(ringOpeningCounter++);
					}
					SMILES+="ccc";
					ringOpeningCounter--;
					for (int j = 1; j <= ringsOnPlane; j++) {
						SMILES+="c";
						SMILES+="c"+ringClosure(ringOpeningCounter--);
					}
					for (int j = 1; j < ringsAbovePlane; j++) {
						SMILES+="c"+ringClosure(ringOpeningCounter--);
						SMILES+="c";
					}

					SMILES+="c12";
				} else if (classOfHydrocarbonFRSystem.equals("polyalene")){
					if (multiplierValue <5){
						throw new PostProcessingException("Invalid polyalene");
					}
					SMILES= "c1";
					for (int j = 3; j < multiplierValue; j++) {
						SMILES+="c";
					}
					SMILES+= "c2";
					for (int j = 3; j <= multiplierValue; j++) {
						SMILES+="c";
					}
					SMILES+="c12";
				} else if (classOfHydrocarbonFRSystem.equals("polyphenylene")){
					if (multiplierValue <2){
						throw new PostProcessingException("Invalid polyphenylene");
					}
					SMILES= "c1cccc2";
					for (int j = 1; j < multiplierValue; j++) {
						SMILES+="c3ccccc3";
					}
					SMILES+= "c12";
				} else if (classOfHydrocarbonFRSystem.equals("polynaphthylene")){
					if (multiplierValue <3){
						throw new PostProcessingException("Invalid polynaphthylene");
					}
					SMILES= "c1cccc2cc3";
					for (int j = 1; j < multiplierValue; j++) {
						SMILES+="c4cc5ccccc5cc4";
					}
					SMILES+= "c3cc12";
				} else if (classOfHydrocarbonFRSystem.equals("polyhelicene")){
					if (multiplierValue <6){
						throw new PostProcessingException("Invalid polyhelicene");
					}
					SMILES= "c1c";
					int ringOpeningCounter=2;
					for (int j = 1; j < multiplierValue; j++) {
						SMILES+="ccc" + ringClosure(ringOpeningCounter++);
					}
					SMILES+= "cccc";
					ringOpeningCounter--;
					for (int j = 2; j < multiplierValue; j++) {
						SMILES+="c" + ringClosure(ringOpeningCounter--);
					}
					SMILES+= "c12";
				}

				else{
					throw new PostProcessingException("Unknown semi-trivially named hydrocarbon fused ring system");
				}

				newGroup.addAttribute(new Attribute("value", SMILES));
				newGroup.addAttribute(new Attribute("valType", "SMILES"));
				newGroup.addAttribute(new Attribute("labels", "fusedRing"));
				newGroup.addAttribute(new Attribute("type", "ring"));
				newGroup.addAttribute(new Attribute("subType", "hydrocarbonFusedRingSystem"));
				newGroup.appendChild(multiplier.getValue() + hydrocarbonFRSystem.getValue());
				hydrocarbonFRSystem.getParent().replaceChild(hydrocarbonFRSystem, newGroup);
				multiplier.detach();
			}
			else{
				throw new PostProcessingException("Invalid semi-trivially named hydrocarbon fused ring system");
			}
		}
	}

	/**
	 * If the integer given is > 9 return %ringClosure else just returns ringClosure
	 * @param ringClosure
	 * @return
	 */
	private String ringClosure(int ringClosure) {
		if (ringClosure >9){
			return "%" +Integer.toString(ringClosure);
		}
		else{
			return Integer.toString(ringClosure);
		}
	}

	/**
	 * Prepares spiro string for processing
	 * @param text - string with spiro
	 * @return array with number of carbons in each group and associated index of spiro atom
	 */
	private int[][] getSpiroGroups(String text)
	{
		int n1 = text.indexOf('[');
		int n2 = text.indexOf(']');
		text = text.substring(n1+1, n2);
		String[] sGroups = matchDot.split(text);
		int gNum = sGroups.length;

		int[][] iGroups = new int[gNum][2]; // array of groups where number of elements and super string present

		for (int i=0; i<sGroups.length; i++)
		{
			String[] elements = matchNonDigit.split(sGroups[i]);
			if (elements.length >1)//a "superscripted" number is present
			{
				iGroups[i][0] = Integer.parseInt(elements[0]);
				iGroups[i][1] = Integer.parseInt(elements[1]);
			}
			else
			{
				iGroups[i][0] = Integer.parseInt(sGroups[i]);
				iGroups[i][1] = -1;
			}
		}

		return iGroups;
	}

	/**
	 * finds atom index in smile string corresponding to atom index in a given structure
	 * @param smile string to search in
	 * @param index index of the atom in given structure
	 * @return correspondent atom index in smile string
	 */
	private String findCSmileIndex(String smile, int index)
	{
		int cnt = 0;
		int pos = -1;
		int i;
		String sIndex;

		for (i=0; i<smile.length(); i++)
		{
			if (smile.charAt(i) == 'C') cnt++;
			if (cnt==index) { pos=i; break;}
		}

		pos++;

		if (smile.charAt(pos)=='%')
		{
			sIndex = "%";
			pos++;
			while (smile.charAt(pos)>='0' && smile.charAt(pos)<='9' && pos<smile.length())
			{
				sIndex += smile.charAt(pos);
				pos++;
			}

		}
		else sIndex = "" + smile.charAt(pos); // Can be non integer char

		if (sIndex.equals("0")) sIndex="1"; // as the first ring can have 2 indices, we need the second
		return sIndex;
	}

	/**Looks (multiplier)cyclo/spiro/cyclo tags before chain
	 * and replaces them with a group with appropriate SMILES
	 * Note that only simple spiro tags are handled at this stage i.e. not dispiro
	 * @param group A group which is potentially a chain
	 * @throws PostProcessingException
	 */
	private void processRings(Element group) throws PostProcessingException {
		Element previous = (Element)XOMTools.getPreviousSibling(group);
		if(previous != null) {
			if(previous.getLocalName().equals("spiro")){
				String text = previous.getValue();
				int[][] groups = getSpiroGroups(text);
				int curIndex = 2;

				Element multiplier =(Element)XOMTools.getPreviousSibling(previous);
				int numberOfSpiros = 1;
				if (multiplier != null && multiplier.getLocalName().equals("multiplier")){
					numberOfSpiros = Integer.parseInt(multiplier.getAttributeValue("value"));
					multiplier.detach();
				}
				int numOfOpenedBrackets = 1;

				String SMILES = "C0" + StringTools.multiplyString("C", groups[0][0]) + "01(";

				// for those molecules where no superstrings compare prefix number with curIndex.
				for (int i=1; i<groups.length; i++)
				{
					if (groups[i][1] >= 0)
					{
						String smileIndex = findCSmileIndex( SMILES, groups[i][1] );

						int pos = SMILES.indexOf(smileIndex);

						// we already had this index twice
						// the molecule has atom connecting more than 1 ring
						if (SMILES.indexOf(smileIndex, pos+1)>=0)
						{
							// insert extra index
							SMILES = SMILES.substring(0, pos+1) + curIndex + SMILES.substring(pos+1);

							// add ring in a new brackets
							String stIndex = ""+curIndex;
							if (curIndex>9) stIndex += "%" + stIndex; // check if is more than 9
							SMILES += "(" + StringTools.multiplyString("C", groups[i][0]) + stIndex + ")";
							curIndex++;
						}
						else
						{
							SMILES += StringTools.multiplyString("C", groups[i][0]) + smileIndex + ")";
						}

					}
					else if (numOfOpenedBrackets >= numberOfSpiros)
					{
						SMILES += StringTools.multiplyString("C", groups[i][0]);

						// take the number before bracket as index for smile
						// we can open more brackets, this considered in prev if
						curIndex--;
						if (curIndex>9)	SMILES += "%";
						SMILES += curIndex + ")";

						// from here start to decrease index for the following
					}
					else
					{
						SMILES += StringTools.multiplyString("C", groups[i][0]);

						if (curIndex>9)
							SMILES += "C%" + curIndex++ + "(";
						else
							SMILES += "C" + curIndex++ + "(";
						numOfOpenedBrackets++;
					}
				}
				group.addAttribute(new Attribute("value", SMILES));
				group.addAttribute(new Attribute("valType", "SMILES"));
				group.getAttribute("type").setValue("ring");
				previous.detach();
			} else if(previous.getLocalName().equals("vonBaeyer")) {
				String vonBaeyerBracket = previous.getValue();
				Element multiplier =(Element)XOMTools.getPreviousSibling(previous);
				int numberOfRings=Integer.parseInt(multiplier.getAttributeValue("value"));
				multiplier.detach();

				int alkylChainLength;
				LinkedList<String> elementSymbolArray = new LinkedList<String>();
				if (group.getAttributeValue("valType").equals("chain")){
					alkylChainLength=Integer.parseInt(group.getAttributeValue("value"));
					for (int i = 0; i < alkylChainLength; i++) {
						elementSymbolArray.add("C");
					}
				}
				else if (group.getAttributeValue("valType").equals("SMILES")){
					String smiles =group.getAttributeValue("value");
					char[] smilesArray =smiles.toCharArray();
					for (int i = 0; i < smilesArray.length; i++) {//only able to interpret the SMILES that should be in an unmodified unbranched chain
						char currentChar =smilesArray[i];
						if (currentChar == '['){
							if ( smilesArray[i +2]==']'){
								elementSymbolArray.add("[" +String.valueOf(smilesArray[i+1]) +"]");
								i=i+2;
							}
							else{
								elementSymbolArray.add("[" + String.valueOf(smilesArray[i+1]) +String.valueOf(smilesArray[i+2]) +"]");
								i=i+3;
							}
						}
						else{
							elementSymbolArray.add(String.valueOf(currentChar));
						}
					}
					alkylChainLength=elementSymbolArray.size();
				}
				else{
					throw new PostProcessingException("unexpected group valType: " + group.getAttributeValue("valType"));
				}


				int totalLengthOfBridges=0;
				int bridgeLabelsUsed=3;//start labelling from 3 upwards
				//3 and 4 will be the atoms on each end of one secondary bridge, 5 and 6 for the next etc.

				ArrayList<HashMap<String, Integer>> bridges = new ArrayList<HashMap<String, Integer>>();
				HashMap<Integer, ArrayList<Integer>> bridgeLocations = new HashMap<Integer, ArrayList<Integer>>(alkylChainLength);
				Matcher m = matchVonBaeyer.matcher(vonBaeyerBracket);
				while(m.find()) {
					String[] lengthOfBridgeArray=m.group(0).split(",");
					HashMap<String, Integer> bridge = new HashMap<String, Integer>();
					int bridgeLength =0;
					if (lengthOfBridgeArray.length > 1){//this is a secondary bridge (chain start/end locations have been specified)

						String coordinatesStr1;
						String coordinatesStr2 =lengthOfBridgeArray[1].replaceAll("\\D", "");
						String[] tempArray = lengthOfBridgeArray[0].split("\\D+");

						if (tempArray.length ==1){
							//there is some ambiguity as it has not been made obvious which number/s are supposed superscripted to be the superscripted locant
							//so we assume that it is more likely that it will be referring to an atom of label >10
							//rather than a secondary bridge of length > 10
							char[] tempCharArray = lengthOfBridgeArray[0].toCharArray();
							if (tempCharArray.length ==2){
								bridgeLength= Character.getNumericValue(tempCharArray[0]);
								coordinatesStr1= Character.toString(tempCharArray[1]);
							}
							else if (tempCharArray.length ==3){
								bridgeLength= Character.getNumericValue(tempCharArray[0]);
								coordinatesStr1=Character.toString(tempCharArray[1]) +Character.toString(tempCharArray[2]);
							}
							else if (tempCharArray.length ==4){
								bridgeLength = Integer.parseInt(Character.toString(tempCharArray[0]) +Character.toString(tempCharArray[1]));
								coordinatesStr1 = Character.toString(tempCharArray[2]) +Character.toString(tempCharArray[3]);
							}
							else{
								throw new PostProcessingException("Unsupported Von Baeyer locant description: " + m.group(0) );
							}
						}
						else{//bracket or other delimiter detected, no ambiguity!
							bridgeLength= Integer.parseInt(tempArray[0]);
							coordinatesStr1= tempArray[1];
						}

						bridge.put("Bridge Length", bridgeLength );
						int coordinates1=Integer.parseInt(coordinatesStr1);
						int coordinates2=Integer.parseInt(coordinatesStr2);
						if (coordinates1 > alkylChainLength || coordinates2 > alkylChainLength){
							throw new PostProcessingException("Indicated bridge position is not on chain: " +coordinates1 +"," +coordinates2);
						}
						if (coordinates2>coordinates1){//makes sure that bridges are built from highest coord to lowest
							int swap =coordinates1;
							coordinates1=coordinates2;
							coordinates2=swap;
						}
						if (bridgeLocations.get(coordinates1)==null){
							bridgeLocations.put(coordinates1, new ArrayList<Integer>());
						}
						if (bridgeLocations.get(coordinates2)==null){
							bridgeLocations.put(coordinates2, new ArrayList<Integer>());
						}
						bridgeLocations.get(coordinates1).add(bridgeLabelsUsed);
						bridge.put("AtomId_Larger_Label", bridgeLabelsUsed);
						bridgeLabelsUsed++;
						if (bridgeLength==0){//0 length bridge, hence want atoms with the same labels so they can join together without a bridge
							bridgeLocations.get(coordinates2).add(bridgeLabelsUsed -1);
							bridge.put("AtomId_Smaller_Label", bridgeLabelsUsed -1);
						}
						else{
							bridgeLocations.get(coordinates2).add(bridgeLabelsUsed);
							bridge.put("AtomId_Smaller_Label", bridgeLabelsUsed);
						}
						bridgeLabelsUsed++;

						bridge.put("AtomId_Larger", coordinates1);
						bridge.put("AtomId_Smaller", coordinates2);
					}
					else{
						bridgeLength= Integer.parseInt(lengthOfBridgeArray[0]);
						bridge.put("Bridge Length", bridgeLength);
					}
					totalLengthOfBridges += bridgeLength;
					bridges.add(bridge);
				}
				if (totalLengthOfBridges + 2 !=alkylChainLength ){
					throw new PostProcessingException("Disagreement between lengths of bridges and alkyl chain length");
				}
				if (numberOfRings +1 != bridges.size()){
					throw new PostProcessingException("Disagreement between number of rings and number of bridges");
				}

				String SMILES="";
				int atomCounter=1;
				int bridgeCounter=1;
				//add standard bridges
				for (HashMap<String, Integer> bridge : bridges) {
					if (bridgeCounter==1){
						SMILES += elementSymbolArray.removeFirst() +"1";
						if (bridgeLocations.get(atomCounter)!=null){
							for (Integer bridgeAtomLabel : bridgeLocations.get(atomCounter)) {
								SMILES += ringClosure(bridgeAtomLabel);
							}
						}
						SMILES += "(";
					}
					int bridgeLength =bridge.get("Bridge Length");

					for (int i = 0; i < bridgeLength; i++) {
						atomCounter++;
						SMILES +=elementSymbolArray.removeFirst();
						if (bridgeLocations.get(atomCounter)!=null){
							for (Integer bridgeAtomLabel : bridgeLocations.get(atomCounter)) {
								SMILES +=ringClosure(bridgeAtomLabel);
							}
						}
					}
					if (bridgeCounter==1){
						atomCounter++;
						SMILES += elementSymbolArray.removeFirst() +"2";
						if (bridgeLocations.get(atomCounter)!=null){
							for (Integer bridgeAtomLabel : bridgeLocations.get(atomCounter)) {
								SMILES +=ringClosure(bridgeAtomLabel);
							}
						}
					}
					if (bridgeCounter==2){
						SMILES += "1)";
					}
					if (bridgeCounter==3){
						SMILES += "2";
					}
					bridgeCounter++;
					if (bridgeCounter >3){break;}
				}

				//create list of secondary bridges that need to be added
				//0 length bridges and the 3 main bridges are dropped
				ArrayList<HashMap<String, Integer>> secondaryBridges = new ArrayList<HashMap<String, Integer>>();
				for (HashMap<String, Integer> bridge : bridges) {
					if(bridge.get("AtomId_Larger")!=null && bridge.get("Bridge Length")!=0){
						secondaryBridges.add(bridge);
					}
				}

				Comparator<HashMap<String, Integer>> sortBridges= new VonBaeyerSecondaryBridgeSort();
				Collections.sort(secondaryBridges, sortBridges);

				ArrayList<HashMap<String, Integer>> dependantSecondaryBridges;
				//add secondary bridges, recursively add dependent secondary bridges
				do{
					dependantSecondaryBridges = new ArrayList<HashMap<String, Integer>>();
					for (HashMap<String, Integer> bridge : secondaryBridges) {
						int bridgeLength =bridge.get("Bridge Length");
						if (bridge.get("AtomId_Larger") > atomCounter){
							dependantSecondaryBridges.add(bridge);
							continue;
						}
						SMILES+=".";
						for (int i = 0; i < bridgeLength; i++) {
							atomCounter++;
							SMILES +=elementSymbolArray.removeFirst();
							if (i==0){SMILES+=ringClosure(bridge.get("AtomId_Larger_Label"));}
							if (bridgeLocations.get(atomCounter)!=null){
								for (Integer bridgeAtomLabel : bridgeLocations.get(atomCounter)) {
									SMILES += ringClosure(bridgeAtomLabel);
								}
							}
						}
						SMILES+= ringClosure(bridge.get("AtomId_Smaller_Label"));
					}
					if (dependantSecondaryBridges.size() >0 && dependantSecondaryBridges.size()==secondaryBridges.size()){
						throw new PostProcessingException("Unable to resolve all dependant bridges!!!");
					}
					secondaryBridges=dependantSecondaryBridges;
				}
				while(dependantSecondaryBridges.size() > 0);

				group.addAttribute(new Attribute("value", SMILES));
				group.addAttribute(new Attribute("valType", "SMILES"));
				group.getAttribute("type").setValue("ring");
				previous.detach();
			}
			else if(previous.getLocalName().equals("cyclo")) {
				if (!group.getAttributeValue("subType").equals("heteroStem")){
					int chainlen = Integer.parseInt(group.getAttributeValue("value"));
					if (chainlen < 3){
						throw new PostProcessingException("Alkane chain too small to create a cyclo alkane: " + chainlen);
					}
					Element next = (Element)XOMTools.getNextSibling(group);
					int groupCount = ((Element)group.getParent()).getChildElements("group").size();
					Boolean conjugate =false;
					if (next!=null && next.getLocalName().equals("unsaturator")){
						if (groupCount >1 && next.getAttributeValue("value").equals("2")){
							conjugate=true;
							next.detach();
						}
					}
					else if (groupCount >1){
						conjugate =true;
					}

					String SMILES;
					if (conjugate){
						//will have conjugated double bonds as is dictated by fusion nomenclature
						SMILES = "c1" + StringTools.multiplyString("c", chainlen -1) + "1";
					}
					else{
						SMILES = "C1" + StringTools.multiplyString("C", chainlen - 1) + "1";
					}
					group.addAttribute(new Attribute("value", SMILES));
					group.addAttribute(new Attribute("valType", "SMILES"));
				}
				else{
					String smiles=group.getAttributeValue("value");
					smiles+="1";
					if (Character.isUpperCase(smiles.charAt(1))){//element is 1 letter long
						smiles= smiles.substring(0,1) +"1" + smiles.substring(1);
					}
					else{
						smiles= smiles.substring(0,2) +"1" + smiles.substring(2);
					}
					group.getAttribute("value").setValue(smiles);
				}
				group.getAttribute("type").setValue("ring");
				previous.detach();
			}
		}
	}

	/**Handles special cases in IUPAC nomenclature.
	 * Benzyl etc.
	 * @param state
	 * @param group The group to look for irregularities in.
	 */
	private void handleIrregularities(Element group, BuildState state) throws PostProcessingException {
		String groupValue =group.getValue();
		/* Benzyl, benzyloxy etc. Add a methylene */
		if(groupValue.equals("benz")) {
			Element possibleSuffix = XOMTools.getNextSiblingIgnoringCertainElements(group, new String[]{"locant", "multiplier"});
			if (possibleSuffix !=null && possibleSuffix.getLocalName().equals("suffix")) {
				group.getAttribute("value").setValue("Cc1ccccc1");
				group.getAttribute("valType").setValue("SMILES");
				group.addAttribute(new Attribute("labels", "alpha/1/2/3/4/5/6"));
			}
		}

		if(groupValue.equals("thiophen")) {//thiophenol is phenol with an O replaced with S not thiophene with a hydroxy
			Element possibleSuffix = (Element) XOMTools.getNextSibling(group);
			if (possibleSuffix !=null && possibleSuffix.getLocalName().equals("suffix")) {
				if (possibleSuffix.getValue().equals("ol")){
					throw new PostProcessingException("thiophenol has been incorrectly interpreted as thiophen, ol instead of thio, phen, ol");
				}
			}
		}

		if (groupValue.equals("ethylene")) {
			Element previous = (Element)XOMTools.getPreviousSibling(group);
			if (previous!=null && previous.getLocalName().equals("multiplier")){
				int multiplierValue = Integer.parseInt(previous.getAttributeValue("value"));
				Element possibleRoot =(Element) XOMTools.getNextSibling(group.getParent());
				if (possibleRoot==null && state.wordRule.equals("glycol")){//e.g. dodecaethylene glycol
					String smiles ="CC";
					for (int i = 1; i < multiplierValue; i++) {
						smiles+="OCC";
					}
					group.getAttribute("outIDs").setValue("1," +Integer.toString(3*(multiplierValue-1) +2));
					group.getAttribute("value").setValue(smiles);
					previous.detach();
				}
				else if (possibleRoot!=null && possibleRoot.getLocalName().equals("root")){
					Elements children = possibleRoot.getChildElements();
					if (children.size()==2){
						Element amineMultiplier =children.get(0);
						Element amine =children.get(1);
						if (amineMultiplier.getLocalName().equals("multiplier") && amine.getValue().equals("amin")){//e.g. Triethylenetetramine
							if (Integer.parseInt(amineMultiplier.getAttributeValue("value"))!=multiplierValue +1){
								throw new PostProcessingException("Invalid polyethylene amine!");
							}
							String smiles ="";
							for (int i = 0; i < multiplierValue; i++) {
								smiles+="NCC";
							}
							smiles+="N";
							group.removeAttribute(group.getAttribute("outIDs"));
							group.getAttribute("value").setValue(smiles);
							previous.detach();
							possibleRoot.detach();
							((Element)group.getParent()).setLocalName("root");
						}
					}
				}
			}
		}

		//anthrone, phenanthrone and xanthone have the one at position 9 by default
		if (groupValue.equals("anthr") || groupValue.equals("phenanthr") || groupValue.equals("xanth")) {
			Element possibleLocant = (Element) XOMTools.getPreviousSibling(group);
			if (possibleLocant==null || !possibleLocant.getLocalName().equals("locant")){//only need to give one a locant of 9 if no locant currently present
				Element possibleOne =(Element) XOMTools.getNextSibling(group);
				if (possibleOne!=null && possibleOne.getValue().equals("one")){
					Element newLocant =new Element("locant");
					newLocant.appendChild("9(10H)");//Rule C-315.2
					XOMTools.insertBefore(possibleOne, newLocant);
				}
			}
		}
	}

	/**
	 * Converts hydro substituents to properties of the next group of type ring
	 * @param group
	 * @throws PostProcessingException
	 */
	private void processHydroSubstituents(Element group) throws PostProcessingException {
		if (group.getValue().equals("hydro")){
			Element hydroSubstituent =(Element) group.getParent();
			Element targetRing =null;
			Node nextSubOrRootOrBracket = XOMTools.getNextSibling(hydroSubstituent);
			//first check adjacent substituent/root. If this is locantless then we can assume the hydro is acting as a nondetachable prefix
			Element potentialRing =((Element)nextSubOrRootOrBracket).getFirstChildElement("group");
			if (potentialRing!=null && potentialRing.getAttributeValue("type").equals("ring")){
				Element possibleLocant =(Element) XOMTools.getPreviousSibling(potentialRing, "locant");
				if (possibleLocant !=null){
					if (potentialRing.getAttribute("frontLocantsExpected")!=null){//check whether the group was expecting a locant e.g. 2-furyl
						String locantValue =StringTools.removeDashIfPresent(possibleLocant.getValue());
						String[] expectedLocants = matchComma.split(potentialRing.getAttributeValue("frontLocantsExpected"));
						for (String expectedLocant : expectedLocants) {
							if (locantValue.equals(expectedLocant)){
								targetRing =potentialRing;
								break;
							}
						}
					}
					//check whether the group is a HW system e.g. 1,3-thiazole
					if (potentialRing.getAttributeValue("subType").equals("hantzschWidman")){
						String locantValue =StringTools.removeDashIfPresent(possibleLocant.getValue());
						int locants = matchComma.split(locantValue).length;
						int heteroCount = 0;
						Element currentElem =  (Element) XOMTools.getNextSibling(possibleLocant);
						while(!currentElem.equals(potentialRing)){
							if(currentElem.getLocalName().equals("heteroatom")) {
								heteroCount++;
							} else if (currentElem.getLocalName().equals("multiplier")){
								heteroCount += Integer.parseInt(currentElem.getAttributeValue("value")) -1;
							}
							currentElem = (Element)XOMTools.getNextSibling(currentElem);
						}
						if (heteroCount==locants){//number of locants must match number
							targetRing =potentialRing;
						}
					}
				}
				else{
					targetRing =potentialRing;
				}
			}

			//that didn't match so the hydro appears to be a detachable prefix. detachable prefixes attach in preference to the rightmost applicable group so search any remaining substituents/roots from right to left
			if (targetRing ==null){
				Element nextSubOrRootOrBracketfromLast = (Element) hydroSubstituent.getParent().getChild(hydroSubstituent.getParent().getChildCount()-1);//the last sibling
				while (!nextSubOrRootOrBracketfromLast.equals(hydroSubstituent)){
					potentialRing = nextSubOrRootOrBracketfromLast.getFirstChildElement("group");
					if (potentialRing!=null && potentialRing.getAttributeValue("type").equals("ring")){
						targetRing =potentialRing;
						break;
					}
					else{
						nextSubOrRootOrBracketfromLast = (Element) XOMTools.getPreviousSibling(nextSubOrRootOrBracketfromLast);
					}
				}
			}
			if (targetRing ==null){
				throw new PostProcessingException("Cannot find ring for hydro substituent to apply to");
			}
			//move the children of the hydro substituent, and change the group tag into a hydro tag
			Elements children =hydroSubstituent.getChildElements();
			for (int i = children.size()-1; i >=0 ; i--) {
				Element child =children.get(i);
				if (!child.getLocalName().equals("hyphen")){
					child.detach();
					if (child.getLocalName().equals("group")){
						child =new Element("hydro");
						child.appendChild("hydro");
					}
					targetRing.getParent().insertChild(child, 0);
				}
			}
			hydroSubstituent.detach();
		}
	}

	/**Finds matching open and close brackets, and places the
	 * elements contained within in a big &lt;bracket&gt; element.
	 *
	 * @param substituentsAndRoot: The substituent/root elements at the current level of the tree
	 * @return Whether the method did something, and so needs to be called again.
	 * @throws PostProcessingException
	 */
	private boolean findAndStructureBrackets(List<Element> substituentsAndRoot) throws PostProcessingException {
		int blevel = 0;
		Element openBracket = null;
		Element closeBracket = null;
		for (Element sub : substituentsAndRoot) {
			Elements children = sub.getChildElements();
			for(int i=0; i<children.size(); i++) {
				Element child = children.get(i);
				if(child.getLocalName().equals("openbracket")) {
					if(openBracket == null) {
						openBracket = child;
					}
					blevel++;
				} else if (child.getLocalName().equals("closebracket")) {
					blevel--;
					if(blevel == 0) {
						closeBracket = child;
						Element bracket = structureBrackets(openBracket, closeBracket);
						while(findAndStructureBrackets(XOMTools.getDescendantElementsWithTagName(bracket, "substituent")));
						return true;
					}
				}
			}
		}
		if (blevel!=0){
			throw new PostProcessingException("Matching closing bracket not found!");
		}
		return false;
	}

	/**Places the elements in substituents containing/between an open and close bracket
	 * in a &lt;bracket&gt; tag.
	 *
	 * @param openBracket The open bracket element
	 * @param closeBracket The close bracket element
	 * @return The bracket element thus created.
	 */
	private Element structureBrackets(Element openBracket, Element closeBracket) {
		Element bracket = new Element("bracket");
		XOMTools.insertBefore(openBracket.getParent(), bracket);
		/* Pick up everything in the substituent before the bracket*/
		while(!openBracket.getParent().getChild(0).equals(openBracket)) {
			Node n = openBracket.getParent().getChild(0);
			n.detach();
			bracket.appendChild(n);
		}
		/* Pick up all nodes from the one with the open bracket,
		 * to the one with the close bracket, inclusive.
		 */
		Node currentNode = openBracket.getParent();
		while(!currentNode.equals(closeBracket.getParent())) {
			Node nextNode = XOMTools.getNextSibling(currentNode);
			currentNode.detach();
			bracket.appendChild(currentNode);
			currentNode = nextNode;
		}
		currentNode.detach();
		bracket.appendChild(currentNode);
		/* Pick up nodes after the close bracket */
		currentNode = XOMTools.getNextSibling(closeBracket);
		while(currentNode != null) {
			Node nextNode = XOMTools.getNextSibling(currentNode);
			currentNode.detach();
			bracket.appendChild(currentNode);
			currentNode = nextNode;
		}
		openBracket.detach();
		closeBracket.detach();

		return bracket;
	}

	/**
	 * Moves substituents into new words if it appears that a space was omitted in the input name
	 * @param state
	 * @param elem
	 */
	private void addOmittedSpaces(BuildState state, Element elem) {
		if (state.wordRule.equals("divalentLiteralFunctionalGroup")){
			List<Element> substituentWords = XOMTools.getChildElementsWithTagNameAndAttribute(elem, "word", "type", "substituent");
			if (substituentWords.size()==1){//potentially been "wrongly" interpreted e.g. ethylmethyl ketone is more likely to mean ethyl methyl ketone
				Elements children  =substituentWords.get(0).getChildElements();
				if (children.size()==2){
					Element firstChildOfFirstSubstituent =(Element)children.get(0).getChild(0);
					//rule out correct usage e.g. diethyl ether and locanted substituents e.g. 2-methylpropyl ether
					if (!firstChildOfFirstSubstituent.getLocalName().equals("locant") && !firstChildOfFirstSubstituent.getLocalName().equals("multiplier")){
						Element subToMove =children.get(1);
						subToMove.detach();
						Element newWord =new Element("word");
						newWord.addAttribute(new Attribute("type", "substituent"));
						newWord.appendChild(subToMove);
						XOMTools.insertAfter(substituentWords.get(0), newWord);
					}
				}
			}
		}
	}

}

