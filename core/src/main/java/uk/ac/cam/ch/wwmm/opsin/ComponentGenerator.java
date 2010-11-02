package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import uk.ac.cam.ch.wwmm.opsin.WordRules.WordRule;
import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

import nu.xom.Attribute;
import nu.xom.Element;
import nu.xom.Elements;
import nu.xom.Node;

/**Does destructive procedural parsing on parser results.
 *
 * @author ptc24/dl387
 *
 */
class ComponentGenerator {

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
	private final Pattern matchNumberLocantsOnlyFusionBracket = Pattern.compile("\\[\\d+(,\\d+)*\\]");
	private final Pattern matchIndicatedHydrogen =Pattern.compile("(\\d+[a-z]?'*)H");
	private final Pattern matchIndicatedHydrogenBracket =Pattern.compile("[\\[\\(\\{][^\\[\\(\\{]*H[\\]\\)\\}]");
	private final Pattern matchCommaOrDot =Pattern.compile("[\\.,]");
	private final Pattern matchAnnulene = Pattern.compile("[\\[\\(\\{]([1-9]\\d*)[\\]\\)\\}]annulen");
	private final String elementSymbols ="(?:He|Li|Be|B|C|N|O|F|Ne|Na|Mg|Al|Si|P|S|Cl|Ar|K|Ca|Sc|Ti|V|Cr|Mn|Fe|Co|Ni|Cu|Zn|Ga|Ge|As|Se|Br|Kr|Rb|Sr|Y|Zr|Nb|Mo|Tc|Ru|Rh|Pd|Ag|Cd|In|Sn|Sb|Te|I|Xe|Cs|Ba|La|Ce|Pr|Nd|Pm|Sm|Eu|Gd|Tb|Dy|Ho|Er|Tm|Yb|Lu|Hf|Ta|W|Re|Os|Ir|Pt|Au|Hg|Tl|Pb|Bi|Po|At|Rn|Fr|Ra|Ac|Th|Pa|U|Np|Pu|Am|Cm|Bk|Cf|Es|Fm|Md|No|Lr|Rf|Db|Sg|Bh|Hs|Mt|Ds)";
	private final Pattern matchStereochemistry = Pattern.compile("(.*?)(SR|RS|[RSEZrsez])");
	private final Pattern matchStar = Pattern.compile("\\*");
	private final Pattern matchRS = Pattern.compile("[RSrs]");
	private final Pattern matchEZ = Pattern.compile("[EZez]");
	private final Pattern matchLambdaConvention = Pattern.compile("(\\S+)?lambda\\D*(\\d+)\\D*");
	private final Pattern matchComma =Pattern.compile(",");
	private final Pattern matchSemiColon =Pattern.compile(";");
	private final Pattern matchHdigit =Pattern.compile("H\\d");
	private final Pattern matchNonDigit =Pattern.compile("\\D+");
	private final Pattern matchSuperscriptedLocant = Pattern.compile("(" + elementSymbols +"'*).*(\\d+[a-z]?'*).*");
	private final Pattern matchIUPAC2004ElementLocant = Pattern.compile("(\\d+'*)-(" + elementSymbols +"'*)");
	private final Pattern matchElementSymbol = Pattern.compile("[A-Z][a-z]?");
	private final Pattern matchInlineSuffixesThatAreAlsoGroups = Pattern.compile("carbonyl|oxy|sulfenyl|sulfinyl|sulfonyl|selenenyl|seleninyl|selenonyl|tellurenyl|tellurinyl|telluronyl");

	/** The master method, processes a parse result destructively adding semantic information by processing the various micro syntaxes .
	 *
	 * @param moleculeEl The element to process.
	 * @param state
	 * @return
	 * @throws Exception
	 */
	void process(Element moleculeEl, BuildState state) throws Exception {
		List<Element> substituentsAndRoot = XOMTools.getDescendantElementsWithTagNames(moleculeEl, new String[]{SUBSTITUENT_EL, ROOT_EL});

		for (Element subOrRoot: substituentsAndRoot) {
			/* Throws exceptions for occurrences that are ambiguous and this parse has picked the incorrect interpretation */
			resolveAmbiguities(subOrRoot);

			processLocants(subOrRoot);
			convertOrthoMetaParaToLocants(subOrRoot);
			formAlkaneStemsFromComponents(subOrRoot);
			processAlkaneStemModifications(subOrRoot);//e.g. tert-butyl
			processHeterogenousHydrides(subOrRoot);//e.g. tetraphosphane, disiloxane
			processIndicatedHydrogens(subOrRoot);
			processStereochemistry(subOrRoot);
			processInfixes(subOrRoot);
			processSuffixPrefixes(subOrRoot);
			processLambdaConvention(subOrRoot);
		}
		List<Element> groups =  XOMTools.getDescendantElementsWithTagName(moleculeEl, GROUP_EL);

		
		/* Converts open/close bracket elements to bracket elements and
		 *  places the elements inbetween within the newly created bracket */
		while(findAndStructureBrackets(substituentsAndRoot));
		
		for (Element subOrRoot: substituentsAndRoot) {
			processHydroCarbonRings(subOrRoot);
			handleSuffixIrregularities(subOrRoot);//handles quinone -->dioxo
		}
		for (Element group : groups) {
			detectAlkaneFusedRingBridges(group);
			processRings(group);//processes cyclo, von baeyer and spiro tokens
			handleGroupIrregularities(group);//handles benzyl, diethylene glycol, phenanthrone and other awkward bits of nomenclature
		}

		addOmittedSpaces(moleculeEl);//e.g. change ethylmethyl ether to ethyl methyl ether
	}

	/**
	 * Resolves common ambiguities e.g. tetradeca being 4x10carbon chain rather than 14carbon chain
	 * @param subOrRoot
	 * @throws ComponentGenerationException
	 */
	private void resolveAmbiguities(Element subOrRoot) throws ComponentGenerationException {
		List<Element> multipliers = XOMTools.getChildElementsWithTagName(subOrRoot, MULTIPLIER_EL);
		for (Element apparentMultiplier : multipliers) {
			if (!BASIC_TYPE_VAL.equals(apparentMultiplier.getAttributeValue(TYPE_ATR)) && !VONBAEYER_TYPE_VAL.equals(apparentMultiplier.getAttributeValue(TYPE_ATR))){
				continue;
			}
			int multiplierNum = Integer.parseInt(apparentMultiplier.getAttributeValue(VALUE_ATR));
			Element nextEl = (Element)XOMTools.getNextSibling(apparentMultiplier);
			if (multiplierNum >=3){//detects ambiguous use of things like tetradeca
				if(nextEl !=null){
					if (nextEl.getLocalName().equals(ALKANESTEMCOMPONENT)){//can ignore the trivial alkanes as ambiguity does not exist for them
						int alkaneChainLength = Integer.parseInt(nextEl.getAttributeValue(VALUE_ATR));
						if (alkaneChainLength >=10 && alkaneChainLength > multiplierNum){
							Element isThisALocant =(Element)XOMTools.getPreviousSibling(apparentMultiplier);
							if (isThisALocant == null ||
									!isThisALocant.getLocalName().equals(LOCANT_EL) ||
									matchComma.split(isThisALocant.getValue()).length != multiplierNum){
								throw new ComponentGenerationException(apparentMultiplier.getValue() + nextEl.getValue() +" should not have been lexed as two tokens!");
							}
						}
					}
				}
			}

			if (multiplierNum >=4 && nextEl !=null && nextEl.getLocalName().equals(HYDROCARBONFUSEDRINGSYSTEM_EL)&& nextEl.getValue().equals("phen")){//deals with tetra phen yl vs tetraphen yl
				Element possibleSuffix = (Element) XOMTools.getNextSibling(nextEl);
				if (possibleSuffix!=null){//null if not used as substituent
					String multiplierAndGroup =apparentMultiplier.getValue() + nextEl.getValue();
					if (possibleSuffix.getValue().equals("yl")){
						throw new ComponentGenerationException(multiplierAndGroup +" should not have been lexed as one token!");
					}
					Element isThisALocant =(Element)XOMTools.getPreviousSibling(apparentMultiplier);
					if (isThisALocant != null && isThisALocant.getLocalName().equals(LOCANT_EL) && matchComma.split(isThisALocant.getValue()).length == multiplierNum){
						throw new ComponentGenerationException(multiplierAndGroup +" should not have been lexed as one token!");
					}
				}
			}
			if (multiplierNum > 4 && !apparentMultiplier.getValue().endsWith("a")){//disambiguate pent oxy and the like. Assume it means pentanoxy rather than 5 oxys
				if (nextEl !=null && nextEl.getLocalName().equals(GROUP_EL)&& matchInlineSuffixesThatAreAlsoGroups.matcher(nextEl.getValue()).matches()){
					throw new ComponentGenerationException(apparentMultiplier.getValue() + nextEl.getValue() +" should have been lexed as [alkane stem, inline suffix], not [multiplier, group]!");
				}
			}
		}

		List<Element> fusions = XOMTools.getChildElementsWithTagName(subOrRoot, FUSION_EL);
		for (Element fusion : fusions) {
			String fusionText = fusion.getValue();
			if (matchNumberLocantsOnlyFusionBracket.matcher(fusionText).matches()){
				Element possibleHWRing =(Element) XOMTools.getNextSiblingIgnoringCertainElements(fusion, new String[]{MULTIPLIER_EL, HETEROATOM_EL});
				if (possibleHWRing !=null && HANTZSCHWIDMAN_SUBTYPE_VAL.equals(possibleHWRing.getAttributeValue(SUBTYPE_ATR))){
					int heteroCount = 0;
					int multiplierValue = 1;
					Element currentElem = (Element) XOMTools.getNextSibling(fusion);
					while(currentElem != null && !currentElem.getLocalName().equals(GROUP_EL)){
						if(currentElem.getLocalName().equals(HETEROATOM_EL)) {
							heteroCount+=multiplierValue;
							multiplierValue =1;
						} else if (currentElem.getLocalName().equals(MULTIPLIER_EL)){
							multiplierValue = Integer.parseInt(currentElem.getAttributeValue(VALUE_ATR));
						}
						currentElem = (Element)XOMTools.getNextSibling(currentElem);
					}
					String[] locants = matchComma.split(fusionText.substring(1, fusionText.length()-1));
					if (locants.length == heteroCount){
						boolean foundLocantNotInHwSystem =false;
						for (String locant : locants) {
							if (Integer.parseInt(locant) > (possibleHWRing.getAttributeValue(VALUE_ATR).length()-2)){
								foundLocantNotInHwSystem =true;
							}
						}
						if (!foundLocantNotInHwSystem){
							throw new ComponentGenerationException("This fusion bracket is in fact more likely to be a description of the locants of a HW ring");
						}
					}
				}
			}
		}
	}
	
	
	/**
	 * Removes hyphens from the end of locants if present
	 * Looks for locants of the form number-letter and converts them to letternumber
	 * e.g. 1-N becomes N1. 1-N is the IUPAC 2004 recommendation, N1 is the previous recommendation
	 * Strips indicated hydrogen out of locants
	 * 
	 * @param subOrRoot
	 * @throws ComponentGenerationException 
	 */
	private void processLocants(Element subOrRoot) throws ComponentGenerationException {
		List<Element> locants = XOMTools.getChildElementsWithTagName(subOrRoot, LOCANT_EL);
		for (Element locant : locants) {
			String[] individualLocantText = matchComma.split(StringTools.removeDashIfPresent(locant.getValue()));
			for (int i = 0; i < individualLocantText.length; i++) {
				String locantText =individualLocantText[i];
				if (locantText.contains("-")){//this checks should avoid having to do the regex match in all cases as locants shouldn't contain -
					Matcher m= matchIUPAC2004ElementLocant.matcher(locantText);
					if (m.matches()){
						individualLocantText[i] = m.group(2) +m.group(1);
					}
					else{
						throw new ComponentGenerationException("Unexpected hyphen in locantText");
					}
				}
				else if (Character.isLetter(locantText.charAt(0))){
					Matcher m =  matchSuperscriptedLocant.matcher(locantText);//remove indications of superscript as the fact a locant is superscripted can be determined from context e.g. N~1~ ->N1
					if (m.matches()){
						individualLocantText[i] = m.group(1) +m.group(2);
					}
					else if (locantText.length()>=3 && Character.isLetter(locantText.charAt(1)) && Character.isLetter(locantText.charAt(2))){//convert greeks to lower case
						individualLocantText[i] = locantText.toLowerCase();
					}
				}
			}
			String locantText = StringTools.arrayToString(individualLocantText, ",");
			//If the indicatedHydrogen has been specified create a tag for it and remove it from the list of locants
			//e.g. 1(9H),5,7 -->indicatedHydrogen tag value (9H) and 1,5,7
			//can get as complicated as 1,2(2H,7H)
			Matcher matches =matchIndicatedHydrogen.matcher(locantText);
			if (matches.find()){
				do {
					Element indicatedHydrogenElement=new Element(INDICATEDHYDROGEN_EL);
					indicatedHydrogenElement.addAttribute(new Attribute(LOCANT_ATR, matches.group(1)));
					XOMTools.insertBefore(locant, indicatedHydrogenElement);
				}
				while (matches.find());
				locant.addAttribute(new Attribute(TYPE_ATR, ADDEDHYDROGENLOCANT_TYPE_VAL));
				locantText =matchIndicatedHydrogenBracket.matcher(locantText).replaceAll("");
			}
			XOMTools.setTextChild(locant, locantText);

			Element afterLocants = (Element)XOMTools.getNextSibling(locant);
			if(afterLocants == null){
				throw new ComponentGenerationException("Nothing after locant tag: " + locant.toXML());
			}
		}
	}

	/**Converts ortho/meta/para into locants
	 * Depending on context para, for example, will either become para or 1,para
	 *
	 * @param subOrRoot
	 * @throws ComponentGenerationException
	 */
	private void convertOrthoMetaParaToLocants(Element subOrRoot) throws ComponentGenerationException{
		List<Element> ompLocants = XOMTools.getChildElementsWithTagName(subOrRoot, ORTHOMETAPARA_EL);
		for (Element ompLocant : ompLocants) {
			String locantText = ompLocant.getValue();
			String firstChar = locantText.substring(0, 1);
			Element afterOmpLocant = (Element)XOMTools.getNextSibling(ompLocant);
			ompLocant.setLocalName(LOCANT_EL);
			ompLocant.removeChildren();
			ompLocant.addAttribute(new Attribute(TYPE_ATR, ORTHOMETAPARA_TYPE_VAL));
			if(afterOmpLocant.getLocalName().equals(MULTIPLIER_EL) || (afterOmpLocant.getAttribute(OUTIDS_ATR)!=null && matchComma.split(afterOmpLocant.getAttributeValue(OUTIDS_ATR)).length>1) ) {
				if ("o".equalsIgnoreCase(firstChar)){
					ompLocant.appendChild("1,ortho");
				}
				else if ("m".equalsIgnoreCase(firstChar)){
					ompLocant.appendChild("1,meta");
				}
				else if ("p".equalsIgnoreCase(firstChar)){
					ompLocant.appendChild("1,para");
				}
				else{
					throw new ComponentGenerationException(locantText + " was not identified as being either ortho, meta or para but according to the chemical grammar it should of been");
				}
			}
			else{
				if ("o".equalsIgnoreCase(firstChar)){
					ompLocant.appendChild("ortho");
				}
				else if ("m".equalsIgnoreCase(firstChar)){
					ompLocant.appendChild("meta");
				}
				else if ("p".equalsIgnoreCase(firstChar)){
					ompLocant.appendChild("para");
				}
				else{
					throw new ComponentGenerationException(locantText + " was not identified as being either ortho, meta or para but according to the chemical grammar it should of been");
				}
			}
		}
	}
	
	/**
	 * Processes adjacent alkane stem component elements into a single alkaneStem group element with the appropriate SMILES
	 * e.g. dodecane would be "do" value=2 and "dec" value=10 -->alkaneStem with 12 carbons
	 * 
	 * @param subOrRoot
	 * @throws ComponentGenerationException 
	 */
	private void formAlkaneStemsFromComponents(Element subOrRoot) throws ComponentGenerationException {
		LinkedList<Element> alkaneStemComponents =new LinkedList<Element>(XOMTools.getChildElementsWithTagName(subOrRoot, ALKANESTEMCOMPONENT));
		while(!alkaneStemComponents.isEmpty()){
			Element alkaneStemComponent = alkaneStemComponents.removeFirst();
			int alkaneChainLength =0;
			StringBuilder alkaneName = new StringBuilder();
			alkaneChainLength += Integer.parseInt(alkaneStemComponent.getAttributeValue(VALUE_ATR));
			alkaneName.append(alkaneStemComponent.getValue());
			while (!alkaneStemComponents.isEmpty() && XOMTools.getNextSibling(alkaneStemComponent)==alkaneStemComponents.get(0)) {
				alkaneStemComponent.detach();
				alkaneStemComponent = alkaneStemComponents.removeFirst();
				alkaneChainLength += Integer.parseInt(alkaneStemComponent.getAttributeValue(VALUE_ATR));
				alkaneName.append(alkaneStemComponent.getValue());
			}
			Element alkaneStem = new Element(GROUP_EL);
			alkaneStem.appendChild(alkaneName.toString());
			alkaneStem.addAttribute(new Attribute(TYPE_ATR, CHAIN_TYPE_VAL));
			alkaneStem.addAttribute(new Attribute(SUBTYPE_ATR, ALKANESTEM_SUBTYPE_VAL));
			alkaneStem.addAttribute(new Attribute(VALTYPE_ATR, SMILES_VALTYPE_VAL));
			alkaneStem.addAttribute(new Attribute(VALUE_ATR, StringTools.multiplyString("C", alkaneChainLength)));
			alkaneStem.addAttribute(new Attribute(USABLEASJOINER_ATR, "yes"));
			StringBuilder labels = new StringBuilder();
			for (int i=1; i<alkaneChainLength; i++) {
				labels.append(i);
				labels.append("/");
			}
			labels.append(alkaneChainLength);
			alkaneStem.addAttribute(new Attribute(LABELS_ATR, labels.toString()));
			XOMTools.insertAfter(alkaneStemComponent, alkaneStem);
			alkaneStemComponent.detach();
		}
	}
	
	/**
	 * Applies the traditional alkane modifiers: iso, tert, sec, neo by modifying the alkane chain's SMILES
	 * 
	 * @param subOrRoot
	 * @throws ComponentGenerationException 
	 */
	private void processAlkaneStemModifications(Element subOrRoot) throws ComponentGenerationException {
		Elements alkaneStemModifiers = subOrRoot.getChildElements(ALKANESTEMMODIFIER_EL);
		for(int i=0;i<alkaneStemModifiers.size();i++) {
			Element alkaneStemModifier =alkaneStemModifiers.get(i);
			Element alkane = (Element) XOMTools.getNextSibling(alkaneStemModifier);
			if (alkane ==null || !CHAIN_TYPE_VAL.equals(alkane.getAttributeValue(TYPE_ATR))
					|| !ALKANESTEM_SUBTYPE_VAL.equals(alkane.getAttributeValue(SUBTYPE_ATR))){
				throw new ComponentGenerationException("OPSIN Bug: AlkaneStem not found after alkaneStemModifier");
			}
			String type;
			if (alkaneStemModifier.getAttribute(VALUE_ATR)!=null){
				type = alkaneStemModifier.getAttributeValue(VALUE_ATR);//identified by token;
			}
			else{
				if (alkaneStemModifier.getValue().equals("n-")){
					type="normal";
				}
				else if (alkaneStemModifier.getValue().equals("i-")){
					type="iso";
				}
				else if (alkaneStemModifier.getValue().equals("s-")){
					type="sec";
				}
				else{
					throw new ComponentGenerationException("Unrecognised alkaneStem modifier");
				}
			}
			alkaneStemModifier.detach();
			if (type.equals("normal")){
				continue;//do nothing
			}
			int chainLength = alkane.getAttributeValue(VALUE_ATR).length();
			boolean suffixPresent = subOrRoot.getChildElements(SUFFIX_EL).size() > 0;
			String smiles;
			if (type.equals("tert")){
				if (chainLength <4){
					throw new ComponentGenerationException("ChainLength to small for tert modifier, required minLength 4. Found: " +chainLength);
				}
				if (chainLength >8){
					throw new ComponentGenerationException("Interpretation of tert on an alkane chain of length: " + chainLength +" is ambiguous");
				}
				if (chainLength ==8){
					smiles = "C(C)(C)CC(C)(C)C";
				}
				else{
					smiles ="C(C)(C)C" + StringTools.multiplyString("C", chainLength-4);
				}
			}
			else if (type.equals("iso")){
				if (chainLength <3){
					throw new ComponentGenerationException("ChainLength to small for iso modifier, required minLength 3. Found: " +chainLength);
				}
				if (chainLength==3 && !suffixPresent){
					throw new ComponentGenerationException("iso has no meaning without a suffix on an alkane chain of length 3");
				}
				smiles =StringTools.multiplyString("C", chainLength-3) +"C(C)C";
			}
			else if (type.equals("sec")){
				if (chainLength <3){
					throw new ComponentGenerationException("ChainLength to small for sec modifier, required minLength 3. Found: " +chainLength);
				}
				if (!suffixPresent){
					throw new ComponentGenerationException("sec has no meaning without a suffix on an alkane chain");
				}
				smiles ="C(C)C" + StringTools.multiplyString("C", chainLength-3);
			}
			else if (type.equals("neo")){
				if (chainLength <5){
					throw new ComponentGenerationException("ChainLength to small for neo modifier, required minLength 5. Found: " +chainLength);
				}
				smiles = StringTools.multiplyString("C", chainLength-5) + "CC(C)(C)C";
			}
			else{
				throw new ComponentGenerationException("Unrecognised alkaneStem modifier");
			}
			alkane.getAttribute(VALUE_ATR).setValue(smiles);
			alkane.removeAttribute(alkane.getAttribute(USABLEASJOINER_ATR));
			alkane.getAttribute(LABELS_ATR).setValue(NONE_LABELS_VAL);
		}
	}

	/**Form heterogeneous hydrides/substituents
	 * These are chains of one heteroatom or alternating heteroatoms and are expressed using SMILES
	 * They are typically treated in an analogous way to alkanes
	 * @param subOrRoot The root/substituents
	 * @throws ComponentGenerationException 
	 */
	private void processHeterogenousHydrides(Element subOrRoot) throws ComponentGenerationException  {
		List<Element> multipliers = XOMTools.getChildElementsWithTagName(subOrRoot, MULTIPLIER_EL);
		for (int i = 0; i < multipliers.size(); i++) {
			Element m = multipliers.get(i);
			if (m.getAttributeValue(TYPE_ATR).equals(GROUP_TYPE_VAL)){
				continue;
			}
			int mvalue = Integer.parseInt(m.getAttributeValue(VALUE_ATR));
			Element multipliedElem = (Element)XOMTools.getNextSibling(m);

			if(multipliedElem.getLocalName().equals(GROUP_EL) &&
					multipliedElem.getAttribute(SUBTYPE_ATR)!=null &&
					multipliedElem.getAttributeValue(SUBTYPE_ATR).equals(HETEROSTEM_SUBTYPE_VAL)) {
				
				Element possiblyALocant = (Element)XOMTools.getPreviousSibling(m);//detect rare case where multiplier does not mean form a chain of heteroatoms e.g. something like 1,2-disulfanylpropane
				if(possiblyALocant !=null && possiblyALocant.getLocalName().equals(LOCANT_EL)&& mvalue==matchComma.split(possiblyALocant.getValue()).length){
					Element suffix =(Element) XOMTools.getNextSibling(multipliedElem, SUFFIX_EL);
					if (suffix !=null && suffix.getAttributeValue(TYPE_ATR).equals(INLINE_TYPE_VAL)){
						Element possibleMultiplier = (Element) XOMTools.getPreviousSibling(suffix);
						if (!possibleMultiplier.getLocalName().equals(MULTIPLIER_EL)){//NOT something like 3,3'-diselane-1,2-diyl
							continue;
						}
					}
				}
				
				//chain of heteroatoms
				String smiles=multipliedElem.getAttributeValue(VALUE_ATR);
				multipliedElem.getAttribute(VALUE_ATR).setValue(StringTools.multiplyString(smiles, mvalue));
				m.detach();
				multipliers.remove(i--);
			}
		}
		for (Element m : multipliers) {
			Element multipliedElem = (Element)XOMTools.getNextSibling(m);
			if(multipliedElem.getLocalName().equals(HETEROATOM_EL)){
				Element possiblyAnotherHeteroAtom = (Element)XOMTools.getNextSibling(multipliedElem);
				if (possiblyAnotherHeteroAtom !=null && possiblyAnotherHeteroAtom.getLocalName().equals(HETEROATOM_EL)){
					Element possiblyAnUnsaturator = XOMTools.getNextSiblingIgnoringCertainElements(possiblyAnotherHeteroAtom, new String[]{LOCANT_EL, MULTIPLIER_EL});//typically ane but can be ene or yne e.g. triphosphaza-1,3-diene
					if (possiblyAnUnsaturator !=null && possiblyAnUnsaturator.getLocalName().equals(UNSATURATOR_EL)){
						//chain of alternating heteroatoms
						if (possiblyAnUnsaturator.getAttributeValue(VALUE_ATR).equals("1")){
							checkForAmbiguityWithHWring(multipliedElem.getAttributeValue(VALUE_ATR), possiblyAnotherHeteroAtom.getAttributeValue(VALUE_ATR));
						}
						int mvalue = Integer.parseInt(m.getAttributeValue(VALUE_ATR));
						String smiles="";
						Element possiblyARingFormingEl = (Element)XOMTools.getPreviousSibling(m);
						boolean heteroatomChainWillFormARing =false;
						if (possiblyARingFormingEl!=null && (possiblyARingFormingEl.getLocalName().equals(CYCLO_EL) || possiblyARingFormingEl.getLocalName().equals(VONBAEYER_EL) || possiblyARingFormingEl.getLocalName().equals(SPIRO_EL))){
							heteroatomChainWillFormARing=true;
							//will be cyclised later.
							for (int j = 0; j < mvalue; j++) {
								smiles+=possiblyAnotherHeteroAtom.getAttributeValue(VALUE_ATR);
								smiles+=multipliedElem.getAttributeValue(VALUE_ATR);
							}
						}
						else{
							for (int j = 0; j < mvalue -1; j++) {
								smiles+=multipliedElem.getAttributeValue(VALUE_ATR);
								smiles+=possiblyAnotherHeteroAtom.getAttributeValue(VALUE_ATR);
							}
							smiles+=multipliedElem.getAttributeValue(VALUE_ATR);
						}
						smiles = matchHdigit.matcher(smiles).replaceAll("H?");//hydrogen count will be determined by standard valency
						multipliedElem.detach();

						Element addedGroup=new Element(GROUP_EL);
						addedGroup.addAttribute(new Attribute(VALUE_ATR, smiles));
						addedGroup.addAttribute(new Attribute(VALTYPE_ATR, SMILES_VALTYPE_VAL));
						addedGroup.addAttribute(new Attribute(TYPE_ATR, CHAIN_TYPE_VAL));
						addedGroup.addAttribute(new Attribute(SUBTYPE_ATR, HETEROSTEM_SUBTYPE_VAL));
						if (!heteroatomChainWillFormARing){
							addedGroup.addAttribute(new Attribute(USABLEASJOINER_ATR, "yes"));
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

	/**
	 * Checks that the first heteroatom is lower priority than the second
	 * If it is higher priority than the second then the ordering is that which is expected for a Hantzch-widman ring
	 * @param firstHeteroatom
	 * @param secondHeteroatom
	 * @throws ComponentGenerationException 
	 */
	private void checkForAmbiguityWithHWring(String firstHeteroAtomSMILES, String secondHeteroAtomSMILES) throws ComponentGenerationException {
		Matcher m = matchElementSymbol.matcher(firstHeteroAtomSMILES);
		if (!m.find()){
			throw new ComponentGenerationException("Failed to extract element from heteroatom");
		}
		String atom1Element = m.group();
		
		m = matchElementSymbol.matcher(secondHeteroAtomSMILES);
		if (!m.find()){
			throw new ComponentGenerationException("Failed to extract element from heteroatom");
		}
		String atom2Element = m.group();
		if (AtomProperties.elementToHwPriority.get(atom1Element)> AtomProperties.elementToHwPriority.get(atom2Element)){
			throw new ComponentGenerationException("Hantzch-widman ring misparsed as a heterogeneous hydride with alternating atoms");
		}
	}

	/** Handle 1H- in 1H-pyrrole etc.
	 *
	 * @param elem The substituent/root to looks for indicated hydrogens in.
	 */
	private void processIndicatedHydrogens(Element elem) {
		Elements hydrogens = elem.getChildElements(HYDROGEN_EL);
		for(int i=0;i<hydrogens.size();i++) {
			Element hydrogen = hydrogens.get(i);
			String txt = hydrogen.getValue();
			String[] hydrogenLocants =matchComma.split(txt);
            for (String hydrogenLocant : hydrogenLocants) {//TODO this should have an else throw exception clause? maybe should employ removedashifpresent
                if (hydrogenLocant.endsWith("H-")) {
                    Element newHydrogenElement = new Element(HYDROGEN_EL);
                    newHydrogenElement.addAttribute(new Attribute(LOCANT_ATR, hydrogenLocant.substring(0, hydrogenLocant.length() - 2)));
                    XOMTools.insertAfter(hydrogen, newHydrogenElement);
                } else if (hydrogenLocant.endsWith("H")) {
                    Element newHydrogenElement = new Element(HYDROGEN_EL);
                    newHydrogenElement.addAttribute(new Attribute(LOCANT_ATR, hydrogenLocant.substring(0, hydrogenLocant.length() - 1)));
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
	 * @throws ComponentGenerationException
	 */
	private void processStereochemistry(Element elem) throws ComponentGenerationException {
		Elements stereoChemistryElements = elem.getChildElements(STEREOCHEMISTRY_EL);
		for(int i=0;i<stereoChemistryElements.size();i++) {
			Element stereoChemistryElement = stereoChemistryElements.get(i);
			if (stereoChemistryElement.getAttributeValue(TYPE_ATR).equals(STEREOCHEMISTRYBRACKET_TYPE_VAL)){
				String txt = stereoChemistryElement.getValue();
				if (txt.startsWith("rel-")){
					txt = txt.substring(4);
				}
				Matcher starMatcher = matchStar.matcher(txt);
				txt = starMatcher.replaceAll("");
				if (!txt.startsWith("rac-")){
					txt =txt.substring(1, txt.length()-1);//remove opening and closing bracket.
					String[] stereoChemistryDescriptors = matchComma.split(txt);
                    for (String stereoChemistryDescriptor : stereoChemistryDescriptors) {
                        if (stereoChemistryDescriptor.length() > 1) {
                            Matcher m = matchStereochemistry.matcher(stereoChemistryDescriptor);
                            if (m.matches()){
                            	if (!m.group(2).equals("RS") && !m.group(2).equals("SR")){
	                                Element stereoChemEl = new Element(STEREOCHEMISTRY_EL);
	                                stereoChemEl.addAttribute(new Attribute(LOCANT_ATR, m.group(1)));
	                                stereoChemEl.addAttribute(new Attribute(VALUE_ATR, m.group(2).toUpperCase()));
	                                stereoChemEl.appendChild(stereoChemistryDescriptor);
	                                XOMTools.insertAfter(stereoChemistryElement, stereoChemEl);
	                                if (matchRS.matcher(m.group(2)).matches()) {
	                                    stereoChemEl.addAttribute(new Attribute(TYPE_ATR, R_OR_S_TYPE_VAL));
	                                } else {
	                                    stereoChemEl.addAttribute(new Attribute(TYPE_ATR, E_OR_Z_TYPE_VAL));
	                                }
                            	}
                            } else {
                                throw new ComponentGenerationException("Malformed stereochemistry element: " + stereoChemistryElement.getValue());
                            }
                        } else {
                            Element stereoChemEl = new Element(STEREOCHEMISTRY_EL);
                            stereoChemEl.addAttribute(new Attribute(VALUE_ATR, stereoChemistryDescriptor.toUpperCase()));
                            stereoChemEl.appendChild(stereoChemistryDescriptor);
                            XOMTools.insertAfter(stereoChemistryElement, stereoChemEl);
                            if (matchRS.matcher(stereoChemistryDescriptor).matches()) {
                                stereoChemEl.addAttribute(new Attribute(TYPE_ATR, R_OR_S_TYPE_VAL));
                            } else if (matchEZ.matcher(stereoChemistryDescriptor).matches()) {
                                stereoChemEl.addAttribute(new Attribute(TYPE_ATR, E_OR_Z_TYPE_VAL));
                            } else {
                                throw new ComponentGenerationException("Malformed stereochemistry element: " + stereoChemistryElement.getValue());
                            }
                        }
                    }
				}
				stereoChemistryElement.detach();
			}
			else if (stereoChemistryElement.getAttributeValue(TYPE_ATR).equals(CISORTRANS_TYPE_VAL)){//assign a locant if one is directly before the cis/trans
				Element possibleLocant = (Element) XOMTools.getPrevious(stereoChemistryElement);
				if (possibleLocant !=null && possibleLocant.getLocalName().equals(LOCANT_EL) && matchComma.split(possibleLocant.getValue()).length==1){
					stereoChemistryElement.addAttribute(new Attribute(LOCANT_ATR, possibleLocant.getValue()));
					possibleLocant.detach();
				}
			}
		}
	}

	/**
	 * Looks for "suffixPrefix" and assigns their value them as an attribute of an adjacent suffix
	 * @param subOrRoot
	 * @throws ComponentGenerationException
	 */
	private void processSuffixPrefixes(Element subOrRoot) throws ComponentGenerationException {
		List<Element> suffixPrefixes =  XOMTools.getChildElementsWithTagName(subOrRoot, SUFFIXPREFIX_EL);
		for (Element suffixPrefix : suffixPrefixes) {
			Element suffix = (Element) XOMTools.getNextSibling(suffixPrefix);
			if (suffix==null || ! suffix.getLocalName().equals(SUFFIX_EL)){
				throw new ComponentGenerationException("OPSIN bug: suffix not found after suffixPrefix: " + suffixPrefix.getValue());
			}
			suffix.addAttribute(new Attribute(SUFFIXPREFIX_ATR, suffixPrefix.getAttributeValue(VALUE_ATR)));
			suffixPrefix.detach();
		}
	}

	/**
	 * Looks for infixes and assigns them to the next suffix using a semicolon delimited infix attribute
	 * If the infix/suffix block has been bracketed e.g (dithioate) then the infix is multiplied out
	 * If preceded by a suffixPrefix e.g. sulfono infixes are also multiplied out
	 * If a multiplier is present and neither of these cases are met then it is ambiguous as to whether the multiplier is referring to the infix or the infixed suffix
	 * This ambiguity is resolved in processInfixFunctionalReplacementNomenclature by looking at the structure of the suffix to be modified
	 * @param subOrRoot
	 * @throws ComponentGenerationException
	 */
	private void processInfixes(Element subOrRoot) throws ComponentGenerationException {
		List<Element> infixes = XOMTools.getChildElementsWithTagName(subOrRoot, INFIX_EL);
		for (Element infix : infixes) {
			Element suffix = XOMTools.getNextSiblingIgnoringCertainElements(infix, new String[]{INFIX_EL, SUFFIXPREFIX_EL});
			if (suffix ==null || !suffix.getLocalName().equals(SUFFIX_EL)){
				throw new ComponentGenerationException("No suffix found next next to infix: "+ infix.getValue());
			}
			List<String> currentInfixInformation;
			if (suffix.getAttribute(INFIX_ATR)==null){
				suffix.addAttribute(new Attribute(INFIX_ATR, ""));
				currentInfixInformation = new ArrayList<String>();
			}
			else{
				currentInfixInformation = StringTools.arrayToList(matchSemiColon.split(suffix.getAttributeValue(INFIX_ATR)));
			}
			String infixValue =infix.getAttributeValue(VALUE_ATR);
			currentInfixInformation.add(infixValue);
			Element possibleMultiplier = (Element) XOMTools.getPreviousSibling(infix);
			Element possibleBracket;
			boolean suffixPrefixPresent =false;
			if (possibleMultiplier.getLocalName().equals(MULTIPLIER_EL)){
				Element possibleSuffixPrefix = XOMTools.getPreviousSiblingIgnoringCertainElements(infix, new String[]{MULTIPLIER_EL, INFIX_EL});
				if (possibleSuffixPrefix!=null && possibleSuffixPrefix.getLocalName().equals(SUFFIXPREFIX_EL)){
					suffixPrefixPresent =true;
				}
				possibleBracket  = (Element) XOMTools.getPreviousSibling(possibleMultiplier);
			}
			else{
				possibleBracket=possibleMultiplier;
				possibleMultiplier=null;
				infix.detach();
			}
			if (possibleBracket.getLocalName().equals(STRUCTURALOPENBRACKET_EL)){
				Element bracket = (Element) XOMTools.getNextSibling(suffix);
				if (!bracket.getLocalName().equals(STRUCTURALCLOSEBRACKET_EL)){
					throw new ComponentGenerationException("Matching closing bracket not found around infix/suffix block");
				}
				if (possibleMultiplier!=null){
					int multiplierVal = Integer.parseInt(possibleMultiplier.getAttributeValue(VALUE_ATR));
					for (int i = 1; i < multiplierVal; i++) {
						currentInfixInformation.add(infixValue);
					}
					possibleMultiplier.detach();
					infix.detach();
				}
				possibleBracket.detach();
				bracket.detach();
			}
			else if (possibleMultiplier!=null && suffixPrefixPresent){
				int multiplierVal = Integer.parseInt(possibleMultiplier.getAttributeValue(VALUE_ATR));
				for (int i = 1; i < multiplierVal; i++) {
					currentInfixInformation.add(infixValue);
				}
				possibleMultiplier.detach();
				infix.detach();
			}
			suffix.getAttribute(INFIX_ATR).setValue(StringTools.stringListToString(currentInfixInformation, ";"));
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
	 * @throws ComponentGenerationException
	 */
	private void processLambdaConvention(Element subOrRoot) throws ComponentGenerationException {
		List<Element> lambdaConventionEls = XOMTools.getChildElementsWithTagName(subOrRoot, LAMBDACONVENTION_EL);
		boolean fusedRingPresent = false;
		if (lambdaConventionEls.size()>0){
			if (subOrRoot.getChildElements(GROUP_EL).size()>1){
				fusedRingPresent = true;
			}
		}
		for (Element lambdaConventionEl : lambdaConventionEls) {
			boolean frontLocantsExpected =false;//Is the lambdaConvention el followed by benz/benzo of a fused ring system (these have front locants which correspond to the final fused rings numbering) or by a polycylicspiro system
			String[] lambdaValues = matchComma.split(StringTools.removeDashIfPresent(lambdaConventionEl.getValue()));
			Element possibleHeteroatomOrMultiplier = (Element) XOMTools.getNextSibling(lambdaConventionEl);
			int heteroCount = 0;
			int multiplierValue = 1;
			while(possibleHeteroatomOrMultiplier != null){
				if(possibleHeteroatomOrMultiplier.getLocalName().equals(HETEROATOM_EL)) {
					heteroCount+=multiplierValue;
					multiplierValue =1;
				} else if (possibleHeteroatomOrMultiplier.getLocalName().equals(MULTIPLIER_EL)){
					multiplierValue = Integer.parseInt(possibleHeteroatomOrMultiplier.getAttributeValue(VALUE_ATR));
				}
				else{
					break;
				}
				possibleHeteroatomOrMultiplier = (Element)XOMTools.getNextSibling(possibleHeteroatomOrMultiplier);
			}
			boolean assignLambdasToHeteroAtoms =false;
			if (lambdaValues.length==heteroCount){//heteroatom and number of locants +lambdas must match
				if (fusedRingPresent && possibleHeteroatomOrMultiplier!=null && possibleHeteroatomOrMultiplier.getLocalName().equals(GROUP_EL) && possibleHeteroatomOrMultiplier.getAttributeValue(SUBTYPE_ATR).equals(HANTZSCHWIDMAN_SUBTYPE_VAL)){
					//You must not set the locants of a HW system which forms a component of a fused ring system. The locant specified corresponds to the complete fused ring system.
				}
				else{
					assignLambdasToHeteroAtoms =true;
				}
			}
			else if(heteroCount==0 && XOMTools.getNextSibling(lambdaConventionEl).equals(possibleHeteroatomOrMultiplier) &&
					 possibleHeteroatomOrMultiplier!=null &&
						(fusedRingPresent && possibleHeteroatomOrMultiplier.getLocalName().equals(GROUP_EL) &&
						(possibleHeteroatomOrMultiplier.getValue().equals("benzo") || possibleHeteroatomOrMultiplier.getValue().equals("benz"))
						&& !((Element)XOMTools.getNextSibling(possibleHeteroatomOrMultiplier)).getLocalName().equals(FUSION_EL)) ||
						(possibleHeteroatomOrMultiplier.getLocalName().equals(POLYCYCLICSPIRO_EL) && 
								(possibleHeteroatomOrMultiplier.getAttributeValue(VALUE_ATR).equals("spirobi")|| possibleHeteroatomOrMultiplier.getAttributeValue(VALUE_ATR).equals("spiroter")))){
				frontLocantsExpected = true;//a benzo fused ring e.g. 1lambda4,3-benzothiazole or a symmetrical poly cyclic spiro system
			}
			List<Element> heteroAtoms = new ArrayList<Element>();//contains the heteroatoms to apply the lambda values too. Can be empty if the values are applied to a group directly rather than to a heteroatom
			if (assignLambdasToHeteroAtoms){//populate heteroAtoms, multiplied heteroatoms are multiplied out
				Element multiplier = null;
				Element heteroatomOrMultiplier = (Element) XOMTools.getNextSibling(lambdaConventionEl);
				while(heteroatomOrMultiplier != null){
					if(heteroatomOrMultiplier.getLocalName().equals(HETEROATOM_EL)) {
						heteroAtoms.add(heteroatomOrMultiplier);
						if (multiplier!=null){
							for (int i = 1; i < Integer.parseInt(multiplier.getAttributeValue(VALUE_ATR)); i++) {
								Element newHeteroAtom = new Element(heteroatomOrMultiplier);
								XOMTools.insertBefore(heteroatomOrMultiplier, newHeteroAtom);
								heteroAtoms.add(newHeteroAtom);
							}
							multiplier.detach();
							multiplier=null;
						}
					} else if (heteroatomOrMultiplier.getLocalName().equals(MULTIPLIER_EL)){
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
					Attribute valencyChange = new Attribute(LAMBDA_ATR, m.group(2));
					Attribute locantAtr = null;
					if (m.group(1)!=null){
						locantAtr = new Attribute(LOCANT_ATR, m.group(1));
					}
					if (frontLocantsExpected){
						if (m.group(1)==null){
							throw new ComponentGenerationException("Locant not found for lambda convention before a benzo fused ring system");
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
						Element newLambda = new Element(LAMBDACONVENTION_EL);
						newLambda.addAttribute(valencyChange);
						if (locantAtr!=null){
							newLambda.addAttribute(locantAtr);
						}
						XOMTools.insertBefore(lambdaConventionEl, newLambda);
					}
				}
				else{//just a locant e.g 1,3lambda5
					if (!assignLambdasToHeteroAtoms){
						if (!frontLocantsExpected){
							throw new ComponentGenerationException("Lambda convention not specified for locant: " + lambdaValue);
						}
					}
					else{
						Element heteroAtom = heteroAtoms.get(i);
						heteroAtom.addAttribute(new Attribute(LOCANT_ATR, lambdaValue));
					}
				}
			}
			if (!frontLocantsExpected){
				lambdaConventionEl.detach();
			}
			else{
				lambdaConventionEl.setLocalName(LOCANT_EL);
				XOMTools.setTextChild(lambdaConventionEl, StringTools.arrayToString(lambdaValues, ","));
			}
		}
	}

	/**Finds matching open and close brackets, and places the
	 * elements contained within in a big &lt;bracket&gt; element.
	 *
	 * @param substituentsAndRoot: The substituent/root elements at the current level of the tree
	 * @return Whether the method did something, and so needs to be called again.
	 * @throws ComponentGenerationException
	 */
	private boolean findAndStructureBrackets(List<Element> substituentsAndRoot) throws ComponentGenerationException {
		int blevel = 0;
		Element openBracket = null;
		Element closeBracket = null;
		for (Element sub : substituentsAndRoot) {
			Elements children = sub.getChildElements();
			for(int i=0; i<children.size(); i++) {
				Element child = children.get(i);
				if(child.getLocalName().equals(OPENBRACKET_EL)) {
					if(openBracket == null) {
						openBracket = child;
					}
					blevel++;
				} else if (child.getLocalName().equals(CLOSEBRACKET_EL)) {
					blevel--;
					if(blevel == 0) {
						closeBracket = child;
						Element bracket = structureBrackets(openBracket, closeBracket);
						while(findAndStructureBrackets(XOMTools.getDescendantElementsWithTagName(bracket, SUBSTITUENT_EL)));
						return true;
					}
				}
			}
		}
		if (blevel != 0){
			throw new ComponentGenerationException("Brackets do not match!");
		}
		return false;
	}

	/**Places the elements in substituents containing/between an open and close bracket
	 * in a &lt;bracket&gt; tag.
	 *
	 * @param openBracket The open bracket element
	 * @param closeBracket The close bracket element
	 * @return The bracket element thus created.
	 * @throws ComponentGenerationException 
	 */
	private Element structureBrackets(Element openBracket, Element closeBracket) throws ComponentGenerationException {
		Element bracket = new Element(BRACKET_EL);
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
			if (currentNode==null){
				throw new ComponentGenerationException("Brackets within a word do not match!");
			}
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

	/**Looks for annulen/polyacene/polyaphene/polyalene/polyphenylene/polynaphthylene/polyhelicene tags and replaces them with a group with appropriate SMILES.
	 * @param subOrRoot The subOrRoot to look for tags in
	 * @throws ComponentGenerationException
	 */
	private void processHydroCarbonRings(Element subOrRoot) throws ComponentGenerationException {
		List<Element> annulens = XOMTools.getChildElementsWithTagName(subOrRoot, ANNULEN_EL);
		for (Element annulen : annulens) {
			String annulenValue =annulen.getValue();
	        Matcher match = matchAnnulene.matcher(annulenValue);
	        match.matches();
	        if (match.groupCount() !=1){
	        	throw new ComponentGenerationException("Invalid annulen tag");
	        }

	        int annulenSize=Integer.valueOf(match.group(1));
	        if (annulenSize <3){
	        	throw new ComponentGenerationException("Invalid annulen tag");
	        }

			//build [annulenSize]annulene ring as SMILES
			String SMILES = "c1" +StringTools.multiplyString("c", annulenSize -1);
			SMILES += "1";

			Element group =new Element(GROUP_EL);
			group.addAttribute(new Attribute(VALUE_ATR, SMILES));
			group.addAttribute(new Attribute(VALTYPE_ATR, SMILES_VALTYPE_VAL));
			group.addAttribute(new Attribute(TYPE_ATR, RING_TYPE_VAL));
			group.addAttribute(new Attribute(SUBTYPE_ATR, ARYLGROUP_SUBTYPE_VAL));
			group.appendChild(annulenValue);
			annulen.getParent().replaceChild(annulen, group);
		}

		List<Element> hydrocarbonFRSystems = XOMTools.getChildElementsWithTagName(subOrRoot, HYDROCARBONFUSEDRINGSYSTEM_EL);
		for (Element hydrocarbonFRSystem : hydrocarbonFRSystems) {
			Element multiplier = (Element)XOMTools.getPreviousSibling(hydrocarbonFRSystem);
			if(multiplier != null && multiplier.getLocalName().equals(MULTIPLIER_EL)) {
				int multiplierValue =Integer.parseInt(multiplier.getAttributeValue(VALUE_ATR));
				String classOfHydrocarbonFRSystem =hydrocarbonFRSystem.getAttributeValue(VALUE_ATR);
				Element newGroup =new Element(GROUP_EL);
				String SMILES="";
				if (classOfHydrocarbonFRSystem.equals("polyacene")){
					if (multiplierValue <=3){
						throw new ComponentGenerationException("Invalid polyacene");
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
						throw new ComponentGenerationException("Invalid polyaphene");
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
						throw new ComponentGenerationException("Invalid polyalene");
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
						throw new ComponentGenerationException("Invalid polyphenylene");
					}
					SMILES= "c1cccc2";
					for (int j = 1; j < multiplierValue; j++) {
						SMILES+="c3ccccc3";
					}
					SMILES+= "c12";
				} else if (classOfHydrocarbonFRSystem.equals("polynaphthylene")){
					if (multiplierValue <3){
						throw new ComponentGenerationException("Invalid polynaphthylene");
					}
					SMILES= "c1cccc2cc3";
					for (int j = 1; j < multiplierValue; j++) {
						SMILES+="c4cc5ccccc5cc4";
					}
					SMILES+= "c3cc12";
				} else if (classOfHydrocarbonFRSystem.equals("polyhelicene")){
					if (multiplierValue <6){
						throw new ComponentGenerationException("Invalid polyhelicene");
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
					throw new ComponentGenerationException("Unknown semi-trivially named hydrocarbon fused ring system");
				}

				newGroup.addAttribute(new Attribute(VALUE_ATR, SMILES));
				newGroup.addAttribute(new Attribute(VALTYPE_ATR, SMILES_VALTYPE_VAL));
				newGroup.addAttribute(new Attribute(LABELS_ATR, FUSEDRING_LABELS_VAL));
				newGroup.addAttribute(new Attribute(TYPE_ATR, RING_TYPE_VAL));
				newGroup.addAttribute(new Attribute(SUBTYPE_ATR, HYDROCARBONFUSEDRINGSYSTEM_EL));
				newGroup.appendChild(multiplier.getValue() + hydrocarbonFRSystem.getValue());
				hydrocarbonFRSystem.getParent().replaceChild(hydrocarbonFRSystem, newGroup);
				multiplier.detach();
			}
			else{
				throw new ComponentGenerationException("Invalid semi-trivially named hydrocarbon fused ring system");
			}
		}
	}
	
	/**
	 * Handles irregular suffixes. Quinone only currently
	 * @param suffixes
	 * @throws ComponentGenerationException 
	 */
	private void handleSuffixIrregularities(Element subOrRoot) throws ComponentGenerationException {
		List<Element> suffixes = XOMTools.getChildElementsWithTagName(subOrRoot, SUFFIX_EL);
		for (Element suffix : suffixes) {
			String suffixValue = suffix.getValue();
			if (suffixValue.equals("ic") || suffixValue.equals("ous")){
				Node next = XOMTools.getNext(suffix);
				if (next == null){
					throw new ComponentGenerationException("\"acid\" not found after " +suffixValue);
				}
			}
			// convert quinone to dione
			else if (suffixValue.equals("quinone")){
				suffix.removeAttribute(suffix.getAttribute(ADDITIONALVALUE_ATR));
				XOMTools.setTextChild(suffix, "one");
				Element multiplier = new Element(MULTIPLIER_EL);
				multiplier.addAttribute(new Attribute(VALUE_ATR, "2"));
				multiplier.appendChild("di");
				XOMTools.insertBefore(suffix, multiplier);
			}
		}
	}

	/**
	 * Looks for alkaneStems followed by a bridge forming 'o' and makes them fused ring bridge elements
	 * @param group
	 */
	private void detectAlkaneFusedRingBridges(Element group) {
		if (ALKANESTEM_SUBTYPE_VAL.equals(group.getAttributeValue(SUBTYPE_ATR))){
			Element possibleBridgeFormer = (Element) XOMTools.getNextSiblingIgnoringCertainElements(group, new String[]{UNSATURATOR_EL});
			if(possibleBridgeFormer!=null && possibleBridgeFormer.getLocalName().equals(BRIDGEFORMINGO_EL)){
				possibleBridgeFormer.detach();
				group.setLocalName(FUSEDRINGBRIDGE_EL);
			}
		}
	}

	/**Looks (multiplier)cyclo/spiro/cyclo tags before chain
	 * and replaces them with a group with appropriate SMILES
	 * Note that only simple spiro tags are handled at this stage i.e. not dispiro
	 * @param group A group which is potentially a chain
	 * @throws ComponentGenerationException
	 */
	private void processRings(Element group) throws ComponentGenerationException {
		Element previous = (Element)XOMTools.getPreviousSibling(group);
		if(previous != null) {
			String previousElType = previous.getLocalName();
			if(previousElType.equals(SPIRO_EL)){
				processSpiroSystem(group, previous);
			} else if(previousElType.equals(VONBAEYER_EL)) {
				processVonBaeyerSystem(group, previous);
			}
			else if(previousElType.equals(CYCLO_EL)) {
				processCyclisedChain(group, previous);
			}
		}
	}

	/**
	 * Processes a spiro descriptor element.
	 * This modifies the provided chainGroup into the spiro system by replacing the value of the chain group with appropriate SMILES
	 * @param chainGroup
	 * @param spiroEl
	 * @throws ComponentGenerationException 
	 * @throws NumberFormatException 
	 */
	private void processSpiroSystem(Element chainGroup, Element spiroEl) throws NumberFormatException, ComponentGenerationException {
		int[][] spiroDescriptors = getSpiroDescriptors(StringTools.removeDashIfPresent(spiroEl.getValue()));

		Element multiplier =(Element)XOMTools.getPreviousSibling(spiroEl);
		int numberOfSpiros = 1;
		if (multiplier != null && multiplier.getLocalName().equals(MULTIPLIER_EL)){
			numberOfSpiros = Integer.parseInt(multiplier.getAttributeValue(VALUE_ATR));
			multiplier.detach();
		}
		int numberOfCarbonInDescriptors =0;
		for (int[] spiroDescriptor : spiroDescriptors) {
			numberOfCarbonInDescriptors += spiroDescriptor[0];
		}
		numberOfCarbonInDescriptors += numberOfSpiros;
		if (numberOfCarbonInDescriptors != chainGroup.getAttributeValue(VALUE_ATR).length()){
			throw new ComponentGenerationException("Disagreement between number of atoms in spiro descriptor: " + numberOfCarbonInDescriptors +" and number of atoms in chain: " + chainGroup.getAttributeValue(VALUE_ATR).length());
		}

		int numOfOpenedBrackets = 1;
		int curIndex = 2;
		String smiles = "C0" + StringTools.multiplyString("C", spiroDescriptors[0][0]) + "10(";

		// for those molecules where no superstrings compare prefix number with curIndex.
		for (int i=1; i< spiroDescriptors.length; i++) {
			if (spiroDescriptors[i][1] >= 0) {
				int ringOpeningPos = findIndexOfRingOpenings( smiles, spiroDescriptors[i][1] );
				String ringOpeningLabel = String.valueOf(smiles.charAt(ringOpeningPos));
				ringOpeningPos++;
				if (ringOpeningLabel.equals("%")){
					while (smiles.charAt(ringOpeningPos)>='0' && smiles.charAt(ringOpeningPos)<='9' && ringOpeningPos < smiles.length()) {
						ringOpeningLabel += smiles.charAt(ringOpeningPos);
						ringOpeningPos++;
					}
				}
				if (smiles.indexOf("C" + ringOpeningLabel, ringOpeningPos)>=0) {
					// this ring opening has already been closed
					// i.e. this atom connects more than one ring in a spiro fusion
					
					// insert extra ring opening
					smiles = smiles.substring(0, ringOpeningPos) + ringClosure(curIndex) + smiles.substring(ringOpeningPos);

					// add ring in new brackets
					smiles += "(" + StringTools.multiplyString("C", spiroDescriptors[i][0]) + ringClosure(curIndex) + ")";
					curIndex++;
				}
				else {
					smiles += StringTools.multiplyString("C", spiroDescriptors[i][0]) + ringOpeningLabel + ")";
				}
			}
			else if (numOfOpenedBrackets >= numberOfSpiros) {
				smiles += StringTools.multiplyString("C", spiroDescriptors[i][0]);

				// take the number before bracket as index for smiles
				// we can open more brackets, this considered in prev if
				curIndex--;
				smiles += ringClosure(curIndex) + ")";

				// from here start to decrease index for the following
			}
			else {
				smiles += StringTools.multiplyString("C", spiroDescriptors[i][0]);
				smiles += "C" + ringClosure(curIndex++) + "(";
				numOfOpenedBrackets++;
			}
		}
		chainGroup.getAttribute(VALUE_ATR).setValue(smiles);
		chainGroup.getAttribute(TYPE_ATR).setValue(RING_TYPE_VAL);
		if (chainGroup.getAttribute(USABLEASJOINER_ATR) !=null){
			chainGroup.removeAttribute(chainGroup.getAttribute(USABLEASJOINER_ATR));
		}
		spiroEl.detach();
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
	 * @param text - string with spiro e.g. spiro[2.2]
	 * @return array with number of carbons in each group and associated index of spiro atom
	 */
	private int[][] getSpiroDescriptors(String text) {
		if (text.indexOf("-")==5){
			text= text.substring(7, text.length()-1);//cut off spiro-[ and terminal ]
		}
		else{
			text= text.substring(6, text.length()-1);//cut off spiro[ and terminal ]
		}
		
		String[] spiroDescriptorStrings = matchCommaOrDot.split(text);
	
		int[][] spiroDescriptors = new int[spiroDescriptorStrings.length][2]; // array of descriptors where number of elements and super string present
	
		for (int i=0; i < spiroDescriptorStrings.length; i++) {
			String[] elements = matchNonDigit.split(spiroDescriptorStrings[i]);
			if (elements.length >1) {//a "superscripted" number is present
				spiroDescriptors[i][0] = Integer.parseInt(elements[0]);
				String superScriptedNumber ="";
				for (int j = 1; j < elements.length; j++){//may be more than one non digit as there are many ways of indicating superscripts
					superScriptedNumber += elements[j];
				}
				spiroDescriptors[i][1] = Integer.parseInt(superScriptedNumber);
			}
			else {
				spiroDescriptors[i][0] = Integer.parseInt(spiroDescriptorStrings[i]);
				spiroDescriptors[i][1] = -1;
			}
		}
	
		return spiroDescriptors;
	}

	/**
	 * Finds the the carbon atom with the given locant in the provided SMILES
	 * Returns the next index which is expected to correspond to the atom's ring opening/s
	 * @param smiles string to search in
	 * @param locant locant of the atom in given structure
	 * @return index of ring openings
	 * @throws ComponentGenerationException 
	 */
	private Integer findIndexOfRingOpenings(String smiles, int locant) throws ComponentGenerationException{
		int count = 0;
		int pos = -1;
		for (int i=0; i<smiles.length(); i++){
			if (smiles.charAt(i) == 'C') {
				count++;
			}
			if (count==locant) {
				pos=i;
				break;
			}
		}
		if (pos == -1){
			throw new ComponentGenerationException("Unable to find atom corresponding to number indicated by superscript in spiro descriptor");
		}
		pos++;
		return pos;
	}

	/**
	 * Given an element corresponding to an alkane or other systematic chain and the preceding vonBaeyerBracket element:
	 * Generates the SMILES of the von baeyer system and assigns this to the chain Element
	 * Checks are done on the von baeyer multiplier and chain length
	 * The multiplier and vonBaeyerBracket are detached
	 * @param chainEl
	 * @param vonBaeyerBracketEl
	 * @throws ComponentGenerationException
	 */
	private void processVonBaeyerSystem(Element chainEl, Element vonBaeyerBracketEl) throws ComponentGenerationException {
		String vonBaeyerBracket = StringTools.removeDashIfPresent(vonBaeyerBracketEl.getValue());
		Element multiplier =(Element)XOMTools.getPreviousSibling(vonBaeyerBracketEl);
		int numberOfRings=Integer.parseInt(multiplier.getAttributeValue(VALUE_ATR));
		multiplier.detach();

		int alkylChainLength;
		LinkedList<String> elementSymbolArray = new LinkedList<String>();
		if (chainEl.getAttributeValue(VALTYPE_ATR).equals(SMILES_VALTYPE_VAL)){
			String smiles =chainEl.getAttributeValue(VALUE_ATR);
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
			throw new ComponentGenerationException("unexpected group valType: " + chainEl.getAttributeValue(VALTYPE_ATR));
		}


		int totalLengthOfBridges=0;
		int bridgeLabelsUsed=3;//start labelling from 3 upwards
		//3 and 4 will be the atoms on each end of one secondary bridge, 5 and 6 for the next etc.

		ArrayList<HashMap<String, Integer>> bridges = new ArrayList<HashMap<String, Integer>>();
		HashMap<Integer, ArrayList<Integer>> bridgeLocations = new HashMap<Integer, ArrayList<Integer>>(alkylChainLength);
		if (vonBaeyerBracket.indexOf("-")==5){
			vonBaeyerBracket = vonBaeyerBracket.substring(7, vonBaeyerBracket.length()-1);//cut off cyclo-[ and terminal ]
		}
		else{
			vonBaeyerBracket = vonBaeyerBracket.substring(6, vonBaeyerBracket.length()-1);//cut off cyclo[ and terminal ]
		}
		String[] bridgeDescriptors = matchCommaOrDot.split(vonBaeyerBracket);//the bridgelengths and positions for secondary bridges
		//all bridges from past the first 3 are secondary bridges and require specification of bridge position which will be partially in the subsequent position in the array
		for (int i = 0; i < bridgeDescriptors.length; i++) {
			String bridgeDescriptor = bridgeDescriptors[i];
			HashMap<String, Integer> bridge = new HashMap<String, Integer>();
			int bridgeLength =0;
			if (i > 2){//this is a secondary bridge (chain start/end locations should have been specified)
				i++;
				String coordinatesStr1;
				String coordinatesStr2 = bridgeDescriptors[i].replaceAll("\\D", "");
				String[] tempArray = bridgeDescriptor.split("\\D+");

				if (tempArray.length ==1){
					//there is some ambiguity as it has not been made obvious which number/s are supposed to be the superscripted locant
					//so we assume that it is more likely that it will be referring to an atom of label >10
					//rather than a secondary bridge of length > 10
					char[] tempCharArray = bridgeDescriptor.toCharArray();
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
						throw new ComponentGenerationException("Unsupported Von Baeyer locant description: " + bridgeDescriptor );
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
					throw new ComponentGenerationException("Indicated bridge position is not on chain: " +coordinates1 +"," +coordinates2);
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
				bridgeLength= Integer.parseInt(bridgeDescriptor);
				bridge.put("Bridge Length", bridgeLength);
			}
			totalLengthOfBridges += bridgeLength;
			bridges.add(bridge);
		}
		if (totalLengthOfBridges + 2 !=alkylChainLength ){
			throw new ComponentGenerationException("Disagreement between lengths of bridges and alkyl chain length");
		}
		if (numberOfRings +1 != bridges.size()){
			throw new ComponentGenerationException("Disagreement between number of rings and number of bridges");
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
				throw new ComponentGenerationException("Unable to resolve all dependant bridges!!!");
			}
			secondaryBridges=dependantSecondaryBridges;
		}
		while(dependantSecondaryBridges.size() > 0);

		chainEl.getAttribute(VALUE_ATR).setValue(SMILES);
		chainEl.getAttribute(TYPE_ATR).setValue(RING_TYPE_VAL);
		if (chainEl.getAttribute(USABLEASJOINER_ATR) !=null){
			chainEl.removeAttribute(chainEl.getAttribute(USABLEASJOINER_ATR));
		}
		vonBaeyerBracketEl.detach();
	}

	/**
	 * Converts a chain group into a ring.
	 * The chain group can either be an alkane or heteroatom chain
	 * @param chainGroup
	 * @param cycloEl
	 * @throws ComponentGenerationException
	 */
	private void processCyclisedChain(Element chainGroup, Element cycloEl) throws ComponentGenerationException {
		String smiles=chainGroup.getAttributeValue(VALUE_ATR);
		int chainlen =0;
		for (int i = smiles.length() -1 ; i >=0; i--) {
			if (Character.isUpperCase(smiles.charAt(i)) && smiles.charAt(i) !='H'){
				chainlen++;
			}
	    }
		if (chainlen < 3){
			throw new ComponentGenerationException("Heteroatom chain too small to create a ring: " + chainlen);
		}
		smiles+="1";
		if (smiles.charAt(0)=='['){
			int closeBracketIndex = smiles.indexOf(']');
			smiles= smiles.substring(0, closeBracketIndex +1) +"1" + smiles.substring(closeBracketIndex +1);
		}
		else{
			if (Character.getType(smiles.charAt(1)) == Character.LOWERCASE_LETTER){//element is 2 letters long
				smiles= smiles.substring(0,2) +"1" + smiles.substring(2);
			}
			else{
				smiles= smiles.substring(0,1) +"1" + smiles.substring(1);
			}
		}
		chainGroup.getAttribute(VALUE_ATR).setValue(smiles);
		if (chainlen==6){//6 membered rings have ortho/meta/para positions
			if (chainGroup.getAttribute(LABELS_ATR)!=null){
				chainGroup.getAttribute(LABELS_ATR).setValue("1/2,ortho/3,meta/4,para/5/6");
			}
			else{
				chainGroup.addAttribute(new Attribute(LABELS_ATR, "1/2,ortho/3,meta/4,para/5/6"));
			}
		}
		chainGroup.getAttribute(TYPE_ATR).setValue(RING_TYPE_VAL);
		if (chainGroup.getAttribute(USABLEASJOINER_ATR) !=null){
			chainGroup.removeAttribute(chainGroup.getAttribute(USABLEASJOINER_ATR));
		}
		cycloEl.detach();
	}

	/**Handles special cases in IUPAC nomenclature.
	 * Benzyl etc.
	 * @param group The group to look for irregularities in.
	 */
	private void handleGroupIrregularities(Element group) throws ComponentGenerationException {
		String groupValue =group.getValue();
		if(groupValue.equals("thiophen")) {//thiophenol is phenol with an O replaced with S not thiophene with a hydroxy
			Element possibleSuffix = (Element) XOMTools.getNextSibling(group);
			if (!"e".equals(group.getAttributeValue(SUBSEQUENTUNSEMANTICTOKEN_EL)) && possibleSuffix !=null && possibleSuffix.getLocalName().equals(SUFFIX_EL)) {
				if (possibleSuffix.getValue().startsWith("ol")){
					throw new ComponentGenerationException("thiophenol has been incorrectly interpreted as thiophen, ol instead of thio, phenol");
				}
			}
		}
		else if (groupValue.equals("methylene")) {//e.g. 3,4-methylenedioxyphenyl
			Element nextSub = (Element) XOMTools.getNextSibling(group.getParent());
			if (nextSub !=null && nextSub.getLocalName().equals(SUBSTITUENT_EL) && XOMTools.getNextSibling(group)==null 
					&& (XOMTools.getPreviousSibling(group)==null || !((Element)XOMTools.getPreviousSibling(group)).getLocalName().equals(MULTIPLIER_EL))){//not trimethylenedioxy
				Elements children = nextSub.getChildElements();
				if (children.size() >=2 && children.get(0).getValue().equals("di")&& children.get(1).getValue().equals("oxy")){
					XOMTools.setTextChild(group, "methylenedioxy");
					group.getAttribute(VALUE_ATR).setValue("C(O)O");
					group.getAttribute(VALTYPE_ATR).setValue(SMILES_VALTYPE_VAL);
					group.getAttribute(OUTIDS_ATR).setValue("2,3");
					group.getAttribute(SUBTYPE_ATR).setValue(EPOXYLIKE_SUBTYPE_VAL);
					if (group.getAttribute(LABELS_ATR)!=null){
						group.getAttribute(LABELS_ATR).setValue(NONE_LABELS_VAL);
					}
					else{
						group.addAttribute(new Attribute(LABELS_ATR, NONE_SUBTYPE_VAL));
					}
					nextSub.detach();
					for (int i = children.size() -1 ; i >=2; i--) {
						children.get(i).detach();
						XOMTools.insertAfter(group, children.get(i));
					}
				}
			}
		}
		else if (groupValue.equals("ethylene")) {
			Element previous = (Element)XOMTools.getPreviousSibling(group);
			if (previous!=null && previous.getLocalName().equals(MULTIPLIER_EL)){
				int multiplierValue = Integer.parseInt(previous.getAttributeValue(VALUE_ATR));
				Element possibleRoot =(Element) XOMTools.getNextSibling(group.getParent());
				if (possibleRoot==null && OpsinTools.getParentWordRule(group).getAttributeValue(WORDRULE_ATR).equals(WordRule.glycol.toString())){//e.g. dodecaethylene glycol
					String smiles ="CC";
					for (int i = 1; i < multiplierValue; i++) {
						smiles+="OCC";
					}
					group.getAttribute(OUTIDS_ATR).setValue("1," +Integer.toString(3*(multiplierValue-1) +2));
					group.getAttribute(VALUE_ATR).setValue(smiles);
					previous.detach();
				}
				else if (possibleRoot!=null && possibleRoot.getLocalName().equals(ROOT_EL)){
					Elements children = possibleRoot.getChildElements();
					if (children.size()==2){
						Element amineMultiplier =children.get(0);
						Element amine =children.get(1);
						if (amineMultiplier.getLocalName().equals(MULTIPLIER_EL) && amine.getValue().equals("amin")){//e.g. Triethylenetetramine
							if (Integer.parseInt(amineMultiplier.getAttributeValue(VALUE_ATR))!=multiplierValue +1){
								throw new ComponentGenerationException("Invalid polyethylene amine!");
							}
							String smiles ="";
							for (int i = 0; i < multiplierValue; i++) {
								smiles+="NCC";
							}
							smiles+="N";
							group.removeAttribute(group.getAttribute(OUTIDS_ATR));
							group.getAttribute(VALUE_ATR).setValue(smiles);
							previous.detach();
							possibleRoot.detach();
							((Element)group.getParent()).setLocalName(ROOT_EL);
						}
					}
				}
			}
		}
		else if (groupValue.equals("propylene")) {
			Element previous = (Element)XOMTools.getPreviousSibling(group);
			if (previous!=null && previous.getLocalName().equals(MULTIPLIER_EL)){
				int multiplierValue = Integer.parseInt(previous.getAttributeValue(VALUE_ATR));
				Element possibleRoot =(Element) XOMTools.getNextSibling(group.getParent());
				if (possibleRoot==null && OpsinTools.getParentWordRule(group).getAttributeValue(WORDRULE_ATR).equals(WordRule.glycol.toString())){//e.g. dodecaethylene glycol
					String smiles ="CCC";
					for (int i = 1; i < multiplierValue; i++) {
						smiles+="OC(C)C";
					}
					group.getAttribute(OUTIDS_ATR).setValue("2," +Integer.toString(4*(multiplierValue-1) +3));
					group.getAttribute(VALUE_ATR).setValue(smiles);
					if (group.getAttribute(LABELS_ATR)!=null){
						group.getAttribute(LABELS_ATR).setValue(NONE_LABELS_VAL);
					}
					else{
						group.addAttribute(new Attribute(LABELS_ATR, NONE_LABELS_VAL));
					}
					previous.detach();
				}
			}
		}

		//anthrone, phenanthrone and xanthone have the one at position 9 by default
		else if (groupValue.equals("anthr") || groupValue.equals("phenanthr") || groupValue.equals("xanth")|| groupValue.equals("xanthen")) {
			Element possibleLocant = (Element) XOMTools.getPreviousSibling(group);
			if (possibleLocant==null || !possibleLocant.getLocalName().equals(LOCANT_EL)){//only need to give one a locant of 9 if no locant currently present
				Element possibleOne =(Element) XOMTools.getNextSibling(group);
				if (possibleOne!=null && possibleOne.getValue().equals("one")){
					//Rule C-315.2
					Element newLocant =new Element(LOCANT_EL);
					newLocant.appendChild("9");
					XOMTools.insertBefore(possibleOne, newLocant);
					Element newIindicatedHydrogen = new Element(INDICATEDHYDROGEN_EL);
					newIindicatedHydrogen.addAttribute(new Attribute(LOCANT_ATR, "10"));
					XOMTools.insertBefore(newLocant, newIindicatedHydrogen);
				}
			}
		}
		else if (groupValue.equals("phospho")){//is this the organic meaning (P(=O)=O) or biochemical meaning (P(=O)(O)O)
			Element substituent = (Element) group.getParent();
			Element nextSubstituent = (Element) XOMTools.getNextSibling(substituent);
			if (nextSubstituent !=null){
				Element nextGroup = nextSubstituent.getFirstChildElement(GROUP_EL);
				if (nextGroup !=null && (nextGroup.getAttributeValue(TYPE_ATR).equals(AMINOACID_TYPE_VAL)||BIOCHEMICAL_SUBTYPE_VAL.equals(nextGroup.getAttributeValue(SUBTYPE_ATR)))){
					group.getAttribute(VALUE_ATR).setValue("-P(=O)(O)O");
					group.addAttribute(new Attribute(USABLEASJOINER_ATR, "yes"));
				}
			}
			
		}
		else if (groupValue.equals("cyste")){//ambiguity between cysteine and cysteic acid
			if (group.getAttributeValue(SUBTYPE_ATR).equals(ENDININE_SUBTYPE_VAL)){//cysteine
				Element ine = (Element) XOMTools.getNextSibling(group);
				if (!ine.getAttributeValue(VALUE_ATR).equals("ine")){
					throw new ComponentGenerationException("This is a cysteic acid derivative, not a cysteine derivative");
				}
			}
		}
		else if (groupValue.equals("hydrogen")){
			Element hydrogenParentEl = (Element) group.getParent();
			Element nextSubOrRoot = (Element) XOMTools.getNextSibling(hydrogenParentEl);
			if (nextSubOrRoot!=null){
				Element possibleSuitableAteGroup = (Element) nextSubOrRoot.getChild(0);
				if (!possibleSuitableAteGroup.getLocalName().equals(GROUP_EL) || !NONCARBOXYLICACID_TYPE_VAL.equals(possibleSuitableAteGroup.getAttributeValue(TYPE_ATR))){
					throw new ComponentGenerationException("Hydrogen is not meant as a substituent in this context!");
				}
				Element possibleMultiplier = (Element) XOMTools.getPreviousSibling(group);
				String multiplier = "1";
				if (possibleMultiplier!=null && possibleMultiplier.getLocalName().equals(MULTIPLIER_EL)){
					multiplier = possibleMultiplier.getAttributeValue(VALUE_ATR);
					possibleMultiplier.detach();
				}
				possibleSuitableAteGroup.addAttribute(new Attribute("numberOfFunctionalAtomsToRemove", multiplier));
				group.detach();
				Elements childrenToMove = hydrogenParentEl.getChildElements();
				for (int i = childrenToMove.size() -1 ; i >=0; i--) {
					childrenToMove.get(i).detach();
					nextSubOrRoot.insertChild(childrenToMove.get(i), 0);
				}
				hydrogenParentEl.detach();
			}
		}
		else if (groupValue.equals("acryl")){
			if (SIMPLESUBSTITUENT_SUBTYPE_VAL.equals(group.getAttributeValue(SUBTYPE_ATR))){
				Element nextEl = (Element) XOMTools.getNext(group);
				if (nextEl!=null && nextEl.getValue().equals("amid")){
					throw new ComponentGenerationException("amide in acrylamide is not [NH2-]");
				}
			}
		}
		else if (groupValue.equals("azo") || groupValue.equals("azoxy") || groupValue.equals("nno-azoxy") || groupValue.equals("non-azoxy") || groupValue.equals("onn-azoxy")){
			Element enclosingSub = (Element) group.getParent();
			Element next = (Element) XOMTools.getNextSiblingIgnoringCertainElements(enclosingSub, new String[]{HYPHEN_EL});
			if (next==null && XOMTools.getPreviousSibling(enclosingSub)==null){//e.g. [(E)-NNO-azoxy]benzene
				next = (Element) XOMTools.getNextSiblingIgnoringCertainElements((Element) enclosingSub.getParent(), new String[]{HYPHEN_EL});
			}
			if (next!=null && next.getLocalName().equals(ROOT_EL)){
				if (!(((Element)next.getChild(0)).getLocalName().equals(MULTIPLIER_EL))){
					List<Element> suffixes = XOMTools.getChildElementsWithTagName(next, SUFFIX_EL);
					if (suffixes.size()==0){//only case without locants is handled so far. suffixes only apply to one of the fragments rather than both!!!
						Element newMultiplier = new Element(MULTIPLIER_EL);
						newMultiplier.addAttribute(new Attribute(VALUE_ATR, "2"));
						next.insertChild(newMultiplier, 0);
						Element interSubstituentHyphen = (Element) XOMTools.getPrevious(group);
						if (interSubstituentHyphen!=null && !interSubstituentHyphen.getLocalName().equals(HYPHEN_EL)){//prevent implicit bracketting
							XOMTools.insertAfter(interSubstituentHyphen, new Element(HYPHEN_EL));
						}
					}
				}
			}
		}
		else if (groupValue.equals("coenzyme a") || groupValue.equals("coa")){
			Element enclosingSubOrRoot = (Element) group.getParent();
			Element previous = (Element) XOMTools.getPreviousSibling(enclosingSubOrRoot);
			if (previous!=null){
				List<Element> groups = XOMTools.getDescendantElementsWithTagName(previous, GROUP_EL);
				if (groups.size()>0){
					Element possibleDiAcid = groups.get(groups.size()-1);
					if (possibleDiAcid.getAttribute(SUFFIXAPPLIESTO_ATR)!=null && ACIDSTEM_TYPE_VAL.equals(possibleDiAcid.getAttributeValue(TYPE_ATR))){
						Element suffix = (Element) XOMTools.getNextSibling(possibleDiAcid, SUFFIX_EL);
						if (suffix.getAttribute(ADDITIONALVALUE_ATR)==null){
							suffix.addAttribute(new Attribute(ADDITIONALVALUE_ATR, "ic"));
						}
					}
				}
			}
			//locanted substitution onto Coenzyme A is rarely intended, so put it in a bracket to disfavour it
			Element newBracket = new Element(BRACKET_EL);
			XOMTools.insertAfter(enclosingSubOrRoot, newBracket);
			enclosingSubOrRoot.detach();
			newBracket.appendChild(enclosingSubOrRoot);
		}
	}
	
	/**
	 * Moves substituents into new words if it appears that a space was omitted in the input name
	 * @param elem
	 */
	private void addOmittedSpaces(Element elem)  {
		List<Element> wordRules = XOMTools.getDescendantElementsWithTagName(elem, WORDRULE_EL);
		for (Element wordRule : wordRules) {
			if (WordRule.valueOf(wordRule.getAttributeValue(WORDRULE_ATR)) == WordRule.divalentFunctionalGroup){
				List<Element> substituentWords = XOMTools.getChildElementsWithTagNameAndAttribute(wordRule, WORD_EL, TYPE_ATR, SUBSTITUENT_TYPE_VAL);
				if (substituentWords.size()==1){//potentially been "wrongly" interpreted e.g. ethylmethyl ketone is more likely to mean ethyl methyl ketone
					Elements children  =substituentWords.get(0).getChildElements();
					if (children.size()==2){
						Element firstChildOfFirstSubstituent =(Element)children.get(0).getChild(0);
						//rule out correct usage e.g. diethyl ether and locanted substituents e.g. 2-methylpropyl ether
						if (!firstChildOfFirstSubstituent.getLocalName().equals(LOCANT_EL) && !firstChildOfFirstSubstituent.getLocalName().equals(MULTIPLIER_EL)){
							Element subToMove =children.get(1);
							subToMove.detach();
							Element newWord =new Element(WORD_EL);
							newWord.addAttribute(new Attribute(TYPE_ATR, SUBSTITUENT_TYPE_VAL));
							newWord.appendChild(subToMove);
							XOMTools.insertAfter(substituentWords.get(0), newWord);
						}
					}
				}
			}
		}
	}

}

