package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Deque;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;
import static uk.ac.cam.ch.wwmm.opsin.OpsinTools.*;

/**Does destructive procedural parsing on parser results.
 *
 * @author ptc24
 * @author dl387
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
	private static class VonBaeyerSecondaryBridgeSort implements Comparator<HashMap<String, Integer>> {

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
	private final static Pattern matchNumberLocantsOnlyFusionBracket = Pattern.compile("\\[\\d+(,\\d+)*\\]");
	private final static Pattern matchCommaOrDot =Pattern.compile("[\\.,]");
	private final static Pattern matchAnnulene = Pattern.compile("[\\[\\(\\{]([1-9]\\d*)[\\]\\)\\}]annulen");
	private final static String elementSymbols ="(?:He|Li|Be|B|C|N|O|F|Ne|Na|Mg|Al|Si|P|S|Cl|Ar|K|Ca|Sc|Ti|V|Cr|Mn|Fe|Co|Ni|Cu|Zn|Ga|Ge|As|Se|Br|Kr|Rb|Sr|Y|Zr|Nb|Mo|Tc|Ru|Rh|Pd|Ag|Cd|In|Sn|Sb|Te|I|Xe|Cs|Ba|La|Ce|Pr|Nd|Pm|Sm|Eu|Gd|Tb|Dy|Ho|Er|Tm|Yb|Lu|Hf|Ta|W|Re|Os|Ir|Pt|Au|Hg|Tl|Pb|Po|At|Rn|Fr|Ra|Ac|Th|Pa|U|Np|Pu|Am|Cm|Bk|Cf|Es|Fm|Md|No|Lr|Rf|Db|Sg|Bh|Hs|Mt|Ds)";
	private final static Pattern matchStereochemistry = Pattern.compile("(.*?)(SR|RS|[RSEZrsezabx]|[cC][iI][sS]|[tT][rR][aA][nN][sS]|[aA][lL][pP][hH][aA]|[bB][eE][tT][aA]|[xX][iI]|[eE][xX][oO]|[eE][nN][dD][oO]|[sS][yY][nN]|[aA][nN][tT][iI])");
	private final static Pattern matchStar = Pattern.compile("\\^?\\*");
	private final static Pattern matchRS = Pattern.compile("[RSrs]");
	private final static Pattern matchEZ = Pattern.compile("[EZez]");
	private final static Pattern matchAlphaBetaStereochem = Pattern.compile("a|b|x|[aA][lL][pP][hH][aA]|[bB][eE][tT][aA]|[xX][iI]");
	private final static Pattern matchCisTrans = Pattern.compile("[cC][iI][sS]|[tT][rR][aA][nN][sS]");
	private final static Pattern matchEndoExoSynAnti = Pattern.compile("[eE][xX][oO]|[eE][nN][dD][oO]|[sS][yY][nN]|[aA][nN][tT][iI]");
	private final static Pattern matchLambdaConvention = Pattern.compile("(\\S+)?lambda\\D*(\\d+)\\D*", Pattern.CASE_INSENSITIVE);
	private final static Pattern matchHdigit =Pattern.compile("H\\d");
	private final static Pattern matchDigit =Pattern.compile("\\d+");
	private final static Pattern matchNonDigit =Pattern.compile("\\D+");
	private final static Pattern matchSuperscriptedLocant = Pattern.compile("(" + elementSymbols +"'*)[\\^\\[\\(\\{~]*(?:[sS][uU][pP][ ]?)?([^\\^\\[\\(\\{~\\]\\)\\}]+)[^\\[\\(\\{]*");
	private final static Pattern matchIUPAC2004ElementLocant = Pattern.compile("(\\d+'*)-(" + elementSymbols +"'*)(.*)");
	private final static Pattern matchBracketAtEndOfLocant = Pattern.compile("-?[\\[\\(\\{](.*)[\\]\\)\\}]$");
	private final static Pattern matchGreek = Pattern.compile("alpha|beta|gamma|delta|epsilon|zeta|eta|omega", Pattern.CASE_INSENSITIVE);
	private final static Pattern matchInlineSuffixesThatAreAlsoGroups = Pattern.compile("carbonyl|oxy|sulfenyl|sulfinyl|sulfonyl|selenenyl|seleninyl|selenonyl|tellurenyl|tellurinyl|telluronyl");

	
	private final NameToStructureConfig n2sConfig;
	
	ComponentGenerator(NameToStructureConfig n2sConfig) {
		this.n2sConfig = n2sConfig;
	}

	/**
	 * Processes a parse result destructively adding semantic information by processing the various micro syntaxes.
	 * @param parse 
	 * @throws ComponentGenerationException 
	 */
	void processParse(Element parse) throws ComponentGenerationException {
		List<Element> substituentsAndRoot = OpsinTools.getDescendantElementsWithTagNames(parse, new String[]{SUBSTITUENT_EL, ROOT_EL});

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
		List<Element> groups =  OpsinTools.getDescendantElementsWithTagName(parse, GROUP_EL);

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
	}

	/**
	 * Resolves common ambiguities e.g. tetradeca being 4x10carbon chain rather than 14carbon chain
	 * @param subOrRoot
	 * @throws ComponentGenerationException
	 */
	static void resolveAmbiguities(Element subOrRoot) throws ComponentGenerationException {
		List<Element> multipliers = subOrRoot.getChildElements(MULTIPLIER_EL);
		for (Element apparentMultiplier : multipliers) {
			if (!BASIC_TYPE_VAL.equals(apparentMultiplier.getAttributeValue(TYPE_ATR)) && !VONBAEYER_TYPE_VAL.equals(apparentMultiplier.getAttributeValue(TYPE_ATR))){
				continue;
			}
			int multiplierNum = Integer.parseInt(apparentMultiplier.getAttributeValue(VALUE_ATR));
			Element nextEl = OpsinTools.getNextSibling(apparentMultiplier);
			if (multiplierNum >=3){//detects ambiguous use of things like tetradeca
				if(nextEl !=null){
					if (nextEl.getName().equals(ALKANESTEMCOMPONENT)){//can ignore the trivial alkanes as ambiguity does not exist for them
						int alkaneChainLength = Integer.parseInt(nextEl.getAttributeValue(VALUE_ATR));
						if (alkaneChainLength >=10 && alkaneChainLength > multiplierNum){
							Element isThisALocant = OpsinTools.getPreviousSibling(apparentMultiplier);
							if (isThisALocant == null ||
									!isThisALocant.getName().equals(LOCANT_EL) ||
									MATCH_COMMA.split(isThisALocant.getValue()).length != multiplierNum){
								throw new ComponentGenerationException(apparentMultiplier.getValue() + nextEl.getValue() +" should not have been lexed as two tokens!");
							}
						}
					}
				}
			}

			if (multiplierNum >=4 && nextEl !=null && nextEl.getName().equals(HYDROCARBONFUSEDRINGSYSTEM_EL) && nextEl.getValue().equals("phen") && !"e".equals(nextEl.getAttributeValue(SUBSEQUENTUNSEMANTICTOKEN_ATR))){//deals with tetra phenyl vs tetraphen yl
				Element possibleLocantOrMultiplierOrSuffix = OpsinTools.getNextSibling(nextEl);
				if (possibleLocantOrMultiplierOrSuffix!=null){//null if not used as substituent
					if (possibleLocantOrMultiplierOrSuffix.getName().equals(SUFFIX_EL)){//for phen the aryl substituent, expect an adjacent suffix e.g. phenyl, phenoxy
						Element isThisALocant = OpsinTools.getPreviousSibling(apparentMultiplier);
						if (isThisALocant == null || !isThisALocant.getName().equals(LOCANT_EL) || MATCH_COMMA.split(isThisALocant.getValue()).length != 1){
							String multiplierAndGroup =apparentMultiplier.getValue() + nextEl.getValue();
							throw new ComponentGenerationException(multiplierAndGroup +" should not have been lexed as one token!");
						}
					}
				}
			}
			if (multiplierNum > 4 && !apparentMultiplier.getValue().endsWith("a")){//disambiguate pent oxy and the like. Assume it means pentanoxy rather than 5 oxys
				if (nextEl !=null && nextEl.getName().equals(GROUP_EL)&& matchInlineSuffixesThatAreAlsoGroups.matcher(nextEl.getValue()).matches()){
					throw new ComponentGenerationException(apparentMultiplier.getValue() + nextEl.getValue() +" should have been lexed as [alkane stem, inline suffix], not [multiplier, group]!");
				}
			}
		}

		List<Element> fusions = subOrRoot.getChildElements(FUSION_EL);
		for (Element fusion : fusions) {
			String fusionText = fusion.getValue();
			if (matchNumberLocantsOnlyFusionBracket.matcher(fusionText).matches()){
				Element possibleHWRing = OpsinTools.getNextSiblingIgnoringCertainElements(fusion, new String[]{MULTIPLIER_EL, HETEROATOM_EL});
				if (possibleHWRing !=null && HANTZSCHWIDMAN_SUBTYPE_VAL.equals(possibleHWRing.getAttributeValue(SUBTYPE_ATR))){
					int heteroCount = 0;
					int multiplierValue = 1;
					Element currentElem = OpsinTools.getNextSibling(fusion);
					while(currentElem != null && !currentElem.getName().equals(GROUP_EL)){
						if(currentElem.getName().equals(HETEROATOM_EL)) {
							heteroCount+=multiplierValue;
							multiplierValue =1;
						} else if (currentElem.getName().equals(MULTIPLIER_EL)){
							multiplierValue = Integer.parseInt(currentElem.getAttributeValue(VALUE_ATR));
						}
						currentElem = OpsinTools.getNextSibling(currentElem);
					}
					String[] locants = MATCH_COMMA.split(fusionText.substring(1, fusionText.length()-1));
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
	 * Strips indication of superscript
	 * Strips added hydrogen out of locants
	 * Strips stereochemistry out of locants
	 * Normalises case on greeks to lower case
	 * 
	 * @param subOrRoot
	 * @throws ComponentGenerationException 
	 */
	static void processLocants(Element subOrRoot) throws ComponentGenerationException {
		List<Element> locants = subOrRoot.getChildElements(LOCANT_EL);
		for (Element locant : locants) {
			List<String> individualLocants = splitIntoIndividualLocants(StringTools.removeDashIfPresent(locant.getValue()));
			for (int i = 0; i < individualLocants.size(); i++) {
				String locantText =individualLocants.get(i);
				
				if (locantText.contains("-")){//avoids this regex being invoked typically
					//rearranges locant to the older equivalent form
					Matcher m= matchIUPAC2004ElementLocant.matcher(locantText);
					if (m.matches()){
						locantText = m.group(2) + m.group(1) + m.group(3);
					}
				}
				
				if (Character.isLetter(locantText.charAt(0))){
					//remove indications of superscript as the fact a locant is superscripted can be determined from context e.g. N~1~ ->N1
					Matcher m =  matchSuperscriptedLocant.matcher(locantText);
					if (m.lookingAt()){
						String replacementString = m.group(1) +m.group(2);
						locantText = m.replaceFirst(replacementString);
					}
					if (locantText.length()>=3){
						//convert greeks to lower case
						m =  matchGreek.matcher(locantText);
						while (m.find()) {
							locantText = locantText.substring(0, m.start()) + m.group().toLowerCase() + locantText.substring(m.end());
						}
					}
				}
				char lastChar = locantText.charAt(locantText.length()-1);
				if(lastChar == ')' || lastChar == ']' || lastChar == '}') {				
					//stereochemistry or added hydrogen that result from the application of this locant as a locant for a substituent may be included in brackets after the locant
					
					Matcher m = matchBracketAtEndOfLocant.matcher(locantText);
					if (m.find()){
						String brackettedText = m.group(1);
						if (StringTools.endsWithCaseInsensitive(brackettedText, "H")){
							locantText = m.replaceFirst("");//strip the bracket from the locantText
							//create individual tags for added hydrogen. Examples of bracketed text include "9H" or "2H,7H"
							String[] addedHydrogens = MATCH_COMMA.split(brackettedText);
							for (String addedHydrogen : addedHydrogens) {
								Element addedHydrogenElement=new TokenEl(ADDEDHYDROGEN_EL);
								addedHydrogenElement.addAttribute(new Attribute(LOCANT_ATR, addedHydrogen.substring(0, addedHydrogen.length()-1)));
								OpsinTools.insertBefore(locant, addedHydrogenElement);
							}
							if (locant.getAttribute(TYPE_ATR)==null){
								locant.addAttribute(new Attribute(TYPE_ATR, ADDEDHYDROGENLOCANT_TYPE_VAL));//this locant must not be used as an indirect locant
							}
						}
						else if (StringTools.endsWithCaseInsensitive(brackettedText, "R") || StringTools.endsWithCaseInsensitive(brackettedText, "S")){
							locantText = m.replaceFirst("");//strip the bracket from the locantText
							String rs = brackettedText;
							Element newStereoChemEl = new TokenEl(STEREOCHEMISTRY_EL, "(" + locantText +rs+")");
							newStereoChemEl.addAttribute(new Attribute(TYPE_ATR, STEREOCHEMISTRYBRACKET_TYPE_VAL));
							OpsinTools.insertBefore(locant, newStereoChemEl);
						}
						else if (matchDigit.matcher(brackettedText).matches()){
							//compounds locant e.g. 1(10). Leave as is, it will be handled by the function that handles unsaturation
						}
						else{
							throw new ComponentGenerationException("OPSIN bug: malformed locant text");
						}
					}
					else{
						throw new ComponentGenerationException("OPSIN bug: malformed locant text");
					}

				}
				individualLocants.set(i, locantText);
			}

			locant.setValue(StringTools.stringListToString(individualLocants, ","));

			Element afterLocants = OpsinTools.getNextSibling(locant);
			if(afterLocants == null){
				throw new ComponentGenerationException("Nothing after locant tag: " + locant.toXML());
			}
			
			if (individualLocants.size()==1){
				ifCarbohydrateLocantConvertToAminoAcidStyleLocant(locant);
			}
		}
	}

	/**
	 * Looks for locants of the form (<anotherLocant>)?<multiplier><locant>
	 * and converts them to <modifiedlocant><multiplier>
	 * e.g. 2,4,6 tri O- -->O2,O4,O6 tri
	 * @param locant
	 */
	private static void ifCarbohydrateLocantConvertToAminoAcidStyleLocant(Element locant) {
		if (MATCH_ELEMENT_SYMBOL.matcher(locant.getValue()).matches()){
			Element possibleMultiplier = OpsinTools.getPreviousSibling(locant);
			if (possibleMultiplier!=null && possibleMultiplier.getName().equals(MULTIPLIER_EL)){
				int multiplierValue = Integer.parseInt(possibleMultiplier.getAttributeValue(VALUE_ATR));
				Element possibleOtherLocant = OpsinTools.getPreviousSibling(possibleMultiplier);
				if (possibleOtherLocant!=null){
					String[] locantValues = MATCH_COMMA.split(possibleOtherLocant.getValue());
					if (locantValues.length == Integer.parseInt(possibleMultiplier.getAttributeValue(VALUE_ATR))){
						for (int i = 0; i < locantValues.length; i++) {
							locantValues[i] = locant.getValue() + locantValues[i];
						}
						possibleOtherLocant.setValue(StringTools.arrayToString(locantValues, ","));
						locant.detach();
					}
				}
				else{
					StringBuilder sb = new StringBuilder();
					for (int i = 0; i < multiplierValue-1; i++) {
						sb.append(locant.getValue());
						sb.append(StringTools.multiplyString("'", i));
						sb.append(',');
					}
					sb.append(locant.getValue());
					sb.append(StringTools.multiplyString("'", multiplierValue-1));
					Element newLocant = new TokenEl(LOCANT_EL, sb.toString());
					OpsinTools.insertBefore(possibleMultiplier, newLocant);
					locant.detach();
				}
			}
		}
	}

	/**
	 * Takes a string of locants and splits on commas, but taking into account brackets
	 * e.g. 1,2(1H,2H),3 becomes [1][2(1H,2H)][3]
	 * @param locantString
	 * @return
	 */
	private static List<String> splitIntoIndividualLocants(String locantString) {
		List<String> individualLocants = new ArrayList<String>();
		char[] charArray = locantString.toCharArray();
		boolean inBracket =false;
		int indiceOfLastMatch =0;
		for (int i = 0; i < charArray.length; i++) {
			char c = charArray[i];
			if (c==','){
				if (!inBracket){
					individualLocants.add(locantString.substring(indiceOfLastMatch, i));
					indiceOfLastMatch = i+1;
				}
			}
			else if (c == '(' || c == '[' || c == '{') {
				inBracket =true;
			}
			else if(c == ')' || c == ']' || c == '}') {
				inBracket =false;
			}
		}
		individualLocants.add(locantString.substring(indiceOfLastMatch, charArray.length));
		return individualLocants;
	}

	/**Converts ortho/meta/para into locants
	 * Depending on context para, for example, will either become para or 1,para
	 *
	 * @param subOrRoot
	 * @throws ComponentGenerationException
	 */
	private void convertOrthoMetaParaToLocants(Element subOrRoot) throws ComponentGenerationException{
		List<Element> ompLocants = subOrRoot.getChildElements(ORTHOMETAPARA_EL);
		for (Element ompLocant : ompLocants) {
			String locantText = ompLocant.getValue();
			String firstChar = locantText.substring(0, 1);
			Element afterOmpLocant = OpsinTools.getNextSibling(ompLocant);
			ompLocant.setName(LOCANT_EL);
			ompLocant.addAttribute(new Attribute(TYPE_ATR, ORTHOMETAPARA_TYPE_VAL));
			if(afterOmpLocant.getName().equals(MULTIPLIER_EL) && afterOmpLocant.getAttributeValue(VALUE_ATR).equals("2") ||
					(afterOmpLocant.getAttribute(OUTIDS_ATR)!=null && MATCH_COMMA.split(afterOmpLocant.getAttributeValue(OUTIDS_ATR)).length>1) ) {
				if ("o".equalsIgnoreCase(firstChar)){
					ompLocant.setValue("1,ortho");
				}
				else if ("m".equalsIgnoreCase(firstChar)){
					ompLocant.setValue("1,meta");
				}
				else if ("p".equalsIgnoreCase(firstChar)){
					ompLocant.setValue("1,para");
				}
				else{
					throw new ComponentGenerationException(locantText + " was not identified as being either ortho, meta or para but according to the chemical grammar it should of been");
				}
			}
			else{
				if ("o".equalsIgnoreCase(firstChar)){
					ompLocant.setValue("ortho");
				}
				else if ("m".equalsIgnoreCase(firstChar)){
					ompLocant.setValue("meta");
				}
				else if ("p".equalsIgnoreCase(firstChar)){
					ompLocant.setValue("para");
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
	 */
	private void formAlkaneStemsFromComponents(Element subOrRoot) {
		Deque<Element> alkaneStemComponents =new ArrayDeque<Element>(subOrRoot.getChildElements(ALKANESTEMCOMPONENT));
		while(!alkaneStemComponents.isEmpty()){
			Element alkaneStemComponent = alkaneStemComponents.removeFirst();
			int alkaneChainLength =0;
			StringBuilder alkaneName = new StringBuilder();
			alkaneChainLength += Integer.parseInt(alkaneStemComponent.getAttributeValue(VALUE_ATR));
			alkaneName.append(alkaneStemComponent.getValue());
			while (!alkaneStemComponents.isEmpty() && OpsinTools.getNextSibling(alkaneStemComponent)==alkaneStemComponents.getFirst()) {
				alkaneStemComponent.detach();
				alkaneStemComponent = alkaneStemComponents.removeFirst();
				alkaneChainLength += Integer.parseInt(alkaneStemComponent.getAttributeValue(VALUE_ATR));
				alkaneName.append(alkaneStemComponent.getValue());
			}
			Element alkaneStem = new TokenEl(GROUP_EL, alkaneName.toString());
			alkaneStem.addAttribute(new Attribute(TYPE_ATR, CHAIN_TYPE_VAL));
			alkaneStem.addAttribute(new Attribute(SUBTYPE_ATR, ALKANESTEM_SUBTYPE_VAL));
			alkaneStem.addAttribute(new Attribute(VALUE_ATR, StringTools.multiplyString("C", alkaneChainLength)));
			alkaneStem.addAttribute(new Attribute(USABLEASJOINER_ATR, "yes"));
			StringBuilder labels = new StringBuilder();
			for (int i=1; i<alkaneChainLength; i++) {
				labels.append(i);
				labels.append("/");
			}
			labels.append(alkaneChainLength);
			alkaneStem.addAttribute(new Attribute(LABELS_ATR, labels.toString()));
			OpsinTools.insertAfter(alkaneStemComponent, alkaneStem);
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
		List<Element> alkaneStemModifiers = subOrRoot.getChildElements(ALKANESTEMMODIFIER_EL);
		for(int i=0;i<alkaneStemModifiers.size();i++) {
			Element alkaneStemModifier =alkaneStemModifiers.get(i);
			Element alkane = OpsinTools.getNextSibling(alkaneStemModifier);
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
			int chainLength = alkane.getAttributeValue(VALUE_ATR).length();
			String smiles;
			String labels = NONE_LABELS_VAL;
			if (type.equals("normal")){
				//normal behaviour is default so don't need to do anything
				//n-methyl and n-ethyl contain redundant information and are probably intended to mean N-methyl/N-ethyl
				if ((chainLength==1 || chainLength ==2) && alkaneStemModifier.getValue().equals("n-")){
					Element locant = new TokenEl(LOCANT_EL, "N");
					OpsinTools.insertBefore(alkane, locant);
				}
				continue;
			}
			else if (type.equals("tert")){
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
				boolean suffixPresent = subOrRoot.getChildElements(SUFFIX_EL).size() > 0;
				if (chainLength==3 && !suffixPresent){
					throw new ComponentGenerationException("iso has no meaning without a suffix on an alkane chain of length 3");
				}
				smiles =StringTools.multiplyString("C", chainLength-3) +"C(C)C";
				StringBuilder sb = new StringBuilder();
				for (int c = 1; c <= chainLength - 2; c++) {
					sb.append(c);
					sb.append('/');
				}
				sb.append('/');
				labels = sb.toString();
			}
			else if (type.equals("sec")){
				if (chainLength <3){
					throw new ComponentGenerationException("ChainLength to small for sec modifier, required minLength 3. Found: " +chainLength);
				}
				boolean suffixPresent = subOrRoot.getChildElements(SUFFIX_EL).size() > 0;
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
			alkane.getAttribute(LABELS_ATR).setValue(labels);
		}
	}

	/**Form heterogeneous hydrides/substituents
	 * These are chains of one heteroatom or alternating heteroatoms and are expressed using SMILES
	 * They are typically treated in an analogous way to alkanes
	 * @param subOrRoot The root/substituents
	 * @throws ComponentGenerationException 
	 */
	private void processHeterogenousHydrides(Element subOrRoot) throws ComponentGenerationException  {
		List<Element> multipliers = subOrRoot.getChildElements(MULTIPLIER_EL);
		for (int i = 0; i < multipliers.size(); i++) {
			Element m = multipliers.get(i);
			if (m.getAttributeValue(TYPE_ATR).equals(GROUP_TYPE_VAL)){
				continue;
			}
			Element multipliedElem = OpsinTools.getNextSibling(m);

			if(multipliedElem.getName().equals(GROUP_EL) &&
					multipliedElem.getAttribute(SUBTYPE_ATR)!=null &&
					multipliedElem.getAttributeValue(SUBTYPE_ATR).equals(HETEROSTEM_SUBTYPE_VAL)) {
				int mvalue = Integer.parseInt(m.getAttributeValue(VALUE_ATR));
				
				Element possiblyALocant = OpsinTools.getPreviousSibling(m);//detect rare case where multiplier does not mean form a chain of heteroatoms e.g. something like 1,2-disulfanylpropane
				if(possiblyALocant !=null && possiblyALocant.getName().equals(LOCANT_EL)&& mvalue==MATCH_COMMA.split(possiblyALocant.getValue()).length){
					Element suffix = OpsinTools.getNextSibling(multipliedElem, SUFFIX_EL);
					if (suffix !=null && suffix.getAttributeValue(TYPE_ATR).equals(INLINE_TYPE_VAL)){
						Element possibleMultiplier = OpsinTools.getPreviousSibling(suffix);
						if (!possibleMultiplier.getName().equals(MULTIPLIER_EL)){//NOT something like 3,3'-diselane-1,2-diyl
							continue;
						}
					}
				}
				
				//chain of heteroatoms
				String heteroatomSmiles=multipliedElem.getAttributeValue(VALUE_ATR);
				if (heteroatomSmiles.equals("B") && OpsinTools.getPreviousSibling(m)==null){
					Element possibleUnsaturator = OpsinTools.getNextSibling(multipliedElem);
					if (possibleUnsaturator !=null && possibleUnsaturator.getName().equals(UNSATURATOR_EL) && possibleUnsaturator.getAttributeValue(VALUE_ATR).equals("1")){
						throw new ComponentGenerationException("Polyboranes are not currently supported");
					}
				}
				String smiles = StringTools.multiplyString(heteroatomSmiles, mvalue);
				multipliedElem.getAttribute(VALUE_ATR).setValue(smiles);
				m.detach();
				multipliers.remove(i--);
			}
		}
		for (Element m : multipliers) {
			if (m.getAttributeValue(TYPE_ATR).equals(GROUP_TYPE_VAL)){
				continue;
			}
			Element multipliedElem = OpsinTools.getNextSibling(m);
			if(multipliedElem.getName().equals(HETEROATOM_EL)){
				Element possiblyAnotherHeteroAtom = OpsinTools.getNextSibling(multipliedElem);
				if (possiblyAnotherHeteroAtom !=null && possiblyAnotherHeteroAtom.getName().equals(HETEROATOM_EL)){
					Element possiblyAnUnsaturator = OpsinTools.getNextSiblingIgnoringCertainElements(possiblyAnotherHeteroAtom, new String[]{LOCANT_EL, MULTIPLIER_EL});//typically ane but can be ene or yne e.g. triphosphaza-1,3-diene
					if (possiblyAnUnsaturator !=null && possiblyAnUnsaturator.getName().equals(UNSATURATOR_EL)){
						StringBuilder newGroupName = new StringBuilder(m.getValue());
						newGroupName.append(multipliedElem.getValue());
						newGroupName.append(possiblyAnotherHeteroAtom.getValue());
						//chain of alternating heteroatoms
						if (possiblyAnUnsaturator.getAttributeValue(VALUE_ATR).equals("1")){
							checkForAmbiguityWithHWring(multipliedElem.getAttributeValue(VALUE_ATR), possiblyAnotherHeteroAtom.getAttributeValue(VALUE_ATR));
						}
						int mvalue = Integer.parseInt(m.getAttributeValue(VALUE_ATR));
						StringBuilder smilesSB= new StringBuilder();
						Element possiblyARingFormingEl = OpsinTools.getPreviousSibling(m);
						boolean heteroatomChainWillFormARing = false;
						if (possiblyARingFormingEl!=null && (possiblyARingFormingEl.getName().equals(CYCLO_EL) || possiblyARingFormingEl.getName().equals(VONBAEYER_EL) || possiblyARingFormingEl.getName().equals(SPIRO_EL))){
							heteroatomChainWillFormARing=true;
							//will be cyclised later.
							for (int j = 0; j < mvalue; j++) {
								smilesSB.append(possiblyAnotherHeteroAtom.getAttributeValue(VALUE_ATR));
								smilesSB.append(multipliedElem.getAttributeValue(VALUE_ATR));
							}
						}
						else{
							for (int j = 0; j < mvalue -1; j++) {
								smilesSB.append(multipliedElem.getAttributeValue(VALUE_ATR));
								smilesSB.append(possiblyAnotherHeteroAtom.getAttributeValue(VALUE_ATR));
							}
							smilesSB.append(multipliedElem.getAttributeValue(VALUE_ATR));
						}
						String smiles =smilesSB.toString();
						smiles = matchHdigit.matcher(smiles).replaceAll("H?");//hydrogen count will be determined by standard valency
						multipliedElem.detach();

						Element addedGroup = new TokenEl(GROUP_EL, newGroupName.toString());
						addedGroup.addAttribute(new Attribute(VALUE_ATR, smiles));
						addedGroup.addAttribute(new Attribute(TYPE_ATR, CHAIN_TYPE_VAL));
						addedGroup.addAttribute(new Attribute(SUBTYPE_ATR, HETEROSTEM_SUBTYPE_VAL));
						if (!heteroatomChainWillFormARing){
							addedGroup.addAttribute(new Attribute(USABLEASJOINER_ATR, "yes"));
						}
						OpsinTools.insertAfter(possiblyAnotherHeteroAtom, addedGroup);

						possiblyAnotherHeteroAtom.detach();
						m.detach();
					}
					else if (possiblyAnUnsaturator!=null && possiblyAnUnsaturator.getValue().equals("an") && HANTZSCHWIDMAN_SUBTYPE_VAL.equals(possiblyAnUnsaturator.getAttributeValue(SUBTYPE_ATR))){
						//check for HWring that should be interpreted as a heterogenous hydride
						boolean foundLocantIndicatingHwRingHeteroatomPositions =false;//allow formally incorrect HW ring systems if they have locants
						Element possibleLocant = OpsinTools.getPreviousSibling(m);
						if (possibleLocant !=null && possibleLocant.getName().equals(LOCANT_EL)){
							int expected = Integer.parseInt(m.getAttributeValue(VALUE_ATR)) + 1;
							if (expected == MATCH_COMMA.split(possibleLocant.getValue()).length){
								foundLocantIndicatingHwRingHeteroatomPositions = true;
							}
						}
						if (!foundLocantIndicatingHwRingHeteroatomPositions){
							checkForAmbiguityWithHeterogenousHydride(multipliedElem.getAttributeValue(VALUE_ATR), possiblyAnotherHeteroAtom.getAttributeValue(VALUE_ATR));
						}
					}
				}
			}
		}
	}

	/**
	 * Throws an exception if the given heteroatoms could be part of a valid Hantzch-widman ring
	 * For this to be true the first heteroatom must be higher priority than the second
	 * and the second must be compatible with a HW ane stem
	 * @param firstHeteroAtomSMILES
	 * @param secondHeteroAtomSMILES
	 * @throws ComponentGenerationException 
	 */
	private void checkForAmbiguityWithHWring(String firstHeteroAtomSMILES, String secondHeteroAtomSMILES) throws ComponentGenerationException {
		Matcher m = MATCH_ELEMENT_SYMBOL.matcher(firstHeteroAtomSMILES);
		if (!m.find()){
			throw new ComponentGenerationException("Failed to extract element from heteroatom");
		}
		String atom1Element = m.group();
		
		m = MATCH_ELEMENT_SYMBOL.matcher(secondHeteroAtomSMILES);
		if (!m.find()){
			throw new ComponentGenerationException("Failed to extract element from heteroatom");
		}
		String atom2Element = m.group();
		if (AtomProperties.elementToHwPriority.get(atom1Element) > AtomProperties.elementToHwPriority.get(atom2Element)){
			if (atom2Element.equals("O") || atom2Element.equals("S")  || atom2Element.equals("Se") || atom2Element.equals("Te") 
					|| atom2Element.equals("Bi")  || atom2Element.equals("Hg")){
				if (!hasSiorGeorSnorPb(atom1Element, atom2Element)){
					throw new ComponentGenerationException("Hantzch-widman ring misparsed as a heterogeneous hydride with alternating atoms");
				}
			}
		}
	}

	/**
	 * Are either of the elements Si/Ge/Sn/Pb
	 * @param atom1Element
	 * @param atom2Element
	 * @return
	 */
	private boolean hasSiorGeorSnorPb(String atom1Element, String atom2Element) {
		return (atom1Element.equals("Si") || atom1Element.equals("Ge") || atom1Element.equals("Sn") ||atom1Element.equals("Pb")
				|| atom2Element.equals("Si") || atom2Element.equals("Ge") || atom2Element.equals("Sn") ||atom2Element.equals("Pb"));
	}
	
	/**
	 * Throws an exception if the given heteroatoms could be part of a heterogenous hydride
	 * For this to be true the second heteroatom must be higher priority than the first
	 * @param firstHeteroAtomSMILES
	 * @param secondHeteroAtomSMILES
	 * @throws ComponentGenerationException 
	 */
	private void checkForAmbiguityWithHeterogenousHydride(String firstHeteroAtomSMILES, String secondHeteroAtomSMILES) throws ComponentGenerationException {
		Matcher m = MATCH_ELEMENT_SYMBOL.matcher(firstHeteroAtomSMILES);
		if (!m.find()){
			throw new ComponentGenerationException("Failed to extract element from heteroatom");
		}
		String atom1Element = m.group();
		
		m = MATCH_ELEMENT_SYMBOL.matcher(secondHeteroAtomSMILES);
		if (!m.find()){
			throw new ComponentGenerationException("Failed to extract element from heteroatom");
		}
		String atom2Element = m.group();
		if (AtomProperties.elementToHwPriority.get(atom2Element) > AtomProperties.elementToHwPriority.get(atom1Element)){
			throw new ComponentGenerationException("heterogeneous hydride with alternating atoms misparsed as a Hantzch-widman ring");
		}
	}

	/** Handle indicated hydrogen  e.g. 1H- in 1H-pyrrole
	 *
	 * @param subOrRoot The substituent/root to looks for indicated hydrogens in.
	 * @throws ComponentGenerationException 
	 */
	private void processIndicatedHydrogens(Element subOrRoot) throws ComponentGenerationException {
		List<Element> indicatedHydrogens = subOrRoot.getChildElements(INDICATEDHYDROGEN_EL);
		for (Element indicatedHydrogenGroup : indicatedHydrogens) {
			String txt = StringTools.removeDashIfPresent(indicatedHydrogenGroup.getValue());
			if (!StringTools.endsWithCaseInsensitive(txt, "h")){//remove brackets if they are present
				txt = txt.substring(1, txt.length()-1);
			}
			String[] hydrogenLocants =MATCH_COMMA.split(txt);
            for (String hydrogenLocant : hydrogenLocants) {
                if (StringTools.endsWithCaseInsensitive(hydrogenLocant, "h")) {
                    Element indicatedHydrogenEl = new TokenEl(INDICATEDHYDROGEN_EL);
                    indicatedHydrogenEl.addAttribute(new Attribute(LOCANT_ATR, hydrogenLocant.substring(0, hydrogenLocant.length() - 1)));
                    OpsinTools.insertBefore(indicatedHydrogenGroup, indicatedHydrogenEl);
                }
                else{
                	throw new ComponentGenerationException("OPSIN Bug: malformed indicated hydrogen element!");
                }
            }
			indicatedHydrogenGroup.detach();
		}
	}

	/** Handles stereoChemistry in brackets: R/Z/E/Z/a/alpha/b/beta and cis/trans
	 *  Will assign a locant to a stereoChemistry element if one was specified/available
	 *
	 * @param subOrRoot The substituent/root to looks for stereoChemistry in.
	 * @throws ComponentGenerationException
	 */
	void processStereochemistry(Element subOrRoot) throws ComponentGenerationException {
		List<Element> stereoChemistryElements = subOrRoot.getChildElements(STEREOCHEMISTRY_EL);
		for (Element stereoChemistryElement : stereoChemistryElements) {
			if (stereoChemistryElement.getAttributeValue(TYPE_ATR).equals(STEREOCHEMISTRYBRACKET_TYPE_VAL)){
				processStereochemistryBracket(stereoChemistryElement);
			}
			else if (stereoChemistryElement.getAttributeValue(TYPE_ATR).equals(CISORTRANS_TYPE_VAL)){
				assignLocantUsingPreviousElementIfPresent(stereoChemistryElement);//assign a locant if one is directly before the cis/trans
			}
			else if (stereoChemistryElement.getAttributeValue(TYPE_ATR).equals(E_OR_Z_TYPE_VAL)){
				stereoChemistryElement.addAttribute(new Attribute(VALUE_ATR, stereoChemistryElement.getValue().toUpperCase()));
				assignLocantUsingPreviousElementIfPresent(stereoChemistryElement);//assign a locant if one is directly before the E/Z
			}
			else if (stereoChemistryElement.getAttributeValue(TYPE_ATR).equals(ENDO_EXO_SYN_ANTI_TYPE_VAL)){
				processLocantAssigningForEndoExoSynAnti(stereoChemistryElement);//assign a locant if one is directly before the endo/exo/syn/anti. Don't neccesarily detach it
			}
			else if (stereoChemistryElement.getAttributeValue(TYPE_ATR).equals(ALPHA_OR_BETA_TYPE_VAL)){
				processUnbracketedAlphaBetaStereochemistry(stereoChemistryElement);
			}
			else if (stereoChemistryElement.getAttributeValue(TYPE_ATR).equals(RELATIVECISTRANS_TYPE_VAL)){
				processRelativeCisTrans(stereoChemistryElement);
			}
		}
	}

	private void processStereochemistryBracket(Element stereoChemistryElement) throws ComponentGenerationException {
		String txt = stereoChemistryElement.getValue();
		if (StringTools.startsWithCaseInsensitive(txt, "rel-")){
			txt = txt.substring(4);
		}
		txt = StringTools.removeDashIfPresent(txt);
		Matcher starMatcher = matchStar.matcher(txt);
		txt = starMatcher.replaceAll("");
		if (!StringTools.startsWithCaseInsensitive(txt, "rac") && txt.length() > 0){//if txt is just "rel-" then it will be length 0 at this point
			List<String> stereoChemistryDescriptors = splitStereoBracketIntoDescriptors(txt);
		    for (String stereoChemistryDescriptor : stereoChemistryDescriptors) {
		        Matcher m = matchStereochemistry.matcher(stereoChemistryDescriptor);
		        if (m.matches()){
		        	if (!m.group(2).equals("RS") && !m.group(2).equals("SR")){
		                Element stereoChemEl = new TokenEl(STEREOCHEMISTRY_EL, stereoChemistryDescriptor);
		                String locantVal = m.group(1);
		                if (locantVal.length() > 0){
		                    stereoChemEl.addAttribute(new Attribute(LOCANT_ATR, StringTools.removeDashIfPresent(locantVal)));
		                }
		                OpsinTools.insertBefore(stereoChemistryElement, stereoChemEl);
		                if (matchRS.matcher(m.group(2)).matches()) {
		                    stereoChemEl.addAttribute(new Attribute(TYPE_ATR, R_OR_S_TYPE_VAL));
		                    stereoChemEl.addAttribute(new Attribute(VALUE_ATR, m.group(2).toUpperCase()));
		                } else if (matchEZ.matcher(m.group(2)).matches()) {
		                    stereoChemEl.addAttribute(new Attribute(TYPE_ATR, E_OR_Z_TYPE_VAL));
		                    stereoChemEl.addAttribute(new Attribute(VALUE_ATR, m.group(2).toUpperCase()));
		                } else if (matchAlphaBetaStereochem.matcher(m.group(2)).matches()){
		                	stereoChemEl.addAttribute(new Attribute(TYPE_ATR, ALPHA_OR_BETA_TYPE_VAL));
		                	if (Character.toLowerCase(m.group(2).charAt(0)) == 'a'){
		                    	stereoChemEl.addAttribute(new Attribute(VALUE_ATR, "alpha"));
		                	}
		                	else if (Character.toLowerCase(m.group(2).charAt(0)) == 'b'){
		                    	stereoChemEl.addAttribute(new Attribute(VALUE_ATR, "beta"));
		                	}
		                 	else if (Character.toLowerCase(m.group(2).charAt(0)) == 'x'){
		                    	stereoChemEl.addAttribute(new Attribute(VALUE_ATR, "xi"));
		                	}
		                	else{
		                		throw new ComponentGenerationException("Malformed alpha/beta stereochemistry element: " + stereoChemistryElement.getValue());
		                	}
		        	 	} else if (matchCisTrans.matcher(m.group(2)).matches()) {
		                    stereoChemEl.addAttribute(new Attribute(TYPE_ATR, CISORTRANS_TYPE_VAL));
		                    stereoChemEl.addAttribute(new Attribute(VALUE_ATR, m.group(2).toLowerCase()));
		        	 	} else if (matchEndoExoSynAnti.matcher(m.group(2)).matches()) {
		                    stereoChemEl.addAttribute(new Attribute(TYPE_ATR, ENDO_EXO_SYN_ANTI_TYPE_VAL));
		                    stereoChemEl.addAttribute(new Attribute(VALUE_ATR, m.group(2).toLowerCase()));
		        	 	} else {
		                    throw new ComponentGenerationException("Malformed stereochemistry element: " + stereoChemistryElement.getValue());
		                }
		        	 		
		        	}
		        } else {
		            throw new ComponentGenerationException("Malformed stereochemistry element: " + stereoChemistryElement.getValue());
		        }
		    }
		}
		stereoChemistryElement.detach();
	}

	private List<String> splitStereoBracketIntoDescriptors(String stereoBracket) {
		List<String> stereoDescriptors = new ArrayList<String>();
		StringBuilder sb = new StringBuilder();
		//ignore first and last character (opening and closing bracket)
		for (int i = 1, l = stereoBracket.length() - 1; i < l; i++) {
			char ch = stereoBracket.charAt(i);
			if (ch ==','){
				stereoDescriptors.add(sb.toString());
				sb.setLength(0);
			}
			else if (ch == '-'){
				if (matchStereochemistry.matcher(sb.toString()).matches()){
					//delimiter between stereochemistry
					stereoDescriptors.add(sb.toString());
					sb.setLength(0);
				}
				else{
					//locanted stereochemistry term
					sb.append(ch);
				}
			}
			else{
				sb.append(ch);
			}
		}
		stereoDescriptors.add(sb.toString());
		return stereoDescriptors;
	}

	private void assignLocantUsingPreviousElementIfPresent(Element stereoChemistryElement) {
		Element possibleLocant = OpsinTools.getPrevious(stereoChemistryElement);
		if (possibleLocant !=null && possibleLocant.getName().equals(LOCANT_EL) && MATCH_COMMA.split(possibleLocant.getValue()).length==1){
			stereoChemistryElement.addAttribute(new Attribute(LOCANT_ATR, possibleLocant.getValue()));
			possibleLocant.detach();
		}
	}
	
	private void processLocantAssigningForEndoExoSynAnti(Element stereoChemistryElement) {
		Element possibleLocant = OpsinTools.getPrevious(stereoChemistryElement);
		if (possibleLocant !=null && possibleLocant.getName().equals(LOCANT_EL) && MATCH_COMMA.split(possibleLocant.getValue()).length==1){
			stereoChemistryElement.addAttribute(new Attribute(LOCANT_ATR, possibleLocant.getValue()));
			Element group = OpsinTools.getNextSibling(stereoChemistryElement, GROUP_EL);
			if (group != null && 
					(CYCLICUNSATURABLEHYDROCARBON_SUBTYPE_VAL.equals(group.getAttributeValue(SUBTYPE_ATR))
						|| OpsinTools.getPreviousSibling(group).getName().equals(VONBAEYER_EL))){
				//detach locant only if we're sure it has no other meaning
				//typically locants in front of endo/exo/syn/anti also indicate the position of a susbtituent/suffix e.g. 3-exo-amino
				possibleLocant.detach();
			}
		}
	}

	private void processUnbracketedAlphaBetaStereochemistry(Element stereoChemistryElement) throws ComponentGenerationException {
		String txt = StringTools.removeDashIfPresent(stereoChemistryElement.getValue());
		String[] stereoChemistryDescriptors = MATCH_COMMA.split(txt);
		List<String> locants = new ArrayList<String>();
		boolean createLocantsEl =false;
		for (String stereoChemistryDescriptor : stereoChemistryDescriptors) {
			Matcher digitMatcher  = matchDigit.matcher(stereoChemistryDescriptor);
			if (digitMatcher.lookingAt()){
				String locant = digitMatcher.group();
				String possibleAlphaBeta = digitMatcher.replaceAll("");
		        locants.add(locant);
				Matcher alphaBetaMatcher = matchAlphaBetaStereochem.matcher(possibleAlphaBeta);
				if (alphaBetaMatcher.matches()){
		            Element stereoChemEl = new TokenEl(STEREOCHEMISTRY_EL, stereoChemistryDescriptor);
		            stereoChemEl.addAttribute(new Attribute(LOCANT_ATR, locant));
		            OpsinTools.insertBefore(stereoChemistryElement, stereoChemEl);
		           	stereoChemEl.addAttribute(new Attribute(TYPE_ATR, ALPHA_OR_BETA_TYPE_VAL));
		        	if (Character.toLowerCase(possibleAlphaBeta.charAt(0)) == 'a'){
		            	stereoChemEl.addAttribute(new Attribute(VALUE_ATR, "alpha"));
		        	}
		        	else if (Character.toLowerCase(possibleAlphaBeta.charAt(0)) == 'b'){
		            	stereoChemEl.addAttribute(new Attribute(VALUE_ATR, "beta"));
		        	}
		         	else if (Character.toLowerCase(possibleAlphaBeta.charAt(0)) == 'x'){
		            	stereoChemEl.addAttribute(new Attribute(VALUE_ATR, "xi"));
		        	}
		        	else{
		        		throw new ComponentGenerationException("Malformed alpha/beta stereochemistry element: " + stereoChemistryElement.getValue());
		        	}
		        }
				else{
					createLocantsEl =true;
				}
			}
		}
		if (!createLocantsEl){
			//create locants unless a group supporting alpha/beta stereochem is within this substituent/root
			createLocantsEl =true;
			List<Element> groups = OpsinTools.getNextSiblingsOfType(stereoChemistryElement, GROUP_EL);
			for (Element group : groups) {
				if (group.getAttributeValue(ALPHABETACLOCKWISEATOMORDERING_ATR)!=null){
					createLocantsEl=false;
					break;
				}
			}
		}
		
		if (createLocantsEl){
			Element newLocantEl = new TokenEl(LOCANT_EL, StringTools.stringListToString(locants, ","));
			OpsinTools.insertAfter(stereoChemistryElement, newLocantEl);
		}
		stereoChemistryElement.detach();
	}
	
	private void processRelativeCisTrans(Element stereoChemistryElement) {
		String value = StringTools.removeDashIfPresent(stereoChemistryElement.getValue());
		StringBuilder sb = new StringBuilder();
		String[] terms = MATCH_COMMA.split(value);
		for (String term : terms) {
			if (term.startsWith("c-")|| term.startsWith("t-") || term.startsWith("r-")){
				if (sb.length() > 0){
					sb.append(',');
				}
				sb.append(term.substring(2));
			}
			else{
				throw new RuntimeException("Malformed relativeCisTrans element");
			}
		}
		Element locantEl = new TokenEl(LOCANT_EL, sb.toString());
		OpsinTools.insertAfter(stereoChemistryElement, locantEl);
	}

	/**
	 * Looks for "suffixPrefix" and assigns their value them as an attribute of an adjacent suffix
	 * @param subOrRoot
	 * @throws ComponentGenerationException
	 */
	private void processSuffixPrefixes(Element subOrRoot) throws ComponentGenerationException {
		List<Element> suffixPrefixes =  subOrRoot.getChildElements(SUFFIXPREFIX_EL);
		for (Element suffixPrefix : suffixPrefixes) {
			Element suffix = OpsinTools.getNextSibling(suffixPrefix);
			if (suffix==null || ! suffix.getName().equals(SUFFIX_EL)){
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
		List<Element> infixes = subOrRoot.getChildElements(INFIX_EL);
		for (Element infix : infixes) {
			Element suffix = OpsinTools.getNextSiblingIgnoringCertainElements(infix, new String[]{INFIX_EL, SUFFIXPREFIX_EL, MULTIPLIER_EL});
			if (suffix ==null || !suffix.getName().equals(SUFFIX_EL)){
				throw new ComponentGenerationException("No suffix found next next to infix: "+ infix.getValue());
			}
			List<String> currentInfixInformation;
			if (suffix.getAttribute(INFIX_ATR)==null){
				suffix.addAttribute(new Attribute(INFIX_ATR, ""));
				currentInfixInformation = new ArrayList<String>();
			}
			else{
				currentInfixInformation = StringTools.arrayToList(MATCH_SEMICOLON.split(suffix.getAttributeValue(INFIX_ATR)));
			}
			String infixValue =infix.getAttributeValue(VALUE_ATR);
			currentInfixInformation.add(infixValue);
			Element possibleMultiplier = OpsinTools.getPreviousSibling(infix);
			Element possibleBracket;
			boolean multiplierKnownToIndicateInfixMultiplicationPresent =false;
			if (possibleMultiplier.getName().equals(MULTIPLIER_EL)){
				//suffix prefix present so multiplier must indicate infix replacement
				Element possibleSuffixPrefix = OpsinTools.getPreviousSiblingIgnoringCertainElements(infix, new String[]{MULTIPLIER_EL, INFIX_EL});
				if (possibleSuffixPrefix!=null && possibleSuffixPrefix.getName().equals(SUFFIXPREFIX_EL)){
					multiplierKnownToIndicateInfixMultiplicationPresent =true;
				}
				Element elementBeforeMultiplier = OpsinTools.getPreviousSibling(possibleMultiplier);
				//double multiplier indicates multiple suffixes which all have their infix multiplied
				//if currentInfixInformation contains more than 1 entry it contains information from an infix from before the multiplier so the interpretation of the multiplier as a suffix multiplier is impossible
				if (elementBeforeMultiplier.getName().equals(MULTIPLIER_EL) || currentInfixInformation.size() > 1){
					multiplierKnownToIndicateInfixMultiplicationPresent =true;
				}
				possibleBracket = elementBeforeMultiplier;
			}
			else{
				possibleBracket=possibleMultiplier;
				possibleMultiplier=null;
				infix.detach();
			}
			if (possibleBracket.getName().equals(STRUCTURALOPENBRACKET_EL)){
				Element bracket = OpsinTools.getNextSibling(suffix);
				if (!bracket.getName().equals(STRUCTURALCLOSEBRACKET_EL)){
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
			else if (multiplierKnownToIndicateInfixMultiplicationPresent){//multiplier unambiguously means multiplication of the infix
				int multiplierVal = Integer.parseInt(possibleMultiplier.getAttributeValue(VALUE_ATR));
				for (int i = 1; i < multiplierVal; i++) {
					currentInfixInformation.add(infixValue);
				}
				possibleMultiplier.detach();
				infix.detach();
			}
			else if (possibleMultiplier!=null && GROUP_TYPE_VAL.equals(possibleMultiplier.getAttributeValue(TYPE_ATR))){//e.g. ethanbisthioic acid == ethanbis(thioic acid)
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
		List<Element> lambdaConventionEls = subOrRoot.getChildElements(LAMBDACONVENTION_EL);
		boolean fusedRingPresent = false;
		if (lambdaConventionEls.size()>0){
			if (subOrRoot.getChildElements(GROUP_EL).size()>1){
				fusedRingPresent = true;
			}
		}
		for (Element lambdaConventionEl : lambdaConventionEls) {
			boolean frontLocantsExpected =false;//Is the lambdaConvention el followed by benz/benzo of a fused ring system (these have front locants which correspond to the final fused rings numbering) or by a polycylicspiro system
			String[] lambdaValues = MATCH_COMMA.split(StringTools.removeDashIfPresent(lambdaConventionEl.getValue()));
			Element possibleHeteroatomOrMultiplier = OpsinTools.getNextSibling(lambdaConventionEl);
			int heteroCount = 0;
			int multiplierValue = 1;
			while(possibleHeteroatomOrMultiplier != null){
				if(possibleHeteroatomOrMultiplier.getName().equals(HETEROATOM_EL)) {
					heteroCount+=multiplierValue;
					multiplierValue =1;
				} else if (possibleHeteroatomOrMultiplier.getName().equals(MULTIPLIER_EL)){
					multiplierValue = Integer.parseInt(possibleHeteroatomOrMultiplier.getAttributeValue(VALUE_ATR));
				}
				else{
					break;
				}
				possibleHeteroatomOrMultiplier = OpsinTools.getNextSibling(possibleHeteroatomOrMultiplier);
			}
			boolean assignLambdasToHeteroAtoms =false;
			if (lambdaValues.length==heteroCount){//heteroatom and number of locants +lambdas must match
				if (fusedRingPresent && possibleHeteroatomOrMultiplier!=null && possibleHeteroatomOrMultiplier.getName().equals(GROUP_EL) && possibleHeteroatomOrMultiplier.getAttributeValue(SUBTYPE_ATR).equals(HANTZSCHWIDMAN_SUBTYPE_VAL)){
					//You must not set the locants of a HW system which forms a component of a fused ring system. The locant specified corresponds to the complete fused ring system.
				}
				else{
					assignLambdasToHeteroAtoms =true;
				}
			}
			else if (possibleHeteroatomOrMultiplier!=null && ((heteroCount==0 && OpsinTools.getNextSibling(lambdaConventionEl).equals(possibleHeteroatomOrMultiplier) &&
					fusedRingPresent && possibleHeteroatomOrMultiplier.getName().equals(GROUP_EL) &&
					(possibleHeteroatomOrMultiplier.getValue().equals("benzo") || possibleHeteroatomOrMultiplier.getValue().equals("benz"))
					&& !OpsinTools.getNextSibling(possibleHeteroatomOrMultiplier).getName().equals(FUSION_EL)
					&& !OpsinTools.getNextSibling(possibleHeteroatomOrMultiplier).getName().equals(LOCANT_EL))
					|| (possibleHeteroatomOrMultiplier.getName().equals(POLYCYCLICSPIRO_EL) && 
							(possibleHeteroatomOrMultiplier.getAttributeValue(VALUE_ATR).equals("spirobi")|| possibleHeteroatomOrMultiplier.getAttributeValue(VALUE_ATR).equals("spiroter"))))){
				frontLocantsExpected = true;//a benzo fused ring e.g. 1lambda4,3-benzothiazole or a symmetrical poly cyclic spiro system
			}
			List<Element> heteroAtoms = new ArrayList<Element>();//contains the heteroatoms to apply the lambda values too. Can be empty if the values are applied to a group directly rather than to a heteroatom
			if (assignLambdasToHeteroAtoms){//populate heteroAtoms, multiplied heteroatoms are multiplied out
				Element multiplier = null;
				Element heteroatomOrMultiplier = OpsinTools.getNextSibling(lambdaConventionEl);
				while(heteroatomOrMultiplier != null){
					if(heteroatomOrMultiplier.getName().equals(HETEROATOM_EL)) {
						heteroAtoms.add(heteroatomOrMultiplier);
						if (multiplier!=null){
							for (int i = 1; i < Integer.parseInt(multiplier.getAttributeValue(VALUE_ATR)); i++) {
								Element newHeteroAtom = heteroatomOrMultiplier.copy();
								OpsinTools.insertBefore(heteroatomOrMultiplier, newHeteroAtom);
								heteroAtoms.add(newHeteroAtom);
							}
							multiplier.detach();
							multiplier=null;
						}
					} else if (heteroatomOrMultiplier.getName().equals(MULTIPLIER_EL)){
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
					heteroatomOrMultiplier = OpsinTools.getNextSibling(heteroatomOrMultiplier);
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
						Element newLambda = new TokenEl(LAMBDACONVENTION_EL);
						newLambda.addAttribute(valencyChange);
						if (locantAtr!=null){
							newLambda.addAttribute(locantAtr);
						}
						OpsinTools.insertBefore(lambdaConventionEl, newLambda);
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
				lambdaConventionEl.setName(LOCANT_EL);
				lambdaConventionEl.setValue(StringTools.arrayToString(lambdaValues, ","));
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
			List<Element> children = sub.getChildElements();
			for(int i=0; i<children.size(); i++) {
				Element child = children.get(i);
				if(child.getName().equals(OPENBRACKET_EL)) {
					if(openBracket == null) {
						openBracket = child;
					}
					blevel++;
				} else if (child.getName().equals(CLOSEBRACKET_EL)) {
					blevel--;
					if(blevel == 0) {
						closeBracket = child;
						Element bracket = structureBrackets(openBracket, closeBracket);
						while(findAndStructureBrackets(OpsinTools.getDescendantElementsWithTagName(bracket, SUBSTITUENT_EL)));
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
		Element bracket = new GroupingEl(BRACKET_EL);
		OpsinTools.insertBefore(openBracket.getParent(), bracket);
		/* Pick up everything in the substituent before the bracket*/
		while(!openBracket.getParent().getChild(0).equals(openBracket)) {
			Element n = openBracket.getParent().getChild(0);
			n.detach();
			bracket.appendChild(n);
		}
		/* Pick up all elements from the one with the open bracket,
		 * to the one with the close bracket, inclusive.
		 */
		Element currentEl = openBracket.getParent();
		while(!currentEl.equals(closeBracket.getParent())) {
			Element nextEl = OpsinTools.getNextSibling(currentEl);
			currentEl.detach();
			bracket.appendChild(currentEl);
			currentEl = nextEl;
			if (currentEl==null){
				throw new ComponentGenerationException("Brackets within a word do not match!");
			}
		}
		currentEl.detach();
		bracket.appendChild(currentEl);
		/* Pick up elements after the close bracket */
		currentEl = OpsinTools.getNextSibling(closeBracket);
		while(currentEl != null) {
			Element nextEl = OpsinTools.getNextSibling(currentEl);
			currentEl.detach();
			bracket.appendChild(currentEl);
			currentEl = nextEl;
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
		List<Element> annulens = subOrRoot.getChildElements(ANNULEN_EL);
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

			Element group =new TokenEl(GROUP_EL, annulenValue);
			group.addAttribute(new Attribute(VALUE_ATR, SMILES));
			group.addAttribute(new Attribute(TYPE_ATR, RING_TYPE_VAL));
			group.addAttribute(new Attribute(SUBTYPE_ATR, ARYLGROUP_SUBTYPE_VAL));
			annulen.getParent().replaceChild(annulen, group);
		}

		List<Element> hydrocarbonFRSystems = subOrRoot.getChildElements(HYDROCARBONFUSEDRINGSYSTEM_EL);
		for (Element hydrocarbonFRSystem : hydrocarbonFRSystems) {
			Element multiplier = OpsinTools.getPreviousSibling(hydrocarbonFRSystem);
			if(multiplier != null && multiplier.getName().equals(MULTIPLIER_EL)) {
				int multiplierValue =Integer.parseInt(multiplier.getAttributeValue(VALUE_ATR));
				String classOfHydrocarbonFRSystem =hydrocarbonFRSystem.getAttributeValue(VALUE_ATR);
				StringBuilder smilesSB= new StringBuilder();
				if (classOfHydrocarbonFRSystem.equals("polyacene")){
					if (multiplierValue <=3){
						throw new ComponentGenerationException("Invalid polyacene");
					}
					smilesSB.append("c1ccc");
					for (int j = 2; j <= multiplierValue; j++) {
						smilesSB.append("c");
						smilesSB.append(ringClosure(j));
						smilesSB.append("c");
					}
					smilesSB.append("ccc");
					for (int j = multiplierValue; j >2; j--) {
						smilesSB.append("c");
						smilesSB.append(ringClosure(j));
						smilesSB.append("c");
					}
					smilesSB.append("c12");
				}else if (classOfHydrocarbonFRSystem.equals("polyaphene")){
					if (multiplierValue <=3){
						throw new ComponentGenerationException("Invalid polyaphene");
					}
					smilesSB.append("c1ccc");

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
						smilesSB.append("c");
						smilesSB.append(ringClosure(ringOpeningCounter++));
						smilesSB.append("c");
					}

					for (int j = 1; j <= ringsOnPlane; j++) {
						smilesSB.append("cc");
						smilesSB.append(ringClosure(ringOpeningCounter++));
					}
					smilesSB.append("ccc");
					ringOpeningCounter--;
					for (int j = 1; j <= ringsOnPlane; j++) {
						smilesSB.append("cc");
						smilesSB.append(ringClosure(ringOpeningCounter--));
					}
					for (int j = 1; j < ringsAbovePlane; j++) {
						smilesSB.append("c");
						smilesSB.append(ringClosure(ringOpeningCounter--));
						smilesSB.append("c");
					}

					smilesSB.append("c12");
				} else if (classOfHydrocarbonFRSystem.equals("polyalene")){
					if (multiplierValue <5){
						throw new ComponentGenerationException("Invalid polyalene");
					}
					smilesSB.append("c1");
					for (int j = 3; j < multiplierValue; j++) {
						smilesSB.append("c");
					}
					smilesSB.append("c2");
					for (int j = 3; j <= multiplierValue; j++) {
						smilesSB.append("c");
					}
					smilesSB.append("c12");
				} else if (classOfHydrocarbonFRSystem.equals("polyphenylene")){
					if (multiplierValue <2){
						throw new ComponentGenerationException("Invalid polyphenylene");
					}
					smilesSB.append("c1cccc2");
					for (int j = 1; j < multiplierValue; j++) {
						smilesSB.append("c3ccccc3");
					}
					smilesSB.append("c12");
				} else if (classOfHydrocarbonFRSystem.equals("polynaphthylene")){
					if (multiplierValue <3){
						throw new ComponentGenerationException("Invalid polynaphthylene");
					}
					smilesSB.append("c1cccc2cc3");
					for (int j = 1; j < multiplierValue; j++) {
						smilesSB.append("c4cc5ccccc5cc4");
					}
					smilesSB.append("c3cc12");
				} else if (classOfHydrocarbonFRSystem.equals("polyhelicene")){
					if (multiplierValue <4){
						throw new ComponentGenerationException("Invalid polyhelicene");
					}
					smilesSB.append("c1c");
					int ringOpeningCounter=2;
					for (int j = 1; j < multiplierValue; j++) {
						smilesSB.append("ccc");
						smilesSB.append(ringClosure(ringOpeningCounter++));
					}
					smilesSB.append("cccc");
					ringOpeningCounter--;
					for (int j = 2; j < multiplierValue; j++) {
						smilesSB.append("c");
						smilesSB.append(ringClosure(ringOpeningCounter--));
					}
					smilesSB.append("c12");
				}

				else{
					throw new ComponentGenerationException("Unknown semi-trivially named hydrocarbon fused ring system");
				}
				Element newGroup =new TokenEl(GROUP_EL, multiplier.getValue() + hydrocarbonFRSystem.getValue());
				newGroup.addAttribute(new Attribute(VALUE_ATR, smilesSB.toString()));
				newGroup.addAttribute(new Attribute(LABELS_ATR, FUSEDRING_LABELS_VAL));
				newGroup.addAttribute(new Attribute(TYPE_ATR, RING_TYPE_VAL));
				newGroup.addAttribute(new Attribute(SUBTYPE_ATR, HYDROCARBONFUSEDRINGSYSTEM_EL));
				hydrocarbonFRSystem.getParent().replaceChild(hydrocarbonFRSystem, newGroup);
				multiplier.detach();
			}
			else{
				throw new ComponentGenerationException("Invalid semi-trivially named hydrocarbon fused ring system");
			}
		}
	}
	
	/**
	 * Handles irregular suffixes. e.g. Quinone and ylene
	 * @param subOrRoot
	 * @throws ComponentGenerationException 
	 */
	private void handleSuffixIrregularities(Element subOrRoot) throws ComponentGenerationException {
		List<Element> suffixes = subOrRoot.getChildElements(SUFFIX_EL);
		for (Element suffix : suffixes) {
			String suffixValue = suffix.getValue();
			if (suffixValue.equals("ic") || suffixValue.equals("ous")){
				if (!n2sConfig.allowInterpretationOfAcidsWithoutTheWordAcid()) {
					Element next = OpsinTools.getNext(suffix);
					if (next == null){
						throw new ComponentGenerationException("\"acid\" not found after " +suffixValue);
					}
				}
			}
			// convert quinone to dione
			else if (suffixValue.equals("quinone") || suffixValue.equals("quinon")){
				suffix.removeAttribute(suffix.getAttribute(ADDITIONALVALUE_ATR));
				suffix.setValue("one");
				Element multiplier = OpsinTools.getPreviousSibling(suffix);
				if (multiplier.getName().equals(MULTIPLIER_EL)){
					Attribute multVal = multiplier.getAttribute(VALUE_ATR);
					int newMultiplier = Integer.parseInt(multVal.getValue()) * 2;
					multVal.setValue(String.valueOf(newMultiplier));
				}
				else{
					multiplier = new TokenEl(MULTIPLIER_EL, "di");
					multiplier.addAttribute(new Attribute(VALUE_ATR, "2"));
					OpsinTools.insertBefore(suffix, multiplier);
				}
			}
			else if (suffixValue.equals("ylene") || suffixValue.equals("ylen")){
				suffix.removeAttribute(suffix.getAttribute(ADDITIONALVALUE_ATR));
				suffix.setValue("yl");
				Element alk = OpsinTools.getPreviousSibling(suffix, GROUP_EL);
				if (alk.getAttribute(USABLEASJOINER_ATR)!=null){
					alk.removeAttribute(alk.getAttribute(USABLEASJOINER_ATR));
				}
				Element multiplier = new TokenEl(MULTIPLIER_EL, "di");
				multiplier.addAttribute(new Attribute(VALUE_ATR, "2"));
				OpsinTools.insertBefore(suffix, multiplier);
			}
			else if (suffixValue.equals("ylium") &&//disambiguate between ylium the charge modifying suffix and ylium the acylium suffix
					"acylium".equals(suffix.getAttributeValue(VALUE_ATR)) &&
					suffix.getAttribute(SUFFIXPREFIX_ATR)==null &&
					suffix.getAttribute(INFIX_ATR)==null){
				Element group = OpsinTools.getPreviousSibling(suffix, GROUP_EL);
				if (group==null || (!ACIDSTEM_TYPE_VAL.equals(group.getAttributeValue(TYPE_ATR)) &&
						!CHALCOGENACIDSTEM_TYPE_VAL.equals(group.getAttributeValue(TYPE_ATR)) &&
								!NONCARBOXYLICACID_TYPE_VAL.equals(group.getAttributeValue(TYPE_ATR)))){
					Element beforeSuffix = OpsinTools.getPreviousSibling(suffix);
					String o = beforeSuffix.getAttributeValue(SUBSEQUENTUNSEMANTICTOKEN_ATR);
					if (o ==null || !StringTools.endsWithCaseInsensitive(o, "o")){
						if (group!=null && ARYLSUBSTITUENT_SUBTYPE_VAL.equals(group.getAttributeValue(SUBTYPE_ATR))){
							//contracted form for removal of hydride e.g. 9-Anthrylium
							suffix.getAttribute(VALUE_ATR).setValue("ylium");
							suffix.getAttribute(TYPE_ATR).setValue(CHARGE_TYPE_VAL);
							suffix.removeAttribute(suffix.getAttribute(SUBTYPE_ATR));
						}
						else{
							throw new ComponentGenerationException("ylium is intended to be the removal of H- in this context not the formation of an acylium ion");
						}
					}
				}
			}
		}
	}

	/**
	 * Looks for alkaneStems followed by a bridge forming 'o' and makes them fused ring bridge elements
	 * @param group
	 */
	private void detectAlkaneFusedRingBridges(Element group) {
		if (ALKANESTEM_SUBTYPE_VAL.equals(group.getAttributeValue(SUBTYPE_ATR))){
			Element possibleBridgeFormer = OpsinTools.getNextSiblingIgnoringCertainElements(group, new String[]{UNSATURATOR_EL});
			if(possibleBridgeFormer != null && possibleBridgeFormer.getName().equals(BRIDGEFORMINGO_EL)){
				possibleBridgeFormer.detach();
				group.setName(FUSEDRINGBRIDGE_EL);
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
		Element previous = OpsinTools.getPreviousSibling(group);
		if(previous != null) {
			String previousElType = previous.getName();
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

		Element multiplier = OpsinTools.getPreviousSibling(spiroEl);
		int numberOfSpiros = 1;
		if (multiplier != null && multiplier.getName().equals(MULTIPLIER_EL)){
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
				StringBuilder superScriptedNumber = new StringBuilder();
				for (int j = 1; j < elements.length; j++){//may be more than one non digit as there are many ways of indicating superscripts
					superScriptedNumber.append(elements[j]);
				}
				spiroDescriptors[i][1] = Integer.parseInt(superScriptedNumber.toString());
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
		Element multiplier = OpsinTools.getPreviousSibling(vonBaeyerBracketEl);
		int numberOfRings=Integer.parseInt(multiplier.getAttributeValue(VALUE_ATR));
		multiplier.detach();

		int alkylChainLength;
		Deque<String> elementSymbolArray = new ArrayDeque<String>();
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

		StringBuilder smilesSB = new StringBuilder();
		int atomCounter=1;
		int bridgeCounter=1;
		//add standard bridges
		for (HashMap<String, Integer> bridge : bridges) {
			if (bridgeCounter==1){
				smilesSB.append(elementSymbolArray.removeFirst());
				smilesSB.append("1");
				if (bridgeLocations.get(atomCounter)!=null){
					for (Integer bridgeAtomLabel : bridgeLocations.get(atomCounter)) {
						smilesSB.append(ringClosure(bridgeAtomLabel));
					}
				}
				smilesSB.append("(");
			}
			int bridgeLength =bridge.get("Bridge Length");

			for (int i = 0; i < bridgeLength; i++) {
				atomCounter++;
				smilesSB.append(elementSymbolArray.removeFirst());
				if (bridgeLocations.get(atomCounter)!=null){
					for (Integer bridgeAtomLabel : bridgeLocations.get(atomCounter)) {
						smilesSB.append(ringClosure(bridgeAtomLabel));
					}
				}
			}
			if (bridgeCounter==1){
				atomCounter++;
				smilesSB.append(elementSymbolArray.removeFirst());
				smilesSB.append("2");
				if (bridgeLocations.get(atomCounter)!=null){
					for (Integer bridgeAtomLabel : bridgeLocations.get(atomCounter)) {
						smilesSB.append(ringClosure(bridgeAtomLabel));
					}
				}
			}
			if (bridgeCounter==2){
				smilesSB.append("1)");
			}
			if (bridgeCounter==3){
				smilesSB.append("2");
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
				smilesSB.append(".");
				for (int i = 0; i < bridgeLength; i++) {
					atomCounter++;
					smilesSB.append(elementSymbolArray.removeFirst());
					if (i==0){
						smilesSB.append(ringClosure(bridge.get("AtomId_Larger_Label")));
					}
					if (bridgeLocations.get(atomCounter)!=null){
						for (Integer bridgeAtomLabel : bridgeLocations.get(atomCounter)) {
							smilesSB.append(ringClosure(bridgeAtomLabel));
						}
					}
				}
				smilesSB.append(ringClosure(bridge.get("AtomId_Smaller_Label")));
			}
			if (dependantSecondaryBridges.size() >0 && dependantSecondaryBridges.size()==secondaryBridges.size()){
				throw new ComponentGenerationException("Unable to resolve all dependant bridges!!!");
			}
			secondaryBridges=dependantSecondaryBridges;
		}
		while(dependantSecondaryBridges.size() > 0);

		chainEl.getAttribute(VALUE_ATR).setValue(smilesSB.toString());
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
     * @throws ComponentGenerationException
	 */
	private void handleGroupIrregularities(Element group) throws ComponentGenerationException {
		String groupValue =group.getValue();
		
		if (!n2sConfig.allowInterpretationOfAcidsWithoutTheWordAcid()) {
			if (group.getAttribute(FUNCTIONALIDS_ATR) !=null && (groupValue.endsWith("ic") || groupValue.endsWith("ous"))){
				Element next = OpsinTools.getNext(group);
				if (next == null){
					throw new ComponentGenerationException("\"acid\" not found after " +groupValue);
				}
			}
		}
		
		if(groupValue.equals("thiophen") || groupValue.equals("selenophen") || groupValue.equals("tellurophen")) {//thiophenol is generally phenol with an O replaced with S not thiophene with a hydroxy
			Element possibleSuffix = OpsinTools.getNextSibling(group);
			if (!"e".equals(group.getAttributeValue(SUBSEQUENTUNSEMANTICTOKEN_ATR)) && possibleSuffix !=null && possibleSuffix.getName().equals(SUFFIX_EL)) {
				if (possibleSuffix.getValue().startsWith("ol")){
					Element isThisALocant = OpsinTools.getPreviousSibling(group);
					if (isThisALocant == null || !isThisALocant.getName().equals(LOCANT_EL) || MATCH_COMMA.split(isThisALocant.getValue()).length != 1){
						throw new ComponentGenerationException(groupValue + "ol has been incorrectly interpreted as "+ groupValue+", ol instead of phenol with the oxgen replaced");
					}
				}
			}
		}
		else if (groupValue.equals("methylene") || groupValue.equals("methylen")) {//e.g. 3,4-methylenedioxyphenyl
			Element nextSub = OpsinTools.getNextSibling(group.getParent());
			if (nextSub !=null && nextSub.getName().equals(SUBSTITUENT_EL) && OpsinTools.getNextSibling(group)==null 
					&& (OpsinTools.getPreviousSibling(group)==null || !OpsinTools.getPreviousSibling(group).getName().equals(MULTIPLIER_EL))){//not trimethylenedioxy
				List<Element> children = nextSub.getChildElements();
				if (children.size() >=2 && children.get(0).getValue().equals("di")&& children.get(1).getValue().equals("oxy")){
					group.setValue(groupValue + "dioxy");
					group.getAttribute(VALUE_ATR).setValue("C(O)O");
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
						OpsinTools.insertAfter(group, children.get(i));
					}
				}
			}
		}
		else if (groupValue.equals("ethylene") || groupValue.equals("ethylen")) {
			Element previous = OpsinTools.getPreviousSibling(group);
			if (previous != null && previous.getName().equals(MULTIPLIER_EL)){
				int multiplierValue = Integer.parseInt(previous.getAttributeValue(VALUE_ATR));
				Element possibleRoot = OpsinTools.getNextSibling(group.getParent());
				if (possibleRoot==null && OpsinTools.getParentWordRule(group).getAttributeValue(WORDRULE_ATR).equals(WordRule.glycol.toString())){//e.g. dodecaethylene glycol
					StringBuilder smiles = new StringBuilder("CC");
					for (int i = 1; i < multiplierValue; i++) {
						smiles.append("OCC");
					}
					group.getAttribute(OUTIDS_ATR).setValue("1," +Integer.toString(3*(multiplierValue-1) +2));
					group.getAttribute(VALUE_ATR).setValue(smiles.toString());
					previous.detach();
					if (group.getAttribute(LABELS_ATR)!=null){//use numeric numbering
						group.getAttribute(LABELS_ATR).setValue(NUMERIC_LABELS_VAL);
					}
					else{
						group.addAttribute(new Attribute(LABELS_ATR, NUMERIC_LABELS_VAL));
					}
				}
				else if (possibleRoot!=null && possibleRoot.getName().equals(ROOT_EL)){
					List<Element> children = possibleRoot.getChildElements();
					if (children.size()==2){
						Element amineMultiplier =children.get(0);
						Element amine =children.get(1);
						if (amineMultiplier.getName().equals(MULTIPLIER_EL) && (amine.getValue().equals("amine") || amine.getValue().equals("amin"))){//e.g. Triethylenetetramine
							if (Integer.parseInt(amineMultiplier.getAttributeValue(VALUE_ATR))!=multiplierValue +1){
								throw new ComponentGenerationException("Invalid polyethylene amine!");
							}
							StringBuilder smiles = new StringBuilder();
							for (int i = 0; i < multiplierValue; i++) {
								smiles.append("NCC");
							}
							smiles.append("N");
							group.removeAttribute(group.getAttribute(OUTIDS_ATR));
							group.getAttribute(VALUE_ATR).setValue(smiles.toString());
							previous.detach();
							possibleRoot.detach();
							group.getParent().setName(ROOT_EL);
							if (group.getAttribute(LABELS_ATR)!=null){//use numeric numbering
								group.getAttribute(LABELS_ATR).setValue(NUMERIC_LABELS_VAL);
							}
							else{
								group.addAttribute(new Attribute(LABELS_ATR, NUMERIC_LABELS_VAL));
							}
						}
					}
				}
			}
			else{
				Element nextSub = OpsinTools.getNextSibling(group.getParent());
				if (nextSub !=null && nextSub.getName().equals(SUBSTITUENT_EL) && OpsinTools.getNextSibling(group)==null){
					List<Element> children = nextSub.getChildElements();
					if (children.size() >=2 && children.get(0).getValue().equals("di")&& children.get(1).getValue().equals("oxy")){
						group.setValue(groupValue + "dioxy");
						group.getAttribute(VALUE_ATR).setValue("C(O)CO");
						group.getAttribute(OUTIDS_ATR).setValue("2,4");
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
							OpsinTools.insertAfter(group, children.get(i));
						}
					}
				}
			}
		}
		else if (groupValue.equals("propylene") || groupValue.equals("propylen")) {
			Element previous = OpsinTools.getPreviousSibling(group);
			if (previous!=null && previous.getName().equals(MULTIPLIER_EL)){
				int multiplierValue = Integer.parseInt(previous.getAttributeValue(VALUE_ATR));
				Element possibleRoot = OpsinTools.getNextSibling(group.getParent());
				if (possibleRoot==null && OpsinTools.getParentWordRule(group).getAttributeValue(WORDRULE_ATR).equals(WordRule.glycol.toString())){//e.g. dodecaethylene glycol
					StringBuilder smiles =new StringBuilder("CCC");
					for (int i = 1; i < multiplierValue; i++) {
						smiles.append("OC(C)C");
					}
					group.getAttribute(OUTIDS_ATR).setValue("2," +Integer.toString(4*(multiplierValue-1) +3));
					group.getAttribute(VALUE_ATR).setValue(smiles.toString());
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

		//acridone (not codified), anthrone, phenanthrone and xanthone have the one at position 9 by default
		else if (groupValue.equals("anthr") || groupValue.equals("phenanthr") || groupValue.equals("acrid") ||
				groupValue.equals("xanth") || groupValue.equals("thioxanth") || groupValue.equals("selenoxanth")|| groupValue.equals("telluroxanth")|| groupValue.equals("xanthen")) {
			Element possibleLocant = OpsinTools.getPreviousSibling(group);
			if (possibleLocant==null || !possibleLocant.getName().equals(LOCANT_EL)){//only need to give one a locant of 9 if no locant currently present
				Element possibleSuffix = OpsinTools.getNextSibling(group);
				if (possibleSuffix!=null && "one".equals(possibleSuffix.getAttributeValue(VALUE_ATR))){
					//Rule C-315.2
					Element newLocant =new TokenEl(LOCANT_EL, "9");
					OpsinTools.insertBefore(possibleSuffix, newLocant);
					Element newAddedHydrogen = new TokenEl(ADDEDHYDROGEN_EL);
					newAddedHydrogen.addAttribute(new Attribute(LOCANT_ATR, "10"));
					OpsinTools.insertBefore(newLocant, newAddedHydrogen);
				}
				else if (possibleSuffix!=null && possibleSuffix.getName().equals(SUFFIX_EL) &&
					groupValue.equals("xanth") || groupValue.equals("thioxanth") || groupValue.equals("selenoxanth")|| groupValue.equals("telluroxanth")){
					//diasambiguate between xanthate/xanthic acid and xanthene
					String suffixVal = possibleSuffix.getAttributeValue(VALUE_ATR);
					if (suffixVal.equals("ic") || suffixVal.equals("ate")){
						throw new ComponentGenerationException(groupValue + possibleSuffix.getValue() +" is not a derivative of xanthene");
					}
				}
			}
		}
		else if (groupValue.equals("phospho")){//is this the organic meaning (P(=O)=O) or biochemical meaning (P(=O)(O)O)
			Element substituent = group.getParent();
			Element nextGroup = OpsinTools.getNextGroup(substituent);
			if (nextGroup != null){
				String type = nextGroup.getAttributeValue(TYPE_ATR);
				String subType = nextGroup.getAttributeValue(SUBTYPE_ATR);
				if (type.equals(AMINOACID_TYPE_VAL) || 
						type.equals(CARBOHYDRATE_TYPE_VAL) ||
						BIOCHEMICAL_SUBTYPE_VAL.equals(subType) ||
						(YLFORACYL_SUBTYPE_VAL.equals(subType) &&
								("glycol".equals(nextGroup.getValue()) || "diglycol".equals(nextGroup.getValue()))
						)
					) {
					group.getAttribute(VALUE_ATR).setValue("-P(=O)(O)O");
					group.addAttribute(new Attribute(USABLEASJOINER_ATR, "yes"));
				}
			}
		}
		else if (groupValue.equals("aspart") || groupValue.equals("glutam")){//aspartyl and glutamyl typically mean alpha-aspartyl/alpha-glutamyl
			if (group.getAttributeValue(SUBTYPE_ATR).equals(ENDINIC_SUBTYPE_VAL)){
				Element yl = OpsinTools.getNextSibling(group);
				if (yl.getAttributeValue(VALUE_ATR).equals("yl")){
					group.removeAttribute(group.getAttribute(SUFFIXAPPLIESTO_ATR));
					if (groupValue.equals("aspart")){
						group.getAttribute("labels").setValue("/2,alpha/3,beta/4,gamma///1/");
						group.getAttribute(VALUE_ATR).setValue("N[C@@H](CC(O)=O)C=O");
						group.addAttribute(new Attribute(OUTIDS_ATR, "7"));
					}
					else {
						group.getAttribute("labels").setValue("/2,alpha/3,beta/4,gamma/5,delta///1/");
						group.getAttribute(VALUE_ATR).setValue("N[C@@H](CCC(O)=O)C=O");
						group.addAttribute(new Attribute(OUTIDS_ATR, "8"));
					}
					yl.detach();
				}
			}
		}
		else if (groupValue.equals("hydrogen")){
			Element hydrogenParentEl = group.getParent();
			Element nextSubOrRoot = OpsinTools.getNextSibling(hydrogenParentEl);
			if (nextSubOrRoot!=null){
				Element possibleSuitableAteGroup = nextSubOrRoot.getChild(0);
				if (!possibleSuitableAteGroup.getName().equals(GROUP_EL) || !NONCARBOXYLICACID_TYPE_VAL.equals(possibleSuitableAteGroup.getAttributeValue(TYPE_ATR))){
					throw new ComponentGenerationException("Hydrogen is not meant as a substituent in this context!");
				}
				Element possibleMultiplier = OpsinTools.getPreviousSibling(group);
				String multiplier = "1";
				if (possibleMultiplier!=null && possibleMultiplier.getName().equals(MULTIPLIER_EL)){
					multiplier = possibleMultiplier.getAttributeValue(VALUE_ATR);
					possibleMultiplier.detach();
				}
				possibleSuitableAteGroup.addAttribute(new Attribute(NUMBEROFFUNCTIONALATOMSTOREMOVE_ATR, multiplier));
				group.detach();
				List<Element> childrenToMove = hydrogenParentEl.getChildElements();
				for (int i = childrenToMove.size() -1 ; i >=0; i--) {
					childrenToMove.get(i).detach();
					nextSubOrRoot.insertChild(childrenToMove.get(i), 0);
				}
				hydrogenParentEl.detach();
			}
		}
		else if (groupValue.equals("acryl")){
			if (SIMPLESUBSTITUENT_SUBTYPE_VAL.equals(group.getAttributeValue(SUBTYPE_ATR))){
				Element nextEl = OpsinTools.getNext(group);
				if (nextEl!=null && nextEl.getValue().equals("amid")){
					throw new ComponentGenerationException("amide in acrylamide is not [NH2-]");
				}
			}
		}
		else if (groupValue.equals("azo") || groupValue.equals("azoxy") || groupValue.equals("nno-azoxy") || groupValue.equals("non-azoxy") || groupValue.equals("onn-azoxy") || groupValue.equals("diazoamino") || groupValue.equals("hydrazo") ){
			Element enclosingSub = group.getParent();
			Element next = OpsinTools.getNextSiblingIgnoringCertainElements(enclosingSub, new String[]{HYPHEN_EL});
			if (next==null && OpsinTools.getPreviousSibling(enclosingSub) == null){//e.g. [(E)-NNO-azoxy]benzene
				next = OpsinTools.getNextSiblingIgnoringCertainElements(enclosingSub.getParent(), new String[]{HYPHEN_EL});
			}
			if (next!=null && next.getName().equals(ROOT_EL)){
				if (!(next.getChild(0).getName().equals(MULTIPLIER_EL))){
					List<Element> suffixes = next.getChildElements(SUFFIX_EL);
					if (suffixes.size()==0){//only case without locants is handled so far. suffixes only apply to one of the fragments rather than both!!!
						Element newMultiplier = new TokenEl(MULTIPLIER_EL);
						newMultiplier.addAttribute(new Attribute(VALUE_ATR, "2"));
						next.insertChild(newMultiplier, 0);
						Element interSubstituentHyphen = OpsinTools.getPrevious(group);
						if (interSubstituentHyphen!=null && !interSubstituentHyphen.getName().equals(HYPHEN_EL)){//prevent implicit bracketting
							OpsinTools.insertAfter(interSubstituentHyphen, new TokenEl(HYPHEN_EL));
						}
					}
				}
			}
		}
		else if (groupValue.equals("coenzyme a") || groupValue.equals("coa")){
			Element enclosingSubOrRoot = group.getParent();
			Element previous = OpsinTools.getPreviousSibling(enclosingSubOrRoot);
			if (previous!=null){
				List<Element> groups = OpsinTools.getDescendantElementsWithTagName(previous, GROUP_EL);
				if (groups.size()>0){
					Element possibleAcid = groups.get(groups.size()-1);
					if (ACIDSTEM_TYPE_VAL.equals(possibleAcid.getAttributeValue(TYPE_ATR))){
						if (possibleAcid.getAttribute(SUFFIXAPPLIESTO_ATR)!=null){//multi acid. yl should be one oyl and the rest carboxylic acids
							Element suffix = OpsinTools.getNextSibling(possibleAcid, SUFFIX_EL);
							if (suffix.getAttribute(ADDITIONALVALUE_ATR)==null){
								suffix.addAttribute(new Attribute(ADDITIONALVALUE_ATR, "ic"));
							}
						}
						String subType = possibleAcid.getAttributeValue(SUBTYPE_ATR);
						if (subType.equals(YLFORYL_SUBTYPE_VAL) || subType.equals(YLFORNOTHING_SUBTYPE_VAL)){
							possibleAcid.getAttribute(SUBTYPE_ATR).setValue(YLFORACYL_SUBTYPE_VAL);//yl always  means an acyl when next to coenzyme A
						}
					}
				}
			}
			//locanted substitution onto Coenzyme A is rarely intended, so put it in a bracket to disfavour it
			Element newBracket = new GroupingEl(BRACKET_EL);
			OpsinTools.insertAfter(enclosingSubOrRoot, newBracket);
			enclosingSubOrRoot.detach();
			newBracket.appendChild(enclosingSubOrRoot);
		}
		else if (groupValue.equals("sphinganine") || groupValue.equals("icosasphinganine") || groupValue.equals("eicosasphinganine") || groupValue.equals("phytosphingosine") || groupValue.equals("sphingosine")){
			Element enclosingSubOrRoot = group.getParent();
			Element previous = OpsinTools.getPreviousSibling(enclosingSubOrRoot);
			if (previous!=null){
				List<Element> groups = OpsinTools.getDescendantElementsWithTagName(previous, GROUP_EL);
				if (groups.size()>0){
					Element possibleAcid = groups.get(groups.size()-1);
					if (ALKANESTEM_SUBTYPE_VAL.equals(possibleAcid.getAttributeValue(SUBTYPE_ATR))){
						List<Element> inlineSuffixes = OpsinTools.getChildElementsWithTagNameAndAttribute(possibleAcid.getParent(), SUFFIX_EL, TYPE_ATR, INLINE_TYPE_VAL);
						if (inlineSuffixes.size()==1 && inlineSuffixes.get(0).getAttributeValue(VALUE_ATR).equals("yl")){
							inlineSuffixes.get(0).getAttribute(VALUE_ATR).setValue("oyl");//yl on a systematic acid next to a fatty acid means acyl
							//c.f. Nomenclature of Lipids 1976, Appendix A, note a
						}
					}
				}
			}
		}
		else if (groupValue.equals("sel")){
			//check that it is not "selenium"
			if (HETEROSTEM_SUBTYPE_VAL.equals(group.getAttributeValue(SUBTYPE_ATR)) && group.getAttribute(SUBSEQUENTUNSEMANTICTOKEN_ATR) ==null){
				Element unsaturator = OpsinTools.getNextSibling(group);
				if (unsaturator !=null && unsaturator.getName().equals(UNSATURATOR_EL) && unsaturator.getValue().equals("en") && group.getAttribute(SUBSEQUENTUNSEMANTICTOKEN_ATR) ==null){
					Element ium = OpsinTools.getNextSibling(unsaturator);
					if (ium !=null && ium.getName().equals(SUFFIX_EL) && ium.getValue().equals("ium")){
						throw new ComponentGenerationException("<multiplier>selenium does not indicate a chain of selenium atoms with a double bond and a positive charge");
					}
				}
			}
		}
		else if ((groupValue.equals("keto") || groupValue.equals("aldehydo")) && SIMPLESUBSTITUENT_SUBTYPE_VAL.equals(group.getAttributeValue(SUBTYPE_ATR))){
			//check for case where this is specifying the open chain form of a ketose/aldose
			Element previousEl = OpsinTools.getPreviousSibling(group);
			if (previousEl ==null || !previousEl.getName().equals(LOCANT_EL) || groupValue.equals("aldehydo")){
				Element parentSubstituent = group.getParent();
				Element nextSubOrRoot = OpsinTools.getNextSibling(parentSubstituent);
				Element parentOfCarbohydate = nextSubOrRoot;
				Element carbohydrate = null;
				while (parentOfCarbohydate != null){
					Element possibleCarbohydrate = parentOfCarbohydate.getFirstChildElement(GROUP_EL);
					if (possibleCarbohydrate !=null && possibleCarbohydrate.getAttributeValue(TYPE_ATR).equals(CARBOHYDRATE_TYPE_VAL)){
						carbohydrate = possibleCarbohydrate;
						break;
					}
					parentOfCarbohydate = OpsinTools.getNextSibling(parentOfCarbohydate);
				}
				if (carbohydrate != null) {
					if (parentOfCarbohydate.getChildElements(CARBOHYDRATERINGSIZE_EL).size() > 0){
						throw new ComponentGenerationException("Carbohydrate has a specified ring size but " + groupValue + " indicates the open chain form!");
					}
					group.detach();
					List<Element> childrenToMove = parentSubstituent.getChildElements();
					for (int i = childrenToMove.size() -1 ; i >=0; i--) {
						Element el = childrenToMove.get(i);
						if (!el.getName().equals(HYPHEN_EL)){
							el.detach();
							nextSubOrRoot.insertChild(el, 0);
						}
					}
					parentSubstituent.detach();
					String carbohydrateAdditionValue = carbohydrate.getAttributeValue(ADDITIONALVALUE_ATR);
					//OPSIN assumes a few trival names are more likely to describe the cyclic form. additonalValue contains the SMILES for the acyclic form
					if (carbohydrateAdditionValue != null){
						if (carbohydrateAdditionValue.equals("n/a")){
							throw new ComponentGenerationException(carbohydrate.getValue() + " can only describe the cyclic form  but " + groupValue + " indicates the open chain form!");
						}
						carbohydrate.getAttribute(VALUE_ATR).setValue(carbohydrateAdditionValue);
					}
				}
				else if (groupValue.equals("aldehydo")){
					throw new ComponentGenerationException("aldehydo is only a valid prefix when it precedes a carbohydrate!");
				}
			}
		}
		else if (groupValue.equals("bor") || groupValue.equals("antimon") 
				|| groupValue.equals("arsen") || groupValue.equals("phosphor") || groupValue.equals("phosphate")
				|| groupValue.equals("silicicacid") || groupValue.equals("silicic acid")
				|| groupValue.equals("silicate")){//fluoroboric acid/fluoroborate are trivial rather than systematic; tetra(fooyl)borate is inorganic
			Element suffix = null;
			Boolean isAcid = null;
			if (groupValue.endsWith("acid")){
				if (OpsinTools.getNext(group) == null){
					isAcid = true;
				}
			}
			else if (groupValue.endsWith("ate")){
				if (OpsinTools.getNext(group) == null){
					isAcid = false;
				}
			}
			else{
				suffix = OpsinTools.getNextSibling(group);
				if (suffix != null && suffix.getName().equals(SUFFIX_EL) &&
						suffix.getAttribute(INFIX_ATR) == null && OpsinTools.getNext(suffix) == null){
					String suffixValue = suffix.getAttributeValue(VALUE_ATR);
					if (suffixValue.equals("ic")){
						isAcid = true;
					}
					else if (suffixValue.equals("ate")){
						isAcid = false;
					}
				}
			}
			if (isAcid != null){//check for inorganic interpretation
				Element substituent = OpsinTools.getPreviousSibling(group.getParent());
				if (substituent !=null && (substituent.getName().equals(SUBSTITUENT_EL) || substituent.getName().equals(BRACKET_EL))){
					List<Element> children = substituent.getChildElements();
					Element firstChild = children.get(0);
					boolean matched = false;
					if (children.size() ==1 && firstChild.getName().equals(GROUP_EL) && (firstChild.getValue().equals("fluoro") || firstChild.getValue().equals("fluor"))){
						if (groupValue.equals("bor")) {
							group.getAttribute(VALUE_ATR).setValue(isAcid ? "F[B-](F)(F)F.[H+]" : "F[B-](F)(F)F");
							matched = true;
						}
						else if (groupValue.equals("antimon")) {
							group.getAttribute(VALUE_ATR).setValue(isAcid ? "F[Sb-](F)(F)(F)(F)F.[H+]" : "F[Sb-](F)(F)(F)(F)F");
							matched = true;
						}
						else if (groupValue.startsWith("silicic")) {
							group.getAttribute(VALUE_ATR).setValue(isAcid ? "F[Si|6-2](F)(F)(F)(F)F.[H+].[H+]" : "F[Si|6-2](F)(F)(F)(F)F");
							matched = true;
						}
						if (matched) {
							substituent.detach();
						}
					}
					else if (firstChild.getName().equals(MULTIPLIER_EL)) {
						String multiplierVal = firstChild.getAttributeValue(VALUE_ATR);
						
						if (groupValue.equals("bor")){
							if (multiplierVal.equals("4") || (multiplierVal.equals("3") && OpsinTools.getPreviousSibling(substituent) != null)) {
								//tri case allows organotrifluoroborates
								group.getAttribute(VALUE_ATR).setValue(isAcid ? "[B-].[H+]" :"[B-]");
								matched = true;
							}
						}
						else if (groupValue.equals("antimon") && multiplierVal.equals("6")) {
							group.getAttribute(VALUE_ATR).setValue(isAcid ? "[Sb-].[H+]" :"[Sb-]");
							matched = true;
						}
						else if (groupValue.equals("arsen") && multiplierVal.equals("6")) {
							group.getAttribute(VALUE_ATR).setValue(isAcid ? "[As-].[H+]" :"[As-]");
							matched = true;
						}
						else if (groupValue.startsWith("phosph") && multiplierVal.equals("6")) {
							group.getAttribute(VALUE_ATR).setValue(isAcid ? "[P-].[H+]" :"[P-]");
							matched = true;
						}
						else if (groupValue.startsWith("silic") && multiplierVal.equals("6")) {
							group.getAttribute(VALUE_ATR).setValue(isAcid ? "[Si|6-2].[H+].[H+]" :"[Si|6-2]");
							matched = true;
						}
					}
					if (matched) {
						group.getAttribute(TYPE_ATR).setValue(SIMPLEGROUP_TYPE_VAL);
						group.getAttribute(SUBTYPE_ATR).setValue(SIMPLEGROUP_SUBTYPE_VAL);
						
						Attribute usableAsJoiner = group.getAttribute(USABLEASJOINER_ATR);
						if (usableAsJoiner != null){
							group.removeAttribute(usableAsJoiner);
						}
						Attribute acceptsAdditiveBonds = group.getAttribute(ACCEPTSADDITIVEBONDS_ATR);
						if (acceptsAdditiveBonds != null){
							group.removeAttribute(acceptsAdditiveBonds);
						}
						Attribute functionalIds = group.getAttribute(FUNCTIONALIDS_ATR);
						if (functionalIds != null){
							group.removeAttribute(functionalIds);
						}
						
						if (suffix != null){
							suffix.detach();
						}
					}
				}
			}
		}
	}
}

