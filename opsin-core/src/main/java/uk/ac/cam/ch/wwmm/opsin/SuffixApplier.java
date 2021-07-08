package uk.ac.cam.ch.wwmm.opsin;

import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import uk.ac.cam.ch.wwmm.opsin.IsotopeSpecificationParser.IsotopeSpecification;


class SuffixApplier {
	
	private final BuildState state;
	private final SuffixRules suffixRules;

	SuffixApplier(BuildState state, SuffixRules suffixRules) {
		this.state = state;
		this.suffixRules = suffixRules;
	}

	/**
	 * Does suffixApplicability.xml have an entry for this group type? 
	 * @param groupType
	 * @return
	 */
	boolean isGroupTypeWithSpecificSuffixRules(String groupType){
		return suffixRules.isGroupTypeWithSpecificSuffixRules(groupType);
	}
	

	/**Process the effects of suffixes upon a fragment. 
	 * Unlocanted non-terminal suffixes are not attached yet. All other suffix effects are performed
	 * @param group The group element for the fragment to which the suffixes will be added
	 * @param suffixes The suffix elements for a fragment.
	 * @throws StructureBuildingException If the suffixes can't be resolved properly.
	 * @throws ComponentGenerationException
	 */
	void resolveSuffixes(Element group, List<Element> suffixes) throws StructureBuildingException, ComponentGenerationException {
		Fragment frag = group.getFrag();
		List<Atom> atomList = frag.getAtomList();//this instance of atomList will not change even once suffixes are merged into the fragment
		String groupType = frag.getType();
		String subgroupType = frag.getSubType();
		String suffixTypeToUse = isGroupTypeWithSpecificSuffixRules(groupType) ? groupType : STANDARDGROUP_TYPE_VAL;
		
		List<Fragment> associatedSuffixFrags = state.xmlSuffixMap.get(group);
		if (associatedSuffixFrags != null) {//null for non-final group in polycyclic spiro systems
			associatedSuffixFrags.clear();
		}
		Map<String, List<Element>> suffixValToSuffixes = new LinkedHashMap<>();//effectively undoes the effect of multiplying out suffixes
		for (Element suffix : suffixes) {
			String suffixValue = suffix.getAttributeValue(VALUE_ATR);
			List<Element> suffixesWithThisVal = suffixValToSuffixes.get(suffixValue);
			if (suffixesWithThisVal == null) {
				suffixesWithThisVal = new ArrayList<>();
				suffixValToSuffixes.put(suffixValue, suffixesWithThisVal);
			}
			suffixesWithThisVal.add(suffix);
			
			//Apply isotopes to suffixes if present
			if (suffix.getFrag() != null) {
				//boughton system applies to preceding suffix
				//iupac system applies to following suffix
				Element boughtonIsotopeSpecification = OpsinTools.getNextSibling(suffix);
				if (boughtonIsotopeSpecification != null && boughtonIsotopeSpecification.getName().equals(ISOTOPESPECIFICATION_EL)) {
					if (BOUGHTONSYSTEM_TYPE_VAL.equals(boughtonIsotopeSpecification.getAttributeValue(TYPE_ATR))) {
						applyIsotopeToSuffix(suffix.getFrag(), boughtonIsotopeSpecification, false);
					}
					else {
						throw new RuntimeException("Unexpected isotope specification after suffix");
					}
				}
				Element iupacIsotopeSpecification = OpsinTools.getPreviousSibling(suffix);
				while (iupacIsotopeSpecification != null && iupacIsotopeSpecification.getName().equals(ISOTOPESPECIFICATION_EL) &&
							IUPACSYSTEM_TYPE_VAL.equals(iupacIsotopeSpecification.getAttributeValue(TYPE_ATR))) {
					Element next = OpsinTools.getPreviousSibling(iupacIsotopeSpecification);
					applyIsotopeToSuffix(suffix.getFrag(), iupacIsotopeSpecification, true);
					iupacIsotopeSpecification = next;
				}
			}
		}
		
		boolean reDetectCycles = false;
		List<Fragment> fragsToMerge = new ArrayList<>();
		for (Entry<String, List<Element>> entry : suffixValToSuffixes.entrySet()) {
			String suffixValue = entry.getKey();
			List<Element> suffixesWithThisVal = entry.getValue();
			List<Atom> possibleAtomsToAttachSuffixTo = null;
			List<SuffixRule> rulesToApply = suffixRules.getSuffixRuleTags(suffixTypeToUse, suffixValue, subgroupType);
			for (int suffixIndex = 0; suffixIndex < suffixesWithThisVal.size(); suffixIndex++) {
				Element suffix = suffixesWithThisVal.get(suffixIndex);
				Fragment suffixFrag = null;
				for (SuffixRule suffixRule : rulesToApply) {
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
								reDetectCycles = true;
							}
							else{
								int bondOrderRequired = firstAtomInSuffix.getIncomingValency();
								Atom fragAtomToUse = getFragAtomToUse(frag, suffix, suffixTypeToUse);
								if (fragAtomToUse == null) {
									if (possibleAtomsToAttachSuffixTo == null) {
										int substitutionsRequired = suffixesWithThisVal.size();
										possibleAtomsToAttachSuffixTo = FragmentTools.findnAtomsForSubstitution(frag, atomList.get(0), substitutionsRequired, bondOrderRequired, true);
										if (possibleAtomsToAttachSuffixTo == null) {
											throw new StructureBuildingException("No suitable atom found to attach " + suffixValue + " suffix");
										}
										for (Atom atom : possibleAtomsToAttachSuffixTo) {
											if (FragmentTools.isCharacteristicAtom(atom)){
												throw new StructureBuildingException("No suitable atom found to attach suffix");
											}
										}
										if ("yes".equals(suffixRule.getAttributeValue(SUFFIXRULES_KETONELOCANT_ATR)) && !atomList.get(0).getAtomIsInACycle()) {
											List<Atom> proKetoneAtoms = getProKetonePositions(possibleAtomsToAttachSuffixTo);
											//Note that names like "ethanone" are allowable as the fragment may subsequently be substituted to form an actual ketone 
											if (proKetoneAtoms.size() >= substitutionsRequired) {
												possibleAtomsToAttachSuffixTo = proKetoneAtoms;
											}
										}
										if (!(substitutionsRequired == 1 && (ALKANESTEM_SUBTYPE_VAL.equals(frag.getSubType()) || HETEROSTEM_SUBTYPE_VAL.equals(frag.getSubType())) && possibleAtomsToAttachSuffixTo.get(0).equals(frag.getFirstAtom()))) {
											if (AmbiguityChecker.isSubstitutionAmbiguous(possibleAtomsToAttachSuffixTo, substitutionsRequired)) {
												state.addIsAmbiguous("Addition of " + suffixValue +" suffix to: " + group.getValue());
											}
										}
									}
									fragAtomToUse = possibleAtomsToAttachSuffixTo.get(suffixIndex);
								}
		
								//create a new bond and associate it with the suffixfrag and both atoms. Remember the suffixFrag has not been imported into the frag yet
								List<Bond> bonds = new ArrayList<>(firstAtomInSuffix.getBonds());
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
					fragsToMerge.add(suffixFrag);
					suffix.setFrag(null);
				}
			}
		}
		for (Fragment suffixFrag : fragsToMerge) {
			state.fragManager.removeAtomAndAssociatedBonds(suffixFrag.getFirstAtom());//the dummy R atom
			Set<String> suffixLocants = new HashSet<>(suffixFrag.getLocants());
			for (String suffixLocant : suffixLocants) {
				if (Character.isDigit(suffixLocant.charAt(0))){//check that numeric locants do not conflict with the parent fragment e.g. hydrazide 2' with biphenyl 2'
					if (frag.hasLocant(suffixLocant)){
						suffixFrag.getAtomByLocant(suffixLocant).removeLocant(suffixLocant);
					}
				}
			}
			state.fragManager.incorporateFragment(suffixFrag, frag);
		}
		if (reDetectCycles) {
			CycleDetector.assignWhetherAtomsAreInCycles(frag);
		}

	}
	
	private void applyIsotopeToSuffix(Fragment frag, Element isotopeSpecification, boolean mustBeApplied) throws StructureBuildingException {
		IsotopeSpecification isotopeSpec = IsotopeSpecificationParser.parseIsotopeSpecification(isotopeSpecification);
		ChemEl chemEl = isotopeSpec.getChemEl();
		int isotope = isotopeSpec.getIsotope();
		int multiplier = isotopeSpec.getMultiplier();
		String[] locants = isotopeSpec.getLocants();
		if (locants != null && !mustBeApplied) {
			//locanted boughton isotope probably applies to the group rather than the suffix
			return;
		}
		if (locants == null) {
			List<Atom> atoms = frag.getAtomList();
			atoms.remove(0);
			if (chemEl == ChemEl.H) {
				List<Atom> parentAtomsToApplyTo = FragmentTools.findnAtomsForSubstitution(atoms, null, multiplier, 1, true);
				if (parentAtomsToApplyTo == null) {
					if (mustBeApplied) {
						throw new StructureBuildingException("Failed to find sufficient hydrogen atoms for unlocanted hydrogen isotope replacement");
					}
					else {
						return;
					}
				}
				if (AmbiguityChecker.isSubstitutionAmbiguous(parentAtomsToApplyTo, multiplier)) {
					state.addIsAmbiguous("Position of hydrogen isotope on " + frag.getTokenEl().getValue());
				}
				for (int j = 0; j < multiplier; j++) {
					Atom atomWithHydrogenIsotope = parentAtomsToApplyTo.get(j);
					Atom hydrogen = state.fragManager.createAtom(isotopeSpec.getChemEl(), frag);
					hydrogen.setIsotope(isotope);
					state.fragManager.createBond(atomWithHydrogenIsotope, hydrogen, 1);
				}
			}
			else {
				List<Atom> parentAtomsToApplyTo = new ArrayList<>();
				for (Atom atom : atoms) {
					if (atom.getElement() == chemEl) {
						parentAtomsToApplyTo.add(atom);
					}
				}
				if (parentAtomsToApplyTo.size() < multiplier) {
					if(mustBeApplied) {
						throw new StructureBuildingException("Failed to find sufficient atoms for " + chemEl.toString() + " isotope replacement");
					}
					else {
						return;
					}
				}
				if (AmbiguityChecker.isSubstitutionAmbiguous(parentAtomsToApplyTo, multiplier)) {
					state.addIsAmbiguous("Position of isotope on " + frag.getTokenEl().getValue());
				}
				for (int j = 0; j < multiplier; j++) {
					parentAtomsToApplyTo.get(j).setIsotope(isotope);
				}
			}
		}
		else {
			if (chemEl == ChemEl.H) {
				for (int j = 0; j < locants.length; j++) {
					Atom atomWithHydrogenIsotope = frag.getAtomByLocantOrThrow(locants[j]);
					Atom hydrogen = state.fragManager.createAtom(isotopeSpec.getChemEl(), frag);
					hydrogen.setIsotope(isotope);
					state.fragManager.createBond(atomWithHydrogenIsotope, hydrogen, 1);
				}
			}
			else {
				for (int j = 0; j < locants.length; j++) {
					Atom atom = frag.getAtomByLocantOrThrow(locants[j]);
					if (chemEl != atom.getElement()) {
						throw new StructureBuildingException("The atom at locant: " + locants[j]  + " was not a " + chemEl.toString() );
					}
					atom.setIsotope(isotope);
				}
			}
		}
		isotopeSpecification.detach();
	}
	

	/**
	 * Return the subset of atoms that are "pro-ketone"
	 * i.e. a [CD2](C)C
	 * @param atoms
	 * @return
	 */
	private List<Atom> getProKetonePositions(List<Atom> atoms) {
		List<Atom> proKetonePositions = new ArrayList<>();
		for (Atom atom : atoms) {
			List<Bond> bonds = atom.getBonds();
			if (bonds.size() == 2 && 
					bonds.get(0).getOrder() == 1 &&
					bonds.get(1).getOrder() == 1 &&
					bonds.get(0).getOtherAtom(atom).getElement() == ChemEl.C &&
					bonds.get(1).getOtherAtom(atom).getElement() == ChemEl.C) {
				proKetonePositions.add(atom);
			}
		}
		return proKetonePositions;
	}
	
	private void processCycleFormingSuffix(Fragment suffixFrag, Fragment suffixableFragment, Element suffix) throws StructureBuildingException, ComponentGenerationException {
		List<Atom> rAtoms = new ArrayList<>();
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
			String[] locants = locant.split(",");
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
			String[] locantIds = locantId.split(",");
			if (locantIds.length !=2){
				throw new ComponentGenerationException("OPSIN bug: Should be exactly 2 locants associated with a cyclic suffix");
			}
			parentAtom1 = suffixableFragment.getAtomByIDOrThrow(Integer.parseInt(locantIds[0]));
			parentAtom2 = suffixableFragment.getAtomByIDOrThrow(Integer.parseInt(locantIds[1]));
		}
		else{
			int chainLength = suffixableFragment.getChainLength();
			if (chainLength > 1 && chainLength == suffixableFragment.getAtomCount()){
				parentAtom1 = suffixableFragment.getAtomByLocantOrThrow("1");
				parentAtom2 = suffixableFragment.getAtomByLocantOrThrow(String.valueOf(chainLength));
			}
			else{
				List<Atom> hydroxyAtoms = FragmentTools.findHydroxyGroups(suffixableFragment);
				if (hydroxyAtoms.size() == 1 && suffixableFragment.getAtomByLocant("1") != null){
					parentAtom1 = suffixableFragment.getAtomByLocantOrThrow("1");
					parentAtom2 = hydroxyAtoms.get(0);
				}
				else{
					throw new ComponentGenerationException("cycle forming suffix: " + suffix.getValue() +" should be locanted!");
				}
			}
		}
		if (parentAtom1.equals(parentAtom2)){
			throw new ComponentGenerationException("cycle forming suffix: " + suffix.getValue() +" attempted to form a cycle involving the same atom twice!");
		}
		
		if (suffixableFragment.getType().equals(CARBOHYDRATE_TYPE_VAL)){
			FragmentTools.removeTerminalOxygen(state, parentAtom1, 2);
			FragmentTools.removeTerminalOxygen(state, parentAtom1, 1);
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
	
	private Atom getFragAtomToUse(Fragment frag, Element suffix, String suffixTypeToUse) throws StructureBuildingException {
		String locant = suffix.getAttributeValue(LOCANT_ATR);
		if (locant != null) {
			return frag.getAtomByLocantOrThrow(locant);
		}
		String locantId = suffix.getAttributeValue(LOCANTID_ATR);
		if (locantId != null) {
			return frag.getAtomByIDOrThrow(Integer.parseInt(locantId));
		}
		String defaultLocantId = suffix.getAttributeValue(DEFAULTLOCANTID_ATR);
		if (defaultLocantId != null) {
			return frag.getAtomByIDOrThrow(Integer.parseInt(defaultLocantId));
		}
		else if (suffixTypeToUse.equals(ACIDSTEM_TYPE_VAL) || suffixTypeToUse.equals(NONCARBOXYLICACID_TYPE_VAL) || suffixTypeToUse.equals(CHALCOGENACIDSTEM_TYPE_VAL)) {//means that e.g. sulfonyl, has an explicit outAtom
			return frag.getFirstAtom();
		}
		return null;
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
		//List of atoms that can accept this charge while remaining in a reasonable valency
		List<Atom> nitrogens = new ArrayList<>();//most likely
		List<Atom> otherHeteroatoms = new ArrayList<>();//plausible
		List<Atom> carbonsAtoms = new ArrayList<>();//rare
		List<Atom> chargedAtoms = new ArrayList<>();//very rare
		if (atomList.isEmpty()) {
			throw new RuntimeException("OPSIN Bug: List of atoms to add charge suffix to was empty");
		}
		for (Atom a : atomList) {
			ChemEl chemEl = a.getElement();
			Integer[] stableValencies = ValencyChecker.getPossibleValencies(chemEl, a.getCharge() + chargeChange);
			if (stableValencies == null) {//unstable valency so seems unlikely
				continue;
			}
			int resultantExpectedValency = (a.getLambdaConventionValency() ==null ? ValencyChecker.getDefaultValency(chemEl) : a.getLambdaConventionValency()) + a.getProtonsExplicitlyAddedOrRemoved() + protonChange;
			
			if (!Arrays.asList(stableValencies).contains(resultantExpectedValency)) {
				//unstable valency so seems unlikely
				continue;
			}
			if (protonChange < 0) {
				int substitableHydrogen = StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(a);
				if (a.hasSpareValency() && !a.getFrag().getIndicatedHydrogen().contains(a)) {
					substitableHydrogen--;
				}
				if (substitableHydrogen < 1) {
					//no hydrogens so operation can't remove one!
					continue;
				}
			}
			if (a.getCharge() == 0) {
				if (chemEl == ChemEl.N) {
					nitrogens.add(a);
				}
				else if (chemEl != ChemEl.C) {
					otherHeteroatoms.add(a);
				}
				else {
					carbonsAtoms.add(a);
				}
			}
			else {
				chargedAtoms.add(a);
			}
		}
		List<Atom> listFromWhichToChoose;
		if (!nitrogens.isEmpty()) {
			listFromWhichToChoose = nitrogens;
			if (AMINOACID_TYPE_VAL.equals(atomList.get(0).getFrag().getType())) {
				//By convention treat names like lysinium as unambiguous (prefer alpha nitrogen)
				if (listFromWhichToChoose.contains(atomList.get(0))){
					listFromWhichToChoose = new ArrayList<>();
					listFromWhichToChoose.add(atomList.get(0));
				}
			}
		}
		else if (!otherHeteroatoms.isEmpty()) {
			listFromWhichToChoose = otherHeteroatoms;
		}
		else if (!carbonsAtoms.isEmpty()) {
			listFromWhichToChoose = carbonsAtoms;
		}
		else if (!chargedAtoms.isEmpty()) {
			listFromWhichToChoose = chargedAtoms;
		}
		else {
			listFromWhichToChoose = atomList;
		}

		Atom chosenAtom =  listFromWhichToChoose.get(0);
		if (!AmbiguityChecker.allAtomsEquivalent(listFromWhichToChoose)) {
			state.addIsAmbiguous("Addition of charge suffix to: " + chosenAtom.getFrag().getTokenEl().getValue());
		}

		chosenAtom.addChargeAndProtons(chargeChange, protonChange);
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
						List<String> tempLocants1 = new ArrayList<>(atom.getLocants());
						List<String> tempLocants2 = new ArrayList<>(atomToSwapWith.getLocants());
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
	 * Creates bonds between the parentAtom and the atoms connected to the R atoms.
	 * Removes bonds to the R atom
	 * @param parentAtom
	 * @param suffixRAtom
	 */
	private void makeBondsToSuffix(Atom parentAtom, Atom suffixRAtom) {
		List<Bond> bonds = new ArrayList<>(suffixRAtom.getBonds());
		for (Bond bondToSuffix : bonds) {
			Atom suffixAtom = bondToSuffix.getOtherAtom(suffixRAtom);
			state.fragManager.createBond(parentAtom, suffixAtom, bondToSuffix.getOrder());
			state.fragManager.removeBond(bondToSuffix);
		}
	}

	List<SuffixRule> getSuffixRuleTags(String suffixTypeToUse, String suffixValue, String subgroupType) throws ComponentGenerationException {
		return suffixRules.getSuffixRuleTags(suffixTypeToUse, suffixValue, subgroupType);
	}
}
