package uk.ac.cam.ch.wwmm.opsin;

import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.ELEMENTARYATOM_SUBTYPE_VAL;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

class AmbiguityChecker {

	static boolean isSubstitutionAmbiguous(List<Atom> substitutableAtoms, int numberToBeSubstituted) {
		if (substitutableAtoms.size() == 0) {
			throw new IllegalArgumentException("OPSIN Bug: Must provide at least one substituable atom");
		}
		if (substitutableAtoms.size() < numberToBeSubstituted) {
			throw new IllegalArgumentException("OPSIN Bug: substitutableAtoms must be >= numberToBeSubstituted");
		}
		if (substitutableAtoms.size() == numberToBeSubstituted){
			return false;
		}
		if (allAtomsConnectToDefaultInAtom(substitutableAtoms, numberToBeSubstituted)) {
			return false;
		}
		Set<Atom> uniqueAtoms = new HashSet<Atom>(substitutableAtoms);
		if (uniqueAtoms.size() == 1) {
			return false;
		}
		if (allAtomsEquivalent(uniqueAtoms) && (numberToBeSubstituted == 1 || numberToBeSubstituted == substitutableAtoms.size() - 1)){
			return false;
		}
		return true;
	}
	
	static boolean allAtomsEquivalent(Collection<Atom> atoms) {
		StereoAnalyser analyser = analyseRelevantAtomsAndBonds(atoms);
		Set<String> uniqueEnvironments = new HashSet<String>();
		for (Atom a : atoms) {
			uniqueEnvironments.add(getAtomEnviron(analyser, a));
		}
		return uniqueEnvironments.size() == 1;
	}

	static boolean allBondsEquivalent(Collection<Bond> bonds) {
		Set<Atom> relevantAtoms = new HashSet<Atom>();
		for (Bond b : bonds) {
			relevantAtoms.add(b.getFromAtom());
			relevantAtoms.add(b.getToAtom());
		}
		StereoAnalyser analyser = analyseRelevantAtomsAndBonds(relevantAtoms);
		Set<String> uniqueBonds = new HashSet<String>();
		for (Bond b : bonds) {
			uniqueBonds.add(bondToCanonicalEnvironString(analyser, b));
		}
		return uniqueBonds.size() == 1;
	}

	private static String bondToCanonicalEnvironString(StereoAnalyser analyser, Bond b) {
		String s1 = getAtomEnviron(analyser, b.getFromAtom());
		String s2 = getAtomEnviron(analyser, b.getToAtom());
		if (s1.compareTo(s2) > 0){
			return s1 + s2;
		}
		else {
			return s2 + s1;
		}
	}

	static String getAtomEnviron(StereoAnalyser analyser, Atom a) {
		Integer env = analyser.getAtomEnvironmentNumber(a);
		if (env == null) {
			throw new RuntimeException("OPSIN Bug: Atom was not part of ambiguity analysis");
		}
		//"identical" atoms may be distinguished by bonds yet to be formed, hence split by outvalency
		// e.g. [PH3] vs [PH3]=
		return env + "\t" + a.getOutValency();
	}

	private static boolean allAtomsConnectToDefaultInAtom(List<Atom> substitutableAtoms, int numberToBeSubstituted) {
		Atom defaultInAtom = substitutableAtoms.get(0).getFrag().getDefaultInAtom();
		if (defaultInAtom != null) {
			for (int i = 0; i < numberToBeSubstituted; i++) {
				if (!substitutableAtoms.get(i).equals(defaultInAtom)) {
					return false;
				}
			}
			return true;
		}
		return false;
	}

	static StereoAnalyser analyseRelevantAtomsAndBonds(Collection<Atom> startingAtoms) {
		Set<Atom> atoms = new HashSet<Atom>();
		Set<Bond> bonds = new HashSet<Bond>();
		Deque<Atom> stack = new ArrayDeque<Atom>(startingAtoms);
		while (!stack.isEmpty()) {
			Atom a = stack.removeLast();
			if (!atoms.contains(a)) {
				atoms.add(a);
				for (Bond b : a.getBonds()) {
					bonds.add(b);
					stack.add(b.getOtherAtom(a));
				}
			}
		}
		
		List<Atom> ghostHydrogens = new ArrayList<Atom>();
		for (Atom atom : atoms) {
			if (atom.getFrag().getSubType().equals(ELEMENTARYATOM_SUBTYPE_VAL)){//these do not have implicit hydrogen e.g. phosphorus is literally just a phosphorus atom
				continue;
			}
			int explicitHydrogensToAdd = StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(atom);
			for (int i = 0; i < explicitHydrogensToAdd; i++) {
				Atom ghostHydrogen = new Atom(ChemEl.H);
				Bond b = new Bond(ghostHydrogen, atom, 1);
				atom.addBond(b);
				ghostHydrogen.addBond(b);
				ghostHydrogens.add(ghostHydrogen);
			}
		}
		atoms.addAll(ghostHydrogens);
		StereoAnalyser analyzer = new StereoAnalyser(atoms, bonds);
		for (Atom ghostHydrogen : ghostHydrogens) {
			Bond b = ghostHydrogen.getFirstBond();
			b.getOtherAtom(ghostHydrogen).removeBond(b);
		}
		return analyzer;
	}

	static List<Atom> useAtomEnvironmentsToGivePlausibleSubstitution(List<Atom> substitutableAtoms, int numberToBeSubstituted) {
		if (substitutableAtoms.size() == 0) {
			throw new IllegalArgumentException("OPSIN Bug: Must provide at least one substituable atom");
		}
		if (substitutableAtoms.size() < numberToBeSubstituted) {
			throw new IllegalArgumentException("OPSIN Bug: substitutableAtoms must be >= numberToBeSubstituted");
		}
		if (substitutableAtoms.size() == numberToBeSubstituted){
			return substitutableAtoms;
		}

		List<Atom> preferredAtoms = findPlausibleSubstitutionPatternUsingSymmmetry(substitutableAtoms, numberToBeSubstituted);
		if (preferredAtoms != null){
			return preferredAtoms;
		}
		return findPlausibleSubstitutionPatternUsingLocalEnvironment(substitutableAtoms, numberToBeSubstituted);
	}

	private static List<Atom> findPlausibleSubstitutionPatternUsingSymmmetry(List<Atom> substitutableAtoms, int numberToBeSubstituted) {
		//cf. octaethylporphyrin (8 identical atoms capable of substitution)
		StereoAnalyser analyser = analyseRelevantAtomsAndBonds(new HashSet<Atom>(substitutableAtoms));
		Map<String, List<Atom>> atomsInEachEnvironment = new HashMap<String, List<Atom>>();
		for (Atom a : substitutableAtoms) {
			String env = getAtomEnviron(analyser, a);
			List<Atom> atomsInEnvironment = atomsInEachEnvironment.get(env);
			if (atomsInEnvironment == null) {
				atomsInEnvironment = new ArrayList<Atom>();
				atomsInEachEnvironment.put(env, atomsInEnvironment);
			}
			atomsInEnvironment.add(a);
		}
		List<Atom> preferredAtoms = null;
		for (List<Atom> atoms : atomsInEachEnvironment.values()) {
			if (atoms.size() == numberToBeSubstituted){
				if (preferredAtoms != null){
					return null;
				}
				preferredAtoms = atoms;
			}
		}
		if (preferredAtoms == null) {
			//check for environments with double the required atoms where this means each atom can support two substitutions c.f. cyclohexane
			for (List<Atom> atoms : atomsInEachEnvironment.values()) {
				if (atoms.size() == (numberToBeSubstituted * 2)){
					Set<Atom> uniquified = new LinkedHashSet<Atom>(atoms);//retain deterministic atom ordering
					if (uniquified.size() == numberToBeSubstituted) {
						if (preferredAtoms != null){
							return null;
						}
						preferredAtoms = new ArrayList<Atom>(uniquified);
					}
				}
			}
		}
		return preferredAtoms;
	}
	
	private static List<Atom> findPlausibleSubstitutionPatternUsingLocalEnvironment(List<Atom> substitutableAtoms, int numberToBeSubstituted) {
		//cf. pentachlorotoluene (5 sp2 carbons vs sp3 methyl)
		Map<String, List<Atom>> atomsInEachLocalEnvironment = new HashMap<String, List<Atom>>();
		for (Atom a : substitutableAtoms) {
			int valency = a.determineValency(true);	
			int currentValency = a.getIncomingValency() + a.getOutValency();
			int numOfBonds = (valency - currentValency) + a.getBondCount();//distinguish sp2 and sp3 atoms
			String s = a.getElement().toString() +"\t" + valency + "\t" + numOfBonds + "\t" + a.hasSpareValency();
			List<Atom> atomsInEnvironment = atomsInEachLocalEnvironment.get(s);
			if (atomsInEnvironment == null) {
				atomsInEnvironment = new ArrayList<Atom>();
				atomsInEachLocalEnvironment.put(s, atomsInEnvironment);
			}
			atomsInEnvironment.add(a);
		}
		List<Atom> preferredAtoms = null;
		for (List<Atom> atoms : atomsInEachLocalEnvironment.values()) {
			if (atoms.size() == numberToBeSubstituted){
				if (preferredAtoms != null){
					return null;
				}
				preferredAtoms = atoms;
			}
		}
		return preferredAtoms;
	}
}
