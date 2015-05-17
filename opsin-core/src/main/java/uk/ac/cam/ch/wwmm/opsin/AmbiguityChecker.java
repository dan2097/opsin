package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
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
		StereoAnalyser analyser = analyzeRelevantAtomsAndBonds(atoms);
		Set<String> uniqueEnvironments = new HashSet<String>();
		for (Atom a : atoms) {
			Integer env = analyser.getAtomEnvironmentNumber(a);
			if (env == null){
				throw new RuntimeException("OPSIN Bug: Atom was not part of ambiguity analysis");
			}
			uniqueEnvironments.add(env + "\t" + a.getOutValency());
		}
		return uniqueEnvironments.size() == 1;
	}

	static boolean allBondsEquivalent(Collection<Bond> bonds) {
		Set<Atom> relevantAtoms = new HashSet<Atom>();
		for (Bond b : bonds) {
			relevantAtoms.add(b.getFromAtom());
			relevantAtoms.add(b.getToAtom());
		}
		StereoAnalyser analyser = analyzeRelevantAtomsAndBonds(relevantAtoms);
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

	private static String getAtomEnviron(StereoAnalyser analyser, Atom a) {
		Integer env = analyser.getAtomEnvironmentNumber(a);
		if (env == null){
			throw new RuntimeException("OPSIN Bug: Atom was not part of ambiguity analysis");
		}
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

	static StereoAnalyser analyzeRelevantAtomsAndBonds(Collection<Atom> startingAtoms) {
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
		return new StereoAnalyser(atoms, bonds);
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
		StereoAnalyser analyzer = analyzeRelevantAtomsAndBonds(new HashSet<Atom>(substitutableAtoms));
		Map<String, List<Atom>> atomsInEachEnvironment = new HashMap<String, List<Atom>>();
		for (Atom a : substitutableAtoms) {
			Integer env = analyzer.getAtomEnvironmentNumber(a);
			if (env == null) {
				throw new RuntimeException("OPSIN Bug: Atom was not part of ambiguity analysis");
			}
			String s = env + "\t" + a.getOutValency();
			List<Atom> atomsInEnvironment = atomsInEachEnvironment.get(s);
			if (atomsInEnvironment == null) {
				atomsInEnvironment = new ArrayList<Atom>();
				atomsInEachEnvironment.put(s, atomsInEnvironment);
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
		return preferredAtoms;
	}
	
	private static List<Atom> findPlausibleSubstitutionPatternUsingLocalEnvironment(List<Atom> substitutableAtoms, int numberToBeSubstituted) {
		Map<String, List<Atom>> atomsInEachLocalEnvironment = new HashMap<String, List<Atom>>();
		for (Atom a : substitutableAtoms) {
			String s = a.getElement().toString() +"\t" + a.determineValency(true) +"\t" + a.hasSpareValency();
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
