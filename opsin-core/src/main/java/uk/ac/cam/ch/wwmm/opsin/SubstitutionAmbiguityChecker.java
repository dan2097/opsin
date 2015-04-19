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

class SubstitutionAmbiguityChecker {

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
		StereoAnalyser analyzer = analyzeRelevantAtomsAndBonds(uniqueAtoms);
		Set<String> uniqueEnvironments = new HashSet<String>();
		for (Atom a : substitutableAtoms) {
			Integer env = analyzer.getAtomEnvironmentNumber(a);
			if (env == null){
				throw new RuntimeException("OPSIN Bug: Atom was not part of ambiguity analysis");
			}
			uniqueEnvironments.add(env + "\t" + a.getOutValency());
		}
		if (uniqueEnvironments.size() == 1 && (numberToBeSubstituted == 1 || numberToBeSubstituted == substitutableAtoms.size() - 1)){
			return false;
		}
		return true;
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
