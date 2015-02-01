package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

class SubstitutionAmbiguityChecker {

	public static boolean isSubstitutionAmbiguous(List<Atom> substitutableAtoms, int numberToBeSubstituted) {
		if (substitutableAtoms.size() == 0) {
			throw new IllegalArgumentException("Must provide at least one substituable atom");
		}
		if (substitutableAtoms.size() == numberToBeSubstituted){
			return false;
		}
		if (numberToBeSubstituted ==1 && substitutableAtoms.get(0).equals(substitutableAtoms.get(0).getFrag().getDefaultInAtom())) {
			return false;
		}
		Set<Atom> uniqueAtoms = new HashSet<Atom>(substitutableAtoms);
		if (uniqueAtoms.size() == 1) {
			return false;
		}
		
		StereoAnalyser analyzer = analyzeRelevantAtomsAndBonds(uniqueAtoms);
		Set<Integer> uniqueEnvironments = new HashSet<Integer>();
		for (Atom a : substitutableAtoms) {
			Integer env = analyzer.getAtomEnvironmentNumber(a);
			if (env == null){
				throw new RuntimeException("OPSIN Bug: Atom was not part of ambiguity analysis");
			}
			uniqueEnvironments.add(env);
		}
		if (uniqueEnvironments.size() == 1 && (numberToBeSubstituted == 1 || numberToBeSubstituted == substitutableAtoms.size() - 1)){
			return false;
		}
		return true;
	}

	private static StereoAnalyser analyzeRelevantAtomsAndBonds(Set<Atom> startingAtoms) {
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
}
