package uk.ac.cam.ch.wwmm.opsin;

import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Replaces the superatoms of a phane simplified skeleton with amplificants (Phane Nomenclature Part I, PhI-2.3; P-26.3).
 * e.g. 1,4(1,4)-dibenzenacyclohexaphane
 */
class PhaneAmplifier {

	private static final Pattern MATCH_PLAIN_INTEGER = Pattern.compile("(?<![\\d^])\\d+(?![\\d^])");
	private static final Pattern MATCH_AMPLIFICATION = Pattern.compile("([1-9][0-9]*(?:,[1-9][0-9]*)*)\\(([^()]+)\\)");

	private static class Amplification {
		final int superatom;
		final String[] attachmentLocants;
		Fragment amplificant;

		Amplification(int superatom, String[] attachmentLocants) {
			this.superatom = superatom;
			this.attachmentLocants = attachmentLocants;
		}
	}

	static void processPhaneAmplification(BuildState state, Element subOrRoot) throws StructureBuildingException {
		Element skeletonGroup = null;
		for (Element group : subOrRoot.getChildElements(GROUP_EL)) {
			if (PHANE_SUBTYPE_VAL.equals(group.getAttributeValue(SUBTYPE_ATR))) {
				skeletonGroup = group;
			}
		}
		if (skeletonGroup == null) {
			return;
		}
		Fragment skeleton = skeletonGroup.getFrag();
		int nodeCount = skeleton.getAtomCount();

		List<Amplification> amplifications = new ArrayList<>();
		List<Element> consumed = new ArrayList<>();
		List<Amplification> pending = null;
		Integer pendingMultiplier = null;
		for (Element el : new ArrayList<>(subOrRoot.getChildElements())) {
			if (el == skeletonGroup) {
				break;
			}
			switch (el.getName()) {
			case PHANELOCANTS_EL:
				pending = parseAmplifications(el.getValue(), nodeCount);
				consumed.add(el);
				break;
			case MULTIPLIER_EL:
				pendingMultiplier = Integer.parseInt(el.getAttributeValue(VALUE_ATR));
				consumed.add(el);
				break;
			case GROUP_EL:
				if (pending == null) {
					throw new StructureBuildingException("Amplificant " + el.getValue() + " has no locants indicating the superatoms it replaces");
				}
				if (pendingMultiplier != null && pendingMultiplier != pending.size() || pendingMultiplier == null && pending.size() != 1) {
					throw new StructureBuildingException("Disagreement between the multiplier and the number of superatoms of the amplificant");
				}
				Fragment template = el.getFrag();
				for (int i = 0; i < pending.size(); i++) {
					pending.get(i).amplificant = i == 0 ? template : state.fragManager.copyFragment(template);
				}
				amplifications.addAll(pending);
				consumed.add(el);
				pending = null;
				pendingMultiplier = null;
				break;
			case STRUCTURALOPENBRACKET_EL:
			case STRUCTURALCLOSEBRACKET_EL:
				consumed.add(el);
				break;
			default:
				break;
			}
		}
		if (pending != null) {
			throw new StructureBuildingException("Amplification locants are not followed by an amplificant");
		}
		if (amplifications.isEmpty()) {
			return;
		}

		Atom[] nodes = new Atom[nodeCount + 1];
		for (int i = 1; i <= nodeCount; i++) {
			nodes[i] = skeleton.getAtomByLocantOrThrow(Integer.toString(i));
		}
		List<List<Integer>> neighbours = new ArrayList<>();
		neighbours.add(Collections.<Integer>emptyList());
		for (int i = 1; i <= nodeCount; i++) {
			List<Integer> locants = new ArrayList<>();
			for (Atom neighbour : nodes[i].getAtomNeighbours()) {
				locants.add(Integer.parseInt(neighbour.getFirstLocant()));
			}
			Collections.sort(locants);
			neighbours.add(locants);
		}

		Amplification[] bySuperatom = new Amplification[nodeCount + 1];
		for (Amplification amplification : amplifications) {
			if (bySuperatom[amplification.superatom] != null) {
				throw new StructureBuildingException("Skeleton node " + amplification.superatom + " is amplified more than once");
			}
			bySuperatom[amplification.superatom] = amplification;
			if (amplification.attachmentLocants.length != neighbours.get(amplification.superatom).size()) {
				throw new StructureBuildingException("Superatom " + amplification.superatom + " is bonded to " + neighbours.get(amplification.superatom).size() +
						" skeleton nodes but " + amplification.attachmentLocants.length + " attachment locants were given");
			}
		}

		for (Amplification amplification : amplifications) {
			nodes[amplification.superatom].clearLocants();
		}
		Map<String, List<Atom>> concatenatedLocantOwners = new LinkedHashMap<>();
		Atom[][] attachmentAtoms = new Atom[nodeCount + 1][];
		for (Amplification amplification : amplifications) {
			Fragment amplificant = amplification.amplificant;
			Atom[] attachments = new Atom[amplification.attachmentLocants.length];
			for (int i = 0; i < attachments.length; i++) {
				attachments[i] = amplificant.getAtomByLocantOrThrow(amplification.attachmentLocants[i]);
			}
			attachmentAtoms[amplification.superatom] = attachments;
			for (Atom attachment : attachments) {
				int skeletonBonds = 0;
				for (Atom other : attachments) {
					if (other == attachment) {
						skeletonBonds++;
					}
				}
				//a mancude carbon needs one valency for its double bond so can form at most three sigma bonds
				if (attachment.getElement() == ChemEl.C && attachment.hasSpareValency() && attachment.getIncomingValency() + skeletonBonds > 3) {
					throw new StructureBuildingException("The mancude amplificant atom " + attachment.getFirstLocant() + " cannot form the requested number of bonds to the simplified skeleton");
				}
			}
			for (Atom atom : amplificant) {
				List<String> locants = new ArrayList<>(atom.getLocants());
				atom.clearLocants();
				for (String locant : locants) {
					atom.addLocant(amplification.superatom + "^" + locant);
					if (Character.isDigit(locant.charAt(0))) {
						String concatenated = amplification.superatom + locant;
						List<Atom> owners = concatenatedLocantOwners.get(concatenated);
						if (owners == null) {
							owners = new ArrayList<>();
							concatenatedLocantOwners.put(concatenated, owners);
						}
						owners.add(atom);
					}
				}
			}
		}
		Set<String> ambiguousLocants = new HashSet<>();
		for (Map.Entry<String, List<Atom>> entry : concatenatedLocantOwners.entrySet()) {
			boolean skeletonAtomHasLocant = isRemainingSkeletonLocant(entry.getKey(), nodeCount, bySuperatom);
			if (entry.getValue().size() == 1 && !skeletonAtomHasLocant) {
				entry.getValue().get(0).addLocant(entry.getKey());
			}
			else if (skeletonAtomHasLocant || entry.getValue().size() > 1) {
				ambiguousLocants.add(entry.getKey());
			}
		}
		if (!ambiguousLocants.isEmpty()) {
			warnIfAmbiguousLocantsAreUsed(state, subOrRoot, ambiguousLocants);
		}
		for (Amplification amplification : amplifications) {
			state.fragManager.incorporateFragment(amplification.amplificant, skeleton);
		}

		for (int i = 1; i <= nodeCount; i++) {
			for (int j : neighbours.get(i)) {
				if (j < i || bySuperatom[i] == null && bySuperatom[j] == null) {
					continue;
				}
				Atom from = getBondingAtom(i, j, nodes, neighbours, attachmentAtoms);
				Atom to = getBondingAtom(j, i, nodes, neighbours, attachmentAtoms);
				state.fragManager.createBond(from, to, 1);
			}
		}
		for (Amplification amplification : amplifications) {
			state.fragManager.removeAtomAndAssociatedBonds(nodes[amplification.superatom]);
		}
		for (Element el : consumed) {
			el.detach();
		}
	}

	/**
	 * A superscript written inline (e.g. 14 for 1^4) cannot be distinguished from a skeleton locant when both exist
	 */
	private static void warnIfAmbiguousLocantsAreUsed(BuildState state, Element subOrRoot, Set<String> ambiguousLocants) {
		Element word = OpsinTools.getParentWordRule(subOrRoot);
		Element scope = word != null ? word : subOrRoot;
		Set<String> used = new HashSet<>();
		for (Element el : OpsinTools.getDescendantElementsWithTagName(scope, LOCANT_EL)) {
			if (el.getValue().indexOf('^') >= 0) {
				return;//superscripts are written explicitly so plain integers are skeleton locants
			}
			collectPlainIntegers(el.getValue(), used);
		}
		for (Element el : OpsinTools.getDescendantElementsWithTagNames(scope, new String[]{INDICATEDHYDROGEN_EL, ADDEDHYDROGEN_EL, HYDRO_EL, HETEROATOM_EL, UNSATURATOR_EL, SUFFIX_EL, STEREOCHEMISTRY_EL})) {
			String locant = el.getAttributeValue(LOCANT_ATR);
			if (locant != null) {
				if (locant.indexOf('^') >= 0) {
					return;
				}
				collectPlainIntegers(locant, used);
			}
		}
		for (String locant : ambiguousLocants) {
			if (used.contains(locant)) {
				state.addIsAmbiguous("The locant " + locant + " could refer to an atom of the simplified skeleton or to an atom of an amplificant (superscript locants should be written e.g. 1^4)");
			}
		}
	}

	private static void collectPlainIntegers(String text, Set<String> used) {
		Matcher m = MATCH_PLAIN_INTEGER.matcher(text);
		while (m.find()) {
			used.add(m.group());
		}
	}

	/**
	 * Attachment locants are cited in the order of the locants of the skeleton nodes they are bonded to (lowest first)
	 */
	private static Atom getBondingAtom(int node, int neighbour, Atom[] nodes, List<List<Integer>> neighbours, Atom[][] attachmentAtoms) {
		if (attachmentAtoms[node] == null) {
			return nodes[node];
		}
		return attachmentAtoms[node][neighbours.get(node).indexOf(neighbour)];
	}

	private static boolean isRemainingSkeletonLocant(String locant, int nodeCount, Amplification[] bySuperatom) {
		if (!locant.chars().allMatch(Character::isDigit) || locant.length() > 4) {
			return false;
		}
		int number = Integer.parseInt(locant);
		return number <= nodeCount && bySuperatom[number] == null;
	}

	private static List<Amplification> parseAmplifications(String text, int nodeCount) throws StructureBuildingException {
		List<Amplification> amplifications = new ArrayList<>();
		Matcher m = MATCH_AMPLIFICATION.matcher(StringTools.removeDashIfPresent(text));
		while (m.find()) {
			String[] attachments = m.group(2).split(",");
			for (String locant : m.group(1).split(",")) {
				int superatom = Integer.parseInt(locant);
				if (superatom > nodeCount) {
					throw new StructureBuildingException("Superatom locant " + superatom + " exceeds the number of skeleton nodes: " + nodeCount);
				}
				amplifications.add(new Amplification(superatom, attachments));
			}
		}
		return amplifications;
	}
}
