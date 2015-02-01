package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * Used to pass the current configuration and FragmentManager around
 * The currentWordRule can be mutated to keep track of what the parent wordRule is at the given time
 *
 * @author dl387
 *
 */
class BuildState {

	final FragmentManager fragManager;
	final HashMap<Element, List<Fragment>> xmlSuffixMap;
	final NameToStructureConfig n2sConfig;
	private final List<String> warningMessages = new ArrayList<String>();
	
	WordRule currentWordRule = null;
	
	
	BuildState(NameToStructureConfig n2sConfig) {
		this.n2sConfig = n2sConfig;
		IDManager idManager = new IDManager();
		fragManager = new FragmentManager(new SMILESFragmentBuilder(idManager), idManager);
		xmlSuffixMap = new HashMap<Element, List<Fragment>>();
	}

	List<String> getWarningMessages() {
		return warningMessages;
	}
	
	String getWarningMessage() {
		int numWarnings = warningMessages.size();
		if (numWarnings == 0) {
			return null;
		}
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < numWarnings; i++) {
			sb.append(warningMessages.get(i));
			if (i + 1 > numWarnings){
				sb.append(System.getProperty("line.separator"));
			}
		}
		return sb.toString();
	}

	void addWarningMessage(String warningMessage) {
		if (warningMessage == null) {
			warningMessage = "null";
		}
		warningMessages.add(warningMessage);
	}

	void checkForAmbiguity(List<Atom> substituentPoints, int numberOfSubstitutionsRequired) {
		boolean isAmbiguous = SubstitutionAmbiguityChecker.isSubstitutionAmbiguous(substituentPoints, numberOfSubstitutionsRequired);
		if (isAmbiguous){
			addWarningMessage("Name appears to be ambiguous");
		}
	}
}
