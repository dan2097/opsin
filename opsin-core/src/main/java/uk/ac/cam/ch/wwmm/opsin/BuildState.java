package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import uk.ac.cam.ch.wwmm.opsin.OpsinWarning.OpsinWarningType;

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
	private final List<OpsinWarning> warnings = new ArrayList<OpsinWarning>();
	
	WordRule currentWordRule = null;

	BuildState(NameToStructureConfig n2sConfig) {
		this.n2sConfig = n2sConfig;
		IDManager idManager = new IDManager();
		fragManager = new FragmentManager(new SMILESFragmentBuilder(idManager), idManager);
		xmlSuffixMap = new HashMap<Element, List<Fragment>>();
	}

	List<OpsinWarning> getWarnings() {
		return warnings;
	}
	
	void addWarning(OpsinWarningType type, String message) {
		warnings.add(new OpsinWarning(type, message));
	}
	
	void addIsAmbiguous() {
		warnings.add(new OpsinWarning(OpsinWarningType.APPEARS_AMBIGUOUS, ""));
	}
	
	void addIsAmbiguous(String message) {
		warnings.add(new OpsinWarning(OpsinWarningType.APPEARS_AMBIGUOUS, message));
	}

	boolean checkForAmbiguity(List<Atom> substituentPoints, int numberOfSubstitutionsRequired) {
		boolean isAmbiguous = AmbiguityChecker.isSubstitutionAmbiguous(substituentPoints, numberOfSubstitutionsRequired);
		if (isAmbiguous){
			addIsAmbiguous();
		}
		return isAmbiguous;
	}


}
