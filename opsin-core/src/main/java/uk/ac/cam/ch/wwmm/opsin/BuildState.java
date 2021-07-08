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
	// counter is used for DL- racemic stereochemistry in oligomers, we place each one in a separate racemic group,
	// there is implicitly one group in-case the input has a combination of (RS)- and then DL-
	int numRacGrps = 1;
	private final List<OpsinWarning> warnings = new ArrayList<>();
	
	WordRule currentWordRule = null;

	BuildState(NameToStructureConfig n2sConfig) {
		this.n2sConfig = n2sConfig;
		IDManager idManager = new IDManager();
		fragManager = new FragmentManager(new SMILESFragmentBuilder(idManager), idManager);
		xmlSuffixMap = new HashMap<>();
	}

	List<OpsinWarning> getWarnings() {
		return warnings;
	}
	
	void addWarning(OpsinWarningType type, String message) {
		warnings.add(new OpsinWarning(type, message));
	}
	
	void addIsAmbiguous(String message) {
		warnings.add(new OpsinWarning(OpsinWarningType.APPEARS_AMBIGUOUS, message));
	}
}
