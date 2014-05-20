package uk.ac.cam.ch.wwmm.opsin;

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
	
	WordRule currentWordRule = null;
	
	private String warningMessage = null;
	
	String getWarningMessage() {
		return warningMessage;
	}

	void addWarningMessage(String warningMessage) {
		if (warningMessage == null){
			this.warningMessage = warningMessage;
		}
		else{
			this.warningMessage += ("\n" + warningMessage);
		}
	}

	BuildState(NameToStructureConfig n2sConfig) {
		this.n2sConfig = n2sConfig;
		IDManager idManager = new IDManager();
		fragManager = new FragmentManager(new SMILESFragmentBuilder(idManager), idManager);
		xmlSuffixMap = new HashMap<Element, List<Fragment>>();
	}
}
