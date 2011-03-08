package uk.ac.cam.ch.wwmm.opsin;

/**
 * Allows OPSIN to be configured e.g. enable processing of radicals
 * Example usage:
 * NameToStructureConfig n2sConfig = new NameToStructureConfig();
 * n2sconfig.setAllowRadicals(true);
 * nts.parseChemicalName(chemicalName, n2sConfig) 
 * where nts is an instance of NameToStructure
 * @author dl387
 *
 */
public class NameToStructureConfig implements Cloneable {
	
	// Fields set with default values
	private boolean allowRadicals = false;
	private boolean detailedFailureAnalysis = false;
//	private boolean slackSpaceHandling;
//	private boolean substituentAbbreviations;
//	private boolean ignoreStereochemistry;
//	private boolean ignoreCurrentlyUninterpretableStereochemistry;
	

	/**
	 * Constructs a NameToStructureConfig with default settings:
	 * allowRadicals = false
	 * detailedFailureAnalysis = false
	 */
	public NameToStructureConfig() {
	}


	/**
	 * Are radicals allowed?  e.g. should fragments such as phenyl be interpretable
	 * @return
	 */
	public boolean isAllowRadicals() {
		return allowRadicals;
	}

	/**
	 * Sets whether radicals allowed? e.g. should fragments such as phenyl be interpretable
	 * @return
	 */
	public void setAllowRadicals(boolean allowRadicals) {
		this.allowRadicals = allowRadicals;
	}

	/**
	 * Should OPSIN attempt reverse parsing to more accurately determine why parsing failed
	 * @return
	 */
	public boolean isDetailedFailureAnalysis() {
		return detailedFailureAnalysis;
	}

	/**
	 * Sets whether OPSIN should attempt reverse parsing to more accurately determine why parsing failed
	 * @return
	 */
	public void setDetailedFailureAnalysis(boolean detailedFailureAnalysis) {
		this.detailedFailureAnalysis = detailedFailureAnalysis;
	}

//	boolean isSlackSpaceHandling() {
//		return slackSpaceHandling;
//	}
//
//	void setSlackSpaceHandling(boolean slackSpaceHandling) {
//		this.slackSpaceHandling = slackSpaceHandling;
//	}

//	boolean isSubstituentAbbreviations() {
//		return substituentAbbreviations;
//	}
//
//	void setSubstituentAbbreviations(boolean substituentAbbreviations) {
//		this.substituentAbbreviations = substituentAbbreviations;
//	}
//
//	boolean isIgnoreStereochemistry() {
//		return ignoreStereochemistry;
//	}
//
//	void setIgnoreStereochemistry(boolean ignoreStereochemistry) {
//		this.ignoreStereochemistry = ignoreStereochemistry;
//	}
//
//	boolean isIgnoreCurrentlyUninterpretableStereochemistry() {
//		return ignoreCurrentlyUninterpretableStereochemistry;
//	}
//
//	void setIgnoreCurrentlyUninterpretableStereochemistry(boolean ignoreCurrentlyUninterpretableStereochemistry) {
//		this.ignoreCurrentlyUninterpretableStereochemistry = ignoreCurrentlyUninterpretableStereochemistry;
//	}

	/**
	 * Constructs a NameToStructureConfig with default settings:
	 * allowRadicals = false
	 * detailedFailureAnalysis = false
	 */
	public static NameToStructureConfig getDefaultConfigInstance() {
		return new NameToStructureConfig();
	}
	
//	/**
//	 * Returns a NameToStructureConfig with the following settings:
//	 * allowRadicals = false
//	 * slackSpaceHandling = false
//	 * substituentAbbreviations = false
//	 * ignoreStereochemistry = false
//	 * ignoreCurrentlyUninterpretableStereochemistry = false
//	 */
//	public static NameToStructureConfig getStrictConfigInstance() {
//		NameToStructureConfig n2sConfig = new NameToStructureConfig();
//		n2sConfig.allowRadicals = false;
//		n2sConfig.slackSpaceHandling = true;
//		n2sConfig.substituentAbbreviations = false;
//		n2sConfig.ignoreStereochemistry = false;
//		n2sConfig.ignoreCurrentlyUninterpretableStereochemistry = true;
//		return n2sConfig;
//	}
	
	
	@Override
	public NameToStructureConfig clone() {
		try {
			NameToStructureConfig copy = (NameToStructureConfig) super.clone();
			return copy;
		} catch (CloneNotSupportedException e) {
			// Can only be thrown if we *don't* implement Cloneable, which we do...
			throw new Error("Impossible!", e);
		}
	}
	
}
