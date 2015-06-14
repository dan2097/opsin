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
	private boolean outputRadicalsAsWildCardAtoms = false;
	private boolean detailedFailureAnalysis = false;
	private boolean interpretAcidsWithoutTheWordAcid = false;
	private boolean warnRatherThanFailOnUninterpretableStereochemistry = false;

	/**
	 * Constructs a NameToStructureConfig with default settings:
	 * allowRadicals = false
	 * outputRadicalsAsWildCardAtoms = false
	 * detailedFailureAnalysis = false
	 * interpretAcidsWithoutTheWordAcid = false
	 * warnRatherThanFailOnUninterpretableStereochemistry = false
	 */
	public NameToStructureConfig() {
	}


	/**
	 * Are radicals allowed?  e.g. should fragments such as phenyl be interpretable
	 * @return whether radicals are allowed
	 */
	public boolean isAllowRadicals() {
		return allowRadicals;
	}

	/**
	 * Sets whether radicals allowed? e.g. should fragments such as phenyl be interpretable
	 */
	public void setAllowRadicals(boolean allowRadicals) {
		this.allowRadicals = allowRadicals;
	}
	
	/**
	 * Are radicals output as wildcard atoms e.g. [*]CC for ethyl
	 * @return whether radicals are output using explicit wildcard atoms
	 */
	public boolean isOutputRadicalsAsWildCardAtoms() {
		return outputRadicalsAsWildCardAtoms;
	}

	/**
	 * Should radicals be output as wildcard atoms e.g. [*]CC for ethyl (as opposed to [CH2]C)<br>
	 * Note that if this is set to true InChIs cannot be generated
	 * @param outputRadicalsAsWildCardAtoms
	 */
	public void setOutputRadicalsAsWildCardAtoms(boolean outputRadicalsAsWildCardAtoms) {
		this.outputRadicalsAsWildCardAtoms = outputRadicalsAsWildCardAtoms;
	}

	/**
	 * Should OPSIN attempt reverse parsing to more accurately determine why parsing failed
	 * @return whether a more precise cause of failure should be determined if parsing fails
	 */
	public boolean isDetailedFailureAnalysis() {
		return detailedFailureAnalysis;
	}

	/**
	 * Sets whether OPSIN should attempt reverse parsing to more accurately determine why parsing failed
	 */
	public void setDetailedFailureAnalysis(boolean detailedFailureAnalysis) {
		this.detailedFailureAnalysis = detailedFailureAnalysis;
	}

	/**
	 * Are acids without the word "acid" interpretable e.g. should "acetic" be interpretable
	 * @return whether acids without the word "acid" should be interpretable
	 */
	public boolean allowInterpretationOfAcidsWithoutTheWordAcid() {
		return interpretAcidsWithoutTheWordAcid;
	}


	/**
	 * Sets whether acids without the word "acid" interpretable e.g. should "acetic" be interpretable
	 * @param interpretAcidsWithoutTheWordAcid
	 */
	public void setInterpretAcidsWithoutTheWordAcid(boolean interpretAcidsWithoutTheWordAcid) {
		this.interpretAcidsWithoutTheWordAcid = interpretAcidsWithoutTheWordAcid;
	}

	/**
	 * If OPSIN cannot understand the stereochemistry in a name should OPSIN's result be a warning
	 * and structure with incomplete stereochemistry, or should failure be returned (Default)
	 * @return whether ignored stereochemistry is a warning (rather than a failure)
	 */
	public boolean warnRatherThanFailOnUninterpretableStereochemistry() {
		return warnRatherThanFailOnUninterpretableStereochemistry;
	}


	/**
	 * Sets whether if OPSIN cannot understand the stereochemistry in a name whether OPSIN's result should be a warning
	 * and structure with incomplete stereochemistry, or should failure be returned (Default)
	 * @param warnRatherThanFailOnUninterpretableStereochemistry
	 */
	public void setWarnRatherThanFailOnUninterpretableStereochemistry(boolean warnRatherThanFailOnUninterpretableStereochemistry) {
		this.warnRatherThanFailOnUninterpretableStereochemistry = warnRatherThanFailOnUninterpretableStereochemistry;
	}


	/**
	 * Constructs a NameToStructureConfig with default settings:
	 * allowRadicals = false
	 * outputRadicalsAsWildCardAtoms = false
	 * detailedFailureAnalysis = false
	 * interpretAcidsWithoutTheWordAcid = false
	 * warnRatherThanFailOnUninterpretableStereochemistry = false
	 */
	public static NameToStructureConfig getDefaultConfigInstance() {
		return new NameToStructureConfig();
	}
	
	@Override
	public NameToStructureConfig clone() {
		try {
			return (NameToStructureConfig) super.clone();
		} catch (CloneNotSupportedException e) {
			// Can only be thrown if we *don't* implement Cloneable, which we do...
			throw new Error("Impossible!", e);
		}
	}
	
}
