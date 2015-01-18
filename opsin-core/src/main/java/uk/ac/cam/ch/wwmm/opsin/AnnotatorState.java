package uk.ac.cam.ch.wwmm.opsin;


/**
 * Contains the state needed during finite-state parsing
 * From this the tokens string and their semantics can be generated
 * @author Daniel
 *
 */
class AnnotatorState {

	/** The current state of the DFA. */
	private final int state;
	/** The annotation so far. */
	private final char annot;
	
	/** The index of the first char in the chemical name that has yet to be tokenised */
	private final int posInName;
	
	private final boolean isCaseSensitive;
	
	private final AnnotatorState previousAs;

	
	AnnotatorState(int state, char annot, int posInName, boolean isCaseSensitive, AnnotatorState previousAs) {
		this.state = state;
		this.annot = annot;
		this.posInName = posInName;
		this.isCaseSensitive = isCaseSensitive;
		this.previousAs = previousAs;
	}

	/**
	 * The current state in the DFA
	 * @return
	 */
	int getState() {
		return state;
	}

	/**
	 * The annotation that was consumed to transition to this state
	 * @return
	 */
	char getAnnot() {
		return annot;
	}

	/**
	 * The index of the first char in the chemical name that has yet to be tokenised (at the point of creating this AnnotatorState)
	 * @return
	 */
	int getPosInName() {
		return posInName;
	}

	/**
	 * Where the corresponding token is case sensitive
	 * @return
	 */
	boolean isCaseSensitive() {
		return isCaseSensitive;
	}

	/**
	 * The last annotator state for the previous token (or null if this is the first)
	 * @return
	 */
	AnnotatorState getPreviousAs() {
		return previousAs;
	}
}
