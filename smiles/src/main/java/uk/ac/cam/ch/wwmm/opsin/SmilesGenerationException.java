package uk.ac.cam.ch.wwmm.opsin;

/**Thrown if a problem is encountered when converting a name to SMILES
 *
 * @author dl387
 *
 */
public class SmilesGenerationException extends Exception {
	private static final long serialVersionUID = 1L;

	SmilesGenerationException() {
		super();
		// TODO Auto-generated constructor stub
	}

	SmilesGenerationException(String message) {
		super(message);
		// TODO Auto-generated constructor stub
	}

	SmilesGenerationException(String message, Throwable cause) {
		super(message, cause);
		// TODO Auto-generated constructor stub
	}

	SmilesGenerationException(Throwable cause) {
		super(cause);
		// TODO Auto-generated constructor stub
	}
}
