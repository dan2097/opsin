package uk.ac.cam.ch.wwmm.opsin;

/**Thrown if a problem is encountered when converting a name to InChI
 *
 * @author dl387
 *
 */
public class InchiGenerationException extends Exception {
	private static final long serialVersionUID = 1L;

	InchiGenerationException() {
		super();
		// TODO Auto-generated constructor stub
	}

	InchiGenerationException(String message) {
		super(message);
		// TODO Auto-generated constructor stub
	}

	InchiGenerationException(String message, Throwable cause) {
		super(message, cause);
		// TODO Auto-generated constructor stub
	}

	InchiGenerationException(Throwable cause) {
		super(cause);
		// TODO Auto-generated constructor stub
	}
}
