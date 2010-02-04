package uk.ac.cam.ch.wwmm.opsin;

/**Thrown if a problem is encountered when converting a name to a depiction
 *
 * @author dl387
 *
 */
public class DepictionGenerationException extends Exception {
	private static final long serialVersionUID = 1L;

	DepictionGenerationException() {
		super();
		// TODO Auto-generated constructor stub
	}

	DepictionGenerationException(String message) {
		super(message);
		// TODO Auto-generated constructor stub
	}

	DepictionGenerationException(String message, Throwable cause) {
		super(message, cause);
		// TODO Auto-generated constructor stub
	}

	DepictionGenerationException(Throwable cause) {
		super(cause);
		// TODO Auto-generated constructor stub
	}
}
