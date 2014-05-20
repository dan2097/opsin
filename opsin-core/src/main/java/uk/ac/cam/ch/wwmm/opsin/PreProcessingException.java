package uk.ac.cam.ch.wwmm.opsin;

/**Thrown during preprocessing.
 *
 * @author dl387
 *
 */
class PreProcessingException extends Exception {

	private static final long serialVersionUID = 1L;

	PreProcessingException() {
		super();
	}

	PreProcessingException(String message) {
		super(message);
	}

	PreProcessingException(String message, Throwable cause) {
		super(message, cause);
	}

	PreProcessingException(Throwable cause) {
		super(cause);
	}

}
