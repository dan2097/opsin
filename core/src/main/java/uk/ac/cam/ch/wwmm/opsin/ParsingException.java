package uk.ac.cam.ch.wwmm.opsin;

/**Thrown during finite-state parsing.
 *
 * @author ptc24
 *
 */
class ParsingException extends Exception {

	private static final long serialVersionUID = 1L;

	ParsingException() {
		super();
	}

	ParsingException(String message) {
		super(message);
	}

	ParsingException(String message, Throwable cause) {
		super(message, cause);
	}

	ParsingException(Throwable cause) {
		super(cause);
	}

}
