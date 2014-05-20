package uk.ac.cam.ch.wwmm.opsin;

/**Thrown during component generation.
 *
 * @author ptc24
 *
 */
class ComponentGenerationException extends Exception {

	private static final long serialVersionUID = 1L;

	ComponentGenerationException() {
		super();
	}

	ComponentGenerationException(String message) {
		super(message);
	}

	ComponentGenerationException(String message, Throwable cause) {
		super(message, cause);
	}

	ComponentGenerationException(Throwable cause) {
		super(cause);
	}

}
