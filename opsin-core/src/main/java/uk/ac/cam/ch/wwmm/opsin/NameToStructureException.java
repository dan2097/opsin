package uk.ac.cam.ch.wwmm.opsin;

/**Thrown if OPSIN failed to initialise
 *
 * @author ptc24
 *
 */
public class NameToStructureException extends RuntimeException {

	private static final long serialVersionUID = 1L;

	NameToStructureException() {
		super();
	}

	NameToStructureException(String message) {
		super(message);
	}

	NameToStructureException(String message, Throwable cause) {
		super(message, cause);
	}

	NameToStructureException(Throwable cause) {
		super(cause);
	}

}
