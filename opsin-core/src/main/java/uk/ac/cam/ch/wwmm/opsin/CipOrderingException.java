package uk.ac.cam.ch.wwmm.opsin;

/**Thrown if the ordering of ligands can now be determined by OPSIN's implementation of the CIP rules.
 * This could be due to a limitation of the implementation or ligands actually being indistinguishable
 *
 * @author dl387
 *
 */
public class CipOrderingException extends RuntimeException {

	private static final long serialVersionUID = 1L;

	CipOrderingException() {
		super();
	}

	CipOrderingException(String message) {
		super(message);
	}

	CipOrderingException(String message, Throwable cause) {
		super(message, cause);
	}

	CipOrderingException(Throwable cause) {
		super(cause);
	}

}
