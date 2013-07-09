package uk.ac.cam.ch.wwmm.opsin;


/**
 * Thrown if stereochemistry cannot be applied to a structure
 * @author Daniel
 *
 */
class StereochemistryException extends StructureBuildingException {

	private static final long serialVersionUID = 1L;

	StereochemistryException() {
		super();
	}

	StereochemistryException(String message) {
		super(message);
	}

	StereochemistryException(String message, Throwable cause) {
		super(message, cause);
	}

	StereochemistryException(Throwable cause) {
		super(cause);
	}
}