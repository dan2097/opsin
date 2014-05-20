package uk.ac.cam.ch.wwmm.opsin;

/**Thrown during assembly of the structure
 *
 * @author ptc24
 *
 */
class StructureBuildingException extends Exception {

	private static final long serialVersionUID = 1L;

	StructureBuildingException() {
		super();
	}

	StructureBuildingException(String message) {
		super(message);
	}

	StructureBuildingException(String message, Throwable cause) {
		super(message, cause);
	}

	StructureBuildingException(Throwable cause) {
		super(cause);
	}

}
