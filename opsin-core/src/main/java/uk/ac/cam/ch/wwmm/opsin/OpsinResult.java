package uk.ac.cam.ch.wwmm.opsin;

import org.apache.log4j.Logger;

/**
 * Holds the structure OPSIN has generated from a name
 * Additionally holds a status code for whether name interpretation was successful
 * @author dl387
 *
 */
public class OpsinResult {
	private static final Logger LOG = Logger.getLogger(OpsinResult.class);
	private final Fragment structure;
	private final OPSIN_RESULT_STATUS status;
	private final String message;
	private final String chemicalName;

	/**
	 * Whether parsing the chemical name was successful, encountered problems or was unsuccessful.<br>
	 * If the result is not {@link OPSIN_RESULT_STATUS#FAILURE} then a structure has been generated
	 * @author dl387
	 *
	 */
	public enum OPSIN_RESULT_STATUS{
		/**
		 * OPSIN successfully interpreted the name
		 */
		SUCCESS,
		/**
		 * OPSIN interpreted the name but detected a potential problem e.g. could not interpret stereochemistry<br>
		 * Currently, by default, WARNING is not used as stereochemistry failures are treated as failures<br>
		 * In the future, ambiguous chemical names may produce WARNING
		 */
		WARNING,
		/**
		 * OPSIN failed to interpret the name
		 */
		FAILURE
	}

	OpsinResult(Fragment frag, OPSIN_RESULT_STATUS status, String message, String chemicalName){
		this.structure = frag;
		this.status = status;
		this.message = message;
		this.chemicalName = chemicalName;
	}

	Fragment getStructure() {
		return structure;
	}

	/**
	 * Returns an enum with values SUCCESS, WARNING and FAILURE
	 * Currently warning is never used
	 * @return OPSIN_RESULT_STATUS status
	 */
	public OPSIN_RESULT_STATUS getStatus() {
		return status;
	}

	/**
	 * Returns a message explaining why generation of a molecule from the name failed
	 * This string will be blank when no problems were encountered
	 * @return String message
	 */
	public String getMessage() {
		return message;
	}
	
	/**
	 * Returns the chemical name that this OpsinResult was generated frm
	 * @return String chemicalName
	 */
	public String getChemicalName() {
		return chemicalName;
	}

	/**
	 * Generates the CML corresponding to the molecule described by the name
	 * If name generation failed i.e. the OPSIN_RESULT_STATUS is FAILURE then null is returned
	 * @return String cml
	 */
	public String getCml() {
		if (structure != null){
			try{
				return CMLWriter.generateCml(structure, chemicalName);
			}
			catch (Exception e) {
				LOG.debug("CML generation failed", e);
			}
		}
		return null;
	}
	
	/**
	 * Generates the CML corresponding to the molecule described by the name
	 * If name generation failed i.e. the OPSIN_RESULT_STATUS is FAILURE then null is returned
	 * The CML is indented
	 * @return String cml
	 */
	public String getPrettyPrintedCml() {
		if (structure != null){
			try{
				return CMLWriter.generateIndentedCml(structure, chemicalName);
			}
			catch (Exception e) {
				LOG.debug("CML generation failed", e);
			}
		}
		return null;
	}

	/**
	 * Generates the SMILES corresponding to the molecule described by the name
	 * If name generation failed i.e. the OPSIN_RESULT_STATUS is FAILURE then null is returned
	 * @return String smiles
	 */
	public String getSmiles() {
		if (structure != null){
			try{
				return new SMILESWriter(structure).generateSmiles();
			}
			catch (Exception e) {
				LOG.debug("SMILES generation failed", e);
			}
		}
		return null;
	}
	
	
}
