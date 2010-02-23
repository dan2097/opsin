package uk.ac.cam.ch.wwmm.opsin;

import nu.xom.Element;

public class OpsinResult {

	private final Fragment structure;
	private final OPSIN_RESULT_STATUS status;
	private final String message;
	private final String chemicalName;
	private Element cml = null;

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
	 * Lazily evaluates and returns the CML corresponding to the molecule described by the name
	 * If name generation failed i.e. the OPSIN_RESULT_STATUS is FAILURE then null is returned
	 * @return Element cml
	 */
	public synchronized Element getCml() {
		if (cml ==null && structure!=null){
			try{
				cml = structure.toCMLMolecule(chemicalName);
			}
			catch (Exception e) {
				//e.printStackTrace();//TODO use log4j
				cml = null;
			}
		}
		return cml;
	}

	enum OPSIN_RESULT_STATUS{
		SUCCESS,
		WARNING,
		FAILURE
	}
	
	OpsinResult(Fragment frag, OPSIN_RESULT_STATUS status, String message, String chemicalName){
		this.structure = frag;
		this.status = status;
		this.message = message;
		this.chemicalName = chemicalName;
	}
	
	
}
