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

	public OPSIN_RESULT_STATUS getStatus() {
		return status;
	}

	public String getMessage() {
		return message;
	}
	
	public String getChemicalName() {
		return chemicalName;
	}

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
