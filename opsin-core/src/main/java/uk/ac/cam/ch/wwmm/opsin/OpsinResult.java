package uk.ac.cam.ch.wwmm.opsin;

import java.util.Collections;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.cam.ch.wwmm.opsin.OpsinWarning.OpsinWarningType;

/**
 * Holds the structure OPSIN has generated from a name
 * Additionally holds a status code for whether name interpretation was successful
 * @author dl387
 *
 */
public class OpsinResult {
	private static final Logger LOG = LogManager.getLogger(OpsinResult.class);
	private final Fragment structure;
	private final OPSIN_RESULT_STATUS status;
	private final String message;
	private final String chemicalName;
	private final List<OpsinWarning> warnings;

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
	
	OpsinResult(Fragment frag, OPSIN_RESULT_STATUS status, List<OpsinWarning> warnings, String chemicalName) {
		this.structure = frag;
		this.status = status;
		StringBuilder sb = new StringBuilder();
		for (int i = 0, l = warnings.size(); i < l; i++) {
			OpsinWarning warning = warnings.get(i);
			sb.append(warning.getType().toString());
			sb.append(": ");
			sb.append(warning.getMessage());
			if (i + 1 < l){
				sb.append("; ");
			}
		}
		this.message = sb.toString();
		this.chemicalName = chemicalName;
		this.warnings = warnings;
	}

	OpsinResult(Fragment frag, OPSIN_RESULT_STATUS status, String message, String chemicalName) {
		this.structure = frag;
		this.status = status;
		this.message = message;
		this.chemicalName = chemicalName;
		this.warnings = Collections.emptyList();
	}

	Fragment getStructure() {
		return structure;
	}

	/**
	 * Returns an enum indicating whether interpreting the chemical name was successful
	 * If an issue was identified but a chemical structure could be still be deduced the status is {@link OPSIN_RESULT_STATUS#WARNING}
	 * @return {@link OPSIN_RESULT_STATUS} status
	 */
	public OPSIN_RESULT_STATUS getStatus() {
		return status;
	}

	/**
	 * Returns a message explaining why generation of a molecule from the name failed
	 * This string will be blank when no problems were encountered
	 * @return String explaining problems encountered
	 */
	public String getMessage() {
		return message;
	}
	
	/**
	 * Returns the chemical name that this OpsinResult was generated from
	 * @return String containing the original chemical name
	 */
	public String getChemicalName() {
		return chemicalName;
	}

	/**
	 * Generates the CML corresponding to the molecule described by the name
	 * If name generation failed i.e. the OPSIN_RESULT_STATUS is FAILURE then null is returned
	 * @return Chemical Markup Language as a String
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
	 * @return Idented Chemical Markup Language as a String
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
	 * @return SMILES as a String
	 */
	public String getSmiles() {
		return getSmiles(SmilesOptions.DEFAULT);
	}

	/**
	 * Generates the SMILES corresponding to the molecule described by the name
	 * If name generation failed i.e. the OPSIN_RESULT_STATUS is FAILURE then null is returned.
	 * <br>
	 * The <code>options</code> parameter is used to control the output by a set of binary flags. This is
	 * primarily used to control the output layers in ChemAxon Extended SMILES (CXSMILES).
	 * <pre>
	 * // only the include the enhanced stereo layers
	 * result.getSmiles(SmilesOptions.CXSMILES_ENHANCED_STEREO);
	 * // only the include the enhanced stereo and polymer layers
	 * result.getSmiles(SmilesOptions.CXSMILES_ENHANCED_STEREO +
	 *                  SmilesOptions.CXSMILES_POLYMERS);
	 * </pre>
	 *
	 * @param options binary flags of {@link SmilesOptions} (default: {@link SmilesOptions#DEFAULT}))
	 * @return SMILES as a String
	 * @see SmilesOptions
	 */
	public String getSmiles(int options) {
		if (structure != null){
			try{
				return SMILESWriter.generateSmiles(structure, options);
			}
			catch (Exception e) {
				LOG.debug("SMILES generation failed", e);
			}
		}
		return null;
	}

	
	/**
	 * Experimental function that generates the extended SMILES corresponding to the molecule described by the name
	 * If name generation failed i.e. the OPSIN_RESULT_STATUS is FAILURE then null is returned
	 * If the molecule doesn't utilise any features made possible by extended SMILES this is equivalent to {@link #getSmiles()}
	 * @return Extended SMILES as a String
	 */
	public String getExtendedSmiles() {
		if (structure != null){
			try{
				return SMILESWriter.generateSmiles(structure, SmilesOptions.CXSMILES);
			}
			catch (Exception e) {
				LOG.debug("Extended SMILES generation failed", e);
			}
		}
		return null;
	}

	/**
	 * A list of warnings encountered when the result was {@link OPSIN_RESULT_STATUS#WARNING}<br>
	 * This list of warnings is immutable
	 * @return A list of {@link OpsinWarning}
	 */
	public List<OpsinWarning> getWarnings() {
		return Collections.unmodifiableList(warnings);
	}
	
	/**
	 * Convenience method to check if one of the associated OPSIN warnings was {@link OpsinWarningType#APPEARS_AMBIGUOUS}
	 * @return true if name appears to be ambiguous
	 */
	public boolean nameAppearsToBeAmbiguous() {
		for (OpsinWarning warning : warnings) {
			if (warning.getType() == OpsinWarningType.APPEARS_AMBIGUOUS) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Convenience method to check if one of the associated OPSIN warnings was {@link OpsinWarningType#STEREOCHEMISTRY_IGNORED}
	 * @return true if stereochemistry was ignored to interpret the name
	 */
	public boolean stereochemistryIgnored() {
		for (OpsinWarning warning : warnings) {
			if (warning.getType() == OpsinWarningType.STEREOCHEMISTRY_IGNORED) {
				return true;
			}
		}
		return false;
	}

}
