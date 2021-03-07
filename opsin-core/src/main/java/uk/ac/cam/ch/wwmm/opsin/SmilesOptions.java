package uk.ac.cam.ch.wwmm.opsin;

/**
 * Options to control SMILES generation.
 * These options can be provided to the {@link OpsinResult#getSmiles(int)} method to control what is included in
 * the generated SMILES. The main use here is to control generation of ChemAxon Extended SMILES (CXSMILES) that supports
 * features beyond plain SMILES.
 * @see <a href="https://docs.chemaxon.com/display/docs/chemaxon-extended-smiles-and-smarts-cxsmiles-and-cxsmarts.md">ChemAxon Extended SMILES and SMARTS - CXSMILES and CXSMARTS</a>
 */
public interface SmilesOptions {
	/**
	 * Default SMILES generation, as Daylight intended.
	 */
	int DEFAULT                  = 0x0;
	/**
	 * Include atom labels in CXSMILES.
	 */
	int CXSMILES_ATOM_LABELS     = 0x1;
	/**
	 * Include atom values in CXSMILES - these are the locant values.
	 */
	int CXSMILES_ATOM_VALUES     = 0x2;
	/**
	 * Include repeat brackets in the CXSMILES layers for polymers.
	 */
	int CXSMILES_POLYMERS        = 0x4;
	/**
	 * Include racemic, relative, and absolute enhanced stereochemistry in the CXSMILES layers.
	 */
	int CXSMILES_ENHANCED_STEREO = 0x8;
	/**
	 * Include all CXSMILES layers that are relevant. This option is equivalent to turning on all CXSMILES features.s
	 */
	int CXSMILES                 = CXSMILES_ATOM_LABELS +
			                       CXSMILES_ATOM_VALUES +
			                       CXSMILES_POLYMERS +
			                       CXSMILES_ATOM_VALUES;
}
