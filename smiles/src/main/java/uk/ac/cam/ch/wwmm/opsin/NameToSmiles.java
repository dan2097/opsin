package uk.ac.cam.ch.wwmm.opsin;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.List;

import sea36.chem.base.ChemicalElement;
import sea36.chem.core.CMLAtom;
import sea36.chem.core.CMLAtomArray;
import sea36.chem.core.CMLAtomParity;
import sea36.chem.core.CMLMolecule;
import sea36.chem.io.smiles.SmilesWriter;

public class NameToSmiles {
	private NameToStructure n2s;
	public NameToSmiles() throws NameToStructureException {
		n2s = NameToStructure.getInstance();
	}

	/**Parses a chemical name, returning an unambiguous SMILES representation of the molecule.
	 *
	 * @param name The chemical name to parse.
	 * @param verbose Whether to print lots of debugging information to stdin and stderr or not.
	 * @return A SMILES string, containing the parsed molecule
	 * @throws SmilesGenerationException Thrown if parsing fails or SMILES could not be generated
	 */
	public String parseToSmiles(String name, boolean verbose) throws SmilesGenerationException {
		OpsinResult result = n2s.parseChemicalName(name, verbose);
		return convertResultToSMILES(result, verbose);
	}
	
	/**
	 * Converts an OPSIN result to SMILES. An exception is thrown if the conversion fails
	 * @param result
	 * @param verbose Whether to print lots of debugging information to stdin and stderr or not.
	 * @return String SMILES
	 * @throws SmilesGenerationException Thrown if conversion failed
	 */
	public static String convertResultToSMILES(OpsinResult result, boolean verbose) throws SmilesGenerationException{
		if (result.getStructure() !=null){
			String smiles = null;
			try{
				smiles = opsinFragmentToSmiles(result.getStructure(), verbose);
			}
			catch (Exception e) {
				if (verbose){
					e.printStackTrace();
				}
				throw new SmilesGenerationException("SMILES generation failed! This probably indicates a bug in the SMILES writer", e);
			}
			if (smiles ==null){
				throw new SmilesGenerationException("SMILES generation failed! This probably indicates a bug in the SMILES writer");
			}
			if(verbose) System.out.println(smiles);
			return smiles;
		}
		throw new SmilesGenerationException(result.getMessage());
	}

	private static String opsinFragmentToSmiles(Fragment frag, boolean verbose){
		OpsinToChemKitWrapper chemKitWrapper = new OpsinToChemKitWrapper(frag);
		CMLMolecule mol = chemKitWrapper.getChemKitMolecule();
		CMLAtomArray atoms = mol.getAtomArray();
		for (CMLAtom atom : atoms) {
            if (CMLAtomParity.findAtomParity(atom) ==null){
				List<CMLAtom> neighbours = mol.getNeighbourList(atom);
				int hydrogen = 0;
				for (CMLAtom neighbour : neighbours) {
					if (neighbour.getChemicalElement().equals(ChemicalElement.H) ){
						mol.removeAtomAndConnectedBonds(neighbour);
						hydrogen++;
					}
				}
				atom.setHydrogenCount(hydrogen);
            }
		}
		SmilesWriter sw = new SmilesWriter();
		return sw.writeMolecule(mol);
	}

	/**Run OPSIN as a standalone component for SMILES generation
	 *
	 * @param args
	 * @throws Exception
	 */
	public static void main(String [] args) throws Exception {
		NameToSmiles ntsmi = new NameToSmiles();
		boolean end = false;
		BufferedReader stdinReader = new BufferedReader(new InputStreamReader(System.in));
		System.out.println("OPSIN Prealpha: enter chemical name:");
		while(!end) {
			String name = stdinReader.readLine();
			if(name == null) {
				System.err.println("Disconnected!");
				end = true;
			} else if(name.equals("END")) {
				end = true;
			} else {
				String output = ntsmi.parseToSmiles(name, false);
				if(output == null) {
					System.out.println("Did not parse.");
					System.out.flush();
				} else {
					System.out.println(output);
					System.out.flush();
				}
			}
		}
	}
}
