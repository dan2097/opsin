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
	public NameToSmiles() throws Exception {
		n2s = NameToStructure.getInstance();
	}

	/**Parses a chemical name, returning an unambiguous SMILES representation of the molecule.
	 *
	 * @param name The chemical name to parse.
	 * @param verbose Whether to print lots of debugging information to stdin and stderr or not.
	 * @return A SMILES string, containing the parsed molecule, or null if the molecule would not parse.
	 */
	public String parseToSmiles(String name, boolean verbose) {
		Fragment frag = n2s.parseToOpsinFragment(name, verbose);
		if (frag == null){
			return null;
		}
		else{
			String smiles = null;
			try{
				smiles = opsinFragmentToSmiles(frag, verbose);
			}
			catch (Exception e) {
				if (verbose){
					e.printStackTrace();
				}
				return null;
			}
			if (smiles ==null){return null;}//smiles generation failed
			if(verbose) System.out.println(smiles);
			return smiles;
		}
	}
	

	private String opsinFragmentToSmiles(Fragment frag, boolean verbose){
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
