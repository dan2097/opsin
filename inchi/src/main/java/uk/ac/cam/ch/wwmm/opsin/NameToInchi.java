package uk.ac.cam.ch.wwmm.opsin;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

import net.sf.jniinchi.INCHI_BOND_TYPE;
import net.sf.jniinchi.INCHI_OPTION;
import net.sf.jniinchi.INCHI_PARITY;
import net.sf.jniinchi.INCHI_RET;
import net.sf.jniinchi.JniInchiAtom;
import net.sf.jniinchi.JniInchiBond;
import net.sf.jniinchi.JniInchiException;
import net.sf.jniinchi.JniInchiInput;
import net.sf.jniinchi.JniInchiOutput;
import net.sf.jniinchi.JniInchiStereo0D;
import net.sf.jniinchi.JniInchiWrapper;
import nu.xom.Element;

public class NameToInchi {

	private NameToStructure n2s;
	public NameToInchi() throws NameToStructureException {
		n2s = NameToStructure.getInstance();
	}

	/**Parses a chemical name, returning an unambiguous InChI representation of the molecule.
	 *
	 * @param name The chemical name to parse.
	 * @param verbose Whether to print lots of debugging information to stdin and stderr or not.
	 * @return An InChI string, containing the parsed molecule, or null if the molecule would not parse.
	 */
	public String parseToInchi(String name, boolean verbose) {
		OpsinResult result = n2s.parseChemicalName(name, verbose);
		return convertResultToInChI(result, verbose);
	}
	
	/**
	 * Converts an OPSIN result to InChI. Null is returned if this conversion fails
	 * @param result
	 * @param verbose Whether to print lots of debugging information to stdin and stderr or not.
	 * @return String InChI
	 */
	public static String convertResultToInChI(OpsinResult result, boolean verbose){
		if (result.getStructure() !=null){
			String inchi = null;
			try{
				inchi = opsinFragmentToInchi(result.getStructure(), verbose);
			}
			catch (Exception e) {
				if (verbose){
					e.printStackTrace();
				}
				return null;
			}
			if (inchi ==null){return null;}//inchi generation failed
			if(verbose) System.out.println(inchi);
			return inchi;
		}
		return null;
	}

	private static String opsinFragmentToInchi(Fragment frag, boolean verbose) throws JniInchiException{
		HashMap<Integer, JniInchiAtom> opsinIdAtomMap = new HashMap<Integer, JniInchiAtom>();
		List<INCHI_OPTION> options = new ArrayList<INCHI_OPTION>();
		options.add(INCHI_OPTION.FixedH);
		JniInchiInput input = new JniInchiInput(options);
		List<Atom> atomList =frag.getAtomList();
		// Generate atoms
		for (Atom atom : atomList) {
			JniInchiAtom jAtom = input.addAtom(new JniInchiAtom(0.0, 0.0, 0.0, atom.getElement()));
			jAtom.setCharge(atom.getCharge());
			jAtom.setImplicitH(0);
			opsinIdAtomMap.put(atom.getID(), jAtom);
		}
		Set<Bond> bondList =frag.getBondSet();
		for (Bond bond : bondList) {
			input.addBond(new JniInchiBond(opsinIdAtomMap.get(bond.getFrom()), opsinIdAtomMap.get(bond.getTo()), INCHI_BOND_TYPE.getValue(bond.getOrder())));
		}

		for (Atom atom : atomList) {//add atomParities
			AtomParity atomParity =atom.getAtomParity();
        	if (atomParity != null){
        		Atom[] atomRefs4 = atomParity.getAtomRefs4();
				int[] atomRefs4AsInt = new int[4];
				for (int i = 0; i < atomRefs4.length; i++) {
					atomRefs4AsInt[i] = atomRefs4[i].getID();
				}
				INCHI_PARITY parity =INCHI_PARITY.UNKNOWN;
				if (atomParity.getParity() > 0){
					parity =INCHI_PARITY.EVEN;
				}
				else if (atomParity.getParity() < 0){
					parity =INCHI_PARITY.ODD;
				}
				input.addStereo0D(JniInchiStereo0D.createNewTetrahedralStereo0D(opsinIdAtomMap.get(atom.getID()), opsinIdAtomMap.get(atomRefs4AsInt[0]), opsinIdAtomMap.get(atomRefs4AsInt[1]), opsinIdAtomMap.get(atomRefs4AsInt[2]), opsinIdAtomMap.get(atomRefs4AsInt[3]), parity));
			}
        }

		for (Bond bond : bondList) {//add bondStereos
			Element bondStereoEl =bond.getBondStereoElement();
			if (bondStereoEl != null){
				String[] atomRefs4 = bondStereoEl.getAttributeValue("atomRefs4").split(" ");
				int[] atomRefs4AsInt = new int[4];
				for (int i = 0; i < atomRefs4.length; i++) {
					atomRefs4AsInt[i] = Integer.parseInt(atomRefs4[i].substring(1));//cut off starting a
				}
				if ("C".equals(bondStereoEl.getValue())){
					input.addStereo0D(JniInchiStereo0D.createNewDoublebondStereo0D(opsinIdAtomMap.get(atomRefs4AsInt[0]), opsinIdAtomMap.get(atomRefs4AsInt[1]), opsinIdAtomMap.get(atomRefs4AsInt[2]), opsinIdAtomMap.get(atomRefs4AsInt[3]), INCHI_PARITY.ODD));
				}
				else if ("T".equals(bondStereoEl.getValue())){
					input.addStereo0D(JniInchiStereo0D.createNewDoublebondStereo0D(opsinIdAtomMap.get(atomRefs4AsInt[0]), opsinIdAtomMap.get(atomRefs4AsInt[1]), opsinIdAtomMap.get(atomRefs4AsInt[2]), opsinIdAtomMap.get(atomRefs4AsInt[3]), INCHI_PARITY.EVEN));
				}
			}
        }
		JniInchiOutput output = JniInchiWrapper.getInchi(input);
		if (output ==null){
			return null;
		}
    	INCHI_RET ret = output.getReturnStatus();
    	if (verbose){
    		System.out.println("Inchi generation status: " + ret);
    		if (!INCHI_RET.OKAY.equals(ret)){
    			System.out.println(output.getMessage());
    		}
    	}
    	if (!INCHI_RET.OKAY.equals(ret) && !INCHI_RET.WARNING.equals(ret)){
    		return null;
    	}
    	return output.getInchi();
	}


	/**Run OPSIN as a standalone component for Inchi generation
	 *
	 * @param args
	 * @throws Exception
	 */
	public static void main(String [] args) throws Exception {
		NameToInchi nti = new NameToInchi();
		boolean end = false;
		BufferedReader stdinReader = new BufferedReader(new InputStreamReader(System.in));
		System.err.println("OPSIN Prealpha: enter chemical name:");
		while(!end) {
			String name = stdinReader.readLine();
			if(name == null || name.equals("END")) {
				end = true;
			} else {
				String output = nti.parseToInchi(name, false);
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
