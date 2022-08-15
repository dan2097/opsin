package uk.ac.cam.ch.wwmm.opsin;

import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import uk.ac.cam.ch.wwmm.opsin.Bond.SMILES_BOND_DIRECTION;
import uk.ac.cam.ch.wwmm.opsin.BondStereo.BondStereoValue;

/** A builder for fragments specified as SMILES. A slightly custom SMILES dialect is used.
 * It includes all common features of SMILES and a few useful extensions:
 * | is used within a square bracketed element to directly set valency e.g. [P|5]. This is the same as using the lambda convention
 * sb/te are allowed (aromatic antimony/tellurium):
 * H? e.g. [SeH?] is used to indicate that the atom should use the default valency. It is equivalent to not using square brackets for organic atoms
 *
 * Allowed:
 * Organic elements B,C,N,O,P,S,F,Cl,Br,I (square brackets not required)
 * Aromatic elements c,n,o,p,s (square brackets not required) si,as,se,sb,te (square brackets required) Note that the inclusion of si/sb/te are an unofficial extension
 * =, # for bond orders
 * . for disconnection
 * (, ) for branching
 * [, ] for placing inorganic elements within and specifying charge. Allowed: [Al3+] or [Al+++]
 * 012345679 - ring closures
 * %10 %99 - more ring closures (%100 is ring closure %10 and 0 as in normal SMILES)
 * / and \ to set double bond stereochemistry to cis/trans
 * @ and @@ to set tetrahedral stereochemistry as in SMILES.
 * Hx where x is a digit is used to sort of set the hydrogen. In actuality the valency of the atom is derived and a valency hint added to the atom
 * This valency hint is the minimum valency that atom may be in. H? as an extension gives you the lowest acceptable valency.
 * |3 |5 etc. can be used to set the valency of an atom e.g.  [Se|2]
 *
 * Also, an = or # at the start of the string indicates that the group attaches to its parent group via a double or triple bond.
 *
 * A -,=,# on the end indicates that in the absence of locants, other groups attach to
 * *it* via the atom at the end of the string, not at the start of the string with -,=,# meaning single,double or triple bond
 * This behaviour is overridden for certain suffixes to give different meanings to the atom the -,=,# is referring to
 *
 * @author ptc24
 * @author dl387
 *
 */
class SMILESFragmentBuilder {

	/**A "struct" to hold information on the parsing stack
	 *
	 * @author ptc24
	 *
	 */
	private static class StackFrame {
		/**The Atom currently under consideration.*/
		Atom atom;

		/**The order of the bond about to be formed.*/
		int bondOrder;

		/**Whether the bond is a \ or / bond for use in determining cis/trans.*/
		SMILES_BOND_DIRECTION slash = null;

		/**The index of a dummy atom in the atom's stereochemistry atomrefs4*/
		Integer indexOfDummyAtom = null;

		/**Creates a stack frame with given parameters.
		 *
		 * @param a An atom or null
		 * @param bondOrderVal The value for bondOrder.
		 */
		StackFrame(Atom a, int bondOrderVal) {
			atom = a;
			bondOrder = bondOrderVal;
		}

		/**Creates a copy of an existing StackFrame.
		 *
		 * @param sf The stackframe to copy.
		 */
		StackFrame(StackFrame sf) {
			atom = sf.atom;
			bondOrder = sf.bondOrder;
		}
	}

	/**Ring opening dummy atom, used as a placeholder in stereochemistry atomrefs4*/
	private static final Atom ringOpeningDummyAtom = new Atom(ChemEl.R);

	/**Organic Atoms.*/
	private static final Set<String> organicAtoms = new HashSet<>();
	/**Aromatic Atoms.*/
	private static final Set<String> aromaticAtoms = new HashSet<>();

	static {
		organicAtoms.add("B");
		organicAtoms.add("C");
		organicAtoms.add("N");
		organicAtoms.add("O");
		organicAtoms.add("P");
		organicAtoms.add("S");
		organicAtoms.add("F");
		organicAtoms.add("Cl");
		organicAtoms.add("Br");
		organicAtoms.add("I");

		aromaticAtoms.add("c");
		aromaticAtoms.add("n");
		aromaticAtoms.add("o");
		aromaticAtoms.add("p");
		aromaticAtoms.add("s");
		aromaticAtoms.add("si");
		aromaticAtoms.add("as");
		aromaticAtoms.add("se");
		aromaticAtoms.add("sb");
		aromaticAtoms.add("te");
	}
	
	private final IDManager idManager;
	
	SMILESFragmentBuilder(IDManager idManager) {
		this.idManager = idManager;
	}

	private class ParserInstance {
		private final Deque<StackFrame> stack = new ArrayDeque<>();
		private final Map<String, StackFrame> ringClosures = new HashMap<>();
		
		private final String smiles;
		private final int endOfSmiles;
		private final Fragment fragment;
		private final int firstAtomOutValency;
		private final int lastAtomOutValency;
		
		private int i;

		ParserInstance(String smiles, Fragment fragment) {
			this.smiles = smiles;
			this.fragment = fragment;

			int lastIndex = smiles.length();
			
			char firstChar = smiles.charAt(0);//used by OPSIN to specify the valency with which this fragment connects
			if (firstChar == '-') {
				this.firstAtomOutValency = 1;
				this.i = 1;
			}
			else if (firstChar == '=') {
				this.firstAtomOutValency = 2;
				this.i = 1;
			}
			else if (firstChar == '#') {
				this.firstAtomOutValency = 3;
				this.i = 1;
			}
			else {
				this.firstAtomOutValency = -1;
				this.i = 0;
			}
			
			char lastChar = smiles.charAt(lastIndex - 1);//used by OPSIN to specify the valency with which this fragment connects and to indicate it connects via the last atom in the SMILES
			if (lastChar == '-') {
				this.lastAtomOutValency = 1;
				this.endOfSmiles = lastIndex - 1;
			}
			else if (lastChar == '=') {
				this.lastAtomOutValency = 2;
				this.endOfSmiles = lastIndex - 1;
			}
			else if (lastChar == '#') {
				this.lastAtomOutValency = 3;
				this.endOfSmiles = lastIndex - 1;
			}
			else {
				this.lastAtomOutValency = -1;
				this.endOfSmiles = lastIndex;
			}
		}
		
		void parseSmiles() throws StructureBuildingException {
			stack.add(new StackFrame(null, 1));
			for (; i < endOfSmiles; i++) {
				char ch = smiles.charAt(i);
				switch (ch) {
				case '(':
					stack.add(new StackFrame(stack.getLast()));
					break;
				case ')':
					stack.removeLast();
					break;
				case '-':
					stack.getLast().bondOrder = 1;
					break;
				case '=':
					if (stack.getLast().bondOrder != 1){
						throw new StructureBuildingException("= in unexpected position: bond order already defined!");
					}
					stack.getLast().bondOrder = 2;
					break;
				case '#':
					if (stack.getLast().bondOrder != 1){
						throw new StructureBuildingException("# in unexpected position: bond order already defined!");
					}
					stack.getLast().bondOrder = 3;
					break;
				case '/':
					if (stack.getLast().slash != null){
						throw new StructureBuildingException("/ in unexpected position: bond configuration already defined!");
					}
					stack.getLast().slash = SMILES_BOND_DIRECTION.RSLASH;
					break;
				case '\\':
					if (stack.getLast().slash != null){
						throw new StructureBuildingException("\\ in unexpected position: bond configuration already defined!");
					}
					stack.getLast().slash = SMILES_BOND_DIRECTION.LSLASH;
					break;
				case '.':
					stack.getLast().atom = null;
					break;
				case 'a':
				case 'b':
				case 'c':
				case 'd':
				case 'e':
				case 'f':
				case 'g':
				case 'h':
				case 'i':
				case 'j':
				case 'k':
				case 'l':
				case 'm':
				case 'n':
				case 'o':
				case 'p':
				case 'q':
				case 'r':
				case 's':
				case 't':
				case 'u':
				case 'v':
				case 'w':
				case 'x':
				case 'y':
				case 'z':
				case 'A':
				case 'B':
				case 'C':
				case 'D':
				case 'E':
				case 'F':
				case 'G':
				case 'H':
				case 'I':
				case 'J':
				case 'K':
				case 'L':
				case 'M':
				case 'N':
				case 'O':
				case 'P':
				case 'Q':
				case 'R':
				case 'S':
				case 'T':
				case 'U':
				case 'V':
				case 'W':
				case 'X':
				case 'Y':
				case 'Z':
				case '*':
					processOrganicAtom(ch);
					break;
				case '[':
					processBracketedAtom();
					break;
				case '0':
				case '1':
				case '2':
				case '3':
				case '4':
				case '5':
				case '6':
				case '7':
				case '8':
				case '9':
				case '%':
					processRingOpeningOrClosure(ch);
					break;
				default: 
					throw new StructureBuildingException(ch + " is in an unexpected position. Check this is not a mistake and that this feature of SMILES is supported by OPSIN's SMILES parser");
				}
			}
			if (!ringClosures.isEmpty()){
				throw new StructureBuildingException("Unmatched ring opening");
			}
			
			if (firstAtomOutValency > 0) {
				fragment.addOutAtom(fragment.getFirstAtom(), firstAtomOutValency, true);
			}

			if (lastAtomOutValency > 0) {
				//note that in something like C(=O)- this would be the carbon not the oxygen
				fragment.addOutAtom(getInscopeAtom(), lastAtomOutValency, true);
			}
		}

		/**
		 * An organic atom e.g. 'C', 'Cl', 'c' etc.
		 * @param ch
		 * @throws StructureBuildingException
		 */
		private void processOrganicAtom(char ch) throws StructureBuildingException {
			String elementType = String.valueOf(ch);
			boolean spareValency = false;
			if(is_A_to_Z(ch)) {//normal atoms
				if(i + 1 < endOfSmiles && is_a_to_z(smiles.charAt(i + 1)) && organicAtoms.contains(smiles.substring(i, i + 2))) {
					elementType = smiles.substring(i, i + 2);
					i++;
				}
				else if (!organicAtoms.contains(elementType)){
					throw new StructureBuildingException(elementType + " is not an organic Element. If it is actually an element it should be in square brackets");
				}
			}
			else if(is_a_to_z(ch)) {//aromatic atoms
				if (!aromaticAtoms.contains(elementType)){
					throw new StructureBuildingException(elementType + " is not an aromatic Element. If it is actually an element it should not be in lower case");
				}
				elementType = String.valueOf((char)(ch - 32));
				spareValency = true;
			}
			else if (ch == '*') {
				elementType = "R";
			}
			Atom atom = createAtom(elementType, fragment);
			atom.setSpareValency(spareValency);
			fragment.addAtom(atom);
		
			StackFrame currentFrame = stack.getLast();
			if(currentFrame.atom != null) {
				Bond b = createBond(currentFrame.atom, atom, currentFrame.bondOrder);
				if (currentFrame.slash != null){
					b.setSmilesStereochemistry(currentFrame.slash);
					currentFrame.slash = null;
				}
				if (currentFrame.atom.getAtomParity() != null){
					addAtomToAtomParity(currentFrame.atom.getAtomParity(), atom);
				}
			}
			currentFrame.atom = atom;
			currentFrame.bondOrder = 1;
		}

		/**
		 * square brackets- contain non-organic atoms or where required to set properties such as charge/chirality etc.
		 * e.g. [Na+]
		 * @throws StructureBuildingException
		 */
		private void processBracketedAtom() throws StructureBuildingException {
			i++;
			int indexOfRightSquareBracket = smiles.indexOf(']', i);
			if (indexOfRightSquareBracket == -1) {
				throw new StructureBuildingException("[ without matching \"]\"");
			}
			// isotope
			String isotope = "";
			while(is_0_to_9(smiles.charAt(i))) {
				isotope += smiles.charAt(i);
				i++;
			}

			char ch;
			if (i < indexOfRightSquareBracket){
				ch = smiles.charAt(i);
				i++;
			}
			else{
				throw new StructureBuildingException("No element found in square brackets");
			}
			// elementType
			String elementType = String.valueOf(ch);
			boolean spareValency = false;
			if(is_A_to_Z(ch)) {//normal atoms
				if(is_a_to_z(smiles.charAt(i))) {
					elementType += smiles.charAt(i);
					i++;
				}
			}
			else if(is_a_to_z(ch)) {//aromatic atoms
				if(is_a_to_z(smiles.charAt(i))) {
					if (aromaticAtoms.contains(elementType + smiles.charAt(i))){
						elementType = String.valueOf((char)(ch - 32)) + smiles.charAt(i);
						i++;
					}
					else{
						throw new StructureBuildingException(elementType + smiles.charAt(i) + " is not an aromatic Element. If it is actually an element it should not be in lower case");
					}
				}
				else{
					if (!aromaticAtoms.contains(elementType)){
						throw new StructureBuildingException(elementType + " is not an aromatic Element.");
					}
					elementType = String.valueOf((char)(ch - 32));
				}
				spareValency = true;
			}
			else if (elementType.equals("*")){
				elementType = "R";
			}
			else{
				throw new StructureBuildingException(elementType + " is not a valid element type!");
			}
			Atom atom = createAtom(elementType, fragment);
			atom.setSpareValency(spareValency);
			if (isotope.length() > 0){
				atom.setIsotope(Integer.parseInt(isotope));
			}
			fragment.addAtom(atom);
			StackFrame currentFrame = stack.getLast();
			if(currentFrame.atom != null) {
				Bond b = createBond(currentFrame.atom, atom, currentFrame.bondOrder);
				if (currentFrame.slash != null){
					b.setSmilesStereochemistry(currentFrame.slash);
					currentFrame.slash = null;
				}
				if (currentFrame.atom.getAtomParity() != null){
					addAtomToAtomParity(currentFrame.atom.getAtomParity(), atom);
				}
			}
			Atom previousAtom = currentFrame.atom;//needed for setting atomParity elements up
			currentFrame.atom = atom;
			currentFrame.bondOrder = 1;

			Integer hydrogenCount = 0;
			int charge = 0;
			Boolean chiralitySet = false;
			for (; i < indexOfRightSquareBracket; i++) {
				ch = smiles.charAt(i);
				if(ch == '@') {// chirality-sets atom parity
					if (chiralitySet){
						throw new StructureBuildingException("Atom parity appeared to be specified twice for an atom in a square bracket!");
					}
					processTetrahedralStereochemistry(atom, previousAtom, fragment.getAtomCount() == 1);
					chiralitySet = true;
				}
				else if (ch == 'H'){// hydrogenCount
					if (hydrogenCount == null || hydrogenCount != 0){
						throw new StructureBuildingException("Hydrogen count appeared to be specified twice for an atom in a square bracket!");
					}
					if (smiles.charAt(i + 1) == '?'){
						//extension to allow standard valency (as determined by the group in the periodic table) to dictate hydrogens
						i++;
						hydrogenCount = null;
					}
					else{
						String hydrogenCountString ="";
						while(is_0_to_9(smiles.charAt(i + 1))) {
							hydrogenCountString += smiles.charAt(i + 1);
							i++;
						}
						if (hydrogenCountString.length() == 0){
							hydrogenCount = 1;
						}
						else{
							hydrogenCount = Integer.parseInt(hydrogenCountString);
						}
						if (atom.hasSpareValency()) {
							if ((!elementType.equals("C") && !elementType.equals("Si")) || hydrogenCount >=2){
								fragment.addIndicatedHydrogen(atom);
							}
						}
					}
				}
				else if(ch == '+' || ch == '-') {// formalCharge
					if (charge != 0){
						throw new StructureBuildingException("Charge appeared to be specified twice for an atom in a square bracket!");
					}
					charge = (ch == '+') ? 1 : -1;
					String changeChargeStr = "";
					int changeCharge = 1;
					while(is_0_to_9(smiles.charAt(i + 1))) {//e.g. [C+2]
						changeChargeStr += smiles.charAt(i + 1);
						i++;
					}
					if (changeChargeStr.length() == 0){
						while(i + 1 < indexOfRightSquareBracket){//e.g. [C++]
							ch = smiles.charAt(i + 1);
							if (ch == '+'){
								if (charge != 1){
									throw new StructureBuildingException("Atom has both positive and negative charges specified!");//e.g. [C+-]
								}
							}
							else if (ch == '-'){
								if (charge != -1){
									throw new StructureBuildingException("Atom has both negative and positive charges specified!");
								}
							}
							else{
								break;
							}
							changeCharge++;
							i++;
						}
					}
					changeCharge = changeChargeStr.length() == 0 ? changeCharge : Integer.parseInt(changeChargeStr);
					atom.setCharge(charge * changeCharge);
				}
				else if(ch == '|') {
					StringBuilder lambda = new StringBuilder();
					while(i < endOfSmiles && is_0_to_9(smiles.charAt(i + 1))) {
						lambda.append(smiles.charAt(i + 1));
						i++;
					}
					atom.setLambdaConventionValency(Integer.parseInt(lambda.toString()));
				}
				else{
					throw new StructureBuildingException("Unexpected character found in square bracket");
				}
			}
			atom.setProperty(Atom.SMILES_HYDROGEN_COUNT, hydrogenCount);
		}

		/**
		 * Adds an atomParity element to the given atom using the information at the current index
		 * @param atom
		 * @param previousAtom
		 * @param isFirstAtom 
		 */
		private void processTetrahedralStereochemistry(Atom atom, Atom previousAtom, boolean isFirstAtom){
			Boolean chiralityClockwise = false;
			if (smiles.charAt(i + 1) == '@'){
				chiralityClockwise = true;
				i++;
			}
			Atom[] atomRefs4 = new Atom[4];
			AtomParity atomParity = new AtomParity(atomRefs4, chiralityClockwise ? 1 : -1);
			int index =0;
			if (previousAtom != null){
				atomRefs4[index] = previousAtom;
				index++;
			}
			else if (isFirstAtom && firstAtomOutValency == 1) {
				atomRefs4[index] = AtomParity.deoxyHydrogen;
				index++;
			}
			if (smiles.charAt(i + 1) == 'H'){
				atomRefs4[index] = AtomParity.hydrogen;
				//this character will also be checked by the hydrogen count check, hence don't increment i
			}
			atom.setAtomParity(atomParity);
		}
		
		/**
		 * Process ring openings and closings e.g. the two 1s in c1ccccc1
		 * @param ch
		 * @throws StructureBuildingException
		 */
		private void processRingOpeningOrClosure(char ch) throws StructureBuildingException {
			String closure = String.valueOf(ch);
			if(ch == '%') {
				if (i + 2 < endOfSmiles && is_0_to_9(smiles.charAt(i + 1)) && is_0_to_9(smiles.charAt(i + 2))) {
					closure = smiles.substring(i + 1, i + 3);
					i +=2;
				}
				else{
					throw new StructureBuildingException("A ring opening indice after a % must be two digits long");
				}
			}
			if(ringClosures.containsKey(closure)) {
				processRingClosure(closure);
			} else {
				if (getInscopeAtom() == null){
					throw new StructureBuildingException("A ring opening has appeared before any atom!");
				}
				processRingOpening(closure);
			}
		}

		private void processRingOpening(String closure) throws StructureBuildingException {
			StackFrame currentFrame = stack.getLast();
			StackFrame sf = new StackFrame(currentFrame);
			if (currentFrame.slash != null){
				sf.slash = currentFrame.slash;
				currentFrame.slash = null;
			}
			AtomParity atomParity = sf.atom.getAtomParity();
			if (atomParity != null){//replace ringclosureX with actual reference to id when it is known
				sf.indexOfDummyAtom = addAtomToAtomParity(atomParity, ringOpeningDummyAtom);
			}
			ringClosures.put(closure, sf);
			currentFrame.bondOrder = 1;
		}

		private void processRingClosure(String closure) throws StructureBuildingException {
			StackFrame sf = ringClosures.remove(closure);
			StackFrame currentFrame = stack.getLast();
			int bondOrder = 1;
			if(sf.bondOrder > 1) {
				if(currentFrame.bondOrder > 1 && sf.bondOrder != currentFrame.bondOrder){
					throw new StructureBuildingException("ring closure has two different bond orders specified!");
				}
				bondOrder = sf.bondOrder;
			} else if(currentFrame.bondOrder > 1) {
				bondOrder = currentFrame.bondOrder;
			}
			Bond b;
			if (currentFrame.slash != null) {
				//stereochemistry specified on ring closure
				//special case e.g. CC1=C/F.O\1  Bond is done from the O to the the C due to the presence of the \
				b = createBond(currentFrame.atom, sf.atom, bondOrder);
				b.setSmilesStereochemistry(currentFrame.slash);
				if(sf.slash != null && sf.slash.equals(currentFrame.slash)) {//specified twice check for contradiction
					throw new StructureBuildingException("Contradictory double bond stereoconfiguration");
				}
				currentFrame.slash = null;
			}
			else {
				b = createBond(sf.atom, currentFrame.atom, bondOrder);
				if (sf.slash != null) {
					//stereochemistry specified on ring opening
					b.setSmilesStereochemistry(sf.slash);
				}
			}

			AtomParity currentAtomParity = currentFrame.atom.getAtomParity();
			if (currentAtomParity != null) {
				addAtomToAtomParity(currentAtomParity, sf.atom);
			}
			
			AtomParity closureAtomParity = sf.atom.getAtomParity();
			if (closureAtomParity != null) {//replace dummy atom with actual atom e.g. N[C@@H]1C.F1 where the 1 initially holds a dummy atom before being replaced with the F atom
				Atom[] atomRefs4 = closureAtomParity.getAtomRefs4();
				if (sf.indexOfDummyAtom == null) {
					throw new RuntimeException("OPSIN Bug: Index of dummy atom representing ring closure atom not set");
				}
				atomRefs4[sf.indexOfDummyAtom] = currentFrame.atom;
			}
			currentFrame.bondOrder = 1;
		}

		/**
		 * Adds an atom at the first non-null position in the atomParity's atomRefs4
		 * @param atomParity
		 * @param atom
		 * @return Returns the index of the atom in the atomParity's atomRefs4
		 * @throws StructureBuildingException
		 */
		private int addAtomToAtomParity(AtomParity atomParity, Atom atom) throws StructureBuildingException {
			Atom[] atomRefs4 = atomParity.getAtomRefs4();
			boolean setAtom = false;
			int i = 0;
			for (; i < atomRefs4.length; i++) {
				if (atomRefs4[i] == null){
					atomRefs4[i] = atom;
					setAtom = true;
					break;
				}
			}
			if (!setAtom){
				throw new StructureBuildingException("Tetrahedral stereocentre specified in SMILES appears to involve more than 4 atoms");
			}
			return i;
		}
		
		/**
		 * For non-empty SMILES will return the atom at the top of the stack i.e. the one that will be bonded to next if the SMILES continued
		 * (only valid during execution of and after {@link ParserInstance#parseSmiles()} has been called)
		 * @return
		 */
		Atom getInscopeAtom(){
			return stack.getLast().atom;
		}
	}
	
	/**
	 * Build a Fragment based on a SMILES string.
	 * The type/subType of the Fragment are the empty String
	 * The fragment has no locants
	 *
	 * @param smiles The SMILES string to build from.
	 * @return The built fragment.
	 * @throws StructureBuildingException
	 */
	Fragment build(String smiles) throws StructureBuildingException {
		return build(smiles, "", NONE_LABELS_VAL);
	}
	
	/**
	 * Build a Fragment based on a SMILES string.
	 * @param smiles The SMILES string to build from.
	 * @param type The type of the fragment retrieved when calling {@link Fragment#getType()}
	 * @param labelMapping A string indicating which locants to assign to each atom. Can be a slash delimited list, "numeric", "fusedRing" or "none"/""
	 * @return
	 * @throws StructureBuildingException
	 */
	Fragment build(String smiles, String type, String labelMapping) throws StructureBuildingException {
		return build(smiles, new Fragment(type), labelMapping);
	}

	/**
	 * Build a Fragment based on a SMILES string.
	 * @param smiles The SMILES string to build from.
	 * @param tokenEl The corresponding tokenEl
	 * @param labelMapping A string indicating which locants to assign to each atom. Can be a slash delimited list, "numeric", "fusedRing" or "none"/""
	 * @return Fragment The built fragment.
	 * @throws StructureBuildingException
	 */
	Fragment build(String smiles, Element tokenEl, String labelMapping) throws StructureBuildingException {
		if (tokenEl == null){
			throw new IllegalArgumentException("tokenEl is null. FragmentManager's DUMMY_TOKEN should be used instead");
		}
		return build(smiles, new Fragment(tokenEl), labelMapping);
	}
	
	private Fragment build(String smiles, Fragment fragment, String labelMapping) throws StructureBuildingException {	
		if (smiles == null) {
			throw new IllegalArgumentException("SMILES specified is null");
		}
		if (labelMapping == null) {
			throw new IllegalArgumentException("labelMapping is null use \"none\" if you do not want any numbering or \"numeric\" if you would like default numbering");
		}
		if (smiles.isEmpty()){
			return fragment;
		}
		ParserInstance instance = new ParserInstance(smiles, fragment);
		instance.parseSmiles();
		
		List<Atom> atomList = fragment.getAtomList();
		processLabelling(labelMapping, atomList);

		verifyAndTakeIntoAccountLonePairsInAtomParities(atomList);
		addBondStereoElements(fragment);
		

		for (Atom atom : atomList) {
			if (atom.getProperty(Atom.SMILES_HYDROGEN_COUNT) != null && atom.getLambdaConventionValency() == null){
				setupAtomValency(atom);
			}
		}
		CycleDetector.assignWhetherAtomsAreInCycles(fragment);
		return fragment;
	}

	private void processLabelling(String labelMapping, List<Atom> atomList) throws StructureBuildingException {
		if (labelMapping.equals(NONE_LABELS_VAL) || labelMapping.length() == 0) {
			return;
		}
		if (labelMapping.equals(NUMERIC_LABELS_VAL)) {
			int atomNumber = 1;
			for (Atom atom : atomList) {
				atom.addLocant(Integer.toString(atomNumber++));
			}
		}
		else if(labelMapping.equals(FUSEDRING_LABELS_VAL)) {//fragment is a fusedring with atoms in the correct order for fused ring numbering
			//this will do stuff like changing labels from 1,2,3,4,5,6,7,8,9,10->1,2,3,4,4a,5,6,7,8,8a
			FragmentTools.relabelLocantsAsFusedRingSystem(atomList);
		}
		else{
			String[] labelMap = labelMapping.split("/", -1);//place slash delimited labels into an array
			int numOfAtoms = atomList.size();
			if (labelMap.length != numOfAtoms){
				throw new StructureBuildingException("Group numbering has been invalidly defined in resource file: labels: " +labelMap.length + ", atoms: " + numOfAtoms );
			}
			for (int i = 0; i < numOfAtoms; i++) {
				String labels[] = labelMap[i].split(",");
				for (String label : labels) {
					if (label.length() > 0) {
						atomList.get(i).addLocant(label);
					}
				}
			}
		}
	}

	private void verifyAndTakeIntoAccountLonePairsInAtomParities(List<Atom> atomList) throws StructureBuildingException {
		for (Atom atom : atomList) {
			AtomParity atomParity = atom.getAtomParity();
			if (atomParity != null){
				Atom[] atomRefs4 = atomParity.getAtomRefs4();
				int nullAtoms = 0;
				int hydrogen = 0;
				for (Atom atomRefs4Atom : atomRefs4) {
					if (atomRefs4Atom == null){
						nullAtoms++;
					}
					else if (atomRefs4Atom.equals(AtomParity.hydrogen)){
						hydrogen++;
					}
				}
				if (nullAtoms != 0){
					if (nullAtoms ==1 && hydrogen==0 && 
							(atom.getElement() == ChemEl.N || atom.getElement() == ChemEl.S || atom.getElement() == ChemEl.Se)){//special case where lone pair is part of the tetrahedron
						if (atomList.indexOf(atomRefs4[0]) < atomList.indexOf(atom)){//is there an atom in the SMILES in front of the stereocentre?
							atomRefs4[3] = atomRefs4[2];
							atomRefs4[2] = atomRefs4[1];
							atomRefs4[1] = atom;
						}
						else{
							atomRefs4[3] = atomRefs4[2];
							atomRefs4[2] = atomRefs4[1];
							atomRefs4[1] = atomRefs4[0];
							atomRefs4[0] = atom;
						}
					}
					else{
						throw new StructureBuildingException("SMILES is malformed. Tetrahedral stereochemistry defined on a non tetrahedral centre");
					}
				}
			}
		}
	}

	private void addBondStereoElements(Fragment currentFrag) throws StructureBuildingException {
		Set<Bond> bonds = currentFrag.getBondSet();
		for (Bond centralBond : bonds) {//identify cases of E/Z stereochemistry and add appropriate bondstereo tags
			if (centralBond.getOrder() == 2) {
				List<Bond> fromAtomBonds = centralBond.getFromAtom().getBonds();
				for (Bond preceedingBond : fromAtomBonds) {
					if (preceedingBond.getSmilesStereochemistry() != null) {
						List<Bond> toAtomBonds = centralBond.getToAtom().getBonds();
						for (Bond followingBond : toAtomBonds) {
							if (followingBond.getSmilesStereochemistry() != null) {//now found a double bond surrounded by two bonds with slashs
								boolean upFirst;
								boolean upSecond;
								Atom atom2 = centralBond.getFromAtom();
								Atom atom3 = centralBond.getToAtom();
								Atom atom1 = preceedingBond.getOtherAtom(atom2);
								Atom atom4 = followingBond.getOtherAtom(atom3);
								if (preceedingBond.getSmilesStereochemistry() == SMILES_BOND_DIRECTION.LSLASH) {
									upFirst = preceedingBond.getToAtom() == atom2;//in normally constructed SMILES this will be the case but you could write C(/F)=C/F instead of F\C=C/F
								}
								else if (preceedingBond.getSmilesStereochemistry() == SMILES_BOND_DIRECTION.RSLASH) {
									upFirst = preceedingBond.getToAtom() != atom2;
								}
								else{
									throw new StructureBuildingException(preceedingBond.getSmilesStereochemistry() + " is not a slash!");
								}

								if (followingBond.getSmilesStereochemistry() == SMILES_BOND_DIRECTION.LSLASH) {
									upSecond = followingBond.getFromAtom() != atom3;
								}
								else if (followingBond.getSmilesStereochemistry() == SMILES_BOND_DIRECTION.RSLASH) {
									upSecond = followingBond.getFromAtom() == atom3;
								}
								else{
									throw new StructureBuildingException(followingBond.getSmilesStereochemistry() + " is not a slash!");
								}
								BondStereoValue cisTrans = upFirst == upSecond ? BondStereoValue.CIS : BondStereoValue.TRANS;
								if (centralBond.getBondStereo() != null) {
									//double bond has redundant specification e.g. C/C=C\\1/NC1 hence need to check it is consistent
									Atom[] atomRefs4 = centralBond.getBondStereo().getAtomRefs4();
									if (atomRefs4[0].equals(atom1) || atomRefs4[3].equals(atom4)) {
										if (centralBond.getBondStereo().getBondStereoValue().equals(cisTrans)){
											throw new StructureBuildingException("Contradictory double bond stereoconfiguration");
										}
									}
									else{
										if (!centralBond.getBondStereo().getBondStereoValue().equals(cisTrans)){
											throw new StructureBuildingException("Contradictory double bond stereoconfiguration");
										}
									}
								}
								else{
									Atom[] atomRefs4= new Atom[4];
									atomRefs4[0] = atom1;
									atomRefs4[1] = atom2;
									atomRefs4[2] = atom3;
									atomRefs4[3] = atom4;
									centralBond.setBondStereoElement(atomRefs4, cisTrans);
								}
							}
						}
					}
				}
			}
		}
		for (Bond bond : bonds) {
			bond.setSmilesStereochemistry(null);
		}
	}
	
	/**
	 * Utilises the atom's hydrogen count as set by the SMILES as well as incoming valency to determine the atom's valency
	 * If the atom is charged whether protons have been added or removed will also need to be determined
	 * @param atom
	 * @throws StructureBuildingException 
	 */
	private void setupAtomValency(Atom atom) throws StructureBuildingException {
		int hydrogenCount = atom.getProperty(Atom.SMILES_HYDROGEN_COUNT);
		int incomingValency = atom.getIncomingValency() + hydrogenCount +atom.getOutValency();
		int charge = atom.getCharge();
		int absoluteCharge =Math.abs(charge);
		ChemEl chemEl = atom.getElement();
		if (atom.hasSpareValency()) {
			Integer hwValency = ValencyChecker.getHWValency(chemEl);
			if (hwValency == null || absoluteCharge > 1) {
				throw new StructureBuildingException(chemEl +" is not expected to be aromatic!");
			}
			if (absoluteCharge != 0) {
				Integer[] possibleVal = ValencyChecker.getPossibleValencies(chemEl, charge);
				if (possibleVal != null && possibleVal.length > 0) {
					hwValency = possibleVal[0];
				}
				else {
					throw new StructureBuildingException(chemEl +" with charge " + charge + " is not expected to be aromatic!");
				}
			}
			if (incomingValency < hwValency){
				incomingValency++;
			}
		}
		Integer defaultVal = ValencyChecker.getDefaultValency(chemEl);
		if (defaultVal !=null){//s or p block element
			if (defaultVal != incomingValency || charge !=0) {
				if (Math.abs(incomingValency - defaultVal) == absoluteCharge) {
					atom.setProtonsExplicitlyAddedOrRemoved(incomingValency - defaultVal);
				}
				else{
					Integer[] unchargedStableValencies = ValencyChecker.getPossibleValencies(chemEl, 0);
					boolean hasPlausibleValency =false;
					for (Integer unchargedStableValency : unchargedStableValencies) {
						if (Math.abs(incomingValency - unchargedStableValency)==Math.abs(charge)){
							atom.setProtonsExplicitlyAddedOrRemoved(incomingValency - unchargedStableValency);
							//we strictly set the valency if a charge is specified but are more loose about things if uncharged e.g. allow penta substituted phosphine
							if (charge != 0) {
								atom.setLambdaConventionValency(unchargedStableValency);
							}
							else{
								atom.setMinimumValency(incomingValency);
							}
							hasPlausibleValency=true;
							break;
						}
					}
					if (!hasPlausibleValency){//could be something like [Sn] which would be expected to be attached to later
						atom.setMinimumValency(incomingValency);
					}
				}
			}
		}
		else{
			if (hydrogenCount > 0){//make hydrogen explicit
				Fragment frag =atom.getFrag();
				for (int i = 0; i < hydrogenCount; i++) {
					Atom hydrogen = createAtom(ChemEl.H, frag);
					createBond(atom, hydrogen, 1);
				}
			}
		}
	}
	
	
	/**
	 * Create a new Atom of the given element belonging to the given fragment
	 * @param elementSymbol
	 * @param frag
	 * @return Atom
	 */
	private Atom createAtom(String elementSymbol, Fragment frag) {
		return createAtom(ChemEl.valueOf(elementSymbol), frag);
	}

	/**
	 * Create a new Atom of the given element belonging to the given fragment
	 * @param chemEl
	 * @param frag
	 * @return Atom
	 */
	private Atom createAtom(ChemEl chemEl, Fragment frag) {
		Atom a = new Atom(idManager.getNextID(), chemEl, frag);
		frag.addAtom(a);
		return a;
	}
	
	/**
	 * Create a new bond between two atoms.
	 * The bond is associated with these atoms.
	 * @param fromAtom
	 * @param toAtom
	 * @param bondOrder
	 * @return Bond
	 */
	private Bond createBond(Atom fromAtom, Atom toAtom, int bondOrder) {
		Bond b = new Bond(fromAtom, toAtom, bondOrder);
		fromAtom.addBond(b);
		toAtom.addBond(b);
		fromAtom.getFrag().addBond(b);
		return b;
	}
	
	private boolean is_A_to_Z(char ch) {
		return ch >= 'A' && ch <= 'Z';
	}
	
	private boolean is_a_to_z(char ch) {
		return ch >= 'a' && ch <= 'z';
	}
	
	private boolean is_0_to_9(char ch){
		return ch >= '0' && ch <= '9'; 
	}

}
