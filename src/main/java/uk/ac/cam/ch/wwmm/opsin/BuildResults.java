package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;


/**A "struct" to hold the results of fragment building.*/
class BuildResults {
	/**Holds the ID or IDs of the atoms that will be used to connect to other fragments.
	 * Also holds the order of the bond that will be created when it is used (valency) e.g.  Eg. chloro 1, oxo 2
	 * setExplicitly says whether the outID absolutely definitely refers to that id or not.
	 * e.g. propyl is stored as prop-1-yl with this set to false while prop-2-yl has it set to true
	 * If at the end of building any IDs are unused R groups will be added to indicate radicals
	 * These OutIDs are the same objects as are present in the fragments*/
	private LinkedList<OutID> outIDs;
	
	/**The ID or IDs of the atoms that may be used to from things like esters*/
	LinkedList<Integer> functionalIDs;
	
	/**A list of fragments that are represented by this class. The first fragment in this list is the "main Fragment"*/
	LinkedHashSet<Fragment> fragments;

	/**A BuildResults is constructed from a given fragment.*/
	BuildResults(Fragment mainFrag) {
		outIDs = new LinkedList<OutID>();
		List<OutID> givenOutIDs =mainFrag.getOutIDs();
		for (OutID outID : givenOutIDs) {
			outIDs.add(outID);
		}
		functionalIDs = new LinkedList<Integer>();
		List<Integer> givenOutFunctionalIDs =mainFrag.getFunctionalIDs();
		for (Integer ID : givenOutFunctionalIDs) {
			functionalIDs.add(new Integer(ID));
		}
		fragments=new LinkedHashSet<Fragment>();
		fragments.add(mainFrag);
	}

	/**
	 * Returns the fragment which this buildResults was originally built around
	 * @return
	 * @throws StructureBuildingException 
	 */
	Fragment getMainFragment() throws StructureBuildingException{
		for (Fragment fragment : fragments) {
			return fragment;
		}
		throw new StructureBuildingException("No fragments in BuildResults");
	}
	
	/**
	 * Returns the atom corresponding to position i in the outID list
	 * @param i index
	 * @return atom
	 * @throws StructureBuildingException 
	 */
	Atom getOutAtom(int i) throws StructureBuildingException {
		int outID =outIDs.get(i).id;
		for (Fragment fragment : fragments) {
			Atom outAtom =fragment.getAtomByID(outID);
			if (outAtom!=null){
				return outAtom;
			}
		}
		throw new StructureBuildingException("No fragment contained this outID: " + outID);
	}
	
	OutID getOutID(int i){
		return outIDs.get(i);
	}
	
	int getOutIDCount(){
		return outIDs.size();
	}
	
	void removeOutID(int i) throws StructureBuildingException{
		OutID outID =outIDs.get(i);
		if (outID.frag!=null){
			outID.frag.removeOutID(outID);
		}
		outIDs.remove(i);
	}
	
	void removeAllOutIDs() throws StructureBuildingException{
		for (int i = outIDs.size() -1; i >=0 ; i--) {
			removeOutID(i);
		}
	}
	
	/**
	 * Returns the atom corresponding to position i in the outFunctionalIDs list
	 * @param i index
	 * @return atom
	 * @throws StructureBuildingException 
	 */
	Atom getFunctionalOutAtom(int i) throws StructureBuildingException {
		int functionalOutID =functionalIDs.get(i);
		for (Fragment fragment : fragments) {
			Atom functionalOutAtom =fragment.getAtomByID(functionalOutID);
			if (functionalOutAtom!=null){
				return functionalOutAtom;
			}
		}
		throw new StructureBuildingException("No fragment contained this functionalOutID: " + functionalOutID);
	}
	
	/**
	 * Returns the first OutId
	 * @return atom
	 */
	OutID getFirstOutID(){
		return outIDs.get(0);
	}
	
	/**
	 * Adds the given outID and associates it with the order of bond it should form and whether it was set explicitly
	 * @param id
	 * @param valency
	 * @param setExplicitly
	 */
	void addOutID(int id, int valency, boolean setExplicitly) {
		outIDs.add(new OutID(id, valency, setExplicitly, null, null));
	}

	/**
	 * Combines a new buildResults into this buildResults.
	 * OutIDs
	 * @param resolveWordOrBracket
	 */
	void mergeBuildResults(BuildResults newBuildResults) {
		fragments.addAll(newBuildResults.fragments);
		outIDs.addAll(newBuildResults.outIDs);
		functionalIDs.addAll(newBuildResults.functionalIDs);
	}

	/**
	 * Returns the atom corresponding to position i in the outFunctionalIDs list
	 * If not set explicitly and atom would violate valency or break aromaticity another is looked for
	 * @param i index
	 * @return atom
	 * @throws StructureBuildingException 
	 */
	Atom getOutAtomTakingIntoAccountWhetherSetExplicitly(int i) throws StructureBuildingException {
		int outID =outIDs.get(i).id;
		if (outIDs.get(i).setExplicitly==true){
			for (Fragment fragment : fragments) {
				Atom outAtom =fragment.getAtomByID(outID);
				if (outAtom!=null){
					return outAtom;
				}
			}
		}
		else{
			for (Fragment fragment : fragments) {
				Atom outAtom =fragment.getAtomByID(outID);
				if (outAtom!=null){
					return fragment.getAtomByIdOrNextSuitableAtomOrThrow(outID, outIDs.get(i).valency);
				}
			}
		}
		throw new StructureBuildingException("No fragment contained this outID: " + outID);
	}
	
	/**
	 * Returns the atom corresponding to the given id assuming the atom the id corresponds to is within the list of fragment in this Buildresults
	 * @param i index
	 * @return atom
	 * @throws StructureBuildingException 
	 */
	Atom getAtomById(int id) throws StructureBuildingException {
		for (Fragment fragment : fragments) {
			Atom outAtom =fragment.getAtomByID(id);
			if (outAtom!=null){
				return outAtom;
			}
		}
		throw new StructureBuildingException("No fragment contained this id: " + id);
	}

	/**
	 * Performs a subStructure search. The input is SMILES. SMARTS queries are not currently supported.
	 * Returns a list of atomLists.
	 * There is one atomlist for each successful match. The atomlist contains the atoms matched. The order is the same as that specified in the SMILES query.
	 * @param String SMILES
	 * @param FragmentManage fm
	 * @throws StructureBuildingException 
	 */
	List<List<Atom>> subStructureSearch(String smiles, FragmentManager fm) throws StructureBuildingException {
		if (smiles.contains(".")){
			throw new StructureBuildingException("Substructure SMILES may not contain disconnections");
		}
		if (smiles.contains("(")){
			throw new StructureBuildingException("Substructure SMILES may not contain branching");
		}
		List<List<Atom>> results =new ArrayList<List<Atom>>();
		Fragment query =fm.buildSMILES(smiles);
		fm.removeFragment(query);
		for (Fragment frag : fragments) {
			List<Atom> atomList =frag.getAtomList();
			for (Atom atom : atomList) {
				List<List<Atom>> resultsWithThisStartingAtom = recursivelyMatch(query, 0, atom, new ArrayList<Atom>(), new ArrayList<List<Atom>>());
				if (resultsWithThisStartingAtom.size() >0){
					 results.addAll(resultsWithThisStartingAtom);
				}
			}
		}
		return results;
	}

	private boolean checkAtomsMatch(Atom a1, Atom a2){
		if (!a1.getElement().equals(a2.getElement())){
			return false;
		}
		if (a1.getCharge() !=a2.getCharge()){
			return false;
		}
		
		return true;
	}
	
	/**
	 * Recurisively substructure matches
	 * @param query he query fragment generated from the original query smiles
	 * @param i The atom position in the query
	 * @param atom The atomToGoToNext
	 * @param result An atomlist containing the results so far from this instance of the method
	 * @param results A list of atomlist containing succesful substructure matches found so far
	 * @return results 
	 * @throws StructureBuildingException 
	 */
	private List<List<Atom>> recursivelyMatch(Fragment query, int i, Atom atom,
			List<Atom> result, ArrayList<List<Atom>> results) throws StructureBuildingException {
		Atom queryAtom =query.getAtomList().get(i);
		if (checkAtomsMatch(atom, queryAtom)){
			result.add(atom);
			if (i==query.getAtomList().size()-1){
				results.add(result);
				return results;
			}
			
			List<Bond> bondsToNeighboursInQuery =queryAtom.getBonds();
			int expectedBondOrder =0;
			for(Bond b :  bondsToNeighboursInQuery) {
				if(b.getFrom() ==queryAtom.getID() +1 || b.getTo() ==queryAtom.getID() +1) {
					expectedBondOrder=b.getOrder();
				}
			}
			
			List<Bond> bondsToNeighbours =atom.getBonds();

			int ID =atom.getID();
			for(Bond b :  bondsToNeighbours) {
				if (b.getOrder()==expectedBondOrder){
					Integer toID =null;
					if(b.getFrom() != ID) {
						toID =b.getFrom();
					} else if(b.getTo() != ID) {
						toID =b.getTo();
					}
					recursivelyMatch(query, i +1, atom.getFrag().getAtomByIDOrThrow(toID), new ArrayList<Atom>(result), results);
				}
			}
		}
		return results;
	}
}
