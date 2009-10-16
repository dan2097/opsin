package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;

class Ring {
	private ArrayList<Atom>  atomSet;
	private ArrayList<Bond>  bondSet;
	private ArrayList<Atom>  cyclicAtomSet;
	private ArrayList<Bond>  cyclicBondSet;
	private ArrayList<Ring>  neighbors = new ArrayList<Ring>();
	
	private int nFusedBonds = 0;
	
	private int size = -1;
	
	public Ring(ArrayList<Atom> atomSet, ArrayList<Bond> bondSet) throws StructureBuildingException {
		this.atomSet = atomSet;
		this.bondSet = bondSet;
		if (atomSet.size() != bondSet.size()) {
			throw new StructureBuildingException("atomSet and bondset different sizes");
		}
	}
	
	public Ring(ArrayList<Bond> bondSet) throws StructureBuildingException {
		if (bondSet==null || bondSet.size()<=0) throw new StructureBuildingException("Bond set is empty");
		this.bondSet = bondSet;			
		this.atomSet = new ArrayList<Atom>();
		
		for(Bond bond: bondSet)
		{
			Atom atom1 = bond.getFromAtom();				
			if (!atomSet.contains(atom1)) atomSet.add(atom1);
				
			Atom atom2 = bond.getToAtom();
			if (!atomSet.contains(atom2)) atomSet.add(atom2);
			
		}
		
		if (atomSet.size() != bondSet.size()) {
			throw new StructureBuildingException("atomSet and bondset different sizes. Ring(bond)");
		}
		
	}
	
	public ArrayList<Bond>  getBondSet() {
		return bondSet;
	}
	
	public ArrayList<Atom>  getAtomSet() {
		return atomSet;
	}
	
	public int size() {
		if (size<0)	size = bondSet.size();
		return size;
	}
	
	public int getNOFusedBonds() {
		return nFusedBonds;
	}
	
	public void setNOFusedBonds(int fusedBonds) {
		this.nFusedBonds = fusedBonds;
	}
	
	public void incrNOFFusedBonds() {
		this.nFusedBonds++;
	}
	
	public ArrayList<Bond> getFusedBonds()
	{
		ArrayList<Bond> bonds = new ArrayList<Bond>();
		
		for (Bond bond : bondSet) {
			if (bond.getFusedRings().size()>0) bonds.add(bond);
		}
		
		return bonds;
	}
	
	public int getBondNumber(Bond bond)
	{
		int i=0;
		for(Bond ibond : cyclicBondSet)
		{			
			if (ibond == bond) return i;
			i++;
		}
		return -1;
	}
	
	public ArrayList<Bond> getCyclicBondSet() 
	{
		return cyclicBondSet;
	}
	
	public ArrayList<Atom> getCyclicAtomSet() 
	{
		return cyclicAtomSet;
	}
	
	public ArrayList<Ring> getNeighbors() {
		return neighbors;
	}
	
	public void addNeighbor(Ring ring) {
		this.neighbors.add(ring);
	}
	/**
	 * Stores atoms and bonds in the order defined by atom and bond
	 * @param stBond - the first bond
	 * @param stAtom - the atom defining in which direction to go
	 */	
	public void makeCyclicSets(Bond stBond, Atom stAtom)
	{
		if (cyclicBondSet==null){
			cyclicBondSet = new ArrayList<Bond>();
			cyclicAtomSet = new ArrayList<Atom>();
			
			Bond bond = stBond;
			Atom atom = stAtom;
			cyclicBondSet.add(bond);
			cyclicAtomSet.add(atom);
			
			for (int i=0; i<size()-1; i++)
			{								
				for(Bond bond2 : bondSet)
				{
					if (cyclicBondSet.contains(bond2)) continue;					
					if (bond2.getFromAtom() == atom) 
					{
						cyclicBondSet.add(bond2);						
						atom = bond2.getToAtom();
						cyclicAtomSet.add(atom);
					}
					else if (bond2.getToAtom() == atom)
					{
						cyclicBondSet.add(bond2); 
						atom = bond2.getFromAtom();
						cyclicAtomSet.add(atom);
					}
				}
			}
		}		
	}
	
	public String toString()
	{
		String st = "";
		ArrayList<Atom> toPrint;
		
		if (cyclicAtomSet != null) toPrint = cyclicAtomSet;
		else toPrint = atomSet;
		for (Atom atom : toPrint) 
		{
			st+=atom.getID() + " ";
		}
		return st;
		//return size + "";
	}
		
}
