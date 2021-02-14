package uk.ac.cam.ch.wwmm.opsin;

import java.io.ByteArrayOutputStream;
import java.io.UnsupportedEncodingException;
import java.util.List;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.ctc.wstx.api.WstxOutputProperties;
import com.ctc.wstx.stax.WstxOutputFactory;

class CMLWriter {
	/**
	 * CML Elements/Attributes/NameSpace
	 */
	static final String CML_NAMESPACE = "http://www.xml-cml.org/schema";

	private static final XMLOutputFactory factory = new WstxOutputFactory();
	static {
		factory.setProperty(WstxOutputProperties.P_OUTPUT_ESCAPE_CR, false);
	}
	
	/**The XML writer*/
	private final XMLStreamWriter writer;

	/**
	 * Creates a CML writer for the given fragment
	 * @param writer 

	 */
	CMLWriter(XMLStreamWriter writer) {
		this.writer = writer;
	}
	
	static String generateCml(Fragment structure, String chemicalName) {
		return generateCml(structure, chemicalName, false);
	}
	
	static String generateIndentedCml(Fragment structure, String chemicalName) {
		return generateCml(structure, chemicalName, true);
	}
	
	private static String generateCml(Fragment structure, String chemicalName, boolean indent) {
		ByteArrayOutputStream out = new ByteArrayOutputStream();
		try {
			XMLStreamWriter xmlWriter = factory.createXMLStreamWriter(out, "UTF-8");
			if (indent) {
				xmlWriter = new IndentingXMLStreamWriter(xmlWriter, 2);
			}
			CMLWriter cmlWriter = new CMLWriter(xmlWriter);
			cmlWriter.writeCmlStart();
			cmlWriter.writeMolecule(structure, chemicalName, 1);
			cmlWriter.writeCmlEnd();
			xmlWriter.close();
		} catch (XMLStreamException e) {
			throw new RuntimeException(e);
		}
		try {
			return out.toString("UTF-8");
		} catch (UnsupportedEncodingException e) {
			throw new RuntimeException("JVM doesn't support UTF-8...but it should do!");
		}
	}

	void writeCmlStart(){
		try {
			writer.writeStartElement("cml");
			writer.writeDefaultNamespace(CML_NAMESPACE);
			writer.writeAttribute("convention", "conventions:molecular");
			writer.writeNamespace("conventions", "http://www.xml-cml.org/convention/");
			writer.writeNamespace("cmlDict", "http://www.xml-cml.org/dictionary/cml/");
			writer.writeNamespace("nameDict", "http://www.xml-cml.org/dictionary/cml/name/");
		} catch (XMLStreamException e) {
			throw new RuntimeException(e);
		}
	}

	void writeCmlEnd(){
		try {		
			writer.writeEndElement();
			writer.flush();
		} catch (XMLStreamException e) {
			throw new RuntimeException(e);
		}
	}

	void writeMolecule(Fragment structure, String chemicalName, int id) throws XMLStreamException {
		writer.writeStartElement("molecule");
		writer.writeAttribute("id", "m" + id);

		writer.writeStartElement("name");
		writer.writeAttribute("dictRef", "nameDict:unknown");
		writer.writeCharacters(chemicalName);
		writer.writeEndElement();
		
		if (structure != null) {
			writer.writeStartElement("atomArray");
			for(Atom atom : structure.getAtomList()) {
				writeAtom(atom);
			}
			writer.writeEndElement();
			
			writer.writeStartElement("bondArray");
			for(Bond bond : structure.getBondSet()) {
				writeBond(bond);
			}
			writer.writeEndElement();
		}
		
		writer.writeEndElement();
	}
	
	private void writeAtom(Atom atom) throws XMLStreamException {
		writer.writeStartElement("atom");
		writer.writeAttribute("id", "a" + Integer.toString(atom.getID()));
		writer.writeAttribute("elementType", atom.getElement().toString());
		if(atom.getCharge() != 0){
			writer.writeAttribute("formalCharge", Integer.toString(atom.getCharge()));
		}
		if(atom.getIsotope() != null){
			writer.writeAttribute("isotopeNumber", Integer.toString(atom.getIsotope()));
		}
		if (atom.getElement() != ChemEl.H){
			int hydrogenCount =0;
			List<Atom> neighbours = atom.getAtomNeighbours();
			for (Atom neighbour : neighbours) {
				if (neighbour.getElement() == ChemEl.H){
					hydrogenCount++;
				}
			}
			if (hydrogenCount==0){//prevent adding of implicit hydrogen
				writer.writeAttribute("hydrogenCount", "0");
			}
		}
		AtomParity atomParity = atom.getAtomParity();
		if(atomParity != null && atomParity.getStereoGroup() != StereoGroup.Rac){
			writeAtomParity(atomParity);
		}
		for(String locant : atom.getLocants()) {
			writer.writeStartElement("label");
			writer.writeAttribute("value", locant);
			writer.writeAttribute("dictRef", "cmlDict:locant");
			writer.writeEndElement();
		}
		writer.writeEndElement();
	}

	private void writeAtomParity(AtomParity atomParity) throws XMLStreamException {
		writer.writeStartElement("atomParity");
		writeAtomRefs4(atomParity.getAtomRefs4());
		writer.writeCharacters(Integer.toString(atomParity.getParity()));
		writer.writeEndElement();
	}

	private void writeBond(Bond bond) throws XMLStreamException {
		writer.writeStartElement("bond");
		writer.writeAttribute("id", "a" + Integer.toString(bond.getFrom()) + "_a" + Integer.toString(bond.getTo()));
		writer.writeAttribute("atomRefs2", "a" + Integer.toString(bond.getFrom()) + " a" + Integer.toString(bond.getTo()));
		switch (bond.getOrder()) {
		case 1:
			writer.writeAttribute("order", "S");
			break;
		case 2:
			writer.writeAttribute("order", "D");
			break;
		case 3:
			writer.writeAttribute("order", "T");
			break;
		default:
			writer.writeAttribute("order", "unknown");
			break;
		}
		BondStereo bondStereo = bond.getBondStereo();
		if (bondStereo != null){
			writeBondStereo(bondStereo);
		}
		writer.writeEndElement();
	}

	private void writeBondStereo(BondStereo bondStereo) throws XMLStreamException {
		writer.writeStartElement("bondStereo");
		writeAtomRefs4(bondStereo.getAtomRefs4());
		writer.writeCharacters(bondStereo.getBondStereoValue().toString());
		writer.writeEndElement();
	}

	private void writeAtomRefs4(Atom[] atomRefs4) throws XMLStreamException {
		StringBuilder atomRefsSb = new StringBuilder();
		for(int i = 0; i< atomRefs4.length - 1; i++) {
			atomRefsSb.append('a');
			atomRefsSb.append(atomRefs4[i].getID());
			atomRefsSb.append(' ');
		}
		atomRefsSb.append('a');
		atomRefsSb.append(atomRefs4[atomRefs4.length - 1].getID());
		writer.writeAttribute("atomRefs4", atomRefsSb.toString());
	}

}
