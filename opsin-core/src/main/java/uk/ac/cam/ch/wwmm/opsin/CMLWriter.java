package uk.ac.cam.ch.wwmm.opsin;

import java.io.ByteArrayOutputStream;
import java.io.UnsupportedEncodingException;
import java.util.List;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import org.codehaus.stax2.util.StreamWriterDelegate;

import com.ctc.wstx.stax.WstxOutputFactory;

public class CMLWriter {
	/**
	 * CML Elements/Attributes/NameSpace
	 */
	private static final String CML_NAMESPACE = "http://www.xml-cml.org/schema";

	private static final XMLOutputFactory factory = new WstxOutputFactory();
	
	/**The structure to be converted to CML*/
	private final Fragment structure;
	
	/**The name of the structure*/
	private final String chemicalName;
	
	/**The XML writer*/
	private final XMLStreamWriter writer;

	/**
	 * Creates a CML writer for the given fragment
	 * @param writer 
	 * @param structure
	 * @param chemicalName
	 */
	private CMLWriter(XMLStreamWriter writer, Fragment structure, String chemicalName) {
		this.structure = structure;
		this.chemicalName = chemicalName;
		this.writer = writer;
	}
	
	static String generateCml(Fragment structure, String chemicalName) {
		ByteArrayOutputStream out = new ByteArrayOutputStream();
		try {
			XMLStreamWriter writer = factory.createXMLStreamWriter(out, "UTF-8");
			new CMLWriter(writer, structure, chemicalName).writeCml();
		} catch (XMLStreamException e) {
			throw new RuntimeException(e);
		}
		try {
			return out.toString("UTF-8");
		} catch (UnsupportedEncodingException e) {
			throw new RuntimeException("JVM doesn't support UTF-8...but it should do!");
		}
	}
	
	static String generateCml(Fragment structure, String chemicalName, int indent){
		ByteArrayOutputStream out = new ByteArrayOutputStream();
		try {
			XMLStreamWriter writer = factory.createXMLStreamWriter(out, "UTF-8");
			writer = new IndentingXMLStreamWriter(writer, indent);
			new CMLWriter(writer, structure, chemicalName).writeCml();
		} catch (XMLStreamException e) {
			throw new RuntimeException(e);
		}
		try {
			return out.toString("UTF-8");
		} catch (UnsupportedEncodingException e) {
			throw new RuntimeException("JVM doesn't support UTF-8...but it should do!");
		}
	}

	private void writeCml(){
		try {
			writer.writeStartDocument();
			writer.writeStartElement("cml");
			writer.writeDefaultNamespace(CML_NAMESPACE);
			writer.writeAttribute("convention", "conventions:molecular");
			writer.writeNamespace("conventions", "http://www.xml-cml.org/convention/");
			writer.writeNamespace("cmlDict", "http://www.xml-cml.org/dictionary/cml/");
			writer.writeNamespace("nameDict", "http://www.xml-cml.org/dictionary/cml/name/");
			
			writer.writeStartElement("molecule");
			writeMolecule();
			writer.writeEndElement();
			
			writer.writeEndElement();
			writer.writeEndDocument();
			writer.flush();
			writer.close();
		} catch (XMLStreamException e) {
			throw new RuntimeException(e);
		}
	}

	private void writeMolecule() throws XMLStreamException {
		writer.writeAttribute("id", "m1");
		writer.writeStartElement("name");
		writer.writeAttribute("dictRef", "nameDict:unknown");
		writer.writeCharacters(chemicalName);
		writer.writeEndElement();
		
		writer.writeStartElement("atomArray");
		for(Atom atom : structure.getAtomList()) {
			writeAtom(atom);
		}
		writer.writeEndElement();
		
		writer.writeStartElement("bondArray");
		for(Bond bond : structure.getBondSet()) {
			writeBond(bond);
		}
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
		if(atomParity != null){
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
	
	/**
	 * This only overrides the commands actually used by the CmlWriter i.e. it isn't general
	 */
	private static class IndentingXMLStreamWriter extends StreamWriterDelegate {
		private final int indentSize;
		private int depth = 0;
		private boolean atStartOfNewline = false;
				
		IndentingXMLStreamWriter(XMLStreamWriter writer, int indentSize) {
			super(writer);
			this.indentSize = indentSize;
		}

		@Override
		public void writeStartElement(String arg0) throws XMLStreamException {
			if (!atStartOfNewline){
				super.writeCharacters("\n");
			}
			super.writeCharacters(StringTools.multiplyString(" ", depth * indentSize));
			super.writeStartElement(arg0);
			atStartOfNewline = false;
			depth++;
		}
		
		@Override
		public void writeEndElement() throws XMLStreamException {
			depth--;
			if (atStartOfNewline) {
				super.writeCharacters(StringTools.multiplyString(" ", depth * indentSize));
			}
			super.writeEndElement();
			super.writeCharacters("\n");
			atStartOfNewline = true;
		}
		
		@Override
		public void writeCharacters(String arg0) throws XMLStreamException {
			super.writeCharacters(arg0);
			atStartOfNewline = false;
		}
		
	}
}
