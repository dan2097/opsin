package uk.ac.cam.ch.wwmm.opsin;

import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import org.codehaus.stax2.util.StreamWriterDelegate;

/**
 * This only overrides the commands actually used by the CmlWriter i.e. it isn't general
 */
class IndentingXMLStreamWriter extends StreamWriterDelegate {

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
			super.writeCharacters(OpsinTools.NEWLINE);
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
		super.writeCharacters(OpsinTools.NEWLINE);
		atStartOfNewline = true;
	}
	
	@Override
	public void writeCharacters(String arg0) throws XMLStreamException {
		super.writeCharacters(arg0);
		atStartOfNewline = false;
	}
	
}
