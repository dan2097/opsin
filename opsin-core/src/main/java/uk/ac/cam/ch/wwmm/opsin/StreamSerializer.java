package uk.ac.cam.ch.wwmm.opsin;

import java.io.IOException;
import java.io.OutputStream;

import nu.xom.Element;
import nu.xom.Serializer;

public class StreamSerializer extends Serializer {

    public StreamSerializer(OutputStream out) {
        super(out);
    }

    @Override
    public void write(Element element) throws IOException {
        super.write(element);
    }

    @Override
    public void writeXMLDeclaration() throws IOException {
        super.writeXMLDeclaration();
    }

    @Override
    public void writeEndTag(Element element) throws IOException {
        super.writeEndTag(element);
    }

    @Override
    public void writeStartTag(Element element) throws IOException {
        super.writeStartTag(element);
    }
}

