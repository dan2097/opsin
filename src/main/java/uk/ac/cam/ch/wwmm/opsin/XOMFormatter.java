package uk.ac.cam.ch.wwmm.opsin;

import java.io.ByteArrayOutputStream;
import java.io.IOException;

import nu.xom.Document;
import nu.xom.Element;
import nu.xom.Serializer;

/**Turns a XOM Element into a pretty indented string.
 * 
 * @author ptc24
 *
 */
public class XOMFormatter {

	ByteArrayOutputStream outStream = new ByteArrayOutputStream();
	Serializer serializer;	

	/**Sets up a new XOMFormatter.
	 * 
	 */
	public XOMFormatter() {
		super();
		try {
			serializer = new Serializer(outStream, "ISO-8859-1");
			serializer.setIndent(4);
			serializer.setMaxLength(300);
		    }
		    catch (IOException ex) {
		       System.err.println(ex); 
		    }
	}

	/**Converts an Element to an indented string.
	 * 
	 * @param elem The Element to convert to a string.
	 * @return The string.
	 */
	public String elemToString(Element elem) {
		try {
			// Grrr protected methods grrr
			outStream.reset();
			// Put the element in a document...
			serializer.write(new Document(new Element(elem)));
			// Then return the document, destroying the evidence
			// that it ever was a document.
			return outStream.toString().substring(45);
		} catch (IOException ex) {
			ex.printStackTrace();
		}
		return null;
	}

}
