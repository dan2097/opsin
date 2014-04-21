package uk.ac.cam.ch.wwmm.opsin;

import java.io.ByteArrayOutputStream;
import java.io.IOException;

import nu.xom.Element;
import nu.xom.Serializer;
/**
 * Turns a XOM Element into a pretty indented string.
 * @author ptc24
 * @author dl387
 *
 */
public class XOMFormatter extends Serializer {

	private final ByteArrayOutputStream outStream;

	/**
	 * Sets up a new XOMFormatter.
	 */
	public XOMFormatter() {
		this(new ByteArrayOutputStream());
	}
	
	private XOMFormatter(ByteArrayOutputStream outStream) {
		super(outStream);
		this.outStream = outStream;
		setIndent(4);
		setMaxLength(300);
	}

	/**Converts an Element to an indented string.
	 *
	 * @param elem The Element to convert to a string.
	 * @return The string.
	 */
	public String elemToString(Element elem) {
		try {
			outStream.reset();
			write(elem);
			flush();
			return outStream.toString();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return null;
	}
}