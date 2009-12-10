package uk.ac.cam.ch.wwmm.opsin;

import java.net.URI;
import java.net.URISyntaxException;
import java.net.URL;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import nu.xom.Document;
import nu.xom.Elements;

import org.junit.Test;
import org.xml.sax.ErrorHandler;
import org.xml.sax.SAXException;
import org.xml.sax.SAXParseException;

public class DtdTest {
	private final String RESOURCE_LOCATION = "uk/ac/cam/ch/wwmm/opsin/resources/";
	private final ResourceGetter resourceGetter = new ResourceGetter(RESOURCE_LOCATION);

	@Test
	public void testTokenFiles() throws Exception {
		Document tokenFileDoc = resourceGetter.getXMLDocument("index.xml");
		Elements tokenFiles = tokenFileDoc.getRootElement().getChildElements();
		for (int i = 0; i < tokenFiles.size(); i++) {
			validate(getUriForFile(tokenFiles.get(i).getValue()));
		}
	}
	
	@Test
	public void testRegexes() throws Exception {
		validate(getUriForFile("regexes.xml"));
	}
	
	@Test
	public void testRegexTokens() throws Exception {
		validate(getUriForFile("regexTokens.xml"));
	}
	
	@Test
	public void testSuffixApplicability() throws Exception {
		validate(getUriForFile("suffixApplicability.xml"));
	}
	
	@Test
	public void testSuffixRules() throws Exception {
		validate(getUriForFile("suffixRules.xml"));
	}
	
	@Test
	public void testWordRules() throws Exception {
		validate(getUriForFile("wordRules.xml"));
	}
	
	public static void validate(URI uri) throws Exception {
		System.out.println("Validating:"+ uri);
		DocumentBuilderFactory f = DocumentBuilderFactory.newInstance();
		f.setValidating(true);
		DocumentBuilder b = f.newDocumentBuilder();
		MyErrorHandler h = new MyErrorHandler();
		b.setErrorHandler(h);
		try {
			b.parse(uri.toString());
		} catch (SAXException e) {
			if (h.error != null) {
				AssertionError ae = new AssertionError("XML Validation error: "+uri.toString());
				ae.initCause(h.error);
				throw ae;
			}
		}
	}

	static class MyErrorHandler implements ErrorHandler {

		private SAXParseException error;
		
		public void error(SAXParseException exception) throws SAXException {
			this.error = exception;
			throw new SAXException("Error");
		}

		public void fatalError(SAXParseException exception) throws SAXException {
			this.error = exception;
			throw new SAXException("Error");
		}

		public void warning(SAXParseException exception) throws SAXException {
			this.error = exception;
			throw new SAXException("Error");
		}
		
	}
	
	private URI getUriForFile (String fileName) throws URISyntaxException {
		ClassLoader l = getClass().getClassLoader();
		URL url = l.getResource(RESOURCE_LOCATION + fileName);
		if (url ==null) {return null;}
		return url.toURI();
	}

}
