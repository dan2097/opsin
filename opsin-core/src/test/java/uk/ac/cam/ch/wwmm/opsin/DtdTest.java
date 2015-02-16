package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.Assert.assertTrue;

import java.net.URI;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.HashSet;
import java.util.Set;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.stream.XMLStreamConstants;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;

import org.junit.Test;
import org.xml.sax.ErrorHandler;
import org.xml.sax.SAXException;
import org.xml.sax.SAXParseException;

public class DtdTest {
	private final static String RESOURCE_LOCATION = "uk/ac/cam/ch/wwmm/opsin/resources/";
	private final ResourceGetter resourceGetter = new ResourceGetter(RESOURCE_LOCATION);

	@Test
	public void testTokenFiles() throws Exception {
		XMLStreamReader reader = resourceGetter.getXMLStreamReader("index.xml");
		while (reader.hasNext()) {
			if (reader.next() == XMLStreamConstants.START_ELEMENT &&
					reader.getLocalName().equals("tokenFile")) {
				validate(getUriForFile(reader.getElementText()));
			}
		}
		reader.close();
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
	
	@Test
	public void testTokenFilesValueValidity() throws Exception {
		XMLStreamReader indexReader = resourceGetter.getXMLStreamReader("index.xml");
		while (indexReader.hasNext()) {
			if (indexReader.next() == XMLStreamConstants.START_ELEMENT &&
					indexReader.getLocalName().equals("tokenFile")) {
				XMLStreamReader tokenReader = resourceGetter.getXMLStreamReader(indexReader.getElementText());
				while (tokenReader.hasNext()) {
					if (tokenReader.next() == XMLStreamConstants.START_ELEMENT) {
						String tagName = tokenReader.getLocalName();
						if (tagName.equals("tokenLists")) {
							while (tokenReader.hasNext()) {
								switch (tokenReader.next()) {
								case XMLStreamConstants.START_ELEMENT:
									if (tokenReader.getLocalName().equals("tokenList")) {
										validateTokenList(tokenReader);
									}
									break;
								}
							}
						}
						else if (tagName.equals("tokenList")) {
							validateTokenList(tokenReader);
						}
					}
				}
			}
		}
		indexReader.close();
	}
	
	private void validateTokenList(XMLStreamReader reader) throws XMLStreamException {
		Set<String> terms = new HashSet<String>();
		while (reader.hasNext()) {
			switch (reader.next()) {
			case XMLStreamConstants.START_ELEMENT:
				if (reader.getLocalName().equals("token")) {
					String tokenString = reader.getElementText();
					assertTrue(tokenString +" occurred more than once in a tokenList",!terms.contains(tokenString));
					terms.add(tokenString);
					char[] characters = tokenString.toCharArray();
					for (char c : characters) {
						assertTrue("Non ascii character found in token: " + tokenString + OpsinTools.NEWLINE + "An ASCII replacement should be used!" ,(int)c < 128);
						assertTrue("Capital letter found in token: " + tokenString + OpsinTools.NEWLINE + "Only lower case letters should be used!" , !(c >='A' && c <='Z'));
					}
				}
				break;
			case XMLStreamConstants.END_ELEMENT:
				if (reader.getLocalName().equals("tokenList")) {
					return;
				}
				break;
			}
		}
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
				System.out.println(h.error);
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
