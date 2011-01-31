package uk.ac.cam.ch.wwmm.opsin;

import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import nu.xom.Document;
import nu.xom.Element;
import nu.xom.Elements;

import org.junit.Test;
import org.xml.sax.ErrorHandler;
import org.xml.sax.SAXException;
import org.xml.sax.SAXParseException;
import static junit.framework.Assert.*;

public class DtdTest {
	private final static String RESOURCE_LOCATION = "uk/ac/cam/ch/wwmm/opsin/resources/";
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
	
	@Test
	public void testTokenFilesValueValidity() throws IOException {
		Document tokenFileDoc = resourceGetter.getXMLDocument("index.xml");
		Elements tokenFiles = tokenFileDoc.getRootElement().getChildElements();
		for (int i = 0; i < tokenFiles.size(); i++) {
			Element rootElement = resourceGetter.getXMLDocument(tokenFiles.get(i).getValue()).getRootElement();
			List<Element> tokenLists =new ArrayList<Element>();
			if (rootElement.getLocalName().equals("tokenLists")){//support for xml files with one "tokenList" or multiple "tokenList" under a "tokenLists" element
				Elements children =rootElement.getChildElements();
				for (int j = 0; j <children.size(); j++) {
					tokenLists.add(children.get(j));
				}
			}
			else{
				tokenLists.add(rootElement);
			}
			for (Element tokenList : tokenLists) {
				Elements tokenElements = tokenList.getChildElements("token");
				for(int j=0;j<tokenElements.size();j++) {
					String tokenString = tokenElements.get(j).getValue();
					assertTrue("The following token is less than 2 characters long!: " +tokenString, tokenString.length()>=2);
					assertEquals("The following token contains upper case characters!: " +tokenString,tokenString.toLowerCase(), tokenString);
				}
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
