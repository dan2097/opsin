package uk.ac.cam.ch.wwmm.opsin;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.Reader;
import java.io.UnsupportedEncodingException;
import java.net.URL;

import org.xml.sax.SAXException;
import org.xml.sax.SAXNotRecognizedException;
import org.xml.sax.SAXNotSupportedException;
import org.xml.sax.XMLReader;
import org.xml.sax.helpers.XMLReaderFactory;

import nu.xom.Builder;
import nu.xom.Document;

/**
 * Handles I/O:
 * Gets resource files from packages which is useful for including data from the JAR file.
 * Provides OutputStreams for the serialisation of automata.
 * This class has its roots in the resourceGetter in OSCAR.
 *
 * @author ptc24
 * @author dl387
 *
 */
class ResourceGetter {

	private final String resourcePath;
	private String workingDirectory;

	/**
	 * Sets up a resourceGetter to get resources from a particular path.
	 *  /-separated - e.g. uk.ac.ch.cam.wwmm.opsin.resources should be
	 *  /uk/ac/cam/ch/wwmm/opsin/resources/
	 *
	 * @param resourcePath The /-separated resource path.
	 */
	public ResourceGetter(String resourcePath) {
		if(resourcePath.startsWith("/")) {
			resourcePath = resourcePath.substring(1);
		}
		this.resourcePath = resourcePath;
		try {
			workingDirectory =new File("").getCanonicalPath();
		} catch (IOException e) {
			//Automata will not be serialisable
			workingDirectory = null;
		}
	}

	/**Fetches a data file from resourcePath,
	 * and parses it to an XML Document.
	 *
	 * @param name The name of the file to parse.
	 * @return The parsed document.
	 */
	public Document getXMLDocument(String name) {
		try {
			if (workingDirectory != null){
				File f = getFile(name);
				if(f != null) {
					return new Builder().build(f);
				}
			}
			ClassLoader l = getClass().getClassLoader();
			XMLReader xmlReader;
			try{
				xmlReader = XMLReaderFactory.createXMLReader();
			}
			catch (SAXException e) {
				throw new RuntimeException("No XML Reader could be initialised!", e);
			}
			try{
				xmlReader.setFeature("http://apache.org/xml/features/nonvalidating/load-external-dtd", false);
			}
			catch (SAXNotSupportedException e) {
				throw new RuntimeException("Your system's default XML Reader does not support disabling DTD loading! Maybe try updating your version of java?", e);
			} catch (SAXNotRecognizedException e) {
				throw new RuntimeException("Your system's default XML Reader has not recognised the DTD loading feature! Maybe try updating your version of java?", e);
			}
			Builder xomBuilder = new Builder(xmlReader);
			URL url = l.getResource(resourcePath + name);
			if (url == null){
				throw new RuntimeException("URL for resource: " + resourcePath + name + " is invalid");
			}
			return xomBuilder.build(url.openStream());
		} catch (Exception e) {
			throw new RuntimeException("Could not get resource file: " + name, e);
		}
	}

	private File getFile(String name) {
		File f = new File(getResDir(), name);
		if(f.isDirectory()){
			return null;
		}
		if(f.exists()){
			return f;
		}
		return null;
	}

	private File getResDir() {
		File resourcesTop = new File(workingDirectory, "resources");
		return new File(resourcesTop, resourcePath);
	}

	/**Fetches a data file from resourcePath, and returns the entire contents
	 * as a string.
	 *
	 * @param name The file to fetch.
	 * @return The string.
	 */
	public String getString(String name){
		try {
			return readText(new InputStreamReader(getStream(name), "UTF-8"));
		} catch (UnsupportedEncodingException e) {
			throw new RuntimeException("Java VM is broken; UTF-8 should be supported");
		} catch (IOException e) {
			throw new RuntimeException("Error reading stream from:" + name);
		}
	}

	/**Fetches a data file from resourcePath as an InputStream.
	 *
	 * @param name The name of the file to get an InputStream of.
	 * @return An InputStream corresponding to the file.
	 */
	public InputStream getStream(String name) {
		if(name == null){
			name="";
		}
		try {
			if (workingDirectory!=null){
				File f = getFile(name);
				if(f != null) {
					return new FileInputStream(f);
				}
			}
			ClassLoader l = getClass().getClassLoader();
			URL url = l.getResource(resourcePath + name);
			if (url == null){
				throw new RuntimeException("URL for resource: " + resourcePath + name + " is invalid");
			}
			return url.openStream();
		} catch (Exception e) {
			throw new RuntimeException("Could not get resource file: " + name, e);
		}
	}

	/**Reads a text file into a single string.
	 *
	 * @param r The Reader to read the text file.
	 * @return The string.
	 * @throws IOException 
	 */
	private String readText(Reader r) throws IOException {
		BufferedReader br = new BufferedReader(r);
		StringBuffer sb = new StringBuffer();
		while(br.ready()){
			sb.append((char)br.read());
		}
		br.close();
		return sb.toString();
	}

	/**Sets up an output stream to which a resource file can be written; this
	 * resource file will be in a subdirectory of the resources directory in
	 * the workspace.
	 *
	 * @param name The name of the file to write.
	 * @return The output stream.
	 */
	public OutputStream getOutputStream(String name) {
		try{
			File f = getFileForWriting(name);
			return new FileOutputStream(f);
		}
		catch (Exception e) {
			throw new RuntimeException("Failed to create outputstream", e);
		}
	}

	private File getFileForWriting(String name) throws IOException {
		File resourcesTop = new File(workingDirectory, "resources");
		File resDir = new File(resourcesTop, resourcePath);
		if(!resDir.exists()){
			if (!resDir.mkdirs()){
				throw new IOException("Failed to generate requested directories to create: " + name);
			}
		}
		return new File(resDir, name);
	}
}
