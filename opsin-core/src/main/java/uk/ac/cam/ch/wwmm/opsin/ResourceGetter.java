package uk.ac.cam.ch.wwmm.opsin;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.UnsupportedEncodingException;
import java.net.URL;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;

import org.apache.commons.io.IOUtils;
import org.codehaus.stax2.XMLInputFactory2;

import com.ctc.wstx.stax.WstxInputFactory;

/**
 * Handles I/O:
 * Gets resource files from packages which is useful for including data from the JAR file.
 * Provides OutputStreams for the serialisation of automata.
 *
 * @author ptc24
 * @author dl387
 *
 */
class ResourceGetter {

	private final String resourcePath;
	private final String workingDirectory;
	private final XMLInputFactory xmlInputFactory;

	/**
	 * Sets up a resourceGetter to get resources from a particular path.
	 *  /-separated - e.g. uk.ac.ch.cam.wwmm.opsin.resources should be
	 *  /uk/ac/cam/ch/wwmm/opsin/resources/
	 *
	 * @param resourcePath The /-separated resource path.
	 */
	ResourceGetter(String resourcePath) {
		if(resourcePath.startsWith("/")) {
			resourcePath = resourcePath.substring(1);
		}
		this.resourcePath = resourcePath;
		String workingDirectory;
		try {
			workingDirectory = new File(".").getCanonicalPath();//works on linux unlike using the system property
		} catch (IOException e) {
			//Automata will not be serialisable
			workingDirectory = null;
		}
		this.workingDirectory = workingDirectory;

		xmlInputFactory = new WstxInputFactory();
		xmlInputFactory.setProperty(XMLInputFactory.SUPPORT_DTD, false);
		xmlInputFactory.setProperty(XMLInputFactory2.P_AUTO_CLOSE_INPUT, true);
	}
	
	/**
	 * Gets the resourcePath used to initialise this ResourceGetter
	 * @return
	 */
	String getResourcePath() {
		return resourcePath;
	}
	
	/**Fetches a data file from resourcePath,
	 * and returns an XML stream reader for it
	 *
	 * @param name The name of the file to parse.
	 * @return An XMLStreamReader
	 * @throws IOException 
	 */
	XMLStreamReader getXMLDocument2(String name) throws IOException {
		if(name == null){
			throw new IllegalArgumentException("Input to function was null");
		}
		try {
			if (workingDirectory != null){
				File f = getFile(name);
				if(f != null) {
					return xmlInputFactory.createXMLStreamReader(new FileInputStream(f));
				}
			}
			ClassLoader l = getClass().getClassLoader();
			URL url = l.getResource(resourcePath + name);
			if (url == null){
				throw new IOException("URL for resource: " + resourcePath + name + " is invalid");
			}
			return xmlInputFactory.createXMLStreamReader(url.openStream());
		} catch (XMLStreamException e) {
			throw new IOException("Validity exception occurred while reading the XML file with name:" +name, e);
		}
	}

	private File getFile(String name) {
		File f = new File(getResDir(), name);
		if(f.isFile()){
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
	 * @return The contents of the file as a string or "" if an IOException occurred
	 */
	String getFileContentsAsString(String name){
		if(name == null){
			throw new IllegalArgumentException("Input to function was null");
		}
		InputStreamReader is = null;
		try {
			try {
				is = new InputStreamReader(getInputstreamFromFileName(name), "UTF-8");
				return IOUtils.toString(is);
			} catch (UnsupportedEncodingException e) {
				throw new RuntimeException("Java VM is broken; UTF-8 should be supported", e);
			}
			finally{
				IOUtils.closeQuietly(is);
			}
		} catch (IOException e) {
			return "";
		}
	}

	/**Fetches a data file from the working directory or resourcePath as an InputStream.
	 *
	 * @param name The name of the file to get an InputStream of.
	 * @return An InputStream corresponding to the file.
	 * @throws IOException 
	 */
	InputStream getInputstreamFromFileName(String name) throws IOException {
		if(name == null){
			throw new IllegalArgumentException("Input to function was null");
		}
		if (workingDirectory!=null){
			File f = getFile(name);
			if(f != null) {
				return new FileInputStream(f);
			}
		}
		ClassLoader l = getClass().getClassLoader();
		URL url = l.getResource(resourcePath + name);
		if (url == null){
			throw new IOException("URL for resource: " + resourcePath + name + " is invalid");
		}
		return url.openStream();
	}

	/**Sets up an output stream to which a resource file can be written; this
	 * resource file will be in a subdirectory of the resources directory in
	 * the working directory.
	 *
	 * @param name The name of the file to write.
	 * @return The output stream.
	 * @throws IOException 
	 */
	OutputStream getOutputStream(String name) throws IOException {
		if(name == null){
			throw new IllegalArgumentException("Input to function was null");
		}
		File f = getFileForWriting(name);
		return new FileOutputStream(f);
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
