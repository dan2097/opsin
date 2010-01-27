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
import java.net.URL;

import org.xml.sax.XMLReader;
import org.xml.sax.helpers.XMLReaderFactory;

import nu.xom.Builder;
import nu.xom.Document;

/**Gets resource files from packages. Useful for including data in JAR files.
 * This is derived from the resourceGetter in OSCAR but with any extraneous functionality cut out.
 * It does NOT use OSCAR3props
 *
 * @author ptc24/d387
 *
 */
final class ResourceGetter {

	private String resourcePath;
	private String workspace;

	/**
	 * Sets up a resourceGetter to get resources from a particular path.
	 *  /-separated - e.g. uk.ac.ch.cam.wwmm.opsin.resources should be
	 *  /uk/ac/cam/ch/wwmm/opsin/resources/
	 *
	 * @param resourcePath The /-separated resource path.
	 * @throws Exception 
	 */
	public ResourceGetter(String resourcePath) {
		if(resourcePath.startsWith("/")) resourcePath = resourcePath.substring(1);
		this.resourcePath = resourcePath;
		try {
			workspace =new File("").getCanonicalPath();
		} catch (IOException e) {
			System.err.println("Unable to determine working directory");
			workspace = null;
		}
	}

	private File getResDir() {
		File resourcesTop = new File(workspace, "resources");
		return new File(resourcesTop, resourcePath);
	}

	private File getFile(String name) {
		File f = new File(getResDir(), name);
		if(f.isDirectory()) return null;
		if(f.exists()) return f;
		return null;
	}

	private File getFileForWriting(String name) {
		File resourcesTop = new File(workspace, "resources");
		File resDir = new File(resourcesTop, resourcePath);
		if(!resDir.exists()) resDir.mkdirs();
		return new File(resDir, name);
	}

	/**Sets up an output stream to which a resource file can be written; this
	 * resource file will be in a subdirectory of the resources directory in
	 * the workspace.
	 *
	 * @param name The name of the file to write.
	 * @return The output stream.
	 * @throws Exception
	 */
	public OutputStream getOutputStream(String name) throws Exception {
		File f = getFileForWriting(name);
		return new FileOutputStream(f);
	}

	/**Fetches a data file from resourcePath,
	 * and parses it to an XML Document.
	 *
	 * @param name The name of the file to parse.
	 * @return The parsed document.
	 * @throws Exception If the document can't be found, or can't parse, or is malformed/invalid.
	 */
	public Document getXMLDocument(String name) throws Exception {
		try {
			if (workspace != null){
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
			catch (Exception e) {
				e.printStackTrace();
				throw new Exception("No XML Reader could be initialised!");
			}
			try{
				xmlReader.setFeature("http://apache.org/xml/features/nonvalidating/load-external-dtd", false);
			}
			catch (Exception e) {
				e.printStackTrace();
				throw new Exception("Your system's default XML Reader does not support disabling DTD loading! Maybe try updating your version of java?");
			}
			Builder xomBuilder = new Builder(xmlReader);
			URL url = l.getResource(resourcePath + name);
			if (url == null){
				throw new Exception("URL for resource: " + resourcePath + name + " is invalid");
			}
			return xomBuilder.build(url.openStream());
		} catch (Exception e) {
			e.printStackTrace();
			throw new Exception("Could not get resource file: " + name);
		}
	}

	/**Fetches a data file from resourcePath as an InputStream.
	 *
	 * @param name The name of the file to get an InputStream of.
	 * @return An InputStream corresponding to the file.
	 * @throws Exception If the resource file couldn't be found.
	 */
	public InputStream getStream(String name) throws Exception {
		if(name == null) name="";
		try {
			if (workspace!=null){
				File f = getFile(name);
				if(f != null) {
					return new FileInputStream(f);
				}
			}
			ClassLoader l = getClass().getClassLoader();
			URL url = l.getResource(resourcePath + name);
			if (url == null){
				throw new Exception("URL for resource: " + resourcePath + name + " is invalid");
			}
			return url.openStream();
		} catch (Exception e) {
			e.printStackTrace();
			throw new Exception("Could not get resource file: " + name);
		}
	}

	/**Fetches a data file from resourcePath, and returns the entire contents
	 * as a string.
	 *
	 * @param name The file to fetch.
	 * @return The string.
	 * @throws Exception
	 */
	public String getString(String name) throws Exception {
		return readText(new InputStreamReader(getStream(name), "UTF-8"));
	}

	/**Reads a text file into a single string.
	 *
	 * @param r The Reader to read the text file.
	 * @return The string.
	 * @throws Exception
	 */
	public String readText(Reader r) throws Exception {
		BufferedReader br = new BufferedReader(r);
		StringBuffer sb = new StringBuffer();
		while(br.ready()) sb.append((char)br.read());
		br.close();
		return sb.toString();
	}

}
