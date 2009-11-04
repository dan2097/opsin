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

import nu.xom.Builder;
import nu.xom.Document;

/**Gets resource files from packages. Useful for including data in JAR files.
 * Cut down version of resourceGetter found in OSCAR
 * 
 * @author ptc24
 *
 */
public final class ResourceGetter {

	String resourcePath;
	String workspace;
	String resourcePrefix = "none";
	private boolean skipFiles = false; 
	
	/**Sets up a resourceGetter to get resources from a particular path.
	 *  /-separated - e.g. uk.ac.ch.cam.wwmm.ptclib.files.resources should be
	 *  /uk/ac/cam/ch/wwmm/ptclib/files/resources/
	 * 
	 * @param resourcePath The /-separated resource path.
	 */
	public ResourceGetter(String resourcePath) {
		if(resourcePath.startsWith("/")) resourcePath = resourcePath.substring(1);
		this.resourcePath = resourcePath;
		String currentWorkspace;
		try {
			currentWorkspace =new File("").getCanonicalPath();
		} catch (IOException e) {
			throw new RuntimeException("Unable to determine working directory");
		}
		this.workspace =currentWorkspace;
	}


	/**
	 * Sets up a resourceGetter to get resources from a particular path.
	 *  /-separated - e.g. uk.ac.ch.cam.wwmm.opsin.resources should be
	 *  /uk/ac/cam/ch/wwmm/opsin/resources/
	 * @param resourcePath The /-separated resource path.
	 * @param skipFiles Whether or not to skip reading files from the oscar3 workspace
	 * @param workspace
	 * @param resourcePrefix
	 */
	public ResourceGetter(String resourcePath, boolean skipFiles, String workspace, String resourcePrefix) {
		if(resourcePath.startsWith("/")) resourcePath = resourcePath.substring(1);
		this.resourcePath = resourcePath;
		this.skipFiles = skipFiles;
		this.workspace = workspace;
		this.resourcePrefix = resourcePrefix; 
		
		//workspace = (String) c.getField("workspace").get(oscar3Props);
		//resourcePrefix = (String) c.getField("resourcePrefix").get(oscar3Props);
	}
	
	private File getResDir() {
		if(skipFiles) return null;
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
		if(skipFiles) return null;
		File resourcesTop = new File(workspace, "resources");
		File resDir = new File(resourcesTop, resourcePath);
		if(!resDir.exists()) resDir.mkdirs();
		File f = new File(resDir, name);
		return f;
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
		if(skipFiles) return null;
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
			File f = getFile(name);
			if(f != null) {
				return new Builder().build(f);
			} else {
				ClassLoader l = getClass().getClassLoader();
				if(!skipFiles && !"none".equals(resourcePrefix)) {
					URL url = l.getResource(resourcePrefix + resourcePath + name);        					
					try {
						Document d = new Builder().build(url.toString());						
						if(d != null) return d;
					} catch (Exception e) {
						// Squelching the exceptions that come from failing to find a file here
					}
				}
				URL url = l.getResource(resourcePath + name);        
				Document d = new Builder().build(url.toString());
				return d;
			}
		} catch (Exception e) {
			e.printStackTrace();
			throw new Exception("Could not get resource file: " + name);
		}
	}

	/**Fetches a data file from resourcePath as an InputStream.
	 * 
	 * @param name The name of the file to get an InputStream of.
	 * @return An InputStream corresponding to the file.
	 * @throws Exception If the resouce file couldn't be found.
	 */
	public InputStream getStream(String name) throws Exception {
		if(name == null) name="";
		try {
			File f = getFile(name);
			if(f != null) {
				return new FileInputStream(f);
			} else {
				ClassLoader l = getClass().getClassLoader();
				if(!skipFiles && !"none".equals(resourcePrefix)) {
					URL url = l.getResource(resourcePrefix + resourcePath + name);        					
					try {
						if(url != null) {
							InputStream i = url.openStream();
							if(i != null) return i;
						}
					} catch (Exception e) {
						// Squelching the exceptions that come from failing to find a file here
					}
				}
				URL url = l.getResource(resourcePath + name);
				if(url == null) return null;
				InputStream i = url.openStream();
				return i;
			}
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
