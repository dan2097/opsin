package uk.ac.cam.ch.wwmm.opsin;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.net.JarURLConnection;
import java.net.URL;
import java.net.URLConnection;
import java.net.URLDecoder;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.jar.JarEntry;

import nu.xom.Builder;
import nu.xom.Document;
import sun.net.www.protocol.file.FileURLConnection;
import uk.ac.cam.ch.wwmm.ptclib.io.FileTools;
import uk.ac.cam.ch.wwmm.ptclib.string.StringTools;

/**Gets resource files from packages. Useful for incuding data in JAR files.
 * 
 * @author ptc24
 *
 */
public final class ResourceGetter {

	//private static Builder builder = new Builder();
	String resourcePath;
	String workspace;
	String resourcePrefix;
	
	/**
	 * Sets the current location of OPSIN. Paths to resources will be relative to this location.
	 * @param workspace
	 */
	public void setWorkspace(String workspace) {
		this.workspace = workspace;
	}

	/**
	 * Sets a prefix to be added to the location of OPSIN resources. By default there is none.
	 * @param resourcePrefix
	 */
	public void setResourcePrefix(String resourcePrefix) {
		this.resourcePrefix = resourcePrefix;
	}

	private boolean skipFiles; 
	
	/**Sets up a resourceGetter to get resources from a particular path.
	 *  /-separated - e.g. uk.ac.ch.cam.wwmm.ptclib.files.resources should be
	 *  /uk/ac/cam/ch/wwmm/ptclib/files/resources/
	 * 
	 * @param resourcePath The /-separated resource path.
	 */
	public ResourceGetter(String resourcePath) {
		this(resourcePath, false);
	}
	
	/**Sets up a resourceGetter to get resources from a particular path.
	 *  /-separated - e.g. uk.ac.ch.cam.wwmm.ptclib.files.resources should be
	 *  /uk/ac/cam/ch/wwmm/ptclib/files/resources/
	 * 
	 * @param resourcePath The /-separated resource path.
	 * @param skipFiles Whether or not to skip reading files from the oscar3 workspace
	 */
	@SuppressWarnings("unchecked")
	public ResourceGetter(String resourcePath, boolean skipFiles) {
		this.skipFiles = skipFiles;
		if(resourcePath.startsWith("/")) resourcePath = resourcePath.substring(1);
		this.resourcePath = resourcePath;
		
		//checks for the prescence of Oscar3Props by reflection. If it is found workspace/resourcePrefix will derive from it instead of local values
		Method oscar3PropsGetInstance;

		// TODO Booo! Stealth cyclic dependency!
		try {
			Class oscar3PropsClass = Class.forName("uk.ac.cam.ch.wwmm.oscar3.Oscar3Props");
			oscar3PropsGetInstance =oscar3PropsClass.getMethod("getInstance");
		} catch(Exception e){
			oscar3PropsGetInstance=null;
		}

		if (oscar3PropsGetInstance!=null){//OSCAR3 props appears to be present
			try {
				Object oscar3Props =oscar3PropsGetInstance.invoke(null);
				Class c =oscar3Props.getClass();
				try {
					workspace = (String) c.getField("workspace").get(oscar3Props);
				} catch (Exception e){
					throw new RuntimeException("Error attempting to read workspace field in oscar3Props");
				}
				try {
					resourcePrefix = (String) c.getField("resourcePrefix").get(oscar3Props);
				}catch (Exception e){
					throw new RuntimeException("Error attempting to read workspace field in oscar3Props");
				}
			} catch (InvocationTargetException e) {
				throw new RuntimeException(e);
			} catch (IllegalAccessException e) {
				throw new RuntimeException(e);
			}
		} else {
			try {
				workspace =new File("").getCanonicalPath();
			} catch (IOException e) {
				throw new RuntimeException("Unable to determine working directory");
			}
			resourcePrefix="none";
		}
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
	
	/**Fetches a data file from resourcePath, and writes it as the given file. 
	 * 
	 * @param name The resource to write.
	 * @param file The file to write it to.
	 * @throws Exception If the files cannot be read or written.
	 */
	public void writeToFile(String name, File file) throws Exception {
		FileOutputStream fos = new FileOutputStream(file);
		FileTools.pipeStreamToStream(getStream(name), fos);
		fos.close();
	}

	/**Copies the contents of the resourcePath into a new directory, recursively.
	 * WARNING: the resources directory must not contain files with no dot in them,
	 * as the presence/absence of a dot is taken to indicate whether or not a particular
	 * resource is a directory or not.
	 * 
	 * @param file
	 * @throws Exception
	 */
	public void writeDirRecursive(File file) throws Exception {
		if(!file.exists()) file.mkdirs();
		InputStream is = getStream("");
		BufferedReader br = new BufferedReader(new InputStreamReader(is));
		String line = br.readLine();
		while(line != null) {
			if(line.contains(".")) {
				writeToFile(line, new File(file, line));
			} else {
				ResourceGetter subRg = new ResourceGetter(resourcePath + line + "/");
				subRg.writeDirRecursive(new File(file, line));
			}
			line = br.readLine();
		}
	}
	
	/**Fetches a data file from resourcePath as an InputStream, removes comments starting with \s#, and
	 * returns each line in a list.
	 * 
	 * @param name The name of the file to get an InputStream of.
	 * @return A List of Strings corresponding to the file.
	 * @throws Exception If the resouce file couldn't be found.
	 */	
	public List<String> getStrings(String name) throws Exception {
		return getStrings(name, true);
	}

	private List<String> getFilesFromClasspath() throws Exception {
		List<String> files = new ArrayList<String>();
		ClassLoader l = Thread.currentThread().getContextClassLoader();
		Enumeration<URL> urls = l.getResources(resourcePath);
		while(urls.hasMoreElements()) {
			URL url = urls.nextElement();
			URLConnection conn = url.openConnection();
			if(conn instanceof JarURLConnection) {
				JarURLConnection jconn = (JarURLConnection)conn;
				if(jconn.getJarEntry().isDirectory()) {
					Enumeration<JarEntry> entries = jconn.getJarFile().entries();
					while(entries.hasMoreElements()) {
						JarEntry entry = entries.nextElement();
						String name = entry.getName();
						if(name.startsWith(resourcePath)) {
							String after = name.substring(resourcePath.length());
							if(after.length() == 0) {
								
							} else if(!after.contains("/")) {
								files.add(after);
							} else if(after.matches("[^/]+/")) {
								files.add(after.substring(0, after.length()-1));
							}
						}
					}		
				}
			} else if(conn instanceof FileURLConnection) {
				File f = new File(URLDecoder.decode(url.getPath(), "UTF-8"));
				if(f.exists() && f.isDirectory()) {
					File [] ff = f.listFiles();
					for(int i=0;i<ff.length;i++) {
						files.add(ff[i].getName());
					}
				}
			}
		}
		return files;
	}
	
	/**Gets a list of files that are available for this resourceGetter.
	 * 
	 * @return The available files.
	 * @throws Exception
	 */
	public List<String> getFiles() throws Exception {
		Set<String> seen = new LinkedHashSet<String>();
		try {
			seen.addAll(getFilesFromClasspath());
		} catch (Exception e) {
			e.printStackTrace();
		}
		if(resourcePath.equals("/") || resourcePath.equals("")) {
			seen.add("uk");
		}
		File resDir = getResDir();
		if(resDir != null && resDir.exists() && resDir.isDirectory()) {
			seen.addAll(StringTools.arrayToList(resDir.list()));
		}
		return new ArrayList<String>(seen);
	}
	
	/**Fetches a data file from resourcePath as an InputStream, removes comments starting with \s#, and
	 * returns each line in a list.
	 * 
	 * @param name The name of the file to get an InputStream of.
	 * @param UTF8 Whether to load the strings in UTF8
	 * @return A List of Strings corresponding to the file.
	 * @throws Exception If the resouce file couldn't be found.
	 */	
	public List<String> getStrings(String name, boolean UTF8) throws Exception {
		List<String> results = new ArrayList<String>();
    	InputStream is = getStream(name);
    	InputStreamReader isr;
		if(UTF8) {
    		isr = new InputStreamReader(is, "UTF-8");
    	} else {
    		isr = new InputStreamReader(is);
    	}
		BufferedReader br = new BufferedReader(isr);
    	String line = br.readLine();
    	while(line != null) {
    		line = line.split("\\s*#")[0];
    		if(line.length() == 0) {
        		line = br.readLine();
    			continue;
    		}
    		results.add(line);
    		line = br.readLine();
    	}
    	return results;
	}

	/**Fetches a data file from resourcePath as an InputStream, removes comments starting with \s#, and
	 * returns each line in a set.
	 * 
	 * @param name The name of the file to get an InputStream of.
	 * @return A Set of Strings corresponding to the file.
	 * @throws Exception If the resouce file couldn't be found.
	 */	
	public Set<String> getStringSet(String name) throws Exception {
		Set<String> results = new HashSet<String>();
    	BufferedReader br = new BufferedReader(new InputStreamReader(getStream(name), "UTF-8"));
    	String line = br.readLine();
    	while(line != null) {
    		line = line.split("\\s*#")[0];
    		if(line.length() == 0) {
        		line = br.readLine();
    			continue;
    		}
    		results.add(line);
    		line = br.readLine();
    	}
    	return results;
	}
	
	/**Fetches a data file from resourcePath, and returns the entire contents
	 * as a string.
	 * 
	 * @param name The file to fetch.
	 * @return The string.
	 * @throws Exception
	 */
	public String getString(String name) throws Exception {
		return FileTools.readText(new InputStreamReader(getStream(name), "UTF-8"));
	}
	
}
