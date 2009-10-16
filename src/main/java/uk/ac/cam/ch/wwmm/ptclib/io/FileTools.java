package uk.ac.cam.ch.wwmm.ptclib.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Reader;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.regex.Pattern;

/** Generic routines for dealing with files.
 * 
 * @author ptc24
 *
 */
public final class FileTools {

	/*public static File getCorrespondingFile(File f, File originalRootDir, File newRootDir) {
		File rf = new File(newRootDir.getAbsolutePath(),
				f.getAbsolutePath().substring(originalRootDir.getAbsolutePath().length()));
		rf.getParentFile().mkdirs();
		return rf;
	}*/

	/*public static File getCorrespondingFile(File f, File originalRootDir, File newRootDir, String extension) {
		String withOldExtension = f.getAbsolutePath().substring(originalRootDir.getAbsolutePath().length());
		String withoutExtension = withOldExtension.replaceAll("\\.[^\\.]+$", "");
		String withNewExtension = withoutExtension + extension;
		File rf = new File(newRootDir.getAbsolutePath(), withNewExtension);
		rf.getParentFile().mkdirs();
		return rf;
	}*/
	
	/*public static File getFileWithNewExtension(File f, String extension) {
		return new File(f.getAbsolutePath().replaceAll("\\.[^\\.]+$", "") + extension);
	}*/

	/**Gets a list of files, with a given name, that are in a directory,
	 * or subdirectories of that directory, recursively.
	 * 
	 * @param directory The directory to search.
	 * @param name The filename
	 * @return The files.
	 */
	public static List<File> getFilesFromDirectoryByName(File directory, String name) {
		List<File> foundFiles = new ArrayList<File>();
		File[] list = directory.listFiles();
		for(int i=0;i<list.length;i++) {
			if(list[i].getName().equals(name)) {
				foundFiles.add(list[i]);
			} else if(list[i].isDirectory()) {
				foundFiles.addAll(getFilesFromDirectoryByName(list[i], name));
			}
		}
		return foundFiles;
	}

	/**Gets a list of files, with a given suffix, that are in a directory,
	 * or subdirectories of that directory, recursively.
	 * 
	 * @param directory The directory to search.
	 * @param suffix The suffix.
	 * @return The files.
	 */
	public static List<File> getFilesFromDirectoryBySuffix(File directory, String suffix) {
		List<File> foundFiles = new ArrayList<File>();
		File[] list = directory.listFiles();
		for(int i=0;i<list.length;i++) {
			if(list[i].getAbsolutePath().endsWith(suffix)) {
				foundFiles.add(list[i]);
			} else if(list[i].isDirectory()) {
				foundFiles.addAll(getFilesFromDirectoryBySuffix(list[i], suffix));
			}
		}
		return foundFiles;
	}
	
	/**Gets a list of files, matching a given regex, that are in a directory,
	 * or subdirectories of that directory, recursively.
	 * 
	 * @param directory The directory to search.
	 * @param regex The regex.
	 * @return The files.
	 */
	public static List<File> getFilesFromDirectoryByRegex(File directory, String regex) {
		Pattern pattern = Pattern.compile(regex);
		List<File> foundFiles = new ArrayList<File>();
		File[] list = directory.listFiles();
		for(int i=0;i<list.length;i++) {
			if(pattern.matcher(list[i].getName()).matches()) {
				foundFiles.add(list[i]);
			} else if(list[i].isDirectory()) {
				foundFiles.addAll(getFilesFromDirectoryBySuffix(list[i], regex));
			}
		}
		return foundFiles;
	}

	/**Reads a text file into a single string.
	 * 
	 * @param r The Reader to read the text file.
	 * @return The string.
	 * @throws Exception
	 */
	public static String readText(Reader r) throws Exception {
		BufferedReader br = new BufferedReader(r);
		StringBuffer sb = new StringBuffer();
		while(br.ready()) sb.append((char)br.read());
		br.close();
		return sb.toString();
	}
	
	/**Reads a text file into a single string. UTF-8 is assumed.
	 * 
	 * @param file The file.
	 * @return The string.
	 * @throws Exception
	 */
	public static String readTextFile(File file) throws Exception {
		return readText(new InputStreamReader(new FileInputStream(file), "UTF-8"));
	}
	
	/**Sends the contents of one stream to another.
	 * 
	 * @param is The input stream.
	 * @param os The output stream.
	 * @throws IOException
	 */
	public static void pipeStreamToStream(InputStream is, OutputStream os) throws IOException {
		for(int i=is.read();i!=-1;i=is.read()) {
			os.write(i);
		}
	}
	
	/**Recursively delete a directory and all its contents.
	 * 
	 * @param dir The directory to delete.
	 */
	public static void deleteDir(File dir) {
        if (dir.isDirectory()) {
            String[] children = dir.list();
            for (int i=0; i<children.length; i++) {
                deleteDir(new File(dir, children[i]));
            }
        }
        dir.delete();
    }

	/**Gets the lines of a file as a list of strings, removing comments. 
	 * Assumes UTF-8. Comments appear as whitespace followed by a hash, 
	 * followed by anything.
	 * 
	 * @param is The input stream to read from.
	 * @return The list of strings.
	 * @throws Exception
	 */
	public static List<String> getStrings(InputStream is) throws Exception {
		return getStrings(is, true);
	}
	
	/**Gets the lines of a file as a list of strings, optionally removing 
	 * comments. Assumes UTF-8. Comments appear as whitespace followed by a
	 * hash, followed by anything.
	 * 
	 * @param is The input stream to read from.
	 * @param removeComments Whether to remove comments.
	 * @return The list of strings.
	 * @throws Exception
	 */
	public static List<String> getStrings(InputStream is, boolean removeComments) throws Exception {
    	InputStreamReader isr = new InputStreamReader(is, "UTF-8");
		BufferedReader br = new BufferedReader(isr);
    	String line = br.readLine();
    	List<String> results = new ArrayList<String>();
    	while(line != null) {
    		if(removeComments) line = line.split("\\s*#")[0];
    		if(line.length() == 0) {
        		line = br.readLine();
    			continue;
    		}
    		results.add(line);
    		line = br.readLine();
    	}
    	return results;
	}
}
