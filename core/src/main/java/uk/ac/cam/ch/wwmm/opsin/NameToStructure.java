package uk.ac.cam.ch.wwmm.opsin;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;


import nu.xom.Document;
import nu.xom.Element;
import nu.xom.Serializer;

/** The "master" class, to turn a name into a structure.
 *
 * @author dl387
 * @author ptc24
 */
public class NameToStructure {
	
	private static final Logger LOG = Logger.getLogger(NameToStructure.class);

	class SortParses implements Comparator<Element>{
		public int compare(Element el1, Element el2){
			int elementsInEl1 = XOMTools.countDescendantElements(el1);
			int elementsInEl2 = XOMTools.countDescendantElements(el2);
			if ( elementsInEl1> elementsInEl2){
				return 1;
			}
			else if (elementsInEl1 < elementsInEl2){
				return -1;
			}
			else{
				return 0;
			}
		}
	}

	/**Identifies name type and rejects a few special cases.*/
	private PreProcessor preProcessor;

	/**Does finite-state non-destructive parsing on chemical names.*/
	private Parser parser;

	/**Does destructive procedural parsing on parser results.*/
	private PostProcessor postProcessor;

	/**Does structure-aware destructive procedural parsing on parser results.*/
	private PreStructureBuilder preStructureBuilder;

	/** A builder for fragments specified as SMILES */
	private SMILESFragmentBuilder sBuilder;

	/** A builder for fragments specified as references to a CML data file */
	private CMLFragmentBuilder cmlBuilder;

	/**Constructs the CML molecule from the postProcessor results.*/
	private StructureBuilder structureBuilder;

	private static NameToStructure myInstance;

	public static void reinitialise() throws Exception {
		myInstance = null;
		getInstance();
	}

	public static NameToStructure getInstance() throws Exception {
		if(myInstance == null) myInstance = new NameToStructure();
		return myInstance;
	}


	/**Initialises the name-to-structure converter.
	 *
	 * @throws Exception If the converter cannot be initialised, most likely due to bad or missing data files.
	 */
	private NameToStructure() throws NameToStructureException {
		LOG.setLevel(Level.ALL);//FIXME remove this
		LOG.info("Initialising OPSIN... ");
		try {
			/*Initialise all of OPSIN's classes. Some classes are injected as dependencies into subsequent classes*/

			preProcessor = new PreProcessor();
			//Allows retrieving of OPSIN resources
			ResourceGetter resourceGetter = new ResourceGetter("uk/ac/cam/ch/wwmm/opsin/resources/");
			TokenManager tokenManager = new TokenManager(resourceGetter);
			WordRules wordRules = new WordRules(resourceGetter);
			ParseRules parseRules = new ParseRules(tokenManager, resourceGetter);
			parser = new Parser(wordRules, parseRules, tokenManager);

			postProcessor = new PostProcessor(tokenManager);

			sBuilder = new SMILESFragmentBuilder();
			cmlBuilder = new CMLFragmentBuilder(resourceGetter);
			structureBuilder = new StructureBuilder();

			FusedRingBuilder fusedRingBuilder = new FusedRingBuilder();
			preStructureBuilder = new PreStructureBuilder(fusedRingBuilder, resourceGetter);

		} catch (Exception e) {
			throw new NameToStructureException(e.getMessage(), e);
		}
		LOG.info("OPSIN initialised");
	}

	public Element parseToCML(String name) {
		return parseToCML(name, false);
	}

	/**Parses a chemical name, returning an unambiguous CML representation of the molecule.
	 *
	 * @param name The chemical name to parse.
	 * @param verbose Whether to print lots of debugging information to stdin and stderr or not.
	 * @return A CML element, containing the parsed molecule, or null if the molecule would not parse.
	 */
	public Element parseToCML(String name, boolean verbose) {
		Fragment frag = parseToOpsinFragment(name, verbose);
		if (frag == null){
			return null;
		}
		else{
			Element cml = null;
			try{
				cml = frag.toCMLMolecule(name);
			}
			catch (Exception e) {
				if (verbose){
					e.printStackTrace();
				}
				return null;
			}
			if(verbose) System.out.println(new XOMFormatter().elemToString(cml));
			return cml;
		}
	}

	/**Parses a chemical name, returning an OPSIN fragment which represents the molecule.
	 * This is null if the name cannot be interpreted
	 *
	 * @param name The chemical name to parse.
	 * @param verbose Whether to print lots of debugging information to stdin and stderr or not.
	 * @return An OPSIN fragment containing the parsed molecule, or null if the molecule would not parse.
	 */
	public synchronized Fragment parseToOpsinFragment(String name, boolean verbose) {
		if (name==null){
			throw new IllegalArgumentException("String given for name was null");
		}
		try {
			name = preProcessor.preProcess(name);
			if (name ==null){return null;}//not a specific chemical name e.g. amine/carboxylic acid or a blank string
			if(verbose) System.out.println(name);
			List<Element> parses = parser.parse(name);
			//if(verbose) for(Element parse : parses) System.out.println(new XOMFormatter().elemToString(parse));
			Collections.sort(parses, new SortParses());//fewer tokens preferred
			Fragment frag = null;
			for(Element parse : parses) {
				try {
					if (verbose && LOG.isDebugEnabled()){
						LOG.debug(new XOMFormatter().elemToString(parse));
					}
					BuildState state = new BuildState(sBuilder, cmlBuilder);
					postProcessor.postProcess(parse, state);
					if (verbose && LOG.isDebugEnabled()){
						LOG.debug(new XOMFormatter().elemToString(parse));
					}
					preStructureBuilder.postProcess(state, parse);
					if (verbose && LOG.isDebugEnabled()){
						LOG.debug(new XOMFormatter().elemToString(parse));
					}
					frag = structureBuilder.buildFragment(state, parse);
					if (verbose && LOG.isDebugEnabled()){
						LOG.debug(new XOMFormatter().elemToString(parse));
					}
					break;
				} catch (Exception e) {
					if (verbose && LOG.isDebugEnabled()){
						LOG.warn(e.getMessage(),e);
					}
				}
			}
			return frag;
		} catch (Exception e) {
			if(verbose) e.printStackTrace();
			return null;
		}
	}

	/**Run OPSIN as a standalone component.
	 *
	 * @param args
	 * @throws Exception
	 */
	public static void main(String [] args) throws Exception {
		NameToStructure nts = new NameToStructure();
		Serializer serializer = new Serializer(System.out);
		serializer.setIndent(2);
		boolean end = false;
		BufferedReader stdinReader = new BufferedReader(new InputStreamReader(System.in));
		System.out.println("OPSIN Prealpha: enter chemical name:");
		while(!end) {
			String name = stdinReader.readLine();
			if(name == null) {
				System.err.println("Disconnected!");
				end = true;
			} else if(name.equals("END")) {
				end = true;
			} else {
				Element output = nts.parseToCML(name);
				if(output == null) {
					System.out.println("Did not parse.");
					System.out.flush();
				} else {
					serializer.write(new Document(output));
					System.out.flush();
				}
			}
		}
	}
}
