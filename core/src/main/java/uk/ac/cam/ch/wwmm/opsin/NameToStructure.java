package uk.ac.cam.ch.wwmm.opsin;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import uk.ac.cam.ch.wwmm.opsin.OpsinResult.OPSIN_RESULT_STATUS;

import nu.xom.Attribute;
import nu.xom.Element;

/** The "master" class, to turn a name into a structure.
 *
 * @author dl387
 * @author ptc24
 */
public class NameToStructure {
	
	private static final Logger LOG = Logger.getLogger(NameToStructure.class);

	/**
	 * Prefer less childless elements e.g. benzal beats benz al
	 * Prefer less elements e.g. <acryl(acidStem)amide(suffix)> beats <acryl(substituent)><amide(group)>
	 */
	private class SortParses implements Comparator<Element>{
		public int compare(Element el1, Element el2){
			int[] counts1 = XOMTools.countNumberOfElementsAndNumberOfChildLessElements(el1);
			int[] counts2 = XOMTools.countNumberOfElementsAndNumberOfChildLessElements(el2);
			int childLessElementsInEl1 = counts1[1];
			int childLessElementsInEl2 = counts2[1];
			if ( childLessElementsInEl1> childLessElementsInEl2){
				return 1;
			}
			else if (childLessElementsInEl1 < childLessElementsInEl2){
				return -1;
			}

			int elementsInEl1 = counts1[0];
			int elementsInEl2  = counts2[0];
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

	/**Does finite-state non-destructive parsing on chemical names.*/
	private Parser parser;
	
	/**Applies OPSIN's grammar to tokenise and assign meanings tokens.*/
	private ParseRules parseRules;

	/**Does destructive procedural parsing on parser results.*/
	private ComponentGenerator componentGenerator;

	/**Does structure-aware destructive procedural parsing on parser results.*/
	private ComponentProcessor componentProcessor;

	/** A builder for fragments specified as SMILES */
	private SMILESFragmentBuilder sBuilder;

	/** A builder for fragments specified as references to a CML data file */
	private CMLFragmentBuilder cmlBuilder;

	/**Constructs a single fragment from the result of the component generation and processing stages.*/
	private StructureBuilder structureBuilder;

	private static NameToStructure NTS_INSTANCE;

	public static synchronized NameToStructure getInstance() throws NameToStructureException {
		if (NTS_INSTANCE ==null){
			NTS_INSTANCE = new NameToStructure();
		}
		return NTS_INSTANCE;
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

			//Allows retrieving of OPSIN resources
			ResourceGetter resourceGetter = new ResourceGetter("uk/ac/cam/ch/wwmm/opsin/resources/");
			ResourceManager resourceManager = new ResourceManager(resourceGetter);
			WordRules wordRules = new WordRules(resourceGetter);
			parseRules = new ParseRules(resourceManager);
			Tokeniser tokeniser = new Tokeniser(parseRules);
			parser = new Parser(wordRules, tokeniser, resourceManager);

			componentGenerator = new ComponentGenerator();

			sBuilder = new SMILESFragmentBuilder();
			cmlBuilder = new CMLFragmentBuilder(resourceGetter);
			structureBuilder = new StructureBuilder();

			componentProcessor = new ComponentProcessor(resourceGetter);

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
		OpsinResult result = parseChemicalName(name, verbose);
		if(result.getCml() != null && verbose) System.out.println(new XOMFormatter().elemToString(result.getCml()));
		return result.getCml();
	}

	/**Parses a chemical name, returning an OpsinResult which represents the molecule.
	 * This object contains in the status whether the name was parsed successfully
	 * A message which may contain additional information if the status was warning/failure
	 * A CML representation of the structure
	 *
	 * @param name The chemical name to parse.
	 * @param verbose Whether to print lots of debugging information to stdin and stderr or not.
	 * @return OpsinResult
	 */
	public OpsinResult parseChemicalName(String name, boolean verbose) {
		NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
		n2sConfig.setVerbose(verbose);
		return parseChemicalName(name, n2sConfig);
	}

	/**Parses a chemical name, returning an OpsinResult which represents the molecule.
	 * This object contains in the status whether the name was parsed successfully
	 * A message which may contain additional information if the status was warning/failure
	 * A CML representation of the structure
	 *
	 * @param name The chemical name to parse.
	 * @param n2sConfig Options to control how OPSIN interprets the name.
	 * @return OpsinResult
	 */
	public OpsinResult parseChemicalName(String name, NameToStructureConfig n2sConfig) {
		n2sConfig =  n2sConfig.clone();//avoid n2sconfig being modified mid name processing
		if (name==null){
			throw new IllegalArgumentException("String given for name was null");
		}
		boolean verbose = n2sConfig.isVerbose();
		String message = "";
		try {
			if(verbose) System.out.println(name);
			String modifiedName = PreProcessor.preProcess(name);
			List<Element> parses = parser.parse(n2sConfig, modifiedName);
			//if(verbose) for(Element parse : parses) System.out.println(new XOMFormatter().elemToString(parse));
			Collections.sort(parses, new SortParses());//fewer tokens preferred
			Fragment frag = null;
			for(Element parse : parses) {
				try {
					if (verbose && LOG.isDebugEnabled()){
						LOG.debug(new XOMFormatter().elemToString(parse));
					}
					BuildState state = new BuildState(n2sConfig, sBuilder, cmlBuilder);
					componentGenerator.process(parse, state);
					if (verbose && LOG.isDebugEnabled()){
						LOG.debug(new XOMFormatter().elemToString(parse));
					}
					componentProcessor.process(state, parse);
					if (verbose && LOG.isDebugEnabled()){
						LOG.debug(new XOMFormatter().elemToString(parse));
					}
					frag = structureBuilder.buildFragment(state, parse);
					if (verbose && LOG.isDebugEnabled()){
						LOG.debug(new XOMFormatter().elemToString(parse));
					}
					break;
				} catch (Exception e) {
					if (message.equals("")){
						message += e.getMessage();
					}
					if (verbose && LOG.isDebugEnabled()){
						LOG.warn(e.getMessage(),e);
					}
				}
			}
			return new OpsinResult(frag, frag != null ? OPSIN_RESULT_STATUS.SUCCESS : OPSIN_RESULT_STATUS.FAILURE, message, name);
		} catch (Exception e) {
			message += e.getMessage();
			if(verbose) e.printStackTrace();
			return new OpsinResult(null, OPSIN_RESULT_STATUS.FAILURE, message, name);
		}
	}
	
	/**
	 * Returns an OPSIN parser
	 * This can be used to determine whether a word can be interpreted as being part of a chemical name.
	 * Just because a word can be split into tokens does not mean the word constitutes a valid chemical name
	 * e.g. ester is interpretable but is not in itself a chemical name
	 * @return
	 * @throws NameToStructureException
	 */
	public static ParseRules getOpsinParser() throws NameToStructureException{
		NameToStructure n2s = NameToStructure.getInstance();
		return n2s.parseRules;
	}

	/**Run OPSIN as a standalone component.
	 *
	 * @param args
	 * @throws Exception
	 */
	public static void main(String [] args) throws Exception {
		NameToStructure nts = NameToStructure.getInstance();
		StreamSerializer serializer = new StreamSerializer(System.out);
		serializer.setIndent(2);
		boolean end = false;
		BufferedReader stdinReader = new BufferedReader(new InputStreamReader(System.in));
		System.err.println("OPSIN Prealpha: enter chemical name:");
		serializer.writeXMLDeclaration();
		Element cml = new Element("cml");
		cml.addAttribute(new Attribute("convention","cmlDict:cmllite"));
		cml.addNamespaceDeclaration("cmlDict", "http://www.xml-cml.org/dictionary/cml/");
		cml.addNamespaceDeclaration("nameDict", "http://www.xml-cml.org/dictionary/cml/name/");
		cml.setNamespaceURI("http://www.xml-cml.org/schema");
		serializer.writeStartTag(cml);
		int id =1;
		while(!end) {
			String name = stdinReader.readLine();
			if(name == null || name.equals("END")) {
				end = true;
			} else {
				OpsinResult result = nts.parseChemicalName(name, false);
				Element output = result.getCml();
				if(output == null) {
					System.err.println(result.getMessage());
					Element uninterpretableMolecule = new Element("molecule");
					uninterpretableMolecule.addAttribute(new Attribute("id", "m" + id++));
					Element nameEl = new Element("name");
					nameEl.appendChild(name);
					nameEl.addAttribute(new Attribute("dictRef", "nameDict:unknown"));
					uninterpretableMolecule.appendChild(nameEl);
					XOMTools.setNamespaceURIRecursively(uninterpretableMolecule, "http://www.xml-cml.org/schema");
					serializer.write(uninterpretableMolecule);
					serializer.flush();
				} else {
					Element molecule = XOMTools.getChildElementsWithTagName(output, "molecule").get(0);
					molecule.getAttribute("id").setValue("m" + id++);
					serializer.write(molecule);
					serializer.flush();
				}
			}
		}
		serializer.writeEndTag(cml);
		serializer.flush();
	}
}
