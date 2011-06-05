package uk.ac.cam.ch.wwmm.opsin;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.lang.reflect.Method;
import java.util.Collections;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.cli.UnrecognizedOptionException;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import uk.ac.cam.ch.wwmm.opsin.OpsinResult.OPSIN_RESULT_STATUS;

import nu.xom.Attribute;
import nu.xom.Element;

/** The "master" class, to turn a name into a structure.
 *
 * @author ptc24
 * @author dl387
 */
public class NameToStructure {
	
	private static final Logger LOG = Logger.getLogger(NameToStructure.class);

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
	 * @throws NameToStructureException If the converter cannot be initialised, most likely due to bad or missing data files.
	 */
	private NameToStructure() throws NameToStructureException {
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

	/**
	 * Convenience method for converting a name to CML with OPSIN's default options
	 * @param name The chemical name to parse.
	 * @return A CML element, containing the parsed molecule, or null if the name was uninterpretable.
	 */
	public Element parseToCML(String name) {
		OpsinResult result = parseChemicalName(name);
		Element cml = result.getCml();
		if(cml != null && LOG.isDebugEnabled()){
			LOG.debug(new XOMFormatter().elemToString(result.getCml()));
		}
		return cml;
	}

	/**
	 * Convenience method for converting a name to SMILES with OPSIN's default options
	 * @param name The chemical name to parse.
	 * @return A SMILES string describing the parsed molecule, or null if the name was uninterpretable.
	 */
	public String parseToSmiles(String name) {
		OpsinResult result = parseChemicalName(name);
		String smiles = result.getSmiles();
		LOG.debug(smiles);
		return smiles;
	}

	/**Parses a chemical name, returning an OpsinResult which represents the molecule.
	 * This object contains in the status whether the name was parsed successfully
	 * A message which may contain additional information if the status was warning/failure
	 * The OpsinResult has methods to generate a SMILES or CML representation
	 * For InChI, the OpsinResult should be given to the NameToInchi class
	 *
	 * @param name The chemical name to parse.
	 * @return OpsinResult
	 */
	public OpsinResult parseChemicalName(String name) {
		NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
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
		String message = "";
		try {
			LOG.debug(name);
			String modifiedName = PreProcessor.preProcess(name);
			List<Element> parses = parser.parse(n2sConfig, modifiedName);
			//if(LOG.isDebugEnabled()) for(Element parse : parses) LOG.debug(new XOMFormatter().elemToString(parse));
			Collections.sort(parses, new SortParses());//fewer tokens preferred
			Fragment frag = null;
			for(Element parse : parses) {
				try {
					if (LOG.isDebugEnabled()){
						LOG.debug(new XOMFormatter().elemToString(parse));
					}
					BuildState state = new BuildState(n2sConfig, sBuilder, cmlBuilder);
					componentGenerator.process(parse);
					if (LOG.isDebugEnabled()){
						LOG.debug(new XOMFormatter().elemToString(parse));
					}
					componentProcessor.process(state, parse);
					if (LOG.isDebugEnabled()){
						LOG.debug(new XOMFormatter().elemToString(parse));
					}
					frag = structureBuilder.buildFragment(state, parse);
					if (LOG.isDebugEnabled()){
						LOG.debug(new XOMFormatter().elemToString(parse));
					}
					break;
				} catch (Exception e) {
					if (message.equals("")){
						message += e.getMessage();
					}
					if (LOG.isDebugEnabled()){
						LOG.debug(e.getMessage(),e);
					}
				}
			}
			return new OpsinResult(frag, frag != null ? OPSIN_RESULT_STATUS.SUCCESS : OPSIN_RESULT_STATUS.FAILURE, message, name);
		} catch (Exception e) {
			message += e.getMessage();
			if(LOG.isDebugEnabled()) {
				LOG.debug(e.getMessage(),e);
			}
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

	/**Run OPSIN as a command-line application.
	 *
	 * @param args
	 * @throws Exception
	 */
	public static void main(String [] args) throws Exception {
		Options options = buildCommandLineOptions();
		CommandLineParser parser = new PosixParser();
		CommandLine cmd = null;
		try{
			cmd =  parser.parse(options, args);
		}
		catch (UnrecognizedOptionException e) {
			System.err.println(e.getMessage());
			System.exit(1);
		}
		if (cmd.hasOption("h")){
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("java -jar opsin-[version]-jar-with-dependencies.jar [options]\n" +
					"OPSIN accepts new line delimited names either interactively or in batch and can output CML, SMILES or InChI/StdInChI\n" +
					"For batch use direct a new line seperated list of names to the program and direct the stdout to an output files e.g. " +
					"java -jar opsin.jar -osmi < inputFile.name > outputFile.smiles", options);
			System.exit(0);
		}
		if(cmd.hasOption("v")){
			Logger.getLogger("uk.ac.cam.ch.wwmm.opsin").setLevel(Level.DEBUG);
		}
		NameToStructure nts = NameToStructure.getInstance();
		NameToStructureConfig n2sconfig = generateOpsinConfigObjectFromCmd(cmd);

		System.err.println("Welcome to OPSIN 1.0, use -h for help. Enter a chemical name:");
		String outputType = cmd.getOptionValue("o", "cml");
		if (outputType.equalsIgnoreCase("cml")){
			interactiveCmlOutput(nts, n2sconfig);
		}
		else if (outputType.equalsIgnoreCase("smi") || outputType.equalsIgnoreCase("smiles")){
			interactiveSmilesOutput(nts, n2sconfig);
		}
		else if (outputType.equalsIgnoreCase("inchi")){
			interactiveInchiOutput(nts, n2sconfig, false);
		}
		else if (outputType.equalsIgnoreCase("stdinchi")){
			interactiveInchiOutput(nts, n2sconfig, true);
		}
		else{
			System.err.println("Unrecognised output format: " + outputType);
			System.err.println("Expected output types are \"cml\", \"smi\", \"inchi\" and \"stdinchi\"");
			System.exit(1);
		}
	}

	private static Options buildCommandLineOptions() throws ParseException {
		Options options = new Options();
		OptionBuilder.withArgName("o");
		OptionBuilder.withLongOpt("output");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Sets OPSIN's output format (default cml)\n" +
				"Allowed values are:\n" +
				"cml for Chemical Markup Language\n" +
				"smi for SMILES\n" +
				"inchi for InChI\n" +
				"stdinchi for StdInChI");
		options.addOption(OptionBuilder.create("o"));
		options.addOption("h", "help", false, "Displays the allowed command line flags");
		options.addOption("v", "verbose", false, "Enables debugging");
		options.addOption("r", "allowRadicals", false, "Enables interpretation of radicals");
		options.addOption("f", "detailedFailureAnalysis", false, "Enables reverse parsing to more accurately determine why parsing failed");
		return options;
	}
	
	/**
	 * Uses the command line parameters to configure a new NameToStructureConfig
	 * @param cmd
	 * @return
	 */
	private static NameToStructureConfig generateOpsinConfigObjectFromCmd(CommandLine cmd) {
		NameToStructureConfig n2sconfig = new NameToStructureConfig();
		n2sconfig.setDetailedFailureAnalysis(cmd.hasOption("f"));
		n2sconfig.setAllowRadicals(cmd.hasOption("r"));
		return n2sconfig;
	}

	private static void interactiveCmlOutput(NameToStructure nts, NameToStructureConfig n2sconfig) throws IOException, NameToStructureException {
		BufferedReader stdinReader = new BufferedReader(new InputStreamReader(System.in));
		StreamSerializer serializer = new StreamSerializer(System.out);
		serializer.setIndent(2);
		serializer.writeXMLDeclaration();
		Element cml = new Element("cml", XmlDeclarations.CML_NAMESPACE);
		cml.addAttribute(new Attribute("convention","conventions:molecular"));
		cml.addNamespaceDeclaration("conventions", "http://www.xml-cml.org/convention/");
		cml.addNamespaceDeclaration("cmlDict", "http://www.xml-cml.org/dictionary/cml/");
		cml.addNamespaceDeclaration("nameDict", "http://www.xml-cml.org/dictionary/cml/name/");
		serializer.writeStartTag(cml);
		int id =1;
		String name = stdinReader.readLine();
		while(name !=null) {
			OpsinResult result = nts.parseChemicalName(name, n2sconfig);
			Element output = result.getCml();
			if(output == null) {
				System.err.println(result.getMessage());
				Element uninterpretableMolecule = new Element("molecule", XmlDeclarations.CML_NAMESPACE);
				uninterpretableMolecule.addAttribute(new Attribute("id", "m" + id++));
				Element nameEl = new Element("name", XmlDeclarations.CML_NAMESPACE);
				nameEl.appendChild(name);
				nameEl.addAttribute(new Attribute("dictRef", "nameDict:unknown"));
				uninterpretableMolecule.appendChild(nameEl);
				serializer.write(uninterpretableMolecule);
				serializer.flush();
			} else {
				Element molecule = XOMTools.getChildElementsWithTagName(output, "molecule").get(0);
				molecule.getAttribute("id").setValue("m" + id++);
				serializer.write(molecule);
				serializer.flush();
			}
			name = stdinReader.readLine();
		}
		serializer.writeEndTag(cml);
		serializer.flush();
	}
	
	private static void interactiveSmilesOutput(NameToStructure nts, NameToStructureConfig n2sconfig) throws IOException, NameToStructureException {
		BufferedReader stdinReader = new BufferedReader(new InputStreamReader(System.in));
		String name = stdinReader.readLine();
		while(name !=null) {
			OpsinResult result = nts.parseChemicalName(name, n2sconfig);
			String output = result.getSmiles();
			if(output == null) {
				System.err.println(result.getMessage());
				System.out.println("");
				System.out.flush();
			} else {
				System.out.println(output);
				System.out.flush();
			}
			name = stdinReader.readLine();
		}
	}
	
	@SuppressWarnings("unchecked")
	private static void interactiveInchiOutput(NameToStructure nts, NameToStructureConfig n2sconfig, boolean produceStdInChI) throws Exception {
		BufferedReader stdinReader = new BufferedReader(new InputStreamReader(System.in));
		Class c;
		try {
			c = Class.forName("uk.ac.cam.ch.wwmm.opsin.NameToInchi");
		} catch (ClassNotFoundException e) {
			System.err.println("Could not initialise NameToInChI module. Is it on your classpath?");
			throw new RuntimeException(e);
		}
		Method m;
		if (produceStdInChI){
			m = c.getMethod("convertResultToStdInChI", new Class[]{OpsinResult.class});
		}
		else{
			m = c.getMethod("convertResultToInChI", new Class[]{OpsinResult.class});
		}

		String name = stdinReader.readLine();
		while(name !=null) {
			OpsinResult result = nts.parseChemicalName(name, n2sconfig);
			String output = (String) m.invoke(null, new Object[]{result});
			if(output == null) {
				System.err.println(result.getMessage());
				System.out.println("");
				System.out.flush();
			} else {
				System.out.println(output);
				System.out.flush();
			}
			name = stdinReader.readLine();
		}
	}
}
