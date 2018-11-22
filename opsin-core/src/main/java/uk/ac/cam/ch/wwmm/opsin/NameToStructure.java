package uk.ac.cam.ch.wwmm.opsin;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.lang.reflect.Method;
import java.util.Collections;
import java.util.List;
import java.util.Properties;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Option.Builder;
import org.apache.commons.io.IOUtils;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.UnrecognizedOptionException;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import com.ctc.wstx.api.WstxOutputProperties;
import com.ctc.wstx.stax.WstxOutputFactory;

import uk.ac.cam.ch.wwmm.opsin.OpsinResult.OPSIN_RESULT_STATUS;

/** The "master" class, to turn a name into a structure.
 *
 * @author ptc24
 * @author dl387
 */
public class NameToStructure {
	
	private static final Logger LOG = Logger.getLogger(NameToStructure.class);
	
	/**Applies OPSIN's grammar to tokenise and assign meaning to tokens*/
	private ParseRules parseRules;

	/**Parses a chemical name into one (or more in the case of ambiguity) parse trees*/
	private Parser parser;
	
	/**Which suffixes apply to what and what their effects are*/
	private SuffixRules suffixRules;

	private static NameToStructure NTS_INSTANCE;

	public static synchronized NameToStructure getInstance() {
		if (NTS_INSTANCE == null) {
			NTS_INSTANCE = new NameToStructure();
		}
		return NTS_INSTANCE;
	}
	
	/**
	 * Returns the version of the OPSIN library
	 * @return Version number String
	 */
	public static String getVersion() {
		try {
			InputStream is = NameToStructure.class.getResourceAsStream("opsinbuild.props");
			try {
				Properties props = new Properties();
				props.load(is);
				return props.getProperty("version");
			}
			finally {
				IOUtils.closeQuietly(is);
			}
		}
		catch (Exception e) {
			return null;
		}
	}

	/**Initialises the name-to-structure converter.
	 *
	 * @throws NameToStructureException If the converter cannot be initialised, most likely due to bad or missing data files.
	 */
	private NameToStructure() {
		LOG.debug("Initialising OPSIN... ");
		try {
			/*Initialise all of OPSIN's classes. Some classes are injected as dependencies into subsequent classes*/

			//Allows retrieving of OPSIN resources
			ResourceGetter resourceGetter = new ResourceGetter("uk/ac/cam/ch/wwmm/opsin/resources/");
			ResourceManager resourceManager = new ResourceManager(resourceGetter);
			WordRules wordRules = new WordRules(resourceGetter);
			parseRules = new ParseRules(resourceManager);
			Tokeniser tokeniser = new Tokeniser(parseRules);
			parser = new Parser(wordRules, tokeniser, resourceManager);
			suffixRules = new SuffixRules(resourceGetter);
		} catch (Exception e) {
			throw new NameToStructureException(e.getMessage(), e);
		}
		LOG.debug("OPSIN initialised");
	}

	/**
	 * Convenience method for converting a name to CML with OPSIN's default options
	 * @param name The chemical name to parse.
	 * @return A CML element, containing the parsed molecule, or null if the name was uninterpretable.
	 */
	public String parseToCML(String name) {
		OpsinResult result = parseChemicalName(name);
		String cml = result.getCml();
		if(cml != null && LOG.isDebugEnabled()){
			LOG.debug(cml);
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
	 * CML and SMILES representations may be retrieved directly from the object
	 * InChI may be generate using NameToInchi
	 *
	 * @param name The chemical name to parse.
	 * @param n2sConfig Options to control how OPSIN interprets the name.
	 * @return OpsinResult
	 */
	public OpsinResult parseChemicalName(String name, NameToStructureConfig n2sConfig) {
		if (name == null){
			throw new IllegalArgumentException("String given for name was null");
		}
		n2sConfig = n2sConfig.clone();//avoid n2sconfig being modified mid name processing

		List<Element> parses;
		try {
			LOG.debug(name);
			String modifiedName = PreProcessor.preProcess(name);
			parses = parser.parse(n2sConfig, modifiedName);
			Collections.sort(parses, new SortParses());//fewer tokens preferred
		} catch (Exception e) {
			if(LOG.isDebugEnabled()) {
				LOG.debug(e.getMessage(), e);
			}
			String message = e.getMessage() != null ? e.getMessage() : "exception with null message";
			return new OpsinResult(null, OPSIN_RESULT_STATUS.FAILURE, message, name);
		}
		String reasonForFailure = "";
		Fragment fragGeneratedWithWarning = null;
		List<OpsinWarning> warnings = Collections.emptyList();
		for(Element parse : parses) {
			try {
				if (LOG.isDebugEnabled()) {
					LOG.debug(parse.toXML());
				}
				//Performs XML manipulation e.g. nesting bracketing, processing some nomenclatures
				new ComponentGenerator(n2sConfig).processParse(parse);
				if (LOG.isDebugEnabled()) {
					LOG.debug(parse.toXML());
				}
				BuildState state = new BuildState(n2sConfig);
				//Converts the XML to fragments (handles many different nomenclatueres for describing structure). Assigns locants 
				new ComponentProcessor(state, new SuffixApplier(state, suffixRules)).processParse(parse);
				if (LOG.isDebugEnabled()) {
					LOG.debug(parse.toXML());
				}
				//Constructs a single fragment from the fragments generated by the ComponentProcessor. Applies stereochemistry
				Fragment frag = new StructureBuilder(state).buildFragment(parse);
				if (LOG.isDebugEnabled()) {
					LOG.debug(parse.toXML());
				}
				if (state.getWarnings().size() == 0) {
					return new OpsinResult(frag, OPSIN_RESULT_STATUS.SUCCESS, "", name);
				}
				if (fragGeneratedWithWarning == null) {
					//record first frag that had a warning but try other parses as they may work without a warning
					fragGeneratedWithWarning = frag;
					warnings = state.getWarnings();
				}
			} catch (Exception e) {
				if (reasonForFailure.length() == 0) {
					reasonForFailure = e.getMessage() != null ? e.getMessage() : "exception with null message";
				}
				if (LOG.isDebugEnabled()) {
					LOG.debug(e.getMessage(), e);
				}
			}
		}
		if (fragGeneratedWithWarning != null) {
			return new OpsinResult(fragGeneratedWithWarning, OPSIN_RESULT_STATUS.WARNING, warnings, name);
		}
		return new OpsinResult(null, OPSIN_RESULT_STATUS.FAILURE, reasonForFailure, name);
	}
	
	/**
	 * Returns an OPSIN parser
	 * This can be used to determine whether a word can be interpreted as being part of a chemical name.
	 * Just because a word can be split into tokens does not mean the word constitutes a valid chemical name
	 * e.g. ester is interpretable but is not in itself a chemical name
	 * @return Opsin parser for recognition/parsing of a chemical word
	 */
	public static ParseRules getOpsinParser() {
		NameToStructure n2s = NameToStructure.getInstance();
		return n2s.parseRules;
	}
	
	private enum InchiType{
		inchiWithFixedH,
		stdInchi,
		stdInchiKey
	}

	/**Run OPSIN as a command-line application.
	 *
	 * @param args
	 * @throws Exception
	 */
	public static void main(String [] args) throws Exception {
		Options options = buildCommandLineOptions();
		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = null;
		try{
			cmd = parser.parse(options, args);
		}
		catch (UnrecognizedOptionException e) {
			System.err.println(e.getMessage());
			System.exit(1);
		}
		if (cmd.hasOption("h")){
			displayUsage(options);
		}
		if(cmd.hasOption("v")){
			Logger.getLogger("uk.ac.cam.ch.wwmm.opsin").setLevel(Level.DEBUG);
		}

		NameToStructureConfig n2sconfig = generateOpsinConfigObjectFromCmd(cmd);
		
		InputStream input;
		OutputStream output;
		String[] unparsedArgs = cmd.getArgs();
		if (unparsedArgs.length == 0){
			input = System.in;
			output = System.out;
		}
		else if (unparsedArgs.length == 1){
			input = new FileInputStream(new File(unparsedArgs[0]));
			output = System.out;
		}
		else if (unparsedArgs.length == 2){
			input = new FileInputStream(new File(unparsedArgs[0]));
			output = new FileOutputStream(new File(unparsedArgs[1]));
		}
		else {
			input = null;
			output = null;
			displayUsage(options);
		}

		System.err.println("Run the jar using the -h flag for help. Enter a chemical name to begin:");
		String outputType = cmd.getOptionValue("o", "smi");
		boolean outputName = cmd.hasOption("n");
		if (outputType.equalsIgnoreCase("cml")) {
			interactiveCmlOutput(input, output, n2sconfig);
		}
		else if (outputType.equalsIgnoreCase("smi") || outputType.equalsIgnoreCase("smiles")) {
			interactiveSmilesOutput(input, output, n2sconfig, false, outputName);
		}
		else if (outputType.equalsIgnoreCase("inchi")) {
			interactiveInchiOutput(input, output, n2sconfig, InchiType.inchiWithFixedH, outputName);
		}
		else if (outputType.equalsIgnoreCase("stdinchi")) {
			interactiveInchiOutput(input, output, n2sconfig, InchiType.stdInchi, outputName);
		}
		else if (outputType.equalsIgnoreCase("stdinchikey")) {
			interactiveInchiOutput(input, output, n2sconfig, InchiType.stdInchiKey, outputName);
		}
		else if (outputType.equalsIgnoreCase("extendedsmi") || outputType.equalsIgnoreCase("extendedsmiles") || 
				outputType.equalsIgnoreCase("cxsmi") || outputType.equalsIgnoreCase("cxsmiles")) {
			interactiveSmilesOutput(input, output, n2sconfig, true, outputName);
		}
		else{
			System.err.println("Unrecognised output format: " + outputType);
			System.err.println("Expected output types are \"cml\", \"smi\", \"inchi\", \"stdinchi\" and \"stdinchikey\"");
			System.exit(1);
		}
		if (unparsedArgs.length == 1) {
			input.close();
		}
		else if (unparsedArgs.length == 2) {
			input.close();
			output.close();
		}
	}

	private static void displayUsage(Options options) {
		HelpFormatter formatter = new HelpFormatter();
		String version = getVersion();
		formatter.printHelp("java -jar opsin-" + (version != null ? version : "[version]") + "-jar-with-dependencies.jar [options] [inputfile] [outputfile]" + OpsinTools.NEWLINE +
				"OPSIN converts systematic chemical names to CML, SMILES or InChI/StdInChI/StdInChIKey" + OpsinTools.NEWLINE +
				"Names should be new line delimited and may be read from stdin (default) or a file and output to stdout (default) or a file", options);
		System.exit(0);
	}

	private static Options buildCommandLineOptions() {
		Options options = new Options();
		Builder outputBuilder = Option.builder("o");
		outputBuilder.longOpt("output");
		outputBuilder.hasArg();
		outputBuilder.argName("format");
		StringBuilder outputOptionsDesc = new StringBuilder();
		outputOptionsDesc.append("Sets OPSIN's output format (default smi)").append(OpsinTools.NEWLINE);
		outputOptionsDesc.append("Allowed values are:").append(OpsinTools.NEWLINE);
		outputOptionsDesc.append("cml for Chemical Markup Language").append(OpsinTools.NEWLINE);
		outputOptionsDesc.append("smi for SMILES").append(OpsinTools.NEWLINE);
		outputOptionsDesc.append("extendedsmi for Extended SMILES").append(OpsinTools.NEWLINE);
		outputOptionsDesc.append("inchi for InChI (with FixedH)").append(OpsinTools.NEWLINE);
		outputOptionsDesc.append("stdinchi for StdInChI").append(OpsinTools.NEWLINE);
		outputOptionsDesc.append("stdinchikey for StdInChIKey");
		outputBuilder.desc(outputOptionsDesc.toString());
		options.addOption(outputBuilder.build());
		options.addOption("h", "help", false, "Displays the allowed command line flags");
		options.addOption("v", "verbose", false, "Enables debugging");
		
		options.addOption("a", "allowAcidsWithoutAcid", false, "Allows interpretation of acids without the word acid e.g. \"acetic\"");
		options.addOption("f", "detailedFailureAnalysis", false, "Enables reverse parsing to more accurately determine why parsing failed");
		options.addOption("n", "name", false, "Include name in SMILES/InChI output (tab delimited)");
		options.addOption("r", "allowRadicals", false, "Enables interpretation of radicals");
		options.addOption("s", "allowUninterpretableStereo", false, "Allows stereochemistry uninterpretable by OPSIN to be ignored");
		options.addOption("w", "wildcardRadicals", false, "Radicals are output as wildcard atoms");
		return options;
	}
	
	/**
	 * Uses the command line parameters to configure a new NameToStructureConfig
	 * @param cmd
	 * @return The configured NameToStructureConfig
	 */
	private static NameToStructureConfig generateOpsinConfigObjectFromCmd(CommandLine cmd) {
		NameToStructureConfig n2sconfig = new NameToStructureConfig();
		n2sconfig.setInterpretAcidsWithoutTheWordAcid(cmd.hasOption("a"));
		n2sconfig.setDetailedFailureAnalysis(cmd.hasOption("f"));
		n2sconfig.setAllowRadicals(cmd.hasOption("r"));
		n2sconfig.setWarnRatherThanFailOnUninterpretableStereochemistry(cmd.hasOption("s"));
		n2sconfig.setOutputRadicalsAsWildCardAtoms(cmd.hasOption("w"));
		return n2sconfig;
	}

	private static void interactiveCmlOutput(InputStream input, OutputStream out, NameToStructureConfig n2sconfig) throws IOException, XMLStreamException {
		NameToStructure nts = NameToStructure.getInstance();
		BufferedReader inputReader = new BufferedReader(new InputStreamReader(input, "UTF-8"));
		XMLOutputFactory factory = new WstxOutputFactory();
		factory.setProperty(WstxOutputProperties.P_OUTPUT_ESCAPE_CR, false);
		XMLStreamWriter writer = factory.createXMLStreamWriter(out, "UTF-8");
		writer = new IndentingXMLStreamWriter(writer, 2);
		writer.writeStartDocument();
		CMLWriter cmlWriter = new CMLWriter(writer);
		cmlWriter.writeCmlStart();
		int id = 1;
		String line;
		while((line =inputReader.readLine()) != null) {
			int splitPoint = line.indexOf('\t');
			String name = splitPoint >=0 ? line.substring(0, splitPoint) : line;
			OpsinResult result = nts.parseChemicalName(name, n2sconfig);
			Fragment structure = result.getStructure();
			cmlWriter.writeMolecule(structure, name, id++);
			writer.flush();
			if(structure == null) {
				System.err.println(result.getMessage());
			}
		}
		cmlWriter.writeCmlEnd();
		writer.writeEndDocument();
		writer.flush();
		writer.close();
	}
	
	private static void interactiveSmilesOutput(InputStream input, OutputStream out, NameToStructureConfig n2sconfig, boolean extendedSmiles, boolean outputName) throws IOException {
		NameToStructure nts = NameToStructure.getInstance();
		BufferedReader inputReader = new BufferedReader(new InputStreamReader(input, "UTF-8"));
		BufferedWriter outputWriter = new BufferedWriter(new OutputStreamWriter(out, "UTF-8"));
		String line;
		while((line =inputReader.readLine()) != null) {
			int splitPoint = line.indexOf('\t');
			String name = splitPoint >=0 ? line.substring(0, splitPoint) : line;
			OpsinResult result = nts.parseChemicalName(name, n2sconfig);
			String output = extendedSmiles ? result.getExtendedSmiles() : result.getSmiles();
			if(output == null) {
				System.err.println(result.getMessage());
			} else {
				outputWriter.write(output);
			}
			if (outputName) {
				outputWriter.write('\t');
				outputWriter.write(line);
			}
			outputWriter.newLine();
			outputWriter.flush();
		}
	}

	private static void interactiveInchiOutput(InputStream input, OutputStream out, NameToStructureConfig n2sconfig, InchiType inchiType, boolean outputName) throws Exception {
		NameToStructure nts = NameToStructure.getInstance();
		BufferedReader inputReader = new BufferedReader(new InputStreamReader(input, "UTF-8"));
		BufferedWriter outputWriter = new BufferedWriter(new OutputStreamWriter(out, "UTF-8"));
		Class<?> c;
		try {
			c = Class.forName("uk.ac.cam.ch.wwmm.opsin.NameToInchi");
		} catch (ClassNotFoundException e) {
			System.err.println("Could not initialise NameToInChI module. Is it on your classpath?");
			throw new RuntimeException(e);
		}
		Method m;
		switch (inchiType) {
		case inchiWithFixedH:
			m = c.getMethod("convertResultToInChI", new Class[]{OpsinResult.class});
			break;
		case stdInchi:
			m = c.getMethod("convertResultToStdInChI", new Class[]{OpsinResult.class});
			break;
		case stdInchiKey:
			m = c.getMethod("convertResultToStdInChIKey", new Class[]{OpsinResult.class});
			break;
		default :
			throw new IllegalArgumentException("Unexepected enum value: " + inchiType);
		}

		String line;
		while((line =inputReader.readLine()) != null) {
			int splitPoint = line.indexOf('\t');
			String name = splitPoint >=0 ? line.substring(0, splitPoint) : line;
			OpsinResult result = nts.parseChemicalName(name, n2sconfig);
			String output = (String) m.invoke(null, result);
			if(output == null) {
				System.err.println(result.getMessage());
			} else {
				outputWriter.write(output);
			}
			if (outputName) {
				outputWriter.write('\t');
				outputWriter.write(line);
			}
			outputWriter.newLine();
			outputWriter.flush();
		}
	}
}
