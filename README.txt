OPSIN - Open Parser for Systematic IUPAC Nomenclature
version 1.3.0 (see ReleaseNotes.txt for what's new in this version)

Daniel Lowe(Current maintainer), Dr. Peter Corbett and Prof. Peter Murray-Rust
We are thankful for contributions from Albina Asadulina and Rich Apodaca

Contact address: dl387@cam.ac.uk
Source code: http://bitbucket.org/dan2097/opsin/
Web interface and informational site: http://opsin.ch.cam.ac.uk/
JavaDoc is available via OPSIN's generated web site: http://apidoc.ch.cam.ac.uk/opsin/

This is a Java(1.5+) library for IUPAC name-to-structure conversion offering high recall and precision on general organic nomenclature.

##################################################

OPSIN is available as a standalone JAR from Bitbucket, http://bitbucket.org/dan2097/opsin/downloads
It is also available as a dependency for use with Apache Maven.
opsin-1.3.0-jar-with-dependencies.jar includes CML (Chemical Markup Language), SMILES, and InChI (IUPAC International Chemical Identifier) output and all dependendencies
The main classes are uk.ac.cam.ch.wwmm.opsin.NameToStructure for CML and SMILES and uk.ac.cam.ch.wwmm.opsin.NameToInchi for InChI.

To use OPSIN as a library add opsin-1.3.0-jar-with-dependencies.jar to your classpath.

If you are using Maven then add:
		<dependency>
			 <groupId>uk.ac.cam.ch.opsin</groupId>
			 <artifactId>opsin-core</artifactId>
			 <version>1.3.0</version>
		</dependency>
	If you need just CML or SMILES output support

	or
		<dependency>
			 <groupId>uk.ac.cam.ch.opsin</groupId>
			 <artifactId>opsin-inchi</artifactId>
			 <version>1.3.0</version>
		</dependency>

	if you also need InChI output support.

##################################################

Using OPSIN as a command-line utility:
	java -jar opsin-1.3.0-jar-with-dependencies.jar will give you a command-line interface to convert names to CML
	java -jar opsin-1.3.0-jar-with-dependencies.jar -h will give you information on more advanced usage

As well as interactive input on the command-line the command-line interface will accept a piped input of newline seperated chemical names.
e.g. java -jar opsin-1.3.0-jar-with-dependencies.jar < input.name > output.cml
e.g. java -jar opsin-1.3.0-jar-with-dependencies.jar -osmi < input.name > output.smiles

Using OPSIN as a library:

1) Create an OPSIN instance, by calling the following static method

NameToStructure n2s = NameToStructure.getInstance();

2) (Optional) create a NameToStructureConfig (this can be set to, for example, allow radicals to convert to structures)
NameToStructureConfig n2sconfig = new NameToStructureConfig();

3) Give OPSIN a name and recieve an OpsinResult
OpsinResult result = n2s.parseChemicalName("acetonitrile", n2sconfig);

4) Check whether conversion was successful by checking result.getStatus()

5) Retrieve the structure as CML/SMILES/InChI:
Element cml = result.getCml();
String smiles = result.getSmiles();
String inchi = NameToInchi.convertResultToInChI(result);

Convenience methods exist to convert directly from names to CML, SMILES and InChI

Converting a name to a structure will typically take of the order of 5-10ms making OPSIN suitable for batch processing of large sets of names.

##################################################

The workings of OPSIN are more fully described in:

Chemical Name to Structure: OPSIN, an Open Source Solution
Daniel M. Lowe, Peter T. Corbett, Peter Murray-Rust, Robert C. Glen
Journal of Chemical Information and Modeling 2011 51 (3), 739-753

If you use OPSIN to produce results for publication, then it would be great if you could cite us.

The following list broadly summarises what OPSIN can currently do and what will be worked on in the future.

Supported nomenclature includes:
alkanes/alkenes/alkynes/heteroatom chains e.g. hexane, hex-1-ene, tetrasiloxane and their cyclic analogues e.g. cyclopropane
All IUPAC 1993 recommended rings
Trivial acids
Hantzsch-Widman e.g. 1,3-oxazole
Spiro systems
All von Baeyer rings e.g. bicyclo[2.2.2]octane
Hydro/dehydro e.g. 2,3-dihydropyridine
Indicated hydrogen e.g. 1H-benzoimidazole
Heteroatom replacement
Specification of charge e.g. ium/ide/ylium/uide
Multiplicative nomenclature e.g. ethylenediaminetetraacetic acid
Conjunctive nomenclature e.g. cyclohexaneethanol
Fused ring systems e.g. imidazo[4,5-d]pyridine
Ring assemblies e.g. biphenyl
Most prefix and infix functional replacement nomenclature
The following functional classes: acids, acetals, alcohols, amides, anhydrides, azides, bromides, chlorides, cyanates, cyanides,
esters, di/tri/tetra esters, ethers, fluorides, fulminates, glycols, glycol ethers, hemiacetals, hemiketal, hydrazones, hydroperoxides,
hydrazides, imides, iodides, isocyanates, isocyanides, isoselenocyanates, isothiocyanates, ketals, ketones, lactams, lactims,
lactones, selenocyanates, thiocyanates, selenols, thiols, mercaptans, oxides, oximes, peroxides, selenides, selenones, selenoxides,
selones, selenoketones, selenosemicarbazones, semicarbazones, sulfides, sulfones, sulfoxides, sultams, sultims, sultines, sultones,
tellurides, telluroketones, tellurosemicarbazones, tellurones, telluroxides, thioketones and thiosemicarbazones
Greek letters
Lambda convention
E/Z/R/S stereochemistry
cis/trans indicating relative stereochemistry on rings and as a synonym of E/Z
Amino Acids and derivatives
Structure-based polymer names e.g. poly(2,2'-diamino-5-hexadecylbiphenyl-3,3'-diyl)
Simple bridge prefixes e.g. methano
Nucleosides, nucleotides and their esters
Steroids including alpha/beta stereochemistry
Specification of oxidation numbers and charge on elements
Perhalogeno terms
Deoxy
Open-chain and simple cyclised carbohydrates
Stoichiometry ratios and mixture indicators
Simple CAS names including inverted CAS names

Currently UNsupported nomenclature includes:
Other less common stereochemical terms
Carbohydrate derivatives e.g. glycosides
Most natural Products other than steroids
Natural product specific nomenclature operations
Multiplied, unsaturated or composite bridge prefixes e.g. epoxymethano

##################################################

Good Luck and let us know if you have problems, comments or suggestions!
You can contact us by posting a message on Bitbucket or you can email me directly (dl387@cam.ac.uk)