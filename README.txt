OPSIN - Open Parser for Systematic IUPAC Nomenclature
version 0.6.1 (see ReleaseNotes.txt for what's new in this version)

Daniel Lowe(Current maintainer), Dr. Peter Corbett and Prof. Peter Murray-Rust

Contact address: dl387@cam.ac.uk

This is a library for IUPAC name-to-structure conversion.
Currently it should be considered to be under development although the interface for using it will remain constant.
OPSIN was formerly a component of OSCAR3 but is now a wholly standalone library.

##################################################

The easiest way to use OPSIN is to use the standalone jar available from sourceforge.
java -jar opsin-0.6.1-jar-with-dependencies.jar will give you a command line interface to convert names to CML (Chemical Markup Language)
opsin-0.6.1-jar-with-dependencies.jar includes InChI and CML output and all dependendencies
core-0.6.1.jar does not include dependencies and InChI output support
inchi-0.6.1.jar does not include dependencies, depends on JNI-InChI and core, and adds InChI suppport.
The main classes are uk.ac.cam.ch.wwmm.opsin.NameToStructure for CML
and uk.ac.cam.ch.wwmm.opsin.NameToInchi for InChI

To use OPSIN as a library within Java add opsin-0.6.1.jar to your classpath then:

1) Learn about XOM (http://xom.nu), the XML processing framework used
   by OPSIN
2) Create an OPSIN instance, by calling the following static method

NameToStructure nameToStructure = NameToStructure.getInstance();

3) Get CML (as XOM Elements):

Element cmlElement = nameToStructure.parseToCML("acetonitrile");

4) Whatever you like. Maybe print it out, thus:

System.out.println(cmlElement.toXML());

parseToCML will typically take 5-10ms to convert a name to CML making OPSIN suitable for use on a large number of names.

CML can, if desired, be converted to other format such as SD, SMILES, InChI etc. by toolkits such as CDK, OpenBabel and JUMBO.
(NOTE: if you want InChI the most efficient way to generate it is to use the InChI module and the corresponding parseToInChI method)

##################################################

The workings of OPSIN are more fully described in:

Peter Corbett, Peter Murray-Rust High-throughput identification of
chemistry in life science texts. Proceedings of Computational Life
Sciences (CompLife) 2006, Cambridge, UK, pp. 107-118.

The following lists broadly summarise what OPSIN can currently do and what will be worked on in the future.

Supported nomenclature includes:
alkanes/alkenes/alkynes/heteroatom chains e.g. hexane, hex-1-ene, tetrasiloxane and their cyclic analogues e.g. cyclopropane
All IUPAC 1993 recommended rings
Trivial acids
Hantzsch-Widman e.g. 1,3-oxazole
Spiro systems (using Von baeyer brackets)
All von Baeyer rings e.g. bicyclo[2.2.2]octane
Hydro e.g. 2,3-dihydropyridine
Indicated hydrogen e.g. 1H-benzoimidazole
Heteroatom replacement
Specification of charge e.g. ium/ide
Multiplicative nomenclature e.g. ethylenediaminetetraacetic acid
Fused ring systems e.g. imidazo[4,5-d]pyridine. For a small number of fused ring systems numbering cannot be determined
Ring assemblies e.g. biphenyl
Most prefix and infix functional replacement nomenclature
The following functional classes: amide, anhydrides, esters, diesters, glycols, acids, azides, bromides, chlorides, cyanates, cyanides, fluorides, fulminates, hydrazones, hydroperoxides, imides,
iodides, isocyanates, isocyanides, isoselenocyanates, isothiocyanates, selenocyanates, thiocyanates, alcohols, selenols, thiols, ethers, ketones, oxides, oximes, peroxides, selenides, 
selenones, selenoxides, selones, selenoketones, selenosemicarbazone, semicarbazones, sulfides, sulfones, sulfoxides, tellurides, telluroketones, tellurosemicarbazones, tellurones, telluroxides, thioketones and thiosemicarbazones
Greek letters
Lambda convention
E/Z/R/S stereochemistry
Amino Acids and derivatives

Currently UNsupported nomenclature includes:
Other less common stereochemical terms
Carbohydrates
Steroids
Nucleic acids
Bridged rings 
Fused ring systems involving non 6-membered rings which are not in a "chain" cannot be numbered e.g. indeno[2,1-c]pyridine can be numbered, benzo[cd]indole cannot
Some conjunctive operations e.g. cyclohexaneethanol
The following functional classes: Hydrazides, lactones, sultams, lactams, sultims, lactims, acetals, hemiacetals and ketals

##################################################

Good Luck and let us know if you have problems, comments or suggestions!
You can contact us by posting a message on SourceForge or you can email me directly (dl387@cam.ac.uk)