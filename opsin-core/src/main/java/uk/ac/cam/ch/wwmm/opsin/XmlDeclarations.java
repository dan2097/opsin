package uk.ac.cam.ch.wwmm.opsin;

/**
 * Contains static final strings corresponding to XML element names and attributes employed by OPSIN
 * @author dl387
 *
 */
class XmlDeclarations {
	
	//TODO are all these types and subtypes actually a good idea considering the vast majority are never used?
	
	/*
	 * The container XML elements. These are generated by OPSIN 
	 */
	/**Define a scope for determining what group a substituent should bond to*/
	static final String BRACKET_EL ="bracket";

	/**Contains a functional group or class. These terms typically effect the chosen wordRule for the name*/
	static final String FUNCTIONALTERM_EL ="functionalTerm";
	
	/**The top most element in OPSIN's parse tree. As a name can describe multiple molecules the same is confusingly true of this element*/
	static final String MOLECULE_EL ="molecule";
	
	/**Contains a substituent. A substituent will after the ComponentProcessor contain one group*/
	static final String SUBSTITUENT_EL = "substituent";

	/**Contains a root group(the rightmost in a word). A root will after the ComponentProcessor contain one group*/
	static final String ROOT_EL ="root";

	/**Contains brackets/substituents/root. Generally these correspond to words in the original chemical name (unless erroneous/omitted spaces were present)*/
	static final String WORD_EL ="word";
	
	/**Contains words/wordRules. The value of the wordRule indicates how the StructureBuilder should process its children*/
	static final String WORDRULE_EL ="wordRule";
	

	/*
	 * The token XML elements. These are generally produced by the parser from the tokenised chemical name
	 * Some are produced by OPSIN in the ComponentGenerator/ComponentProcessor
	 */
	
	/**Adds a hydrogen to an unsaturated system, this is hydrogen that is added due to a suffix and is expressed in a locant e.g. 1(2H) */
	static final String ADDEDHYDROGEN_EL ="addedHydrogen";
	
	/**A component of an alkaneStem e.g. [octa][hexaconta][tetract]ane will have three alkaneStemComponents*/
	static final String ALKANESTEMCOMPONENT ="alkaneStemComponent";
	
	/**Something like tert/iso/sec Modifies an alkaneStem in the ComponentGenerator*/
	static final String ALKANESTEMMODIFIER_EL ="alkaneStemModifier";

	/**An annulene. Converted to a group by the ComponentGenerator*/
	static final String ANNULEN_EL ="annulen";
	
	/**A bridge described in SMILES for used on rings*/
	static final String FUSEDRINGBRIDGE_EL ="fusedRingBridge";
	
	/**An O that indicates that the preceding alkaneStem is in fact a bridge*/
	static final String BRIDGEFORMINGO_EL ="bridgeFormingO";

	/**A locant indicating the positions for a glycosidic linkage. The first locant will point to an alpha carbon
	 * Also used to indicate joining of nucleosyl groups*/
	static final String BIOCHEMICALLINKAGE_EL ="biochemicalLinkage";	

	/**Indicates the size of the ring in a carbohydrate e.g. furanose = 5*/
	static final String CARBOHYDRATERINGSIZE_EL ="carbohydrateRingSize";

	/**A charge specifier e.g. (2+). Value is the charge to set something to*/
	static final String CHARGESPECIFIER_EL ="chargeSpecifier";
	
	/**Created by the ComponentProcessor. Something like the acetic acid in benzene-1,3,5-triacetic acid*/
	static final String CONJUNCTIVESUFFIXGROUP_EL ="conjunctiveSuffixGroup";
	
	/**Used by the ComponentGenerator to group elements into bracket elements*/
	static final String CLOSEBRACKET_EL ="closebracket";
	
	/**Used by the ComponentGenerator to modify alkanes into cycloalkanes*/
	static final String CYCLO_EL ="cyclo";
	
	/** A delta used to indicate the position of a double bond in older nomenclature*/
	static final String DELTA_EL ="delta";
	
	/** Used in amino acid and carbohydrate nomenclature to indicate stereochemistry*/
	static final String DLSTEREOCHEMISTRY_EL ="dlStereochemistry";
	
	/**A fractional multiplier e.g. hemi*/
	static final String FRACTIONALMULTIPLIER_EL ="fractionalMultiplier";

	/**A functional Class such as acid. Does not correspond to a fragment*/
	static final String FUNCTIONALCLASS_EL ="functionalClass";
	
	/**A functional group such as alcohol or sulfone. Describes a fragment*/
	static final String FUNCTIONALGROUP_EL ="functionalGroup";
	
	/**Currently just poly or oligo for polymers*/
	static final String FUNCTIONALMODIFIER_EL ="functionalModifier";

	/**A fusion bracket. Used in fusion nomenclature*/
	static final String FUSION_EL ="fusion";

	/**Define a scope for determining what group a substituent should bond to*/
	static final String GROUP_EL ="group";
	
	/**A heteroatom. Could be part of a Hantzsch Widman ring or a replacement prefix*/
	static final String HETEROATOM_EL ="heteroatom";

	/**Adds a hydrogen to an unsaturated system (hydro/perhydro)*/
	static final String HYDRO_EL ="hydro";

	/**One of the systematic hydrocarbon fused ring series e.g. tetralene, pentalene. Converted to a group by the ComponentGenerator*/
	static final String HYDROCARBONFUSEDRINGSYSTEM_EL ="hydrocarbonFusedRingSystem";

	/**Adds a hydrogen to an unsaturated system to indicate what atoms are saturated in a system where not all atoms with spare valency can form double bonds e.g.  e.g. 2H-pyran*/
	static final String INDICATEDHYDROGEN_EL ="indicatedHydrogen";
	
	/**Was the word salt encountered indicating that a salt was expected? */
	static final String ISSALT_ATR = "isSalt";
	
	/**A hyphen between two substituents. Used as hint that the two substituents do not join together*/
	static final String HYPHEN_EL ="hyphen";
	
	/**ine as in the end of an aminoAcid. Has no meaning*/
	static final String INE_EL ="ine";
	
	/**An infix. This performs functionalReplacement on a suffix*/
	static final String INFIX_EL ="infix";

	/**Indicates that a heteroatom or atom should be in a specific valency*/
	static final String LAMBDACONVENTION_EL ="lambdaConvention";
	
	/**A locant e.g. where a substituent should attach*/
	static final String LOCANT_EL ="locant";
	
	/**Used by the ComponentGenerator to group elements into bracket elements*/
	static final String OPENBRACKET_EL ="openbracket";

	/**otho/meta/para Converted to a locant by the ComponentProcessor*/
	static final String ORTHOMETAPARA_EL ="orthoMetaPara";

	/**Describes the number of spiro centres in a poly cyclic spiro system*/
	static final String POLYCYCLICSPIRO_EL ="polyCyclicSpiro";
	
	/**A locant indicating through which atoms a multiplied parent in multiplicative nomenclature is connected*/
	static final String MULTIPLICATIVELOCANT_EL ="multiplicativeLocant";
	
	/**A multiplier e.g. indicating multiplication of a heteroatom or substituent*/
	static final String MULTIPLIER_EL ="multiplier";
	
	/**e.g. (III), Specifies the oxidation number of an atom. Value is the oxidation number to set something to*/
	static final String OXIDATIONNUMBERSPECIFIER_EL ="oxidationNumberSpecifier";

	/**Used amongst other things to indicate how the rings of a ring assembly should be assembled*/
	static final String COLONORSEMICOLONDELIMITEDLOCANT_EL ="colonOrSemiColonDelimitedLocant";

	/**Used to indicate how many rings are in a ring assembly*/
	static final String RINGASSEMBLYMULTIPLIER_EL ="ringAssemblyMultiplier";
	
	/**A spiro system. Converted to a group by the ComponentGenerator*/
	static final String SPIRO_EL ="spiro";
	
	/**A locant that seperates components of a spiro system*/
	static final String SPIROLOCANT_EL ="spiroLocant";
	
	/**Something like R/S/E/Z. Indicates stereochemical configuration*/
	static final String STEREOCHEMISTRY_EL ="stereoChemistry";

	/**Present in complicated nomenclature e.g. ring assemblies to avoid ambiguity*/
	static final String STRUCTURALCLOSEBRACKET_EL ="structuralCloseBracket";

	/**Present in complicated nomenclature to avoid ambiguity*/
	static final String STRUCTURALOPENBRACKET_EL ="structuralOpenBracket";
	
	/**Indicates replacement of a group by hydrogen e.g. deoxy means replace OH with H*/
	static final String SUBTRACTIVEPREFIX_EL ="subtractivePrefix";
	
	/**A suffix e.g. amide, al, yl etc.*/
	static final String SUFFIX_EL ="suffix";
	
	/**Something like sulfon/carbo/carbox that modifies a following suffix*/
	static final String SUFFIXPREFIX_EL ="suffixPrefix";
	
	/**ene/yne, indicated that a double/triple bond should be formed at a saturated location*/
	static final String UNSATURATOR_EL ="unsaturator";
	
	/**A vonBaeyer system. Converted to a group by the ComponentGenerator*/
	static final String VONBAEYER_EL ="vonBaeyer";

	/*
	 * The token XML attributes. These are generally produced by the parser from the tokenised chemical name
	 * Some are produced by OPSIN in the ComponentGenerator/ComponentProcessor
	 */

	static final String VALUE_ATR ="value";
	static final String LABELS_ATR = "labels";
	static final String FUSEDRINGNUMBERING_ATR = "fusedRingNumbering";
	static final String DEFAULTINLOCANT_ATR = "defaultInLocant";
	static final String DEFAULTINID_ATR = "defaultInID";
	static final String OUTIDS_ATR = "outIDs";
	
	static final String ALPHABETACLOCKWISEATOMORDERING_ATR="alphaBetaClockWiseAtomOrdering";	
	static final String ACCEPTSADDITIVEBONDS_ATR = "acceptsAdditiveBonds";
	
	/**Works like a locant but refers to the atoms OPSIN id. Will be overridden by the locant/locantId attribute*/
    static final String DEFAULTLOCANTID_ATR = "defaultLocantID";

	static final String IMINOLIKE_ATR = "iminoLike";
	
	/**The functional replacement specified by an infix to be performed on the suffix*/
	static final String INFIX_ATR = "infix";

	/**Indicates that an element has been multiplied. Prevents badly assigning indirect locants*/
	static final String MULTIPLIED_ATR = "multiplied";

	/**A comma separated list of relative ids at which to add functionalAtoms*/
	static final String FUNCTIONALIDS_ATR = "functionalIDs";
	static final String ADDGROUP_ATR = "addGroup";
	static final String ADDHETEROATOM_ATR = "addHeteroAtom";
	static final String ADDBOND_ATR = "addBond";
	
	/**Can the substituent be implicitly bracketted to a previous substitutent e.g. methylaminobenzene --> (methylamino)benzene as amino has this atr*/
    static final String USABLEASJOINER_ATR = "usableAsAJoiner";
    
    /**A comma separated list of locants that are expected in front of a group for either xylene-like nomenclature or as indirect locants*/
    static final String FRONTLOCANTSEXPECTED_ATR = "frontLocantsExpected";

    /** Used as a fudge for some hydrogen esters e.g. dihydrogenphosphate*/
    static final String NUMBEROFFUNCTIONALATOMSTOREMOVE_ATR = "numberOfFunctionalAtomsToRemove";
    
    /**A comma seperated list of relatives ids indicating where to add suffix/es*/
    static final String SUFFIXAPPLIESTO_ATR = "suffixAppliesTo";

    /**A relatives id indicating at what position to attach a suffix to by default*/
    static final String SUFFIXAPPLIESTOBYDEFAULT_ATR = "suffixAppliesToByDefault";
    static final String COMMONOXIDATIONSTATESANDMAX_ATR = "commonOxidationStatesAndMax";
    static final String ADDITIONALVALUE_ATR = "additionalValue";
    static final String LOCANT_ATR = "locant";

	/**Works like a locant but refers to the atoms OPSIN id*/
    static final String LOCANTID_ATR = "locantID";
    
    
    static final String TYPE_ATR = "type";
    static final String SUBTYPE_ATR = "subType";

	/**Defines the locants for which a radical will connect to another group in multiplicative nomenclature e.g. in 2,2'-methylenedipyridine the 2,2' become inlocants of the pyridine*/
    static final String INLOCANTS_ATR = "inLocants"; 
    
	/**Determined by the {@link ComponentProcessor}. True if a fragment has more than two radical positions e.g. ethan-1,2-diyl not ethanylidene*/
    static final String ISAMULTIRADICAL_ATR = "isAMultiRadical"; 
    
	/**Added to a heteroatom or LAMBDACONVENTION_EL to indicate the desired valency*/
    static final String LAMBDA_ATR = "lambda";
    
	/**Indicates how many times a bracket/substituent should be multiplied*/
	static final String MULTIPLIER_ATR ="multiplier";
	
	/** The name that was inputted into OPSIN's parser. Attribute of molecule */
    static final String NAME_ATR = "name";
	
	/**Indicates that a substituent/bracket has been processed by StructureBuildingMethods*/
	static final String RESOLVED_ATR ="resolved";
	
	/**Indicates that the natural enantiomer is the opposite that is expected for the class of compound e.g. a natural L sugar*/
	static final String NATURALENTISOPPOSITE_ATR ="naturalEntIsOpposite";
	
	/**Placed on a word rule if explicit stoichiometry has been provided. Value is always an integer */
    static final String STOICHIOMETRY_ATR = "stoichiometry";
    
	/** Holds the value of any tokens for which XML was not generated by the parser e.g. an optional e. Multiple elided tokens will be concatenated*/
	static final String SUBSEQUENTUNSEMANTICTOKEN_ATR ="subsequentUnsemanticToken";
    
	/**Added by the ComponentGenerator to a suffix*/
    static final String SUFFIXPREFIX_ATR = "suffixPrefix";
    
	/**The wordRule that a wordRule element corresponds to*/
	static final String WORDRULE_ATR ="wordRule";

	/*
	 * The values the type attribute can take
	 * Type is expected to be present at minimum on all group elements
	 */
	/**A term like amide or hydrazide that replaces a functional hydroxy group*/
	static final String ACIDREPLACINGFUNCTIONALGROUP_TYPE_VAL ="acidReplacingFunctionalGroup";
	
	/**A trivial carboxylic acid. These by default do not have their acid groups which are then added on using suffixes*/
	static final String ACIDSTEM_TYPE_VAL ="acidStem";
	
	/**This stereochemistry element conveys alpha/beta stereochemistry*/
	static final String ALPHA_OR_BETA_TYPE_VAL ="alphaOrBeta";
	
	/**An aminoAcid. These by default do not have their acid groups which are then added on using suffixes. Notably these suffixes do NOT correspond to tokens in the input chemical name!*/
	static final String AMINOACID_TYPE_VAL ="aminoAcid";
	
	/**A subtractive prefix that removes a terminal chalcogen and forms an intramolecular bridge to another*/
	static final String ANHYDRO_TYPE_VAL ="anhydro";
	
	/**This stereochemistry element conveys axial stereochemistry
	 * These indicate the postion of groups are an axis/plane/helix. This is expressed by the descriptors: M, P, Ra, Sa, Rp, Sp*/
	static final String AXIAL_TYPE_VAL ="axial";
    
	/**A normal multiplier e.g. di*/
	static final String BASIC_TYPE_VAL ="basic";
	
	/**A locant enclosed in square brackets e.g. [5]*/
	static final String BRACKETEDLOCANT_TYPE_VAL ="bracketedLocant";
	
	/**This stereochemistry element specifies stereochemistry in a carbohydrate e.g. gluco is  r/l/r/r (position of hydroxy in a fischer projection)*/
	static final String CARBOHYDRATECONFIGURATIONPREFIX_TYPE_VAL ="carbohydrateConfigurationalPrefix";
	
	/**Groups formed in accordance with carbohydrate nomenclature */
	static final String CARBOHYDRATE_TYPE_VAL ="carbohydrate";
	
	/**Indicates the group should be acyclic*/
	static final String CHAIN_TYPE_VAL ="chain";
	
	/**This suffix modifies charge*/
	static final String CHARGE_TYPE_VAL ="charge";
	
	/**This stereochemistry element conveys cis/trans stereochemistry*/
	static final String CISORTRANS_TYPE_VAL ="cisOrTrans";
	
	/**This stereochemistry element conveys R/S stereochemistry*/
	static final String R_OR_S_TYPE_VAL ="RorS";
	
	/**This stereochemistry element conveys E/Z stereochemistry*/
	static final String E_OR_Z_TYPE_VAL ="EorZ";
	
	/**This group is a sulfur/selenium/tellurium acid with the acidic hydroxy missing*/
	static final String CHALCOGENACIDSTEM_TYPE_VAL ="chalcogenAcidStem";
	
	/**A subtractive prefix that removes a hydrogen to covert a hydroxy into a carbonyl or convert a bond to a double/triple bond*/
	static final String DEHYDRO_TYPE_VAL ="dehydro";
	
	/**A subtractive prefix that removes a terminal hydroxy like atom*/
	static final String DEOXY_TYPE_VAL ="deoxy";
	
	/**A functional group describing a divalent group*/
	static final String DIVALENTGROUP_TYPE_VAL ="diValentGroup";
	
	/**This stereochemistry element conveys endo/exo/syn/anti stereochemistry
	 * These indicate relative orientation of groups attached to non-bridgehead atoms in a bicyclo[x.y.z]alkane (x >= y > z > 0)*/
	static final String ENDO_EXO_SYN_ANTI_TYPE_VAL ="endoExoSynAnti";
	
	/**A group that is functional class e.g. O for anhydride*/
	static final String FUNCTIONALCLASS_TYPE_VAL ="functionalClass";
	
	/**A multiplier for groups of terms e.g. bis*/
	static final String GROUP_TYPE_VAL ="group";
	
	/**An implicit bracket. Implicit brackets are added where a bracket is needed to give the intended meaning*/
	static final String IMPLICIT_TYPE_VAL ="implicit";
	
	/**This suffix adds a radical to the preceding group e.g. yl, oyl*/
	static final String INLINE_TYPE_VAL ="inline";
	
	/**This functional group is monovalent e.g. alcohol*/
	static final String MONOVALENTGROUP_TYPE_VAL ="monoValentGroup";
	
	/**This functional group is monovalent and describes a specific compound e.g. cyanide*/
	static final String MONOVALENTSTANDALONEGROUP_TYPE_VAL ="monoValentStandaloneGroup";
	
	/**A non carboxylic acid e.g. phosphoric*/
	static final String NONCARBOXYLICACID_TYPE_VAL ="nonCarboxylicAcid";
	
	/**This stereochemistry element describes the direction that plane polarised light is rotated*/
	static final String OPTICALROTATION_TYPE_VAL ="opticalRotation";

	/**Indicates the locant was made from an ortho/meta/para term*/
	static final String ORTHOMETAPARA_TYPE_VAL ="orthoMetaPara";
	
	/**This stereochemistry element conveys relative cis/trans stereochemistry e.g. r-1, c-2, t-3*/
	static final String RELATIVECISTRANS_TYPE_VAL ="relativeCisTrans";
	
	/**Indicates the group should be, at least in part, cyclic*/
	static final String RING_TYPE_VAL ="ring";
	
	/**Indicates a group that does not allow suffixes*/
	static final String SIMPLEGROUP_TYPE_VAL ="simpleGroup";
	
	/**Groups that do not have any special rules for suffix handling*/
	static final String STANDARDGROUP_TYPE_VAL ="standardGroup";
	
	/**A bracket containing R/S/E/Z descriptors*/
	static final String STEREOCHEMISTRYBRACKET_TYPE_VAL ="stereochemistryBracket";
	
	/**Indicates a group that is a substituent*/
	static final String SUBSTITUENT_TYPE_VAL ="substituent";
	
	/**A locant that also indicated the addition of hydrogen e.g.2(1H); not used to locant onto another group*/
	static final String ADDEDHYDROGENLOCANT_TYPE_VAL ="addedHydrogenLocant";
	
	/**Indicates a group that is a suffix*/
	static final String SUFFIX_TYPE_VAL ="suffix";
	
	/**A suffix that does not add a radical, hence will be present only on the root group */
	static final String ROOT_TYPE_VAL ="root";
	
	/**A multiplier for a Von Baeyer system e.g. bi in bicyclo*/
	static final String VONBAEYER_TYPE_VAL ="VonBaeyer";
	
	
	/*
	 * The values the subType attribute can take
	 * subType is expected to be present at minimum on all group elements
	 */
	
	/**The stem of an alkane e.g. "eth" */
	static final String ALKANESTEM_SUBTYPE_VAL ="alkaneStem";
	/**An anhydride functional term e.g. "thioanhydride"*/
	static final String ANHYDRIDE_SUBTYPE_VAL ="anhydride";
	/**An aryl group e.g. "pyridin" */
	static final String ARYLGROUP_SUBTYPE_VAL ="arylGroup";
	/**An aryl subsituent or stem e.g. "phenyl", "styr" */
	static final String ARYLSUBSTITUENT_SUBTYPE_VAL ="arylSubstituent";
	//FIXME ideally carbohydrates and nucleotides/nucleosides/natural products should have a common type or subtype
	/**Nucleotides/nucleosides/natural products. Carbohydrates can be detected by {@link XmlDeclarations#CARBOHYDRATE_TYPE_VAL}*/
	static final String BIOCHEMICAL_SUBTYPE_VAL ="biochemical";
	/**A trivial carbohydrate stem for an aldose e.g. "galact"*/
	static final String CARBOHYDRATESTEMALDOSE_SUBTYPE_VAL ="carbohydrateStemAldose";
	/**A trivial carbohydrate stem for a ketose e.g. "fruct"*/
	static final String CARBOHYDRATESTEMKETOSE_SUBTYPE_VAL ="carbohydrateStemKetose";
	/**A suffix that forms a cycle e.g. imide, lactam, sultam*/
	static final String CYCLEFORMER_SUBTYPE_VAL ="cycleformer";
	/**A hydrocarbon stem that is typically followed by an unsaturator e.g. "adamant" */
	static final String CYCLICUNSATURABLEHYDROCARBON_SUBTYPE_VAL ="cyclicUnsaturableHydrocarbon";
	/**Replacmenet terms that are not substituents e.g.  amido/hydrazido/imido/nitrido*/
	static final String DEDICATEDFUNCTIONALREPLACEMENTPREFIX_SUBTYPE_VAL = "dedicatedFunctionalReplacementPrefix";
	/**An atom e.g. "lithium" */
	static final String ELEMENTARYATOM_SUBTYPE_VAL ="elementaryAtom";
	/**An amino acid that ends in an e.g. tryptoph */
	static final String ENDINAN_SUBTYPE_VAL ="endInAn";
	/**An amino acid that ends in ic e.g. aspart */
	static final String ENDINIC_SUBTYPE_VAL ="endInIc";
	/**An amino acid that ends in ine e.g. alan */
	static final String ENDININE_SUBTYPE_VAL ="endInIne";
	/**A substituent that is expected to form a bridge e.g. "epoxy", "epiimino" */
	static final String EPOXYLIKE_SUBTYPE_VAL ="epoxyLike";
	/**A ring that will be fused onto another ring e.g. "benzo", "pyrido", "pyridino" */
	static final String FUSIONRING_SUBTYPE_VAL ="fusionRing";
	/**A group that can be suffixed e.g. "hydrazin" */
	static final String GROUPSTEM_SUBTYPE_VAL ="groupStem";
	/**A halide or pseudo halide e.g. "bromo", "cyano". Can be functional replacment terms when preceding certain non-carboxylic acids */
	static final String HALIDEORPSEUDOHALIDE_SUBTYPE_VAL = "halideOrPseudoHalide";
	/**The stem of a hantzch Widman ring sytem e.g. "an", "ol", "olidin"  */
	static final String HANTZSCHWIDMAN_SUBTYPE_VAL ="hantzschWidman";
	/**A heteroatom hydride e.g. "az" "sulf" (will be followed by an unsaturator, may be preceded by a multiplier to form the heteroatom equivalent of alkanes)*/
	static final String HETEROSTEM_SUBTYPE_VAL ="heteroStem";
	/**A group with no special properties Similar to: {@link XmlDeclarations#NONE_SUBTYPE_VAL}*/
	static final String SIMPLEGROUP_SUBTYPE_VAL ="simpleGroup";
	/**A substituent which intrinsically forms multiple bonds e.g. "siloxane", "thio" */
	static final String MULTIRADICALSUBSTITUENT_SUBTYPE_VAL ="multiRadicalSubstituent";
	/**A non-carboxylic acid which cannot form a substituent e.g. "bor" */
	static final String NOACYL_SUBTYPE_VAL ="noAcyl";
	/**A group with no special properties Similar to: {@link XmlDeclarations#SIMPLEGROUP_SUBTYPE_VAL}*/
	static final String NONE_SUBTYPE_VAL ="none";
	/**oxido/sulfido/selenido/tellurido These are handled similarly to oxide e.g. might give -[O-] or =O*/
	static final String OXIDOLIKE_SUBTYPE_VAL ="oxidoLike";
	/**A term indicating replacement of all substitutable hydrogens by a halogen e.g. "perchloro" */
	static final String PERHALOGENO_SUBTYPE_VAL ="perhalogeno";
	/** phospho and other very related substituents. Strongly prefer forming bonds to hydroxy groups */
	static final String PHOSPHO_SUBTYPE_VAL ="phospho";
	/** A component of a salt e.g "hydrate", "2HCl" */
	static final String SALTCOMPONENT_SUBTYPE_VAL ="saltComponent";
	/**A substitutent with no suffix e.g. "amino" */
	static final String SIMPLESUBSTITUENT_SUBTYPE_VAL ="simpleSubstituent";
	/**A substituent expecting a suffix e.g."bor" "vin" */
	static final String SUBSTITUENT_SUBTYPE_VAL ="substituent";
	/**A group representing a straight chain carbohydrate of a certain length with undefined stereochemistry e.g. hex in hexose */
	static final String SYSTEMATICCARBOHYDRATESTEMALDOSE_SUBTYPE_VAL ="systematicCarbohydrateStemAldose";
	/**A group representing a straight chain carbohydrate of a certain length with undefined stereochemistry e.g. hex in hex-2-ulose */
	static final String SYSTEMATICCARBOHYDRATESTEMKETOSE_SUBTYPE_VAL ="systematicCarbohydrateStemKetose";
	/**A suffix that attaches to the end of a chain e.g. "aldehyde", "ic acid" */
	static final String TERMINAL_SUBTYPE_VAL ="terminal";
	/**An acid that when suffixed with yl gives an acyl group e.g. "acet" */
	static final String YLFORACYL_SUBTYPE_VAL ="ylForAcyl";
	/**An acid that has undefined meaning when suffixed with yl */
	static final String YLFORNOTHING_SUBTYPE_VAL ="ylForNothing";
	/**An acid that when suffixed with yl gives an alkyl group e.g. "laur" */
	static final String YLFORYL_SUBTYPE_VAL ="ylForYl";
	
	/**Requests that no labelling should be applied */
	static final String NONE_LABELS_VAL ="none";//TODO no labels attribute should probably mean no labelling

	/**Requests that labelling be done like a fused ring. It is assumed that the order of the atoms is locant 1 as the first atom*/
	static final String FUSEDRING_LABELS_VAL ="fusedRing";
	
	/**Requests that labelling be 1, 2, 3 etc. It is assumed that the order of the atoms is locant 1 as the first atom*/
	static final String NUMERIC_LABELS_VAL ="numeric";

	/** InLocants have not been specified */
	static final String INLOCANTS_DEFAULT = "default";

	/**
	 * See suffixRules.dtd
	 */
	static final String SUFFIXRULES_RULE_EL = "rule";
	static final String SUFFIXRULES_VALUE_ATR = "value";
	static final String SUFFIXRULES_SMILES_ATR = "SMILES";
	static final String SUFFIXRULES_LABELS_ATR = "labels";
	static final String SUFFIXRULES_FUNCTIONALIDS_ATR = "functionalIDs";
	static final String SUFFIXRULES_OUTIDS_ATR = "outIDs";
	static final String SUFFIXRULES_KETONELOCANT_ATR = "ketoneLocant";
	static final String SUFFIXRULES_ORDER_ATR = "order";
	static final String SUFFIXRULES_OUTVALENCY_ATR = "outValency";
	static final String SUFFIXRULES_CHARGE_ATR = "charge";
	static final String SUFFIXRULES_PROTONS_ATR = "protons";
	static final String SUFFIXRULES_ELEMENT_ATR = "element";
	
	/**
	 * See suffixApplicability.dtd
	 */
	static final String SUFFIXAPPLICABILITY_GROUPTYPE_EL = "groupType";
	static final String SUFFIXAPPLICABILITY_SUFFIX_EL = "suffix";
	static final String SUFFIXAPPLICABILITY_TYPE_ATR = "type";
	static final String SUFFIXAPPLICABILITY_VALUE_ATR = "value";
	static final String SUFFIXAPPLICABILITY_SUBTYPE_ATR = "subType";
}
