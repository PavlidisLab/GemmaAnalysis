package chibi.gemmaanalysis;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;

import cern.colt.list.DoubleArrayList;

//import ubic.gemma.model.analysis.expression.coexpression.GeneCoexpressionAnalysis;
//import ubic.gemma.model.association.coexpression.Gene2GeneCoexpression;
//import ubic.gemma.model.association.coexpression.Gene2GeneCoexpressionService;
//import ubic.gemma.model.association.coexpression.HumanGeneCoExpression;
//import ubic.gemma.model.association.coexpression.MouseGeneCoExpression;
//import ubic.gemma.model.association.coexpression.OtherGeneCoExpression;
//import ubic.gemma.model.association.coexpression.RatGeneCoExpression;
import ubic.gemma.model.analysis.expression.coexpression.CoexpressionCollectionValueObject;
import ubic.gemma.model.analysis.expression.coexpression.CoexpressionValueObject;
import ubic.gemma.model.association.coexpression.Probe2ProbeCoexpressionService;
import ubic.gemma.model.expression.arrayDesign.ArrayDesign;
import ubic.gemma.model.expression.arrayDesign.ArrayDesignService;
import ubic.gemma.model.expression.bioAssayData.DoubleVectorValueObject;
import ubic.gemma.model.expression.bioAssayData.ProcessedExpressionDataVectorService;
import ubic.gemma.model.expression.bioAssayData.ProcessedExpressionDataVectorDao.RankMethod;
import ubic.gemma.model.expression.designElement.CompositeSequence;
import ubic.gemma.model.expression.designElement.CompositeSequenceService;
import ubic.gemma.model.expression.designElement.DesignElement;
import ubic.gemma.model.expression.experiment.BioAssaySet;
import ubic.gemma.model.expression.experiment.ExpressionExperiment;
import ubic.gemma.model.expression.experiment.ExpressionExperimentService;
import ubic.gemma.model.genome.Gene;
import ubic.gemma.model.genome.PhysicalLocation;
//import ubic.gemma.model.genome.PhysicalLocationService;
import ubic.gemma.model.genome.Taxon;
import ubic.gemma.model.genome.TaxonService;
import ubic.gemma.model.genome.biosequence.BioSequence;
import ubic.gemma.model.genome.biosequence.SequenceType;
import ubic.gemma.model.genome.gene.GeneService;
import ubic.gemma.model.genome.sequenceAnalysis.BlatAssociation;
import ubic.gemma.model.genome.sequenceAnalysis.BlatAssociationService;
import ubic.gemma.util.AbstractSpringAwareCLI;

import ubic.gemma.analysis.expression.coexpression.ProbeLinkCoexpressionAnalyzer;
//import ubic.gemma.datastructure.matrix.ExpressionDataDoubleMatrix;
//import ubic.gemma.datastructure.matrix.ExpressionDataMatrixRowElement;

public class PARMapperAnalyzeCLI extends AbstractSpringAwareCLI {

	private Taxon taxon;

	// Input files and output directories
	// These can and should be overrided by arguments -i, -p, and -o
	private String inFile = "/home/hmokada/scratch/human-par-gene-relations.subset.txt";
	private String pargeneexpFile = null;// = "/home/hmokada/scratch/allpargeneexps_dontuse.txt";
	private String ExperimentListFile = null;// = "/home/hmokada/scratch/expressionexperimentsDoubleStranded.txt";
	private String outFileDir = "/home/hmokada/scratch/outputfiles";
	
	// All output files will be appended to the stem, it is overridden by argument -s
	private String stem = "human-par-gene-relations";

	// Options that determine analysis is to be performed
	private boolean checkPAR = false;					// 
	private boolean checkGeneRank = false;				// 
	private boolean checkPARprobe = false;				// 
	private boolean checkGeneCorank = false;			// 
	private boolean checkAllCoexp = false;				// 
	private boolean checkAllCoexpInDB = false;			// 
	private boolean checkAllPairsCoexp = false;			// 
	
	private boolean checkExperimentTypes = false;		// 
	private boolean checkUniqueProbeMappings = false;	// 


	private static Map headerLookup;
//	private static String[] headers;

	private static Collection<String[]> records;

	private static Map<Long, long[]> parToCoExps;
	private static Map parFileEntries;

	// Services used
	GeneService parService;
	TaxonService taxonService;
	// ExpressionExperimentSetService expressionExperimentSetService;
	ExpressionExperimentService expressionExperimentService;
	ProcessedExpressionDataVectorService processedExpressionDataVectorService;
	CompositeSequenceService compositeSequenceService;
	ArrayDesignService arrayDesignService;
	BlatAssociationService blatAssociationService;
	// Gene2GeneCoexpressionService gene2GeneCoexpressionService;
	Probe2ProbeCoexpressionService probe2ProbeCoexpressionService;

	@SuppressWarnings("static-access")
	@Override
	protected void buildOptions() {
		// TODO Auto-generated method stub

		Option taxonOption = OptionBuilder.hasArg().isRequired()
				.withDescription("taxon name").withDescription("taxon to use")
				.withLongOpt("taxon").create('t');
		addOption(taxonOption);

		Option inFileOption = OptionBuilder.hasArg().withDescription(
				"infile name").withDescription("file to read").withLongOpt(
				"inFile").create('i');
		addOption(inFileOption);

		Option outFileDirOption = OptionBuilder.hasArg().withDescription(
				"outdirectory name").withDescription("directory to write to")
				.withLongOpt("outFileDir").create('o');
		addOption(outFileDirOption);

		Option pargeneexpFileOption = OptionBuilder.hasArg().withDescription(
				"PAR-gene-experiment map filename")
				.withDescription("file of Par/Gene/Exps to read").withLongOpt(
				"pargeneexpFile").create('p');
		addOption(pargeneexpFileOption);
		
		Option ExperimentListFileOption = OptionBuilder.hasArg().withDescription(
				"experiment list filename")
				.withDescription("file of expression experiments to read").withLongOpt(
				"pargeneexpFile").create('e');
		addOption(ExperimentListFileOption);
		
		Option outputStemOption = OptionBuilder.hasArg().withDescription(
				"output file stem name")
				.withDescription("output file stem name").withLongOpt(
				"stem").create('s');
		addOption(outputStemOption);

		Option checkParOption = OptionBuilder.create("checkPAR");
		addOption(checkParOption);

		Option checkPARprobeOption = OptionBuilder.create("checkPARprobe");
		addOption(checkPARprobeOption);

		Option checkGeneCorankOption = OptionBuilder.create("checkGeneCorank");
		addOption(checkGeneCorankOption);

		Option checkGeneRankOption = OptionBuilder.create("checkGeneRank");
		addOption(checkGeneRankOption);

		Option checkExperimentTypesOption = OptionBuilder
				.create("checkExperimentTypes");
		addOption(checkExperimentTypesOption);

		Option checkUniqueProbeMappingsOption = OptionBuilder
				.create("checkUniqueProbeMappings");
		addOption(checkUniqueProbeMappingsOption);
		
		Option checkAllCoexpOption = OptionBuilder
		.create("checkAllCoexp");
		addOption(checkAllCoexpOption);
		
		Option checkAllCoexpInDBOption = OptionBuilder
		.create("checkAllCoexpInDB");
		addOption(checkAllCoexpInDBOption);

		Option checkAllPairsCoexpOption = OptionBuilder
		.create("checkAllPairsCoexp");
		addOption(checkAllPairsCoexpOption);

		
	}

	@Override
	protected void processOptions() {
		super.processOptions();

		this.taxonService = (TaxonService) this.getBean("taxonService");
		this.parService = (GeneService) this.getBean("geneService");

		/*
		 * Add processing o
		 */

		if (hasOption("i")) {
			inFile = getOptionValue("i");

			// if ( == null) {
			// log.error("ERROR: Cannot find file " + inFile);
			// }
		}

		if (hasOption("o")) {
			outFileDir = getOptionValue("o");

			// if ( == null) {
			// log.error("ERROR: Cannot find directory " + outFileDir);
			// }
		}

		if (hasOption("p")) {
			pargeneexpFile = getOptionValue("p");
		}
		
		if (hasOption("e")) {
			ExperimentListFile = getOptionValue("e");
		}
		
		if (hasOption("s")) {
			stem = getOptionValue("s");
		}

		if (hasOption("t")) {
			taxon = taxonService.findByCommonName(getOptionValue("t"));

			if (taxon == null) {
				log.error("ERROR: Cannot find taxon " + getOptionValue("t"));
			}
		}

		this.checkPAR = this.hasOption("checkPAR");
		this.checkPARprobe = this.hasOption("checkPARprobe");
		this.checkGeneCorank = this.hasOption("checkGeneCorank");
		this.checkGeneRank = this.hasOption("checkGeneRank");

		this.checkExperimentTypes = this.hasOption("checkExperimentTypes");

		this.checkUniqueProbeMappings = this
				.hasOption("checkUniqueProbeMappings");
		this.checkAllCoexp = this.hasOption("checkAllCoexp");
		this.checkAllCoexpInDB = this.hasOption("checkAllCoexpInDB");
		this.checkAllPairsCoexp = this.hasOption("checkAllPairsCoexp");
		
	}

	@Override
	protected Exception doWork(String[] args) {
		// TODO Auto-generated method stub

		Exception exc = processCommandLine("test", args);

		this.parService = (GeneService) this.getBean("geneService");
		this.taxonService = (TaxonService) this.getBean("taxonService");
		// this.expressionExperimentSetService = (
		// ExpressionExperimentSetService ) this.getBean(
		// "expressionExperimentSetService" );
		this.expressionExperimentService = (ExpressionExperimentService) this
				.getBean("expressionExperimentService");
		this.processedExpressionDataVectorService = (ProcessedExpressionDataVectorService) this
				.getBean("processedExpressionDataVectorService");
		this.compositeSequenceService = (CompositeSequenceService) this
				.getBean("compositeSequenceService");
		this.arrayDesignService = (ArrayDesignService) this
				.getBean("arrayDesignService");
		this.blatAssociationService = (BlatAssociationService) this
				.getBean("blatAssociationService");
		// this.gene2GeneCoexpressionService = ( Gene2GeneCoexpressionService )
		// this.getBean( "gene2GeneCoexpressionService" );
		this.probe2ProbeCoexpressionService = (Probe2ProbeCoexpressionService) this
				.getBean("probe2ProbeCoexpressionService");

		// CompositeSequenceService.getGenes or getGenesWithSpecificity.
		/*
		 * //CompositeSequence cz = compositeSequenceService.load(new
		 * Long(2365242)); Collection<CompositeSequence> cccz = new ArrayList<CompositeSequence>();
		 * //cccz.add(cz); cccz.add(compositeSequenceService.load(new
		 * Long(2365242))); cccz.add(compositeSequenceService.load(new
		 * Long(406118))); Map<CompositeSequence, Map<PhysicalLocation,
		 * Collection<BlatAssociation>>> ddd =
		 * compositeSequenceService.getGenesWithSpecificity(cccz);
		 * //compositeSequenceService.get
		 * 
		 * for (CompositeSequence cs : ddd.keySet()) { System.out.println(cs);
		 * for (PhysicalLocation pl : ddd.get(cs).keySet()) {
		 * System.out.println("\t"+pl); for (BlatAssociation ba :
		 * ddd.get(cs).get(pl)) { blatAssociationService.thaw(ba);
		 * System.out.println("\t\t"+ba); } } }
		 */
		
		
		
		

		// Load up experiments - all experiments or from a list
		Collection<ExpressionExperiment> eeCol;
		if (ExperimentListFile == null) {
			eeCol = expressionExperimentService.findByTaxon(taxon);
		} else {
			eeCol = loadExpressionExperimentsByFile(ExperimentListFile);
		}
		System.out.println("Experiments loaded:" + eeCol.size());
		
		
		// prints out experiment information (exits program when finished)
		if (checkExperimentTypes) {
			printExperimentTypes(eeCol);
			return null;
		}
		
		
		// log.info(taxonService.loadAll());
		// For my mac
		// String inFile =
		// "/Users/mokada/development/PARs/data/human-par-gene-relations.txt";
		// String outFileDir =
		// "/Users/mokada/development/PARs/data/outputfiles";
		int batchSize = 50;

		PrintStream pxx = null; // declare a print stream object
		PrintStream pxe = null;
		PrintStream pex = null;
		PrintStream pee = null;

		PrintStream ppr = null;

		PrintStream pgp_xx = null; // for gene vs. PAR rankings (within same
									// experiments)
		PrintStream pgp_xe = null;
		PrintStream pgp_ex = null;
		PrintStream pgp_ee = null;

		PrintStream pco = null;	//correlation between par/gene pairs
		PrintStream pcA = null; //correlation between par and all others in Link analysis stored in DB
		PrintStream pap = null; //correlation between all par/gene pairs - calculates all

		try {
			// Connect print stream to the output stream
			pxx = new PrintStream(new FileOutputStream(outFileDir
					+ "/" + stem + ".output.xx.txt"));
			pxe = new PrintStream(new FileOutputStream(outFileDir
					+ "/" + stem + ".output.xe.txt"));
			pex = new PrintStream(new FileOutputStream(outFileDir
					+ "/" + stem + ".output.ex.txt"));
			pee = new PrintStream(new FileOutputStream(outFileDir
					+ "/" + stem + ".output.ee.txt"));

			ppr = new PrintStream(new FileOutputStream(outFileDir
					+ "/" + stem + ".output.probelevel.txt"));

			pgp_xx = new PrintStream(new FileOutputStream(outFileDir
					+ "/" + stem + ".output.geneVsPar.xx.txt"));
			pgp_xe = new PrintStream(new FileOutputStream(outFileDir
					+ "/" + stem + ".output.geneVsPar.xe.txt"));
			pgp_ex = new PrintStream(new FileOutputStream(outFileDir
					+ "/" + stem + ".output.geneVsPar.ex.txt"));
			pgp_ee = new PrintStream(new FileOutputStream(outFileDir
					+ "/" + stem + ".output.geneVsPar.ee.txt"));

			pco = new PrintStream(new FileOutputStream(outFileDir
					+ "/" + stem + ".output.genecoexpression.txt"));
			pcA = new PrintStream(new FileOutputStream(outFileDir
					+ "/" + stem + ".output.genecoexpressionAll.txt"));
			pap = new PrintStream(new FileOutputStream(outFileDir
					+ "/" + stem + ".output.genecoexpressionAllPairs.txt"));

		} catch (Exception e) {
			System.err.println("Error writing to file");
			System.exit(0);
		}

		// allow a user to enter a filename through the command line
		// if (0 < args.length && args[0] != null) {
		// inFile = args[0];
		// } *** replaced with java options

		System.out.println("Reading file: " + inFile);
		readPARFile(inFile);
		if (pargeneexpFile != null) readpargeneFile(pargeneexpFile);

		// Establish the column numbers

		int ParIDIdx = getIndex("ParID");
		int ParNameIdx = getIndex("ParName");
		int ChromIdx = getIndex("Chrom");
		int NucIdx = getIndex("Nuc");
		int GeneIdIdx = getIndex("GeneId");
		int GeneSymbolIdx = getIndex("GeneSymbol");
		int DistanceIdx = getIndex("Distance");
		int GeneContainsParIdx = getIndex("GeneContainsPar");
		int SameStrandIdx = getIndex("SameStrand");

		// Output column ordering
		String outputHeader = "ParID,GeneId,ParName,Chrom,Nuc,GeneSymbol,Distance,GeneContainsPar,SameStrand,NumExperiments,NumSamples,Rank";
		pxx.println(outputHeader);
		pxe.println(outputHeader);
		pex.println(outputHeader);
		pee.println(outputHeader);

		// for probe level rannkings
		ppr.println("ParID,GeneId,ParName,Chrom,Nuc,GeneSymbol,Distance,GeneContainsPar,SameStrand,CompSeqId,NumExperiments,RankEMean_MethMean,RankEMean_MethMax,RankEMax_MethMean,RankEMax_MethMax");

		// for Gene vs Par rannkings
		pgp_xx.println("ParID,GeneId,ParName,Chrom,Nuc,GeneSymbol,Distance,GeneContainsPar,SameStrand,ParNumExperiments,ParNumSamples,ParRank,GeneNumExperiments,GeneNumSamples,GeneRank");
		pgp_xe.println("ParID,GeneId,ParName,Chrom,Nuc,GeneSymbol,Distance,GeneContainsPar,SameStrand,ParNumExperiments,ParNumSamples,ParRank,GeneNumExperiments,GeneNumSamples,GeneRank");
		pgp_ex.println("ParID,GeneId,ParName,Chrom,Nuc,GeneSymbol,Distance,GeneContainsPar,SameStrand,ParNumExperiments,ParNumSamples,ParRank,GeneNumExperiments,GeneNumSamples,GeneRank");
		pgp_ee.println("ParID,GeneId,ParName,Chrom,Nuc,GeneSymbol,Distance,GeneContainsPar,SameStrand,ParNumExperiments,ParNumSamples,ParRank,GeneNumExperiments,GeneNumSamples,GeneRank");
		
		pco.println("ParID,GeneId,ParName,Chrom,Nuc,GeneSymbol,Distance,GeneContainsPar,SameStrand,PosCorr,NegCorr");
		pcA.println("ParID,GeneId,ParName,Chrom,Nuc,GeneSymbol,Distance,GeneContainsPar,SameStrand,CoexpGeneId,PosCorr,NegCorr");
		pap.println("ParID,GeneId,ParName,Chrom,Nuc,GeneSymbol,Distance,GeneContainsPar,SameStrand,ExperimentId,NumParMapProbe,NumGenesMapProbe,VectorSize,ParProbe,GeneProbe,Pearson");

		// System.out.println("\n\nExperiments by taxon\n\n");
		// taxon = taxonService.findByCommonName( "human" );
		// Add to options

		

		// Get the distribution of ranks for ALL human genes!!!!
		// ///////////////////////////////////////////////// Start
		// mess!!!!!!!!!!
		/*
		 * ppr.println("GeneId,GeneName,DesignElmId,numOfExperiments,valExpMean_RankMean,valExpMean_RankMax,valExpMax_RankMean,valExpMax_RankMax");
		 * 
		 * 
		 * 
		 * Collection<Gene> allGenes = this.parService.loadKnownGenes(taxon);
		 * System.out.println("All known genes: "+ allGenes.size());
		 * 
		 * ppr.println("#All genes: "+ allGenes.size()); for (Gene g : allGenes) {
		 * ppr.println("#"+ g.getId() + "\t" + g.getName()); }
		 * 
		 * //if (true) return null; Iterator geneItr = allGenes.iterator(); int
		 * batchg = 0;
		 * 
		 * //for (Gene g : allGenes) { while (geneItr.hasNext() //&& batchg < 2 ) {
		 * Collection<Gene> pars = new ArrayList<Gene>(); int count =
		 * batchSize;
		 *  // work with a small batch while (geneItr.hasNext() && 0 < count) {
		 * 
		 * Gene g = (Gene) geneItr.next();
		 *  // TODO: place as log if (g == null) { System.out.println("Gene
		 * doesn't exist: "+ g.getId()); continue; }
		 * 
		 * pars.add(g);
		 * 
		 * count--; // count down the number of genes in a batch
		 *  }
		 * 
		 * 
		 * ppr.println("#Batch Number: "+ batchg); batchg++; ppr.flush();
		 * 
		 * 
		 * outputAll_probelevel(eeCol, pars, ppr);
		 *  }
		 * 
		 * 
		 * 
		 * 
		 * ppr.close();
		 * 
		 * 
		 * if (true) return null;
		 */
		// ///////////////////////////////////////////////// End mess!!!!!!!!!!

		Iterator recordItr = records.iterator();
		int batch = 0;

		String[] record;
		while (recordItr.hasNext()) {

			parFileEntries = new HashMap<Long, String>(batchSize);

			HashMap<Gene, Gene> parToGene = new HashMap<Gene, Gene>();

			Collection<Gene> pars = new ArrayList<Gene>();
			Collection<Gene> genes = new ArrayList<Gene>();
			int count = batchSize;

			// work with a small batch
			while (recordItr.hasNext() && 0 < count) {

				record = (String[]) recordItr.next();
				// System.out.println(record);

				// accessing data elements
				int ParID = Integer.parseInt(record[ParIDIdx]);
				int GeneId = Integer.parseInt(record[GeneIdIdx]);

				// it's not necessary to parse thees to int/boolean as it'll
				// printed off to file anyway
				String ParName = record[ParNameIdx];
				String Chrom = record[ChromIdx];
				// int Nuc = Integer.parseInt( record[NucIdx] );
				String Nuc = record[NucIdx];
				String GeneSymbol = record[GeneSymbolIdx];
				// int Distance = Integer.parseInt( record[DistanceIdx] );
				// boolean GeneContainsPar = Boolean.parseBoolean(
				// record[GeneContainsParIdx] );
				// boolean SameStrand = Boolean.parseBoolean(
				// record[SameStrandIdx] );
				String Distance = record[DistanceIdx];
				String GeneContainsPar = record[GeneContainsParIdx];
				String SameStrand = record[SameStrandIdx];

				String allEntries = ParID + "," + GeneId + "," + ParName + ","
						+ Chrom + "," + Nuc + "," + GeneSymbol + "," + Distance
						+ "," + GeneContainsPar + "," + SameStrand;

				// HashMap<Integer, String> parFileEntries = new
				// HashMap<Integer, String>(batchSize);

				parFileEntries.put(Long.valueOf(ParID), allEntries); // .put(ParID,
																		// record);

				Gene par = this.parService.load(Long.valueOf(ParID)); // this.parService.load(ParID);
				Gene g = this.parService.load(GeneId);

				// TODO: place as log
				if (par == null || g == null) {
					System.out.println("PAR or Gene doesn't exist: " + par
							+ "\t" + g);
					continue;
				}

				pars.add(par);
				genes.add(g);

				parToGene.put(par, g);

				count--; // count down the number of genes in a batch

			}

			// print results to file

			
			pxx.println("Batch Number: "+ batch);
			pxe.println("Batch Number: "+ batch);
			pex.println("Batch Number: "+ batch);
			pee.println("Batch Number: "+ batch);
			
			ppr.println("Batch Number: "+ batch);
			
			pgp_xx.println("Batch Number: "+ batch);
			pgp_xe.println("Batch Number: "+ batch);
			pgp_ex.println("Batch Number: "+ batch);
			pgp_ee.println("Batch Number: "+ batch);
			
			pco.println("Batch Number: " + batch);
			pcA.println("Batch Number: " + batch);
			pap.println("Batch Number: " + batch);
			
			// ParID,GeneId,ParName,Chrom,Nuc,GeneSymbol,Distance,GeneContainsPar,SameStrand,ParNumExperiments,ParNumSamples,ParRank,GeneNumExperiments,GeneNumSamples,GeneRank

			batch++;

			pxx.flush();
			pxe.flush();
			pex.flush();
			pee.flush();

			ppr.flush();

			pgp_xx.flush();
			pgp_xe.flush();
			pgp_ex.flush();
			pgp_ee.flush();

			pco.flush();
			pcA.flush();
			pap.flush();

			if (checkPARprobe) {
				outputRankingsProbelevel(eeCol, pars, ppr);
			}

			if (checkPAR) {
				outputRankings(eeCol, pars, pxx, pxe, "max");
				outputRankings(eeCol, pars, pex, pee, "mean");
			}
			
			if (checkGeneRank) {
				outputRankings(eeCol, genes, pxx, pxe, "max");
				outputRankings(eeCol, genes, pex, pee, "mean");
			}

			if (checkGeneCorank) {
				outputPARGeneCorank(eeCol, pars, genes, pgp_xx, pgp_xe, "max");
				outputPARGeneCorank(eeCol, pars, genes, pgp_ex, pgp_ee, "mean");

				// not working
				// outputAll_geneCorank(eeCol, pars, genes, pgp_xx);

				// outputAll_geneCorank(eeCol, parToGene, pxx, pxe, "max");
				// outputAll_geneCorank(eeCol, parToGene, pex, pee, "mean");
			}

			
			if (checkAllCoexpInDB) {
				//if (parToCoExps != null) {
				if (pargeneexpFile == null) {
					System.out.println("Error, need to supply pargeneexp file");
					System.exit(0);
				}
				//check all coexpressions given a par or look at each par-gene combination
				outputPARGeneCoexpressionLinks(pars, genes, pco);
				
			}
			
			if (checkAllCoexp) {
				outputAllGeneCoexpressionLinks(pars, eeCol, 1, pcA);
			}
			
			if (checkAllPairsCoexp) {
				if (pargeneexpFile == null) {
					System.out.println("Error, need to supply pargeneexp file");
					System.exit(0);
				}
					
				outputAllPARGeneCoexpressionCalculated(pars, genes, pap);
			}

		}
		pxx.close();
		pxe.close();
		pex.close();
		pee.close();

		ppr.close();
		
		
		pgp_xx.close();
		
		pgp_xe.close();
		pgp_ex.close();
		pgp_ee.close();

		pco.close();
		pcA.close();
		pap.close();

		

		return null;
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		PARMapperAnalyzeCLI p = new PARMapperAnalyzeCLI();
		Exception e = p.doWork(args);
		if (e != null) {
			throw new RuntimeException(e);
		}
	}

	/*
	 * This method finds the array types for the expressions given It might be
	 * important to know which array designs make use of double stranded probes
	 */
	private void printExperimentTypes(Collection<ExpressionExperiment> eeCol) {

		// Printing out array information
		for (ExpressionExperiment e : eeCol) {
			System.out.println("#" + e.getId() + ": " + e.getName());
			Collection<ArrayDesign> ad = expressionExperimentService
					.getArrayDesignsUsed(e);
			for (ArrayDesign a : ad) {
				System.out.println("#\t" + a.getId() + ": " + a.getName());
				// Collection<CompositeSequence> csc =
				// a.getCompositeSequences();
				// Iterator csItr = csc.iterator();
				// CompositeSequence cs = (CompositeSequence) csItr.next();
				// CompositeSequence cs = (CompositeSequence)
				// a.getCompositeSequences().iterator().next();

				// If everything was perfect and easy
				// SequenceType st = ((CompositeSequence)
				// a.getCompositeSequences().iterator().next()).getBiologicalCharacteristic().getType();

				SequenceType st = getSequenceType(a);

				System.out.println("#\t\t" + st);

				System.out.println(e.getId() + "\t" + a.getId() + "\t" + st);
				// System.out.println("\t"+st.toString());
				// System.out.println("\t"+st.getValue());

				// System.out.println("\t"+a.getTechnologyType());
				// for (String name : a.getTechnologyType().names()) {
				// System.out.println("\t\t."+name);
				// }
			}

		}
	}

	private SequenceType getSequenceType(ArrayDesign a) {
		arrayDesignService.thawLite(a);

		Collection<CompositeSequence> csc = a.getCompositeSequences();
		Iterator cscItr = csc.iterator();
		if (cscItr == null) {
			System.out.println("***No iterator: " + csc);
			return null;
		}

		if (!cscItr.hasNext()) {
			System.out.println("***No iterator next obj: " + csc);
			return null;
		}

		while (cscItr.hasNext()) {
			CompositeSequence cs = (CompositeSequence) cscItr.next();

			if (cs == null) {
				// System.out.println("***No compseq: "+cscItr);
				continue;
			}
			BioSequence bs = cs.getBiologicalCharacteristic();
			if (bs == null) {
				// System.out.println("***No bio sequence: "+cs);
				continue;
			}
			SequenceType st = bs.getType();
			if (st == null) {
				// System.out.println("***No sequencetype: "+bs);
				continue;
			}

			return st;

		}
		System.out.println("***No probes have sequence type: " + csc);
		return null;

	}


	// Retrieves all coexpression data for pars and genes that are in the same
	// experiment
	private void outputPARGeneCoexpressionLinks(Collection<Gene> pars,
		Collection<Gene> genes, PrintStream pco) {
		
		ProbeLinkCoexpressionAnalyzer pca = new ProbeLinkCoexpressionAnalyzer();
		pca.setGeneService(parService);
		pca.setProbe2ProbeCoexpressionService(probe2ProbeCoexpressionService);

		
		Iterator<Gene> pItr = pars.iterator();
		Iterator<Gene> gItr = genes.iterator();

		// for (Gene par: pars) {
		while (pItr.hasNext()) {
			Gene par = pItr.next();
			Gene gene = gItr.next();

			//bad
			par  = parService.load(1445595);
			gene = parService.load(127958);

			//good
//			par  = parService.load(1445513);
//			gene = parService.load(442313);
			
			// if (par)
			long[] expIds = parToCoExps.get(par.getId());

			if (expIds == null) {
				//System.out.println("No information: " + par.getId());
				continue;
			}
			
			// change expression exp objects to bioassaysets
			Collection<BioAssaySet> eexps = new ArrayList<BioAssaySet>();
			for (int i = 0; i < expIds.length; i++) {
				eexps.add((BioAssaySet) expressionExperimentService.load(expIds[i]));
			}

			Collection<Gene> pargene = new ArrayList<Gene>();
			pargene.add(par);
			pargene.add(gene);

			if (pargene == null) {
				System.out.println("No pargene");
				continue;
			}
			if (eexps == null) {
				System.out.println("No exps");
				continue;
			}


			System.out.println("Linkanalysis " + par.getId() + "\t" + gene.getId());
			
			Map<Gene, CoexpressionCollectionValueObject> coexp;
			
//			try {
				//Map<Gene, CoexpressionCollectionValueObject> 
				coexp = pca.linkAnalysis(pargene, eexps, 1, false, true, 0);
//			} catch (java.lang.IndexOutOfBoundsException e) {
//				e.printStackTrace();
//				System.out.println("Error performing linkAnalysis "
//						+ e.getMessage());
//				continue;
//			}
			

			// System.out.print("Par "+ par.getId()+"\t");

			for (Gene g : coexp.keySet()) {
				CoexpressionCollectionValueObject ccvo = coexp.get(g);
				
				Collection<CoexpressionValueObject> cvos = ccvo.getAllGeneCoexpressionData(0);
				
				if (cvos.size() >  0) {
					// found evidence of coexpression
				}
				
				System.out.println("GeneID: " + g.getId()+" Size: " + cvos.size());

				
				for (CoexpressionValueObject cvo : cvos) {
					//System.out.println("\t\tTS: " + cvo);
					
					System.out.println(parFileEntries.get(par.getId())
							+ "," + g.getId()
							+ "," + cvo.getPositiveScore()
							+ "," + cvo.getNegativeScore());
					
					
					pco.println(parFileEntries.get(par.getId())
							+ "," + cvo.getPositiveScore()
							+ "," + cvo.getNegativeScore());
					
					
				}
				// System.exit(0);

			}

		}
	}

	
	
	
	
	
	
	
	
	
	// Retrieves all coexpression data for a par and all genes that are in the same
	// experiment above a coexpression threshold
	//
	// This method is different from the previous method outgenecoexpressions examining par-gene pairs.
	// Here a given par is compared with all genes containing a link in the database (implying that they
	// already have some coexpression)
	private void outputAllGeneCoexpressionLinks(Collection<Gene> pars, Collection<ExpressionExperiment> ccCol, int stringency, PrintStream pco) {
		

		
		Iterator<Gene> pItr = pars.iterator();

		// for (Gene par: pars) {
		while (pItr.hasNext()) {
			Gene par = pItr.next();

			// if (par)
			//long[] expIds = parToCoExps.get(par.getId());
//
//			if (expIds == null) {
//				//System.out.println("No information: " + par.getId());
//				continue;
//			}

//			Collection<ExpressionExperiment> exps = new ArrayList<ExpressionExperiment>();
			Collection<BioAssaySet> exps2 = new ArrayList<BioAssaySet>();
			
			for (ExpressionExperiment ee : ccCol) {
				exps2.add(ee);
			}

//			for (int i = 0; i < eez.length; i++) {
				// gene2GeneCoexpressionService.
				// GeneCoexpressionAnalysis gca;
				// Gene2GeneCoexpression g = getNewGGCOInstance();
				// g.setFirstGene(par);
				// g.setSecondGene(gene);
				// Double d = g.getPvalue();
				// if (d==null) d = new Double(-5);
//				exps.add(expressionExperimentService.load(expIds[i]));
//				exps2.add((BioAssaySet) expressionExperimentService
//						.load(expIds[i]));
//			}
			// System.out.println();

			Collection<Gene> pargene = new ArrayList<Gene>();
			pargene.add(par);

			if (pargene == null) {
				System.out.println("No pargene");
				continue;
			}
			if (ccCol == null) {
				System.out.println("No exps");
				continue;
			}


//			System.out.println("Linkanalysis " + par.getId());
			
			Map<Gene, CoexpressionCollectionValueObject> coexp;
			
			try {
				ProbeLinkCoexpressionAnalyzer pca = new ProbeLinkCoexpressionAnalyzer();
				pca.setGeneService(parService);
				pca.setProbe2ProbeCoexpressionService(probe2ProbeCoexpressionService);
				//Map<Gene, CoexpressionCollectionValueObject> 
				coexp = pca.linkAnalysis(pargene, exps2, stringency, false, false, 0);
			} catch (java.lang.IndexOutOfBoundsException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				System.out.println("Error performing linkAnalysis "
						+ e.getMessage());
				continue;
			}
			
			
			if (coexp == null || coexp.keySet().size() < 1) {
				System.out.println("No coexpression :(");
				continue;
			}
			
			
			// System.out.print("Par "+ par.getId()+"\t");
			
			//Collection<Integer> sizeCount = new ArrayList<Integer>();
			//int zeroCount = 0;
			//int maxSize = 0;
			//int[] sizes = new int[100];
			
			//for (Gene g : coexp.keySet()) {
			
				//CoexpressionCollectionValueObject ccvo = coexp.get(g);
			CoexpressionCollectionValueObject ccvo = coexp.get(par);
			Collection<CoexpressionValueObject> cvos = ccvo.getAllGeneCoexpressionData(0);
			
			int zeroCount = 0;
			String output_half = (String) parFileEntries.get(par.getId());
			
			int csize = cvos.size();
			
			if (csize ==  0) {
//System.out.println("Zerocount!!");
				zeroCount++;
				continue;
			}
			
			// keep track of how many coexpressions are non-zero
			//sizeCount.add(new Integer(csize));
			//int csizeTmp = csize;
			//if (99 < csizeTmp) csizeTmp = 99;
			//sizes[csizeTmp]++;
			//if (maxSize < csize) maxSize = csize;
			
			
//				System.out.println("GeneID: " + g.getId()+" Size: " + csize);
			
			// for (BioAssaySet ee : exps2) {
			// System.out.println("\tExpid: " + ee.getId());
			// }
			
			for (CoexpressionValueObject cvo : cvos) {
				//System.out.println("\t\tTS: " + cvo);
				String output = output_half
							+ "," + cvo.getGeneId()
							+ "," + csize
							+ "," + cvo.getPositiveScore()
							+ "," + cvo.getNegativeScore();
				System.out.println(output);
				pco.println(output);
				
				
			}
			
			// print out the numbers
			System.out.print("Zerocount: " + par.getId()+"\t"+par.getId()+"\t"+"0:"+zeroCount);
			//for (int i=1; i<99; i++) {
			//	if (sizes[i] == 0) continue;
			//	System.out.print("\t"+i+":"+sizes[i]);
			//}
			//System.out.println();
			
			//}

		}
	}

	
	
	
	
	
	
	
	
	// Looks at all par/gene pairs to see if there are any experiemnts with both
	// This outputs A LOT of information to standard out!!
	// Note: the stdout can be parsed by pargeneExplist.pl which sorts the
	// par/gene/exp combo
	// into lists.
	private void outputPARGeneCorank(Collection<ExpressionExperiment> eeCol,
			Collection<Gene> pars, Collection<Gene> genes, PrintStream px,
			PrintStream pe, String rankMethodStr) {

		RankMethod method;

		// todo: make these boolean
		if (rankMethodStr.equalsIgnoreCase("max")) {
			method = RankMethod.max;
		} else {
			method = RankMethod.mean;
		}

		Collection<Gene> allProbesGenes = new ArrayList<Gene>();
		allProbesGenes.addAll(pars);
		allProbesGenes.addAll(genes);

		// use getRanks function to map genes to experiments
		Map<ExpressionExperiment, Map<Gene, Collection<Double>>> expressRankings = processedExpressionDataVectorService
				.getRanks(eeCol, allProbesGenes, method);
		// Map<ExpressionExperiment, Map<Gene, Collection<Double>>>
		// expressRankings =
		// processedExpressionDataVectorService.getRanks(eeCol, pars, method);

		// need to map experiments to genes to calculate max/mean across
		// experiments
		Map<Gene, Collection<Double>> allGeneRankings = new HashMap<Gene, Collection<Double>>();
		Map<Gene, Collection<Double>> allParRankings = new HashMap<Gene, Collection<Double>>();

		Map<Gene, int[]> geneCounts = new HashMap<Gene, int[]>();

		Map<Gene, Collection<ExpressionExperiment>> numOfExpsPars =  new HashMap<Gene, Collection<ExpressionExperiment>>();
		Map<Gene, Collection<ExpressionExperiment>> numOfExpsGenes = new HashMap<Gene, Collection<ExpressionExperiment>>();

		if (expressRankings.size() < 1) {
			System.out.println("No experiments for this batch:");
			Iterator<Gene> gItr = genes.iterator();
			Iterator<Gene> pItr = pars.iterator();
			Gene g, p;

			while (gItr.hasNext()) {
				g = gItr.next();
				p = pItr.next();
				// System.out.println(+g.getId()+"\t"+p.getId());
			}

			return;
		}

		for (ExpressionExperiment ee : expressRankings.keySet()) {
			System.out.println("Experiment: " + ee.getId());
			Map<Gene, Collection<Double>> expressRankingsExperiment = expressRankings
					.get(ee);

			Iterator<Gene> gItr = genes.iterator();
			Iterator<Gene> pItr = pars.iterator();

			Gene g, p;

			// for ( Gene g : expressRankings.get(ee).keySet() ) {
			while (pItr.hasNext()) {
				g = gItr.next();
				p = pItr.next();
				System.out.print("Checking: " + g.getId() + "\t" + p.getId()
						+ "\t");
				if (!expressRankingsExperiment.containsKey(p))
					System.out.print("NoPar\t");
				if (!expressRankingsExperiment.containsKey(g))
					System.out.print("NoGene\t");
				// System.out.print("Checking:
				// "+g.getId()+"\t"+p.getId()+"\t"+expressRankingsExperiment.get(g).size()+"\t"+expressRankingsExperiment.get(p).size());

				// skip if both gene and PAR not in experiment
				if ((!expressRankingsExperiment.containsKey(g))
						|| (!expressRankingsExperiment.containsKey(p))) {
					System.out.println("skipping, no info");
					continue;
				}

				Collection<Double> candidateGeneRanks = new ArrayList<Double>();
				Collection<Double> candidatePARRanks = new ArrayList<Double>();
				for (Double d : expressRankingsExperiment.get(g)) {
					// subtract this entry from the list of probes, and if
					// necessary, the entire experiment
					if (d == null)
						continue;

					candidateGeneRanks.add(d);
				}

				for (Double d : expressRankingsExperiment.get(p)) {
					// subtract this entry from the list of probes, and if
					// necessary, the entire experiment
					if (d == null)
						continue;

					candidatePARRanks.add(d);
				}

				// skip if both gene and PAR not in experiment
				if (candidateGeneRanks.isEmpty() || candidatePARRanks.isEmpty()) {
					System.out.println("skipping again, at least one is empty");
					continue;
				}

				if (!allGeneRankings.containsKey(p))
					allGeneRankings.put(p, new ArrayList<Double>());
				if (!allParRankings.containsKey(p))
					allParRankings.put(p, new ArrayList<Double>());

				// note that PAR is key for both
				allGeneRankings.get(p).addAll(candidateGeneRanks);
				allParRankings.get(p).addAll(candidatePARRanks);

				// note how many experiments are performed per PAR and gene
				if (!numOfExpsPars.containsKey(g)) {
					numOfExpsPars.put(p, new ArrayList<ExpressionExperiment>());
					numOfExpsGenes
							.put(g, new ArrayList<ExpressionExperiment>());
				}
				numOfExpsPars.get(p).add(ee);
				numOfExpsGenes.get(g).add(ee);

				System.out.println("----Exists!!!----");

			}
		}

		// this time, combine the entries for PARs and genes

		Iterator<Gene> gItr = genes.iterator();
		Iterator<Gene> pItr = pars.iterator();

		Gene g, p;

		// for ( Gene g : expressRankings.get(ee).keySet() ) {
		while (gItr.hasNext()) {
			// for ( Gene p : allGeneRankings.keySet() ) {

			g = gItr.next();
			p = pItr.next();

			// weed out all those that don't have experiments
			if (!allParRankings.containsKey(p)) {
				System.out.println("No Gene information for " + p.getId());
				continue;
			}

			double maxrankPar = 0;
			double meanrankPar = 0;

			double maxrankGene = 0;
			double meanrankGene = 0;

			int numSamplesPar = allParRankings.get(p).size();
			int numSamplesGene = allGeneRankings.get(p).size();

			int numExpsPar = numOfExpsPars.get(p).size();
			int numExpsGene = numOfExpsGenes.get(g).size();

			for (Double d : allParRankings.get(p)) {
				double rank = d.doubleValue();
				if (maxrankPar < rank)
					maxrankPar = rank;
				meanrankPar += rank;
			}
			meanrankPar = meanrankPar / numSamplesPar;

			for (Double d : allGeneRankings.get(p)) {
				double rank = d.doubleValue();
				if (maxrankGene < rank)
					maxrankGene = rank;
				meanrankGene += rank;
			}
			meanrankGene = meanrankGene / numSamplesGene;

			px.println(parFileEntries.get(p.getId()) + "," + numExpsPar + ","
					+ numSamplesPar + "," + maxrankPar + "," + numExpsGene
					+ "," + numSamplesGene + "," + maxrankGene);
			pe.println(parFileEntries.get(p.getId()) + "," + numExpsPar + ","
					+ numSamplesPar + "," + meanrankPar + "," + numExpsGene
					+ "," + numSamplesGene + "," + meanrankGene);

			// px.println(g.getId()+","+numExperiments+","+numSamples+","+maxrank);
			// pe.println(g.getId()+","+numExperiments+","+numSamples+","+meanrank);

		}

	}

	private void outputAll_geneCorank_old(
			Collection<ExpressionExperiment> eeCol,
			HashMap<Gene, Gene> parToGene, PrintStream px, PrintStream pe,
			String rankMethodStr) {

		Collection<Gene> pars = parToGene.keySet();
		Collection<Gene> combined = parToGene.values();
		combined.addAll(pars);

		RankMethod method;

		// todo: make these boolean
		if (rankMethodStr.equalsIgnoreCase("max")) {
			method = RankMethod.max;
		} else {
			method = RankMethod.mean;
		}

		// use getRanks function to map genes to experiments
		Map<ExpressionExperiment, Map<Gene, Collection<Double>>> expressRankings = processedExpressionDataVectorService
				.getRanks(eeCol, pars, method);

		// need to map experiments to genes to calculate max/mean across
		// experiments
		HashMap<Gene, double[]> allGeneRankings = new HashMap<Gene, double[]>();
		HashMap<Gene, int[]> geneCounts = new HashMap<Gene, int[]>();

		for (ExpressionExperiment ee : expressRankings.keySet()) {
			for (Gene p : expressRankings.get(ee).keySet()) {

				double[] maxAndAveRank = new double[2];
				maxAndAveRank[0] = 0; // max
				maxAndAveRank[1] = 0; // mean

				int numNullRanks = 0; // keep track of how many null rankings
										// there are
				boolean allNull = true;

				for (Double d : expressRankings.get(ee).get(p)) {

					// subtract this entry from the list of probes, and if
					// necessary, the entire experiment
					if (d == null) {
						numNullRanks++;
						// System.out.println("Null rank value "+p.getId()+",
						// numNullRanks: "+ numNullRanks);

					} else {

						// calculate the maximum and average rankings
						maxAndAveRank[1] += d;
						if (maxAndAveRank[0] < d)
							maxAndAveRank[0] = d;

						allNull = false;
					}

				}

				// check if genes for this experiment are null-ranked. if so,
				// skip
				if (allNull) {
					// System.out.println("Null: No exps for this gene
					// "+p.getId()+", exp: "+ ee.getId());
					continue;
				}

				maxAndAveRank[1] = maxAndAveRank[1]
						/ (expressRankings.get(ee).get(p).size() - numNullRanks);

				// update hash of expression levels
				if (allGeneRankings.containsKey(p)) {
					if (allGeneRankings.get(p)[1] < maxAndAveRank[0])
						allGeneRankings.get(p)[1] = maxAndAveRank[0];
					allGeneRankings.get(p)[1] += maxAndAveRank[1];

					geneCounts.get(p)[0]++;
					geneCounts.get(p)[1] += expressRankings.get(ee).get(p)
							.size();
				}
				// new gene, must create new entities in hash
				else {
					allGeneRankings.put(p, maxAndAveRank);

					// number of experiments and number samples for this gene
					int[] counts = new int[2];
					counts[0] = 1;
					counts[1] = expressRankings.get(ee).get(p).size();

					geneCounts.put(p, counts);
				}

			}
		}

		for (Gene g : allGeneRankings.keySet()) {

			int numExperiments = geneCounts.get(g)[0];
			int numSamples = geneCounts.get(g)[1];

			double maxrank = allGeneRankings.get(g)[0];
			double meanrank = allGeneRankings.get(g)[1] / numExperiments;

			px.println(parFileEntries.get(g.getId()) + "," + numExperiments
					+ "," + numSamples + "," + maxrank);
			pe.println(parFileEntries.get(g.getId()) + "," + numExperiments
					+ "," + numSamples + "," + meanrank);
		}

	}

	private void outputAll_geneCorank_probelevel(
			Collection<ExpressionExperiment> eeCol, Collection<Gene> pars,
			Collection<Gene> genes, PrintStream p) {

		Collection<Gene> allProbesGenes = new ArrayList<Gene>();
		allProbesGenes.addAll(pars);
		allProbesGenes.addAll(genes);

		/*
		 * for ( ExpressionExperiment ee : expressRankings.keySet() ) { Map<Gene,
		 * Collection<Double>> expressRankingsExperiment =
		 * expressRankings.get(ee);
		 * 
		 * Iterator<Gene> gItr = genes.iterator(); Iterator<Gene> pItr =
		 * pars.iterator();
		 * 
		 * Gene g, p;
		 * 
		 * 
		 * //for ( Gene g : expressRankings.get(ee).keySet() ) { while
		 * (gItr.hasNext()) { g = gItr.next(); p = pItr.next();
		 */

		// use getRanks function to map genes to experiments
		Map<ExpressionExperiment, Map<Gene, Map<DesignElement, Double[]>>> expressRankings = processedExpressionDataVectorService
				.getRanksProbes(eeCol, allProbesGenes);
		// Map<ExpressionExperiment, Map<Gene, Map<DesignElement,Double[]>>>
		// expressRankings = new HashMap<ExpressionExperiment, Map<Gene,
		// Map<DesignElement,Double[]>>>();

		// need to map experiments to genes to calculate max/mean across
		// experiments
		Map<Gene, Double[]> allGeneRankings = new HashMap<Gene, Double[]>();
		Map<Gene, int[]> geneCounts = new HashMap<Gene, int[]>();

		// Map<Gene, Map<DesignElement, Collection<Double>>>
		// allProbeRankingsMeanrank = new HashMap<Gene,
		// Map<DesignElement,Collection<Double>>>();
		// Map<Gene, Map<DesignElement, Collection<Double>>>
		// allProbeRankingsMaxrank = new HashMap<Gene,
		// Map<DesignElement,Collection<Double>>>();

		// Contains all the data for Genes and Pars
		Map<Gene, Map<DesignElement, Collection<Double[]>>> allGeneProbeRankings = new HashMap<Gene, Map<DesignElement, Collection<Double[]>>>();
		Map<Gene, Map<DesignElement, Collection<Double[]>>> allPARProbeRankings = new HashMap<Gene, Map<DesignElement, Collection<Double[]>>>();

		for (ExpressionExperiment ee : expressRankings.keySet()) {
			Map<Gene, Map<DesignElement, Double[]>> expressRankingsGene = expressRankings
					.get(ee);

			Iterator<Gene> gItr = genes.iterator();
			Iterator<Gene> pItr = pars.iterator();

			Gene gene, par;

			// for ( Gene g : expressRankingsGene.keySet() ) {
			while (gItr.hasNext()) {
				gene = gItr.next();
				par = pItr.next();

				if ((!expressRankingsGene.containsKey(gene))
						|| (!expressRankingsGene.containsKey(par)))
					continue;

				Map<DesignElement, Double[]> expressRankingsGeneCompseq = expressRankingsGene
						.get(gene);
				Map<DesignElement, Double[]> expressRankingsPARCompseq = expressRankingsGene
						.get(par);

				// take all the samples for that gene in this experiment
				for (DesignElement de : expressRankingsGeneCompseq.keySet()) {

					// Double d = expressRankings.get(ee).get(g).get(de);
					Double[] d = expressRankingsGeneCompseq.get(de);

					// subtract this entry from the list of probes, and if
					// necessary, the entire experiment
					if (d == null) {
						// numNullRanks++;
						// System.out.println("Null rank value "+gene.getId()+",
						// with probeID: "+ de.getId());

					} else {

						if (!allGeneProbeRankings.containsKey(par)) {
							allGeneProbeRankings
									.put(
											par,
											new HashMap<DesignElement, Collection<Double[]>>());

							// allProbeRankings.put(g, new
							// HashMap<DesignElement,Collection<Double>>());
							// allProbeRankingsMaxrank.put( g, new
							// HashMap<DesignElement,Collection<Double>>());
						}

						if (!allGeneProbeRankings.get(par).containsKey(de)) {
							allGeneProbeRankings.get(par).put(de,
									new ArrayList<Double[]>());
							// allProbeRankingsMaxrank.get( g).put(de, new
							// ArrayList<Double>());
						}

						allGeneProbeRankings.get(par).get(de).add(d);
						// allProbeRankingsMaxrank.get( g).get(de).add(d[1]);

						// calculate the maximum and average rankings
						// maxAndAveRank[1] += d;
						// if (maxAndAveRank[0] < d) maxAndAveRank[0] = d;

						// allNull = false;
					}

				}

				// //////////////////////////////////// repeat for pars, doesn't
				// make sense to do this at probe
				// level!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				// //////////////////////////////////////////////////

				// check if genes for this experiment are null-ranked. if so,
				// skip
				// if (allNull) {
				// System.out.println("Null: No exps for this gene
				// "+g.getId()+", exp: "+ ee.getId());
				// continue;
				// }

				// maxAndAveRank[1] = maxAndAveRank[1] /
				// (expressRankings.get(ee).get(g).size() - numNullRanks);

				// update hash of expression levels
				// if (allGeneRankings.containsKey(g)) {
				// if (allGeneRankings.get(g)[1] < maxAndAveRank[0])
				// allGeneRankings.get(g)[1] = maxAndAveRank[0];
				// allGeneRankings.get(g)[1] += maxAndAveRank[1];
				//	
				// geneCounts.get(g)[0]++;
				// geneCounts.get(g)[1] +=
				// expressRankings.get(ee).get(g).size();
				// }
				// new gene, must create new entities in hash
				// else {
				// allGeneRankings.put(g, maxAndAveRank);
				//	
				//	
				// // number of experiments and number samples for this gene
				// int[] counts = new int[2];
				// counts[0] = 1;
				// counts[1] = expressRankings.get(ee).get(g).size();
				//	
				// geneCounts.put(g, counts);
				// }

			}
		}

		// allProbeRankings.get(g).put(de, new ArrayList<Double>());
		for (Gene g : allGeneProbeRankings.keySet()) {

			Map<DesignElement, Collection<Double[]>> allProbeRankingsGenes = allGeneProbeRankings
					.get(g);

			// int numExperiments = geneCounts.get(g)[0];
			// int numSamples = geneCounts.get(g)[1];
			//
			// double maxrank = allGeneRankings.get(g)[0];
			// double meanrank = allGeneRankings.get(g)[1] / numExperiments;
			//
			// //px.println(parFileEntries.get(g.getId())+","+numExperiments+","+numSamples+","+maxrank);
			// //pe.println(parFileEntries.get(g.getId())+","+numExperiments+","+numSamples+","+meanrank);

			for (DesignElement de : allProbeRankingsGenes.keySet()) {

				// element 0 is the mean, element 1 is the max ranking

				// from list of probe expressios, do stats

				int numOfExperiments = allProbeRankingsGenes.get(de).size();

				double valExpMean_RankMean = 0;
				double valExpMean_RankMax = 0;
				double valExpMax_RankMean = 0;
				double valExpMax_RankMax = 0;

				for (Double[] entry : allProbeRankingsGenes.get(de)) {

					// if (entry[0] == null || entry[1] == null) continue;

					double valRankMean = entry[0].doubleValue();

					valExpMean_RankMean += valRankMean;
					if (valExpMax_RankMean < valRankMean)
						valExpMax_RankMean = valRankMean;

					double valRankMax = entry[1].doubleValue();

					valExpMean_RankMax += valRankMax;
					if (valExpMax_RankMax < valRankMax) {
						valExpMax_RankMax = valRankMax;
					}

				}

				valExpMean_RankMean = valExpMean_RankMean / numOfExperiments;
				valExpMean_RankMax = valExpMean_RankMax / numOfExperiments;

				p.println(parFileEntries.get(g.getId()) + "," + de.getId()
						+ "," + numOfExperiments + "," + valExpMean_RankMean
						+ "," + valExpMean_RankMax + "," + valExpMax_RankMean
						+ "," + valExpMax_RankMax);

				// px.println(g.getId()+","+numExperiments+","+numSamples+","+maxrank);
				// pe.println(g.getId()+","+numExperiments+","+numSamples+","+meanrank);
			}

		}

	}

	private void outputRankingsProbelevel(Collection<ExpressionExperiment> eeCol,
			Collection<Gene> pars, PrintStream p) {

		// System.out.println("reading probes");
		// for (ExpressionExperiment e : eeCol) {System.out.println(e.getId());}
		// System.out.println("Genes");
		// for (Gene e : pars) {System.out.println(e.getId());}


		// RankMethod method;
		//
		// // todo: make these boolean
		// if (rankMethodStr.equalsIgnoreCase("max")) {
		// method = RankMethod.max;
		// } else {
		// method = RankMethod.mean;
		// }

		// use getRanks function to map genes to experiments
		Map<ExpressionExperiment, Map<Gene, Map<DesignElement, Double[]>>> expressRankings = processedExpressionDataVectorService
				.getRanksProbes(eeCol, pars);
		// Map<ExpressionExperiment, Map<Gene, Map<DesignElement,Double[]>>>
		// expressRankings = new HashMap<ExpressionExperiment, Map<Gene,
		// Map<DesignElement,Double[]>>>();
		// Map<ExpressionExperiment, Map<Gene, Collection<Double>>>
		// expressRankings =
		// processedExpressionDataVectorService.getRanks(eeCol, pars,
		// RankMethod.mean);

		// Map<ExpressionExperiment, Map<Gene, Collection<Double>>>
		// expressRankings =
		// processedExpressionDataVectorService.getRanks(eeCol, pars, method);

		// need to map experiments to genes to calculate max/mean across
		// experiments
		Map<Gene, Double[]> allGeneRankings = new HashMap<Gene, Double[]>();
		Map<Gene, int[]> geneCounts = new HashMap<Gene, int[]>();

		// Map<Gene, Map<DesignElement, Collection<Double>>>
		// allProbeRankingsMeanrank = new HashMap<Gene,
		// Map<DesignElement,Collection<Double>>>();
		// Map<Gene, Map<DesignElement, Collection<Double>>>
		// allProbeRankingsMaxrank = new HashMap<Gene,
		// Map<DesignElement,Collection<Double>>>();
		Map<Gene, Map<DesignElement, Collection<Double[]>>> allProbeRankings = new HashMap<Gene, Map<DesignElement, Collection<Double[]>>>();
		// Map<Gene, Map<DesignElement,Integer>> probeCounts = new HashMap<Gene,
		// Map<DesignElement,Integer>>();

		// if (expressRankings.keySet().isEmpty()) {
		// System.out.println("Batch is empty");
		// }
		for (ExpressionExperiment ee : expressRankings.keySet()) {
			// System.out.println(ee.getId());
			Map<Gene, Map<DesignElement, Double[]>> expressRankingsGene = expressRankings
					.get(ee);
			for (Gene g : expressRankingsGene.keySet()) {
				// System.out.println(ee.getId()+"t"+g.getId());
				Map<DesignElement, Double[]> expressRankingsGeneCompseq = expressRankingsGene
						.get(g);

				Collection<DesignElement> cde;

				// Add filter here to obtain only unique probes
				if (checkUniqueProbeMappings) {
					// System.out.println("checking probe mappings");
					Collection<CompositeSequence> csCasted = new ArrayList<CompositeSequence>();

					// cast all composite sequences as design elements
					for (DesignElement de : expressRankingsGeneCompseq.keySet()) {
						// System.out.println("casting de's: " + de.getId());
						csCasted.add((CompositeSequence) de);
					}

					// use the getGenesWithSpecificity service to obtain all
					// physical locations - and other information not needed
					Map<CompositeSequence, Map<PhysicalLocation, Collection<BlatAssociation>>> gws = compositeSequenceService
							.getGenesWithSpecificity(csCasted);
					cde = new ArrayList<DesignElement>();

					for (CompositeSequence cs : gws.keySet()) {

						int numLocations = gws.get(cs).keySet().size();
						if (numLocations == 1) {
							// System.out.println("CS unique! "+cs.getId());
							cde.add(cs);
						} else {
							// System.out.println("CS is not unique!
							// "+cs.getId() +"\tTotal Locations: "+
							// numLocations);
						}

						// System.out.println(cs);
						// for (PhysicalLocation pl : gws.get(cs).keySet()) {
						// //System.out.println("\t"+pl);
						// for (BlatAssociation ba : gws.get(cs).get(pl)) {
						// blatAssociationService.thaw(ba);
						// System.out.println("\t\t"+ba);
						// }
						// }
					}

					// cccz.add(cz);
					// cccz.add(compositeSequenceService.load(new
					// Long(2365242)));
					// cccz.add(compositeSequenceService.load(new
					// Long(406118)));
					// Map<CompositeSequence, Map<PhysicalLocation,
					// Collection<BlatAssociation>>> ddd =
					// compositeSequenceService.getGenesWithSpecificity(cccz);
					// compositeSequenceService.get

					// for (CompositeSequence cs : ddd.keySet()) {
					// System.out.println(cs);
					// for (PhysicalLocation pl : ddd.get(cs).keySet()) {
					// System.out.println("\t"+pl);
					// for (BlatAssociation ba : ddd.get(cs).get(pl)) {
					// blatAssociationService.thaw(ba);
					// System.out.println("\t\t"+ba);
					// }
					// }
					// }

				}

				else {
					cde = expressRankingsGeneCompseq.keySet();
				}

				// double[] maxAndAveRank = new double[2];
				// maxAndAveRank[0] = 0; // max
				// maxAndAveRank[1] = 0; // mean

				// int numNullRanks = 0; // keep track of how many null rankings
				// there are
				// boolean allNull = true;

				// for ( Double d : expressRankings.get(ee).get(g) ) {
				// for ( DesignElement de :
				// expressRankings.get(ee).get(g).keySet() ) {
				for (DesignElement de : cde) {
					// System.out.println("Checking DE: " + de.getId());
					// if (checkUniqueProbeMappings) {
					// CompositeSequence cse = (CompositeSequence) de;
					//						
					// Map<CompositeSequence, Map<PhysicalLocation,
					// Collection<BlatAssociation>>> ddd =
					// compositeSequenceService.getGenesWithSpecificity(cse);
					// }

					// Double d = expressRankings.get(ee).get(g).get(de);
					Double[] d = expressRankingsGeneCompseq.get(de);

					// subtract this entry from the list of probes, and if
					// necessary, the entire experiment
					if (d == null) {
						// numNullRanks++;
						// System.out.println("Null rank value "+g.getId()+",
						// with probeID: "+ de.getId());

					} else {

						// System.out.println(ee.getId()+"\t"+g.getId()+"\t"+d[0]);
						if (!allProbeRankings.containsKey(g)) {
							allProbeRankings
									.put(
											g,
											new HashMap<DesignElement, Collection<Double[]>>());

							// allProbeRankings.put(g, new
							// HashMap<DesignElement,Collection<Double>>());
							// allProbeRankingsMaxrank.put( g, new
							// HashMap<DesignElement,Collection<Double>>());
						}

						if (!allProbeRankings.get(g).containsKey(de)) {
							allProbeRankings.get(g).put(de,
									new ArrayList<Double[]>());
							// allProbeRankingsMaxrank.get( g).put(de, new
							// ArrayList<Double>());
						}

						allProbeRankings.get(g).get(de).add(d);
						// allProbeRankingsMaxrank.get( g).get(de).add(d[1]);

						// calculate the maximum and average rankings
						// maxAndAveRank[1] += d;
						// if (maxAndAveRank[0] < d) maxAndAveRank[0] = d;

						// allNull = false;
					}

				}

				// check if genes for this experiment are null-ranked. if so,
				// skip
				// if (allNull) {
				// System.out.println("Null: No exps for this gene
				// "+g.getId()+", exp: "+ ee.getId());
				// continue;
				// }

				// maxAndAveRank[1] = maxAndAveRank[1] /
				// (expressRankings.get(ee).get(g).size() - numNullRanks);

				// update hash of expression levels
				// if (allGeneRankings.containsKey(g)) {
				// if (allGeneRankings.get(g)[1] < maxAndAveRank[0])
				// allGeneRankings.get(g)[1] = maxAndAveRank[0];
				// allGeneRankings.get(g)[1] += maxAndAveRank[1];
				//	
				// geneCounts.get(g)[0]++;
				// geneCounts.get(g)[1] +=
				// expressRankings.get(ee).get(g).size();
				// }
				// new gene, must create new entities in hash
				// else {
				// allGeneRankings.put(g, maxAndAveRank);
				//	
				//	
				// // number of experiments and number samples for this gene
				// int[] counts = new int[2];
				// counts[0] = 1;
				// counts[1] = expressRankings.get(ee).get(g).size();
				//	
				// geneCounts.put(g, counts);
				// }

			}
		}

		// allProbeRankings.get(g).put(de, new ArrayList<Double>());
		for (Gene g : allProbeRankings.keySet()) {

			Map<DesignElement, Collection<Double[]>> allProbeRankingsGenes = allProbeRankings
					.get(g);

			// int numExperiments = geneCounts.get(g)[0];
			// int numSamples = geneCounts.get(g)[1];
			//
			// double maxrank = allGeneRankings.get(g)[0];
			// double meanrank = allGeneRankings.get(g)[1] / numExperiments;
			//
			// //px.println(parFileEntries.get(g.getId())+","+numExperiments+","+numSamples+","+maxrank);
			// //pe.println(parFileEntries.get(g.getId())+","+numExperiments+","+numSamples+","+meanrank);

			for (DesignElement de : allProbeRankingsGenes.keySet()) {

				// element 0 is the mean, element 1 is the max ranking

				// from list of probe expressios, do stats

				int numOfExperiments = allProbeRankingsGenes.get(de).size();

				double valExpMean_RankMean = 0;
				double valExpMean_RankMax = 0;
				double valExpMax_RankMean = 0;
				double valExpMax_RankMax = 0;

				for (Double[] entry : allProbeRankingsGenes.get(de)) {

					// if (entry[0] == null || entry[1] == null) continue;

					double valRankMean = entry[0].doubleValue();

					valExpMean_RankMean += valRankMean;
					if (valExpMax_RankMean < valRankMean)
						valExpMax_RankMean = valRankMean;

					double valRankMax = entry[1].doubleValue();

					valExpMean_RankMax += valRankMax;
					if (valExpMax_RankMax < valRankMax) {
						valExpMax_RankMax = valRankMax;
					}

				}

				valExpMean_RankMean = valExpMean_RankMean / numOfExperiments;
				valExpMean_RankMax = valExpMean_RankMax / numOfExperiments;

				// / TEMPORARY TO GET RANKS FOR ALL GENES!!!
				// p.println(g.getId()+","+g.getName()+","+de.getId()+","+numOfExperiments+","+valExpMean_RankMean+","+valExpMean_RankMax+","+valExpMax_RankMean+","+valExpMax_RankMax);
				p.println(parFileEntries.get(g.getId()) + "," + de.getId()
						+ "," + numOfExperiments + "," + valExpMean_RankMean
						+ "," + valExpMean_RankMax + "," + valExpMax_RankMean
						+ "," + valExpMax_RankMax);

				// px.println(g.getId()+","+numExperiments+","+numSamples+","+maxrank);
				// pe.println(g.getId()+","+numExperiments+","+numSamples+","+meanrank);
			}

		}

	}

	// the original - look at PARs/Genes as a whole, not distinguishing probes
	private void outputRankings(Collection<ExpressionExperiment> eeCol,
			Collection<Gene> pars, PrintStream px, PrintStream pe,
			String rankMethodStr) {

		RankMethod method;

		// todo: make these boolean
		if (rankMethodStr.equalsIgnoreCase("max")) {
			method = RankMethod.max;
		} else {
			method = RankMethod.mean;
		}

		// use getRanks function to map genes to experiments
		Map<ExpressionExperiment, Map<Gene, Collection<Double>>> expressRankings = processedExpressionDataVectorService
				.getRanks(eeCol, pars, method);

		// need to map experiments to genes to calculate max/mean across
		// experiments
		HashMap<Gene, double[]> allGeneRankings = new HashMap<Gene, double[]>();
		HashMap<Gene, int[]> geneCounts = new HashMap<Gene, int[]>();

		for (ExpressionExperiment ee : expressRankings.keySet()) {
			for (Gene g : expressRankings.get(ee).keySet()) {

				double[] maxAndAveRank = new double[2];
				maxAndAveRank[0] = 0; // max
				maxAndAveRank[1] = 0; // mean

				int numNullRanks = 0; // keep track of how many null rankings
										// there are
				boolean allNull = true;

				for (Double d : expressRankings.get(ee).get(g)) {

					// subtract this entry from the list of probes, and if
					// necessary, the entire experiment
					if (d == null) {
						numNullRanks++;
						// System.out.println("Null rank value "+g.getId()+",
						// numNullRanks: "+ numNullRanks);

					} else {

						// calculate the maximum and average rankings
						maxAndAveRank[1] += d;
						if (maxAndAveRank[0] < d)
							maxAndAveRank[0] = d;

						allNull = false;
					}

				}

				// check if genes for this experiment are null-ranked. if so,
				// skip
				if (allNull) {
					// System.out.println("Null: No exps for this gene
					// "+g.getId()+", exp: "+ ee.getId());
					continue;
				}

				maxAndAveRank[1] = maxAndAveRank[1]
						/ (expressRankings.get(ee).get(g).size() - numNullRanks);

				// update hash of expression levels
				if (allGeneRankings.containsKey(g)) {
					if (allGeneRankings.get(g)[0] < maxAndAveRank[0])
						allGeneRankings.get(g)[0] = maxAndAveRank[0];
					allGeneRankings.get(g)[1] += maxAndAveRank[1];

					geneCounts.get(g)[0]++;
					geneCounts.get(g)[1] += expressRankings.get(ee).get(g)
							.size();
				}
				// new gene, must create new entities in hash
				else {
					allGeneRankings.put(g, maxAndAveRank);

					// number of experiments and number samples for this gene
					int[] counts = new int[2];
					counts[0] = 1;
					counts[1] = expressRankings.get(ee).get(g).size();

					geneCounts.put(g, counts);
				}

			}
		}

		for (Gene g : allGeneRankings.keySet()) {

			int numExperiments = geneCounts.get(g)[0];
			int numSamples = geneCounts.get(g)[1];

			double maxrank = allGeneRankings.get(g)[0];
			double meanrank = allGeneRankings.get(g)[1] / numExperiments;

			String label = (String) parFileEntries.get(g.getId());
			if (label == null)
				label = g.getId().toString();

			px.println(label + "," + numExperiments + "," + numSamples + ","
					+ maxrank);
			pe.println(label + "," + numExperiments + "," + numSamples + ","
					+ meanrank);

			// px.println(g.getId()+","+numExperiments+","+numSamples+","+maxrank);
			// pe.println(g.getId()+","+numExperiments+","+numSamples+","+meanrank);

		}

	}

	private static HashMap getIndices(String header) {
		String[] labels = header.trim().split("\t");
		HashMap<String, Integer> hash = new HashMap<String, Integer>(
				labels.length);

		for (int i = 0; i < labels.length; i++) {
			hash.put(labels[i], i);
		}

		return hash;
	}

	private static int getIndex(String header) {
		return ((Integer) headerLookup.get(header)).intValue();
	}
	
	/*
	 * This method reads the PAR file, the main source of information for this CLI. 
	 * The PAR file is created by PARMapper.java
	 * 
	 * The file format should be tabe delimited file like the following:
	 * ParID	ParName	Chrom	Nuc	GeneId	GeneSymbol	Distance	GeneContainsPar	SameStrand
	 * 1124999	R49730.par.8.239068.239461	8	239068	3887948	LOC100131718	462	false	true
	 * 1125323	AA431402.par.9.35913976.35914288	9	35913976	400013	LOC158376	12294	false	false
	 * 1126681	W90363.par.2.215298022.215298505	2	215298022	12860	BARD1	3001	false	true
	 * 
	 */
	private static void readPARFile(String inFile) {

		BufferedReader in;
		Collection<String[]> fRecords = new ArrayList<String[]>();
		HashMap fHash = null;
		String[] fHeaders = null;

		String header;

		try {
			in = new BufferedReader(new FileReader(inFile));
			String line;
			header = in.readLine();
			fHash = getIndices(header);
			fHeaders = header.trim().split("\t");

			while ((line = in.readLine()) != null) {
				if (line.startsWith("#"))
					continue;
				String[] s = line.trim().split("\t");
				fRecords.add(s);
			}

			in.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.out.println("File " + inFile + " not found - "
					+ e.getMessage());
			System.exit(0);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.out.println("File " + inFile + " reads a bit wonky - "
					+ e.getMessage());
			System.exit(0);
		}

		headerLookup = fHash;
//		headers = fHeaders;
		records = fRecords;
	}

	// reads a par/gene/experiments file and saves the results to hash
	// parToCoExps
	private static void readpargeneFile(String inFile) {

		parToCoExps = new HashMap<Long, long[]>();

		BufferedReader in;
//		Collection<String[]> fRecords = new ArrayList<String[]>();

		String header;

		try {
			in = new BufferedReader(new FileReader(inFile));
			String line;
			header = in.readLine();

			while ((line = in.readLine()) != null) {
				if (line.startsWith("#"))
					continue;
				String[] s = line.trim().split("\t");
				// fRecords.add(s);
				String[] expIdsStr = s[2].split(",");
				long[] expIds = new long[expIdsStr.length];

				for (int i = 0; i < expIdsStr.length; i++) {
					expIds[i] = Long.parseLong(expIdsStr[i]);
				}

				parToCoExps.put(new Long(Long.parseLong(s[0])), expIds);
				// parToCoExps.put(s[0], expIds);
			}

			in.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.out.println("File " + inFile + " not found - "
					+ e.getMessage());
			System.exit(0);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.out.println("File " + inFile + " reads a bit wonky - "
					+ e.getMessage());
			System.exit(0);
		}

	}
	
	/*
	 * This method outputs the pairwise correlations between pars and genes
	 */
	private void outputAllPARGeneCoexpressionCalculated(Collection<Gene> pars, Collection<Gene> genes, PrintStream pap) {
		Iterator<Gene> pItr = pars.iterator();
		Iterator<Gene> gItr = genes.iterator();

		while (pItr.hasNext()) {
			Gene par = pItr.next();
			Gene gene = gItr.next();

			long[] expIds = parToCoExps.get(par.getId());

			if (expIds == null) {
				//System.out.println("No information: " + par.getId());
				continue;
			}

			Collection<ExpressionExperiment> exps = new ArrayList<ExpressionExperiment>();

			for (int i = 0; i < expIds.length; i++) {
				exps.add(expressionExperimentService.load(expIds[i]));
			}

			Collection<Gene> pargene = new ArrayList<Gene>();
			pargene.add(par);
			pargene.add(gene);

			if (pargene == null) {
				System.out.println("No pargene");
				continue;
			}
			if (exps == null) {
				System.out.println("No exps");
				continue;
			}
			
			Map<Long, Collection<DoubleVectorValueObject>> expParMap  = new HashMap<Long, Collection<DoubleVectorValueObject>>();
			Map<Long, Collection<DoubleVectorValueObject>> expGeneMap = new HashMap<Long, Collection<DoubleVectorValueObject>>();
			
//System.out.println("Calling API");
			Collection<DoubleVectorValueObject> dvvos = processedExpressionDataVectorService.getProcessedDataArrays(exps, pargene);
//System.out.println("Finished API");
			
			// for ever probe, assign it to either par or gene for the experiment
			for (DoubleVectorValueObject dvvo : dvvos) {
				ExpressionExperiment ee = dvvo.getExpressionExperiment();
				expressionExperimentService.thawLite(ee);
				Long eeId = new Long(ee.getId());
				ee = null;
				
//System.out.println(ee.getId());
				//try {
				//	expressionExperimentService.thaw(ee);
				//} catch (ubic.gemma.model.expression.experiment.ExpressionExperimentServiceException e) {
				//	System.out.println("Thaw failed for " + ee);
				//	continue;
				//}
//System.out.println("Dvvo" + dvvo.toString());
				
				Collection<Gene> genelist = dvvo.getGenes();
				
				boolean isParGene = false;
				for (Gene g: genelist) {
//System.out.println("\t" + g.getId() + "\t" + g.getName());
					long id = g.getId();
					if (id == gene.getId()) {
						//geneData.put(dvvo.getDesignElement(), dvvo.getData());
						if (!expGeneMap.containsKey(eeId)) {
							expGeneMap.put(eeId, new ArrayList<DoubleVectorValueObject>());
						}
						expGeneMap.get(eeId).add(dvvo);
						isParGene = true;
						break;
					} else if (id == par.getId()) {
						//parData.put(dvvo.getDesignElement(), dvvo.getData());
						if (!expParMap.containsKey(eeId)) {
							expParMap.put(eeId, new ArrayList<DoubleVectorValueObject>());
						}
						expParMap.get(eeId).add(dvvo);
						isParGene = true;
						break;
					}
				}
				if (isParGene == false) {
					System.out.println("Probe not assigned, skipping: "+ dvvo.getDesignElement()+ " (Experiment: "
							+ dvvo.getExpressionExperiment().getId() + ")");
					continue;
				}
				
				
				
				
			}
			
			for (Long eeId: expParMap.keySet()) {
				
				
				
				if (! expGeneMap.containsKey(eeId)) {
					System.out.println("No experiment data for gene "+gene.getId()+" in this experiment: "+eeId);
					continue;
				}
				
				//iterate through the pars
				for (DoubleVectorValueObject parDvvo : expParMap.get(eeId)) {
					DesignElement pd = parDvvo.getDesignElement();
					
					//double[] pdata = parData.get(pd);
					double[] pdata = parDvvo.getData();
					
					int parProbeMappings  = parDvvo.getGenes().size();
					
					
					//iterate through the genes
					for (DoubleVectorValueObject geneDvvo : expGeneMap.get(eeId)) {
						DesignElement gd = geneDvvo.getDesignElement();
						
						//double[] gdata = geneData.get(gd);
						double[] gdata = geneDvvo.getData();
						
						int geneProbeMappings = geneDvvo.getGenes().size();
						
						// check if vectors the same length
						if (pdata.length != gdata.length) {
							System.out.println("Error: vector sizes do not match: ParProbeID:"
									+pd.getId()+"\tGeneProbeID:"+gd.getId()
									+"\tSizes:"+ pdata.length +"\t"+gdata.length);
							continue;
						}
						
						int numbadIndices = 0;
						boolean[] badIndices = new boolean[pdata.length];
						for (int i=0; i<pdata.length; i++) {
							badIndices[i] = false;
						}
						
						// check data integrity - remove missing values
						//boolean integrity = true;
						for (int i=0; i<pdata.length; i++) {
							if ((new Double(pdata[i])).isNaN() || (new Double(gdata[i])).isNaN()) {
								badIndices[i] = true;
								numbadIndices++;
								//integrity = false;
								//break;
							}
						}
						
						if (pdata.length - numbadIndices < 3) {
							System.out.println("Not enough samples ("
									+ (pdata.length - numbadIndices)
									+") to calc coexpression: ParProbeID:"
									+pd.getId()+"\tGeneProbeID:"+gd.getId());
						}
						
						if (0 < numbadIndices) {
							double[] newpdata = new double[pdata.length - numbadIndices];
							double[] newgdata = new double[pdata.length - numbadIndices];
							
							int j=0;
							for (int i=0; i<pdata.length; i++) {
								//System.out.println(i+"--"+j+"__"+badIndices[i]);
								if (badIndices[i]) {
									//System.out.println("skipping...");
									//continue;
								} else {
									newpdata[j] = pdata[i];
									newgdata[j] = gdata[i];
									j++;
								}
							}
							
							pdata = newpdata;
							gdata = newgdata;
						}
						
						
						
						
						// do co-expression with these two pairs!
						DoubleArrayList pdal = new DoubleArrayList(pdata);
						DoubleArrayList gdal = new DoubleArrayList(gdata);
						double corr = ubic.basecode.math.DescriptiveWithMissing.correlation(pdal, gdal);
						
						pap.println(parFileEntries.get(par.getId()) 
								+ "," + eeId
								+ "," + parProbeMappings
								+ "," + geneProbeMappings
								+ "," + pdata.length
								+ "," + pd.getId()
								+ "," + gd.getId()
								+ "," + corr
								);


						
					}
					
					// just to save some memory when finished with this
					expParMap.put(eeId, null);
					
					//System.out.println();
				}
			}
		}
		
		
		

	}
	
	
	
	/*
	 * This method outputs the pairwise correlations between pars and genes
	 * this was found to be too slow when calling on getProcessedDataArrays for each experiment.  It should
	 * input all experiments instead of feeding them one by one
	 */
	private void out_pargene_allcoexpressions_orig_slow(Collection<Gene> pars, Collection<Gene> genes, PrintStream pap) {
		Iterator<Gene> pItr = pars.iterator();
		Iterator<Gene> gItr = genes.iterator();

		while (pItr.hasNext()) {
			Gene par = pItr.next();
			Gene gene = gItr.next();

			long[] expIds = parToCoExps.get(par.getId());

			if (expIds == null) {
				//System.out.println("No information: " + par.getId());
				continue;
			}

			Collection<ExpressionExperiment> exps = new ArrayList<ExpressionExperiment>();
			//Collection<BioAssaySet> exps2 = new ArrayList<BioAssaySet>();

			for (int i = 0; i < expIds.length; i++) {
				// gene2GeneCoexpressionService.
				// GeneCoexpressionAnalysis gca;
				// Gene2GeneCoexpression g = getNewGGCOInstance();
				// g.setFirstGene(par);
				// g.setSecondGene(gene);
				// Double d = g.getPvalue();
				// if (d==null) d = new Double(-5);
				exps.add(expressionExperimentService.load(expIds[i]));
				//exps2.add((BioAssaySet) expressionExperimentService
				//		.load(expIds[i]));
			}
			// System.out.println();

			Collection<Gene> pargene = new ArrayList<Gene>();
			pargene.add(par);
			pargene.add(gene);

			if (pargene == null) {
				System.out.println("No pargene");
				continue;
			}
			if (exps == null) {
				System.out.println("No exps");
				continue;
			}
			
//			Collection<ExpressionExperiment> ees = new ArrayList<ExpressionExperiment>();
			//Collection<Gene> gs = new ArrayList<Gene>();
			
			//Gene pp = parService.load((long) 1126681);
			//Gene gg = parService.load((long) 12860);
			
			// just a random par/gene pair
			//ees.add(expressionExperimentService.load((long) 441));
			//ees.add(expressionExperimentService.load((long) 442));
			//gs.add(pp);
			//gs.add(gg);
			
			// has link data
			//ees.add(expressionExperimentService.load((long) 1));
			//ees.add(expressionExperimentService.load((long) 4));
			//ees.add(expressionExperimentService.load((long) 111));
			//ees.add(expressionExperimentService.load((long) 468));
			//ees.add(expressionExperimentService.load((long) 440));
			//gs.add(parService.load((long) 1445507));
			//gs.add(parService.load((long) 442313));
			// ---> weird because these probes maps to MANY other genes (like 10 to 30+)
			
			
			// just a random par/gene pair
			//ees.add(expressionExperimentService.load((long) 167));
			//gs.add(parService.load((long) 1186837));
			//gs.add(parService.load((long) 4641158));
			
			for (ExpressionExperiment ee: exps) {
				Collection<ExpressionExperiment> ees = new ArrayList<ExpressionExperiment>();
				ees.add(ee);
				
				
				Collection<DoubleVectorValueObject> dvvos = processedExpressionDataVectorService.getProcessedDataArrays(ees, pargene);
				
				Map<DesignElement, double[]> parData  = new HashMap<DesignElement, double[]>();
				Map<DesignElement, double[]> geneData = new HashMap<DesignElement, double[]>();
				
				Map<DesignElement, Integer> probemappingCount = new HashMap<DesignElement, Integer>();
				
				for (DoubleVectorValueObject dvvo : dvvos) {
					//System.out.println(">>");
					//System.out.println("Probe: " + dvvo.getDesignElement().getName() +" - "+ dvvo.getDesignElement().getId());
					Collection<Gene> genelist = dvvo.getGenes();
					
					//System.out.println("Genes:");
					boolean isParGene = false;
					for (Gene g: genelist) {
						//System.out.println("\t" + g.getId() + "\t" + g.getName());
						long id = g.getId();
						if (id == gene.getId()) {
							geneData.put(dvvo.getDesignElement(), dvvo.getData());
							isParGene = true;
							break;
						} else if (id == par.getId()) {
							parData.put(dvvo.getDesignElement(), dvvo.getData());
							isParGene = true;
							break;
						}
					}
					if (isParGene == false) {
						System.out.println("Probe not assigned, skipping: "+ dvvo.getDesignElement()+ " (Experiment: "
								+ dvvo.getExpressionExperiment().getId() + ")");
						continue;
					}
					probemappingCount.put(dvvo.getDesignElement(), new Integer(dvvo.getGenes().size()));
				}
				
				for (DesignElement pd: parData.keySet()) {
					double[] pdata = parData.get(pd);
					
					for (DesignElement gd: geneData.keySet()) {
						double[] gdata = geneData.get(gd);
						
						// check if vectors the same length
						if (pdata.length != gdata.length) {
							System.out.println("Error: vector sizes do not match: ParProbeID:"
									+pd.getId()+"\tGeneProbeID:"+gd.getId()
									+"\tSizes:"+ pdata.length +"\t"+gdata.length);
							continue;
						}
						
						int numbadIndices = 0;
						boolean[] badIndices = new boolean[pdata.length];
						for (int i=0; i<pdata.length; i++) {
							badIndices[i] = false;
						}
						
						// check data integrity - remove missing values
						//boolean integrity = true;
						for (int i=0; i<pdata.length; i++) {
							if ((new Double(pdata[i])).isNaN() || (new Double(gdata[i])).isNaN()) {
								badIndices[i] = true;
								numbadIndices++;
								//integrity = false;
								//break;
							}
						}
						
						if (pdata.length - numbadIndices < 3) {
							System.out.println("Not enough samples ("
									+ (pdata.length - numbadIndices)
									+") to calc coexpression: ParProbeID:"
									+pd.getId()+"\tGeneProbeID:"+gd.getId());
						}
						
						if (0 < numbadIndices) {
							double[] newpdata = new double[pdata.length - numbadIndices];
							double[] newgdata = new double[pdata.length - numbadIndices];
							
							int j=0;
							for (int i=0; i<pdata.length; i++) {
								//System.out.println(i+"--"+j+"__"+badIndices[i]);
								if (badIndices[i]) {
									//System.out.println("skipping...");
									//continue;
								} else {
									newpdata[j] = pdata[i];
									newgdata[j] = gdata[i];
									j++;
								}
							}
							/*
							System.out.println("Old");
							for (int i=0; i<pdata.length; i++) {
								System.out.print("\t"+pdata[i]);
							}
							System.out.println();
							for (int i=0; i<pdata.length; i++) {
								System.out.print("\t"+gdata[i]);
							}
							System.out.println();
							for (int i=0; i<pdata.length; i++) {
								System.out.print("\t"+badIndices[i]);
							}
							System.out.println();
							
							System.out.println("New");
							for (int i=0; i<newpdata.length; i++) {
								System.out.print("\t"+newpdata[i]);
							}
							System.out.println();
							for (int i=0; i<newpdata.length; i++) {
								System.out.print("\t"+newgdata[i]);
							}
							System.out.println();
							System.exit(0);*/
							pdata = newpdata;
							gdata = newgdata;
						}
						
						
						//if (integrity != true) {
						//	System.out.println("Dataset cannot be used, has missing values: ParProbeID:"
						//		+pd.getId()+"\tGeneProbeID:"+gd.getId());
						//	continue;
						//}
						
						
						// do co-expression with these two pairs!
						DoubleArrayList pdal = new DoubleArrayList(pdata);
						DoubleArrayList gdal = new DoubleArrayList(gdata);
						double corr = ubic.basecode.math.DescriptiveWithMissing.correlation(pdal, gdal);
						
						//System.out.println(pd.getId() + "\t" + gd.getId()+ "\t" + corr);
						//NumGenesMapProbe,ParProbe,GeneProbe,Pearson");
						int parProbeMappings  = probemappingCount.get(pd).intValue();
						int geneProbeMappings = probemappingCount.get(gd).intValue();
						
						pap.println(parFileEntries.get(par.getId()) 
								+ "," + ee.getId()
								+ "," + parProbeMappings
								+ "," + geneProbeMappings
								+ "," + pdata.length
								+ "," + pd.getId()
								+ "," + gd.getId()
								+ "," + corr
								);


						
					}
					
					
					//System.out.println();
				}
			}
		}
		
		
		

	}
	
	private Collection<ExpressionExperiment> loadExpressionExperimentsByFile (String ExperimentListFile) {

		Collection<ExpressionExperiment> eeCol = new ArrayList<ExpressionExperiment>();
		Collection<Long> eeids = new ArrayList<Long>();
		
		BufferedReader in;

		try {
			in = new BufferedReader(new FileReader(ExperimentListFile));
			String line;
			
			while ((line = in.readLine()) != null) {
				if (line.startsWith("#"))
					continue;
				//String[] s = line.trim().split("\t");
				String expIDstr = line.trim();
				
				eeids.add(Long.decode(expIDstr));
			}
			
			eeCol = expressionExperimentService.loadMultiple(eeids);
			
			in.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.out.println("File " + ExperimentListFile + " not found - "
					+ e.getMessage());
			System.exit(0);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.out.println("File " + ExperimentListFile + " reads a bit wonky - "
					+ e.getMessage());
			System.exit(0);
		}
		
		return eeCol;
	
	}

}
