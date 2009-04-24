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
import java.util.HashSet;
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
import ubic.gemma.model.common.auditAndSecurity.AuditEvent;
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
import ubic.gemma.model.expression.experiment.ExpressionExperimentValueObject;
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
	private boolean checkTroubledExperiments = false;	//
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

		Option taxonOption = OptionBuilder.hasArg().isRequired()
				.withDescription("taxon name").withDescription("taxon to use")
				.withLongOpt("taxon").create('t');
		addOption(taxonOption);

		Option inFileOption = OptionBuilder.hasArg()
				.withDescription("infile name").withDescription("file to read").withLongOpt(
				"inFile").create('i');
		addOption(inFileOption);

		Option outFileDirOption = OptionBuilder.hasArg()
				.withDescription("outdirectory name").withDescription("directory to write to")
				.withLongOpt("outFileDir").create('o');
		addOption(outFileDirOption);

		Option pargeneexpFileOption = OptionBuilder.hasArg()
				.withDescription("PAR-gene-experiment map filename")
				.withDescription("file of Par/Gene/Exps to read").withLongOpt(
				"pargeneexpFile").create('p');
		addOption(pargeneexpFileOption);
		
		Option ExperimentListFileOption = OptionBuilder.hasArg()
				.withDescription("experiment list filename")
				.withDescription("file of expression experiments to read").withLongOpt(
				"pargeneexpFile").create('e');
		addOption(ExperimentListFileOption);
		
		Option outputStemOption = OptionBuilder.hasArg()
				.withDescription("output file stem name")
				.withDescription("output file stem name").withLongOpt(
				"stem").create('s');
		addOption(outputStemOption);

		Option checkParOption = OptionBuilder
				.withDescription("PAR rankings")
				.create("checkPAR");
		addOption(checkParOption);

		Option checkPARprobeOption = OptionBuilder
				.withDescription("PAR probe level rankings")
				.create("checkPARprobe");
		addOption(checkPARprobeOption);

		Option checkGeneCorankOption = OptionBuilder
				.withDescription("PAR and Gene rankings")
				.create("checkGeneCorank");
		addOption(checkGeneCorankOption);

		Option checkGeneRankOption = OptionBuilder
				.withDescription("Gene rankings")
				.create("checkGeneRank");
		addOption(checkGeneRankOption);

		Option checkExperimentTypesOption = OptionBuilder
				.withDescription("Prints detailed information on experiments")
				.create("checkExperimentTypes");
		addOption(checkExperimentTypesOption);
		
		Option checkTroubledExperimentsOption = OptionBuilder
				.withDescription("Determines which experiments are troubled")
				.create("checkTroubledExperiments");
		addOption(checkTroubledExperimentsOption);

		Option checkUniqueProbeMappingsOption = OptionBuilder
				.withDescription("[Optional: used with PARprobe]")
				.create("checkUniqueProbeMappings");
		addOption(checkUniqueProbeMappingsOption);
		
		Option checkAllCoexpOption = OptionBuilder
				.withDescription("Co-expression links for all PARs")
				.create("checkAllCoexp");
		addOption(checkAllCoexpOption);
		
		Option checkAllCoexpInDBOption = OptionBuilder
				.withDescription("Co-expression links in DB for PAR/gene pairs")
				.create("checkAllCoexpInDB");
		addOption(checkAllCoexpInDBOption);

		Option checkAllPairsCoexpOption = OptionBuilder
				.withDescription("All co-expression values for PAR/gene pairs")
				.create("checkAllPairsCoexp");
		addOption(checkAllPairsCoexpOption);

		
	}

	@Override
	protected void processOptions() {
		super.processOptions();

		this.taxonService = (TaxonService) this.getBean("taxonService");
		this.parService = (GeneService) this.getBean("geneService");

		/*
		 * Process arguments
		 */

		// File IO
		if (hasOption("i")) {
			inFile = getOptionValue("i");
		}

		if (hasOption("p")) {
			pargeneexpFile = getOptionValue("p");
		}
		
		if (hasOption("e")) {
			ExperimentListFile = getOptionValue("e");
		}
		
		if (hasOption("o")) {
			outFileDir = getOptionValue("o");
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

		
		this.checkPARprobe = this.hasOption("checkPARprobe");
		this.checkUniqueProbeMappings = this.hasOption("checkUniqueProbeMappings");
		this.checkPAR = this.hasOption("checkPAR");
		this.checkGeneCorank = this.hasOption("checkGeneCorank");
		this.checkGeneRank = this.hasOption("checkGeneRank");
		this.checkAllCoexp = this.hasOption("checkAllCoexp");
		this.checkAllCoexpInDB = this.hasOption("checkAllCoexpInDB");
		this.checkAllPairsCoexp = this.hasOption("checkAllPairsCoexp");
		
		this.checkExperimentTypes = this.hasOption("checkExperimentTypes");
		this.checkTroubledExperiments = this.hasOption("checkTroubledExperiments");

	}

	@Override
	protected Exception doWork(String[] args) {
		// TODO Auto-generated method stub

		Exception exc = processCommandLine("test", args);

		
		// get services used
		this.parService = (GeneService) this.getBean("geneService");
		this.taxonService = (TaxonService) this.getBean("taxonService");
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
		this.probe2ProbeCoexpressionService = (Probe2ProbeCoexpressionService) this
				.getBean("probe2ProbeCoexpressionService");


		
		// Load up experiments - all experiments or from a list
		Collection<ExpressionExperiment> eeCol;
		if (ExperimentListFile == null) {
			eeCol = expressionExperimentService.findByTaxon(taxon);
		} else {
			eeCol = loadExpressionExperimentsByFile(ExperimentListFile);
		}
		System.out.println("Experiments loaded:" + eeCol.size());
		
		
		// test for troubled experiments
		if (checkTroubledExperiments) {
			
			// use expressionexperiment value objects which seem to have the same id as expressionexperiments
			Collection<Long> ees = new HashSet<Long>();
			for (ExpressionExperiment ee : eeCol) {
				ees.add(ee.getId());
			}

	        //int size = ees.size();
	        Map<Long, AuditEvent> trouble = expressionExperimentService.getLastTroubleEvent( ees );
			
	        System.out.println("Called trouble event");
	        
	        for (Long l: ees) {
	        	AuditEvent ae = trouble.get(l);
	        	if (ae == null) {
	        		System.out.println(l +"\tok");
	        	} else {
	        		System.out.println(l +"\ttrouble");
	        	}
	        }
	        
			return null;
		}
		
		
		// prints out experiment information (exits program when finished)
		if (checkExperimentTypes) {
			printExperimentTypes(eeCol);
			return null;
		}
		
		// Read PAR file
		System.out.println("Reading file: " + inFile);
		readPARFile(inFile);
		if (pargeneexpFile != null) readpargeneFile(pargeneexpFile);
		
		
		// Creates output files
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

		
		// Open the file print streams and create the columns
		try {
			// Connect print stream to the output stream
			if (checkPAR || checkGeneRank) {
				pxx = new PrintStream(new FileOutputStream(outFileDir
						+ "/" + stem + ".output.xx.txt"));
				pxe = new PrintStream(new FileOutputStream(outFileDir
						+ "/" + stem + ".output.xe.txt"));
				pex = new PrintStream(new FileOutputStream(outFileDir
						+ "/" + stem + ".output.ex.txt"));
				pee = new PrintStream(new FileOutputStream(outFileDir
						+ "/" + stem + ".output.ee.txt"));
				
				String outputHeader = "ParID,GeneId,ParName,Chrom,Nuc,GeneSymbol,Distance,GeneContainsPar,SameStrand,NumExperiments,NumSamples,Rank";
				pxx.println(outputHeader);
				pxe.println(outputHeader);
				pex.println(outputHeader);
				pee.println(outputHeader);

			}
			
			if (checkPARprobe) {
				//for probe level rannkings
				ppr = new PrintStream(new FileOutputStream(outFileDir
						+ "/" + stem + ".output.probelevel.txt"));
				
				ppr.println("ParID,GeneId,ParName,Chrom,Nuc,GeneSymbol,Distance,GeneContainsPar,SameStrand,CompSeqId,NumExperiments,RankEMean_MethMean,RankEMean_MethMax,RankEMax_MethMean,RankEMax_MethMax");
			}
			
			if (checkGeneCorank) {
				// for Gene vs Par rannkings
				pgp_xx = new PrintStream(new FileOutputStream(outFileDir
						+ "/" + stem + ".output.geneVsPar.xx.txt"));
				pgp_xe = new PrintStream(new FileOutputStream(outFileDir
						+ "/" + stem + ".output.geneVsPar.xe.txt"));
				pgp_ex = new PrintStream(new FileOutputStream(outFileDir
						+ "/" + stem + ".output.geneVsPar.ex.txt"));
				pgp_ee = new PrintStream(new FileOutputStream(outFileDir
						+ "/" + stem + ".output.geneVsPar.ee.txt"));
				
				pgp_xx.println("ParID,GeneId,ParName,Chrom,Nuc,GeneSymbol,Distance,GeneContainsPar,SameStrand,ParNumExperiments,ParNumSamples,ParRank,GeneNumExperiments,GeneNumSamples,GeneRank");
				pgp_xe.println("ParID,GeneId,ParName,Chrom,Nuc,GeneSymbol,Distance,GeneContainsPar,SameStrand,ParNumExperiments,ParNumSamples,ParRank,GeneNumExperiments,GeneNumSamples,GeneRank");
				pgp_ex.println("ParID,GeneId,ParName,Chrom,Nuc,GeneSymbol,Distance,GeneContainsPar,SameStrand,ParNumExperiments,ParNumSamples,ParRank,GeneNumExperiments,GeneNumSamples,GeneRank");
				pgp_ee.println("ParID,GeneId,ParName,Chrom,Nuc,GeneSymbol,Distance,GeneContainsPar,SameStrand,ParNumExperiments,ParNumSamples,ParRank,GeneNumExperiments,GeneNumSamples,GeneRank");
			}
			
			if (checkAllCoexpInDB) {
				pco = new PrintStream(new FileOutputStream(outFileDir
						+ "/" + stem + ".output.genecoexpression.txt"));
				
				pco.println("ParID,GeneId,ParName,Chrom,Nuc,GeneSymbol,Distance,GeneContainsPar,SameStrand,PosCorr,NegCorr");
			}
			
			if (checkAllCoexp) {
				pcA = new PrintStream(new FileOutputStream(outFileDir
						+ "/" + stem + ".output.genecoexpressionAll.txt"));
				
				pcA.println("ParID,GeneId,ParName,Chrom,Nuc,GeneSymbol,Distance,GeneContainsPar,SameStrand,CoexpGeneId,PosCorr,NegCorr");
			}
			
			if (checkAllPairsCoexp) {
				pap = new PrintStream(new FileOutputStream(outFileDir
						+ "/" + stem + ".output.genecoexpressionAllPairs.txt"));
				
				pap.println("ParID,GeneId,ParName,Chrom,Nuc,GeneSymbol,Distance,GeneContainsPar,SameStrand,ExperimentId,NumParMapProbe,NumGenesMapProbe,VectorSize,ParProbe,GeneProbe,Pearson");
			}
			
			
		} catch (Exception e) {
			System.err.println("Error writing to file");
			System.exit(0);
		}



		

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

		
		
		// Iterate through the gene list in batches
		
		Iterator recordItr = records.iterator();
		int batch = 0;

		String[] record;
		while (recordItr.hasNext()) {

			// A map that displays the original input PAR file information for a PAR
			// eg. parFileEntries(1124999) = "1124999,R49730.par.8.239068.239461,8,239068, ... false,true"
			parFileEntries = new HashMap<Long, String>(batchSize);

			HashMap<Gene, Gene> parToGene = new HashMap<Gene, Gene>();

			Collection<Gene> pars = new ArrayList<Gene>();
			Collection<Gene> genes = new ArrayList<Gene>();
			int count = batchSize;

			// work with a small batch
			while (recordItr.hasNext() && 0 < count) {

				record = (String[]) recordItr.next();

				// accessing data elements
				long ParID = Long.decode(record[ParIDIdx]);
				long GeneId = Long.decode(record[GeneIdIdx]);

				String ParName = record[ParNameIdx];
				String Chrom = record[ChromIdx];
				String Nuc = record[NucIdx];
				String GeneSymbol = record[GeneSymbolIdx];
				String Distance = record[DistanceIdx];
				String GeneContainsPar = record[GeneContainsParIdx];
				String SameStrand = record[SameStrandIdx];

				String allEntries = ParID + "," + GeneId + "," + ParName + ","
						+ Chrom + "," + Nuc + "," + GeneSymbol + "," + Distance
						+ "," + GeneContainsPar + "," + SameStrand;


				parFileEntries.put(ParID, allEntries);

				Gene par = this.parService.load(ParID);
				Gene g = this.parService.load(GeneId);

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

			// print batch numbers to output file
			/*
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
			*/
			
			
			batch++;

			

						
			
			/*
			 * Do experiments, depends on which switch is selected
			 */
			
			// PAR rankings at probe level
			if (checkPARprobe) {
				ppr.flush();
				outputRankingsProbelevel(eeCol, pars, ppr);
			}

			// PAR rankings
			if (checkPAR) {
				pxx.flush();
				pxe.flush();
				pex.flush();
				pee.flush();
				outputRankings(eeCol, pars, pxx, pxe, "max");
				outputRankings(eeCol, pars, pex, pee, "mean");
			}
			
			// Gene rankings
			if (checkGeneRank) {
				pxx.flush();
				pxe.flush();
				pex.flush();
				pee.flush();
				outputRankings(eeCol, genes, pxx, pxe, "max");
				outputRankings(eeCol, genes, pex, pee, "mean");
			}

			// PAR and Gene rankings together
			if (checkGeneCorank) {
				pgp_xx.flush();
				pgp_xe.flush();
				pgp_ex.flush();
				pgp_ee.flush();
				outputPARGeneCorank(eeCol, pars, genes, pgp_xx, pgp_xe, "max");
				outputPARGeneCorank(eeCol, pars, genes, pgp_ex, pgp_ee, "mean");
			}

			// PAR-gene co-expression Links found in database
			if (checkAllCoexpInDB) {
				//if (parToCoExps != null) {
				if (pargeneexpFile == null) {
					System.out.println("Error, need to supply pargeneexp file");
					System.exit(0);
				}
				//check all coexpressions given a par or look at each par-gene combination
				pco.flush();
				outputPARGeneCoexpressionLinks(pars, genes, pco);
			}
			
			
			// All co-expression between PAR and all genes with Links found in database
			if (checkAllCoexp) {
				pcA.flush();
				outputAllGeneCoexpressionLinks(pars, eeCol, 1, pcA);
			}
			
			// Calcualate all PAR/Gene pairs for all experiments where they co-exist
			if (checkAllPairsCoexp) {
				if (pargeneexpFile == null) {
					System.out.println("Error, need to supply pargeneexp file");
					System.exit(0);
				}
				pap.flush();
				outputAllPARGeneCoexpressionCalculated(pars, genes, pap);
			}

		}
		
		
		
		if (pxx != null) pxx.close();
		if (pxe != null) pxe.close();
		if (pex != null) pex.close();
		if (pee != null) pee.close();

		if (ppr != null) ppr.close();
		
		
		if (pgp_xx != null) pgp_xx.close();
		if (pgp_xe != null) pgp_xe.close();
		if (pgp_ex != null) pgp_ex.close();
		if (pgp_ee != null) pgp_ee.close();

		if (pco != null) pco.close();
		if (pcA != null) pcA.close();
		if (pap != null) pap.close();

		

		return null;
	}

	
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		PARMapperAnalyzeCLI p = new PARMapperAnalyzeCLI();
		Exception e = p.doWork(args);
		if (e != null) {
			throw new RuntimeException(e);
		}
	}

	/*
	 * This method finds the array types for the expressions given It might be
	 * important to know which array designs make use of double stranded probes
	 * 
	 * This is an example of what is printed to standard out
	 * 
	 * #399: Molecular profiles (HG-U95A) of dystrophin-deficient and normal human muscle
	 * #       8: Affymetrix GeneChip Human Genome U95 Version [1 or 2] Set HG-U95A
	 * #               AFFY_COLLAPSED
	 * 399     8       AFFY_COLLAPSED
	 * #403: Molecular heterogeneity in acute renal allograft rejection
	 * #       197: LC-17
	 * #               EST
	 * 403     197     EST
	 * #406: Prostaglandin J2 alters gene expression patterns and 26S proteasome...
	 * #       162: CAG Human 19k array v1.0
	 * #               mRNA
	 * 406     162     mRNA
	 * #409: Molecular portraits of human breast tumors
	 * #       40: SVC
	 * #               EST
	 * 409     40      EST
	 * #437: Serum stimulation of fibroblasts
	 * #       225: Stanford Human 10k prints 1-3
	 * #               DNA
	 * 437     225     DNA
	 * 
	 * To determine which are double stranded, a separate process will have to
	 * be used to extract the names and the types.
	 * 
	 */
	private void printExperimentTypes(Collection<ExpressionExperiment> eeCol) {

		// Printing out array information
		for (ExpressionExperiment e : eeCol) {
			System.out.println("#" + e.getId() + ": " + e.getName());
			Collection<ArrayDesign> ad = expressionExperimentService
					.getArrayDesignsUsed(e);
			for (ArrayDesign a : ad) {
				System.out.println("#\t" + a.getId() + ": " + a.getName());

				SequenceType st = getSequenceType(a);

				System.out.println("#\t\t" + st);
				System.out.println(e.getId() + "\t" + a.getId() + "\t" + st);
			}

		}
	}
	
	// called by printExperimentTypes to get the array type
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

		while (pItr.hasNext()) {
			Gene par = pItr.next();
			Gene gene = gItr.next();

			long[] expIds = parToCoExps.get(par.getId());

			// no information - skip
			if (expIds == null) {
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


			//System.out.println("Linkanalysis " + par.getId() + "\t" + gene.getId());
			
			// search Gemma for links
			Map<Gene, CoexpressionCollectionValueObject> coexp = pca.linkAnalysis(pargene, eexps, 1, false, true, 0);
			

			for (Gene g : coexp.keySet()) {
				CoexpressionCollectionValueObject ccvo = coexp.get(g);
				Collection<CoexpressionValueObject> cvos = ccvo.getAllGeneCoexpressionData(0);
				
				//System.out.println("GeneID: " + g.getId()+" Size: " + cvos.size());

				
				for (CoexpressionValueObject cvo : cvos) {
					//System.out.println("\t\tTS: " + cvo);
					
					//System.out.println(parFileEntries.get(par.getId())
					//		+ "," + g.getId()
					//		+ "," + cvo.getPositiveScore()
					//		+ "," + cvo.getNegativeScore());
					
					
					pco.println(parFileEntries.get(par.getId())
							+ "," + cvo.getPositiveScore()
							+ "," + cvo.getNegativeScore());
					
					
				}

			}

		}
	}
	
	
	
	// Retrieves all coexpression data for a par and all genes that are in the same
	// experiment above a coexpression threshold
	//
	// This method is different from the previous method outgenecoexpressions examining par-gene pairs.
	// Here a given par is compared with all genes containing a link in the database (implying that they
	// already have some coexpression)
	private void outputAllGeneCoexpressionLinks(
			Collection<Gene> pars, Collection<ExpressionExperiment> ccCol, int stringency, PrintStream pco) {
		
		Iterator<Gene> pItr = pars.iterator();

		// for (Gene par: pars) {
		while (pItr.hasNext()) {
			Gene par = pItr.next();


			Collection<BioAssaySet> exps = new ArrayList<BioAssaySet>();
			
			for (ExpressionExperiment ee : ccCol) {
				exps.add(ee);
			}


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


			Map<Gene, CoexpressionCollectionValueObject> coexp;
			
			try {
				ProbeLinkCoexpressionAnalyzer pca = new ProbeLinkCoexpressionAnalyzer();
				pca.setGeneService(parService);
				pca.setProbe2ProbeCoexpressionService(probe2ProbeCoexpressionService);
				coexp = pca.linkAnalysis(pargene, exps, stringency, false, false, 0);
			} catch (java.lang.IndexOutOfBoundsException e) {
				e.printStackTrace();
				System.out.println("Error performing linkAnalysis "
						+ e.getMessage());
				continue;
			}
			
			
			if (coexp == null || coexp.keySet().size() < 1) {
				System.out.println("No coexpression :(");
				continue;
			}
			
			CoexpressionCollectionValueObject ccvo = coexp.get(par);
			Collection<CoexpressionValueObject> cvos = ccvo.getAllGeneCoexpressionData(0);
			
			int zeroCount = 0;
			String output_half = (String) parFileEntries.get(par.getId());
			
			int csize = cvos.size();
			if (csize ==  0) {
				zeroCount++;
				continue;
			}
			
			
			for (CoexpressionValueObject cvo : cvos) {
				String output = output_half
							+ "," + cvo.getGeneId()
							+ "," + csize
							+ "," + cvo.getPositiveScore()
							+ "," + cvo.getNegativeScore();
				System.out.println(output);
				pco.println(output);
				
				
			}
			
			// print out the number of zeros
			System.out.println("Zerocount: " + par.getId()+"\t"+par.getId()+"\t"+"0:"+zeroCount);
			
		}
	}

	
	
	
	/* Looks at all par/gene pairs to see if there are any experiemnts with both
	 * and extracts the rank information for both
	 * 
	 * To get the PAR-gene-experiment map, parse the output to standard out
	 * This outputs A LOT of information to standard out!!
	 * Note: the stdout can be parsed by pargeneExplist.pl which sorts the
	 * par/gene/exp combo into lists.
	 */
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

		// need to map experiments to genes to calculate max/mean across experiments
		Map<Gene, Collection<Double>> allGeneRankings = new HashMap<Gene, Collection<Double>>();
		Map<Gene, Collection<Double>> allParRankings = new HashMap<Gene, Collection<Double>>();

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

		while (gItr.hasNext()) {

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

		}

	}

	
	/*
	 * outputs the ranks for PARs at the probe level.
	 * 
	 * Note, this method calls getRanksProbes, which is not in the Gemma
	 * CVS code base.
	 * 
	 */
	private void outputRankingsProbelevel(Collection<ExpressionExperiment> eeCol,
			Collection<Gene> pars, PrintStream p) {


		// use getRanksProbes function to map genes to experiments
		Map<ExpressionExperiment, Map<Gene, Map<DesignElement, Double[]>>> expressRankings = processedExpressionDataVectorService
				.getRanksProbes(eeCol, pars);

		// need to map experiments to genes to calculate max/mean across experiments
		Map<Gene, Map<DesignElement, Collection<Double[]>>> allProbeRankings = 
			new HashMap<Gene, Map<DesignElement, Collection<Double[]>>>();

		
		// go through the information that getRanksProbes returned
		for (ExpressionExperiment ee : expressRankings.keySet()) {
			Map<Gene, Map<DesignElement, Double[]>> expressRankingsGene = expressRankings.get(ee);
			for (Gene g : expressRankingsGene.keySet()) {
				Map<DesignElement, Double[]> expressRankingsGeneCompseq = expressRankingsGene.get(g);

				Collection<DesignElement> cde;

				// Add filter here to obtain only unique probes - if user requests it
				if (checkUniqueProbeMappings) {
					Collection<CompositeSequence> csCasted = new ArrayList<CompositeSequence>();

					// cast all composite sequences as design elements
					for (DesignElement de : expressRankingsGeneCompseq.keySet()) {
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
							// CS is unique
							cde.add(cs);
						} else {
							// CS non-unique
						}
					}

				}

				// don't filter for only unique probes
				else {
					cde = expressRankingsGeneCompseq.keySet();
				}

				
				for (DesignElement de : cde) {
					Double[] d = expressRankingsGeneCompseq.get(de);

					// subtract this entry from the list of probes, and if
					// necessary, the entire experiment
					if (d == null) {
						// Null rank value - just skip
					} else {

						if (!allProbeRankings.containsKey(g)) {
							allProbeRankings.put(
											g,
											new HashMap<DesignElement, Collection<Double[]>>());

						}

						if (!allProbeRankings.get(g).containsKey(de)) {
							allProbeRankings.get(g).put(de,
									new ArrayList<Double[]>());
						}

						allProbeRankings.get(g).get(de).add(d);

					}
				}


			}
		}

		// go through the data, calculate the rank information and print out to file
		for (Gene g : allProbeRankings.keySet()) {

			Map<DesignElement, Collection<Double[]>> allProbeRankingsGenes = allProbeRankings
					.get(g);

			for (DesignElement de : allProbeRankingsGenes.keySet()) {

				// element 0 is the mean, element 1 is the max ranking
				// from list of probe expressios, do stats
				int numOfExperiments = allProbeRankingsGenes.get(de).size();

				double valExpMean_RankMean = 0;
				double valExpMean_RankMax = 0;
				double valExpMax_RankMean = 0;
				double valExpMax_RankMax = 0;

				// calculate the mean rannkings
				for (Double[] entry : allProbeRankingsGenes.get(de)) {

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

			}

		}

	}

	/*
	 * the original method - get the rankings for PARs/Genes as a whole,
	 * not distinguishing probes
	 */
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

				int numNullRanks = 0; // keep track of how many null rankings there are
				boolean allNull = true;

				for (Double d : expressRankings.get(ee).get(g)) {

					// subtract this entry from the list of probes, and if
					// necessary, the entire experiment
					if (d == null) {
						numNullRanks++;

					} else {

						// calculate the maximum and average rankings
						maxAndAveRank[1] += d;
						if (maxAndAveRank[0] < d)
							maxAndAveRank[0] = d;

						allNull = false;
					}

				}

				// check if genes for this experiment are null-ranked. if so, skip
				if (allNull) {
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


		}

	}

	/*
	 * Return a map for column name to index number
	 */
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

		String header;

		try {
			in = new BufferedReader(new FileReader(inFile));
			String line;
			header = in.readLine();
			fHash = getIndices(header);

			while ((line = in.readLine()) != null) {
				if (line.startsWith("#"))
					continue;
				String[] s = line.trim().split("\t");
				fRecords.add(s);
			}

			in.close();

		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.out.println("File " + inFile + " not found - "
					+ e.getMessage());
			System.exit(0);
		} catch (IOException e) {
			e.printStackTrace();
			System.out.println("File " + inFile + " reads a bit wonky - "
					+ e.getMessage());
			System.exit(0);
		}

		headerLookup = fHash;
		records = fRecords;
	}

	/*
	 * Reads a par/gene/experiments file and saves the results to hash
	 * this file maps each PAR to a gene and a list of experiments that has both
	 * 
	 * The format should by the like the following
	 * ParID	GeneID	Experiments
	 * 1126681	12860	441,442,443,445,519,521
	 * 1130907	166867	443
	 * 1134105	7379	442
	 * 1134841	14121	442
	 * 1135031	119682	441,442,519,521
	 */
	// parToCoExps
	private static void readpargeneFile(String inFile) {

		parToCoExps = new HashMap<Long, long[]>();

		BufferedReader in;

		String header;

		try {
			in = new BufferedReader(new FileReader(inFile));
			String line;
			header = in.readLine();

			while ((line = in.readLine()) != null) {
				if (line.startsWith("#"))
					continue;
				String[] s = line.trim().split("\t");
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
			e.printStackTrace();
			System.out.println("File " + inFile + " not found - "
					+ e.getMessage());
			System.exit(0);
		} catch (IOException e) {
			e.printStackTrace();
			System.out.println("File " + inFile + " reads a bit wonky - "
					+ e.getMessage());
			System.exit(0);
		}

	}
	
	/*
	 * UNSTABLE
	 * 
	 * This method outputs the pairwise correlations between pars and genes
	 * The correlations are calculated here using Pearson.
	 * 
	 * This loads all experiments at once, therefore it calls the database
	 * less but this method crashes very frequently!!!!!
	 * 
	 * Note: this method requires that a PAR/gene/experiment map file be used
	 * Note: this method takes a LONG time to finish (order of 10+ days) and
	 * it might run out of memory.  May need to split into other jobs or resume
	 * from previous runs.
	 */
	private void outputAllPARGeneCoexpressionCalculated_unstable(Collection<Gene> pars, Collection<Gene> genes, PrintStream pap) {
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
			
			Collection<DoubleVectorValueObject> dvvos = processedExpressionDataVectorService.getProcessedDataArrays(exps, pargene);
			
			// for every probe, assign it to either par or gene for the experiment
			for (DoubleVectorValueObject dvvo : dvvos) {
				ExpressionExperiment ee = dvvo.getExpressionExperiment();
				expressionExperimentService.thawLite(ee);
				Long eeId = new Long(ee.getId());
				ee = null;
				
				Collection<Gene> genelist = dvvo.getGenes();
				
				// go through the gene list to determine if PAR or gene and put into lists
				// the data returned by getProcessedDataArrays does not tell us this directly
				boolean isParGene = false;
				for (Gene g: genelist) {
					long id = g.getId();
					if (id == gene.getId()) {
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
			
			// go through list of PAR/gene 
			for (Long eeId: expParMap.keySet()) {
				
				
				
				if (! expGeneMap.containsKey(eeId)) {
					System.out.println("No experiment data for gene "+gene.getId()+" in this experiment: "+eeId);
					continue;
				}
				
				//iterate through the par probes
				for (DoubleVectorValueObject parDvvo : expParMap.get(eeId)) {
					DesignElement pd = parDvvo.getDesignElement();
					
					double[] pdata = parDvvo.getData();
					int parProbeMappings  = parDvvo.getGenes().size();
					
					
					//iterate through the gene probes
					for (DoubleVectorValueObject geneDvvo : expGeneMap.get(eeId)) {
						DesignElement gd = geneDvvo.getDesignElement();
						
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
						
						// check data integrity - tag if missing value in column
						for (int i=0; i<pdata.length; i++) {
							if ((new Double(pdata[i])).isNaN() || (new Double(gdata[i])).isNaN()) {
								badIndices[i] = true;
								numbadIndices++;
							}
						}
						
						// make sure there are at least 3 samples in vector
						if (pdata.length - numbadIndices < 3) {
							System.out.println("Not enough samples ("
									+ (pdata.length - numbadIndices)
									+") to calc coexpression: ParProbeID:"
									+pd.getId()+"\tGeneProbeID:"+gd.getId());
						}
						
						// remove bad vector entries from both if either is labeled as bad
						if (0 < numbadIndices) {
							double[] newpdata = new double[pdata.length - numbadIndices];
							double[] newgdata = new double[pdata.length - numbadIndices];
							
							int j=0;
							for (int i=0; i<pdata.length; i++) {
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
					
				}
			}
		}
		
		
		

	}
	
	
	/*
	 * STABLE
	 * 
	 * This method outputs the pairwise correlations between pars and genes
	 * The correlations are calculated here using Pearson.
	 * 
	 * This loads all the data one experiment at a time.  It is slower than the
	 * other method but it is much more stable.
	 * 
	 * Note: this method requires that a PAR/gene/experiment map file be used
	 * Note: this method takes a LONG time to finish (order of 10+ days) and
	 * it might run out of memory.  May need to split into other jobs or resume
	 * from previous runs.
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
			
			// for each experiment, call getProcessedDataArrays
			for (ExpressionExperiment ee: exps) {
				Collection<ExpressionExperiment> ees = new ArrayList<ExpressionExperiment>();
				ees.add(ee);
				
				
				Collection<DoubleVectorValueObject> dvvos = processedExpressionDataVectorService.getProcessedDataArrays(ees, pargene);
				
				Map<DesignElement, double[]> parData  = new HashMap<DesignElement, double[]>();
				Map<DesignElement, double[]> geneData = new HashMap<DesignElement, double[]>();
				
				Map<DesignElement, Integer> probemappingCount = new HashMap<DesignElement, Integer>();
				
				// separate data into PAR or gene - the returned object does not
				// differenciate which probes belong to PAR or gene
				for (DoubleVectorValueObject dvvo : dvvos) {
					Collection<Gene> genelist = dvvo.getGenes();
					
					boolean isParGene = false;
					for (Gene g: genelist) {
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
						for (int i=0; i<pdata.length; i++) {
							if ((new Double(pdata[i])).isNaN() || (new Double(gdata[i])).isNaN()) {
								badIndices[i] = true;
								numbadIndices++;
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
								if (badIndices[i]) {
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
					
					
				}
			}
		}
		
		
		

	}
	
	/*
	 * Loads experiments from a list on a file
	 */
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
			e.printStackTrace();
			System.out.println("File " + ExperimentListFile + " not found - "
					+ e.getMessage());
			System.exit(0);
		} catch (IOException e) {
			e.printStackTrace();
			System.out.println("File " + ExperimentListFile + " reads a bit wonky - "
					+ e.getMessage());
			System.exit(0);
		}
		
		return eeCol;
	
	}

}
