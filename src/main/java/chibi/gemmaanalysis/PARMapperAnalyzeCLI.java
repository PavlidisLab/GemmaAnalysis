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

import ubic.gemma.model.expression.bioAssayData.ProcessedExpressionDataVectorService;
import ubic.gemma.model.expression.bioAssayData.ProcessedExpressionDataVectorDao.RankMethod;
import ubic.gemma.model.expression.designElement.CompositeSequenceService;
import ubic.gemma.model.expression.designElement.DesignElement;
import ubic.gemma.model.expression.experiment.ExpressionExperiment;
import ubic.gemma.model.expression.experiment.ExpressionExperimentService;
import ubic.gemma.model.genome.Gene;
import ubic.gemma.model.genome.PhysicalLocationService;
import ubic.gemma.model.genome.Taxon;
import ubic.gemma.model.genome.TaxonService;
import ubic.gemma.model.genome.gene.GeneService;
import ubic.gemma.util.AbstractSpringAwareCLI;

public class PARMapperAnalyzeCLI extends AbstractSpringAwareCLI {
	
	
	private Taxon taxon;
	private String inFile = "/home/hmokada/scratch/human-par-gene-relations.subset.txt";
	private String outFileDir = "/home/hmokada/scratch/outputfiles";
	
	private boolean checkPAR = false;
	private boolean checkPARprobe = false;
	private boolean checkGeneCorank = false;
	private boolean checkGeneRank = false;
	
	private static Map headerLookup;
	private static String[] headers;
	private static Collection<String[]> records;
	private static Map parFileEntries;

	GeneService parService;
	TaxonService taxonService;
	//ExpressionExperimentSetService expressionExperimentSetService;
	ExpressionExperimentService expressionExperimentService;
	ProcessedExpressionDataVectorService processedExpressionDataVectorService;
	//CompositeSequenceService compositeSequenceService;
	
	@Override
	protected void buildOptions() {
		// TODO Auto-generated method stub

	}
	
	@Override
	protected void processOptions() {
		super.processOptions();

		this.taxonService = (TaxonService) this.getBean("taxonService");
		this.parService = (GeneService) this.getBean("geneService");

		/*
		 * Add processing o
		 */
		
		if (hasOption("inFile")) {
			inFile = getOptionValue("inFile");
			
			//if ( == null) {
			//	log.error("ERROR: Cannot find file " + inFile);
			//}
		}
		
		if (hasOption("outFileDir")) {
			outFileDir = getOptionValue("outFileDir");
			
			//if ( == null) {
			//	log.error("ERROR: Cannot find directory " + outFileDir);
			//}
		}
		
		this.checkPAR 		= this.hasOption("checkPAR");
		this.checkPARprobe	= this.hasOption("checkPARprobe");
		this.checkGeneCorank= this.hasOption("checkGeneCorank");
		this.checkGeneRank	= this.hasOption("checkGeneRank");
		
	}

	@Override
	protected Exception doWork(String[] args) {
		// TODO Auto-generated method stub
		
		Exception exc = processCommandLine( "test", args );
		
		this.parService = ( GeneService ) this.getBean( "geneService" );
		this.taxonService = ( TaxonService ) this.getBean( "taxonService" );
		//this.expressionExperimentSetService = ( ExpressionExperimentSetService ) this.getBean( "expressionExperimentSetService" );
		this.expressionExperimentService = ( ExpressionExperimentService ) this.getBean( "expressionExperimentService" );
		this.processedExpressionDataVectorService = ( ProcessedExpressionDataVectorService ) this.getBean( "processedExpressionDataVectorService" );
		//this.compositeSequenceService = ( CompositeSequenceService ) this.getBean( "compositeSequenceService" );
		
		
		
//CompositeSequenceService.getGenes or getGenesWithSpecificity.

		
		//log.info(taxonService.loadAll());
				
		// For my mac
//		String inFile =  "/Users/mokada/development/PARs/data/human-par-gene-relations.txt";
//		String outFileDir = "/Users/mokada/development/PARs/data/outputfiles";
		
		int batchSize = 50;
		
		
		PrintStream pxx = null; // declare a print stream object
		PrintStream pxe = null;
		PrintStream pex = null;
		PrintStream pee = null;
		
		PrintStream ppr = null;
		
		
		PrintStream pgp_xx = null; // for gene vs. PAR rankings (within same experiments)
		PrintStream pgp_xe = null;
		PrintStream pgp_ex = null;
		PrintStream pgp_ee = null;
		
		
		try
		{
			// Connect print stream to the output stream
			pxx = new PrintStream(new FileOutputStream(outFileDir + "/human-par-gene-relations.output.xx.txt"));
			pxe = new PrintStream(new FileOutputStream(outFileDir + "/human-par-gene-relations.output.xe.txt"));
			pex = new PrintStream(new FileOutputStream(outFileDir + "/human-par-gene-relations.output.ex.txt"));
			pee = new PrintStream(new FileOutputStream(outFileDir + "/human-par-gene-relations.output.ee.txt"));
			
			ppr = new PrintStream(new FileOutputStream(outFileDir + "/human-par-gene-relations.output.probelevel.txt"));
			
			pgp_xx = new PrintStream(new FileOutputStream(outFileDir + "/human-par-gene-relations.output.geneVsPar.xx.txt"));
			pgp_xe = new PrintStream(new FileOutputStream(outFileDir + "/human-par-gene-relations.output.geneVsPar.xe.txt"));
			pgp_ex = new PrintStream(new FileOutputStream(outFileDir + "/human-par-gene-relations.output.geneVsPar.ex.txt"));
			pgp_ee = new PrintStream(new FileOutputStream(outFileDir + "/human-par-gene-relations.output.geneVsPar.ee.txt"));
		} catch (Exception e) {
			System.err.println("Error writing to file");
			System.exit(0);
		}
		
		// allow a user to enter a filename through the command line
		if (0 < args.length && args[0] != null) {
			inFile = args[0];
		}
		
		
		System.out.println("Reading file: " + inFile);
		readPARFile(inFile);
		
		
		// Establish the column numbers
		
		int ParIDIdx		= getIndex("ParID");
		int ParNameIdx		= getIndex("ParName");
		int ChromIdx		= getIndex("Chrom");
		int NucIdx			= getIndex("Nuc");
		int GeneIdIdx		= getIndex("GeneId");
		int GeneSymbolIdx	= getIndex("GeneSymbol");
		int DistanceIdx		= getIndex("Distance");
		int GeneContainsParIdx	= getIndex("GeneContainsPar");
		int SameStrandIdx	= getIndex("SameStrand");
		
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
		
		//System.out.println("\n\nExperiments by taxon\n\n");
		//Taxon taxon = taxonService.findByCommonName( "human" );
		
		Collection<ExpressionExperiment> eeCol = expressionExperimentService.findByTaxon(taxon);
		System.out.println(eeCol.size());
		
		
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
//System.out.println(record);
				
				// accessing data elements
				int ParID			= Integer.parseInt( record[ParIDIdx] );
				int GeneId			= Integer.parseInt( record[GeneIdIdx] );
				
				// it's not necessary to parse thees to int/boolean as it'll printed off to file anyway
				String ParName		= record[ ParNameIdx];
				String Chrom		= record[ ChromIdx];
				//int Nuc				= Integer.parseInt( record[NucIdx] );
				String Nuc		= record[ NucIdx];
				String GeneSymbol	= record[ GeneSymbolIdx ];
				//int Distance		= Integer.parseInt( record[DistanceIdx] );
				//boolean GeneContainsPar	= Boolean.parseBoolean( record[GeneContainsParIdx] );
				//boolean SameStrand		= Boolean.parseBoolean( record[SameStrandIdx] );
				String Distance		= record[DistanceIdx];
				String GeneContainsPar	= record[GeneContainsParIdx];
				String SameStrand		= record[SameStrandIdx];
				
				String allEntries = ParID + "," 
				+ GeneId + "," 
				+ ParName + "," 
				+ Chrom + "," 
				+ Nuc + "," 
				+ GeneSymbol + "," 
				+ Distance + "," 
				+ GeneContainsPar + "," 
				+ SameStrand;
				
				
				//HashMap<Integer, String> parFileEntries = new HashMap<Integer, String>(batchSize);
				
				
				parFileEntries.put(Long.valueOf(ParID), allEntries); //.put(ParID, record);
				
				
				Gene par = this.parService.load(ParID);
				Gene g = this.parService.load(GeneId);
				
				
				// TODO: place as log
				if (par == null || g == null) {
					System.out.println("PAR or Gene doesn't exist: "+ par.getId() +"\t"+ g.getId());
					continue;
				}
				
				pars.add(par);
				genes.add(g);
				
				parToGene.put(par, g);
				
				
				count--;  // count down the number of genes in a batch
				
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
			
			
//ParID,GeneId,ParName,Chrom,Nuc,GeneSymbol,Distance,GeneContainsPar,SameStrand,ParNumExperiments,ParNumSamples,ParRank,GeneNumExperiments,GeneNumSamples,GeneRank
			
			
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
			
			if (checkPARprobe) {
				outputAll_probelevel(eeCol, pars, ppr);
			}
			
			if (checkPAR) {
				outputAll(eeCol, pars, pxx, pxe, "max");	
				outputAll(eeCol, pars, pex, pee, "mean");
			}
			
			if (checkGeneCorank) {
				outputAll_geneCorank(eeCol, pars, genes, pgp_xx, pgp_xe, "max");
				outputAll_geneCorank(eeCol, pars, genes, pgp_ex, pgp_ee, "mean");
				
				// not working
				//outputAll_geneCorank(eeCol, pars, genes, pgp_xx);
				
				//outputAll_geneCorank(eeCol, parToGene, pxx, pxe, "max");
				//outputAll_geneCorank(eeCol, parToGene, pex, pee, "mean");
			}
			
			if (checkGeneRank) {
				outputAll(eeCol, genes, pxx, pxe, "max");
				outputAll(eeCol, genes, pex, pee, "mean");
			}
			

			
			
			
		}
		pxx.close();
		pxe.close();
		pex.close();
		pee.close();
		
		ppr.close();

		
		return null;
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		PARMapperAnalyzeCLI p = new PARMapperAnalyzeCLI();
        Exception e = p.doWork( args );
        if ( e != null ) {
            throw new RuntimeException( e );
        }
	}
	
	
	
	private void outputAll_geneCorank(Collection<ExpressionExperiment> eeCol, Collection<Gene> pars, Collection<Gene> genes, PrintStream px, PrintStream pe, String rankMethodStr) {
		
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
		Map<ExpressionExperiment, Map<Gene, Collection<Double>>> expressRankings = processedExpressionDataVectorService.getRanks(eeCol, allProbesGenes, method);
		//Map<ExpressionExperiment, Map<Gene, Collection<Double>>> expressRankings = processedExpressionDataVectorService.getRanks(eeCol, pars, method);
		
		
//for (Gene g : genes) {System.out.println("Gene "+g.getId());}
//for (Gene p : genes) {System.out.println("Par  "+p.getId());}
//for (Gene p : allProbesGenes) {System.out.println("Allp "+p.getId());}
//for (ExpressionExperiment e : expressRankings.keySet()) {System.out.println("Experiment "+e.getId());}
//if (true) return;
		
		// need to map experiments to genes to calculate max/mean across experiments
		Map<Gene, Collection<Double>> allGeneRankings = new HashMap<Gene, Collection<Double>>();
		Map<Gene, Collection<Double>> allParRankings  = new HashMap<Gene, Collection<Double>>();
		
		Map<Gene, int[]> geneCounts = new HashMap<Gene, int[]>();
		
		Map<Gene, Collection<ExpressionExperiment>> numOfExpsPars  = new HashMap<Gene, Collection<ExpressionExperiment>>();
		Map<Gene, Collection<ExpressionExperiment>> numOfExpsGenes = new HashMap<Gene, Collection<ExpressionExperiment>>();
		
		for ( ExpressionExperiment ee : expressRankings.keySet() ) {
System.out.println("Experiment: "+ee.getId());
			Map<Gene, Collection<Double>> expressRankingsExperiment = expressRankings.get(ee);
			
			Iterator<Gene> gItr = genes.iterator();
			Iterator<Gene> pItr = pars.iterator();
			
			Gene g, p;
			
			
			//for ( Gene g : expressRankings.get(ee).keySet() ) {
			while (pItr.hasNext()) {
				g = gItr.next();
				p = pItr.next();
System.out.print("Checking: "+g.getId()+"\t"+p.getId()+"\t");
if (! expressRankingsExperiment.containsKey(p)) System.out.print("NoPar\t"); 
if (! expressRankingsExperiment.containsKey(g)) System.out.print("NoGene\t"); 
//System.out.print("Checking: "+g.getId()+"\t"+p.getId()+"\t"+expressRankingsExperiment.get(g).size()+"\t"+expressRankingsExperiment.get(p).size());
				
				// skip if both gene and PAR not in experiment
				if ((! expressRankingsExperiment.containsKey(g)) || (! expressRankingsExperiment.containsKey(p))) {
System.out.println("skipping, no info");
					continue;
				}
				
				
				Collection<Double> candidateGeneRanks = new ArrayList<Double>();
				Collection<Double> candidatePARRanks  = new ArrayList<Double>();
				for ( Double d : expressRankingsExperiment.get(g) ) {
					// subtract this entry from the list of probes, and if necessary, the entire experiment
					if (d == null)
						continue;
					
					candidateGeneRanks.add(d);
				}

				for ( Double d : expressRankingsExperiment.get(p) ) {
					// subtract this entry from the list of probes, and if necessary, the entire experiment
					if (d == null)
						continue;
					
					candidatePARRanks.add(d);
				}
				
				
				//skip if both gene and PAR not in experiment
				if (candidateGeneRanks.isEmpty() || candidatePARRanks.isEmpty()) {
System.out.println("skipping again, at least one is empty");
					continue;
				}
				
				
				// note that PAR is key for both
				allGeneRankings.get(p).addAll(candidateGeneRanks);
				allParRankings.get( p).addAll(candidatePARRanks);
				
				// note how many experiments are performed per PAR and gene
				if (! numOfExpsPars.containsKey(g)) {
					numOfExpsPars.put( g, new ArrayList<ExpressionExperiment>());
					numOfExpsGenes.put(g, new ArrayList<ExpressionExperiment>());
				}
				numOfExpsPars.get( g).add(ee);
				numOfExpsGenes.get(g).add(ee);

				
System.out.println("----Exists!!!----");

			}
		}
		
		
		
		////////////////////////////////////////////////////////////
		// WARNING: There are no experiments with both PAR and nearest
		// gene, so this part of the code is untested!!!!!
		////////////////////////////////////////////////////////////
				
		
		// this time, combine the entries for PARs and genes
		
		Iterator<Gene> gItr = genes.iterator();
		Iterator<Gene> pItr = pars.iterator();
		
		Gene g, p;
		
		
		//for ( Gene g : expressRankings.get(ee).keySet() ) {
		while (gItr.hasNext()) {
		//for ( Gene p : allGeneRankings.keySet() ) {
			
			g = gItr.next();
			p = pItr.next();
			
			// weed out all those that don't have experiments
			if (! allParRankings.containsKey(p)) {
System.out.println("No Gene information for "+p.getId());
				continue;
			}
			
			
			
			double maxrankPar = 0;
			double meanrankPar = 0;
			
			double maxrankGene = 0;
			double meanrankGene = 0;
			
			int numSamplesPar  = allParRankings.get( p).size();
			int numSamplesGene = allGeneRankings.get(p).size();
			
			int numExpsPar  = numOfExpsPars.get( p).size();
			int numExpsGene = numOfExpsGenes.get(g).size();
			
			for (Double d : allParRankings.get(p)) {
				double rank = d.doubleValue();
				if (maxrankPar < rank)
					maxrankPar = rank;
				meanrankPar += rank;;
			}
			meanrankPar = meanrankPar / numSamplesPar;
			
			for (Double d : allGeneRankings.get(p)) {
				double rank = d.doubleValue();
				if (maxrankGene < rank)
					maxrankGene = rank;
				meanrankGene += rank;;
			}
			meanrankGene = meanrankGene / numSamplesGene;
			
			
			
			px.println(parFileEntries.get(p.getId())+","+
					numExpsPar+","+
					numSamplesPar+","+
					maxrankPar+","+
					numExpsGene+","+
					numSamplesGene+","+
					maxrankGene);
			pe.println(parFileEntries.get(p.getId())+","+
					numExpsPar+","+
					numSamplesPar+","+
					meanrankPar+","+			
					numExpsGene+","+
					numSamplesGene+","+
					meanrankGene);
			
			//px.println(g.getId()+","+numExperiments+","+numSamples+","+maxrank);
			//pe.println(g.getId()+","+numExperiments+","+numSamples+","+meanrank);
			
		}
		
		
		
		
		
	}

	
	
	
	
	
	private void outputAll_geneCorank_old(Collection<ExpressionExperiment> eeCol, HashMap<Gene, Gene> parToGene, PrintStream px, PrintStream pe, String rankMethodStr) {
		
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
		Map<ExpressionExperiment, Map<Gene, Collection<Double>>> expressRankings = processedExpressionDataVectorService.getRanks(eeCol, pars, method);
		
		// need to map experiments to genes to calculate max/mean across experiments
		HashMap<Gene, double[]> allGeneRankings = new HashMap<Gene, double[]>();
		HashMap<Gene, int[]> geneCounts = new HashMap<Gene, int[]>();
		
		for ( ExpressionExperiment ee : expressRankings.keySet() ) {
			for ( Gene p : expressRankings.get(ee).keySet() ) {
				
				double[] maxAndAveRank = new double[2];
				maxAndAveRank[0] = 0;  // max
				maxAndAveRank[1] = 0;  // mean
								
				int numNullRanks = 0;  // keep track of how many null rankings there are 
				boolean allNull = true;
				
				for ( Double d : expressRankings.get(ee).get(p) ) {
					
					// subtract this entry from the list of probes, and if necessary, the entire experiment
					if (d == null) {
						numNullRanks++;
						System.out.println("Null rank value "+p.getId()+", numNullRanks: "+ numNullRanks);
						
					} else {
						
						// calculate the maximum and average rankings
						maxAndAveRank[1] += d;
						if (maxAndAveRank[0] < d) maxAndAveRank[0] = d;
						
						allNull = false;
					}
					
					
				}
				
				
				// check if genes for this experiment are null-ranked. if so, skip
				if (allNull) {
					System.out.println("Null: No exps for this gene "+p.getId()+", exp: "+ ee.getId());
					continue;
				}
				
				
				maxAndAveRank[1] = maxAndAveRank[1] / (expressRankings.get(ee).get(p).size() - numNullRanks);
				
				
				//update hash of expression levels
				if (allGeneRankings.containsKey(p)) {
					if (allGeneRankings.get(p)[1] < maxAndAveRank[0])
						allGeneRankings.get(p)[1] = maxAndAveRank[0];
					allGeneRankings.get(p)[1] += maxAndAveRank[1];
					
					geneCounts.get(p)[0]++;
					geneCounts.get(p)[1] += expressRankings.get(ee).get(p).size();
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
		
		
		for ( Gene g : allGeneRankings.keySet() ) {
			
			int numExperiments = geneCounts.get(g)[0];
			int numSamples = geneCounts.get(g)[1];
			
			double maxrank  = allGeneRankings.get(g)[0];
			double meanrank = allGeneRankings.get(g)[1] / numExperiments;
			
			px.println(parFileEntries.get(g.getId())+","+numExperiments+","+numSamples+","+maxrank);
			pe.println(parFileEntries.get(g.getId())+","+numExperiments+","+numSamples+","+meanrank);
		}
		
		
		
		
		
	}

	
	
	private void outputAll_geneCorank_probelevel(Collection<ExpressionExperiment> eeCol, Collection<Gene> pars, Collection<Gene> genes, PrintStream p) {
		
		Collection<Gene> allProbesGenes = new ArrayList<Gene>();
		allProbesGenes.addAll(pars);
		allProbesGenes.addAll(genes);
		
		
		/*
		for ( ExpressionExperiment ee : expressRankings.keySet() ) {
			Map<Gene, Collection<Double>> expressRankingsExperiment = expressRankings.get(ee);
			
			Iterator<Gene> gItr = genes.iterator();
			Iterator<Gene> pItr = pars.iterator();
			
			Gene g, p;
			
			
			//for ( Gene g : expressRankings.get(ee).keySet() ) {
			while (gItr.hasNext()) {
				g = gItr.next();
				p = pItr.next();
		*/

		
		// use getRanks function to map genes to experiments
		Map<ExpressionExperiment, Map<Gene, Map<DesignElement,Double[]>>> expressRankings = processedExpressionDataVectorService.getRanksProbes(eeCol, allProbesGenes);
//		Map<ExpressionExperiment, Map<Gene, Map<DesignElement,Double[]>>> expressRankings = new HashMap<ExpressionExperiment, Map<Gene, Map<DesignElement,Double[]>>>();
		
		// need to map experiments to genes to calculate max/mean across experiments
		Map<Gene, Double[]> allGeneRankings = new HashMap<Gene, Double[]>();
		Map<Gene, int[]> geneCounts = new HashMap<Gene, int[]>();
		
//		Map<Gene, Map<DesignElement, Collection<Double>>> allProbeRankingsMeanrank = new HashMap<Gene, Map<DesignElement,Collection<Double>>>();
//		Map<Gene, Map<DesignElement, Collection<Double>>> allProbeRankingsMaxrank = new HashMap<Gene, Map<DesignElement,Collection<Double>>>();
		
		// Contains all the data for Genes and Pars
		Map<Gene, Map<DesignElement, Collection<Double[]>>> allGeneProbeRankings = new HashMap<Gene, Map<DesignElement,Collection<Double[]>>>();
		Map<Gene, Map<DesignElement, Collection<Double[]>>> allPARProbeRankings  = new HashMap<Gene, Map<DesignElement,Collection<Double[]>>>();
		
		
		for ( ExpressionExperiment ee : expressRankings.keySet() ) {
			Map<Gene, Map<DesignElement,Double[]>> expressRankingsGene = expressRankings.get(ee);
			
			Iterator<Gene> gItr = genes.iterator();
			Iterator<Gene> pItr = pars.iterator();
			
			Gene gene, par;
			
			//for ( Gene g : expressRankingsGene.keySet() ) {
			while (gItr.hasNext()) {
				gene = gItr.next();
				par = pItr.next();
				
				if ((! expressRankingsGene.containsKey(gene)) || (! expressRankingsGene.containsKey(par)))
					continue;
				
				Map<DesignElement, Double[]> expressRankingsGeneCompseq = expressRankingsGene.get(gene);
				Map<DesignElement, Double[]> expressRankingsPARCompseq = expressRankingsGene.get(par);
				
				// take all the samples for that gene in this experiment
				for ( DesignElement de : expressRankingsGeneCompseq.keySet() ) {
					
					//Double d = expressRankings.get(ee).get(g).get(de);
					Double[] d = expressRankingsGeneCompseq.get(de);
					
					// subtract this entry from the list of probes, and if necessary, the entire experiment
					if (d == null) {
						//numNullRanks++;
						System.out.println("Null rank value "+gene.getId()+", with probeID: "+ de.getId());
						
					} else {
						
						if (! allGeneProbeRankings.containsKey(par)) {
							allGeneProbeRankings.put(par, new HashMap<DesignElement,Collection<Double[]>>());
							
							//allProbeRankings.put(g, new HashMap<DesignElement,Collection<Double>>());
							//allProbeRankingsMaxrank.put( g, new HashMap<DesignElement,Collection<Double>>());
						}
						
						if (! allGeneProbeRankings.get(par).containsKey(de)) {
							allGeneProbeRankings.get(par).put(de, new ArrayList<Double[]>());
							//allProbeRankingsMaxrank.get( g).put(de, new ArrayList<Double>());
						}
						
						allGeneProbeRankings.get(par).get(de).add(d);
						//allProbeRankingsMaxrank.get( g).get(de).add(d[1]);
						
						// calculate the maximum and average rankings
						//maxAndAveRank[1] += d;
						//if (maxAndAveRank[0] < d) maxAndAveRank[0] = d;
						
						//allNull = false;
					}
					
					
				}
				
				
				////////////////////////////////////// repeat for pars, doesn't make sense to do this at probe level!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				////////////////////////////////////////////////////
				
				
				// check if genes for this experiment are null-ranked. if so, skip
				//if (allNull) {
				//	System.out.println("Null: No exps for this gene "+g.getId()+", exp: "+ ee.getId());
				//	continue;
				//}
				
				
				//maxAndAveRank[1] = maxAndAveRank[1] / (expressRankings.get(ee).get(g).size() - numNullRanks);
				
				
				//update hash of expression levels
				//if (allGeneRankings.containsKey(g)) {
				//	if (allGeneRankings.get(g)[1] < maxAndAveRank[0])
				//		allGeneRankings.get(g)[1] = maxAndAveRank[0];
				//	allGeneRankings.get(g)[1] += maxAndAveRank[1];
				//	
				//	geneCounts.get(g)[0]++;
				//	geneCounts.get(g)[1] += expressRankings.get(ee).get(g).size();
				//} 
				// new gene, must create new entities in hash
				//else {
				//	allGeneRankings.put(g, maxAndAveRank);
				//	
				//	
				//	// number of experiments and number samples for this gene
				//	int[] counts = new int[2];
				//	counts[0] = 1;
				//	counts[1] = expressRankings.get(ee).get(g).size();
				//	
				//	geneCounts.put(g, counts);
				//}
				
			}
		}
		
		
		//allProbeRankings.get(g).put(de, new ArrayList<Double>());
		for ( Gene g : allGeneProbeRankings.keySet() ) {
			
			Map<DesignElement,Collection<Double[]>> allProbeRankingsGenes = allGeneProbeRankings.get(g);
			
			//int numExperiments = geneCounts.get(g)[0];
			//int numSamples = geneCounts.get(g)[1];
			//
			//double maxrank  = allGeneRankings.get(g)[0];
			//double meanrank = allGeneRankings.get(g)[1] / numExperiments;
			//
			////px.println(parFileEntries.get(g.getId())+","+numExperiments+","+numSamples+","+maxrank);
			////pe.println(parFileEntries.get(g.getId())+","+numExperiments+","+numSamples+","+meanrank);
			
			for (DesignElement de : allProbeRankingsGenes.keySet()) {
				
				
				// element 0 is the mean, element 1 is the max ranking

				// from list of probe expressios, do stats
				
				int numOfExperiments = allProbeRankingsGenes.get(de).size();
				
				double valExpMean_RankMean = 0;
				double valExpMean_RankMax = 0;
				double valExpMax_RankMean = 0;
				double valExpMax_RankMax = 0;
				
				
				for (Double[] entry : allProbeRankingsGenes.get(de)) {
					
//					if (entry[0] == null || entry[1] == null) continue;
					
					double valRankMean = entry[0].doubleValue();
					
					valExpMean_RankMean += valRankMean;
					if (valExpMax_RankMean < valRankMean)
						valExpMax_RankMean = valRankMean;
					
					double valRankMax  = entry[1].doubleValue();
					
					valExpMean_RankMax += valRankMax;
					if (valExpMax_RankMax < valRankMax) {
						valExpMax_RankMax = valRankMax;
					}
					
				}
				
				
				
				valExpMean_RankMean = valExpMean_RankMean / numOfExperiments;
				valExpMean_RankMax = valExpMean_RankMax / numOfExperiments;

				
				
				p.println(parFileEntries.get(g.getId())+","+de.getId()+","+numOfExperiments+","+valExpMean_RankMean+","+valExpMean_RankMax+","+valExpMax_RankMean+","+valExpMax_RankMax);

				
				//px.println(g.getId()+","+numExperiments+","+numSamples+","+maxrank);
				//pe.println(g.getId()+","+numExperiments+","+numSamples+","+meanrank);
			}
			
			
		}
		
		
		
		
		
	}
	
	
	
	
	
	private void outputAll_probelevel(Collection<ExpressionExperiment> eeCol, Collection<Gene> pars, PrintStream p) {
		
System.out.println("reading probes");
for (ExpressionExperiment e : eeCol) {System.out.println(e.getId());}
System.out.println("Genes");
for (Gene e : pars) {System.out.println(e.getId());}

//WHY THE F ISN'T THIS WORKING??  FOR SOME REASON THEY AREN'T RETURNING ANYTHING FROM THE GETRANKS METHOD!!!

		//RankMethod method;
		//
		//// todo: make these boolean
		//if (rankMethodStr.equalsIgnoreCase("max")) {
		//	method = RankMethod.max;
		//} else {
		//	method = RankMethod.mean;
		//}
		
		
		
		// use getRanks function to map genes to experiments
		Map<ExpressionExperiment, Map<Gene, Map<DesignElement,Double[]>>> expressRankings = processedExpressionDataVectorService.getRanksProbes(eeCol, pars);
//		Map<ExpressionExperiment, Map<Gene, Map<DesignElement,Double[]>>> expressRankings = new HashMap<ExpressionExperiment, Map<Gene, Map<DesignElement,Double[]>>>();
		//Map<ExpressionExperiment, Map<Gene, Collection<Double>>> expressRankings = processedExpressionDataVectorService.getRanks(eeCol, pars, RankMethod.mean);
		
		
		//Map<ExpressionExperiment, Map<Gene, Collection<Double>>> expressRankings = processedExpressionDataVectorService.getRanks(eeCol, pars, method);
		
		// need to map experiments to genes to calculate max/mean across experiments
		Map<Gene, Double[]> allGeneRankings = new HashMap<Gene, Double[]>();
		Map<Gene, int[]> geneCounts = new HashMap<Gene, int[]>();
		
//		Map<Gene, Map<DesignElement, Collection<Double>>> allProbeRankingsMeanrank = new HashMap<Gene, Map<DesignElement,Collection<Double>>>();
//		Map<Gene, Map<DesignElement, Collection<Double>>> allProbeRankingsMaxrank = new HashMap<Gene, Map<DesignElement,Collection<Double>>>();
		Map<Gene, Map<DesignElement, Collection<Double[]>>> allProbeRankings = new HashMap<Gene, Map<DesignElement,Collection<Double[]>>>();
		//Map<Gene, Map<DesignElement,Integer>> probeCounts = new HashMap<Gene, Map<DesignElement,Integer>>();
		
		
		for ( ExpressionExperiment ee : expressRankings.keySet() ) {
System.out.println(ee.getId());
			Map<Gene, Map<DesignElement,Double[]>> expressRankingsGene = expressRankings.get(ee);
			for ( Gene g : expressRankingsGene.keySet() ) {
System.out.println(ee.getId()+"t"+g.getId());
				Map<DesignElement, Double[]> expressRankingsGeneCompseq = expressRankingsGene.get(g);
				
				//double[] maxAndAveRank = new double[2];
				//maxAndAveRank[0] = 0;  // max
				//maxAndAveRank[1] = 0;  // mean
				
				//int numNullRanks = 0;  // keep track of how many null rankings there are 
				//boolean allNull = true;
				
				//for ( Double d : expressRankings.get(ee).get(g) ) {
				//for ( DesignElement de : expressRankings.get(ee).get(g).keySet() ) {
				for ( DesignElement de : expressRankingsGeneCompseq.keySet() ) {
					
					//Double d = expressRankings.get(ee).get(g).get(de);
					Double[] d = expressRankingsGeneCompseq.get(de);
					
					// subtract this entry from the list of probes, and if necessary, the entire experiment
					if (d == null) {
						//numNullRanks++;
						System.out.println("Null rank value "+g.getId()+", with probeID: "+ de.getId());
						
					} else {
						
System.out.println(ee.getId()+"\t"+g.getId()+"\t"+d[0]);
						if (! allProbeRankings.containsKey(g)) {
							allProbeRankings.put(g, new HashMap<DesignElement,Collection<Double[]>>());
							
							//allProbeRankings.put(g, new HashMap<DesignElement,Collection<Double>>());
							//allProbeRankingsMaxrank.put( g, new HashMap<DesignElement,Collection<Double>>());
						}
						
						if (! allProbeRankings.get(g).containsKey(de)) {
							allProbeRankings.get(g).put(de, new ArrayList<Double[]>());
							//allProbeRankingsMaxrank.get( g).put(de, new ArrayList<Double>());
						}
						
						allProbeRankings.get(g).get(de).add(d);
						//allProbeRankingsMaxrank.get( g).get(de).add(d[1]);
						
						// calculate the maximum and average rankings
						//maxAndAveRank[1] += d;
						//if (maxAndAveRank[0] < d) maxAndAveRank[0] = d;
						
						//allNull = false;
					}
					
					
				}
				
				
				// check if genes for this experiment are null-ranked. if so, skip
				//if (allNull) {
				//	System.out.println("Null: No exps for this gene "+g.getId()+", exp: "+ ee.getId());
				//	continue;
				//}
				
				
				//maxAndAveRank[1] = maxAndAveRank[1] / (expressRankings.get(ee).get(g).size() - numNullRanks);
				
				
				//update hash of expression levels
				//if (allGeneRankings.containsKey(g)) {
				//	if (allGeneRankings.get(g)[1] < maxAndAveRank[0])
				//		allGeneRankings.get(g)[1] = maxAndAveRank[0];
				//	allGeneRankings.get(g)[1] += maxAndAveRank[1];
				//	
				//	geneCounts.get(g)[0]++;
				//	geneCounts.get(g)[1] += expressRankings.get(ee).get(g).size();
				//} 
				// new gene, must create new entities in hash
				//else {
				//	allGeneRankings.put(g, maxAndAveRank);
				//	
				//	
				//	// number of experiments and number samples for this gene
				//	int[] counts = new int[2];
				//	counts[0] = 1;
				//	counts[1] = expressRankings.get(ee).get(g).size();
				//	
				//	geneCounts.put(g, counts);
				//}
				
			}
		}
		
		
		//allProbeRankings.get(g).put(de, new ArrayList<Double>());
		for ( Gene g : allProbeRankings.keySet() ) {
			
			Map<DesignElement,Collection<Double[]>> allProbeRankingsGenes = allProbeRankings.get(g);
			
			//int numExperiments = geneCounts.get(g)[0];
			//int numSamples = geneCounts.get(g)[1];
			//
			//double maxrank  = allGeneRankings.get(g)[0];
			//double meanrank = allGeneRankings.get(g)[1] / numExperiments;
			//
			////px.println(parFileEntries.get(g.getId())+","+numExperiments+","+numSamples+","+maxrank);
			////pe.println(parFileEntries.get(g.getId())+","+numExperiments+","+numSamples+","+meanrank);
			
			for (DesignElement de : allProbeRankingsGenes.keySet()) {
				
				
				// element 0 is the mean, element 1 is the max ranking

				// from list of probe expressios, do stats
				
				int numOfExperiments = allProbeRankingsGenes.get(de).size();
				
				double valExpMean_RankMean = 0;
				double valExpMean_RankMax = 0;
				double valExpMax_RankMean = 0;
				double valExpMax_RankMax = 0;
				
				
				for (Double[] entry : allProbeRankingsGenes.get(de)) {
					
//					if (entry[0] == null || entry[1] == null) continue;
					
					double valRankMean = entry[0].doubleValue();
					
					valExpMean_RankMean += valRankMean;
					if (valExpMax_RankMean < valRankMean)
						valExpMax_RankMean = valRankMean;
					
					double valRankMax  = entry[1].doubleValue();
					
					valExpMean_RankMax += valRankMax;
					if (valExpMax_RankMax < valRankMax) {
						valExpMax_RankMax = valRankMax;
					}
					
				}
				
				
				
				valExpMean_RankMean = valExpMean_RankMean / numOfExperiments;
				valExpMean_RankMax = valExpMean_RankMax / numOfExperiments;

				
				
				p.println(parFileEntries.get(g.getId())+","+de.getId()+","+numOfExperiments+","+valExpMean_RankMean+","+valExpMean_RankMax+","+valExpMax_RankMean+","+valExpMax_RankMax);

				
				//px.println(g.getId()+","+numExperiments+","+numSamples+","+maxrank);
				//pe.println(g.getId()+","+numExperiments+","+numSamples+","+meanrank);
			}
			
			
		}
		
		
		
		
	}


	
	
	
	// the original - look at PARs as a whole, not distinguishing probes
	private void outputAll(Collection<ExpressionExperiment> eeCol, Collection<Gene> pars, PrintStream px, PrintStream pe, String rankMethodStr) {
		
		RankMethod method;
		
		// todo: make these boolean
		if (rankMethodStr.equalsIgnoreCase("max")) {
			method = RankMethod.max;
		} else {
			method = RankMethod.mean;
		}
		
		
		
		// use getRanks function to map genes to experiments
		Map<ExpressionExperiment, Map<Gene, Collection<Double>>> expressRankings = processedExpressionDataVectorService.getRanks(eeCol, pars, method);
		
		// need to map experiments to genes to calculate max/mean across experiments
		HashMap<Gene, double[]> allGeneRankings = new HashMap<Gene, double[]>();
		HashMap<Gene, int[]> geneCounts = new HashMap<Gene, int[]>();
		
		for ( ExpressionExperiment ee : expressRankings.keySet() ) {
			for ( Gene g : expressRankings.get(ee).keySet() ) {
				
				double[] maxAndAveRank = new double[2];
				maxAndAveRank[0] = 0;  // max
				maxAndAveRank[1] = 0;  // mean
								
				int numNullRanks = 0;  // keep track of how many null rankings there are 
				boolean allNull = true;
				
				for ( Double d : expressRankings.get(ee).get(g) ) {
					
					// subtract this entry from the list of probes, and if necessary, the entire experiment
					if (d == null) {
						numNullRanks++;
						System.out.println("Null rank value "+g.getId()+", numNullRanks: "+ numNullRanks);
						
					} else {
						
						// calculate the maximum and average rankings
						maxAndAveRank[1] += d;
						if (maxAndAveRank[0] < d) maxAndAveRank[0] = d;
						
						allNull = false;
					}
					
					
				}
				
				
				// check if genes for this experiment are null-ranked. if so, skip
				if (allNull) {
					System.out.println("Null: No exps for this gene "+g.getId()+", exp: "+ ee.getId());
					continue;
				}
				
				
				maxAndAveRank[1] = maxAndAveRank[1] / (expressRankings.get(ee).get(g).size() - numNullRanks);
				
				
				//update hash of expression levels
				if (allGeneRankings.containsKey(g)) {
					if (allGeneRankings.get(g)[0] < maxAndAveRank[0])
						allGeneRankings.get(g)[0] = maxAndAveRank[0];
					allGeneRankings.get(g)[1] += maxAndAveRank[1];
					
					geneCounts.get(g)[0]++;
					geneCounts.get(g)[1] += expressRankings.get(ee).get(g).size();
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
		
		
		for ( Gene g : allGeneRankings.keySet() ) {
			
			int numExperiments = geneCounts.get(g)[0];
			int numSamples = geneCounts.get(g)[1];
			
			double maxrank  = allGeneRankings.get(g)[0];
			double meanrank = allGeneRankings.get(g)[1] / numExperiments;
			
			px.println(parFileEntries.get(g.getId())+","+numExperiments+","+numSamples+","+maxrank);
			pe.println(parFileEntries.get(g.getId())+","+numExperiments+","+numSamples+","+meanrank);
			
			//px.println(g.getId()+","+numExperiments+","+numSamples+","+maxrank);
			//pe.println(g.getId()+","+numExperiments+","+numSamples+","+meanrank);
			
		}
		
		
		
		
		
	}

	
	
	
	
	private static HashMap getIndices (String header) {
		String[] labels = header.trim().split("\t");
		HashMap<String, Integer> hash = new HashMap<String, Integer>(labels.length);
		
		for (int i=0; i<labels.length; i++) {
			hash.put(labels[i], i);
		}
		
		return hash;
	}
	
	private static int getIndex( String header) {
		return ((Integer) headerLookup.get(header)).intValue();
	}
	
	private static void readPARFile (String inFile) {
		
		BufferedReader in;
		Collection<String[]> fRecords = new ArrayList<String []>();
		HashMap fHash = null;
		String[] fHeaders = null;
		
		String header;
		
		try {
			in = new BufferedReader( new FileReader( inFile ) );
			String line;
			header = in.readLine();
			fHash = getIndices(header);
			fHeaders = header.trim().split("\t");
			
			while ( ( line = in.readLine() ) != null ) {
			    if ( line.startsWith( "#" ) ) continue;
			    String[] s = line.trim().split("\t");
			    fRecords.add(s);
			}
			
			in.close();
			
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.out.println("File "+ inFile +" not found - " + e.getMessage());
			System.exit(0);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.out.println("File "+ inFile +" reads a bit wonky - " + e.getMessage());
			System.exit(0);
		}
		
		headerLookup = fHash;
		headers = fHeaders;
		records = fRecords;
	}
	


}
