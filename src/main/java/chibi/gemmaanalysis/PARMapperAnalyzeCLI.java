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
import ubic.gemma.model.expression.experiment.ExpressionExperiment;
import ubic.gemma.model.expression.experiment.ExpressionExperimentService;
import ubic.gemma.model.genome.Gene;
import ubic.gemma.model.genome.Taxon;
import ubic.gemma.model.genome.TaxonService;
import ubic.gemma.model.genome.gene.GeneService;
import ubic.gemma.util.AbstractSpringAwareCLI;

public class PARMapperAnalyzeCLI extends AbstractSpringAwareCLI {
	
	
	
	private static HashMap headerLookup;
	private static String[] headers;
	private static Collection<String[]> records;
	private static HashMap parFileEntries;

	GeneService parService;
	TaxonService taxonService;
	//ExpressionExperimentSetService expressionExperimentSetService;
	ExpressionExperimentService expressionExperimentService;
	ProcessedExpressionDataVectorService processedExpressionDataVectorService;
	
	
	
	@Override
	protected void buildOptions() {
		// TODO Auto-generated method stub

	}

	@Override
	protected Exception doWork(String[] args) {
		// TODO Auto-generated method stub
		
		this.parService = ( GeneService ) this.getBean( "geneService" );
		this.taxonService = ( TaxonService ) this.getBean( "taxonService" );
		//this.expressionExperimentSetService = ( ExpressionExperimentSetService ) this.getBean( "expressionExperimentSetService" );
		this.expressionExperimentService = ( ExpressionExperimentService ) this.getBean( "expressionExperimentService" );
		this.processedExpressionDataVectorService = ( ProcessedExpressionDataVectorService ) this.getBean( "processedExpressionDataVectorService" );
		
		
		//log.info(taxonService.loadAll());
		
		// For my mac
//		String inFile =  "/Users/mokada/development/PARs/data/human-par-gene-relations.txt";
//		String outFileDir = "/Users/mokada/development/PARs/data";
		
		String inFile = "/home/hmokada/scratch/human-par-gene-relations.txt";
		String outFileDir = "/home/hmokada/scratch";
		
		int batchSize = 100;
		
		
		PrintStream pxx = null; // declare a print stream object
		PrintStream pxe = null;
		PrintStream pex = null;
		PrintStream pee = null;
		
		try
		{
			// Connect print stream to the output stream
			pxx = new PrintStream(new FileOutputStream(outFileDir + "/human-par-gene-relations.output.xx.txt"));
			pxe = new PrintStream(new FileOutputStream(outFileDir + "/human-par-gene-relations.output.xe.txt"));
			pex = new PrintStream(new FileOutputStream(outFileDir + "/human-par-gene-relations.output.ex.txt"));
			pee = new PrintStream(new FileOutputStream(outFileDir + "/human-par-gene-relations.output.ee.txt"));

		} catch (Exception e) {
			System.err.println("Error writing to file");
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
		String outputHeader = "ParID,ParName,Chrom,Nuc,GeneId,GeneSymbol,Distance,GeneContainsPar,SameStrand,NumExperiments,NumSamples,Rank";
		pxx.println(outputHeader);
		pxe.println(outputHeader);
		pex.println(outputHeader);
		pee.println(outputHeader);

		
		//System.out.println("\n\nExperiments by taxon\n\n");
		Taxon taxon_human = taxonService.findByCommonName( "human" );
		
		Collection<ExpressionExperiment> eeCol = expressionExperimentService.findByTaxon(taxon_human);
		System.out.println(eeCol.size());
		
		
		Iterator recordItr = records.iterator();
		int batch = 0;
		
		String[] record;
		while (recordItr.hasNext()) {
			
			parFileEntries = new HashMap<Long, String>(batchSize);
			
			
			Collection<Gene> pars = new ArrayList<Gene>();
			Collection<Gene> genes = new ArrayList<Gene>();
			int count = batchSize;
			
			// work with a small batch
			while (recordItr.hasNext() && 0 < count) {
				
				record = (String[]) recordItr.next();

				
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
				
				
				// todo: place as log
				if (par == null || g == null) {
					System.out.println("PAR or Gene doesn't exist: "+ par.getId() +"\t"+ g.getId());
					continue;
				}
				
				pars.add(par);
				genes.add(g);
				
				count--;  // count down the number of genes in a batch
				
			}
			
			
			
			
			// print results to file
			
			pxx.println("Batch Number: "+ batch);
			pxe.println("Batch Number: "+ batch);
			pex.println("Batch Number: "+ batch);
			pee.println("Batch Number: "+ batch);
			
			batch++;
			
			pxx.flush();
			pxe.flush();
			pex.flush();
			pee.flush();
			
			
			outputAll(eeCol, pars, pxx, pxe, "max");
			outputAll(eeCol, pars, pex, pee, "mean");
			
			
			
		}
		pxx.close();
		pxe.close();
		pex.close();
		pee.close();
		
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
	
	
	
	
	
	
	private void outputAll(Collection<ExpressionExperiment> eeCol, Collection<Gene> genes, PrintStream px, PrintStream pe, String rankMethodStr) {
		
		RankMethod method;
		
		// todo: make these boolean
		if (rankMethodStr.equalsIgnoreCase("max")) {
			method = RankMethod.max;
		} else {
			method = RankMethod.mean;
		}
		
		
		
		// use getRanks function to map genes to experiments
		Map<ExpressionExperiment, Map<Gene, Collection<Double>>> expressRankings = processedExpressionDataVectorService.getRanks(eeCol, genes, method);
		
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
					if (allGeneRankings.get(g)[1] < maxAndAveRank[0])
						allGeneRankings.get(g)[1] = maxAndAveRank[0];
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
