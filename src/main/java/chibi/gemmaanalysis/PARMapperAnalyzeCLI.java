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

import ubic.gemma.model.analysis.expression.ExpressionExperimentSet;
import ubic.gemma.model.analysis.expression.ExpressionExperimentSetService;
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
		
		Exception exc = processCommandLine( "test", args );
		
		// TODO Auto-generated method stub
//			ExpressionExperimentService eeService = (ExpressionExperimentService)getBean("expressionExperimentService");
		this.parService = ( GeneService ) this.getBean( "geneService" );
		this.taxonService = ( TaxonService ) this.getBean( "taxonService" );
		//this.expressionExperimentSetService = ( ExpressionExperimentSetService ) this.getBean( "expressionExperimentSetService" );
		this.expressionExperimentService = ( ExpressionExperimentService ) this.getBean( "expressionExperimentService" );
		this.processedExpressionDataVectorService = ( ProcessedExpressionDataVectorService ) this.getBean( "processedExpressionDataVectorService" );
		
		
		Taxon taxon = taxonService.findByCommonName( "mouse" );
		long taxonId = 1;
		
		log.info(taxonService.loadAll());
		log.info(taxonService.load(taxonId));
		log.info(taxonService.find(taxon));
		
		//int taxonId = 1; // human=1
		//String taxonLabel = "beef";//taxonService.toString(); 
		
		
		//Gene g1 = this.parService.load(889481); //System.out.println("Gene 889481\t"+this.parService.getCompositeSequencesById(new Integer(889481).longValue()));
		//System.out.println("Gene 889481\t"+g1+"\t|\t"+g1.getName()+"\t|\t"+g1.getDescription()+"\t|\t"+g1.getNcbiId()+"\t|\t"+g1.getTaxon());
		//System.out.println("Gene 3609073\t"+this.parService.getCompositeSequencesById(new Integer(3609073).longValue()));
		
		//System.out.println("Reading file");
		System.out.println(taxon.toString());
		//System.out.println("Reading the file");
		
		
		
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

			//p.println ("This is written to a file");
		} catch (Exception e) {
			System.err.println("Error writing to file");
		}
		
		
		
		
		
		if (0 < args.length && args[0] != null) {
			inFile = args[0];
		}
		
		System.out.println("Reading file: " + inFile);
		
		
		
		//HashMap headerLookup;
		//String headers
		
		readPARFile(inFile);
		
		
		/*for (int i=0; i<records.size(); i++) {
			System.out.println(records[i]);
		}*/
		
		
		
		
		int ParIDIdx		= getIndex("ParID");
		int ParNameIdx		= getIndex("ParName");
		int ChromIdx		= getIndex("Chrom");
		int NucIdx			= getIndex("Nuc");
		int GeneIdIdx		= getIndex("GeneId");
		int GeneSymbolIdx	= getIndex("GeneSymbol");
		int DistanceIdx		= getIndex("Distance");
		int GeneContainsParIdx	= getIndex("GeneContainsPar");
		int SameStrandIdx	= getIndex("SameStrand");
		
		
		String outputHeader = "ParID,ParName,Chrom,Nuc,GeneId,GeneSymbol,Distance,GeneContainsPar,SameStrand,NumExperiments,NumSamples,Rank";
		pxx.println(outputHeader);
		pxe.println(outputHeader);
		pex.println(outputHeader);
		pee.println(outputHeader);

		
		//System.out.println("\n\nExperiments by taxon\n\n");
		Taxon taxon_human = taxonService.findByCommonName( "human" );
		//System.out.println(taxon_human);
		
		Collection<ExpressionExperiment> eeCol = expressionExperimentService.findByTaxon(taxon_human);
		System.out.println(eeCol.size());
		/*for ( ExpressionExperiment ees : eeCol ) {
			System.out.println(ees.getDescription()
					+"\t|\t"+ ees.toString()
					//+"\t|\t"+ ees.getTaxon()  //why does this give me errors???
					+"\t|\t"+ ees.getId()
					);
		}*/
		
		
		
		
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
					System.out.println("PAR or Gene doesn't exist: "+ g);
					continue;
				}
				
				pars.add(par);
				genes.add(g);
				/*
				System.out.println("Gene "
						//+g
						+"\t|\t"+GeneId
						//+"\t|\t"+g.getPhysicalLocation()  // gives errors
						+"\t|\t"+g.getName()
						+"\t|\t"+g.getOfficialSymbol()
						//+"\t|\t"+g.getDescription()
						//+"\t|\t"+g.getNcbiId()
						//+"\t|\t"+g.getTaxon()
						);
				
				System.out.println("PAR  "
						//+g
						+"\t|\t"+ParID
						//+"\t|\t"+par.getPhysicalLocation()  // gives errors
						+"\t|\t"+par.getName()
						+"\t|\t"+par.getOfficialSymbol()
						//+"\t|\t"+g.getDescription()
						//+"\t|\t"+g.getNcbiId()
						//+"\t|\t"+g.getTaxon()
						);
				*/
				count--;
				
			}
			
			
			
			
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
			
			
			//RankMethod method = RankMethod.max;
//			outputAll(eeCol, genes, pxx, "max", "max");
			
			//method = RankMethod.mean;
//			outputAll(eeCol, genes, pex, "mean", "max");

			//RankMethod method = RankMethod.max;
//			outputAll(eeCol, genes, pxe, "max", "mean");
			
			//method = RankMethod.mean;
//			outputAll(eeCol, genes, pee, "mean", "mean");
			
			
			/*Map<ExpressionExperiment, Map<Gene, Collection<Double>>> res = processedExpressionDataVectorService.getRanks(eeCol, genes, method);
			
			for ( ExpressionExperiment ee : res.keySet() ) {
				//p.println(">> Exp: " + ee.getName());
				for ( Gene g : res.get(ee).keySet() ) {
					//p.println("> Gene: " + g.getName());
					for ( Double d : res.get(ee).get(g) ) {
						//p.print(d+"\t");
					}
					//p.println();
				}
				//p.println();
			}*/
			
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
		//HashMap<Gene, Map<ExpressionExperiment, Collection<Double>>> allGeneRankings = new HashMap<Gene, Map<ExpressionExperiment, Collection<Double>>>();
		HashMap<Gene, double[]> allGeneRankings = new HashMap<Gene, double[]>();
		HashMap<Gene, int[]> geneCounts = new HashMap<Gene, int[]>();
		
		for ( ExpressionExperiment ee : expressRankings.keySet() ) {
			//p.println(">> Exp: " + ee.getName());
			for ( Gene g : expressRankings.get(ee).keySet() ) {
				//p.println("> Gene: " + g.getName());
				
				double[] maxAndAveRank = new double[2];
				maxAndAveRank[0] = 0;  // max
				maxAndAveRank[1] = 0;  // mean
				//double averageRank = 0;
				//double maxRank = 0;
				
				//if (expressRankings.get(ee).get(g).size() != 1)
				//	System.out.println("### "+ expressRankings.get(ee).get(g).size());
				
				for ( Double d : expressRankings.get(ee).get(g) ) {
					//p.print(d+"\t");
					//p.println(g.getId()+","+ee.getId()+","+d+","+rankMethodStr);
					
					
					// TODO: subtract this entry from the list of probes, and if necessary, the entire experiment
					if (d == null) {
						//d = 0.0;
						//System.out.println("Null");
						
						
					}
					
					maxAndAveRank[1] += d;
					if (maxAndAveRank[0] < d) maxAndAveRank[0] = d;
					
				}
				
				maxAndAveRank[1] = maxAndAveRank[1] / expressRankings.get(ee).get(g).size();
				
				
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
				
				//p.println();
			}
			//p.println();
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

	
	
	private void outputAll_old(Collection<ExpressionExperiment> eeCol, Collection<Gene> genes, PrintStream p, String rankMethodStr, String probeRankStr) {
		
		RankMethod method;
		
		// todo: make these boolean
		if (rankMethodStr.equalsIgnoreCase("max")) {
			method = RankMethod.max;
		} else {
			method = RankMethod.mean;
		}
		
		Map<ExpressionExperiment, Map<Gene, Collection<Double>>> res = processedExpressionDataVectorService.getRanks(eeCol, genes, method);
		
		int numExperiments = res.size();
		
		for ( ExpressionExperiment ee : res.keySet() ) {
			//p.println(">> Exp: " + ee.getName());
			for ( Gene g : res.get(ee).keySet() ) {
				//p.println("> Gene: " + g.getName());
				for ( Double d : res.get(ee).get(g) ) {
					//p.print(d+"\t");
					p.println(g.getId()+","+ee.getId()+","+d+","+rankMethodStr);
				}
				//p.println();
			}
			//p.println();
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
