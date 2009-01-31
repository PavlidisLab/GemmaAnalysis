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
		
		
		
		String inFile = "/home/hmokada/scratch/human-par-gene-relations.subset.txt";
		String outFile = "/home/hmokada/scratch/human-par-gene-relations.output.subset.txt";
		int batchSize = 100;
		
		
		PrintStream p = null; // declare a print stream object
		try
		{
			// Connect print stream to the output stream
			p = new PrintStream(new FileOutputStream(outFile));
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
		
		
		
		
		
		
		
		
		

		
		System.out.println("\n\nExperiments by taxon\n\n");
		Taxon taxon_human = taxonService.findByCommonName( "human" );
		
		System.out.println(taxon_human);
		
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
			
			Collection<Gene> genes = new ArrayList<Gene>();
			int count = batchSize;
			
			// work with a small batch
			while (recordItr.hasNext() && 0 < count) {
				
				record = (String[]) recordItr.next();

				
				// accessing data elements
				int ParID			= Integer.parseInt( record[ParIDIdx] );
				String ParName		= record[ ParNameIdx];
				String Chrom		= record[ ChromIdx];
				int Nuc				= Integer.parseInt( record[NucIdx] );
				int GeneId			= Integer.parseInt( record[GeneIdIdx] );
				String GeneSymbol	= record[ GeneSymbolIdx ];
				
				Gene g = this.parService.load(GeneId);
				
				if (g == null) {
					System.out.println("Gene doesn't exist: "+ g);
					continue;
				}
				
				genes.add(g);
				System.out.println("Gene "
						//+g
						+"\t|\t"+GeneId
						+"\t|\t"+g.getName()
						+"\t|\t"+g.getOfficialSymbol()
						//+"\t|\t"+g.getDescription()
						//+"\t|\t"+g.getNcbiId()
						//+"\t|\t"+g.getTaxon()
						);
				
				count--;
				
			}
			
			
			
			
			
			
			
			
			p.println("Batch Number: "+ batch++);
			//RankMethod method = RankMethod.max;
			outputAll(eeCol, genes, p, "max");
			
			//method = RankMethod.mean;
			outputAll(eeCol, genes, p, "mean");

			
			
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
		p.close();
		
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
	
	
	
	
	
	
	private void outputAll(Collection<ExpressionExperiment> eeCol, Collection<Gene> genes, PrintStream p, String methodStr) {
		
		
		RankMethod method;
		
		if (methodStr.equalsIgnoreCase("max")) {
			method = RankMethod.max;
		} else {
			method = RankMethod.mean;
		}
		
		Map<ExpressionExperiment, Map<Gene, Collection<Double>>> res = processedExpressionDataVectorService.getRanks(eeCol, genes, method);
		
		for ( ExpressionExperiment ee : res.keySet() ) {
			//p.println(">> Exp: " + ee.getName());
			for ( Gene g : res.get(ee).keySet() ) {
				//p.println("> Gene: " + g.getName());
				for ( Double d : res.get(ee).get(g) ) {
					//p.print(d+"\t");
					p.println(g.getId()+","+ee.getId()+","+d+","+methodStr);
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
