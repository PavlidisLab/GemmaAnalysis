package chibi.gemmaanalysis;

//import ubic.gemma.model.expression.experiment.ExpressionExperimentService;
import ubic.gemma.model.genome.Taxon;
import ubic.gemma.model.genome.TaxonService;
import ubic.gemma.model.genome.gene.GeneService;
import ubic.gemma.util.AbstractSpringAwareCLI;


/**
 * Tester file to get acquanted with Gemma.  Sorry for the mess.
 * 
 * @author mark
 * @version $Id$
 */
public class PARMapperTestFileTester extends AbstractSpringAwareCLI {
	
	
	GeneService parService;
	TaxonService taxonService;
	
	@Override
	protected void buildOptions() {
		// TODO Auto-generated method stub

	}

	@Override
	protected Exception doWork(String[] args) {
		 Exception exc = processCommandLine( "test", args );
		 
		// TODO Auto-generated method stub
//		ExpressionExperimentService eeService = (ExpressionExperimentService)getBean("expressionExperimentService");
		this.parService = ( GeneService ) this.getBean( "geneService" );
 		this.taxonService = ( TaxonService ) this.getBean( "taxonService" );
 		
 		Taxon taxon = taxonService.findByCommonName( "mouse" );
 		long taxonId = 1;
 		
 		log.info(taxonService.loadAll());
 		log.info(taxonService.load(taxonId));
 		log.info(taxonService.find(taxon));
 		
//		int taxonId = 1; // human=1
//		String taxonLabel = "beef";//taxonService.toString(); 
		
		System.out.println("Reading file");
		System.out.println(taxon.toString());
		//System.out.println("Reading the file");
		
		return null;
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		PARMapperTestFileTester p = new PARMapperTestFileTester();
        Exception e = p.doWork( args );
        if ( e != null ) {
            throw new RuntimeException( e );
        }
	}

}
