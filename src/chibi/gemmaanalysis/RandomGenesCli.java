package chibi.gemmaanalysis;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;

import ubic.gemma.analysis.expression.coexpression.GeneCoexpressionService;
import ubic.gemma.model.association.coexpression.Gene2GeneCoexpressionService;
import ubic.gemma.model.genome.Gene;
import ubic.gemma.model.genome.TaxonService;
import ubic.gemma.model.genome.gene.GeneService;
import ubic.gemma.util.AbstractSpringAwareCLI;
import ubic.gemma.util.AbstractCLI.ErrorCode;

public class RandomGenesCli extends AbstractSpringAwareCLI{

    private GeneService geneService;
    private Integer listSize;
    private String fileName;
    /**
     * @param args
     */
    public static void main( String[] args ) {
        RandomGenesCli cli = new RandomGenesCli();
        cli.buildOptions();
        
        try {
            Exception ex = cli.doWork( args );
            if ( ex != null ) {
                ex.printStackTrace();
            }
        } catch ( Exception e ) {
            throw new RuntimeException( e );
        }
       
    }
    @SuppressWarnings("static-access")
    @Override
    protected void buildOptions() {
        
        Option geneFileOption = OptionBuilder.hasArg().withArgName( "Gene List File Name" ).withDescription(
        "A text file that contains a list of gene IDs, with one gene ID on each line" ).withLongOpt(
        "geneFile" ).create( 'g' );
        addOption( geneFileOption );
    
        Option numberOption = OptionBuilder.hasArg().withArgName( "NumberOfGenes" ).withDescription(
        "The number of genes to use." ).withLongOpt( "number" ).create( 'n' );
        addOption( numberOption );
       
    }
    @Override
    protected Exception doWork( String[] args ) {
        Exception err = processCommandLine( "Get Random Genes ", args );
        if ( err != null ) return err;

        try {
            getRandomGenes(fileName, listSize); //also prints it to file
            //printToFile(genes);
        } catch ( IOException e ) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        return null;
    }
    
    protected void processOptions() {
        super.processOptions();

        //initialize bean
        geneService = ( GeneService ) getBean( "geneService" );
               
        //get list size
        if ( this.hasOption( 'n' ) ) {           
                listSize = Integer.parseInt( getOptionValue('n') );         
            }
        
        //get filename
        if (this.hasOption( 'g' )){
            fileName = getOptionValue('g');
        }
    }
    
    private void getRandomGenes(String inFile, int number) throws IOException{
        Object[] rawLinesInFile = readGeneListFile( inFile ).toArray();
        SecureRandom randomSeed;
        Collection<Gene> genes = new HashSet<Gene>();
        int geneCount = 0;
        String line;
        //boolean added = false;
        try {
            
            BufferedWriter out = new BufferedWriter( new FileWriter( "random"+number+".genes.txt" ));
     
          randomSeed = SecureRandom.getInstance( "SHA1PRNG" );

          while ( geneCount < number ) {

              randomSeed.generateSeed( 10 );
              int randomInt = randomSeed.nextInt( rawLinesInFile.length );

              line = ( String ) rawLinesInFile[randomInt];

              Gene gene = geneService.load( Long.parseLong( line ) );
              if ( gene == null ) {
                  log.error( "ERROR: Cannot find genes for ID: " + line );
                  continue;
              }
              geneService.thaw( gene );

              // because we are using a HashSet, no duplicates are allowed - which is what we want
              log.info( "Index used: "+randomInt );
//              added = genes.add( gene );
//              if (added)
             out.write(gene.getOfficialSymbol()+"\n");
                       
             geneCount++;
          }
          out.close(); 
      } catch ( NoSuchAlgorithmException e ) {
        // TODO Auto-generated catch block
        e.printStackTrace();
    }
         
      
    }
      
      /**
       * Read in a list of genes; calls helper method, <code>readGeneListFile()</code>
       * 
       * @param inFile - file name to read
       * @param taxon
       * @return collection of genes
       * @throws IOException
       */
      private Collection<String> readGeneListFile( String inFile ) throws IOException {
          log.info( "Reading " + inFile );
          Collection<String> lines = new ArrayList<String>();
          BufferedReader in = new BufferedReader( new FileReader( inFile ) );
          String line;
          while ( ( line = in.readLine() ) != null ) {
              // if ( line.startsWith( "#" ) ) continue;
              String s = line.trim();
              lines.add( s );
          }
          in.close();
          return lines;
      }
      
//      private void printToFile(Collection<Gene> genes){
//          try {
//              
//              BufferedWriter out = new BufferedWriter( new FileWriter( "random1000.genes.txt" ));
//
//              for ( Gene gene : genes ) {
//
//                  out.write(gene.getOfficialSymbol()+"\n");
//
//              }
//              out.close();
//          } catch ( IOException e ) {
//
//              e.printStackTrace();
//          }
//             
//          
//      }

}
