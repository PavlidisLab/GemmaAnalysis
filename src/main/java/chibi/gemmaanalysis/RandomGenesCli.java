package chibi.gemmaanalysis;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Random;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;

import ubic.gemma.core.genome.gene.service.GeneService;
import ubic.gemma.core.util.AbstractSpringAwareCLI;
import ubic.gemma.model.genome.Gene;

public class RandomGenesCli extends AbstractSpringAwareCLI {

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

    private GeneService geneService;
    private Integer listSize;

    private String fileName;

    /*
     * (non-Javadoc)
     *
     * @see ubic.gemma.util.AbstractCLI#getCommandName()
     */
    @Override
    public String getCommandName() {
        return null;
    }

    @SuppressWarnings("static-access")
    @Override
    protected void buildOptions() {

        OptionBuilder.hasArg();
        OptionBuilder.withArgName( "Gene List File Name" );
        OptionBuilder
                .withDescription( "A text file that contains a list of gene IDs, with one gene ID on each line" );
        OptionBuilder
                .withLongOpt( "geneFile" );
        Option geneFileOption = OptionBuilder.create( 'g' );
        addOption( geneFileOption );

        OptionBuilder.hasArg();
        OptionBuilder.withArgName( "NumberOfGenes" );
        OptionBuilder
                .withDescription( "The number of genes to use." );
        OptionBuilder.withLongOpt( "number" );
        Option numberOption = OptionBuilder.create( 'n' );
        addOption( numberOption );

    }

    @Override
    protected Exception doWork( String[] args ) {
        Exception err = processCommandLine( args );
        if ( err != null ) return err;

        try {
            getRandomGenes( fileName, listSize ); // also prints it to file
            // printToFile(genes);
        } catch ( IOException e ) {
            e.printStackTrace();
        }
        return null;
    }

    @Override
    protected void processOptions() {
        super.processOptions();

        // initialize bean
        geneService = getBean( GeneService.class );

        // get list size
        if ( this.hasOption( 'n' ) ) {
            listSize = Integer.parseInt( getOptionValue( 'n' ) );
        }

        // get filename
        if ( this.hasOption( 'g' ) ) {
            fileName = getOptionValue( 'g' );
        }
    }

    private void getRandomGenes( String inFile, int number ) throws IOException {
        Object[] rawLinesInFile = readGeneListFile( inFile ).toArray();
        int geneCount = 0;
        String line;
        try (BufferedWriter out = new BufferedWriter( new FileWriter( "random" + number + ".genes.txt" ) );) {

            Random r = new Random();

            while ( geneCount < number ) {

                int randomInt = r.nextInt( rawLinesInFile.length );

                line = ( String ) rawLinesInFile[randomInt];

                Gene gene = geneService.load( Long.parseLong( line ) );
                if ( gene == null ) {
                    log.error( "ERROR: Cannot find genes for ID: " + line );
                    continue;
                }
                gene = geneService.thaw( gene );

                out.write( gene.getOfficialSymbol() + "\n" );

                geneCount++;
            }
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
        Collection<String> lines = new ArrayList<>();
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

}
