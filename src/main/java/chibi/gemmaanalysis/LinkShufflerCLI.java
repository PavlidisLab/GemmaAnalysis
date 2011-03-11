package chibi.gemmaanalysis;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.lang.time.StopWatch;

import ubic.gemma.apps.ExpressionExperimentManipulatingCLI;
import ubic.gemma.model.genome.Gene;

/**
 * Used to count up links and to generate the link(gene pair) background distribution, which could be used to estimate
 * the false positive rates under different levels of confirmation. When shuffling, there are two steps to finish this
 * process. The first step is to prepare the working table. Then the analysis is done using the working table.
 * 
 * @author xwan
 * @author raymond (refactor)
 */
public class LinkShufflerCLI extends ExpressionExperimentManipulatingCLI {
    /*
     * How many shuffled runs to do. 2 or 3 is enough to get a quick idea of what the results will look like; 100 is
     * better for final analysis.
     */
    private int numIterations = 2;

    /*
     * Print out the link details for shuffled data sets. This is only useful for debugging (big files).
     */
    private boolean doShuffledOutput = false;

    /*
     * If false, just do shuffling. This is primarily for debugging.
     */
    private boolean doRealAnalysis = false;

    private String outFileName;

    @Override
    protected Exception doWork( String[] args ) {
        Exception err = processCommandLine( "Shuffle Links ", args );
        if ( err != null ) {
            return err;
        }

        LinkStatisticsService lss = ( LinkStatisticsService ) getBean( "linkStatisticsService" );

        Collection<Gene> genes = geneService.loadKnownGenes( taxon );

        LinkConfirmationStatistics confStats = null;

        if ( doRealAnalysis ) { // Currently this is really just for debugging
            // purposes, though reading in from a
            // file might be useful.
            log.info( "Doing real analysis" );
            LinkStatistics realStats = lss.analyze( expressionExperiments, genes, taxon, false, true );
            log.info( realStats.getTotalLinkCount() + " gene links in total" );
            confStats = realStats.getLinkConfirmationStats();

            try {
                Writer linksOut = new BufferedWriter( new FileWriter( new File( "link-data.txt" ) ) );
                realStats.writeLinks( linksOut, 0 );
            } catch ( IOException e ) {
                return e;
            }
        }

        List<LinkConfirmationStatistics> shuffleRuns = new ArrayList<LinkConfirmationStatistics>();
        System.out.println( "Running shuffled runs" );
        for ( int i = 0; i < numIterations; i++ ) {
            System.out.println( "*** Iteration " + i + " ****" );

            LinkStatistics sr = lss.analyze( expressionExperiments, genes, taxon, true, true );
            log.info( sr.getTotalLinkCount() + " gene links in total" );

            shuffleRuns.add( sr.getLinkConfirmationStats() );

            if ( doShuffledOutput ) {
                try {
                    Writer linksOut = new BufferedWriter(
                            new FileWriter( new File( "shuffled-link-data-" + i + ".txt" ) ) );
                    sr.writeLinks( linksOut, 2 );
                } catch ( IOException e ) {
                    return e;
                }
            }

        }

        if ( numIterations > 0 ) {
            try {
                Writer out = new PrintWriter( outFileName );
                lss.writeStats( out, confStats, shuffleRuns );
            } catch ( IOException e ) {
                return e;
            }
        }

        return null;
    }

    public static void main( String[] args ) {
        LinkShufflerCLI shuffle = new LinkShufflerCLI();
        StopWatch watch = new StopWatch();
        watch.start();
        Exception ex = shuffle.doWork( args );
        if ( ex != null ) {
            ex.printStackTrace();
        }
        watch.stop();
        log.info( "Finished: " + watch );
    }

    @SuppressWarnings("static-access")
    @Override
    protected void buildOptions() {
        super.buildOptions();

        Option iterationNum = OptionBuilder.hasArg().withArgName( " Shuffling iterations " ).withDescription(
                " The number of shuffling iterations (default = 2) " ).withLongOpt( "iterationNum" ).create( 'i' );
        addOption( iterationNum );

        Option outputShuffledLinks = OptionBuilder.withArgName( "Shuffled output" ).withDescription(
                "Print out link details for shuffled data sets" ).withLongOpt( "outputShuffledData" ).create( 'l' );
        addOption( outputShuffledLinks );

        Option realAnalysis = OptionBuilder.withArgName( "Real analysis (unshuffled)" ).withDescription(
                "Perform a real link analysis and output it to link-data.txt" ).withLongOpt( "realAnalysis" ).create(
                'r' );
        addOption( realAnalysis );

        Option outputFile = OptionBuilder.hasArg().isRequired().withArgName( "Output File" ).withDescription(
                "Output shuffle statistics to this file" ).withLongOpt( "outFile" ).create( 'o' );
        addOption( outputFile );

        /*
         * Not supported
         */
        // Option linkStringency = OptionBuilder.hasArg().withArgName( "Link
        // support threshold (stringency)" )
        // .withDescription( "Link Stringency " ).withLongOpt( "linkStringency"
        // ).create( 'l' );
        // addOption( linkStringency );
    }

    @Override
    protected void processOptions() {
        super.processOptions();

        if ( hasOption( 'i' ) ) {
            this.numIterations = getIntegerOptionValue( 'i' );
        }

        if ( hasOption( 'l' ) ) {
            this.doShuffledOutput = true;
        }

        if ( hasOption( 'r' ) ) {
            this.doRealAnalysis = true;
        }

        outFileName = getOptionValue( 'o' );

    }

    @Override
    protected String[] getAdditionalSpringConfigLocations() {
        return new String[] { "classpath*:chibi/beans.xml" };
    }

}
