package chibi.gemmaanalysis;

import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.Collection;
import java.util.Map;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.lang3.time.StopWatch;

import chibi.gemmaanalysis.CoexpressionAnalysisService.CoexpressionMatrices;
import ubic.basecode.dataStructure.matrix.DenseDouble3dMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.io.writer.MatrixWriter;
import ubic.gemma.core.analysis.preprocess.filter.FilterConfig;
import ubic.gemma.core.genome.gene.service.GeneService;
import ubic.gemma.core.ontology.providers.GeneOntologyService;
import ubic.gemma.model.expression.experiment.BioAssaySet;
import ubic.gemma.model.expression.experiment.ExpressionExperiment;
import ubic.gemma.model.genome.Gene;
import ubic.gemma.model.genome.Taxon;
import ubic.gemma.persistence.service.expression.experiment.ExpressionExperimentService;

/**
 * Calculate the effect size
 *
 * @author xwan
 * @author raymond
 */
public class EffectSizeCalculationCli extends AbstractGeneCoexpressionManipulatingCLI {
    public static final int DEFAULT_STRINGENCY = 3;

    public static void main( String[] args ) {
        EffectSizeCalculationCli analysis = new EffectSizeCalculationCli();
        StopWatch watch = new StopWatch();
        watch.start();
        log.info( "Starting Effect Size Analysis" );
        Exception exc = analysis.doWork( args );
        if ( exc != null ) {
            log.error( exc.getMessage() );
        }
        log.info( "Finished analysis in " + watch.getTime() / 1000 + " seconds" );
    }

    private String goTerm;

    private String outFilePrefix;

    private CoexpressionAnalysisService coexpressionAnalysisService;

    private GeneOntologyService goService;

    public EffectSizeCalculationCli() {
        super();
    }

    /*
     * (non-Javadoc)
     *
     * @see ubic.gemma.util.AbstractCLI#getCommandName()
     */
    @Override
    public String getCommandName() {
        // TODO Auto-generated method stub
        return null;
    }

    @SuppressWarnings("static-access")
    @Override
    protected void buildOptions() {
        super.buildOptions();
        OptionBuilder.hasArg();
        OptionBuilder.withArgName( "GOTerm" );
        OptionBuilder.withDescription( "Target GO term" );
        OptionBuilder
                .withLongOpt( "GOTerm" );
        Option goOption = OptionBuilder.create( 'g' );
        addOption( goOption );
        OptionBuilder.hasArg();
        OptionBuilder.isRequired();
        OptionBuilder.withArgName( "outFilePrefix" );
        OptionBuilder
                .withDescription( "File prefix for saving the correlation data" );
        OptionBuilder.withLongOpt( "outFilePrefix" );
        Option outputFileOption = OptionBuilder
                .create( 'o' );
        addOption( outputFileOption );
        OptionBuilder.hasArg();
        OptionBuilder.withArgName( "stringency" );
        OptionBuilder
                .withDescription( "Vote count stringency for link selection" );
        OptionBuilder.withLongOpt( "stringency" );
        Option stringencyOption = OptionBuilder.create( 'r' );
        addOption( stringencyOption );
    }

    @SuppressWarnings({ "rawtypes", "unchecked" })
    @Override
    protected Exception doWork( String[] args ) {
        Exception exc = processCommandLine( args );
        if ( exc != null ) {
            return exc;
        }
        StopWatch watch = new StopWatch();
        watch.start();

        Collection<Gene> queryGenes, targetGenes;
        try {

            queryGenes = getQueryGenes();
            targetGenes = getTargetGenes();
        } catch ( IOException e ) {
            return e;
        }
        if ( goTerm != null ) {
            while ( !goService.isReady() ) {
                try {
                    Thread.sleep( 1000 );
                } catch ( InterruptedException e ) {
                }
            }
            targetGenes.addAll( goService.getGenes( goTerm, this.getTaxon() ) );
        }

        if ( targetGenes.size() == 0 || queryGenes.size() == 0 ) {
            return new Exception( "No genes in query/target" );
        }

        FilterConfig filterConfig = new FilterConfig();
        CoexpressionMatrices matrices = coexpressionAnalysisService.calculateCoexpressionMatrices(
                getExpressionExperiments(), queryGenes, targetGenes, filterConfig, null );
        DenseDouble3dMatrix<Gene, Gene, BioAssaySet> correlationMatrix = matrices.getCorrelationMatrix();
        DenseDouble3dMatrix<Gene, Gene, BioAssaySet> sampleSizeMatrix = matrices.getSampleSizeMatrix();
        DoubleMatrix<Gene, Gene> effectSizeMatrix = coexpressionAnalysisService.calculateEffectSizeMatrix(
                correlationMatrix, sampleSizeMatrix );

        // create 2D correlation heat map
        // ColorMatrix dataColorMatrix = new ColorMatrix( correlationMatrix2D );
        // dataColorMatrix.setColorMap( ColorMap.GREENRED_COLORMAP );
        // JMatrixDisplay dataMatrixDisplay = new JMatrixDisplay( dataColorMatrix );
        // String figureFileName = outFilePrefix + ".corr.png";

        // create row/col name maps
        Map<ExpressionExperiment, String> eeNameMap = matrices.getEeNameMap();

        DecimalFormat formatter = ( DecimalFormat ) NumberFormat.getNumberInstance();
        formatter.applyPattern( "0.0000" );
        DecimalFormatSymbols symbols = formatter.getDecimalFormatSymbols();
        symbols.setNaN( "" );
        formatter.setDecimalFormatSymbols( symbols );
        String topLeft = "GenePair";
        try {
            MatrixWriter out = new MatrixWriter( outFilePrefix + ".corr.txt", formatter );
            out.setColNameMap( eeNameMap );
            out.setTopLeft( topLeft );
            out.writeMatrix( correlationMatrix, true );

            out = new MatrixWriter( outFilePrefix + ".effect_size.txt", formatter );
            out.setColNameMap( eeNameMap );
            out.setTopLeft( topLeft );
            out.writeMatrix( effectSizeMatrix, true );

            // dataMatrixDisplay.saveImage( figureFileName, true );
            // log.info( "Saved correlation image to " + figureFileName );
        } catch ( IOException e ) {
            return e;
        }

        return null;
    }

    protected void initBeans() {
        coexpressionAnalysisService = this.getBean( CoexpressionAnalysisService.class );

    }

    @Override
    protected void processOptions() {
        super.processOptions();
        if ( hasOption( 'g' ) ) {
            this.goTerm = getOptionValue( 'g' );
        }
        if ( hasOption( 'g' ) ) {
            this.goTerm = getOptionValue( 'g' );
        }

        if ( hasOption( 'o' ) ) {
            this.outFilePrefix = getOptionValue( 'o' );
        }

        initBeans();
    }
}