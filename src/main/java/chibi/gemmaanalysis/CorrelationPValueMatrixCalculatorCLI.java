package chibi.gemmaanalysis;

import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Locale;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;

import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.io.reader.DoubleMatrixReader;
import ubic.basecode.io.writer.MatrixWriter;
import ubic.gemma.core.apps.ExpressionExperimentManipulatingCLI;

/**
 * CLI for reading in a max correlation matrix to calculate p values from correlation histograms
 *
 * @author raymond
 */
public class CorrelationPValueMatrixCalculatorCLI extends ExpressionExperimentManipulatingCLI {

    public static void main( String[] args ) {
        CorrelationPValueMatrixCalculatorCLI cli = new CorrelationPValueMatrixCalculatorCLI();
        Exception exc = cli.doWork( args );
        if ( exc != null ) log.error( exc.getMessage() );

    }

    private String inFile;

    private String outFile;

    private CoexpressionAnalysisService coexpService;

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

    @Override
    protected void buildOptions() {
        super.buildOptions();
        Option inFileo = OptionBuilder.create( 'i' );
        addOption( inFileo );
        Option outFileo = OptionBuilder.create( 'o' );
        addOption( outFileo );
    }

    @Override
    protected Exception doWork( String[] args ) {
        processCommandLine( args );

        //        String taxonName = "human";
        //        taxon = Taxon.Factory.newInstance();
        //        taxon.setCommonName( taxonName );
        //        taxon = taxonService.find( taxon );

        DecimalFormat formatter = ( DecimalFormat ) NumberFormat.getNumberInstance( Locale.US );
        // formatter.applyPattern("0.0000");

        DoubleMatrixReader in = new DoubleMatrixReader();
        try {
            DoubleMatrix<String, String> matrix = in.read( inFile );
            DoubleMatrix<String, String> pMatrix = coexpService.calculateMaxCorrelationPValueMatrix( matrix, 0,
                    getExpressionExperiments() );
            MatrixWriter<String, String> out = new MatrixWriter<>( outFile, formatter );
            out.writeMatrix( pMatrix, true );
        } catch ( IOException e ) {
            return e;
        }

        return null;
    }

    @Override
    protected void processOptions() {
        super.processOptions();
        inFile = getOptionValue( 'i' );
        outFile = getOptionValue( 'o' );

        coexpService = getBean( CoexpressionAnalysisService.class );
    }

}
