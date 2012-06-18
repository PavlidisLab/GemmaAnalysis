/*
 * The Gemma project
 * 
 * Copyright (c) 2007 University of British Columbia
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
package chibi.gemmaanalysis;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.Collection;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.lang.time.StopWatch;

import ubic.basecode.dataStructure.matrix.DenseDouble3dMatrix;
import ubic.basecode.io.writer.MatrixWriter;
import ubic.gemma.analysis.preprocess.filter.FilterConfig;
import ubic.gemma.model.expression.experiment.BioAssaySet;
import ubic.gemma.model.expression.experiment.ExpressionExperiment;
import ubic.gemma.expression.experiment.service.ExpressionExperimentService;
import ubic.gemma.model.genome.Gene;
import chibi.gemmaanalysis.CoexpressionAnalysisService.CoexpressionMatrices;
import chibi.gemmaanalysis.CoexpressionAnalysisService.CorrelationMethod;

/**
 * @author raymond
 * @version $Id$
 */
public class CorrelationAnalysisCLI extends AbstractGeneCoexpressionManipulatingCLI {
    private String outFilePrefix;

    private CoexpressionAnalysisService coexpressionAnalysisService;

    private FilterConfig filterConfig;

    public CorrelationAnalysisCLI() {
        super();
        filterConfig = new FilterConfig();
    }

    @SuppressWarnings("static-access")
    @Override
    protected void buildOptions() {
        super.buildOptions();
        Option outputFileOption = OptionBuilder.hasArg().isRequired().withArgName( "File prefix" )
                .withDescription( "File prefix for saving the output" ).withLongOpt( "outFilePrefix" ).create( 'o' );
        addOption( outputFileOption );

        Option kMaxOption = OptionBuilder.hasArg().withArgName( "k" ).withDescription( "Select the kth largest value" )
                .withType( Integer.class ).withLongOpt( "kValue" ).create( 'k' );
        addOption( kMaxOption );
    }

    @Override
    protected void processOptions() {
        super.processOptions();
        if ( hasOption( 'o' ) ) {
            this.outFilePrefix = getOptionValue( 'o' );
        }

        initBeans();
    }

    protected void initBeans() {
        coexpressionAnalysisService = this.getBean( CoexpressionAnalysisService.class );
        eeService = this.getBean( ExpressionExperimentService.class );
    }

    @Override
    protected Exception doWork( String[] args ) {
        Exception exc = processCommandLine( "CorrelationAnalysis", args );
        if ( exc != null ) {
            return exc;
        }

        Collection<Gene> queryGenes, targetGenes;
        try {
            queryGenes = getQueryGenes();
            targetGenes = getTargetGenes();
        } catch ( IOException e ) {
            return e;
        }

        // calculate matrices
        CoexpressionMatrices matrices = coexpressionAnalysisService.calculateCoexpressionMatrices(
                expressionExperiments, queryGenes, targetGenes, filterConfig, CorrelationMethod.SPEARMAN );
        DenseDouble3dMatrix<Gene, Gene, BioAssaySet> correlationMatrix = matrices.getCorrelationMatrix();
        // DenseDoubleMatrix3DNamed sampleSizeMatrix = matrices
        // .getSampleSizeMatrix();

        // DoubleMatrixNamed maxCorrelationMatrix = coexpressionAnalysisService
        // .getMaxCorrelationMatrix(correlationMatrix, kMax);
        // DoubleMatrixNamed pValMatrix = coexpressionAnalysisService
        // .calculateMaxCorrelationPValueMatrix(maxCorrelationMatrix,
        // kMax, ees);
        // DoubleMatrixNamed effectSizeMatrix = coexpressionAnalysisService
        // .calculateEffectSizeMatrix(correlationMatrix, sampleSizeMatrix);

        // get row/col name maps
        Map<Gene, String> geneNameMap = matrices.getGeneNameMap();
        Map<ExpressionExperiment, String> eeNameMap = matrices.getEeNameMap();

        DecimalFormat formatter = ( DecimalFormat ) NumberFormat.getNumberInstance( Locale.US );
        formatter.applyPattern( "0.0000" );
        DecimalFormatSymbols symbols = formatter.getDecimalFormatSymbols();
        symbols.setNaN( "NaN" );
        formatter.setDecimalFormatSymbols( symbols );

        try {
            MatrixWriter<Gene, Gene> matrixOut;
            matrixOut = new MatrixWriter<Gene, Gene>( outFilePrefix + ".corr.txt", formatter );
            matrixOut.setSliceNameMap( eeNameMap );
            matrixOut.setRowNameMap( geneNameMap );
            matrixOut.setColNameMap( geneNameMap );
            matrixOut.writeMatrix( correlationMatrix, false );

            PrintWriter out = new PrintWriter( new FileWriter( outFilePrefix + ".corr.row_names.txt" ) );
            List<Gene> rows = correlationMatrix.getRowNames();
            for ( Gene row : rows ) {
                out.println( row );
            }
            out.close();

            out = new PrintWriter( new FileWriter( outFilePrefix + ".corr.col_names.txt" ) );
            Collection<BioAssaySet> cols = correlationMatrix.getSliceNames();
            for ( BioAssaySet bas : cols ) {
                ExpressionExperiment ee = ( ExpressionExperiment ) bas;
                out.println( ee.getShortName() );
            }
            out.close();

            // out = new MatrixWriter(outFilePrefix + ".max_corr.txt", formatter,
            // geneNameMap, geneNameMap);
            // out.writeMatrix(maxCorrelationMatrix, true);
            // out.close();
            //
            // out = new MatrixWriter(outFilePrefix + ".max_corr.pVal.txt",
            // formatter, geneNameMap, geneNameMap);
            // out.writeMatrix(pValMatrix, true);
            // out.close();
            //
            // out = new MatrixWriter(outFilePrefix + ".effect_size.txt",
            // formatter, geneNameMap, geneNameMap);
            // out.writeMatrix(effectSizeMatrix, true);
            // out.close();

        } catch ( IOException e ) {
            return e;
        }

        return null;
    }

    public static void main( String[] args ) {
        CorrelationAnalysisCLI analysis = new CorrelationAnalysisCLI();
        StopWatch watch = new StopWatch();
        watch.start();
        log.info( "Starting Correlation Analysis" );
        Exception exc = analysis.doWork( args );
        if ( exc != null ) {
            log.error( exc.getMessage() );
        }
        log.info( "Finished analysis in " + watch );
    }
}
