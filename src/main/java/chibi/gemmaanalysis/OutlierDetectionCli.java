/*
 * The GemmaAnalysis project
 * 
 * Copyright (c) 2011 University of British Columbia
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

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.Collection;
import java.util.Date;

import ubic.gemma.analysis.preprocess.OutlierDetails;
import ubic.gemma.analysis.preprocess.OutlierDetectionService;
import ubic.gemma.apps.ExpressionExperimentManipulatingCLI;
import ubic.gemma.model.expression.bioAssay.BioAssay;
import ubic.gemma.model.expression.experiment.BioAssaySet;
import ubic.gemma.model.expression.experiment.ExpressionExperiment;
import ubic.gemma.util.DateUtil;

/**
 * Checks experiments for outliers
 * 
 * @author paul
 * @version $Id$
 */
public class OutlierDetectionCli extends ExpressionExperimentManipulatingCLI {

    private static final String REGRESSION_OPT = "regression";

    /**
     * @param args
     */
    public static void main( String[] args ) {
        OutlierDetectionCli c = new OutlierDetectionCli();
        c.doWork( args );

    }

    OutlierDetectionService outlierDetectionService;

    Writer summaryFile;

    boolean useRegression = false;

    double thresholdFraction = 0.9;

    int quantile = 15;

    /**
     * @param fileName
     * @return
     * @throws IOException
     */
    private Writer initOutputFile( String fileName ) throws IOException {
        File f = new File( fileName );
        if ( f.exists() ) {
            f.delete();
        }
        f.createNewFile();
        log.info( "New file: " + f.getAbsolutePath() );
        return new FileWriter( f );
    }

    @Override
    protected void buildOptions() {
        super.buildOptions();

        this.addOption( REGRESSION_OPT, false, "Set to use regression to attempt to "
                + "subtract the effects of experimental factors before looking for outliers." );

        this.addOption( "quantile", true,
                "Set the quantile of correlations that would be considered 'suspiciously low'" );

        this.addOption( "thresholdFraction", true, "Set the fraction of samples which must have "
                + "correlations < quantile for a sample to be considered an outlier." );
    }

    /*
     * (non-Javadoc)
     * 
     * @see ubic.gemma.util.AbstractCLI#doWork(java.lang.String[])
     */
    @Override
    protected Exception doWork( String[] args ) {
        Exception error = super.processCommandLine( "outlier detection", args );
        if ( error != null ) return error;

        this.outlierDetectionService = this.getBean( OutlierDetectionService.class );

        try {
            summaryFile = initOutputFile( "outlier.summary.txt" );
            summaryFile.write( "# Created by Gemma OutlierDetectionCli\n" );
            summaryFile.write( "# " + DateUtil.convertDateToString( new Date() ) + "\n" );
            summaryFile.write( "# Quantile = " + quantile + "\n" );
            summaryFile.write( "# Threshold fraction = " + String.format( "%.2f", this.thresholdFraction ) + "\n" );
            summaryFile.write( "# Regression = " + this.useRegression + "\n" );
            summaryFile.write( "EEID\tEENAME\tBAID\tBANAME\tSCORE\tTHRESHOLDCORR\n" );

            for ( BioAssaySet bas : this.expressionExperiments ) {
                if ( !( bas instanceof ExpressionExperiment ) ) {
                    continue;
                }
                ExpressionExperiment ee = eeService.thawLite( ( ExpressionExperiment ) bas );
                processExperiment( ee );
            }
            summaryFile.close();

        } catch ( Exception e ) {
            log.error( e, e );
            return e;
        }

        summarizeProcessing();

        return null;
    }

    /**
     * @param ee
     */
    protected void processExperiment( ExpressionExperiment ee ) {
        try {
            StringBuilder sb = new StringBuilder();

            Collection<OutlierDetails> outliers = this.outlierDetectionService.identifyOutliers( ee,
                    this.useRegression, this.quantile, this.thresholdFraction );

            String baseS = ee.getId() + "\t" + ee.getShortName();

            for ( OutlierDetails outlier : outliers ) {
                BioAssay bioAssay = outlier.getBioAssay();
                double score = outlier.getScore();
                String line = baseS + "\t" + bioAssay.getId() + "\t" + bioAssay.getName() + "\t"
                        + String.format( "%.2f\t%.2f", score, outlier.getThresholdCorrelation() );
                sb.append( line + "\n" );
                log.info( "Outlier: " + line );
            }

            summaryFile.write( sb.toString() );
            summaryFile.flush();

            successObjects.add( ee );
        } catch ( Exception e ) {
            log.error( e, e );
            errorObjects.add( ee + e.getMessage() );
        } finally {
        }
    }

    @Override
    protected void processOptions() {
        super.processOptions();
        if ( this.hasOption( REGRESSION_OPT ) ) {
            this.useRegression = true;
        }
    }

}
