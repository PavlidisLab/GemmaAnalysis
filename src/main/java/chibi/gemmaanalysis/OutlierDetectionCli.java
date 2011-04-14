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

import ubic.gemma.analysis.preprocess.OutlierDetectionService;
import ubic.gemma.apps.DifferentialExpressionAnalysisCli;
import ubic.gemma.model.expression.bioAssay.BioAssay;
import ubic.gemma.model.expression.experiment.BioAssaySet;
import ubic.gemma.model.expression.experiment.ExpressionExperiment;

/**
 * Checks experiments for outliers
 * 
 * @author paul
 * @version $Id$
 */
public class OutlierDetectionCli extends DifferentialExpressionAnalysisCli {

    /**
     * @param args
     */
    public static void main( String[] args ) {
        OutlierDetectionCli c = new OutlierDetectionCli();
        c.doWork( args );

    }

    OutlierDetectionService outlierDetectionService;

    Writer summaryFile;

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

    /*
     * (non-Javadoc)
     * 
     * @see ubic.gemma.util.AbstractCLI#doWork(java.lang.String[])
     */
    @Override
    protected Exception doWork( String[] args ) {
        Exception error = super.processCommandLine( "outlier detection", args );
        if ( error != null ) return error;

        this.outlierDetectionService = ( OutlierDetectionService ) this.getBean( "outlierDetectionService" );

        try {
            summaryFile = initOutputFile( "outlier.summary.txt" );
            summaryFile.write( "EEID\tEENAME\tBAID\tBANAME\n" );

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

    /*
     * (non-Javadoc)
     * 
     * @seeubic.gemma.apps.DifferentialExpressionAnalysisCli#processExperiment(ubic.gemma.model.expression.experiment.
     * ExpressionExperiment)
     */
    @Override
    protected void processExperiment( ExpressionExperiment ee ) {
        try {
            StringBuilder sb = new StringBuilder();

            Collection<BioAssay> outliers = this.outlierDetectionService.identifyOutliers( ee );

            String baseS = ee.getId() + "\t" + ee.getShortName();

            // TODO we might get more details here (like the correlation)
            for ( BioAssay bioAssay : outliers ) {
                String line = baseS + "\t" + bioAssay.getId() + "\t" + bioAssay.getName();
                sb.append( line + "\n" );
                log.info( line );
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

}
