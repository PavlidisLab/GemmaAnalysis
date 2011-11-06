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

import cern.colt.list.DoubleArrayList;

import ubic.basecode.math.DescriptiveWithMissing;
import ubic.gemma.apps.DifferentialExpressionAnalysisCli;
import ubic.gemma.model.common.quantitationtype.QuantitationType;
import ubic.gemma.model.common.quantitationtype.ScaleType;
import ubic.gemma.model.expression.bioAssayData.DoubleVectorValueObject;
import ubic.gemma.model.expression.bioAssayData.ProcessedExpressionDataVectorService;
import ubic.gemma.model.expression.experiment.BioAssaySet;
import ubic.gemma.model.expression.experiment.ExpressionExperiment;

/**
 * Performs multiple differential expression analyses under different conditions: Without including a batch covariate;
 * with including it; and repeating those, after batch correction
 * 
 * @author paul
 * @version $Id$
 */
public class ScalingCheckCli extends DifferentialExpressionAnalysisCli {

    /**
     * @param args
     */
    public static void main( String[] args ) {
        ScalingCheckCli c = new ScalingCheckCli();
        c.doWork( args );

    }

    ProcessedExpressionDataVectorService processedExpressionDataVectorService;

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
        Exception error = super.processCommandLine( "batch diff ex test", args );
        if ( error != null ) return error;
        this.processedExpressionDataVectorService = ( ProcessedExpressionDataVectorService ) this
                .getBean( "processedExpressionDataVectorService" );

        try {
            summaryFile = initOutputFile( "scale.info.txt" );
            summaryFile.write( "State\tEEID\tEENAME\tMIN\tMED\tMAX\n" );

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

        return null;
    }

    private boolean onLogScale( Collection<DoubleVectorValueObject> vectos ) throws IOException {

        if ( vectos.isEmpty() ) throw new IllegalArgumentException( "No vectors!" );

        Long id = null;
        String shortName = null;
        QuantitationType quantitationType = null;
        DoubleArrayList r = new DoubleArrayList();

        for ( DoubleVectorValueObject v : vectos ) {

            quantitationType = v.getQuantitationType();

            // take it at its word.
            if ( quantitationType.getScale() != null ) {
                id = v.getExpressionExperiment().getId();
                shortName = v.getExpressionExperiment().getShortName();
                if ( quantitationType.getScale().equals( ScaleType.LOG2 ) ) {
                    summaryFile.write( "SUPPOSEDLY LOG2 SCALED\t" + id + "\t" + shortName + "\n" );
                } else if ( quantitationType.getScale().equals( ScaleType.LOG10 ) ) {
                    summaryFile.write( "SUPPOSEDLY LOG10 SCALED\t" + id + "\t" + shortName + "\n" );

                } else if ( quantitationType.getScale().equals( ScaleType.LOGBASEUNKNOWN ) ) {
                    summaryFile.write( "SUPPOSEDLY LOGBASEUNKNOWN SCALED\t" + id + "\t" + shortName + "\n" );
                }
                return true;
            }

            for ( int j = 0; j < v.getData().length; j++ ) {
                double a = v.getData()[j];
                r.add( a );
            }

        }

        double max = DescriptiveWithMissing.max( r );
        double min = DescriptiveWithMissing.min( r );
        double median = DescriptiveWithMissing.median( r );

        if ( quantitationType == null ) {
            throw new IllegalStateException( "QT was null" );
        }

        String m = "\t" + quantitationType.getName() + "\t" + quantitationType.getDescription() + "\t"
                + String.format( "%.2f", min ) + "\t" + String.format( "%.2f", median ) + "\t"
                + String.format( "%.2f", max );

        if ( max - min < 10 ) {
            summaryFile.write( "POSSIBLY LOG SCALED\t" + id + "\t" + shortName + m + "\n" );
            log.info( "Range is narrow, could be log scaled" );
            return true;
        }

        if ( max > 50 ) {
            log.info( "Data has large values, doesn't look log transformed: " + max );
            summaryFile.write( "NOT LOG SCALED\t" + id + "\t" + shortName + m + "\n" );
            return false;
        }
        if ( min < 1 ) {
            log.info( "Data has very small values, doesn't look log transformed: " + min );
            summaryFile.write( "NOT LOG SCALED\t" + id + "\t" + shortName + m + "\n" );
            return false;
        }

        log.info( "Can't rule out possibility of log scale" );
        summaryFile.write( "POSSIBLY LOG SCALED\t" + id + "\t" + shortName + m + "\n" );
        return true;

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

            /*
             * Extract data
             */
            Collection<DoubleVectorValueObject> vectos = processedExpressionDataVectorService.getProcessedDataArrays(
                    ee, 100, false );

            onLogScale( vectos );

        } catch ( Exception e ) {
            log.error( e, e );
        } finally {
        }
    }
}
