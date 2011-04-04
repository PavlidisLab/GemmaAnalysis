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
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

import ubic.gemma.analysis.expression.diff.DifferentialExpressionAnalysisConfig;
import ubic.gemma.analysis.expression.diff.LinearModelAnalyzer;
import ubic.gemma.analysis.preprocess.batcheffects.ExpressionExperimentBatchCorrectionService;
import ubic.gemma.apps.DifferentialExpressionAnalysisCli;
import ubic.gemma.datastructure.matrix.ExpressionDataDoubleMatrix;
import ubic.gemma.model.analysis.expression.diff.DifferentialExpressionAnalysis;
import ubic.gemma.model.analysis.expression.diff.DifferentialExpressionAnalysisResult;
import ubic.gemma.model.analysis.expression.diff.ExpressionAnalysisResultSet;
import ubic.gemma.model.analysis.expression.diff.ProbeAnalysisResult;
import ubic.gemma.model.expression.bioAssayData.ProcessedExpressionDataVector;
import ubic.gemma.model.expression.bioAssayData.ProcessedExpressionDataVectorService;
import ubic.gemma.model.expression.designElement.CompositeSequence;
import ubic.gemma.model.expression.experiment.BioAssaySet;
import ubic.gemma.model.expression.experiment.ExperimentalFactor;
import ubic.gemma.model.expression.experiment.ExpressionExperiment;

/**
 * Performs multiple differential expression analyses under different conditions: Without including a batch covariate;
 * with including it; and repeating those, after batch correction
 * 
 * @author paul
 * @version $Id$
 */
public class BatchDiffExCli extends DifferentialExpressionAnalysisCli {
    private static final int LOGGING_FREQ = 20000;

    /**
     * @param args
     */
    public static void main( String[] args ) {
        BatchDiffExCli c = new BatchDiffExCli();
        c.doWork( args );

    }

    ExpressionExperimentBatchCorrectionService expressionExperimentBatchCorrectionService;

    LinearModelAnalyzer lma;

    ProcessedExpressionDataVectorService processedExpressionDataVectorService;

    /**
     * @param revisedResultDetails
     * @param ef
     * @param r
     * @param c
     * @return c
     */
    private int tally( Map<CompositeSequence, Map<ExperimentalFactor, Double>> revisedResultDetails,
            ExperimentalFactor ef, DifferentialExpressionAnalysisResult r, int c ) {
        Double pval = r.getCorrectedPvalue();
        if ( pval != null && pval < 0.01 ) {
            c++;
        }
        /*
         * Map of probe -> factor -> pval
         */
        CompositeSequence probe = ( ( ProbeAnalysisResult ) r ).getProbe();
        if ( !revisedResultDetails.containsKey( probe ) ) {
            revisedResultDetails.put( probe, new HashMap<ExperimentalFactor, Double>() );
        }
        revisedResultDetails.get( probe ).put( ef, pval );
        return c;
    }

    Writer summaryFile;

    /*
     * (non-Javadoc)
     * 
     * @see ubic.gemma.util.AbstractCLI#doWork(java.lang.String[])
     */
    @Override
    protected Exception doWork( String[] args ) {
        Exception error = super.processCommandLine( "batch diff ex test", args );
        if ( error != null ) return error;

        this.expressionExperimentBatchCorrectionService = ( ExpressionExperimentBatchCorrectionService ) this
                .getBean( "expressionExperimentBatchCorrectionService" );
        this.lma = ( LinearModelAnalyzer ) this.getBean( "genericAncovaAnalyzer" );
        this.processedExpressionDataVectorService = ( ProcessedExpressionDataVectorService ) this
                .getBean( "processedExpressionDataVectorService" );

        try {
            summaryFile = initOutputFile( "batch.proc.summary.txt" );
            summaryFile.write( "State\tEEID\tEENAME\tEFID\tEFNAME\tNUM\tNUMDIFF\n" );

            for ( BioAssaySet bas : this.expressionExperiments ) {
                if ( !( bas instanceof ExpressionExperiment ) ) {
                    continue;
                }
                ExpressionExperiment ee = eeService.thawLite( ( ExpressionExperiment ) bas );
                log.info( "Processing: " + ee );

                processExperiment( ee );
            }
            summaryFile.close();

        } catch ( Exception e ) {
            log.error( e, e );
            return e;
        }
        return null;
    }

    @Override
    protected void processExperiment( ExpressionExperiment ee ) {
        Writer detailFile = null;
        try {

            Collection<ExperimentalFactor> experimentalFactors = ee.getExperimentalDesign().getExperimentalFactors();

            ExperimentalFactor batchFactor = expressionExperimentBatchCorrectionService.getBatchFactor( ee );
            if ( null == batchFactor ) {
                this.errorObjects.add( "No batch factor: " + ee.getShortName() );
                return;
            }

            if ( experimentalFactors.size() < 2 ) {
                // need at least two factors, one of which has to be the batch.
                this.errorObjects.add( "Too few factors: " + ee.getShortName() );
                return;
            }

            if ( ee.getBioAssays().size() < 8 ) {
                this.errorObjects.add( "Too small (" + ee.getBioAssays().size() + " samples): " + ee.getShortName() );
                return;
            }

            if ( experimentalFactors.size() > 4 ) {
                /*
                 * This could be modified to select just a few factors, at random ... but that's probably
                 */
                this.errorObjects.add( "Too many factors (" + experimentalFactors.size() + " factors): "
                        + ee.getShortName() );
                return;
            }

            /* TODO use this, or skip it... we have this information elsewhere already */
            expressionExperimentBatchCorrectionService.checkBatchEffectSeverity( ee );
            expressionExperimentBatchCorrectionService.checkCorrectability( ee );

            /*
             * Extract data
             */
            Collection<ProcessedExpressionDataVector> vectos = processedExpressionDataVectorService
                    .getProcessedDataVectors( ee );
            ExpressionDataDoubleMatrix mat = new ExpressionDataDoubleMatrix( vectos );

            /*
             * first do an analysis without batch; this is our baseline. Let's ignore interactions to keep things
             * simple.
             */
            Collection<ExperimentalFactor> factors2 = new HashSet<ExperimentalFactor>();
            for ( ExperimentalFactor ef : experimentalFactors ) {
                if ( ef.equals( batchFactor ) ) continue;
                factors2.add( ef );
            }
            int j = 0;
            DifferentialExpressionAnalysisConfig configWithoutBatch = new DifferentialExpressionAnalysisConfig();
            configWithoutBatch.setFactorsToInclude( factors2 );
            DifferentialExpressionAnalysis beforeResults = lma.run( ee, mat, configWithoutBatch ).iterator().next();
            Map<CompositeSequence, Map<ExperimentalFactor, Double>> beforeResultDetails = new HashMap<CompositeSequence, Map<ExperimentalFactor, Double>>();
            for ( ExpressionAnalysisResultSet brs : beforeResults.getResultSets() ) {
                assert brs.getExperimentalFactors().size() == 1;
                ExperimentalFactor ef = brs.getExperimentalFactors().iterator().next();
                Collection<DifferentialExpressionAnalysisResult> results = brs.getResults();
                int c = 0;
                for ( DifferentialExpressionAnalysisResult r : results ) {
                    c = tally( beforeResultDetails, ef, r, c );
                    if ( ++j % LOGGING_FREQ == 0 ) {
                        log.info( j + " processed" );
                    }
                }
                summaryFile.write( "Before\t" + ee.getId() + "\t" + ee.getShortName() + "\t" + ef.getId() + "\t"
                        + ef.getName() + "\t" + results.size() + "\t" + c + "\n" );
            }

            /*
             * Then do it with batch.
             */
            Collection<ExperimentalFactor> factors = experimentalFactors;
            assert factors.contains( batchFactor );
            DifferentialExpressionAnalysisConfig configIncludingBatch = new DifferentialExpressionAnalysisConfig();
            configIncludingBatch.setFactorsToInclude( factors );
            DifferentialExpressionAnalysis withBatchEffectResults = lma.run( ee, mat, configIncludingBatch ).iterator()
                    .next();

            /*
             * Determine how many genes are diff ex wrt batch. The other factors are tracked; this shows how we would do
             * if we tried to simply directly include batch in the model
             */
            Map<CompositeSequence, Map<ExperimentalFactor, Double>> batchEffectDetails = new HashMap<CompositeSequence, Map<ExperimentalFactor, Double>>();

            for ( ExpressionAnalysisResultSet brs : withBatchEffectResults.getResultSets() ) {
                assert brs.getExperimentalFactors().size() == 1;
                ExperimentalFactor ef = brs.getExperimentalFactors().iterator().next();
                Collection<DifferentialExpressionAnalysisResult> results = brs.getResults();

                int c = 0;
                for ( DifferentialExpressionAnalysisResult r : results ) {
                    c = tally( batchEffectDetails, ef, r, c );

                    if ( ++j % LOGGING_FREQ == 0 ) {
                        log.info( j + " processed" );
                    }

                }
                summaryFile.write( "Batch\t" + ee.getId() + "\t" + ee.getShortName() + "\t" + ef.getId() + "\t"
                        + ef.getName() + "\t" + results.size() + "\t" + c + "\n" );
            }

            /*
             * Correct for batch effects
             */
            log.info( "ComBat-ing" );
            ExpressionDataDoubleMatrix comBat = expressionExperimentBatchCorrectionService.comBat( ee, mat );
            assert comBat != null;

            /*
             * Check if we have removed the batch effect: there should be no diff ex wrt batch. This is just a sanity
             * check, really. The other factors are tracked just for completeness.
             */
            DifferentialExpressionAnalysis revisedResultWithBatch = lma.run( ee, comBat, configIncludingBatch )
                    .iterator().next();

            Map<CompositeSequence, Map<ExperimentalFactor, Double>> batchEffectAfterCorrDetails = new HashMap<CompositeSequence, Map<ExperimentalFactor, Double>>();

            for ( ExpressionAnalysisResultSet brs : revisedResultWithBatch.getResultSets() ) {
                assert brs.getExperimentalFactors().size() == 1;
                ExperimentalFactor ef = brs.getExperimentalFactors().iterator().next();
                Collection<DifferentialExpressionAnalysisResult> results = brs.getResults();
                int c = 0;
                for ( DifferentialExpressionAnalysisResult r : results ) {
                    c = tally( batchEffectAfterCorrDetails, ef, r, c );

                    if ( ++j % LOGGING_FREQ == 0 ) {
                        log.info( j + " processed" );
                    }

                }
                summaryFile.write( "BatchAftCorr\t" + ee.getId() + "\t" + ee.getShortName() + "\t" + ef.getId() + "\t"
                        + ef.getName() + "\t" + results.size() + "\t" + c + "\n" );
            }

            /*
             * Now without batch as a factor, which is what we really want.
             */
            DifferentialExpressionAnalysis revisedResult = lma.run( ee, comBat, configWithoutBatch ).iterator().next();
            Map<CompositeSequence, Map<ExperimentalFactor, Double>> revisedResultDetails = new HashMap<CompositeSequence, Map<ExperimentalFactor, Double>>();
            for ( ExpressionAnalysisResultSet brs : revisedResult.getResultSets() ) {
                assert brs.getExperimentalFactors().size() == 1;
                ExperimentalFactor ef = brs.getExperimentalFactors().iterator().next();
                Collection<DifferentialExpressionAnalysisResult> results = brs.getResults();
                int c = 0;
                for ( DifferentialExpressionAnalysisResult r : results ) {
                    c = tally( revisedResultDetails, ef, r, c );

                    if ( ++j % LOGGING_FREQ == 0 ) {
                        log.info( j + " processed" );
                    }

                }
                summaryFile.write( "After\t" + ee.getId() + "\t" + ee.getShortName() + "\t" + ef.getId() + "\t"
                        + ef.getName() + "\t" + results.size() + "\t" + c + "\n" );

            }

            /*
             * Print out a summary
             */
            detailFile = initOutputFile( "batch.proc.detail." + ee.getId() + "."
                    + ee.getShortName().replaceAll( "[\\W\\s]+", "_" ) + ".txt" );

            detailFile
                    .write( "EEID\tEENAME\tEFID\tEFNAME\tPROBEID\tPROBENAME\tBEFOREQVAL\tBATCHQVAL\tBATAFTERQVAL\tAFTERQVAL\n" );

            for ( CompositeSequence c : beforeResultDetails.keySet() ) {
                for ( ExperimentalFactor ef : factors ) {
                    detailFile.write( ee.getId() + "\t" + ee.getShortName() + "\t" + ef.getId() + "\t" + ef.getName()
                            + "\t" + c.getId() + "\t" + c.getName() + "\t" );

                    Double bpval = beforeResultDetails.get( c ).get( ef ); // will be null for 'batch'

                    Double batpval = batchEffectDetails.get( c ).get( ef ); // when batch was included.

                    Double batapval = batchEffectAfterCorrDetails.get( c ).get( ef ); // when batch was included.

                    Double aftpval = revisedResultDetails.get( c ).get( ef ); // will be null for 'batch'

                    detailFile.write( String.format( "%.4g\t%.4g\t%.4g\t%.4g\n", bpval, batpval, batapval, aftpval ) );

                }
            }
            detailFile.close();
            summaryFile.flush();
            successObjects.add( ee );
        } catch ( Exception e ) {
            log.error( e, e );
            errorObjects.add( ee + e.getMessage() );
        } finally {
            if ( detailFile != null ) {
                try {
                    detailFile.close();
                } catch ( IOException e ) {
                    log.error( e, e );
                }
            }
        }
    }

    /**
     * @param fileName
     * @return
     * @throws IOException
     */
    protected Writer initOutputFile( String fileName ) throws IOException {
        File f = new File( fileName );
        if ( f.exists() ) {
            f.delete();
        }
        f.createNewFile();
        log.info( "New file: " + f.getAbsolutePath() );
        return new FileWriter( f );
    }

}