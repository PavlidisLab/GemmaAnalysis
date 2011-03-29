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

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

import ubic.gemma.analysis.expression.diff.DifferentialExpressionAnalysisConfig;
import ubic.gemma.analysis.expression.diff.GenericAncovaAnalyzer;
import ubic.gemma.analysis.expression.diff.LinearModelAnalyzer;
import ubic.gemma.analysis.preprocess.batcheffects.BatchInfoPopulationService;
import ubic.gemma.analysis.preprocess.batcheffects.ExpressionExperimentBatchCorrectionService;
import ubic.gemma.apps.DifferentialExpressionAnalysisCli;
import ubic.gemma.datastructure.matrix.ExpressionDataDoubleMatrix;
import ubic.gemma.model.analysis.expression.diff.DifferentialExpressionAnalysis;
import ubic.gemma.model.analysis.expression.diff.DifferentialExpressionAnalysisResult;
import ubic.gemma.model.analysis.expression.diff.ExpressionAnalysisResultSet;
import ubic.gemma.model.expression.bioAssayData.ProcessedExpressionDataVector;
import ubic.gemma.model.expression.bioAssayData.ProcessedExpressionDataVectorService;
import ubic.gemma.model.expression.experiment.BioAssaySet;
import ubic.gemma.model.expression.experiment.ExperimentalFactor;
import ubic.gemma.model.expression.experiment.ExpressionExperiment;

/**
 * TODO Document Me
 * 
 * @author paul
 * @version $Id$
 */
public class BatchDiffExCli extends DifferentialExpressionAnalysisCli {
    ExpressionExperimentBatchCorrectionService expressionExperimentBatchCorrectionService;

    LinearModelAnalyzer lma;

    ProcessedExpressionDataVectorService processedExpressionDataVectorService;

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

        for ( BioAssaySet bas : this.expressionExperiments ) {
            if ( !( bas instanceof ExpressionExperiment ) ) {
                continue;
            }
            ExpressionExperiment ee = eeService.thawLite( ( ExpressionExperiment ) bas );
            log.info( "Processing: " + ee );

            /*
             * Check if it has a batch factor.
             */
            ExperimentalFactor batchFactor = expressionExperimentBatchCorrectionService.getBatchFactor( ee );
            if ( null == batchFactor ) {
                // this.errorObjects.add( "No batch factor: " + ee.getShortName() );
                continue;
            }

            /* TODO use this */
            expressionExperimentBatchCorrectionService.checkBatchEffectSeverity( ee );
            expressionExperimentBatchCorrectionService.checkCorrectability( ee );

            /*
             * Extract data
             */
            Collection<ProcessedExpressionDataVector> vectos = processedExpressionDataVectorService
                    .getProcessedDataVectors( ee );
            ExpressionDataDoubleMatrix mat = new ExpressionDataDoubleMatrix( vectos );

            /*
             * first do an analysis without batch. Let's ignore interactions to keep things simple.
             */
            Collection<ExperimentalFactor> factors2 = new HashSet<ExperimentalFactor>();
            for ( ExperimentalFactor ef : ee.getExperimentalDesign().getExperimentalFactors() ) {
                if ( ef.equals( batchFactor ) ) continue;
                factors2.add( ef );
            }
            DifferentialExpressionAnalysisConfig config2 = new DifferentialExpressionAnalysisConfig();
            config2.setFactorsToInclude( factors2 );
            DifferentialExpressionAnalysis beforeResults = lma.run( ee, mat, config2 ).iterator().next();
            Map beforeResultCount = new HashMap();
            for ( ExpressionAnalysisResultSet brs : beforeResults.getResultSets() ) {
                assert brs.getExperimentalFactors().size() == 1;
                ExperimentalFactor ef = brs.getExperimentalFactors().iterator().next();
                Collection<DifferentialExpressionAnalysisResult> results = brs.getResults();
                int c = 0;
                for ( DifferentialExpressionAnalysisResult r : results ) {
                    Double pval = r.getCorrectedPvalue();
                    if ( pval < 0.01 ) {
                        c++;
                    }
                }
                System.out.println( "Before\t" + ef.getId() + "\t" + ef.getName() + "\t" + c );
                beforeResultCount.put( ef, c );
            }

            /*
             * Then do it with batch.
             */
            Collection<ExperimentalFactor> factors = ee.getExperimentalDesign().getExperimentalFactors();
            assert factors.contains( batchFactor );
            DifferentialExpressionAnalysisConfig configIncludingBatch = new DifferentialExpressionAnalysisConfig();
            configIncludingBatch.setFactorsToInclude( factors );
            DifferentialExpressionAnalysis withBatchEffectResults = lma.run( ee, mat, configIncludingBatch ).iterator().next();

            /*
             * Determine how many genes are diff ex wrt batch.
             */
            int batchEffectCount = 0;
            for ( ExpressionAnalysisResultSet brs : withBatchEffectResults.getResultSets() ) {
                assert brs.getExperimentalFactors().size() == 1;
                ExperimentalFactor ef = brs.getExperimentalFactors().iterator().next();
                Collection<DifferentialExpressionAnalysisResult> results = brs.getResults();
                if ( ef.equals( batchFactor ) ) {
                    for ( DifferentialExpressionAnalysisResult r : results ) {
                        Double pval = r.getCorrectedPvalue();
                        if ( pval < 0.01 ) {
                            batchEffectCount++;
                        }
                    }
                    System.out.println( "Batch\t" + ef.getId() + "\t" + ef.getName() + "\t" + batchEffectCount );
                    break;
                }
                log.info( "Factor: " + ef );
            }

            /*
             * Correct for batch effects and repeat the analysis
             */
            ExpressionDataDoubleMatrix comBat = expressionExperimentBatchCorrectionService.comBat( ee, mat );
            assert comBat != null;
            DifferentialExpressionAnalysis revisedResult = lma.run( ee, comBat, config2 ).iterator().next();

            /*
             * Enumer # of probes for each factor ...
             */
            Map revisedResultCount = new HashMap();
            for ( ExpressionAnalysisResultSet brs : revisedResult.getResultSets() ) {
                assert brs.getExperimentalFactors().size() == 1;
                ExperimentalFactor ef = brs.getExperimentalFactors().iterator().next();
                Collection<DifferentialExpressionAnalysisResult> results = brs.getResults();
                int c = 0;
                for ( DifferentialExpressionAnalysisResult r : results ) {
                    Double pval = r.getCorrectedPvalue();
                    if ( pval < 0.01 ) {
                        c++;
                    }
                }
                System.out.println( "After\t" + ef.getId() + "\t" + ef.getName() + "\t" + c );
                revisedResultCount.put( ef, c );
            }

            /*
             * Print summary of results
             */

            successObjects.add( ee );
        }

        return null;
    }

    /**
     * @param args
     */
    public static void main( String[] args ) {
        BatchDiffExCli c = new BatchDiffExCli();
        c.doWork( args );

    }

}
