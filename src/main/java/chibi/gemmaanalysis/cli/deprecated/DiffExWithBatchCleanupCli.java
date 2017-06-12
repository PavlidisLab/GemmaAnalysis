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

package chibi.gemmaanalysis.cli.deprecated;

import java.util.Collection;

import ubic.gemma.core.analysis.expression.diff.DifferentialExpressionAnalyzerService;
import ubic.gemma.core.analysis.report.ExpressionExperimentReportService;
import ubic.gemma.core.apps.ExpressionExperimentManipulatingCLI;
import ubic.gemma.model.analysis.expression.diff.DifferentialExpressionAnalysis;
import ubic.gemma.model.analysis.expression.diff.ExpressionAnalysisResultSet;
import ubic.gemma.model.expression.experiment.BioAssaySet;
import ubic.gemma.model.expression.experiment.ExperimentalFactor;
import ubic.gemma.model.expression.experiment.ExpressionExperiment;
import ubic.gemma.persistence.service.analysis.expression.diff.DifferentialExpressionAnalysisService;
import ubic.gemma.persistence.service.expression.experiment.ExperimentalFactorService;

/**
 * One-off; dDelete analyses which include 'batch' as a factor.
 *
 * @author paul
 * @version $Id: DiffExWithBatchCleanupCli.java,v 1.15 2015/11/12 19:37:11 paul Exp $
 */
public class DiffExWithBatchCleanupCli extends ExpressionExperimentManipulatingCLI {

    public static void main( String[] args ) {
        DiffExWithBatchCleanupCli c = new DiffExWithBatchCleanupCli();
        c.doWork( args );
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

    /*
     * (non-Javadoc)
     *
     * @see ubic.gemma.util.AbstractCLI#doWork(java.lang.String[])
     */
    @Override
    protected Exception doWork( String[] args ) {

        processCommandLine( args );

        DifferentialExpressionAnalysisService differentialExpressionAnalysisService = this
                .getBean( DifferentialExpressionAnalysisService.class );

        DifferentialExpressionAnalyzerService differentialExpressionAnalyzerService = this
                .getBean( DifferentialExpressionAnalyzerService.class );

        ExpressionExperimentReportService expressionExperimentReportService = this
                .getBean( ExpressionExperimentReportService.class );
        for ( BioAssaySet ee : expressionExperiments ) {
            if ( !( ee instanceof ExpressionExperiment ) ) {
                continue;
            }

            ExpressionExperiment expressionExperiment = ( ExpressionExperiment ) ee;

            log.info( "Processing: " + ee );

            try {

                Collection<DifferentialExpressionAnalysis> diffAnalyses = differentialExpressionAnalysisService
                        .findByInvestigation( expressionExperiment );

                if ( diffAnalyses.isEmpty() ) continue;

                differentialExpressionAnalysisService.thaw( diffAnalyses );

                for ( DifferentialExpressionAnalysis existingAnalysis : diffAnalyses ) {

                    for ( ExpressionAnalysisResultSet resultSet : existingAnalysis.getResultSets() ) {
                        if ( resultSet.getExperimentalFactors().size() > 1 ) {
                            continue;
                        }
                        ExperimentalFactor factor = resultSet.getExperimentalFactors().iterator().next();
                        if ( factor.getName().equals( ExperimentalFactorService.BATCH_FACTOR_NAME ) ) {
                            log.info( "Deleting analysis with batch factor, Id=" + existingAnalysis.getId() );
                            differentialExpressionAnalyzerService.deleteAnalysis( expressionExperiment,
                                    existingAnalysis );
                            expressionExperimentReportService.generateSummary( expressionExperiment.getId() );
                            this.successObjects.add( ee );
                        }
                    }
                }

            } catch ( Exception e ) {
                this.errorObjects.add( ee + ": " + e.getMessage() );
                log.error( e, e );
            }
        }
        super.summarizeProcessing();
        return null;
    }

}
