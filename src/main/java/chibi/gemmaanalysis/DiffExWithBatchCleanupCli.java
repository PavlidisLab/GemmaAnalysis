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

import ubic.gemma.analysis.expression.diff.DifferentialExpressionAnalyzerService;
import ubic.gemma.apps.ExpressionExperimentManipulatingCLI;
import ubic.gemma.model.analysis.expression.diff.DifferentialExpressionAnalysis;
import ubic.gemma.model.analysis.expression.diff.DifferentialExpressionAnalysisService;
import ubic.gemma.model.analysis.expression.diff.ExpressionAnalysisResultSet;
import ubic.gemma.model.expression.experiment.BioAssaySet;
import ubic.gemma.model.expression.experiment.ExperimentalFactor;
import ubic.gemma.model.expression.experiment.ExpressionExperiment;

/**
 * Delete analyses which include 'batch' as a factor.
 * 
 * @author paul
 * @version $Id$
 */
public class DiffExWithBatchCleanupCli extends ExpressionExperimentManipulatingCLI {

    /*
     * (non-Javadoc)
     * 
     * @see ubic.gemma.util.AbstractCLI#doWork(java.lang.String[])
     */
    @Override
    protected Exception doWork( String[] args ) {
        DifferentialExpressionAnalyzerService ds = ( DifferentialExpressionAnalyzerService ) this
                .getBean( "differentialExpressionAnalyzerService" );
        DifferentialExpressionAnalysisService differentialExpressionAnalysisService = ( DifferentialExpressionAnalysisService ) this
                .getBean( "differentialExpressionAnalysisService" );

        for ( BioAssaySet ee : expressionExperiments ) {
            if ( !( ee instanceof ExpressionExperiment ) ) {
                continue;
            }

            ExpressionExperiment expressionExperiment = ( ExpressionExperiment ) ee;

            log.info( "Processing: " + ee );

            Collection<DifferentialExpressionAnalysis> diffAnalyses = differentialExpressionAnalysisService
                    .findByInvestigation( expressionExperiment );

            differentialExpressionAnalysisService.thaw( diffAnalyses );

            a: for ( DifferentialExpressionAnalysis existingAnalysis : diffAnalyses ) {

                for ( ExpressionAnalysisResultSet resultSet : existingAnalysis.getResultSets() ) {
                    if ( resultSet.getExperimentalFactors().size() > 1 ) {
                        continue a;
                    }
                    ExperimentalFactor factor = resultSet.getExperimentalFactors().iterator().next();
                    if ( factor.getName().equals( "batch" ) ) {
                        log.info( "Found analysis with batch factor, Id=" + existingAnalysis.getId() );
                        ds.deleteOldAnalysis( expressionExperiment, existingAnalysis );

                    }
                }

            }

        }

        return null;
    }

}
