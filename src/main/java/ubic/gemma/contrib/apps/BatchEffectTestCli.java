/*
 * The Gemma project
 *
 * Copyright (c) 2011 University of British Columbia
 *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with
 * the License. You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on
 * an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the
 * specific language governing permissions and limitations under the License.
 */
package ubic.gemma.contrib.apps;

import org.springframework.beans.factory.annotation.Autowired;
import ubic.gemma.core.analysis.preprocess.batcheffects.BatchConfound;
import ubic.gemma.core.analysis.preprocess.batcheffects.BatchConfoundUtils;
import ubic.gemma.core.analysis.preprocess.svd.SVDResult;
import ubic.gemma.core.analysis.preprocess.svd.SVDService;
import ubic.gemma.core.apps.ExpressionExperimentManipulatingCLI;
import ubic.gemma.model.expression.experiment.BioAssaySet;
import ubic.gemma.model.expression.experiment.ExperimentalFactor;
import ubic.gemma.model.expression.experiment.ExpressionExperiment;

import java.util.Collection;
import java.util.Map;

/**
 * For bulk processing of batch-info-fetching.
 *
 * @author paul
 * @version $Id: BatchEffectTestCli.java,v 1.6 2015/11/12 19:37:12 paul Exp $
 */
public class BatchEffectTestCli extends ExpressionExperimentManipulatingCLI {

    @Autowired
    private SVDService svdService;

    @Override
    public String getCommandName() {
        return "batchEffectTest";
    }

    @Override
    public String getShortDesc() {
        return "Test for batch effects";
    }

    @Override
    protected void doWork() {

        for ( BioAssaySet bas : expressionExperiments ) {
            if ( bas instanceof ExpressionExperiment ) {
                ExpressionExperiment ee = ( ExpressionExperiment ) bas;
                ee = eeService.thawLite( ee );
                log.info( "Processing: " + ee );

                try {

                    boolean success = false;

                    pcaFactorTest( svdService, ee );

                    Collection<BatchConfound> results = BatchConfoundUtils.test( ee );

                    for ( BatchConfound r : results ) {
                        System.out.println( r );
                    }

                    success = true;
                    if ( success ) {
                        this.addSuccessObject( bas.toString(), "" );
                    } else {
                        this.addErrorObject( bas.toString() + ": No dates found", "" );

                    }

                } catch ( Exception e ) {
                    log.error( e, e );
                    this.addErrorObject( bas + ": " + e.getMessage(), "" );
                }

            }
        }

    }

    /**
     * Just extracts information from the SVD and prints it out.
     */
    private void pcaFactorTest( SVDService svdService, ExpressionExperiment ee ) {
        SVDResult svdo = svdService.getSvdFactorAnalysis( ee.getId() );
        /*
         * Compare PCs to batches.
         */
        Map<Integer, Double> dateCorrelations = svdo.getDateCorrelations();
        Map<Integer, Double> datePvals = svdo.getDatePVals();
        Map<Integer, Map<ExperimentalFactor, Double>> factorCorrelations = svdo.getFactorCorrelations();
        Map<Integer, Map<ExperimentalFactor, Double>> factorPvals = svdo.getFactorPVals();

        for ( Integer cmp : dateCorrelations.keySet() ) {
            Double dateCorr = dateCorrelations.get( cmp );
            Double datePval = datePvals.get( cmp );
            System.out.println( "PCA\t" + ee.getId() + "\t" + ee.getShortName() + "\tRunDate\tRunDate\tPC" + cmp + "\t"
                    + String.format( "%.2f", dateCorr ) + "\t" + String.format( "%.2g", datePval ) );
        }

        for ( Integer cmp : factorCorrelations.keySet() ) {
            for ( ExperimentalFactor ef : factorCorrelations.get( cmp ).keySet() ) {
                Double factorCorr = factorCorrelations.get( cmp ).get( ef );
                Double factorPval = factorPvals.get( cmp ).get( ef );
                System.out.println( "PCA\t" + ee.getId() + "\t" + ee.getShortName() + "\t" + ef + "\t"
                        + ef.getName() + "\tPC" + cmp + "\t" + String.format( "%.2f", factorCorr )
                        + "\t" + String.format( "%.2g", factorPval ) );
            }
        }
    }

}
