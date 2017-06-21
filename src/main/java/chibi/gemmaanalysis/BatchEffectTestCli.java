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
package chibi.gemmaanalysis;

import java.util.Collection;
import java.util.Map;

import ubic.gemma.core.analysis.preprocess.batcheffects.BatchConfound;
import ubic.gemma.core.analysis.preprocess.batcheffects.BatchConfoundValueObject;
import ubic.gemma.core.analysis.preprocess.svd.SVDService;
import ubic.gemma.core.analysis.preprocess.svd.SVDValueObject;
import ubic.gemma.core.apps.ExpressionExperimentManipulatingCLI;
import ubic.gemma.core.apps.GemmaCLI.CommandGroup;
import ubic.gemma.model.expression.experiment.BioAssaySet;
import ubic.gemma.model.expression.experiment.ExperimentalFactor;
import ubic.gemma.model.expression.experiment.ExpressionExperiment;
import ubic.gemma.persistence.util.EntityUtils;

/**
 * For bulk processing of batch-info-fetching.
 *
 * @author paul
 * @version $Id: BatchEffectTestCli.java,v 1.6 2015/11/12 19:37:12 paul Exp $
 */
public class BatchEffectTestCli extends ExpressionExperimentManipulatingCLI {

    public static void main( String[] args ) {
        BatchEffectTestCli b = new BatchEffectTestCli();
        b.doWork( args );
    }

    @Override
    public CommandGroup getCommandGroup() {
        return CommandGroup.EXPERIMENT;
    }

    /*
     * (non-Javadoc)
     *
     * @see ubic.gemma.util.AbstractCLI#getCommandName()
     */
    @Override
    public String getCommandName() {
        return "batchEffectTest";
    }

    @Override
    public String getShortDesc() {
        return "Test for batch effects";
    }

    /*
     * (non-Javadoc)
     *
     * @see ubic.gemma.apps.ExpressionExperimentManipulatingCLI#buildOptions()
     */
    @Override
    protected void buildOptions() {
        super.buildOptions();
    }

    /*
     * (non-Javadoc)
     *
     * @see ubic.gemma.util.AbstractCLI#doWork(java.lang.String[])
     */
    @Override
    protected Exception doWork( String[] args ) {

        Exception ex = super.processCommandLine( args );
        if ( ex != null ) return ex;
        SVDService svdService = this.getBean( SVDService.class );

        for ( BioAssaySet bas : this.expressionExperiments ) {
            if ( bas instanceof ExpressionExperiment ) {
                ExpressionExperiment ee = ( ExpressionExperiment ) bas;
                eeService.thawLite( ee );
                log.info( "Processing: " + ee );

                try {

                    boolean success = false;

                    pcaFactorTest( svdService, ee );

                    Collection<BatchConfoundValueObject> results = BatchConfound.test( ee );

                    for ( BatchConfoundValueObject r : results ) {
                        System.out.println( r );
                    }

                    success = true;
                    if ( success ) {
                        this.successObjects.add( bas.toString() );
                    } else {
                        this.errorObjects.add( bas.toString() + ": No dates found" );

                    }

                } catch ( Exception e ) {
                    log.error( e, e );
                    this.errorObjects.add( bas + ": " + e.getMessage() );
                }

            }
        }

        summarizeProcessing();
        return null;
    }

    /**
     * Just extracts information from the SVD and prints it out.
     *
     * @param svdService
     * @param ee
     */
    private void pcaFactorTest( SVDService svdService, ExpressionExperiment ee ) {
        SVDValueObject svdo = svdService.getSvdFactorAnalysis( ee.getId() );
        Map<Long, ExperimentalFactor> efMap = EntityUtils
                .getIdMap( ee.getExperimentalDesign().getExperimentalFactors() );
        /*
         * Compare PCs to batches.
         */
        Map<Integer, Double> dateCorrelations = svdo.getDateCorrelations();
        Map<Integer, Double> datePvals = svdo.getDatePvals();
        Map<Integer, Map<Long, Double>> factorCorrelations = svdo.getFactorCorrelations();
        Map<Integer, Map<Long, Double>> factorPvals = svdo.getFactorPvals();

        for ( Integer cmp : dateCorrelations.keySet() ) {
            Double dateCorr = dateCorrelations.get( cmp );
            Double datePval = datePvals.get( cmp );
            System.out.println( "PCA\t" + ee.getId() + "\t" + ee.getShortName() + "\tRunDate\tRunDate\tPC" + cmp + "\t"
                    + String.format( "%.2f", dateCorr ) + "\t" + String.format( "%.2g", datePval ) );
        }

        for ( Integer cmp : factorCorrelations.keySet() ) {
            for ( Long efId : factorCorrelations.get( cmp ).keySet() ) {
                Double factorCorr = factorCorrelations.get( cmp ).get( efId );
                Double factorPval = factorPvals.get( cmp ).get( efId );
                System.out.println( "PCA\t" + ee.getId() + "\t" + ee.getShortName() + "\t" + efId + "\t"
                        + efMap.get( efId ).getName() + "\tPC" + cmp + "\t" + String.format( "%.2f", factorCorr )
                        + "\t" + String.format( "%.2g", factorPval ) );
            }
        }
    }

}
