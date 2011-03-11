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

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math.distribution.ChiSquaredDistribution;
import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;
import org.apache.commons.math.stat.inference.ChiSquareTest;
import org.apache.commons.math.stat.inference.ChiSquareTestImpl;

import ubic.gemma.analysis.preprocess.svd.SVDService;
import ubic.gemma.analysis.preprocess.svd.SVDServiceImpl;
import ubic.gemma.analysis.preprocess.svd.SVDValueObject;
import ubic.gemma.apps.ExpressionExperimentManipulatingCLI;
import ubic.gemma.model.expression.bioAssay.BioAssay;
import ubic.gemma.model.expression.biomaterial.BioMaterial;
import ubic.gemma.model.expression.experiment.BioAssaySet;
import ubic.gemma.model.expression.experiment.ExperimentalFactor;
import ubic.gemma.model.expression.experiment.ExpressionExperiment;
import ubic.gemma.model.expression.experiment.FactorValue;
import ubic.gemma.util.EntityUtils;

/**
 * For bulk processing of batch-info-fetching.
 * 
 * @author paul
 * @version $Id$
 */
public class BatchEffectTestCli extends ExpressionExperimentManipulatingCLI {

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

        Exception ex = super.processCommandLine( "BatchEffectPopulation", args );
        if ( ex != null ) return ex;
        SVDService svdService = ( SVDService ) this.getBean( "sVDService" );

        for ( BioAssaySet bas : this.expressionExperiments ) {
            if ( bas instanceof ExpressionExperiment ) {
                ExpressionExperiment ee = eeService.thawLite( ( ExpressionExperiment ) bas );
                log.info( "Processing: " + ee );

                try {
                    Map<ExperimentalFactor, Map<Long, Double>> bioMaterialFactorMap = getBioMaterialFactorMap( ee );

                    Map<Long, ExperimentalFactor> efMap = EntityUtils.getIdMap( ee.getExperimentalDesign()
                            .getExperimentalFactors() );

                    boolean success = false;

                    SVDValueObject svdo = svdService.svdFactorAnalysis( ee );

                    /*
                     * Compare PCs to batches.
                     */
                    Map<Integer, Double> dateCorrelations = svdo.getDateCorrelations();
                    Map<Integer, Map<Long, Double>> factorCorrelations = svdo.getFactorCorrelations();

                    for ( Integer cmp : dateCorrelations.keySet() ) {
                        Double dateCorr = dateCorrelations.get( cmp );
                        System.out.println( "PCA\t" + ee.getId() + "\t" + ee.getShortName() + "\tRunDate\tRunDate\t"
                                + cmp + "\t" + String.format( "%.2f", dateCorr ) );
                    }

                    for ( Integer cmp : factorCorrelations.keySet() ) {
                        for ( Long efId : factorCorrelations.get( cmp ).keySet() ) {
                            Double factorCorr = factorCorrelations.get( cmp ).get( efId );
                            System.out.println( "PCA\t" + ee.getId() + "\t" + ee.getShortName() + "\t" + efId + "\t"
                                    + efMap.get( efId ).getName() + "\tPC" + cmp + "\t"
                                    + String.format( "%.2f", factorCorr ) );
                        }
                    }

                    Long[] svdBioMaterials = svdo.getBioMaterialIds();
                    // double[] batchFactorValues = new double[svdBioMaterials.length];

                    Map<Long, Long> batchMembership = new HashMap<Long, Long>();
                    ExperimentalFactor batchFactor = null;
                    Map<Long, Integer> batchIndexes = new HashMap<Long, Integer>();
                    for ( ExperimentalFactor ef : ee.getExperimentalDesign().getExperimentalFactors() ) {
                        if ( ef.getName().equalsIgnoreCase( "batch" ) ) {
                            batchFactor = ef;

                            int index = 0;
                            for ( FactorValue fv : ef.getFactorValues() ) {
                                batchIndexes.put( fv.getId(), index++ );
                            }

                            Map<Long, Double> bmToFv = bioMaterialFactorMap.get( ef );
                            for ( int j = 0; j < svdBioMaterials.length; j++ ) {
                                batchMembership.put( svdBioMaterials[j], bmToFv.get( svdBioMaterials[j] ).longValue() );
                            }
                            break;
                        }
                    }

                    if ( batchFactor == null ) {
                        continue;
                    }

                    /*
                     * Compare other factors to batches to look for confounds.
                     */

                    for ( ExperimentalFactor ef : ee.getExperimentalDesign().getExperimentalFactors() ) {

                        if ( ef.equals( batchFactor ) ) continue;
                        if ( SVDServiceImpl.isContinuous( ef ) ) {
                            // // // TODO
                        } else {

                            Map<Long, Integer> factorValueIndexes = new HashMap<Long, Integer>();
                            int index = 0;
                            for ( FactorValue fv : ef.getFactorValues() ) {
                                factorValueIndexes.put( fv.getId(), index++ );
                            }
                            Map<Long, Long> factorValueMembership = new HashMap<Long, Long>();
                            Map<Long, Double> bmToFv = bioMaterialFactorMap.get( ef );
                            for ( int j = 0; j < svdBioMaterials.length; j++ ) {
                                factorValueMembership.put( svdBioMaterials[j], bmToFv.get( svdBioMaterials[j] )
                                        .longValue() );
                            }

                            long[][] counts = new long[batchFactor.getFactorValues().size()][ef.getFactorValues()
                                    .size()];

                            for ( int i = 0; i < batchIndexes.size(); i++ ) {
                                for ( int j = 0; j < factorValueIndexes.size(); j++ ) {
                                    counts[i][j] = 0;
                                }
                            }

                            for ( Long bm : svdBioMaterials ) {
                                long fv = factorValueMembership.get( bm );
                                long batch = batchMembership.get( bm );
                                int batchIndex = batchIndexes.get( batch );
                                int factorIndex = factorValueIndexes.get( fv );
                                counts[batchIndex][factorIndex]++;
                            }

                            ChiSquareTest cst = new ChiSquareTestImpl();
                            double chiSquare = cst.chiSquare( counts );
                            double df = ( ( double ) counts.length - 1 ) * ( ( double ) counts[0].length - 1 );
                            ChiSquaredDistribution distribution = new ChiSquaredDistributionImpl( df );
                            double p = 1.0 - distribution.cumulativeProbability( chiSquare );

                            System.out.println( "ChiSq\t" + ee.getId() + "\t" + ee.getShortName() + "\t" + ef.getId()
                                    + "\t" + ef.getName() + "\t" + String.format( "%.2f", chiSquare ) + "\t"
                                    + ( int ) df + "\t" + String.format( "%.2g", p ) );
                        }
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

    private Map<ExperimentalFactor, Map<Long, Double>> getBioMaterialFactorMap( ExpressionExperiment ee ) {
        Map<ExperimentalFactor, Map<Long, Double>> bioMaterialFactorMap = new HashMap<ExperimentalFactor, Map<Long, Double>>();

        // code copied from SVDServiceImpl
        for ( BioAssay bioAssay : ee.getBioAssays() ) {
            for ( BioMaterial bm : bioAssay.getSamplesUsed() ) {
                for ( FactorValue fv : bm.getFactorValues() ) {
                    ExperimentalFactor experimentalFactor = fv.getExperimentalFactor();
                    if ( !bioMaterialFactorMap.containsKey( experimentalFactor ) ) {
                        bioMaterialFactorMap.put( experimentalFactor, new HashMap<Long, Double>() );
                    }

                    double valueToStore;
                    if ( fv.getMeasurement() != null ) {
                        try {
                            valueToStore = Double.parseDouble( fv.getMeasurement().getValue() );
                        } catch ( NumberFormatException e ) {
                            log.warn( "Measurement wasn't a number for " + fv );
                            valueToStore = Double.NaN;
                        }

                    } else {
                        /*
                         * This is a hack. We're storing the ID but as a double.
                         */
                        valueToStore = fv.getId().doubleValue();
                    }
                    bioMaterialFactorMap.get( experimentalFactor ).put( bm.getId(), valueToStore );
                }

            }
        }
        return bioMaterialFactorMap;
    }

    public static void main( String[] args ) {
        BatchEffectTestCli b = new BatchEffectTestCli();
        b.doWork( args );
    }

}
