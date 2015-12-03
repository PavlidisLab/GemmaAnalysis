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

import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.collections.Transformer;
import org.apache.commons.lang3.StringUtils;

import ubic.gemma.analysis.expression.diff.DiffExAnalyzer;
import ubic.gemma.analysis.expression.diff.DifferentialExpressionAnalysisConfig;
import ubic.gemma.analysis.expression.diff.LinearModelAnalyzer;
import ubic.gemma.analysis.preprocess.batcheffects.ExpressionExperimentBatchCorrectionService;
import ubic.gemma.apps.DifferentialExpressionAnalysisCli;
import ubic.gemma.datastructure.matrix.ExpressionDataDoubleMatrix;
import ubic.gemma.datastructure.matrix.MatrixWriter;
import ubic.gemma.model.analysis.expression.diff.DifferentialExpressionAnalysis;
import ubic.gemma.model.analysis.expression.diff.DifferentialExpressionAnalysisResult;
import ubic.gemma.model.analysis.expression.diff.ExpressionAnalysisResultSet;
import ubic.gemma.model.expression.arrayDesign.ArrayDesign;
import ubic.gemma.model.expression.arrayDesign.ArrayDesignService;
import ubic.gemma.model.expression.bioAssayData.ProcessedExpressionDataVector;
import ubic.gemma.model.expression.bioAssayData.ProcessedExpressionDataVectorService;
import ubic.gemma.model.expression.designElement.CompositeSequence;
import ubic.gemma.model.expression.designElement.CompositeSequenceService;
import ubic.gemma.model.expression.experiment.BioAssaySet;
import ubic.gemma.model.expression.experiment.ExperimentalFactor;
import ubic.gemma.model.expression.experiment.ExpressionExperiment;
import ubic.gemma.model.genome.Gene;
import ubic.gemma.util.EntityUtils;

/**
 * Performs multiple differential expression analyses under different conditions: Without including a batch covariate;
 * with including it; and repeating those, after batch correction
 * 
 * @author paul
 * @version $Id$
 */
public class BatchDiffExCli extends DifferentialExpressionAnalysisCli {
    private static final int LOGGING_FREQ = 20000;

    @Override
    public String getCommandName() {
        return null;
    }

    @Override
    public String getShortDesc() {
        return "Performs multiple differential expression analyses under different conditions: Without including a batch covariate; "
                + "with including it; and repeating those, after batch correction";
    }

    /**
     * @param args
     */
    public static void main( String[] args ) {
        BatchDiffExCli c = new BatchDiffExCli();
        c.doWork( args );

    }

    private ExpressionExperimentBatchCorrectionService expressionExperimentBatchCorrectionService;

    private DiffExAnalyzer lma;

    private ArrayDesignService arrayDesignService;

    private CompositeSequenceService compositeSequenceService;

    private ProcessedExpressionDataVectorService processedExpressionDataVectorService;

    private Collection<ArrayDesign> seenArrays = new HashSet<ArrayDesign>();

    private Map<CompositeSequence, Collection<Gene>> genes = new HashMap<CompositeSequence, Collection<Gene>>();

    Transformer geneSymbolTransformer = new Transformer() {
        @Override
        public Object transform( Object input ) {
            return ( ( Gene ) input ).getOfficialSymbol();
        }
    };

    /**
     * This only affects the summaries that are output.
     */
    private double summaryQvalThreshold = 0.01;

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

    /**
     * @param mat
     * @param filename
     * @throws IOException
     */
    private void saveData( ExpressionDataDoubleMatrix mat, String filename ) throws IOException {
        MatrixWriter mw = new MatrixWriter();
        try (FileWriter fw = new FileWriter( new File( filename ) );) {
            mw.write( fw, mat, null, true, true );
        }
    }

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

        if ( pval != null && pval < summaryQvalThreshold ) {
            c++;
        }
        /*
         * Map of probe -> factor -> pval
         */
        CompositeSequence probe = r.getProbe();
        if ( !revisedResultDetails.containsKey( probe ) ) {
            revisedResultDetails.put( probe, new HashMap<ExperimentalFactor, Double>() );
        }
        revisedResultDetails.get( probe ).put( ef, pval );
        return c;
    }

    /*
     * (non-Javadoc)
     * 
     * @see ubic.gemma.util.AbstractCLI#doWork(java.lang.String[])
     */
    @Override
    protected Exception doWork( String[] args ) {
        Exception error = super.processCommandLine( args );
        if ( error != null ) return error;
        this.expressionExperimentBatchCorrectionService = this
                .getBean( ExpressionExperimentBatchCorrectionService.class );
        this.lma = new LinearModelAnalyzer();
        this.processedExpressionDataVectorService = this.getBean( ProcessedExpressionDataVectorService.class );
        this.compositeSequenceService = this.getBean( CompositeSequenceService.class );

        arrayDesignService = this.getBean( ArrayDesignService.class );

        try {
            summaryFile = initOutputFile( "batch.proc.summary.txt" );
            summaryFile.write( "State\tEEID\tEENAME\tEFID\tEFNAME\tNUM\tNUMDIFF\n" );

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
        String fileprefix = ee.getId() + "." + ee.getShortName().replaceAll( "[\\W\\s]+", "_" );

        try (Writer detailFile = initOutputFile( "batch.proc.detail." + fileprefix + ".txt" );) {

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

            if ( experimentalFactors.size() > 10 ) {
                /*
                 * This could be modified to select just a few factors, at random ... but that's probably
                 */
                this.errorObjects.add( "Too many factors (" + experimentalFactors.size()
                        + " factors, including 'batch'): " + ee.getShortName() );
                return;
            }

            /* TODO use this, or skip it... we have this information elsewhere already */
            // expressionExperimentBatchCorrectionService.checkBatchEffectSeverity( ee );

            boolean correctable = expressionExperimentBatchCorrectionService.checkCorrectability( ee );
            if ( !correctable ) {

                /*
                 * Note that later on we can still end up with a model that is not of full rank, so combat will fail.
                 * This can sometimes be ameliorated by dropping covariates.
                 */
                this.errorObjects
                        .add( "Batch effect is not correctable; possibly contains batches with only one sample: "
                                + ee.getShortName() );

                return;
            }

            log.info( "Processing: " + ee );

            /*
             * Extract data
             */
            Collection<ProcessedExpressionDataVector> vectos = processedExpressionDataVectorService
                    .getProcessedDataVectors( ee );

            ExpressionDataDoubleMatrix mat = new ExpressionDataDoubleMatrix( vectos );

            /*
             * TODO for some data sets we should re-normalize?
             */

            StringBuilder summaryBuf = new StringBuilder();

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
            configWithoutBatch.setQvalueThreshold( null );
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
                summaryBuf.append( "Before\t" + ee.getId() + "\t" + ee.getShortName() + "\t" + ef.getId() + "\t"
                        + ef.getName() + "\t" + results.size() + "\t" + c + "\n" );
            }

            /*
             * Then do it with batch.
             */
            Collection<ExperimentalFactor> factors = experimentalFactors;
            assert factors.contains( batchFactor );
            DifferentialExpressionAnalysisConfig configIncludingBatch = new DifferentialExpressionAnalysisConfig();
            configIncludingBatch.setQvalueThreshold( summaryQvalThreshold );

            configIncludingBatch.setFactorsToInclude( factors );

            DifferentialExpressionAnalysis withBatchEffectResults = lma.run( ee, mat, configIncludingBatch ).iterator()
                    .next();

            /*
             * Determine how many genes are diff ex wrt batch. The other factors are tracked; this shows how we would do
             * if we tried to simply directly include batch in the model
             */
            Map<CompositeSequence, Map<ExperimentalFactor, Double>> batchEffectDetails = new HashMap<>();

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
                summaryBuf.append( "Batch\t" + ee.getId() + "\t" + ee.getShortName() + "\t" + ef.getId() + "\t"
                        + ef.getName() + "\t" + results.size() + "\t" + c + "\n" );
            }

            /*
             * Correct for batch effects; covariates which are "unimportant" will be dropped.
             */
            log.info( "ComBat-ing" );
            boolean parametric = true;
            double importanceThreshold = 0.01;
            ExpressionDataDoubleMatrix comBat = expressionExperimentBatchCorrectionService.comBat( mat, parametric,
                    importanceThreshold );
            assert comBat != null;

            /*
             * Check if we have removed the batch effect: there should be no diff ex wrt batch. This is just a sanity
             * check, really. The other factors are tracked just for completeness. Note that Combat log transforms the
             * data if necessary, but transforms it back.
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
                summaryBuf.append( "BatchAftCorr\t" + ee.getId() + "\t" + ee.getShortName() + "\t" + ef.getId() + "\t"
                        + ef.getName() + "\t" + results.size() + "\t" + c + "\n" );
            }

            /*
             * Now without batch as a factor, which is what we really want.
             */
            boolean hasNonNulls = false;
            DifferentialExpressionAnalysis revisedResult = lma.run( ee, comBat, configWithoutBatch ).iterator().next();
            Map<CompositeSequence, Map<ExperimentalFactor, Double>> revisedResultDetails = new HashMap<CompositeSequence, Map<ExperimentalFactor, Double>>();
            for ( ExpressionAnalysisResultSet brs : revisedResult.getResultSets() ) {
                assert brs.getExperimentalFactors().size() == 1;
                ExperimentalFactor ef = brs.getExperimentalFactors().iterator().next();
                Collection<DifferentialExpressionAnalysisResult> results = brs.getResults();
                int c = 0;
                for ( DifferentialExpressionAnalysisResult r : results ) {

                    if ( r.getCorrectedPvalue() != null && !Double.isNaN( r.getCorrectedPvalue() ) ) {
                        hasNonNulls = true;
                    }

                    c = tally( revisedResultDetails, ef, r, c );

                    if ( ++j % LOGGING_FREQ == 0 ) {
                        log.info( j + " processed" );
                    }

                }
                summaryBuf.append( "After\t" + ee.getId() + "\t" + ee.getShortName() + "\t" + ef.getId() + "\t"
                        + ef.getName() + "\t" + results.size() + "\t" + c + "\n" );

            }

            if ( !hasNonNulls ) {
                // this means something went wrong ... somewhere. Possibly the model cannot be fit.
                errorObjects.add( "No valid pvalues after correction: " + ee.getShortName() );
                return;
            }

            /*
             * Print out details
             */

            detailFile
                    .write( "EEID\tEENAME\tEFID\tEFNAME\tPROBEID\tPROBENAME\tGENESYMBS\tGENEIDS\tBEFOREQVAL\tBATCHQVAL\tBATAFTERQVAL\tAFTERQVAL\n" );

            getGeneAnnotations( ee );

            for ( CompositeSequence c : beforeResultDetails.keySet() ) {

                // Get the gene information
                String geneSymbs = "";
                String geneIds = "";
                if ( genes.containsKey( c ) ) {
                    Collection<Gene> g = new HashSet<Gene>();
                    g.addAll( genes.get( c ) );
                    CollectionUtils.transform( g, geneSymbolTransformer );
                    geneSymbs = StringUtils.join( g, "|" );
                    geneIds = StringUtils.join( EntityUtils.getIds( genes.get( c ) ), "|" );
                }

                for ( ExperimentalFactor ef : factors ) {

                    detailFile.write( ee.getId() + "\t" + ee.getShortName() + "\t" + ef.getId() + "\t" + ef.getName()
                            + "\t" + c.getId() + "\t" + c.getName() + "\t" + geneSymbs + "\t" + geneIds + "\t" );

                    Double bpval = beforeResultDetails.get( c ).get( ef ); // will be null for 'batch'

                    Double batpval = batchEffectDetails.get( c ).get( ef ); // when batch was included.

                    Double batapval = batchEffectAfterCorrDetails.get( c ).get( ef ); // when batch was included.

                    Double aftpval = revisedResultDetails.get( c ).get( ef ); // will be null for 'batch'

                    detailFile.write( String.format( "%.4g\t%.4g\t%.4g\t%.4g\n", bpval, batpval, batapval, aftpval ) );

                }
            }
            detailFile.close();

            summaryFile.write( summaryBuf.toString() );
            summaryFile.flush();

            String rawDataFileName = fileprefix + ".originaldata.txt";
            saveData( mat, rawDataFileName );
            String correctedDataFileName = fileprefix + ".correcteddata.txt";
            saveData( comBat, correctedDataFileName );

            successObjects.add( ee );
        } catch ( Exception e ) {
            log.error( e, e );
            errorObjects.add( ee + e.getMessage() );
        }
    }

    /**
     * Gets annotations for this experiment, if the array designs have not already been seen in this run.
     * 
     * @param ee
     */
    private void getGeneAnnotations( ExpressionExperiment ee ) {
        Collection<ArrayDesign> arrayDesigns = this.eeService.getArrayDesignsUsed( ee );
        for ( ArrayDesign ad : arrayDesigns ) {
            if ( seenArrays.contains( ad ) ) continue;
            ad = this.arrayDesignService.thaw( ad );
            genes.putAll( compositeSequenceService.getGenes( ad.getCompositeSequences() ) );
            seenArrays.add( ad );
        }
    }
}
