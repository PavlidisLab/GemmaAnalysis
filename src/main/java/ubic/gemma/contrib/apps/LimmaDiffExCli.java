/*
 * The GemmaAnalysis project
 *
 * Copyright (c) 2017 University of British Columbia
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

package ubic.gemma.contrib.apps;

import cern.colt.list.DoubleArrayList;
import org.apache.commons.lang3.StringUtils;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.beans.factory.annotation.Value;
import ubic.basecode.math.Distance;
import ubic.gemma.apps.DifferentialExpressionAnalysisCli;
import ubic.gemma.core.analysis.expression.diff.DiffExAnalyzer;
import ubic.gemma.core.analysis.expression.diff.DifferentialExpressionAnalysisConfig;
import ubic.gemma.core.datastructure.matrix.ExpressionDataDoubleMatrix;
import ubic.gemma.model.analysis.expression.diff.ContrastResult;
import ubic.gemma.model.analysis.expression.diff.DifferentialExpressionAnalysis;
import ubic.gemma.model.analysis.expression.diff.DifferentialExpressionAnalysisResult;
import ubic.gemma.model.analysis.expression.diff.ExpressionAnalysisResultSet;
import ubic.gemma.model.expression.arrayDesign.ArrayDesign;
import ubic.gemma.model.expression.bioAssayData.ProcessedExpressionDataVector;
import ubic.gemma.model.expression.designElement.CompositeSequence;
import ubic.gemma.model.expression.experiment.ExperimentalDesignUtils;
import ubic.gemma.model.expression.experiment.ExperimentalFactor;
import ubic.gemma.model.expression.experiment.ExpressionExperiment;
import ubic.gemma.model.expression.experiment.FactorValue;
import ubic.gemma.model.genome.Gene;
import ubic.gemma.persistence.service.expression.arrayDesign.ArrayDesignService;
import ubic.gemma.persistence.service.expression.bioAssayData.ProcessedExpressionDataVectorService;
import ubic.gemma.persistence.service.expression.designElement.CompositeSequenceService;
import ubic.gemma.persistence.util.IdentifiableUtils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Performs differential expression analyses with and without ebayes. Based on BachDiffExCli
 *
 * @author paul
 */
public class LimmaDiffExCli extends DifferentialExpressionAnalysisCli {
    private static final int LOGGING_FREQ = 20000;

    /**
     * This only affects the summaries that are output.
     */
    private static final double summaryQvalThreshold = 0.01;

    @Autowired
    private DiffExAnalyzer lma;

    @Autowired
    private ArrayDesignService arrayDesignService;

    @Autowired
    private CompositeSequenceService compositeSequenceService;

    @Autowired
    private ProcessedExpressionDataVectorService processedExpressionDataVectorService;

    @Value("${gemma.download.path}")
    private Path downloadPath;

    private final Collection<ArrayDesign> seenArrays = new HashSet<>();

    private final Map<CompositeSequence, Collection<Gene>> genes = new HashMap<>();

    private Writer summaryFile;

    @Override
    public String getCommandName() {
        return "limmacompare";
    }

    @Override
    public String getShortDesc() {
        return "Performs multiple differential expression analyses with and without ebayes, generate comparison stats";
    }

    @Override
    protected void doAuthenticatedWork() throws Exception {
        try ( Writer summaryFile = initOutputFile( "limma.proc.summary.txt" ) ) {
            summaryFile.write( "State\tEEID\tEENAME\tEFID\tEFNAME\tNUM\tNUMDIFF\n" );
            this.summaryFile = summaryFile;
            super.doAuthenticatedWork();
        } finally {
            this.summaryFile = null;
        }
    }

    @Override
    protected void processExpressionExperiment( ExpressionExperiment ee ) {
        String fileprefix = ee.getId() + "." + ee.getShortName().replaceAll( "[\\W\\s]+", "_" );

        try ( Writer detailFile = initOutputFile( "ebayes.proc.detail." + fileprefix + ".txt" ) ) {

            ee = eeService.thawLite( ee );


            Collection<ExperimentalFactor> experimentalFactors = ee.getExperimentalDesign().getExperimentalFactors();

            if ( experimentalFactors.size() > 10 ) {
                /*
                 * This could be modified to select just a few factors, at random ... but that's probably
                 */
                this.addErrorObject( ee, "Too many factors (" + experimentalFactors.size()
                        + " factors : " + ee.getShortName() );
                return;
            }

            log.info( "===== Processing: " + ee );

            /*
             * Extract data FIXME make this ONE STEP to getting the data matrix.
             */
            Collection<ProcessedExpressionDataVector> vectos = processedExpressionDataVectorService
                    .getProcessedDataVectors( ee );
            processedExpressionDataVectorService.thaw( vectos );
            ExpressionDataDoubleMatrix mat = new ExpressionDataDoubleMatrix( vectos );

            StringBuilder summaryBuf = new StringBuilder();

            /*
             * first do an analysis without ebayes; this is our baseline. Let's ignore interactions to keep things
             * simple.
             */
            Collection<ExperimentalFactor> factorsToAnalyze = new HashSet<>();
            for ( ExperimentalFactor ef : experimentalFactors ) {
                if ( ExperimentalDesignUtils.isBatchFactor( ef ) ) continue;
                factorsToAnalyze.add( ef );
            }
            int j = 0;
            DifferentialExpressionAnalysisConfig config1 = new DifferentialExpressionAnalysisConfig();
            config1.addFactorsToInclude( factorsToAnalyze );
            config1.setMakeArchiveFile( false );
            config1.setModerateStatistics( false ); // <----
            log.info( "=== Nobayes === " );
            Collection<DifferentialExpressionAnalysis> deas = lma.run( ee, mat, config1 );
            if ( deas.isEmpty() ) {
                log.error( "No differential expression results obtained, moving on" );
                addErrorObject( ee, "No differential expression results obtained" );
                return;
            }
            DifferentialExpressionAnalysis beforeResults = deas.iterator().next();
            Map<CompositeSequence, Map<ExperimentalFactor, Double>> beforeResultDetails = new HashMap<>();
            for ( ExpressionAnalysisResultSet brs : beforeResults.getResultSets() ) {
                ExperimentalFactor ef = brs.getExperimentalFactors().iterator().next();
                Collection<DifferentialExpressionAnalysisResult> results = brs.getResults();
                int c = 0;
                for ( DifferentialExpressionAnalysisResult r : results ) {
                    c = tally( beforeResultDetails, ef, r, c );
                    if ( ++j % LOGGING_FREQ == 0 ) {
                        log.info( j + " processed" );
                    }
                }
                summaryBuf.append( "Nobayes\t" + ee.getId() + "\t" + ee.getShortName() + "\t" + ef.getId() + "\t"
                        + ef.getName() + "\t" + results.size() + "\t" + c + "\n" );
            }

            /*
             * Then do it with ebayes.
             */
            DifferentialExpressionAnalysisConfig config2 = new DifferentialExpressionAnalysisConfig();
            config2.setModerateStatistics( true ); // <----
            config2.addFactorsToInclude( factorsToAnalyze );
            config2.setMakeArchiveFile( false );
            log.info( "=== With ebayes ===" );
            Collection<DifferentialExpressionAnalysis> deas2 = lma.run( ee, mat, config2 );
            if ( deas2.isEmpty() ) {
                log.error( "No differential expression results obtained with eBayes, moving on" );
                addErrorObject( ee, "No differential expression results obtained with eBayes" );
                return;
            }
            DifferentialExpressionAnalysis eBayesResults = deas2.iterator()
                    .next();

            Map<CompositeSequence, Map<ExperimentalFactor, Double>> afterDetails = new HashMap<>();

            for ( ExpressionAnalysisResultSet brs : eBayesResults.getResultSets() ) {
                ExperimentalFactor ef = brs.getExperimentalFactors().iterator().next();
                Collection<DifferentialExpressionAnalysisResult> results = brs.getResults();
                int c = 0;
                for ( DifferentialExpressionAnalysisResult r : results ) {
                    c = tally( afterDetails, ef, r, c );
                    if ( ++j % LOGGING_FREQ == 0 ) {
                        log.info( j + " processed" );
                    }

                }
                summaryBuf.append( "Ebayes\t" ).append( ee.getId() ).append( "\t" )
                        .append( ee.getShortName() ).append( "\t" )
                        .append( ef.getId() ).append( "\t" ).append( ef.getName() ).append( "\t" )
                        .append( results.size() ).append( "\t" ).append( c ).append( "\n" );
            }

            List<String> correlations = compare( beforeResults, eBayesResults );
            summaryBuf.append( "Correlations\t" ).append( ee.getId() ).append( "\t" )
                    .append( ee.getShortName() ).append( "\t" )
                    .append( StringUtils.join( correlations, " " ) ).append( "\n" );

            /*
             * Print out details
             */

            detailFile
                    .write( "EEID\tEENAME\tEFID\tEFNAME\tPROBEID\tPROBENAME\tGENESYMBS\tGENEIDS\tNOBAYESQVAL\tEBAYESQVAL\n" );

            getGeneAnnotations( ee );

            for ( CompositeSequence c : beforeResultDetails.keySet() ) {

                // Get the gene information
                String geneSymbs = "";
                String geneIds = "";
                if ( genes.containsKey( c ) ) {
                    LinkedHashSet<Gene> g = new LinkedHashSet<>( genes.get( c ) );
                    geneSymbs = g.stream().map( Gene::getOfficialSymbol ).collect( Collectors.joining( "|" ) );
                    geneIds = StringUtils.join( IdentifiableUtils.getIds( g ), "|" );
                }

                for ( ExperimentalFactor ef : factorsToAnalyze ) {

                    detailFile.write( ee.getId() + "\t" + ee.getShortName() + "\t" + ef.getId() + "\t" + ef.getName()
                            + "\t" + c.getId() + "\t" + c.getName() + "\t" + geneSymbs + "\t" + geneIds + "\t" );

                    Double bpval = beforeResultDetails.get( c ).get( ef ); // will be null for 'batch'

                    Double batpval = afterDetails.get( c ).get( ef ); // when batch was included.

                    detailFile.write( String.format( "%.4g\t%.4g\n", bpval, batpval ) );

                }
            }
            detailFile.close();

            summaryFile.write( summaryBuf.toString() );
            summaryFile.flush();
            addSuccessObject( ee, "" );
        } catch ( Exception e ) {
            log.error( e, e );
            addErrorObject( ee, e.getMessage() );
        }
        log.info( "==== Completed processing: " + ee );
    }

    /**
     * @return vector of rank correlations.
     */
    private List<String> compare( DifferentialExpressionAnalysis beforeResults, DifferentialExpressionAnalysis eBayesResults ) {

        Collection<ExpressionAnalysisResultSet> beforeResultSets = beforeResults.getResultSets();
        Collection<ExpressionAnalysisResultSet> eBayesResultSets = eBayesResults.getResultSets();

        /*
         * Collect vectors of values and compute rank correlations. Other comparisions can be done as well...
         */

        Map<ExperimentalFactor, Map<CompositeSequence, Double[]>> apv = new HashMap<>(); // ANOVA effects
        Map<FactorValue, Map<CompositeSequence, Double[]>> cpv = new HashMap<>(); // contrasts

        for ( ExpressionAnalysisResultSet brs : beforeResultSets ) {
            // Note: ignoring interactions.
            if ( brs.getExperimentalFactors().size() > 1 ) {
                continue;
            }
            ExperimentalFactor ef = brs.getExperimentalFactors().iterator().next();
            if ( !apv.containsKey( ef ) ) apv.put( ef, new HashMap<CompositeSequence, Double[]>() );

            for ( DifferentialExpressionAnalysisResult r : brs.getResults() ) {

                Double pvalue2 = r.getPvalue();
                CompositeSequence probe = r.getProbe();
                if ( !apv.get( ef ).containsKey( probe ) ) apv.get( ef ).put( probe, new Double[2] );
                apv.get( ef ).get( probe )[0] = pvalue2;

                for ( ContrastResult cr : r.getContrasts() ) {
                    FactorValue factorValue = cr.getFactorValue();
                    Double pvalue = cr.getPvalue();

                    if ( !cpv.containsKey( factorValue ) )
                        cpv.put( factorValue, new HashMap<CompositeSequence, Double[]>() );
                    if ( !cpv.get( factorValue ).containsKey( probe ) )
                        cpv.get( factorValue ).put( probe, new Double[2] );

                    cpv.get( factorValue ).get( probe )[0] = pvalue;

                }
            }
        }

        boolean warned = false;
        for ( ExpressionAnalysisResultSet brs : eBayesResultSets ) {
            if ( brs.getExperimentalFactors().size() > 1 ) {
                continue;
            }
            ExperimentalFactor ef = brs.getExperimentalFactors().iterator().next();

            for ( DifferentialExpressionAnalysisResult r : brs.getResults() ) {

                Double pvalue2 = r.getPvalue();
                CompositeSequence probe = r.getProbe();
                assert apv.get( ef ) != null;
                if ( apv.get( ef ).get( probe ) == null ) {
                    // If it's missing from the nobayes analysis we put in NaN. Not sure why that happens.
                    if ( !warned ) {
                        log.warn( "No nobayes result for " + probe + ", further warnings suppressed" );
                        warned = true;
                    }
                    apv.get( ef ).put( probe, new Double[2] );
                    apv.get( ef ).get( probe )[0] = Double.NaN;
                }

                apv.get( ef ).get( probe )[1] = pvalue2;

                for ( ContrastResult cr : r.getContrasts() ) {
                    FactorValue factorValue = cr.getFactorValue();
                    Double pvalue = cr.getPvalue();
                    if ( cpv.get( factorValue ).get( probe ) == null ) {
                        cpv.get( factorValue ).put( probe, new Double[2] );
                        cpv.get( factorValue ).get( probe )[0] = Double.NaN;
                    }
                    cpv.get( factorValue ).get( probe )[1] = pvalue;

                }
            }
        }

        List<String> r = new ArrayList<>();

        // unroll
        // ANOVA effects
        for ( ExperimentalFactor ef : apv.keySet() ) {
            int n = apv.get( ef ).size();
            DoubleArrayList a = new DoubleArrayList( n );
            DoubleArrayList b = new DoubleArrayList( n );

            for ( CompositeSequence cs : apv.get( ef ).keySet() ) {
                Double[] pvs = apv.get( ef ).get( cs );
                a.add( pvs[0] );
                b.add( pvs[1] );
            }

            String corr = String.format( "%.5f", Distance.spearmanRankCorrelation( a, b ) );
            r.add( corr );
        }

        // Contrasts
        for ( FactorValue fv : cpv.keySet() ) {
            int n = cpv.get( fv ).size();
            DoubleArrayList a = new DoubleArrayList( n );
            DoubleArrayList b = new DoubleArrayList( n );
            for ( CompositeSequence cs : cpv.get( fv ).keySet() ) {
                Double[] pvs = cpv.get( fv ).get( cs );
                a.add( pvs[0] == null ? Double.NaN : pvs[0] );
                b.add( pvs[1] == null ? Double.NaN : pvs[1] );
            }

            String corr = String.format( "%.5f", Distance.spearmanRankCorrelation( a, b ) );
            r.add( corr );
        }
        return r;
    }

    /**
     * Gets annotations for this experiment, if the array designs have not already been seen in this run.
     */
    private void getGeneAnnotations( ExpressionExperiment ee ) {
        Collection<ArrayDesign> arrayDesigns = eeService.getArrayDesignsUsed( ee );
        for ( ArrayDesign ad : arrayDesigns ) {
            if ( seenArrays.contains( ad ) ) continue;
            ad = this.arrayDesignService.thaw( ad );
            genes.putAll( compositeSequenceService.getGenes( ad.getCompositeSequences() ) );
            seenArrays.add( ad );
        }
    }

    private Writer initOutputFile( String fileName ) throws IOException {
        File f = downloadPath.resolve( fileName ).toFile();
        if ( f.exists() ) {
            f.delete();
        }
        f.createNewFile();
        log.info( "New file: " + f.getAbsolutePath() );
        return new FileWriter( f );
    }

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
}
