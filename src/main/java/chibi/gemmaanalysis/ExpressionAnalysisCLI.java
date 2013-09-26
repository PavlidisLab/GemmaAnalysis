/**
 * 
 */
package chibi.gemmaanalysis;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.lang3.time.StopWatch;

import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.io.writer.MatrixWriter;
import ubic.gemma.expression.experiment.service.ExpressionExperimentService;
import ubic.gemma.model.expression.arrayDesign.ArrayDesign;
import ubic.gemma.model.expression.arrayDesign.ArrayDesignService;
import ubic.gemma.model.expression.bioAssayData.ProcessedExpressionDataVector;
import ubic.gemma.model.expression.bioAssayData.ProcessedExpressionDataVectorService;
import ubic.gemma.model.expression.designElement.CompositeSequence;
import ubic.gemma.model.expression.designElement.CompositeSequenceService;
import ubic.gemma.model.expression.experiment.BioAssaySet;
import ubic.gemma.model.expression.experiment.ExpressionExperiment;
import ubic.gemma.model.genome.Gene;

/**
 * Create a relative expression level (dedv rank) matrix for a list of genes
 * 
 * @author raymond
 */
public class ExpressionAnalysisCLI extends AbstractGeneCoexpressionManipulatingCLI {
    /**
     * @param args
     */
    public static void main( String[] args ) {
        ExpressionAnalysisCLI analysis = new ExpressionAnalysisCLI();
        StopWatch watch = new StopWatch();
        watch.start();

        log.info( "Starting expression analysis" );
        Exception e = analysis.doWork( args );
        if ( e != null ) log.error( e.getMessage() );
        watch.stop();
        log.info( "Finished expression analysis in " + watch );
    }

    private String outFilePrefix;

    private ArrayDesignService adService;

    public static final double DEFAULT_FILTER_THRESHOLD = 0.8;

    ProcessedExpressionDataVectorService processedExpressionDataVectorService;

    CompositeSequenceService compositeSequenceService;

    /*
     * (non-Javadoc)
     * 
     * @see ubic.gemma.util.AbstractCLI#buildOptions()
     */
    @Override
    protected void buildOptions() {
        super.buildOptions();

        Option outFileOption = OptionBuilder.create( 'o' );
        addOption( outFileOption );

        Option filterOption = OptionBuilder.create( "threshold" );
        addOption( filterOption );
    }

    /*
     * (non-Javadoc)
     * 
     * @see ubic.gemma.util.AbstractCLI#doWork(java.lang.String[])
     */
    @Override
    protected Exception doWork( String[] args ) {
        Exception e = processCommandLine( "ExpressionAnalysis", args );
        if ( e != null ) return e;

        Collection<Gene> genes;

        log.info( "Getting genes" );
        genes = geneService.loadAll( taxon );
        log.info( "Loaded " + genes.size() + " genes" );

        DoubleMatrix<Gene, ExpressionExperiment> rankMatrix = getRankMatrix( genes, expressionExperiments );
        // rankMatrix = filterRankmatrix(rankMatrix);

        // gene names
        Collection<Gene> rowGenes = rankMatrix.getRowNames();
        try (PrintWriter out = new PrintWriter( new FileWriter( outFilePrefix + ".row_names.txt" ) );) {
            for ( Gene gene : rowGenes ) {
                String s = gene.getOfficialSymbol();
                if ( s == null ) s = gene.getId().toString();
                out.println( s );
            }
        } catch ( IOException exc ) {
            return exc;
        }

        // expression experiment names
        Collection<ExpressionExperiment> colEes = rankMatrix.getColNames();
        try (PrintWriter out = new PrintWriter( new FileWriter( outFilePrefix + ".col_names.txt" ) );) {
            for ( ExpressionExperiment ee : colEes ) {
                out.println( ee.getShortName() );
            }
        } catch ( IOException exc ) {
            return exc;
        }

        DecimalFormat formatter = ( DecimalFormat ) NumberFormat.getNumberInstance( Locale.US );
        formatter.applyPattern( "0.0000" );
        formatter.getDecimalFormatSymbols().setNaN( "NaN" );
        try {
            MatrixWriter<Gene, ExpressionExperiment> out = new MatrixWriter<Gene, ExpressionExperiment>( outFilePrefix
                    + ".txt", formatter );
            out.writeMatrix( rankMatrix, false );
        } catch ( IOException exc ) {
            return exc;
        }

        return null;
    }

    /**
     * 
     */
    protected void initBeans() {
        eeService = getBean( ExpressionExperimentService.class );
        adService = getBean( ArrayDesignService.class );
        processedExpressionDataVectorService = this.getBean( ProcessedExpressionDataVectorService.class );
    }

    @Override
    protected void processOptions() {
        super.processOptions();
        if ( hasOption( 'o' ) ) {
            outFilePrefix = getOptionValue( 'o' );
        }

        initBeans();
    }

    private Map<Long, Collection<Long>> getCs2GeneMap( Collection<Long> csIds ) {
        Map<CompositeSequence, Collection<Gene>> genes = compositeSequenceService.getGenes( compositeSequenceService
                .loadMultiple( csIds ) );
        Map<Long, Collection<Long>> result = new HashMap<Long, Collection<Long>>();
        for ( CompositeSequence cs : genes.keySet() ) {
            result.put( cs.getId(), new HashSet<Long>() );
            for ( Gene g : genes.get( cs ) ) {
                result.get( cs.getId() ).add( g.getId() );
            }
        }
        return result;
    }

    /**
     * @param genes
     * @param ees
     * @return
     */
    private DenseDoubleMatrix<Gene, ExpressionExperiment> getRankMatrix( Collection<Gene> genes,
            Collection<BioAssaySet> ees ) {
        DenseDoubleMatrix<Gene, ExpressionExperiment> matrix = new DenseDoubleMatrix<Gene, ExpressionExperiment>(
                genes.size(), ees.size() );
        for ( int i = 0; i < matrix.rows(); i++ ) {
            for ( int j = 0; j < matrix.columns(); j++ ) {
                matrix.set( i, j, Double.NaN );
            }
        }
        // name rows + cols
        for ( Gene gene : genes ) {
            matrix.addRowName( gene );
        }
        for ( BioAssaySet bas : ees ) {
            ExpressionExperiment ee = ( ExpressionExperiment ) bas;
            matrix.addColumnName( ee );
        }

        int eeCount = 1;
        for ( BioAssaySet bas : ees ) {
            ExpressionExperiment ee = ( ExpressionExperiment ) bas;
            int col = matrix.getColIndexByName( ee );
            log.info( "Processing " + ee.getShortName() + " (" + eeCount++ + " of " + ees.size() + ")" );
            Collection<ArrayDesign> ads = eeService.getArrayDesignsUsed( ee );
            Collection<CompositeSequence> css = new HashSet<CompositeSequence>();
            for ( ArrayDesign ad : ads ) {
                css.addAll( adService.getCompositeSequences( ad ) );
            }

            Collection<ProcessedExpressionDataVector> dedvs = processedExpressionDataVectorService
                    .getProcessedDataVectors( ee );

            // get cs2gene map
            Collection<ArrayDesign> ADs = eeService.getArrayDesignsUsed( ee );
            Collection<Long> csIds = new HashSet<Long>();
            for ( ArrayDesign AD : ADs ) {
                Collection<CompositeSequence> CSs = adService.getCompositeSequences( AD );
                for ( CompositeSequence CS : CSs ) {
                    csIds.add( CS.getId() );
                }
            }
            // FIXME this use dto return only known genes.
            Map<Long, Collection<Long>> cs2geneMap = getCs2GeneMap( csIds );

            // invert dedv2geneMap
            Map<Long, Collection<ProcessedExpressionDataVector>> gene2dedvMap = new HashMap<Long, Collection<ProcessedExpressionDataVector>>();
            for ( ProcessedExpressionDataVector dedv : dedvs ) {
                Collection<Long> c = cs2geneMap.get( dedv.getDesignElement().getId() );
                for ( Long gene : c ) {
                    Collection<ProcessedExpressionDataVector> vs = gene2dedvMap.get( gene );
                    if ( vs == null ) {
                        vs = new HashSet<ProcessedExpressionDataVector>();
                        gene2dedvMap.put( gene, vs );
                    }
                    vs.add( dedv );
                }

            }
            log.info( "Loaded design element data vectors" );

            // construct the rank matrix
            int rankCount = 0;
            for ( Gene gene : genes ) {
                int row = matrix.getRowIndexByName( gene );
                Double rank;
                List<Double> ranks = new ArrayList<Double>();
                Collection<ProcessedExpressionDataVector> vs = gene2dedvMap.get( gene.getId() );
                if ( vs == null ) continue;
                for ( ProcessedExpressionDataVector dedv : vs ) {
                    ranks.add( dedv.getRankByMean() );
                }
                if ( ranks.size() < 1 ) continue;

                // take the median rank
                Collections.sort( ranks );
                rank = ranks.get( ranks.size() / 2 );
                if ( rank == null ) continue;
                matrix.set( row, col, rank );
                rankCount++;
            }
            log.info( "Saved " + rankCount + " gene ranks" );
        }

        return matrix;
    }

}
