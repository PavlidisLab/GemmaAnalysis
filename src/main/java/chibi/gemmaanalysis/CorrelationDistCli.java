/*
 * The Gemma project
 * 
 * Copyright (c) 2007 University of British Columbia
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
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.lang.time.StopWatch;

import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.graphics.ColorMap;
import ubic.basecode.graphics.ColorMatrix;
import ubic.basecode.graphics.MatrixDisplay;
import ubic.basecode.io.ByteArrayConverter;
import ubic.basecode.math.CorrelationStats;
import ubic.gemma.apps.ExpressionExperimentManipulatingCLI;
import ubic.gemma.model.common.quantitationtype.QuantitationType;
import ubic.gemma.model.expression.arrayDesign.ArrayDesign;
import ubic.gemma.model.expression.arrayDesign.ArrayDesignService;
import ubic.gemma.model.expression.bioAssayData.DesignElementDataVector;
import ubic.gemma.model.expression.bioAssayData.ProcessedExpressionDataVector;
import ubic.gemma.model.expression.bioAssayData.ProcessedExpressionDataVectorService;
import ubic.gemma.model.expression.designElement.CompositeSequence;
import ubic.gemma.model.expression.designElement.CompositeSequenceService;
import ubic.gemma.model.expression.designElement.DesignElement;
import ubic.gemma.model.expression.experiment.BioAssaySet;
import ubic.gemma.model.expression.experiment.ExpressionExperiment;
import ubic.gemma.model.genome.Gene;
import cern.colt.list.DoubleArrayList;

/**
 * @author raymond,xwan
 * @version $Id$
 */
public class CorrelationDistCli extends ExpressionExperimentManipulatingCLI {

    private ArrayDesignService adService = null;

    private int binNum = 100;
    private int[][] histogram = null;
    private Map<ExpressionExperiment, Integer> eeIndexMap = null;
    private Collection<ExpressionExperiment> noLinkEEs = null;
    private ByteArrayConverter bac = new ByteArrayConverter();

    @SuppressWarnings("static-access")
    @Override
    protected void buildOptions() {

        Option binNumOption = OptionBuilder.hasArg().withArgName( "Bin Num for Histogram" ).withDescription(
                "Bin Num for Histogram" ).withLongOpt( "binNum" ).create( 'b' );
        addOption( binNumOption );
    }

    @Override
    protected void processOptions() {
        super.processOptions();

        if ( hasOption( 'b' ) ) {
            this.binNum = Integer.valueOf( getOptionValue( 'b' ) );
        }
        adService = ( ArrayDesignService ) this.getBean( "arrayDesignService" );
        noLinkEEs = new HashSet<ExpressionExperiment>();
        processedExpressionDataVectorService = ( ProcessedExpressionDataVectorService ) this
                .getBean( "processedExpressionDataVectorService" );
        compositeSequenceService = ( CompositeSequenceService ) this.getBean( "compositeSequenceService" );
    }

    /**
     * @param needles
     * @param geneIds
     * @return
     */
    boolean known( Collection<Long> needles, Collection<Long> geneIds ) {
        boolean res = false;
        for ( Long id : needles ) {
            if ( geneIds.contains( id ) ) return true;
        }
        return res;
    }

    /**
     * @param source
     * @param target
     * @return
     */
    private double medianEPCorrelation( Collection<ExpressionProfile> source, Collection<ExpressionProfile> target ) {
        DoubleArrayList data = new DoubleArrayList();
        for ( ExpressionProfile ep1 : source ) {
            for ( ExpressionProfile ep2 : target ) {
                if ( ep1.val.length == ep2.val.length
                        && ep1.val.length > GeneEffectSizeCoExpressionAnalyzer.MINIMUM_SAMPLE ) {
                    data.add( CorrelationStats.correl( ep1.val, ep2.val ) );
                }
            }
        }
        data.sort();
        return ( data.size() > 0 ) ? data.get( data.size() / 2 ) : 0.0;
    }

    /**
     * @param genes
     * @return
     */
    private Object[] shuffling( Object[] genes ) {
        Object[] shuffledGenes = new Object[genes.length];
        System.arraycopy( genes, 0, shuffledGenes, 0, genes.length );
        Random random = new Random();
        for ( int i = genes.length - 1; i >= 0; i-- ) {
            int pos = random.nextInt( i + 1 );
            Object tmp = shuffledGenes[pos];
            shuffledGenes[pos] = shuffledGenes[i];
            shuffledGenes[i] = tmp;
        }
        return shuffledGenes;
    }

    ProcessedExpressionDataVectorService processedExpressionDataVectorService;

    /**
     * @param ee
     * @param cs2knowngenes
     * @return
     */
    @SuppressWarnings("unchecked")
    Collection<Double> calculateCorrs( ExpressionExperiment ee, Map<Long, Collection<Long>> cs2knowngenes ) {
        ArrayList<Double> corrs = new ArrayList<Double>();

        Collection<ProcessedExpressionDataVector> processedDataVectors = processedExpressionDataVectorService
                .getProcessedDataVectors( ee );

        // get cs2gene map
        Collection<ArrayDesign> ADs = eeService.getArrayDesignsUsed( ee );
        Collection<Long> csIds = new HashSet<Long>();
        for ( ArrayDesign AD : ADs ) {
            Collection<CompositeSequence> CSs = adService.loadCompositeSequences( AD );
            for ( CompositeSequence CS : CSs ) {
                csIds.add( CS.getId() );
            }
        }
        Map<Long, Collection<Long>> cs2geneMap = getCs2GeneMap( csIds );

        Map<Long, Collection<ExpressionProfile>> geneID2EPs = new HashMap<Long, Collection<ExpressionProfile>>();
        for ( ProcessedExpressionDataVector dedv : processedDataVectors ) {
            Collection<Long> geneIds = cs2geneMap.get( dedv.getDesignElement().getId() );
            for ( Long id : geneIds ) {
                Collection<ExpressionProfile> eps = geneID2EPs.get( id );
                if ( eps == null ) {
                    eps = new HashSet<ExpressionProfile>();
                    geneID2EPs.put( id, eps );
                }
                eps.add( new ExpressionProfile( dedv ) );
            }
        }
        Object[] geneIds = geneID2EPs.keySet().toArray();
        for ( int i = 0; i < 100; i++ ) {
            Object[] shuffledGeneIds = shuffling( geneIds );
            for ( int j = 0; j < shuffledGeneIds.length; j++ ) {
                Collection<ExpressionProfile> source = geneID2EPs.get( geneIds[j] );
                Collection<ExpressionProfile> target = geneID2EPs.get( shuffledGeneIds[j] );
                if ( source != null && target != null ) {
                    double corr = medianEPCorrelation( source, target );
                    corrs.add( corr );
                }
            }
        }
        return corrs;
    }

    CompositeSequenceService compositeSequenceService;

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
     * @param ee
     * @param geneIds
     */
    @SuppressWarnings("unchecked")
    void fillHistogram( ExpressionExperiment ee, Collection<Long> geneIds ) {
        int halfBin = binNum / 2;
        Collection<Long> csIds = new HashSet<Long>();
        Collection<DesignElement> allCSs = new HashSet<DesignElement>();
        Collection<ArrayDesign> ads = eeService.getArrayDesignsUsed( ee );
        for ( ArrayDesign ad : ads ) {
            allCSs.addAll( adService.loadCompositeSequences( ad ) );
        }
        for ( DesignElement cs : allCSs ) {
            csIds.add( cs.getId() );
        }

        Map<Long, Collection<Long>> cs2genes = getCs2GeneMap( csIds );
        Map<Long, Collection<Long>> cs2knowngenes = new HashMap<Long, Collection<Long>>();
        for ( Long csId : cs2genes.keySet() ) {
            Collection<Long> mappedGeneIds = cs2genes.get( csId );
            if ( known( mappedGeneIds, geneIds ) ) cs2knowngenes.put( csId, mappedGeneIds );
        }
        Collection<Double> corrs = calculateCorrs( ee, cs2knowngenes );
        int eeIndex = eeIndexMap.get( ee );
        for ( Double corr : corrs ) {
            int bin = Math.min( ( int ) ( ( 1.0 + corr ) * halfBin ), binNum - 1 );
            histogram[eeIndex][bin]++;
        }
    }

    /**
     * 
     */
    void saveHistogram() {
        try {
            FileWriter out = new FileWriter( new File( "correlationDist.txt" ) );
            List<String> rowLabels = new ArrayList<String>();
            List<String> colLabels = new ArrayList<String>();
            for ( int i = 0; i < binNum; i++ ) {
                out.write( "\t" + i );
                colLabels.add( Integer.toString( i ) );
            }
            out.write( "\n" );
            // double culmulative = 0.0;
            int culmulatives[] = new int[histogram.length];
            for ( int i = 0; i < histogram.length; i++ ) {
                for ( int j = 0; j < binNum; j++ ) {
                    culmulatives[i] = culmulatives[i] + histogram[i][j];
                }
            }
            for ( ExpressionExperiment ee : eeIndexMap.keySet() ) {
                if ( noLinkEEs.contains( ee ) ) continue;
                rowLabels.add( ee.getShortName() );
            }
            double data[][] = new double[histogram.length - noLinkEEs.size()][binNum];
            int dataIndex = 0;
            for ( ExpressionExperiment ee : eeIndexMap.keySet() ) {
                if ( noLinkEEs.contains( ee ) ) continue;
                out.write( eeService.getTaxon( ee.getId() ).getCommonName() + ee.getShortName() );
                int eeIndex = eeIndexMap.get( ee );
                for ( int j = 0; j < binNum; j++ ) {
                    data[dataIndex][j] = ( double ) histogram[eeIndex][j] / ( double ) culmulatives[eeIndex];
                    out.write( "\t" + data[dataIndex][j] );
                }
                out.write( "\n" );
                log.info( ee.getShortName() + "---->" + culmulatives[eeIndex] );
                dataIndex++;
            }
            DoubleMatrix<String, String> dataMatrix = new DenseDoubleMatrix<String, String>( data );
            dataMatrix.setRowNames( rowLabels );
            dataMatrix.setColumnNames( colLabels );

            ColorMatrix dataColorMatrix = new ColorMatrix( dataMatrix );
            // dataColorMatrix.setColorMap( ColorMap.GREENRED_COLORMAP );
            dataColorMatrix.setColorMap( ColorMap.BLACKBODY_COLORMAP );
            MatrixDisplay dataMatrixDisplay = new MatrixDisplay( dataColorMatrix );
            dataMatrixDisplay.saveImage( "correlationDist.png", true );

            out.write( "\n" );
            out.close();
        } catch ( Exception e ) {
            e.printStackTrace();
        }
    }

    @SuppressWarnings("unchecked")
    @Override
    protected Exception doWork( String[] args ) {
        Exception err = processCommandLine( "Correlation Distribution ", args );
        if ( err != null ) {
            return err;
        }

        try {

            histogram = new int[expressionExperiments.size()][binNum];
            eeIndexMap = new HashMap<ExpressionExperiment, Integer>();
            int index = 0;
            for ( BioAssaySet bas : expressionExperiments ) {
                ExpressionExperiment ee = ( ExpressionExperiment ) bas;
                eeIndexMap.put( ee, index );
                index++;
            }
            Collection<Gene> genes = geneService.loadKnownGenes( taxon );
            Collection<Long> geneIds = new HashSet<Long>();
            for ( Gene gene : genes )
                geneIds.add( gene.getId() );
            for ( BioAssaySet bas : expressionExperiments ) {
                ExpressionExperiment ee = ( ExpressionExperiment ) bas;
                fillHistogram( ee, geneIds );
            }
            saveHistogram();
            return null;
        } catch ( Exception e ) {
            return e;
        }
    }

    /**
     */
    public class ExpressionProfile {
        DesignElementDataVector dedv = null;

        double[] val = null;

        long id;

        /**
         * Construct an ExpressionProfile from the specified DesignElementDataVector
         * 
         * @param dedv - vector to convert
         */
        public ExpressionProfile( DesignElementDataVector dedv ) {
            this.dedv = dedv;
            this.id = dedv.getId();
            byte[] bytes = dedv.getData();
            val = bac.byteArrayToDoubles( bytes );
        }

        /**
         * Get the ID of the vector
         * 
         * @return the vector ID
         */
        public long getId() {
            return this.id;
        }
    }

    /**
     * @param args
     */
    public static void main( String[] args ) {
        CorrelationDistCli corrDist = new CorrelationDistCli();
        StopWatch watch = new StopWatch();
        watch.start();
        try {
            Exception ex = corrDist.doWork( args );
            if ( ex != null ) {
                ex.printStackTrace();
            }
            watch.stop();
            log.info( watch.getTime() / 1000 );
        } catch ( Exception e ) {
            throw new RuntimeException( e );
        }
    }

}
