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
package chibi.gemmaanalysis.cli.deprecated;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.time.StopWatch;

import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.graphics.ColorMap;
import ubic.basecode.graphics.ColorMatrix;
import ubic.basecode.graphics.MatrixDisplay;
import ubic.gemma.expression.experiment.service.ExpressionExperimentService;
import ubic.gemma.genome.gene.service.GeneService;
import ubic.gemma.genome.taxon.service.TaxonService;
import ubic.gemma.model.association.coexpression.Probe2ProbeCoexpressionService;
import ubic.gemma.model.association.coexpression.ProbeLink;
import ubic.gemma.model.expression.designElement.CompositeSequence;
import ubic.gemma.model.expression.designElement.CompositeSequenceService;
import ubic.gemma.model.expression.experiment.ExpressionExperiment;
import ubic.gemma.model.genome.Gene;
import ubic.gemma.model.genome.Taxon;
import ubic.gemma.ontology.providers.GeneOntologyService;
import ubic.gemma.util.AbstractSpringAwareCLI;

/**
 * @author xwan
 * @version $Id$
 * @deprecated because we don't know what this is for. see LinkGOStatsCli.
 */
@Deprecated
public class LinkGOAnalysisCli extends AbstractSpringAwareCLI {

    private final static int GO_MAXIMUM_COUNT = 50;
    private final static int ITERATION_NUM = 1;

    /**
     * @param args
     */
    public static void main( String[] args ) {
        LinkGOAnalysisCli goAnalysis = new LinkGOAnalysisCli();
        StopWatch watch = new StopWatch();
        watch.start();
        try {
            Exception ex = goAnalysis.doWork( args );
            if ( ex != null ) {
                ex.printStackTrace();
            }
            watch.stop();
            log.info( watch.getTime() / 1000 );
        } catch ( Exception e ) {
            throw new RuntimeException( e );
        }
    }

    private Probe2ProbeCoexpressionService p2pService = null;
    private ExpressionExperimentService eeService = null;
    private GeneService geneService = null;
    private GeneOntologyService goService;

    private TaxonService taxonService;
    private String taxonName = "mouse";
    private String eeNameFile = null;
    private Map<ExpressionExperiment, Integer> eeIndexMap = new HashMap<ExpressionExperiment, Integer>();

    private Map<Long, Gene> geneMap = new HashMap<Long, Gene>();
    private int[][] realStats = null;
    private int[][] simulatedStats = null;
    private Collection<ExpressionExperiment> noLinkEEs = null;

    private Collection<Gene> coveredGenes = null;

    CompositeSequenceService compositeSequenceService;

    @SuppressWarnings("static-access")
    @Override
    protected void buildOptions() {
        Option taxonOption = OptionBuilder.hasArg().isRequired().withArgName( "Taxon" )
                .withDescription( "the taxon name" ).withLongOpt( "Taxon" ).create( 't' );
        addOption( taxonOption );
        Option eeNameFileo = OptionBuilder.hasArg().withArgName( "File having Expression Experiment Names" )
                .withDescription( "File having Expression Experiment Names" ).withLongOpt( "eeFileName" ).create( 'f' );
        addOption( eeNameFileo );

    }

    @Override
    protected Exception doWork( String[] args ) {
        Exception err = processCommandLine( "Shuffle Links ", args );
        if ( err != null ) {
            return err;
        }

        Taxon taxon = taxonService.findByCommonName( taxonName );
        compositeSequenceService = this.getBean( CompositeSequenceService.class );
        Collection<ExpressionExperiment> ees = eeService.findByTaxon( taxon );
        Collection<ExpressionExperiment> eeCandidates = getCandidateEE( this.eeNameFile, ees );
        Collection<Gene> allGenes = geneService.loadAll( taxon );

        log.info( "Load " + allGenes.size() + " genes" );

        for ( Gene gene : allGenes ) {
            geneMap.put( gene.getId(), gene );
        }
        int index = 0;
        for ( ExpressionExperiment ee : eeCandidates ) {
            eeIndexMap.put( ee, index );
            index++;
        }
        realStats = new int[index][GO_MAXIMUM_COUNT];
        simulatedStats = new int[index][GO_MAXIMUM_COUNT];

        for ( ExpressionExperiment ee : eeCandidates ) {
            log.info( "Shuffling " + ee.getShortName() );
            Collection<ProbeLink> links = p2pService.getProbeCoExpression( ee, this.taxonName, false );
            coveredGenes = new HashSet<Gene>();

            if ( links == null || links.size() == 0 ) {
                noLinkEEs.add( ee );
                continue;
            }
            StopWatch watch = new StopWatch();
            watch.start();
            counting( realStats[eeIndexMap.get( ee )], links );
            watch.stop();
            log.info( "Time spent at first time for GO search " + watch.getTime() );
            Gene[] genes = new Gene[coveredGenes.size()];
            for ( Gene gene : coveredGenes )
                genes[index++] = gene;

            for ( int i = 0; i < ITERATION_NUM; i++ ) {
                shuffleLinks( links );
                watch.reset();
                watch.start();
                counting( simulatedStats[eeIndexMap.get( ee )], links );
                watch.stop();
                log.info( "Time spent at second time for GO search " + watch.getTime() );

            }
        }
        output();
        return null;
    }

    @Override
    protected void processOptions() {
        super.processOptions();
        if ( hasOption( 't' ) ) {
            this.taxonName = getOptionValue( 't' );
        }
        if ( hasOption( 'f' ) ) {
            this.eeNameFile = getOptionValue( 'f' );
        }

        p2pService = this.getBean( Probe2ProbeCoexpressionService.class );
        eeService = this.getBean( ExpressionExperimentService.class );
        geneService = this.getBean( GeneService.class );
        noLinkEEs = new HashSet<ExpressionExperiment>();
        goService = this.getBean( GeneOntologyService.class );
        taxonService = this.getBean( TaxonService.class );
    }

    private void counting( int[] stats, Collection<ProbeLink> links ) {
        Collection<Long> csIds = new HashSet<Long>();
        for ( ProbeLink link : links ) {
            csIds.add( link.getFirstDesignElementId() );
            csIds.add( link.getSecondDesignElementId() );
        }
        Map<Long, Collection<Long>> cs2genes = getCs2GeneMap( csIds );
        for ( ProbeLink link : links ) {
            if ( link.getFirstDesignElementId() == link.getSecondDesignElementId() ) continue;
            Collection<Long> firstGeneIds = cs2genes.get( link.getFirstDesignElementId() );
            Collection<Long> secondGeneIds = cs2genes.get( link.getSecondDesignElementId() );
            if ( firstGeneIds == null || secondGeneIds == null ) {
                continue;
            }

            for ( Long firstGeneId : firstGeneIds ) {
                for ( Long secondGeneId : secondGeneIds ) {
                    firstGeneId = firstGeneIds.iterator().next();
                    secondGeneId = secondGeneIds.iterator().next();
                    Gene gene1 = geneMap.get( firstGeneId );
                    Gene gene2 = geneMap.get( secondGeneId );
                    if ( gene1 == null || gene2 == null ) {
                        log.info( "Wrong setting for gene" + firstGeneId + "\t" + secondGeneId );
                        continue;
                    }
                    int goOverlap = goService.calculateGoTermOverlap( gene1, gene2 ).size();
                    if ( goOverlap >= GO_MAXIMUM_COUNT )
                        stats[GO_MAXIMUM_COUNT - 1]++;
                    else
                        stats[goOverlap]++;

                }
            }
        }
    }

    private Collection<ExpressionExperiment> getCandidateEE( String fileName, Collection<ExpressionExperiment> ees ) {
        if ( fileName == null ) return ees;
        Collection<ExpressionExperiment> candidates = new HashSet<ExpressionExperiment>();
        Collection<String> eeNames = new HashSet<String>();
        try {
            InputStream is = new FileInputStream( fileName );
            BufferedReader br = new BufferedReader( new InputStreamReader( is ) );
            String shortName = null;
            while ( ( shortName = br.readLine() ) != null ) {
                if ( StringUtils.isBlank( shortName ) ) continue;
                eeNames.add( shortName.trim().toUpperCase() );
            }
        } catch ( Exception e ) {
            e.printStackTrace();
            return candidates;
        }
        for ( ExpressionExperiment ee : ees ) {
            String shortName = ee.getShortName();
            if ( eeNames.contains( shortName.trim().toUpperCase() ) ) candidates.add( ee );
        }
        return candidates;
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

    private void output() {
        try {
            FileWriter out = new FileWriter( new File( "analysis.txt" ) );
            out.write( "Real GO Stats:\n" );
            for ( int i = 0; i < realStats.length; i++ ) {
                for ( int j = 0; j < realStats[i].length; j++ ) {
                    out.write( realStats[i][j] + "\t" );
                }
                out.write( "\n" );
            }
            out.write( "Simulate GO Stats:\n" );
            for ( int i = 0; i < realStats.length; i++ ) {
                for ( int j = 0; j < realStats[i].length; j++ ) {
                    out.write( realStats[i][j] + "\t" );
                }
                out.write( "\n" );
            }

            output( realStats, "realGODist.png" );
            output( simulatedStats, "simulateGODist.png" );

            out.close();
        } catch ( Exception e ) {
            e.printStackTrace();
        }
    }

    private void output( int[][] stats, String imageName ) throws Exception {
        List<String> rowLabels = new ArrayList<String>();
        List<String> colLabels = new ArrayList<String>();
        for ( int i = 0; i < GO_MAXIMUM_COUNT; i++ ) {
            colLabels.add( Integer.toString( i ) );
        }
        // double culmulative = 0.0;
        int culmulatives[] = new int[stats.length];
        for ( int i = 0; i < stats.length; i++ ) {
            for ( int j = 0; j < GO_MAXIMUM_COUNT; j++ ) {
                culmulatives[i] = culmulatives[i] + stats[i][j];
            }
        }
        for ( ExpressionExperiment ee : eeIndexMap.keySet() ) {
            if ( noLinkEEs.contains( ee ) ) continue;
            rowLabels.add( ee.getShortName() );
        }
        double data[][] = new double[stats.length - noLinkEEs.size()][GO_MAXIMUM_COUNT];
        int dataIndex = 0;
        for ( ExpressionExperiment ee : eeIndexMap.keySet() ) {
            if ( noLinkEEs.contains( ee ) ) continue;
            int eeIndex = eeIndexMap.get( ee );
            for ( int j = 0; j < GO_MAXIMUM_COUNT; j++ ) {
                data[dataIndex][j] = ( double ) stats[eeIndex][j] / ( double ) culmulatives[eeIndex];
            }
            dataIndex++;
        }
        DoubleMatrix<String, String> dataMatrix = new DenseDoubleMatrix<String, String>( data );
        dataMatrix.setRowNames( rowLabels );
        dataMatrix.setColumnNames( colLabels );

        ColorMatrix<String, String> dataColorMatrix = new ColorMatrix<String, String>( dataMatrix );
        // dataColorMatrix.setColorMap( ColorMap.GREENRED_COLORMAP );
        dataColorMatrix.setColorMap( ColorMap.BLACKBODY_COLORMAP );
        MatrixDisplay<String, String> dataMatrixDisplay = new MatrixDisplay<String, String>( dataColorMatrix );
        dataMatrixDisplay.saveImage( imageName, true, false );

    }

    private void shuffleLinks( Collection<ProbeLink> links ) {
        // Do shuffling
        Random random = new Random();
        Object[] linksInArray = links.toArray();
        for ( int i = linksInArray.length - 1; i >= 0; i-- ) {
            int pos = random.nextInt( i + 1 );
            Long tmpId = ( ( ProbeLink ) linksInArray[pos] ).getSecondDesignElementId();
            ( ( ProbeLink ) linksInArray[pos] ).setSecondDesignElementId( ( ( ProbeLink ) linksInArray[i] )
                    .getSecondDesignElementId() );
            ( ( ProbeLink ) linksInArray[i] ).setSecondDesignElementId( tmpId );
        }
    }

}
