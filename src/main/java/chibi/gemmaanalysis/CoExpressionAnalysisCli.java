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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.lang.time.StopWatch;

import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.dataStructure.matrix.StringMatrix;
import ubic.basecode.graphics.ColorMap;
import ubic.basecode.graphics.ColorMatrix;
import ubic.basecode.graphics.MatrixDisplay;
import ubic.basecode.io.reader.StringMatrixReader;
import chibi.gemmaanalysis.GeneEffectSizeCoExpressionAnalyzer;
import ubic.gemma.analysis.expression.coexpression.ProbeLinkCoexpressionAnalyzer;
import ubic.gemma.model.analysis.expression.coexpression.CoexpressionCollectionValueObject;
import ubic.gemma.model.expression.arrayDesign.ArrayDesign;
import ubic.gemma.model.expression.arrayDesign.ArrayDesignService;
import ubic.gemma.model.expression.bioAssayData.DesignElementDataVector;
import ubic.gemma.model.expression.bioAssayData.DoubleVectorValueObject;
import ubic.gemma.model.expression.bioAssayData.ProcessedExpressionDataVector;
import ubic.gemma.model.expression.bioAssayData.ProcessedExpressionDataVectorService;
import ubic.gemma.model.expression.designElement.CompositeSequence;
import ubic.gemma.model.expression.designElement.CompositeSequenceService;
import ubic.gemma.model.expression.experiment.ExpressionExperiment;
import ubic.gemma.model.expression.experiment.ExpressionExperimentService;
import ubic.gemma.model.genome.Gene;
import ubic.gemma.model.genome.Taxon;
import ubic.gemma.model.genome.TaxonService;
import ubic.gemma.genome.gene.service.GeneService;
import ubic.gemma.util.AbstractSpringAwareCLI;

/**
 * "Global" coexpression analysis: for a set of coexpressed genes, get all the data for all the genes it is coexpressed
 * with and construct a summary matrix.
 * 
 * @author xwan
 * @version $Id$
 */
public class CoExpressionAnalysisCli extends AbstractSpringAwareCLI {

    private ProcessedExpressionDataVectorService dedvService;
    private ExpressionExperimentService eeService;
    private GeneService geneService;
    private String geneList = null;
    private String taxonName = null;
    private int stringency = 3;
    private static String DIVIDOR = "-----";
    private String experimentListFile;
    private String eeExcludeFile;

    @SuppressWarnings("static-access")
    @Override
    protected void buildOptions() {
        Option geneFileOption = OptionBuilder.hasArg().isRequired().withArgName( "geneFile" )
                .withDescription( "Short names of the genes to analyze" ).withLongOpt( "geneFile" ).create( 'g' );
        addOption( geneFileOption );
        Option taxonOption = OptionBuilder.hasArg().isRequired().withArgName( "Taxon" )
                .withDescription( "the taxon of the genes to analyze" ).withLongOpt( "Taxon" ).create( 't' );
        addOption( taxonOption );

        Option stringencyFileOption = OptionBuilder.hasArg().withArgName( "stringency" )
                .withDescription( "The minimum support for links to be selected (Default 3)" )
                .withLongOpt( "stringency" ).create( 's' );
        addOption( stringencyFileOption );

        Option eeListOption = OptionBuilder.hasArg().withArgName( "Expression experiment list file" )
                .withDescription( "File with list of short names of expression experiments to use" )
                .withLongOpt( "eeListfile" ).create( 'f' );
        addOption( eeListOption );

        Option eeExcludeListFile = OptionBuilder.hasArg().withArgName( "Expression experiment exclude list file" )
                .withDescription( "File with list of short names of expression experiments to exclude" )
                .withLongOpt( "eeExcludefile" ).create( 'x' );
        addOption( eeExcludeListFile );
    }

    /**
     * 
     */
    @Override
    protected void processOptions() {
        super.processOptions();
        if ( hasOption( 'g' ) ) {
            this.geneList = getOptionValue( 'g' );
        }
        if ( hasOption( 't' ) ) {
            this.taxonName = getOptionValue( 't' );
        }

        if ( hasOption( 's' ) ) {
            this.stringency = Integer.parseInt( getOptionValue( 's' ) );
        }
        if ( hasOption( 'f' ) ) {
            this.experimentListFile = getOptionValue( 'f' );
        }
        if ( hasOption( 'x' ) ) {
            this.eeExcludeFile = getOptionValue( 'x' );
        }
        dedvService = ( ProcessedExpressionDataVectorService ) this.getBean( "processedExpressionDataVectorService" );
        eeService = ( ExpressionExperimentService ) this.getBean( "expressionExperimentService" );
        geneService = ( GeneService ) this.getBean( "geneService" );
        probeLinkCoexpressionAnalyzer = ( ProbeLinkCoexpressionAnalyzer ) this
                .getBean( "probeLinkCoexpressionAnalyzer" );
        processedExpressionDataVectorService = ( ProcessedExpressionDataVectorService ) this
                .getBean( "processedExpressionDataVectorService" );
        adService = ( ArrayDesignService ) getBean( "arrayDesignService" );
    }

    ProbeLinkCoexpressionAnalyzer probeLinkCoexpressionAnalyzer;

    /**
     * @param geneService
     * @param geneNames
     * @param taxon
     * @return
     */
    private Collection<Gene> getGenes( Object[] geneNames, Taxon taxon ) {
        HashSet<Gene> genes = new HashSet<Gene>();
        for ( int i = 0; i < geneNames.length; i++ ) {
            Gene gene = getGene( ( String ) geneNames[i], taxon );
            if ( gene != null ) genes.add( gene );
        }
        return genes;
    }

    /**
     * @param geneService
     * @param geneName
     * @param taxon
     * @return
     */
    private Gene getGene( String geneName, Taxon taxon ) {
        Gene gene = Gene.Factory.newInstance();
        gene.setOfficialSymbol( geneName.trim() );
        gene.setTaxon( taxon );
        gene = geneService.find( gene );
        if ( gene == null ) {
            log.info( "Can't Load gene " + geneName );
        }
        return gene;
    }

    /**
     * @param name
     * @return
     */
    private Taxon getTaxon( String name ) {
        Taxon taxon = Taxon.Factory.newInstance();
        taxon.setCommonName( name );
        TaxonService taxonService = ( TaxonService ) this.getBean( "taxonService" );
        taxon = taxonService.find( taxon );
        if ( taxon == null ) {
            log.info( "NO Taxon found!" );
        }
        return taxon;
    }

    private ProcessedExpressionDataVectorService processedExpressionDataVectorService;

    Collection<Gene> getCoExpressedGenes( Collection<Gene> queryGenes ) {
        Set<Gene> coExpressedGenes = new HashSet<Gene>();
        Collection<Long> geneIds = new HashSet<Long>();
        for ( Gene gene : queryGenes ) {
            log.info( "Get co-expressed genes for " + gene.getName() );
            CoexpressionCollectionValueObject coexpressed = probeLinkCoexpressionAnalyzer.linkAnalysis( gene, null,
                    this.stringency, 0 );
            Map<Long, Collection<Long>> geneEEMap = coexpressed.getKnownGeneCoexpression()
                    .getExpressionExperimentsWithSpecificProbeForCoexpressedGenes();
            for ( Long geneId : geneEEMap.keySet() ) {
                Collection<Long> ees = geneEEMap.get( geneId );
                if ( ees.size() >= this.stringency ) geneIds.add( geneId );
            }
        }
        log.info( " Got " + geneIds.size() + " genes for the CoExpression analysis" );
        if ( geneIds.size() > 0 ) coExpressedGenes.addAll( geneService.loadMultiple( geneIds ) );
        return coExpressedGenes;
    }

    private ArrayDesignService adService = null;

    /**
     * @param cs2gene
     * @param ee
     * @return
     */
    private Map<DoubleVectorValueObject, Collection<Long>> getDesignElementDataVector( ExpressionExperiment ee ) {

        Collection<ProcessedExpressionDataVector> dedvs = processedExpressionDataVectorService
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

        Map<DoubleVectorValueObject, Collection<Long>> dedv2genes = new HashMap<DoubleVectorValueObject, Collection<Long>>();
        for ( DesignElementDataVector dedv : dedvs ) {
            dedv2genes.put( new DoubleVectorValueObject( dedv ), cs2geneMap.get( dedv.getDesignElement().getId() ) );
        }
        return dedv2genes;
    }

    // For small number of queryGenes.
    Map<DoubleVectorValueObject, Collection<Long>> getDedv2GenesMap( Collection<Gene> queryGenes,
            Collection<Gene> coExpressedGenes, Collection<ExpressionExperiment> allEEs ) {
        StopWatch qWatch = new StopWatch();
        qWatch.start();
        Map<DoubleVectorValueObject, Collection<Long>> dedv2queryGenes = new HashMap<DoubleVectorValueObject, Collection<Long>>();
        Map<DoubleVectorValueObject, Collection<Long>> dedv2coExpressedGenes = new HashMap<DoubleVectorValueObject, Collection<Long>>();
        // ArrayList<Object> saved = new ArrayList<Object>();
        log.info( "Start the Query for " + queryGenes.size() + " query genes" );
        Collection<DoubleVectorValueObject> processedDataArrays = dedvService.getProcessedDataArrays( allEEs,
                queryGenes );
        for ( DoubleVectorValueObject doubleVectorValueObject : processedDataArrays ) {
            Collection<Long> gids = new HashSet<Long>();
            for ( Gene g : doubleVectorValueObject.getGenes() ) {
                gids.add( g.getId() );
            }
            dedv2queryGenes.put( doubleVectorValueObject, gids );
        }

        Collection<ExpressionExperiment> ees = new HashSet<ExpressionExperiment>();
        Collection<Long> eeIds = new HashSet<Long>();
        log.info( "loading designElementDataVector from expression experiments" );
        int count = 0;
        for ( DoubleVectorValueObject dedv : dedv2queryGenes.keySet() ) {
            ExpressionExperiment ee = dedv.getExpressionExperiment();

            if ( !eeIds.contains( ee.getId() ) ) {
                Map<DoubleVectorValueObject, Collection<Long>> dedvs = getDesignElementDataVector( ee );
                dedv2coExpressedGenes.putAll( dedvs );
                count = count + dedvs.keySet().size();
                eeIds.add( ee.getId() );
            }
        }
        for ( ExpressionExperiment ee : allEEs ) {
            if ( eeIds.contains( ee.getId() ) ) {
                ees.add( ee );
            }
        }
        allEEs.clear();
        allEEs.addAll( ees );
        log.info( "Get " + allEEs.size() + " expression experiments for analysis" );
        dedv2coExpressedGenes.putAll( dedv2queryGenes );
        qWatch.stop();
        log.info( "Query takes " + qWatch.getTime() + " to get " + dedv2coExpressedGenes.keySet().size() + ":" + count
                + " DEDVs for " + ( coExpressedGenes.size() + queryGenes.size() ) + " genes" );
        return dedv2coExpressedGenes;
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

    /*
     * (non-Javadoc)
     * 
     * @see ubic.gemma.util.AbstractCLI#doWork(java.lang.String[])
     */
    @SuppressWarnings("unchecked")
    @Override
    protected Exception doWork( String[] args ) {
        Exception err = processCommandLine( "CoExpression Analysis", args );
        if ( err != null ) {
            return err;
        }

        Taxon taxon = getTaxon();
        Collection<ExpressionExperiment> allEEs = null;
        if ( this.experimentListFile != null ) {
            try {
                allEEs = readExpressionExperimentListFile( this.experimentListFile );
            } catch ( IOException e ) {
                return e;
            }
        } else if ( taxon != null ) {
            allEEs = eeService.findByTaxon( taxon );
            if ( this.eeExcludeFile != null ) {
                try {
                    Collection<ExpressionExperiment> eesToExclude = readExpressionExperimentListFile( this.eeExcludeFile );
                    allEEs.removeAll( eesToExclude );
                } catch ( IOException e ) {
                    return e;
                }
            }
        } else {
            log.error( "You must provide either the taxon or a list of expression experiments in a file" );
            bail( ErrorCode.MISSING_OPTION );
        }

        Collection<String> queryGeneNames = new HashSet<String>();
        Collection<String> coExpressedGeneNames = new HashSet<String>();
        boolean readingQueryGene = true;
        /*
         * The gene input file could contain query genes and co-expressed genes divided by DIVIDOR; If user doesn't
         * provide the co-expressed genes, then use the service to find the co-expressed genes in database.
         */
        try {
            InputStream is = new FileInputStream( this.geneList );
            BufferedReader br = new BufferedReader( new InputStreamReader( is ) );
            String shortName = null;
            while ( ( shortName = br.readLine() ) != null ) {
                if ( StringUtils.isBlank( shortName ) ) continue;
                if ( shortName.trim().contains( DIVIDOR ) ) {
                    readingQueryGene = false;
                    continue;
                }
                if ( readingQueryGene ) {
                    queryGeneNames.add( shortName.trim() );
                } else {
                    coExpressedGeneNames.add( shortName.trim() );
                }
            }
        } catch ( Exception e ) {
            return e;
        }

        if ( queryGeneNames.isEmpty() ) {
            log.info( "No gene is read from the input file" );
            return null;
        }
        Collection<Gene> queryGenes = this.getGenes( queryGeneNames.toArray(), taxon );
        if ( queryGenes.isEmpty() ) {
            log.info( "Can't load any of genes" + queryGeneNames );
            return null;
        }
        Collection<Gene> coExpressedGenes = null;
        if ( !coExpressedGeneNames.isEmpty() ) {
            coExpressedGenes = this.getGenes( coExpressedGeneNames.toArray(), taxon );
        } else {
            // coexpressed genes using vote count
            coExpressedGenes = this.getCoExpressedGenes( queryGenes );
        }
        log.info( "Start the Query for " + queryGenes.size() + " genes" );
        Map<DoubleVectorValueObject, Collection<Long>> dedv2genes = getDedv2GenesMap( queryGenes, coExpressedGenes,
                allEEs );

        if ( dedv2genes.isEmpty() || queryGenes.isEmpty() || coExpressedGenes.isEmpty() || allEEs == null
                || allEEs.isEmpty() ) return null;
        GeneEffectSizeCoExpressionAnalyzer coExpression = new GeneEffectSizeCoExpressionAnalyzer( queryGenes,
                coExpressedGenes, new HashSet( allEEs ) );

        coExpression.setDedv2Genes( dedv2genes );
        coExpression.setExpressionExperimentService( eeService );
        coExpression.analyze( dedv2genes.keySet() );

        for ( Gene gene : queryGenes ) {
            try {
                makeClusterGrams( coExpression, gene );
            } catch ( Exception e ) {
                return e;
            }
        }
        coExpression.calculateMatrixEffectSize();
        return null;
    }

    /**
     * @return
     */
    private Taxon getTaxon() {
        Taxon taxon = getTaxon( this.taxonName );
        if ( taxon == null ) {
            log.error( "No taxon is found " + this.taxonName );
            bail( ErrorCode.INVALID_OPTION );
        }
        return taxon;
    }

    /**
     * @param coExpression
     * @throws FileNotFoundException
     * @throws IOException
     * @throws InterruptedException
     */
    private void makeClusterGrams( GeneEffectSizeCoExpressionAnalyzer coExpression, Gene inputGene )
            throws FileNotFoundException, IOException, InterruptedException {
        String filebaseName = inputGene.getOfficialSymbol() + "_coexp";
        // Generate the data file for Cluster3
        PrintStream output = new PrintStream( new FileOutputStream( new File( filebaseName ) ) );
        double presencePercent = 0.8;
        coExpression.output( output, presencePercent, inputGene );
        output.close();

        // Running Cluster3 to geneate .cdt file
        Runtime rt = Runtime.getRuntime();
        Process clearOldFiles = rt.exec( "rm *.cdt -f" );
        clearOldFiles.waitFor();

        String clusterCmd = "eisen-cluster";
        String commonOptions = "-g 2 -e 2 -m a"; // pearson, average linkage.
        String cmdToRun = clusterCmd + " -f " + filebaseName + " " + commonOptions;
        log.info( "Running: " + cmdToRun );
        Process cluster = rt.exec( cmdToRun );
        cluster.waitFor();

        DoubleMatrix<String, String> dataMatrix = getClusteredMatrix( filebaseName );

        // Get the rank Matrix
        DoubleMatrix<String, String> rankMatrix = coExpression.getRankMatrix( dataMatrix );

        // generate the png figures
        ColorMatrix dataColorMatrix = new ColorMatrix( dataMatrix );
        dataColorMatrix.setColorMap( ColorMap.GREENRED_COLORMAP );
        ColorMatrix rankColorMatrix = new ColorMatrix( rankMatrix );
        rankColorMatrix.setColorMap( ColorMap.GREENRED_COLORMAP );

        MatrixDisplay dataMatrixDisplay = new MatrixDisplay( dataColorMatrix );
        MatrixDisplay rankMatrixDisplay = new MatrixDisplay( rankColorMatrix );

        dataMatrixDisplay.saveImage( filebaseName + ".png", true, false );
        rankMatrixDisplay.saveImage( filebaseName + ".ranks.png", true, false );
    }

    /**
     * @return
     * @throws IOException
     */
    private DoubleMatrix<String, String> getClusteredMatrix( String baseName ) throws IOException {
        // Read the generated file into a String Matrix
        StringMatrixReader mReader = new StringMatrixReader();

        String CDTMatrixFile = baseName;
        StringMatrix<String, String> cdtMatrix = mReader.read( CDTMatrixFile + ".cdt" );

        // Read String Matrix and convert into DenseDoubleMatrix
        int extra_rows = 2, extra_cols = 3;
        double[][] data = new double[cdtMatrix.rows() - extra_rows][];
        List<String> rowLabels = new ArrayList<String>();
        List<String> colLabels = new ArrayList<String>();

        List<String> colNames = cdtMatrix.getColNames();
        for ( int i = extra_cols; i < colNames.size(); i++ )
            colLabels.add( colNames.get( i ) );

        int rowIndex = 0;
        for ( int i = extra_rows; i < cdtMatrix.rows(); i++ ) {
            Object row[] = cdtMatrix.getRow( i );
            rowLabels.add( ( String ) row[0] );
            data[rowIndex] = new double[row.length - extra_cols];
            for ( int j = extra_cols; j < row.length; j++ )
                try {
                    data[rowIndex][j - extra_cols] = Double.valueOf( ( String ) row[j] );
                } catch ( Exception e ) {
                    data[rowIndex][j - extra_cols] = Double.NaN;
                    continue;
                }
            rowIndex++;
        }
        DoubleMatrix<String, String> dataMatrix = new DenseDoubleMatrix<String, String>( data );
        dataMatrix.setRowNames( rowLabels );
        dataMatrix.setColumnNames( colLabels );
        return dataMatrix;
    }

    /**
     * @param args
     */
    public static void main( String[] args ) {
        CoExpressionAnalysisCli analysis = new CoExpressionAnalysisCli();
        StopWatch watch = new StopWatch();
        watch.start();
        try {
            Exception ex = analysis.doWork( args );
            if ( ex != null ) {
                ex.printStackTrace();
            }
            watch.stop();
            log.info( "Time elapsed: " + watch.getTime() );
            System.exit( 0 );
        } catch ( Exception e ) {
            throw new RuntimeException( e );
        }
    }

    /**
     * FIXME duplicated code from AbstractedGeneExpressionExperimentManipulatingCli
     * 
     * @param fileName
     * @return
     * @throws IOException
     */
    private Collection<String> readExpressionExperimentListFileToStrings( String fileName ) throws IOException {
        Collection<String> eeNames = new HashSet<String>();
        BufferedReader in = new BufferedReader( new FileReader( fileName ) );
        while ( in.ready() ) {
            String eeName = in.readLine().trim();
            if ( eeName.startsWith( "#" ) ) {
                continue;
            }
            eeNames.add( eeName );
        }
        return eeNames;
    }

    /**
     * Load expression experiments based on a list of short names in a file. FIXME duplicated code from
     * AbstractedGeneExpressionExperimentManipulatingCli
     * 
     * @param fileName
     * @return
     * @throws IOException
     */
    public Collection<ExpressionExperiment> readExpressionExperimentListFile( String fileName ) throws IOException {
        Collection<ExpressionExperiment> ees = new HashSet<ExpressionExperiment>();
        for ( String eeName : readExpressionExperimentListFileToStrings( fileName ) ) {
            ExpressionExperiment ee = eeService.findByShortName( eeName );
            if ( ee == null ) {
                log.error( "No experiment " + eeName + " found" );
                continue;
            }
            ees.add( ee );
        }
        return ees;
    }

}
