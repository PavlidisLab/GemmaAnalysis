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
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.io.Serializable;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.lang.time.StopWatch;

import cern.jet.stat.Descriptive;
import cern.colt.list.DoubleArrayList;

import ubic.basecode.dataStructure.matrix.DoubleMatrixNamed;
import ubic.basecode.dataStructure.matrix.SparseRaggedDoubleMatrix2DNamed;
import ubic.basecode.math.RandomChooser;
import ubic.gemma.model.association.Gene2GOAssociationService;
import ubic.gemma.model.genome.Gene;
import ubic.gemma.model.genome.Taxon;
import ubic.gemma.model.genome.TaxonService;
import ubic.gemma.model.genome.gene.GeneService;
import ubic.gemma.model.genome.PredictedGene;
import ubic.gemma.model.genome.ProbeAlignedRegion;
import ubic.gemma.ontology.GeneOntologyService;
import ubic.gemma.ontology.GoMetric;
import ubic.gemma.ontology.OntologyTerm;
import ubic.gemma.ontology.GoMetric.Metric;
import ubic.gemma.util.AbstractSpringAwareCLI;
import ubic.gemma.util.ConfigUtils;
import ubic.gemma.model.expression.arrayDesign.*;
import ubic.gemma.model.expression.designElement.*;

/**
 * @author meeta
 * @version $Id$
 */
public class ComputeGoOverlapCli extends AbstractSpringAwareCLI {

    private static final String HASH_MAP_RETURN = "HashMapReturn";

    /*
     * (non-Javadoc)
     * 
     * @see ubic.gemma.util.AbstractCLI#buildOptions()
     */

    private static final String GO_PROB_MAP = "GoProbMap";

    private static final String HOME_DIR = ConfigUtils.getString( "gemma.appdata.home" );

    private static final String RANDOM_SUBSET = "RandomSubset1K";

    private static final String VECTOR_MATRIX = "VectorMatrix";

    private static final String GENE_CACHE = "geneCache";

    private static Integer SET_SIZE = 0;

    public static void main( String[] args ) {
        ComputeGoOverlapCli p = new ComputeGoOverlapCli();
        try {
            Exception ex = p.doWork( args );
            if ( ex != null ) {
                ex.printStackTrace();
            }
        } catch ( Exception e ) {
            throw new RuntimeException( e );
        }
    }

    private int firstGeneColumn = 0;

    private int secondGeneColumn = 1;

    private String arrayDesign = "";// holds inputted string indicating array

    // design

    private Taxon taxon;

    // A list of service beans
    private GeneService geneService;

    private Gene2GOAssociationService gene2GOAssociationService;

    private GeneOntologyService ontologyEntryService;

    private ArrayDesignService arrayDesignService;

    private CompositeSequenceService css;

    private TaxonService taxonService;

    private GoMetric goMetric;

    private Map<Long, Collection<String>> geneGoMap = new HashMap<Long, Collection<String>>();

    private Map<String, Gene> geneCache = new HashMap<String, Gene>();

    private Map<String, Integer> rootMap = new HashMap<String, Integer>();

    Map<String, Integer> GOcountMap = new HashMap<String, Integer>();

    Map<String, Double> GOProbMap = new HashMap<String, Double>();

    DoubleMatrixNamed<Long, String> geneVectorMatrix = new SparseRaggedDoubleMatrix2DNamed<Long, String>();

    private String filePath = "";

    private Metric metric = GoMetric.Metric.simple;

    private boolean max = false;

    private boolean weight = false;

    private boolean random = false;

    private boolean arraySelected = false;// true when array design option is

    // inputted

    private String outFile;

    // INCLUDE PARTOF OR CHANGE STRINGENCY
    private boolean partOf = true;

    private String process = "http://purl.org/obo/owl/GO#GO_0008150";

    private String function = "http://purl.org/obo/owl/GO#GO_0003674";

    private String component = "http://purl.org/obo/owl/GO#GO_0005575";

    private int pCount = 0;

    private int fCount = 0;

    private int cCount = 0;

    @Override
    @SuppressWarnings("static-access")
    protected void buildOptions() {
        Option goMetricOption = OptionBuilder.hasArg().withArgName( "Choice of GO Metric" ).withDescription(
                "resnik, lin, jiang, percent, cosine, kappa; default = simple" ).withLongOpt( "metric" ).create( 'm' );
        addOption( goMetricOption );

        Option maxOption = OptionBuilder.hasArg().withArgName( "Choice of using MAX calculation" ).withDescription(
                "MAX" ).withLongOpt( "max" ).create( 'x' );
        addOption( maxOption );

        Option weightedOption = OptionBuilder.hasArg().withArgName( "Choice of using weighted matrix" )
                .withDescription( "weight" ).withLongOpt( "weight" ).create( 'w' );
        addOption( weightedOption );

        Option dataOption = OptionBuilder.hasArg().withArgName(
                "Choice of generating random gene pairs OR Input data file" ).withDescription( "dataType" )
                .isRequired().create( 'd' );
        addOption( dataOption );

        Option taxonOption = OptionBuilder.hasArg().withArgName( "Choice of taxon" ).withDescription(
                "human, rat, mouse" ).isRequired().create( 't' );
        addOption( taxonOption );

        Option firstGeneOption = OptionBuilder.hasArg().withArgName( "colindex" ).withDescription(
                "Column index with gene1 (starting from 0; default=0)" ).create( "g1col" );
        Option secondGeneOption = OptionBuilder.hasArg().withArgName( "colindex" ).withDescription(
                "Column index with gene2 (starting from 0; default=1)" ).create( "g2col" );
        addOption( firstGeneOption );
        addOption( secondGeneOption );

        Option arrayDesignOption = OptionBuilder.hasArg().withArgName( "Array Design" ).withDescription(
                "Complete Name of Microarray Design " ).create( "array" );
        addOption( arrayDesignOption );

        Option outFileOption = OptionBuilder.hasArg().isRequired().withArgName( "outFile" ).withDescription(
                "Write output to this file" ).create( 'o' );
        addOption( outFileOption );

    }

    /**
     * @param taxon
     */
    private void cacheGeneGoInformationToDisk( Taxon taxon ) {
        Collection<Gene> mouseGenes = geneService.loadKnownGenes( taxon );

        populateGeneGoMap( mouseGenes );
        saveCacheToDisk( ( HashMap ) geneGoMap, HASH_MAP_RETURN );
    }

    /**
     * 
     */
    private void computeTermProbabilities() {
        File f2 = new File( HOME_DIR + File.separatorChar + GO_PROB_MAP );
        if ( f2.exists() ) {
            GOProbMap = ( HashMap<String, Double> ) getCacheFromDisk( f2 );
            log.info( "Found probability file!" );
        }

        else {
            log.info( "Calculating probabilities... " );

            GOcountMap = goMetric.getTermOccurrence( geneGoMap );
            makeRootMap( GOcountMap.keySet() );
            GOProbMap = makeProbeMap( GOcountMap );

            this.saveCacheToDisk( ( HashMap ) GOProbMap, GO_PROB_MAP );
        }
    }

    private void createGene2TermMatrix() {

        File f3 = new File( HOME_DIR + File.separatorChar + VECTOR_MATRIX );
        if ( f3.exists() ) {
            geneVectorMatrix = ( DoubleMatrixNamed<Long, String> ) getCacheFromDisk( f3 );
            log.info( "Found vector matrix file!" );
        }

        else {
            log.info( "Creating sparse matrix... " );

            geneVectorMatrix = goMetric.createVectorMatrix( geneGoMap, weight );

            this.saveCacheToDisk( ( DoubleMatrixNamed ) geneVectorMatrix, VECTOR_MATRIX );
        }

    }

    /*
     * (non-Javadoc)
     * 
     * @see ubic.gemma.util.AbstractCLI#doWork(java.lang.String[])
     */
    @SuppressWarnings("unchecked")
    @Override
    protected Exception doWork( String[] args ) {

        Exception err = processCommandLine( "Compute Go Overlap ", args );
        if ( err != null ) return err;

        // populateGeneGoMapForTaxon();

        Collection<GenePair> genePairs = getLinks();

        populateGeneGoMapBasedOnGenesInLinks( genePairs );

        Map<GenePair, Double> scoreMap = scorePairs( genePairs );

        writeGoSimilarityResults( scoreMap );
        return null;

    }

    public Serializable getCacheFromDisk( File f ) {
        Serializable returnObject = null;
        try {
            if ( f.exists() ) {
                FileInputStream fis = new FileInputStream( f );
                ObjectInputStream ois = new ObjectInputStream( fis );
                returnObject = ( Serializable ) ois.readObject();
                ois.close();
                fis.close();
            }
        } catch ( Throwable e ) {
            return null;
        }
        return returnObject;
    }

    /**
     * @return
     */
    private Collection<GenePair> getLinks() {
        Collection<GenePair> genePairs = new HashSet<GenePair>();

        if ( random ) {
            genePairs = loadRandomPairs();
        } else if ( StringUtils.isNotBlank( filePath ) ) {
            try {
                File links = new File( filePath );
                genePairs = loadLinks( links, taxon );
                log.info( "Done loading data..." );
            } catch ( IOException e ) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }

        } else {
            log.error( "What should I do?" );
            bail( ErrorCode.INVALID_OPTION );
        }
        return genePairs;
    }

    /**
     * @param take a collection of genes and size of susbset
     * @return a collection of random gene pairs that have GO annotations
     */
    private Collection<GenePair> getRandomPairs( int size, Collection<Gene> genes ) {

        Collection<GenePair> subsetPairs = new HashSet<GenePair>();
        int i = 0;

        while ( i < size ) {
            List<Gene> twoGenes = new ArrayList( RandomChooser.chooseRandomSubset( 2, genes ) );

            if ( twoGenes.size() != 2 ) {
                log.warn( "A pair consists of two objects. More than two is unacceptable. Fix it!!" );
            }
            Collections.sort( twoGenes, new GeneComparator() );

            GenePair genePair = new GenePair();

            Gene g1 = twoGenes.get( 0 );
            Gene g2 = twoGenes.get( 1 );

            genePair.addFirstGene( g1 );
            genePair.addSecondGene( g2 );

            // Collections.sort( genePair );

            if ( subsetPairs.contains( genePair ) ) continue;

            Iterator<Gene> iterator = twoGenes.iterator();
            boolean noGoTerms = false;
            while ( iterator.hasNext() ) {
                if ( !geneGoMap.containsKey( iterator.next().getId() ) ) {
                    noGoTerms = true;
                    break;
                }
            }

            if ( noGoTerms ) continue;
            subsetPairs.add( genePair );
            log.info( "Added pair to subset!" );
            i++;
        }

        return subsetPairs;
    }

    /**
     * @param take a collection of probes and size of susbset
     * @return a collection of random gene pairs that may or may not have GO annotations
     */
    private Collection<GenePair> getRandomPairsFromProbes( int size, Collection<CompositeSequence> probes ) {

        Collection<GenePair> subsetPairs = new HashSet<GenePair>();
        int i = 0;
        GenePair genePair;

        while ( i < size ) {
            List<CompositeSequence> twoProbes = new ArrayList( RandomChooser.chooseRandomSubset( 2, probes ) );

            if ( twoProbes.size() != 2 ) {
                log.warn( "A pair consists of two objects. More than two is unacceptable. Fix it!!" );
            }
            Collections.sort( twoProbes, new ProbeComparator() );

            genePair = new GenePair();

            Collection<Gene> firstGenes = css.getGenes( twoProbes.get( 0 ) );
            Collection<Gene> secondGenes = css.getGenes( twoProbes.get( 1 ) );

            for ( Gene firstGene : firstGenes ) {
                if ( !( firstGene instanceof PredictedGene ) && !( firstGene instanceof ProbeAlignedRegion ) ) {// known
                    // gene
                    genePair.addFirstGene( firstGene );
                }
            }

            for ( Gene secondGene : secondGenes ) {
                if ( !( secondGene instanceof PredictedGene ) && !( secondGene instanceof ProbeAlignedRegion ) ) {// known
                    // gene
                    genePair.addSecondGene( secondGene );
                }
            }

            if ( genePair.getFirstGenes().isEmpty() || genePair.getSecondGenes().isEmpty() ) continue;

            if ( subsetPairs.contains( genePair ) ) continue;// genePair already added

            subsetPairs.add( genePair );
            log.info( "Added pair to subset!" );
            i++;
        }

        return subsetPairs;
    }

    private Collection<OntologyTerm> getTermOverlap( Gene g, Gene coexpG ) {

        Collection<String> masterGO = geneGoMap.get( g.getId() );
        Collection<String> coExpGO = geneGoMap.get( coexpG.getId() );
        Collection<OntologyTerm> overlapTerms = new HashSet<OntologyTerm>();

        if ( ( coExpGO == null ) || coExpGO.isEmpty() ) return null;

        if ( ( masterGO == null ) || masterGO.isEmpty() ) return null;

        for ( String ontologyEntry : masterGO ) {
            if ( ontologyEntry.equalsIgnoreCase( process ) || ontologyEntry.equalsIgnoreCase( function )
                    || ontologyEntry.equalsIgnoreCase( component ) ) continue;
            for ( String ontologyEntryC : coExpGO ) {

                if ( ontologyEntry.equalsIgnoreCase( process ) || ontologyEntry.equalsIgnoreCase( function )
                        || ontologyEntry.equalsIgnoreCase( component ) ) continue;

                if ( ontologyEntry.equalsIgnoreCase( ontologyEntryC ) )
                    overlapTerms.add( GeneOntologyService.getTermForURI( ontologyEntry ) );
            }
        }

        return overlapTerms;
    }

    protected void initBeans() {
        taxonService = ( TaxonService ) getBean( "taxonService" );
        geneService = ( GeneService ) getBean( "geneService" );
        gene2GOAssociationService = ( Gene2GOAssociationService ) getBean( "gene2GOAssociationService" );
        ontologyEntryService = ( GeneOntologyService ) getBean( "geneOntologyService" );
        goMetric = ( GoMetric ) getBean( "goMetric" );
        arrayDesignService = ( ArrayDesignService ) getBean( "arrayDesignService" );
        css = ( CompositeSequenceService ) getBean( "compositeSequenceService" );

    }

    /**
     * 
     */
    private void initGO() {
        /*
         * Initialize the Gene Ontology.
         */
        this.ontologyEntryService.init( true );

        while ( !ontologyEntryService.isReady() ) {
            log.info( "waiting for ontology load.." );
            try {
                Thread.sleep( 10000 );
                if ( !ontologyEntryService.isRunning() ) break;
            } catch ( InterruptedException e ) {
                throw new RuntimeException( e );
            }
        }

        if ( !ontologyEntryService.isReady() ) {
            throw new RuntimeException( "Gene Ontology was not loaded successfully." );
        }
    }

    /**
     * Opens a file for writing anda adds the header.
     * 
     * @param fileName if Null, output will be written to standard output.
     * @throws IOException
     */
    protected Writer initOutputFile( String fileName ) throws IOException {

        Writer writer;
        if ( StringUtils.isBlank( fileName ) ) {
            log.info( "Output to stdout" );
            writer = new PrintWriter( System.out );
        } else {

            // write into file
            log.info( "Creating new annotation file " + fileName + " \n" );

            File f = new File( fileName );

            if ( f.exists() ) {
                log.warn( "Will overwrite existing file " + f );
                f.delete();
            }

            f.createNewFile();
            writer = new FileWriter( f );
        }

        writer.write( "Gene1\tGene2\tScore\tG1GOTerms\tG2GOTerms\tTermOverlap\tGOTerms\n" );

        return writer;
    }

    private Set<GenePair> loadLinks( File f, Taxon taxon ) throws IOException {

        log.info( "Loading data from " + f );
        BufferedReader in = new BufferedReader( new FileReader( f ) );

        Set<GenePair> geneMap = new HashSet<GenePair>();

        String line;
        boolean alreadyWarned = false;
        while ( ( line = in.readLine() ) != null ) {
            line = line.trim();
            if ( line.startsWith( "#" ) ) {
                continue;
            }

            String[] fields = StringUtils.split( line );

            if ( fields.length < 2 ) {
                if ( !alreadyWarned ) {
                    log.warn( "Bad field on line: " + line + " (subsequent errors suppressed)" );
                    alreadyWarned = true;
                }
                continue;
            }

            String g1 = fields[firstGeneColumn];
            String g2 = fields[secondGeneColumn];
            // skip any self links.
            // if ( g1.equals( g2 ) ) continue;

            String[] gene1Strings = StringUtils.split( g1, "," );
            String[] gene2Strings = StringUtils.split( g2, "," );

            GenePair genePair = new GenePair();
            Collection<Gene> genes;
            for ( String gene1string : gene1Strings ) {

                Gene gene1 = null;
                if ( geneCache.containsKey( gene1string ) ) {
                    gene1 = geneCache.get( gene1string );
                } else {
                    genes = geneService.findByOfficialSymbol( gene1string );
                    for ( Gene gene : genes ) {
                        if ( gene.getTaxon().equals( taxon ) && ontologyEntryService.getGOTerms( gene ).size() > 0 ) {
                            geneCache.put( gene1string, gene );
                            gene1 = gene;
                            break;
                        }
                    }
                }

                if ( gene1 == null ) continue;

                genePair.addFirstGene( gene1 );

                for ( String gene2string : gene2Strings ) {
                    if ( gene1string.equals( gene2string ) ) {
                        continue;
                    }

                    Gene gene2 = null;

                    if ( geneCache.containsKey( gene2string ) ) {
                        gene2 = geneCache.get( gene2string );
                    } else {
                        genes = geneService.findByOfficialSymbol( gene2string );
                        for ( Gene gene : genes ) {
                            if ( gene.getTaxon().equals( taxon ) && ontologyEntryService.getGOTerms( gene ).size() > 0 ) {
                                geneCache.put( gene2string, gene );
                                gene2 = gene;
                                break;
                            }
                        }
                    }

                    if ( gene2 == null ) continue;

                    genePair.addSecondGene( gene2 );

                    // Collections.sort( genePair, new GeneComparator() );
                    geneMap.add( genePair );

                    if ( geneMap.size() % 50000 == 0 ) {
                        log.info( "Loaded " + geneMap.size() + " links" );
                    }

                }
            }

            // compute the median of gooverlaps and do something with the
            // result.

        }
        log.info( "Loaded " + geneMap.size() + " links" );
        saveCacheToDisk( ( HashMap ) geneCache, GENE_CACHE );
        return geneMap;
    }

    /**
     * @return
     */
    private Collection<GenePair> loadRandomPairs() {
        Collection<GenePair> randomPairs = null;
        Gene tempGene = null;
        if ( arraySelected ) {
            ArrayDesign ad = arrayDesignService.findByShortName( arrayDesign );
            arrayDesignService.thawLite( ad );
            if ( ad == null ) {
                System.out.println( "Array design " + arrayDesign + " not found" );
                System.exit( 0 );
            }
            Collection<CompositeSequence> probes = ad.getCompositeSequences();
            Collection<Gene> adGenes = new ArrayList<Gene>();
            css.thaw( probes );
            randomPairs = getRandomPairsFromProbes( SET_SIZE, probes );
        } else {
            try {
                File f3 = new File( HOME_DIR + File.separatorChar + RANDOM_SUBSET );
                if ( f3.exists() ) {
                    randomPairs = ( HashSet<GenePair> ) loadLinks( f3, this.taxon );
                    log.info( "Found cached subset file!" );
                }

                else {
                    Collection<Gene> allKnownGenes = geneService.loadKnownGenes( taxon );
                    Collection<Gene> allGOGenes = new HashSet<Gene>();
                    populateGeneGoMap( allKnownGenes );
                    //saveCacheToDisk( ( Serializable ) geneGoMap, HASH_MAP_RETURN );

                    for ( Long g : geneGoMap.keySet() ) {
                        Gene gene = geneService.load( g );
                        allGOGenes.add( gene );
                    }

                    randomPairs = getRandomPairs( SET_SIZE, allGOGenes );

                    Writer w = new FileWriter( f3 );
                    this.writeLinks( randomPairs, w );
                }
            } catch ( IOException e ) {
                throw new RuntimeException( e );
            }
        }
        return randomPairs;
    }

    private Map<String, Double> makeProbeMap( Map<String, Integer> GOcountMap ) {

        for ( String uri : GOcountMap.keySet() ) {
            int total = 0;
            log.info( "Counting children for " + uri );
            int count = goMetric.getChildrenOccurrence( GOcountMap, uri );
            if ( rootMap.get( uri ) == 1 ) total = pCount;
            if ( rootMap.get( uri ) == 2 ) total = fCount;
            if ( rootMap.get( uri ) == 3 ) total = cCount;

            GOProbMap.put( uri, ( double ) count / total );
        }

        return GOProbMap;
    }

    /**
     * @param take a collection of GOTerm URIs
     * @return Identify the root of each term and put it in the rootMap
     */
    private void makeRootMap( Collection<String> terms ) {

        Collection<String> remove = new HashSet<String>();

        for ( String t : terms ) {
            Collection<OntologyTerm> parents = ontologyEntryService.getAllParents( GeneOntologyService
                    .getTermForURI( t ), partOf );

            for ( OntologyTerm p : parents ) {
                if ( p.getUri().equalsIgnoreCase( process ) ) {
                    rootMap.put( t, 1 );
                    pCount += GOcountMap.get( t );
                    break;
                }
                if ( p.getUri().equalsIgnoreCase( function ) ) {
                    rootMap.put( t, 2 );
                    fCount += GOcountMap.get( t );
                    break;
                }
                if ( p.getUri().equalsIgnoreCase( component ) ) {
                    rootMap.put( t, 3 );
                    cCount += GOcountMap.get( t );
                    break;
                }
            }
            if ( !( rootMap.containsKey( t ) ) ) {
                log.warn( "Couldn't get root for term: " + t );
                remove.add( t );
            }
        }

        for ( String s : remove ) {
            GOcountMap.remove( s );
        }
    }

    /**
     * @param genes
     */
    private void populateGeneGoMap( Collection<Gene> genes ) {
        for ( Gene gene : genes ) {
            Collection<OntologyTerm> GOTerms = ontologyEntryService.getGOTerms( gene, partOf );

            if ( GOTerms == null || GOTerms.isEmpty() ) continue;
            log.info( "Got go terms for " + gene.getName() );

            Collection<String> termString = new HashSet<String>();
            for ( OntologyTerm oe : GOTerms ) {
                termString.add( oe.getUri() );
            }
            geneGoMap.put( gene.getId(), termString );
        }
    }

    /**
     * @param genePairs
     */
    private void populateGeneGoMapBasedOnGenesInLinks( Collection<GenePair> genePairs ) {
        Collection<Gene> allGenes = new HashSet<Gene>();

        if ( metric.equals( GoMetric.Metric.simple ) && geneGoMap.isEmpty() ) {
            for ( GenePair pair : genePairs ) {
                allGenes.addAll( pair.getGenes() );
            }
            populateGeneGoMap( allGenes );
        } else {
            populateGeneGoMapForTaxon();
        }
    }

    /**
     * 
     */
    private void populateGeneGoMapForTaxon() {

        if ( geneGoMap.isEmpty() ) {
            log.info( "Checking for Gene2GO Map file..." );
            File f = new File( HOME_DIR + File.separatorChar + HASH_MAP_RETURN );
            if ( f.exists() ) {
                geneGoMap = ( Map<Long, Collection<String>> ) getCacheFromDisk( f );
                log.info( "Found file!" );
            } else {
                cacheGeneGoInformationToDisk( taxon );
            }
        }

        if ( metric.equals( GoMetric.Metric.resnik ) || metric.equals( GoMetric.Metric.lin )
                || metric.equals( GoMetric.Metric.jiang ) ) {
            log.info( "Computing term probabilities for all GOterms ..." );
            computeTermProbabilities();
        }

        if ( metric.equals( GoMetric.Metric.cosine ) || metric.equals( GoMetric.Metric.kappa ) ) {
            createGene2TermMatrix();
        }

    }

    @Override
    protected void processOptions() {
        super.processOptions();
        initBeans();
        String commonName = getOptionValue( 't' );
        if ( StringUtils.isBlank( commonName ) ) {
            System.out.println( "MUST enter a valid taxon!" );
            System.exit( 0 );
        }
        if ( !StringUtils.equalsIgnoreCase( "mouse", commonName ) && !StringUtils.equalsIgnoreCase( "rat", commonName )
                && !StringUtils.equalsIgnoreCase( "human", commonName ) ) {
            System.out.println( "MUST enter a valid taxon!" );
            System.exit( 0 );
        }

        this.taxon = taxonService.findByCommonName( commonName );
        if ( taxon == null ) {
            System.out.println( "Taxon " + commonName + " not found." );
            System.exit( 0 );
        }

        if ( hasOption( 'd' ) ) {
            String input = getOptionValue( 'd' );
            if ( this.hasOption( "array" ) ) {
                this.arrayDesign = getOptionValue( "array" );
                if ( StringUtils.isNumeric( input ) ) {
                    this.random = true;
                    this.arraySelected = true;
                    SET_SIZE = Integer.parseInt( input );

                } else {
                    System.out
                            .println( input
                                    + " is not a numeric value, dataset size must be a numeric value when -array option is present" );
                    System.exit( 0 );
                }
            } else {
                if ( StringUtils.isNumeric( input ) ) {
                    SET_SIZE = Integer.parseInt( input );
                    this.random = true;
                    System.out.println( "Will create a set of " + SET_SIZE + " random gene pairs!" );
                } else {
                    File f = new File( input );
                    if ( f.canRead() ) {
                        this.filePath = input;
                        this.random = false;
                    } else {
                        System.out.println( input
                                + "is NOT a valid filename! You MUST enter a valid filename OR size of dataset" );
                        System.exit( 0 );
                    }
                }
            }
        }

        if ( hasOption( 'x' ) ) {
            String max = getOptionValue( 'x' );
            if ( max.equalsIgnoreCase( "MAX" ) )
                this.max = true;
            else
                this.max = false;
        }

        if ( hasOption( 'w' ) ) {
            String weight = getOptionValue( 'w' );
            if ( weight.equalsIgnoreCase( "weight" ) )
                this.weight = true;
            else
                this.weight = false;
        }

        if ( hasOption( 'm' ) ) {
            String metricName = getOptionValue( 'm' );
            if ( metricName.equalsIgnoreCase( "resnik" ) )
                this.metric = GoMetric.Metric.resnik;
            else if ( metricName.equalsIgnoreCase( "lin" ) )
                this.metric = GoMetric.Metric.lin;
            else if ( metricName.equalsIgnoreCase( "jiang" ) )
                this.metric = GoMetric.Metric.jiang;
            else if ( metricName.equalsIgnoreCase( "percent" ) )
                this.metric = GoMetric.Metric.percent;
            else if ( metricName.equalsIgnoreCase( "cosine" ) )
                this.metric = GoMetric.Metric.cosine;
            else if ( metricName.equalsIgnoreCase( "kappa" ) )
                this.metric = GoMetric.Metric.kappa;
            else {
                this.metric = GoMetric.Metric.simple;
                this.max = false;
            }
        }

        if ( this.hasOption( "g1col" ) ) {
            this.firstGeneColumn = getIntegerOptionValue( "g1col" );
        }

        if ( this.hasOption( "g2col" ) ) {
            this.secondGeneColumn = getIntegerOptionValue( "g2col" );
        }

        outFile = getOptionValue( 'o' );

        initGO();
    }

    public void saveCacheToDisk( Serializable toSave, String filename ) {

        log.info( "Generating file ... " );

        try {
            // remove file first
            File f = new File( HOME_DIR + File.separatorChar + filename );
            if ( f.exists() ) {
                f.delete();
            }
            FileOutputStream fos = new FileOutputStream( f );
            ObjectOutputStream oos = new ObjectOutputStream( fos );
            oos.writeObject( toSave );
            oos.flush();
            oos.close();
        } catch ( Throwable e ) {
            log.error( "Cannot write to file." );
            return;
        }
        log.info( "Done making report." );
    }

    /**
     * @param genePairs
     * @return
     */
    private Map<GenePair, Double> scorePairs( Collection<GenePair> genePairs ) {
        Map<GenePair, Double> scoreMap = new HashMap<GenePair, Double>();
        log.info( genePairs.size() + " pairs to score" );
        for ( GenePair pair : genePairs ) {

            List<Gene> genes1 = pair.getFirstGenes();
            List<Gene> genes2 = pair.getSecondGenes();

            List<Double> scores = new ArrayList<Double>();

            for ( Gene g1 : genes1 ) {
                for ( Gene g2 : genes2 ) {
                    if ( g1.getId() == g2.getId() ) continue;
                    double score = 0.0;

                    if ( metric.equals( GoMetric.Metric.cosine ) || metric.equals( GoMetric.Metric.kappa ) ) {
                        score = goMetric.computeMatrixSimilarity( g1, g2, geneVectorMatrix, metric );
                    } else if ( max ) {
                        log.info( "getting MAX scores for " + metric );
                        score = goMetric.computeMaxSimilarity( g1, g2, GOProbMap, metric );
                    } else {
                        score = goMetric.computeSimilarity( g1, g2, GOProbMap, metric );
                    }
                    scores.add( score );
                }
            }
            DoubleArrayList dlist = new DoubleArrayList();
            dlist.addAllOf( scores );

            double median = Descriptive.median( dlist );
            log.info( "Adding pair " + pair + " with score " + median );
            scoreMap.put( pair, median );
        }
        return scoreMap;
    }

    /**
     * @param scoreMap
     * @param overallWatch
     */
    private void writeGoSimilarityResults( Map<GenePair, Double> scoreMap ) {
        StopWatch overallWatch = new StopWatch();
        overallWatch.start();
        try {
            Writer write = initOutputFile( outFile );

            for ( GenePair pair : scoreMap.keySet() ) {
                List<Gene> firstGenes = pair.getFirstGenes();
                List<Gene> secondGenes = pair.getSecondGenes();

                String pairString = pair.toString();
                Double score = scoreMap.get( pair );
                if ( score == null ) {
                    log.warn( "Score is null for " + pair );
                    continue;
                }
                if ( firstGenes.size() == 1 && secondGenes.size() == 1 ) {
                    Gene mGene = firstGenes.iterator().next();
                    Gene cGene = secondGenes.iterator().next();
                    int masterGOTerms;
                    int coExpGOTerms;

                    if ( !geneGoMap.containsKey( mGene.getId() ) )
                        masterGOTerms = 0;
                    else
                        masterGOTerms = ( geneGoMap.get( mGene.getId() ) ).size();

                    if ( !geneGoMap.containsKey( cGene.getId() ) )
                        coExpGOTerms = 0;
                    else
                        coExpGOTerms = ( geneGoMap.get( cGene.getId() ) ).size();

                    Collection<OntologyTerm> goTerms = getTermOverlap( mGene, cGene );
                    writeOverlapLine( write, pairString, score, goTerms, masterGOTerms, coExpGOTerms );
                } else {
                    writeOverlapLine( write, pairString, score, null, 0, 0 );
                }

            }
            overallWatch.stop();
            log.info( "Compute GoOverlap takes " + overallWatch.getTime() + "ms" );

            // printResults(masterTermCountMap);
        } catch ( IOException ioe ) {
            log.error( "Couldn't write to file: " + ioe );

        }
    }

    private void writeLinks( Collection<GenePair> pairs, Writer writer ) {
        /*
         * print the pairs out in the same format that we can read in.
         */
    }

    /**
     * @param writer
     * @param pairString
     * @param overlap
     * @param goTerms
     * @param masterGOTerms
     * @param coExpGOTerms
     * @throws IOException
     */
    protected void writeOverlapLine( Writer writer, String pairString, double score, Collection<OntologyTerm> goTerms,
            int masterGOTerms, int coExpGOTerms ) throws IOException {

        if ( log.isDebugEnabled() ) log.debug( "Generating line for annotation file \n" );

        writer.write( pairString + "\t" + score + "\t" + masterGOTerms + "\t" + coExpGOTerms + "\t" );

        if ( goTerms == null || goTerms.isEmpty() ) {
            writer.write( "\n" );
            writer.flush();
            return;
        }

        boolean wrote = false;

        for ( OntologyTerm oe : goTerms ) {
            if ( oe == null ) continue;
            if ( wrote )
                writer.write( "|" + GeneOntologyService.asRegularGoId( oe ) );
            else
                writer.write( GeneOntologyService.asRegularGoId( oe ) );
            wrote = true;
        }

        writer.write( "\n" );
        writer.flush();
    }

    class GeneComparator implements Comparator {

        public int compare( Object o1, Object o2 ) {

            if ( o1 instanceof Gene && o2 instanceof Gene ) {
                Long g1 = ( ( Gene ) o1 ).getId();
                Long g2 = ( ( Gene ) o2 ).getId();

                if ( g1 > g2 )
                    return 1;
                else if ( g1 < g2 ) return -1;
            }
            return 0;
        }

    }

    class ProbeComparator implements Comparator {

        public int compare( Object o1, Object o2 ) {

            if ( o1 instanceof CompositeSequence && o2 instanceof CompositeSequence ) {
                Long g1 = ( ( CompositeSequence ) o1 ).getId();
                Long g2 = ( ( CompositeSequence ) o2 ).getId();

                if ( g1 > g2 )
                    return 1;
                else if ( g1 < g2 ) return -1;
            }
            return 0;
        }

    }

    class GenePair extends ArrayList<List<Gene>> implements Comparable<GenePair> {

        public GenePair() {
            super( 0 );
            this.add( new ArrayList<Gene>() );
            this.add( new ArrayList<Gene>() );
        }

        public void addFirstGene( Gene g ) {
            List<Gene> firstGenes = this.get( 0 );
            for ( Gene firstGene : firstGenes ) {
                if ( firstGene.getId() == g.getId() ) {
                    return;
                }
            }
            this.get( 0 ).add( g );
        }

        public void addSecondGene( Gene g ) {
            List<Gene> secondGenes = this.get( 0 );
            for ( Gene secondGene : secondGenes ) {
                if ( secondGene.getId() == g.getId() ) {
                    return;
                }
            }
            this.get( 1 ).add( g );
        }

        /*
         * (non-Javadoc)
         * 
         * @see java.lang.Comparable#compareTo(java.lang.Object)
         */
        public int compareTo( GenePair o2 ) {
            return this.getFirstGenes().iterator().next().getName().compareTo(
                    o2.getFirstGenes().iterator().next().getName() );
        }

        /*
         * (non-Javadoc)
         * 
         * @see java.lang.Object#equals(java.lang.Object)
         */
        @Override
        public boolean equals( Object obj ) {
            if ( this == obj ) return true;
            if ( obj == null ) return false;
            if ( getClass() != obj.getClass() ) return false;
            final GenePair other = ( GenePair ) obj;
            for ( Gene g : other.getFirstGenes() ) {
                for ( Gene fg : this.getFirstGenes() ) {
                    if ( fg.getId() == g.getId() ) {
                        return true;
                    }
                }
            }
            for ( Gene g : other.getSecondGenes() ) {
                for ( Gene sg : this.getSecondGenes() ) {
                    if ( sg.getId() == g.getId() ) {
                        return true;
                    }
                }
            }
            return false;
            /*
             * for ( Gene g : other.getFirstGenes() ) { if ( !this.getFirstGenes().contains( g ) ) return false; } for (
             * Gene g : other.getSecondGenes() ) { if ( !this.getSecondGenes().contains( g ) ) return false; } return
             * true;
             */
        }

        /**
         * @return the firstGenes
         */
        public List<Gene> getFirstGenes() {
            return this.get( 0 );
        }

        public Set<Gene> getGenes() {
            Set<Gene> result = new HashSet<Gene>();
            result.addAll( getFirstGenes() );
            result.addAll( getSecondGenes() );
            return result;
        }

        /**
         * @return the secondGenes
         */
        public List<Gene> getSecondGenes() {
            return this.get( 1 );
        }

        /*
         * (non-Javadoc)
         * 
         * @see java.lang.Object#hashCode()
         */
        @Override
        public int hashCode() {
            final int PRIME = 31;
            int result = 1;
            for ( Gene g : getFirstGenes() ) {
                result = PRIME * result + ( ( g == null ) ? 0 : g.hashCode() );
            }
            for ( Gene g : getSecondGenes() ) {
                result = PRIME * result + ( ( g == null ) ? 0 : g.hashCode() );
            }
            return result;
        }

        /**
         * @param firstGenes the firstGenes to set
         */
        public void setFirstGenes( List<Gene> firstGenes ) {
            this.set( 0, firstGenes );
        }

        /**
         * @param secondGenes the secondGenes to set
         */
        public void setSecondGenes( List<Gene> secondGenes ) {
            this.set( 1, secondGenes );
        }

        @Override
        public String toString() {
            StringBuilder buf = new StringBuilder();
            Collections.sort( this.get( 0 ), new GeneComparator() );
            Collections.sort( this.get( 1 ), new GeneComparator() );
            for ( Iterator<Gene> it = this.get( 0 ).iterator(); it.hasNext(); ) {
                buf.append( it.next().getName() );
                if ( it.hasNext() ) buf.append( "," );
            }
            buf.append( "\t" );
            for ( Iterator<Gene> it = this.get( 1 ).iterator(); it.hasNext(); ) {
                buf.append( it.next().getName() );
                if ( it.hasNext() ) buf.append( "," );
            }
            return buf.toString();
        }

    }

}
