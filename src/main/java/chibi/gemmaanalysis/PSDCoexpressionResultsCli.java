package chibi.gemmaanalysis;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.gemma.analysis.expression.coexpression.CoexpressionMetaValueObject;
import ubic.gemma.analysis.expression.coexpression.CoexpressionValueObjectExt;
import ubic.gemma.analysis.expression.coexpression.GeneCoexpressionService;
import ubic.gemma.model.analysis.expression.coexpression.GeneCoexpressionAnalysisService;
import ubic.gemma.model.association.coexpression.Gene2GeneCoexpressionService;
import ubic.gemma.model.genome.Gene;
import ubic.gemma.model.genome.Taxon;
import ubic.gemma.model.genome.TaxonService;
import ubic.gemma.model.genome.gene.GeneService;
import ubic.gemma.util.AbstractSpringAwareCLI;

public class PSDCoexpressionResultsCli extends AbstractSpringAwareCLI {

    private static final int DEFAULT_STRINGINCY = 2;

    private int stringency;

    // private Collection<Gene> geneList;
    private String fileName;

    private boolean allGenes;

    private Collection<String> rawGeneList;

    private Taxon taxon;

    private boolean details;

    private boolean queryGenesOnly;

    private static Log log = LogFactory.getLog( PSDCoexpressionResultsCli.class );

    protected GeneService geneService;

    protected Gene2GeneCoexpressionService gene2GeneCoexpressionService;

    protected TaxonService taxonService;

    protected Long analysisID;

    protected GeneCoexpressionService geneCoexpressionService;

    protected GeneCoexpressionAnalysisService geneCoexpressionAnalysisService;

    public static void main( String args[] ) {
        PSDCoexpressionResultsCli run = new PSDCoexpressionResultsCli();
        try {
            Exception ex = run.doWork( args );
            if ( ex != null ) {
                ex.printStackTrace();
            }
        } catch ( Exception e ) {
            throw new RuntimeException( e );
        }
    }

    @SuppressWarnings("static-access")
    @Override
    protected void buildOptions() {
        Option geneFileOption = OptionBuilder
                .hasArg()
                .isRequired()
                .withArgName( "Gene List File Name" )
                .withDescription( "A text file that contains a list of gene symbols, with one gene symbol on each line" )
                .withLongOpt( "geneFile" ).create( 'g' );
        addOption( geneFileOption );

        Option stringencyOption = OptionBuilder.hasArg().withArgName( "Stringency" ).withDescription(
                "The stringency value: Defaults to " + DEFAULT_STRINGINCY ).withLongOpt( "stringency" ).create( 's' );
        addOption( stringencyOption );

        Option taxonOption = OptionBuilder.hasArg().isRequired().withArgName( "Taxon" ).withDescription(
                "The name of the taxon." ).withLongOpt( "taxon" ).create( 't' );
        addOption( taxonOption );

        Option analysisOption = OptionBuilder.hasArg().withArgName( "CoexpressionAnalysis" ).withDescription(
                "Gene coexpression analysis ID to use." ).withLongOpt( "analysis" ).create( 'a' );
        addOption( analysisOption );

        Option queryOption = OptionBuilder.hasArg().withArgName( "QueryGenesOnly" ).withDescription(
                "Link analysis on query genes only?" ).withLongOpt( "query" ).create( 'q' );
        addOption( queryOption );

        Option randomOption = OptionBuilder.hasArg().withArgName( "RandomFileOrCmd" ).withDescription(
                "all = all known genes in taxon will be used" ).withLongOpt( "random" ).create( 'r' );
        addOption( randomOption );

        Option detailOption = OptionBuilder.hasArg().withArgName( "ExtraDetails" ).withDescription(
                "true = extra details, false = just coexpressed genes" ).withLongOpt( "detail" ).create( 'd' );
        addOption( detailOption );
    }

    @SuppressWarnings("unchecked")
    @Override
    protected Exception doWork( String[] args ) {
        Exception err = processCommandLine( "Gene 2 Gene Coexpression Results ", args );
        if ( err != null ) return err;

        try {
            Collection<Gene> genes = null;
            if ( allGenes == true )
                genes = geneService.loadKnownGenes( taxon );
            else
                genes = getGenes( rawGeneList, taxon );

            // Batch up the genes to reduce memory footprint

            int count = 0;
            int CHUNK_SIZE = 100;

            Collection<Gene> batch = new HashSet<Gene>();

            for ( Gene g : genes ) {
                batch.add( g );
                count++;
                if ( count % CHUNK_SIZE == 0 ) {
                    outputCoexpressionResults( batch, taxon, stringency, queryGenesOnly );
                    batch.clear();
                }
            }

            if ( batch.size() > 0 ) {
                outputCoexpressionResults( batch, taxon, stringency, queryGenesOnly );
            }

        } catch ( IOException e ) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

        return null;
    }

    @SuppressWarnings("unchecked")
    @Override
    protected void processOptions() {
        super.processOptions();

        geneService = ( GeneService ) getBean( "geneService" );
        taxonService = ( TaxonService ) getBean( "taxonService" );
        gene2GeneCoexpressionService = ( Gene2GeneCoexpressionService ) getBean( "gene2GeneCoexpressionService" );
        geneCoexpressionService = ( GeneCoexpressionService ) getBean( "geneCoexpressionService" );
        geneCoexpressionAnalysisService = ( GeneCoexpressionAnalysisService ) getBean( "geneCoexpressionAnalysisService" );

        // taxon option - t
        String taxonName = "";
        if ( hasOption( 't' ) ) {
            taxonName = this.getOptionValue( 't' );
            taxon = taxonService.findByCommonName( taxonName );
            if ( taxon == null ) {
                log.error( "ERROR: Cannot find taxon " + taxonName );
            }
        }
        // check to make sure taxon is also given if gene file is given
        if ( this.hasOption( 'g' ) ) {
            if ( !this.hasOption( 't' ) ) {
                log.info( "You must provide the taxon if you provide a gene file" );
                bail( ErrorCode.MISSING_ARGUMENT );
                return;
            }
        }

        // random option - g
        // if "all" is given, then use all genes in a given taxon
        if ( this.hasOption( 'g' ) ) {
            String rOption = this.getOptionValue( 'g' );
            assert taxon != null;
            allGenes = false; // default
            if ( rOption.equals( "all" ) ) // geneList = geneService.loadKnownGenes( taxon );
                allGenes = true;

            // no random option
            else {
                assert taxon != null;
                try {
                    // geneList = getGeneList( this.getOptionValue( 'g' ), taxon );
                    fileName = this.getOptionValue( 'g' );
                    readGeneListFile( fileName );
                } catch ( IOException e ) {
                    throw new RuntimeException( e );
                }
            }
        }

        // stringency option - s
        // default stringency = 2
        stringency = DEFAULT_STRINGINCY;
        if ( this.hasOption( 's' ) ) {
            stringency = Integer.parseInt( this.getOptionValue( 's' ) );
        }

        // link analysis on genelist only?
        queryGenesOnly = false; // default
        if ( this.hasOption( 'q' ) ) {
            String query = this.getOptionValue( 'q' );
            if ( query.equals( "true" ) || query.equals( "1" ) ) queryGenesOnly = true;
        }

        // extra details options - d
        details = false; // default
        if ( this.hasOption( 'd' ) ) {
            String d = this.getOptionValue( 'd' );
            if ( d.equals( "true" ) || d.equals( "1" ) ) details = true;
        }

        // analysis options - a
        if ( this.hasOption( 'a' ) ) {
            String analysis = this.getOptionValue( 'a' );
            analysisID = Long.parseLong( analysis );
        } else { // default cases
            if ( taxonName.equals( "mouse" ) )
                analysisID = 708L; // default, if mouse, is all mouse analysis
            else
                analysisID = 717L; // default is all human analysis
        }

    }

    /**
     * Read in a list of genes; calls helper method, <code>readGeneListFile()</code>
     * 
     * @param inFile - file name to read
     * @param taxon
     * @return collection of genes
     * @throws IOException
     */
    private void readGeneListFile( String inFile ) throws IOException {
        log.info( "Reading " + inFile );
        Collection<String> lines = new ArrayList<String>();
        BufferedReader in = new BufferedReader( new FileReader( inFile ) );
        String line;
        while ( ( line = in.readLine() ) != null ) {
            // if ( line.startsWith( "#" ) ) continue;
            String s = line.trim();
            lines.add( s );
        }
        in.close();
        // to be used later
        rawGeneList = lines;

    }

    /**
     * Helper method to read in the lines in a text file and store in a collection
     * 
     * @param inFile
     * @return - ArrayList storing the lines in a text file
     * @throws IOException
     */
    private Collection<Gene> getGenes( Collection<String> rawGeneList, Taxon taxon ) throws IOException {
        // Collection<String> rawLinesInFile = readGeneListFile( inFile );
        Collection<Gene> genes = new ArrayList<Gene>();

        Iterator linesIt = rawGeneList.iterator();
        while ( linesIt.hasNext() ) {
            String s = ( String ) linesIt.next();
            Gene gene = geneService.findByOfficialSymbol( s, taxon );
            if ( gene == null ) {
                log.error( "ERROR: Cannot find genes for " + s );
                continue;
            }
            geneService.thaw( gene );
            genes.add( gene );

        }
        // geneService.thawLite( genes );
        return genes;
    }

    // /**
    // * Uses precomputed file of IDs; calls helper method, <code>readGeneListFile()</code>
    // *
    // * @param inFile - Precomputed file of IDs
    // * @param i - number of random genes to use
    // * @return - list of 'i' random number of genes from precomputed list of IDs
    // * @throws IOException
    // */
    // private Collection<Gene> getRandomKnownGenes( String inFile, int i ) throws IOException {
    // Object[] rawLinesInFile = readGeneListFile( inFile ).toArray();
    // Collection<Gene> genes = new HashSet<Gene>();
    // String line;
    // int geneCount = 0;
    // SecureRandom randomSeed;
    // try {
    // randomSeed = SecureRandom.getInstance( "SHA1PRNG" );
    //
    // while ( geneCount < 1000 ) {
    //
    // randomSeed.generateSeed( 10 );
    // int randomInt = randomSeed.nextInt( rawLinesInFile.length );
    //
    // line = ( String ) rawLinesInFile[randomInt];
    //
    // Gene gene = geneService.load( Long.parseLong( line ) );
    // if ( gene == null ) {
    // log.error( "ERROR: Cannot find genes for ID: " + line );
    // continue;
    // }
    // geneService.thaw( gene );
    //
    // // because we are using a HashSet, no duplicates are allowed - which is what we want
    // genes.add( gene );
    // geneCount++;
    // }
    //
    // } catch ( NoSuchAlgorithmException e ) {
    // // TODO Auto-generated catch block
    // e.printStackTrace();
    // }
    //
    // return genes;
    // }

    @SuppressWarnings("unchecked")
    private void outputCoexpressionResults( Collection<Gene> geneList, Taxon taxon, int stringency,
            boolean queryGenesOnly ) throws IOException {
        log.debug( "Running Gene Coexpression..." );

        CoexpressionMetaValueObject cmvo = geneCoexpressionService.coexpressionSearch( analysisID, geneList,
                stringency, 0, queryGenesOnly );
        Collection<CoexpressionValueObjectExt> cvoExtCol = cmvo.getKnownGeneResults();
        log.debug( "Printing results..." );
        printCoexpressionResults( cvoExtCol );
    }

    private void printCoexpressionResults( Collection<CoexpressionValueObjectExt> cvoExtCol ) {
        createColumnHeadings();

        for ( CoexpressionValueObjectExt cvo : cvoExtCol )
            System.out.println( coexpressedResult( cvo ) );

    }

    // @SuppressWarnings("unchecked")
    // private Long getCannedAnalysisID( Taxon taxon2 ) {
    // Collection<GeneCoexpressionAnalysis> analysisCol = geneCoexpressionAnalysisService.findByTaxon( taxon );
    // GeneCoexpressionAnalysis analysis2Use = null;
    // // use the 1st canned analysis that isn't virtual for the given taxon (should be the all"Taxon" analysis)
    // for ( GeneCoexpressionAnalysis analysis : analysisCol ) {
    // if ( analysis instanceof GeneCoexpressionVirtualAnalysis )
    // continue;
    // else {
    // analysis2Use = analysis;
    // break;
    // }
    // }
    // // use the brain analysis
    // // for (GeneCoexpressionAnalysis analysis : analysisCol){
    // // if (analysis instanceof GeneCoexpressionVirtualAnalysis && analysis.getId() == 737){
    // // analysis2Use = analysis;
    // // break;
    // // }
    // // else
    // // continue;
    // // }
    //
    // return analysis2Use.getId();
    // }

    private void createColumnHeadings() {
        // output column headings for main coexpression results file
        if ( details == true )
            System.out.println( "Query_Gene\tCoexpressed_Gene\tSupport\tSign" );
        else
            System.out.println( "Query_Gene\tCoexpressed_Gene" );
    }

    /**
     * map of queryGene to nicely formated result - coexpressed gene, # of datasets, etc
     * 
     * @param cvoExtCol
     * @return mapping of the query gene to all expressed
     */
    // private Map<String, Collection<String>> organizeCoexpressionValueObjectResults(
    // Collection<CoexpressionValueObjectExt> cvoExtCol ) {
    // Map<String, Collection<String>> results = new HashMap<String, Collection<String>>();
    // for ( CoexpressionValueObjectExt cvo : cvoExtCol ) {
    // // see if query gene official symbol is already a key in the map
    // String queryGene = cvo.getQueryGene().getOfficialSymbol();
    // // Long foundGene = cvo.getFoundGene().getId();
    // String coexpressedResult = coexpressedResult( cvo );
    //
    // if ( results.containsKey( queryGene ) )
    // // yes, in map, then add the found gene object to the Collection value corresponding to the key (query
    // // gene)
    // results.get( queryGene ).add( coexpressedResult );
    // else {
    // // no, not in map, then add the query gene as new key and insert the first found gene for the Collection
    // // value
    // results.put( queryGene, new HashSet<String>() );
    // results.get( queryGene ).add( coexpressedResult );
    // }
    //
    // }
    // return results;
    // }
    private String coexpressedResult( CoexpressionValueObjectExt cvo ) {
        String queryGene = cvo.getQueryGene().getOfficialSymbol();
        String foundGene = cvo.getFoundGene().getOfficialSymbol();
        // Integer numDSTested = cvo.getNumTestedIn();
        // StringBuilder buf = new StringBuilder();
        // String[] fields;
        if ( details == true ) {
            return cvo.toString();
            // // only positive correlation
            // if ( cvo.getPosLinks() > 0 && cvo.getNegLinks() <= 0 ) {
            // buf.append( cvo.getPosLinks() + "\t" );
            // }
            // // only negative correlation
            // else if ( cvo.getNegLinks() > 0 && cvo.getPosLinks() <= 0 ) {
            // buf.append( "\t" + cvo.getNegLinks() );
            // }
            // // negative and positive
            // else
            // buf.append( cvo.getPosLinks() + "\t" + cvo.getNegLinks() );
            // fields = new String[] { cvo.getQueryGene().getOfficialSymbol(), foundGene, numDSTested.toString(),
            // buf.toString() };
            // return StringUtils.join( fields, "\t" );
        } else {
            return queryGene + "\t" + foundGene;
        }
    }

    /**
     * @param coexpressionList.keySet()
     */
    // private void printCoexpressionResults( Map<String, Collection<String>> coexpressionList ) {
    //
    // for ( String queryGene : coexpressionList.keySet() ) {
    // for ( String result : coexpressionList.get( queryGene ) )
    // System.out.println( result );
    //
    // }
    //
    // }
    /**
     * prints out, in a file, the degree count for each query gene. Degree count = number of coexpressed genes
     * 
     * @param coexpressionList
     * @throws IOException
     * @throws InterruptedException
     */
    // private void printDegreeCount( Map<String, Collection<String>> coexpressionList ) throws IOException {
    // try {
    // String dir = "../";
    // String query = "";
    // if ( queryGenesOnly == true ) query = "-qgOnly";
    // String outFile = getOptionValue( 'g' ) + "Stringency" + getOptionValue( 's' ) + query + "-Summary.txt";
    // BufferedWriter out = new BufferedWriter( new FileWriter( dir + outFile, true ) );
    // Collection<String> foundResults;
    // Formatter fmt;
    //
    // Set<String> querySet = coexpressionList.keySet();
    // for ( String queryGene : querySet ) {
    //
    // foundResults = coexpressionList.get( queryGene );
    //
    // // cut out all the details; reassigns foundResults; original foundResults - garbaged collected?
    // if ( details == true ) foundResults = trimToFoundGenesOnly( foundResults );
    //
    // out.write( queryGene + "\t" );
    // // total degree
    // int degree = foundResults.size();
    // out.write( degree + "\t" );
    //
    // // # of coexpressed genes within query gene set
    // // make shallow copy
    // Collection<String> intersection = new TreeSet<String>( foundResults );
    //
    // intersection.retainAll( rawGeneList );
    // out.write( intersection.size() + "\t" );
    //
    // // # of coexpressed genes not within query gene set
    // Collection<String> difference = new TreeSet<String>( foundResults );
    //
    // difference.removeAll( rawGeneList );
    // out.write( difference.size() + "\t" );
    //
    // // output results
    // fmt = new Formatter();
    // String percentQueryGenes = ( fmt.format( "%.4f", ( intersection.size() / ( double ) degree ) * 100 ) )
    // .toString();
    // out.write( percentQueryGenes + "\n" );
    //
    // }
    //
    // out.close();
    // } catch ( IOException e ) {
    //
    // e.printStackTrace();
    // }
    // }
    /**
     * makes new copy array with trimmed results
     * 
     * @param foundResults
     * @return
     */
    // private Collection<String> trimToFoundGenesOnly( Collection<String> foundResults ) {
    // Collection<String> trimmedResult = new TreeSet<String>();
    // for ( String result : foundResults )
    // trimmedResult.add( StringUtils.substringBefore( result, "\t" ) );
    // return trimmedResult;
    // }
    // private boolean isNumber( String string ) {
    // for ( int i = 0; i < string.length(); i++ ) {
    // if ( !Character.isDigit( string.charAt( i ) ) ) return false;
    // }
    // return true;
    // }
}
