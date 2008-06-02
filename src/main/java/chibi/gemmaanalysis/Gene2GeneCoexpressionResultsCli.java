package chibi.gemmaanalysis;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.gemma.analysis.expression.coexpression.CoexpressionMetaValueObject;
import ubic.gemma.analysis.expression.coexpression.CoexpressionValueObjectExt;
import ubic.gemma.analysis.expression.coexpression.GeneCoexpressionService;
import ubic.gemma.model.association.coexpression.Gene2GeneCoexpressionService;
import ubic.gemma.model.genome.Gene;
import ubic.gemma.model.genome.Taxon;
import ubic.gemma.model.genome.TaxonService;
import ubic.gemma.model.genome.gene.GeneService;
import ubic.gemma.util.AbstractSpringAwareCLI;

public class Gene2GeneCoexpressionResultsCli extends AbstractSpringAwareCLI {

    private static final int DEFAULT_STRINGINCY = 2;

    private int stringency;

    private Collection<Gene> geneList;

    private Taxon taxon;

    private int SUBSETSIZE;
    
    private boolean queryGenesOnly;

    private static Log log = LogFactory.getLog( Gene2GeneCoexpressionResultsCli.class );

    protected GeneService geneService;

    protected Gene2GeneCoexpressionService gene2GeneCoexpressionService;

    protected TaxonService taxonService;

    protected GeneCoexpressionService geneCoexpressionService;

    public static void main( String args[] ) {
        Gene2GeneCoexpressionResultsCli run = new Gene2GeneCoexpressionResultsCli();
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
        Option geneFileOption = OptionBuilder.hasArg().withArgName( "Gene List File Name" ).withDescription(
                "A text file that contains a list of gene symbols, with one gene symbol on each line" ).withLongOpt(
                "geneFile" ).create( 'g' );
        addOption( geneFileOption );

        Option stringencyOption = OptionBuilder.hasArg().withArgName( "Stringency" ).withDescription(
                "The stringency value: Defaults to " + DEFAULT_STRINGINCY ).withLongOpt( "stringency" ).create( 's' );
        addOption( stringencyOption );

        Option taxonOption = OptionBuilder.hasArg().isRequired().withArgName( "Taxon" ).withDescription(
                "The name of the taxon." ).withLongOpt( "taxon" ).create( 't' );
        addOption( taxonOption );
        
        Option queryOption = OptionBuilder.hasArg().withArgName( "QueryGenesOnly" ).withDescription(
        "Link analysis on query genes only?" ).withLongOpt( "query" ).create( 'q' );
        addOption( queryOption );

        Option subsetOption = OptionBuilder.hasArg().withArgName( "GeneSubsetSize" ).withDescription(
                "The size of the gene subset to use.  true/(1) or false/(0)" ).withLongOpt( "subsetSize" ).create( 'z' );
        addOption( subsetOption );

        Option randomOption = OptionBuilder.hasArg().withArgName( "RandomFileOrCmd" ).withDescription(
                "all = all known genes in taxon will be used" ).withLongOpt( "random" ).create( 'r' );
        addOption( randomOption );
    }

    @Override
    protected Exception doWork( String[] args ) {
        Exception err = processCommandLine( "Gene 2 Gene Coexpression Results ", args );
        if ( err != null ) return err;

        try {
            outputCoexpressionResults( geneList, taxon, stringency, queryGenesOnly );
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

        // taxon option - t
        if ( hasOption( 't' ) ) {
            String taxonName = this.getOptionValue( 't' );
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

        // random option - r
        // if "all" is given, then use all genes in a given taxon        
        if ( this.hasOption( 'r' ) ) {
            String rOption = this.getOptionValue( 'r' );
            assert taxon != null;
            if ( rOption.equals( "all" ) ) geneList = geneService.loadKnownGenes( taxon );
        } 
        //no random option
        else {
            assert taxon != null;
            try {
                geneList = getGeneList( this.getOptionValue( 'g' ), taxon );
            } catch ( IOException e ) {
                throw new RuntimeException( e );
            }
        }

        // stringency option - s
        // default stringency = 2
        stringency = DEFAULT_STRINGINCY;
        if ( this.hasOption( 's' ) ) {
            stringency = Integer.parseInt( this.getOptionValue( 's' ) );
        }

        // subset gene list size to use
        SUBSETSIZE = 25; // default
        if ( this.hasOption( 'z' ) ) 
            SUBSETSIZE = Integer.parseInt( this.getOptionValue( 'z' ) );
        
        //link analysis on genelist only?
        queryGenesOnly = false; // default
        if ( this.hasOption( 'q' ) ) {
            String query = this.getOptionValue( 'q' );
            if (query.equals( "true") || query.equals( "1") )
                queryGenesOnly =  true;           
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
    private Collection<String> readGeneListFile( String inFile ) throws IOException {
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
        return lines;
    }

    /**
     * Helper method to read in the lines in a text file and store in a collection
     * 
     * @param inFile
     * @return - ArrayList storing the lines in a text file
     * @throws IOException
     */
    private Collection<Gene> getGeneList( String inFile, Taxon taxon ) throws IOException {
        Collection<String> rawLinesInFile = readGeneListFile( inFile );
        Collection<Gene> genes = new ArrayList<Gene>();

        Iterator linesIt = rawLinesInFile.iterator();
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
        return genes;
    }

    //    /**
    //     * Uses precomputed file of IDs; calls helper method, <code>readGeneListFile()</code>
    //     * 
    //     * @param inFile - Precomputed file of IDs
    //     * @param i - number of random genes to use
    //     * @return - list of 'i' random number of genes from precomputed list of IDs
    //     * @throws IOException
    //     */
    //    private Collection<Gene> getRandomKnownGenes( String inFile, int i ) throws IOException {
    //        Object[] rawLinesInFile = readGeneListFile( inFile ).toArray();
    //        Collection<Gene> genes = new HashSet<Gene>();
    //        String line;
    //        int geneCount = 0;
    //        SecureRandom randomSeed;
    //        try {
    //            randomSeed = SecureRandom.getInstance( "SHA1PRNG" );
    //
    //            while ( geneCount < 1000 ) {
    //
    //                randomSeed.generateSeed( 10 );
    //                int randomInt = randomSeed.nextInt( rawLinesInFile.length );
    //
    //                line = ( String ) rawLinesInFile[randomInt];
    //
    //                Gene gene = geneService.load( Long.parseLong( line ) );
    //                if ( gene == null ) {
    //                    log.error( "ERROR: Cannot find genes for ID: " + line );
    //                    continue;
    //                }
    //                geneService.thaw( gene );
    //
    //                // because we are using a HashSet, no duplicates are allowed - which is what we want
    //                genes.add( gene );
    //                geneCount++;
    //            }
    //
    //        } catch ( NoSuchAlgorithmException e ) {
    //            // TODO Auto-generated catch block
    //            e.printStackTrace();
    //        }
    //
    //        return genes;
    //    }

    @SuppressWarnings("unchecked")
    private void outputCoexpressionResults( Collection<Gene> geneList, Taxon taxon, int stringency, boolean queryGenesOnly ) throws IOException {
        int geneCount = 0;
       
       createColumnHeadings();
       
        while ( geneCount < geneList.size() ) {

            // use a subset of genes
            Collection<Gene> geneSubset = new ArrayList<Gene>( SUBSETSIZE );
            Object[] geneArray = geneList.toArray();
            for ( int i = 0; i < SUBSETSIZE; i++ ) {
                if ( geneCount == geneList.size() ) break;
                geneSubset.add( ( Gene ) geneArray[geneCount] );
                geneCount++;
            }
            // use subset list of genes, size SUBSETSIZE
            CoexpressionMetaValueObject cmvo = geneCoexpressionService.getCannedAnalysisResults( 717L, geneSubset,
                    stringency, 0, queryGenesOnly );
            Collection<CoexpressionValueObjectExt> cvoExtCol = cmvo.getKnownGeneResults();

            Map<String, Collection<String>> coexpressionList = organizeCoexpressionValueObjectResults( cvoExtCol );

            printDegreeCount( coexpressionList );

            printCoexpressionResults( coexpressionList );

            log.info( "Genes analyzed so far: " + geneCount );
        }
        log.info( "Total number of genes analyzed" + geneCount );

    }

    private void createColumnHeadings() throws IOException {
 System.out.println( "Query_Gene\tCoexpressed_Gene\tDatasets_Tested\t+Correlation\t-Correlation" );
        
        String dir = "../";
        String query = "";
        if (queryGenesOnly == true)
              query = "-qgOnly";
        String outFile = getOptionValue('g')+"Stringency" + getOptionValue( 's' )+ query+"-Summary.txt";
        BufferedWriter out = new BufferedWriter( new FileWriter( dir + outFile, true ) );
        out.write( "Query_Gene\tTotal_Degree\n" );
        out.close();
    }

    /**
     * map of queryGene to coexpressed gene, # of datasets, etc
     * @param cvoExtCol
     * @return mapping of the query gene to all expressed
     */
    private Map<String, Collection<String>> organizeCoexpressionValueObjectResults(
            Collection<CoexpressionValueObjectExt> cvoExtCol ) {
        Map<String, Collection<String>> results = new TreeMap<String, Collection<String>>();
        for ( CoexpressionValueObjectExt cvo : cvoExtCol ) {
            // see if query gene official symbol is already a key in the map
            String queryGene = cvo.getQueryGene().getOfficialSymbol();
            String coexpressedResult = coexpressedResult( cvo );

            if ( results.containsKey( queryGene ) )
                // yes, in map, then add the found gene object to the Collection value corresponding to the key (query
                // gene)
                results.get( queryGene ).add( coexpressedResult );
            else {
                // no, not in map, then add the query gene as new key and insert the first found gene for the Collection
                // value
                results.put( queryGene, new TreeSet<String>() );
                results.get( queryGene ).add( coexpressedResult );
            }

        }
        return results;
    }

    private String coexpressedResult( CoexpressionValueObjectExt cvo ) {
        String queryGene = cvo.getQueryGene().getOfficialSymbol();
        String foundGene = cvo.getFoundGene().getOfficialSymbol();
        Integer numDSTested = cvo.getNumTestedIn();
        StringBuilder buf = new StringBuilder();

        //only positive correlation
        if ( cvo.getPosLinks() > 0 && cvo.getNegLinks() <= 0 ) {
            buf.append( cvo.getPosLinks() + "\t" );
        }
        //only negative correlation
        else if ( cvo.getNegLinks() > 0 && cvo.getPosLinks() <= 0 ) {
            buf.append( "\t" + cvo.getNegLinks() );
        }
        //negative and positive
        else
            buf.append( cvo.getPosLinks() + "\t" + cvo.getNegLinks() );
        String[] fields = new String[] { queryGene, foundGene, numDSTested.toString(), buf.toString() };
        return StringUtils.join( fields, "\t" );

    }

    /**
     * 
     * @param coexpressionList
     */
    private void printCoexpressionResults( Map<String, Collection<String>> coexpressionList ) {
        
        for ( String queryGene : coexpressionList.keySet() ) {
            // String outFile = "hub-coexpression-canned-results-queryGenesOnly-false.txt";
            // try {
            // BufferedWriter out = new BufferedWriter(new FileWriter(outFile));
            // out.write( "\n"+queryGene +"\n");
            // out.write( "=============================================================" );
            for ( String coexpressedResult : coexpressionList.get( queryGene ) ) {
                System.out.println( coexpressedResult );
                // out.write( foundGene.getOfficialSymbol()+"\n" );
            }
            // out.close();

            // } catch ( IOException e ) {
            // // TODO Auto-generated catch block
            // e.printStackTrace();
            // }
        }
    }

    /**
     * prints out, in a file, the degree count for each query gene.  Degree count = number of coexpressed genes
     * @param coexpressionList
     * @throws IOException
     */
    private void printDegreeCount( Map<String, Collection<String>> coexpressionList ) throws IOException {
        try {
            String dir = "../";
            String query = "";
            if (queryGenesOnly == true)
                  query = "-qgOnly";
            String outFile = getOptionValue('g')+"Stringency" + getOptionValue( 's' )+ query+"-Summary.txt";
            BufferedWriter out = new BufferedWriter( new FileWriter( dir + outFile, true ) );

            for ( String queryGene : coexpressionList.keySet() ) {

                int degree = coexpressionList.get( queryGene ).size();
                out.write(queryGene + "\t" + degree +"\n" );

            }
            out.close();
        } catch ( IOException e ) {

            e.printStackTrace();
        }
    }

    //    private boolean isNumber( String string ) {
    //        for ( int i = 0; i < string.length(); i++ ) {
    //            if ( !Character.isDigit( string.charAt( i ) ) ) return false;
    //        }
    //        return true;
    //    }
}
