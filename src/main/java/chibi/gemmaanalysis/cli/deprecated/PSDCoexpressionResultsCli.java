package chibi.gemmaanalysis.cli.deprecated;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;

import ubic.gemma.analysis.expression.coexpression.CoexpressionValueObjectExt;
import ubic.gemma.analysis.expression.coexpression.GeneCoexpressionService;
import ubic.gemma.model.analysis.expression.coexpression.GeneCoexpressionAnalysisService;
import ubic.gemma.model.association.coexpression.Gene2GeneCoexpressionService;
import ubic.gemma.model.genome.Gene;
import ubic.gemma.model.genome.Taxon;
import ubic.gemma.genome.gene.service.GeneService;
import ubic.gemma.genome.taxon.service.TaxonService;
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

        Option stringencyOption = OptionBuilder.hasArg().withArgName( "Stringency" )
                .withDescription( "The stringency value: Defaults to " + DEFAULT_STRINGINCY )
                .withLongOpt( "stringency" ).create( 's' );
        addOption( stringencyOption );

        Option taxonOption = OptionBuilder.hasArg().isRequired().withArgName( "Taxon" )
                .withDescription( "The name of the taxon." ).withLongOpt( "taxon" ).create( 't' );
        addOption( taxonOption );

        Option analysisOption = OptionBuilder.hasArg().withArgName( "CoexpressionAnalysis" )
                .withDescription( "Gene coexpression analysis ID to use." ).withLongOpt( "analysis" ).create( 'a' );
        addOption( analysisOption );

        Option queryOption = OptionBuilder.hasArg().withArgName( "QueryGenesOnly" )
                .withDescription( "Link analysis on query genes only?" ).withLongOpt( "query" ).create( 'q' );
        addOption( queryOption );

        Option randomOption = OptionBuilder.hasArg().withArgName( "RandomFileOrCmd" )
                .withDescription( "all = all known genes in taxon will be used" ).withLongOpt( "random" ).create( 'r' );
        addOption( randomOption );

        Option detailOption = OptionBuilder.hasArg().withArgName( "ExtraDetails" )
                .withDescription( "true = extra details, false = just coexpressed genes" ).withLongOpt( "detail" )
                .create( 'd' );
        addOption( detailOption );
    }

    @Override
    protected Exception doWork( String[] args ) {
        Exception err = processCommandLine( "Gene 2 Gene Coexpression Results ", args );
        if ( err != null ) return err;

        try {
            Collection<Gene> genes = null;
            if ( allGenes == true )
                genes = geneService.loadKnownGenes( taxon );
            else
                genes = getGenes();

            // Batch up the genes to reduce memory footprint

            int count = 0;
            int CHUNK_SIZE = 100;

            Collection<Gene> batch = new HashSet<Gene>();

            for ( Gene g : genes ) {
                batch.add( g );
                count++;
                if ( count % CHUNK_SIZE == 0 ) {
                    outputCoexpressionResults( batch );
                    batch.clear();
                }
            }

            if ( batch.size() > 0 ) {
                outputCoexpressionResults( batch );
            }

        } catch ( IOException e ) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

        return null;
    }

    @Override
    protected void processOptions() {
        super.processOptions();

        geneService = getBean( GeneService.class );
        taxonService = getBean( TaxonService.class );
        gene2GeneCoexpressionService = getBean( Gene2GeneCoexpressionService.class );
        geneCoexpressionService = getBean( GeneCoexpressionService.class );
        geneCoexpressionAnalysisService = getBean( GeneCoexpressionAnalysisService.class );

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
    private Collection<Gene> getGenes() throws IOException {
        // Collection<String> rawLinesInFile = readGeneListFile( inFile );
        Collection<Gene> genes = new ArrayList<Gene>();

        Iterator<String> linesIt = rawGeneList.iterator();
        while ( linesIt.hasNext() ) {
            String s = linesIt.next();
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

    private void outputCoexpressionResults( Collection<Gene> geneList ) {
        log.debug( "Running Gene Coexpression..." );

        Collection<CoexpressionValueObjectExt> cmvo = geneCoexpressionService.coexpressionSearchQuick( analysisID,
                geneList, stringency, 0, queryGenesOnly, true );
        log.debug( "Printing results..." );
        printCoexpressionResults( cmvo );
    }

    private void printCoexpressionResults( Collection<CoexpressionValueObjectExt> cvoExtCol ) {
        createColumnHeadings();

        for ( CoexpressionValueObjectExt cvo : cvoExtCol )
            System.out.println( coexpressedResult( cvo ) );

    }

    private void createColumnHeadings() {
        // output column headings for main coexpression results file
        if ( details == true )
            System.out.println( "Query_Gene\tCoexpressed_Gene\tSupport\tSign" );
        else
            System.out.println( "Query_Gene\tCoexpressed_Gene" );
    }

    private String coexpressedResult( CoexpressionValueObjectExt cvo ) {
        String queryGene = cvo.getQueryGene().getOfficialSymbol();
        String foundGene = cvo.getFoundGene().getOfficialSymbol();
        // Integer numDSTested = cvo.getNumTestedIn();
        // StringBuilder buf = new StringBuilder();
        // String[] fields;
        if ( details == true ) {
            return cvo.toString();

        }
        return queryGene + "\t" + foundGene;

    }

}
