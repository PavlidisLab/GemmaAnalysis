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

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.lang.time.StopWatch;
import org.springframework.util.StringUtils;

import ubic.gemma.analysis.service.ExpressionDataMatrixService;
import ubic.gemma.apps.ExpressionExperimentManipulatingCLI;
import ubic.gemma.datastructure.matrix.ExpressionDataDoubleMatrix;
import ubic.gemma.expression.experiment.service.ExpressionExperimentService;
import ubic.gemma.model.expression.arrayDesign.ArrayDesign;
import ubic.gemma.model.expression.arrayDesign.ArrayDesignService;
import ubic.gemma.model.expression.biomaterial.BioMaterial;
import ubic.gemma.model.expression.designElement.CompositeSequence;
import ubic.gemma.model.expression.designElement.CompositeSequenceService;
import ubic.gemma.model.expression.experiment.BioAssaySet;
import ubic.gemma.model.expression.experiment.ExpressionExperiment;
import ubic.gemma.model.genome.Gene;

/**
 * Writes the expression level for a select group of genes for each sample.
 * 
 * @author ptan
 * @version $Id$
 */
public class GeneExpressionWriterCLI extends ExpressionExperimentManipulatingCLI {

    public static void main( String[] args ) {
        GeneExpressionWriterCLI cli = new GeneExpressionWriterCLI();
        Exception e = cli.doWork( args );
        if ( e != null ) log.error( e.getMessage() );
    }

    private ArrayDesignService adService;

    private ExpressionDataMatrixService expressionDataMatrixService;

    private String outFilePrefix;
    private String[] queryGeneSymbols;

    private String queryGeneFile;

    private CompositeSequenceService csService;

    /*
     * (non-Javadoc)
     * 
     * @see ubic.gemma.util.AbstractCLI#getCommandName()
     */
    @Override
    public String getCommandName() {
        return "writeGeneData";
    }

    /**
     * Get a gene to composite sequence map // FIXME This corresponds to an existing service method? Yes, it does -PP
     * 
     * @param css
     * @return gene to composite sequences map
     */
    public Map<Gene, Collection<CompositeSequence>> getGene2CsMap( Collection<CompositeSequence> css ) {
        Map<CompositeSequence, Collection<Gene>> cs2gene = csService.getGenes( css );
        // filter for specific cs 2 gene
        for ( Iterator<Map.Entry<CompositeSequence, Collection<Gene>>> it = cs2gene.entrySet().iterator(); it.hasNext(); ) {
            Map.Entry<CompositeSequence, Collection<Gene>> entry = it.next();
            Collection<Gene> genes = entry.getValue();
            if ( genes.size() > 1 ) it.remove();
        }

        Map<Gene, Collection<CompositeSequence>> gene2css = new HashMap<Gene, Collection<CompositeSequence>>();
        for ( Map.Entry<CompositeSequence, Collection<Gene>> entry : cs2gene.entrySet() ) {
            CompositeSequence cs = entry.getKey();
            Collection<Gene> genes = entry.getValue();
            for ( Gene gene : genes ) {
                Collection<CompositeSequence> c = gene2css.get( gene );
                if ( c == null ) {
                    c = new HashSet<CompositeSequence>();
                    gene2css.put( gene, c );
                }
                c.add( cs );
            }
        }
        return gene2css;
    }

    public Collection<Gene> getQueryGenes() throws IOException {
        Collection<Gene> genes = new HashSet<Gene>();
        if ( queryGeneFile != null ) genes.addAll( readGeneListFile( queryGeneFile, taxon ) );
        if ( queryGeneSymbols != null ) {
            for ( int i = 0; i < queryGeneSymbols.length; i++ ) {
                genes.add( findGeneByOfficialSymbol( queryGeneSymbols[i], taxon ) );
            }
        }

        return genes;
    }

    @Override
    protected void buildOptions() {
        super.buildOptions();

        addForceOption( null );

        OptionBuilder.hasArg();
        OptionBuilder.withDescription( "Query file containing list of gene official symbols" );
        OptionBuilder.withArgName( "File name" );
        OptionBuilder.withLongOpt( "queryGeneFile" );
        Option queryGeneFileOption = OptionBuilder.create();
        addOption( queryGeneFileOption );
        OptionBuilder.hasArgs();
        OptionBuilder.withArgName( "Gene symbol(s)" );
        OptionBuilder.withDescription( "The query gene symbol(s), comma separated" );
        OptionBuilder.withValueSeparator( ',' );
        OptionBuilder.withLongOpt( "queryGene" );
        Option queryGeneOption = OptionBuilder.create();
        addOption( queryGeneOption );

        OptionBuilder.hasArg();
        OptionBuilder.withArgName( "outfile" );
        OptionBuilder.withDescription( "Output filename prefix" );
        OptionBuilder.withLongOpt( "outfile" );
        // Option outFilePrefixOption = OptionBuilder.create( 'o' );
        Option outFilePrefixOption = OptionBuilder.create( 'o' );
        addOption( outFilePrefixOption );

    }

    @Override
    protected Exception doWork( String[] args ) {

        StopWatch timer = new StopWatch();
        timer.start();

        String HEADER = "HEADER";

        processCommandLine( args );

        Collection<Gene> genes;
        try {
            genes = getQueryGenes();
        } catch ( IOException e ) {
            return e;
        }

        StringBuffer sb = new StringBuffer();
        for ( Gene gene : genes ) {
            sb.append( gene.toString() );
            sb.append( ", " );
        }
        log.info( genes.size() + " genes: " + sb.toString() );

        sb = new StringBuffer();
        for ( BioAssaySet bas : this.expressionExperiments ) {
            sb.append( bas.toString() );
            sb.append( ", " );
        }
        log.info( this.expressionExperiments.size() + " experiments: " + sb.toString() );

        Map<Object, StringBuffer> outBuffs = new HashMap<>();

        String fileName = outFilePrefix + ".txt";

        log.info( "Output filename " + fileName );

        for ( BioAssaySet bas : this.expressionExperiments ) {
            ExpressionExperiment ee = ( ExpressionExperiment ) bas;
            Collection<ArrayDesign> ads = eeService.getArrayDesignsUsed( ee );
            Collection<CompositeSequence> css = new HashSet<CompositeSequence>();
            for ( ArrayDesign ad : ads ) {
                css.addAll( adService.getCompositeSequences( ad ) );
            }

            log.info( "====================================" );
            log.info( "Experiment " + ee + "; Array Design " + ads + "; Composite sequences " + css.size() );

            Map<Gene, Collection<CompositeSequence>> gene2css = getGene2CsMap( css );
            // ExpressionDataDoubleMatrix dataMatrix = expressionDataMatrixService.getFilteredMatrix( ee, filterConfig
            // );
            ExpressionDataDoubleMatrix dataMatrix = expressionDataMatrixService.getProcessedExpressionDataMatrix( ee );

            // store values inside a buffer
            if ( !outBuffs.containsKey( HEADER ) ) {
                StringBuffer hb = new StringBuffer();
                hb.append( "Gene\tProbe\tID" );
                outBuffs.put( HEADER, hb );
            }
            for ( BioMaterial bm : dataMatrix.getMatrix().getColNames() ) {
                // bm.getFactorValues()
                // bm.getCharacteristics()
                // outBuffs.get( HEADER )
                // .append(
                // "\t" + bm.getName() + "."
                // + StringUtils.collectionToDelimitedString( bm.getFactorValues(), "," ) );
                outBuffs.get( HEADER ).append(
                        "\t" + ee.getShortName() + "." + bm.getName() + "."
                                + StringUtils.collectionToDelimitedString( bm.getCharacteristics(), "," ) + "."
                                + StringUtils.collectionToDelimitedString( bm.getFactorValues(), "," ) );
            }

            for ( Gene gene : genes ) {
                log.debug( " Getting component sequence for gene " + gene );
                Collection<CompositeSequence> c = gene2css.get( gene );

                if ( c == null ) {
                    log.error( "No composite sequences found for gene " + gene );
                    continue;
                }

                for ( CompositeSequence cs : c ) {

                    Double[] row = dataMatrix.getRow( cs );
                    if ( row == null ) {
                        log.error( "Cannot get data from data matrix for " + gene.getOfficialSymbol() + " ("
                                + cs.getName() + ")" );
                        // continue;
                        row = new Double[dataMatrix.getMatrix().columns()];
                    }

                    if ( !outBuffs.containsKey( cs ) ) {
                        StringBuffer gb = new StringBuffer();
                        gb.append( gene.getOfficialSymbol() + "\t" + cs.getName() + "\t" + cs.getId() );
                        outBuffs.put( cs, gb );
                    }

                    StringBuffer buf = new StringBuffer();
                    for ( Double d : row ) {
                        if ( d == null )
                            buf.append( "NA" );
                        else
                            buf.append( d );
                        buf.append( "\t" );
                    }
                    buf.deleteCharAt( buf.length() - 1 );
                    outBuffs.get( cs ).append( "\t" + buf.toString() );
                }
            }

        }

        // Output to file
        try (PrintWriter out = new PrintWriter( new FileWriter( fileName ) );) {
            out.println( outBuffs.get( HEADER ) );
            for ( Object key : outBuffs.keySet() ) {
                if ( key.equals( HEADER ) ) {
                    continue;
                }
                out.println( outBuffs.get( key ) );
            }
        } catch ( IOException e ) {
            return e;
        }

        log.info( "Done. Wrote " + genes.size() + " genes and " + ( outBuffs.keySet().size() - 1 )
                + " composite sequences in " + this.expressionExperiments.size() + " experiments which took "
                + timer.getTime() + " ms. Output file " + fileName );

        return null;
    }

    protected Map<String, String> getGeneIdPair2NameMap( Collection<Gene> queryGenes, Collection<Gene> targetGenes ) {
        Map<String, String> map = new HashMap<String, String>();
        for ( Gene qGene : queryGenes ) {
            String qName = ( qGene.getOfficialSymbol() != null ) ? qGene.getOfficialSymbol() : qGene.getId().toString();
            for ( Gene tGene : targetGenes ) {
                String tName = ( tGene.getOfficialSymbol() != null ) ? tGene.getOfficialSymbol() : tGene.getId()
                        .toString();
                map.put( qGene.getId() + ":" + tGene.getId(), qName + ":" + tName );
            }
        }
        return map;
    }

    @Override
    protected void processOptions() {
        super.processOptions();
        if ( hasOption( "queryGeneFile" ) ) queryGeneFile = getOptionValue( "queryGeneFile" );
        if ( hasOption( "queryGene" ) ) queryGeneSymbols = getOptionValues( "queryGene" );

        outFilePrefix = getOptionValue( 'o' );

        log.info( "outFilePrefix " + outFilePrefix );

        adService = getBean( ArrayDesignService.class );
        eeService = getBean( ExpressionExperimentService.class );
        expressionDataMatrixService = getBean( ExpressionDataMatrixService.class );
        csService = getBean( CompositeSequenceService.class );
    }

}
