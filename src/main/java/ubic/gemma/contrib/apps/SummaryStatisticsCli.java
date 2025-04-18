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
package ubic.gemma.contrib.apps;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.springframework.beans.factory.annotation.Autowired;
import ubic.basecode.dataStructure.matrix.CompressedSparseDoubleMatrix;
import ubic.gemma.cli.util.AbstractCLI;
import ubic.gemma.model.expression.arrayDesign.ArrayDesign;
import ubic.gemma.model.expression.designElement.CompositeSequence;
import ubic.gemma.model.expression.experiment.ExpressionExperiment;
import ubic.gemma.model.genome.Gene;
import ubic.gemma.model.genome.Taxon;
import ubic.gemma.model.genome.biosequence.BioSequence;
import ubic.gemma.persistence.service.expression.arrayDesign.ArrayDesignService;
import ubic.gemma.persistence.service.expression.designElement.CompositeSequenceService;
import ubic.gemma.persistence.service.expression.experiment.ExpressionExperimentService;
import ubic.gemma.persistence.service.genome.gene.GeneService;
import ubic.gemma.persistence.service.genome.taxon.TaxonService;

import javax.annotation.Nullable;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigInteger;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

/**
 * Computing different statistics about the database to assist in computing probabilities
 *
 * @author pavlidis
 * @version $Id: SummaryStatistics.java,v 1.14 2015/11/12 19:37:11 paul Exp $
 */
public class SummaryStatisticsCli extends AbstractCLI {

    private static final int MAX_EXPS = 5;

    private static final int MAX_GENES = 100000;

    @Autowired
    private ExpressionExperimentService expressionExperimentService;
    @Autowired
    private GeneService geneService;
    @Autowired
    private TaxonService taxonService;
    @Autowired
    private CompositeSequenceService compositeSequenceService;
    @Autowired
    private ArrayDesignService adService;

    private String taxonName;
    private String outFileName;

    @Override
    public String getCommandName() {
        return "summaryStatistics";
    }

    @Override
    protected void buildOptions( Options options ) {
        options.addOption( Option.builder( "o" ).hasArg().argName( "outFile" ).desc( "Output file" ).build() );
        options.addOption( Option.builder( "t" ).hasArg().argName( "taxon" ).desc( "Taxon common name (e.g., human)" ).build() );
    }

    @Override
    protected void processOptions( CommandLine c ) {
        if ( c.hasOption( 't' ) ) {
            this.taxonName = c.getOptionValue( 't' );
        }
        if ( c.hasOption( 'o' ) ) {
            this.outFileName = c.getOptionValue( 'o' );
        }
    }

    @Override
    protected void doWork() {
        Taxon taxon = taxonService.findByCommonName( taxonName );
        genesPerProbe( taxon );
    }

    /**
     * For each gene, count how many expression experiments it appears in.
     */
    private void geneOccurrenceDistributions( Taxon taxon ) {

        Map<Long, Integer> counts = new HashMap<>();

        // get all expression experiments

        Collection<ExpressionExperiment> eeColl = expressionExperimentService.loadAll();

        int i = 0;
        for ( ExpressionExperiment experiment : eeColl ) {
            if ( i > MAX_EXPS ) break;
            Taxon eeTax = expressionExperimentService.getTaxon( experiment );
            if ( eeTax == null || !eeTax.equals( taxon ) ) continue;
            Collection<ArrayDesign> ads = expressionExperimentService.getArrayDesignsUsed( experiment );

            // only count each gene once per data set.
            Collection<Long> seenids = new HashSet<>();

            for ( ArrayDesign design : ads ) {
                log.info( i + " " + design );
                Collection<Object[]> vals = compositeSequenceService.getRawSummary( design, -1 );
                log.info( "Got " + vals.size() + " reports" );
                for ( Object[] objects : vals ) {

                    BigInteger geneidi = ( BigInteger ) objects[10];
                    if ( geneidi == null ) {
                        continue;
                    }
                    Long geneid = geneidi.longValue();

                    if ( seenids.contains( geneid ) ) continue;

                    if ( counts.get( geneid ) == null ) {
                        counts.put( geneid, 0 );
                    }
                    counts.put( geneid, counts.get( geneid ) + 1 );
                    seenids.add( geneid );
                }
            }
            i++;
        }

        for ( Long l : counts.keySet() ) {
            System.out.println( l + "\t" + counts.get( l ) );
        }
    }

    /**
     * For each pair of genes, count how many expression experiments both appear in.
     */
    private void genePairOccurrenceDistributions( Taxon taxon ) {

        Collection<ExpressionExperiment> eeColl = expressionExperimentService.loadAll();

        CompressedSparseDoubleMatrix<Long, Long> mat = new CompressedSparseDoubleMatrix<>( MAX_GENES,
                MAX_GENES );

        int numEEs = 0;
        for ( ExpressionExperiment experiment : eeColl ) {
            if ( numEEs > MAX_EXPS ) break;
            Taxon eeTax = expressionExperimentService.getTaxon( experiment );
            if ( eeTax == null || !eeTax.equals( taxon ) ) continue;
            Collection<ArrayDesign> ads = expressionExperimentService.getArrayDesignsUsed( experiment );

            // only count each gene once per data set.
            Collection<Long> seenids = new HashSet<>();

            for ( ArrayDesign design : ads ) {

                Collection<Object[]> vals = compositeSequenceService.getRawSummary( design, -1 );
                log.info( numEEs + " " + design + "Got " + vals.size() + " reports" );

                for ( Object[] objects : vals ) {

                    BigInteger geneidi = ( BigInteger ) objects[10];
                    if ( geneidi == null ) {
                        continue;
                    }
                    Long geneid = geneidi.longValue();

                    if ( seenids.contains( geneid ) ) continue;

                    if ( !mat.containsRowName( geneid ) ) {
                        mat.addRowName( geneid );
                    }

                    int outerIndex = mat.getRowIndexByName( geneid );

                    int j = 0;
                    for ( Object[] ojbB : vals ) {

                        BigInteger geneBidi = ( BigInteger ) ojbB[10];
                        if ( geneBidi == null || geneBidi.equals( geneidi ) ) {
                            continue;
                        }
                        Long geneBid = geneBidi.longValue();
                        if ( seenids.contains( geneBid ) ) continue;
                        int innerIndex;
                        if ( !mat.containsColumnName( geneBid ) ) {
                            mat.addColumnName( geneBid );
                            innerIndex = mat.getColIndexByName( geneBid );
                            mat.set( outerIndex, innerIndex, 0.0 ); // initialize
                        }

                        innerIndex = mat.getColIndexByName( geneBid );
                        mat.set( outerIndex, innerIndex, mat.get( outerIndex, innerIndex ) + 1 );

                        if ( mat.columns() > MAX_GENES ) {
                            log.warn( "Too many genes!" );
                            break;
                        }
                        j++;
                        if ( j > 1000 ) break;
                    }
                    seenids.add( geneid );

                    if ( mat.rows() > MAX_GENES ) {
                        break;
                    }

                }

            }
            numEEs++;
        }

        // print the histogram.
        int[] counts = new int[MAX_EXPS + 1];
        for ( Long outer : mat.getRowNames() ) {
            double[] row = mat.getRowByName( outer );
            for ( double d : row ) {
                counts[( int ) d]++;
            }
        }

        for ( int j = 0; j < counts.length; j++ ) {
            System.out.println( j + "\t" + counts[j] );
        }
    }

    /**
     * For each composites sequence, count how many genes there are. This does not take into account multiple occurrence
     * os the same sequence in different probes!
     *
     * @param taxon
     */
    public void genesPerProbe( Taxon taxon ) {
        Collection<ArrayDesign> allAds = adService.loadAll();
        Collection<ArrayDesign> ads = new HashSet<>();
        for ( ArrayDesign ad : allAds ) {
            Taxon t = ad.getPrimaryTaxon();
            if ( t != null && t.equals( taxon ) ) {
                ads.add( ad );
            }
        }

        Map<ArrayDesign, Map<Integer, Integer>> countMap = new HashMap<>();
        Collection<Long> seenSeqs = new HashSet<>();
        int count = 0;
        for ( ArrayDesign design : ads ) {
            log.info( design + " : " + ++count + " of " + ads.size() );
            adService.thawLite( design );
            Map<Integer, Integer> counts = new HashMap<>();

            int i = 0;
            for ( CompositeSequence cs : design.getCompositeSequences() ) {

                BioSequence bs = cs.getBiologicalCharacteristic();

                if ( bs == null ) continue; // these don't count.

                if ( seenSeqs.contains( bs ) ) continue;

                Integer numGenes = 0;
                Collection<Gene> genes = compositeSequenceService.getGenes( cs );

                if ( genes == null )
                    numGenes = 0;
                else
                    numGenes = genes.size();

                if ( !counts.containsKey( numGenes ) ) {
                    counts.put( numGenes, 1 );
                } else {
                    counts.put( numGenes, counts.get( numGenes ) + 1 );
                }
                if ( ++i % 1000 == 0 ) {
                    log.info( "Processed " + i + " compositeSequences" );
                }

                seenSeqs.add( bs.getId() );
            }
            countMap.put( design, counts );

        }

        try {
            printGenesPerProbeCountMap( countMap );
        } catch ( IOException e ) {
            e.printStackTrace();
        }
    }

    /**
     * For each gene, count how many microarray probes there are.
     *
     * @param taxon
     */
    private void probesPerGene( Taxon taxon ) {
        Collection<Gene> genes = geneService.loadAll( taxon );
        Map<Long, Integer> counts = new HashMap<>();
        int i = 0;
        for ( Gene gene : genes ) {
            Collection<CompositeSequence> compositeSequences = geneService.getCompositeSequencesById( gene.getId() );
            counts.put( gene.getId(), compositeSequences.size() );
            if ( ++i % 1000 == 0 ) {
                log.info( "Processed " + i + " genes" );
            }
        }
        for ( Long l : counts.keySet() ) {
            System.out.println( l + "\t" + counts.get( l ) );
        }
    }

    private void printGenesPerProbeCountMap( Map<ArrayDesign, Map<Integer, Integer>> countMap ) throws IOException {
        PrintWriter out;
        if ( outFileName == null ) {
            out = new PrintWriter( System.out );
        } else {
            out = new PrintWriter( new FileWriter( outFileName ) );
        }

        // get max number of genes
        int maxNumGenes = 0;
        for ( Map<Integer, Integer> counts : countMap.values() ) {
            for ( Integer n : counts.keySet() ) {
                int count = counts.get( n );
                if ( count > 0 && n > maxNumGenes ) maxNumGenes = n;
            }
        }

        StringBuffer buf = new StringBuffer( "Count\t" );
        for ( ArrayDesign ad : countMap.keySet() ) {
            buf.append( ad.getShortName() + "\t" );
        }
        buf.deleteCharAt( buf.length() - 1 );
        out.println( buf );

        for ( int numGenes = 0; numGenes <= maxNumGenes; numGenes++ ) {
            buf = new StringBuffer();
            buf.append( numGenes + "\t" );
            for ( ArrayDesign ad : countMap.keySet() ) {
                Map<Integer, Integer> counts = countMap.get( ad );
                if ( counts.get( numGenes ) != null )
                    buf.append( counts.get( numGenes ) );
                else
                    buf.append( "0" );
                buf.append( "\t" );
            }
            buf.deleteCharAt( buf.length() - 1 );
            out.println( buf );
        }
        out.close();
    }
}
