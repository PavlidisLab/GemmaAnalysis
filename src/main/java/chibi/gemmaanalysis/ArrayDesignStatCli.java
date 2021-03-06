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
import java.util.Map;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.lang.time.StopWatch;

import cern.colt.Arrays;
import ubic.gemma.core.apps.ArrayDesignSequenceManipulatingCli;
import ubic.gemma.core.genome.gene.service.GeneService;
import ubic.gemma.model.common.auditAndSecurity.curation.CurationDetails;
import ubic.gemma.model.expression.arrayDesign.ArrayDesign;
import ubic.gemma.model.expression.designElement.CompositeSequence;
import ubic.gemma.model.genome.Gene;
import ubic.gemma.model.genome.Taxon;
import ubic.gemma.persistence.service.expression.arrayDesign.ArrayDesignService;
import ubic.gemma.persistence.service.expression.designElement.CompositeSequenceService;
import ubic.gemma.persistence.service.genome.taxon.TaxonService;

/**
 * CLI for ArrayDesignMapSummaryService
 *
 * @author xwan
 */
public class ArrayDesignStatCli extends ArrayDesignSequenceManipulatingCli {

    private final static int MAXIMUM_COUNT = 10;
    private final static String DEFAULT_OUT_FILE = "arraydesignsummary.txt";
    //private StatusService statusService;
    private final static String NA = "NA";

    public static void main( String[] args ) {
        ArrayDesignStatCli s = new ArrayDesignStatCli();
        try {
            Exception ex = s.doWork( args );
            if ( ex != null ) {
                ex.printStackTrace();
            }
        } catch ( Exception e ) {
            throw new RuntimeException( e );
        }
    }

    private ArrayDesignService adService;
    private GeneService geneService;
    private Collection<Long> geneIds = new HashSet<>();
    private CompositeSequenceService compositeSequenceService;
    private TaxonService taxonService;
    private String outFile;

    private Taxon taxon;

    /*
     * (non-Javadoc)
     *
     * @see ubic.gemma.util.AbstractCLI#getCommandName()
     */
    @Override
    public String getCommandName() {
        return "platformStats";
    }

    Map<Long, Collection<Long>> getGeneId2CSIdsMap( Map<Long, Collection<Long>> csId2geneIds ) {
        Map<Long, Collection<Long>> geneId2csIds = new HashMap<>();
        for ( Long csId : csId2geneIds.keySet() ) {
            Collection<Long> gids = csId2geneIds.get( csId );
            for ( Long geneId : gids ) {
                Collection<Long> csIds = geneId2csIds.get( geneId );
                if ( csIds == null ) {
                    csIds = new HashSet<>();
                    geneId2csIds.put( geneId, csIds );
                }
                csIds.add( csId );
            }
        }
        return geneId2csIds;
    }

    int[] getStats( Map<Long, Collection<Long>> dataMap, boolean geneIdKey ) {
        int[] res = new int[MAXIMUM_COUNT];
        for ( Long id : dataMap.keySet() ) {
            if ( geneIdKey ) {
                if ( !geneIds.contains( id ) ) continue;
            }
            int size = dataMap.get( id ).size();
            if ( !geneIdKey ) {
                size = 0;
                Collection<Long> ids = dataMap.get( id );
                for ( Long geneId : ids ) {
                    if ( geneIds.contains( geneId ) ) size++;
                }
            }
            if ( size == 0 ) continue;
            if ( size > MAXIMUM_COUNT ) size = MAXIMUM_COUNT;
            res[size - 1]++;
        }
        return res;
    }

    @Override
    @SuppressWarnings("static-access")
    protected void buildOptions() {
        super.buildOptions();

        OptionBuilder.hasArg();
        OptionBuilder.withArgName( "outfile" );
        OptionBuilder.withDescription( "TSV output filename" );
        OptionBuilder
                .withLongOpt( "outfile" );
        Option expOption = OptionBuilder.create( 'o' );
        addOption( expOption );

        OptionBuilder.hasArg();
        OptionBuilder.withDescription( "taxon name" );
        OptionBuilder
                .withDescription( "Taxon of the expression experiments and genes" );
        OptionBuilder.withLongOpt( "taxon" );
        Option taxonOption = OptionBuilder.create( 't' );
        addOption( taxonOption );
    }

    /*
     * (non-Javadoc)
     *
     * @see ubic.gemma.util.AbstractCLI#doWork(java.lang.String[])
     */
    @Override
    protected Exception doWork( String[] args ) {
        Collection<String> failedAds = new ArrayList<>();
        Exception err = processCommandLine( args );
        if ( err != null ) return err;

        Map<Taxon, Collection<ArrayDesign>> taxon2arraydesign = new HashMap<>();
        Collection<Long> adIds = new HashSet<>();
        for ( ArrayDesign ad : this.getArrayDesignsToProcess() ) {

            Taxon taxon = ad.getPrimaryTaxon();
            if ( taxon == null ) {
                System.err.println( "ArrayDesign " + ad.getName() + " doesn't have a taxon" );
                continue;
            }
            taxon = taxonService.load( taxon.getId() );

            if ( taxon != null && taxon.getCommonName() == null ) {
                log.warn( ad.getShortName() + " taxon common name is null" );
                failedAds.add( ad.getShortName() );
                continue;
            }

            // filter out taxon
            if ( this.taxon != null && this.taxon.getCommonName() != null
                    && !taxon.getCommonName().equalsIgnoreCase( this.taxon.getCommonName() ) ) {
                continue;
            }

            adIds.add( ad.getId() );

            Collection<ArrayDesign> ads = null;
            ads = taxon2arraydesign.get( taxon );
            if ( ads == null ) {
                ads = new HashSet<>();
                taxon2arraydesign.put( taxon, ads );
            }
            ads.add( ad );
        }
        Map<Long, Boolean> isMerged = adService.isMerged( adIds );
        Map<Long, Boolean> isSubsumed = adService.isSubsumed( adIds );
        StopWatch timer = new StopWatch();
        timer.start();
        int lineCount = 0;
        try (FileWriter out = new FileWriter( new File( this.outFile ) );) {
            String header = "taxon\tshortName\tname\tisTroubled\texperiments\tmergees\tsubsumes\tsubsumedBy\tgenes\tprobes\tcsWithGenes\tcsBioSequences\tcsBlatResults\tP2G_0";
            out.write( header );
            for ( int i = 1; i <= MAXIMUM_COUNT; i++ )
                out.write( "\tP2G_" + i );
            for ( int i = 1; i <= MAXIMUM_COUNT; i++ )
                out.write( "\tG2P_" + i );
            out.write( "\n" );
            System.err.print( header + "\n" );
            for ( Taxon t : taxon2arraydesign.keySet() ) {
                Collection<ArrayDesign> ads = taxon2arraydesign.get( t );
                Collection<Gene> allGenes = geneService.getGenesByTaxon( t );
                for ( Gene gene : allGenes ) {
                    geneIds.add( gene.getId() );

                }
                for ( ArrayDesign ad : ads ) {

                    try {
                        CurationDetails status = ad.getCurationDetails();
                        String isTroubled = status != null ? Boolean.toString( status.getTroubled().booleanValue() )
                                : NA;
                        getArrayDesignService().thawLite( ad );
                        long mergees = ad.getMergees().size();
                        long subsumes = ad.getSubsumedArrayDesigns().size();
                        String subsumedBy = ad.getSubsumingArrayDesign() != null ? ad.getSubsumingArrayDesign()
                                .getShortName() : NA;
                        long numEEs = this.getArrayDesignService().getExpressionExperiments( ad ).size();
                        // boolean merged = isMerged.get( ad.getId() );
                        // if ( merged ) continue;
                        // boolean subsumed = isSubsumed.get( ad.getId() );
                        // if ( subsumed ) continue;
                        long numProbes = getArrayDesignService().getCompositeSequenceCount( ad ).longValue();
                        long numCsBioSequences = getArrayDesignService().numCompositeSequenceWithBioSequences( ad );
                        long numCsBlatResults = getArrayDesignService().numCompositeSequenceWithBlatResults( ad );
                        long numCsGenes = getArrayDesignService().numCompositeSequenceWithGenes( ad );
                        long numGenes = getArrayDesignService().numGenes( ad );
                        Collection<CompositeSequence> allCSs = getArrayDesignService().getCompositeSequences( ad );
                        Collection<Long> csIds = new HashSet<>();
                        for ( CompositeSequence cs : allCSs )
                            csIds.add( cs.getId() );
                        // FIXME this used to provide only known genes.
                        Map<Long, Collection<Long>> csId2geneIds = this.getCs2GeneMap( csIds );
                        Map<Long, Collection<Long>> geneId2csIds = getGeneId2CSIdsMap( csId2geneIds );
                        int[] csStats = getStats( csId2geneIds, false );
                        int[] geneStats = getStats( geneId2csIds, true );
                        int cs2NoneGene = allCSs.size() - csId2geneIds.keySet().size();
                        String line = taxon.getCommonName() + "\t" + ad.getShortName() + "\t" + ad.getName() + "\t"
                                + isTroubled + "\t" + numEEs + "\t" + mergees + "\t" + subsumes + "\t" + subsumedBy
                                + "\t" + numGenes + "\t" + numProbes + "\t" + numCsGenes + "\t" + numCsBioSequences
                                + "\t" + numCsBlatResults + "\t" + cs2NoneGene;
                        out.write( line );
                        for ( int i = 0; i < MAXIMUM_COUNT; i++ ) {
                            out.write( "\t" + csStats[i] );
                        }
                        for ( int i = 0; i < MAXIMUM_COUNT; i++ ) {
                            out.write( "\t" + geneStats[i] );
                        }
                        out.write( "\n" );
                        System.err.print( line + "\n" );

                        lineCount++;
                    } catch ( Exception e ) {
                        log.error( e, e );
                        failedAds.add( ad.getShortName() );
                        continue;
                    }
                }
            }
            out.close();
            log.info( "Skipped " + failedAds.size() + " array designs : " + Arrays.toString( failedAds.toArray() ) );
            log.info( "Finished running in " + timer.getTime() + " ms." );
            log.info( "Wrote " + lineCount + " lines to " + outFile );

        } catch ( Exception e ) {
            return e;
        }
        return null;
    }

    @Override
    protected void processOptions() {
        super.processOptions();
        // FIXME: add HTML output option.

        this.adService = this.getBean( ArrayDesignService.class );
        this.compositeSequenceService = this.getBean( CompositeSequenceService.class );
        this.geneService = this.getBean( GeneService.class );
        this.taxonService = this.getBean( TaxonService.class );

        if ( hasOption( 'o' ) ) {
            this.outFile = getOptionValue( 'o' );
            log.info( "Output will be written to " + outFile );
        } else {
            this.outFile = DEFAULT_OUT_FILE;
        }

        if ( hasOption( 't' ) ) {
            String taxonName = getOptionValue( 't' );
            this.taxon = taxonService.findByCommonName( taxonName );
            if ( this.taxon == null ) {
                log.error( "ERROR: Cannot find taxon " + taxonName );
            }
        }

    }

    private Map<Long, Collection<Long>> getCs2GeneMap( Collection<Long> csIds ) {
        Map<CompositeSequence, Collection<Gene>> genes = compositeSequenceService.getGenes( compositeSequenceService
                .load( csIds ) );
        Map<Long, Collection<Long>> result = new HashMap<>();
        for ( CompositeSequence cs : genes.keySet() ) {
            result.put( cs.getId(), new HashSet<Long>() );
            for ( Gene g : genes.get( cs ) ) {
                result.get( cs.getId() ).add( g.getId() );
            }
        }
        return result;
    }

}
