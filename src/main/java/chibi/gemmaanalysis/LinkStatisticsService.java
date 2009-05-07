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

import java.io.IOException;
import java.io.Writer;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.apache.commons.lang.time.StopWatch;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.gemma.model.analysis.expression.coexpression.ProbeCoexpressionAnalysisService;
import ubic.gemma.model.association.coexpression.Probe2ProbeCoexpressionService;
import ubic.gemma.model.association.coexpression.Probe2ProbeCoexpressionDaoImpl.ProbeLink;
import ubic.gemma.model.expression.designElement.CompositeSequence;
import ubic.gemma.model.expression.designElement.CompositeSequenceService;
import ubic.gemma.model.expression.experiment.BioAssaySet;
import ubic.gemma.model.expression.experiment.ExpressionExperiment;
import ubic.gemma.model.genome.Gene;
import ubic.gemma.model.genome.Taxon;

/**
 * Methods for analyzing links from the database.
 * 
 * @author paul
 * @author xwan
 * @version $Id$
 */
public class LinkStatisticsService {

    private static Log log = LogFactory.getLog( LinkStatisticsService.class.getName() );

    private Probe2ProbeCoexpressionService p2pService = null;

    private ProbeCoexpressionAnalysisService probeCoexpressionAnalysisService;

    private CompositeSequenceService compositeSequenceService;

    /**
     * Tally and summarize the coexpression links stored for the given data sets among the given genes.
     * 
     * @param ees ExpressionExperiments to use
     * @param genes Genes to consider
     * @param taxon for the experiments (convenient, this shouldn't really be needed)
     * @return LinkStatistic object holding the results.
     */
    public LinkStatistics analyze( Collection<BioAssaySet> ees, Collection<Gene> genes, Taxon taxon,
            boolean shuffleLinks, boolean filterNonSpecific ) {
        LinkStatistics stats = new LinkStatistics( ees, genes );
        int numLinks = 0;
        for ( BioAssaySet bas : ees ) {
            ExpressionExperiment ee = ( ExpressionExperiment ) bas;
            numLinks += tallyLinks( stats, ee, taxon, shuffleLinks, filterNonSpecific );
        }
        log.info( numLinks + " gene links in total for " + ees.size() + " expression experiments " );
        return stats;
    }

    /**
     * Tabulate links in an ExpressionExperiment, adding them to the tally.
     * 
     * @param stats object to hold the results
     * @param ee ExpressionExperiment to analyze.
     * @param taxon (should be able to figure this out from the ee, but it's easier)
     * @param shuffleLinks
     * @param filterNonSpecific
     * @return the number of links added
     */
    @SuppressWarnings("unchecked")
    private int tallyLinks( LinkStatistics stats, ExpressionExperiment ee, Taxon taxon, boolean shuffleLinks,
            boolean filterNonSpecific ) {
        assert ee != null;

        // FIXME if not shuffling, don't use the working table, so we can 'get on with it' without worrying about
        // creating that table first.
        log.info( "Loading links for  " + ee );

        Map<Long, Collection<Long>> assayedProbe2Genes = getRelevantProbe2GeneIds( ee, stats, filterNonSpecific );

        Collection<ProbeLink> links = p2pService.getProbeCoExpression( ee, taxon.getCommonName(), filterNonSpecific );
        Collection<ProbeLink> filteredLinks = filterLinks( stats, links, assayedProbe2Genes, filterNonSpecific );
        if ( filteredLinks == null || filteredLinks.size() == 0 ) return 0;

        Collection<GeneLink> geneLinks = getGeneLinks( filteredLinks, stats, assayedProbe2Genes );

        if ( shuffleLinks ) {
            log.info( "Shuffling links for " + ee );
            /*
             * We start with a list of all the probes for genes that are in our matrix and which are expressed.
             */

            /* for shuffling at the probe level */

            shuffleProbeLinks( stats, ee, filteredLinks, assayedProbe2Genes.keySet() );
            geneLinks = getGeneLinks( filteredLinks, stats, assayedProbe2Genes );

            /* shuffle at gene level */
            // geneLinks = shuffleGeneLinks( geneLinks, assayedProbe2Genes );
        }

        return stats.addLinks( geneLinks, ee );
    }

    /**
     * Limit the links to those for genes that were included in the initial setup. Typically this means that links
     * between PARs or PARs and Genes can be removed. This step is important when shuffling, so that the probes used for
     * shuffling match those used in the real links.
     * 
     * @param stats
     * @param links
     * @param cs2genes
     * @param filterNonSpecific
     * @return filtered collection of links.
     */
    private Collection<ProbeLink> filterLinks( LinkStatistics stats, Collection<ProbeLink> links,
            Map<Long, Collection<Long>> cs2genes, boolean filterNonSpecific ) {

        if ( links.size() == 0 ) return links;

        Collection<Long> probeIdsToKeepLinksFor = filterProbes( stats, cs2genes, filterNonSpecific );

        Collection<ProbeLink> filteredLinks = new HashSet<ProbeLink>();
        for ( ProbeLink pl : links ) {
            if ( probeIdsToKeepLinksFor.contains( pl.getFirstDesignElementId() )
                    && probeIdsToKeepLinksFor.contains( pl.getSecondDesignElementId() ) ) {
                filteredLinks.add( pl );
            }
        }

        log.info( filteredLinks.size() + "/" + links.size()
                + " links retained after removing links for excluded genes." );
        return filteredLinks;
    }

    /**
     * Filter probes to include just those used in the matrix.
     * 
     * @param stats
     * @param cs2genes map of probes -> genes for probes you want to filter.
     * @param filterNonSpecific.
     * @return Collection of probes
     */
    private Collection<Long> filterProbes( LinkStatistics stats, Map<Long, Collection<Long>> cs2genes,
            boolean filterNonSpecific ) {
        Collection<Long> matrixGenes = stats.getGeneIds();

        // when probe is for multiple genes, we keep even if only one of the genes is on our query list - unless
        // filternonspecific.
        Collection<Long> probeIdsToKeepLinksFor = new HashSet<Long>();
        for ( Long probeId : cs2genes.keySet() ) {
            Collection<Long> genesForProbesInLinks = cs2genes.get( probeId );

            /*
             * Important if we dont' filter out non-specific probes here, we should not filter them out when we first
             * fetch links.
             */
            if ( filterNonSpecific && genesForProbesInLinks.size() > 1 ) continue;

            for ( Long geneForLink : genesForProbesInLinks ) {
                if ( matrixGenes.contains( geneForLink ) ) {
                    probeIdsToKeepLinksFor.add( probeId );
                    break; // we know we have to keep this probe.
                }
            }
        }
        log.info( "Filtered irrelevant probes: kept " + probeIdsToKeepLinksFor.size() + "/" + cs2genes.keySet().size()
                + " probes" );
        return probeIdsToKeepLinksFor;
    }

    /**
     * @param csIds
     * @return
     */
    private Map<Long, Collection<Long>> getCs2GeneMapFromProbeIds( Collection<Long> csIds ) {
        Map<Long, Collection<Long>> result = new HashMap<Long, Collection<Long>>();

        if ( csIds.size() == 0 ) return result;

        // note that this includes PARs etc.
        Map<CompositeSequence, Collection<Gene>> genes = compositeSequenceService.getGenes( compositeSequenceService
                .loadMultiple( csIds ) );

        for ( CompositeSequence cs : genes.keySet() ) {
            result.put( cs.getId(), new HashSet<Long>() );
            for ( Gene g : genes.get( cs ) ) {
                result.get( cs.getId() ).add( g.getId() );
            }
        }
        return result;
    }

    /**
     * Given probe links, convert to gene links.
     * 
     * @param links
     * @return collection of GeneLink objects.
     */
    private Collection<GeneLink> getGeneLinks( Collection<ProbeLink> links, LinkStatistics stats,
            Map<Long, Collection<Long>> cs2genes ) {
        log.info( "Converting " + links.size() + " probe links to gene links ..." );
        Collection<GeneLink> result = new HashSet<GeneLink>();

        /*
         * Once again, ignore links that don't show up in the original query
         */
        Collection<Long> geneIds = stats.getGeneIds();

        int missing = 0;
        int count = 0;
        for ( ProbeLink link : links ) {
            Collection<Long> firstGeneIds = cs2genes.get( link.getFirstDesignElementId() );
            Collection<Long> secondGeneIds = cs2genes.get( link.getSecondDesignElementId() );
            if ( firstGeneIds == null || secondGeneIds == null ) {
                ++missing;
                log.debug( "No gene found for one/both of CS:  " + link.getFirstDesignElementId() + ","
                        + link.getSecondDesignElementId() );
                continue;
            }

            for ( Long firstGeneId : firstGeneIds ) {
                if ( !geneIds.contains( firstGeneId ) ) continue;
                for ( Long secondGeneId : secondGeneIds ) {
                    if ( !geneIds.contains( secondGeneId ) ) continue;
                    result.add( new GeneLink( firstGeneId, secondGeneId, link.getScore() ) );
                }
            }
            if ( ++count % 5e5 == 0 ) {
                log.info( count + " links converted" );
            }
        }

        if ( missing > 0 ) {
            log.warn( missing + " links had probes with no genes" );
        }

        log.info( links.size() + " probe links --> " + result.size() + " gene links" );

        return result;
    }

    /**
     * Used for shuffling, to choose the universe of probes to consider 'random links' among. This must be chosen
     * carefully to get an accurate null distribution.
     * 
     * @param ee
     * @param stats
     * @param filterNonSpecific
     * @return map of probes to genes which are 1) assayed and 2) have gene mappings (that is, alignments to the genome;
     *         this includes non-refseq genes etc.) and 3) used as inputs for the 'real' analysis we are comparing to.
     */
    private Map<Long, Collection<Long>> getRelevantProbe2GeneIds( ExpressionExperiment ee, LinkStatistics stats,
            boolean filterNonSpecific ) {

        Collection<CompositeSequence> probesAssayed = probeCoexpressionAnalysisService.getAssayedProbes( ee );

        List<Long> assayedProbeIds = new ArrayList<Long>();
        for ( CompositeSequence cs : probesAssayed ) {
            assayedProbeIds.add( cs.getId() );
        }
        log.info( assayedProbeIds.size() + " probes assayed " );

        /*
         * Further filter the list to only include probes that have alignments. This is a map of CS to genes (this is
         * probably not strictly needed as the next step also looks at genes.)
         */
        Map<Long, Collection<Long>> geneMap = getCs2GeneMapFromProbeIds( assayedProbeIds );
        log.info( geneMap.size() + " probes with alignments" );

        /*
         * Further filter to include only probes that are also for genes we used as potential inputs.
         */
        Collection<Long> probesToKeep = filterProbes( stats, geneMap, filterNonSpecific );

        log.info( probesToKeep.size() + " probes after filtering." );

        // return map so we can still have access to the genes. Remove all 'extraneous' genes.
        Collection<Long> geneIdsToKeep = stats.getGeneIds();
        Map<Long, Collection<Long>> finalGeneMap = new HashMap<Long, Collection<Long>>();
        for ( Long p : probesToKeep ) {
            Collection<Long> geneIdsToClean = geneMap.get( p );
            geneIdsToClean.retainAll( geneIdsToKeep );
            finalGeneMap.put( p, geneIdsToClean );
        }

        // return new ArrayList<Long>( probesToKeep );
        return finalGeneMap;
    }

    /**
     * @param ees
     * @param taxonName
     * @param filterNonSpecific
     */
    public void prepareDatabase( Collection<BioAssaySet> ees, String taxonName, boolean filterNonSpecific ) {
        log.info( "Creating working table for link analysis" );
        StopWatch watch = new StopWatch();
        watch.start();
        p2pService.prepareForShuffling( ees, taxonName, filterNonSpecific );
        watch.stop();
        log.info( "Done, spent " + watch.getTime() / 1000 + "s preparing the database" );
    }

    /**
     * @param compositeSequenceService the compositeSequenceService to set
     */
    public void setCompositeSequenceService( CompositeSequenceService compositeSequenceService ) {
        this.compositeSequenceService = compositeSequenceService;
    }

    public void setP2pService( Probe2ProbeCoexpressionService service ) {
        p2pService = service;
    }

    public void setProbeCoexpressionAnalysisService( ProbeCoexpressionAnalysisService probeCoexpressionAnalysisService ) {
        this.probeCoexpressionAnalysisService = probeCoexpressionAnalysisService;
    }

    /**
     * @param geneLinks unshuffled links.
     * @param probeGeneMap
     * @return
     */
    protected Collection<GeneLink> shuffleGeneLinks( Collection<GeneLink> geneLinks,
            Map<Long, Collection<Long>> probeGeneMap ) {
        Collection<Long> geneIds = new HashSet<Long>();
        for ( Long probeId : probeGeneMap.keySet() ) {
            geneIds.addAll( probeGeneMap.get( probeId ) );
        }

        List<Long> geneIdList = new ArrayList<Long>( geneIds );
        List<Long> shuffledGeneIdList = new ArrayList<Long>( geneIdList );

        Collections.shuffle( shuffledGeneIdList );

        log.debug( geneIdList.size() + " genes to shuffle" );
        // assign each gene ID to a random one
        Map<Long, Long> shuffleMap = new HashMap<Long, Long>();
        for ( int i = 0; i < geneIdList.size(); i++ ) {
            shuffleMap.put( geneIdList.get( i ), shuffledGeneIdList.get( i ) );
        }

        // keep track of genes in ee
        Set<Long> usedGenes = new HashSet<Long>();
        for ( GeneLink gl : geneLinks ) {
            if ( !shuffleMap.containsKey( gl.getFirstGene() ) ) {
                log.error( "Shuffled Gene->Probe map did not contain gene used in links: " + gl.getFirstGene() );
                continue;
            }
            if ( !shuffleMap.containsKey( gl.getSecondGene() ) ) {
                log.error( "Shuffled Gene->Probe map did not contain gene used in links: " + gl.getSecondGene() );
                continue;
            }

            usedGenes.add( gl.getFirstGene() );
            usedGenes.add( gl.getSecondGene() );

            // reassign gene links
            gl.setFirstGene( shuffleMap.get( gl.getFirstGene() ) );
            gl.setSecondGene( shuffleMap.get( gl.getSecondGene() ) );
        }

        log.info( "Gene links used a total of " + usedGenes.size() + " genes, we shuffled using " + geneIdList.size()
                + " available genes" );

        log.debug( geneLinks.size() + " links after shuffling" );
        return geneLinks;
    }

    /**
     * Given links, shuffle the probes they are between. The way this is done is designed to reflect the null
     * distribution of link selection in the data set, maintaining the same number of links per probe.
     * 
     * @param stats
     * @param ee
     * @param probeUniverse set probe ids, including only probes to consider.
     * @param linksToShuffle
     * @return
     */
    private void shuffleProbeLinks( LinkStatistics stats, ExpressionExperiment ee,
            Collection<ProbeLink> linksToShuffle, Collection<Long> probeUniverse ) {
        log.info( "Shuffling links for " + ee.getShortName() );

        /*
         * Make a copy so we can make a shuffled mapping
         */
        List<Long> assayedProbeIdsforShuffle = new ArrayList<Long>();
        List<Long> assayedProbeIds = new ArrayList<Long>( probeUniverse );
        for ( Long id : assayedProbeIds ) {
            assayedProbeIdsforShuffle.add( id );
        }

        /*
         * Shuffle the probes and create map. We exchange old probes for random new ones.
         */
        Collections.shuffle( assayedProbeIdsforShuffle );
        // shuffleList( assayedProbeIdsforShuffle );
        Map<Long, Long> shuffleMap = new HashMap<Long, Long>();
        for ( int i = 0; i < assayedProbeIdsforShuffle.size(); i++ ) {
            Long old = assayedProbeIds.get( i );
            Long shuffled = assayedProbeIdsforShuffle.get( i );
            shuffleMap.put( old, shuffled );
        }

        // Now replace the in the probes in the links with the shuffled ones, from the (potentially) expanded list.. Now
        // the links are among random probes, but number per probe is the same.
        for ( ProbeLink p : linksToShuffle ) {

            // Sanity checks.
            if ( !shuffleMap.containsKey( p.getFirstDesignElementId() ) ) {
                throw new IllegalStateException( "Probe " + p.getFirstDesignElementId()
                        + " was used in the links but doesn't show up in the 'assayedProbes'" );
            }
            if ( !shuffleMap.containsKey( p.getSecondDesignElementId() ) ) {
                throw new IllegalStateException( "Probe " + p.getSecondDesignElementId()
                        + " was used in the links but doesn't show up in the 'assayedProbes'" );
            }

            // Replace with the shuffled replacement.
            p.setFirstDesignElementId( shuffleMap.get( p.getFirstDesignElementId() ) );
            p.setSecondDesignElementId( shuffleMap.get( p.getSecondDesignElementId() ) );
        }

        // this step might be redundant.
        shuffleProbeLinks( new ArrayList<ProbeLink>( linksToShuffle ) );

    }

    /**
     * Do shuffling at probe level
     * 
     * @param links
     */
    private void shuffleProbeLinks( List<ProbeLink> links ) {
        Random random = new Random();
        int i = links.size();
        while ( --i > 0 ) {
            int k = random.nextInt( i + 1 );
            Long tmpId = links.get( i ).getSecondDesignElementId();
            links.get( i ).setSecondDesignElementId( links.get( k ).getSecondDesignElementId() );
            links.get( k ).setSecondDesignElementId( tmpId );
        }
    }

    /**
     * @param out
     * @param realStats
     * @param shuffleStats
     * @param maxSupport
     * @throws IOException
     */
    private void writeShuffleStats( Writer out, LinkConfirmationStatistics realStats,
            Collection<LinkConfirmationStatistics> shuffleStats, int maxSupport ) throws IOException {
        NumberFormat nf = NumberFormat.getInstance();
        nf.setMaximumFractionDigits( 3 );

        if ( realStats != null ) {
            out.write( "ShuffleMean" );
            double[] falsePositiveRates = new double[maxSupport + 1];
            for ( LinkConfirmationStatistics shuffleStat : shuffleStats ) {

                for ( int j = 1; j <= maxSupport; j++ ) {
                    if ( realStats.getRepCount( j ) != 0 ) {
                        falsePositiveRates[j] = falsePositiveRates[j] + ( double ) shuffleStat.getRepCount( j )
                                / ( double ) realStats.getRepCount( j );
                    }
                }
            }
            for ( int j = 1; j <= maxSupport; j++ ) {
                out.write( "\t" + nf.format( falsePositiveRates[j] / shuffleStats.size() ) );
            }
            out.write( "\n" );
        }

        // FP rates for each run.
        int i = 1;
        for ( LinkConfirmationStatistics shuffleStat : shuffleStats ) {
            out.write( "ShuffleRun_" + i );
            for ( int j = 1; j <= maxSupport; j++ ) {
                out.write( "\t" + shuffleStat.getRepCount( j ) );
                // if ( realStats.getCumulativeRepCount( j ) != 0 ) {
                // out.write( nf.format( ( double ) shuffleStat.getCumulativeRepCount( j )
                // / ( double ) realStats.getCumulativeRepCount( j ) ) );
                // }
            }
            ++i;
            out.write( "\n" );
        }
    }

    /**
     * @param out
     * @param realStats
     * @param shuffleStats
     */
    public void writeStats( Writer out, LinkConfirmationStatistics realStats,
            Collection<LinkConfirmationStatistics> shuffleStats ) {

        // Determine the largest number of data sets any link was seen in.
        int maxSupport = -1;
        if ( realStats != null ) {
            maxSupport = realStats.getMaxLinkSupport();
        } else {
            for ( LinkConfirmationStatistics statistics : shuffleStats ) {
                int s = statistics.getMaxLinkSupport();
                if ( s > maxSupport ) maxSupport = s;
            }
        }

        try {
            NumberFormat nf = NumberFormat.getInstance();
            nf.setMaximumFractionDigits( 3 );

            // Header
            out.write( "Support" );
            for ( int j = 1; j <= maxSupport; j++ )
                out.write( "\t" + j /* + "+" */);
            out.write( "\n" );

            // True links
            if ( realStats != null ) {
                out.write( "RealLinks" );
                for ( int j = 1; j <= maxSupport; j++ ) {
                    out.write( "\t" + realStats.getRepCount( j ) );
                }
                out.write( "\n" );
            }

            if ( shuffleStats != null && shuffleStats.size() > 0 ) {
                // number of shuffled links / # of real links
                writeShuffleStats( out, realStats, shuffleStats, maxSupport );
            }
            out.close();
        } catch ( Exception e ) {
            throw new RuntimeException( e );
        }
    }
}
