package chibi.gemmaanalysis;

import java.util.Collection;

import ubic.gemma.model.genome.Gene;
import ubic.gemma.model.genome.PhysicalLocation;
import ubic.gemma.model.genome.PhysicalLocationService;
import ubic.gemma.model.genome.ProbeAlignedRegion;
import ubic.gemma.model.genome.RelativeLocationData;
import ubic.gemma.model.genome.Taxon;
import ubic.gemma.model.genome.TaxonService;
import ubic.gemma.model.genome.gene.GeneService;
import ubic.gemma.util.AbstractSpringAwareCLI;

/**
 * Collect statistics about how far PARs are from known genes.
 * 
 * @author paul
 * @version $Id$
 */
public class PARMapper extends AbstractSpringAwareCLI {

    GeneService parService;
    TaxonService taxonService;
    PhysicalLocationService physicalLocationservice;

    @Override
    protected void buildOptions() {
    }

    @Override
    protected Exception doWork( String[] args ) {
        Exception exc = processCommandLine( "PARMapper", args );
        if ( exc != null ) return exc;

        this.taxonService = ( TaxonService ) this.getBean( "taxonService" );
        this.parService = ( GeneService ) this.getBean( "geneService" );
        this.physicalLocationservice = ( PhysicalLocationService ) this.getBean( "physicalLocationService" );

        /*
         * FIXME don't hard-code this.
         */
        Taxon taxon = taxonService.findByCommonName( "mouse" );
        if ( taxon == null ) {
            throw new IllegalArgumentException();
        }
        Collection<ProbeAlignedRegion> pars = parService.loadProbeAlignedRegions( taxon );
        log.info( pars.size() + " " + taxon.getCommonName() + " PARS" );

        System.out.println( "ParID\tParName\tChrom\tNuc\tGeneId\tGeneSymbol\tDistance\tGeneContainsPar\tSameStrand" );

        // test case
        // ProbeAlignedRegion par = ( ProbeAlignedRegion ) parService.load( 1450985 );

        for ( ProbeAlignedRegion par : pars ) {
            parService.thaw( par );

            PhysicalLocation loc = par.getPhysicalLocation();

            physicalLocationservice.thaw( loc );

            RelativeLocationData nearest = parService.findNearest( loc );
            if ( nearest != null ) {
                Gene gene = nearest.getNearestGene();

                System.out
                        .println( par.getId() + "\t" + par.getName() + "\t" + loc.getChromosome().getName() + "\t"
                                + loc.getNucleotide() + "\t" + gene.getId() + "\t" + gene.getOfficialSymbol() + "\t"
                                + nearest.getRange() + "\t" + nearest.isContainedWithinGene() + "\t"
                                + nearest.isOnSameStrand() );
            }
        }
        return null;
    }

    /**
     * @param args
     */
    public static void main( String[] args ) {
        PARMapper p = new PARMapper();
        Exception e = p.doWork( args );
        if ( e != null ) {
            throw new RuntimeException( e );
        }

    }

}
