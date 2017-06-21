/*
 * The GemmaAnalysis project
 *
 * Copyright (c) 2010 University of British Columbia
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
package chibi.gemmaanalysis.cli.deprecated;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;

import ubic.gemma.core.genome.gene.service.GeneService;
import ubic.gemma.core.util.AbstractSpringAwareCLI;
import ubic.gemma.model.genome.Taxon;
import ubic.gemma.persistence.service.genome.PhysicalLocationService;
import ubic.gemma.persistence.service.genome.taxon.TaxonService;

/**
 * Collect statistics about how far PARs are from known genes.
 *
 * @author paul
 * @version $Id: PARMapper.java,v 1.6 2015/11/12 19:37:13 paul Exp $
 * @deprecated
 */
@Deprecated
public class PARMapper extends AbstractSpringAwareCLI {

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

    GeneService parService;

    TaxonService taxonService;

    Taxon taxon;

    boolean useStrand = true;

    PhysicalLocationService physicalLocationservice;

    /*
     * (non-Javadoc)
     *
     * @see ubic.gemma.util.AbstractCLI#getCommandName()
     */
    @Override
    public String getCommandName() {
        return null;
    }

    @SuppressWarnings("static-access")
    @Override
    protected void buildOptions() {

        OptionBuilder.hasArg();
        OptionBuilder.isRequired();
        OptionBuilder.withDescription( "taxon name" );
        OptionBuilder
                .withDescription( "taxon to use" );
        OptionBuilder.withLongOpt( "taxon" );
        Option taxonOption = OptionBuilder.create( 't' );
        addOption( taxonOption );

        Option useStrandOption = OptionBuilder.create( "useStrand" );
        addOption( useStrandOption );

    }

    @Override
    protected Exception doWork( String[] args ) {
        Exception exc = processCommandLine( args );
        if ( exc != null ) return exc;

        this.taxonService = this.getBean( TaxonService.class );
        this.parService = this.getBean( GeneService.class );
        this.physicalLocationservice = this.getBean( PhysicalLocationService.class );

        /*
         * Add processing o
         */

        if ( hasOption( 't' ) ) {
            String taxonName = getOptionValue( 't' );
            taxon = taxonService.findByCommonName( taxonName );
            if ( taxon == null ) {
                log.error( "ERROR: Cannot find taxon " + taxonName );
            }
        }

        this.useStrand = this.hasOption( "useStrand" );

        // Collection<ProbeAlignedRegion> pars = parService.loadProbeAlignedRegions( taxon );
        // log.info( pars.size() + " " + taxon.getCommonName() + " PARS" );
        //
        // System.out.println( "ParID\tParName\tChrom\tNuc\tGeneId\tGeneSymbol\tDistance\tGeneContainsPar\tSameStrand"
        // );

        // test case
        // ProbeAlignedRegion par = ( ProbeAlignedRegion ) parService.load(
        // 1450985 );

        // for ( ProbeAlignedRegion par : pars ) {
        // par = ( ProbeAlignedRegion ) parService.thaw( par );
        //
        // PhysicalLocation loc = par.getPhysicalLocation();
        //
        // physicalLocationservice.thaw( loc );
        //
        // RelativeLocationData nearest = parService.findNearest( loc, this.useStrand );
        // if ( nearest != null ) {
        // Gene gene = nearest.getNearestGene();
        //
        // System.out
        // .println( par.getId() + "\t" + par.getName() + "\t" + loc.getChromosome().getName() + "\t"
        // + loc.getNucleotide() + "\t" + gene.getId() + "\t" + gene.getOfficialSymbol() + "\t"
        // + nearest.getRange() + "\t" + nearest.isContainedWithinGene() + "\t"
        // + nearest.isOnSameStrand() );
        // }
        // }
        return null;
    }

}
