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
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.lang.time.StopWatch;

import ubic.gemma.externalDb.GoldenPathSequenceAnalysis;
import ubic.gemma.model.expression.arrayDesign.ArrayDesign;
import ubic.gemma.model.expression.arrayDesign.ArrayDesignService;
import ubic.gemma.model.expression.designElement.CompositeSequence;
import ubic.gemma.model.expression.designElement.CompositeSequenceService;
import ubic.gemma.model.genome.Gene;
import ubic.gemma.model.genome.PhysicalLocation;
import ubic.gemma.model.genome.PhysicalLocationService;
import ubic.gemma.model.genome.Taxon;
import ubic.gemma.model.genome.gene.GeneProduct;
import ubic.gemma.model.genome.gene.GeneService;
import ubic.gemma.model.genome.sequenceAnalysis.BlatAssociation;
import ubic.gemma.model.genome.sequenceAnalysis.BlatAssociationService;
import ubic.gemma.model.genome.sequenceAnalysis.BlatResult;
import ubic.gemma.util.AbstractSpringAwareCLI;

/**
 * @author xwan
 * @version $Id$
 */
public class MicroRNAFinderCli extends AbstractSpringAwareCLI {

    private String arrayDesignName = null;
    private String outFileName = null;
    private Collection<GeneProduct> miRNAs = new HashSet<GeneProduct>();
    PhysicalLocationService plService;

    @Override
    protected void processOptions() {
        super.processOptions();
        if ( hasOption( 'a' ) ) {
            this.arrayDesignName = getOptionValue( 'a' );
        }
        if ( hasOption( 'o' ) ) {
            this.outFileName = getOptionValue( 'o' );
        }

    }

    @SuppressWarnings("static-access")
    @Override
    protected void buildOptions() {
        Option ADOption = OptionBuilder.hasArg().isRequired().withArgName( "arrayDesign" ).withDescription(
                "Array Design Short Name (GPLXXX) " ).withLongOpt( "arrayDesign" ).create( 'a' );
        addOption( ADOption );
        Option OutOption = OptionBuilder.hasArg().isRequired().withArgName( "outputFile" ).withDescription(
                "The name of the file to save the output " ).withLongOpt( "outputFile" ).create( 'o' );
        addOption( OutOption );

    }

    /**
     * @param targetLocation
     * @return
     */
    private Collection<GeneProduct> checkMappedRNAs( PhysicalLocation targetLocation ) {
        Collection<GeneProduct> returnedRNAs = new HashSet<GeneProduct>();
        if ( targetLocation == null ) {
            return returnedRNAs;
        }

        for ( GeneProduct miRNA : miRNAs ) {
            // if(!miRNA.getName().equals("mmu-let-7b")) continue;
            PhysicalLocation location = miRNA.getPhysicalLocation();
            plService.thaw( location );
            if ( targetLocation.computeOverlap( location ) >= location.getNucleotideLength() )
                returnedRNAs.add( miRNA );
        }
        return returnedRNAs;

    }

    @SuppressWarnings("unchecked")
    @Override
    protected Exception doWork( String[] args ) {
        Exception err = processCommandLine( "MicroRNAFinder", args );
        if ( err != null ) {
            return err;
        }
        ArrayDesignService adService = ( ArrayDesignService ) this.getBean( "arrayDesignService" );
        GeneService geneService = ( GeneService ) this.getBean( "geneService" );
        BlatAssociationService blatAssociationService = ( BlatAssociationService ) this
                .getBean( "blatAssociationService" );
        plService = ( PhysicalLocationService ) this.getBean( "physicalLocationService" );

        CompositeSequenceService compositeSequenceService = ( CompositeSequenceService ) this
                .getBean( "compositeSequenceService" );

        ArrayDesign arrayDesign = adService.findByShortName( this.arrayDesignName );
        if ( arrayDesign == null ) {
            System.err.println( " Array Design " + this.arrayDesignName + " doesn't exist" );
            return null;
        }

        Taxon taxon = arrayDesign.getPrimaryTaxon();

        log.info( "Loading microRNAs for " + taxon.getCommonName() );
        Collection<Gene> genes = geneService.loadMicroRNAs( taxon );

        for ( Gene gene : genes ) {
            geneService.thaw( gene );
            miRNAs.addAll( gene.getProducts() );
        }

        log.info( "Got " + miRNAs.size() + " microRNAs" );

        Map<CompositeSequence, Collection<GeneProduct>> results = new HashMap<CompositeSequence, Collection<GeneProduct>>();

        try {

            GoldenPathSequenceAnalysis analysis = new GoldenPathSequenceAnalysis( taxon );
            Collection<CompositeSequence> allCSs = adService.loadCompositeSequences( arrayDesign );

            compositeSequenceService.thaw( allCSs );

            int count = 0;
            for ( CompositeSequence cs : allCSs ) {
                // if(!cs.getName().equals("1440357_at")) continue;

                Collection bs2gps = cs.getBiologicalCharacteristic().getBioSequence2GeneProduct();
                blatAssociationService.thaw( bs2gps );

                Collection<GeneProduct> mappedRNAs = new HashSet<GeneProduct>();

                for ( Object object : bs2gps ) {

                    BlatAssociation blatAssociation = ( BlatAssociation ) object;
                    blatAssociationService.thaw( blatAssociation );
                    GeneProduct geneProduct = blatAssociation.getGeneProduct();
                    PhysicalLocation targetLocation = geneProduct.getPhysicalLocation();

                    if ( targetLocation != null ) {
                        plService.thaw( targetLocation );
                        mappedRNAs.addAll( checkMappedRNAs( targetLocation ) );
                    }

                    /*
                     * Re-check (??)
                     */
                    BlatResult blatResult = blatAssociation.getBlatResult();
                    if ( blatResult == null ) continue;

                    long start = blatResult.getTargetStart();
                    long end = blatResult.getTargetEnd();

                    String chromosome = "chr" + blatResult.getTargetChromosome().getName();

                    log.debug( chromosome );

                    mappedRNAs.addAll( analysis.findMicroRNAGenesByLocation( chromosome, start, end, blatResult
                            .getStrand() ) );

                    Collection<Gene> alignedGenes = analysis.findRNAs( chromosome, start, end, blatResult.getStrand() );
                    for ( Gene gene : alignedGenes ) {
                        for ( GeneProduct gp : gene.getProducts() ) {
                            PhysicalLocation loc = gp.getPhysicalLocation();
                            mappedRNAs.addAll( checkMappedRNAs( loc ) );
                        }
                    }
                }

                if ( mappedRNAs.size() > 0 ) results.put( cs, mappedRNAs );

                if ( ++count % 1000 == 0 ) {
                    log.info( "Processed " + count + " probes" );
                }
            }

            /*
             * Output phase
             */
            PrintStream output = new PrintStream( new FileOutputStream( new File( this.outFileName ) ) );
            for ( CompositeSequence cs : results.keySet() ) {
                output.print( cs.getName() + "\t" );
                Collection<GeneProduct> mappedmiRNAs = results.get( cs );
                for ( GeneProduct miRNA : mappedmiRNAs ) {
                    output.print( miRNA.getName() + "," );
                }
                output.println();
            }
            output.close();
        } catch ( Exception e ) {
            return e;
        }
        return null;
    }

    /**
     * @param args
     */
    public static void main( String[] args ) {
        MicroRNAFinderCli finder = new MicroRNAFinderCli();
        StopWatch watch = new StopWatch();
        watch.start();
        try {
            Exception ex = finder.doWork( args );
            if ( ex != null ) {
                ex.printStackTrace();
            }
            watch.stop();
            log.info( "Elapsed time: " + watch.getTime() / 1000 + " seconds" );
        } catch ( Exception e ) {
            throw new RuntimeException( e );
        }
    }

}
