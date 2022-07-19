/*
 * The GemmaAnalysis project
 * 
 * Copyright (c) 2018 University of British Columbia
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

package ubic.gemma.core.apps;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.lang3.StringUtils;
import ubic.gemma.core.apps.GemmaCLI.CommandGroup;
import ubic.gemma.core.genome.gene.service.GeneService;
import ubic.gemma.core.genome.gene.service.GeneSetService;
import ubic.gemma.core.loader.genome.gene.ncbi.NcbiGeneHistoryParser;
import ubic.gemma.core.util.AbstractCLIContextCLI;
import ubic.gemma.model.association.BioSequence2GeneProduct;
import ubic.gemma.model.association.Gene2GOAssociation;
import ubic.gemma.model.association.phenotype.PhenotypeAssociation;
import ubic.gemma.model.common.description.DatabaseEntry;
import ubic.gemma.model.expression.designElement.CompositeSequence;
import ubic.gemma.model.genome.Gene;
import ubic.gemma.model.genome.Taxon;
import ubic.gemma.model.genome.biosequence.BioSequence;
import ubic.gemma.model.genome.gene.GeneProduct;
import ubic.gemma.model.genome.gene.GeneSet;
import ubic.gemma.model.genome.gene.GeneSetMember;
import ubic.gemma.model.genome.sequenceAnalysis.AnnotationAssociation;
import ubic.gemma.model.genome.sequenceAnalysis.BlatAssociation;
import ubic.gemma.persistence.service.association.Gene2GOAssociationService;
import ubic.gemma.persistence.service.association.phenotype.service.PhenotypeAssociationService;
import ubic.gemma.persistence.service.expression.designElement.CompositeSequenceService;
import ubic.gemma.persistence.service.genome.biosequence.BioSequenceService;
import ubic.gemma.persistence.service.genome.gene.GeneProductService;
import ubic.gemma.persistence.service.genome.sequenceAnalysis.AnnotationAssociationService;
import ubic.gemma.persistence.service.genome.sequenceAnalysis.BlatAssociationService;
import ubic.gemma.persistence.service.genome.taxon.TaxonService;
import ubic.gemma.persistence.util.Settings;

import java.io.File;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

/**
 * Identify, diagnose and optionally fix instances where genes are duplicated (due to changes in NCBI IDs)
 * 
 * @author paul
 */
public class GeneDuplicateResolveCli extends AbstractCLIContextCLI {

    private String filePath = Settings.getDownloadPath();

    private String taxonCommonName = null;

    private boolean doFix = false;

    /*
     * (non-Javadoc)
     * 
     * @see ubic.gemma.core.util.AbstractCLI#getCommandName()
     */
    @Override
    public String getCommandName() {
        return "geneDupResolve";
    }

    /*
     * (non-Javadoc)
     * 
     * @see ubic.gemma.core.util.AbstractCLI#buildOptions()
     */
    @Override
    protected void buildOptions( Options options ) {

        Option pathOption = Option.builder( "f" ).hasArg().argName( "Input File Path" )
                .desc( "Optional path to the directory where gene_history.gz is; default is in configed download directory" ).longOpt( "file" )
                .build();

        options.addOption( pathOption );

        options.addOption( "t", "taxon", true, "Specific taxon for which to update genes" );

        options.addOption( "fix", "Fix problems if possible; otherwise just log them" );

    }

    GeneService gs;
    GeneProductService gps;
    TaxonService ts;
    BlatAssociationService blatService;
    AnnotationAssociationService aaService;
    CompositeSequenceService csService;
    BioSequenceService bsService;
    GeneSetService gsService;

    /*
     * (non-Javadoc)
     * 
     * @see ubic.gemma.core.util.AbstractCLI#doWork(java.lang.String[])
     */
    @Override
    protected void doWork() {

        gs = this.getBean( GeneService.class );
        gps = this.getBean( GeneProductService.class );
        ts = this.getBean( TaxonService.class );
        blatService = this.getBean( BlatAssociationService.class );
        aaService = this.getBean( AnnotationAssociationService.class );
        csService = this.getBean( CompositeSequenceService.class );
        gsService = this.getBean( GeneSetService.class );
        bsService = this.getBean( BioSequenceService.class );

        PhenotypeAssociationService paService = this.getBean( PhenotypeAssociationService.class );

        try {
            String geneHistoryFile = filePath + File.separatorChar + "gene_history.gz";

            NcbiGeneHistoryParser history = new NcbiGeneHistoryParser();
            history.parse( new File( geneHistoryFile ) );

            Collection<Gene> genes;

            if ( this.taxonCommonName != null )
                genes = gs.loadAll( ts.findByCommonName( taxonCommonName ) );
            else
                genes = gs.loadAll();

            Map<Taxon, Map<String, Collection<Gene>>> perTaxonDuplicates = new HashMap<>();

            for ( Gene g : genes ) {

                Taxon tax = g.getTaxon();

                if ( !perTaxonDuplicates.containsKey( tax ) ) {
                    perTaxonDuplicates.put( tax, new HashMap<String, Collection<Gene>>() );
                }

                String sym = g.getOfficialSymbol();

                Map<String, Collection<Gene>> c1 = perTaxonDuplicates.get( tax );

                if ( !c1.containsKey( sym ) ) {
                    c1.put( sym, new HashSet<Gene>() );
                }
                c1.get( sym ).add( g );

            }

            Collection<String> couldNotFix = new HashSet<>();
            Collection<String> genesWithRemovedProbeAssociations = new HashSet<>();
            Collection<Gene> genesWithNoProducts = new HashSet<>();
            for ( Taxon t : perTaxonDuplicates.keySet() ) {

                for ( String sym : perTaxonDuplicates.get( t ).keySet() ) {
                    Collection<Gene> duplicates = perTaxonDuplicates.get( t ).get( sym );

                    if ( duplicates.size() == 1 ) continue;

                    Integer discontinuedId = null;
                    String discontinuedIdForSymbolString = history.discontinuedIdForSymbol( sym, t.getNcbiId() );
                    if ( StringUtils.isNotBlank( discontinuedIdForSymbolString ) ) {
                        discontinuedId = Integer.parseInt( discontinuedIdForSymbolString );
                    }

                    System.err.println( sym + ": " + duplicates.size() + " duplicate symbols for " + t.getCommonName() );
                    Collection<Gene> discontinued = new HashSet<>();
                    Gene useThisOne = null;
                    //  Gene maxNcbiIDGene = null; // was experimenting with using most recent ID as the current for unresolved cases, but this isn't reliable.
                    boolean canFix = false;
                    //  long maximumNcbiID = 0;

                    for ( Gene d : duplicates ) {

                        //                        if ( d.getNcbiGeneId() > maximumNcbiID ) {
                        //                            maximumNcbiID = d.getNcbiGeneId();
                        //                            maxNcbiIDGene = d;
                        //                        }

                        if ( discontinuedId != null ) {
                            if ( d.getNcbiGeneId().equals( discontinuedId ) ) {
                                discontinued.add( d );
                                System.err.println( sym + "\t" + d + "\tDiscontinued" );
                            } else {
                                if ( useThisOne != null ) {
                                    // more than one possiblility - would have to fix by hand
                                    canFix = false;
                                    System.err.println( sym + "\t" + d + "\tIn use; Multiple non-discontinued matches: " + discontinuedId );

                                } else {
                                    canFix = true;
                                    useThisOne = d;
                                    System.err.println( sym + "\t" + d + "\tUse this; Another is discontinued: " + discontinuedId );

                                }
                            }
                        } else {
                            canFix = false;
                            System.err.println( sym + "\t" + d + "\tNone discontinued" );
                        }
                    }

                    /*
                     * If some are discontinued, we know what to do. Otherwise, do we just delete the oldest NCBI IDs?
                     */
                    //  if ( useThisOne == null ) {
                    // useThisOne = maxNcbiIDGene;
                    //  canFix = true;
                    // Unfortunately, this has no general solution. Example: https://www.ncbi.nlm.nih.gov/gene/?term=317241%5Buid%5D and https://www.ncbi.nlm.nih.gov/gene/?term=108348168
                    //  System.err.println( "If fixing, will delete genes with older NCBI ids than " + maxNcbiIDGene );
                    //  } else {
                    System.err.println( "If fixing, will only retain " + useThisOne );
                    //  }

                    /*
                     * If we switch existing assocations to a "real" gene, then we have to deal with:
                     * GENE2GO_ASSOCIATION (but easily updated), PHENOTYPE_ASSOCIATION, GENE_SET_MEMBER,
                     * 
                     * 
                     * Then, probe-gene mappings should probably be retained - but this is at the gene product level,
                     * and if the products are
                     * being deleted, this is not a big deal.
                     * 
                     */
                    if ( doFix && canFix && useThisOne != null ) {
                        // Collection<GeneProduct> products = useThisOne.getProducts();

                        //                        GeneProduct soleProduct = null;
                        //                        GeneProduct annotationAssociationProduct = null;
                        //                        if ( products.size() == 1 ) {
                        //                            soleProduct = products.iterator().next();
                        //                            annotationAssociationProduct = soleProduct;
                        //                        } else {
                        //                            annotationAssociationProduct = products.iterator().next();
                        //                        }

                        useThisOne = gs.thaw( useThisOne );

                        for ( Gene d : discontinued ) {

                            System.err.println( ">> Updating/removing problems with discontinued gene " + d );

                            d = gs.thaw( d );

                            fixGOAssociations( d );

                            fixPhenotypes( paService, discontinued, useThisOne, d );

                            fixGeneSets( discontinued, useThisOne, d );

                            fixPlatformAssociations( genesWithRemovedProbeAssociations, genesWithNoProducts, sym, useThisOne, d );

                            //  gene products can share database entries with biosequences - deleting the database entry in such cases is not ok

                            if ( !d.getProducts().isEmpty() ) {
                                System.err.println( "Will be deleting " + d.getProducts().size() + " gene products from " + d );
                            }
                            for ( GeneProduct gp : d.getProducts() ) {
                                gp.getAccessions().clear();
                                gps.update( gp );
                            }

                            // else no products.
                            //   gps.remove( d.getProducts() ); // do via cascade.

                            try {
                                gs.remove( gs.load( d.getId() ) ); // avoid staleness.
                                System.err.println( "Deleted " + d );
                            } catch ( Exception e3 ) {
                                System.err.println( "Error while deleting: " + d + "; " + e3.getMessage() );
                            }
                        }
                    } else {
                        couldNotFix.add( sym );
                        System.err.println( ">>> Could not fix: " + sym );
                    }

                }
                System.err.println();
            }
            if ( !couldNotFix.isEmpty() ) {
                System.err.println( ">>>The following could not be fixed" );
                for ( String string : couldNotFix ) {
                    System.err.println( "Could not fix:\t" + string );
                }
            }

            if ( !genesWithRemovedProbeAssociations.isEmpty() ) {
                System.err.println( genesWithRemovedProbeAssociations.size() + " genes had probe (non-annotation based) assocations removed:" );
                System.err.println( StringUtils.join( genesWithRemovedProbeAssociations, "," ) );
            }

            if ( !genesWithNoProducts.isEmpty() ) {
                System.err.println( "These 'to use' genes had no products" );
                System.err.println( StringUtils.join( genesWithNoProducts, "\n" ) );
            }

        } catch ( Exception e1 ) {

        }

    }

    /**
     * @param d
     */
    void fixGOAssociations( Gene d ) {
        // go associations: just delete them
        Gene2GOAssociationService g2goService = this.getBean( Gene2GOAssociationService.class );
        Collection<Gene2GOAssociation> goassocs = g2goService.findAssociationByGene( d );
        if ( !goassocs.isEmpty() )
            System.err.println( "Removing " + goassocs.size() + " GO associations for discontinued gene " + d );
        g2goService.remove( goassocs );
    }

    /**
     * @param paService
     * @param discontinued
     * @param useThisOne
     * @param d
     */
    void fixPhenotypes( PhenotypeAssociationService paService, Collection<Gene> discontinued, Gene useThisOne, Gene d ) {
        // phenotypes: replace
        Collection<PhenotypeAssociation> pass = paService.findPhenotypeAssociationForGeneId( d.getId() );
        if ( !pass.isEmpty() ) System.err.println(
                "Updating " + pass.size() + " phenotype associations for " + discontinued + " ---> switch to " + useThisOne );
        for ( PhenotypeAssociation pas : pass ) {
            pas.setGene( useThisOne );
            paService.update( pas );
        }
    }

    /**
     * @param discontinued
     * @param useThisOne
     * @param d
     */
    void fixGeneSets( Collection<Gene> discontinued, Gene useThisOne, Gene d ) {
        // gene set: replace
        Collection<GeneSet> geneSets = gsService.findByGene( d );
        if ( !geneSets.isEmpty() ) System.err.println(
                "Updating " + geneSets.size() + " gene set associations for " + discontinued + " ---> switch to " + useThisOne );
        for ( GeneSet gset : geneSets ) {
            gsService.thaw( gset );
            GeneSetMember toRemoveFromSet = null;
            for ( GeneSetMember gsm : gset.getMembers() ) {
                if ( gsm.getGene().equals( d ) ) {
                    toRemoveFromSet = gsm;
                    break;// should only be one...
                }
            }
            if ( toRemoveFromSet != null ) {
                GeneSetMember newmember = GeneSetMember.Factory.newInstance( toRemoveFromSet.getScore(), useThisOne );
                gset.getMembers().add( newmember );
                gset.getMembers().remove( toRemoveFromSet );
                gsService.update( gset );
                assert !gsService.findByGene( d ).contains( gset );
                assert gsService.findByGene( useThisOne ).contains( gset );
                System.err.println( "Updated " + gset );
            }
        }
    }

    /**
     * @param genesWithRemovedProbeAssociations
     * @param genesWithNoProducts
     * @param sym
     * @param useThisOne
     * @param d
     */
    void fixPlatformAssociations( Collection<String> genesWithRemovedProbeAssociations, Collection<Gene> genesWithNoProducts, String sym,
            Gene useThisOne, Gene d ) {
        /*
         * Platforms - remove the associations. We will eventually regenerate them on the new
         * gene if they are valid. Note that that currently geneDao.remove takes care of these
         * already, this is just clearer and allows tracking
         */
        Collection<CompositeSequence> cses = csService.findByGene( d );
        if ( !cses.isEmpty() ) {
            System.err.println( "Removing or updating any associations to " + cses.size() + " platform element for " + d );
        }

        Collection<BlatAssociation> blatAssocRemove = new HashSet<>();
        Collection<AnnotationAssociation> annotationAssocRemove = new HashSet<>();

        boolean replacementGeneHasProducts = !useThisOne.getProducts().isEmpty();
        if ( !replacementGeneHasProducts ) {
            // yikes. The gene we're supposed to use has no products. Pseudogenes?
            System.err.println( "Warning: replacement gene has no products: " + useThisOne );
            genesWithNoProducts.add( useThisOne );
        }

        for ( CompositeSequence cs : cses ) {
            cs = csService.thaw( cs );
            for ( BioSequence2GeneProduct b2gp : cs.getBiologicalCharacteristic().getBioSequence2GeneProduct() ) {
                if ( b2gp.getGeneProduct().getGene().equals( d ) ) {
                    if ( b2gp instanceof AnnotationAssociation ) {

                        if ( replacementGeneHasProducts ) {
                            // then maybe we can use as a replacement.
                            Collection<AnnotationAssociation> find = aaService.find( useThisOne );
                            if ( find.isEmpty() ) {
                                // then we don't have a redundant one
                                // update the existing (bad) association to use this gene.
                                GeneProduct arbitraryProduct = useThisOne.getProducts().iterator().next();
                                b2gp.setGeneProduct( arbitraryProduct );

                                Collection<DatabaseEntry> accessions = arbitraryProduct.getAccessions();
                                BioSequence bs = bsService.findByAccession( accessions.iterator().next() );
                                b2gp.setBioSequence( bs );

                                aaService.update( ( AnnotationAssociation ) b2gp );
                            } else {
                                // No existing association to use; we'll delete this one, and then a new one will be made when we update the generic platform.
                                annotationAssocRemove.add( ( AnnotationAssociation ) b2gp );
                            }
                        } else {
                            // not much can be done to salvage this
                            annotationAssocRemove.add( ( AnnotationAssociation ) b2gp );
                        }
                    } else {
                        blatAssocRemove.add( ( BlatAssociation ) b2gp );
                        genesWithRemovedProbeAssociations.add( sym );
                    }
                }
            }
        }

        if ( !blatAssocRemove.isEmpty() ) {
            System.err.println( "Removing " + blatAssocRemove.size() + " blat associations for " + d );
            blatService.remove( blatAssocRemove );
        }
        if ( !annotationAssocRemove.isEmpty() ) {
            System.err.println( "Removing " + annotationAssocRemove.size() + " unneeded annotation associations for " + d );
            aaService.remove( annotationAssocRemove );
        }
    }

    @Override
    protected void processOptions( CommandLine c ) {

        if ( c.hasOption( 'f' ) ) {
            this.filePath = c.getOptionValue( 'f' );
        }
        if ( c.hasOption( "t" ) ) {
            this.taxonCommonName = c.getOptionValue( "t" );
        }
        if ( c.hasOption( "fix" ) ) {
            this.doFix = true;
        }
    }

    /*
     * (non-Javadoc)
     * 
     * @see ubic.gemma.core.util.AbstractCLIContextCLI#getCommandGroup()
     */
    @Override
    public CommandGroup getCommandGroup() {
        return CommandGroup.MISC;
    }

}
