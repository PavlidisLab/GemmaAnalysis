/*
 * The Gemma project
 *
 * Copyright (c) 2011 University of British Columbia
 *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with
 * the License. You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on
 * an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the
 * specific language governing permissions and limitations under the License.
 */
package ubic.gemma.contrib.apps;

import org.springframework.beans.factory.annotation.Autowired;
import ubic.basecode.ontology.model.OntologyTerm;
import ubic.basecode.ontology.providers.ExperimentalFactorOntologyService;
import ubic.basecode.ontology.search.OntologySearchException;
import ubic.gemma.core.analysis.expression.diff.BaselineSelection;
import ubic.gemma.core.analysis.preprocess.batcheffects.BatchConfound;
import ubic.gemma.core.analysis.preprocess.batcheffects.BatchConfoundUtils;
import ubic.gemma.core.analysis.preprocess.svd.SVDService;
import ubic.gemma.core.analysis.preprocess.svd.SVDValueObject;
import ubic.gemma.core.apps.ExpressionExperimentManipulatingCLI;
import ubic.gemma.core.apps.GemmaCLI.CommandGroup;
import ubic.gemma.core.ontology.OntologyService;
import ubic.gemma.core.ontology.OntologyUtils;
import ubic.gemma.model.common.description.Characteristic;
import ubic.gemma.model.expression.experiment.*;
import ubic.gemma.persistence.service.expression.experiment.ExperimentalDesignService;
import ubic.gemma.persistence.util.EntityUtils;

import java.util.*;

/**
 *
 * @author paul
 */
public class FactorValueCharacteristicAnalysis extends ExpressionExperimentManipulatingCLI {

    @Override
    public CommandGroup getCommandGroup() {
        return CommandGroup.EXPERIMENT;
    }

    /*
     * (non-Javadoc)
     *
     * @see ubic.gemma.util.AbstractCLI#getCommandName()
     */
    @Override
    public String getCommandName() {
        return "fvAnalysis";
    }

    @Override
    public String getShortDesc() {
        return "Analyze factor value characteristics for grouping oppotunities";
    }


    private OntologyService ontologyService;

    /*
     * (non-Javadoc)
     *
     * @see ubic.gemma.util.AbstractCLI#doWork(java.lang.String[])
     */
    @Override
    protected void doWork() {

        ExperimentalDesignService eds = this.getBean( ExperimentalDesignService.class ); // autowire not available here

        this.ontologyService = this.getBean( OntologyService.class );

        try {
            OntologyUtils.ensureInitialized( this.getBean( ExperimentalFactorOntologyService.class ) );
        } catch ( InterruptedException e ) {
            throw new RuntimeException( e );
        }

        int total = 0;
        int nosolution = 0;

        for ( BioAssaySet bas : this.getExpressionExperiments() ) {
            if ( bas instanceof ExpressionExperiment ) {
                ExpressionExperiment ee = ( ExpressionExperiment ) bas;
                this.getEeService().thawLite( ee );
                log.info( "Processing: " + ee );
                System.out.println( "================" );
                System.out.println( ee );

                try {

                    boolean success = false;

                    // retrieve the experimental factors
                    ExperimentalDesign ed = eds.loadWithExperimentalFactors( ee.getExperimentalDesign().getId() );
                    for ( ExperimentalFactor factor : ed.getExperimentalFactors() ) {
                        for ( FactorValue fv : factor.getFactorValues() ) {

                            if ( fv.getCharacteristics().size() <= 1 ) continue;

                            System.out.println( "--------------------" );
                            System.out.println( fv );

                            fv.getCharacteristics().forEach( c -> {
                                System.out.println( " - c - " + c.getCategory() + ": " + c.getValue() + " " + ( c.getValueUri() == null ? "[free text]" : c.getValueUri() ) );
                            } );


                            /*
                             * If there two characteristics, and one is Dose and one is Treatment, print it out like Treatment has_dose Dose
                             */
                            if ( fv.getCharacteristics().size() == 2 ) {

                                total++;

                                Characteristic[] cs = fv.getCharacteristics().toArray( new Characteristic[0] );
                                Characteristic cs1 = cs[0];
                                Characteristic cs2 = cs[1];

                                // drug - dose
                                if ( cs1.getCategory() != null && cs1.getCategory().equals( "dose" ) ) {
                                    if ( cs2.getCategory() != null && cs2.getCategory().equals( "treatment" ) ) {
                                        System.out.println( "Treatment -> " + cs2.getValue() + " [" + ( cs2.getValueUri() == null ? "free text" : cs2.getValueUri() ) + "]"
                                                + " has_dose " + cs1.getValue() + " [" + ( cs1.getValueUri() == null ? "free text" : cs1.getValueUri() ) + "]" );
                                    }
                                } else if ( cs2.getCategory() != null && cs2.getCategory().equals( "dose" ) ) {
                                    if ( cs1.getCategory() != null && cs1.getCategory().equals( "treatment" ) ) {
                                        System.out.println( "Treatment -> " + cs1.getValue() + " [" + ( cs1.getValueUri() == null ? "free text" : cs1.getValueUri() ) + "]"
                                                + " has_dose " + cs2.getValue() + " [" + ( cs2.getValueUri() == null ? "free text" : cs2.getValueUri() ) + "]" );
                                    }

                                    // gene - `genotype
                                } else if ( cs1.getCategory() != null && cs1.getCategory().equals( "genotype" ) && cs1.getValueUri() != null && cs1.getValueUri().contains( "ncbi_gene" ) ) {
                                    if ( cs2.getCategory() != null && cs2.getCategory().equals( "genotype" ) ) {
                                        System.out.println( "Genotype -> " + cs1.getValue() + " [" + ( cs1.getValueUri() == null ? "free text" : cs1.getValueUri() ) + "]"
                                                + " has_genotype " + cs2.getValue() + " [" + ( cs2.getValueUri() == null ? "free text" : cs2.getValueUri() ) + "]" );
                                    }
                                } else if ( cs2.getCategory() != null && cs2.getCategory().equals( "genotype" ) && cs2.getValueUri() != null && cs2.getValueUri().contains( "ncbi_gene" ) ) {
                                    if ( cs1.getCategory() != null && cs1.getCategory().equals( "genotype" ) ) {
                                        System.out.println( "Genotype -> " + cs2.getValue() + " [" + ( cs2.getValueUri() == null ? "free text" : cs2.getValueUri() ) + "]"
                                                + " has_genotype " + cs1.getValue() + " [" + ( cs1.getValueUri() == null ? "free text" : cs1.getValueUri() ) + "]" );
                                    }
                                }


                                // developmental stages described by free text
                                else if ( cs1.getCategory() != null && cs1.getCategory().equals( "developmental stage" ) && isFormalDevelopmentalStage( cs1 ) ) {
                                    if ( cs2.getCategory() != null && cs2.getCategory().equals( "developmental stage" ) && cs2.getValueUri() == null ) {
                                        System.out.println( "Stage -> " + cs1.getValue() + " [" + ( cs1.getValueUri() == null ? "free text" : cs1.getValueUri() ) + "]"
                                                + " has_embryo_stage " + cs2.getValue() + " [" + ( cs2.getValueUri() == null ? "free text" : cs2.getValueUri() ) + "]" );

                                    }
                                } else if ( cs2.getCategory() != null && cs2.getCategory().equals( "developmental stage" ) && isFormalDevelopmentalStage( cs2 ) ) {
                                    if ( cs1.getCategory() != null && cs1.getCategory().equals( "developmental stage" ) && cs1.getValueUri() == null ) {
                                        System.out.println( "Stage -> " + cs2.getValue() + " [" + ( cs2.getValueUri() == null ? "free text" : cs2.getValueUri() ) + "]"
                                                + " has_embryo_stage " + cs1.getValue() + " [" + ( cs1.getValueUri() == null ? "free text" : cs1.getValueUri() ) + "]" );
                                    }
                                }

                                // an UBERON term and a location modifier
                                else if ( cs1.getCategory() != null && cs1.getCategory().equals( "organism part" ) && cs2.getValueUri() != null && cs2.getValueUri().contains( "PATO_" ) ) {
                                    System.out.println( "OrganismPart -> " + cs1.getValue() + " [" + ( cs1.getValueUri() == null ? "free text" : cs1.getValueUri() ) + "]"
                                            + " has_location[?] " + cs2.getValue() + " [" + ( cs2.getValueUri() == null ? "free text" : cs2.getValueUri() ) + "]" );
                                } else if ( cs2.getCategory() != null && cs2.getCategory().equals( "organism part" ) && cs1.getValueUri() != null && cs1.getValueUri().contains( "PATO_" ) ) {
                                    System.out.println( "Organism Part location -> " + cs2.getValue() + " [" + ( cs2.getValueUri() == null ? "free text" : cs2.getValueUri() ) + "]"
                                            + " has_location[?] " + cs1.getValue() + " [" + ( cs1.getValueUri() == null ? "free text" : cs1.getValueUri() ) + "]" );
                                }

                                // fall through to generic case where one value is an ontology term and the other is free text. We assume the free text is the modifier.
                                else if ( cs1.getValueUri() != null && cs2.getValueUri() == null ) {
                                    System.out.println( "Generic -> " + cs1.getValue() + " [" + cs1.getValueUri() + "]"
                                            + " has_modifier " + cs2.getValue() );

                                } else if ( cs2.getValueUri() != null && cs1.getValueUri() == null ) {

                                    System.out.println( "Generic -> " + cs2.getValue() + " [" + cs2.getValueUri() + "]"
                                            + " has_modifier " + cs1.getValue() );
                                } else {
                                    System.out.println( "*** NO SOLUTION *** " );
                                    nosolution++;
                                }
                            }


                            System.out.println( "--------------------" );

                        }
                    }

//
//                    success = true;
//                    if ( success ) {
//                        this.addSuccessObject( bas.toString(), "" );
//                    } else {
//                        this.addErrorObject( bas.toString() + ": No dates found", "" );
//
//                    }


                } catch ( Exception e ) {
                    log.error( e, e );
                    this.addErrorObject( bas + ": " + e.getMessage(), "" );
                }

            }

        }
        System.err.println( "Total two-characteristic = " + total + " of which " + ( total - nosolution ) + " had a solution (" + String.format( "%.2f", 100.0 * ( 1.0 - ( double ) nosolution / total ) ) + "%)" );

    }


    private List<String> formalDevStageUris = Arrays.asList(
            /* FIXME can probably use
                http://purl.obolibrary.org/obo/UBERON_0000105 (life cycle stage)
                and http://www.ebi.ac.uk/efo/EFO_0000399 (developmental stage) children */
            /* embryo stage */ "http://purl.obolibrary.org/obo/UBERON_0000068",
            /*neonate */ "http://purl.obolibrary.org/obo/UBERON_0007221",
            /*infant*/ "http://www.ebi.ac.uk/efo/EFO_0001355" );


    private List<String> formalLocationUris = Arrays.asList(
            /*dorsal*/ "http://www.ebi.ac.uk/efo/EFO_0001656",
            /* ventral */ "http://www.ebi.ac.uk/efo/EFO_0001662" );

    private boolean isFormalLocation( Characteristic cs ) throws Exception {
        if ( cs.getValueUri() == null ) return false;
        String patophys = "http://purl.obolibrary.org/obo/PATO_0001241"; /* physical object quality*/
        String efoanatmod = "http://www.ebi.ac.uk/efo/EFO_0001646"; /* anatomical modifier */

        Collection<OntologyTerm> terms = getChildTerms( patophys );
        terms.addAll( getChildTerms( efoanatmod ) );

        /* get the valueUris from the terms */
        Set<String> uris = new HashSet<>();
        for ( OntologyTerm t : terms ) {
            uris.add( t.getUri() );
        }

        return uris.contains( cs.getValueUri() );

    }

    /**
     * FIXME put this in OntologyService probably
     * FIXME are these cached?
     * @param uri
     * @return
     * @throws OntologySearchException
     */
    private Collection<OntologyTerm> getChildTerms( String uri ) throws OntologySearchException {
        Collection<OntologyTerm> terms = ontologyService.findTerms( uri );
        if ( terms.size() > 1 ) throw new IllegalStateException( "More than one term found for " + uri );
        if ( terms.size() == 0 ) throw new IllegalStateException( "No term found for " + uri );
        terms = ontologyService.getChildren( terms, false, false );
        return terms;
    }

    private boolean isControlCondition( Characteristic cs ) {
        return BaselineSelection.isBaselineCondition( cs );
    }

    private boolean isFormalDevelopmentalStage( Characteristic cs ) {

        if ( cs.getValueUri() == null ) return false;

        return formalDevStageUris.contains( cs.getValueUri() );
    }


}
