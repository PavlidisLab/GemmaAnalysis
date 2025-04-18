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

import org.apache.commons.lang3.StringUtils;
import org.springframework.beans.factory.annotation.Autowired;
import ubic.gemma.apps.ExpressionExperimentManipulatingCLI;
import ubic.gemma.core.analysis.expression.diff.BaselineSelection;
import ubic.gemma.model.common.description.Characteristic;
import ubic.gemma.model.expression.experiment.*;
import ubic.gemma.persistence.service.expression.experiment.ExperimentalDesignService;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 *
 * @author paul
 */
public class FactorValueCharacteristicAnalysis extends ExpressionExperimentManipulatingCLI {

    @Autowired
    private ExperimentalDesignService eds;

    Collection<RemappingInfo> results = new HashSet<>();
    Collection<RemappingInfo> unresolved = new HashSet<>();
    private int numSingletons = 0;
    private int numMultiplex = 0;
    private int numDuplex = 0;
    private int nosolution = 0;
    private int experimentsExamined = 0; // not counting ones that have no experimental design
    private int numExperiments;

    @Override
    public String getCommandName() {
        return "fvAnalysis";
    }

    @Override
    public String getShortDesc() {
        return "Analyze factor value characteristics for grouping oppotunities";
    }

    @Override
    protected void doAuthenticatedWork() throws Exception {

        log.info( "Starting examining experiments ..." );

        // this.ontologyService = this.getBean( OntologyService.class );

//        try {
//            System.err.println( "Initializing ontologies ..." );
//
//            log.info( "Initializing EFO" );
//            OntologyUtils.ensureInitialized( this.getBean( ExperimentalFactorOntologyService.class ) );
//            log.info( "EFO initialized" );
//            System.err.println( "Initialized EFO" );
//
//            log.info( "Initializing MONDO" );
//            OntologyUtils.ensureInitialized( this.getBean( MondoOntologyService.class ) );
//            log.info( "MONDO initialized" );
//            System.err.println( "Initialized MONDO" );
//        } catch ( InterruptedException e ) {
//            throw new RuntimeException( e );
//        }

        super.doAuthenticatedWork();

        System.err.println( "Total experiments examined = " + experimentsExamined );
        System.err.println( "Total two-characteristic factor values = " + numDuplex + " of which " + ( numDuplex - nosolution ) + " had a solution (" + String.format( "%.2f", 100.0 * ( 1.0 - ( double ) nosolution / numDuplex ) ) + "%)" );
        System.err.println( "Total singleton factor values = " + numSingletons );
        System.err.println( "Total multiplex factor values = " + numMultiplex );

        /*
        Iterate over results and call tabularize() on each ri, and print to a file.
         */
        try {
            String filename = "remapping-" + new Date().toString().replace( " ", "_" ).replace( ":", "_" ) + ".tsv";
            BufferedWriter writer = Files.newBufferedWriter( Paths.get( filename ) );
            writer.write( "Experiment\tExperiment ID\tFactorValue\tFVID\tSummary\tC1 ID\tC2 ID\tSubjectCategory\tSubject\tPredicate\tObjectCategory\tObject\n" );
            for ( RemappingInfo ri : results ) {
                String c = ri.tabularize();
                if ( c == null ) continue;
                writer.write( c + "\n" );
            }
            writer.close();
            System.err.println( "Wrote results to " + Paths.get( filename ) );
        } catch ( IOException e ) {
            throw new RuntimeException( e );
        }


        try {
            String filename = "unresolved-" + new Date().toString().replace( " ", "_" ).replace( ":", "_" ) + ".tsv";
            BufferedWriter writer = Files.newBufferedWriter( Paths.get( filename ) );
            writer.write( "Experiment\tExperiment ID\tFactorValue\tFVID\tSummary\tC1 ID\tC2 ID\tSubjectCategory\tSubject\tPredicate\tObjectCategory\tObject\n" );
            for ( RemappingInfo ri : unresolved ) {
                String c = ri.tabularize();
                if ( c == null ) continue;
                writer.write( c + "\n" );
            }
            writer.close();
            System.err.println( "Wrote unresolved to " + Paths.get( filename ) );
        } catch ( IOException e ) {
            throw new RuntimeException( e );
        }
    }

    @Override
    protected void processBioAssaySets( Collection<BioAssaySet> expressionExperiments ) {
        numExperiments = expressionExperiments.size();
        super.processBioAssaySets( expressionExperiments );
    }

    @Override
    protected void processExpressionExperiment( ExpressionExperiment ee ) throws Exception {
        if ( experimentsExamined > 0 && experimentsExamined % 500 == 0 ) {
            log.info( "Processed " + experimentsExamined + " experiments / " + numExperiments );
        }

        try {

            boolean success = false;

            // retrieve the experimental factors
            ExperimentalDesign ed = eds.loadWithExperimentalFactors( ee.getExperimentalDesign().getId() );

            if ( !ed.getExperimentalFactors().isEmpty() ) {
                experimentsExamined++;
                // log.info( "Processing: " + ee );
                // System.out.println( "================" );
            }

            for ( ExperimentalFactor factor : ed.getExperimentalFactors() ) {
                for ( FactorValue fv : factor.getFactorValues() ) {

                    if ( fv.getCharacteristics().size() <= 1 ) {
                        numSingletons++;
                        continue;
                    }

                    boolean solved = false;
                    RemappingInfo ri = new RemappingInfo( ee, fv ); // note: for some complex cases we will recreate this


                    /*
                     * Handle cases like If there two characteristics, and one is Dose and one is Treatment, print it out like Treatment has_dose Dose
                     */
                    if ( fv.getCharacteristics().size() == 2 ) {

                        numDuplex++;

                        Characteristic[] cs = fv.getCharacteristics().toArray( new Characteristic[0] );
                        Characteristic cs1 = cs[0];
                        Characteristic cs2 = cs[1];

                        // drug - dose
                        if ( cs1.getCategory() != null && cs1.getCategory().equals( "dose" ) && cs2.getCategory() != null && cs2.getCategory().equals( "treatment" ) ) {
                            ri.summary = "Treatment-dose";
                            ri.subject = cs2;
                            ri.object = cs1;
                            ri.predicate = "has_dose";
                            solved = true;
                        } else if ( cs2.getCategory() != null && cs2.getCategory().equals( "dose" ) && cs1.getCategory() != null && cs1.getCategory().equals( "treatment" ) ) {
                            ri.summary = "Treatment-dose";
                            ri.subject = cs1;
                            ri.object = cs2;
                            ri.predicate = "has_dose";
                            solved = true;

                            /* case where we have a dose but the other thing isn't obviously a drug - we just treat it as the modifier anyway */
                        } else if ( cs1.getCategory() != null && cs1.getCategory().equals( "dose" ) ) {
                            ri.summary = "Generic dosing ";
                            ri.subject = cs2;
                            ri.object = cs1;
                            ri.predicate = "has_dose";
                            solved = true;

                        } else if ( cs2.getCategory() != null && cs2.getCategory().equals( "dose" ) ) {
                            ri.summary = "Generic dosing ";
                            ri.subject = cs1;
                            ri.object = cs2;
                            ri.predicate = "has_dose";
                            solved = true;


                            // gene - `genotype -- assumine the genotype is a genetic manipulation or some details of the genotype.
                        } else if ( isGene( cs1 ) &&
                                cs2.getCategory() != null && cs2.getCategory().equals( "genotype" ) ) {
                            ri.summary = "Genotype ";
                            ri.subject = cs1;
                            ri.object = cs2;
                            ri.predicate = "has_genotype";
                            solved = true;


                        } else if ( isGene( cs2 ) &&
                                cs1.getCategory() != null && cs1.getCategory().equals( "genotype" ) ) {
                            ri.summary = "Genotype ";
                            ri.subject = cs2;
                            ri.object = cs1;
                            ri.predicate = "has_genotype";
                            solved = true;

                        }


                                /*
                                A cell type and a cell line. We assume this means "this cell type from this cell line"
                                 */
                        else if ( cs1.getCategory() != null && cs1.getCategory().equals( "cell type" ) && cs2.getCategory() != null && cs2.getCategory().equals( "cell line" ) ) {
                            ri.summary = "Cell type+line";
                            ri.subject = cs1;
                            ri.object = cs2;
                            ri.predicate = "from_cell_line";
                            solved = true;

                        } else if ( cs2.getCategory() != null && cs2.getCategory().equals( "cell type" ) && cs1.getCategory() != null && cs1.getCategory().equals( "cell line" ) ) {
                            ri.summary = "Cell type+line";
                            ri.subject = cs2;
                            ri.object = cs1;
                            ri.predicate = "from_cell_line";
                            solved = true;

                        }


                        // developmental stages described by free text
                        else if ( cs1.getCategory() != null && cs1.getCategory().equals( "developmental stage" ) && isDevelopmentalStage( cs1 ) &&
                                cs2.getCategory() != null && cs2.getCategory().equals( "developmental stage" ) && cs2.getValueUri() == null ) {
                            ri.summary = "Stage";
                            ri.subject = cs1;
                            ri.object = cs2;
                            ri.predicate = "has_stage";
                            solved = true;


                        } else if ( cs2.getCategory() != null && cs2.getCategory().equals( "developmental stage" ) && isDevelopmentalStage( cs2 ) &&
                                cs1.getCategory() != null && cs1.getCategory().equals( "developmental stage" ) && cs1.getValueUri() == null ) {
                            ri.summary = "Stage";
                            ri.subject = cs2;
                            ri.object = cs1;
                            ri.predicate = "has_stage";
                            solved = true;
                        }

                        // an UBERON term and a location modifier
                        else if ( cs1.getCategory() != null && cs1.getCategory().equals( "organism part" ) && isLocation( cs2 ) ) {
                            ri.summary = "OrganismPart location";
                            ri.subject = cs1;
                            ri.object = cs2;
                            ri.predicate = "has_location";
                            solved = true;

                        } else if ( cs2.getCategory() != null && cs2.getCategory().equals( "organism part" ) && isLocation( cs1 ) ) {
                            ri.summary = "OrganismPart location";
                            ri.subject = cs2;
                            ri.object = cs1;
                            ri.predicate = "has_location";
                            solved = true;
                        }



                                /*  disease in an organism part (location)
                                *
                                * FactorValue 134421: phenotype:kidney | disease related to solid organ transplantation |
 - c - organism part: kidney http://purl.obolibrary.org/obo/UBERON_0002113
 - c - disease: disease related to solid organ transplantation http://purl.obolibrary.org/obo/MONDO_0700221 */

                        else if ( cs1.getCategory() != null && cs1.getCategory().equals( "organism part" ) && cs2.getCategory() != null && isDisease( cs2 ) ) {
                            ri.summary = "Disease in part";
                            ri.subject = cs2;
                            ri.object = cs1;
                            ri.predicate = "has_location";
                            solved = true;

                        } else if ( cs2.getCategory() != null && cs2.getCategory().equals( "organism part" ) && cs1.getCategory() != null && isDisease( cs1 ) ) {
                            ri.summary = "Disease in part";
                            ri.subject = cs1;
                            ri.object = cs2;
                            ri.predicate = "has_location";
                            solved = true;

                        }


                        /* disease + stage */
                        else if ( cs1.getCategory() != null && cs1.getCategory().equals( "disease staging" ) && isDisease( cs2 ) ) {
                            ri.summary = "Disease stage";
                            ri.subject = cs2;
                            ri.object = cs1;
                            ri.predicate = "has_stage";
                            solved = true;

                        } else if ( cs2.getCategory() != null && cs2.getCategory().equals( "disease staging" ) && isDisease( cs1 ) ) {
                            ri.summary = "Disease stage";
                            ri.subject = cs1;
                            ri.object = cs2;
                            ri.predicate = "has_stage";
                            solved = true;

                        }





                        /* disease + modifier */
                        else if ( isDiseaseModifier( cs1 ) && isDisease( cs2 ) ) {
                            ri.summary = "Disease with modifier";
                            ri.subject = cs2;
                            ri.object = cs1;
                            ri.predicate = "has_modifier";
                            solved = true;

                        } else if ( isDiseaseModifier( cs2 ) && isDisease( cs1 ) ) {
                            ri.summary = "Disease with modifier";
                            ri.subject = cs1;
                            ri.object = cs2;
                            ri.predicate = "has_modifier";
                            solved = true;

                        }


                        /* cell type from a particular organism part */
                                /*
                                 - c - organism part: astrocyte http://purl.obolibrary.org/obo/CL_0000127
                                 - c - organism part: prefrontal cortex http://purl.obolibrary.org/obo/UBERON_0000451

                                 */
                        else if ( cs1.getCategory() != null && cs1.getCategory().equals( "organism part" ) && cs2.getCategory() != null && cs2.getCategory().equals( "cell type" ) ) {
                            ri.summary = "Cell type from part";
                            ri.subject = cs2;
                            ri.object = cs1;
                            ri.predicate = "from_location";
                            solved = true;


                        } else if ( cs2.getCategory() != null && cs2.getCategory().equals( "organism part" ) && cs1.getCategory() != null && cs1.getCategory().equals( "cell type" ) ) {
                            ri.summary = "Cell type from part";
                            ri.subject = cs1;
                            ri.object = cs2;
                            ri.predicate = "from_location";
                            solved = true;

                        }


                                /*  drug -> delivery (children of OBI "adding a material entity into a target"[http://purl.obolibrary.org/obo/OBI_0000274]?)
                                    FactorValue 137132: treatment:sodium metaarsenite | intraperitoneal injection |
                                    - c - treatment: sodium metaarsenite http://purl.obolibrary.org/obo/CHEBI_29678
                                    - c - treatment: intraperitoneal injection http://purl.obolibrary.org/obo/OBI_0000281
                                  */
                        else if ( isDrug( cs2 ) && isDelivery( cs1 ) ) {
                            ri.summary = "Drug with delivery";
                            ri.subject = cs2;
                            ri.object = cs1;
                            ri.predicate = "has_delivery";
                            solved = true;

                        } else if ( isDrug( cs1 ) && isDelivery( cs2 ) ) {
                            ri.summary = "Drug with delivery";
                            ri.subject = cs1;
                            ri.object = cs2;
                            ri.predicate = "has_delivery";
                            solved = true;

                        }


                        // control + reference subject role etc. We assume that 'reference subject role' (or whatever) has_role 'control', I guess.
                                /*
                                FactorValue 144634: treatment:control | reference subject role |
 - c - treatment: control http://www.ebi.ac.uk/efo/EFO_0001461
 - c - treatment: reference subject role http://purl.obolibrary.org/obo/OBI_0000220

 FactorValue 147812: genotype:wild type genotype | control |
 - c - genotype: wild type genotype http://www.ebi.ac.uk/efo/EFO_0005168
 - c - genotype: control http://www.ebi.ac.uk/efo/EFO_0001461
                                 */
                        else if ( cs2.getValueUri() != null && cs2.getValueUri().equals( "http://www.ebi.ac.uk/efo/EFO_0001461" /*control*/ ) ) {
                            ri.summary = "Control";
                            ri.subject = cs1;
                            ri.object = cs2;
                            ri.predicate = "has_role[?]";
                            solved = true;

                        } else if ( cs1.getValueUri() != null && cs1.getValueUri().equals( "http://www.ebi.ac.uk/efo/EFO_0001461" /*control*/ ) ) {
                            ri.summary = "Control";
                            ri.subject = cs2;
                            ri.object = cs1;
                            ri.predicate = "has_role[?]";
                            solved = true;
                        }

                                /*
                                when we have two cell types, and one of them is either iPSC derived cell line or ESC derived cell line
                                 */
                        else if ( cs2.getCategory() != null && cs1.getCategory() != null && cs1.getCategory().equals( "cell type" ) && cs2.getCategory().equals( "cell type" ) ) {
                            if ( cs1.getValueUri() != null && ( cs1.getValueUri().equals( "http://www.ebi.ac.uk/efo/EFO_0005738" ) || cs1.getValueUri().equals( "http://www.ebi.ac.uk/efo/EFO_0005740" ) ) ) {
                                ri.summary = "Derived cells";
                                ri.subject = cs1;
                                ri.object = cs2;
                                ri.predicate = "cell_derived_from";
                                solved = true;

                            } else if ( cs2.getValueUri() != null && ( cs2.getValueUri().equals( "http://www.ebi.ac.uk/efo/EFO_0005738" ) || cs2.getValueUri().equals( "http://www.ebi.ac.uk/efo/EFO_0005740" ) ) ) {
                                ri.summary = "Derived cells";
                                ri.subject = cs2;
                                ri.object = cs1;
                                ri.predicate = "cell_derived_from";
                                solved = true;

                            }
                        }



                                   /*  chemical has_role reference substance role
                                --------------------
                                FactorValue 165331: treatment:reference substance role | DMSO |
                                    - c - treatment: reference substance role http://purl.obolibrary.org/obo/OBI_0000025
                                 - c - treatment: DMSO http://purl.obolibrary.org/obo/CHEBI_2826
                                 and so on:

                                 FactorValue 186321: phenotype:asynchronous | reference subject role |
                                - c - phenotype: asynchronous http://purl.obolibrary.org/obo/PATO_0000688
                                   - c - phenotype: reference subject role http://purl.obolibrary.org/obo/OBI_0000220


                                iment Id=23348 Name=Transcriptomic analysis of corticosteroids-treated hiPSC-derived RPE cells Short Name=GSE172478
                                FactorValue 187865: treatment:ethanol | reference substance role |
                                - c - treatment: ethanol http://purl.obolibrary.org/obo/CHEBI_16236
                                    - c - treatment: reference substance role http://purl.obolibrary.org/obo/OBI_0000025
                                 */
                        else if ( isControlCondition( cs1 ) ) {
                            ri.summary = "Control";
                            ri.subject = cs2;
                            ri.object = cs1;
                            ri.predicate = "has_role";
                            solved = true;

                        } else if ( isControlCondition( cs2 ) ) {
                            ri.summary = "Control";
                            ri.subject = cs1;
                            ri.object = cs2;
                            ri.predicate = "has_role";
                            solved = true;

                        }

                        /* something and a timepoint where we don't really know what the other thing is*/
                        else if ( isTimepoint( cs1 ) && !isTimepoint( cs2 ) ) {
                            ri.summary = "Timepoint";
                            ri.object = cs1;
                            ri.subject = cs2;
                            ri.predicate = "has_timepoint";
                            solved = true;
                        } else if ( isTimepoint( cs2 ) && !isTimepoint( cs1 ) ) {
                            ri.summary = "Timepoint";
                            ri.object = cs2;
                            ri.subject = cs1;
                            ri.predicate = "has_timepoint";
                            solved = true;

                        }


                        // somewhat generic case of where there is a physical object property or occurrence type of term, which we treat as the modifier.
                        // The exception here is when one is free text, is might be a dose or details like for "radiation exposure"
                        else if ( isQualityProperty( cs1 ) ) {

                            if ( cs2.getValueUri() == null ) {
                                // then print is as cs1 cs2
                                ri.summary = "Property";
                                ri.subject = cs1;
                                ri.object = cs2;
                                ri.predicate = "has_property[?]";
                                solved = true;

                            } else {
                                // also not quite clear - cs2 is free text.
                                ri.summary = "Physical property";
                                ri.subject = cs2;
                                ri.object = cs1;
                                ri.predicate = "has_property[?]";
                                solved = true;

                            }

                        } else if ( isQualityProperty( cs2 ) ) {
                            if ( cs1.getValueUri() == null ) {
                                ri.summary = "Property";
                                ri.subject = cs2;
                                ri.object = cs1;
                                ri.predicate = "has_property[?]";
                                solved = true;

                            } else {
                                ri.summary = "Physical property";
                                ri.subject = cs1;
                                ri.object = cs2;
                                ri.predicate = "has_property[?]";
                                solved = true;

                            }
                        } // FIXME also deal with children of http://purl.obolibrary.org/obo/PATO_0000069 - deviation (from_normal)



                                /*
                                If it's just two drugs, we can't do anything, it's just how it is.
                                 */
                        else if ( isDrug( cs1 ) && isDrug( cs2 ) ) {
                            // we count this as solved.
                            ri.summary = "Two unrelated drugs";
                            ri.subject = cs1;
                            ri.object = cs2;
                            ri.predicate = null;
                            solved = true;

                        }


                        // fall through to generic case where one value is an ontology term and the other is free text. We assume the free text is the modifier.
                        else if ( cs1.getValueUri() != null && cs2.getValueUri() == null ) {
                            ri.summary = "Generic";
                            ri.subject = cs1;
                            ri.object = cs2;
                            ri.predicate = "has_modifier[?]";
                            solved = true;


                        } else if ( cs2.getValueUri() != null && cs1.getValueUri() == null ) {
                            ri.summary = "Generic";
                            ri.subject = cs2;
                            ri.object = cs1;
                            ri.predicate = "has_modifier[?]";
                            solved = true;


                        } else {
                            solved = false;
                            nosolution++;
                        }

                    } else /* more than 2 characteristics */ {
                        numMultiplex++;

                        // case of one genetic modifier, and > 1 gene: assume modifier applies to each.
                        // which case, we need to duplicate the modifier
                        int numGenes = 0;
                        boolean hasgeneticmodifier = false;
                        for ( Characteristic c : fv.getCharacteristics() ) {
                            if ( isGene( c ) ) {
                                numGenes++;
                            } else if ( isGeneticManipulation( c ) ) {
                                hasgeneticmodifier = true;
                            }
                        }
                        if ( hasgeneticmodifier && numGenes == fv.getCharacteristics().size() - 1 ) {
                            // then we have a genetic modifier and a bunch of genes
                            // we need to duplicate the modifier
                            Characteristic geneticmodifier = null;
                            for ( Characteristic c : fv.getCharacteristics() ) {
                                if ( isGeneticManipulation( c ) ) {
                                    geneticmodifier = c;
                                    break;
                                }
                            }
                            if ( geneticmodifier == null ) {
                                throw new IllegalStateException( "Didn't find genetic modifier for: " + ee + " " + fv );
                            }
                            for ( Characteristic c : fv.getCharacteristics() ) {
                                if ( c.equals( geneticmodifier ) ) continue;
                                ri = new RemappingInfo( ee, fv );
                                ri.subject = c;
                                ri.predicate = "has_genetic_modifier";
                                ri.object = geneticmodifier; // FIXME: we need to create a new copy of the modifier Characteristic!!
                                ri.summary = "Genetic modifier [multigene]";
                                results.add( ri );
                            }
                            log.info( "Resolve multi-gene genotype" );
                            solved = true;
                        }

                        // look for fusion_gene and two genes
                        Characteristic cf = isGeneFusion( fv );
                        if ( !solved && fv.getCharacteristics().size() == 3 && cf != null ) {
                            boolean firstGene = true;
                            ri.summary = "Gene fusion";
                            ri.predicate = "fusion_gene";
                            Characteristic[] characteristics = fv.getCharacteristics().toArray( new Characteristic[0] );
                            for ( int i = 0; i < characteristics.length; i++ ) {
                                Characteristic c = characteristics[i];
                                if ( c.equals( cf ) ) continue;
                                if ( isGene( c ) ) {
                                    if ( ri.object == null )
                                        ri.object = c; // arbitrarily pick one gene to be the object.
                                    else ri.subject = c;
                                }
                            }

                            if ( ri.subject == null || ri.object == null || ri.predicate == null ) {
                                throw new IllegalStateException( "Didn't find all three parts of the fusion for: " + ee + " " + fv );
                            }

                            log.info( "Resolved gene fusion" );
                            solved = true;
                        }
                        //nosolution++; let's not count these yet

                        //solved = false;


                    } // end inspections of the characteristics


                    if ( !solved ) {
                        System.out.println( "--- Unsolved --" );
                        System.out.println( ri );
                        unresolved.add( ri );
                    } else {
                        results.add( ri );
                    }
                }
            }


        } catch ( Exception e ) {
            log.error( e, e );
            this.addErrorObject( ee + ": " + e.getMessage(), "" );
        }
    }

    private static boolean isTimepoint( Characteristic cs1 ) {
        return cs1.getCategory() != null && cs1.getCategory().equals( "timepoint" );
    }

    private static boolean isGene( Characteristic c ) {
        return c.getValueUri() != null && c.getValueUri().contains( "ncbi_gene" );
    }

    private Characteristic isGeneFusion( FactorValue fv ) {
        for ( Characteristic c : fv.getCharacteristics() ) {
            if ( c.getValueUri() != null && c.getValueUri().equals( "http://purl.obolibrary.org/obo/SO_0000287" ) ) {
                return c;
            }
        }
        return null;
    }


    private boolean isDrug( Characteristic cs ) {
        return cs.getValueUri() != null && cs.getValueUri().contains( "CHEBI_" );
    }


    private boolean isDisease( Characteristic cs ) {
        return ( cs.getValueUri() != null && cs.getValueUri().contains( "MONDO_" ) ) || cs.getCategory().toLowerCase().contains( "disease" ); /* includes disease model */
    }

    private boolean isDelivery( Characteristic cs ) {
        return ( cs.getValueUri() != null && addingMaterialentityToTargetTerms.contains( cs.getValueUri() ) );
    }


    private boolean isQualityProperty( Characteristic cs ) {
        if ( cs.getValueUri() == null ) return false;
        return qualityTerms.contains( cs.getValueUri() ) || deviationFromNormalTerms.contains( cs.getValueUri() );
    }

    // e.g. children of http://purl.obolibrary.org/obo/MONDO_0021125 - disease characteristic
    private boolean isDiseaseModifier( Characteristic cs ) {

        if ( cs.getValueUri() == null ) return false;
        return diseaseModifierTerms.contains( cs.getValueUri() );
    }


    private boolean isDevelopmentalStage( Characteristic cs ) {
        return stageTerms.contains( cs.getValueUri() );
    }

    private boolean isLocation( Characteristic cs ) {
        if ( cs.getValueUri() == null ) return false;
        return locationTerms.contains( cs.getValueUri() ) || cs.getValueUri().contains( "PATO_" ) /*FIXME DON"T DO THIS */;

    }

    private boolean isControlCondition( Characteristic cs ) {
        return BaselineSelection.isBaselineCondition( cs );
    }

    private boolean isGeneticManipulation( Characteristic cs ) {
        // e.g. knockdown, knockout, overexpression
        return cs.getValueUri() != null && geneticTerms.contains( cs.getValueUri() );
    }

    private final Set<String> geneticTerms = new HashSet<>( Arrays.asList(
            /*Dominant negative mutation */ "http://gemma.msl.ubc.ca/ont/TGEMO_00009",
            /*Constitutive active mutation */ "http://gemma.msl.ubc.ca/ont/TGEMO_00008",
            /*Knockdown */ "http://gemma.msl.ubc.ca/ont/TGEMO_00007",
            /*Overexpression */ "http://gemma.msl.ubc.ca/ont/TGEMO_00004",
            /*Double-copy overexpression */ "http://gemma.msl.ubc.ca/ont/TGEMO_00006",
            /*Single-copy overexpression */ "http://gemma.msl.ubc.ca/ont/TGEMO_00005",
            /*Heterozygous */ "http://gemma.msl.ubc.ca/ont/TGEMO_00003",
            /*Homozygous negative */ "http://gemma.msl.ubc.ca/ont/TGEMO_00001",
            /*TGEMO Genetics/Mutations */ "http://gemma.msl.ubc.ca/ont/TGEMO_00002",
            /*Rescue by external protein */ "http://gemma.msl.ubc.ca/ont/TGEMO_00102"
    ) );

    private final Set<String> stageTerms = new HashSet<>( Arrays.asList(
            /*  http://purl.obolibrary.org/obo/UBERON_0000105 (life cycle stage)
                and http://www.ebi.ac.uk/efo/EFO_0000399 (developmental stage) children */
            /* embryo stage */ "http://purl.obolibrary.org/obo/UBERON_0000068",
            /*neonate */ "http://purl.obolibrary.org/obo/UBERON_0007221",
            /*infant*/ "http://www.ebi.ac.uk/efo/EFO_0001355",/* fruit ripening stage */ "http://purl.obolibrary.org/obo/PO_0025502",
            /* fertilized egg stage */ "http://www.ebi.ac.uk/efo/EFO_0001322",
            /* postpartum */ "http://www.ebi.ac.uk/efo/EFO_0008562",
            /* coleorhiza emergence stage */ "http://purl.obolibrary.org/obo/PO_0025475",
            /* floral transition */ "http://www.ebi.ac.uk/efo/EFO_0002683",
            /* cleavage 64-cell */ "http://www.ebi.ac.uk/efo/EFO_0001288",
            /* animal zygote */ "http://purl.obolibrary.org/obo/CL_0000365",
            /* gastrula germ-ring */ "http://www.ebi.ac.uk/efo/EFO_0001295",
            /* radicle emergence */ "http://purl.obolibrary.org/obo/PO_0007015",
            /* embryonic stage 1 */ "http://www.ebi.ac.uk/efo/EFO_0005860",
            /* embryonic stage 14 */ "http://www.ebi.ac.uk/efo/EFO_0005874",
            /* mouse embryo stage */ "http://www.ebi.ac.uk/efo/EFO_0005857",
            /* LP.08 eight leaves visible stage */ "http://purl.obolibrary.org/obo/PO_0007095",
            /* third instar larva stage */ "http://www.ebi.ac.uk/efo/EFO_0002682",
            /* embryonic day 10.5 */ "http://www.ebi.ac.uk/efo/EFO_0007643",
            /* coleoptile emergence stage */ "http://purl.obolibrary.org/obo/PO_0007045",
            /* embryonic day 9 */ "http://www.ebi.ac.uk/efo/EFO_0007644",
            /* blastula oblong */ "http://www.ebi.ac.uk/efo/EFO_0001281",
            /* embryonic day 12 */ "http://www.ebi.ac.uk/efo/EFO_0007640",
            /* LP.13 thirteen leaves visible stage */ "http://purl.obolibrary.org/obo/PO_0007083",
            /* segmentation stage */ "http://www.ebi.ac.uk/efo/EFO_0001315",
            /* FL.00 first flower(s) open stage */ "http://purl.obolibrary.org/obo/PO_0007026",
            /* 3 inflorescence detectable stage */ "http://purl.obolibrary.org/obo/PO_0007047",
            /* Platyhelminthes life stage */ "http://www.ebi.ac.uk/efo/EFO_0007711",
            /* embryonic stage 2 */ "http://www.ebi.ac.uk/efo/EFO_0005861",
            /* embryonic stage 13 */ "http://www.ebi.ac.uk/efo/EFO_0005873",
            /* LP.04 four leaves visible stage */ "http://purl.obolibrary.org/obo/PO_0007115",
            /* fertilized ovule stage */ "http://purl.obolibrary.org/obo/PO_0004505",
            /* blastula 256-cell */ "http://www.ebi.ac.uk/efo/EFO_0001276",
            /* pharyngula prim-5 */ "http://www.ebi.ac.uk/efo/EFO_0001314",
            /* seed maturation stage */ "http://purl.obolibrary.org/obo/PO_0007632",
            /* late rosette growth stage */ "http://purl.obolibrary.org/obo/PO_0007076",
            /* rosette growth stage */ "http://purl.obolibrary.org/obo/PO_0007113",
            /* SE.00 stem elongation begins stage */ "http://purl.obolibrary.org/obo/PO_0007079",
            /* embryonic day 14.5 */ "http://www.ebi.ac.uk/efo/EFO_0002565",
            /* post dauer stage */ "http://www.ebi.ac.uk/efo/EFO_0005511",
            /* Theiler stage 27 */ "http://www.ebi.ac.uk/efo/EFO_0004391",
            /* Theiler stage 7 */ "http://www.ebi.ac.uk/efo/EFO_0004399",
            /* drosophila developmental stage */ "http://www.ebi.ac.uk/efo/EFO_0005651",
            /* pupal stage */ "http://purl.obolibrary.org/obo/UBERON_0000070",
            /* segmentation 1-4 somites */ "http://www.ebi.ac.uk/efo/EFO_0001316",
            /* LP.12 twelve leaves visible stage */ "http://purl.obolibrary.org/obo/PO_0007064",
            /* whole plant fruit ripening stage */ "http://purl.obolibrary.org/obo/PO_0007010",
            /* LP.20 twenty or more leaves whorls visible stage */ "http://purl.obolibrary.org/obo/PO_0007082",
            /* L2 larva */ "http://www.ebi.ac.uk/efo/EFO_0002719",
            /* 2 formation of axillary shoot stage */ "http://purl.obolibrary.org/obo/PO_0007073",
            /* beginning of whole plant fruit ripening stage */ "http://purl.obolibrary.org/obo/PO_0007036",
            /* cleavage stage */ "http://purl.obolibrary.org/obo/UBERON_0000107",
            /* embryonic stage 16 */ "http://www.ebi.ac.uk/efo/EFO_0005876",
            /* pharyngula stage */ "http://purl.obolibrary.org/obo/UBERON_0004707",
            /* SE.99 maximum stem length reached stage */ "http://purl.obolibrary.org/obo/PO_0007109",
            /* L4 larva */ "http://www.ebi.ac.uk/efo/EFO_0005509",
            /* cleavage 16-cell */ "http://www.ebi.ac.uk/efo/EFO_0001284",
            /* Theiler stage 20 */ "http://www.ebi.ac.uk/efo/EFO_0004410",
            /* Theiler stage 21 */ "http://www.ebi.ac.uk/efo/EFO_0002584",
            /* pharyngula prim-15 */ "http://www.ebi.ac.uk/efo/EFO_0001312",
            /* gastrula 50%-epiboly */ "http://www.ebi.ac.uk/efo/EFO_0001291",
            /* flower development stage */ "http://purl.obolibrary.org/obo/PO_0007615",
            /* dorsal closure stage */ "http://purl.obolibrary.org/obo/FBdv_00005331",
            /* FL.02 1/2 of flowers open stage */ "http://purl.obolibrary.org/obo/PO_0007053",
            /* Theiler stage 10 */ "http://www.ebi.ac.uk/efo/EFO_0004402",
            /* mouse postnatal */ "http://www.ebi.ac.uk/efo/EFO_0004390",
            /* embryonic stage 15 */ "http://www.ebi.ac.uk/efo/EFO_0005875",
            /* LP.06 six leaves visible stage */ "http://purl.obolibrary.org/obo/PO_0007123",
            /* cleavage 32-cell */ "http://www.ebi.ac.uk/efo/EFO_0001286",
            /* larval day 21-29 */ "http://www.ebi.ac.uk/efo/EFO_0001308",
            /* sporophyte development stage */ "http://purl.obolibrary.org/obo/PO_0028002",
            /* early rosette growth stage */ "http://purl.obolibrary.org/obo/PO_0007081",
            /* Theiler stage 22 */ "http://www.ebi.ac.uk/efo/EFO_0002585",
            /* Theiler stage 5 */ "http://www.ebi.ac.uk/efo/EFO_0004397",
            /* ring stage trophozoite */ "http://www.ebi.ac.uk/efo/EFO_0002590",
            /* blastula dome */ "http://www.ebi.ac.uk/efo/EFO_0001279",
            /* FL.04 end of flowering stage */ "http://purl.obolibrary.org/obo/PO_0007024",
            /* juvenile days 45-89 */ "http://www.ebi.ac.uk/efo/EFO_0001302",
            /* embryonic stage 5 */ "http://www.ebi.ac.uk/efo/EFO_0005864",
            /* hatching pec-fin */ "http://www.ebi.ac.uk/efo/EFO_0001299",
            /* embryonic stage 10 */ "http://www.ebi.ac.uk/efo/EFO_0005870",
            /* sexually immature stage */ "http://purl.obolibrary.org/obo/UBERON_0000112",
            /* cleavage 2-cell */ "http://www.ebi.ac.uk/efo/EFO_0001285",
            /* seedling development stage */ "http://purl.obolibrary.org/obo/PO_0007131",
            /* whole plant fruit formation stage 70% to final size */ "http://purl.obolibrary.org/obo/PO_0007027",
            /* Theiler stage 19 */ "http://www.ebi.ac.uk/efo/EFO_0004409",
            /* plant embryo development stage */ "http://purl.obolibrary.org/obo/PO_0007631",
            /* Theiler stage 6 */ "http://www.ebi.ac.uk/efo/EFO_0004398",
            /* Caenorhabditis elegans larval stage */ "http://www.ebi.ac.uk/efo/EFO_0007694",
            /* floral organ formation stage */ "http://purl.obolibrary.org/obo/PO_0025585",
            /* LP.18 eighteen leaves visible stage */ "http://purl.obolibrary.org/obo/PO_0007072",
            /* juvenile days 30-44 */ "http://www.ebi.ac.uk/efo/EFO_0001301",
            /* embryonic stage 6 */ "http://www.ebi.ac.uk/efo/EFO_0005865",
            /* gynoecium development stage */ "http://purl.obolibrary.org/obo/PO_0007606",
            /* whole plant fruit ripening complete stage */ "http://purl.obolibrary.org/obo/PO_0007038",
            /* 0 seed germination stage */ "http://purl.obolibrary.org/obo/PO_0007057",
            /* 2.03 main shoot and axillary shoots visible at three nodes stage */ "http://purl.obolibrary.org/obo/PO_0007080",
            /* sepal primordium visible stage */ "http://purl.obolibrary.org/obo/PO_0007607",
            /* elongating embryo Ce */ "http://www.ebi.ac.uk/efo/EFO_0005506",
            /* inflorescence emergence stage */ "http://purl.obolibrary.org/obo/PO_0007041",
            /* Theiler stage 3 */ "http://www.ebi.ac.uk/efo/EFO_0004395",
            /* fruit size up to 10% stage */ "http://purl.obolibrary.org/obo/PO_0025504",
            /* embryonic day 18.5 */ "http://www.ebi.ac.uk/efo/EFO_0002570",
            /* embryonic day 13 */ "http://www.ebi.ac.uk/efo/EFO_0007642",
            /* 3-fold embryo Ce */ "http://www.ebi.ac.uk/efo/EFO_0005499",
            /* cleavage 4-cell */ "http://www.ebi.ac.uk/efo/EFO_0001287",
            /* neonate */ "http://www.ebi.ac.uk/efo/EFO_0001372",
            /* sporophyte senescent stage */ "http://purl.obolibrary.org/obo/PO_0007017",
            /* newly molted young adult hermaphrodite */ "http://www.ebi.ac.uk/efo/EFO_0005498",
            /* petal differentiation and expansion stage */ "http://purl.obolibrary.org/obo/PO_0007611",
            /* embryonic stage 3 */ "http://www.ebi.ac.uk/efo/EFO_0005862",
            /* embryonic stage 12 */ "http://www.ebi.ac.uk/efo/EFO_0005872",
            /* LP.01 one leaf visible stage */ "http://purl.obolibrary.org/obo/PO_0007094",
            /* gastrulating embryo Ce */ "http://www.ebi.ac.uk/efo/EFO_0005503",
            /* infant */ "http://www.ebi.ac.uk/efo/EFO_0001355",
            /* LP.07 seven leaves visible stage */ "http://purl.obolibrary.org/obo/PO_0007063",
            /* embryonic day 13.5 */ "http://www.ebi.ac.uk/efo/EFO_0002564",
            /* stamen primordium visible stage */ "http://purl.obolibrary.org/obo/PO_0007613",
            /* seed imbibition stage */ "http://purl.obolibrary.org/obo/PO_0007022",
            /* Theiler stage 4 */ "http://www.ebi.ac.uk/efo/EFO_0004396",
            /* LP.17 seventeen leaves visible stage */ "http://purl.obolibrary.org/obo/PO_0007067",
            /* embryonic day 8.25 */ "http://www.ebi.ac.uk/efo/EFO_0002561",
            /* 1 main shoot growth stage */ "http://purl.obolibrary.org/obo/PO_0007112",
            /* mid rosette growth stage */ "http://purl.obolibrary.org/obo/PO_0007068",
            /* dauer larva */ "http://www.ebi.ac.uk/efo/EFO_0005507",
            /* gastrula bud */ "http://www.ebi.ac.uk/efo/EFO_0001294",
            /* calyx development stage */ "http://purl.obolibrary.org/obo/PO_0007603",
            /* embryonic stage 4 */ "http://www.ebi.ac.uk/efo/EFO_0005863",
            /* embryonic stage 11 */ "http://www.ebi.ac.uk/efo/EFO_0005871",
            /* SE.02 two nodes or internodes visible stage */ "http://purl.obolibrary.org/obo/PO_0007117",
            /* rosette growth complete stage */ "http://purl.obolibrary.org/obo/PO_0007078",
            /* adult */ "http://www.ebi.ac.uk/efo/EFO_0001272",
            /* Theiler stage 1 */ "http://www.ebi.ac.uk/efo/EFO_0004393",
            /* larval day 4 */ "http://www.ebi.ac.uk/efo/EFO_0001304",
            /* plasmodium parasite stage */ "http://www.ebi.ac.uk/efo/EFO_0002544",
            /* blastula high */ "http://www.ebi.ac.uk/efo/EFO_0001280",
            /* embryo stage */ "http://www.ebi.ac.uk/efo/EFO_0007725",
            /* whole plant fruit formation stage 30 to 50% */ "http://purl.obolibrary.org/obo/PO_0007029",
            /* larval stage */ "http://purl.obolibrary.org/obo/UBERON_0000069",
            /* Theiler stage 15 */ "http://www.ebi.ac.uk/efo/EFO_0004406",
            /* embryonic stage 9 */ "http://www.ebi.ac.uk/efo/EFO_0005869",
            /* embryonic day 11.5 */ "http://www.ebi.ac.uk/efo/EFO_0002562",
            /* mid whole plant fruit ripening stage */ "http://purl.obolibrary.org/obo/PO_0007031",
            /* fully-elongated embryo Ce */ "http://www.ebi.ac.uk/efo/EFO_0005504",
            /* Theiler stage 2 */ "http://www.ebi.ac.uk/efo/EFO_0004394",
            /* mouse prenatal */ "http://www.ebi.ac.uk/efo/EFO_0002543",
            /* Theiler stage 17 */ "http://www.ebi.ac.uk/efo/EFO_0002583",
            /* segmentation 20-25 somites */ "http://www.ebi.ac.uk/efo/EFO_0001320",
            /* cercarium */ "http://www.ebi.ac.uk/efo/EFO_0007712",
            /* late cleavage stage embryo Ce */ "http://www.ebi.ac.uk/efo/EFO_0005502",
            /* juvenile stage */ "http://purl.obolibrary.org/obo/UBERON_0034919",
            /* cotyledon emergence stage */ "http://purl.obolibrary.org/obo/PO_0007049",
            /* Theiler stage 16 */ "http://www.ebi.ac.uk/efo/EFO_0004407",
            /* LP.11 eleven leaves visible stage */ "http://purl.obolibrary.org/obo/PO_0007116",
            /* Theiler stage 28 */ "http://www.ebi.ac.uk/efo/EFO_0002588",
            /* developing seed stage */ "http://purl.obolibrary.org/obo/PO_0004506",
            /* embryo Ce */ "http://purl.obolibrary.org/obo/WBls_0000003",
            /* LP.02 two leaves visible stage */ "http://purl.obolibrary.org/obo/PO_0007098",
            /* segmentation 14-19 somites */ "http://www.ebi.ac.uk/efo/EFO_0001319",
            /* 3 hr schistosomulum */ "http://www.ebi.ac.uk/efo/EFO_0007713",
            /* plant structure development stage */ "http://purl.obolibrary.org/obo/PO_0009012",
            /* L3 larva */ "http://www.ebi.ac.uk/efo/EFO_0002720",
            /* embryonic day 9.5 */ "http://www.ebi.ac.uk/efo/EFO_0007641",
            /* cleavage 8-cell */ "http://www.ebi.ac.uk/efo/EFO_0001289",
            /* erythrocytic schizont */ "http://www.ebi.ac.uk/efo/EFO_0002591",
            /* embryonic stage 7 */ "http://www.ebi.ac.uk/efo/EFO_0005866",
            /* L1 larva */ "http://www.ebi.ac.uk/efo/EFO_0005508",
            /* LP.14 fourteen leaves visible stage */ "http://purl.obolibrary.org/obo/PO_0007085",
            /* fruit size 30 to 50% stage */ "http://purl.obolibrary.org/obo/PO_0025506",
            /* androecium development stage */ "http://purl.obolibrary.org/obo/PO_0007605",
            /* larval protruding mouth */ "http://www.ebi.ac.uk/efo/EFO_0007696",
            /* Platyhelminthes adult */ "http://www.ebi.ac.uk/efo/EFO_0007715",
            /* enclosing embryo Ce */ "http://www.ebi.ac.uk/efo/EFO_0005505",
            /* leaf production stage */ "http://purl.obolibrary.org/obo/PO_0007133",
            /* embryonic day 17.5 */ "http://www.ebi.ac.uk/efo/EFO_0002568",
            /* hatching stage */ "http://www.ebi.ac.uk/efo/EFO_0001298",
            /* IL.00 inflorescence just visible stage */ "http://purl.obolibrary.org/obo/PO_0007006",
            /* proliferating embryo Ce */ "http://www.ebi.ac.uk/efo/EFO_0005501",
            /* inflorescence development stage */ "http://purl.obolibrary.org/obo/PO_0001083",
            /* floral organ differentiation stage */ "http://purl.obolibrary.org/obo/PO_0007600",
            /* blastula stage */ "http://purl.obolibrary.org/obo/UBERON_0000108",
            /* 4-cell embryo Ce */ "http://www.ebi.ac.uk/efo/EFO_0005500",
            /* embryonic stage 8 */ "http://www.ebi.ac.uk/efo/EFO_0005868",
            /* blastula 128-cell */ "http://www.ebi.ac.uk/efo/EFO_0001273",
            /* booting stage */ "http://purl.obolibrary.org/obo/PO_0007014",
            /* corolla development stage */ "http://purl.obolibrary.org/obo/PO_0007604",
            /* embryonic day 12.5 */ "http://www.ebi.ac.uk/efo/EFO_0002563",
            /* stem elongation stage */ "http://purl.obolibrary.org/obo/PO_0007089",
            /* Theiler stage 18 */ "http://www.ebi.ac.uk/efo/EFO_0004408",
            /* LP.10 ten leaves visible stage */ "http://purl.obolibrary.org/obo/PO_0007103",
            /* Theiler stage 8 */ "http://www.ebi.ac.uk/efo/EFO_0004400",
            /* FL.01 1/4 of flowers open stage */ "http://www.ebi.ac.uk/efo/EFO_0005793",
            /* larval day 7-13 */ "http://www.ebi.ac.uk/efo/EFO_0001309",
            /* embryonic day 18 */ "http://www.ebi.ac.uk/efo/EFO_0002569",
            /* LP.09 nine leaves visible */ "http://purl.obolibrary.org/obo/PO_0007101",
            /* fruit development stage */ "http://purl.obolibrary.org/obo/PO_0001002",
            /* hepatic schizont */ "http://www.ebi.ac.uk/efo/EFO_0002589",
            /* blastula 1k-cell */ "http://www.ebi.ac.uk/efo/EFO_0001274",
            /* LP.03 three leaves visible stage */ "http://purl.obolibrary.org/obo/PO_0007106",
            /* fruit formation stage */ "http://purl.obolibrary.org/obo/PO_0025501",
            /* late embryonic stage */ "http://purl.obolibrary.org/obo/UBERON_0007220",
            /* early whole plant fruit ripening stage */ "http://purl.obolibrary.org/obo/PO_0007001",
            /* Theiler stage 11 */ "http://www.ebi.ac.uk/efo/EFO_0002582",
            /* Drosophila embryo stage */ "http://www.ebi.ac.uk/efo/EFO_0005859",
            /* pharyngula prim-25 */ "http://www.ebi.ac.uk/efo/EFO_0001313",
            /* embryonic day 15.5 */ "http://www.ebi.ac.uk/efo/EFO_0002566",
            /* Danio rerio larval stage */ "http://www.ebi.ac.uk/efo/EFO_0007695",
            /* Theiler stage 9 */ "http://www.ebi.ac.uk/efo/EFO_0004401",
            /* fruit size 50 to 70% stage */ "http://purl.obolibrary.org/obo/PO_0025507",
            /* floral organ meristem development stage */ "http://purl.obolibrary.org/obo/PO_0007601",
            /* blastula sphere */ "http://www.ebi.ac.uk/efo/EFO_0001283",
            /* LP.19 nineteen leaves visible stage */ "http://purl.obolibrary.org/obo/PO_0007120",
            /* LP.16 sixteen leaves visible stage */ "http://purl.obolibrary.org/obo/PO_0007119",
            /* gastrula stage */ "http://purl.obolibrary.org/obo/UBERON_0000109",
            /* sporophyte vegetative stage */ "http://purl.obolibrary.org/obo/PO_0007134",
            /* larval day 14-20 */ "http://www.ebi.ac.uk/efo/EFO_0001307",
            /* Theiler stage 12 */ "http://www.ebi.ac.uk/efo/EFO_0004403",
            /* postnatal */ "http://www.ebi.ac.uk/efo/EFO_0002948",
            /* L2d-dauer molt */ "http://www.ebi.ac.uk/efo/EFO_0005510",
            /* embryonic stage 17 */ "http://www.ebi.ac.uk/efo/EFO_0005877",
            /* 24 hr schistosomulum */ "http://www.ebi.ac.uk/efo/EFO_0007714",
            /* LP.05 five leaves visible stage */ "http://purl.obolibrary.org/obo/PO_0007065",
            /* gastrula shield */ "http://www.ebi.ac.uk/efo/EFO_0001297",
            /* Theiler stage 24 */ "http://www.ebi.ac.uk/efo/EFO_0002586",
            /* Theiler stage 23 */ "http://www.ebi.ac.uk/efo/EFO_0004411",
            /* LP.15 fifteen leaves visible stage */ "http://purl.obolibrary.org/obo/PO_0007104",
            /* larval day 6 */ "http://www.ebi.ac.uk/efo/EFO_0001306",
            /* hatching long-pec */ "http://www.ebi.ac.uk/efo/EFO_0002718",
            /* gastrula 75%-epiboly */ "http://www.ebi.ac.uk/efo/EFO_0001292",
            /* pharyngula high-pec */ "http://www.ebi.ac.uk/efo/EFO_0001311",
            /* trophozoite */ "http://www.ebi.ac.uk/efo/EFO_0002592",
            /* Theiler stage 13 */ "http://www.ebi.ac.uk/efo/EFO_0004404",
            /* segmentation 26+ somites */ "http://www.ebi.ac.uk/efo/EFO_0001321",
            /* fruit size 70% to final size stage */ "http://purl.obolibrary.org/obo/PO_0025508",
            /* whole plant fruit formation stage */ "http://purl.obolibrary.org/obo/PO_0007042",
            /* segmentation 5-9 somites */ "http://www.ebi.ac.uk/efo/EFO_0001317",
            /* blastula 30%-epiboly */ "http://www.ebi.ac.uk/efo/EFO_0001277",
            /* gastrula 90%-epiboly */ "http://www.ebi.ac.uk/efo/EFO_0001293",
            /* segmentation 10-13 somites */ "http://www.ebi.ac.uk/efo/EFO_0001318",
            /* postmenopausal */ "http://www.ebi.ac.uk/efo/EFO_0002721",
            /* larval day 5 */ "http://www.ebi.ac.uk/efo/EFO_0001305",
            /* hypocotyl emergence stage */ "http://purl.obolibrary.org/obo/PO_0007043",
            /* whole plant flowering stage */ "http://purl.obolibrary.org/obo/PO_0007016",
            /* C. elegans embryo stage */ "http://www.ebi.ac.uk/efo/EFO_0005858",
            /* Theiler stage 14 */ "http://www.ebi.ac.uk/efo/EFO_0004405",
            /* gastrula 80%-epiboly */ "http://www.ebi.ac.uk/efo/EFO_0002560",
            /* blastula 512-cell */ "http://www.ebi.ac.uk/efo/EFO_0001278",
            /* Theiler stage 26 */ "http://www.ebi.ac.uk/efo/EFO_0002587",
            /* Theiler stage 25 */ "http://www.ebi.ac.uk/efo/EFO_0004412",
            /* seed development stage */ "http://purl.obolibrary.org/obo/PO_0001170",
            /* embryonic day 16.5 */ "http://www.ebi.ac.uk/efo/EFO_0002567" ) );


    private final Set<String> qualityTerms = new HashSet<>( Arrays.asList( /* maybe this is too broad */
            /* occurrence */ "http://purl.obolibrary.org/obo/PATO_0000057", // occurrence terms
            /* decreased occurrence */ "http://purl.obolibrary.org/obo/PATO_0002052",
            /* arrested */ "http://purl.obolibrary.org/obo/PATO_0000297",
            /* asynchronous */ "http://purl.obolibrary.org/obo/PATO_0000688",
            /* non-progressive */ "http://purl.obolibrary.org/obo/PATO_0002026",
            /* progressive */ "http://purl.obolibrary.org/obo/PATO_0001818",
            /* normal occurrence */ "http://purl.obolibrary.org/obo/PATO_0045086",
            /* recurrent */ "http://purl.obolibrary.org/obo/PATO_0000427",
            /* refractory */ "http://purl.obolibrary.org/obo/PATO_0002631",
            /* repetitive */ "http://purl.obolibrary.org/obo/PATO_0000441",
            /* alternation */ "http://purl.obolibrary.org/obo/PATO_0000999",
            /* increased occurrence */ "http://purl.obolibrary.org/obo/PATO_0002051",
            /* continuous */ "http://purl.obolibrary.org/obo/PATO_0000689",
            /* episodic */ "http://purl.obolibrary.org/obo/PATO_0002630",
            /* synchronous */ "http://purl.obolibrary.org/obo/PATO_0000695",
            /* sporadic */ "http://purl.obolibrary.org/obo/PATO_0000428",
            /* discontinuous */ "http://purl.obolibrary.org/obo/PATO_0000690",
            /* heterochronic */ "http://purl.obolibrary.org/obo/PATO_0000692",
            /* disrupted */ "http://purl.obolibrary.org/obo/PATO_0001507",
            /* adjacent tissue */ "http://gemma.msl.ubc.ca/ont/TGEMO_00058", // from TGEMO
            /* physical object quality */ "http://purl.obolibrary.org/obo/PATO_0001241", // physical object quality eterms
            /* has number of */ "http://purl.obolibrary.org/obo/PATO_0001555",
            /* altered number of */ "http://purl.obolibrary.org/obo/PATO_0002083",
            /* has extra parts of type */ "http://purl.obolibrary.org/obo/PATO_0002002",
            /* lacks parts or has fewer parts of type */ "http://purl.obolibrary.org/obo/PATO_0001999",
            /* has fewer parts of type */ "http://purl.obolibrary.org/obo/PATO_0002001",
            /* lacks all parts of type */ "http://purl.obolibrary.org/obo/PATO_0002000",
            /* has normal numbers of parts of type */ "http://purl.obolibrary.org/obo/PATO_0001905",
            /* structurally discontinuous */ "http://purl.obolibrary.org/obo/PATO_0040026",
            /* organismal quality */ "http://purl.obolibrary.org/obo/PATO_0001995",
            /* arboreal */ "http://purl.obolibrary.org/obo/PATO_0095004",
            /* pathogenicity */ "http://purl.obolibrary.org/obo/PATO_0040003",
            /* viability */ "http://purl.obolibrary.org/obo/PATO_0000169",
            /* dead */ "http://purl.obolibrary.org/obo/PATO_0001422",
            /* semi-lethal (sensu genetics) */ "http://purl.obolibrary.org/obo/PATO_0001768",
            /* brood viability */ "http://purl.obolibrary.org/obo/PATO_0001497",
            /* alive */ "http://purl.obolibrary.org/obo/PATO_0001421",
            /* viable */ "http://purl.obolibrary.org/obo/PATO_0000719",
            /* lethal (sensu genetics) */ "http://purl.obolibrary.org/obo/PATO_0000718",
            /* decayed */ "http://purl.obolibrary.org/obo/PATO_0001432",
            /* semi-viable */ "http://purl.obolibrary.org/obo/PATO_0001770",
            /* immortal */ "http://purl.obolibrary.org/obo/PATO_0001991",
            /* sexually dimorphic */ "http://purl.obolibrary.org/obo/PATO_0002451",
            /* sex-specific */ "http://purl.obolibrary.org/obo/PATO_0060001",
            /* female-specific */ "http://purl.obolibrary.org/obo/PATO_0060003",
            /* male-specific */ "http://purl.obolibrary.org/obo/PATO_0060002",
            /* cellularity */ "http://purl.obolibrary.org/obo/PATO_0001992",
            /* unicellular */ "http://purl.obolibrary.org/obo/PATO_0001994",
            /* multicellular */ "http://purl.obolibrary.org/obo/PATO_0001993",
            /* biological sex */ "http://purl.obolibrary.org/obo/PATO_0000047",
            /* genotypic sex */ "http://purl.obolibrary.org/obo/PATO_0020000",
            /* female genotypic sex */ "http://purl.obolibrary.org/obo/PATO_0020002",
            /* male genotypic sex */ "http://purl.obolibrary.org/obo/PATO_0020001",
            /* mating type */ "http://purl.obolibrary.org/obo/PATO_0001895",
            /* yeast mating type */ "http://purl.obolibrary.org/obo/PATO_0001337",
            /* Schizosaccharomyces pombe mating type */ "http://purl.obolibrary.org/obo/PATO_0001343",
            /* h plus */ "http://purl.obolibrary.org/obo/PATO_0001346",
            /* h minus */ "http://purl.obolibrary.org/obo/PATO_0001345",
            /* Saccharomyces cerevisiae mating type */ "http://purl.obolibrary.org/obo/PATO_0001342",
            /* alpha mating type (yeast) */ "http://purl.obolibrary.org/obo/PATO_0001344",
            /* a mating type (yeast) */ "http://purl.obolibrary.org/obo/PATO_0001341",
            /* bacterial mating type */ "http://purl.obolibrary.org/obo/PATO_0001335",
            /* F minus mating type */ "http://purl.obolibrary.org/obo/PATO_0001348",
            /* high frequency recombinant */ "http://purl.obolibrary.org/obo/PATO_0001349",
            /* F mating type */ "http://purl.obolibrary.org/obo/PATO_0001347",
            /* phenotypic sex */ "http://purl.obolibrary.org/obo/PATO_0001894",
            /* female with DSD */ "http://purl.obolibrary.org/obo/PATO_0040056",
            /* monandrous */ "http://purl.obolibrary.org/obo/PATO_0040051",
            /* diandrous */ "http://purl.obolibrary.org/obo/PATO_0040050",
            /* male with DSD */ "http://purl.obolibrary.org/obo/PATO_0040049",
            /* pseudohermaphrodite */ "http://purl.obolibrary.org/obo/PATO_0001827",
            /* female pseudohermaphrodite */ "http://purl.obolibrary.org/obo/PATO_0001829",
            /* male pseudohermaphrodite */ "http://purl.obolibrary.org/obo/PATO_0001828",
            /* hermaphrodite */ "http://purl.obolibrary.org/obo/PATO_0001340",
            /* sequential hermaphrodite */ "http://purl.obolibrary.org/obo/PATO_0070000",
            /* protandrous */ "http://purl.obolibrary.org/obo/PATO_0040053",
            /* protandrous hermaphroditism */ "http://purl.obolibrary.org/obo/PATO_0040055",
            /* protogynous */ "http://purl.obolibrary.org/obo/PATO_0040052",
            /* protogynous hermaphroditism */ "http://purl.obolibrary.org/obo/PATO_0040054",
            /* male */ "http://purl.obolibrary.org/obo/PATO_0000384",
            /* castrated male */ "http://purl.obolibrary.org/obo/PATO_0002367",
            /* intact male */ "http://purl.obolibrary.org/obo/PATO_0002366",
            /* female */ "http://purl.obolibrary.org/obo/PATO_0000383",
            /* spayed female */ "http://purl.obolibrary.org/obo/PATO_0040020",
            /* intact female */ "http://purl.obolibrary.org/obo/PATO_0002365",
            /* organismal area to mass ratio */ "http://purl.obolibrary.org/obo/PATO_0095001",
            /* resistance to */ "http://purl.obolibrary.org/obo/PATO_0001046",
            /* resistant to */ "http://purl.obolibrary.org/obo/PATO_0001178",
            /* normal resistance to */ "http://purl.obolibrary.org/obo/PATO_0045047",
            /* decreased resistance to */ "http://purl.obolibrary.org/obo/PATO_0001651",
            /* increased resistance to */ "http://purl.obolibrary.org/obo/PATO_0001650",
            /* tolerant to */ "http://purl.obolibrary.org/obo/PATO_0000515",
            /* normal tolerance to */ "http://purl.obolibrary.org/obo/PATO_0045070",
            /* decreased tolerance to */ "http://purl.obolibrary.org/obo/PATO_0002394",
            /* increased tolerance to */ "http://purl.obolibrary.org/obo/PATO_0002393",
            /* anecic */ "http://purl.obolibrary.org/obo/PATO_0095009",
            /* homeostatic */ "http://purl.obolibrary.org/obo/PATO_0002541",
            /* shedability */ "http://purl.obolibrary.org/obo/PATO_0001729",
            /* herbaceous */ "http://purl.obolibrary.org/obo/PATO_0002352",
            /* semi-deciduous(plant) */ "http://purl.obolibrary.org/obo/PATO_0001734",
            /* evergreen (plant) */ "http://purl.obolibrary.org/obo/PATO_0001733",
            /* non-deciduous (any body part) */ "http://purl.obolibrary.org/obo/PATO_0001732",
            /* deciduous (plant) */ "http://purl.obolibrary.org/obo/PATO_0001731",
            /* deciduous (generic) */ "http://purl.obolibrary.org/obo/PATO_0001730",
            /* sensitivity toward */ "http://purl.obolibrary.org/obo/PATO_0000085",
            /* sensitivity to irradiation */ "http://purl.obolibrary.org/obo/PATO_0001806",
            /* normal sensitivity to irradiation */ "http://purl.obolibrary.org/obo/PATO_0045048",
            /* increased sensitivity to irradiation */ "http://purl.obolibrary.org/obo/PATO_0001808",
            /* decreased sensitivity to irradiation */ "http://purl.obolibrary.org/obo/PATO_0001807",
            /* sensitivity to oxygen */ "http://purl.obolibrary.org/obo/PATO_0001454",
            /* anaerobic */ "http://purl.obolibrary.org/obo/PATO_0001456",
            /* aerobic */ "http://purl.obolibrary.org/obo/PATO_0001455",
            /* susceptibility toward */ "http://purl.obolibrary.org/obo/PATO_0001043",
            /* insusceptible toward */ "http://purl.obolibrary.org/obo/PATO_0001153",
            /* susceptible toward */ "http://purl.obolibrary.org/obo/PATO_0001152",
            /* normal susceptibility toward */ "http://purl.obolibrary.org/obo/PATO_0045066",
            /* decreased susceptibility toward */ "http://purl.obolibrary.org/obo/PATO_0001670",
            /* increased susceptibility toward */ "http://purl.obolibrary.org/obo/PATO_0001669",
            /* photosensitivity */ "http://purl.obolibrary.org/obo/PATO_0000927",
            /* phototoxic */ "http://purl.obolibrary.org/obo/PATO_0001803",
            /* photosensitive */ "http://purl.obolibrary.org/obo/PATO_0000547",
            /* normal photosensitivity */ "http://purl.obolibrary.org/obo/PATO_0045038",
            /* increased photosensitivity */ "http://purl.obolibrary.org/obo/PATO_0001698",
            /* decreased photosensitivity */ "http://purl.obolibrary.org/obo/PATO_0001697",
            /* photoinsensitive */ "http://purl.obolibrary.org/obo/PATO_0000546",
            /* sensitive toward */ "http://purl.obolibrary.org/obo/PATO_0000516",
            /* normal sensitivity toward */ "http://purl.obolibrary.org/obo/PATO_0045049",
            /* decreased sensitivity toward */ "http://purl.obolibrary.org/obo/PATO_0001550",
            /* increased sensitivity toward */ "http://purl.obolibrary.org/obo/PATO_0001549",
            /* insensitive toward */ "http://purl.obolibrary.org/obo/PATO_0000513",
            /* threshold */ "http://purl.obolibrary.org/obo/PATO_0000152",
            /* normal threshold */ "http://purl.obolibrary.org/obo/PATO_0045069",
            /* decreased threshold */ "http://purl.obolibrary.org/obo/PATO_0000708",
            /* increased threshold */ "http://purl.obolibrary.org/obo/PATO_0000706",
            /* trophic quality */ "http://purl.obolibrary.org/obo/PATO_0000056",
            /* prototrophic */ "http://purl.obolibrary.org/obo/PATO_0000423",
            /* auxotrophic */ "http://purl.obolibrary.org/obo/PATO_0000422",
            /* fossorial */ "http://purl.obolibrary.org/obo/PATO_0095003",
            /* cursorial */ "http://purl.obolibrary.org/obo/PATO_0095002",
            /* maturity */ "http://purl.obolibrary.org/obo/PATO_0000261",
            /* prepupal */ "http://purl.obolibrary.org/obo/PATO_0001188",
            /* adolescent */ "http://purl.obolibrary.org/obo/PATO_0001189",
            /* juvenile */ "http://purl.obolibrary.org/obo/PATO_0001190",
            /* immature */ "http://purl.obolibrary.org/obo/PATO_0001501",
            /* prepubescent */ "http://purl.obolibrary.org/obo/PATO_0001186",
            /* larval */ "http://purl.obolibrary.org/obo/PATO_0001185",
            /* neonatal */ "http://purl.obolibrary.org/obo/PATO_0002206",
            /* pupal */ "http://purl.obolibrary.org/obo/PATO_0001187",
            /* pubescent */ "http://purl.obolibrary.org/obo/PATO_0000455",
            /* mature */ "http://purl.obolibrary.org/obo/PATO_0001701",
            /* response to */ "http://purl.obolibrary.org/obo/PATO_0000077",
            /* unresponsive to */ "http://purl.obolibrary.org/obo/PATO_0000488",
            /* responsive to */ "http://purl.obolibrary.org/obo/PATO_0000487",
            /* hyporesponsive to */ "http://purl.obolibrary.org/obo/PATO_0001194",
            /* hyperresponsive to */ "http://purl.obolibrary.org/obo/PATO_0001192",
            /* reproductive quality */ "http://purl.obolibrary.org/obo/PATO_0001434",
            /* parous */ "http://purl.obolibrary.org/obo/PATO_0040028",
            /* primiparous */ "http://purl.obolibrary.org/obo/PATO_0002371",
            /* multiparous */ "http://purl.obolibrary.org/obo/PATO_0002369",
            /* grand multiparous */ "http://purl.obolibrary.org/obo/PATO_0002372",
            /* parity */ "http://purl.obolibrary.org/obo/PATO_0002370",
            /* nulliparous */ "http://purl.obolibrary.org/obo/PATO_0002368",
            /* brood quality */ "http://purl.obolibrary.org/obo/PATO_0001496",
            /* brood size */ "http://purl.obolibrary.org/obo/PATO_0000276",
            /* fertility */ "http://purl.obolibrary.org/obo/PATO_0000274",
            /* lack of fertility in offspring */ "http://purl.obolibrary.org/obo/PATO_0001862",
            /* semi-fertile */ "http://purl.obolibrary.org/obo/PATO_0001767",
            /* male semi-fertile */ "http://purl.obolibrary.org/obo/PATO_0001761",
            /* female semi-fertile */ "http://purl.obolibrary.org/obo/PATO_0001760",
            /* sterile */ "http://purl.obolibrary.org/obo/PATO_0000956",
            /* female sterile */ "http://purl.obolibrary.org/obo/PATO_0000892",
            /* ovariectomized female */ "http://purl.obolibrary.org/obo/PATO_0002380",
            /* ovariohysterectomized female */ "http://purl.obolibrary.org/obo/PATO_0002379",
            /* male sterile */ "http://purl.obolibrary.org/obo/PATO_0000890",
            /* fertile */ "http://purl.obolibrary.org/obo/PATO_0000955",
            /* normal fertility */ "http://purl.obolibrary.org/obo/PATO_0045020",
            /* increased fertility */ "http://purl.obolibrary.org/obo/PATO_0001835",
            /* decreased fertility */ "http://purl.obolibrary.org/obo/PATO_0001834",
            /* male fertility */ "http://purl.obolibrary.org/obo/PATO_0000279",
            /* male semi-sterile */ "http://purl.obolibrary.org/obo/PATO_0001762",
            /* male fertile */ "http://purl.obolibrary.org/obo/PATO_0000891",
            /* normal male fertility */ "http://purl.obolibrary.org/obo/PATO_0045028",
            /* decreased male fertility */ "http://purl.obolibrary.org/obo/PATO_0001833",
            /* increased male fertility */ "http://purl.obolibrary.org/obo/PATO_0001832",
            /* female fertility */ "http://purl.obolibrary.org/obo/PATO_0000277",
            /* female semi-sterile */ "http://purl.obolibrary.org/obo/PATO_0001763",
            /* female fertile */ "http://purl.obolibrary.org/obo/PATO_0000888",
            /* normal female fertility */ "http://purl.obolibrary.org/obo/PATO_0045018",
            /* increased female fertility */ "http://purl.obolibrary.org/obo/PATO_0001831",
            /* decreased female fertility */ "http://purl.obolibrary.org/obo/PATO_0001830",
            /* fecundity */ "http://purl.obolibrary.org/obo/PATO_0000273",
            /* normal fecundity */ "http://purl.obolibrary.org/obo/PATO_0045017",
            /* decreased fecundity */ "http://purl.obolibrary.org/obo/PATO_0001696",
            /* increased fecundity */ "http://purl.obolibrary.org/obo/PATO_0001695",
            /* behavioral quality */ "http://purl.obolibrary.org/obo/PATO_0000186",
            /* handedness */ "http://purl.obolibrary.org/obo/PATO_0002201",
            /* ambidextrous handedness */ "http://purl.obolibrary.org/obo/PATO_0002204",
            /* right handedness */ "http://purl.obolibrary.org/obo/PATO_0002203",
            /* left handedness */ "http://purl.obolibrary.org/obo/PATO_0002202",
            /* receptivity */ "http://purl.obolibrary.org/obo/PATO_0001719",
            /* male receptivity */ "http://purl.obolibrary.org/obo/PATO_0001721",
            /* normal male receptivity */ "http://purl.obolibrary.org/obo/PATO_0045029",
            /* decreased male receptivity */ "http://purl.obolibrary.org/obo/PATO_0001726",
            /* increased male receptivity */ "http://purl.obolibrary.org/obo/PATO_0001725",
            /* female receptivity */ "http://purl.obolibrary.org/obo/PATO_0001720",
            /* normal female receptivity */ "http://purl.obolibrary.org/obo/PATO_0045019",
            /* decreased female receptivity */ "http://purl.obolibrary.org/obo/PATO_0001724",
            /* increased female receptivity */ "http://purl.obolibrary.org/obo/PATO_0001723",
            /* discrimination */ "http://purl.obolibrary.org/obo/PATO_0000189",
            /* discriminate */ "http://purl.obolibrary.org/obo/PATO_0001319",
            /* indiscriminate */ "http://purl.obolibrary.org/obo/PATO_0001318",
            /* preference */ "http://purl.obolibrary.org/obo/PATO_0000773",
            /* indifference */ "http://purl.obolibrary.org/obo/PATO_0000772",
            /* aversion */ "http://purl.obolibrary.org/obo/PATO_0000771",
            /* behavioural activity */ "http://purl.obolibrary.org/obo/PATO_0002265",
            /* behavioural active */ "http://purl.obolibrary.org/obo/PATO_0001707",
            /* normal behavioural activity */ "http://purl.obolibrary.org/obo/PATO_0045007",
            /* decreased behavioural activity */ "http://purl.obolibrary.org/obo/PATO_0000761",
            /* increased behavioural activity */ "http://purl.obolibrary.org/obo/PATO_0000760",
            /* behavioural inactive */ "http://purl.obolibrary.org/obo/PATO_0001706",
            /* lethargic */ "http://purl.obolibrary.org/obo/PATO_0001418",
            /* compatibility */ "http://purl.obolibrary.org/obo/PATO_0000021",
            /* incompatible */ "http://purl.obolibrary.org/obo/PATO_0000345",
            /* compatible */ "http://purl.obolibrary.org/obo/PATO_0000344",
            /* movement behavioral quality */ "http://purl.obolibrary.org/obo/PATO_0002076",
            /* circling */ "http://purl.obolibrary.org/obo/PATO_0001911",
            /* partially paralysed */ "http://purl.obolibrary.org/obo/PATO_0001858",
            /* paralysed */ "http://purl.obolibrary.org/obo/PATO_0000763",
            /* passive */ "http://purl.obolibrary.org/obo/PATO_0000764",
            /* balance */ "http://purl.obolibrary.org/obo/PATO_0000185",
            /* unbalanced */ "http://purl.obolibrary.org/obo/PATO_0000758",
            /* balanced */ "http://purl.obolibrary.org/obo/PATO_0000757",
            /* coordination */ "http://purl.obolibrary.org/obo/PATO_0000188",
            /* uncoordinated */ "http://purl.obolibrary.org/obo/PATO_0000770",
            /* coordinated */ "http://purl.obolibrary.org/obo/PATO_0000769",
            /* normal coordination */ "http://purl.obolibrary.org/obo/PATO_0045012",
            /* decreased coordination */ "http://purl.obolibrary.org/obo/PATO_0001860",
            /* increased coordination */ "http://purl.obolibrary.org/obo/PATO_0001859",
            /* bang sensitive */ "http://purl.obolibrary.org/obo/PATO_0000759",
            /* migratory */ "http://purl.obolibrary.org/obo/PATO_0095008",
            /* mosaicism */ "http://purl.obolibrary.org/obo/PATO_0002048",
            /* virulence */ "http://purl.obolibrary.org/obo/PATO_0002146",
            /* normal virulence */ "http://purl.obolibrary.org/obo/PATO_0045080",
            /* increased virulence */ "http://purl.obolibrary.org/obo/PATO_0002148",
            /* reduced virulence */ "http://purl.obolibrary.org/obo/PATO_0002147",
            /* increased object quality */ "http://purl.obolibrary.org/obo/PATO_0002305",
            /* increased adhesivity */ "http://purl.obolibrary.org/obo/PATO_0002333",
            /* increased radiopacity */ "http://purl.obolibrary.org/obo/PATO_0002144",
            /* increased fragility */ "http://purl.obolibrary.org/obo/PATO_0002055",
            /* increased magnetism */ "http://purl.obolibrary.org/obo/PATO_0001683",
            /* increased waxiness */ "http://purl.obolibrary.org/obo/PATO_0002382",
            /* increased stability */ "http://purl.obolibrary.org/obo/PATO_0015027",
            /* high-arched */ "http://purl.obolibrary.org/obo/PATO_0002162",
            /* increased position */ "http://purl.obolibrary.org/obo/PATO_0001475",
            /* increased angle to */ "http://purl.obolibrary.org/obo/PATO_0002327",
            /* increased elevation */ "http://purl.obolibrary.org/obo/PATO_0001688",
            /* increased distribution */ "http://purl.obolibrary.org/obo/PATO_0001671",
            /* increased distance */ "http://purl.obolibrary.org/obo/PATO_0000374",
            /* increased tonicity */ "http://purl.obolibrary.org/obo/PATO_0001618",
            /* increased mass density */ "http://purl.obolibrary.org/obo/PATO_0001788",
            /* ivory */ "http://purl.obolibrary.org/obo/PATO_0002149",
            /* increased osmolarity */ "http://purl.obolibrary.org/obo/PATO_0001657",
            /* increased affinity */ "http://purl.obolibrary.org/obo/PATO_0002071",
            /* increased elasticity */ "http://purl.obolibrary.org/obo/PATO_0002287",
            /* increased coiling */ "http://purl.obolibrary.org/obo/PATO_0001795",
            /* increased radioactivity */ "http://purl.obolibrary.org/obo/PATO_0001742",
            /* increased strength */ "http://purl.obolibrary.org/obo/PATO_0001778",
            /* decreased fatigability */ "http://purl.obolibrary.org/obo/PATO_0001817",
            /* increased velocity */ "http://purl.obolibrary.org/obo/PATO_0002471",
            /* increased wetness */ "http://purl.obolibrary.org/obo/PATO_0001825",
            /* hard */ "http://purl.obolibrary.org/obo/PATO_0000386",
            /* scirrhous */ "http://purl.obolibrary.org/obo/PATO_0002127",
            /* increased age */ "http://purl.obolibrary.org/obo/PATO_0001764",
            /* increased viscosity */ "http://purl.obolibrary.org/obo/PATO_0001693",
            /* increased pigmentation */ "http://purl.obolibrary.org/obo/PATO_0002250",
            /* increased efficiency */ "http://purl.obolibrary.org/obo/PATO_0001676",
            /* increased fluorescence */ "http://purl.obolibrary.org/obo/PATO_0001926",
            /* increased force */ "http://purl.obolibrary.org/obo/PATO_0002245",
            /* increased weight */ "http://purl.obolibrary.org/obo/PATO_0000582",
            /* increased avidity */ "http://purl.obolibrary.org/obo/PATO_0002074",
            /* increased life span */ "http://purl.obolibrary.org/obo/PATO_0001603",
            /* increased mobility */ "http://purl.obolibrary.org/obo/PATO_0002282",
            /* increased concentration */ "http://purl.obolibrary.org/obo/PATO_0001162",
            /* increased solubility */ "http://purl.obolibrary.org/obo/PATO_0001663",
            /* increased speed */ "http://purl.obolibrary.org/obo/PATO_0000303",
            /* hyperplastic */ "http://purl.obolibrary.org/obo/PATO_0000644",
            /* increased flexibility */ "http://purl.obolibrary.org/obo/PATO_0001776",
            /* increased osmolality */ "http://purl.obolibrary.org/obo/PATO_0002029",
            /* increased mass */ "http://purl.obolibrary.org/obo/PATO_0001563",
            /* increased permeability */ "http://purl.obolibrary.org/obo/PATO_0001577",
            /* increased acidity */ "http://purl.obolibrary.org/obo/PATO_0001844",
            /* increased size */ "http://purl.obolibrary.org/obo/PATO_0000586",
            /* gigantic */ "http://purl.obolibrary.org/obo/PATO_0001940",
            /* increased area */ "http://purl.obolibrary.org/obo/PATO_0002057",
            /* increased height */ "http://purl.obolibrary.org/obo/PATO_0000570",
            /* distended */ "http://purl.obolibrary.org/obo/PATO_0001602",
            /* increased width */ "http://purl.obolibrary.org/obo/PATO_0000600",
            /* increased width and length */ "http://purl.obolibrary.org/obo/PATO_0040031",
            /* increased thickness */ "http://purl.obolibrary.org/obo/PATO_0000591",
            /* increased volume */ "http://purl.obolibrary.org/obo/PATO_0000595",
            /* ballooning */ "http://purl.obolibrary.org/obo/PATO_0002093",
            /* hypertrophic */ "http://purl.obolibrary.org/obo/PATO_0000584",
            /* increased diameter */ "http://purl.obolibrary.org/obo/PATO_0001714",
            /* increased anterior-posterior diameter */ "http://purl.obolibrary.org/obo/PATO_0002043",
            /* dilated */ "http://purl.obolibrary.org/obo/PATO_0001571",
            /* increased depth */ "http://purl.obolibrary.org/obo/PATO_0001596",
            /* increased length */ "http://purl.obolibrary.org/obo/PATO_0000573",
            /* increased perimeter */ "http://purl.obolibrary.org/obo/PATO_0001712",
            /* increased circumference */ "http://purl.obolibrary.org/obo/PATO_0001898",
            /* increased contractility */ "http://purl.obolibrary.org/obo/PATO_0001580",
            /* increased humidity */ "http://purl.obolibrary.org/obo/PATO_0015010",
            /* increased odor */ "http://purl.obolibrary.org/obo/PATO_0001893",
            /* having extra function */ "http://purl.obolibrary.org/obo/PATO_0001559",
            /* increased functionality */ "http://purl.obolibrary.org/obo/PATO_0001625",
            /* increased fluid flow */ "http://purl.obolibrary.org/obo/PATO_0001839",
            /* increased cellular motility */ "http://purl.obolibrary.org/obo/PATO_0002298",
            /* increased porosity */ "http://purl.obolibrary.org/obo/PATO_0015024",
            /* increased tendency */ "http://purl.obolibrary.org/obo/PATO_0002361",
            /* increased phosphorylation */ "http://purl.obolibrary.org/obo/PATO_0002221",
            /* increased degree of illumination */ "http://purl.obolibrary.org/obo/PATO_0015014",
            /* increased pressure */ "http://purl.obolibrary.org/obo/PATO_0001576",
            /* increased curvature */ "http://purl.obolibrary.org/obo/PATO_0001592",
            /* increased contamination */ "http://purl.obolibrary.org/obo/PATO_0015032",
            /* increased turgor */ "http://purl.obolibrary.org/obo/PATO_0001622",
            /* increased combustibility */ "http://purl.obolibrary.org/obo/PATO_0015022",
            /* increased temperature */ "http://purl.obolibrary.org/obo/PATO_0001305",
            /* severe increased temperature */ "http://purl.obolibrary.org/obo/PATO_0001317",
            /* moderate increased temperature */ "http://purl.obolibrary.org/obo/PATO_0001316",
            /* mild increased temperature */ "http://purl.obolibrary.org/obo/PATO_0001315",
            /* disposition */ "http://purl.obolibrary.org/obo/PATO_0001727",
            /* tendency */ "http://purl.obolibrary.org/obo/PATO_0002360",
            /* normal tendency */ "http://purl.obolibrary.org/obo/PATO_0045068",
            /* decreased tendency */ "http://purl.obolibrary.org/obo/PATO_0002362",
            /* multi-cellular organismal disposition */ "http://purl.obolibrary.org/obo/PATO_0001728",
            /* normal object quality */ "http://purl.obolibrary.org/obo/PATO_0045001",
            /* normal cellular motility */ "http://purl.obolibrary.org/obo/PATO_0045008",
            /* normal velocity */ "http://purl.obolibrary.org/obo/PATO_0045079",
            /* normal osmolarity */ "http://purl.obolibrary.org/obo/PATO_0045035",
            /* normal magnetism */ "http://purl.obolibrary.org/obo/PATO_0045027",
            /* normal fluid flow */ "http://purl.obolibrary.org/obo/PATO_0045022",
            /* normal distance */ "http://purl.obolibrary.org/obo/PATO_0045014",
            /* normal position */ "http://purl.obolibrary.org/obo/PATO_0045040",
            /* normal elevation */ "http://purl.obolibrary.org/obo/PATO_0045043",
            /* normal distribution */ "http://purl.obolibrary.org/obo/PATO_0045042",
            /* normal angle to */ "http://purl.obolibrary.org/obo/PATO_0045041",
            /* normal phosphorylation */ "http://purl.obolibrary.org/obo/PATO_0045037",
            /* normal turgor */ "http://purl.obolibrary.org/obo/PATO_0045072",
            /* normal force */ "http://purl.obolibrary.org/obo/PATO_0045024",
            /* normal radioactivity */ "http://purl.obolibrary.org/obo/PATO_0045045",
            /* normal elasticity */ "http://purl.obolibrary.org/obo/PATO_0045016",
            /* normal viscosity */ "http://purl.obolibrary.org/obo/PATO_0045081",
            /* normal contractility */ "http://purl.obolibrary.org/obo/PATO_0045011",
            /* normal pigmentation */ "http://purl.obolibrary.org/obo/PATO_0045039",
            /* normal adhesivity */ "http://purl.obolibrary.org/obo/PATO_0045003",
            /* normal life span */ "http://purl.obolibrary.org/obo/PATO_0045026",
            /* normal curvature */ "http://purl.obolibrary.org/obo/PATO_0045013",
            /* normal tonicity */ "http://purl.obolibrary.org/obo/PATO_0045071",
            /* normal mass */ "http://purl.obolibrary.org/obo/PATO_0045030",
            /* normal pressure */ "http://purl.obolibrary.org/obo/PATO_0045044",
            /* normal age */ "http://purl.obolibrary.org/obo/PATO_0045005",
            /* normal wetness */ "http://purl.obolibrary.org/obo/PATO_0045083",
            /* normal speed */ "http://purl.obolibrary.org/obo/PATO_0045063",
            /* normal efficiency */ "http://purl.obolibrary.org/obo/PATO_0045015",
            /* normal size */ "http://purl.obolibrary.org/obo/PATO_0045050",
            /* normal width */ "http://purl.obolibrary.org/obo/PATO_0045061",
            /* normal volume */ "http://purl.obolibrary.org/obo/PATO_0045060",
            /* normal thickness */ "http://purl.obolibrary.org/obo/PATO_0045059",
            /* normal diameter */ "http://purl.obolibrary.org/obo/PATO_0045055",
            /* normal anterior-posterior diameter */ "http://purl.obolibrary.org/obo/PATO_0045056",
            /* normal length */ "http://purl.obolibrary.org/obo/PATO_0045054",
            /* normal perimeter */ "http://purl.obolibrary.org/obo/PATO_0045057",
            /* normal circumference */ "http://purl.obolibrary.org/obo/PATO_0045058",
            /* normal height */ "http://purl.obolibrary.org/obo/PATO_0045053",
            /* normal depth */ "http://purl.obolibrary.org/obo/PATO_0045052",
            /* normal area */ "http://purl.obolibrary.org/obo/PATO_0045051",
            /* normal acidity */ "http://purl.obolibrary.org/obo/PATO_0045002",
            /* normal mobility */ "http://purl.obolibrary.org/obo/PATO_0045032",
            /* normal radiopacity */ "http://purl.obolibrary.org/obo/PATO_0045046",
            /* normal affinity */ "http://purl.obolibrary.org/obo/PATO_0045004",
            /* normal waxiness */ "http://purl.obolibrary.org/obo/PATO_0045082",
            /* normal solubility */ "http://purl.obolibrary.org/obo/PATO_0045062",
            /* normal osmolality */ "http://purl.obolibrary.org/obo/PATO_0045034",
            /* normal flexibility */ "http://purl.obolibrary.org/obo/PATO_0045021",
            /* normal permeability */ "http://purl.obolibrary.org/obo/PATO_0045036",
            /* normal mass density */ "http://purl.obolibrary.org/obo/PATO_0045031",
            /* normal avidity */ "http://purl.obolibrary.org/obo/PATO_0045006",
            /* normal coiling */ "http://purl.obolibrary.org/obo/PATO_0045009",
            /* normal temperature */ "http://purl.obolibrary.org/obo/PATO_0045067",
            /* normal fluorescence */ "http://purl.obolibrary.org/obo/PATO_0045023",
            /* normal concentration */ "http://purl.obolibrary.org/obo/PATO_0045010",
            /* normal odor */ "http://purl.obolibrary.org/obo/PATO_0045033",
            /* normal strength */ "http://purl.obolibrary.org/obo/PATO_0045064",
            /* normal fatigability */ "http://purl.obolibrary.org/obo/PATO_0045065",
            /* normal fragility */ "http://purl.obolibrary.org/obo/PATO_0045025",
            /* physical quality */ "http://purl.obolibrary.org/obo/PATO_0001018",
            /* conductivity */ "http://purl.obolibrary.org/obo/PATO_0001585",
            /* electrical conductivity */ "http://purl.obolibrary.org/obo/PATO_0001757",
            /* nerve conductivity */ "http://purl.obolibrary.org/obo/PATO_0001758",
            /* heat conductivity */ "http://purl.obolibrary.org/obo/PATO_0001756",
            /* neuron projection quality */ "http://purl.obolibrary.org/obo/PATO_0070033",
            /* intratelencephalic projecting */ "http://purl.obolibrary.org/obo/PATO_0070034",
            /* corticomedulla projecting */ "http://purl.obolibrary.org/obo/PATO_0070031",
            /* near projecting */ "http://purl.obolibrary.org/obo/PATO_0070030",
            /* corticothalamic projecting */ "http://purl.obolibrary.org/obo/PATO_0070029",
            /* extratelencephalic projecting */ "http://purl.obolibrary.org/obo/PATO_0070028",
            /* potability */ "http://purl.obolibrary.org/obo/PATO_0025000",
            /* potable */ "http://purl.obolibrary.org/obo/PATO_0025002",
            /* non-potable */ "http://purl.obolibrary.org/obo/PATO_0025001",
            /* direction */ "http://purl.obolibrary.org/obo/PATO_0000039",
            /* back */ "http://purl.obolibrary.org/obo/PATO_0001901",
            /* bi-directional */ "http://purl.obolibrary.org/obo/PATO_0001903",
            /* left */ "http://purl.obolibrary.org/obo/PATO_0000366",
            /* front */ "http://purl.obolibrary.org/obo/PATO_0001900",
            /* unidirectional */ "http://purl.obolibrary.org/obo/PATO_0001902",
            /* cardinal direction */ "http://purl.obolibrary.org/obo/PATO_0045090",
            /* west */ "http://purl.obolibrary.org/obo/PATO_0045094",
            /* south */ "http://purl.obolibrary.org/obo/PATO_0045093",
            /* east */ "http://purl.obolibrary.org/obo/PATO_0045092",
            /* north */ "http://purl.obolibrary.org/obo/PATO_0045091",
            /* down */ "http://purl.obolibrary.org/obo/PATO_0000365",
            /* circling direction */ "http://purl.obolibrary.org/obo/PATO_0001904",
            /* up */ "http://purl.obolibrary.org/obo/PATO_0000370",
            /* right */ "http://purl.obolibrary.org/obo/PATO_0000367",
            /* strength */ "http://purl.obolibrary.org/obo/PATO_0001230",
            /* fatigability */ "http://purl.obolibrary.org/obo/PATO_0001815",
            /* increased fatigability */ "http://purl.obolibrary.org/obo/PATO_0001816",
            /* decreased strength */ "http://purl.obolibrary.org/obo/PATO_0001779",
            /* elasticity */ "http://purl.obolibrary.org/obo/PATO_0001031",
            /* inelastic */ "http://purl.obolibrary.org/obo/PATO_0001172",
            /* elastic */ "http://purl.obolibrary.org/obo/PATO_0001171",
            /* decreased elasticity */ "http://purl.obolibrary.org/obo/PATO_0002288",
            /* catalytic activity */ "http://purl.obolibrary.org/obo/PATO_0001414",
            /* time */ "http://purl.obolibrary.org/obo/PATO_0000165",
            /* life span */ "http://purl.obolibrary.org/obo/PATO_0000050",
            /* decreased life span */ "http://purl.obolibrary.org/obo/PATO_0001604",
            /* age */ "http://purl.obolibrary.org/obo/PATO_0000011",
            /* young */ "http://purl.obolibrary.org/obo/PATO_0000309",
            /* old */ "http://purl.obolibrary.org/obo/PATO_0000308",
            /* senescent */ "http://purl.obolibrary.org/obo/PATO_0001487",
            /* chronological age */ "http://purl.obolibrary.org/obo/PATO_0001486",
            /* decreased age */ "http://purl.obolibrary.org/obo/PATO_0001765",
            /* impulse */ "http://purl.obolibrary.org/obo/PATO_0001022",
            /* thermoplastic */ "http://purl.obolibrary.org/obo/PATO_0040070",
            /* adhesivity */ "http://purl.obolibrary.org/obo/PATO_0001431",
            /* non-adhesive */ "http://purl.obolibrary.org/obo/PATO_0002332",
            /* adhesive */ "http://purl.obolibrary.org/obo/PATO_0002331",
            /* decreased adhesivity */ "http://purl.obolibrary.org/obo/PATO_0002334",
            /* cellular adhesivity */ "http://purl.obolibrary.org/obo/PATO_0001531",
            /* heterophilic */ "http://purl.obolibrary.org/obo/PATO_0001441",
            /* homophilic */ "http://purl.obolibrary.org/obo/PATO_0001440",
            /* edibility */ "http://purl.obolibrary.org/obo/PATO_0002138",
            /* inedible */ "http://purl.obolibrary.org/obo/PATO_0002140",
            /* edible */ "http://purl.obolibrary.org/obo/PATO_0002139",
            /* protruding */ "http://purl.obolibrary.org/obo/PATO_0001598",
            /* prominent */ "http://purl.obolibrary.org/obo/PATO_0001482",
            /* umbonate */ "http://purl.obolibrary.org/obo/PATO_0001358",
            /* exserted */ "http://purl.obolibrary.org/obo/PATO_0000623",
            /* papillary */ "http://purl.obolibrary.org/obo/PATO_0002108",
            /* caudate */ "http://purl.obolibrary.org/obo/PATO_0001880",
            /* protruding into */ "http://purl.obolibrary.org/obo/PATO_0001645",
            /* herniated into */ "http://purl.obolibrary.org/obo/PATO_0002163",
            /* protruding out of */ "http://purl.obolibrary.org/obo/PATO_0001646",
            /* herniated out of */ "http://purl.obolibrary.org/obo/PATO_0002167",
            /* rhizoidal */ "http://purl.obolibrary.org/obo/PATO_0001364",
            /* herniated */ "http://purl.obolibrary.org/obo/PATO_0000643",
            /* hydrophobicity */ "http://purl.obolibrary.org/obo/PATO_0001884",
            /* hydrophobic */ "http://purl.obolibrary.org/obo/PATO_0001885",
            /* wavelength */ "http://purl.obolibrary.org/obo/PATO_0001242",
            /* sound wavelength */ "http://purl.obolibrary.org/obo/PATO_0001523",
            /* magnetism */ "http://purl.obolibrary.org/obo/PATO_0001682",
            /* non-magnetic */ "http://purl.obolibrary.org/obo/PATO_0001686",
            /* magnetic */ "http://purl.obolibrary.org/obo/PATO_0001685",
            /* decreased magnetism */ "http://purl.obolibrary.org/obo/PATO_0001684",
            /* temperature */ "http://purl.obolibrary.org/obo/PATO_0000146",
            /* decreased temperature */ "http://purl.obolibrary.org/obo/PATO_0001306",
            /* frozen */ "http://purl.obolibrary.org/obo/PATO_0001985",
            /* wetness */ "http://purl.obolibrary.org/obo/PATO_0001822",
            /* dry */ "http://purl.obolibrary.org/obo/PATO_0001824",
            /* wet */ "http://purl.obolibrary.org/obo/PATO_0001823",
            /* decreased wetness */ "http://purl.obolibrary.org/obo/PATO_0001826",
            /* prominence */ "http://purl.obolibrary.org/obo/PATO_0015007",
            /* force */ "http://purl.obolibrary.org/obo/PATO_0001035",
            /* decreased force */ "http://purl.obolibrary.org/obo/PATO_0002246",
            /* decreased weight */ "http://purl.obolibrary.org/obo/PATO_0000583",
            /* weight */ "http://purl.obolibrary.org/obo/PATO_0000128",
            /* position */ "http://purl.obolibrary.org/obo/PATO_0000140",
            /* lateralized */ "http://purl.obolibrary.org/obo/PATO_0000626",
            /* concealed */ "http://purl.obolibrary.org/obo/PATO_0002508",
            /* orientation */ "http://purl.obolibrary.org/obo/PATO_0000133",
            /* posterolateral orientation */ "http://purl.obolibrary.org/obo/PATO_0002503",
            /* dorsal orientation */ "http://purl.obolibrary.org/obo/PATO_0002495",
            /* medial orientation */ "http://purl.obolibrary.org/obo/PATO_0002498",
            /* dorsomedial orientation */ "http://purl.obolibrary.org/obo/PATO_0040013",
            /* angle */ "http://purl.obolibrary.org/obo/PATO_0002326",
            /* right angle to */ "http://purl.obolibrary.org/obo/PATO_0001321",
            /* convex angle to */ "http://purl.obolibrary.org/obo/PATO_0001053",
            /* external angle */ "http://purl.obolibrary.org/obo/PATO_0002025",
            /* internal angle */ "http://purl.obolibrary.org/obo/PATO_0001054",
            /* reflex angle to */ "http://purl.obolibrary.org/obo/PATO_0001055",
            /* straight angle to */ "http://purl.obolibrary.org/obo/PATO_0001322",
            /* obtuse angle to */ "http://purl.obolibrary.org/obo/PATO_0001052",
            /* decreased angle to */ "http://purl.obolibrary.org/obo/PATO_0002328",
            /* acute angle to */ "http://purl.obolibrary.org/obo/PATO_0001051",
            /* mesiolateral orientation */ "http://purl.obolibrary.org/obo/PATO_0002524",
            /* transverse orientation */ "http://purl.obolibrary.org/obo/PATO_0002502",
            /* distal orientation */ "http://purl.obolibrary.org/obo/PATO_0002494",
            /* displaced to */ "http://purl.obolibrary.org/obo/PATO_0002168",
            /* adduction */ "http://purl.obolibrary.org/obo/PATO_0002133",
            /* abduction */ "http://purl.obolibrary.org/obo/PATO_0002131",
            /* anterodorsal orientation */ "http://purl.obolibrary.org/obo/PATO_0002491",
            /* anterodistal orientation */ "http://purl.obolibrary.org/obo/PATO_0002526",
            /* posterodorsal orientation */ "http://purl.obolibrary.org/obo/PATO_0002504",
            /* dorsolateral orientation */ "http://purl.obolibrary.org/obo/PATO_0002496",
            /* parasagittal orientation */ "http://purl.obolibrary.org/obo/PATO_0005023",
            /* anteromedial orientation */ "http://purl.obolibrary.org/obo/PATO_0002514",
            /* anterior orientation */ "http://purl.obolibrary.org/obo/PATO_0002490",
            /* posteroventral orientation */ "http://purl.obolibrary.org/obo/PATO_0040005",
            /* ventrolaterally orientation */ "http://purl.obolibrary.org/obo/PATO_0002500",
            /* anterolateral orientation */ "http://purl.obolibrary.org/obo/PATO_0002492",
            /* longitudinal orientation */ "http://purl.obolibrary.org/obo/PATO_0005024",
            /* posterodistal orientation */ "http://purl.obolibrary.org/obo/PATO_0002527",
            /* lateral orientation */ "http://purl.obolibrary.org/obo/PATO_0002497",
            /* distomedial orientation */ "http://purl.obolibrary.org/obo/PATO_0040023",
            /* ventral orientation */ "http://purl.obolibrary.org/obo/PATO_0002501",
            /* anteroventral orientation */ "http://purl.obolibrary.org/obo/PATO_0002493",
            /* oriented towards */ "http://purl.obolibrary.org/obo/PATO_0002448",
            /* oblique orientation */ "http://purl.obolibrary.org/obo/PATO_0002481",
            /* posterior orientation */ "http://purl.obolibrary.org/obo/PATO_0002499",
            /* alignment */ "http://purl.obolibrary.org/obo/PATO_0001652",
            /* misaligned with */ "http://purl.obolibrary.org/obo/PATO_0001654",
            /* misaligned away from */ "http://purl.obolibrary.org/obo/PATO_0002174",
            /* misaligned towards */ "http://purl.obolibrary.org/obo/PATO_0002173",
            /* aligned with */ "http://purl.obolibrary.org/obo/PATO_0001653",
            /* adjacent to */ "http://purl.obolibrary.org/obo/PATO_0002259",
            /* distal to */ "http://purl.obolibrary.org/obo/PATO_0001234",
            /* diagonal to */ "http://purl.obolibrary.org/obo/PATO_0002392",
            /* opposite */ "http://purl.obolibrary.org/obo/PATO_0001933",
            /* stacked */ "http://purl.obolibrary.org/obo/PATO_0002077",
            /* located in */ "http://purl.obolibrary.org/obo/PATO_0002261",
            /* adaxial to */ "http://purl.obolibrary.org/obo/PATO_0002047",
            /* exposed */ "http://purl.obolibrary.org/obo/PATO_0002425",
            /* proximal to */ "http://purl.obolibrary.org/obo/PATO_0001195",
            /* dispersed */ "http://purl.obolibrary.org/obo/PATO_0001630",
            /* oblique to */ "http://purl.obolibrary.org/obo/PATO_0002513",
            /* misrouted */ "http://purl.obolibrary.org/obo/PATO_0000629",
            /* subterminal */ "http://purl.obolibrary.org/obo/PATO_0002475",
            /* medial to */ "http://purl.obolibrary.org/obo/PATO_0001191",
            /* surrounding */ "http://purl.obolibrary.org/obo/PATO_0001772",
            /* abutting */ "http://purl.obolibrary.org/obo/PATO_0002435",
            /* vertical to */ "http://purl.obolibrary.org/obo/PATO_0002523",
            /* basal to */ "http://purl.obolibrary.org/obo/PATO_0002349",
            /* ipsilateral to */ "http://purl.obolibrary.org/obo/PATO_0002035",
            /* far from */ "http://purl.obolibrary.org/obo/PATO_0002233",
            /* dorsal to */ "http://purl.obolibrary.org/obo/PATO_0001233",
            /* centered */ "http://purl.obolibrary.org/obo/PATO_0002268",
            /* lateral to */ "http://purl.obolibrary.org/obo/PATO_0001193",
            /* erect */ "http://purl.obolibrary.org/obo/PATO_0000622",
            /* semi erect */ "http://purl.obolibrary.org/obo/PATO_0002260",
            /* external to */ "http://purl.obolibrary.org/obo/PATO_0002483",
            /* horizontal */ "http://purl.obolibrary.org/obo/PATO_0001855",
            /* anterior to */ "http://purl.obolibrary.org/obo/PATO_0001632",
            /* anteroventral to */ "http://purl.obolibrary.org/obo/PATO_0001917",
            /* anterodorsal to */ "http://purl.obolibrary.org/obo/PATO_0001915",
            /* decreased position */ "http://purl.obolibrary.org/obo/PATO_0001476",
            /* decreased elevation */ "http://purl.obolibrary.org/obo/PATO_0001689",
            /* decreased distribution */ "http://purl.obolibrary.org/obo/PATO_0001672",
            /* perpendicular to */ "http://purl.obolibrary.org/obo/PATO_0002434",
            /* peripheral */ "http://purl.obolibrary.org/obo/PATO_0002107",
            /* transmural */ "http://purl.obolibrary.org/obo/PATO_0002417",
            /* position originates from */ "http://purl.obolibrary.org/obo/PATO_0040002",
            /* concealed by */ "http://purl.obolibrary.org/obo/PATO_0002516",
            /* crowded */ "http://purl.obolibrary.org/obo/PATO_0000619",
            /* parallel to */ "http://purl.obolibrary.org/obo/PATO_0002272",
            /* terminal */ "http://purl.obolibrary.org/obo/PATO_0002476",
            /* ventral to */ "http://purl.obolibrary.org/obo/PATO_0001196",
            /* offset */ "http://purl.obolibrary.org/obo/PATO_0002436",
            /* elevation */ "http://purl.obolibrary.org/obo/PATO_0001687",
            /* increased elevation relative to */ "http://purl.obolibrary.org/obo/PATO_0002515",
            /* sunken */ "http://purl.obolibrary.org/obo/PATO_0002416",
            /* cauline to */ "http://purl.obolibrary.org/obo/PATO_0002350",
            /* extends to */ "http://purl.obolibrary.org/obo/PATO_0002463",
            /* prostrate */ "http://purl.obolibrary.org/obo/PATO_0000631",
            /* procumbent */ "http://purl.obolibrary.org/obo/PATO_0002389",
            /* decumbent */ "http://purl.obolibrary.org/obo/PATO_0002343",
            /* axial to */ "http://purl.obolibrary.org/obo/PATO_0002036",
            /* left side of */ "http://purl.obolibrary.org/obo/PATO_0001792",
            /* flush */ "http://purl.obolibrary.org/obo/PATO_0002518",
            /* level with */ "http://purl.obolibrary.org/obo/PATO_0002443",
            /* distichous */ "http://purl.obolibrary.org/obo/PATO_0001952",
            /* bilateral */ "http://purl.obolibrary.org/obo/PATO_0000618",
            /* extends beyond */ "http://purl.obolibrary.org/obo/PATO_0002464",
            /* abaxial to */ "http://purl.obolibrary.org/obo/PATO_0002046",
            /* superficial */ "http://purl.obolibrary.org/obo/PATO_0001665",
            /* posterior to */ "http://purl.obolibrary.org/obo/PATO_0001633",
            /* posteromedial to */ "http://purl.obolibrary.org/obo/PATO_0002449",
            /* posteroventral to */ "http://purl.obolibrary.org/obo/PATO_0001918",
            /* posterodorsal to */ "http://purl.obolibrary.org/obo/PATO_0001916",
            /* divergent from */ "http://purl.obolibrary.org/obo/PATO_0002424",
            /* retracted */ "http://purl.obolibrary.org/obo/PATO_0001477",
            /* positional polarity */ "http://purl.obolibrary.org/obo/PATO_0001769",
            /* mediolateral polarity */ "http://purl.obolibrary.org/obo/PATO_0002373",
            /* anterior-posterior polarity */ "http://purl.obolibrary.org/obo/PATO_0002024",
            /* apical-basal polarity */ "http://purl.obolibrary.org/obo/PATO_0002023",
            /* dorsal-ventral polarity */ "http://purl.obolibrary.org/obo/PATO_0001775",
            /* right side of */ "http://purl.obolibrary.org/obo/PATO_0001793",
            /* continuous with */ "http://purl.obolibrary.org/obo/PATO_0005011",
            /* surrounded by */ "http://purl.obolibrary.org/obo/PATO_0005016",
            /* inserted into */ "http://purl.obolibrary.org/obo/PATO_0000624",
            /* confluent with */ "http://purl.obolibrary.org/obo/PATO_0002512",
            /* displaced */ "http://purl.obolibrary.org/obo/PATO_0002181",
            /* dislocated */ "http://purl.obolibrary.org/obo/PATO_0001852",
            /* partially dislocated */ "http://purl.obolibrary.org/obo/PATO_0002157",
            /* mislocalised */ "http://purl.obolibrary.org/obo/PATO_0000628",
            /* mislocalized adaxially */ "http://purl.obolibrary.org/obo/PATO_0002396",
            /* mislocalized abaxially */ "http://purl.obolibrary.org/obo/PATO_0002395",
            /* mislocalised radially */ "http://purl.obolibrary.org/obo/PATO_0002178",
            /* mislocalised ventrally */ "http://purl.obolibrary.org/obo/PATO_0001920",
            /* mislocalised dorsally */ "http://purl.obolibrary.org/obo/PATO_0001919",
            /* mislocalised proximally */ "http://purl.obolibrary.org/obo/PATO_0002179",
            /* mislocalised posteriorly */ "http://purl.obolibrary.org/obo/PATO_0001922",
            /* mislocalised medially */ "http://purl.obolibrary.org/obo/PATO_0001924",
            /* mislocalised anteriorly */ "http://purl.obolibrary.org/obo/PATO_0001921",
            /* mislocalised laterally */ "http://purl.obolibrary.org/obo/PATO_0001923",
            /* consistency */ "http://purl.obolibrary.org/obo/PATO_0000037",
            /* work */ "http://purl.obolibrary.org/obo/PATO_0001026",
            /* energy */ "http://purl.obolibrary.org/obo/PATO_0001021",
            /* contractility */ "http://purl.obolibrary.org/obo/PATO_0001579",
            /* non-contractile */ "http://purl.obolibrary.org/obo/PATO_0001691",
            /* contractile */ "http://purl.obolibrary.org/obo/PATO_0001690",
            /* decreased contractility */ "http://purl.obolibrary.org/obo/PATO_0001581",
            /* tonicity */ "http://purl.obolibrary.org/obo/PATO_0001439",
            /* atonicity */ "http://purl.obolibrary.org/obo/PATO_0001813",
            /* decreased tonicity */ "http://purl.obolibrary.org/obo/PATO_0001619",
            /* dystonicity */ "http://purl.obolibrary.org/obo/PATO_0001814",
            /* proportionality to */ "http://purl.obolibrary.org/obo/PATO_0001470",
            /* increased proportionality to */ "http://purl.obolibrary.org/obo/PATO_0040043",
            /* decreased proportionality to */ "http://purl.obolibrary.org/obo/PATO_0040042",
            /* hydrophilicity */ "http://purl.obolibrary.org/obo/PATO_0001886",
            /* hydrophilic */ "http://purl.obolibrary.org/obo/PATO_0001887",
            /* tension */ "http://purl.obolibrary.org/obo/PATO_0002284",
            /* buoyancy */ "http://purl.obolibrary.org/obo/PATO_0001420",
            /* distance */ "http://purl.obolibrary.org/obo/PATO_0000040",
            /* insertion depth */ "http://purl.obolibrary.org/obo/PATO_0002207",
            /* decreased distance */ "http://purl.obolibrary.org/obo/PATO_0000375",
            /* strain */ "http://purl.obolibrary.org/obo/PATO_0001034",
            /* surface tension */ "http://purl.obolibrary.org/obo/PATO_0001461",
            /* movement quality */ "http://purl.obolibrary.org/obo/PATO_0001906",
            /* velocity */ "http://purl.obolibrary.org/obo/PATO_0002242",
            /* decreased velocity */ "http://purl.obolibrary.org/obo/PATO_0002472",
            /* angular velocity */ "http://purl.obolibrary.org/obo/PATO_0001413",
            /* flow rate */ "http://purl.obolibrary.org/obo/PATO_0001574",
            /* mass flow rate */ "http://purl.obolibrary.org/obo/PATO_0002244",
            /* fluid flow rate */ "http://purl.obolibrary.org/obo/PATO_0002243",
            /* decreased fluid flow */ "http://purl.obolibrary.org/obo/PATO_0001838",
            /* flux */ "http://purl.obolibrary.org/obo/PATO_0001030",
            /* acceleration */ "http://purl.obolibrary.org/obo/PATO_0001028",
            /* angular acceleration */ "http://purl.obolibrary.org/obo/PATO_0001350",
            /* momentum */ "http://purl.obolibrary.org/obo/PATO_0001023",
            /* speed */ "http://purl.obolibrary.org/obo/PATO_0000008",
            /* sound speed */ "http://purl.obolibrary.org/obo/PATO_0001522",
            /* decreased speed */ "http://purl.obolibrary.org/obo/PATO_0000304",
            /* efficiency */ "http://purl.obolibrary.org/obo/PATO_0001029",
            /* efficient */ "http://purl.obolibrary.org/obo/PATO_0001678",
            /* decreased efficiency */ "http://purl.obolibrary.org/obo/PATO_0001675",
            /* inefficient */ "http://purl.obolibrary.org/obo/PATO_0001677",
            /* power */ "http://purl.obolibrary.org/obo/PATO_0001024",
            /* shaded */ "http://purl.obolibrary.org/obo/PATO_0040029",
            /* mass density */ "http://purl.obolibrary.org/obo/PATO_0001019",
            /* irregular density */ "http://purl.obolibrary.org/obo/PATO_0002141",
            /* volumetric density */ "http://purl.obolibrary.org/obo/PATO_0001353",
            /* linear density */ "http://purl.obolibrary.org/obo/PATO_0001352",
            /* area density */ "http://purl.obolibrary.org/obo/PATO_0001351",
            /* dense */ "http://purl.obolibrary.org/obo/PATO_0001164",
            /* decreased mass density */ "http://purl.obolibrary.org/obo/PATO_0001790",
            /* combustibility */ "http://purl.obolibrary.org/obo/PATO_0015021",
            /* decreased combustibility */ "http://purl.obolibrary.org/obo/PATO_0015023",
            /* pressure */ "http://purl.obolibrary.org/obo/PATO_0001025",
            /* decreased pressure */ "http://purl.obolibrary.org/obo/PATO_0001575",
            /* mobility */ "http://purl.obolibrary.org/obo/PATO_0000004",
            /* immobile relative to */ "http://purl.obolibrary.org/obo/PATO_0040011",
            /* mobile relative to */ "http://purl.obolibrary.org/obo/PATO_0040010",
            /* immobile */ "http://purl.obolibrary.org/obo/PATO_0000300",
            /* mobile */ "http://purl.obolibrary.org/obo/PATO_0000299",
            /* decreased mobility */ "http://purl.obolibrary.org/obo/PATO_0002283",
            /* quality of a substance */ "http://purl.obolibrary.org/obo/PATO_0002198",
            /* quality of a solid */ "http://purl.obolibrary.org/obo/PATO_0001546",
            /* solid configuration */ "http://purl.obolibrary.org/obo/PATO_0001736",
            /* crystal configuration */ "http://purl.obolibrary.org/obo/PATO_0002066",
            /* flexibility */ "http://purl.obolibrary.org/obo/PATO_0001543",
            /* inflexible */ "http://purl.obolibrary.org/obo/PATO_0001545",
            /* flexible */ "http://purl.obolibrary.org/obo/PATO_0001544",
            /* decreased flexibility */ "http://purl.obolibrary.org/obo/PATO_0001777",
            /* hardness */ "http://purl.obolibrary.org/obo/PATO_0000048",
            /* firm */ "http://purl.obolibrary.org/obo/PATO_0002450",
            /* soft */ "http://purl.obolibrary.org/obo/PATO_0000387",
            /* irradiated */ "http://purl.obolibrary.org/obo/PATO_0001744",
            /* quality of a colloid */ "http://purl.obolibrary.org/obo/PATO_0015017",
            /* quality of an aerosol */ "http://purl.obolibrary.org/obo/PATO_0015018",
            /* meltability */ "http://purl.obolibrary.org/obo/PATO_0002199",
            /* quality of a liquid */ "http://purl.obolibrary.org/obo/PATO_0001548",
            /* miscibility */ "http://purl.obolibrary.org/obo/PATO_0001888",
            /* liquid configuration */ "http://purl.obolibrary.org/obo/PATO_0001735",
            /* viscosity */ "http://purl.obolibrary.org/obo/PATO_0000992",
            /* viscous */ "http://purl.obolibrary.org/obo/PATO_0000998",
            /* decreased viscosity */ "http://purl.obolibrary.org/obo/PATO_0001694",
            /* quality of a gas */ "http://purl.obolibrary.org/obo/PATO_0001547",
            /* quality of a plasma */ "http://purl.obolibrary.org/obo/PATO_0015012",
            /* humidity */ "http://purl.obolibrary.org/obo/PATO_0015009",
            /* decreased humidity */ "http://purl.obolibrary.org/obo/PATO_0015011",
            /* gaseus configuration */ "http://purl.obolibrary.org/obo/PATO_0001737",
            /* quality of a suspension */ "http://purl.obolibrary.org/obo/PATO_0015029",
            /* quality of a fluid */ "http://purl.obolibrary.org/obo/PATO_0080001",
            /* quality of interaction of a substance with electromagnetic radiation */ "http://purl.obolibrary.org/obo/PATO_0070060",
            /* degree of illumination */ "http://purl.obolibrary.org/obo/PATO_0015013",
            /* decreased degree of illumination */ "http://purl.obolibrary.org/obo/PATO_0015015",
            /* radiopacity */ "http://purl.obolibrary.org/obo/PATO_0002136",
            /* decreased radiopacity */ "http://purl.obolibrary.org/obo/PATO_0002145",
            /* radiopaque */ "http://purl.obolibrary.org/obo/PATO_0002137",
            /* optical quality */ "http://purl.obolibrary.org/obo/PATO_0001300",
            /* focus */ "http://purl.obolibrary.org/obo/PATO_0001516",
            /* blurry */ "http://purl.obolibrary.org/obo/PATO_0001518",
            /* focused */ "http://purl.obolibrary.org/obo/PATO_0001517",
            /* chromatic property */ "http://purl.obolibrary.org/obo/PATO_0001301",
            /* color saturation */ "http://purl.obolibrary.org/obo/PATO_0000017",
            /* high saturation */ "http://purl.obolibrary.org/obo/PATO_0001229",
            /* low saturation */ "http://purl.obolibrary.org/obo/PATO_0000328",
            /* color hue */ "http://purl.obolibrary.org/obo/PATO_0000015",
            /* luminous flux */ "http://purl.obolibrary.org/obo/PATO_0001296",
            /* luminance */ "http://purl.obolibrary.org/obo/PATO_0001718",
            /* fluorescence */ "http://purl.obolibrary.org/obo/PATO_0000018",
            /* autofluorescence */ "http://purl.obolibrary.org/obo/PATO_0001868",
            /* phosphorescence */ "http://purl.obolibrary.org/obo/PATO_0001298",
            /* fluorescent */ "http://purl.obolibrary.org/obo/PATO_0001290",
            /* decreased fluorescence */ "http://purl.obolibrary.org/obo/PATO_0001927",
            /* opacity */ "http://purl.obolibrary.org/obo/PATO_0000957",
            /* translucent */ "http://purl.obolibrary.org/obo/PATO_0001354",
            /* transparent */ "http://purl.obolibrary.org/obo/PATO_0000964",
            /* opaque */ "http://purl.obolibrary.org/obo/PATO_0000963",
            /* color brightness */ "http://purl.obolibrary.org/obo/PATO_0000016",
            /* high brightness */ "http://purl.obolibrary.org/obo/PATO_0000665",
            /* low brightness */ "http://purl.obolibrary.org/obo/PATO_0000327",
            /* color */ "http://purl.obolibrary.org/obo/PATO_0000014",
            /* brown */ "http://purl.obolibrary.org/obo/PATO_0000952",
            /* yellow brown */ "http://purl.obolibrary.org/obo/PATO_0002411",
            /* light yellow brown */ "http://purl.obolibrary.org/obo/PATO_0002413",
            /* beige */ "http://purl.obolibrary.org/obo/PATO_0002410",
            /* dark yellow brown */ "http://purl.obolibrary.org/obo/PATO_0002412",
            /* bronze */ "http://purl.obolibrary.org/obo/PATO_0002363",
            /* red brown */ "http://purl.obolibrary.org/obo/PATO_0001287",
            /* light red brown */ "http://purl.obolibrary.org/obo/PATO_0001289",
            /* dark red brown */ "http://purl.obolibrary.org/obo/PATO_0001288",
            /* desaturated brown */ "http://purl.obolibrary.org/obo/PATO_0001268",
            /* saturated brown */ "http://purl.obolibrary.org/obo/PATO_0001267",
            /* light brown */ "http://purl.obolibrary.org/obo/PATO_0001246",
            /* dark brown */ "http://purl.obolibrary.org/obo/PATO_0001245",
            /* purple */ "http://purl.obolibrary.org/obo/PATO_0000951",
            /* lilac */ "http://purl.obolibrary.org/obo/PATO_0001943",
            /* desaturated purple */ "http://purl.obolibrary.org/obo/PATO_0001282",
            /* saturated purple */ "http://purl.obolibrary.org/obo/PATO_0001281",
            /* light purple */ "http://purl.obolibrary.org/obo/PATO_0001260",
            /* dark purple */ "http://purl.obolibrary.org/obo/PATO_0001259",
            /* orange */ "http://purl.obolibrary.org/obo/PATO_0000953",
            /* yellow orange */ "http://purl.obolibrary.org/obo/PATO_0001944",
            /* ochre */ "http://purl.obolibrary.org/obo/PATO_0001945",
            /* desaturated orange */ "http://purl.obolibrary.org/obo/PATO_0001278",
            /* saturated orange */ "http://purl.obolibrary.org/obo/PATO_0001277",
            /* dark orange */ "http://purl.obolibrary.org/obo/PATO_0001256",
            /* light orange */ "http://purl.obolibrary.org/obo/PATO_0001255",
            /* white */ "http://purl.obolibrary.org/obo/PATO_0000323",
            /* magenta */ "http://purl.obolibrary.org/obo/PATO_0000321",
            /* desaturated magenta */ "http://purl.obolibrary.org/obo/PATO_0001276",
            /* saturated magenta */ "http://purl.obolibrary.org/obo/PATO_0001275",
            /* dark magenta */ "http://purl.obolibrary.org/obo/PATO_0001254",
            /* light magenta */ "http://purl.obolibrary.org/obo/PATO_0001253",
            /* cyan */ "http://purl.obolibrary.org/obo/PATO_0000319",
            /* desaturated cyan */ "http://purl.obolibrary.org/obo/PATO_0001270",
            /* saturated cyan */ "http://purl.obolibrary.org/obo/PATO_0001269",
            /* dark cyan */ "http://purl.obolibrary.org/obo/PATO_0001248",
            /* light cyan */ "http://purl.obolibrary.org/obo/PATO_0001247",
            /* green */ "http://purl.obolibrary.org/obo/PATO_0000320",
            /* brown green */ "http://purl.obolibrary.org/obo/PATO_0001942",
            /* yellow green */ "http://purl.obolibrary.org/obo/PATO_0001941",
            /* desaturated green */ "http://purl.obolibrary.org/obo/PATO_0001272",
            /* saturated green */ "http://purl.obolibrary.org/obo/PATO_0001271",
            /* light green */ "http://purl.obolibrary.org/obo/PATO_0001250",
            /* dark green */ "http://purl.obolibrary.org/obo/PATO_0001249",
            /* rosy */ "http://purl.obolibrary.org/obo/PATO_0001425",
            /* yellow */ "http://purl.obolibrary.org/obo/PATO_0000324",
            /* desaturated yellow */ "http://purl.obolibrary.org/obo/PATO_0001286",
            /* saturated yellow */ "http://purl.obolibrary.org/obo/PATO_0001285",
            /* light yellow */ "http://purl.obolibrary.org/obo/PATO_0001264",
            /* dark yellow */ "http://purl.obolibrary.org/obo/PATO_0001263",
            /* red */ "http://purl.obolibrary.org/obo/PATO_0000322",
            /* desaturated red */ "http://purl.obolibrary.org/obo/PATO_0001284",
            /* saturated red */ "http://purl.obolibrary.org/obo/PATO_0001283",
            /* light red */ "http://purl.obolibrary.org/obo/PATO_0001262",
            /* dark red */ "http://purl.obolibrary.org/obo/PATO_0001261",
            /* pink */ "http://purl.obolibrary.org/obo/PATO_0000954",
            /* pale pink */ "http://purl.obolibrary.org/obo/PATO_0002020",
            /* dark pale pink */ "http://purl.obolibrary.org/obo/PATO_0001280",
            /* deep pink */ "http://purl.obolibrary.org/obo/PATO_0001258",
            /* light deep pink */ "http://purl.obolibrary.org/obo/PATO_0001257",
            /* violet */ "http://purl.obolibrary.org/obo/PATO_0001424",
            /* dark violet */ "http://purl.obolibrary.org/obo/PATO_0001705",
            /* light violet */ "http://purl.obolibrary.org/obo/PATO_0001704",
            /* desaturated violet */ "http://purl.obolibrary.org/obo/PATO_0001703",
            /* saturated violet */ "http://purl.obolibrary.org/obo/PATO_0001702",
            /* blue */ "http://purl.obolibrary.org/obo/PATO_0000318",
            /* desaturated blue */ "http://purl.obolibrary.org/obo/PATO_0001266",
            /* saturated blue */ "http://purl.obolibrary.org/obo/PATO_0001265",
            /* dark blue */ "http://purl.obolibrary.org/obo/PATO_0001244",
            /* light blue */ "http://purl.obolibrary.org/obo/PATO_0001243",
            /* colored */ "http://purl.obolibrary.org/obo/PATO_0000336",
            /* maroon */ "http://purl.obolibrary.org/obo/PATO_0001426",
            /* grey */ "http://purl.obolibrary.org/obo/PATO_0000950",
            /* light grey */ "http://purl.obolibrary.org/obo/PATO_0001252",
            /* dark grey */ "http://purl.obolibrary.org/obo/PATO_0001251",
            /* black */ "http://purl.obolibrary.org/obo/PATO_0000317",
            /* vermilion */ "http://purl.obolibrary.org/obo/PATO_0001302",
            /* full-spectrum EM radiation quality */ "http://purl.obolibrary.org/obo/PATO_0001292",
            /* radiation emitting quality */ "http://purl.obolibrary.org/obo/PATO_0001299",
            /* emmision wavelength */ "http://purl.obolibrary.org/obo/PATO_0002059",
            /* radiation emitting intensity quality */ "http://purl.obolibrary.org/obo/PATO_0001717",
            /* radiation reflective quality */ "http://purl.obolibrary.org/obo/PATO_0001294",
            /* side scatter */ "http://purl.obolibrary.org/obo/PATO_0002033",
            /* forward scatter */ "http://purl.obolibrary.org/obo/PATO_0002032",
            /* refractivity */ "http://purl.obolibrary.org/obo/PATO_0001372",
            /* refractile */ "http://purl.obolibrary.org/obo/PATO_0001423",
            /* reflectivity */ "http://purl.obolibrary.org/obo/PATO_0001297",
            /* glistening */ "http://purl.obolibrary.org/obo/PATO_0001373",
            /* albedo */ "http://purl.obolibrary.org/obo/PATO_0001295",
            /* absorption quality */ "http://purl.obolibrary.org/obo/PATO_0001293",
            /* absorption wavelength */ "http://purl.obolibrary.org/obo/PATO_0002060",
            /* vaporizability */ "http://purl.obolibrary.org/obo/PATO_0002200",
            /* activity (of a radionuclide) */ "http://purl.obolibrary.org/obo/PATO_0001740",
            /* radioactive */ "http://purl.obolibrary.org/obo/PATO_0001741",
            /* decreased radioactivity */ "http://purl.obolibrary.org/obo/PATO_0001743",
            /* odor */ "http://purl.obolibrary.org/obo/PATO_0000058",
            /* odorous */ "http://purl.obolibrary.org/obo/PATO_0001331",
            /* decreased odor */ "http://purl.obolibrary.org/obo/PATO_0001892",
            /* odorless */ "http://purl.obolibrary.org/obo/PATO_0001208",
            /* sound quality */ "http://purl.obolibrary.org/obo/PATO_0001519",
            /* loud */ "http://purl.obolibrary.org/obo/PATO_0001528",
            /* quiet */ "http://purl.obolibrary.org/obo/PATO_0001527",
            /* mass */ "http://purl.obolibrary.org/obo/PATO_0000125",
            /* molar mass */ "http://purl.obolibrary.org/obo/PATO_0001681",
            /* decreased mass */ "http://purl.obolibrary.org/obo/PATO_0001562",
            /* radiation quality */ "http://purl.obolibrary.org/obo/PATO_0001739",
            /* electric potential */ "http://purl.obolibrary.org/obo/PATO_0001464",
            /* action potential */ "http://purl.obolibrary.org/obo/PATO_0001463",
            /* membrane potential */ "http://purl.obolibrary.org/obo/PATO_0001462",
            /* flavor */ "http://purl.obolibrary.org/obo/PATO_0000043",
            /* flavourless */ "http://purl.obolibrary.org/obo/PATO_0001330",
            /* flavourful */ "http://purl.obolibrary.org/obo/PATO_0001329",
            /* bitter */ "http://purl.obolibrary.org/obo/PATO_0002474",
            /* functionality */ "http://purl.obolibrary.org/obo/PATO_0001509",
            /* activation quality */ "http://purl.obolibrary.org/obo/PATO_0002353",
            /* inactive */ "http://purl.obolibrary.org/obo/PATO_0002355",
            /* active */ "http://purl.obolibrary.org/obo/PATO_0002354",
            /* constitutively active */ "http://purl.obolibrary.org/obo/PATO_0002356",
            /* necessity (continuant) */ "http://purl.obolibrary.org/obo/PATO_0001634",
            /* unnecessary (continuant) */ "http://purl.obolibrary.org/obo/PATO_0001636",
            /* necessary (continuant) */ "http://purl.obolibrary.org/obo/PATO_0001635",
            /* sufficiency */ "http://purl.obolibrary.org/obo/PATO_0001626",
            /* insufficient */ "http://purl.obolibrary.org/obo/PATO_0001628",
            /* sufficient */ "http://purl.obolibrary.org/obo/PATO_0001627",
            /* non-functional */ "http://purl.obolibrary.org/obo/PATO_0001511",
            /* functional */ "http://purl.obolibrary.org/obo/PATO_0001510",
            /* decreased functionality */ "http://purl.obolibrary.org/obo/PATO_0001624",
            /* aplastic/hypoplastic */ "http://purl.obolibrary.org/obo/PATO_0002290",
            /* aplastic */ "http://purl.obolibrary.org/obo/PATO_0001483",
            /* hypoplastic */ "http://purl.obolibrary.org/obo/PATO_0000645",
            /* morphology */ "http://purl.obolibrary.org/obo/PATO_0000051",
            /* differentiated from */ "http://purl.obolibrary.org/obo/PATO_0005006",
            /* texture */ "http://purl.obolibrary.org/obo/PATO_0000150",
            /* foveate */ "http://purl.obolibrary.org/obo/PATO_0002296",
            /* looseness */ "http://purl.obolibrary.org/obo/PATO_0002010",
            /* pilosity */ "http://purl.obolibrary.org/obo/PATO_0000066",
            /* hairy */ "http://purl.obolibrary.org/obo/PATO_0000454",
            /* arachnose */ "http://purl.obolibrary.org/obo/PATO_0002344",
            /* tomentose */ "http://purl.obolibrary.org/obo/PATO_0002341",
            /* hispid */ "http://purl.obolibrary.org/obo/PATO_0002339",
            /* hispidulous */ "http://purl.obolibrary.org/obo/PATO_0002340",
            /* setose */ "http://purl.obolibrary.org/obo/PATO_0002289",
            /* pubescent hair */ "http://purl.obolibrary.org/obo/PATO_0001320",
            /* glabrous */ "http://purl.obolibrary.org/obo/PATO_0000453",
            /* flaky */ "http://purl.obolibrary.org/obo/PATO_0001805",
            /* abrased */ "http://purl.obolibrary.org/obo/PATO_0001849",
            /* scaly */ "http://purl.obolibrary.org/obo/PATO_0001804",
            /* coating */ "http://purl.obolibrary.org/obo/PATO_0002012",
            /* greasy */ "http://purl.obolibrary.org/obo/PATO_0001606",
            /* viscid */ "http://purl.obolibrary.org/obo/PATO_0001370",
            /* rusty */ "http://purl.obolibrary.org/obo/PATO_0070059",
            /* grooved */ "http://purl.obolibrary.org/obo/PATO_0002255",
            /* furrowed */ "http://purl.obolibrary.org/obo/PATO_0001948",
            /* wrinkled */ "http://purl.obolibrary.org/obo/PATO_0001810",
            /* scrobiculate */ "http://purl.obolibrary.org/obo/PATO_0002294",
            /* smooth */ "http://purl.obolibrary.org/obo/PATO_0000701",
            /* plush */ "http://purl.obolibrary.org/obo/PATO_0040004",
            /* blistered */ "http://purl.obolibrary.org/obo/PATO_0001928",
            /* rough */ "http://purl.obolibrary.org/obo/PATO_0000700",
            /* warty */ "http://purl.obolibrary.org/obo/PATO_0001361",
            /* deformed */ "http://purl.obolibrary.org/obo/PATO_0001617",
            /* malformed */ "http://purl.obolibrary.org/obo/PATO_0000646",
            /* ventro-posteriorized */ "http://purl.obolibrary.org/obo/PATO_0030011",
            /* wholly ventro-posteriorized */ "http://purl.obolibrary.org/obo/PATO_0030013",
            /* partially ventro-posteriorized */ "http://purl.obolibrary.org/obo/PATO_0030012",
            /* dorso-anteriorized */ "http://purl.obolibrary.org/obo/PATO_0030008",
            /* wholly dorso-anteriorized */ "http://purl.obolibrary.org/obo/PATO_0030010",
            /* partially dorso-anteriorized */ "http://purl.obolibrary.org/obo/PATO_0030009",
            /* ventralized */ "http://purl.obolibrary.org/obo/PATO_0030003",
            /* partially ventralized */ "http://purl.obolibrary.org/obo/PATO_0030007",
            /* wholly ventralized */ "http://purl.obolibrary.org/obo/PATO_0000636",
            /* posteriorized */ "http://purl.obolibrary.org/obo/PATO_0030002",
            /* partially posteriorized */ "http://purl.obolibrary.org/obo/PATO_0030006",
            /* wholly posteriorized */ "http://purl.obolibrary.org/obo/PATO_0000630",
            /* dorsalized */ "http://purl.obolibrary.org/obo/PATO_0030001",
            /* partially dorsalized */ "http://purl.obolibrary.org/obo/PATO_0030005",
            /* wholly dorsalized */ "http://purl.obolibrary.org/obo/PATO_0000620",
            /* anteriorized */ "http://purl.obolibrary.org/obo/PATO_0030000",
            /* partially anteriorized */ "http://purl.obolibrary.org/obo/PATO_0030004",
            /* wholly anteriorized */ "http://purl.obolibrary.org/obo/PATO_0000615",
            /* monstrous */ "http://purl.obolibrary.org/obo/PATO_0001465",
            /* erythematous */ "http://purl.obolibrary.org/obo/PATO_0040048",
            /* anaplastic */ "http://purl.obolibrary.org/obo/PATO_0002092",
            /* closure */ "http://purl.obolibrary.org/obo/PATO_0000136",
            /* open */ "http://purl.obolibrary.org/obo/PATO_0000610",
            /* closure incomplete */ "http://purl.obolibrary.org/obo/PATO_0000609",
            /* closed */ "http://purl.obolibrary.org/obo/PATO_0000608",
            /* atretic */ "http://purl.obolibrary.org/obo/PATO_0001819",
            /* obstructed */ "http://purl.obolibrary.org/obo/PATO_0000648",
            /* amorphous */ "http://purl.obolibrary.org/obo/PATO_0001332",
            /* transformed to */ "http://purl.obolibrary.org/obo/PATO_0002470",
            /* adenomatous */ "http://purl.obolibrary.org/obo/PATO_0002090",
            /* shape */ "http://purl.obolibrary.org/obo/PATO_0000052",
            /* arrow-shaped */ "http://purl.obolibrary.org/obo/PATO_0001881",
            /* peg-like */ "http://purl.obolibrary.org/obo/PATO_0002535",
            /* anchor-shaped */ "http://purl.obolibrary.org/obo/PATO_0002446",
            /* plug shaped */ "http://purl.obolibrary.org/obo/PATO_0040012",
            /* sloped */ "http://purl.obolibrary.org/obo/PATO_0001481",
            /* sloped downward */ "http://purl.obolibrary.org/obo/PATO_0002143",
            /* coiling */ "http://purl.obolibrary.org/obo/PATO_0001794",
            /* uncoiled */ "http://purl.obolibrary.org/obo/PATO_0000415",
            /* coiled */ "http://purl.obolibrary.org/obo/PATO_0000404",
            /* decreased coiling */ "http://purl.obolibrary.org/obo/PATO_0001796",
            /* digitate */ "http://purl.obolibrary.org/obo/PATO_0001980",
            /* bent */ "http://purl.obolibrary.org/obo/PATO_0000617",
            /* kinked */ "http://purl.obolibrary.org/obo/PATO_0001798",
            /* flared */ "http://purl.obolibrary.org/obo/PATO_0040022",
            /* pinnate */ "http://purl.obolibrary.org/obo/PATO_0000410",
            /* sigmoid */ "http://purl.obolibrary.org/obo/PATO_0001878",
            /* falciform */ "http://purl.obolibrary.org/obo/PATO_0002215",
            /* concavity */ "http://purl.obolibrary.org/obo/PATO_0002005",
            /* concavo-convex */ "http://purl.obolibrary.org/obo/PATO_0002538",
            /* concave */ "http://purl.obolibrary.org/obo/PATO_0001857",
            /* biconcave */ "http://purl.obolibrary.org/obo/PATO_0002039",
            /* invaginated */ "http://purl.obolibrary.org/obo/PATO_0001748",
            /* arched */ "http://purl.obolibrary.org/obo/PATO_0001594",
            /* notched */ "http://purl.obolibrary.org/obo/PATO_0001495",
            /* emarginate */ "http://purl.obolibrary.org/obo/PATO_0002234",
            /* cleft */ "http://purl.obolibrary.org/obo/PATO_0000403",
            /* convex */ "http://purl.obolibrary.org/obo/PATO_0001355",
            /* biconvex */ "http://purl.obolibrary.org/obo/PATO_0002040",
            /* fimbriated */ "http://purl.obolibrary.org/obo/PATO_0002311",
            /* cicatricial */ "http://purl.obolibrary.org/obo/PATO_0002421",
            /* torsioned */ "http://purl.obolibrary.org/obo/PATO_0002445",
            /* band shaped */ "http://purl.obolibrary.org/obo/PATO_0040014",
            /* scute-like */ "http://purl.obolibrary.org/obo/PATO_0002520",
            /* ruffled */ "http://purl.obolibrary.org/obo/PATO_0001799",
            /* pleomorphic */ "http://purl.obolibrary.org/obo/PATO_0001356",
            /* x-shaped */ "http://purl.obolibrary.org/obo/PATO_0002429",
            /* strap-shaped */ "http://purl.obolibrary.org/obo/PATO_0002430",
            /* crown like */ "http://purl.obolibrary.org/obo/PATO_0040006",
            /* radiating */ "http://purl.obolibrary.org/obo/PATO_0005005",
            /* hyponastic */ "http://purl.obolibrary.org/obo/PATO_0002329",
            /* undivided */ "http://purl.obolibrary.org/obo/PATO_0002034",
            /* slender */ "http://purl.obolibrary.org/obo/PATO_0002212",
            /* tendrilous */ "http://purl.obolibrary.org/obo/PATO_0015005",
            /* elongated */ "http://purl.obolibrary.org/obo/PATO_0001154",
            /* telescopic */ "http://purl.obolibrary.org/obo/PATO_0002313",
            /* papillomatous */ "http://purl.obolibrary.org/obo/PATO_0002423",
            /* branchiness */ "http://purl.obolibrary.org/obo/PATO_0002009",
            /* unbranched */ "http://purl.obolibrary.org/obo/PATO_0000414",
            /* branched */ "http://purl.obolibrary.org/obo/PATO_0000402",
            /* tripartite */ "http://purl.obolibrary.org/obo/PATO_0001890",
            /* Y-shaped */ "http://purl.obolibrary.org/obo/PATO_0001201",
            /* T-shaped */ "http://purl.obolibrary.org/obo/PATO_0001200",
            /* increased branchiness */ "http://purl.obolibrary.org/obo/PATO_0002285",
            /* dendritic */ "http://purl.obolibrary.org/obo/PATO_0002045",
            /* bifurcated */ "http://purl.obolibrary.org/obo/PATO_0001784",
            /* brochidodromous */ "http://purl.obolibrary.org/obo/PATO_0001970",
            /* reticulodromous */ "http://purl.obolibrary.org/obo/PATO_0001972",
            /* actinodromous */ "http://purl.obolibrary.org/obo/PATO_0001967",
            /* decreased branchiness */ "http://purl.obolibrary.org/obo/PATO_0002286",
            /* cladodromous */ "http://purl.obolibrary.org/obo/PATO_0001971",
            /* craspedodromous */ "http://purl.obolibrary.org/obo/PATO_0001969",
            /* parallelodromous */ "http://purl.obolibrary.org/obo/PATO_0001968",
            /* plumose */ "http://purl.obolibrary.org/obo/PATO_0005010",
            /* scaphoid */ "http://purl.obolibrary.org/obo/PATO_0002426",
            /* 3-D shape */ "http://purl.obolibrary.org/obo/PATO_0002266",
            /* molariform */ "http://purl.obolibrary.org/obo/PATO_0005009",
            /* brachydont */ "http://purl.obolibrary.org/obo/PATO_0005008",
            /* hypsodont */ "http://purl.obolibrary.org/obo/PATO_0005007",
            /* concave 3-D shape */ "http://purl.obolibrary.org/obo/PATO_0002008",
            /* obclavate */ "http://purl.obolibrary.org/obo/PATO_0002213",
            /* J-shaped */ "http://purl.obolibrary.org/obo/PATO_0015020",
            /* anvil */ "http://purl.obolibrary.org/obo/PATO_0002386",
            /* trumpet-shaped */ "http://purl.obolibrary.org/obo/PATO_0002375",
            /* trough shaped */ "http://purl.obolibrary.org/obo/PATO_0040015",
            /* U-shaped */ "http://purl.obolibrary.org/obo/PATO_0001879",
            /* bowl shaped */ "http://purl.obolibrary.org/obo/PATO_0040009",
            /* limaciform */ "http://purl.obolibrary.org/obo/PATO_0001882",
            /* cup-shaped */ "http://purl.obolibrary.org/obo/PATO_0002227",
            /* cupulate */ "http://purl.obolibrary.org/obo/PATO_0002342",
            /* reniform */ "http://purl.obolibrary.org/obo/PATO_0001871",
            /* dumbbell-shaped */ "http://purl.obolibrary.org/obo/PATO_0001876",
            /* bicornuate */ "http://purl.obolibrary.org/obo/PATO_0002161",
            /* spoon-shaped */ "http://purl.obolibrary.org/obo/PATO_0002208",
            /* clavate */ "http://purl.obolibrary.org/obo/PATO_0001883",
            /* villiform */ "http://purl.obolibrary.org/obo/PATO_0002022",
            /* multipartite */ "http://purl.obolibrary.org/obo/PATO_0002510",
            /* bipartite */ "http://purl.obolibrary.org/obo/PATO_0002533",
            /* quadripartite */ "http://purl.obolibrary.org/obo/PATO_0002447",
            /* snowman-shaped */ "http://purl.obolibrary.org/obo/PATO_0002346",
            /* heart shaped */ "http://purl.obolibrary.org/obo/PATO_0000948",
            /* hourglass-shaped */ "http://purl.obolibrary.org/obo/PATO_0002239",
            /* C-shaped */ "http://purl.obolibrary.org/obo/PATO_0015019",
            /* convex 3-D shape */ "http://purl.obolibrary.org/obo/PATO_0002007",
            /* hemispheroid */ "http://purl.obolibrary.org/obo/PATO_0005000",
            /* spindle-shaped */ "http://purl.obolibrary.org/obo/PATO_0001409",
            /* fusiform */ "http://purl.obolibrary.org/obo/PATO_0002400",
            /* bullet-shaped */ "http://purl.obolibrary.org/obo/PATO_0002087",
            /* teardrop-shaped */ "http://purl.obolibrary.org/obo/PATO_0002240",
            /* spheroid */ "http://purl.obolibrary.org/obo/PATO_0001865",
            /* subspherical */ "http://purl.obolibrary.org/obo/PATO_0005014",
            /* subovoid */ "http://purl.obolibrary.org/obo/PATO_0002537",
            /* bulbous */ "http://purl.obolibrary.org/obo/PATO_0002210",
            /* ovate */ "http://purl.obolibrary.org/obo/PATO_0001891",
            /* obovate */ "http://purl.obolibrary.org/obo/PATO_0001936",
            /* prolate */ "http://purl.obolibrary.org/obo/PATO_0001866",
            /* spherical */ "http://purl.obolibrary.org/obo/PATO_0001499",
            /* oblate */ "http://purl.obolibrary.org/obo/PATO_0000409",
            /* conical */ "http://purl.obolibrary.org/obo/PATO_0002021",
            /* obconical */ "http://purl.obolibrary.org/obo/PATO_0002347",
            /* lemon-shaped */ "http://purl.obolibrary.org/obo/PATO_0002345",
            /* pulvinate */ "http://purl.obolibrary.org/obo/PATO_0001357",
            /* pear shaped */ "http://purl.obolibrary.org/obo/PATO_0005002",
            /* D-shaped */ "http://purl.obolibrary.org/obo/PATO_0002357",
            /* cylindrical */ "http://purl.obolibrary.org/obo/PATO_0001873",
            /* tubular */ "http://purl.obolibrary.org/obo/PATO_0002299",
            /* subcylindrical */ "http://purl.obolibrary.org/obo/PATO_0002226",
            /* fiber shaped */ "http://purl.obolibrary.org/obo/PATO_0002309",
            /* discoid */ "http://purl.obolibrary.org/obo/PATO_0001874",
            /* cuboid */ "http://purl.obolibrary.org/obo/PATO_0001872",
            /* surface feature shape */ "http://purl.obolibrary.org/obo/PATO_0001925",
            /* botryoidal */ "http://purl.obolibrary.org/obo/PATO_0001907",
            /* rugose */ "http://purl.obolibrary.org/obo/PATO_0001359",
            /* areolate */ "http://purl.obolibrary.org/obo/PATO_0002295",
            /* knobbled */ "http://purl.obolibrary.org/obo/PATO_0002427",
            /* sculpted surface */ "http://purl.obolibrary.org/obo/PATO_0002433",
            /* spiny */ "http://purl.obolibrary.org/obo/PATO_0001365",
            /* alobate */ "http://purl.obolibrary.org/obo/PATO_0002506",
            /* ornamentation */ "http://purl.obolibrary.org/obo/PATO_0002440",
            /* unornamented */ "http://purl.obolibrary.org/obo/PATO_0002442",
            /* ornamented */ "http://purl.obolibrary.org/obo/PATO_0002441",
            /* lobate */ "http://purl.obolibrary.org/obo/PATO_0001367",
            /* folded */ "http://purl.obolibrary.org/obo/PATO_0001910",
            /* keel-shaped */ "http://purl.obolibrary.org/obo/PATO_0002522",
            /* curled */ "http://purl.obolibrary.org/obo/PATO_0000405",
            /* convolute */ "http://purl.obolibrary.org/obo/PATO_0001966",
            /* reclinate */ "http://purl.obolibrary.org/obo/PATO_0001965",
            /* circinate */ "http://purl.obolibrary.org/obo/PATO_0001964",
            /* revolute */ "http://purl.obolibrary.org/obo/PATO_0001963",
            /* involute */ "http://purl.obolibrary.org/obo/PATO_0001962",
            /* slit-like */ "http://purl.obolibrary.org/obo/PATO_0002482",
            /* paddle shaped */ "http://purl.obolibrary.org/obo/PATO_0005003",
            /* waisted */ "http://purl.obolibrary.org/obo/PATO_0002431",
            /* keyhole shaped */ "http://purl.obolibrary.org/obo/PATO_0002466",
            /* auriculate */ "http://purl.obolibrary.org/obo/PATO_0001981",
            /* decurrent */ "http://purl.obolibrary.org/obo/PATO_0001984",
            /* corymb-like */ "http://purl.obolibrary.org/obo/PATO_0002455",
            /* leaf-like */ "http://purl.obolibrary.org/obo/PATO_0002457",
            /* punctiform */ "http://purl.obolibrary.org/obo/PATO_0001366",
            /* linear */ "http://purl.obolibrary.org/obo/PATO_0001199",
            /* oblanceolate */ "http://purl.obolibrary.org/obo/PATO_0002330",
            /* subulate */ "http://purl.obolibrary.org/obo/PATO_0001954",
            /* stubby */ "http://purl.obolibrary.org/obo/PATO_0001643",
            /* curvature */ "http://purl.obolibrary.org/obo/PATO_0001591",
            /* flattened */ "http://purl.obolibrary.org/obo/PATO_0002254",
            /* mesiodistally compressed */ "http://purl.obolibrary.org/obo/PATO_0005018",
            /* labiolingually compressed */ "http://purl.obolibrary.org/obo/PATO_0005017",
            /* antero-posteriorly flattened */ "http://purl.obolibrary.org/obo/PATO_0002252",
            /* laterally compressed */ "http://purl.obolibrary.org/obo/PATO_0002054",
            /* dorso-ventrally flattened */ "http://purl.obolibrary.org/obo/PATO_0002053",
            /* flat */ "http://purl.obolibrary.org/obo/PATO_0000407",
            /* platelike */ "http://purl.obolibrary.org/obo/PATO_0002253",
            /* sinuous */ "http://purl.obolibrary.org/obo/PATO_0001989",
            /* rotational curvature */ "http://purl.obolibrary.org/obo/PATO_0001787",
            /* lateral and rotional curvature */ "http://purl.obolibrary.org/obo/PATO_0002049",
            /* curved */ "http://purl.obolibrary.org/obo/PATO_0000406",
            /* antrorse */ "http://purl.obolibrary.org/obo/PATO_0002238",
            /* curved lingually */ "http://purl.obolibrary.org/obo/PATO_0005019",
            /* asymmetrically curved */ "http://purl.obolibrary.org/obo/PATO_0001848",
            /* curved ventral */ "http://purl.obolibrary.org/obo/PATO_0001469",
            /* domed */ "http://purl.obolibrary.org/obo/PATO_0001789",
            /* upturned */ "http://purl.obolibrary.org/obo/PATO_0002031",
            /* curved medial */ "http://purl.obolibrary.org/obo/PATO_0002164",
            /* splayed */ "http://purl.obolibrary.org/obo/PATO_0001785",
            /* splayed ventral */ "http://purl.obolibrary.org/obo/PATO_0002154",
            /* splayed rostral */ "http://purl.obolibrary.org/obo/PATO_0002153",
            /* splayed lateral */ "http://purl.obolibrary.org/obo/PATO_0002152",
            /* splayed dorsal */ "http://purl.obolibrary.org/obo/PATO_0002151",
            /* splayed caudal */ "http://purl.obolibrary.org/obo/PATO_0002150",
            /* curved lateral */ "http://purl.obolibrary.org/obo/PATO_0001649",
            /* retrorse */ "http://purl.obolibrary.org/obo/PATO_0002237",
            /* decreased curvature */ "http://purl.obolibrary.org/obo/PATO_0001593",
            /* curved rostral */ "http://purl.obolibrary.org/obo/PATO_0001466",
            /* recurved */ "http://purl.obolibrary.org/obo/PATO_0002211",
            /* curved dorsal */ "http://purl.obolibrary.org/obo/PATO_0001468",
            /* curved caudal */ "http://purl.obolibrary.org/obo/PATO_0001467",
            /* angular */ "http://purl.obolibrary.org/obo/PATO_0001977",
            /* inverted-V shaped */ "http://purl.obolibrary.org/obo/PATO_0002484",
            /* w-shaped */ "http://purl.obolibrary.org/obo/PATO_0002439",
            /* L-shaped */ "http://purl.obolibrary.org/obo/PATO_0002225",
            /* V-shaped */ "http://purl.obolibrary.org/obo/PATO_0002224",
            /* straight */ "http://purl.obolibrary.org/obo/PATO_0002180",
            /* shortened */ "http://purl.obolibrary.org/obo/PATO_0002364",
            /* brush-like shape */ "http://purl.obolibrary.org/obo/PATO_0002315",
            /* pyramidal */ "http://purl.obolibrary.org/obo/PATO_0002336",
            /* bracket */ "http://purl.obolibrary.org/obo/PATO_0002142",
            /* undulate */ "http://purl.obolibrary.org/obo/PATO_0000967",
            /* spade-shaped */ "http://purl.obolibrary.org/obo/PATO_0002432",
            /* filamentous */ "http://purl.obolibrary.org/obo/PATO_0001360",
            /* prism shaped */ "http://purl.obolibrary.org/obo/PATO_0002465",
            /* inflorescence-like */ "http://purl.obolibrary.org/obo/PATO_0002456",
            /* striated */ "http://purl.obolibrary.org/obo/PATO_0001410",
            /* obliquely striated */ "http://purl.obolibrary.org/obo/PATO_0002479",
            /* transversely striated */ "http://purl.obolibrary.org/obo/PATO_0002478",
            /* striate-angular */ "http://purl.obolibrary.org/obo/PATO_0002385",
            /* ridged */ "http://purl.obolibrary.org/obo/PATO_0002358",
            /* sepal-like */ "http://purl.obolibrary.org/obo/PATO_0002459",
            /* carpel-like */ "http://purl.obolibrary.org/obo/PATO_0002454",
            /* plowshare shaped */ "http://purl.obolibrary.org/obo/PATO_0002534",
            /* plume-shaped */ "http://purl.obolibrary.org/obo/PATO_0015030",
            /* shell shaped */ "http://purl.obolibrary.org/obo/PATO_0040007",
            /* lobed */ "http://purl.obolibrary.org/obo/PATO_0001979",
            /* acinus */ "http://purl.obolibrary.org/obo/PATO_0002378",
            /* trilobed */ "http://purl.obolibrary.org/obo/PATO_0002241",
            /* bilobed */ "http://purl.obolibrary.org/obo/PATO_0002214",
            /* spur shaped */ "http://purl.obolibrary.org/obo/PATO_0002540",
            /* stepped */ "http://purl.obolibrary.org/obo/PATO_0015016",
            /* edge shape */ "http://purl.obolibrary.org/obo/PATO_0002267",
            /* fringed */ "http://purl.obolibrary.org/obo/PATO_0040008",
            /* serration */ "http://purl.obolibrary.org/obo/PATO_0001976",
            /* unserrated */ "http://purl.obolibrary.org/obo/PATO_0001975",
            /* scalloped */ "http://purl.obolibrary.org/obo/PATO_0001889",
            /* erose */ "http://purl.obolibrary.org/obo/PATO_0001368",
            /* serrated */ "http://purl.obolibrary.org/obo/PATO_0001206",
            /* dentated */ "http://purl.obolibrary.org/obo/PATO_0001205",
            /* epinastic */ "http://purl.obolibrary.org/obo/PATO_0000945",
            /* tholiform */ "http://purl.obolibrary.org/obo/PATO_0002335",
            /* fasciated */ "http://purl.obolibrary.org/obo/PATO_0000949",
            /* columnar */ "http://purl.obolibrary.org/obo/PATO_0002063",
            /* drooping */ "http://purl.obolibrary.org/obo/PATO_0002165",
            /* wilty */ "http://purl.obolibrary.org/obo/PATO_0002461",
            /* lathlike */ "http://purl.obolibrary.org/obo/PATO_0002467",
            /* split */ "http://purl.obolibrary.org/obo/PATO_0001786",
            /* multifid */ "http://purl.obolibrary.org/obo/PATO_0002231",
            /* split radially */ "http://purl.obolibrary.org/obo/PATO_0002172",
            /* split bilaterally */ "http://purl.obolibrary.org/obo/PATO_0002171",
            /* split laterally */ "http://purl.obolibrary.org/obo/PATO_0002170",
            /* split medially */ "http://purl.obolibrary.org/obo/PATO_0002169",
            /* interdigitated */ "http://purl.obolibrary.org/obo/PATO_0001960",
            /* broad */ "http://purl.obolibrary.org/obo/PATO_0002359",
            /* boomerang shaped */ "http://purl.obolibrary.org/obo/PATO_0002536",
            /* pin-like */ "http://purl.obolibrary.org/obo/PATO_0002458",
            /* truncated */ "http://purl.obolibrary.org/obo/PATO_0000936",
            /* raised */ "http://purl.obolibrary.org/obo/PATO_0001369",
            /* saddle-shaped */ "http://purl.obolibrary.org/obo/PATO_0002517",
            /* tripodal */ "http://purl.obolibrary.org/obo/PATO_0002428",
            /* cut */ "http://purl.obolibrary.org/obo/PATO_0001978",
            /* triradiate */ "http://purl.obolibrary.org/obo/PATO_0002391",
            /* parallel-sided */ "http://purl.obolibrary.org/obo/PATO_0002485",
            /* sharpness */ "http://purl.obolibrary.org/obo/PATO_0000944",
            /* blunt */ "http://purl.obolibrary.org/obo/PATO_0001950",
            /* retuse */ "http://purl.obolibrary.org/obo/PATO_0001974",
            /* sharp */ "http://purl.obolibrary.org/obo/PATO_0001419",
            /* pointed */ "http://purl.obolibrary.org/obo/PATO_0002258",
            /* incisiform */ "http://purl.obolibrary.org/obo/PATO_0002209",
            /* pointleted */ "http://purl.obolibrary.org/obo/PATO_0001949",
            /* lanceolate */ "http://purl.obolibrary.org/obo/PATO_0001877",
            /* lanceolate-triangular */ "http://purl.obolibrary.org/obo/PATO_0002338",
            /* lance-ovate */ "http://purl.obolibrary.org/obo/PATO_0002337",
            /* tapered */ "http://purl.obolibrary.org/obo/PATO_0001500",
            /* blade-like */ "http://purl.obolibrary.org/obo/PATO_0002235",
            /* acuminate */ "http://purl.obolibrary.org/obo/PATO_0002228",
            /* attenuate */ "http://purl.obolibrary.org/obo/PATO_0001982",
            /* cuspidate */ "http://purl.obolibrary.org/obo/PATO_0001973",
            /* unicuspidate */ "http://purl.obolibrary.org/obo/PATO_0005021",
            /* multicuspidate */ "http://purl.obolibrary.org/obo/PATO_0002257",
            /* biscupidate */ "http://purl.obolibrary.org/obo/PATO_0002281",
            /* tricuspidate */ "http://purl.obolibrary.org/obo/PATO_0002256",
            /* star shaped */ "http://purl.obolibrary.org/obo/PATO_0002065",
            /* cane-like */ "http://purl.obolibrary.org/obo/PATO_0002511",
            /* robust */ "http://purl.obolibrary.org/obo/PATO_0002310",
            /* 2-D shape */ "http://purl.obolibrary.org/obo/PATO_0002006",
            /* triangular */ "http://purl.obolibrary.org/obo/PATO_0001875",
            /* scalene triangular */ "http://purl.obolibrary.org/obo/PATO_0002308",
            /* isosceles triangular */ "http://purl.obolibrary.org/obo/PATO_0002307",
            /* equilateral triangular */ "http://purl.obolibrary.org/obo/PATO_0002306",
            /* subtriangular */ "http://purl.obolibrary.org/obo/PATO_0002230",
            /* cuneate */ "http://purl.obolibrary.org/obo/PATO_0001955",
            /* superelliptic */ "http://purl.obolibrary.org/obo/PATO_0002318",
            /* subelliptical */ "http://purl.obolibrary.org/obo/PATO_0005004",
            /* hyperelliptic */ "http://purl.obolibrary.org/obo/PATO_0002322",
            /* hypoelliptic */ "http://purl.obolibrary.org/obo/PATO_0002321",
            /* diamond shaped */ "http://purl.obolibrary.org/obo/PATO_0002320",
            /* squircle */ "http://purl.obolibrary.org/obo/PATO_0002319",
            /* elliptic */ "http://purl.obolibrary.org/obo/PATO_0000947",
            /* ring shaped */ "http://purl.obolibrary.org/obo/PATO_0002539",
            /* circular */ "http://purl.obolibrary.org/obo/PATO_0000411",
            /* obtuse */ "http://purl.obolibrary.org/obo/PATO_0001935",
            /* orbicular */ "http://purl.obolibrary.org/obo/PATO_0001934",
            /* aliform */ "http://purl.obolibrary.org/obo/PATO_0002236",
            /* crescent-shaped */ "http://purl.obolibrary.org/obo/PATO_0001870",
            /* oblong */ "http://purl.obolibrary.org/obo/PATO_0000946",
            /* spatulate */ "http://purl.obolibrary.org/obo/PATO_0001937",
            /* subcircular */ "http://purl.obolibrary.org/obo/PATO_0002397",
            /* quadrangular */ "http://purl.obolibrary.org/obo/PATO_0001988",
            /* parallelogram */ "http://purl.obolibrary.org/obo/PATO_0002317",
            /* rhomboid */ "http://purl.obolibrary.org/obo/PATO_0001938",
            /* rectangular */ "http://purl.obolibrary.org/obo/PATO_0001867",
            /* square */ "http://purl.obolibrary.org/obo/PATO_0000413",
            /* subrectangular */ "http://purl.obolibrary.org/obo/PATO_0002229",
            /* trapezoid */ "http://purl.obolibrary.org/obo/PATO_0002044",
            /* hexagonal */ "http://purl.obolibrary.org/obo/PATO_0002509",
            /* pentagonal */ "http://purl.obolibrary.org/obo/PATO_0040016",
            /* semicircular */ "http://purl.obolibrary.org/obo/PATO_0002232",
            /* fan-shaped */ "http://purl.obolibrary.org/obo/PATO_0002219",
            /* irregularly shaped */ "http://purl.obolibrary.org/obo/PATO_0005020",
            /* funnel-shaped */ "http://purl.obolibrary.org/obo/PATO_0002521",
            /* structure */ "http://purl.obolibrary.org/obo/PATO_0000141",
            /* composition */ "http://purl.obolibrary.org/obo/PATO_0000025",
            /* cartilaginous */ "http://purl.obolibrary.org/obo/PATO_0001449",
            /* calcified */ "http://purl.obolibrary.org/obo/PATO_0001447",
            /* ossified */ "http://purl.obolibrary.org/obo/PATO_0001448",
            /* poorly ossified */ "http://purl.obolibrary.org/obo/PATO_0002480",
            /* biomaterial purity */ "http://purl.obolibrary.org/obo/PATO_0001339",
            /* hydropic */ "http://purl.obolibrary.org/obo/PATO_0002119",
            /* edematous */ "http://purl.obolibrary.org/obo/PATO_0001450",
            /* hydrocephalic */ "http://purl.obolibrary.org/obo/PATO_0001853",
            /* mixed */ "http://purl.obolibrary.org/obo/PATO_0002122",
            /* contamination */ "http://purl.obolibrary.org/obo/PATO_0015031",
            /* decreased contamination */ "http://purl.obolibrary.org/obo/PATO_0015033",
            /* fibrinoid */ "http://purl.obolibrary.org/obo/PATO_0002115",
            /* water composition */ "http://purl.obolibrary.org/obo/PATO_0001800",
            /* decreased water composition */ "http://purl.obolibrary.org/obo/PATO_0001801",
            /* fibrillary */ "http://purl.obolibrary.org/obo/PATO_0002134",
            /* serous */ "http://purl.obolibrary.org/obo/PATO_0002128",
            /* fibrinopurulent */ "http://purl.obolibrary.org/obo/PATO_0002116",
            /* waxiness */ "http://purl.obolibrary.org/obo/PATO_0002381",
            /* waxy */ "http://purl.obolibrary.org/obo/PATO_0002384",
            /* decreased waxiness */ "http://purl.obolibrary.org/obo/PATO_0002383",
            /* collagenous */ "http://purl.obolibrary.org/obo/PATO_0002462",
            /* granular */ "http://purl.obolibrary.org/obo/PATO_0001759",
            /* amylose composition */ "http://purl.obolibrary.org/obo/PATO_0001539",
            /* non-glutinous */ "http://purl.obolibrary.org/obo/PATO_0001541",
            /* glutinous */ "http://purl.obolibrary.org/obo/PATO_0001540",
            /* mineralized */ "http://purl.obolibrary.org/obo/PATO_0002444",
            /* fibrotic */ "http://purl.obolibrary.org/obo/PATO_0040019",
            /* watery */ "http://purl.obolibrary.org/obo/PATO_0002408",
            /* suppurative */ "http://purl.obolibrary.org/obo/PATO_0002120",
            /* fleshy */ "http://purl.obolibrary.org/obo/PATO_0002351",
            /* osseous */ "http://purl.obolibrary.org/obo/PATO_0002126",
            /* inflamed */ "http://purl.obolibrary.org/obo/PATO_0002104",
            /* mucoid */ "http://purl.obolibrary.org/obo/PATO_0001371",
            /* ligamentous */ "http://purl.obolibrary.org/obo/PATO_0002314",
            /* fatty */ "http://purl.obolibrary.org/obo/PATO_0002114",
            /* keratinous */ "http://purl.obolibrary.org/obo/PATO_0002507",
            /* brittle */ "http://purl.obolibrary.org/obo/PATO_0002477",
            /* constricted */ "http://purl.obolibrary.org/obo/PATO_0001847",
            /* overlap with */ "http://purl.obolibrary.org/obo/PATO_0001590",
            /* overlapped by */ "http://purl.obolibrary.org/obo/PATO_0002489",
            /* overlapping */ "http://purl.obolibrary.org/obo/PATO_0002488",
            /* wholeness */ "http://purl.obolibrary.org/obo/PATO_0001442",
            /* fragmented */ "http://purl.obolibrary.org/obo/PATO_0040047",
            /* whole */ "http://purl.obolibrary.org/obo/PATO_0001446",
            /* broken */ "http://purl.obolibrary.org/obo/PATO_0001444",
            /* ruptured */ "http://purl.obolibrary.org/obo/PATO_0040044",
            /* fractured */ "http://purl.obolibrary.org/obo/PATO_0040034",
            /* avulsion fractured */ "http://purl.obolibrary.org/obo/PATO_0040041",
            /* compression-fractured */ "http://purl.obolibrary.org/obo/PATO_0040040",
            /* spirally fractured */ "http://purl.obolibrary.org/obo/PATO_0040039",
            /* obliquely fractured */ "http://purl.obolibrary.org/obo/PATO_0040038",
            /* linearly fractured */ "http://purl.obolibrary.org/obo/PATO_0040037",
            /* impact-fractured */ "http://purl.obolibrary.org/obo/PATO_0040036",
            /* transversely fractured */ "http://purl.obolibrary.org/obo/PATO_0040035",
            /* partially broken */ "http://purl.obolibrary.org/obo/PATO_0002082",
            /* shattered */ "http://purl.obolibrary.org/obo/PATO_0002081",
            /* broken into two pieces */ "http://purl.obolibrary.org/obo/PATO_0002080",
            /* disassembled */ "http://purl.obolibrary.org/obo/PATO_0001445",
            /* shriveled */ "http://purl.obolibrary.org/obo/PATO_0002460",
            /* degeneration */ "http://purl.obolibrary.org/obo/PATO_0002037",
            /* non-degenerate */ "http://purl.obolibrary.org/obo/PATO_0002038",
            /* degenerate */ "http://purl.obolibrary.org/obo/PATO_0000639",
            /* absence due to degeneration */ "http://purl.obolibrary.org/obo/PATO_0015001",
            /* lysed */ "http://purl.obolibrary.org/obo/PATO_0065001",
            /* spongy */ "http://purl.obolibrary.org/obo/PATO_0001480",
            /* pedunculate */ "http://purl.obolibrary.org/obo/PATO_0002388",
            /* attachment quality */ "http://purl.obolibrary.org/obo/PATO_0001435",
            /* attached to */ "http://purl.obolibrary.org/obo/PATO_0001667",
            /* detached from */ "http://purl.obolibrary.org/obo/PATO_0001453",
            /* pedicellate */ "http://purl.obolibrary.org/obo/PATO_0001438",
            /* sessile (sensu zoology) */ "http://purl.obolibrary.org/obo/PATO_0001437",
            /* sessile (sensu botany) */ "http://purl.obolibrary.org/obo/PATO_0001436",
            /* fused with */ "http://purl.obolibrary.org/obo/PATO_0000642",
            /* disconnected */ "http://purl.obolibrary.org/obo/PATO_0010001",
            /* nodular */ "http://purl.obolibrary.org/obo/PATO_0002125",
            /* subdermal */ "http://purl.obolibrary.org/obo/PATO_0002438",
            /* cancellous */ "http://purl.obolibrary.org/obo/PATO_0002519",
            /* perfoliate */ "http://purl.obolibrary.org/obo/PATO_0001983",
            /* apoptotic */ "http://purl.obolibrary.org/obo/PATO_0000638",
            /* tangled */ "http://purl.obolibrary.org/obo/PATO_0001846",
            /* matted */ "http://purl.obolibrary.org/obo/PATO_0001607",
            /* dysplastic */ "http://purl.obolibrary.org/obo/PATO_0000640",
            /* neoplastic */ "http://purl.obolibrary.org/obo/PATO_0002011",
            /* neoplastic, spontaneous */ "http://purl.obolibrary.org/obo/PATO_0002473",
            /* neoplastic, non-invasive */ "http://purl.obolibrary.org/obo/PATO_0002132",
            /* neoplastic, invasive */ "http://purl.obolibrary.org/obo/PATO_0002129",
            /* neoplastic, deeply invasive */ "http://purl.obolibrary.org/obo/PATO_0002130",
            /* neoplastic, metastatic */ "http://purl.obolibrary.org/obo/PATO_0002098",
            /* neoplastic, malignant */ "http://purl.obolibrary.org/obo/PATO_0002097",
            /* neoplastic, non-malignant */ "http://purl.obolibrary.org/obo/PATO_0002096",
            /* acinar */ "http://purl.obolibrary.org/obo/PATO_0002422",
            /* friability */ "http://purl.obolibrary.org/obo/PATO_0002405",
            /* indurated */ "http://purl.obolibrary.org/obo/PATO_0002407",
            /* friable */ "http://purl.obolibrary.org/obo/PATO_0002406",
            /* accumulation */ "http://purl.obolibrary.org/obo/PATO_0002269",
            /* decreased accumulation */ "http://purl.obolibrary.org/obo/PATO_0002271",
            /* increased accumulation */ "http://purl.obolibrary.org/obo/PATO_0002270",
            /* maximally connected */ "http://purl.obolibrary.org/obo/PATO_0010000",
            /* fenestrated */ "http://purl.obolibrary.org/obo/PATO_0002064",
            /* laminar */ "http://purl.obolibrary.org/obo/PATO_0002124",
            /* structured */ "http://purl.obolibrary.org/obo/PATO_0001411",
            /* complete structure */ "http://purl.obolibrary.org/obo/PATO_0005012",
            /* eroding */ "http://purl.obolibrary.org/obo/PATO_0002453",
            /* distensibility */ "http://purl.obolibrary.org/obo/PATO_0015008",
            /* distensible */ "http://purl.obolibrary.org/obo/PATO_0002468",
            /* scarred */ "http://purl.obolibrary.org/obo/PATO_0001850",
            /* delaminated */ "http://purl.obolibrary.org/obo/PATO_0001514",
            /* fasciculation */ "http://purl.obolibrary.org/obo/PATO_0002013",
            /* defasciculated */ "http://purl.obolibrary.org/obo/PATO_0001959",
            /* fasciculated */ "http://purl.obolibrary.org/obo/PATO_0001861",
            /* necrotic */ "http://purl.obolibrary.org/obo/PATO_0000647",
            /* degree of pigmentation */ "http://purl.obolibrary.org/obo/PATO_0002247",
            /* unpigmented */ "http://purl.obolibrary.org/obo/PATO_0002249",
            /* pigmented */ "http://purl.obolibrary.org/obo/PATO_0002248",
            /* patchy pigmentation */ "http://purl.obolibrary.org/obo/PATO_0065002",
            /* decreased pigmentation */ "http://purl.obolibrary.org/obo/PATO_0002251",
            /* fragility */ "http://purl.obolibrary.org/obo/PATO_0001662",
            /* non-fragile */ "http://purl.obolibrary.org/obo/PATO_0001716",
            /* fragile */ "http://purl.obolibrary.org/obo/PATO_0001362",
            /* decreased fragility */ "http://purl.obolibrary.org/obo/PATO_0002056",
            /* polymeric */ "http://purl.obolibrary.org/obo/PATO_0015006",
            /* organization quality */ "http://purl.obolibrary.org/obo/PATO_0002264",
            /* organized */ "http://purl.obolibrary.org/obo/PATO_0000938",
            /* disorganized */ "http://purl.obolibrary.org/obo/PATO_0000937",
            /* congested */ "http://purl.obolibrary.org/obo/PATO_0001836",
            /* decondensed */ "http://purl.obolibrary.org/obo/PATO_0002452",
            /* hemorrhagic */ "http://purl.obolibrary.org/obo/PATO_0002105",
            /* porosity */ "http://purl.obolibrary.org/obo/PATO_0000973",
            /* decreased porosity */ "http://purl.obolibrary.org/obo/PATO_0015025",
            /* non-porous */ "http://purl.obolibrary.org/obo/PATO_0000985",
            /* porous */ "http://purl.obolibrary.org/obo/PATO_0000984",
            /* complexity */ "http://purl.obolibrary.org/obo/PATO_0001502",
            /* complex */ "http://purl.obolibrary.org/obo/PATO_0001504",
            /* simple */ "http://purl.obolibrary.org/obo/PATO_0001503",
            /* vestigial */ "http://purl.obolibrary.org/obo/PATO_0000588",
            /* demyelinated */ "http://purl.obolibrary.org/obo/PATO_0002218",
            /* unfused from */ "http://purl.obolibrary.org/obo/PATO_0000651",
            /* associated with */ "http://purl.obolibrary.org/obo/PATO_0001668",
            /* infiltrative */ "http://purl.obolibrary.org/obo/PATO_0002103",
            /* permeability */ "http://purl.obolibrary.org/obo/PATO_0000970",
            /* impermeable */ "http://purl.obolibrary.org/obo/PATO_0000983",
            /* permeable */ "http://purl.obolibrary.org/obo/PATO_0000982",
            /* decreased permeability */ "http://purl.obolibrary.org/obo/PATO_0001578",
            /* incomplete structure */ "http://purl.obolibrary.org/obo/PATO_0005013",
            /* permanent */ "http://purl.obolibrary.org/obo/PATO_0002293",
            /* interlocked with */ "http://purl.obolibrary.org/obo/PATO_0002437",
            /* damage */ "http://purl.obolibrary.org/obo/PATO_0001020",
            /* undamaged */ "http://purl.obolibrary.org/obo/PATO_0001168",
            /* damaged */ "http://purl.obolibrary.org/obo/PATO_0001167",
            /* in contact with */ "http://purl.obolibrary.org/obo/PATO_0001961",
            /* articulated with */ "http://purl.obolibrary.org/obo/PATO_0002278",
            /* broadly articulated with */ "http://purl.obolibrary.org/obo/PATO_0002280",
            /* tightly articulated with */ "http://purl.obolibrary.org/obo/PATO_0002279",
            /* stability */ "http://purl.obolibrary.org/obo/PATO_0015026",
            /* decreased stability */ "http://purl.obolibrary.org/obo/PATO_0015028",
            /* segmented */ "http://purl.obolibrary.org/obo/PATO_0002312",
            /* sutured to */ "http://purl.obolibrary.org/obo/PATO_0002469",
            /* separated from */ "http://purl.obolibrary.org/obo/PATO_0001505",
            /* diastatic */ "http://purl.obolibrary.org/obo/PATO_0001506",
            /* swollen */ "http://purl.obolibrary.org/obo/PATO_0001851",
            /* autogenous */ "http://purl.obolibrary.org/obo/PATO_0002316",
            /* structure, cavities */ "http://purl.obolibrary.org/obo/PATO_0002014",
            /* hollow */ "http://purl.obolibrary.org/obo/PATO_0002078",
            /* fluid-filled */ "http://purl.obolibrary.org/obo/PATO_0002409",
            /* unlumenized */ "http://purl.obolibrary.org/obo/PATO_0001896",
            /* inflated */ "http://purl.obolibrary.org/obo/PATO_0002376",
            /* saccular */ "http://purl.obolibrary.org/obo/PATO_0001987",
            /* cystic */ "http://purl.obolibrary.org/obo/PATO_0001673",
            /* polycystic */ "http://purl.obolibrary.org/obo/PATO_0002089",
            /* monocystic */ "http://purl.obolibrary.org/obo/PATO_0002088",
            /* perforate */ "http://purl.obolibrary.org/obo/PATO_0002112",
            /* cribriform */ "http://purl.obolibrary.org/obo/PATO_0002113",
            /* vacuolated */ "http://purl.obolibrary.org/obo/PATO_0000941",
            /* lumenized */ "http://purl.obolibrary.org/obo/PATO_0001897",
            /* uninflated */ "http://purl.obolibrary.org/obo/PATO_0002377",
            /* imperforate */ "http://purl.obolibrary.org/obo/PATO_0001821",
            /* unstructured */ "http://purl.obolibrary.org/obo/PATO_0001412",
            /* ligneous */ "http://purl.obolibrary.org/obo/PATO_0002348",
            /* separating */ "http://purl.obolibrary.org/obo/PATO_0002525",
            /* dissociated from */ "http://purl.obolibrary.org/obo/PATO_0001738",
            /* condensed */ "http://purl.obolibrary.org/obo/PATO_0001485",
            /* turgor */ "http://purl.obolibrary.org/obo/PATO_0001620",
            /* decreased turgor */ "http://purl.obolibrary.org/obo/PATO_0001621",
            /* lesioned */ "http://purl.obolibrary.org/obo/PATO_0040025",
            /* caseous */ "http://purl.obolibrary.org/obo/PATO_0002109",
            /* collapsed */ "http://purl.obolibrary.org/obo/PATO_0001478",
            /* transient */ "http://purl.obolibrary.org/obo/PATO_0002292",
            /* spatial pattern */ "http://purl.obolibrary.org/obo/PATO_0000060",
            /* spatial deviation */ "http://purl.obolibrary.org/obo/PATO_0002175",
            /* deviation towards the medial side */ "http://purl.obolibrary.org/obo/PATO_0002177",
            /* deviation towards the lateral side */ "http://purl.obolibrary.org/obo/PATO_0002176",
            /* distributed */ "http://purl.obolibrary.org/obo/PATO_0001566",
            /* symmetry */ "http://purl.obolibrary.org/obo/PATO_0000965",
            /* isometrical */ "http://purl.obolibrary.org/obo/PATO_0001573",
            /* radial symmetry */ "http://purl.obolibrary.org/obo/PATO_0001325",
            /* actinomorphic */ "http://purl.obolibrary.org/obo/PATO_0001328",
            /* bilateral symmetry */ "http://purl.obolibrary.org/obo/PATO_0001324",
            /* zygomorphic */ "http://purl.obolibrary.org/obo/PATO_0001327",
            /* symmetrical */ "http://purl.obolibrary.org/obo/PATO_0000632",
            /* asymmetrical */ "http://purl.obolibrary.org/obo/PATO_0000616",
            /* regular spatial pattern */ "http://purl.obolibrary.org/obo/PATO_0000440",
            /* decussate */ "http://purl.obolibrary.org/obo/PATO_0001953",
            /* alternate placement */ "http://purl.obolibrary.org/obo/PATO_0001932",
            /* tight */ "http://purl.obolibrary.org/obo/PATO_0001809",
            /* oriented */ "http://purl.obolibrary.org/obo/PATO_0000614",
            /* introverted */ "http://purl.obolibrary.org/obo/PATO_0001856",
            /* rotated */ "http://purl.obolibrary.org/obo/PATO_0001599",
            /* anteromedially rotated */ "http://purl.obolibrary.org/obo/PATO_0002399",
            /* laterally rotated */ "http://purl.obolibrary.org/obo/PATO_0002156",
            /* medially rotated */ "http://purl.obolibrary.org/obo/PATO_0002155",
            /* ventrally rotated */ "http://purl.obolibrary.org/obo/PATO_0001659",
            /* dorsally rotated */ "http://purl.obolibrary.org/obo/PATO_0001658",
            /* anteriorly rotated */ "http://purl.obolibrary.org/obo/PATO_0001601",
            /* posteriorly rotated */ "http://purl.obolibrary.org/obo/PATO_0001600",
            /* everted */ "http://purl.obolibrary.org/obo/PATO_0001597",
            /* anteverted */ "http://purl.obolibrary.org/obo/PATO_0001474",
            /* inverted */ "http://purl.obolibrary.org/obo/PATO_0000625",
            /* heterotaxic */ "http://purl.obolibrary.org/obo/PATO_0040000",
            /* localized */ "http://purl.obolibrary.org/obo/PATO_0000627",
            /* focally extensive */ "http://purl.obolibrary.org/obo/PATO_0002415",
            /* color pattern */ "http://purl.obolibrary.org/obo/PATO_0000019",
            /* high contrast color pattern */ "http://purl.obolibrary.org/obo/PATO_0002275",
            /* multi-colored */ "http://purl.obolibrary.org/obo/PATO_0001533",
            /* netted */ "http://purl.obolibrary.org/obo/PATO_0001947",
            /* motley */ "http://purl.obolibrary.org/obo/PATO_0001534",
            /* mottled */ "http://purl.obolibrary.org/obo/PATO_0002274",
            /* barred */ "http://purl.obolibrary.org/obo/PATO_0002276",
            /* marbled */ "http://purl.obolibrary.org/obo/PATO_0002273",
            /* blotchy */ "http://purl.obolibrary.org/obo/PATO_0000329",
            /* banded */ "http://purl.obolibrary.org/obo/PATO_0001946",
            /* dappled */ "http://purl.obolibrary.org/obo/PATO_0001535",
            /* spotted */ "http://purl.obolibrary.org/obo/PATO_0000333",
            /* mono-colored */ "http://purl.obolibrary.org/obo/PATO_0001532",
            /* variable color */ "http://purl.obolibrary.org/obo/PATO_0001515",
            /* colorless */ "http://purl.obolibrary.org/obo/PATO_0000337",
            /* discolored */ "http://purl.obolibrary.org/obo/PATO_0000331",
            /* uncrowded */ "http://purl.obolibrary.org/obo/PATO_0000633",
            /* stratification */ "http://purl.obolibrary.org/obo/PATO_0002067",
            /* unstratified */ "http://purl.obolibrary.org/obo/PATO_0002069",
            /* stratified */ "http://purl.obolibrary.org/obo/PATO_0002068",
            /* multi-localised */ "http://purl.obolibrary.org/obo/PATO_0001791",
            /* multifocal to coalescing */ "http://purl.obolibrary.org/obo/PATO_0002402",
            /* generalized */ "http://purl.obolibrary.org/obo/PATO_0002403",
            /* random pattern */ "http://purl.obolibrary.org/obo/PATO_0002401",
            /* unilateral */ "http://purl.obolibrary.org/obo/PATO_0000634",
            /* irregular spatial pattern */ "http://purl.obolibrary.org/obo/PATO_0000330",
            /* trabecular */ "http://purl.obolibrary.org/obo/PATO_0002121",
            /* loose */ "http://purl.obolibrary.org/obo/PATO_0001802",
            /* aggregated */ "http://purl.obolibrary.org/obo/PATO_0001629",
            /* patchy */ "http://purl.obolibrary.org/obo/PATO_0001608",
            /* disheveled */ "http://purl.obolibrary.org/obo/PATO_0001605",
            /* punctate */ "http://purl.obolibrary.org/obo/PATO_0001512",
            /* disoriented */ "http://purl.obolibrary.org/obo/PATO_0000613",
            /* sparse */ "http://purl.obolibrary.org/obo/PATO_0001609",
            /* whorled */ "http://purl.obolibrary.org/obo/PATO_0001951",
            /* vertical */ "http://purl.obolibrary.org/obo/PATO_0001854",
            /* unlocalised */ "http://purl.obolibrary.org/obo/PATO_0000635",
            /* undistributed */ "http://purl.obolibrary.org/obo/PATO_0001567",
            /* segmental */ "http://purl.obolibrary.org/obo/PATO_0002404",
            /* size */ "http://purl.obolibrary.org/obo/PATO_0000117",
            /* decreased size */ "http://purl.obolibrary.org/obo/PATO_0000587",
            /* decreased volume */ "http://purl.obolibrary.org/obo/PATO_0000596",
            /* hypotrophic */ "http://purl.obolibrary.org/obo/PATO_0000585",
            /* decreased height */ "http://purl.obolibrary.org/obo/PATO_0000569",
            /* decreased thickness */ "http://purl.obolibrary.org/obo/PATO_0000592",
            /* decreased diameter */ "http://purl.obolibrary.org/obo/PATO_0001715",
            /* decreased anterior-posterior diameter */ "http://purl.obolibrary.org/obo/PATO_0002042",
            /* decreased width */ "http://purl.obolibrary.org/obo/PATO_0000599",
            /* decreased width and length */ "http://purl.obolibrary.org/obo/PATO_0040030",
            /* decreased depth */ "http://purl.obolibrary.org/obo/PATO_0001472",
            /* dwarf-like */ "http://purl.obolibrary.org/obo/PATO_0000969",
            /* atrophied */ "http://purl.obolibrary.org/obo/PATO_0001623",
            /* dystrophic */ "http://purl.obolibrary.org/obo/PATO_0001780",
            /* decreased length */ "http://purl.obolibrary.org/obo/PATO_0000574",
            /* decreased perimeter */ "http://purl.obolibrary.org/obo/PATO_0001713",
            /* decreased circumference */ "http://purl.obolibrary.org/obo/PATO_0001899",
            /* decreased area */ "http://purl.obolibrary.org/obo/PATO_0002058",
            /* 2-D extent */ "http://purl.obolibrary.org/obo/PATO_0001709",
            /* area */ "http://purl.obolibrary.org/obo/PATO_0001323",
            /* 3-D extent */ "http://purl.obolibrary.org/obo/PATO_0001710",
            /* dendritic field size */ "http://purl.obolibrary.org/obo/PATO_0070064",
            /* volume */ "http://purl.obolibrary.org/obo/PATO_0000918",
            /* molar volume */ "http://purl.obolibrary.org/obo/PATO_0001680",
            /* specific volume */ "http://purl.obolibrary.org/obo/PATO_0001679",
            /* 1-D extent */ "http://purl.obolibrary.org/obo/PATO_0001708",
            /* perimeter */ "http://purl.obolibrary.org/obo/PATO_0001711",
            /* circumference */ "http://purl.obolibrary.org/obo/PATO_0001648",
            /* depth */ "http://purl.obolibrary.org/obo/PATO_0001595",
            /* diameter */ "http://purl.obolibrary.org/obo/PATO_0001334",
            /* uniform diameter */ "http://purl.obolibrary.org/obo/PATO_0005022",
            /* anterior-posterior diameter */ "http://purl.obolibrary.org/obo/PATO_0002041",
            /* width */ "http://purl.obolibrary.org/obo/PATO_0000921",
            /* thickness */ "http://purl.obolibrary.org/obo/PATO_0000915",
            /* irregular thickness */ "http://purl.obolibrary.org/obo/PATO_0001781",
            /* length */ "http://purl.obolibrary.org/obo/PATO_0000122",
            /* radius */ "http://purl.obolibrary.org/obo/PATO_0002390",
            /* height */ "http://purl.obolibrary.org/obo/PATO_0000119",
            /* tapered size */ "http://purl.obolibrary.org/obo/PATO_0005015",
            /* cell morphology */ "http://purl.obolibrary.org/obo/PATO_0010006",
            /* monostratified dendrite cell morphology */ "http://purl.obolibrary.org/obo/PATO_0070063",
            /* bistratified dendrite cell morphology */ "http://purl.obolibrary.org/obo/PATO_0070062",
            /* retinal bipolar morphology */ "http://purl.obolibrary.org/obo/PATO_0070042",
            /* bitufted dendrite cell morphology */ "http://purl.obolibrary.org/obo/PATO_0070027",
            /* double bouquet morphology */ "http://purl.obolibrary.org/obo/PATO_0070061",
            /* bitufted cell morphology */ "http://purl.obolibrary.org/obo/PATO_0070012",
            /* cortical bipolar morphology */ "http://purl.obolibrary.org/obo/PATO_0070006",
            /* multipolar neuron morphology */ "http://purl.obolibrary.org/obo/PATO_0070026",
            /* pyramidal family morphology */ "http://purl.obolibrary.org/obo/PATO_0070015",
            /* standard pyramidal morphology */ "http://purl.obolibrary.org/obo/PATO_0070017",
            /* inverted pyramidal morphology */ "http://purl.obolibrary.org/obo/PATO_0070021",
            /* stellate pyramidal morphology */ "http://purl.obolibrary.org/obo/PATO_0070020",
            /* untufted pyramidal morphology */ "http://purl.obolibrary.org/obo/PATO_0070019",
            /* tufted pyramidal morphology */ "http://purl.obolibrary.org/obo/PATO_0070018",
            /* horizontal pyramidal morphology */ "http://purl.obolibrary.org/obo/PATO_0070016",
            /* spiny stellate cell morphology */ "http://purl.obolibrary.org/obo/PATO_0070013",
            /* stellate morphology */ "http://purl.obolibrary.org/obo/PATO_0070010",
            /* Martinotti morphology */ "http://purl.obolibrary.org/obo/PATO_0070007",
            /* fan Martinotti morphology */ "http://purl.obolibrary.org/obo/PATO_0070009",
            /* T Martinotti morphology */ "http://purl.obolibrary.org/obo/PATO_0070008",
            /* basket cell morphology */ "http://purl.obolibrary.org/obo/PATO_0070002",
            /* nest basket cell morphology */ "http://purl.obolibrary.org/obo/PATO_0070005",
            /* large basket cell morphology */ "http://purl.obolibrary.org/obo/PATO_0070004",
            /* small basket cell morphology */ "http://purl.obolibrary.org/obo/PATO_0070003",
            /* neurogliaform morphology */ "http://purl.obolibrary.org/obo/PATO_0070001",
            /* unipolar neuron morphology */ "http://purl.obolibrary.org/obo/PATO_0070025",
            /* chandelier cell morphology */ "http://purl.obolibrary.org/obo/PATO_0070011",
            /* absence of anatomical entity */ "http://purl.obolibrary.org/obo/PATO_0040060",
            /* agenesis */ "http://purl.obolibrary.org/obo/PATO_0002291",
            /* glandular */ "http://purl.obolibrary.org/obo/PATO_0002117",
            /* molecular quality */ "http://purl.obolibrary.org/obo/PATO_0002182",
            /* electric polarizability */ "http://purl.obolibrary.org/obo/PATO_0002189",
            /* polarity */ "http://purl.obolibrary.org/obo/PATO_0002186",
            /* nonpolar polarity */ "http://purl.obolibrary.org/obo/PATO_0002188",
            /* polar polarity */ "http://purl.obolibrary.org/obo/PATO_0002187",
            /* nitrated */ "http://purl.obolibrary.org/obo/PATO_0002217",
            /* specificity to */ "http://purl.obolibrary.org/obo/PATO_0002030",
            /* electric charge */ "http://purl.obolibrary.org/obo/PATO_0002193",
            /* negative charge */ "http://purl.obolibrary.org/obo/PATO_0002196",
            /* positive charge */ "http://purl.obolibrary.org/obo/PATO_0002195",
            /* neutral charge */ "http://purl.obolibrary.org/obo/PATO_0002194",
            /* cyclicity */ "http://purl.obolibrary.org/obo/PATO_0002183",
            /* acyclic cyclicity */ "http://purl.obolibrary.org/obo/PATO_0002185",
            /* cyclic cyclicity */ "http://purl.obolibrary.org/obo/PATO_0002184",
            /* electronegativity */ "http://purl.obolibrary.org/obo/PATO_0002197",
            /* concentration of */ "http://purl.obolibrary.org/obo/PATO_0000033",
            /* osmolality */ "http://purl.obolibrary.org/obo/PATO_0002027",
            /* decreased osmolality */ "http://purl.obolibrary.org/obo/PATO_0002028",
            /* acidity */ "http://purl.obolibrary.org/obo/PATO_0001842",
            /* decreased acidity */ "http://purl.obolibrary.org/obo/PATO_0001843",
            /* medium acidity */ "http://purl.obolibrary.org/obo/PATO_0001428",
            /* alkaline */ "http://purl.obolibrary.org/obo/PATO_0001430",
            /* acidic */ "http://purl.obolibrary.org/obo/PATO_0001429",
            /* compound acidity */ "http://purl.obolibrary.org/obo/PATO_0001427",
            /* catalytic (activity) concentration */ "http://purl.obolibrary.org/obo/PATO_0001674",
            /* osmolarity */ "http://purl.obolibrary.org/obo/PATO_0001655",
            /* salinity */ "http://purl.obolibrary.org/obo/PATO_0085001",
            /* decreased osmolarity */ "http://purl.obolibrary.org/obo/PATO_0001656",
            /* diluted */ "http://purl.obolibrary.org/obo/PATO_0001161",
            /* concentrated */ "http://purl.obolibrary.org/obo/PATO_0001159",
            /* decreased concentration */ "http://purl.obolibrary.org/obo/PATO_0001163",
            /* capacitance */ "http://purl.obolibrary.org/obo/PATO_0002205",
            /* solubility */ "http://purl.obolibrary.org/obo/PATO_0001536",
            /* dissolved */ "http://purl.obolibrary.org/obo/PATO_0001986",
            /* insoluble in */ "http://purl.obolibrary.org/obo/PATO_0001538",
            /* soluble in */ "http://purl.obolibrary.org/obo/PATO_0001537",
            /* decreased solubility */ "http://purl.obolibrary.org/obo/PATO_0001664",
            /* oxidized */ "http://purl.obolibrary.org/obo/PATO_0002223",
            /* avidity */ "http://purl.obolibrary.org/obo/PATO_0002073",
            /* decreased avidity */ "http://purl.obolibrary.org/obo/PATO_0002075",
            /* ubiquinated */ "http://purl.obolibrary.org/obo/PATO_0002216",
            /* aromaticity */ "http://purl.obolibrary.org/obo/PATO_0002190",
            /* non-aromatic */ "http://purl.obolibrary.org/obo/PATO_0002192",
            /* aromatic */ "http://purl.obolibrary.org/obo/PATO_0002191",
            /* affinity */ "http://purl.obolibrary.org/obo/PATO_0002070",
            /* decreased affinity */ "http://purl.obolibrary.org/obo/PATO_0002072",
            /* phosphorylation */ "http://purl.obolibrary.org/obo/PATO_0002262",
            /* dephosphorylated */ "http://purl.obolibrary.org/obo/PATO_0002263",
            /* phosphorylated */ "http://purl.obolibrary.org/obo/PATO_0002220",
            /* decreased phosphorylation */ "http://purl.obolibrary.org/obo/PATO_0002222",
            /* population quality */ "http://purl.obolibrary.org/obo/PATO_0002003",
            /* quality of an ecological community */ "http://purl.obolibrary.org/obo/PCO_0000004",
            /* quality of a population */ "http://purl.obolibrary.org/obo/PCO_0000003",
            /* morbidity */ "http://purl.obolibrary.org/obo/PATO_0001415",
            /* mixed sex */ "http://purl.obolibrary.org/obo/PATO_0001338",
            /* absence of physical object */ "http://purl.obolibrary.org/obo/PATO_0040058",
            /* anatomical structure quality */ "http://purl.obolibrary.org/obo/PATO_0070044",
            /* anatomical histological quality */ "http://purl.obolibrary.org/obo/PATO_0070045",
            /* polychromatophilic */ "http://purl.obolibrary.org/obo/PATO_0070047",
            /* neutrophillic */ "http://purl.obolibrary.org/obo/PATO_0070046",
            /* amphiphilic */ "http://purl.obolibrary.org/obo/PATO_0002420",
            /* acidophilic */ "http://purl.obolibrary.org/obo/PATO_0002418",
            /* basophilic */ "http://purl.obolibrary.org/obo/PATO_0002094",
            /* cellular quality */ "http://purl.obolibrary.org/obo/PATO_0001396",
            /* cellular motility */ "http://purl.obolibrary.org/obo/PATO_0001488",
            /* decreased cellular motility */ "http://purl.obolibrary.org/obo/PATO_0002297",
            /* nucleate quality */ "http://purl.obolibrary.org/obo/PATO_0001404",
            /* nucleated */ "http://purl.obolibrary.org/obo/PATO_0002505",
            /* multinucleate */ "http://purl.obolibrary.org/obo/PATO_0001908",
            /* trinucleate */ "http://purl.obolibrary.org/obo/PATO_0001909",
            /* binucleate */ "http://purl.obolibrary.org/obo/PATO_0001406",
            /* mononucleate */ "http://purl.obolibrary.org/obo/PATO_0001407",
            /* anucleate */ "http://purl.obolibrary.org/obo/PATO_0001405",
            /* proliferative */ "http://purl.obolibrary.org/obo/PATO_0002102",
            /* self-renewal */ "http://purl.obolibrary.org/obo/PATO_0001398",
            /* nuclear/cytoplasmic ratio */ "http://purl.obolibrary.org/obo/PATO_0040071",
            /* ploidy */ "http://purl.obolibrary.org/obo/PATO_0001374",
            /* euploid */ "http://purl.obolibrary.org/obo/PATO_0001393",
            /* doubled haploid */ "http://purl.obolibrary.org/obo/PATO_0040027",
            /* diploid */ "http://purl.obolibrary.org/obo/PATO_0001394",
            /* haplodiploid */ "http://purl.obolibrary.org/obo/PATO_0001395",
            /* polyploid */ "http://purl.obolibrary.org/obo/PATO_0001377",
            /* endopolyploid */ "http://purl.obolibrary.org/obo/PATO_0001392",
            /* hexaploid */ "http://purl.obolibrary.org/obo/PATO_0001384",
            /* pentaploid */ "http://purl.obolibrary.org/obo/PATO_0001383",
            /* tetraploid */ "http://purl.obolibrary.org/obo/PATO_0001382",
            /* triploid */ "http://purl.obolibrary.org/obo/PATO_0001381",
            /* paleopolyploid */ "http://purl.obolibrary.org/obo/PATO_0001380",
            /* allopolyploidy */ "http://purl.obolibrary.org/obo/PATO_0001379",
            /* autopolyploid */ "http://purl.obolibrary.org/obo/PATO_0001378",
            /* monoploid */ "http://purl.obolibrary.org/obo/PATO_0001376",
            /* haploid */ "http://purl.obolibrary.org/obo/PATO_0001375",
            /* aneuploid */ "http://purl.obolibrary.org/obo/PATO_0001385",
            /* trisomy */ "http://purl.obolibrary.org/obo/PATO_0001389",
            /* mosaic trisomy */ "http://purl.obolibrary.org/obo/PATO_0001391",
            /* partial trisomy */ "http://purl.obolibrary.org/obo/PATO_0001390",
            /* disomy */ "http://purl.obolibrary.org/obo/PATO_0001387",
            /* uniparental disomy */ "http://purl.obolibrary.org/obo/PATO_0001388",
            /* monosomy */ "http://purl.obolibrary.org/obo/PATO_0001386",
            /* ciliatedness */ "http://purl.obolibrary.org/obo/PATO_0001408",
            /* cellular potency */ "http://purl.obolibrary.org/obo/PATO_0001397",
            /* metaplastic */ "http://purl.obolibrary.org/obo/PATO_0002101",
            /* undifferentiated */ "http://purl.obolibrary.org/obo/PATO_0002100",
            /* differentiated */ "http://purl.obolibrary.org/obo/PATO_0002099",
            /* moderately well differentiated */ "http://purl.obolibrary.org/obo/PATO_0002111",
            /* well differentiated */ "http://purl.obolibrary.org/obo/PATO_0002110",
            /* poorly differentiated */ "http://purl.obolibrary.org/obo/PATO_0002106",
            /* pluripotent */ "http://purl.obolibrary.org/obo/PATO_0001403",
            /* multipotent */ "http://purl.obolibrary.org/obo/PATO_0001402",
            /* oligopotent */ "http://purl.obolibrary.org/obo/PATO_0001401",
            /* unipotent */ "http://purl.obolibrary.org/obo/PATO_0001400",
            /* totipotent */ "http://purl.obolibrary.org/obo/PATO_0001399",
            /* high nuclear/cytoplasmic ratio */ "http://purl.obolibrary.org/obo/PATO_0040072",
            /* decreased object quality */ "http://purl.obolibrary.org/obo/PATO_0002303"
    ) );

    private final Set<String> locationTerms = new HashSet<>( Arrays.asList(
            /* prothoracic leg disc */ "http://purl.obolibrary.org/obo/FBbt_00001781",
            /* otic placode */ "http://purl.obolibrary.org/obo/UBERON_0003069",
            /* left */ "http://www.ebi.ac.uk/efo/EFO_0001658",
            /* aerial part */ "http://www.ebi.ac.uk/efo/EFO_0001901",
            /* medial */ "http://www.ebi.ac.uk/efo/EFO_0001660",
            /* anterior lateral line placode */ "http://purl.obolibrary.org/obo/UBERON_2001316",
            /* precursor */ "http://www.ebi.ac.uk/efo/EFO_0001651",
            /* neurogenic placode */ "http://purl.obolibrary.org/obo/UBERON_0009955",
            /* gonad primordium */ "http://purl.obolibrary.org/obo/UBERON_0005564",
            /* lateral */ "http://www.ebi.ac.uk/efo/EFO_0001657",
            /* apical */ "http://www.ebi.ac.uk/efo/EFO_0001653",
            /* insect visual anlage */ "http://purl.obolibrary.org/obo/UBERON_6005434",
            /* posterior lateral line placode */ "http://purl.obolibrary.org/obo/UBERON_2001156",
            /* primordium */ "http://purl.obolibrary.org/obo/UBERON_0001048",
            /* insect tracheal primordium */ "http://purl.obolibrary.org/obo/UBERON_6005037",
            /* alveolus */ "http://purl.obolibrary.org/obo/UBERON_0003215",
            /* tissue modifier */ "http://www.ebi.ac.uk/efo/EFO_0001647",
            /* insect clypeo-labral primordium */ "http://purl.obolibrary.org/obo/UBERON_6005538",
            /* ventral */ "http://www.ebi.ac.uk/efo/EFO_0001662",
            /* vagal placode 3 */ "http://purl.obolibrary.org/obo/UBERON_2001299",
            /* vagal placode 1 */ "http://purl.obolibrary.org/obo/UBERON_2001297",
            /* vagal placode 4 */ "http://purl.obolibrary.org/obo/UBERON_2001300",
            /* alveolus of lung */ "http://purl.obolibrary.org/obo/UBERON_0002299",
            /* pancreas primordium */ "http://purl.obolibrary.org/obo/UBERON_0003921",
            /* lens placode */ "http://purl.obolibrary.org/obo/UBERON_0003073",
            /* anlage */ "http://purl.obolibrary.org/obo/UBERON_0007688",
            /* liver primordium */ "http://purl.obolibrary.org/obo/UBERON_0003894",
            /* geometric modifier */ "http://www.ebi.ac.uk/efo/EFO_0001648",
            /* thymus primordium */ "http://purl.obolibrary.org/obo/UBERON_0005562",
            /* caudal */ "http://www.ebi.ac.uk/efo/EFO_0001908",
            /* insect embryonic optic lobe primordium */ "http://purl.obolibrary.org/obo/UBERON_6000186",
            /* dorsal */ "http://www.ebi.ac.uk/efo/EFO_0001656",
            /* eye primordium */ "http://purl.obolibrary.org/obo/UBERON_0003071",
            /* placode */ "http://www.ebi.ac.uk/efo/EFO_0001650",
            /* insect ventral epidermis primordium */ "http://purl.obolibrary.org/obo/UBERON_6005533",
            /* distal */ "http://www.ebi.ac.uk/efo/EFO_0001655",
            /* epibranchial placode */ "http://purl.obolibrary.org/obo/UBERON_0003078",
            /* trigeminal placode complex */ "http://purl.obolibrary.org/obo/UBERON_0003070",
            /* insect visual primordium */ "http://purl.obolibrary.org/obo/UBERON_6001059",
            /* insect dorsal epidermis primordium */ "http://purl.obolibrary.org/obo/UBERON_6005526",
            /* insect trunk mesoderm anlage */ "http://purl.obolibrary.org/obo/UBERON_6005436",
            /* primordial vasculature */ "http://purl.obolibrary.org/obo/UBERON_0014903",
            /* proximal */ "http://www.ebi.ac.uk/efo/EFO_0001661",
            /* vagal placode 2 */ "http://purl.obolibrary.org/obo/UBERON_2001298",
            /* adenohypophyseal placode */ "http://purl.obolibrary.org/obo/UBERON_0009122",
            /* olfactory placode */ "http://purl.obolibrary.org/obo/UBERON_0003050",
            /* basal */ "http://www.ebi.ac.uk/efo/EFO_0001654",
            /* right */ "http://www.ebi.ac.uk/efo/EFO_0001659" ) );


    private final Set<String> deviationFromNormalTerms = new HashSet<>( Arrays.asList(
            /* deviation (from_normal) */ "http://purl.obolibrary.org/obo/PATO_0000069",
            /* decreased quality */ "http://purl.obolibrary.org/obo/PATO_0002301",
            /* decreased intensity */ "http://purl.obolibrary.org/obo/PATO_0001783",
            /* decreased variability of temperature */ "http://purl.obolibrary.org/obo/PATO_0001307",
            /* decreased variability of color */ "http://purl.obolibrary.org/obo/PATO_0001613",
            /* decreased variability of size */ "http://purl.obolibrary.org/obo/PATO_0001957",
            /* decreased variability of rate */ "http://purl.obolibrary.org/obo/PATO_0001588",
            /* decreased object quality */ "http://purl.obolibrary.org/obo/PATO_0002303",
            /* decreased porosity */ "http://purl.obolibrary.org/obo/PATO_0015025",
            /* decreased contractility */ "http://purl.obolibrary.org/obo/PATO_0001581",
            /* decreased flexibility */ "http://purl.obolibrary.org/obo/PATO_0001777",
            /* decreased tendency */ "http://purl.obolibrary.org/obo/PATO_0002362",
            /* decreased threshold */ "http://purl.obolibrary.org/obo/PATO_0000708",
            /* decreased female fertility */ "http://purl.obolibrary.org/obo/PATO_0001830",
            /* soft */ "http://purl.obolibrary.org/obo/PATO_0000387",
            /* decreased humidity */ "http://purl.obolibrary.org/obo/PATO_0015011",
            /* decreased resistance to */ "http://purl.obolibrary.org/obo/PATO_0001651",
            /* decreased viscosity */ "http://purl.obolibrary.org/obo/PATO_0001694",
            /* decreased wetness */ "http://purl.obolibrary.org/obo/PATO_0001826",
            /* decreased acidity */ "http://purl.obolibrary.org/obo/PATO_0001843",
            /* decreased magnetism */ "http://purl.obolibrary.org/obo/PATO_0001684",
            /* decreased affinity */ "http://purl.obolibrary.org/obo/PATO_0002072",
            /* decreased contamination */ "http://purl.obolibrary.org/obo/PATO_0015033",
            /* decreased force */ "http://purl.obolibrary.org/obo/PATO_0002246",
            /* decreased weight */ "http://purl.obolibrary.org/obo/PATO_0000583",
            /* decreased fragility */ "http://purl.obolibrary.org/obo/PATO_0002056",
            /* decreased radioactivity */ "http://purl.obolibrary.org/obo/PATO_0001743",
            /* decreased male fertility */ "http://purl.obolibrary.org/obo/PATO_0001833",
            /* reduced virulence */ "http://purl.obolibrary.org/obo/PATO_0002147",
            /* decreased combustibility */ "http://purl.obolibrary.org/obo/PATO_0015023",
            /* decreased fertility */ "http://purl.obolibrary.org/obo/PATO_0001834",
            /* decreased osmolality */ "http://purl.obolibrary.org/obo/PATO_0002028",
            /* decreased efficiency */ "http://purl.obolibrary.org/obo/PATO_0001675",
            /* decreased mass density */ "http://purl.obolibrary.org/obo/PATO_0001790",
            /* decreased tonicity */ "http://purl.obolibrary.org/obo/PATO_0001619",
            /* decreased age */ "http://purl.obolibrary.org/obo/PATO_0001765",
            /* decreased susceptibility toward */ "http://purl.obolibrary.org/obo/PATO_0001670",
            /* decreased behavioural activity */ "http://purl.obolibrary.org/obo/PATO_0000761",
            /* decreased temperature */ "http://purl.obolibrary.org/obo/PATO_0001306",
            /* decreased coiling */ "http://purl.obolibrary.org/obo/PATO_0001796",
            /* decreased waxiness */ "http://purl.obolibrary.org/obo/PATO_0002383",
            /* has fewer parts of type */ "http://purl.obolibrary.org/obo/PATO_0002001",
            /* decreased life span */ "http://purl.obolibrary.org/obo/PATO_0001604",
            /* decreased mobility */ "http://purl.obolibrary.org/obo/PATO_0002283",
            /* decreased photosensitivity */ "http://purl.obolibrary.org/obo/PATO_0001697",
            /* decreased degree of illumination */ "http://purl.obolibrary.org/obo/PATO_0015015",
            /* decreased permeability */ "http://purl.obolibrary.org/obo/PATO_0001578",
            /* decreased distance */ "http://purl.obolibrary.org/obo/PATO_0000375",
            /* decreased male receptivity */ "http://purl.obolibrary.org/obo/PATO_0001726",
            /* decreased velocity */ "http://purl.obolibrary.org/obo/PATO_0002472",
            /* decreased mass */ "http://purl.obolibrary.org/obo/PATO_0001562",
            /* decreased turgor */ "http://purl.obolibrary.org/obo/PATO_0001621",
            /* decreased speed */ "http://purl.obolibrary.org/obo/PATO_0000304",
            /* hypoplastic */ "http://purl.obolibrary.org/obo/PATO_0000645",
            /* decreased phosphorylation */ "http://purl.obolibrary.org/obo/PATO_0002222",
            /* decreased curvature */ "http://purl.obolibrary.org/obo/PATO_0001593",
            /* decreased coordination */ "http://purl.obolibrary.org/obo/PATO_0001860",
            /* decreased size */ "http://purl.obolibrary.org/obo/PATO_0000587",
            /* decreased volume */ "http://purl.obolibrary.org/obo/PATO_0000596",
            /* hypotrophic */ "http://purl.obolibrary.org/obo/PATO_0000585",
            /* decreased height */ "http://purl.obolibrary.org/obo/PATO_0000569",
            /* decreased thickness */ "http://purl.obolibrary.org/obo/PATO_0000592",
            /* decreased diameter */ "http://purl.obolibrary.org/obo/PATO_0001715",
            /* decreased anterior-posterior diameter */ "http://purl.obolibrary.org/obo/PATO_0002042",
            /* decreased width */ "http://purl.obolibrary.org/obo/PATO_0000599",
            /* decreased width and length */ "http://purl.obolibrary.org/obo/PATO_0040030",
            /* decreased depth */ "http://purl.obolibrary.org/obo/PATO_0001472",
            /* dwarf-like */ "http://purl.obolibrary.org/obo/PATO_0000969",
            /* atrophied */ "http://purl.obolibrary.org/obo/PATO_0001623",
            /* dystrophic */ "http://purl.obolibrary.org/obo/PATO_0001780",
            /* decreased length */ "http://purl.obolibrary.org/obo/PATO_0000574",
            /* decreased perimeter */ "http://purl.obolibrary.org/obo/PATO_0001713",
            /* decreased circumference */ "http://purl.obolibrary.org/obo/PATO_0001899",
            /* decreased area */ "http://purl.obolibrary.org/obo/PATO_0002058",
            /* hyporesponsive to */ "http://purl.obolibrary.org/obo/PATO_0001194",
            /* decreased odor */ "http://purl.obolibrary.org/obo/PATO_0001892",
            /* decreased radiopacity */ "http://purl.obolibrary.org/obo/PATO_0002145",
            /* decreased fluid flow */ "http://purl.obolibrary.org/obo/PATO_0001838",
            /* decreased strength */ "http://purl.obolibrary.org/obo/PATO_0001779",
            /* increased fatigability */ "http://purl.obolibrary.org/obo/PATO_0001816",
            /* decreased cellular motility */ "http://purl.obolibrary.org/obo/PATO_0002297",
            /* decreased stability */ "http://purl.obolibrary.org/obo/PATO_0015028",
            /* decreased adhesivity */ "http://purl.obolibrary.org/obo/PATO_0002334",
            /* decreased concentration */ "http://purl.obolibrary.org/obo/PATO_0001163",
            /* decreased solubility */ "http://purl.obolibrary.org/obo/PATO_0001664",
            /* decreased pigmentation */ "http://purl.obolibrary.org/obo/PATO_0002251",
            /* decreased water composition */ "http://purl.obolibrary.org/obo/PATO_0001801",
            /* decreased sensitivity to irradiation */ "http://purl.obolibrary.org/obo/PATO_0001807",
            /* decreased pressure */ "http://purl.obolibrary.org/obo/PATO_0001575",
            /* decreased fluorescence */ "http://purl.obolibrary.org/obo/PATO_0001927",
            /* decreased elasticity */ "http://purl.obolibrary.org/obo/PATO_0002288",
            /* decreased fecundity */ "http://purl.obolibrary.org/obo/PATO_0001696",
            /* decreased female receptivity */ "http://purl.obolibrary.org/obo/PATO_0001724",
            /* decreased functionality */ "http://purl.obolibrary.org/obo/PATO_0001624",
            /* decreased position */ "http://purl.obolibrary.org/obo/PATO_0001476",
            /* decreased angle to */ "http://purl.obolibrary.org/obo/PATO_0002328",
            /* decreased elevation */ "http://purl.obolibrary.org/obo/PATO_0001689",
            /* decreased distribution */ "http://purl.obolibrary.org/obo/PATO_0001672",
            /* decreased avidity */ "http://purl.obolibrary.org/obo/PATO_0002075",
            /* decreased osmolarity */ "http://purl.obolibrary.org/obo/PATO_0001656",
            /* decreased tolerance to */ "http://purl.obolibrary.org/obo/PATO_0002394",
            /* decreased sensitivity toward */ "http://purl.obolibrary.org/obo/PATO_0001550",
            /* decreased amount */ "http://purl.obolibrary.org/obo/PATO_0001997",
            /* decreased process quality */ "http://purl.obolibrary.org/obo/PATO_0002302",
            /* decreased spatial extent of a process */ "http://purl.obolibrary.org/obo/PATO_0055001",
            /* decreased occurrence */ "http://purl.obolibrary.org/obo/PATO_0002052",
            /* arrested */ "http://purl.obolibrary.org/obo/PATO_0000297",
            /* decreased propagation velocity */ "http://purl.obolibrary.org/obo/PATO_0010004",
            /* decreased rate */ "http://purl.obolibrary.org/obo/PATO_0000911",
            /* decreased rate of continuous process */ "http://purl.obolibrary.org/obo/PATO_0055006",
            /* decreased rate of occurrence */ "http://purl.obolibrary.org/obo/PATO_0055004",
            /* decreased frequency */ "http://purl.obolibrary.org/obo/PATO_0000381",
            /* decreased efficacy */ "http://purl.obolibrary.org/obo/PATO_0015003",
            /* decreased sensitivity of a process */ "http://purl.obolibrary.org/obo/PATO_0001552",
            /* decreased sensitivity of a process to oxygen */ "http://purl.obolibrary.org/obo/PATO_0001554",
            /* decreased linear velocity */ "http://purl.obolibrary.org/obo/PATO_0040033",
            /* decreased duration */ "http://purl.obolibrary.org/obo/PATO_0000499",
            /* decreased duration of temperature */ "http://purl.obolibrary.org/obo/PATO_0001311",
            /* decreased variability */ "http://purl.obolibrary.org/obo/PATO_0001583",
            /* increased quality */ "http://purl.obolibrary.org/obo/PATO_0002300",
            /* increased intensity */ "http://purl.obolibrary.org/obo/PATO_0001782",
            /* increased variability */ "http://purl.obolibrary.org/obo/PATO_0001584",
            /* increased object quality */ "http://purl.obolibrary.org/obo/PATO_0002305",
            /* increased adhesivity */ "http://purl.obolibrary.org/obo/PATO_0002333",
            /* increased radiopacity */ "http://purl.obolibrary.org/obo/PATO_0002144",
            /* increased fragility */ "http://purl.obolibrary.org/obo/PATO_0002055",
            /* increased magnetism */ "http://purl.obolibrary.org/obo/PATO_0001683",
            /* increased waxiness */ "http://purl.obolibrary.org/obo/PATO_0002382",
            /* increased stability */ "http://purl.obolibrary.org/obo/PATO_0015027",
            /* high-arched */ "http://purl.obolibrary.org/obo/PATO_0002162",
            /* increased position */ "http://purl.obolibrary.org/obo/PATO_0001475",
            /* increased angle to */ "http://purl.obolibrary.org/obo/PATO_0002327",
            /* increased elevation */ "http://purl.obolibrary.org/obo/PATO_0001688",
            /* increased distribution */ "http://purl.obolibrary.org/obo/PATO_0001671",
            /* increased distance */ "http://purl.obolibrary.org/obo/PATO_0000374",
            /* increased tonicity */ "http://purl.obolibrary.org/obo/PATO_0001618",
            /* increased susceptibility toward */ "http://purl.obolibrary.org/obo/PATO_0001669",
            /* increased mass density */ "http://purl.obolibrary.org/obo/PATO_0001788",
            /* ivory */ "http://purl.obolibrary.org/obo/PATO_0002149",
            /* increased osmolarity */ "http://purl.obolibrary.org/obo/PATO_0001657",
            /* increased affinity */ "http://purl.obolibrary.org/obo/PATO_0002071",
            /* increased male receptivity */ "http://purl.obolibrary.org/obo/PATO_0001725",
            /* increased elasticity */ "http://purl.obolibrary.org/obo/PATO_0002287",
            /* increased coiling */ "http://purl.obolibrary.org/obo/PATO_0001795",
            /* increased fecundity */ "http://purl.obolibrary.org/obo/PATO_0001695",
            /* increased radioactivity */ "http://purl.obolibrary.org/obo/PATO_0001742",
            /* increased strength */ "http://purl.obolibrary.org/obo/PATO_0001778",
            /* decreased fatigability */ "http://purl.obolibrary.org/obo/PATO_0001817",
            /* increased virulence */ "http://purl.obolibrary.org/obo/PATO_0002148",
            /* increased velocity */ "http://purl.obolibrary.org/obo/PATO_0002471",
            /* increased wetness */ "http://purl.obolibrary.org/obo/PATO_0001825",
            /* hard */ "http://purl.obolibrary.org/obo/PATO_0000386",
            /* scirrhous */ "http://purl.obolibrary.org/obo/PATO_0002127",
            /* increased female fertility */ "http://purl.obolibrary.org/obo/PATO_0001831",
            /* increased age */ "http://purl.obolibrary.org/obo/PATO_0001764",
            /* increased male fertility */ "http://purl.obolibrary.org/obo/PATO_0001832",
            /* increased viscosity */ "http://purl.obolibrary.org/obo/PATO_0001693",
            /* increased pigmentation */ "http://purl.obolibrary.org/obo/PATO_0002250",
            /* increased coordination */ "http://purl.obolibrary.org/obo/PATO_0001859",
            /* increased efficiency */ "http://purl.obolibrary.org/obo/PATO_0001676",
            /* increased fluorescence */ "http://purl.obolibrary.org/obo/PATO_0001926",
            /* increased force */ "http://purl.obolibrary.org/obo/PATO_0002245",
            /* increased weight */ "http://purl.obolibrary.org/obo/PATO_0000582",
            /* increased avidity */ "http://purl.obolibrary.org/obo/PATO_0002074",
            /* increased life span */ "http://purl.obolibrary.org/obo/PATO_0001603",
            /* increased mobility */ "http://purl.obolibrary.org/obo/PATO_0002282",
            /* increased female receptivity */ "http://purl.obolibrary.org/obo/PATO_0001723",
            /* increased concentration */ "http://purl.obolibrary.org/obo/PATO_0001162",
            /* increased solubility */ "http://purl.obolibrary.org/obo/PATO_0001663",
            /* increased speed */ "http://purl.obolibrary.org/obo/PATO_0000303",
            /* hyperplastic */ "http://purl.obolibrary.org/obo/PATO_0000644",
            /* increased flexibility */ "http://purl.obolibrary.org/obo/PATO_0001776",
            /* hyperresponsive to */ "http://purl.obolibrary.org/obo/PATO_0001192",
            /* increased osmolality */ "http://purl.obolibrary.org/obo/PATO_0002029",
            /* increased mass */ "http://purl.obolibrary.org/obo/PATO_0001563",
            /* increased permeability */ "http://purl.obolibrary.org/obo/PATO_0001577",
            /* increased acidity */ "http://purl.obolibrary.org/obo/PATO_0001844",
            /* increased size */ "http://purl.obolibrary.org/obo/PATO_0000586",
            /* gigantic */ "http://purl.obolibrary.org/obo/PATO_0001940",
            /* increased area */ "http://purl.obolibrary.org/obo/PATO_0002057",
            /* increased height */ "http://purl.obolibrary.org/obo/PATO_0000570",
            /* distended */ "http://purl.obolibrary.org/obo/PATO_0001602",
            /* increased width */ "http://purl.obolibrary.org/obo/PATO_0000600",
            /* increased width and length */ "http://purl.obolibrary.org/obo/PATO_0040031",
            /* increased thickness */ "http://purl.obolibrary.org/obo/PATO_0000591",
            /* increased volume */ "http://purl.obolibrary.org/obo/PATO_0000595",
            /* ballooning */ "http://purl.obolibrary.org/obo/PATO_0002093",
            /* hypertrophic */ "http://purl.obolibrary.org/obo/PATO_0000584",
            /* increased diameter */ "http://purl.obolibrary.org/obo/PATO_0001714",
            /* increased anterior-posterior diameter */ "http://purl.obolibrary.org/obo/PATO_0002043",
            /* dilated */ "http://purl.obolibrary.org/obo/PATO_0001571",
            /* increased depth */ "http://purl.obolibrary.org/obo/PATO_0001596",
            /* increased length */ "http://purl.obolibrary.org/obo/PATO_0000573",
            /* increased perimeter */ "http://purl.obolibrary.org/obo/PATO_0001712",
            /* increased circumference */ "http://purl.obolibrary.org/obo/PATO_0001898",
            /* increased contractility */ "http://purl.obolibrary.org/obo/PATO_0001580",
            /* increased tolerance to */ "http://purl.obolibrary.org/obo/PATO_0002393",
            /* increased photosensitivity */ "http://purl.obolibrary.org/obo/PATO_0001698",
            /* increased behavioural activity */ "http://purl.obolibrary.org/obo/PATO_0000760",
            /* increased fertility */ "http://purl.obolibrary.org/obo/PATO_0001835",
            /* increased humidity */ "http://purl.obolibrary.org/obo/PATO_0015010",
            /* increased odor */ "http://purl.obolibrary.org/obo/PATO_0001893",
            /* having extra function */ "http://purl.obolibrary.org/obo/PATO_0001559",
            /* increased functionality */ "http://purl.obolibrary.org/obo/PATO_0001625",
            /* increased threshold */ "http://purl.obolibrary.org/obo/PATO_0000706",
            /* increased fluid flow */ "http://purl.obolibrary.org/obo/PATO_0001839",
            /* increased sensitivity to irradiation */ "http://purl.obolibrary.org/obo/PATO_0001808",
            /* increased cellular motility */ "http://purl.obolibrary.org/obo/PATO_0002298",
            /* has extra parts of type */ "http://purl.obolibrary.org/obo/PATO_0002002",
            /* increased porosity */ "http://purl.obolibrary.org/obo/PATO_0015024",
            /* increased sensitivity toward */ "http://purl.obolibrary.org/obo/PATO_0001549",
            /* increased tendency */ "http://purl.obolibrary.org/obo/PATO_0002361",
            /* increased degree of illumination */ "http://purl.obolibrary.org/obo/PATO_0015014",
            /* increased phosphorylation */ "http://purl.obolibrary.org/obo/PATO_0002221",
            /* increased pressure */ "http://purl.obolibrary.org/obo/PATO_0001576",
            /* increased curvature */ "http://purl.obolibrary.org/obo/PATO_0001592",
            /* increased contamination */ "http://purl.obolibrary.org/obo/PATO_0015032",
            /* increased turgor */ "http://purl.obolibrary.org/obo/PATO_0001622",
            /* increased resistance to */ "http://purl.obolibrary.org/obo/PATO_0001650",
            /* increased combustibility */ "http://purl.obolibrary.org/obo/PATO_0015022",
            /* increased temperature */ "http://purl.obolibrary.org/obo/PATO_0001305",
            /* severe increased temperature */ "http://purl.obolibrary.org/obo/PATO_0001317",
            /* moderate increased temperature */ "http://purl.obolibrary.org/obo/PATO_0001316",
            /* mild increased temperature */ "http://purl.obolibrary.org/obo/PATO_0001315",
            /* increased process quality */ "http://purl.obolibrary.org/obo/PATO_0002304",
            /* increased duration */ "http://purl.obolibrary.org/obo/PATO_0000498",
            /* increased duration of temperature */ "http://purl.obolibrary.org/obo/PATO_0001312",
            /* increased rate */ "http://purl.obolibrary.org/obo/PATO_0000912",
            /* increased rate of continuous process */ "http://purl.obolibrary.org/obo/PATO_0055005",
            /* increased rate of occurrence */ "http://purl.obolibrary.org/obo/PATO_0055003",
            /* increased frequency */ "http://purl.obolibrary.org/obo/PATO_0000380",
            /* increased linear velocity */ "http://purl.obolibrary.org/obo/PATO_0040032",
            /* temporally extended */ "http://purl.obolibrary.org/obo/PATO_0001333",
            /* increased propagation velocity */ "http://purl.obolibrary.org/obo/PATO_0010003",
            /* increased occurrence */ "http://purl.obolibrary.org/obo/PATO_0002051",
            /* increased sensitivity of a process */ "http://purl.obolibrary.org/obo/PATO_0001551",
            /* increased sensitivity of a process to oxygen */ "http://purl.obolibrary.org/obo/PATO_0001553",
            /* increased spatial extent of a process */ "http://purl.obolibrary.org/obo/PATO_0055002",
            /* increased efficacy */ "http://purl.obolibrary.org/obo/PATO_0015004",
            /* increased amount */ "http://purl.obolibrary.org/obo/PATO_0000470",
            /* multiple */ "http://purl.obolibrary.org/obo/PATO_0002118",
            /* duplicated */ "http://purl.obolibrary.org/obo/PATO_0001473",
            /* increased variability of temperature */ "http://purl.obolibrary.org/obo/PATO_0001308",
            /* increased variability of size */ "http://purl.obolibrary.org/obo/PATO_0001958",
            /* increased variability of rate */ "http://purl.obolibrary.org/obo/PATO_0001587",
            /* increased variability of color */ "http://purl.obolibrary.org/obo/PATO_0001612",
            /* abnormal */ "http://purl.obolibrary.org/obo/PATO_0000460",
            /* pathological */ "http://purl.obolibrary.org/obo/PATO_0001869"
    ) );

    private final Set<String> diseaseModifierTerms = new HashSet<>( Arrays.asList(
            /* has an isolated presentation */ "http://purl.obolibrary.org/obo/MONDO_0021128",
            /* classic or non-classic genetic disease presentation */ "http://purl.obolibrary.org/obo/MONDO_0100355",
            /* mosaic */ "http://purl.obolibrary.org/obo/MONDO_0700062",
            /* hereditary vs non-hereditary etiology */ "http://purl.obolibrary.org/obo/MONDO_0021149",
            /* syndromic or isolated */ "http://purl.obolibrary.org/obo/MONDO_0021126",
            /* idiopathic */ "http://purl.obolibrary.org/obo/MONDO_0700005",
            /* non-classic presentation */ "http://purl.obolibrary.org/obo/MONDO_0100357",
            /* non-genetic */ "http://purl.obolibrary.org/obo/MONDO_0021151",
            /* inherited */ "http://purl.obolibrary.org/obo/MONDO_0021152",
            /* idiopathic vs non-idiopathic */ "http://purl.obolibrary.org/obo/MONDO_0700004",
            /* general tumor grading characteristic */ "http://purl.obolibrary.org/obo/MONDO_0024489",
            /* tumor grade 2, general grading system */ "http://purl.obolibrary.org/obo/MONDO_0024492",
            /* non-idiopathic */ "http://purl.obolibrary.org/obo/MONDO_0700006",
            /* tumor grading characteristic */ "http://purl.obolibrary.org/obo/MONDO_0024488",
            /* tumor grade 3, general grading system */ "http://purl.obolibrary.org/obo/MONDO_0024493",
            /* tumor grade 1 or 2, general grading system */ "http://purl.obolibrary.org/obo/MONDO_0024495",
            /* tumor grade 1, general grading system */ "http://purl.obolibrary.org/obo/MONDO_0024491",
            /* complete */ "http://purl.obolibrary.org/obo/MONDO_0700063",
            /* mosaic vs complete */ "http://purl.obolibrary.org/obo/MONDO_0700061",
            /* non-iatrogenic */ "http://purl.obolibrary.org/obo/MONDO_0100427",
            /* tumor grade 3 or 4, general grading system */ "http://purl.obolibrary.org/obo/MONDO_0024497",
            /* iatrogenic */ "http://purl.obolibrary.org/obo/MONDO_0100426",
            /* opportunistic infectious */ "http://purl.obolibrary.org/obo/MONDO_0045035",
            /* infectious disease characteristic */ "http://purl.obolibrary.org/obo/MONDO_0045034",
            /* tumor grade X, general grading system */ "http://purl.obolibrary.org/obo/MONDO_0024490",
            /* restricted to specific location */ "http://purl.obolibrary.org/obo/MONDO_0045042",
            /* rare or common */ "http://purl.obolibrary.org/obo/MONDO_0021135",
            /* tumor grade 4, general grading system */ "http://purl.obolibrary.org/obo/MONDO_0024494",
            /* rare */ "http://purl.obolibrary.org/obo/MONDO_0021136",
            /* locational disease characteristic */ "http://purl.obolibrary.org/obo/MONDO_0045040",
            /* primary infectious */ "http://purl.obolibrary.org/obo/MONDO_0045036",
            /* not rare */ "http://purl.obolibrary.org/obo/MONDO_0021137",
            /* acquired */ "http://purl.obolibrary.org/obo/MONDO_0021141",
            /* classic presentation */ "http://purl.obolibrary.org/obo/MONDO_0100356",
            /* tumor grade 2 or 3, general grading system */ "http://purl.obolibrary.org/obo/MONDO_0024496",
            /* disseminated */ "http://purl.obolibrary.org/obo/MONDO_0022202",
            /* congenital */ "http://purl.obolibrary.org/obo/MONDO_0021140",
            /* has a syndromic presentation */ "http://purl.obolibrary.org/obo/MONDO_0021127",
            /* iatrogenic or non-iatrogenic */ "http://purl.obolibrary.org/obo/MONDO_0100369",
            /* congenital or acquired */ "http://purl.obolibrary.org/obo/MONDO_0021139",
            /* metastatic */ "http://purl.obolibrary.org/obo/PATO_0002098",
            /* metastatic carcinoma */ "http://purl.obolibrary.org/obo/MONDO_0024879",
            /* metastasis */ "http://www.ebi.ac.uk/efo/EFO_0009708",
            /* mild intensity */ "http://purl.obolibrary.org/obo/PATO_0000394",
            /* moderate intensity */ "http://purl.obolibrary.org/obo/PATO_0000395",
            /* increased intensity */ "http://purl.obolibrary.org/obo/PATO_0001782",
            /* severe intensity */ "http://purl.obolibrary.org/obo/PATO_0000396",
            /* decreased intensity */ "http://purl.obolibrary.org/obo/PATO_0001783",
            /* profound intensity */ "http://purl.obolibrary.org/obo/PATO_0002692",
            /* borderline intensity */ "http://purl.obolibrary.org/obo/PATO_0002628",
            /* remittent intensity */ "http://purl.obolibrary.org/obo/PATO_0001841"
    ) );

    private final Set<String> addingMaterialentityToTargetTerms = new HashSet<>( Arrays.asList(
            /* adding a material entity into a target */ "http://purl.obolibrary.org/obo/OBI_0000274",
            /* administering substance in vivo */ "http://purl.obolibrary.org/obo/OBI_0600007",
            /* passive immunization */ "http://purl.obolibrary.org/obo/OBI_0001174",
            /* epitope specific immune intervention */ "http://purl.obolibrary.org/obo/OBI_0001192",
            /* B cell epitope specific immune intervention */ "http://purl.obolibrary.org/obo/OBI_0001487",
            /* T cell epitope specific immune intervention */ "http://purl.obolibrary.org/obo/OBI_0001471",
            /* vaccination */ "http://purl.obolibrary.org/obo/VO_0000002",
            /* antibiotic administration process */ "http://purl.obolibrary.org/obo/OBI_0003057",
            /* intramuscular injection */ "http://purl.obolibrary.org/obo/OBI_0000934",
            /* administering substance in vivo to cause disease */ "http://purl.obolibrary.org/obo/OBI_0003413",
            /* in vivo challenge */ "http://purl.obolibrary.org/obo/OBI_0003456",
            /* tumor challenge */ "http://purl.obolibrary.org/obo/OBI_0003455",
            /* pathogen challenge */ "http://purl.obolibrary.org/obo/OBI_0000712",
            /* intranasal mucosal administration */ "http://purl.obolibrary.org/obo/OBI_0000983",
            /* administering substance in vivo to reduce or prevent disease */ "http://purl.obolibrary.org/obo/OBI_0003414",
            /* subcutaneous injection */ "http://purl.obolibrary.org/obo/OBI_0000954",
            /* oral administration */ "http://purl.obolibrary.org/obo/OBI_0000952",
            /* oral ingestion of pill */ "http://purl.obolibrary.org/obo/OBI_0000837",
            /* intravenous injection */ "http://purl.obolibrary.org/obo/OBI_0000994",
            /* administration in vivo with infectious agent */ "http://purl.obolibrary.org/obo/OBI_1110007",
            /* intradermal injection */ "http://purl.obolibrary.org/obo/OBI_0000942",
            /* intraperitoneal administration */ "http://purl.obolibrary.org/obo/OBI_0000429",
            /* intraperitoneal injection */ "http://purl.obolibrary.org/obo/OBI_0000281",
            /* adding substance to cell culture */ "http://purl.obolibrary.org/obo/OBI_0600000",
            /* experimental infection of cell culture */ "http://purl.obolibrary.org/obo/OBI_1110030",
            /* administration of material to specimen */ "http://purl.obolibrary.org/obo/OBI_0000995",
            /* injection */ "http://purl.obolibrary.org/obo/OBI_0000426",
            /* injection into organ section */ "http://purl.obolibrary.org/obo/OBI_0000431"
    ) );

    private class RemappingInfo {

        private FactorValue fv;
        private ExpressionExperiment ee;

        /* note that it is possible we would have more than one statement, but I'm assuming we have just one for now as
        more complicated cases (>2 characteristics) are going to be hard.
         */
        private Characteristic subject = null;
        private String predicate = null;
        private Characteristic object = null;

        private String summary = null;


        public RemappingInfo( ExpressionExperiment ee, FactorValue fv ) {
            this.ee = ee;
            this.fv = fv;
        }

        @Override
        public int hashCode() {
            return Objects.hash( fv, ee ); // this is unique to each fv
        }


        /**
         pretty-print a characteristic
         */
        private String sc( Characteristic c ) {
            return c.getValue() + " [" + ( c.getValueUri() == null ? "free text" : c.getValueUri() ) + "]";
        }

        private String getString() {
            return ee + "\n" + fv + "\n" + fv.getCharacteristics().stream().map( c ->
                    new String( " - c - " + c.getCategory() + ": " + c.getValue() + " " + ( c.getValueUri() == null ? "[free text]" : c.getValueUri() )
                            + " drug=" + isDrug( c )
                            + " phys=" + isQualityProperty( c )
                            + " deliv=" + isDelivery( c )
                            + " dismodif=" + isDiseaseModifier( c )
                            + " control=" + isControlCondition( c )
                            + " stage=" + isDevelopmentalStage( c )
                            + " time=" + isTimepoint( c )
                            + " loc=" + isLocation( c )
                            + " genet=" + isGeneticManipulation( c )
                            + "\n" ) ).collect( Collectors.joining() );
        }

        public String toString() {
            String result = getString();
            if ( summary != null ) {
                if ( predicate == null ) {
                    result = result
                            + summary + " -> " + sc( subject ) + " " + "no relation" + " " + sc( object );
                } else {
                    result = result
                            + summary + " -> " + sc( subject ) + " " + predicate + " " + sc( object );
                }
            }
            return result;
        }

        private final List<String> customOrder = Arrays.asList( "treatment", "genotype", "cell type", "cell line", "developmental stage", "strain", "organism part", "disease", "disease staging", "disease model", "phenotype", "timepoint", "age", "behavior", "delivery", "growth condition", "control", "reference subject role", "dose" );

        public String tabularize() {
            List<String> fields = new ArrayList<>();
            fields.add( ee.getShortName() );
            fields.add( ee.getId().toString() );
            fields.add( StringUtils.strip( fv.toString() ) );
            fields.add( fv.getId().toString() );

            if ( summary == null || predicate == null || subject == null ) {
                // unresolved
                List<Characteristic> listchars = new ArrayList<>( fv.getCharacteristics() );
                listchars.sort( new Comparator<Characteristic>() {
                    @Override
                    public int compare( Characteristic o1, Characteristic o2 ) {

                        try {
                            Integer o1x = customOrder.indexOf( o1.getCategory().toLowerCase() );
                            Integer o2x = customOrder.indexOf( o2.getCategory().toLowerCase() );

                            if ( o1x != null && o2x != null ) { // should always be true as NPE will be thrown if either is null (or if the category is null)
                                return o1x.compareTo( o2x );
                            }
                        } catch ( NullPointerException e ) {
                            // ignore - it was a different category
                            //System.err.println(o1.getCategory() + " " + o2.getCategory());
                        }

                        // Otherwise: put free text at the end
                        if ( o1.getValueUri() == null ) {
                            return 1;
                        } else if ( o2.getValueUri() == null ) {
                            return 1;
                        } else {
                            return o1.getValue().compareTo( o2.getValue() );
                        }
                    }
                } );
                fields.add( "" + fv.getCharacteristics().size() );
                fields.add( listchars.stream().map( c -> c.getId() + "\t" + ( c.getCategory() == null ? "" : c.getCategory() ) + "\t" + c.getValue() + "\t" + ( c.getValueUri() == null ? "" : c.getValueUri() ) ).collect( Collectors.joining( "\t" ) ) );
            } else {
                // resolved.
                fields.add( summary );
                fields.add( subject.getId().toString() );
                fields.add( object.getId().toString() );
                fields.add( subject.getCategory() );
                fields.add( sc( subject ) );
                fields.add( predicate );
                fields.add( object.getCategory() );
                fields.add( sc( object ) );
            }
            return StringUtils.join( fields, "\t" );
        }
    }
}


