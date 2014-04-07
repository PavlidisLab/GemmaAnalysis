/*
 * The GemmaAnalysis project
 * 
 * Copyright (c) 2014 University of British Columbia
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

import java.util.Collection;
import java.util.List;

import ubic.gemma.apps.ExpressionExperimentManipulatingCLI;
import ubic.gemma.expression.experiment.service.ExperimentalDesignService;
import ubic.gemma.model.common.auditAndSecurity.AuditEvent;
import ubic.gemma.model.common.auditAndSecurity.AuditTrailService;
import ubic.gemma.model.common.auditAndSecurity.Status;
import ubic.gemma.model.common.auditAndSecurity.eventType.TroubleStatusFlagEvent;
import ubic.gemma.model.common.description.BibliographicReference;
import ubic.gemma.model.expression.arrayDesign.ArrayDesign;
import ubic.gemma.model.expression.arrayDesign.ArrayDesignService;
import ubic.gemma.model.expression.experiment.BioAssaySet;
import ubic.gemma.model.expression.experiment.ExperimentalDesign;
import ubic.gemma.model.expression.experiment.ExpressionExperiment;
import ubic.gemma.model.expression.experiment.ExpressionExperimentSubSet;
import ubic.gemma.model.expression.experiment.ExpressionExperimentValueObject;
import ubic.gemma.model.genome.Taxon;

/**
 * TODO Document Me
 * 
 * @author paul
 * @version $Id$
 */
public class ExperimentMetaDataExtractorCli extends ExpressionExperimentManipulatingCLI {

    private ArrayDesignService adService;

    /*
     * (non-Javadoc)
     * 
     * @see ubic.gemma.util.AbstractCLI#doWork(java.lang.String[])
     */
    @Override
    protected Exception doWork( String[] args ) {
        super.processCommandLine( "experiment metadata extract", args );
        adService = getBean( ArrayDesignService.class );
        auditTrailService = getBean( AuditTrailService.class );
        experimentalDesignService = getBean( ExperimentalDesignService.class );
        for ( BioAssaySet bas : super.expressionExperiments ) {

            process( bas );

        }

        return null;
    }

    /**
     * @param bas
     */
    private void process( BioAssaySet bas ) {
        // File 1: For all experiments on Gemma
        // Experiment GSE
        // Taxon
        // Date of upload
        // Date of first curation
        //
        // Platform (GPL#)
        // Two or one channelled
        // Total number profiles Profiles:
        // Number profiles remaining after filtering:
        //
        // Number samples
        // Number of conditions
        // Number of factors
        //
        // Number of replicates per condition
        // suspected number of outliers (via autodetection)
        // Manually removed outlier number
        //
        // Batch effect p-value (NA if none)
        // Marked as trouble anywhere in history?
        //
        // Category (Normal, Exon, RNASeq)
        // Publication status?
        // Publication year
        // Publication journal
        //

        /*
         * Skip subsets.
         */
        if ( bas instanceof ExpressionExperimentSubSet ) return;

        ExpressionExperiment ee = ( ExpressionExperiment ) bas;

        log.info( "Processing: " + ee );

        ee = eeService.thawLite( ee );

        Collection<ArrayDesign> arrayDesignsUsed = eeService.getArrayDesignsUsed( ee );

        for ( ArrayDesign ad : arrayDesignsUsed ) {
            ad = adService.thawLite( ad );
        }

        Taxon t = eeService.getTaxon( ee );

        List<AuditEvent> events = auditEventService.getEvents( ee );
        for ( AuditEvent auditEvent : events ) {
            if ( auditEvent.getEventType() != null ) {
                // this is the first "curation" event? Not clear, because many steps are automated. What should we
                // count?

                // if not the first, then we can see if was 'trouble'.
                if ( auditEvent.getEventType() instanceof TroubleStatusFlagEvent ) {
                    // has trouble at some point...
                }

            }
        }

        Status status = ee.getStatus();

        if ( status.getTroubled() ) {

        }

        ExpressionExperimentValueObject vo = eeService.loadValueObject( ee.getId() );

        ExperimentalDesign experimentalDesign = ee.getExperimentalDesign();

        BibliographicReference primaryPublication = ee.getPrimaryPublication();

    }

    AuditTrailService auditTrailService;
    ExperimentalDesignService experimentalDesignService;
}
