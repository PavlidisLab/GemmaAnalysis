/*
 * The Gemma project
 * 
 * Copyright (c) 2011 University of British Columbia
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

package chibi.gemmaanalysis.cli;

import java.util.Collection;
import java.util.Map;

import ubic.gemma.model.common.auditAndSecurity.AuditEvent;
import ubic.gemma.model.common.auditAndSecurity.StatusService;
import ubic.gemma.model.expression.arrayDesign.ArrayDesign;
import ubic.gemma.model.expression.arrayDesign.ArrayDesignService;
import ubic.gemma.model.expression.experiment.ExpressionExperiment;
import ubic.gemma.expression.experiment.service.ExpressionExperimentService;
import ubic.gemma.util.AbstractSpringAwareCLI;
import ubic.gemma.util.EntityUtils;

/**
 * Initialize the 'status'
 * 
 * @author paul
 * @version $Id$
 */
public class StatusPopulatorCli extends AbstractSpringAwareCLI {

    /*
     * (non-Javadoc)
     * 
     * @see ubic.gemma.util.AbstractCLI#buildOptions()
     */
    @Override
    protected void buildOptions() {
        // nothing do to
    }

    /*
     * (non-Javadoc)
     * 
     * @see ubic.gemma.util.AbstractCLI#doWork(java.lang.String[])
     */
    @Override
    protected Exception doWork( String[] args ) {

        processCommandLine( "Populate status", args );

        ExpressionExperimentService eeService = ( ExpressionExperimentService ) this
                .getBean( "expressionExperimentService" );
        ArrayDesignService adService = ( ArrayDesignService ) this.getBean( "arrayDesignService" );
        StatusService statusService = ( StatusService ) this.getBean( "statusService" );

        /*
         * Array designs trouble
         */
        Collection<ArrayDesign> ads = adService.loadAll();
        Map<Long, AuditEvent> adTroubles = adService.getLastTroubleEvent( EntityUtils.getIds( ads ) );
        Map<Long, ArrayDesign> adidmap = EntityUtils.getIdMap( ads );

        for ( Long id : EntityUtils.getIds( ads ) ) {
            ArrayDesign ad = adidmap.get( id );
            if ( ad == null ) continue;

            if ( adTroubles.containsKey( id ) && !ad.getStatus().getTroubled() ) {
                ad.getStatus().setTroubled( true );
                statusService.update( ad.getStatus() );
            } else if ( !adTroubles.containsKey( id ) && ad.getStatus().getTroubled() ) {
                // flip it.
                ad.getStatus().setTroubled( false );
                statusService.update( ad.getStatus() );
            }
        }

        /*
         * AD validations
         */
        Map<Long, AuditEvent> adValidations = adService.getLastValidationEvent( EntityUtils.getIds( ads ) );
        for ( Long id : EntityUtils.getIds( ads ) ) {
            ArrayDesign ad = adidmap.get( id );
            if ( ad == null ) continue;

            if ( adValidations.containsKey( id ) && !ad.getStatus().getValidated() ) {
                ad.getStatus().setValidated( true );
                statusService.update( ad.getStatus() );
            } else if ( !adTroubles.containsKey( id ) && ad.getStatus().getTroubled() ) {
                ad.getStatus().setValidated( false );
                statusService.update( ad.getStatus() );
            }
        }

        /*
         * Experiments trouble
         */
        Collection<ExpressionExperiment> ees = eeService.loadAll();
        Map<Long, ExpressionExperiment> eeidmap = EntityUtils.getIdMap( ees );
        Map<Long, AuditEvent> eeTroubles = eeService.getLastTroubleEvent( EntityUtils.getIds( ees ) );

        for ( Long id : EntityUtils.getIds( ees ) ) {
            ExpressionExperiment ee = eeidmap.get( id );

            if ( ee == null ) continue;

            if ( eeTroubles.containsKey( id ) && !ee.getStatus().getTroubled() ) {
                ee.getStatus().setTroubled( true );
                statusService.update( ee.getStatus() );

            } else if ( !eeTroubles.containsKey( id ) && ee.getStatus().getTroubled() ) {
                ee.getStatus().setTroubled( false );
                statusService.update( ee.getStatus() );

            }
        }

        /*
         * Experiments validated
         */
        Map<Long, AuditEvent> eeValidations = eeService.getLastValidationEvent( EntityUtils.getIds( ees ) );
        for ( Long id : EntityUtils.getIds( ees ) ) {
            ExpressionExperiment ee = eeidmap.get( id );

            if ( ee == null ) continue;

            if ( eeValidations.containsKey( id ) && !ee.getStatus().getValidated() ) {
                ee.getStatus().setValidated( true );
                statusService.update( ee.getStatus() );
            } else if ( !eeTroubles.containsKey( id ) && ee.getStatus().getValidated() ) {
                ee.getStatus().setValidated( false );
                statusService.update( ee.getStatus() );
            }
        }

        return null;
    }

    public static void main( String[] args ) {
        StatusPopulatorCli c = new StatusPopulatorCli();
        c.doWork( args );
    }

}
