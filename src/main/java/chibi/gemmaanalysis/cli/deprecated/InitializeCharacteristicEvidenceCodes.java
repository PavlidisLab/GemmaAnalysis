/*
 * The Gemma project
 * 
 * Copyright (c) 2006 University of British Columbia
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

import java.util.Collection;
import java.util.Map;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;

import ubic.gemma.model.association.GOEvidenceCode;
import ubic.gemma.model.association.Gene2GOAssociation;
import ubic.gemma.model.association.Gene2GOAssociationImpl;
import ubic.gemma.model.common.description.Characteristic;
import ubic.gemma.model.common.description.CharacteristicService;
import ubic.gemma.model.expression.biomaterial.BioMaterialImpl;
import ubic.gemma.model.expression.experiment.ExperimentalFactor;
import ubic.gemma.model.expression.experiment.ExperimentalFactorService;
import ubic.gemma.model.expression.experiment.ExpressionExperimentImpl;
import ubic.gemma.model.expression.experiment.FactorValueImpl;
import ubic.gemma.util.AbstractSpringAwareCLI;

/**
 * @author luke
 * @version $Id$
 */
public class InitializeCharacteristicEvidenceCodes extends AbstractSpringAwareCLI {

    private boolean justTesting = false;

    @SuppressWarnings("static-access")
    @Override
    protected void buildOptions() {
        Option justTestingOption = OptionBuilder.withLongOpt( "tryout" ).withDescription(
                "Set to run without any database modifications" ).create( 'f' );
        addOption( justTestingOption );
    }

    public static void main( String[] args ) {
        InitializeCharacteristicEvidenceCodes p = new InitializeCharacteristicEvidenceCodes();
        try {
            Exception ex = p.doWork( args );
            if ( ex != null ) {
                ex.printStackTrace();
            }
        } catch ( Exception e ) {
            throw new RuntimeException( e );
        }
    }

    CharacteristicService characteristicService;
    ExperimentalFactorService experimentalFactorService;

    @Override
    protected void processOptions() {
        super.processOptions();
        if ( this.hasOption( 'f' ) ) {
            this.justTesting = true;
            log.info( "TEST MODE: NO DATABASE UPDATES WILL BE PERFORMED" );
        }

        characteristicService = ( CharacteristicService ) this.getBean( "characteristicService" );
        experimentalFactorService = ( ExperimentalFactorService ) this.getBean( "experimentalFactorService" );
    }

    @Override
    protected Exception doWork( String[] args ) {

        Exception err = processCommandLine( "Initialize Characteristic Evidence Codes", args );
        if ( err != null ) return err;

        Map<Characteristic, Object> eeChars = characteristicService.findByParentClass( ExpressionExperimentImpl.class );
        initializeCharacteristicEvidenceCode( eeChars.keySet(), GOEvidenceCode.IC );

        Map<Characteristic, Object> bmChars = characteristicService.findByParentClass( BioMaterialImpl.class );
        initializeCharacteristicEvidenceCode( bmChars.keySet(), GOEvidenceCode.IEA );

        Map<Characteristic, Object> fvChars = characteristicService.findByParentClass( FactorValueImpl.class );
        initializeFactorValueCharacteristics( fvChars.keySet() );

        Map<Characteristic, Object> goChars = characteristicService.findByParentClass( Gene2GOAssociationImpl.class );
        intializeGene2GOAssociationCharacteristics( goChars );

        initializeExperimentalFactorCategoryCharacteristics();

        return null;

    }

    private void initializeCharacteristicEvidenceCode( Collection<Characteristic> chars, GOEvidenceCode ec ) {
        for ( Characteristic c : chars ) {
            c.setEvidenceCode( ec );
            if ( !justTesting ) characteristicService.update( c );
        }
    }

    private void initializeFactorValueCharacteristics( Collection<Characteristic> chars ) {
        for ( Characteristic c : chars ) {
            if ( c.getDescription() != null
                    && ( c.getDescription().startsWith( "GEO" ) || c.getDescription().startsWith( "Converted from GEO" ) ) ) {
                c.setEvidenceCode( GOEvidenceCode.IEA );

                if ( !justTesting ) characteristicService.update( c );
            } else {
                c.setEvidenceCode( GOEvidenceCode.IC );
                if ( !justTesting ) characteristicService.update( c );
            }
        }
    }

    private void intializeGene2GOAssociationCharacteristics( Map<Characteristic, Object> charToParent ) {
        for ( Characteristic c : charToParent.keySet() ) {
            Gene2GOAssociation go = ( Gene2GOAssociation ) charToParent.get( c );
            c.setEvidenceCode( go.getEvidenceCode() );
            if ( !justTesting ) characteristicService.update( c );
        }
    }

    private void initializeExperimentalFactorCategoryCharacteristics() {
        Collection<ExperimentalFactor> experimentalFactors = experimentalFactorService.loadAll();
        for ( ExperimentalFactor ef : experimentalFactors ) {
            Characteristic c = ef.getCategory();
            if ( c == null ) continue;
            if ( c.getDescription() != null
                    && ( c.getDescription().startsWith( "GEO" ) || c.getDescription().startsWith( "Converted from GEO" ) ) ) {
                c.setEvidenceCode( GOEvidenceCode.IEA );
                if ( !justTesting ) characteristicService.update( c );
            } else {
                c.setEvidenceCode( GOEvidenceCode.IC );
                if ( !justTesting ) characteristicService.update( c );
            }
        }
    }

}
