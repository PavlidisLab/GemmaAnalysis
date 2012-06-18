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

package chibi.gemmaanalysis;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.basecode.io.ByteArrayConverter;
import ubic.gemma.model.expression.bioAssayData.DesignElementDataVector;
import ubic.gemma.model.expression.bioAssayData.ProcessedExpressionDataVector;
import ubic.gemma.model.expression.designElement.CompositeSequence;
import ubic.gemma.model.expression.experiment.ExpressionExperiment;

/**
 * @author xiangwan
 * @version $Id$
 */
public class ExpressionDataLoader {

    final String actualExperimentsPath = "C:/TestData/";
    final String analysisResultsPath = "C:/Results/";

    protected ExpressionExperiment experiment = null;

    protected String experimentName = null;

    protected static final Log log = LogFactory.getLog( ExpressionDataLoader.class );

    protected int uniqueItems = 0;

    protected Collection<ProcessedExpressionDataVector> designElementDataVectors = null;

    public ExpressionDataLoader( ExpressionExperiment paraExperiment ) {
        this.experiment = paraExperiment;
        if ( this.experiment != null ) {
            this.getValidDesignmentDataVector();
            this.experimentName = this.experiment.getName();
        }

    }

    private void getValidDesignmentDataVector() {
        this.designElementDataVectors = this.experiment.getProcessedExpressionDataVectors();
    }

    public void writeExpressionDataToFile( String paraFileName ) {
        BufferedWriter writer = null;
        try {
            writer = new BufferedWriter( new FileWriter( this.analysisResultsPath + paraFileName ) );
        } catch ( IOException e ) {
            throw new RuntimeException( "File for output expression data " + this.analysisResultsPath + paraFileName
                    + "could not be opened" );
        }

        try {
            writer.write( "Experiment Name: " + this.experimentName + "\n" );
            writer.write( "Accession: " + this.experiment.getAccession().getAccession() + "\n" );
            writer.write( "Name: " + this.experiment.getName() + "\n" );
            writer.write( "Description: " + this.experiment.getDescription() + "\n" );
            writer.write( "Source: " + this.experiment.getSource() + "\n" );

            for ( DesignElementDataVector dataVector : this.designElementDataVectors ) {
                CompositeSequence designElement = dataVector.getDesignElement();
                String probId = designElement.getName();
                byte[] expressionByteData = dataVector.getData();
                ByteArrayConverter byteConverter = new ByteArrayConverter();
                double[] expressionData = byteConverter.byteArrayToDoubles( expressionByteData );
                writer.write( probId + "\t" );
                for ( int i = 0; i < expressionData.length; i++ )
                    writer.write( expressionData[i] + "\t" );
                writer.write( dataVector.getQuantitationType().getName() + "\t" );
                writer.write( dataVector.getQuantitationType().getRepresentation() + "\t" );
                writer.write( dataVector.getQuantitationType().getScale().getValue() + "\t" );
                writer.write( dataVector.getQuantitationType().getType().getValue() + "\t" );
                writer.write( "\n" );
            }
            writer.close();
        } catch ( IOException e ) {
            log.error( "Error in write data into file" );
        }
    }
}