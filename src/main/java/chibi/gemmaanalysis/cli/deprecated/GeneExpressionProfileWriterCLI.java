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
package chibi.gemmaanalysis.cli.deprecated;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.HashSet;
import java.util.Map;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;

import chibi.gemmaanalysis.AbstractGeneCoexpressionManipulatingCLI;
import chibi.gemmaanalysis.CoexpressionAnalysisService;
import ubic.gemma.core.analysis.preprocess.filter.FilterConfig;
import ubic.gemma.core.apps.GemmaCLI.CommandGroup;
import ubic.gemma.core.datastructure.matrix.ExpressionDataDoubleMatrix;
import ubic.gemma.model.expression.arrayDesign.ArrayDesign;
import ubic.gemma.model.expression.designElement.CompositeSequence;
import ubic.gemma.model.expression.experiment.BioAssaySet;
import ubic.gemma.model.expression.experiment.ExpressionExperiment;
import ubic.gemma.model.genome.Gene;
import ubic.gemma.persistence.service.expression.arrayDesign.ArrayDesignService;

/**
 * @author raymond
 * @version $Id: GeneExpressionProfileWriterCLI.java,v 1.8 2015/11/30 23:50:56 paul Exp $
 * @deprecated this seems redundant with GeneExpressionWriterCLI
 */
@Deprecated
public class GeneExpressionProfileWriterCLI extends AbstractGeneCoexpressionManipulatingCLI {

    public static void main( String[] args ) {
        GeneExpressionProfileWriterCLI cli = new GeneExpressionProfileWriterCLI();
        Exception e = cli.doWork( args );
        if ( e != null ) log.error( e.getMessage() );
    }

    private ArrayDesignService adService;

    private CoexpressionAnalysisService coexpAnalysisService;

    private String outFilePrefix;

    @Override
    public CommandGroup getCommandGroup() {
        return CommandGroup.ANALYSIS;
    }

    /*
     * (non-Javadoc)
     *
     * @see ubic.gemma.util.AbstractCLI#getCommandName()
     */
    @Override
    public String getCommandName() {
        return "null";
    }

    @Override
    protected void buildOptions() {
        super.buildOptions();

        Option outFilePrefixOption = OptionBuilder.create( 'o' );
        addOption( outFilePrefixOption );
    }

    @Override
    protected Exception doWork( String[] args ) {
        processCommandLine( args );

        Collection<Gene> genes;
        try {
            genes = getQueryGenes();
        } catch ( IOException e ) {
            return e;
        }

        FilterConfig filterConfig = new FilterConfig();
        for ( BioAssaySet bas : this.expressionExperiments ) {
            ExpressionExperiment ee = ( ExpressionExperiment ) bas;
            Collection<ArrayDesign> ads = eeService.getArrayDesignsUsed( ee );
            Collection<CompositeSequence> css = new HashSet<>();
            for ( ArrayDesign ad : ads ) {
                css.addAll( adService.getCompositeSequences( ad ) );
            }
            Map<Gene, Collection<CompositeSequence>> gene2css = coexpAnalysisService.getGene2CsMap( css );
            ExpressionDataDoubleMatrix dataMatrix = coexpAnalysisService.getExpressionDataMatrix( ee, filterConfig );
            String fileName = outFilePrefix + "." + ee.getShortName() + ".txt";
            try (PrintWriter out = new PrintWriter( new FileWriter( fileName ) );) {
                for ( Gene gene : genes ) {
                    Collection<CompositeSequence> c = gene2css.get( gene );
                    for ( CompositeSequence cs : c ) {
                        Double[] row = dataMatrix.getRow( cs );
                        if ( row == null ) {
                            log.error( "Cannot get data from data matrix for " + gene.getOfficialSymbol() + " ("
                                    + cs.getName() + ")" );
                            continue;
                        }
                        StringBuffer buf = new StringBuffer();
                        buf.append( gene.getOfficialSymbol() + "\t" + cs.getName() + "\t" );
                        for ( Double d : row ) {
                            if ( d == null )
                                buf.append( "NaN" );
                            else
                                buf.append( d );
                            buf.append( "\t" );
                        }
                        buf.deleteCharAt( buf.length() - 1 );
                        out.println( buf.toString() );
                    }
                }

            } catch ( IOException e ) {
                return e;
            }
        }

        return null;
    }

    @Override
    protected void processOptions() {
        super.processOptions();
        outFilePrefix = getOptionValue( 'o' );

        adService = getBean( ArrayDesignService.class );
        coexpAnalysisService = getBean( CoexpressionAnalysisService.class );
    }

}