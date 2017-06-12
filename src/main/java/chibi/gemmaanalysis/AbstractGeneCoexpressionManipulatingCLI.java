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

import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;

import ubic.gemma.core.apps.ExpressionExperimentManipulatingCLI;
import ubic.gemma.core.apps.GemmaCLI.CommandGroup;
import ubic.gemma.model.genome.Gene;

/**
 * Class for CLIs that manipulate a list of genes
 *
 * @author Raymond
 * @version $Id: AbstractGeneCoexpressionManipulatingCLI.java,v 1.3 2015/11/30 23:21:39 paul Exp $
 */
public abstract class AbstractGeneCoexpressionManipulatingCLI extends ExpressionExperimentManipulatingCLI {

    private String[] queryGeneSymbols;
    private String queryGeneFile;

    private String[] targetGeneSymbols;
    private String targetGeneFile;

    @Override
    public CommandGroup getCommandGroup() {
        return CommandGroup.ANALYSIS;
    }

    public Collection<Gene> getQueryGenes() throws IOException {
        Collection<Gene> genes = new HashSet<>();
        if ( queryGeneFile != null ) genes.addAll( readGeneListFile( queryGeneFile, taxon ) );
        if ( queryGeneSymbols != null ) {
            for ( int i = 0; i < queryGeneSymbols.length; i++ ) {
                genes.add( findGeneByOfficialSymbol( queryGeneSymbols[i], taxon ) );
            }
        }

        return genes;
    }

    public Collection<Gene> getTargetGenes() throws IOException {
        Collection<Gene> genes = new HashSet<>();
        if ( targetGeneFile != null ) genes.addAll( readGeneListFile( targetGeneFile, taxon ) );
        if ( targetGeneSymbols != null ) {
            for ( int i = 0; i < targetGeneSymbols.length; i++ ) {
                genes.add( findGeneByOfficialSymbol( targetGeneSymbols[i], taxon ) );
            }
        }
        return genes;
    }

    @SuppressWarnings("static-access")
    @Override
    protected void buildOptions() {
        super.buildOptions();
        OptionBuilder.hasArg();
        OptionBuilder
                .withDescription( "Query file containing list of gene official symbols" );
        OptionBuilder.withArgName( "File name" );
        OptionBuilder
                .withLongOpt( "queryGeneFile" );
        Option queryGeneFileOption = OptionBuilder.create();
        addOption( queryGeneFileOption );
        OptionBuilder.hasArgs();
        OptionBuilder.withArgName( "Gene symbol(s)" );
        OptionBuilder
                .withDescription( "The query gene(s)" );
        OptionBuilder.withLongOpt( "queryGene" );
        Option queryGeneOption = OptionBuilder.create();
        addOption( queryGeneOption );

        OptionBuilder.hasArg();
        OptionBuilder.withArgName( "File name" );
        OptionBuilder
                .withDescription( "File containing list of target gene official symbols" );
        OptionBuilder
                .withLongOpt( "targetGeneFile" );
        Option targetFileOption = OptionBuilder.create();
        addOption( targetFileOption );
        OptionBuilder.hasArgs();
        OptionBuilder.withArgName( "Gene symbol(s)" );
        OptionBuilder
                .withDescription( "The target gene(s)" );
        OptionBuilder.withLongOpt( "targetGene" );
        Option targetGeneOption = OptionBuilder.create();
        addOption( targetGeneOption );

    }

    protected Map<String, String> getGeneIdPair2NameMap( Collection<Gene> queryGenes, Collection<Gene> targetGenes ) {
        Map<String, String> map = new HashMap<>();
        for ( Gene qGene : queryGenes ) {
            String qName = ( qGene.getOfficialSymbol() != null ) ? qGene.getOfficialSymbol() : qGene.getId().toString();
            for ( Gene tGene : targetGenes ) {
                String tName = ( tGene.getOfficialSymbol() != null ) ? tGene.getOfficialSymbol() : tGene.getId()
                        .toString();
                map.put( qGene.getId() + ":" + tGene.getId(), qName + ":" + tName );
            }
        }
        return map;
    }

    @Override
    protected void processOptions() {
        super.processOptions();
        if ( hasOption( "queryGeneFile" ) ) queryGeneFile = getOptionValue( "queryGeneFile" );
        if ( hasOption( "queryGene" ) ) queryGeneSymbols = getOptionValues( "queryGene" );

        if ( hasOption( "targetGeneFile" ) ) targetGeneFile = getOptionValue( "targetGeneFile" );
        if ( hasOption( "targetGene" ) ) targetGeneSymbols = getOptionValues( "targetGene" );
    }

}
