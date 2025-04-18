/*
 * The Gemma project
 * 
 * Copyright (c) 2013 University of British Columbia
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

package ubic.gemma.script.example

import groovy.cli.commons.CliBuilder
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix
import ubic.gemma.core.analysis.service.ExpressionDataMatrixService
import ubic.gemma.groovy.framework.SpringSupport
import ubic.gemma.model.expression.experiment.ExpressionExperiment
import ubic.gemma.model.genome.Gene
import ubic.gemma.persistence.service.expression.experiment.ExpressionExperimentService
import ubic.gemma.persistence.service.genome.gene.GeneService
import ubic.gemma.persistence.service.genome.taxon.TaxonService

/**
 * Parse params
 */
cli = new CliBuilder()
cli.h(longOpt: 'help', 'Usage: groovy ExprMatrix -togiv')
cli.o(argName: 'file name', longOpt: 'outFile', required: true, args: 1, 'Output results to this file')
cli.t(argName: 'common name', longOpt: 'taxon', required: true, args: 1, 'Taxon of genes to fetch')
cli.g(argName: 'file name', longOpt: 'geneFile', required: false, args: 1, 'File containing list of gene official symbols to load')
cli.i(longOpt: 'filterNonSpecific', 'Filter non-specific probes')
cli.e(argName: 'file name', longOpt: 'eeFile', required: false, args: 1, 'File containing list of data sets to load')

opts = cli.parse(args)
if (!opts) return
if (opts.h) cli.usage()

geneSymbols = (opts.g) ? new File(opts.g as String).readLines() : null
filterNonSpecific = opts.i
if (filterNonSpecific) {
    System.out.println "Filtering non-specific probes"
}

taxonName = opts.t as String
eeNames = (opts.e) ? new File(opts.e as String).readLines() : null

/**
 * Gemma services
 */
sx = new SpringSupport()
taxonService = sx.getBean(TaxonService.class)
geneService = sx.getBean(GeneService.class)
eeService = sx.getBean(ExpressionExperimentService.class)
expressionDataMatrixService = sx.getBean(ExpressionDataMatrixService.class)

taxon = taxonService.findByCommonName(taxonName)

// load gene list or everything for a given taxon
if (geneSymbols != null && geneSymbols.size() > 0) {
    System.out.println "${new Date()}: Attempting to load ${geneSymbols.size()} $taxonName genes..."
    genes = geneSymbols.collect { geneService.findByOfficialSymbol(it, taxon) }
    genes = genes.findAll { it != null }
} else {
    System.out.println "${new Date()}: Loading all known $taxonName genes ..."
    //genes = geneService.loadKnownGenes(taxon);
    genes = geneService.getGenesByTaxon(taxon)
}
System.out.println "${new Date()}: Loaded ${genes.size()} $taxonName genes ..."

// load experiment list or everything for a given taxon
if (eeNames != null && eeNames.size() > 0) {
    System.out.println "${new Date()}: Attempting to load ${eeNames.size()} $taxonName experiments..."
    ees = eeNames.collect { eeService.findByShortName(it) }
    ees = ees.findAll { it != null }
} else {
    System.out.println "${new Date()}: Loading all known $taxonName experiments..."
    ees = eeService.findByTaxon(taxon)
}
System.out.println "${new Date()}: Loaded ${ees.size()} $taxonName experiments."

//TODO
/*genes = genes[1..3]
ees = ees[1..3]
println "Genes ${genes}"
println "EEs ${ees}"
*/

/**
 * Get expression ranks
 */
rssUsed = []
int count = 1
int filterCount = 0
geneRsVal = [:]
ids = new HashSet()

// main call to expressionDataMatrixService to obtain rank results
DenseDoubleMatrix<Gene, ExpressionExperiment> rankMatrix = expressionDataMatrixService.getRankMatrix(
        genes, ees, ProcessedExpressionDataVectorDao.RankMethod.mean)

f = opts.o
outFile = new File(f)
outFile.delete()
println "${new Date()}: Writing results to $f"
fOut = new BufferedWriter(new PrintWriter(outFile))

// header
count = 1
fOut << "Gene\t${ees*.shortName.collect { it }.join('\t')}\n"
//fOut << "Gene\t${ees*.id.collect{"EEID." + it}.join('\t')}\n"
for (gene in rankMatrix.getRowNames()) {
    line = "${gene.officialSymbol}"
    for (ee in ees) {
        val = rankMatrix.getByKeys(gene, ee)
        line += "\t" + ((val != null && !Double.isNaN(val)) ? sprintf("%.4f", val) : "NA")
    }
    line += "\n"
    fOut << line
    count++
}
fOut.close()

System.out.println "${new Date()}: Wrote $count genes across ${ees.size()} experiments"
System.out.println "${new Date()}: Finished"

sx.shutdown()
