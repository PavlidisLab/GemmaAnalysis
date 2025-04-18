#!/usr/bin/groovy
package ubic.gemma.script.example

import ubic.gemma.core.analysis.expression.coexpression.CoexpressionValueObjectExt
import ubic.gemma.core.analysis.expression.coexpression.GeneCoexpressionSearchService
import ubic.gemma.groovy.framework.SpringSupport
import ubic.gemma.model.genome.Gene
import ubic.gemma.persistence.service.expression.experiment.ExpressionExperimentService
import ubic.gemma.persistence.service.genome.gene.GeneService
import ubic.gemma.persistence.service.genome.taxon.TaxonService

// Tool for downloading coexpression links. See Bug 4158.
def cli = new CliBuilder(usage: 'groovy CoexSearch [opts] -tgoesr')
        .arg('o', argName: 'file name', longOpt: 'outFile', required: true, args: 1, 'Output results to this file')
        .arg('t', argName: 'common name', longOpt: 'taxon', required: true, args: 1, 'Taxon of genes to fetch')
        .arg('g', argName: 'file name', longOpt: 'geneFile', required: false, args: 1, 'File containing list of gene official symbols to load')
        .arg('e', argName: 'file name', longOpt: 'eeFile', required: false, args: 1, 'File containing list of data sets to load')
        .arg('s', argName: 'stringency', longOpt: 'stringency', required: false, args: 1, 'Minimum support per coex link')
        .arg('r', argName: 'maxResultsPerGene', longOpt: 'maxResultsPerGene', required: false, args: 1, 'Maximum number of genes per query')
        .arg('q', argName: 'queryGenesOnly', longOpt: 'queryGenesOnly', required: false, args: 1, 'Output query genes only?')

// ------------ Parse Args
println "Parsing args ..."
opts = cli.parse(args)
if (!opts) {
    println "Opts is ${opts}"
    return
}
if (opts.hasOption("h")) cli.usage()

outFile = (opts.hasOption("o")) ? new File(opts.getOptionValue("o")) : null
taxonName = opts.getOptionValue("t")
geneSymbols = (opts.hasOption("g")) ? new File(opts.getOptionValue("g")).readLines() : null
eeNames = (opts.hasOption("e")) ? new File(opts.getOptionValue("e")).readLines() : null
stringency = (opts.hasOption("s")) ? Integer.parseInt(opts.getOptionValue("s")) : 4
maxResultsPerGene = (opts.hasOption("r")) ? Integer.parseInt(opts.getOptionValue("r")) : 500
queryGenesOnly = (opts.hasOption("q")) ? Integer.parseInt(opts.getOptionValue("q")) : false


// ------------ Register beans

sx = new SpringSupport()
geneCoexpressionSearchService = sx.getBean(GeneCoexpressionSearchService.class)
taxonService = sx.getBean(TaxonService.class)
geneService = sx.getBean(GeneService.class)
eeService = sx.getBean(ExpressionExperimentService.class)

// ------------ Get IDs 
println "Loading IDs ..."
taxon = taxonService.findByCommonName(taxonName)

if (geneSymbols != null && geneSymbols.size() > 0) {
    System.out.println "${new Date()}: Attempting to load ${geneSymbols.size()} $taxonName genes..."
    genes = geneSymbols.collect { geneService.findByOfficialSymbol(it, taxon) }
    Collection<Gene> genes = genes.findAll { it != null }
    geneIds = genes.collect { it.getId() }
} else {
    System.out.println "${new Date()}: Loading all known $taxonName genes..."
    //genes = geneService.loadKnownGenes(taxon);
    Collection<Gene> genes = geneService.loadAll(taxon)
    geneIds = genes.collect { it.getId() }
}
System.out.println "${new Date()}: Loaded ${genes.size()} $taxonName genes."

ees = (eeNames == null) ? eeService.loadAll() : eeNames.collect { eeService.findByShortName(it) }
eeIds = ees.collect { it.getId() }

// ------------ Do work
println "Started coex search for ${geneIds.size()} genes across ${eeIds.size()} experiments on ${new Date()}"
println "eeIds.size = ${eeIds.size()}; geneIds.size = ${geneIds.size()}; stringency = ${stringency}; maxResultsPerGene = ${maxResultsPerGene}; queryGenesOnly = ${queryGenesOnly}"
results = geneCoexpressionSearchService.coexpressionSearch(eeIds, geneIds, stringency, maxResultsPerGene, queryGenesOnly).getResults()
if (results == null) {
    println "Error. No results found."
    return
}
println "Coex search completed on ${new Date()} and found ${results.size()} links"

// ------------ Write results
println "${new Date()}: Writing ${results.size()} results to ${outFile.getAbsolutePath()}"
fOut = new BufferedWriter(new PrintWriter(outFile))
header = ["Query Gene", "Query Gene NCBI Id", "Coexpressed Gene", "Coexpressed Gene NCBI Id", "Specificity",
          "Positive Support", "Negative Support", "Datasets tested"]
delim = "\t"
fOut << header.join(delim) + "\n"

print "Starting..."
lineCount = 0
for (CoexpressionValueObjectExt coexresult : results) {
    line = [coexresult.queryGene.officialSymbol, coexresult.queryGene.ncbiId,
            coexresult.foundGene.officialSymbol, coexresult.foundGene.ncbiId,
            coexresult.foundGeneNodeDegree > coexresult.queryGeneNodeDegree ? coexresult.foundGeneNodeDegree : coexresult.queryGeneNodeDegree, coexresult.posSupp, coexresult.negSupp, coexresult.numTestedIn]
    fOut << line.join(delim) + "\n"
    //if( lineCount % 100 == 0) {
    //    print "${lineCount}...";
    //}
    lineCount++
}
println "Done"

fOut.close()

println "Finished on ${new Date()}. Written ${lineCount} to ${outFile.getAbsolutePath()}"

sx.shutdown()
