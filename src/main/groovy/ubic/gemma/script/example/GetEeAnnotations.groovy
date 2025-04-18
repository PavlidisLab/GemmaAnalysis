package ubic.gemma.script.example

import ubic.gemma.persistence.service.expression.experiment.ExpressionExperimentService
import ubic.gemma.persistence.service.genome.taxon.TaxonService

cli = new CliBuilder(usage: 'groovy GetEeAnnotations.groovy -e EXPERIMENT -t TAXON -o FILE')
        .arg('o', argName: 'file name', longOpt: 'outFile prefix', required: true, args: 1, '')
        .arg('t', argName: 'taxon', longOpt: 'taxon', required: true, args: 1, '')
        .arg('e', argName: 'file name', longOpt: 'eeFile', required: false, args: 1, 'File containing list of data sets to load')

opts = cli.parse(args)
if (!opts) return
if (opts.hasOption("h")) cli.usage()
if (!opts.hasOption("t")) cli.usage()

sx = new SpringSupport()
taxonService = sx.getBean(TaxonService.class)
eeService = sx.getBean(ExpressionExperimentService.class)

f = opts.getOptionValue("o") as String
outFile = new File(f)
outFile.delete()

/**
 * Load taxon
 */
taxonName = opts.getOptionValue("t")
taxon = taxonService.findByCommonName(taxonName)

/**
 * Load list of ees from a file or from Gemma (for a given taxon)
 */
//ees = eeService.loadAll()
eeNames = (opts.hasOption("e")) ? new File(opts.getOptionValue("e")).readLines() : null
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

i = 0

f = opts.getOptionValue("o") + "-eeMetaData.txt"
println "${new Date()}: Starting ...."
out = new BufferedWriter(new FileWriter(f))
//out << "ARSID\tName\tDescription\tValues\teeShortName\teeName\tTaxon\tNumSamples\tArrayDesignNames\tArrayDesignSizes\n"
out << "eeShortName\tNumAnnot\tAnnotations\n"
count = 1
for (ee in ees) {
    annots = eeService.getAnnotations(ee.id)
    annotStr = annots.collect { "\"${it.termName}\" : ${it.termUri}" }.join(', ')

    line = "${ee.shortName}\t${annots.size()}\t${annotStr}\n"
    out << line

    if (count++ % 50 == 0) {
        System.out.println "${new Date()}: Processed $count experiments"
    }
}

out.flush()
out.close()

System.out.println "${new Date()}: Done. Output written to \'${(new File(f)).getAbsolutePath()}\'"

sx.shutdown()

