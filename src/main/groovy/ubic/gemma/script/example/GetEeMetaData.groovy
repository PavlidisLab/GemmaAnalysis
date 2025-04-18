package ubic.gemma.script.example

import groovy.cli.commons.CliBuilder
import ubic.gemma.groovy.framework.SpringSupport
import ubic.gemma.persistence.service.expression.experiment.ExpressionExperimentService
import ubic.gemma.persistence.service.genome.gene.GeneService
import ubic.gemma.persistence.service.genome.taxon.TaxonService

cli = new CliBuilder()
cli.h(longOpt: 'help', 'Usage: groovy GetEeMetaData -o')
cli.f(argName: 'file name', longOpt: 'inFile', required: false, args: 1, '')
cli.o(argName: 'file name', longOpt: 'outFile prefix', required: true, args: 1, '')
cli.t(argName: 'taxon', longOpt: 'taxon', required: true, args: 1, '')
cli.e(argName: 'file name', longOpt: 'eeFile', required: false, args: 1, 'File containing list of data sets to load')

opts = cli.parse(args)
if (!opts) return
if (opts.h) cli.usage()
if (!opts.t) cli.usage()

sx = new SpringSupport()
taxonService = sx.getBean(TaxonService.class)
eeService = sx.getBean(ExpressionExperimentService.class)
geneService = sx.getBean(GeneService.class)

f = opts.o as String
outFile = new File(f)
outFile.delete()

/**
 * Load taxon
 */
taxonName = opts.t as String
taxon = taxonService.findByCommonName(taxonName)

/**
 * Load genes
 */
genes = geneService.loadAll(taxon)
System.out.println "${new Date()}: Loaded ${genes.size()} $taxonName genes ...)."
f = opts.o + "-geneMetaData.txt"
println "${new Date()}: Writing results to $f"
//println "Only genes with an ENSEMBL_ID will be written"
fOut = new BufferedWriter(new FileWriter(f))
fOut << "ID\tOFFICIAL_SYMBOL\tENSEMBL_ID\tNCBI_ID\tTAXON\tOFFICIAL_NAME\n"
count = 1
for (gene in genes) {
    //if(gene.ensemblId == null) continue
    line = "${gene.id}\t${gene.officialSymbol}\t${gene.ensemblId}\t${gene.ncbiId}\t${gene.taxon.commonName}\t\"${gene.officialName}\"\n"
    fOut << line
    count++
}
println "${new Date()}: Wrote $count lines"
fOut.close()

/**
 * Load list of ees from a file or from Gemma (for a given taxon)
 */
//ees = eeService.loadAll()
eeNames = (opts.e) ? new File(opts.e).readLines() : null
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

f = opts.o + "-eeMetaData.txt"
println "${new Date()}: Writing results to $f"
out = new BufferedWriter(new FileWriter(f))
out << "ARSID\tName\tDescription\tValues\teeShortName\teeName\tTaxon\tNumSamples\tArrayDesignNames\tArrayDesignSizes\n"
count = 1
for (ee in ees) {

    //taxon = eeService.getTaxon(ee);
    //if (!taxon.getCommonName().equals(opts.getOptionValue('t'))) {
    //	continue;
    //}

    numSamples = eeService.getBioMaterialCount(ee)
    ads = eeService.getArrayDesignsUsed(ee)
    adSizes = ads*.advertisedNumberOfDesignElements
    adNames = ads*.shortName

    ee = eeService.thawLite(ee)
    des = ee.experimentalDesign
    efs = des.experimentalFactors


    //efs = efService.thaw(design.experimentalFactors.iterator().next());
    //efs = ee.experimentalDesign.experimentalFactors;
    //efNames = efs*.name;
    //System.out.println("efNames $efNames")

    for (ef in efs) {
        vals = ef.factorValues*.value

        //out << "${ee.shortName}\t${ee.name}\t${taxon.commonName}\t${numSamples}\t${adNames.join('|')}\t${adSizes.join('|')}\t${efNames.join('|')}\n"
        //line = "${ee.shortName}\t${ee.name}\t${taxon.commonName}\t${numSamples}\t${adNames.join('|')}\t${adSizes.join('|')}\t${efNames.join('|')}\n"
        //line = "${ee.shortName}\t${ee.name}\t${taxon.commonName}\t${numSamples}\t${adNames.join('|')}\t${adSizes.join('|')}\t${des.normalizationDescription}\n";
        line = "ARSID.${ee.id}.${ef.id}\t\"${ef.name}\"\t\"${ef.description}\"\t\"${vals.join('|')}\"\t${ee.shortName}\t\"${ee.name}\"\t${taxon.commonName}\t${numSamples}\t${adNames.join('|')}\t${adSizes.join('|')}\n"
        //System.out.print "$line"
        out << line

        if (count++ % 50 == 0) {
            System.out.println "${new Date()}: Processed $count experimental factors"
        }
    }

    //TODO
    /*if ( count >= 20 ) {
        break;
    }*/
}
out.flush()
out.close()


sx.shutdown()

