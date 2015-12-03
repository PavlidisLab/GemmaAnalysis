package ubic.gemma.script.example

import java.util.concurrent.ForkJoinPool;

import ubic.gemma.script.framework.SpringSupport;
import ubic.gemma.tasks.visualization.DifferentialExpressionVisualizationValueObject.GeneScore;
import ubic.gemma.model.analysis.expression.diff.Direction;
import ubic.gemma.model.analysis.expression.diff.ContrastVO;
import ubic.gemma.model.analysis.expression.diff.ContrastsValueObject;
import ubic.gemma.tasks.visualization.DifferentialExpressionSearchTaskCommand;

CITATION = """\
# Generated by Gemma 
# 
# If you use this file for your research, please cite: 
# Zoubarev, A., et al., Gemma: A resource for the re-use, sharing and meta-analysis of expression profiling data. Bioinformatics, 2012. 
# 
# This functionality is currently in beta. The file format may change in the near future. 
# Fields are separated by tabs and delimited with double quotes 
# 
"""

// flags for outputting results to diff files
P = 0;
CORR_P = 1;
LOG_FC = 2;

cli = new CliBuilder()
cli.h(longOpt: 'help', 'Usage: groovy DiffExprMatrix -togiv')
cli.o(argName: 'file name', longOpt: 'outFile', required: true, args: 1, 'Output results to this file')
cli.t(argName: 'common name', longOpt: 'taxon', required: true, args: 1, 'Taxon of genes to fetch')
cli.g(argName: 'file name', longOpt: 'geneFile', required: false, args: 1, 'File containing list of gene official symbols to load')
cli.i(longOpt: 'filterNonSpecific', 'Filter non-specific probes')
cli.e(argName: 'file name', longOpt: 'eeFile', required: false, args: 1, 'File containing list of data sets to load')

cli.c(longOpt: 'correctedPval', required: false, 'Output corrected p-values')
cli.d(longOpt: 'direction', required: false, 'Output direction')

opts = cli.parse(args)
if (!opts) return;
if (opts.hasOption("h")) cli.usage();

geneSymbols = (opts.hasOption("g"))? new File(opts.getOptionValue("g")).readLines() : null;
filterNonSpecific = opts.i;
if (filterNonSpecific) {
    System.out.println "Filtering non-specific probes";
}
outputDirection = opts.d;
outputCorrectedPval = opts.c;
if (outputDirection) {
    System.out.println "Retrieving effect sizes";
} else if (outputCorrectedPval) {
    System.out.println "Retrieving corrected p-values";
} else {
    System.out.println "Retrieving (uncorrected) p-values";
}

sx = new SpringSupport();
taxonService = sx.getBean("taxonService");
geneService = sx.getBean("geneService");
deaService = sx.getBean("differentialExpressionAnalysisService")
gdeService = sx.getBean("geneDifferentialExpressionService")
csService = sx.getBean("compositeSequenceService")
eeService = sx.getBean("expressionExperimentService")
derService = sx.getBean( "differentialExpressionResultService" )
deSearchTask = sx.getBean("differentialExpressionSearchTask");

taxonName = opts.getOptionValue("t");
taxon = taxonService.findByCommonName(taxonName);

eeNames = (opts.hasOption("e"))? new File(opts.getOptionValue("e")).readLines() : null;

if (geneSymbols != null && geneSymbols.size() > 0) {
    System.out.println "${new Date()}: Attempting to load ${geneSymbols.size()} $taxonName genes...";
    genes = geneSymbols.collect { geneService.findByOfficialSymbol(it, taxon) };
    genes = genes.findAll { it != null };
} else {
    System.out.println "${new Date()}: Loading all known $taxonName genes...";
    //genes = geneService.loadKnownGenes(taxon);
    genes = geneService.getGenesByTaxon(taxon);
}
geneIds = genes.collect { it.id };

System.out.println "${new Date()}: Loaded ${genes.size()} $taxonName genes.";

ees = eeNames.collect { eeService.findByShortName(it) };
eeIds = ees.collect { it.id };

// VOs
genes = geneService.loadValueObjects( geneIds );
ordered = true;
experiments = eeService.loadValueObjects( eeIds, ordered );


//DifferentialExpressionSearchTaskCommand 
taskCommand = new DifferentialExpressionSearchTaskCommand( genes, experiments, "my gene set", "my ee set" );
deSearchTask.setTaskCommand( taskCommand );
taskResult = deSearchTask.execute();
ex = taskResult.getException();
if ( ex != null ) {
	throw new Error(ex);
}

answer = taskResult.getAnswer(); // returns DifferentialExpressionGenesConditionsValueObject

def printMatrix( answer, prop ) {
	cellData = answer.getCellData();
	buf = new StringBuilder();
	buf.append("Gene");

	for( c in answer.getConditions() ) {
		buf.append( "\t" + c.getId() );
	}
	buf.append("\n");
	
	for( g in answer.getGenes() ) {
		buf.append( g.getName() );
		for( c in answer.getConditions() ) {
			cell = cellData.get( c.getId() ).get( g.getId() );
			if ( prop == LOG_FC ) {
				val = cell.getLogFoldChange();
			} else if ( prop == P ) {
				val = cell.getpValue();
			} else if ( prop == CORR_P ) {
				val = cell.getCorrectedPValue();
			} else {
				throw new Error("Property ${prop} is not supported.");
			}
			//if ( val == null ) {
			//	val = "";
			//}
			// null values are printed as "nu"
			buf.append( String.format( "\t%.2f", val ) );
		}
		buf.append( "\n" );
	}
	
	return buf.toString();
}

def printConditions( answer ) {
	buf = new StringBuilder();

	buf.append( "Id\tDatasetShortName\tDatasetName\tContrastFactorValue\tFactorName\tBaselineFactorValue\n" );

	for( c in answer.getConditions() ) {
		buf.append( String.format( "%s\t%s\t%s\t%s\t%s\t%s\n", c.getId(), c.getDatasetShortName(), c.getDatasetName(), c.getContrastFactorValue(), c.getFactorName(), c.getBaselineFactorValue() ) );
	}
	
	return buf.toString();
}

def save( str, fileName ) {
	outFile = new File( fileName );
	outFile.delete();
	fOut = new BufferedWriter(new PrintWriter(outFile));
	fOut << CITATION;
	fOut << str;
	fOut.close();
	println "${new Date()}: Written results to $fileName";
}

println "===========Conditions================";
s = printConditions( answer );
//println "${s}";
save( "# Expression experiment metadata\n" + s, opts.getOptionValue("o") + "-eeMetaData.txt" );

println "===========PValue===========";
s = printMatrix( answer, P );
//println "${s}";
save( "# P-value\n" + s, opts.getOptionValue("o") + "-pvalue.txt" );

println "===========Corrected PValue===========";
s = printMatrix( answer, CORR_P );
//println "${s}";
save( "# Corrected p-value\n" + s, opts.getOptionValue("o") + "-qvalue.txt" );

println "===========Log Fold Change===========";
s = printMatrix( answer, LOG_FC );
//println "${s}";
save( "# Log fold change\n" + s, opts.getOptionValue("o") + "-logFC.txt" );

//println "Answer ${answer}";
