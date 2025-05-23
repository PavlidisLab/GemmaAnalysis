#!/usr/bin/groovy
package ubic.gemma.script.example

import groovy.cli.commons.CliBuilder
import groovy.time.TimeCategory
import groovy.time.TimeDuration
import ubic.gemma.core.analysis.preprocess.MeanVarianceService
import ubic.gemma.model.expression.experiment.ExpressionExperiment
import ubic.gemma.persistence.service.expression.experiment.ExpressionExperimentService

/* Parse arguments */
def cli = new CliBuilder(usage: 'groovy CreateMeanVariance [options] -u *** -p *** -e ***')
cli.p(longOpt: 'password', args: 1, 'password', required: true)
cli.u(longOpt: 'username', args: 1, 'username', required: true)
cli.e(args: 1, 'file that contains EE shortNames, one per line', required: true)
cli.r(longOpt: 'recompute', 'recompute mean-variance?', required: false)
cli.h(longOpt: 'help', 'usage information', required: false)
def opt = cli.parse(args)

username = System.getenv('GEMMA_USERNAME')
password = System.getenv('GEMMA_PASSWORD')
forceRecompute = opt.r
infile = opt.e as String
eeIds = new File(infile).readLines()

/* Do the actual work */
sx = new SpringSupport(username, password)
ees = sx.getBean(ExpressionExperimentService.class)
mvs = sx.getBean(MeanVarianceService.class)

pass = 0
fail = []
for (id in eeIds) {
    ExpressionExperiment ee = ees.findByShortName(id)

    if (ee == null) {
        System.err.println('ERROR: Could not find experiment ' + id)
        fail.add(id)
        continue
    }

    print "Processing mean-variance for experiment " + ee.getShortName() + ' ...'

    try {
        def timeStart = new Date()
        mvs.create(ee, forceRecompute)
        def timeStop = new Date()
        TimeDuration duration = TimeCategory.minus(timeStop, timeStart)
        println 'It took ' + duration
        pass++
    } catch (e) {
        System.err.println('ERROR: Error computing mean-variance')
        e.printStackTrace()
        fail.add(id)
    }
}

println 'Finished processing ' + pass + ' experiments.'

if (fail.size() > 0) {
    println 'Failed to process these experiments ' + fail
}

