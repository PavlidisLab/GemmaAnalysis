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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.lang.time.StopWatch;

import ubic.gemma.apps.ExpressionExperimentManipulatingCLI;
import ubic.gemma.model.genome.Gene;
import ubic.gemma.model.genome.PredictedGene;
import ubic.gemma.model.genome.ProbeAlignedRegion;

/**
 * <p>
 * The reason to make a temporary table; In Gemma, the link
 * tables store each links twice (the duplicate one with firstDesignElement and
 * secondDesignElement switched) to speed up the online co-expression query.
 * Some huge expression experiments give rise to a large amount of links more
 * than 10M. However, the shuffling need to go through all expression
 * experiments one by one to extract all links for each expression experiment
 * and this process is required to repeat many times (default is 100) to get
 * better estimation on the background distribution. Therefore, to speed up the
 * shuffling process, the first step will create a new table to save the links
 * without redundancy. It could also do some filtering (only save links for
 * known genes). Then the next step will do the shuffling on the working table,
 * which runs much faster.
 * 
 * @author xwan
 * @version $Id$
 */
public class LinkShuffleTempTableBuilderCLI extends ExpressionExperimentManipulatingCLI {

	public static void main(String[] args) {
		LinkShuffleTempTableBuilderCLI shuffle = new LinkShuffleTempTableBuilderCLI();
		StopWatch watch = new StopWatch();
		watch.start();
		try {
			Exception ex = shuffle.doWork(args);
			if (ex != null) {
				ex.printStackTrace();
			}
			watch.stop();
			log.info(watch);
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

	private boolean filterNonSpecific = true;

	@SuppressWarnings("static-access")
	@Override
	protected void buildOptions() {
		super.buildOptions();
		Option filterOption = OptionBuilder.withArgName("noFilterNonSpecific")
				.withDescription("Turn off filtering of non-specific probes from table")
				.withLongOpt("noFilter").create("nf");
		addOption(filterOption);

	}

	@SuppressWarnings("unchecked")
	@Override
	protected Exception doWork(String[] args) {
		Exception err = processCommandLine("Shuffle Links ", args);
		if (err != null) {
			return err;
		}

		LinkStatisticsService lss = (LinkStatisticsService) getBean("linkStatisticsService");

		lss.prepareDatabase(expressionExperiments, taxon.getCommonName(),
				filterNonSpecific);
		return null;
	}

	/**
	 * 
	 */
	@Override
	protected void processOptions() {
		super.processOptions();
		
		if (hasOption("nf")) 
		    filterNonSpecific = false;

	}

	protected String[] getAdditionalSpringConfigLocations() {
		return new String[] { "classpath*:chibi/beans.xml" };
	}

}
