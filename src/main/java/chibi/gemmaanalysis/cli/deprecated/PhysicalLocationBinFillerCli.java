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

import java.util.Collection;

import ubic.gemma.core.util.AbstractSpringAwareCLI;
import ubic.gemma.model.genome.PhysicalLocation;
import ubic.gemma.persistence.service.genome.PhysicalLocationDao;
import ubic.gemma.persistence.util.SequenceBinUtils;

/**
 * This is a one-off.
 *
 * @author pavlidis
 * @version $Id: PhysicalLocationBinFillerCli.java,v 1.5 2015/11/12 19:37:12 paul Exp $
 */
@Deprecated
public class PhysicalLocationBinFillerCli extends AbstractSpringAwareCLI {

    public static void main( String[] args ) {
        PhysicalLocationBinFillerCli p = new PhysicalLocationBinFillerCli();
        p.doWork( args );
    }

    /*
     * (non-Javadoc)
     *
     * @see ubic.gemma.util.AbstractCLI#getCommandName()
     */
    @Override
    public String getCommandName() {
        return null;
    }

    @Override
    protected void buildOptions() {
        //
    }

    @Override
    protected Exception doWork( String[] args ) {
        processCommandLine( args );

        PhysicalLocationDao pld = this.getBean( PhysicalLocationDao.class );

        Collection<? extends PhysicalLocation> lcs = pld.loadAll();
        int count = 0;
        for ( PhysicalLocation location : lcs ) {
            if ( location.getNucleotide() == null || location.getNucleotideLength() == null ) continue;
            int bin = SequenceBinUtils.binFromRange( location.getNucleotide().intValue(), location.getNucleotide()
                    .intValue() + location.getNucleotideLength().intValue() );
            location.setBin( bin );
            pld.update( location );
            if ( ++count % 10000 == 0 ) {
                log.info( "Processed " + count );
            }
        }

        return null;
    }
}
