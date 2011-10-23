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

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.lang.time.StopWatch;

import ubic.basecode.io.ByteArrayConverter;
import ubic.basecode.math.DescriptiveWithMissing;
import ubic.gemma.model.common.quantitationtype.QuantitationType;
import ubic.gemma.model.common.quantitationtype.StandardQuantitationType;
import ubic.gemma.model.expression.arrayDesign.ArrayDesign;
import ubic.gemma.model.expression.arrayDesign.ArrayDesignService;
import ubic.gemma.model.expression.bioAssay.BioAssay;
import ubic.gemma.model.expression.bioAssayData.DesignElementDataVector;
import ubic.gemma.model.expression.bioAssayData.DesignElementDataVectorService;
import ubic.gemma.model.expression.bioAssayData.ProcessedExpressionDataVector;
import ubic.gemma.model.expression.designElement.CompositeSequence;
import ubic.gemma.model.expression.designElement.CompositeSequenceService;
import ubic.gemma.model.expression.experiment.ExpressionExperiment;
import ubic.gemma.model.expression.experiment.ExpressionExperimentService;
import ubic.gemma.model.genome.Gene;
import ubic.gemma.util.AbstractSpringAwareCLI;
import cern.colt.list.DoubleArrayList;
import cern.colt.list.ObjectArrayList;
import cern.jet.stat.Descriptive;

/**
 * @author xwan
 * @version $Id$
 */
public class AffyPlatFormAnalysisCli extends AbstractSpringAwareCLI {

    private class SortedElement implements Comparable<SortedElement> {
        private Double max, median, presentAbsentCall;
        private CompositeSequence de;

        public SortedElement( CompositeSequence de, double max, double median, double presentAbsentCall ) {
            this.de = de;
            this.max = max;
            this.median = median;
            this.presentAbsentCall = presentAbsentCall;
        }

        public int compareTo( SortedElement o ) {
            return median.compareTo( o.median );
        }

        public Double getMax() {
            return max;
        }

        public CompositeSequence getDE() {
            return de;
        }

        public Double getPresentAbsentCall() {
            return presentAbsentCall;
        }
    }

    public static final int MIN = 1;
    public static final int MAX = 2;
    public static final int MEDIAN = 3;
    public static final int MEAN = 4;
    public static final int STD = 5;

    private String arrayDesignName = null;
    private String outFileName = null;
    private Map<CompositeSequence, DoubleArrayList> rankData = new HashMap<CompositeSequence, DoubleArrayList>();
    private Map<CompositeSequence, DoubleArrayList> presentAbsentData = new HashMap<CompositeSequence, DoubleArrayList>();
    private DesignElementDataVectorService devService = null;
    private ExpressionExperimentService eeService = null;
    private Map<CompositeSequence, Collection<Gene>> probeToGeneAssociation = null;

    @Override
    protected void processOptions() {
        super.processOptions();
        if ( hasOption( 'a' ) ) {
            this.arrayDesignName = getOptionValue( 'a' );
        }
        if ( hasOption( 'o' ) ) {
            this.outFileName = getOptionValue( 'o' );
        }
    }

    @SuppressWarnings("static-access")
    @Override
    protected void buildOptions() {
        Option ADOption = OptionBuilder.hasArg().isRequired().withArgName( "arrayDesign" )
                .withDescription( "Array Design Short Name (GPLXXX) " ).withLongOpt( "arrayDesign" ).create( 'a' );
        addOption( ADOption );
        Option OutOption = OptionBuilder.hasArg().isRequired().withArgName( "outputFile" )
                .withDescription( "The name of the file to save the output " ).withLongOpt( "outputFile" ).create( 'o' );
        addOption( OutOption );
    }

    private QuantitationType getQuantitationType( ExpressionExperiment ee, StandardQuantitationType requiredQT,
            boolean isPreferedQT ) {
        QuantitationType qtf = null;
        Collection<QuantitationType> eeQT = this.eeService.getQuantitationTypes( ee );
        for ( QuantitationType qt : eeQT ) {
            if ( isPreferedQT ) {
                if ( qt.getIsPreferred() ) {
                    qtf = qt;
                    StandardQuantitationType tmpQT = qt.getType();
                    if ( tmpQT != StandardQuantitationType.AMOUNT ) {
                        log.warn( "Preferred Quantitation Type may not be correct." + ee.getShortName() + ":"
                                + tmpQT.toString() );
                    }
                    break;
                }
            } else {
                if ( qt.getType().equals( requiredQT ) ) {
                    qtf = qt;
                    break;
                }
            }
        }
        if ( qtf == null ) {
            log.info( "Expression Experiment " + ee.getShortName() + " doesn't have required QT " );
        }
        return qtf;
    }

    String processEEForPercentage( ExpressionExperiment ee ) {
        // eeService.thaw( ee );
        QuantitationType qt = this.getQuantitationType( ee, StandardQuantitationType.PRESENTABSENT, false );
        if ( qt == null ) return ( "No usable quantitation type in " + ee.getShortName() );
        log.info( "Load Data for  " + ee.getShortName() );

        Collection<? extends DesignElementDataVector> dataVectors = devService.find( qt );
        if ( dataVectors == null ) return ( "No data vector " + ee.getShortName() );
        ByteArrayConverter bac = new ByteArrayConverter();

        for ( DesignElementDataVector vector : dataVectors ) {
            CompositeSequence de = vector.getDesignElement();
            DoubleArrayList presentAbsentList = this.presentAbsentData.get( de );
            if ( presentAbsentList == null ) {
                // return (" EE data vectors don't match array design for probe " + de.getName());
                continue;
            }
            byte[] bytes = vector.getData();
            // String vals = bac.byteArrayToAsciiString( bytes );
            char[] chars = bac.byteArrayToChars( bytes );
            double presents = 0;
            double total = 0;
            for ( int i = 0; i < chars.length; i++ ) {
                if ( chars[i] == 'P' ) presents++;
                if ( chars[i] != '\t' ) total++;
            }
            presentAbsentList.add( presents / total );
        }
        return null;
    }

    private Map<CompositeSequence, Collection<Gene>> getDevToGeneAssociation(
            Collection<ProcessedExpressionDataVector> datavectors ) {

        Collection<CompositeSequence> cs = new HashSet<CompositeSequence>();
        for ( ProcessedExpressionDataVector designElementDataVector : datavectors ) {
            cs.add( designElementDataVector.getDesignElement() );
        }
        CompositeSequenceService css = ( CompositeSequenceService ) this.getBean( "compositeSequenceService" );
        return css.getGenes( cs );
    }

    String processEE( ExpressionExperiment ee ) {
        // eeService.thaw( ee );
        QuantitationType qt = this.getQuantitationType( ee, null, true );
        if ( qt == null ) return ( "No usable quantitation type in " + ee.getShortName() );
        log.info( "Load Data for  " + ee.getShortName() );

        Collection<ProcessedExpressionDataVector> dataVectors = eeService.getProcessedDataVectors( ee );
        if ( dataVectors == null ) return ( "No data vector " + ee.getShortName() );
        if ( this.probeToGeneAssociation == null ) {
            this.probeToGeneAssociation = this.getDevToGeneAssociation( dataVectors );
        }

        for ( ProcessedExpressionDataVector vector : dataVectors ) {
            CompositeSequence de = vector.getDesignElement();
            DoubleArrayList rankList = this.rankData.get( de );
            if ( rankList == null ) {
                return ( " EE data vectors don't match array design for probe " + de.getName() );
            }
            Double rank = vector.getRankByMean();
            if ( rank != null ) {
                rankList.add( rank.doubleValue() );
            }
        }
        return null;
    }

    private double getStatValue( DoubleArrayList valList, int method ) {
        double value = 0.0;
        switch ( method ) {
            case MIN:
                value = DescriptiveWithMissing.min( valList );
                break;
            case MAX:
                value = DescriptiveWithMissing.max( valList );
                break;
            case MEAN:
                value = DescriptiveWithMissing.mean( valList );
                break;
            case MEDIAN:
                valList.sort();
                value = DescriptiveWithMissing.median( valList );
                break;
            case STD:
                int N = valList.size();
                double sum = DescriptiveWithMissing.sum( valList );
                double ss = DescriptiveWithMissing.sumOfSquares( valList );
                value = Descriptive.standardDeviation( DescriptiveWithMissing.variance( N, sum, ss ) );
                break;
        }
        if ( Double.isNaN( value ) ) value = 0.0;
        return value;
    }

    int getNumberofArraysinEE( ExpressionExperiment ee, ArrayDesign ad ) {
        int numberofArrays = 0;
        eeService.thawLite( ee );
        Collection<BioAssay> bioAssays = ee.getBioAssays();
        for ( BioAssay assay : bioAssays ) {
            ArrayDesign design = assay.getArrayDesignUsed();
            if ( ad.equals( design ) ) {
                numberofArrays++;
            }
        }
        System.err.println( "Got " + numberofArrays );
        return numberofArrays;
    }

    @Override
    protected Exception doWork( String[] args ) {
        Exception err = processCommandLine( "AffYPlatFormAnalysisCli ", args );
        if ( err != null ) {
            return err;
        }
        ArrayDesignService adService = ( ArrayDesignService ) this.getBean( "arrayDesignService" );
        this.eeService = ( ExpressionExperimentService ) this.getBean( "expressionExperimentService" );
        this.devService = ( DesignElementDataVectorService ) this.getBean( "designElementDataVectorService" );

        ArrayDesign arrayDesign = adService.findByShortName( this.arrayDesignName );
        if ( arrayDesign == null ) {
            System.err.println( " Array Design " + this.arrayDesignName + " doesn't exist" );
            return null;
        }

        // adService.thaw(arrayDesign);
        Collection<CompositeSequence> allCSs = adService.loadCompositeSequences( arrayDesign );
        for ( CompositeSequence cs : allCSs ) {
            this.rankData.put( cs, new DoubleArrayList() );
            this.presentAbsentData.put( cs, new DoubleArrayList() );
        }

        Collection<ExpressionExperiment> relatedEEs = adService.getExpressionExperiments( arrayDesign );

        int numberofAllArrays = 0;
        for ( ExpressionExperiment ee : relatedEEs ) {
            System.err.println( ee.getName() );
            if ( this.processEEForPercentage( ee ) != null ) continue;

            if ( this.processEE( ee ) != null ) continue;
            numberofAllArrays = numberofAllArrays + this.getNumberofArraysinEE( ee, arrayDesign );
        }
        log.info( "The total number of all arrays is " + numberofAllArrays );
        ObjectArrayList sortedList = new ObjectArrayList();
        for ( CompositeSequence de : this.rankData.keySet() ) {
            DoubleArrayList rankList = this.rankData.get( de );
            DoubleArrayList presentAbsentList = this.presentAbsentData.get( de );
            if ( rankList.size() > 0 ) {
                SortedElement oneElement = new SortedElement( de, getStatValue( rankList, MAX ), getStatValue(
                        rankList, MEDIAN ), getStatValue( presentAbsentList, MEAN ) );
                sortedList.add( oneElement );
            } else {
                System.err.print( de.getName() );
                System.err.println( " Empty " );
            }
        }
        sortedList.sort();
        try {
            PrintStream output = new PrintStream( new FileOutputStream( new File( this.outFileName ) ) );
            for ( int i = 0; i < sortedList.size(); i++ ) {
                SortedElement oneElement = ( SortedElement ) sortedList.get( i );
                if ( oneElement.getDE().getName().toUpperCase().contains( "AFFY" ) ) continue;
                output.print( oneElement.getDE().getName() );
                output.print( "\t" + oneElement.getMax() );
                double lower_threshould = 0.001;
                if ( oneElement.getPresentAbsentCall() < lower_threshould )
                    output.print( "\t" + lower_threshould );
                else
                    output.print( "\t" + oneElement.getPresentAbsentCall() );
                Collection<Gene> mappedGenes = this.probeToGeneAssociation.get( oneElement.getDE() );
                if ( mappedGenes != null )
                    output.println( "\t" + mappedGenes.size() );
                else
                    output.println( "\t0" );
            }
            output.close();
        } catch ( Exception e ) {
            return e;
        }
        return null;
    }

    /**
     * @param args
     */
    public static void main( String[] args ) {
        AffyPlatFormAnalysisCli analysis = new AffyPlatFormAnalysisCli();
        StopWatch watch = new StopWatch();
        watch.start();
        try {
            Exception ex = analysis.doWork( args );
            if ( ex != null ) {
                ex.printStackTrace();
            }
            watch.stop();
            log.info( "Elapsed time: " + watch.getTime() / 1000 + " seconds" );
        } catch ( Exception e ) {
            throw new RuntimeException( e );
        }
    }

}
