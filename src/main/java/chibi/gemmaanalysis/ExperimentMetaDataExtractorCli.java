/*
 * The GemmaAnalysis project
 *
 * Copyright (c) 2014 University of British Columbia
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
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.io.Writer;
import java.text.DateFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.lang.time.StopWatch;

import gemma.gsec.SecurityService;
import ubic.basecode.dataStructure.CountingMap;
import ubic.gemma.core.analysis.preprocess.OutlierDetails;
import ubic.gemma.core.analysis.preprocess.OutlierDetectionService;
import ubic.gemma.core.analysis.preprocess.batcheffects.BatchEffectDetails;
import ubic.gemma.core.analysis.preprocess.batcheffects.BatchInfoPopulationServiceImpl;
import ubic.gemma.core.analysis.preprocess.filter.FilterConfig;
import ubic.gemma.core.analysis.preprocess.svd.SVDService;
import ubic.gemma.core.analysis.preprocess.svd.SVDValueObject;
import ubic.gemma.core.analysis.service.ExpressionDataMatrixService;
import ubic.gemma.core.analysis.util.ExperimentalDesignUtils;
import ubic.gemma.core.apps.ExpressionExperimentManipulatingCLI;
import ubic.gemma.core.apps.GemmaCLI.CommandGroup;
import ubic.gemma.core.expression.experiment.service.ExperimentalDesignService;
import ubic.gemma.model.common.description.BibliographicReference;
import ubic.gemma.model.common.quantitationtype.QuantitationType;
import ubic.gemma.model.expression.arrayDesign.ArrayDesign;
import ubic.gemma.model.expression.bioAssay.BioAssay;
import ubic.gemma.model.expression.experiment.BioAssaySet;
import ubic.gemma.model.expression.experiment.ExperimentalFactor;
import ubic.gemma.model.expression.experiment.ExpressionExperiment;
import ubic.gemma.model.expression.experiment.ExpressionExperimentSubSet;
import ubic.gemma.model.expression.experiment.ExpressionExperimentValueObject;
import ubic.gemma.persistence.service.common.auditAndSecurity.AuditTrailService;
import ubic.gemma.persistence.util.FactorValueVector;
import ubic.gemma.persistence.util.Settings;

/**
 * Extracts expression experiment meta data such as experimental design, array design, outlier count, and publication
 * into a .txt.gz TSV file. See Bug 3968.
 *
 * @author paul
 * @version $Id: ExperimentMetaDataExtractorCli.java,v 1.16 2015/11/12 19:37:12 paul Exp $
 */
public class ExperimentMetaDataExtractorCli extends ExpressionExperimentManipulatingCLI {

    private static final String EXPERIMENT_META_DATA_BASENAME = "ExperimentMetaData";
    private static final String VIEW_FILE_SUFFIX = ".txt.gz";
    public static final String DEFAULT_VIEW_FILE = Settings.getString( "gemma.appdata.home" ) + File.separatorChar
            + "dataFiles" + File.separatorChar + EXPERIMENT_META_DATA_BASENAME + VIEW_FILE_SUFFIX;
    private static final String NA = "";

    public static void main( String[] args ) {
        ExperimentMetaDataExtractorCli s = new ExperimentMetaDataExtractorCli();
        try {
            Exception ex = s.doWork( args );
            if ( ex != null ) {
                ex.printStackTrace();
            }

            System.exit( 0 );
        } catch ( Exception e ) {
            // throw new RuntimeException( e );
            log.error( e.getMessage(), e );
            System.exit( 1 );
        }
    }

    private OutlierDetectionService outlierDetectionService;
    //   private StatusService statusService;
    private ExperimentalDesignService edService;
    private SecurityService securityService;
    private String viewFile = DEFAULT_VIEW_FILE;
    private SVDService svdService;

    private ExpressionDataMatrixService expressionDataMatrixService;

    public void generateExperimentMetaData( Collection<BioAssaySet> expressionExperiments ) throws IOException {

        File file = getOutputFile( this.viewFile );

        try (Writer writer = new OutputStreamWriter( new GZIPOutputStream( new FileOutputStream( file ) ) );) {

            String[] colNames = { "ShortName", "Taxon", "DateUpload", "IsPublic", "NumPlatform", "Platform", "Channel",
                    "IsExonArray", "QtIsRatio", "QtIsNormalized", "QtScale", "NumProfiles", "NumFilteredProfiles",
                    "NumSamples", "NumConditions", "NumReplicatesPerCondition", "PossibleOutliers", "CuratedOutlier",
                    "IsTroubled", "PubTroubled", "PubYear", "PubJournal", "Batch.PC1.Var", "Batch.PC2.Var",
                    "Batch.PC3.Var", "Batch.PC1.Pval", "Batch.PC2.Pval", "Batch.PC3.Pval", "NumFactors", "FactorNames",
                    "FactorCategories", "NumFactorValues" };
            // log.info( StringUtils.join( colNames, "\t" ) + "\n" );
            writer.write( StringUtils.join( colNames, "\t" ) + "\n" );

            int i = 0;
            Collection<String> failedEEs = new ArrayList<>();

            StopWatch timer = new StopWatch();
            timer.start();

            for ( BioAssaySet bas : expressionExperiments ) {
                /*
                 * Skip subsets
                 */
                if ( bas instanceof ExpressionExperimentSubSet ) return;

                try {

                    ExpressionExperiment ee = ( ExpressionExperiment ) bas;
                    ee = eeService.thawLite( ee );
                    ExpressionExperimentValueObject vo = eeService.loadValueObject( ee.getId() );
                    vo.setIsPublic( !securityService.isPrivate( ee ) );
                    log.info( "Processing (" + ++i + "/" + expressionExperiments.size() + ") : " + ee );

                    BibliographicReference primaryPublication = ee.getPrimaryPublication();

                    Collection<ArrayDesign> arrayDesignsUsed = eeService.getArrayDesignsUsed( ee );

                    Collection<String> arrayDesignIsExon = new ArrayList<>();
                    Collection<String> arrayDesignTechTypes = new ArrayList<>();
                    Collection<String> arrayDesignShortNames = new ArrayList<>();

                    // for multiple platforms e.g. GSE5949
                    for ( ArrayDesign ad : arrayDesignsUsed ) {
                        arrayDesignShortNames.add( ad.getShortName() );
                        arrayDesignTechTypes.add( ad.getTechnologyType().getValue() );
                        arrayDesignIsExon.add( ad.getName().toLowerCase().contains( "exon" ) + "" );
                    }

                    QuantitationType qt = null;
                    for ( QuantitationType q : ee.getQuantitationTypes() ) {
                        if ( q.getIsPreferred().booleanValue() ) {
                            qt = q;
                            break;
                        }
                    }

                    int manualOutlierCount = 0;
                    for ( BioAssay ba : ee.getBioAssays() ) {
                        if ( ba.getIsOutlier().booleanValue() ) {
                            manualOutlierCount++;
                        }
                    }

                    //  ExperimentalDesign experimentalDesign = edService.load( vo.getExperimentalDesign() );

                    // Batch PCs
                    int maxcomp = 3;
                    BatchEffectDetails batchEffectPC1 = null;
                    BatchEffectDetails batchEffectPC2 = null;
                    BatchEffectDetails batchEffectPC3 = null;
                    Collection<BatchEffectDetails> batchEffects = getBatchEffect( ee, maxcomp );
                    Iterator<BatchEffectDetails> batchEffectsIterator;
                    if ( batchEffects == null || batchEffects.size() == 0 ) {
                        log.warn( "No batch effect info" );
                    } else {
                        batchEffectsIterator = batchEffects.iterator();
                        if ( batchEffectsIterator.hasNext() ) {
                            batchEffectPC1 = batchEffectsIterator.next();
                        }
                        if ( batchEffectsIterator.hasNext() ) {
                            batchEffectPC2 = batchEffectsIterator.next();
                        }
                        if ( batchEffectsIterator.hasNext() ) {
                            batchEffectPC3 = batchEffectsIterator.next();
                        }
                    }

                    // eeService.getExperimentsWithOutliers();
                    StopWatch timerOutlier = new StopWatch();
                    timerOutlier.start();
                    // log.info( "Outlier detection service started " + timer.getTime() + "ms" );
                    Collection<OutlierDetails> possibleOutliers = outlierDetectionService.identifyOutliers( ee );
                    if ( timerOutlier.getTime() > 10000 ) {
                        log.info( "Automatic outlier detection took " + timerOutlier.getTime() + "ms" );
                    }
                    // Collection<OutlierDetails> possibleOutliers = null;

                    // samples per condition
                    //      boolean removeBatchFactor = false;
                    Collection<String> samplesPerConditionCount = new ArrayList<>();
                    CountingMap<FactorValueVector> assayCount = ExperimentalDesignUtils.getDesignMatrix( ee );
                    List<FactorValueVector> keys = assayCount.sortedKeyList( true );
                    for ( FactorValueVector key : keys ) {
                        samplesPerConditionCount.add( Integer.toString( assayCount.get( key ).intValue() ) );
                    }

                    // factor names
                    Collection<ExperimentalFactor> factors = ee.getExperimentalDesign().getExperimentalFactors();
                    Collection<String> factorNames = new ArrayList<>();
                    Collection<String> factorCategories = new ArrayList<>();
                    Collection<Integer> factorValues = new ArrayList<>();
                    for ( ExperimentalFactor f : factors ) {
                        factorNames.add( f.getName() );
                        String cat = f.getCategory() != null ? f.getCategory().getCategory() : NA;
                        factorCategories.add( cat );
                        factorValues.add( Integer.valueOf( f.getFactorValues().size() ) );
                    }

                    int filteredProfilesCount = -1;

                    try {
                        FilterConfig filterConfig = new FilterConfig();
                        filterConfig.setIgnoreMinimumSampleThreshold( true );
                        filteredProfilesCount = expressionDataMatrixService.getFilteredMatrix( ee, filterConfig,
                                expressionDataMatrixService.getProcessedExpressionDataVectors( ee ) ).rows();
                    } catch ( Exception e ) {
                        log.error( e.getMessage(), e );
                    }

                    String val[] = {
                            vo.getShortName(),
                            vo.getTaxon(),

                            // FIXME this property is no longer in the VOs
                            //    DateFormat.getDateInstance( DateFormat.MEDIUM ).format( vo.getDateCreated() ),
                            null,

                            Boolean.toString( vo.getIsPublic() ),
                            Integer.toString( arrayDesignsUsed.size() ),
                            StringUtils.join( arrayDesignShortNames, ',' ),
                            StringUtils.join( arrayDesignTechTypes, ',' ), // arrayDesign.getTechnologyType().getValue(),
                            // ONE-COLOR, TWO-COLOR, NONE (RNA-seq
                            // GSE37646), DUAL-MODE (one or two color)
                            StringUtils.join( arrayDesignIsExon, ',' ), // exon GSE28383
                            qt != null ? Boolean.toString( qt.getIsRatio().booleanValue() ) : NA,
                            qt != null ? Boolean.toString( qt.getIsNormalized().booleanValue() ) : NA,
                            qt != null ? qt.getScale().getValue() : NA,
                            Integer.toString( vo.getProcessedExpressionVectorCount().intValue() ), // NumProfiles
                            Integer.toString( filteredProfilesCount ), // NumFilteredProfiles
                            Integer.toString( vo.getBioAssayCount().intValue() ), // NumSamples
                            Integer.toString( assayCount.size() ), // NumConditions
                            StringUtils.join( samplesPerConditionCount, "," ),
                            possibleOutliers != null ? Integer.toString( possibleOutliers.size() ) : NA,
                            Integer.toString( manualOutlierCount ),
                            Boolean.toString( vo.getTroubled() ),
                            primaryPublication != null ? DateFormat.getDateInstance( DateFormat.MEDIUM ).format(
                                    primaryPublication.getPublicationDate() ) : NA,
                            primaryPublication != null ? primaryPublication.getPublication() : NA,

                            batchEffectPC1 != null ? Double.toString( batchEffectPC1.getComponentVarianceProportion() )
                                    : NA,
                            batchEffectPC2 != null ? Double.toString( batchEffectPC2.getComponentVarianceProportion() )
                                    : NA,
                            batchEffectPC3 != null ? Double.toString( batchEffectPC3.getComponentVarianceProportion() )
                                    : NA,

                            batchEffectPC1 != null ? Double.toString( batchEffectPC1.getPvalue() ) : NA,
                            batchEffectPC2 != null ? Double.toString( batchEffectPC2.getPvalue() ) : NA,
                            batchEffectPC3 != null ? Double.toString( batchEffectPC3.getPvalue() ) : NA,

                            // factors
                            Integer.toString( factors.size() ),
                            StringUtils.join( factorNames, "," ),
                            StringUtils.join( factorCategories, "," ),
                            StringUtils.join( factorValues, "," ),

                    };

                    // log.info( StringUtils.join( val, "\t" ) + "\n" );
                    writer.write( StringUtils.join( val, "\t" ) + "\n" );

                } catch ( Exception e ) {
                    failedEEs.add( ( ( ExpressionExperiment ) bas ).getShortName() );
                    StringWriter sw = new StringWriter();
                    e.printStackTrace( new PrintWriter( sw ) );
                    log.error( sw.toString() );
                }
            }

            log.info( "Finished processing " + expressionExperiments.size() + " datasets in " + timer.getTime()
                    + " ms. " );
            log.info( "Writen to " + file );
            log.info( "Number of failed experiment metadata extraction(s): " + failedEEs.size() + " / "
                    + expressionExperiments.size() );

            if ( failedEEs.size() > 0 ) {
                log.info( "Skipped experiments:" );
                log.info( StringUtils.join( failedEEs, "," ) );
            }
        }
    }

    public Collection<BatchEffectDetails> getBatchEffect( ExpressionExperiment ee, int maxcomp ) {
        Collection<BatchEffectDetails> ret = new ArrayList<>();

        for ( ExperimentalFactor ef : ee.getExperimentalDesign().getExperimentalFactors() ) {
            if ( BatchInfoPopulationServiceImpl.isBatchFactor( ef ) ) {
                SVDValueObject svd = svdService.getSvdFactorAnalysis( ee.getId() );
                if ( svd == null ) break;

                for ( Integer component : svd.getFactorPvals().keySet() ) {
                    if ( component.intValue() >= maxcomp ) {
                        break;
                    }
                    Map<Long, Double> cmpEffects = svd.getFactorPvals().get( component );

                    Double pval = cmpEffects.get( ef.getId() );
                    if ( pval != null ) {
                        BatchEffectDetails details = new BatchEffectDetails();
                        details.setPvalue( pval.doubleValue() );
                        details.setComponent( new Integer( component.intValue() + 1 ) );
                        details.setComponentVarianceProportion( svd.getVariances()[component.intValue()].doubleValue() );
                        details.setHasBatchInformation( true );
                        ret.add( details );
                    }

                }
            }
        }
        return ret;
    }

    @Override
    public CommandGroup getCommandGroup() {
        return CommandGroup.METADATA;
    }

    /*
     * (non-Javadoc)
     *
     * @see ubic.gemma.util.AbstractCLI#getCommandName()
     */
    @Override
    public String getCommandName() {
        // TODO Auto-generated method stub
        return null;
    }

    public File getOutputFile( String filename ) {
        String fullFilePath = filename;
        File f = new File( fullFilePath );

        if ( f.exists() ) {
            return f;
        }

        File parentDir = f.getParentFile();
        if ( !parentDir.exists() ) parentDir.mkdirs();
        return f;
    }

    @Override
    @SuppressWarnings("static-access")
    protected void buildOptions() {
        super.buildOptions();

        OptionBuilder.hasArg();
        OptionBuilder.withArgName( "outfile" );
        OptionBuilder.withDescription( "GZipped output filename" );
        OptionBuilder
                .withLongOpt( "outfile" );
        Option expOption = OptionBuilder.create( 'o' );

        addOption( expOption );

        // to keep troubled experiments
        this.addForceOption();
    }

    /*
     * (non-Javadoc)
     *
     * @see ubic.gemma.util.AbstractCLI#doWork(java.lang.String[])
     */
    @Override
    protected Exception doWork( String[] args ) {
        super.processCommandLine( args );
        auditTrailService = getBean( AuditTrailService.class );
        outlierDetectionService = getBean( OutlierDetectionService.class );
        //  statusService = getBean( StatusService.class );
        edService = getBean( ExperimentalDesignService.class );
        securityService = getBean( SecurityService.class );
        svdService = getBean( SVDService.class );
        expressionDataMatrixService = getBean( ExpressionDataMatrixService.class );

        process( super.expressionExperiments );

        return null;
    }

    @Override
    protected void processOptions() {
        super.processOptions();

        if ( hasOption( 'o' ) ) {
            this.viewFile = getOptionValue( 'o' );
            log.info( "GZipped txt output will be written to " + viewFile );
        } else {
            this.viewFile = DEFAULT_VIEW_FILE;
        }

    }

    /**
     * @param bas
     */
    private void process( Collection<BioAssaySet> expressionExperiments ) {
        try {
            generateExperimentMetaData( expressionExperiments );
        } catch ( IOException e ) {
            throw new RuntimeException( e );
        }
    }

}
