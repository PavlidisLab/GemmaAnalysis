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

import gemma.gsec.SecurityService;

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
import java.util.List;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.lang.time.StopWatch;

import ubic.basecode.dataStructure.CountingMap;
import ubic.gemma.analysis.preprocess.OutlierDetails;
import ubic.gemma.analysis.preprocess.OutlierDetectionService;
import ubic.gemma.analysis.preprocess.batcheffects.BatchEffectDetails;
import ubic.gemma.analysis.util.ExperimentalDesignUtils;
import ubic.gemma.apps.ExpressionExperimentManipulatingCLI;
import ubic.gemma.expression.experiment.service.ExperimentalDesignService;
import ubic.gemma.model.common.auditAndSecurity.AuditTrailService;
import ubic.gemma.model.common.auditAndSecurity.Status;
import ubic.gemma.model.common.auditAndSecurity.StatusService;
import ubic.gemma.model.common.description.BibliographicReference;
import ubic.gemma.model.common.quantitationtype.QuantitationType;
import ubic.gemma.model.expression.arrayDesign.ArrayDesign;
import ubic.gemma.model.expression.bioAssay.BioAssay;
import ubic.gemma.model.expression.experiment.BioAssaySet;
import ubic.gemma.model.expression.experiment.ExperimentalDesign;
import ubic.gemma.model.expression.experiment.ExpressionExperiment;
import ubic.gemma.model.expression.experiment.ExpressionExperimentSubSet;
import ubic.gemma.model.expression.experiment.ExpressionExperimentValueObject;
import ubic.gemma.util.FactorValueVector;
import ubic.gemma.util.Settings;

/**
 * Extracts expression experiment meta data such as experimental design, array design, outlier count, and publication
 * into a .txt.gz TSV file. See Bug 3968.
 * 
 * @author paul
 * @version $Id$
 */
public class ExperimentMetaDataExtractorCli extends ExpressionExperimentManipulatingCLI {

    private static final String EXPERIMENT_META_DATA_BASENAME = "ExperimentMetaData";
    private static final String VIEW_FILE_SUFFIX = ".txt.gz";
    public static final String DEFAULT_VIEW_FILE = Settings.getString( "gemma.appdata.home" ) + File.separatorChar
            + "dataFiles" + File.separatorChar + EXPERIMENT_META_DATA_BASENAME + VIEW_FILE_SUFFIX;
    private static final String NA = "";
    private OutlierDetectionService outlierDetectionService;
    private StatusService statusService;
    private ExperimentalDesignService edService;
    private SecurityService securityService;
    private String viewFile = DEFAULT_VIEW_FILE;

    /*
     * (non-Javadoc)
     * 
     * @see ubic.gemma.util.AbstractCLI#doWork(java.lang.String[])
     */
    @Override
    protected Exception doWork( String[] args ) {
        super.processCommandLine( "experiment metadata extract", args );
        auditTrailService = getBean( AuditTrailService.class );
        outlierDetectionService = getBean( OutlierDetectionService.class );
        statusService = getBean( StatusService.class );
        edService = getBean( ExperimentalDesignService.class );
        securityService = getBean( SecurityService.class );

        process( super.expressionExperiments );

        return null;
    }

    @Override
    @SuppressWarnings("static-access")
    protected void buildOptions() {
        super.buildOptions();

        Option expOption = OptionBuilder.hasArg().withArgName( "outfile" ).withDescription( "GZipped output filename" )
                .withLongOpt( "outfile" ).create( 'o' );

        addOption( expOption );

        // to keep troubled experiments
        this.addForceOption();
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

    public void generateExperimentMetaData( Collection<BioAssaySet> expressionExperiments ) throws IOException {

        File file = getOutputFile( this.viewFile );

        try (Writer writer = new OutputStreamWriter( new GZIPOutputStream( new FileOutputStream( file ) ) );) {

            String[] colNames = { "ShortName", "Taxon", "DateUpload", "IsPublic", "Platform", "Channel", "IsExonArray",
                    "QtIsRatio", "QtIsNormalized", "QtScale", "NumProfiles", "NumSamples", "NumFactors",
                    "NumConditions", "NumReplicatesPerCondition", "PossibleOutliers", "CuratedOutlier", "BatchPval",
                    "IsTroubled", "PubTroubled", "PubYear", "PubJournal" };
            // log.debug( StringUtils.join( colNames, "\t" ) + "\n" );
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
                    Status pubStatus = primaryPublication != null ? statusService.getStatus( primaryPublication )
                            : null;
                    Collection<ArrayDesign> arrayDesignsUsed = eeService.getArrayDesignsUsed( ee );
                    if ( arrayDesignsUsed.size() > 1 ) {
                        log.warn( "Multiple array designs found. Only the first array design will be reported." );
                    }
                    ArrayDesign arrayDesign = arrayDesignsUsed.iterator().next();

                    QuantitationType qt = null;
                    for ( QuantitationType q : ee.getQuantitationTypes() ) {
                        if ( q.getIsPreferred() ) {
                            qt = q;
                            break;
                        }
                    }

                    int manualOutlierCount = 0;
                    for ( BioAssay ba : ee.getBioAssays() ) {
                        if ( ba.getIsOutlier() ) {
                            manualOutlierCount++;
                        }
                    }

                    ExperimentalDesign experimentalDesign = edService.load( vo.getExperimentalDesign() );

                    BatchEffectDetails batchEffect = eeService.getBatchEffect( ee );
                    if ( batchEffect == null ) {
                        log.warn( "Null batch effect info" );
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

                    Collection<String> samplesPerConditionCount = new ArrayList<>();
                    // warning: this removes batchEffect factor!
                    CountingMap<FactorValueVector> assayCount = ExperimentalDesignUtils.getDesignMatrix( ee, true );
                    List<FactorValueVector> keys = assayCount.sortedKeyList( true );
                    for ( FactorValueVector key : keys ) {
                        samplesPerConditionCount.add( Integer.toString( assayCount.get( key ) ) );
                    }

                    String val[] = {
                            vo.getShortName(),
                            vo.getTaxon(),
                            DateFormat.getDateInstance( DateFormat.MEDIUM ).format( vo.getDateCreated() ),
                            vo != null ? Boolean.toString( vo.getIsPublic() ) : NA,
                            arrayDesign.getShortName(),
                            arrayDesign.getTechnologyType().getValue(), // ONE-COLOR, TWO-COLOR, NONE (RNA-seq
                                                                        // GSE37646), DUAL-MODE
                                                                        // (one or two color)
                            Boolean.toString( arrayDesign.getName().toLowerCase().contains( "exon" ) ), // exon GSE28383
                            qt != null ? Boolean.toString( qt.getIsRatio() ) : NA,
                            qt != null ? Boolean.toString( qt.getIsNormalized() ) : NA,
                            qt != null ? qt.getScale().getValue() : NA,
                            Integer.toString( vo.getProcessedExpressionVectorCount() ), // NumProfiles
                            Integer.toString( vo.getBioAssayCount() ), // NumSamples
                            batchEffect != null ? Integer
                                    .toString( experimentalDesign.getExperimentalFactors().size() - 1 ) : Integer
                                    .toString( experimentalDesign.getExperimentalFactors().size() ), // NumFactors
                            Integer.toString( assayCount.size() ), // NumConditions
                            StringUtils.join( samplesPerConditionCount, "," ),
                            possibleOutliers != null ? Integer.toString( possibleOutliers.size() ) : NA,
                            Integer.toString( manualOutlierCount ),
                            batchEffect != null ? String.format( "%.2g", batchEffect.getPvalue() ) : NA,
                            Boolean.toString( vo.getTroubled() ),
                            pubStatus != null ? Boolean.toString( pubStatus.getTroubled() ) : NA,
                            primaryPublication != null ? DateFormat.getDateInstance( DateFormat.MEDIUM ).format(
                                    primaryPublication.getPublicationDate() ) : NA,
                            primaryPublication != null ? primaryPublication.getPublication() : NA };

                    // log.debug( StringUtils.join( val, "\t" ) + "\n" );
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

}
