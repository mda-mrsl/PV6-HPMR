/****************************************************************
 *
 * $Source: /pv/CvsTree/pv/gen/src/prg/methods/FLASH/initMeth.c,v $
 *
 * Copyright (c) 2002 - 2010
 * Bruker BioSpin MRI GmbH
 * D-76275 Ettlingen, Germany
 *
 * All Rights Reserved
 *
 * $Id: initMeth.c,v 1.60.2.1 2015/03/23 07:46:19 mawi Exp $
 *
 ****************************************************************/

static const char resid[] = "$Id: initMeth.c,v 1.60.2.1 2015/03/23 07:46:19 mawi Exp $(C) 2002 Bruker BioSpin MRI GmbH";

#define DEBUG		0
#define DB_MODULE	0
#define DB_LINE_NR	0

#include "method.h"

/*:=MPB=:=======================================================*
 *
 * Global Function: initMeth
 *
 * Description: This procedure is implicitly called when this
 *	method is selected.
 *
 * Error History:
 *
 * Interface:							*/

void initMeth()
    /*:=MPE=:=======================================================*/
{
    /* KAM change for hpMR to force 2D, allow as few as possible phase encodes */
    int dimRange[2] = { 2,2 };
    int lowMat[3]   = { 8, 8, 8 };        // This method won't allow fewer than 8
    int upMat[3]    = { 2048, 2048, 512 }; // or more than 2048 points...
    int defaultMat[3] = {128, 128, 64};


    DB_MSG(( "-->FLASH:initMeth\n" ));

    /* which version of toolboxes should be active */

    PTB_VersionRequirement( Yes,20100101,"");

    // KAM add Method version
    strcpy(PVM_MethodVersion.feature, "hpMR");
    strcpy(PVM_MethodVersion.versionId, "1.0");
    strcpy(PVM_MethodVersion.institution, "UT MD Anderson");
    strcpy(PVM_MethodVersion.author, "Keith Michel");
    strcpy(PVM_MethodVersion.date, "20200108");
    // end add

    // KAM add for decoupling
    STB_DecOnOffRange();
    STB_InitDecModule();
    // end add

    /*  initialize local and redirected parameters */

    STB_InitFreqDriftCorr();

    NAveragesRange();
    InversionTimeRange();
    if(ParxRelsParHasValue("PVM_NMovieFrames") == No)
    {
        PVM_NMovieFrames = 1;
    }
    if(ParxRelsParHasValue("PVM_NRepetitions") == No)
    {
        PVM_NRepetitions = 1;
    }

    /* encoding group */
    STB_InitEncoding( 2, dimRange, lowMat, upMat, defaultMat);

    /* Initialisation of rf pulse parameters */


    /*
     * 1: pulses declared in parDefinitions.h
     * in this case - ExcPulse1. We initalise it to default name,
     * 5000.0 Hz, and 30 deg flip angle
     * In addition we declare ExcPulse1Enum ExcPulse1Ampl and ExcPulse1Shape
     * to be handled together with pulse struct ExcPulse1. Only if a double
     * array parameter is provided as shape, the pulse supports calculated
     * shapes.
     */

    STB_InitRFPulse("ExcPulse1",      // name of pulse structure
            "ExcPulse1Enum",  // name of pulse list parameter
            "ExcPulse1Ampl",  // name of pulse amplitude parameter
            "ExcPulse1Shape", // name of pulse shape (for calc. pulses)
            0,                // used for excitation
            "Calculated",     // default shape
            3000.0,           // default bandwidth
            30.0);            // default pulse angle

    // check valid range for this specific pulse see parsRelations.c
    ExcPulse1Range();

    /* KAM Initialize fid excitation pulse */
    STB_InitRFPulse("KAM_FidExcPulse",      // name of pulse structure
            "KAM_FidExcPulseEnum",  // name of pulse list parameter
            "KAM_FidExcPulseAmpl",  // name of pulse amplitude parameter
            "KAM_FidExcPulseShape", // name of pulse shape (for calc. pulses)
            0,                // used for excitation
            "Calculated",     // default shape
            3000.0,           // default bandwidth
            30.0);            // default pulse angle

    // check valid range for this specific pulse see parsRelations.c
    KAM_FidExcPulseRange();

    /* KAM Initialize sat excitation pulse */
    STB_InitRFPulse("KAM_SatExcPulse",      // name of pulse structure
            "KAM_SatExcPulseEnum",  // name of pulse list parameter
            "KAM_SatExcPulseAmpl",  // name of pulse amplitude parameter
            "KAM_SatExcPulseShape", // name of pulse shape (for calc. pulses)
            0,                // used for excitation
            "gauss",     // default shape
            3000.0,           // default bandwidth
            30.0);            // default pulse angle

    // check valid range for this specific pulse see parsRelations.c
    KAM_SatExcPulseRange();

    /* Initialisation of nucleus */
    STB_InitNuclei(1);

    /* Initialisation of geometry parameters */
    STB_InitImageGeometry();

    /* KAM no oblique slices allowed */
    PVM_MajSliceOri = Yes;
    ParxRelsMakeNonEditable("PVM_MajSliceOri");

    /* 1: method specific initialisation */

    if(ParxRelsParHasValue("PVM_RepetitionTime") == No)
        PVM_RepetitionTime = 100.0;
    if(ParxRelsParHasValue("PVM_EchoTime") == No)
        PVM_EchoTime = 4.0;
    if(ParxRelsParHasValue("PVM_DeriveGains") == No)
        PVM_DeriveGains = Yes;

    /* kam add for custom readout */
    if(ParxRelsParHasValue("KAM_NAcqPoints") == No)
        KAM_NAcqPoints = 128;


    /* Initialisation of spoilers */
    MRT_InitSpoiler("ReadSpoiler");
    MRT_InitSpoiler("SliceSpoiler");

    if (ParxRelsParHasValue("PVM_MotionSupOnOff") == 0)
        PVM_MotionSupOnOff = Off;


    /* initialize digitizer parameter */

    STB_InitDigPars();
    EffSWhRange();

    /* not a csi experiment */
    PTB_SetSpectrocopyDims( 0, 0 );

    /* Initialisation of modules */
    STB_InitFatSupModule();
    STB_InitMagTransModule();
    STB_InitFovSatModule();
    STB_InitFlowSaturationModule();
    STB_InitTriggerModule();
    STB_InitTaggingModule();
    STB_InitEvolutionModule();
    STB_InitSelIrModule();
    STB_InitBlBloodModule();
    if(!ParxRelsParHasValue("RFSpoiling"))
        RFSpoiling=No;
    if(!ParxRelsParHasValue("AngioMode"))
        AngioMode=No;
    STB_InitDummyScans(1000.0);

    /* initialize mapshim parameter class */
    STB_InitMapShim();
    STB_InitStartupShims();

    /* initialization of method specific reconstruction */
    if(ParxRelsParHasValue("RecoMethMode") == No)
        RecoMethMode=Default;
    if(ParxRelsParHasValue("WeightingMode") == No)
        WeightingMode=positive_mask;
    GaussBroadRange();
    MaskWeightRange();

    /* Visibility of parameters */
    ParxRelsMakeNonEditable("PVM_EchoPosition");
    ParxRelsMakeNonEditable("PVM_MinEchoTime,PVM_AcquisitionTime");
    ParxRelsMakeNonEditable("EncGradDur");

    /*
     * Once all parameters have initial values, the backbone is called
     * to assure they are consistent
     */

    //backbone();

    /* KAM Dummy scan duration keeps resetting to 1sec. Disable them */
    PVM_DummyScans = 0;
    PVM_DummyScansDur = 0.0;
    ParxRelsHideInEditor("PVM_DummyScans, PVM_DummyScansDur");

    /* KAM retrieve gradient timing resolution, max amplitude and slew rate */
    KAM_GradientTime     = CFG_GradientShapeResolution() * 1e3;
    if( ParxRelsParHasValue("KAM_GradientMaxAmp") == No)
        KAM_GradientMaxAmp   = CFG_MaxGradientStrength();
    if( ParxRelsParHasValue("KAM_GradientMaxSlew") == No)
        KAM_GradientMaxSlew  = KAM_GradientMaxAmp / (CFG_GradientRampTime() + CFG_InterGradientWaitTime());

    /* KAM initialize EPI design parameters */
    if( ParxRelsParHasValue("KAM_EpiType") == No)
        KAM_EpiType = Flyback;
    if( ParxRelsParHasValue("KAM_EpiMatrix") == No)
        KAM_EpiMatrix = 16;
    if( ParxRelsParHasValue("KAM_EpiDwellFactor") == No)
        KAM_EpiDwellFactor = 4;
    if( ParxRelsParHasValue("KAM_EpiPartialFourierY") == No)
        KAM_EpiPartialFourierY = 1;

    /* KAM initialize HyperSense trigger */
    if( ParxRelsParHasValue("KAM_HSTrig") == No)
        KAM_HSTrig = No;

    /* KAM initialize IDEAL echo time shifting parameters */
    if( ParxRelsParHasValue("KAM_IdealRepDelay") == No)
        KAM_IdealRepDelay = 3000.0;
    if( ParxRelsParHasValue("KAM_IdealMinDelay") == No)
        KAM_IdealMinDelay = No;
    if( ParxRelsParHasValue("KAM_IdealEchoShift") == No)
        KAM_IdealEchoShift = 0.0;
    if( ParxRelsParHasValue("KAM_FidEcho") == No)
        KAM_FidEcho = No;

    /* KAM initialize variable flip angle parameters */
    if( ParxRelsParHasValue("KAM_VFAOn") == No)
        KAM_VFAOn = No;
    if( ParxRelsParHasValue("KAM_VFAEffectiveOn") == No)
        KAM_VFAEffectiveOn = No;
    if( ParxRelsParHasValue("KAM_VFAEffectiveDegrees") == No)
        KAM_VFAEffectiveDegrees = 90.0;
    if( ParxRelsParHasValue("KAM_VFADegrees") == No)
    {
        PARX_change_dims("KAM_VFADegrees",1);
        KAM_VFADegrees[0] = 30.0;
    }
    if( ParxRelsParHasValue("KAM_VFAWatts") == No)
    {
        PARX_change_dims("KAM_VFAWatts",1);
        KAM_VFAWatts[0] = 0.0;
    }

    /* KAM initialize variable offset frequency parameters */
    if( ParxRelsParHasValue("KAM_SpecExc") == No)
        KAM_SpecExc = No;
    if( ParxRelsParHasValue("KAM_VOFTxHertz") == No)
    {
        PARX_change_dims("KAM_VOFTxHertz",1);
        KAM_VOFTxHertz[0] = 0.0;
    }
    if( ParxRelsParHasValue("KAM_VOFRxHertz") == No)
    {
        PARX_change_dims("KAM_VOFRxHertz",1);
        KAM_VOFRxHertz[0] = 0.0;
    }
    if( ParxRelsParHasValue("KAM_SatFreq") == No)
    {
        PARX_change_dims("KAM_SatFreq",1);
        KAM_SatFreq[0] = 0.0;
    }

    /* KAM add for hpMR - initialize normal excitation and readout */
    if(ParxRelsParHasValue("KAM_ExcMode") == No)
        KAM_ExcMode=Normal_Exc;
    if(ParxRelsParHasValue("KAM_ReadMode") == No)
        KAM_ReadMode=Normal_Read;
    KAM_ExcModeRel();  // These
    KAM_PrepModeRel(); // call
    KAM_ReadModeRel(); // backbone() !!


    DB_MSG(("<--FLASH:initMeth\n"));

}



/****************************************************************/
/*		E N D   O F   F I L E				*/
/****************************************************************/









