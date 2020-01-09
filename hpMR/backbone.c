
/* ***************************************************************
 *
 * $Source: /pv/CvsTree/pv/gen/src/prg/methods/FLASH/backbone.c,v $
 *
 * Copyright (c) 2009
 * Bruker BioSpin MRI GmbH
 * D-76275 Ettlingen, Germany
 *
 * All Rights Reserved
 *
 *
 * $Id: backbone.c,v 1.50.2.1 2015/03/23 07:43:19 mawi Exp $
 *
 * ***************************************************************/


static const char resid[] = "$Id $ (c) 2007Bruker BioSpin MRI GmbH";

#define DEBUG		0
#define DB_MODULE	0
#define DB_LINE_NR	0

#include "method.h"

/* global variables */
double EchoDelay, EffPulseDur;

void backbone(void)
{
    double minFov[3];
    double minSliceThick;

    DB_MSG(("-->backbone"));

    /* update nuclei parameter group                            */

    STB_UpdateNuclei(No);

    /* KAM enforce hpMR design parameters */
    EnforceHpmrDesignParams();
    KAM_IdealEchoShift = MAX_OF( KAM_IdealEchoShift, 0.0);

    /* update encoding parameter group                          */

    STB_UpdateEncoding(
            NULL,
            SEG_SEQUENTIAL,
            Yes,
            Yes,
            Yes,
            &PVM_EchoPosition);

    /* update parameters controlling the data sampling period   */

    STB_UpdateDigPars(&PVM_EffSWh,
            KAM_NAcqPoints, //PVM_EncMatrix[0],
            &PVM_AntiAlias[0],
            &PVM_AcquisitionTime);

    /* KAM ensure consistent inner loop sizes */
    if( KAM_ReadMode == Normal_Read) {
        KAM_NEchoes = PVM_EncMatrix[1] / KAM_NSlices;
        PVM_EncMatrix[1] = KAM_NSlices * KAM_NEchoes;
        PVM_Matrix[1]    = KAM_NSlices * KAM_NEchoes;
    } else {
        PVM_EncMatrix[1] = KAM_NSlices * KAM_NEchoes;
        PVM_Matrix[1]    = KAM_NSlices * KAM_NEchoes;
        if( KAM_FidEcho == Yes) {
            PVM_EncMatrix[1]++;
            PVM_Matrix[1]++;
        }
        if( KAM_SatEcho == Yes) {
            PVM_EncMatrix[1]++;
            PVM_Matrix[1]++;
        }
    }

    /* KAM handle prep grad logic */
    PARX_change_dims("KAM_PrepGradOnOuter",PVM_NRepetitions);
    for(int i=0; i<PVM_NRepetitions; i++)
    {
        KAM_PrepGradOnOuter[i] = MAX_OF( MIN_OF( KAM_PrepGradOnOuter[i], 1), -1);
    }
    PARX_change_dims("KAM_PrepGradOnInner",KAM_NEchoes);
    for(int i=0; i<KAM_NEchoes; i++)
    {
        KAM_PrepGradOnInner[i] = MAX_OF( MIN_OF( KAM_PrepGradOnInner[i], 1), -1);
    }

    /* KAM handle decoupling module */
    PVM_DecDuration = MAX_OF( MIN_OF(PVM_DecDuration, 100.0), 0.0);
    if( PVM_DecOnOff == On) {
        if( PVM_NumberOfNuclei != 2) {
            STB_InitNuclei(2);
            if( PVM_NumberOfNuclei < 2) {
                /* system config does not support 2nd RF channel */
                PVM_DecOnOff = Off;
            }
        }
    }
    else {
        if( PVM_NumberOfNuclei != 1) {
            STB_InitNuclei(1);
        }
    }
    STB_UpdateDecModule(PVM_Nucleus2, PVM_AcquisitionTime);

    /*KAM Hack for phantom tests*/
    /*
       if( KAM_SpecExc == Yes) {
       PARX_change_dims("KAM_VFADegrees", KAM_NEchoes);
       for(int i=0; i<KAM_NEchoes; i++) {
       KAM_VFADegrees[i] = 90.0;
       }
       }
       */

    /* KAM handle variable flip angles */
    HandleVFA();

    /* KAM handle variable offset frequencies */
    if( KAM_SpecExc == Yes) {
        PARX_change_dims("KAM_VOFTxHertz",KAM_NEchoes);
        for(int i=0; i<KAM_NEchoes; i++) {
            //KAM_VOFTxHertz[i] = -2000.0 + i*4000.0/(KAM_NEchoes-1); // Hack for spsp phantom tests
            //KAM_VOFTxHertz[i] = -200.0 + i*400.0/(KAM_NEchoes-1); // Hack for spsp phantom tests
            KAM_VOFTxHertz[i] = MAX_OF( MIN_OF( KAM_VOFTxHertz[i], 10000.0), -10000.0);
        }
        PARX_change_dims("KAM_VOFRxHertz",KAM_NEchoes);
        for(int i=0; i<KAM_NEchoes; i++) {
            if( KAM_SpecTxRxSeparate == Yes)
                KAM_VOFRxHertz[i] = MAX_OF( MIN_OF( KAM_VOFRxHertz[i], 10000.0), -10000.0);
            else
                KAM_VOFRxHertz[i] = KAM_VOFTxHertz[i];
        }
        PARX_change_dims("KAM_SatFreq",1);
        KAM_SatFreq[0] = MAX_OF( MIN_OF( KAM_SatFreq[0], 10000.0), -10000.0);
    }

    /* update excitation pulse                                  */
    // KAM- this is done in VFA handling, where indicated
    //UpdateRFPulses();

    /* general features of the method */

    PVM_NEchoImages = 1;

    /* KAM no preparation modules allowed */
    PVM_EvolutionOnOff = Off;
    PVM_TaggingOnOff = Off;
    PVM_SelIrOnOff = Off;
    PVM_BlBloodOnOff = Off;
    PVM_FatSupOnOff = Off;
    PVM_MagTransOnOff = Off;
    PVM_FovSatOnOff = Off;
    PVM_InFlowSatOnOff = Off;
    PVM_MotionSupOnOff = Off;

    /* coexistance of modules: */
    if(PVM_SelIrOnOff == On)
    {
        PVM_BlBloodOnOff = Off;
        PVM_NMovieFrames = 1;
    }

    /* set limits for read, phase and slice gradients            */

    ControlGradientLimits(PVM_MajSliceOri);


    /* calculate minima for FOV and slice thickness             */

    UpdateGeometryMinima(minFov,
            &minSliceThick);


    /* update geometry parameters                               */

    int dim=PTB_GetSpatDim();

    // only one package if black-blood module on
    int maxPackages = PVM_BlBloodOnOff == On? 1:0;

    // only one slice per package if 3D
    int maxPerPackage = dim>2? 1:0;

    STB_UpdateImageGeometry(dim,
            PVM_Matrix,
            minFov,
            0, // total slices (no restr)
            maxPackages,
            maxPerPackage,
            minSliceThick,
            1.0); // sliceFovRatio in 3D

    /* update slice spoiler */
    double mindurSlice = 1.5*CFG_GradientRiseTime();
    double spoilerThick = dim>2 ? PVM_SpatResol[2]*PVM_EncZf[2] : PVM_SliceThick;
    MRT_UpdateSpoiler("SliceSpoiler",2.0,ExcSliceGradLim,mindurSlice,PVM_GradCalConst,spoilerThick);

    /* handling of modules */
    STB_UpdateFatSupModule(PVM_Nucleus1, PVM_DeriveGains, spoilerThick);
    STB_UpdateMagTransModule(PVM_DeriveGains);
    STB_UpdateFovSatModule(PVM_Nucleus1, PVM_DeriveGains);
    STB_UpdateFlowSaturationModule(PVM_Nucleus1,PVM_DeriveGains);
    STB_UpdateTriggerModule();
    STB_UpdateTaggingModule(PVM_Fov,PVM_Matrix,PVM_SpatResol[0]*PVM_EncZf[0],PVM_DeriveGains);
    STB_UpdateEvolutionModule();

    /* Black Blood Module   */
    {
        double slabthick = PVM_SPackArrSliceDistance[0] * (PVM_SPackArrNSlices[0]-1)
            + PVM_SliceThick;

        double fixedTime = SliceSpoiler.dur + CFG_AmplifierEnable()
            + ExcPulse1.Length/2.0 + PVM_TaggingModuleTime;

        STB_UpdateBlBloodModule(&slabthick,PVM_SPackArrSliceOffset,1,fixedTime,PVM_SliceThick,PVM_DeriveGains);
    }

    UpdateMovie();
    UpdateRFSpoiling();

    /* Calculate read and slice gradients */
    UpdateReadSliceGradients();

    /* KAM Minimize time course delay */
    if( KAM_IdealMinDelay == Yes)
        KAM_IdealRepDelay = 0.0;

    /* Sequence elements, TE, TR, total time: */
    UpdateSequenceTiming();

    /* Dummy Scans */
    int averages;
    if(PVM_MotionSupOnOff==On || AngioMode ==Yes)
        averages = 1;
    else
        averages = PVM_NAverages;

    STB_UpdateDummyScans(PVM_RepetitionTime,averages);

    /* calculate dephasing and phase-encoding gradients         */
    UpdatePhaseGradients();

    /* calculate frequency offsets                              */
    UpdateFrequencyOffsets();

    /* update mapshim parameter class */
    STB_UpdateMapShim(PVM_Nucleus1,"PVM_SliceGeoObj");
    STB_UpdateStartupShims();

    /* Set ACQP parameters */
    SetBaseLevelParam();

    /* Set RECO parameters */
    SetRecoParam();

    DB_MSG(("<--backbone"));
    return;
}

/*-------------------------------------------------------
 * local utility routines to simplify the backbone
 *------------------------------------------------------*/


void UpdateSequenceTiming()
    /* -------------------------------------------------------
       Adapt sequence elements to the current geometry
       (in this method, only EncGradDur is concerned),
       update TE and TR.
       ReadGrad and ExSliceGrad must be already set.
       ------------------------------------------------------*/
{
    double minEnc2d, minEnc3d, minRephSlice, minDephRead,
           extensionPossible, extensionAllowed, extension;
    int nSlices;

    /* Part of the exctiation pulse to be refocused: */
    // KAM - Not sure if echo time should include 1/2 of
    //    spec-spat pulse or not.... Leaving this for now.
    EffPulseDur = ExcPulse1.Length * (ExcPulse1.Rpfac/100);

    /* KAM spectral-spatial pulses don't require refocusing */
    if( KAM_ExcMode == Normal_Exc)
    {
        minRephSlice = MRT_DephaseTime(EffPulseDur,
                CFG_GradientRiseTime(),
                ExcSliceGrad,
                ExcSliceRephGradLim);
    }
    else
        minRephSlice = 0.0;

    /* KAM spiral trajectories start at k-space center */
    if( KAM_ReadMode == Normal_Read)
    {
        /* Part of the echo to be refocused */
        EchoDelay = PVM_AcquisitionTime * PVM_EchoPosition / 100;

        /* Minimum durations of all phase-gradient pulses */
        minEnc2d     = PTB_GetSpatDim()>1?
            MRT_EncodingTime(PVM_SpatResol[1],
                    PVM_GradCalConst*Phase2DGradLim/100)
            : 0;
        minEnc3d     = PTB_GetSpatDim()>2?
            MRT_EncodingTime(PVM_SpatResol[2],
                    PVM_GradCalConst*Phase3DGradLim/100)
            : 0;

        /* Minimum duration of read dephase */
        minDephRead  = MRT_DephaseTime(EchoDelay,
                CFG_GradientRiseTime(),
                ReadGrad,
                ReadDephGradLim);

        /* In this sequence all phase-gradient pulses are
           simultaneous with duration EncGradDur. We set it first
           to the common (longest) minimum: */
        EncGradDur = MAX_OF( MAX_OF(minEnc2d ,  minEnc3d),
                MAX_OF(minRephSlice ,minDephRead ) );
    }
    else
    {
        /* Spirals start at k-space center, other trajs use KAM_ReadEchoDelay */
        EchoDelay = minEnc2d = minEnc3d = minDephRead = 0.0;
        /* Normal excitations still require slice rephase */
        EncGradDur = minRephSlice;
    }

    /* EncGradDur should also contain one ramp, thus: */
    EncGradDur = RewGradDur = MAX_OF(EncGradDur, CFG_GradientRiseTime());

    DB_MSG(("minimum EncGradDur = %f", EncGradDur));

    /* Update TE with the mimimum EncGradDur */

    UpdateEchoTime();

    /* KAM extension not needed for arbitrary readouts */
    if( KAM_ReadMode == Normal_Read)
    {
        /* If there is some freedom, make EncGradDur longer
           (to avoid unnecessarily strong phase gradients),
           but not longer than allowed by gradient resolution */
        extensionPossible = PVM_EchoTime - PVM_MinEchoTime;
        extensionAllowed  = PTB_GetSpatDim()>1? MRT_MaxEncodingTime(PVM_Fov[1], PVM_GradCalConst)-EncGradDur :20;

        DB_MSG(("ext possible = %f, allowed = %f", extensionPossible, extensionAllowed));

        extension = MIN_OF(extensionPossible, extensionAllowed);
        extension = MAX_OF(extension, 0);

        DB_MSG(("extension = %f",extension));
    }
    else
        extension = 0.0;
    EncGradDur += extension;
    TeFillDelay = PVM_EchoTime - PVM_MinEchoTime - extension;



    /* Other sequence elements, not involved in TE, e.g. spoilers  */
    double mindurRead =  MAX_OF(RewGradDur+CFG_GradientRiseTime(),PVM_DigEndDelOpt);
    MRT_UpdateSpoiler("ReadSpoiler",2.0,Phase2DGradLim,mindurRead,PVM_GradCalConst,PVM_SpatResol[0]*PVM_EncZf[0]);

    /* Update modules dependent on TE, here: inversion recovery: */
    SliceSegDur = minLoopRepetitionTime(); /* depends on TE */
    nSlices = GTB_NumberOfSlices( PVM_NSPacks, PVM_SPackArrNSlices);
    double fixedTime = SliceSpoiler.dur + CFG_AmplifierEnable()
        + ExcPulse1.Length/2.0 + PVM_TaggingModuleTime;
    STB_UpdateSelIrModule(PVM_SliceThick,PVM_SliceOffset,nSlices,&SliceSegDur,0,fixedTime,PVM_DeriveGains);
    if(PVM_SelIrOnOff==On)
        ParxRelsCopyPar("PVM_InversionTime","PVM_SelIrInvTime");
    if(PVM_BlBloodOnOff==On)
        ParxRelsCopyPar("PVM_InversionTime","PVM_BlBloodInvTime");
    if(PVM_BlBloodOnOff==Off&&PVM_SelIrOnOff==Off)
    {
        PVM_InversionTime = 0.0;
        ParxRelsHideInEditor("PVM_InversionTime");
    }

    /* Find min TR and update TR */
    UpdateRepetitionTime();

    /* delay after the slice loop, used only in IR mode to cotrol TR: */
    SliceSegEndDelay = PVM_SelIrOnOff==On? PVM_RepetitionTime - PVM_MinRepetitionTime : 0;

    /* Calculate total experiment time */

    UpdateTotalTime();
}




void UpdateRepetitionTime(void)
{
    int nSlices;

    DB_MSG(("-->UpdateRepetitionTime"));


    nSlices = GTB_NumberOfSlices( PVM_NSPacks, PVM_SPackArrNSlices );

    if(AngioMode == No)
    {
        if(PVM_SelIrOnOff == On)
        {
            PVM_MinRepetitionTime =
                PVM_SelIrModuleTime +
                PVM_TaggingModuleTime +
                nSlices * SliceSegDur;
        }
        else
        {
            /* min TR in a movie: */
            PVM_MinRepetitionTime = nSlices * minLoopRepetitionTime();

            /* if there is no movie, TR also includes some modules: */
            if(PVM_NMovieFrames == 1)
            {
                PVM_MinRepetitionTime +=
                    PVM_BlBloodModuleTime +
                    PVM_SelIrModuleTime +
                    PVM_TaggingModuleTime;
            }
        }
    }
    else
        PVM_MinRepetitionTime = minLoopRepetitionTime();

    PVM_RepetitionTime=MAX_OF(PVM_MinRepetitionTime,PVM_RepetitionTime);

    /*** adapt RewGradDur to desired TR ***/
    double riseT = CFG_GradientRiseTime();
    double extensionAllowed = EncGradDur - RewGradDur;
    double extensionPossible = ReadSpoiler.dur-riseT - RewGradDur;
    double extension = MAX_OF(MIN_OF(extensionPossible, extensionAllowed),0);
    /* 1.) if possible, extend RewGradDur to duration of read spoiler */
    if(extension > 0)
    {
        RewGradDur += extension;
        extensionAllowed = EncGradDur - RewGradDur;
    }
    int slices = AngioMode==Yes? 1:nSlices;
    extensionPossible = (PVM_RepetitionTime - PVM_MinRepetitionTime)/slices;
    extension = MAX_OF(MIN_OF(extensionPossible, extensionAllowed),0);
    /* 2.) if possible, extend RewGradDur to EncGradDur  */
    if(extension > 0)
    {
        RewGradDur += extension;
        MRT_UpdateSpoiler("ReadSpoiler",2.0,Phase2DGradLim,RewGradDur+riseT,PVM_GradCalConst,PVM_SpatResol[0]*PVM_EncZf[0]);
        PVM_MinRepetitionTime += extension*slices;
    }

    DB_MSG(("<--UpdateRepetitionTime"));
    return;
}



/* calculates PVM_ScanTimeStr and TimeForMovieFrames */

void UpdateTotalTime(void)
{
    int dim = PTB_GetSpatDim();
    double TotalTime=0;

    if( dim >1 )
    {
        if(PVM_NMovieFrames > 1)
        {
            /* TR is one movie frame, without prep modules: */
            TimeForMovieFrames = PVM_RepetitionTime * PVM_NMovieFrames + PVM_TaggingModuleTime;
            TotalTime = (PVM_BlBloodModuleTime + TimeForMovieFrames)
                * PVM_EncMatrix[1] * PVM_NAverages;
        }
        else
        {
            if(AngioMode==No)
            {
                /* TR includes prep modules */
                TimeForMovieFrames = PVM_RepetitionTime - PVM_BlBloodModuleTime - PVM_TaggingModuleTime;
                TotalTime = PVM_RepetitionTime * PVM_EncMatrix[1] * PVM_NAverages;
            }
            else
            {
                int nSlices = GTB_NumberOfSlices( PVM_NSPacks, PVM_SPackArrNSlices );
                TotalTime = (PVM_RepetitionTime * PVM_EncMatrix[1] +
                        PVM_BlBloodModuleTime+PVM_SelIrModuleTime+PVM_TaggingModuleTime)*nSlices*PVM_NAverages;
            }
        }

    }

    if( dim >2 )
        TotalTime = TotalTime * PVM_EncMatrix[2];

    /* KAM add for IDEAL time course delay */
    KAM_IdealRepDelay = MAX_OF( KAM_IdealRepDelay, TotalTime);
    KAM_IdealRepFill = KAM_IdealRepDelay - TotalTime;
    TotalTime = KAM_IdealRepDelay;

    /* time for one repetition */
    OneRepTime = TotalTime/1000.0;

    TotalTime = TotalTime * PVM_EvolutionCycles + PVM_EvolutionModuleTime;
    TotalTime = TotalTime * PVM_NRepetitions;

    PVM_ScanTime = TotalTime;
    UT_ScanTimeStr(PVM_ScanTimeStr,TotalTime);
    ParxRelsShowInEditor("PVM_ScanTimeStr");
    ParxRelsMakeNonEditable("PVM_ScanTimeStr");

}

void UpdateGeometryMinima( double *minFov,
        double *minSliceThick)
{
    int dim;


    DB_MSG(("-->UpdateGeometryMinima"));


    dim=PTB_GetSpatDim();

    // KAM only constrain minFov for normal readouts
    if( KAM_ReadMode == Normal_Read) {
        minFov[0]     = PVM_EffSWh /
            (1e-2*ReadDephGradLim * PVM_GradCalConst);
    } else {
        minFov[0] = 1.0;
    }
    *minSliceThick = ExcPulse1.Bandwidth /
        (1e-2*ExcSliceGradLim * PVM_GradCalConst);

    if(dim >= 2)
    {
        minFov[1] = minFov[0]/8;
    }

    if(dim >= 3)
    {

        minFov[2] = *minSliceThick;
    }

    DB_MSG(("<--UpdateGeometryMinima"));
}


void UpdateReadSliceGradients(void)
{
    DB_MSG(("-->UpdateReadSliceGradients"));

    /* KAM this only needs to be calculated for normal readouts */
    if( KAM_ReadMode == Normal_Read) {
        ReadGrad = MRT_ReadGrad(PVM_EffSWh,
                PVM_Fov[0],
                PVM_GradCalConst);
    } else
        ReadGrad = 0.0;

    ExcSliceGrad = MRT_SliceGrad(ExcPulse1.Bandwidth,
            PVM_SliceThick,
            PVM_GradCalConst);

    DB_MSG(("<--UpdateReadSliceGradients"));
}

void UpdatePhaseGradients()
{

    DB_MSG(("-->UpdatePhaseGradients"));

    /* Calculation of phase-encoding,
       dephasing and rephasing gradients.

       (ReadGrad, ExcSliceGrad, EchoDelay, EffPulseDur,
       and EncGradDur must be calculated before)       */

    double rise = CFG_GradientRiseTime();

    /* KAM this only needs to be calculated for normal excitations */
    if( KAM_ExcMode == Normal_Exc)
        ExcSliceRephGrad = MRT_DephaseGrad(EncGradDur, EffPulseDur, rise, ExcSliceGrad);
    else
        ExcSliceRephGrad = 0.0;

    /* KAM these only need to be calculated for normal readouts */
    if( KAM_ReadMode == Normal_Read)
    {
        ReadDephGrad =     MRT_DephaseGrad(EncGradDur, EchoDelay,   rise, ReadGrad);

        Phase2DGrad = PTB_GetSpatDim() > 1 ?
            MRT_PhaseGrad(EncGradDur, PVM_Matrix[1], PVM_Fov[1], PVM_GradCalConst) : 0.0;

        Phase3DGrad = PTB_GetSpatDim() == 3 ?
            MRT_PhaseGrad(EncGradDur, PVM_Matrix[2], PVM_Fov[2], PVM_GradCalConst) : 0.0;

        Rew2DGrad = PTB_GetSpatDim() > 1 ?
            MRT_PhaseGrad(RewGradDur, PVM_Matrix[1], PVM_Fov[1], PVM_GradCalConst) : 0.0;

        Rew3DGrad = PTB_GetSpatDim() == 3 ?
            MRT_PhaseGrad(RewGradDur, PVM_Matrix[2], PVM_Fov[2], PVM_GradCalConst) : 0.0;
    }
    else
        ReadDephGrad = Phase2DGrad = Phase3DGrad = Rew2DGrad = Rew3DGrad = 0.0;

    DB_MSG(("<--UpdatePhaseGradients"));
    return;
}

void UpdateFrequencyOffsets( void )
{
    int spatDim;
    int i,nslices;

    spatDim = PTB_GetSpatDim();
    nslices = GTB_NumberOfSlices(PVM_NSPacks,PVM_SPackArrNSlices);

    MRT_FrequencyOffsetList(nslices,
            PVM_EffReadOffset,
            ReadGrad,
            PVM_GradCalConst,
            PVM_ReadOffsetHz );

    MRT_FrequencyOffsetList(nslices,
            PVM_EffSliceOffset,
            ExcSliceGrad,
            PVM_GradCalConst,
            PVM_SliceOffsetHz );

    if(spatDim == 3)
    {
        for(i=0;i<nslices;i++)
            PVM_EffPhase2Offset[i] = -PVM_EffSliceOffset[i];
    }


}


/*--------------------------------------------------------
 * Routine to update RF pulse parameters
 *-------------------------------------------------------*/

void UpdateRFPulses(void)
{

    /* Updates all parameters that belong to ExcPulse1 pulse structure
       (as initialized by STB_InitRFPulse see initMeth.c)
       */

    STB_UpdateRFPulse("ExcPulse1",1,PVM_DeriveGains,Conventional);

    // KAM add for fid excitation
    STB_UpdateRFPulse("KAM_FidExcPulse",1,PVM_DeriveGains,Conventional);
    // KAM add for sat excitation
    STB_UpdateRFPulse("KAM_SatExcPulse",1,PVM_DeriveGains,Conventional);

    if(PVM_DeriveGains==Yes)
    {
        ParxRelsHideInEditor("ExcPulse1Ampl");
        // KAM add for fid excitation
        ParxRelsHideInEditor("KAM_FidExcPulseAmpl");
        // KAM add for sat excitation
        ParxRelsHideInEditor("KAM_SatExcPulseAmpl");
    }
    else
    {
        ParxRelsShowInEditor("ExcPulse1Ampl");
        // KAM add for fid excitation
        ParxRelsShowInEditor("KAM_FidExcPulseAmpl");
        // KAM add for sat excitation
        ParxRelsShowInEditor("KAM_SatExcPulseAmpl");
    }

    ParxRelsShowInFile("ExcPulse1Ampl");
    // KAM add for fid excitation
    ParxRelsShowInFile("KAM_FidExcPulseAmpl");
    // KAM add for sat excitation
    ParxRelsShowInFile("KAM_SatExcPulseAmpl");

    DB_MSG(("<--UpdateRFPulses"));

    return;
}

/*--------------------------------------------------------
 * Routine to control the visibility of parameters
 *-------------------------------------------------------*/


void ControlGradientLimits(YesNo NotOblique)
{
    DB_MSG(("-->ControlGradientLimits: Obliqueness forbidden: %s",NotOblique==Yes ? "Yes":"No"));


    if(NotOblique==Yes)
    {
        ReadDephGradLim     =
            Phase2DGradLim      =
            Phase3DGradLim      =
            ExcSliceRephGradLim = 100.0;
    }
    else
    {
        /* Gradient limits in % of max. Value 57 (1/sqrt(3))
           is needed for arbitrary oblique slices. */
        ReadDephGradLim     =
            Phase2DGradLim      =
            Phase3DGradLim      =
            ExcSliceRephGradLim = 57.0;
    }

    ReadGradLim        = 100.0;
    ExcSliceGradLim    = 100.0;

    DB_MSG(("-->ControlGradientLimits"));
}

void UpdateMovie(void)
{
    if(AngioMode==Yes)
        PVM_NMovieFrames=1;

    if(PVM_NMovieFrames == 1)
    {
        ParxRelsHideInEditor("TimeForMovieFrames");
    }
    else
    {
        ParxRelsShowInEditor("TimeForMovieFrames");
    }
}




/* Calculates PVM_MinEchoTime and restrict PVM_EchoTIme.
   EffPulseDur, EncGradDur must be set before */
void UpdateEchoTime( void )
{

    double riseTime = CFG_GradientRiseTime(),
           rampTime = CFG_GradientRampTime()+CFG_InterGradientWaitTime();

    DB_MSG(("-->UpdateEchoTime\n"));

    PVM_MinEchoTime =
        KAM_ReadEchoDelay +      //KAM - delay between start of read grads and k-space center
        KAM_ReadAcqDelay/1000 +  //KAM - delay between start of read grads and acq start
        EffPulseDur    +
        rampTime       + //falling ramp of slice grad
        EncGradDur     + //enc. time (ramp+plateau
        rampTime       + //min te-filling (end ramp of encoding)
        riseTime       + //read-on ramp          +
        EchoDelay;

    /*KAM add preparation gradient duration*/
    if( KAM_PrepMode != No_Prep)
        PVM_MinEchoTime += KAM_PrepLength;

    PVM_EchoTime = MAX_OF(PVM_EchoTime, PVM_MinEchoTime);

    DB_MSG(("<--echoTimeRels\n"));
}



double minLoopRepetitionTime(void)
    /* ---------------------------------------------------------
       this function returns the minimum duration of the innermost
       pulse program loop
       ----------------------------------------------------------*/
{
    double minTr, minD0,
           riseTime = CFG_GradientRiseTime(),
           igwt = CFG_InterGradientWaitTime();

    minD0 = 0.01  /* ADC_END */ + igwt;

    minTr =
        0.02                           + /*reload B0*/
        PVM_FatSupModuleTime           +
        PVM_MagTransModuleTime         +
        PVM_FovSatModuleTime           +
        PVM_InFlowSatModuleTime        +
        SliceSpoiler.dur               +
        riseTime                       +
        CFG_AmplifierEnable()          +
        ExcPulse1.Length - EffPulseDur +
        PVM_EchoTime - EchoDelay       +
        PVM_AcquisitionTime            +
        ReadSpoiler.dur                +
        riseTime                       +
        0.02                           +
        minD0;

    /*KAM add IDEAL echo time shifts*/
    minTr += KAM_IdealEchoShift * (KAM_NEchoes-1);

    /*KAM add RF excitation delay*/
    minTr += KAM_ExcRfDelay/1000;

    /*KAM subtract custom readout delays*/
    minTr -= KAM_ReadAcqDelay/1000 + KAM_ReadEchoDelay;

    /*KAM handle custom readout grad duration*/
    double readDur = KAM_ReadLength > PVM_AcquisitionTime+KAM_ReadAcqDelay/1e3 ?
        KAM_ReadLength : PVM_AcquisitionTime+KAM_ReadAcqDelay/1e3;
    minTr += readDur;

    return minTr;
}

void UpdateRFSpoiling(void)
{
    DB_MSG(("-->UpdateRFSpoiling"));

    if(RFSpoiling==No)
    {
        PARX_change_dims("RFPhaseList",1);
        RFPhaseList[0] = 0;
    }
    else
    {
        int max = PVM_EncMatrix[1]+PVM_DummyScans;
        int size = MAX_OF(256,max);
        PARX_change_dims("RFPhaseList",size);
        MRT_RFSpoilPhaseList(117,size,RFPhaseList,Yes);
    }

    DB_MSG(("<--UpdateRFSpoiling"));
}

/*-------------------------------------------------------
 * KAM- additional local routines to simplify the backbone
 *------------------------------------------------------*/

/* This function handles parameters related to variable RF powers */
void HandleVFA( void )
{
    DB_MSG(("-->HandleVFA"));

    int nVFA = KAM_NEchoes;
    char errString[512];
    double maxRfPower;
    if( KAM_FidEcho == Yes) {
        STB_UpdateRFPulse("KAM_FidExcPulse",1,PVM_DeriveGains,Conventional);
        STB_UpdateRFPulse("KAM_SatExcPulse",1,PVM_DeriveGains,Conventional);
    }
    if( KAM_VFAOn == Yes) // Allow user to vary RF powers for ExcPulse1
    {
        PARX_change_dims("KAM_VFADegrees",nVFA);
        PARX_change_dims("KAM_VFAWatts",nVFA);
        ParxRelsShowInEditor("KAM_VFADegrees, KAM_VFAWatts, KAM_VFAEffectiveOn");

        if( KAM_ExcMode == Normal_Exc)
            STB_UpdateRFPulse("ExcPulse1",1,PVM_DeriveGains,Conventional);
        maxRfPower = CFG_MaxRFPower(1, PVM_Nucleus1, ExcPulse1.Pint, ExcPulse1.Length);

        // If RefPow is set, try to set excitation angle to 90 deg
        if( ParxRelsParHasValue("PVM_RefPowCh1") == Yes) {
            // UT_ReportError("Reference power must be set to use variable excitation powers");
            ExcPulse1.Pow = MRT_RFPulsePower(ExcPulse1.Sint, 90.0, ExcPulse1.Length, PVM_RefPowCh1, errString);
            if( *errString)
                UT_ReportError(errString);
        }

        // Power limit check
        ExcPulse1.Pow = MIN_OF(ExcPulse1.Pow, maxRfPower);

        // If RefPow is set and max power is exceeded, use max excitation angle
        if( ParxRelsParHasValue("PVM_RefPowCh1") == Yes) {
            ExcPulse1.Flipangle = MRT_RFPulseAngle(ExcPulse1.Sint, ExcPulse1.Pow, ExcPulse1.Length, PVM_RefPowCh1, errString);
            if( *errString)
                UT_ReportError(errString);
        }
        double maxVFADeg = ExcPulse1.Flipangle;

        if( PVM_DeriveGains == Yes) // User specifies flipangles in degrees
        {
            ParxRelsMakeEditable("KAM_VFAEffectiveOn");
            if( KAM_VFAEffectiveOn == No) // Each flipangle input individually
            {
                ParxRelsMakeEditable("KAM_VFADegrees");
                ParxRelsHideInEditor("KAM_VFAEffectiveDegrees");
            } else {    // Calculate individual flipangles for specified effective flipangle
                ParxRelsMakeNonEditable("KAM_VFADegrees");
                ParxRelsShowInEditor("KAM_VFAEffectiveDegrees");
                KAM_VFAEffectiveDegrees = MAX_OF( MIN_OF( KAM_VFAEffectiveDegrees, 90.0), 0.01);
                double vfaMx = sin(UT_Radians(KAM_VFAEffectiveDegrees)) / sqrt(nVFA);
                double Mz = 1;
                double Mx = 0;
                for( int i=0; i<nVFA; i++)
                {
                    KAM_VFADegrees[i] = UT_Degrees(asin(vfaMx / Mz));
                    Mx = Mz * sin(UT_Radians(KAM_VFADegrees[i]));
                    Mz = Mz * cos(UT_Radians(KAM_VFADegrees[i]));
                }
            }
            for( int i=0; i<nVFA; i++)
            {
                KAM_VFADegrees[i] = MAX_OF( MIN_OF( KAM_VFADegrees[i], maxVFADeg), 0.01);
                //KAM_VFAWatts[i] = ExcPulse1.Pow * pow( KAM_VFADegrees[i] / maxVFADeg, 2.0);
                if( ParxRelsParHasValue("PVM_RefPowCh1") == Yes) {
                    KAM_VFAWatts[i] = MRT_RFPulsePower(ExcPulse1.Sint, KAM_VFADegrees[i], ExcPulse1.Length, PVM_RefPowCh1, errString);
                    if(*errString)
                        UT_ReportError(errString);
                }
                KAM_VFAWatts[i] = MIN_OF( KAM_VFAWatts[i], maxRfPower);
            }
        }
        else // User specifies excitation powers in Watts
        {
            ParxRelsMakeNonEditable("KAM_VFADegrees, KAM_VFAEffectiveDegrees, KAM_VFAEffectiveOn");
            ParxRelsMakeEditable("KAM_VFAWatts");
            for( int i=0; i<nVFA; i++)
            {
                KAM_VFADegrees[i] = -1.0;
                KAM_VFAWatts[i] = MAX_OF( MIN_OF( KAM_VFAWatts[i], maxRfPower), 0.0);
            }
        }
    }
    else // A single excitation power is used for ExcPulse1
    {
        PARX_change_dims("KAM_VFADegrees",1);
        PARX_change_dims("KAM_VFAWatts",1);
        if( KAM_ExcMode == Normal_Exc)
        {
            ParxRelsMakeNonEditable("KAM_VFADegrees, KAM_VFAWatts");
            ParxRelsHideInEditor("KAM_VFADegrees, KAM_VFAWatts, KAM_VFAEffectiveOn, KAM_VFAEffectiveDegrees");
            UpdateRFPulses();
            KAM_VFAWatts[0] = ExcPulse1.Pow;
            KAM_VFADegrees[0] = ExcPulse1.Flipangle;
        } else {
            maxRfPower = CFG_MaxRFPower(1, PVM_Nucleus1, ExcPulse1.Pint, ExcPulse1.Length);
            if( PVM_DeriveGains == Yes ) // Flipangle set in KAM_VFADegrees[0]
            {
                ParxRelsMakeEditable("KAM_VFADegrees");
                ParxRelsMakeNonEditable("KAM_VFAWatts");
                KAM_VFADegrees[0] = MIN_OF(MAX_OF(KAM_VFADegrees[0], 0.01), 90.0);
                KAM_VFAWatts[0] = MRT_RFPulsePower(ExcPulse1.Sint, KAM_VFADegrees[0], ExcPulse1.Length, PVM_RefPowCh1, errString);
                if(*errString)
                    UT_ReportError(errString);
                KAM_VFAWatts[0] = MIN_OF( KAM_VFAWatts[0], maxRfPower);
                KAM_VFADegrees[0] = MRT_RFPulseAngle(ExcPulse1.Sint, KAM_VFAWatts[0], ExcPulse1.Length, PVM_RefPowCh1, errString);
                if( *errString)
                    UT_ReportError(errString);
                ExcPulse1.Pow = KAM_VFAWatts[0];
                ExcPulse1.Flipangle = KAM_VFADegrees[0];
            } else {    // Power set in KAM_VFAWatts[0]
                ParxRelsMakeNonEditable("KAM_VFADegrees");
                ParxRelsMakeEditable("KAM_VFAWatts");
                KAM_VFAWatts[0] = MIN_OF( KAM_VFAWatts[0], maxRfPower);
                ExcPulse1.Pow = KAM_VFAWatts[0];
                ExcPulse1.Flipangle = -1.0;
            }
            ParxRelsHideInEditor("KAM_VFAEffectiveOn, KAM_VFAEffectiveDegrees");
        }
    }

    DB_MSG(("<--HandleVFA"));
}


/* This function ensures consistency with spec-spat excitation and
 * arbitrary readout trajectory designs, and avoids having to
 * repeatedly call the relations for KAM_ExcMode and KAM_ReadMode */
void EnforceHpmrDesignParams( void )
{
    DB_MSG(("-->EnforceHpmrDesignParams"));

    // Excitation Parameters
    if( KAM_ExcMode == Normal_Exc)
    {
        KAM_ExcGradSize = 1;
        PARX_change_dims("KAM_ExcGradShape", KAM_ExcGradSize);
        KAM_ExcGradShape[0] = 0;
        KAM_NSlices = 1;
    }
    // else if( KAM_ExcMode == SpSp_13C_sb_10mm_20190730)
    // {
    //     // Enforce design parameters
    //     // Single slice at isocenter
    //     PVM_NSPacks = 1;
    //     PVM_SPackArrNSlices[0] = 1;
    //     PVM_SPackArrSliceOffset[0] = 0.0;
    //     // Slice thickness
    //     PVM_SliceThick = 10.0;  //mm
    //     KAM_NSlices = 1;
    // }


    // Readout Parameters
    if( KAM_ReadMode == Normal_Read)
    {
        KAM_NAcqPoints = PVM_Matrix[0];
        KAM_IdealEchoShift = 0.0;
        KAM_FidEcho = No;

        KAM_ReadShapeSize = 1;
        PARX_change_dims("KAM_ReadShapeX", KAM_ReadShapeSize);
        KAM_ReadShapeX[0] = 0;
        PARX_change_dims("KAM_ReadShapeY", KAM_ReadShapeSize);
        KAM_ReadShapeY[0] = 0;
        PARX_change_dims("KAM_ReadShapeZ", KAM_ReadShapeSize);
        KAM_ReadShapeZ[0] = 0;
    }
    // else if( KAM_ReadMode == flybackEPI_fov40_mtx16_dw48_13C_20190103)
    // {
    //     PVM_Matrix[0] = 391;
    //     PVM_EffSWh = 20833.333333;
    //     PVM_Fov[0] = 40;
    //     PVM_Fov[1] = 40;
    //     PVM_SPackArrReadOffset[0] = 0;
    //     PVM_SPackArrPhase1Offset[0] = 0;
    // }

    // allow readout to extend beyond encoding gradient
    KAM_NAcqPoints = MAX_OF(KAM_NAcqPoints, PVM_Matrix[0]);
    KAM_NAcqPoints = MIN_OF(KAM_NAcqPoints, 16384);


    DB_MSG(("<--EnforceHpmrDesignParams"));
}
