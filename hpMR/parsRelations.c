/****************************************************************
 *
 * $Source: /pv/CvsTree/pv/gen/src/prg/methods/FLASH/parsRelations.c,v $
 *
 * Copyright (c) 2002 - 2003
 * Bruker BioSpin MRI GmbH
 * D-76275 Ettlingen, Germany
 *
 * All Rights Reserved
 *
 *
 * $Id: parsRelations.c,v 1.73.2.1 2013/12/10 16:21:57 mawi Exp $
 *
 ****************************************************************/

static const char resid[] = "$Id: parsRelations.c,v 1.73.2.1 2013/12/10 16:21:57 mawi Exp $ (C) 2002 Bruker BioSpin MRI GmbH";

#define DEBUG		0
#define DB_MODULE	0
#define DB_LINE_NR	0


#include "method.h"

#include "create2dEpi_output.h"
#include "ovl_toolbox/Utils.h"


/*==========================================================
 *
 *  examples for relations concearning special pulses and
 *  pulselists
 *
 *==========================================================*/



/*===============================================================
 * ExcPulse1EnumRelation
 * Relation of ExcPulse1Enum (a dynamic enmueration parameter which
 * allows to select one of the existing  pulses)
 *===============================================================*/

void ExcPulse1EnumRelation(void)
{
    DB_MSG(("-->ExcPulse1EnumRelation"));

    UT_SetRequest("ExcPulse1Enum");
    backbone();

    DB_MSG(("<--ExcPulse1EnumRelation"));
}

/*===============================================================
 * ExcPulse1AmplRel
 * Relation of ExcPulseAmpl
 * This parameter is used in the setup parameter card to adjust
 * the RF pulse amplitude manually
 *===============================================================*/

void ExcPulse1AmplRel(void)
{
    DB_MSG(("-->ExcPulse1AmplRel"));
    UT_SetRequest("ExcPulse1Ampl");
    HandleRFPulseAmplitude();
    DB_MSG(("-->ExcPulse1AmplRel"));
}

void HandleRFPulseAmplitude(void)
{
    DB_MSG(("-->HandleRFPulseAmplitude"));

    STB_UpdateRFShapeAmplitude("ExcPulse1Ampl",No);
    ATB_SetRFPulse("ExcPulse1","ACQ_RfShapes[0]");

    DB_MSG(("<--HandleRFPulseAmplitude"));
}



/* ===================================================================
 * Relation of ExcPulse
 *
 * All pulses of type PVM_RF_PULSE must have relations like this.
 * However, if you clone this funtion for a different pulse parameter
 * remember to replace the param name in the call to UT_SetRequest!
 *
 * IMPORTANT: this function should not be invoked in the backbone!
 ====================================================================*/

void ExcPulse1Relation(void)
{
    DB_MSG(("-->ExcPulse1Relation"));

    /*
     * Tell the request handling system that the parameter
     * ExcPulse has been edited
     */

    UT_SetRequest("ExcPulse1");

    /* Check the values of ExcPulse */

    ExcPulse1Range();

    /*
     * call the backbone; further handling will take place there
     * (by means of STB_UpdateRFPulse)
     */

    backbone();

    DB_MSG(("<--ExcPulse1Relation"));
}



void ExcPulse1Range(void)
{
    DB_MSG(("-->ExcPulse1Range"));

    // range checker fields to be controlled may be
    // .Length
    // .Bandwidth
    // .Flipangle
    // .Calculated
    // .Sharpness
    // .Flatness
    double dval=ExcPulse1.Flipangle;

    ExcPulse1.Flipangle = MIN_OF(90.0,MAX_OF(dval,1.0));

    DB_MSG(("<--ExcPulseRange"));
}



/*===============================================================
 *
 * Range checking routine for parameter PVM_NAverages
 *
 *==============================================================*/


void NAveragesRange(void)
{
    DB_MSG(("-->NAveragesRange\n"));

    /*
     *  Range check of PVM_NAverages: prevent it to be negative or 0
     */

    if(ParxRelsParHasValue("PVM_NAverages") == No)
    {
        PVM_NAverages = 1;
    }

    if (PVM_NAverages > 1)
    {
        ParxRelsShowInEditor("PVM_MotionSupOnOff");
    }
    else
    {
        ParxRelsHideInEditor("PVM_MotionSupOnOff");
    }

    DB_MSG(("<--NAveragesRange\n"));
}



void NAveragesRels(void)
{

    DB_MSG(("-->NAveragesRels\n"));


    NAveragesRange();

    /*
     *   Averages range check is finished, handle the request by
     *   the method:
     */

    backbone();


    DB_MSG(("<--NAveragesRels\n"));
    return;
}





/* rangechecking and redirected relations of PVM_EffSWh */

void EffSWhRange(void)
{
    DB_MSG(("-->EffSWhRange"));

    if(!ParxRelsParHasValue("PVM_EffSWh"))
    {
        PVM_EffSWh = 50000;
    }
    else
    {
        PVM_EffSWh = MIN_OF(MAX_OF(PVM_EffSWh,2000.0),1000000);
    }

    DB_MSG(("-->EffSWhRange"));
    return;
}

void EffSWhRel(void)
{
    DB_MSG(("-->EffSWhRel"));

    EffSWhRange();
    backbone();

    DB_MSG(("-->EffSWhRel"));
    return;
}

void InversionTimeRels(void)
{
    DB_MSG(("-->InversionTimeRel"));


    InversionTimeRange();

    if(PVM_SelIrOnOff==On)
        PVM_SelIrInvTime = PVM_InversionTime;
    if(PVM_BlBloodOnOff==On)
        PVM_BlBloodInvTime = PVM_InversionTime;

    backbone();
    DB_MSG(("-->InversionTimeRel"));

}

void InversionTimeRange(void)
{
    if(!ParxRelsParHasValue("PVM_InversionTime"))
        PVM_InversionTime = 0.0;

    PVM_InversionTime = MAX_OF(PVM_InversionTime,0.0);
}


/* relations of read/slice spoiler */
void ReadSpoilerRel(void)
{
    DB_MSG(("-->ReadSpoilerRel"));
    UT_SetRequest("ReadSpoiler");
    backbone();
    DB_MSG(("<--ReadSpoilerRel"));
}

void SliceSpoilerRel(void)
{
    DB_MSG(("-->SliceSpoilerRel"));
    UT_SetRequest("SliceSpoiler");
    backbone();
    DB_MSG(("<--SliceSpoilerRel"));
}

void RecoMethModeVisPar(void)
{
    DB_MSG(("-->RecoMethModeVisPar"));

    if(RecoMethMode==SWI)
    {
        ParxRelsShowInEditor("WeightingMode,GaussBroadening");
        if(WeightingMode==phase_image)
            ParxRelsHideInEditor("MaskWeighting");
        else
            ParxRelsShowInEditor("MaskWeighting");
    }
    else
        ParxRelsHideInEditor("WeightingMode,GaussBroadening,MaskWeighting");


    DB_MSG((">--RecoMethModeVisPar"));
}

void MaskModeRel(void)
{
    DB_MSG(("-->MaskModeRel"));

    if(WeightingMode==phase_image)
        ParxRelsHideInEditor("MaskWeighting");
    else
        ParxRelsShowInEditor("MaskWeighting");
    SetRecoParam();

    DB_MSG(("<--MaskModeRel"));
}

void GaussBroadRange(void)
{
    double max;

    DB_MSG(("-->GaussBroadRange"));

    if(ParxRelsParHasValue("GaussBroadening") == No)
        GaussBroadening=1.0;

    max=PTB_GetSpatDim()>1? MAX_OF(PVM_Fov[0],PVM_Fov[1]) : PVM_Fov[0];
    if(PTB_GetSpatDim()==3)
        max=MAX_OF(max,PVM_Fov[2]);

    GaussBroadening = MIN_OF(MAX_OF(0,GaussBroadening),max);

    DB_MSG(("<--GaussBroadRange"));
}

void MaskWeightRange(void)
{
    DB_MSG(("-->MaskWeightRange"));

    if(ParxRelsParHasValue("MaskWeighting") == No)
        MaskWeighting=4.0;

    MaskWeighting = MIN_OF(MAX_OF(0,MaskWeighting),20);

    DB_MSG(("<--MaskWeightRange"));
}

void MyRgInitRel(void)
{

    DB_MSG(("-->MyRgInitRel"));
    if(AngioMode==Yes)
    {
        /* initialize 1 slice per package, upadte without angio mode
           for RG adjustment
           */
        int packs=PVM_NSPacks;
        int *ppersl=PVM_SPackArrNSlices;
        for(int i=0; i<packs;i++)
            ppersl[i]=1;
        AngioMode=No;
        backbone();

    }
    ParxRelsParRelations("PVM_AutoRgInitHandler",Yes);
    DB_MSG(("<--MyRgInitRel"));
}

/* KAM add for hpMR */

// Relations for separate fid excitation
void KAM_FidExcPulseEnumRelation(void)
{
    DB_MSG(("-->KAM_FidExcPulseEnumRelation"));

    UT_SetRequest("KAM_FidExcPulseEnum");
    backbone();

    DB_MSG(("<--KAM_FidExcPulseEnumRelation"));
}

void KAM_FidExcPulseAmplRel(void)
{
    DB_MSG(("-->KAM_FidExcPulseAmplRel"));
    UT_SetRequest("KAM_FidExcPulseAmpl");
    KAM_FidHandleRFPulseAmplitude();
    DB_MSG(("-->KAM_FidExcPulseAmplRel"));
}

void KAM_FidHandleRFPulseAmplitude(void)
{
    DB_MSG(("-->KAM_FidHandleRFPulseAmplitude"));

    STB_UpdateRFShapeAmplitude("KAM_FidExcPulseAmpl",No);
    ATB_SetRFPulse("KAM_FidExcPulse","ACQ_RfShapes[1]");

    DB_MSG(("<--KAM_FidHandleRFPulseAmplitude"));
}

void KAM_FidExcPulseRelation(void)
{
    DB_MSG(("-->KAM_FidExcPulseRelation"));

    /*
     * Tell the request handling system that the parameter
     * ExcPulse has been edited
     */
    UT_SetRequest("KAM_FidExcPulse");

    /* Check the values of KAM_FidExcPulse */
    KAM_FidExcPulseRange();

    /*
     * call the backbone; further handling will take place there
     * (by means of STB_UpdateRFPulse)
     */
    backbone();

    DB_MSG(("<--KAM_FidExcPulseRelation"));
}

void KAM_FidExcPulseRange(void)
{
    DB_MSG(("-->KAM_FidExcPulseRange"));

    // range checker fields to be controlled may be
    // .Length
    // .Bandwidth
    // .Flipangle
    // .Calculated
    // .Sharpness
    // .Flatness
    double dval=KAM_FidExcPulse.Flipangle;
    KAM_FidExcPulse.Flipangle = MIN_OF(90.0,MAX_OF(dval,1.0));

    // Fid excitation pulse can't be longer than other excitation pulses
    double ddval=KAM_FidExcPulse.Length;
    KAM_FidExcPulse.Length = MIN_OF(ExcPulse1.Length,ddval);

    DB_MSG(("<--KAM_FidExcPulseRange"));
}

// Relations for separate fid excitation
void KAM_SatExcPulseEnumRelation(void)
{
    DB_MSG(("-->KAM_SatExcPulseEnumRelation"));

    UT_SetRequest("KAM_SatExcPulseEnum");
    backbone();

    DB_MSG(("<--KAM_SatExcPulseEnumRelation"));
}

void KAM_SatExcPulseAmplRel(void)
{
    DB_MSG(("-->KAM_SatExcPulseAmplRel"));
    UT_SetRequest("KAM_SatExcPulseAmpl");
    KAM_SatHandleRFPulseAmplitude();
    DB_MSG(("-->KAM_SatExcPulseAmplRel"));
}

void KAM_SatHandleRFPulseAmplitude(void)
{
    DB_MSG(("-->KAM_SatHandleRFPulseAmplitude"));

    STB_UpdateRFShapeAmplitude("KAM_SatExcPulseAmpl",No);
    ATB_SetRFPulse("KAM_SatExcPulse","ACQ_RfShapes[1]");

    DB_MSG(("<--KAM_SatHandleRFPulseAmplitude"));
}

void KAM_SatExcPulseRelation(void)
{
    DB_MSG(("-->KAM_SatExcPulseRelation"));

    /*
     * Tell the request handling system that the parameter
     * ExcPulse has been edited
     */
    UT_SetRequest("KAM_SatExcPulse");

    /* Check the values of KAM_SatExcPulse */
    KAM_SatExcPulseRange();

    /*
     * call the backbone; further handling will take place there
     * (by means of STB_UpdateRFPulse)
     */
    backbone();

    DB_MSG(("<--KAM_SatExcPulseRelation"));
}

void KAM_SatExcPulseRange(void)
{
    DB_MSG(("-->KAM_SatExcPulseRange"));

    // range checker fields to be controlled may be
    // .Length
    // .Bandwidth
    // .Flipangle
    // .Calculated
    // .Sharpness
    // .Flatness
    double dval=KAM_SatExcPulse.Flipangle;
    KAM_SatExcPulse.Flipangle = MIN_OF(90.0,MAX_OF(dval,1.0));

    // Sat excitation pulse can't be longer than other excitation pulses
    double ddval=KAM_SatExcPulse.Length;
    KAM_SatExcPulse.Length = MIN_OF(ExcPulse1.Length,ddval);

    DB_MSG(("<--KAM_SatExcPulseRange"));
}


// Relations for excitation
void KAM_ExcModeRel(void)
{
    DB_MSG(("-->KAM_ExcModeRel"));

    UT_SetRequest("KAM_ExcMode");

    if( KAM_SpecExc == Yes)
    {
        // Method will show a single slice at isocenter
        PVM_NSPacks = 1;
        PVM_SPackArrNSlices[0] = 1;
        PVM_SPackArrSliceOffset[0] = 0.0;

        ParxRelsMakeNonEditable("PVM_NSPacks, PVM_SPackArrNSlices, PVM_SPackArrSliceOffset");
        ParxRelsShowInEditor("KAM_SpecTxRxSeparate, KAM_VOFTxHertz, KAM_VOFRxHertz, KAM_SatEcho, KAM_SatFreq");
        if( KAM_SpecTxRxSeparate == Yes)
            ParxRelsMakeNonEditable("KAM_VOFRxHertz");
        else
            ParxRelsMakeEditable("KAM_VOFRxHertz");
    } else {
        KAM_SatEcho = No;
        ParxRelsMakeEditable("KAM_ExcMode, PVM_NSPacks, PVM_SPackArrNSlices, PVM_SPackArrSliceOffset");
        ParxRelsHideInEditor("KAM_SpecTxRxSeparate, KAM_VOFTxHertz, KAM_VOFRxHertz, KAM_SatEcho, KAM_SatFreq");
    }

    if( KAM_FidEcho == Yes) {
        ParxRelsShowInEditor("KAM_FidExcPulse, KAM_FidExcPulseEnum");
    } else {
        ParxRelsHideInEditor("KAM_FidExcPulse, KAM_FidExcPulseEnum");
    }

    if( KAM_SatEcho == Yes) {
        ParxRelsShowInEditor("KAM_SatExcPulse, KAM_SatExcPulseEnum");
    } else {
        ParxRelsHideInEditor("KAM_SatExcPulse, KAM_SatExcPulseEnum");
    }

    if( KAM_ExcMode == Normal_Exc)
    {
        ParxRelsMakeEditable("ExcPulse1, ExcPulse1Enum, PVM_NSPacks, PVM_SPackArrNSlices, PVM_SPackArrSliceOffset, PVM_SliceThick");
        ParxRelsHideInEditor("KAM_NSlices");

        strcpy(KAM_ExcShapeFile,"NULL");
        strcpy(KAM_ExcGradFile,"NULL");
        strcpy(KAM_ExcMetaFile,"NULL");

        KAM_ExcGradSize = 1;
        PARX_change_dims("KAM_ExcGradShape", KAM_ExcGradSize);
        KAM_ExcGradShape[0] = 0;
        KAM_ExcRfDelay = 0.0;
        KAM_ExcGradScale = 0.0;
    }
    else
    {
        ParxRelsMakeNonEditable("ExcPulse1, ExcPulse1Enum, PVM_NSPacks, PVM_SPackArrNSlices, PVM_SPackArrSliceOffset, PVM_SliceThick");
        ParxRelsShowInEditor("KAM_NSlices");
        KAM_ExcRfDelay = 0.0;      //us
        KAM_ExcGradScale = 1.0;

        strcpy(KAM_ExcGradFile, getenv("XWINNMRHOME"));
        strcat(KAM_ExcGradFile, "/prog/curdir/");
        strcat(KAM_ExcGradFile, getenv("USER"));
        strcat(KAM_ExcGradFile, "/ParaVision/exp/lists/gp/");
        strcpy(KAM_ExcMetaFile, getenv("XWINNMRHOME"));
        strcat(KAM_ExcMetaFile, "/prog/curdir/");
        strcat(KAM_ExcMetaFile, getenv("USER"));
        strcat(KAM_ExcMetaFile, "/ParaVision/exp/lists/wave/");

        // if( KAM_ExcMode == SpSp_13C_sb_10mm_20190730)
        // {
        //     PVM_NSPacks = 1;
        //     PVM_SPackArrNSlices[0] = 1;
        //     PVM_SPackArrSliceOffset[0] = 0.0;
        //     PVM_SliceThick = 10.0;     //mm
        //     KAM_ExcRfDelay += 0.0;    //us, added to base delay from above
        //     KAM_ExcGradScale *= 1.0;
        //     KAM_NSlices = 1;

        //     strcpy(KAM_ExcShapeFile,"SpSp_13C_sb_10mm_20190730.exc");
        //     strcat(KAM_ExcGradFile,"SpSp_13C_sb_10mm_20190730.gp");
        //     strcat(KAM_ExcMetaFile,"SpSp_13C_sb_10mm_20190730.meta");
        // }

        // RF parameters
        strcpy(ExcPulse1.Shape, KAM_ExcShapeFile);
        FILE * metaFile;
        metaFile = fopen(KAM_ExcMetaFile,"r");
        // excitation metaFiles contain 3 values:
        //  RF pulse length (ms)
        //  RF pulse shape integral
        //  Number of gradient waveform points
        if( metaFile == NULL)
            UT_ReportError("RF pulse meta file not found.");
        else
            fscanf(metaFile, "%lf %lf %d",
                    &ExcPulse1.Length, &ExcPulse1.Sint, &KAM_ExcGradSize);
        fclose(metaFile);

        double plen = ExcPulse1.Length;
        double glen = KAM_GradientTime * KAM_ExcGradSize / 1e3;
        if( fabs(plen - glen) > 1e-6)
        {
            char errString[512];
            sprintf(errString,
                    "%s \n RF Length: %f ms \n Gradient Length: %f ms \n %s",
                    "Excitation gradient and RF pulse length do not match.",
                    plen, glen,
                    "This pulse may have been designed for a different gradient time resolution.");
            UT_ReportError(errString);
        }

        // Calculate power integral ratio
        /* CFG_RFPulseGetPowIntFac(KAM_ExcShapeFile, &ExcPulse1.Pint); */
        ExcPulse1.Pint = 1.0; // Worst-case scenario for hardware protection

        // TODO: calculate rephasing factor in MATLAB pulse design (TE values are not strictly correct)
        ExcPulse1.Rpfac = 50.0;

        // See EnforceHpmrDesignParams() and HandleVFA() in backbone.c for power handling

        // Gradient parameters
        PARX_change_dims("KAM_ExcGradShape", KAM_ExcGradSize);
        FILE * gradFile;
        gradFile = fopen(KAM_ExcGradFile,"r");
        if( gradFile == NULL)
            UT_ReportError("Excitation gradient file not found.");
        else
        {
            for( int ii=0; ii<KAM_ExcGradSize; ii++)
                fscanf(gradFile, "%lf", &KAM_ExcGradShape[ii]);

            // Double check that gradients obey system limits
            double gradMax, slewMax;
            double gradLim = CFG_MaxGradientStrength();
            double slewLim = CFG_MaxGradientStrength() / (CFG_GradientRampTime() + CFG_InterGradientWaitTime());
            KAM_GradCheck(KAM_ExcGradShape, KAM_ExcGradSize, gradLim, slewLim, &gradMax, &slewMax);

            for( int ii=0; ii<KAM_ExcGradSize; ii++)
            {
                // Convert from mT/m to fraction of max grad
                KAM_ExcGradShape[ii] = KAM_ExcGradScale*KAM_ExcGradShape[ii]/KAM_GradientMaxAmp;
            }
        }
        fclose(gradFile);
    }

    backbone();
    DB_MSG(("<--KAM_ExcModeRel"));
}

// Relations for preparation
void KAM_PrepModeRel(void)
{
    DB_MSG(("-->KAM_PrepModeRel"));

    UT_SetRequest("KAM_PrepMode");
    if( KAM_PrepMode == No_Prep)
    {
        strcpy(KAM_PrepShapeFile,"NULL");

        KAM_PrepShapeSize = 1;
        KAM_PrepLength = 0.0;
        PARX_change_dims("KAM_PrepShapeX", KAM_PrepShapeSize);
        KAM_PrepShapeX[0] = 0;
        PARX_change_dims("KAM_PrepShapeY", KAM_PrepShapeSize);
        KAM_PrepShapeY[0] = 0;
        PARX_change_dims("KAM_PrepShapeZ", KAM_PrepShapeSize);
        KAM_PrepShapeZ[0] = 0;
    }
    else
    {
        strcpy(KAM_PrepShapeFile, getenv("XWINNMRHOME"));
        strcat(KAM_PrepShapeFile, "/prog/curdir/");
        strcat(KAM_PrepShapeFile, getenv("USER"));
        strcat(KAM_PrepShapeFile, "/ParaVision/exp/lists/gp/");

        //if( KAM_PrepMode == diff_13C_iso4pulse_b30)
        //{
        //    KAM_PrepShapeSize = 2008;
        //    KAM_PrepLength = 16.064;
        //    strcat(KAM_PrepShapeFile, "diff_13C_iso4pulse_b30.gp");
        //}

        double plen = KAM_PrepLength;
        double glen = KAM_GradientTime * KAM_PrepShapeSize / 1e3;
        if( fabs(plen - glen) > 1e-6)
        {
            char errString[512];
            sprintf(errString,
                    "%s \n Expected Length: %f ms \n Gradient Length: %f ms \n %s",
                    "Expected and actual preparation gradient lengths do not match.",
                    plen, glen,
                    "This pulse may have been designed for a different gradient time resolution.");
            UT_ReportError(errString);
        }

        // Preparation gradient parameters
        PARX_change_dims("KAM_PrepShapeX", KAM_PrepShapeSize);
        PARX_change_dims("KAM_PrepShapeY", KAM_PrepShapeSize);
        PARX_change_dims("KAM_PrepShapeZ", KAM_PrepShapeSize);
        FILE * gradFile;
        gradFile = fopen(KAM_PrepShapeFile,"r");
        if( gradFile == NULL)
            UT_ReportError("Preparation gradient file not found.");
        else
        {
            for( int ii=0; ii<KAM_PrepShapeSize; ii++)
                fscanf(gradFile, "%lf %lf %lf", &KAM_PrepShapeX[ii], &KAM_PrepShapeY[ii], &KAM_PrepShapeZ[ii]);

            // Double check that gradients obey system limits
            double gradMax, slewMax;
            double gradLim = CFG_MaxGradientStrength();
            double slewLim = CFG_MaxGradientStrength() / (CFG_GradientRampTime() + CFG_InterGradientWaitTime());
            KAM_GradCheck(KAM_PrepShapeX, KAM_PrepShapeSize, gradLim, slewLim, &gradMax, &slewMax);
            KAM_GradCheck(KAM_PrepShapeY, KAM_PrepShapeSize, gradLim, slewLim, &gradMax, &slewMax);
            KAM_GradCheck(KAM_PrepShapeZ, KAM_PrepShapeSize, gradLim, slewLim, &gradMax, &slewMax);

            for( int ii=0; ii<KAM_PrepShapeSize; ii++)
            {
                // Convert from mT/m to fraction of max grad
                KAM_PrepShapeX[ii] = KAM_PrepShapeX[ii]/KAM_GradientMaxAmp;
                KAM_PrepShapeY[ii] = KAM_PrepShapeY[ii]/KAM_GradientMaxAmp;
                KAM_PrepShapeZ[ii] = KAM_PrepShapeZ[ii]/KAM_GradientMaxAmp;
            }
        }
    fclose(gradFile);
}

backbone();
DB_MSG(("<--KAM_PrepModeRel"));
}

// Relations for readout
void KAM_ReadModeRel(void)
{
    DB_MSG(("-->KAM_ReadModeRel"));

    UT_SetRequest("KAM_ReadMode");

    STB_UpdateNuclei(No);
    KAM_GradientMaxAmp   = MAX_OF(MIN_OF(KAM_GradientMaxAmp, CFG_MaxGradientStrength()), 10.0);
    KAM_GradientMaxSlew  = MAX_OF(MIN_OF(KAM_GradientMaxSlew, CFG_MaxGradientStrength() /
                (CFG_GradientRampTime() + CFG_InterGradientWaitTime())), 100.0);
    KAM_GammaMHzT = 42.5764 * CFG_GammaRatio(PVM_Nucleus1);
    KAM_EpiMatrix = MAX_OF(MIN_OF(KAM_EpiMatrix, 256), 4);
    KAM_EpiDwellFactor = MAX_OF(MIN_OF(KAM_EpiDwellFactor, 15), 2);
    if( KAM_EpiDwellFactor % 2 > 0)
        KAM_EpiDwellFactor += 1;
    KAM_EpiPartialFourierY = MAX_OF(MIN_OF(KAM_EpiPartialFourierY, 1.0), 0.5);

    if( KAM_ReadMode == Normal_Read)
    {
        ParxRelsMakeEditable("PVM_EffSWh, PVM_Fov, PVM_SPackArrReadOffset, PVM_SPackArrPhase1Offset");
        ParxRelsMakeNonEditable("KAM_FidEcho, KAM_NEchoes, KAM_NAcqPoints");
        //ParxRelsMakeNonEditable("KAM_FidEcho, KAM_NEchoes, KAM_NAcqPoints, KAM_EpiType, KAM_EpiMatrix, KAM_EpiDwellFactor, KAM_EpiPartialFourierY");

        strcpy(KAM_ReadShapeFile,"NULL");

        KAM_ReadShapeSize = 1;
        KAM_ReadLength = 0.0;
        PARX_change_dims("KAM_ReadShapeX", KAM_ReadShapeSize);
        KAM_ReadShapeX[0] = 0;
        PARX_change_dims("KAM_ReadShapeY", KAM_ReadShapeSize);
        KAM_ReadShapeY[0] = 0;
        PARX_change_dims("KAM_ReadShapeZ", KAM_ReadShapeSize);
        KAM_ReadShapeZ[0] = 0;
        KAM_ReadAcqDelay = 0.0;
        KAM_ReadGradScaleX = 1.0;
        KAM_ReadGradScaleY = 1.0;
        KAM_ReadGradScaleZ = 1.0;
        KAM_ReadEchoDelay = 0.0;
    }
    else if( KAM_ReadMode == create2dEpi)
    {
        ParxRelsMakeEditable("KAM_FidEcho, KAM_NEchoes, KAM_NAcqPoints, KAM_EpiType, KAM_EpiMatrix, KAM_EpiDwellFactor, KAM_EpiPartialFourierY");
        ParxRelsMakeNonEditable("PVM_EffSWh, PVM_SPackArrReadOffset, PVM_SPackArrPhase1Offset");
        KAM_ReadAcqDelay = 0.0;        //us
        KAM_ReadGradScaleX = 1.0;
        KAM_ReadGradScaleY = 1.0;
        KAM_ReadGradScaleZ = 1.0;

        // Enforce square FOV, limit ranges of EPI design parameters
        PVM_Fov[0] = MAX_OF(MIN_OF(PVM_Fov[0], 200.0), 10.0); // [1, 20] cm
        PVM_Fov[1] = PVM_Fov[0];

        // Initialize create2dEpi_codegen and set outputs to default values
        strcpy(KAM_ReadShapeFile,"NULL");
        double GX[65535], GY[65535] = {0.0};
        KAM_ReadShapeSize = 1;
        PARX_change_dims("KAM_ReadShapeX", KAM_ReadShapeSize);
        KAM_ReadShapeX[0] = 0.0;
        PARX_change_dims("KAM_ReadShapeY", KAM_ReadShapeSize);
        KAM_ReadShapeY[0] = 0.0;
        PARX_change_dims("KAM_ReadShapeZ", KAM_ReadShapeSize);
        KAM_ReadShapeZ[0] = 0.0;
        KAM_ReadEchoDelay = 0.0;

        // Run create2dEpi
        int epiType  = (KAM_EpiType == Flyback) ? 0 : 1;
        double fovCm = PVM_Fov[0] / 10;
        double dw    = KAM_GradientTime * KAM_EpiDwellFactor;
        char scanPath[256];
        PvOvlUtilGetExpnoPath(scanPath, 256, "");
        int mtxOut;
        double pfyOut;
        create2dEpi_output(epiType, fovCm, KAM_EpiMatrix, dw,
                KAM_GradientMaxAmp, KAM_GradientMaxSlew, KAM_GradientTime,
                KAM_GammaMHzT, KAM_EpiPartialFourierY, scanPath, GX, GY,
                &KAM_ReadShapeSize, &mtxOut, &KAM_ReadEchoDelay, &pfyOut);
        if( KAM_ReadShapeSize != 1){
            KAM_ReadLength = KAM_ReadShapeSize * KAM_GradientTime / 1e3;
            PVM_EffSWh = 1e6 / dw;
            PVM_Matrix[0] = mtxOut;
            KAM_EpiPartialFourierY = pfyOut;
            //printf("NGRAD = %d\nNACQ = %d\nTE = %g\n", KAM_ReadShapeSize, PVM_Matrix[0], KAM_ReadEchoDelay);
            strcpy(KAM_ReadShapeFile, "create2dEpi");
        } else {
            printf("Error encountered in create2dEpi, resetting Readout Mode to Normal\n");
            KAM_ReadMode = Normal_Read;
            KAM_ReadModeRel();
        }

        // Double check that gradients obey system limits
        double gradMax, slewMax;
        double gradLim = CFG_MaxGradientStrength();
        double slewLim = CFG_MaxGradientStrength() / (CFG_GradientRampTime() + CFG_InterGradientWaitTime());
        KAM_GradCheck(GX, KAM_ReadShapeSize, gradLim, slewLim, &gradMax, &slewMax);
        KAM_GradCheck(GY, KAM_ReadShapeSize, gradLim, slewLim, &gradMax, &slewMax);

        // Collect gradient waveform in KAM_ReadShape vectors
        PARX_change_dims("KAM_ReadShapeX", KAM_ReadShapeSize);
        PARX_change_dims("KAM_ReadShapeY", KAM_ReadShapeSize);
        PARX_change_dims("KAM_ReadShapeZ", KAM_ReadShapeSize);
        for(int ii=0; ii<KAM_ReadShapeSize; ii++) {
            // Convert from mT/m to fraction of max grad
            KAM_ReadShapeX[ii] = KAM_ReadGradScaleX * GX[ii] / KAM_GradientMaxAmp;
            KAM_ReadShapeY[ii] = KAM_ReadGradScaleY * GY[ii] / KAM_GradientMaxAmp;
            KAM_ReadShapeZ[ii] = KAM_ReadGradScaleZ * 0.0 / KAM_GradientMaxAmp;
        }

    }
    else
    {
        ParxRelsMakeNonEditable("PVM_EffSWh, PVM_Fov, PVM_SPackArrReadOffset, PVM_SPackArrPhase1Offset, KAM_EpiType, KAM_EpiMatrix, KAM_EpiDwellFactor, KAM_EpiPartialFourierY");
        ParxRelsMakeEditable("KAM_FidEcho, KAM_NEchoes, KAM_NAcqPoints");
        KAM_ReadAcqDelay = 0.0;        //us
        KAM_ReadGradScaleX = 1.0;
        KAM_ReadGradScaleY = 1.0;
        KAM_ReadGradScaleZ = 1.0;
        KAM_ReadEchoDelay = 0.0;        //ms

        strcpy(KAM_ReadShapeFile, getenv("XWINNMRHOME"));
        strcat(KAM_ReadShapeFile, "/prog/curdir/");
        strcat(KAM_ReadShapeFile, getenv("USER"));
        strcat(KAM_ReadShapeFile, "/ParaVision/exp/lists/gp/");

        // if( KAM_ReadMode == flybackEPI_fov40_mtx16_dw48_13C_20190103)
        // {
        //     KAM_ReadShapeSize = 2345;
        //     KAM_ReadLength = 18.76;
        //     PVM_Matrix[0] = 391;
        //     PVM_EffSWh = 20833.333333;
        //     PVM_Fov[0] = 40;
        //     PVM_Fov[1] = 40;
        //     PVM_SPackArrReadOffset[0] = 0;
        //     PVM_SPackArrPhase1Offset[0] = 0;
        //     KAM_ReadAcqDelay += 18.0;      //us, added to base delay from above
        //     KAM_ReadGradScaleX *= 1.0;
        //     KAM_ReadGradScaleY *= 1.0;
        //     KAM_ReadGradScaleZ *= 1.0;
        //     KAM_ReadEchoDelay += 9.984;   //ms, added to base delay from above

        //     strcat(KAM_ReadShapeFile,"flybackEPI_fov40_mtx16_dw48_13C_20190103.gp");
        // }

        double plen = KAM_ReadLength;
        double glen = KAM_GradientTime * KAM_ReadShapeSize / 1e3;
        if( fabs(plen - glen) > 1e-6)
        {
            char errString[512];
            sprintf(errString,
                    "%s \n Expected Length: %f ms \n Gradient Length: %f ms \n %s",
                    "Expected and actual readout gradient lengths do not match.",
                    plen, glen,
                    "This pulse may have been designed for a different gradient time resolution.");
            UT_ReportError(errString);
        }

        // Readout gradient parameters
        PARX_change_dims("KAM_ReadShapeX", KAM_ReadShapeSize);
        PARX_change_dims("KAM_ReadShapeY", KAM_ReadShapeSize);
        PARX_change_dims("KAM_ReadShapeZ", KAM_ReadShapeSize);
        FILE * gradFile;
        gradFile = fopen(KAM_ReadShapeFile,"r");
        if( gradFile == NULL)
            UT_ReportError("Readout gradient file not found.");
        else
        {
            for( int ii=0; ii<KAM_ReadShapeSize; ii++)
                fscanf(gradFile, "%lf %lf %lf", &KAM_ReadShapeX[ii], &KAM_ReadShapeY[ii], &KAM_ReadShapeZ[ii]);

            // Double check that gradients obey system limits
            double gradMax, slewMax;
            double gradLim = CFG_MaxGradientStrength();
            double slewLim = CFG_MaxGradientStrength() / (CFG_GradientRampTime() + CFG_InterGradientWaitTime());
            KAM_GradCheck(KAM_ReadShapeX, KAM_ReadShapeSize, gradLim, slewLim, &gradMax, &slewMax);
            KAM_GradCheck(KAM_ReadShapeY, KAM_ReadShapeSize, gradLim, slewLim, &gradMax, &slewMax);

            for( int ii=0; ii<KAM_ReadShapeSize; ii++)
            {
                // Convert from mT/m to fraction of max grad
                KAM_ReadShapeX[ii] = KAM_ReadGradScaleX*KAM_ReadShapeX[ii]/KAM_GradientMaxAmp;
                KAM_ReadShapeY[ii] = KAM_ReadGradScaleY*KAM_ReadShapeY[ii]/KAM_GradientMaxAmp;
                KAM_ReadShapeZ[ii] = KAM_ReadGradScaleZ*KAM_ReadShapeZ[ii]/KAM_GradientMaxAmp;
            }
        }
        fclose(gradFile);
    }

    backbone();
    DB_MSG(("<--KAM_ReadModeRel"));
}

void KAM_GradCheck(double* gradWave, int nGrad, double gradLim, double slewLim, double* gradMax, double* slewMax)
{
    DB_MSG(("-->KAM_GradCheck"));

    // Return the max gradient value and slew rate of input 1D gradient waveform
    // If these values exceed the input limitations an error is reported
    // Gradient values are [mT/m], slew rates are [T/m/s]

    if (nGrad > 65535)
        UT_ReportError("Gradient waveform longer than 65535 points");

    double gradMag[65535] = {};
    double slewRate[65535] = {};
    double tmp;

    for (int ii=0; ii<nGrad; ii++)
        gradMag[ii] = abs(gradWave[ii]);
    tmp = UT_MaxOfDoubleArray(nGrad, gradMag);
    if (tmp > gradLim)
        UT_ReportError("Maximum gradient amplitude exceeded");
    *gradMax = tmp;

    for (int ii=1; ii<nGrad; ii++)
        slewRate[ii-1] = 1e3 * abs(gradWave[ii] - gradWave[ii-1]) / KAM_GradientTime;
    tmp = UT_MaxOfDoubleArray(nGrad-1, slewRate);
    if (tmp > slewLim)
        UT_ReportError("Maximum gradient slew rate exceeded");
    *slewMax = tmp;

    DB_MSG(("<--KAM_GradCheck"));
}

/****************************************************************/
/*		E N D   O F   F I L E				*/
/****************************************************************/
