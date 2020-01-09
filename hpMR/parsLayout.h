/****************************************************************
 *
 * $Source: /pv/CvsTree/pv/gen/src/prg/methods/FLASH/parsLayout.h,v $
 *
 * Copyright (c) 1999-2007
 * Bruker BioSpin MRI GmbH
 * D-76275 Ettlingen, Germany
 *
 * All Rights Reserved
 *
 * $Id: parsLayout.h,v 1.42.2.1 2015/03/23 07:47:39 mawi Exp $
 *
 ****************************************************************/

/****************************************************************/
/*	PARAMETER CLASSES				       	*/
/****************************************************************/



/*--------------------------------------------------------------*
 * Definition of the PV class...
 *--------------------------------------------------------------*/

parclass
{
    PVM_EffSWh;
    PVM_EchoPosition;
    EncGradDur;
    PVM_AcquisitionTime;
    ReadSpoiler;
    SliceSpoiler;
    DigitizerPars;
}
attributes
{
    display_name "Sequence Details";
} Sequence_Details;

parclass
{
    ExcPulse1Enum;
    ExcPulse1;
    ExcPulse1Ampl;
    //ExcPulse1Shape;
    KAM_FidExcPulseEnum;
    KAM_FidExcPulse;
    KAM_FidExcPulseAmpl;
    //KAM_FidExcPulseShape;
    KAM_SatExcPulseEnum;
    KAM_SatExcPulse;
    KAM_SatExcPulseAmpl;
    //KAM_SatExcPulseShape;
}
attributes
{
    display_name "RF Pulses";
} RF_Pulses;

parclass
{
    DummyScans_Parameters;
    PVM_FreqDriftYN;

    PVM_NMovieFrames;
    TimeForMovieFrames;

    PVM_EvolutionOnOff;
    Evolution_Parameters;

    PVM_TriggerModule;
    Trigger_Parameters;

    PVM_TaggingOnOff;
    Tagging_Parameters;

    PVM_SelIrOnOff;
    Selective_IR_Parameters;

    PVM_BlBloodOnOff;
    BlackBlood_Parameters;

    PVM_FatSupOnOff;
    Fat_Sup_Parameters;

    PVM_MagTransOnOff;
    Magn_Transfer_Parameters;

    PVM_FovSatOnOff;
    Fov_Sat_Parameters;

    PVM_InFlowSatOnOff;
    Flow_Sat_Parameters;

    PVM_MotionSupOnOff;

    RFSpoiling;
    AngioMode;

} Preparation;

// KAM parameters specific to hpMR
parclass
{
    PVM_MethodVersion;

    KAM_GradientTime;
    KAM_GradientMaxAmp;
    KAM_GradientMaxSlew;
    KAM_GammaMHzT;

    KAM_EpiType;
    KAM_EpiMatrix;
    KAM_EpiDwellFactor;
    KAM_EpiPartialFourierY;

    KAM_AcqMat;
    KAM_NAcqPoints;
    KAM_NEchoes;
    KAM_NSlices;

    KAM_ExcMode;
    KAM_ExcShapeFile;
    KAM_ExcGradFile;
    KAM_ExcMetaFile;
    //KAM_ExcGradShape;
    KAM_ExcGradSize;
    KAM_ExcRfDelay;
    KAM_ExcGradScale;

    KAM_PrepMode;
    KAM_PrepGradOnOuter;
    KAM_PrepGradOnInner;
    KAM_PrepShapeSize;
    KAM_PrepLength;

    KAM_ReadMode;
    KAM_ReadShapeFile;
    //KAM_ReadShapeX;
    //KAM_ReadShapeY;
    //KAM_ReadShapeZ;
    KAM_ReadShapeSize;
    KAM_ReadLength;
    KAM_ReadAcqDelay;
    KAM_ReadGradScaleX;
    KAM_ReadGradScaleY;
    KAM_ReadGradScaleZ;
    KAM_ReadEchoDelay;

    KAM_IdealRepDelay;
    KAM_IdealMinDelay;
    KAM_IdealRepFill;
    KAM_IdealEchoShift;
    KAM_FidEcho;
    KAM_SatEcho;
    KAM_SatFreq;

    KAM_VFAOn;
    KAM_VFAEffectiveOn;
    KAM_VFAEffectiveDegrees;
    KAM_VFADegrees;
    KAM_VFAWatts;

    KAM_SpecExc;
    KAM_SpecTxRxSeparate;
    KAM_VOFTxHertz;
    KAM_VOFRxHertz;

    PVM_DecOnOff;
    Decoupling_Parameters;

    KAM_HSTrig;
}
attributes
{
    display_name "hpMR Parameters";
}hpMRGroup;

parclass
{
    Method;
    PVM_EchoTime;
    PVM_RepetitionTime;
    PVM_NEchoImages;
    PVM_NAverages;
    PVM_NRepetitions;
    PVM_ScanTimeStr;
    PVM_ScanTime;
    PVM_DeriveGains;
    hpMRGroup;
    RF_Pulses;
    Nuclei;
    Encoding;
    Sequence_Details;
    ImageGeometry;
    MapShim;
    StartupShims;
} MethodClass;

// parameters that should be tested after any editing
conflicts
{
    PVM_EchoTime;
    PVM_RepetitionTime;
    PVM_Fov;
    PVM_SliceThick;
};

// parameters for reconstruction
parclass
{
    RecoMethMode;
    WeightingMode;
    MaskWeighting;
    GaussBroadening;
}attributes
{
    display_name "Reconstruction Options";
}MethodRecoGroup;




/****************************************************************/
/*	E N D   O F   F I L E					*/
/****************************************************************/



