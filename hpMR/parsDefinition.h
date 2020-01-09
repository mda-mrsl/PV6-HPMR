/****************************************************************
 *
 * $Source: /pv/CvsTree/pv/gen/src/prg/methods/FLASH/parsDefinition.h,v $
 *
 * Copyright (c) 1999-2003
 * Bruker BioSpin MRI GmbH
 * D-76275 Ettlingen, Germany
 *
 * All Rights Reserved
 *
 * $Id: parsDefinition.h,v 1.30 2012/11/27 12:48:01 anba Exp $
 *
 ****************************************************************/



/****************************************************************/
/* INCLUDE FILES						*/
/****************************************************************/
double parameter OneRepTime;

PVM_SPOILER_TYPE parameter
{
    display_name "Read Spoiler";
    relations ReadSpoilerRel;
}ReadSpoiler;

PVM_SPOILER_TYPE parameter
{
    display_name "Slice Spoiler";
    relations SliceSpoilerRel;
}SliceSpoiler;

PV_PULSE_LIST parameter
{
    display_name "Excitation Pulse Shape";
    relations    ExcPulse1EnumRelation;
}ExcPulse1Enum;


PVM_RF_PULSE parameter
{
    display_name "Excitation Pulse";
    relations    ExcPulse1Relation;
}ExcPulse1;

PVM_RF_PULSE_AMP_TYPE parameter
{
    display_name "RF Pulse Amplitude";
    relations ExcPulse1AmplRel;
}ExcPulse1Ampl;

double parameter
{
    editable false;
}ExcPulse1Shape[];


double parameter
{
    display_name "Inter Slice Delay";
    relations backbone;
    units "ms";
    format "%.2f";
}SliceSegDur;

double parameter SliceSegDelay;
double parameter MinSliceSegDur;
double parameter SliceSegEndDelay;

double parameter
{
    display_name "Time for Movie";
    units "ms";
    format "%.2f";
    relations backbone;
    editable false;
}TimeForMovieFrames;

YesNo parameter
{
    display_name "RF Spoiling";
    relations backbone;
}RFSpoiling;

double parameter RFPhaseList[];

YesNo parameter
{
    display_name "Angio Mode";
    short_description "Inflow contrast is achieved by changing the loop structure in multi-slice experiments.";
    relations backbone;
}AngioMode;

/* ---------------------------------------------------------
 * remaining local method parameters
 * --------------------------------------------------------*/

double parameter
{
    display_name "Read Dephase Gradient";
    format "%f";
    units  "%";
    relations backbone;
}ReadDephGrad;

double parameter
{
    display_name "Max Read Dephase Gradient";
    format "%f";
    units  "%";
    relations backbone;
}ReadDephGradLim;

double parameter
{
    display_name "Read Gradient";
    format "%f";
    units "%";
    relations backbone;
}ReadGrad;

double parameter
{
    display_name "Max Read Gradient";
    format "%f";
    units "%";
    relations backbone;
}ReadGradLim;

double parameter
{
    display_name "2D Phase Gradient";
    format "%f";
    units "%";
    relations backbone;
}Phase2DGrad;

double parameter
{
    display_name "Max. 2D Phase Gradient";
    format "%f";
    units "%";
    relations backbone;
}Phase2DGradLim;


double parameter
{
    display_name "3D Phase Gradient";
    format "%f";
    units "%";
    relations backbone;
}Phase3DGrad;

double parameter
{
    display_name "Max. 3D Phase Gradient";
    format "%f";
    units "%";
    relations backbone;
}Phase3DGradLim;

double parameter
{
    display_name "Exc. Slice Gradient";
    format "%f";
    units  "%";
    relations backbone;
}ExcSliceGrad;


double parameter
{
    display_name "Max. Exc. Slice Gradient";
    format "%f";
    units  "%";
    relations backbone;
}ExcSliceGradLim;

double parameter
{
    display_name "Exc. Slice Reph. Gradient";
    format "%f";
    units  "%";
    relations backbone;
}ExcSliceRephGrad;

double parameter
{
    display_name "Max. Exc. Slice Reph. Gradient";
    format "%f";
    units  "%";
    relations backbone;
}ExcSliceRephGradLim;

double parameter
{
    display_name "Encoding Duration";
    short_description "Duration of encoding gradient.";
    relations backbone;
    units "ms";
    format "%.3f";
}EncGradDur;

double parameter {relations backbone;}Rew2DGrad;
double parameter {relations backbone;}Rew3DGrad;
double parameter {relations backbone;}RewGradDur;
double parameter TeFillDelay; /* placeholder, no relations */
double parameter TrFillDelay;

/* new parameters for SWI Reconstruction */
RecoMeth_MODE parameter
{
    display_name "Reconstruction Mode";
    short_description "Switches between standard and susceptibility-weighted imagig reconstruction.";
    relations SetRecoParam;
}RecoMethMode;

MASK_MODE parameter
{
    display_name "Weighting Mode";
    short_description "Selection of the way phase information will be used to influence image intensity.";
    relations MaskModeRel;
}WeightingMode;

double parameter
{
    display_name "Mask Weighting";
    short_description "Strength of the weighting mask function.";
    relations MaskWeightRange;
    format "%.2f";
}MaskWeighting;

double parameter
{
    display_name "Gauss Broadening";
    short_description "Defines the broadening (smoothing) effect of the Gauss filter.";
    relations GaussBroadRange;
    format "%.2f";
    units "mm";
}GaussBroadening;

/* KAM add for hpMR */

// Gradient system parameters
double parameter
{
    display_name "Gradient timing resolution";
    store true;
    units "us";
    format "%.6f";
    editable false;
}KAM_GradientTime;

double parameter
{
    display_name "Max Gradient Amplitude";
    store true;
    units "mT/m";
    format "%.2f";
    editable true;
    relations KAM_ReadModeRel;
}KAM_GradientMaxAmp;

double parameter
{
    display_name "Max Gradient Slew Rate";
    store true;
    units "T/m/s";
    format "%.2f";
    editable true;
    relations KAM_ReadModeRel;
}KAM_GradientMaxSlew;

double parameter
{
    display_name "Gyromagnetic Ratio";
    store true;
    units "MHz/T";
    format "%.4f";
    editable false;
}KAM_GammaMHzT;


// Parameters needed for create2dEpi (custom snapshot EPI readouts)
Epi_MODE parameter
{
    display_name "EPI Type";
    short_description "EPI waveform design type.";
    store true;
    relations KAM_ReadModeRel;
}KAM_EpiType;

int parameter
{
    display_name "EPI Matrix Size";
    short_description "2D Image Matrix Size (Square).";
    store true;
    relations KAM_ReadModeRel;
}KAM_EpiMatrix;

int parameter
{
    display_name "EPI Dwell Time Factor";
    short_description "Number of gradient points per acquisition point (Must be even integer).";
    store true;
    relations KAM_ReadModeRel;
}KAM_EpiDwellFactor;

double parameter
{
    display_name "EPI Partial Fourier Factor";
    short_description "Fraction of k-space sampled in blipped direction.";
    store true;
    relations KAM_ReadModeRel;
}KAM_EpiPartialFourierY;


// Parameters related to custom readout sizes
int parameter
{
    display_name "hpMR Acquisition Size";
    editable false;
    relations backbone;
}KAM_AcqMat[];

int parameter
{
    display_name "Number of Readout Points";
    editable true;
    relations backbone;
}KAM_NAcqPoints;

int parameter
{
    display_name "Number of Spectral-Spatial Slices";
    editable false;
}KAM_NSlices;

int parameter
{
    display_name "Number of Echoes";
    editable true;
    relations backbone;
}KAM_NEchoes;


// Parameters related to custom excitation
Exc_MODE parameter
{
    display_name "Excitation Mode";
    relations KAM_ExcModeRel;
}KAM_ExcMode;

char parameter
{
    display_name "RF Excitation Pulse File";
    store true;
    editable false;
}KAM_ExcShapeFile[64];

char parameter
{
    display_name "Excitation Gradient File";
    store true;
    editable false;
}KAM_ExcGradFile[256];

char parameter
{
    display_name "RF Excitation Pulse Metadata";
    store true;
    editable false;
}KAM_ExcMetaFile[256];

double parameter
{
    display_name "Excitation Gradient Shape";
    store true;
    editable false;
}KAM_ExcGradShape[];

int parameter
{
    display_name "Excitation Gradient Shape Size";
    store true;
    editable false;
}KAM_ExcGradSize;

double parameter
{
    display_name "RF Excitation Delay";
    short_description "Delay between start of excitation gradient and RF pulse.";
    units "us";
    format "%.2f";
    store true;
    editable false;
}KAM_ExcRfDelay;

double parameter
{
    display_name "Excitation Gradient Scale Factor";
    short_description "Scale factor applied to excitation gradient waveform.";
    format "%.3f";
    store true;
    editable false;
}KAM_ExcGradScale;


// Parameters related to custom preparation gradients
Prep_MODE parameter
{
    display_name "Preparation Mode";
    relations KAM_PrepModeRel;
}KAM_PrepMode;

char parameter
{
    display_name "Preparation Gradient File";
    store true;
    editable false;
}KAM_PrepShapeFile[256];

double parameter
{
    display_name "Preparation Gradient Shape X";
    store true;
    editable false;
}KAM_PrepShapeX[];

double parameter
{
    display_name "Preparation Gradient Shape Y";
    store true;
    editable false;
}KAM_PrepShapeY[];

double parameter
{
    display_name "Preparation Gradient Shape Z";
    store true;
    editable false;
}KAM_PrepShapeZ[];

int parameter
{
    display_name "Preparation Gradient Shape Size";
    short_description "Number of points in preparation gradient.";
    store true;
    editable false;
}KAM_PrepShapeSize;

double parameter
{
    display_name "Preparation Gradient Duration";
    units "ms";
    format "%.6f";
    store true;
    editable false;
}KAM_PrepLength;

int parameter
{
    display_name "Preparation Gradient On - Outer";
    short_description "Outer loop repetitions with preparation gradient waveform active.";
    store true;
    editable true;
}KAM_PrepGradOnOuter[];

int parameter
{
    display_name "Preparation Gradient On - Inner";
    short_description "Inner loop repetitions with preparation gradient waveform active.";
    store true;
    editable true;
}KAM_PrepGradOnInner[];


// Parameters related to custom readout
Read_MODE parameter
{
    display_name "Readout Mode";
    relations KAM_ReadModeRel;
}KAM_ReadMode;

char parameter
{
    display_name "Readout Gradient File";
    store true;
    editable false;
}KAM_ReadShapeFile[256];

double parameter
{
    display_name "Readout Gradient Shape X";
    store true;
    editable false;
}KAM_ReadShapeX[];

double parameter
{
    display_name "Readout Gradient Shape Y";
    store true;
    editable false;
}KAM_ReadShapeY[];

double parameter
{
    display_name "Readout Gradient Shape Z";
    store true;
    editable false;
}KAM_ReadShapeZ[];

int parameter
{
    display_name "Readout Gradient Shape Size";
    short_description "Number of points in readout gradient.";
    store true;
    editable false;
}KAM_ReadShapeSize;

double parameter
{
    display_name "Readout Gradient Duration";
    units "ms";
    format "%.6f";
    store true;
    editable false;
}KAM_ReadLength;

double parameter
{
    display_name "Acquisition Delay";
    short_description "Delay between start of readout gradient and first acquisition point.";
    units "us";
    format "%.2f";
    store true;
    editable false;
}KAM_ReadAcqDelay;

double parameter
{
    display_name "Readout X-Gradient Scale Factor";
    short_description "Scale factor applied to x-gradient of readout waveform.";
    format "%.3f";
    store true;
    editable false;
}KAM_ReadGradScaleX;

double parameter
{
    display_name "Readout Y-Gradient Scale Factor";
    short_description "Scale factor applied to y-gradient of readout waveform.";
    format "%.3f";
    store true;
    editable false;
}KAM_ReadGradScaleY;

double parameter
{
    display_name "Readout Z-Gradient Scale Factor";
    short_description "Scale factor applied to z-gradient of readout waveform.";
    format "%.3f";
    store true;
    editable false;
}KAM_ReadGradScaleZ;

double parameter
{
    display_name "Echo Delay";
    short_description "Delay between start of readout gradient and k-space center.";
    units "ms";
    format "%.3f";
    store true;
    editable false;
}KAM_ReadEchoDelay;


// Parameters related to IDEAL echo time shifting
double parameter
{
    display_name "Time Course Delay";
    short_description "Time between first excitations in subsequent repetitions.";
    relations backbone;
    units "ms";
    format "%.3f";
}KAM_IdealRepDelay;

YesNo parameter
{
    display_name "Minimize Time Course Delay";
    short_description "Minimize time between first excitations in subsequent repetitions.";
    relations backbone;
}KAM_IdealMinDelay;

double parameter
{
    display_name "Time Course Filling";
    short_description "Time course delay minus imaging TRs.";
    editable false;
    units "ms";
    format "%.3f";
}KAM_IdealRepFill;

double parameter
{
    display_name "Echo Time Shift";
    short_description "Delay added to subsequent echo times within a time course for IDEAL phase evolution.";
    relations backbone;
    units "ms";
    format "%.4f";
}KAM_IdealEchoShift;


// Parameters related to variable flip angles
YesNo parameter
{
    display_name "Vary Flip Angles";
    short_description "Use custom flip angle for each excitation in a repetition.";
    relations backbone;
}KAM_VFAOn;

YesNo parameter
{
    display_name "Use Effective Flip Angle";
    short_description "Calculate VFA scheme to provide specified effective flip angle with constant Mx (neglecting T1).";
    relations backbone;
}KAM_VFAEffectiveOn;

double parameter
{
    display_name "Effective Flip Angle";
    short_description "Flip angle imparted by overall VFA scheme with constant Mx (without T1 loss/recovery).";
    relations backbone;
    units "deg";
    format "%.2f";
}KAM_VFAEffectiveDegrees;

double parameter
{
    display_name "Variable Flip Angles";
    short_description "List of variable flip angles.";
    relations backbone;
    units "deg";
    format "%.2f";
}KAM_VFADegrees[];

double parameter
{
    display_name "Variable RF Powers";
    short_description "List of variable RF powers.";
    relations backbone;
    units "W";
    format "%.4f";
}KAM_VFAWatts[];


// Parameters related to fid excitation and readout
YesNo parameter
{
    display_name "Fid Echo";
    short_description "Acquire first echo of each time course without spatial encoding.";
    relations KAM_ExcModeRel;
}KAM_FidEcho;

PV_PULSE_LIST parameter
{
    display_name "Fid Excitation Pulse Shape";
    relations    KAM_FidExcPulseEnumRelation;
}KAM_FidExcPulseEnum;

PVM_RF_PULSE parameter
{
    display_name "Fid Excitation Pulse";
    relations    KAM_FidExcPulseRelation;
}KAM_FidExcPulse;

PVM_RF_PULSE_AMP_TYPE parameter
{
    display_name "Fid RF Pulse Amplitude";
    relations KAM_FidExcPulseAmplRel;
}KAM_FidExcPulseAmpl;

double parameter
{
    editable false;
}KAM_FidExcPulseShape[];


// Parameters related to presaturation excitation and readout
YesNo parameter
{
    display_name "Saturation Echo";
    short_description "Acquire first image of each time course with spectral saturation pulse.";
    relations KAM_ExcModeRel;
}KAM_SatEcho;

PV_PULSE_LIST parameter
{
    display_name "Saturation Pulse Shape";
    relations    KAM_SatExcPulseEnumRelation;
}KAM_SatExcPulseEnum;

PVM_RF_PULSE parameter
{
    display_name "Saturation Pulse";
    relations    KAM_SatExcPulseRelation;
}KAM_SatExcPulse;

PVM_RF_PULSE_AMP_TYPE parameter
{
    display_name "Sat RF Pulse Amplitude";
    relations KAM_SatExcPulseAmplRel;
}KAM_SatExcPulseAmpl;

double parameter
{
    editable false;
}KAM_SatExcPulseShape[];


// Parameters related to variable offset frequencies
YesNo parameter
{
    display_name "Spectrally Selective Excitations";
    short_description "Allows custom setting of transmit and receive frequency offsets.";
    relations KAM_ExcModeRel;
}KAM_SpecExc;

YesNo parameter
{
    display_name "Separate Tx and Rx Offsets";
    short_description "Specify transmit and receive frequency offsets separately.";
    relations KAM_ExcModeRel;
}KAM_SpecTxRxSeparate;

double parameter
{
    display_name "Variable Transmit Offset Frequencies";
    short_description "List of transmit frequency offsets.";
    relations backbone;
    units "Hz";
    format "%.1f";
}KAM_VOFTxHertz[];

double parameter
{
    display_name "Variable Receive Offset Frequencies";
    short_description "List of receive frequency offsets.";
    relations backbone;
    units "Hz";
    format "%.1f";
}KAM_VOFRxHertz[];

double parameter
{
    display_name "Saturation Offset Frequency";
    relations backbone;
    units "Hz";
    format "%.1f";
}KAM_SatFreq[];


// HyperSense trigger
YesNo parameter
{
    display_name "Trigger By HyperSense";
    short_description "Triggers start of acquisition from HyperSense dissolution. CHECK TRIGGER CABLE!";
    relations backbone;
}KAM_HSTrig;


/****************************************************************/
/*	E N D   O F   F I L E					*/
/****************************************************************/

