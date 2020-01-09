/* /opt/PV6.0.1/prog/curdir/kamichel/ParaVision/methods/src/hpMR/initMeth.c */
void initMeth(void);
/* /opt/PV6.0.1/prog/curdir/kamichel/ParaVision/methods/src/hpMR/loadMeth.c */
void loadMeth(const char *);
/* /opt/PV6.0.1/prog/curdir/kamichel/ParaVision/methods/src/hpMR/parsRelations.c */
void ExcPulse1EnumRelation(void);
void ExcPulse1AmplRel(void);
void HandleRFPulseAmplitude(void);
void ExcPulse1Relation(void);
void ExcPulse1Range(void);
void NAveragesRange(void);
void NAveragesRels(void);
void EffSWhRange(void);
void EffSWhRel(void);
void InversionTimeRels(void);
void InversionTimeRange(void);
void ReadSpoilerRel(void);
void SliceSpoilerRel(void);
void RecoMethModeVisPar(void);
void MaskModeRel(void);
void GaussBroadRange(void);
void MaskWeightRange(void);
void MyRgInitRel(void);
void KAM_FidExcPulseEnumRelation(void);
void KAM_FidExcPulseAmplRel(void);
void KAM_FidHandleRFPulseAmplitude(void);
void KAM_FidExcPulseRelation(void);
void KAM_FidExcPulseRange(void);
void KAM_SatExcPulseEnumRelation(void);
void KAM_SatExcPulseAmplRel(void);
void KAM_SatHandleRFPulseAmplitude(void);
void KAM_SatExcPulseRelation(void);
void KAM_SatExcPulseRange(void);
void KAM_ExcModeRel(void);
void KAM_PrepModeRel(void);
void KAM_ReadModeRel(void);
void KAM_GradCheck(double *, int, double, double, double *, double *);
/* /opt/PV6.0.1/prog/curdir/kamichel/ParaVision/methods/src/hpMR/BaseLevelRelations.c */
void SetBaseLevelParam(void);
void SetBasicParameters(void);
void SetFrequencyParameters(void);
void SetGradientParameters(void);
void SetInfoParameters(void);
void SetMachineParameters(void);
void SetPpgParameters(void);
void SetACQ_obj_orderForMovie(void);
/* /opt/PV6.0.1/prog/curdir/kamichel/ParaVision/methods/src/hpMR/RecoRelations.c */
void SetRecoParam(void);
void RecoDerive(void);
void SWIRecoRel(void);
/* /opt/PV6.0.1/prog/curdir/kamichel/ParaVision/methods/src/hpMR/backbone.c */
void backbone(void);
void UpdateSequenceTiming(void);
void UpdateRepetitionTime(void);
void UpdateTotalTime(void);
void UpdateGeometryMinima(double *, double *);
void UpdateReadSliceGradients(void);
void UpdatePhaseGradients(void);
void UpdateFrequencyOffsets(void);
void UpdateRFPulses(void);
void ControlGradientLimits(YesNo);
void UpdateMovie(void);
void UpdateEchoTime(void);
double minLoopRepetitionTime(void);
void UpdateRFSpoiling(void);
void HandleVFA(void);
void EnforceHpmrDesignParams(void);
/* /opt/PV6.0.1/prog/curdir/kamichel/ParaVision/methods/src/hpMR/deriveVisu.c */
void deriveVisu(void);
