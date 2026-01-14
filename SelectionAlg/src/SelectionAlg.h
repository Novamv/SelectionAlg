#ifndef SELECTIONALG_H
#define SELECTIONALG_H

#include "SniperKernel/AlgBase.h"
#include "SniperKernel/ToolBase.h"
#include "EvtNavigator/NavBuffer.h"

#include <vector>
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TTimeStamp.h"

#include "SpmtElecConfigSvc/SpmtElecConfigSvc.h"
#include "OECTagSvc/OECTagSvc.h"
#include "OECTagID/OECTagID.h"

#include "EventTagSvc/EventTagSvc.h"
#include "MuonClassificationTool/MuonClassificationTool.h"
#include "SpalNeutronSelectionTool/SpalNeutronSelectionTool.h"
#include "IBDSelectionTool/IBDSelectionTool.h"

#include "TOF.h"

#include <cstdint>
#include <string>
#include <map>

class RecGeomSvc; 
class CdGeom;

enum class SkipReason {
  None,
  StartOfFile,
  BigGap,
  MissingHeader
};


class SelectionAlg : public AlgBase
{
private :
  
  double ComputeLTOF(double pmtid, double evtx, double evty, double evtz);
  double ComputeSTOF(double pmtid, double evtx, double evty, double evtz);
  
public :
  
  SelectionAlg(const std::string& name);
  
  
  bool initialize();
  bool execute();
  bool finalize();
  bool Book_tree();
  void clearAllTrees();	
  
private :

  //Alg Property
  bool saveElec, saveSim, saveHitInfo;
  int interface;
  std::string recEDMPath, filename;

  JM::NavBuffer* m_buf;

  // Tools
  SpmtElecConfigSvc* m_spmtSvc;
  OECTagSvc* m_tagsvc;
  EventTagSvc* m_eventTagsvc;
  MuonClassificationTool* m_MuClassifier;
  SpalNeutronSelectionTool* m_NeutronClassifier;
  IBDSelectionTool* m_IBDClassifier;

  //OEC Tag list
  uint32_t i_pIBD; 
  uint32_t i_dIBD;


  std::vector<std::string> m_classifiernames;

  // Some useful variables
  TTimeStamp FirstTime, theTime, PreviousTime, tLastMuon;
  TTimeStamp prevCDTime;
  TTimeStamp prevWPTime;
  int m_DelayEvt;

  // skipping
  SkipReason skipReason = SkipReason::None;
  TTimeStamp skipStartTime;
  double dt_skip;
  double dtCD = 0.0;
  double dtWP = 0.0;

  // ------------------ TOF ------------------
  TOFCalculator m_TOF;
  
  // ------------------ PMT POS ------------------
  unsigned int TotalLPMT = 17612;
  unsigned int TotalSPMT = 25600;
  
  std::vector<TVector3> ALL_LPMT_pos;
  std::vector<TVector3> ALL_SPMT_pos;
 

  // ------------------ OUTPUT TREE ------------------

  int m_iEvt;
  TString m_Tag, m_PrevTag;

  uint32_t m_EvtID;
  uint64_t m_AssembleID;
  int m_iRun;
  long long m_Trigger;
  uint64_t m_TimeStamp;

  std::vector<std::string> m_TriggerName;

  double m_ChargeTotLPMT;
  double m_ChargeTotUpLPMT;
  double m_ChargeTotSPMT;
  double m_ChargeTotWP;
  unsigned long m_TimeStampInNanoSec; 

  TTree *m_ntuple1; // calibration tree
  int m_NbHitLPMTCalib;
  int m_NbHitSPMTCalib;
  int m_NbHitWPCalib;
  double m_HitTime_std;
  double m_HitTime_mean;
  std::vector<int> m_PmtIdCalib;
  std::vector<double> m_HitTimeCalib;
  std::vector<double> m_HitTimeCalibTOF;
  std::vector<double> m_ChargeCalib;
  


  TTree *m_ntuple2; // General info of the event
  double runtime = 0.0;
  double effruntime = 0.0;
  int nIBD = 0;
  int nCDMuons = 0;
  int nWPMuons = 0;
  int nCDWPMuons = 0;
  int nNeutrons = 0;
  TString m_fname;
  // std::vector<double> MuonProcessTime;
  // std::vector<double> NeutronProcessTime;
  // std::vector<double> IBDProcessTime;


  TTree *m_ntuple3; // reco tree
  float m_TotalPERec;
  float m_NFiredPMT;
  float m_RecE;
  float m_RecX;
  float m_RecY;
  float m_RecZ;
  float m_T0;
  float m_PosQuality;
  float m_EnergyQuality;
  float m_chisq;
  int m_nTrack;
  std::vector<float> m_MuPE;
  std::vector<float> m_MuQuality;
  std::vector<float> m_MuEntryY;
  std::vector<float> m_MuEntryX;
  std::vector<float> m_MuEntryZ;
  std::vector<float> m_MuEntryTheta;
  std::vector<float> m_MuEntryPhi;
  std::vector<float> m_MuExitX;
  std::vector<float> m_MuExitY;
  std::vector<float> m_MuExitZ;
  // std::vector<float> m_MuDX;
  // std::vector<float> m_MuDY;
  // std::vector<float> m_MuDZ;

  //Oec related
  int m_myOECtag;
  TString m_OECtag;


  double PMT_R;
  double LS_R;
  
  double RfrIndxLS;
  double RfrIndxWR;
  double c;
};

#endif
