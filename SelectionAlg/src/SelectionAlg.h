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

struct PendingIBD {
  int     promptEntry;      // m_iEvt at the time isPrompt() fired
  int     delayEntryOffset; // relative offset returned by getDelayOffset()

  float   pEnergy;
  float   pX, pY, pZ;
  float   pNPE;
  uint64_t pTimeStamp;      // ns
  std::vector<double> pHitTimeTOF;

  int NeutronVeto;
}

struct IBDPair {
  // --- Coincidence ---
  double  dt_ns;      // delay - prompt time in ns 
  float   dR_mm;      // |pVertex - dVertex| in mm

  // --- Prompt ---
  float   pEnergy;
  float   pX, pY, pZ;
  float   pNPE;
  uint64_t pTimeStamp;
  std::vector<double> pHitTimeTOF;

  // --- Delay ---
  float   dEnergy;
  float   dX, dY, dZ;
  float   dNPE;
  uint64_t dTimeStamp;
  std::vector<double> dHitTimeTOF;


  // --- Flags ---
  int     neutronVeto;    // 1 if prompt is within spallation-neutron sphere
  int     runNumber;
  int     promptEntry;
  int     delayEntry;
}

class SelectionAlg : public AlgBase
{
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
  uint32_t i_pBiPo214pair;
  uint32_t i_dBiPo214pair;


  std::vector<std::string> m_classifiernames;

  // Some useful variables
  TTimeStamp FirstTime, theTime, PreviousTime, tLastMuon;
  TTimeStamp prevCDTime;
  TTimeStamp prevWPTime;
  int m_DelayEvt;
  std::vector<PendingIBD> m_pendingIBD;

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


  TVector3 NeutronVertex; 
  TTimeStamp NeutronTime = NULL;

  TTree *m_ntuple1; // event tree
  int m_NeutronVeto;
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


  TTree *m_ntuple3; // Mu reco tree
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

  TTree *m_ibdtree; 
  IBDPair m_pair;

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
