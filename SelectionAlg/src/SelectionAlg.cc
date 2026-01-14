#include "SelectionAlg.h"
#include "TOF.h"

#include "EvtNavigator/NavBuffer.h"
#include "EvtNavigator/EvtNavHelper.h"
#include "SniperKernel/AlgFactory.h"
#include "SniperKernel/ToolFactory.h"

#include "Event/CdVertexRecHeader.h"
#include "Event/CdLpmtCalibHeader.h"
#include "Event/CdSpmtCalibHeader.h"
#include "Event/CdLpmtElecTruthHeader.h"
#include "Event/CdSpmtElecTruthHeader.h"
#include "Event/CdWaveformHeader.h"
#include "Event/WpCalibHeader.h"
#include "Event/WpCalibEvt.h"

#include "Event/SimHeader.h"

#include "Event/CdLpmtElecEvt.h"
#include "Event/CdSpmtElecEvt.h"
#include "Event/CdLpmtElecHeader.h"
#include "Event/CdSpmtElecHeader.h"

#include "Event/CdTriggerHeader.h"
#include "Event/CdTriggerEvt.h"
#include "Event/WpTriggerHeader.h"
#include "Event/WpTriggerEvt.h"
#include "Event/TtTriggerHeader.h"
#include "Event/MMTriggerHeader.h"

#include "Event/CdVertexRecHeader.h"
#include "Event/CdTrackRecHeader.h"

#include "Identifier/Identifier.h"
#include "Identifier/IDService.h"
#include "Identifier/CdID.h"

#include "RootWriter/RootWriter.h"
#include "BufferMemMgr/IDataMemMgr.h"

#include "SpmtElecConfigSvc/SpmtElecConfigSvc.h"

#include "OECTagSvc/OECTagSvc.h"
#include "OECTagID/OECTagID.h"

#include "Geometry/IPMTParamSvc.h"

#include "Event/OecHeader.h"
#include "Event/OecEvt.h"

#include "CLHEP/Vector/LorentzVector.h"

#include "TTree.h"
//#include "Geometry/IRecGeomSvc.hh"
#include <fstream>

DECLARE_ALGORITHM(SelectionAlg);

SelectionAlg::SelectionAlg(const std::string& name)
: AlgBase(name),
  m_iEvt(-1),
  m_buf(0),
  m_spmtSvc(0),
  m_tagsvc(0)
{
	declProp("enableElec", saveElec);
	declProp("interface", interface);
	declProp("enableSim", saveSim);
	declProp("recEDMPath", recEDMPath);
	declProp("saveHitInfo", saveHitInfo);
	declProp("ClassifierNames", m_classifiernames);
	declProp("Filename", filename);

	PMT_R = 35.4; //m
	LS_R = 17.7; //m
	RfrIndxLS = 1.5;
	RfrIndxWR = 1.355;

	c = 299792458.0; //m/s
}

bool SelectionAlg::initialize()
{
	//----------------------------------------------------------------------------
	const std::string compilation_date = __DATE__;
	const std::string compilation_time = __TIME__;
	std::cout <<"##################################################################"<<std::endl
		<<"The source file was compiled on " << compilation_date<< " at " << compilation_time <<std::endl
		<<"##################################################################"<<std::endl;
	//----------------------------------------------------------------------------

	// =======================================================================
	// Loading PMT positions
	// =======================================================================

	SniperPtr<IPMTParamSvc> pmtsvc(getParent(), "PMTParamSvc");

	if (ALL_LPMT_pos.size()==0 && pmtsvc.valid()) {
        TotalLPMT = pmtsvc->get_NTotal_CD_LPMT();

        std::cout << " PMT Information " << std::endl;

		std::cout << "LPMT" <<std::endl;
        for (unsigned int ith = 0; ith < TotalLPMT; ith++)
        {
            TVector3 all_pmtCenter(pmtsvc->getPMTX(ith), pmtsvc->getPMTY(ith), pmtsvc->getPMTZ(ith));
            ALL_LPMT_pos.push_back(all_pmtCenter);
        }
	}

	if (ALL_SPMT_pos.size()==0 && pmtsvc.valid()) {
		TotalSPMT = pmtsvc->get_NTotal_CD_SPMT();

		std::cout << "SPMT" <<std::endl;
        for (unsigned int ith = 0; ith < TotalSPMT; ith++)
        {
            TVector3 all_pmtCenter(pmtsvc->getPMTX(ith+20000), pmtsvc->getPMTY(ith+20000), pmtsvc->getPMTZ(ith+20000));
            ALL_SPMT_pos.push_back(all_pmtCenter);
        }
	}

	// =======================================================================
    // GET EVENT
    // =======================================================================


	LogDebug << "initializing" << std::endl;
	std::cout<<"36"<<std::endl;
	// gDirectory->pwd();

	SniperDataPtr<JM::NavBuffer> navBuf(getParent()/*getRoot()*/,"/Event");
	std::cout<<"Parent="<<getParent()<<"Root="<<getRoot()<<std::endl;//" Parent Name"<<getParentName()<<std::endl;
	if ( navBuf.invalid() ) {
	    LogError << "cannot get the NavBuffer @ /Event" << std::endl;
	    return false;
	}
	m_buf = navBuf.data();
	// gDirectory->pwd();
	SniperPtr<SpmtElecConfigSvc> svc(*getRoot(), "SpmtElecConfigSvc");
	if (svc.invalid()) {
	    LogError << "can't find service SpmtElecConfigSvc" << std::endl;
	    return false;
	}
	m_spmtSvc = svc.data();

	//get EventTag service
	SniperPtr<EventTagSvc> eventTagsvc(getParent(),"EventTagSvc");
	if( eventTagsvc.invalid()) {
	    LogError << "Unable to locate EventTagSvc" << std::endl;
	    return false;
	}
	m_eventTagsvc = eventTagsvc.data();

	
	// ================================================
	// =============  Classifiers  =============
	// ================================================
	
		// Muon
	m_MuClassifier = tool<MuonClassificationTool>("MuonClassificationTool");
	if(m_MuClassifier==NULL){
		LogError << "Failed to retrieve MuonClassificationTool" << std::endl;
		return false;
	}
	else if(!m_MuClassifier->initialize()){
		LogError << "Failed to initilaise MuonClassificationTool" << std::endl;
		return false;
	}

		// Neutron
	m_NeutronClassifier = tool<SpalNeutronSelectionTool>("SpalNeutronSelectionTool");
	if(m_NeutronClassifier==NULL){
		LogError << "Failed to retrieve SpalNeutronTool" << std::endl;
		return false;
	}
	else if(!m_NeutronClassifier->initialize()){
		LogError << "Failed to initilaise SpalNeutronTool" << std::endl;
		return false;
	}
		// IBD
	m_IBDClassifier = tool<IBDSelectionTool>("IBDSelectionTool");
	if(m_IBDClassifier==NULL){
		LogError << "Failed to retrieve IBDSelectionTool" << std::endl;
		return false;
	}
	else if(!m_IBDClassifier->initialize()){
		LogError << "Failed to initilaise IBDSelectionTool" << std::endl;
		return false;
	}


	// ================================================
	// =============  Oec Tag Service  =============
	// ================================================


	SniperPtr<OECTagSvc> tagsvc(getParent(),"OECTagSvc");
	if( tagsvc.invalid()) {
	    LogError << "Unable to locate tagsvc" << std::endl;
	    return false;
	}
	m_tagsvc = tagsvc.data();

	i_pIBD = m_tagsvc->getpTag("InverseBetaDecay");
	i_dIBD = m_tagsvc->getdTag("InverseBetaDecay");



	m_fname = filename;

	Book_tree();
  	return true;
}


bool SelectionAlg::finalize()
{
	LogDebug << "finalizing" << std::endl;

	LogInfo << "Filling RunInfo tree" << std::endl;
	m_ntuple2->Fill();


	m_MuClassifier->finalize();
	m_NeutronClassifier->finalize();
	m_IBDClassifier->finalize();

	return true;
}


bool SelectionAlg::execute()
{
	LogInfo << "executing: " << m_iEvt++ << std::endl;

	// gDirectory->pwd();

//      JM::EvtNavigator* navig = 0;
	JM::SimEvt* simevent = 0;
	JM::CdVertexRecEvt* recevent = 0;
	JM::CdLpmtCalibEvt* calibeventLPMT = 0;
	JM::CdSpmtCalibEvt* calibeventSPMT = 0;
	JM::WpCalibEvt* calibeventWP = 0;
	JM::CdWaveformEvt* elecevent = 0;
	JM::CdLpmtElecTruthEvt *trutheventLPMT = 0;
	JM::CdSpmtElecTruthEvt *trutheventSPMT = 0;
	JM::CdLpmtElecEvt *eventLPMT = 0;
	JM::CdSpmtElecEvt *eventSPMT = 0;
	JM::CdTriggerEvt *triggerevent =0;
	JM::WpTriggerEvt *Wptriggerevent =0;
	JM::OecEvt *oecevt = 0;
	JM::CdTrackRecEvt *trackrecevt = 0;
	// Get the events of different stages
	// calculate block charge and time
	
	LogInfo << "Buffer Size " << m_buf->size() << std::endl;
	auto nav = m_buf->curEvt();
	
	const auto& paths = nav->getPath();
	const auto& refs = nav->getRef();
	
	LogInfo << "Detector type is  " <<nav->getDetectorType()<<std::endl;
	LogInfo << "Start to Explore SmartRef: " << std::endl;
	LogInfo << "Size of paths: " << paths.size() << std::endl;
	LogInfo << "Size of refs: " << refs.size() << std::endl;
	
	for (size_t i = 0; i < paths.size(); ++i) {
		LogInfo << refs[i]<<" -> ref: " << std::endl;
		const std::string& path = paths[i];
		JM::SmartRef* ref = refs[i];
		JM::EventObject* evtobj = ref->GetObject();
		
		LogInfo << " path: " << path
		<< " ref->entry(): " << ref->entry()
		<< " evtobj: " << evtobj;
		
		if (path=="/Event/Sim") {
			auto hdr = dynamic_cast<JM::SimHeader*>(evtobj);
			LogInfo <<i<<" SimHeader: " << hdr;
		}
		LogInfo << std::endl;
	}

	clearAllTrees(); //Clearing all trees variables
	
// ================================================
// =============  Inititalization  =============
// ================================================

	theTime = nav->TimeStamp();

	if(m_iEvt == 0) {
		PreviousTime = theTime;
		skipReason = SkipReason::StartOfFile;
		skipStartTime = theTime;
		FirstTime = theTime;
		dtCD = 0;
		dtWP = 0;
	}

	// =============  Fast skip  =============
	if(skipReason != SkipReason::None){
		
		dt_skip = (theTime.GetSec() - skipStartTime.GetSec()) * 1000000000ULL + (theTime.GetNanoSec() - skipStartTime.GetNanoSec());

		bool skip1 = (skipReason == SkipReason::StartOfFile || skipReason == SkipReason::BigGap) && dt_skip < 1.2e9; //1.2 s skip if Start of file of big gap
		bool skip2 = skipReason == SkipReason::MissingHeader && dt_skip < 5e6; // 5 ms skip if no headers

		if(skip1 || skip2){
			PreviousTime = theTime;
			return true;
		}

		skipReason = SkipReason::None;
	}

	// =============  Compute dt  =============

	double globalTime = (theTime.GetSec() - FirstTime.GetSec())*1000000000ULL + (theTime.GetNanoSec() - FirstTime.GetNanoSec());

	uint64_t dt = (theTime.GetSec() - PreviousTime.GetSec())*1000000000ULL + (theTime.GetNanoSec() - PreviousTime.GetNanoSec());
	double dtLastMuon = (theTime.GetSec() - tLastMuon.GetSec())*1000000000ULL + (theTime.GetNanoSec() - tLastMuon.GetNanoSec());

	LogInfo << "Global time: " << globalTime << std::endl;


	
	// =============  Check	for missing headers  =============
	
	auto calibheaderWP = JM::getHeaderObject<JM::WpCalibHeader>(nav);
	auto calibheaderLPMT = JM::getHeaderObject<JM::CdLpmtCalibHeader>(nav);
	auto recheader = JM::getHeaderObject<JM::CdVertexRecHeader>(nav, recEDMPath);
	
	if(!recheader && !calibheaderLPMT && !calibheaderWP){
		skipReason = SkipReason::MissingHeader;
		skipStartTime = theTime;
		PreviousTime = theTime;
		LogInfo << "New skip reason: Missing Header" << std::endl;
		return true;
	}
	
	// =============  Check	for big gaps  =============
	
	if(calibheaderLPMT){
		if(dtCD == 0) dtCD = (theTime.GetSec() - PreviousTime.GetSec())*1000000000ULL + (theTime.GetNanoSec() - PreviousTime.GetNanoSec());
		else dtCD = (theTime.GetSec() - prevCDTime.GetSec())*1000000000ULL + (theTime.GetNanoSec() - prevCDTime.GetNanoSec());
		
		prevCDTime = theTime;
		LogInfo << "dtCD: " << dtCD << std::endl;

		if(dtCD > 50e6) {
			skipReason = SkipReason::BigGap;
			skipStartTime = theTime;
			LogInfo << "New skip reason: CD gap > 50ms" << std::endl;
			return true;
		}
	}
	if(calibheaderWP){
		if(dtWP == 0) dtWP = (theTime.GetSec() - PreviousTime.GetSec())*1000000000ULL + (theTime.GetNanoSec() - PreviousTime.GetNanoSec());
		else dtWP = (theTime.GetSec() - prevWPTime.GetSec())*1000000000ULL + (theTime.GetNanoSec() - prevWPTime.GetNanoSec());

		prevWPTime = theTime;
		LogInfo << "dtWP: " << dtWP << std::endl;
		
		if(dtWP > 70e6) {
			skipReason = SkipReason::BigGap;
			skipStartTime = theTime;
			LogInfo << "New skip reason: WP gap > 50ms" << std::endl;
			return true;
		}
	}
	
	// =============  Runtime counting  =============
	
	runtime += dt;
	if(dtLastMuon > 5e6){
		effruntime += dt;
	}
	
	PreviousTime = theTime;
	
	LogInfo << "Run Time: " << runtime << std::endl;
	LogInfo << "Effective Run Time: " << effruntime << std::endl;
	
	
// ================================================
// =============  Run Classification  =============
// ================================================

	LogInfo << "dtLastMuon: " << dtLastMuon << std::endl;

	if(m_iEvt == m_DelayEvt) m_Tag = m_eventTagsvc->getTag(nav);
	else if(dtLastMuon > 50e3 && m_MuClassifier->isMuon(nav)){ //50 us dead time
		tLastMuon = theTime;
		m_Tag = m_eventTagsvc->getTag(nav);
		if(m_Tag == "CDMuon") nCDMuons++;
		if(m_Tag == "CDWPMuon") nCDWPMuons++;
		if(m_Tag == "WPMuon") nWPMuons++;
	}
	else if((dtLastMuon > 20e3 && dtLastMuon < 2e6) && m_NeutronClassifier->isSpalNeutron(nav)){
		m_Tag = m_eventTagsvc->getTag(nav);
		nNeutrons++;
	}
	else if(dtLastMuon > 5e6 && m_IBDClassifier->isPrompt(nav)) {
		m_Tag = m_eventTagsvc->getTag(nav);
		m_DelayEvt = m_iEvt + m_IBDClassifier->getDelayOffset();
		nIBD++;
	}
		

	LogInfo << "Current Tag: " << m_Tag << std::endl;
	if(m_Tag == "") return true;

	// nav = m_buf->curEvt(); // in case classification messes up the nav (shouldn't)

// ================================================
// ========  Fill Trees for tagged events  ========
// ================================================

	// ---------------------- Load EDM classes ----------------------
	

	auto simheader = JM::getHeaderObject<JM::SimHeader>(nav);
	if(simheader){
		simevent = (JM::SimEvt*)simheader->event();
		LogInfo << "SimEvent Read in: " << simevent << std::endl;
		LogInfo << "SimEvent Track: " << simevent->getTracksVec().size() << std::endl;
		LogInfo << "SimEvent Hits: " << simevent->getCDHitsVec().size() << std::endl;
	}

	// auto recheader = JM::getHeaderObject<JM::CdVertexRecHeader>(nav, recEDMPath);
	if (recheader) {
	  recevent = recheader->event();
	  LogInfo << "RecEvent Read in: " << recevent << std::endl;
	}
	auto trackheader = JM::getHeaderObject<JM::CdTrackRecHeader>(nav, "/Event/CdTrackRecClassify");
	if(trackheader){
		trackrecevt = trackheader->event();
		LogInfo << "TrackRecEvt Read in: " << trackrecevt << std::endl;
	}
	// auto calibheaderLPMT = JM::getHeaderObject<JM::CdLpmtCalibHeader>(nav);
	if (calibheaderLPMT) {
	  calibeventLPMT = calibheaderLPMT->event();
	  LogInfo << "CalibEventLPMT Read in: " << calibeventLPMT << std::endl;
	}
	auto calibheaderSPMT = JM::getHeaderObject<JM::CdSpmtCalibHeader>(nav);
	if (calibheaderSPMT) {
	  calibeventSPMT = calibheaderSPMT->event();
	  LogInfo << "CalibEventSPMT Read in: " << calibeventSPMT << std::endl;
	}
	// auto calibheaderWP = JM::getHeaderObject<JM::WpCalibHeader>(nav);
	if (calibheaderWP) {
	  calibeventWP = calibheaderWP->event();
	  LogInfo << "CalibEventWP Read in: " << calibeventWP << std::endl;
	}

	auto triggerheader = JM::getHeaderObject<JM::CdTriggerHeader>(nav);
	if(triggerheader){triggerevent = triggerheader->event();
	  LogInfo <<"CD TriggerEvent Read in: " << triggerevent <<std::endl;
	}

	auto OEChdr = JM::getHeaderObject<JM::OecHeader>(nav);
	if(OEChdr){
		oecevt = dynamic_cast<JM::OecEvt*>(OEChdr->event("JM::OecEvt"));
		LogInfo <<"OEC Event Read in: " << oecevt <<std::endl;
	}

	// ---------------------- Set all trees variables ----------------------

	if(oecevt) {
		m_myOECtag = oecevt->getTag();
		if((m_myOECtag & i_pIBD)== i_pIBD) m_OECtag = "Prompt";
        if((m_myOECtag & i_dIBD)== i_dIBD) m_OECtag="Delay";
		//add BiPo tag
	}

	if(triggerevent)
	{
	    const auto& type = triggerevent->triggerType();
	    const auto& pmtFired = triggerevent->nHitMultiplicity();
	    const auto& trigTime = triggerevent->triggerTime();
	    m_Trigger = trigTime.GetNanoSec();

		double dtrig = (theTime.GetSec() - trigTime.GetSec())*1000000000ULL + (theTime.GetNanoSec() - trigTime.GetNanoSec());

		if(type.size() != 0){
			for(auto it = 0; it<type.size(); it++){
				LogInfo<<"Trigger type = "<<type[it]<<std::endl;
				m_TriggerName.push_back(type[it]);
			}
		}
	}

	const auto& timestamp = nav->TimeStamp();
	int RunNumber = nav->RunID();
	uint32_t EventNumber = nav->EventID();
	m_iRun = RunNumber;
	m_EvtID = EventNumber;
	m_AssembleID = nav->AssembleID();
	unsigned long long tempTS =  timestamp.GetSec()*1e9 + timestamp.GetNanoSec() ;
	uint64_t TimeRef_2 =  timestamp.GetSec()*1000000000ULL + timestamp.GetNanoSec() ;
	unsigned long long TimeStampLPMT = (tempTS&0xFFFFFFFF00000000) ;

	m_TimeStamp = TimeRef_2;

	double ChargeTot=0;
	std::vector<int> tempPmtIds;
	std::vector<double> tempHitTimes;
	std::vector<double> tempCharges;

	if(calibeventLPMT){

		// TH1F* hTime = new TH1F("hTime", "hTime", 100, 0, 1000);
		const auto& chhlistLPMT = calibheaderLPMT->event()->calibPMTCol();
		float sum = 0, sum2 = 0;
		int n = 0;
		for (auto chit = chhlistLPMT.begin(); chit!=chhlistLPMT.end(); ++chit){
			auto calib = *chit;
			unsigned int pmtId = calib->pmtId();
			Identifier id = Identifier(pmtId);
			int TruePM=CdID::module(id);
			m_NbHitLPMTCalib+=calib->size();

			for(unsigned int j=0;j<calib->size();j++){
				tempPmtIds.push_back(TruePM);
				tempCharges.push_back(calib->charge(j));

				double time = calib->time(j);
				// hTime->Fill(time);
				sum += time;
				sum2 += time*time;
				++n;

				tempHitTimes.push_back(time);
				ChargeTot+=calib->charge(j);
			}
		}
		m_ChargeTotLPMT = ChargeTot;
		m_HitTime_mean = sum/n;
		m_HitTime_std = std::sqrt((sum2 - sum*sum/n) / (n-1));

		if(m_Tag == "Neutron" && m_HitTime_std < 275){
			m_Tag = "SpalNeutron";
		}
		// m_HitTime_mean = hTime->GetMean();
		// m_HitTime_std = hTime->GetRMS();

		if(m_Tag == "Prompt" || m_Tag == "Delay"){
			m_PmtIdCalib.insert(m_PmtIdCalib.end(), tempPmtIds.begin(), tempPmtIds.end());
			m_HitTimeCalib.insert(m_HitTimeCalib.end(), tempHitTimes.begin(), tempHitTimes.end());
			m_ChargeCalib.insert(m_ChargeCalib.end(), tempCharges.begin(), tempCharges.end());
		}

		// delete hTime;
	}


	ChargeTot=0;
	if(calibeventSPMT)
	{
		const auto& chhlistSPMT = calibheaderSPMT->event()->calibPMTCol();
		for (auto cchit = chhlistSPMT.begin(); cchit!=chhlistSPMT.end(); ++cchit) {
			auto calibSPMT = *cchit;

			unsigned int pmtId = calibSPMT->pmtId();
			Identifier id = Identifier(pmtId);
			int TruePM=CdID::module(id);
			int CircleNo = CdID::circleNo(id);
			TruePM+=20000-17612;

			m_NbHitSPMTCalib+=calibSPMT->size();
			for(unsigned int j=0;j<calibSPMT->size();j++)
			{
				m_PmtIdCalib.push_back(TruePM);
				double timeTOF2 = 0.0;
				//timeTOF2 = calibSPMT->time(j) - ComputeSTOF(TruePM, VtxReco, VtyReco, VtzReco);
				timeTOF2 = calibSPMT->time(j);
				m_HitTimeCalib.push_back(timeTOF2);
				//m_HitTimeCalib.push_back(calibSPMT->time(j));
				double CalibSPMTBackADC = 124.*0.48*calibSPMT->charge(j)+76.;
				// m_ChargeCalib.push_back( CalibSPMTBackADC);
				// ChargeTot+=CalibSPMTBackADC;

				m_ChargeCalib.push_back(calibSPMT->charge(j));
				ChargeTot+=calibSPMT->charge(j);
			}
		}
	}
	m_ChargeTotSPMT = ChargeTot;

	ChargeTot=0;
	if(calibeventWP)
	{
	    const auto& chhlistWP = calibheaderWP->event()->calibPMTCol();
	    for (auto wphit = chhlistWP.begin(); wphit!=chhlistWP.end(); ++wphit) {
			auto calibWP = *wphit;
			m_NbHitWPCalib+=calibWP->size();
			for(unsigned int j=0;j<calibWP->size();j++){
				ChargeTot+=calibWP->charge(j);
		 	}
	    }
	}
	 m_ChargeTotWP = ChargeTot;

	if (recevent)
	{
		const auto& recvertices = recevent->vertices();
		LogInfo << " CdVertexRecEvt: " << std::endl;
		LogInfo << " - number of vertices: " << recevent->nVertices() << std::endl;

		for(auto vertex: recvertices)
		{
			m_TotalPERec = vertex->peSum();
			m_NFiredPMT = vertex->nfiredpmts();
			m_RecE = vertex->energy();
			m_RecX = vertex->x();
			m_RecY = vertex->y();
			m_RecZ = vertex->z();
			m_T0 = vertex->t0();
			m_PosQuality = vertex->positionQuality();
			m_EnergyQuality = vertex->energyQuality();
			m_chisq = vertex->chisq();
		}
		
		TVector3 vertex(m_RecX, m_RecY, m_RecZ);
		for (size_t i = 0; i < m_PmtIdCalib.size(); i++)
		{
			int pmtid = m_PmtIdCalib.at(i);
			m_TOF.setVertex(vertex);
			m_TOF.setInterface(interface);
			if(pmtid < 17612){
				// TOFCalculator TOF(vertex, ALL_LPMT_pos.at(pmtid), interface); //interface is set to -40 m by default
				m_TOF.setPMTVertex(ALL_LPMT_pos.at(pmtid));
				double timeTOF = m_HitTimeCalib.at(i) - m_TOF.CalTOF();
				m_HitTimeCalibTOF.push_back(timeTOF);
			}
			else if (pmtid >= 20000){
				m_TOF.setPMTVertex(ALL_SPMT_pos.at(pmtid - 20000));
				// TOFCalculator TOF(vertex, ALL_SPMT_pos.at(pmtid - 20000), interface);
				double timeTOF = m_HitTimeCalib.at(i) - m_TOF.CalTOF();
				m_HitTimeCalibTOF.push_back(timeTOF);
			}
		}


	}


	if(trackrecevt){
		LogInfo << "Fill Muon Track" << std::endl;
		m_nTrack = trackrecevt->nTracks();
		const auto tracks = trackrecevt->tracks();
		
		LogInfo << "nTracks: " << m_nTrack << " vector size " << tracks.size() << std::endl;
		for(auto it = tracks.begin(); it != tracks.end(); it++){
			m_MuQuality.push_back((*it)->quality());
			m_MuPE.push_back((*it)->peSum());

			const CLHEP::HepLorentzVector entry = (*it)->start();
			// m_MuEntryX = entry.x(); m_MuEntryY = entry.y(); m_MuEntryZ = entry.z();
			m_MuEntryX.push_back(entry.x());
			m_MuEntryY.push_back(entry.y());
			m_MuEntryZ.push_back(entry.z());
			m_MuEntryTheta.push_back(entry.theta());
			m_MuEntryPhi.push_back(entry.phi());

			const CLHEP::HepLorentzVector exit = (*it)->end();
			m_MuExitX.push_back(exit.x());
			m_MuExitY.push_back(exit.y());
			m_MuExitZ.push_back(exit.z());
		}
	}

	// if(calibeventLPMT || calibeventSPMT || calibeventWP || oecevt) //If there is a trigger fill all trees

	if(m_Tag!="" || m_OECtag != ""){
		LogInfo << "Filling event tree" << std::endl;
		m_ntuple1->Fill();
	}

	 return true;

}

bool SelectionAlg::Book_tree()
{
	SniperPtr<RootWriter> svc(*getRoot(),"RootWriter");

	m_ntuple2 = svc->bookTree(*m_par,"Data/RunInfo", "EventInfo");
	// m_ntuple2->Branch("EntryNb", &m_iEvt, "EntryNb/I");
	m_ntuple2->Branch("RunNb", &m_iRun,"Run/I");
	m_ntuple2->Branch("RunTime", &runtime);
	m_ntuple2->Branch("EffectiveRunTime", &effruntime);
	m_ntuple2->Branch("nCDMuon", &nCDMuons);
	m_ntuple2->Branch("nCDWPMuon", &nCDWPMuons);
	m_ntuple2->Branch("nWPMuon", &nWPMuons);
	m_ntuple2->Branch("nNeutron", &nNeutrons);
	m_ntuple2->Branch("nIBD", &nIBD);
	
	
	m_ntuple1 = svc->bookTree(*m_par,"Data/event", "Event Level Tree");
	m_ntuple1->Branch("Filename", &m_fname);
	m_ntuple1->Branch("EntryNb", &m_iEvt, "EntryNb/I");
	m_ntuple1->Branch("Tag", &m_Tag);
	m_ntuple1->Branch("OecTag", &m_OECtag);
	m_ntuple1->Branch("RunNb", &m_iRun,"RunNb/I");
	m_ntuple1->Branch("TriggerName", &m_TriggerName);
	m_ntuple1->Branch("TimeStamp",&m_TimeStamp,"TimeStamp/l");
	m_ntuple1->Branch("ChargeTotLPMT",&m_ChargeTotLPMT,"ChargeTotLPMT/D");
	m_ntuple1->Branch("ChargeTotSPMT",&m_ChargeTotSPMT,"ChargeTotSPMT/D");
	m_ntuple1->Branch("ChargeTotWP",&m_ChargeTotWP,"ChargeTotWP/D");
	m_ntuple1->Branch("NbHitLPMTCalib", &m_NbHitLPMTCalib, "NbHitLPMTCalib/I");
	m_ntuple1->Branch("NbHitSPMTCalib", &m_NbHitSPMTCalib, "NbHitSPMTCalib/I");
	m_ntuple1->Branch("NbHitWPCalib", &m_NbHitWPCalib, "NbHitWPCalib/I");
	m_ntuple1->Branch("HitTimeRMS", &m_HitTime_std);
	m_ntuple1->Branch("HitTimeMean", &m_HitTime_mean);
	m_ntuple1->Branch("HitTimeCalibTOF", &m_HitTimeCalibTOF);
	if(saveHitInfo){
		m_ntuple1->Branch("PmtIDCalib", &m_PmtIdCalib);
		m_ntuple1->Branch("ChargeCalib", &m_ChargeCalib);
	}
	m_ntuple1->Branch("TotalPERec", &m_TotalPERec);
	m_ntuple1->Branch("RecEnergy", &m_RecE, "RecEnergy/F");
	m_ntuple1->Branch("Recx", &m_RecX, "Recx/F");
	m_ntuple1->Branch("Recy", &m_RecY, "Recy/F");
	m_ntuple1->Branch("Recz", &m_RecZ, "Recz/F");
	m_ntuple1->Branch("RecT0", &m_T0, "RecT0/F");
	m_ntuple1->Branch("PositionQuality", &m_PosQuality);
	m_ntuple1->Branch("EnergyQuality", &m_EnergyQuality);
	m_ntuple1->Branch("RecChi2", &m_chisq);
	// m_ntuple1->Branch("nTrack", &m_nTrack);
	// m_ntuple1->Branch("MuPE", &m_MuPE);
	// m_ntuple1->Branch("MuEntryX", &m_MuEntryX);
	// m_ntuple1->Branch("MuEntryY", &m_MuEntryY);
	// m_ntuple1->Branch("MuEntryZ", &m_MuEntryZ);
	// m_ntuple1->Branch("MuEntryTheta", &m_MuEntryTheta);
	// m_ntuple1->Branch("MuEntryPhi", &m_MuEntryPhi);
	// m_ntuple1->Branch("MuExitX", &m_MuExitX);
	// m_ntuple1->Branch("MuExitY", &m_MuExitY);
	// m_ntuple1->Branch("MuExitZ", &m_MuExitZ);
	// // m_ntuple1->Branch("MuDX", &m_MuDX);
	// // m_ntuple1->Branch("MuDY", &m_MuDY);
	// // m_ntuple1->Branch("MuDZ", &m_MuDZ);
	// m_ntuple1->Branch("MuQuality", &m_MuQuality);


	return true;
}
void SelectionAlg::clearAllTrees()
{

	m_Tag = "";

	//Calib
	m_ChargeTotLPMT=0;
	m_ChargeTotUpLPMT=0;
	m_ChargeTotSPMT=0;
	m_ChargeTotWP=0;
	m_TimeStampInNanoSec=0;
	m_NbHitLPMTCalib=0;
	m_NbHitSPMTCalib=0;
	m_NbHitWPCalib=0;
	m_PmtIdCalib.clear();
	m_HitTimeCalib.clear();
	m_HitTimeCalibTOF.clear();
	m_ChargeCalib.clear();

	//Oec
	m_myOECtag = 0;
	m_OECtag = "";
	m_TriggerName.clear();
	m_NFiredPMT=0;

	// Reco
	m_TotalPERec=0;
	m_RecE=0;
	m_RecX=0;
	m_RecY=0;
	m_RecZ=0;
	m_T0=0;
	m_PosQuality = 0.0;
	m_EnergyQuality = 0.0;
	m_chisq = 0.0;
	m_nTrack = 0;
	m_MuPE.clear();
	m_MuEntryX.clear(); 
	m_MuEntryY.clear(); 
	m_MuEntryZ.clear();
	m_MuEntryTheta.clear();
	m_MuEntryPhi.clear();
	m_MuExitX.clear(); 
	m_MuExitY.clear();
	m_MuExitZ.clear();
	// m_MuDX.clear(); 
	// m_MuDY.clear(); 
	// m_MuDZ.clear();
	m_MuQuality.clear();

	PMT_R=0;
	LS_R=0;

	RfrIndxLS=0;
	RfrIndxWR=0;
	c=0;
}
