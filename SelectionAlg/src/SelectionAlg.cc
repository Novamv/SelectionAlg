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

	PMT_R = 19.434; //m
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
	i_pBiPo214pair= m_tagsvc->getpTag("BiPo214Pair");
	i_dBiPo214pair= m_tagsvc->getdTag("BiPo214Pair");


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
	++m_iEvt;
	LogInfo << "executing: " << m_iEvt << std::endl;
	
	TString fullpath = gDirectory->GetFile()->GetName();
	prevname = gSystem->BaseName(fullpath);
	if (prevname != m_fname){
		// gDirectory->pwd();
		m_fname = prevname;
		std::cout << "Current File: " << m_fname << std::endl;
 	}

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
	
	LogDebug << "Detector type is  " <<nav->getDetectorType()<<std::endl;
	LogDebug << "Start to Explore SmartRef: " << std::endl;
	LogDebug << "Size of paths: " << paths.size() << std::endl;
	LogDebug << "Size of refs: " << refs.size() << std::endl;
	
	for (size_t i = 0; i < paths.size(); ++i) {
		LogDebug << refs[i]<<" -> ref: " << std::endl;
		const std::string& path = paths[i];
		JM::SmartRef* ref = refs[i];
		JM::EventObject* evtobj = ref->GetObject();
		
		LogDebug << " path: " << path
		<< " ref->entry(): " << ref->entry()
		<< " evtobj: " << evtobj;
		
		if (path=="/Event/Sim") {
			auto hdr = dynamic_cast<JM::SimHeader*>(evtobj);
			LogDebug <<i<<" SimHeader: " << hdr;
		}
		LogDebug << std::endl;
	}

	clearAllTrees(); //Clearing all trees variables
	
// ================================================
// =============  Inititalization  =============
// ================================================

	theTime = nav->TimeStamp();

	if(m_iEvt == 0) {
		PreviousTime = theTime;
		skipReason = SkipReason::StartOfFile; // Will skip the first 1.2 sec of the file
		skipStartTime = theTime;
		FirstTime = theTime;
		dtCD = 0;
		dtWP = 0;
		
	}

	// =============  Fast skip  =============
	if(skipReason != SkipReason::None){
		
		dt_skip = deltaT_ns(theTime, skipStartTime);

		bool skip1 = (skipReason == SkipReason::StartOfFile || skipReason == SkipReason::BigGap) && dt_skip < 1'200'000'000LL; //1.2 s skip if Start of file of big gap
		bool skip2 = skipReason == SkipReason::MissingHeader && dt_skip < 5'000'000LL; // 5 ms skip if no headers

		if(skip1 || skip2){
			PreviousTime = theTime;
			return true;
		}

		prevCDTime = theTime;
		prevWPTime = theTime;
		dtCD = 0;
		dtWP = 0;
		skipReason = SkipReason::None;
		
	}

	// =============  Compute dt  =============

	int64_t globalTime = deltaT_ns(theTime, FirstTime);
	LogInfo << "Global time: " << globalTime << std::endl;

	int64_t dt 		   = deltaT_ns(theTime, PreviousTime);
	int64_t dtLastMuon = deltaT_ns(theTime, tLastMuon);

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
		if(dtCD == 0) dtCD = deltaT_ns(theTime, PreviousTime);
		else 		  dtCD = deltaT_ns(theTime, prevCDTime);
		
		LogDebug << "dtCD: " << dtCD << std::endl;
		prevCDTime = theTime;

		if(dtCD * 1e-6 > 50.0) {
			skipReason = SkipReason::BigGap;
			skipStartTime = theTime;
			PreviousTime = theTime;
			LogInfo << "New skip reason: CD gap > 50 ms" << std::endl;
			return true;
		}
	}
	if(calibheaderWP){
		if(dtWP == 0) dtWP = deltaT_ns(theTime, PreviousTime);
		else 		  dtWP = deltaT_ns(theTime, prevWPTime);
		
		LogDebug << "dtWP: " << dtWP << std::endl;
		prevWPTime = theTime;
		
		if(dtWP * 1e-6 > 70.0) {
			skipReason = SkipReason::BigGap;
			skipStartTime = theTime;
			PreviousTime = theTime;
			LogInfo << "New skip reason: WP gap > 70 ms" << std::endl;
			return true;
		}
	}
	
	// =============  Runtime counting  =============
	
	runtime += dt;
	if(dtLastMuon * 1e-6 > 5.0){ // 5 ms muon veto 
		effruntime += dt;
	}
	
	PreviousTime = theTime;
	
	LogDebug << "Run Time: " << runtime << std::endl;
	LogDebug << "Effective Run Time: " << effruntime << std::endl;

// ================================================
// =============  Run Classification  =============
// ================================================

	// Quick check for Delay event 
	for (auto& cand : m_pendingIBD){
		int64_t dt = deltaT_ns(theTime, cand.pTimeStamp);
		if(dt > 1'000'000LL) continue;  // will be expired later
		if(m_iEvt == cand.promptEntry + cand.delayEntryOffset) {
			m_Tag = "Delay";
			break;
    	}
	}

	LogInfo << "dtLastMuon: " << dtLastMuon << std::endl;

	if(dtLastMuon * 1e-3 > 50.0 && m_MuClassifier->isMuon(nav)){ //50 us gap time
		tLastMuon = theTime;
		m_Tag = m_eventTagsvc->getTag(nav);
		if(m_Tag == "CDMuon") nCDMuons++;
		if(m_Tag == "CDWPMuon") nCDWPMuons++;
		if(m_Tag == "WPMuon") nWPMuons++;
	}
	else if((dtLastMuon * 1e-3 > 20.0 && dtLastMuon * 1e-6 < 2.0) && m_NeutronClassifier->isSpalNeutron(nav)){
		m_Tag = m_eventTagsvc->getTag(nav);
		nNeutrons++;
	}
	else if(dtLastMuon * 1e-6 > 5.0 && m_IBDClassifier->isPrompt(nav)) { // If IBD out of Muon veto fill prompt information
		m_Tag = m_eventTagsvc->getTag(nav);
		m_DelayEvt = m_iEvt + m_IBDClassifier->getDelayOffset();
		nIBD++;
	}

	LogInfo << "Current Tag: " << m_Tag << std::endl;
	if(m_Tag == "") return true;


// ================================================
// ========  Fill Trees for tagged events  ========
// ================================================

	// ---------------------- Load EDM classes ----------------------
	

	auto simheader = JM::getHeaderObject<JM::SimHeader>(nav);
	if(simheader){
		simevent = (JM::SimEvt*)simheader->event();
		LogDebug << "SimEvent Read in: " << simevent << std::endl;
		LogDebug << "SimEvent Track: " << simevent->getTracksVec().size() << std::endl;
		LogDebug << "SimEvent Hits: " << simevent->getCDHitsVec().size() << std::endl;
	}

	// auto recheader = JM::getHeaderObject<JM::CdVertexRecHeader>(nav, recEDMPath);
	if (recheader) {
	  recevent = recheader->event();
	  LogDebug << "RecEvent Read in: " << recevent << std::endl;
	}
	auto trackheader = JM::getHeaderObject<JM::CdTrackRecHeader>(nav, "/Event/CdTrackRecClassify");
	if(trackheader){
		trackrecevt = trackheader->event();
		LogDebug << "TrackRecEvt Read in: " << trackrecevt << std::endl;
	}
	// auto calibheaderLPMT = JM::getHeaderObject<JM::CdLpmtCalibHeader>(nav);
	if (calibheaderLPMT) {
	  calibeventLPMT = calibheaderLPMT->event();
	  LogDebug << "CalibEventLPMT Read in: " << calibeventLPMT << std::endl;
	}
	auto calibheaderSPMT = JM::getHeaderObject<JM::CdSpmtCalibHeader>(nav);
	if (calibheaderSPMT) {
	  calibeventSPMT = calibheaderSPMT->event();
	  LogDebug << "CalibEventSPMT Read in: " << calibeventSPMT << std::endl;
	}
	// auto calibheaderWP = JM::getHeaderObject<JM::WpCalibHeader>(nav);
	if (calibheaderWP) {
	  calibeventWP = calibheaderWP->event();
	  LogDebug << "CalibEventWP Read in: " << calibeventWP << std::endl;
	}

	auto triggerheader = JM::getHeaderObject<JM::CdTriggerHeader>(nav);
	if(triggerheader){triggerevent = triggerheader->event();
	  LogDebug <<"CD TriggerEvent Read in: " << triggerevent <<std::endl;
	}

	auto OEChdr = JM::getHeaderObject<JM::OecHeader>(nav);
	if(OEChdr){
		oecevt = dynamic_cast<JM::OecEvt*>(OEChdr->event("JM::OecEvt"));
		LogDebug <<"OEC Event Read in: " << oecevt <<std::endl;
	}

	// ---------------------- Set all trees variables ----------------------

	if(oecevt) {
		m_myOECtag = oecevt->getTag();
		if((m_myOECtag & i_pIBD)== i_pIBD) m_OECtag = "PromptIBD";
        if((m_myOECtag & i_dIBD)== i_dIBD) m_OECtag = "DelayIBD";
		if((m_myOECtag & i_pBiPo214pair)== i_pBiPo214pair) m_OECtag = "PromptBiPo";
		if((m_myOECtag & i_dBiPo214pair)== i_dBiPo214pair) m_OECtag = "DelayBiPo";
	}

	bool myIBD = (m_Tag == "Prompt" || m_Tag == "Delay");
	bool OecIBD = (m_OECtag == "PromptIBD" || m_OECtag == "DelayIBD");
	bool OecBiPo = (m_OECtag == "PromptBiPo" || m_OECtag == "DelayBiPo");

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
	m_iRun = RunNumber;

	m_TimeStamp = toKey(timestamp);

	double ChargeTot=0;
	std::vector<int> tempPmtIds;
	std::vector<double> tempHitTimes;
	std::vector<double> tempCharges;

	if(calibeventLPMT){

		// TH1F* hTime = new TH1F("hTime", "hTime", 100, 0, 1000);
		const auto& chhlistLPMT = calibheaderLPMT->event()->calibPMTCol();
		float sum = 0, sum2 = 0;
		float nhit_sum = 0, nhit_sum2 = 0;
		int n = 0, nhit = 0;
		for (auto chit = chhlistLPMT.begin(); chit!=chhlistLPMT.end(); ++chit){
			auto calib = *chit;
			unsigned int pmtId = calib->pmtId();
			Identifier id = Identifier(pmtId);
			int TruePM=CdID::module(id);
			m_NbHitLPMTCalib+=calib->size();

			//for flasher cut
			nhit++;
			nhit_sum2 += calib->size() * calib->size();
			
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

		float nhit_std = std::sqrt((nhit_sum2 - m_NbHitLPMTCalib*m_NbHitLPMTCalib/nhit) / (nhit - 1));

		// If the event is a flasher, skip the event
		bool flasher = std::pow((nhit_std - 0.55)/0.45, 2) + std::pow((m_HitTime_std - 170)/80, 2) >= 1;

		LogInfo << "Is Flasher: " << flasher << std::endl;

		if((m_Tag == "Prompt" || m_Tag == "Delay") && flasher) return true; //only skip flashers if it's an IBD tag

		if(m_Tag == "Neutron" && m_HitTime_std < 275){
			m_eventTagsvc->addTag(nav, "SpalNeutron");
			m_Tag = "SpalNeutron";
		}


		// Fill Hit level information only for IBD tagged events
		if(myIBD || m_OECtag != ""){
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
				double time = calibSPMT->time(j);
				double CalibSPMTBackADC = 124.*0.48*calibSPMT->charge(j)+76.;				
				ChargeTot+=calibSPMT->charge(j);
				
				if(myIBD){
					m_PmtIdCalib.push_back(TruePM);
					m_HitTimeCalib.push_back(time);
					m_ChargeCalib.push_back(calibSPMT->charge(j));
				}
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


		// Correct HitTime with TOF
		TVector3 vertex(m_RecX, m_RecY, m_RecZ);
		for (size_t i = 0; i < m_PmtIdCalib.size(); i++)
		{
			int pmtid = m_PmtIdCalib.at(i);
			m_TOF.setVertex(vertex);
			m_TOF.setInterface(interface);

			if(pmtid < 17612) m_TOF.setPMTVertex(ALL_LPMT_pos.at(pmtid));
			else if (pmtid >= 20000) m_TOF.setPMTVertex(ALL_SPMT_pos.at(pmtid - 20000));

			double timeTOF = m_HitTimeCalib.at(i) - m_TOF.CalTOF();
			m_HitTimeCalibTOF.push_back(timeTOF);
		}

		if(m_Tag == "SpalNeutron"){
			NeutronVertex = vertex;
			NeutronTime = nav->TimeStamp();
		}
		else if(m_Tag == "Prompt"){ // To be removed
			float dR = (NeutronVertex - vertex).Mag();
			const TTimeStamp& pTime = nav->TimeStamp();
			double dt_Neutron = (pTime.GetSec() - NeutronTime.GetSec())*1000000000ULL + (pTime.GetNanoSec() - NeutronTime.GetNanoSec());
			if(dt_Neutron*1e-9 < 1.2 && dR < 4000){
				m_NeutronVeto = 1;
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

// =============  Fill Trees and check pending Prompts  =============

	for(auto it = m_pendingIBD.begin(); it != m_pendingIBD.end(); ) {
		int64_t dt_since_prompt = deltaT_ns(theTime, it->pTimeStamp);
		if(dt_since_prompt > 1'000'000LL) { // 1.0 ms
			LogInfo << "Expired IBD Prompt (entry " << it->promptEntry << ")" << std::endl;
			it = m_pendingIBD.erase(it);
			continue;
		}
		++it;
	}

	if(m_Tag == "Prompt"){
		PendingIBD cand;
		cand.promptEntry 		= m_iEvt;
		cand.delayEntryOffset 	= m_IBDClassifier->getDelayOffset();
		cand.pTimeStamp 		= theTime;
		cand.pHitTimeTOF 		= m_HitTimeCalibTOF;
		cand.pEnergy 			= m_RecE;
		cand.pX 				= m_RecX; 
		cand.pY 				= m_RecY;
		cand.pZ 				= m_RecZ;
		cand.pNPE 				= m_ChargeTotLPMT;
		cand.NeutronVeto		= (int)m_IBDClassifier->isNeutronVetoed(nav);

		m_pendingIBD.push_back(cand);
	}
	else{
		for(auto it = m_pendingIBD.begin(); it != m_pendingIBD.end();){

			int64_t dt_since_prompt = deltaT_ns(theTime, it->pTimeStamp);
			if(m_iEvt == it->promptEntry + it->delayEntryOffset){
				LogInfo << "IBD delay match prompt entry" << std::endl;
				m_Tag = "Delay";

				TVector3 pVtx(it->pX, it->pY, it->pZ);
                TVector3 dVtx(m_RecX, m_RecY, m_RecZ);

				m_pair.dt_ns = dt_since_prompt;
				m_pair.dR_mm = (float)(pVtx - dVtx).Mag();
				
				m_pair.promptEntry = it->promptEntry;
				m_pair.pEnergy     = it->pEnergy;
                m_pair.pX          = it->pX;
                m_pair.pY          = it->pY;
                m_pair.pZ          = it->pZ;
                m_pair.pNPE        = it->pNPE;
                m_pair.pTime	   = toKey(it->pTimeStamp);
				m_pair.pHitTimeTOF = it->pHitTimeTOF;

				m_pair.delayEntry  = m_iEvt;
				m_pair.dEnergy 	   = m_RecE;
				m_pair.dX          = m_RecX;
                m_pair.dY          = m_RecY;
                m_pair.dZ          = m_RecZ;
                m_pair.dNPE        = m_ChargeTotLPMT;
                m_pair.dTime	   = toKey(theTime);
				m_pair.dHitTimeTOF = m_HitTimeCalibTOF;
				
				m_pair.neutronVeto = it->NeutronVeto;

				m_ibdtree->Fill(); 
                LogInfo << "IBD pair filled: dt="
                        << m_pair.dt_ns * 1e-3 << " µs  dR="
                        << m_pair.dR_mm        << " mm" << std::endl;

				it = m_pendingIBD.erase(it);
				continue;
			}

			++it;
		}
	}


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
	m_ntuple2->Branch("RunNb", 			  &m_iRun,"Run/I");
	m_ntuple2->Branch("RunTime", 		  &runtime);
	m_ntuple2->Branch("EffectiveRunTime", &effruntime);
	m_ntuple2->Branch("nCDMuon", 		  &nCDMuons);
	m_ntuple2->Branch("nCDWPMuon", 		  &nCDWPMuons);
	m_ntuple2->Branch("nWPMuon", 		  &nWPMuons);
	m_ntuple2->Branch("nNeutron", 		  &nNeutrons);
	m_ntuple2->Branch("nIBD", 			  &nIBD);
	
	
	m_ntuple1 = svc->bookTree(*m_par,"Data/event", "Event Level Tree");
	m_ntuple1->Branch("Filename", 		 &m_fname);
	m_ntuple1->Branch("EntryNb", 		 &m_iEvt, "EntryNb/I");
	m_ntuple1->Branch("Tag", 			 &m_Tag);
	m_ntuple1->Branch("OecTag", 		 &m_OECtag);
	m_ntuple1->Branch("RunNb", 		 	 &m_iRun,"RunNb/I");
	m_ntuple1->Branch("TriggerName", 	 &m_TriggerName);
	m_ntuple1->Branch("TimeStamp",		 &m_TimeStamp,"TimeStamp/l");
	m_ntuple1->Branch("ChargeTotLPMT",	 &m_ChargeTotLPMT,"ChargeTotLPMT/D");
	m_ntuple1->Branch("ChargeTotSPMT",	 &m_ChargeTotSPMT,"ChargeTotSPMT/D");
	m_ntuple1->Branch("ChargeTotWP",	 &m_ChargeTotWP,"ChargeTotWP/D");
	m_ntuple1->Branch("NbHitLPMTCalib",  &m_NbHitLPMTCalib, "NbHitLPMTCalib/I");
	m_ntuple1->Branch("NbHitSPMTCalib",  &m_NbHitSPMTCalib, "NbHitSPMTCalib/I");
	m_ntuple1->Branch("NbHitWPCalib", 	 &m_NbHitWPCalib, "NbHitWPCalib/I");
	m_ntuple1->Branch("NeutronVeto", 	 &m_NeutronVeto, "NeutronVeto/I");
	m_ntuple1->Branch("HitTimeRMS", 	 &m_HitTime_std);
	m_ntuple1->Branch("HitTimeMean", 	 &m_HitTime_mean);
	m_ntuple1->Branch("HitTimeCalibTOF", &m_HitTimeCalibTOF);
	if(saveHitInfo){
		m_ntuple1->Branch("PmtIDCalib",  &m_PmtIdCalib);
		m_ntuple1->Branch("ChargeCalib", &m_ChargeCalib);
	}
	m_ntuple1->Branch("TotalPERec", 	 &m_TotalPERec);
	m_ntuple1->Branch("RecEnergy", 		 &m_RecE, "RecEnergy/F");
	m_ntuple1->Branch("Recx", 			 &m_RecX, "Recx/F");
	m_ntuple1->Branch("Recy", 			 &m_RecY, "Recy/F");
	m_ntuple1->Branch("Recz",			 &m_RecZ, "Recz/F");
	m_ntuple1->Branch("RecT0", 			 &m_T0, "RecT0/F");
	m_ntuple1->Branch("PositionQuality", &m_PosQuality);
	m_ntuple1->Branch("EnergyQuality", 	 &m_EnergyQuality);
	m_ntuple1->Branch("RecChi2", 	  	 &m_chisq);

	// m_ntuple3 = svc->bookTree(*m_par, "Data/muons", "Muon Reco Tree");
	// m_ntuple3->Branch("EntryNb", &m_iEvt);
	// m_ntuple3->Branch("MuNTrack", &m_nTrack);
	// m_ntuple3->Branch("MuPE", &m_MuPE);
	// m_ntuple3->Branch("MuEntryX", &m_MuEntryX);
	// m_ntuple3->Branch("MuEntryY", &m_MuEntryY);
	// m_ntuple3->Branch("MuEntryZ", &m_MuEntryZ);
	// m_ntuple3->Branch("MuEntryTheta", &m_MuEntryTheta);
	// m_ntuple3->Branch("MuEntryPhi", &m_MuEntryPhi);
	// m_ntuple3->Branch("MuExitX", &m_MuExitX);
	// m_ntuple3->Branch("MuExitY", &m_MuExitY);
	// m_ntuple3->Branch("MuExitZ", &m_MuExitZ);
	// // m_ntuple1->Branch("MuDX", &m_MuDX);
	// // m_ntuple1->Branch("MuDY", &m_MuDY);
	// // m_ntuple1->Branch("MuDZ", &m_MuDZ);
	// m_ntuple3->Branch("MuQuality", &m_MuQuality);

	// ============= IBD Pair Tree =============
	m_ibdtree = svc->bookTree(*m_par, "Data/ibd", "IBD pair events");
	m_ibdtree->Branch("Filename",    &m_fname);
	m_ibdtree->Branch("RunNb", 		 &m_pair.runNumber);
	m_ibdtree->Branch("dt", 		 &m_pair.dt_ns);
	m_ibdtree->Branch("dR", 		 &m_pair.dR_mm);
	//Prompt
	m_ibdtree->Branch("pEntry", 	 &m_pair.promptEntry);
	m_ibdtree->Branch("pEnergy", 	 &m_pair.pEnergy);
	m_ibdtree->Branch("pX", 		 &m_pair.pX);
	m_ibdtree->Branch("pY", 		 &m_pair.pY);
	m_ibdtree->Branch("pZ", 		 &m_pair.pZ);
	m_ibdtree->Branch("pNPE", 		 &m_pair.pNPE);
	m_ibdtree->Branch("pTimeStamp",  &m_pair.pTime, "pTimeStamp/l");
	m_ibdtree->Branch("pHitTimeTOF", &m_pair.pHitTimeTOF);
	//Delay
	m_ibdtree->Branch("dEntry", 	 &m_pair.delayEntry);
	m_ibdtree->Branch("dEnergy",  	 &m_pair.dEnergy);
	m_ibdtree->Branch("dX", 	  	 &m_pair.dX);
	m_ibdtree->Branch("dY", 	  	 &m_pair.dY);
	m_ibdtree->Branch("dZ", 	  	 &m_pair.dZ);
	m_ibdtree->Branch("dNPE",	  	 &m_pair.dNPE);
	m_ibdtree->Branch("dTimeStamp",  &m_pair.dTime, "dTimeStamp/l");
	m_ibdtree->Branch("dHitTimeTOF", &m_pair.dHitTimeTOF);

	m_ibdtree->Branch("NeutronVeto", &m_pair.neutronVeto);

	return true;
}
void SelectionAlg::clearAllTrees()
{

	m_Tag = "";

	m_NeutronVeto = 0;

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
}
