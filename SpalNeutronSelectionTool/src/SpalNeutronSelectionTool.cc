#include "SpalNeutronSelectionTool/SpalNeutronSelectionTool.h"

#include "SniperKernel/ToolFactory.h"
#include "SniperKernel/SniperDataPtr.h"
#include "SniperKernel/SniperPtr.h"


DECLARE_TOOL(SpalNeutronSelectionTool);

SpalNeutronSelectionTool::SpalNeutronSelectionTool(const std::string& name)
    : ToolBase(name), m_buf(NULL), m_calibevt(NULL), m_recevt(NULL)
{
    declProp("recEDMPath", recEDMPath="/Event/CdVertexRec");
    declProp("FiducialVolume", FV_cut=17200); // Fiducial Volume cut (17.2 m default)

    NeutronEnergyCut[0] = 1.5; NeutronEnergyCut[1] = 20; // MeV // Could be a property
}

SpalNeutronSelectionTool::~SpalNeutronSelectionTool(){}

bool SpalNeutronSelectionTool::initialize(){
    LogInfo << "Initializing SpalNeutron classification tool..." << std::endl;

    SniperDataPtr<JM::NavBuffer> navbuf(getParent(), "/Event");
    if(navbuf.invalid()){
        LogError << "Failed to get nav buffer!" << std::endl;
        return false;
    }
    m_buf = navbuf.data();

    // get EventTagSvc
    SniperPtr<EventTagSvc> eventTagSvc(getParent(), "EventTagSvc");
    if (!eventTagSvc.valid()) {
        LogError << "Failed to get EventTagSvc!" << std::endl;
        return false;
    }
    m_eventTagSvc = eventTagSvc.data();

    return true;
}


bool SpalNeutronSelectionTool::finalize(){
    return true;
}



bool SpalNeutronSelectionTool::isSpalNeutron(JM::EvtNavigator* nav){
    
    LogInfo << "Neutron Search" << std::endl;

    if(!nav){
        LogInfo << "EvtNavigator not found" << std::endl;
        return false;
    }

    // OEC
    auto oechdr = JM::getHeaderObject<JM::OecHeader>(nav);
    if(!oechdr) return false;
    m_oecevt = dynamic_cast<JM::OecEvt*>(oechdr->event("JM::OecEvt"));
    if(!m_oecevt) return false;
    const TTimeStamp& ttime = m_oecevt->getTime();

    
    std::string lastMuTag = m_eventTagSvc->getLastMuTag();
    if(lastMuTag != "CDMuon" && lastMuTag != "CDWPMuon") return false;
    JM::EvtNavigator* lastMuNav = m_eventTagSvc->getLastMuNav();
    if(m_buf->find(lastMuNav) == m_buf->end()){
        LogInfo << "Last Muon out of NavBuffer memory" <<std::endl;
        return false;
    }
    auto mu_oechdr = JM::getHeaderObject<JM::OecHeader>(lastMuNav);
    if(!mu_oechdr){
        LogInfo << "Info: No Oec Header" << std::endl;
        return false;

    }
    JM::OecEvt* mu_oecevt = dynamic_cast<JM::OecEvt*>(mu_oechdr->event("JM::OecEvt"));
    const TTimeStamp& btime = mu_oecevt->getTime();

    // ------- Check if within dt window -------

    double dtime = (ttime.GetSec() - btime.GetSec())*1000000000ULL + (ttime.GetNanoSec() - btime.GetNanoSec()); //ns
    if(!(dtime > 20e3 && dtime < 2e6)){ // 20us < dt < 2ms
        return false;
    }

    
    // ------- Load EDM -------
    // Calib
    auto calibhdr = JM::getHeaderObject<JM::CdLpmtCalibHeader>(nav);
    if(!calibhdr) return false; //if no calib header skip
    m_calibevt = calibhdr->event();

    // could add neutron time distribution
    
    
    // Reco
    auto rechdr = JM::getHeaderObject<JM::CdVertexRecHeader>(nav, recEDMPath);
    if(!rechdr){
        LogInfo << "No Rec Header" << std::endl;
        return false;
    }
    m_recevt = rechdr->event();
    
    float recenergy = 0.0;
    const auto& vertices = m_recevt->vertices();

    for(auto vertex: vertices) {
        recenergy = vertex->energy();
        vecNeutron.SetXYZ(vertex->x(), vertex->y(), vertex->z());
    }

    LogInfo << "Neutron Candidate" << std::endl;
    LogInfo << " - energy: " << recenergy << std::endl;
    LogInfo << " - r: " << vecNeutron.Mag() << std::endl;

    // ------- Neutron selection -------
    bool position = vecNeutron.Mag() < FV_cut;
    bool energy = recenergy >= NeutronEnergyCut[0] && recenergy <= NeutronEnergyCut[1];
    if(position && energy){
        m_eventTagSvc->addTag(nav, "Neutron");
        return true;
    }
    
    return false;

}