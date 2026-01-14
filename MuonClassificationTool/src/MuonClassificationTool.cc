#include "MuonClassificationTool/MuonClassificationTool.h"
#include "SniperKernel/ToolFactory.h"
#include "SniperKernel/SniperDataPtr.h"
#include "SniperKernel/SniperPtr.h"

DECLARE_TOOL(MuonClassificationTool);

MuonClassificationTool::MuonClassificationTool(const std::string& name)
    : ToolBase(name), m_buf(NULL), m_calibevt(NULL), m_wpcalibevt(NULL)
{
    // can add properties for veto cut
}
MuonClassificationTool::~MuonClassificationTool(){
}

bool MuonClassificationTool::initialize(){
    LogInfo << "Initializing Muon classification tool..." << std::endl;

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

bool MuonClassificationTool::finalize(){
    return true;
}

bool MuonClassificationTool::isInDeadTime(JM::EvtNavigator* nav, JM::OecEvt* tOecEvt, const std::string& muTag){
    if(!nav || !tOecEvt){
        LogInfo << "EvtNavigator or OecEvt not found" << std::endl;
        return false;
    }

    //GetNavBuffer iterator
    JM::NavBuffer::Iterator navit = m_buf->find(nav);
    const TTimeStamp& ttime = tOecEvt->getTime();

    JM::NavBuffer::Iterator it = navit - 1;
    while(navit != m_buf->begin() && it != m_buf->begin()){
        --it;

        std::string thisTag = m_eventTagSvc->getTag(it->get());

        JM::OecHeader* bOecHdr = JM::getHeaderObject<JM::OecHeader>(it->get());
        if(!bOecHdr) continue;
        JM::OecEvt* bOecEvt = dynamic_cast<JM::OecEvt*>(bOecHdr->event("JM::OecEvt"));
        if(!bOecEvt){
            LogInfo << "Could not load OecEvt" << std::endl;
            continue;
        }
        const TTimeStamp& btime = bOecEvt->getTime();

        double dtime = ((ttime.GetSec() -  btime.GetSec())*1000000000ULL + (ttime.GetNanoSec() - btime.GetNanoSec()));

        if((dtime * 1e-3 > 4 && muTag == "WP") || (dtime * 1e-3 > 50 && muTag == "CD")) { // 4us dead time for WP and 50us dead time for CD
            LogInfo << "Muon candidate not in dead time" << std::endl;
            return false;
        }
        else if( ((thisTag == "WPMuon" || thisTag == "CDWPMuon") && muTag == "WP") || ((thisTag == "CDMuon" || thisTag == "CDWPMuon") && muTag == "CD") ) {
            LogInfo << "Muon candidate in dead time !" << std::endl;
            return true;
        }
    }

    return false;

}

bool MuonClassificationTool::isVetoed(JM::EvtNavigator* nav, JM::OecEvt* tOecEvt) {
    if(!nav || !tOecEvt){
        LogInfo << "EvtNavigator or OecEvt not found" << std::endl;
        return false;
    }

    //Get NavBuffer Iterator
    JM::NavBuffer::Iterator navit = m_buf->find(nav);
    const TTimeStamp& ttime = tOecEvt->getTime();

    JM::NavBuffer::Iterator it = navit;

    while (it != m_buf->begin()) {
        --it;

        std::string thisTag = m_eventTagSvc->getTag(it->get());

        JM::OecHeader* bOecHdr = JM::getHeaderObject<JM::OecHeader>(it->get());
        if(!bOecHdr) continue;
        JM::OecEvt* bOecEvt = dynamic_cast<JM::OecEvt*>(bOecHdr->event("JM::OecEvt"));
        if(!bOecEvt){
            LogInfo << "Could not load OecEvt" << std::endl;
            continue;
        }
        const TTimeStamp& btime = bOecEvt->getTime();

        double dtime = ((ttime.GetSec() - btime.GetSec())*1000000000ULL + (ttime.GetNanoSec() - btime.GetNanoSec()));

        if(dtime * 1e-6 > 2) { // 2ms after last muon
            return false;
        }
        else if (thisTag == "WPMuon" || thisTag == "CDMuon" || thisTag == "CDWPMuon"){
            return true;
        }
    }

    return false;
}

// MyClassificationType MuonClassificationTool::classify(JM::EvtNavigator* nav) {
bool MuonClassificationTool::isMuon(JM::EvtNavigator* nav) {
    LogInfo << "Muon Search" << std::endl;

    if(!nav){
        LogInfo << "EvtNavigator not found" << std::endl;
        return false;
    }

    //Get NavBuffer Iterator
    JM::NavBuffer::Iterator navit = m_buf->find(nav);

    JM::OecHeader* oechdr = JM::getHeaderObject<JM::OecHeader>(nav);
    if(!oechdr) return false;
    LogDebug << "Header Set" << std::endl;
    JM::OecEvt* oecevt = dynamic_cast<JM::OecEvt*>(oechdr->event("JM::OecEvt"));
    if(!oecevt) return false;
    LogDebug << "Oec Evt Set at " << oecevt << std::endl;

    const TTimeStamp& ttime = oecevt->getTime();
    LogDebug << "TimeStamp set" << std::endl;

    float charge = oecevt->getTotalCharge();

    
    JM::WpCalibHeader* wpcalibhdr = JM::getHeaderObject<JM::WpCalibHeader>(nav);
    JM::CdLpmtCalibHeader* cdcalibhdr = JM::getHeaderObject<JM::CdLpmtCalibHeader>(nav);

    LogInfo << "Muon candidate: " << charge << " PE" << std::endl;
    LogInfo << "Wp event " << wpcalibhdr << std::endl;
    LogInfo << "Cd event " << cdcalibhdr << std::endl;

    if(wpcalibhdr && charge > WPMuonChargeCut && !isInDeadTime(nav, oecevt, "WP")){
        LogInfo << "Wp Muon with charge = " << charge << " PE" << std::endl;
        m_eventTagSvc->addTag(nav, "WPMuon");
        return true;
    }

    if(cdcalibhdr && charge > CDMuonChargeCut && !isInDeadTime(nav, oecevt, "CD")){
        //look for WP Muon before (500 ns)
        LogInfo << "Searching WP correlation 500 ns before" << std::endl;
        JM::NavBuffer::Iterator it = navit - 1;
        while(navit != m_buf->begin() && it != m_buf->begin()){
            
            JM::OecHeader* bOecHdr = JM::getHeaderObject<JM::OecHeader>(it->get());
            if(!bOecHdr) continue;
            JM::OecEvt* bOecEvt = dynamic_cast<JM::OecEvt*>(bOecHdr->event("JM::OecEvt"));
            if(!bOecEvt) continue;
            
            auto bWpHdr = JM::getHeaderObject<JM::WpCalibHeader>(it->get()); 
            
            float wpcharge = bOecEvt->getTotalCharge();
            const TTimeStamp& btime = bOecEvt->getTime();
            double dtime = ((ttime.GetSec() -  btime.GetSec())*1000000000ULL + (ttime.GetNanoSec() - btime.GetNanoSec()));
            
            if(dtime > 500) break;
            else if (bWpHdr && wpcharge > 400 ){
                LogInfo << "CdWp Muon with charge = " << charge << " PE" << std::endl;
                
                m_eventTagSvc->addTag(nav, "CDWPMuon");
                return true;
            }
            
            --it;
        }

        //look for WP Muon after (500 ns)
        LogInfo << "Searching WP correlation 500 ns after" << std::endl;
        it = navit + 1;
        while( navit != m_buf->end() && it != m_buf->end()){
            
            auto aOecHdr = JM::getHeaderObject<JM::OecHeader>(it->get());
            if(!aOecHdr) continue;
            JM::OecEvt* aOecEvt = dynamic_cast<JM::OecEvt*>(aOecHdr->event("JM::OecEvt"));
            
            auto aWpHdr = JM::getHeaderObject<JM::WpCalibHeader>(it->get()); 
            
            float wpcharge = aOecEvt->getTotalCharge();
            const TTimeStamp& atime = aOecEvt->getTime();
            double dtime = ((atime.GetSec() -  ttime.GetSec())*1000000000ULL + (atime.GetNanoSec() - ttime.GetNanoSec()));
            
            
            if(dtime > 500) break;
            else if (aWpHdr && wpcharge > 400) {
                LogInfo << "CdWp Muon with charge = " << charge << " PE" << std::endl;
                m_eventTagSvc->addTag(nav, "CDWPMuon");
                return true;
            }
            ++it;
        }
        

        if(!isVetoed(nav, oecevt)){
            LogInfo << "Cd Muon with charge = " << charge << " PE" << std::endl;
            m_eventTagSvc->addTag(nav, "CDMuon");
            return true;
        }
    }

    return false;
}