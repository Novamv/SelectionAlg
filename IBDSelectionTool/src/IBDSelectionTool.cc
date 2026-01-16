#include "IBDSelectionTool/IBDSelectionTool.h"

#include "SniperKernel/ToolFactory.h"
#include "SniperKernel/SniperDataPtr.h"
#include "SniperKernel/SniperPtr.h"


DECLARE_TOOL(IBDSelectionTool);

IBDSelectionTool::IBDSelectionTool(const std::string& name)
    :ToolBase(name), m_oecevt(NULL), m_recevt(NULL)
{
    declProp("recEDMPath", recEDMPath="/Event/CdVertexRec");
    declProp("FiducialVolume", FV_cut); // Fiducial volume cut: 17.2 m default
    
    PromptEnergyCut[0] = 0.7; PromptEnergyCut[1] = 12.0;
    PromptChargeCut[0] = 1500; PromptChargeCut[1] = 21000;
    DelayEnergyCut[0] = 2.0; DelayEnergyCut[1] = 2.5;
    DelayChargeCut[0] = 4000; DelayChargeCut[1] = 6000;

}

IBDSelectionTool::~IBDSelectionTool(){}


bool IBDSelectionTool::initialize() {
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

bool IBDSelectionTool::finalize(){
    return true;
}


bool IBDSelectionTool::isVetoed(JM::OecEvt* pOecEvt, JM::EvtNavigator* munav) {
    const TTimeStamp& ptime = pOecEvt->getTime();

    LogInfo << "Checking Muon veto" << std::endl;

    if(m_buf->find(munav) == m_buf->end()){
        LogInfo << "Last Muon out of NavBuffer memory" <<std::endl;
        return false;
    }
    auto muOechdr = JM::getHeaderObject<JM::OecHeader>(munav);
    if(!muOechdr){
        LogInfo << "No Oec Header" << std::endl;
        return false;
    }
    JM::OecEvt* muOecEvt = dynamic_cast<JM::OecEvt*>(muOechdr->event("JM::OecEvt"));
    const TTimeStamp& muTime = muOecEvt->getTime();
    double dtime = (ptime.GetSec() - muTime.GetSec())*1000000000ULL + (ptime.GetNanoSec() - muTime.GetNanoSec());

    if(dtime*1e-6 > 5){ // 5ms muon veto
        LogInfo << "Out of veto window" << std::endl;
        return false;
    }

    return true;
}

bool IBDSelectionTool::isIsolated(JM::EvtNavigator* pnav, JM::EvtNavigator* dnav, JM::OecEvt* pOecEvt, JM::OecEvt* dOecEvt){

    JM::NavBuffer::Iterator pit = m_buf->find(pnav);
    JM::NavBuffer::Iterator dit = m_buf->find(dnav);

    const TTimeStamp& ptime = pOecEvt->getTime();
    const TTimeStamp& dtime = dOecEvt->getTime();

    // Check before Prompt
    LogInfo << "Info: Checking Isolation Before Prompt" << std::endl;
    for(JM::NavBuffer::Iterator it = pit - 1; it != m_buf->begin(); --it) {

        auto oechdr = JM::getHeaderObject<JM::OecHeader>(it->get());
        if(!oechdr) continue;
        JM::OecEvt* oecevt = dynamic_cast<JM::OecEvt*>(oechdr->event("JM::OecEvt"));
        if(!oecevt){
            LogInfo << "Enable to load OEC Event" << std::endl;
            continue;
        }

        const TTimeStamp& time = oecevt->getTime();
        double dt = (ptime.GetSec() - time.GetSec())*1000000000ULL + (ptime.GetNanoSec() - time.GetNanoSec());

        auto rechdr = JM::getHeaderObject<JM::CdVertexRecHeader>(it->get());
        if(!rechdr) continue;
        JM::CdVertexRecEvt* recevt = rechdr->event();

        const auto& vertex = recevt->getVertex(0);
        float energy = vertex->energy();

        if (dt*1e-6 > 1) break;
        else if(energy > DelayEnergyCut[0]) return false;
    }

    // Check between prompt and delay
    for(JM::NavBuffer::Iterator it = pit + 1; it != m_buf->end(); ++it) {
        auto rechdr = JM::getHeaderObject<JM::CdVertexRecHeader>(it->get());
        if(!rechdr) continue;
        JM::CdVertexRecEvt* recevt = rechdr->event();

        const auto& vertex = recevt->getVertex(0);
        float energy = vertex->energy();

        if (it == dit) break;
        else if(energy > DelayEnergyCut[0]) return false;
    }

    // Check after delay
    for(JM::NavBuffer::Iterator it = dit + 1; it != m_buf->end(); ++it) {
        auto oechdr = JM::getHeaderObject<JM::OecHeader>(it->get());
        if(!oechdr) continue;
        JM::OecEvt* oecevt = dynamic_cast<JM::OecEvt*>(oechdr->event("JM::OecEvt"));
        if(!oecevt){
            LogInfo << "Could not load OecEvt" << std::endl;
        }
        const TTimeStamp& time = oecevt->getTime();
        double dt = (time.GetSec() - dtime.GetSec())*1000000000ULL + (time.GetNanoSec() - dtime.GetNanoSec());

        auto rechdr = JM::getHeaderObject<JM::CdVertexRecHeader>(it->get());
        if(!rechdr) continue;
        JM::CdVertexRecEvt* recevt = rechdr->event();

        const auto& vertex = recevt->getVertex(0);
        float energy = vertex->energy();

        if (dt*1e-6 > 1) break;
        else if(energy > DelayEnergyCut[0]) return false;
    }

    return true;
}


bool IBDSelectionTool::isPrompt(JM::EvtNavigator* nav) {
    LogInfo << "IBD Search" << std::endl;
    if(!nav){
        LogInfo << "EvtNavigator not found" << std::endl;
        return false;
    }

    std::string lastMuTag = m_eventTagSvc->getLastMuTag();
    JM::EvtNavigator* lastMuNav = m_eventTagSvc->getLastMuNav();

    // Load EDM
        // OEC
    auto oechdr = JM::getHeaderObject<JM::OecHeader>(nav);
    if(!oechdr) return false;
    m_oecevt = dynamic_cast<JM::OecEvt*>(oechdr->event("JM::OecEvt"));
    if(!m_oecevt) return false;
    LogInfo << "Last Mu Tag: " << lastMuTag << std::endl;
    if(!lastMuTag.empty() && isVetoed(m_oecevt, lastMuNav)) return false;
    const TTimeStamp& ptime = m_oecevt->getTime();


        // Reco
    auto rechdr = JM::getHeaderObject<JM::CdVertexRecHeader>(nav, recEDMPath);
    if(!rechdr) return false;
    m_recevt = rechdr->event();

    const auto vertex = m_recevt->getVertex(0);
    pEnergy = vertex->energy();
    pCharge = vertex->peSum();
    pVertex.SetXYZ(vertex->x(), vertex->y(), vertex->z());

    bool energy_cut = pEnergy >= PromptEnergyCut[0] && pEnergy <= PromptEnergyCut[1];
    bool charge_cut = pCharge >= PromptChargeCut[0] && pCharge <= PromptChargeCut[1];
    bool position_cut = pVertex.Mag() < FV_cut && !(pVertex.Perp() < 2000 && std::abs(pVertex.Z() < 15500));

    if(energy_cut && position_cut /*&& charge_cut*/){
        LogInfo << "Prompt IBD Candidate: " << std::endl;
        LogInfo << " - energy: " << pEnergy << " r: " << pVertex.Mag() << std::endl;

        JM::NavBuffer::Iterator navit = m_buf->find(nav);
        int offset = 0;
        for(JM::NavBuffer::Iterator it = navit + 1; it != m_buf->end(); ++it) {

            JM::EvtNavigator* dnav = it->get();
            auto dOecHdr = JM::getHeaderObject<JM::OecHeader>(dnav);
            JM::OecEvt* dOecEvt = dynamic_cast<JM::OecEvt*>(dOecHdr->event("JM::OecEvt"));
            // if(isVetoed(dOecEvt, lastMuNav)) return false;
            const TTimeStamp& dtime = dOecEvt->getTime();

            double dt = (dtime.GetSec() - ptime.GetSec())*1000000000ULL + (dtime.GetNanoSec() - ptime.GetNanoSec());
            
            LogInfo << "dtime prompt/delay candidate: " << dt*1e-6 << " ms" << std::endl;
            if(dt*1e-6 > 1) return false; // if dt > 1ms no delay
            offset++;
            
            // Reco
            auto dRecHdr = JM::getHeaderObject<JM::CdVertexRecHeader>(dnav, recEDMPath);
            if(!dRecHdr){
                LogInfo << "No Rec header" << std::endl;
                continue;
            }
            JM::CdVertexRecEvt* dRecEvt = dRecHdr->event();
            
            const auto avertex = dRecEvt->getVertex(0);
            dEnergy = avertex->energy();
            dCharge = avertex->peSum();
            dVertex.SetXYZ(avertex->x(), avertex->y(), avertex->z());
            
            bool energy = dEnergy >= DelayEnergyCut[0] && dEnergy <= DelayEnergyCut[1];
            bool charge = dCharge >= DelayChargeCut[0] && dCharge <= DelayChargeCut[1];
            bool distance = (pVertex - dVertex).Mag() < 1500;
            LogInfo << "Energy: " << dEnergy << " distance: " << (pVertex - dVertex).Mag() << std::endl; 
                        
            if(dt*1e-3 > 5 && distance && energy && charge && isIsolated(nav, dnav, m_oecevt, dOecEvt)){
                LogInfo << "Delay found!" << std::endl;
                LogInfo << " - energy: " << dEnergy << " r: " << dVertex.Mag() << std::endl;

                m_eventTagSvc->addTag(nav, "Prompt");
                m_eventTagSvc->addTag(dnav, "Delay");

                DelayEntryOffset = offset;

                return true;
            }
        }
        
    }

    return false;

}