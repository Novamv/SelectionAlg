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


bool IBDSelectionTool::isVetoed(JM::OecEvt* pOecEvt) {
    
    LogInfo << "Checking Muon veto" << std::endl;
    const TTimeStamp& ptime = pOecEvt->getTime();
    const TTimeStamp& muTime = m_eventTagSvc->getLastMuTime();

    int64_t dtime = deltaT_ns(ptime, muTime);

    LogInfo << "Last Muon dtime from prompt candidate: " << dtime*1e-6 << " ms" << std::endl;

    if(dtime < 0) return false;
    return dtime * 1e-6 <= 5.0;
}


bool IBDSelectionTool::isIsolated(JM::EvtNavigator* pnav, JM::EvtNavigator* dnav, JM::OecEvt* pOecEvt, JM::OecEvt* dOecEvt){

    LogInfo << "Prompt nav " << pnav << " Delay Nav " << dnav << std::endl;

    JM::NavBuffer::Iterator pit = m_buf->find(pnav);
    JM::NavBuffer::Iterator dit = m_buf->find(dnav);

    const TTimeStamp& ptime = pOecEvt->getTime();
    const TTimeStamp& dtime = dOecEvt->getTime();

    // Check before Prompt
    if(pit != m_buf->begin()){
        LogInfo << "Checking Isolation Before Prompt" << std::endl;
        for(JM::NavBuffer::Iterator it = pit - 1; it != m_buf->begin(); --it) {
    
            auto oechdr = JM::getHeaderObject<JM::OecHeader>(it->get());
            if(!oechdr) {
                LogInfo << "No Oec Header" << std::endl;
                continue;
            }
            JM::OecEvt* oecevt = dynamic_cast<JM::OecEvt*>(oechdr->event("JM::OecEvt"));
            if(!oecevt){
                LogInfo << "Unable to load OEC Event" << std::endl;
                continue;
            }
    
            const TTimeStamp& time = oecevt->getTime();
            int64_t dt = deltaT_ns(ptime, time);
    
            auto rechdr = JM::getHeaderObject<JM::CdVertexRecHeader>(it->get(), recEDMPath);
            if(!rechdr) {
                LogInfo << "No Rec Header" << std::endl;
                continue;
            }
            JM::CdVertexRecEvt* recevt = rechdr->event();
            if(!recevt || recevt->nVertices() == 0) continue;
    
            const auto& vertex = recevt->getVertex(0);
            float energy = vertex->energy();
    
            bool inEnergyRange = energy >= DelayEnergyCut[0] && energy <= PromptEnergyCut[1];
            bool isdt = dt*1e-6 > 1;
            LogInfo << "Energy: " << energy << " dtime: " << dt*1e-6 << std::endl;
            LogInfo << "Is in Energy range: " << inEnergyRange << "; Is in dt range: " << !isdt << std::endl;

            if (isdt) break;
            else if(inEnergyRange && !isVetoed(oecevt)) {
                LogInfo << "Multiplicity before!" << std::endl;
                return false;
            }
        }
    }

    // Check between prompt and delay
    LogInfo << "Checking Isolation Between Prompt & Delay" << std::endl;
    for(JM::NavBuffer::Iterator it = pit + 1; it != m_buf->end(); it++) {
        
        if (it == dit) break; //Break if reached delay iterator

        auto rechdr = JM::getHeaderObject<JM::CdVertexRecHeader>(it->get(), recEDMPath);
        if(!rechdr){
            LogInfo << "No Rec Header" << std::endl;
            continue;
        }

        JM::CdVertexRecEvt* recevt = rechdr->event();
        if(!recevt || recevt->nVertices() == 0) continue;

        const auto& vertex = recevt->getVertex(0);
        float energy = vertex->energy();
        
        LogInfo << "Energy: " << energy << std::endl;
        bool inEnergyRange = energy >= DelayEnergyCut[0] && energy <= PromptEnergyCut[1];
        LogInfo << "Is in Energy range: " << inEnergyRange << std::endl;

        
        if(inEnergyRange){
            LogInfo << "Multiplicity in between!" << std::endl;
            return false;
        }
    }

    // Check after delay
    LogInfo << "Checking After Delay" << std::endl;
    if(dit != m_buf->end()){
        for(JM::NavBuffer::Iterator it = dit + 1; it != m_buf->end(); it++) {
            auto oechdr = JM::getHeaderObject<JM::OecHeader>(it->get());
            if(!oechdr) continue;
            JM::OecEvt* oecevt = dynamic_cast<JM::OecEvt*>(oechdr->event("JM::OecEvt"));
            if(!oecevt){
                LogInfo << "Could not load OecEvt" << std::endl;
                continue;
            }
            const TTimeStamp& time = oecevt->getTime();
            int64_t dt = deltaT_ns(time, dtime);
    
            auto rechdr = JM::getHeaderObject<JM::CdVertexRecHeader>(it->get(), recEDMPath);
            if(!rechdr) continue;
            JM::CdVertexRecEvt* recevt = rechdr->event();
            if(!recevt || recevt->nVertices() == 0) continue;
    
            const auto& vertex = recevt->getVertex(0);
            float energy = vertex->energy();

            bool inEnergyRange = energy >= DelayEnergyCut[0] && energy <= PromptEnergyCut[1];
            bool isdt = dt*1e-6 > 1.0;
            
            LogInfo << "Energy: " << energy << " dtime: " << dt*1e-6 << std::endl;
            LogInfo << "Is in Energy range: " << inEnergyRange << "; Is in dt range: " << !isdt << std::endl;
            
            if (dt*1e-6 > 1.0) break;
            else if(inEnergyRange && !isVetoed(oecevt)){
                LogInfo << "Multiplicity after!" << std::endl;
                return false;
            }
        }
    }

    return true;
}


bool IBDSelectionTool::isPrompt(JM::EvtNavigator* nav){
    LogInfo << "IBD Search" << std::endl;
    if(!nav){
        LogInfo << "EvtNavigator not found" << std::endl;
        return false;
    }

    std::string lastMuTag = m_eventTagSvc->getLastMuTag();
    const TTimeStamp lastMuTime = m_eventTagSvc->getLastMuTime();

    // Load EDM
        // OEC
    auto oechdr = JM::getHeaderObject<JM::OecHeader>(nav);
    if(!oechdr) return false;
    m_oecevt = dynamic_cast<JM::OecEvt*>(oechdr->event("JM::OecEvt"));
    if(!m_oecevt) return false;
    
    const TTimeStamp& ptime = m_oecevt->getTime();
    int64_t Mudtime = deltaT_ns(ptime, lastMuTime);

    LogInfo << "Last Mu Tag: " << lastMuTag << std::endl;
    LogInfo << "Prompt dtime from Muon " << Mudtime << std::endl;

    if(!lastMuTag.empty() && isVetoed(m_oecevt)) return false;

        // Reco
    auto rechdr = JM::getHeaderObject<JM::CdVertexRecHeader>(nav, recEDMPath);
    if(!rechdr) return false;
    m_recevt = rechdr->event();
    if(!m_recevt || m_recevt->nVertices() == 0) return false;


    const auto vertex = m_recevt->getVertex(0);
    pEnergy = vertex->energy();
    pCharge = vertex->peSum();
    pVertex.SetXYZ(vertex->x(), vertex->y(), vertex->z());

    bool energy_cut = pEnergy > PromptEnergyCut[0] && pEnergy < PromptEnergyCut[1];
    bool charge_cut = pCharge > PromptChargeCut[0] && pCharge < PromptChargeCut[1];
    bool position_cut = pVertex.Mag() < FV_cut && !(pVertex.Perp() < 2000 && std::abs(pVertex.Z()) < 15500);

    LogInfo << "Prompt IBD Candidate: " << nav << std::endl;
    LogInfo << " - energy: " << pEnergy << " r: " << pVertex.Mag() << " z: " << pVertex.Z() << " rho: " << pVertex.Perp() << std::endl;

    JM::NavBuffer::Iterator navit = m_buf->find(nav);

    if(energy_cut && position_cut && navit != m_buf->end() /*&& charge_cut*/){

        int offset = 0;
        for(JM::NavBuffer::Iterator it = navit + 1; it != m_buf->end(); ++it) {

            JM::EvtNavigator* dnav = it->get();

            auto dOecHdr = JM::getHeaderObject<JM::OecHeader>(dnav);
            if(!dOecHdr) continue;
            JM::OecEvt* dOecEvt = dynamic_cast<JM::OecEvt*>(dOecHdr->event("JM::OecEvt"));
            if(!dOecEvt) continue;
            
            const TTimeStamp& dtime = dOecEvt->getTime();

            int64_t dt = deltaT_ns(dtime, ptime);
            
            LogInfo << "dtime prompt/delay candidate: " << dnav << " dtime: " << dt*1e-6 << " ms" << std::endl;
            // pair them if 5us < dt < 1ms
            if(dt*1e-3 <= 5) continue;
            if(dt*1e-6 > 1) return false; // if dt > 1ms no delay
            offset++;
            
            // Reco
            auto dRecHdr = JM::getHeaderObject<JM::CdVertexRecHeader>(dnav, recEDMPath);
            if(!dRecHdr){
                LogInfo << "No Rec header" << std::endl;
                continue;
            }
            JM::CdVertexRecEvt* dRecEvt = dRecHdr->event();
            if(!dRecEvt || dRecEvt->nVertices() == 0) continue;
            
            const auto avertex = dRecEvt->getVertex(0);
            dEnergy = avertex->energy();
            dCharge = avertex->peSum();
            dVertex.SetXYZ(avertex->x(), avertex->y(), avertex->z());
            
            bool energy = dEnergy > DelayEnergyCut[0] && dEnergy < DelayEnergyCut[1];
            bool charge = dCharge > DelayChargeCut[0] && dCharge < DelayChargeCut[1];
            bool distance = (pVertex - dVertex).Mag() < 1500;
            LogInfo << "Energy: " << dEnergy << " Charge: " << dCharge << " r: " << dVertex.Mag() << " Distance: " << (pVertex - dVertex).Mag() << std::endl; 
                        
            if(dt*1e-3 > 5 && distance && energy /*&& charge*/ && isIsolated(nav, dnav, m_oecevt, dOecEvt)){
                LogInfo << "Delay found!" << std::endl;
                LogInfo << " - energy: " << dEnergy << " r: " << dVertex.Mag() << std::endl;

                m_eventTagSvc->addTag(nav, "Prompt");
                m_eventTagSvc->addTag(dnav, "Delay");

                DelayEntryOffset = offset;
                LogInfo << "Delay Entry Offset: " << DelayEntryOffset << std::endl;

                return true;
            }
        }
        
    }

    return false;
}

bool IBDSelectionTool::isNeutronVetoed(JM::EvtNavigator* nav){
    LogInfo << "Prompt Neutron Veto search" << std::endl;

    TTimeStamp pTime = nav->TimeStamp();

    auto prechdr = JM::getHeaderObject<JM::CdVertexRecHeader>(nav, recEDMPath);
    if(!prechdr) return false;
    m_recevt = prechdr->event();
    if(!m_recevt || m_recevt->nVertices() == 0) return false;

    const auto vertex = m_recevt->getVertex(0); // Assumes 1 vertex per event (might be false)
    TVector3 p_vertex(vertex->x(), vertex->y(), vertex->z());

    JM::NavBuffer::Iterator navit = m_buf->find(nav);
    if(navit == m_buf->end() || navit == m_buf->begin()) return false;

    for(JM::NavBuffer::Iterator it = navit - 1; it != m_buf->begin(); --it) {
        
        JM::EvtNavigator* candidate = it->get();
        
        TTimeStamp nTime = candidate->TimeStamp();
        int64_t dtime = deltaT_ns(pTime, nTime);
        
        if(dtime < 0) continue;
        if(dtime >= 1'200'000'000LL) break; // beyond 1.2 s, stop looking

        std::string tag = m_eventTagSvc->getTag(candidate);
        if(tag != "SpalNeutron") continue;

        //Load neutron vertex
        auto nrechdr = JM::getHeaderObject<JM::CdVertexRecHeader>(candidate, recEDMPath);
        if(!nrechdr) return false;
        JM::CdVertexRecEvt* nrecevt = nrechdr->event();
        if(!nrecevt || nrecevt->nVertices() == 0) continue;

        const auto nVtx = nrecevt->getVertex(0); 
        TVector3 n_vertex(nVtx->x(), nVtx->y(), nVtx->z());
        double dR = (p_vertex - n_vertex).Mag();

        if(dR < 4000) return true;   
    }
    return false;
}


