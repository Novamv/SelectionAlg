#ifndef IBDSELECTION_H
#define IBDSELECTION_H

#include "SniperKernel/ToolBase.h"

#include "EventTagSvc/EventTagSvc.h"
#include "EvtNavigator/NavBuffer.h"
#include "EvtNavigator/EvtNavHelper.h"

#include "Event/OecHeader.h"
#include "Event/CdLpmtCalibHeader.h"
#include "Event/CdVertexRecHeader.h"

#include "TVector3.h"

class IBDSelectionTool: public ToolBase
{
private:
    // EDM
    JM::NavBuffer* m_buf;
    JM::OecEvt* m_oecevt;
    JM::CdVertexRecEvt* m_recevt;

    EventTagSvc* m_eventTagSvc;

    // Property
    std::string recEDMPath;
    float FV_cut;
    float PromptEnergyCut[2];
    float DelayEnergyCut[2];
    int PromptChargeCut[2];
    int DelayChargeCut[2];

    // Variables
    float pEnergy, dEnergy;
    TVector3 pVertex, dVertex;
    int DelayEntryOffset;

    // Methods
    bool isVetoed(JM::OecEvt*, JM::EvtNavigator*);
    bool isIsolated(JM::EvtNavigator*, JM::EvtNavigator*, JM::OecEvt*, JM::OecEvt*);
    
    public:
    IBDSelectionTool(const std::string& name);
    ~IBDSelectionTool();
    
    bool initialize();
    bool finalize();
    
    bool isPrompt(JM::EvtNavigator*);
    int getDelayOffset(){return DelayEntryOffset;};

};


#endif