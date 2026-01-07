#ifndef NEUTRONSPAL_H
#define NEUTRONSPAL_H

#include "SniperKernel/ToolBase.h"

#include "EventTagSvc/EventTagSvc.h"
#include "EvtNavigator/NavBuffer.h"
#include "EvtNavigator/EvtNavHelper.h"

#include "Event/OecHeader.h"
#include "Event/CdLpmtCalibHeader.h"
#include "Event/CdVertexRecHeader.h"

#include "TVector3.h"

class SpalNeutronSelectionTool: public ToolBase
{
private:
    JM::NavBuffer* m_buf;
    JM::CdLpmtCalibEvt* m_calibevt;
    JM::CdVertexRecEvt* m_recevt;
    JM::OecEvt* m_oecevt;

    EventTagSvc* m_eventTagSvc;

    std::string recEDMPath;
    float FV_cut;
    float NeutronEnergyCut[2];
    float NeutronChargeCut[2];

    TVector3 vecNeutron;


public:
    SpalNeutronSelectionTool(const std::string&);
    ~SpalNeutronSelectionTool();

    // Methods Inherited from ToolBase
    bool initialize();
    bool finalize();

    bool isSpalNeutron(JM::EvtNavigator*);
};



#endif