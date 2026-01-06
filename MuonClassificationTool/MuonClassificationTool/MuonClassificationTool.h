#ifndef MUON_HH
#define MUON_HH

#include "SniperKernel/ToolBase.h"

// #include "MyClassification/src/IMyClassification.h"
#include "EventTagSvc/EventTagSvc.h"

#include "EvtNavigator/NavBuffer.h"
#include "EvtNavigator/EvtNavigator.h"
#include "EvtNavigator/EvtNavHelper.h"

#include "Event/CdLpmtCalibHeader.h"
#include "Event/WpCalibHeader.h"
#include "Event/OecHeader.h"


class MuonClassificationTool: public ToolBase /*, public IMyClassificationTool */
{
    private:
        JM::NavBuffer* m_buf;
        JM::CdLpmtCalibEvt* m_calibevt;
        JM::WpCalibEvt* m_wpcalibevt;

        EventTagSvc* m_eventTagSvc;

        float WPMuonChargeCut = 700; //PE
        float CDMuonChargeCut = 30000; //PE

        bool isInDeadTime(JM::EvtNavigator*, JM::OecEvt*, const std::string&);
        bool isVetoed(JM::EvtNavigator*, JM::OecEvt*);
        
    public:
        MuonClassificationTool(const std::string&);
        ~MuonClassificationTool();
        
        //inherited from ToolBase
        virtual bool initialize();
        virtual bool finalize();
        
        //inherited from IMyClassificationTool
        // virtual MyClassificationType classify(JM::EvtNavigator*);

        bool isMuon(JM::EvtNavigator*);
};  


#endif
