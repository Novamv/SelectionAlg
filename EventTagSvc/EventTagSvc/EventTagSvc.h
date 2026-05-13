#ifndef EVENTTAGSVC_HH
#define EVENTTAGSVC_HH

#include "SniperKernel/SvcBase.h"
#include "EvtNavigator/NavBuffer.h"

#include <unordered_map>
#include <string>
#include <vector>

#include "TTimeStamp.h"
#include "TVector3.h"

//----------------- List Of Tags -----------------
// 1. Muon Related: CDMuon, WPMuon, CDWPMuon
// 2. Neutron, SpalNeutron (Neutron is AfterPulse-like while SpalNeutron has the correct Hit time std and mean)
// 3. IBD related: Prompt, Delay
//------------------------------------------------


class EventTagSvc : public SvcBase {
    public:
        EventTagSvc(const std::string& name);
        ~EventTagSvc();

        virtual bool initialize();
        virtual bool finalize();

        void addTag(JM::EvtNavigator* nav, const std::string& tag);
        std::string getTag(JM::EvtNavigator* nav);
        bool hasTag(JM::EvtNavigator* nav);

        JM::EvtNavigator* getLastMuNav(){return m_LastMuNav;};
        std::string getLastMuTag(){return m_LastMuTag;};
        TTimeStamp getLastMuTime(){return m_LastMuTS;};

    private:
        // std::unordered_map<JM::EvtNavigator*, std::string> m_tags;
        std::unordered_map<uint64_t, std::string> m_tags;

        JM::EvtNavigator* m_LastNav {nullptr};
        std::string m_LastTag;
        JM::EvtNavigator* m_LastMuNav {nullptr};
        std::string m_LastMuTag;
        TTimeStamp m_LastMuTS = NULL;
};


#endif