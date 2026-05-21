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
// 2. Neutron, SpalNeutron                      (Neutron is AfterPulse-like while SpalNeutron has the correct Hit time std and mean)
// 3. IBD related: Prompt, Delay
//------------------------------------------------


class EventTagSvc : public SvcBase {
    public:
        EventTagSvc(const std::string& name);
        ~EventTagSvc();

        virtual bool initialize();
        virtual bool finalize();

        // ── Tag Handler ─────────────────────────────────────────
        void addTag(JM::EvtNavigator* nav, const std::string& tag);
        std::string getTag(JM::EvtNavigator* nav);
        bool hasTag(JM::EvtNavigator* nav);

        // ── Fast access to Muon information ─────────────────────────────────────────
        std::string getLastMuTag(){return m_LastMuTag;};
        TTimeStamp getLastMuTime(){return m_LastMuTS;};

    private:
        std::unordered_map<uint64_t, std::string> m_tags;

        TTimeStamp m_LastMuTS{0, 0};
        std::string m_LastTag;
        std::string m_LastMuTag;

        // ── Update Tag list ─────────────────────────────────────────
        void CleanupTags(uint64_t currentTime);
};


#endif