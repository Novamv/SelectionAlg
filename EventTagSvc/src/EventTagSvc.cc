#include "EventTagSvc/EventTagSvc.h"
#include "SniperKernel/SvcFactory.h"
#include "SniperKernel/SniperPtr.h"


DECLARE_SERVICE(EventTagSvc);

EventTagSvc::EventTagSvc(const std::string& name)
    : SvcBase(name) 
{
    // Add property if needed
}

EventTagSvc::~EventTagSvc() {
}

bool EventTagSvc::initialize() {
    LogInfo << "EventTagSvc initialized." << std::endl;
    return true;
}

bool EventTagSvc::finalize() {
    LogInfo << "EventTagSvc finalized." << std::endl;
    return true;
}

void EventTagSvc::CleanupTags(uint64_t currentTime){
    const uint64_t maxAge = 2'000'000'000ULL; // 2 seconds

    while(!m_tagTimes.empty() && currentTime - m_tagTimes.front() > maxAge){
        m_tags.erase(m_tagTimes.front());
        m_tagTimes.pop_front();
    }
}

void EventTagSvc::addTag(JM::EvtNavigator* nav, const std::string& tag) {
    const TTimeStamp& TS(nav->TimeStamp());
    uint64_t ttime = TS.GetSec()*1000000000ULL + TS.GetNanoSec();

    if(m_tags.find(ttime) == m_tags.end()) m_tagTimes.push_back(ttime); // only push if new time
    m_tags[ttime] = tag;

    CleanupTags(ttime);

    if(tag=="CDMuon" || tag=="WPMuon" || tag=="CDWPMuon"){
        m_LastMuTag = tag;
        m_LastMuTS = TS;
    }
    else{
        m_LastTag = tag;
    }
}

std::string EventTagSvc::getTag(JM::EvtNavigator* nav) {
    const TTimeStamp& TS(nav->TimeStamp());
    uint64_t ttime = TS.GetSec()*1000000000ULL + TS.GetNanoSec();
    
    CleanupTags(ttime);

    auto it = m_tags.find(ttime);
    if (it != m_tags.end()) return it->second;
    
    LogInfo << "Tag not found." << std::endl;
    
    return "";
}

bool EventTagSvc::hasTag(JM::EvtNavigator* nav){
    const TTimeStamp& TS(nav->TimeStamp());
    uint64_t ttime = TS.GetSec()*1000000000ULL + TS.GetNanoSec();
    return m_tags.find(ttime) != m_tags.end();
}
