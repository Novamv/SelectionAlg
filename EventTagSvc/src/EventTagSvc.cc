#include "EventTagSvc/EventTagSvc.h"
#include "SniperKernel/SvcFactory.h"
#include "SniperKernel/SniperPtr.h"

#include "EventTagSvc/TimeUtils.h"

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

void EventTagSvc::CleanupTags(const TTimeStamp& currentTime){

    while(!m_tagTimes.empty()){
        int64_t age = deltaT_ns(currentTime, m_tagTimes.front());
        if(age < 0){
            m_tagTimes.pop_front();
            continue;
        }
        if(age <= 2'000'000'000LL) break; // if dt < 2s don't clean

        m_tags.erase(toKey(m_tagTimes.front()));
        m_tagTimes.pop_front();
    }
}

void EventTagSvc::addTag(JM::EvtNavigator* nav, const std::string& tag) {
    const TTimeStamp& TS(nav->TimeStamp());
    int64_t key = toKey(TS);

    if(m_tags.find(key) == m_tags.end()) m_tagTimes.push_back(TS); // only push if new time
    m_tags[key] = tag;

    CleanupTags(TS);

    if(tag == "CDMuon" || tag == "CDWPMuon") m_LastCDMuTS = TS;
    if(tag == "WPMuon" || tag == "CDWPMuon") m_LastWPMuTS = TS;
    if(tag == "CDMuon" || tag == "WPMuon" || tag == "CDWPMuon"){
        m_LastMuTag = tag;
        m_LastMuTS = TS;
    }
    else{
        m_LastTag = tag;
    }
}

std::string EventTagSvc::getTag(JM::EvtNavigator* nav) {
    const TTimeStamp& TS(nav->TimeStamp());
    int64_t key = toKey(TS);

    CleanupTags(TS);

    auto it = m_tags.find(key);
    if (it != m_tags.end()) return it->second;
    
    LogInfo << "Tag not found." << std::endl;
    return "";
}

bool EventTagSvc::hasTag(JM::EvtNavigator* nav){
    const TTimeStamp& TS(nav->TimeStamp());
    int64_t key = toKey(TS);
    return m_tags.find(key) != m_tags.end();
}
