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

void EventTagSvc::addTag(JM::EvtNavigator* nav, const std::string& tag) {
    const TTimeStamp& TS(nav->TimeStamp());
    uint64_t ttime = TS.GetSec()*1000000000ULL + TS.GetNanoSec();
    m_tags[ttime] = tag;

    if(tag=="CDMuon" || tag=="WPMuon" || tag=="CDWPMuon"){
        m_LastMuTag = tag;
        m_LastMuNav = nav;
    }
    else{
        m_LastTag = tag;
        m_LastNav = nav;
    }
}

std::string EventTagSvc::getTag(JM::EvtNavigator* nav) {
    const TTimeStamp& TS(nav->TimeStamp());
    uint64_t ttime = TS.GetSec()*1000000000ULL + TS.GetNanoSec();
    auto it = m_tags.find(ttime);
    if (it != m_tags.end()) {
        return it->second;
    }
    else{
        LogInfo << "Tag not found." << std::endl;
    }
    return "";
}

bool EventTagSvc::hasTag(JM::EvtNavigator* nav){
    const TTimeStamp& TS(nav->TimeStamp());
    uint64_t ttime = TS.GetSec()*1000000000ULL + TS.GetNanoSec();
    return m_tags.find(ttime) != m_tags.end();
}
