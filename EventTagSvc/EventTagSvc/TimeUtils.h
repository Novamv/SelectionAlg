#pragma once
#include "TTimeStamp.h"
#include <cstdint>

inline int64_t toKey(const TTimeStamp& ts) { 
    return (int64_t)(ts.GetSec() * 1000000000LL + ts.GetNanoSec());
}

inline int64_t deltaT_ns(const TTimeStamp& a, const TTimeStamp& b) {
    return (int64_t)(a.GetSec()    - b.GetSec())    * 1000000000LL
         + (int64_t)(a.GetNanoSec() - b.GetNanoSec());
}

inline double deltaT_us(const TTimeStamp& a, const TTimeStamp& b) {
    return deltaT_ns(a, b) * 1e-3;
}

inline double deltaT_ms(const TTimeStamp& a, const TTimeStamp& b) {
    return deltaT_ns(a, b) * 1e-6;
}

inline double deltaT_s(const TTimeStamp& a, const TTimeStamp& b) {
    return deltaT_ns(a, b) * 1e-9;
}