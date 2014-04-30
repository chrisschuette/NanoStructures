#include "utils.h"
#include <algorithm>
#include <time.h>       /* time_t, struct tm, time, localtime, asctime */

namespace utils {
std::string timestamp() {
    time_t rawtime;
    struct tm * timeinfo;

    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    return std::string(asctime (timeinfo));
}

std::string toUpperCase(std::string strToConvert) {
    std::transform(strToConvert.begin(), strToConvert.end(), strToConvert.begin(), ::toupper);
    return strToConvert;
}
}
