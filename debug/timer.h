#ifndef TIMER_H
#define TIMER_H

#include <sys/time.h>

namespace debug {
class Timer
{
public:
    Timer() { gettimeofday(&m_start, NULL); }
    void start()  { gettimeofday(&m_start, NULL); }
    double elapsed() { gettimeofday(&m_end, NULL);
                       return   (m_end.tv_sec +   m_end.tv_usec * 1e-6) -
                                        (m_start.tv_sec + m_start.tv_usec*1e-6); }
protected:
    struct timeval m_start, m_end;

};
}

#endif // TIMER_H
