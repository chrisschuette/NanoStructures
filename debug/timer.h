#ifndef TIMER_H
#define TIMER_H

/**
 * @file
 *
 * @ingroup debug
 *
 * @brief Simple timer class
 */

/*
 * Copyright (c) 2014 Christoph Schuette.
 *
 * The license and distribution terms for this file may be
 * found in the file LICENSE in this distribution
 */

#include <sys/time.h>

namespace debug {

/**
 * \brief A simple timer class
 */
    class Timer {
    public:

        Timer() {
            gettimeofday(&m_start, NULL);
        }

        /**
         * \brief saves the current time as the start time.
         */
        void start() {
            gettimeofday(&m_start, NULL);
        }

        /**
         * \brief determines the elapsed time since start time.
         * @return elapsed time since start in seconds.
         */
        double elapsed() {
            gettimeofday(&m_end, NULL);
            return (m_end.tv_sec + m_end.tv_usec * 1e-6) -
                    (m_start.tv_sec + m_start.tv_usec * 1e-6);
        }
    protected:
        struct timeval m_start, m_end;

    };
}

#endif // TIMER_H
