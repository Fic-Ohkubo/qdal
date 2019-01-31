
#include "../low/chrono.hpp"
#include <cmath>
#include <string>
using namespace std;

const int USECSPERSEC = 1000000;

Chronometer::Chronometer() {
    total_secs = 0;
    total_usecs = 0;
    running = false;
}

void Chronometer::Start() {
    gettimeofday(&this->timer[0], NULL);
    running = true;
}

void Chronometer::Stop() {
    gettimeofday(&this->timer[1], NULL);
    running = false;

    int secs (this->timer[1].tv_sec  - this->timer[0].tv_sec);
    int usecs(this->timer[1].tv_usec - this->timer[0].tv_usec);

    total_secs  +=  secs;
    total_usecs += usecs;
}

void Chronometer::Restart() {
    total_secs = 0;
    total_usecs = 0;
}

double Chronometer::GetTotalTime() const {
    double sec = double(total_usecs);
    sec /= 1000000.;
    sec += double(total_secs);

    return sec;
}

ostream & operator<<(ostream & os, const Chronometer & chrono) {

    int total_usecs = chrono.total_usecs;
    int total_secs  = chrono.total_secs;

    while (total_usecs < 0) {
        total_usecs += USECSPERSEC;
        total_secs  -= 1;
    }
    while (total_usecs >= USECSPERSEC) {
        total_usecs -= USECSPERSEC;
        total_secs  += 1;
    }


    int tot = total_secs;

    int s = tot%60;
    tot /= 60;

    int m = tot%60;
    tot /= 60;

    int h = tot%24;
    tot /= 24;

    int d = tot;

    if (d>0) os << d << " days ";
    if (h>0) os << h << " hours ";
    if (m>0) os << m << " minutes ";
    if (s>0) os << s << " seconds ";
    //os.precision(3);
    os << total_usecs << " microseconds";

    return os;
}
