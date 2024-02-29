#include "sph_int.h"
#include <queue>
#include <set>
#include <array>
#include <math.h>
#include <algorithm>
#include <assert.h>
#include <stdio.h>

struct gkp {
    double x, w1, w2;
};

std::array<gkp, 15> gk15 = {{
    {+0.991455371120813,+0.022935322010529,+0.000000000000000},
    {+0.949107912342759,+0.063092092629979,+0.129484966168870},
    {+0.864864423359769,+0.104790010322250,+0.000000000000000},
    {+0.741531185599394,+0.140653259715525,+0.279705391489277},
    {+0.586087235467691,+0.169004726639267,+0.000000000000000},
    {+0.405845151377397,+0.190350578064785,+0.381830050505119},
    {+0.207784955007898,+0.204432940075298,+0.000000000000000},
    {+0.000000000000000,+0.209482141084728,+0.417959183673469},
    {-0.207784955007898,+0.204432940075298,+0.000000000000000},
    {-0.405845151377397,+0.190350578064785,+0.381830050505119},
    {-0.586087235467691,+0.169004726639267,+0.000000000000000},
    {-0.741531185599394,+0.140653259715525,+0.279705391489277},
    {-0.864864423359769,+0.104790010322250,+0.000000000000000},
    {-0.949107912342759,+0.063092092629979,+0.129484966168870},
    {-0.991455371120813,+0.022935322010529,+0.000000000000000}
}};

struct inter {
    double a;
    double b;
    double I;
    double eps;
    double err1;
    double err2;
    double err;
    bool operator<(const inter& other) const {
        if (err > other.err) return true;
        if (err < other.err) return false;
        if (a > other.a) return true;
        if (a < other.a) return false;
        if (b > other.b) return true;
        if (b < other.b) return false;
        return false;
    }
};

struct thr {
    double x;
    bool start;
    bool operator<(const thr& other) const {
        return x < other.x;
    }
};

retval integrate(const fun1d& f, double a_, double b_, double eps, bool verbose) {
    if (verbose) printf("verbose integrate %lf %lf\n",a_,b_);
    if (fabs(b_-a_) < eps) return {0,0};
    std::set< inter > q;
    double tot_err=0;
    auto add = [&](double a, double b, double eps_) {
        double ab = (a+b)/2;
        double d = (b-a)/2;
        double I1=0, I2=0, err1=0;
        for (const auto& p : gk15) {
            retval rv = f(ab+p.x*d, eps_);
            I1 += rv.val*p.w1;
            I2 += rv.val*p.w2;
            err1 += rv.err*p.w2;
        }
        I1 = I1 * d;
        I2 = I2 * d;
        double err2 = fabs(I1-I2);
        double err = err1 + err2;
        //printf("added off: %le\n", err);
        tot_err += err;
        q.insert({a,b,I2,eps_,err1,err2,err});
    };
    add(a_,b_,eps);
    if (verbose) printf("start error: %le\n", tot_err);
    while (tot_err > eps) {
        if (verbose) printf("current: %.16le (%ld)\n", tot_err,q.size());
        auto it = q.begin();
        inter v = *it;
        q.erase(it);
        tot_err -= v.err;
        if (verbose) printf("errs: %.16le %.16le %.16le, eps: %.16le\n", v.err, v.err1, v.err2, v.eps);
        if (v.err2 > v.err1) {
            if (verbose) printf("Halving interval\n");
            double ab = (v.a+v.b)/2;
            add(v.a, ab, v.eps*2);
            add(ab, v.b, v.eps*2);
        } else {
            if (verbose) printf("Halving eps\n");
            add(v.a, v.b, v.eps/2);
        }
        //if (q.size() > 300) break;
    };
    double val = 0;
    for (const inter& v : q) val += v.I;
    return {val, tot_err};
}

sphere_set::sphere_set(const std::vector<sph>& spheres_)
    : spheres(spheres_) {
        assert(spheres.size() <= MAXN);
};

retval sphere_set::seg_fun(const fun3d& f, double x, double y, double z, double eps) const {
    return f(x,y,z,eps);
}

retval sphere_set::circ_fun(const fun3d& f, double x, double y, double eps) const {
    std::array<thr, 2*MAXN> t;
    auto it = t.begin();
    for (const auto& s : spheres) {
        double b2 = s.r*s.r-(s.x-x)*(s.x-x)-(s.y-y)*(s.y-y);
        if (b2 < 0) continue;
        double b = sqrt(b2);
        *it = {s.z-b, true}; it++;
        *it = {s.z+b, false}; it++;
    }
    std::sort(t.begin(),it);
    retval val = {0.0,0.0};
    int in = 0;
    double a,b;
    for (auto it2 = t.begin(); it2 != it; it2++) {
        if (in == 0) a = it2->x;
        if (it2->start) in++; else in--;
        if (in == 0) {
            b = it2->x;
            val += integrate([&](double z, double eps_) { return seg_fun(f,x,y,z,eps_); },a,b,eps);
        }
    }
    return val;
}

retval sphere_set::sphere_fun(const fun3d& f, double x, double eps) const {
    std::array<thr, 2*MAXN> t;
    auto it = t.begin();
    for (const auto& s : spheres) {
        double b2 = s.r*s.r-(s.x-x)*(s.x-x);
        if (b2 < 0) continue;
        double b = sqrt(b2);
        *it = {s.y-b, true}; it++;
        *it = {s.y+b, false}; it++;
    }
    std::sort(t.begin(),it);
    retval val = {0.0,0.0};
    int in = 0;
    double a,b;
    for (auto it2 = t.begin(); it2 != it; it2++) {
        if (in == 0) a = it2->x;
        if (it2->start) in++; else in--;
        if (in == 0) {
            b = it2->x;
            //printf("%lf - %lf - %lf\n", x, a, b); fflush(stdout);
            val += integrate([&](double y, double eps_) { return circ_fun(f,x,y,eps_); },a,b,eps);
        }
    }
    return val;
}

retval sphere_set::integrate_fun(const fun3d& f, double eps) const {
    std::array<thr, 2*MAXN> t;
    auto it = t.begin();
    for (const auto& s : spheres) {
        double b = s.r;
        *it = {s.x-b, true}; it++;
        *it = {s.x+b, false}; it++;
    }
    std::sort(t.begin(),it);
    retval val = {0.0,0.0};
    int in = 0;
    double a,b;
    for (auto it2 = t.begin(); it2 != it; it2++) {
        if (in == 0) a = it2->x;
        if (it2->start) in++; else in--;
        if (in == 0) {
            b = it2->x;
            val += integrate([&](double x, double eps_) { return sphere_fun(f,x,eps_); },a,b,eps);
        }
    }
    return val;
}
