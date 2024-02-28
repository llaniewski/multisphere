#include <algorithm>
#include <vector>
#include <functional>

#define FUN3D(x__) [&](double x, double y, double z, double eps) -> retval { return {x__,0.0};}

const int MAXN = 10;

struct sph {
    double x,y,z,r;
};

struct retval {
    double val;
    double err;
    inline retval operator+(const retval& other) const {
        return {val+other.val, err+other.err};
    }
    inline retval& operator+=(const retval& other) {
        val += other.val;
        err += other.err;
        return *this;
    }
    inline operator double () {
        return val;
    }
};

typedef std::function< retval(double,double,double,double) > fun3d;
typedef std::function< retval(double,double) > fun1d;

retval integrate(const fun1d& f, double a_, double b_, double eps, bool verbose=false);

class sphere_set {
    const std::vector<sph>& spheres;
    retval seg_fun(const fun3d& f, double x, double y, double z, double eps) const;
    retval circ_fun(const fun3d& f, double x, double y, double eps) const;
    retval sphere_fun(const fun3d& f, double x, double eps) const;
public:
    sphere_set(const std::vector<sph>& spheres_);
    retval integrate_fun(const fun3d& f, double eps=1e-6) const;
};

