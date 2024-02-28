#include <Rcpp.h>
#include "sph_int.h"

using namespace Rcpp;

// [[Rcpp::export]]
List sphere_integrals(DataFrame tab) {
  List ret;
  if ((!tab.containsElementNamed("x")) || 
      (!tab.containsElementNamed("x")) ||
      (!tab.containsElementNamed("x")) ||
      (!tab.containsElementNamed("x"))) {
    return ret;
  }
  std::vector<sph> spheres;
  NumericVector x = tab["x"];
  NumericVector y = tab["y"];
  NumericVector z = tab["z"];
  NumericVector r = tab["r"];
  for (int i=0; i<tab.nrow(); i++) {
    spheres.push_back({x[i],y[i],z[i],r[i]});
  }
  sphere_set myset(spheres);
  double vol = myset.integrate_fun(FUN3D(1));
  ret["volume"] = vol;
  double mx = myset.integrate_fun(FUN3D(x));
  double my = myset.integrate_fun(FUN3D(y));
  double mz = myset.integrate_fun(FUN3D(z));
  ret["center"] = NumericVector({mx,my,mz});
  double Ixx = myset.integrate_fun(FUN3D((x-mx)*(x-mx)));
  double Ixy = myset.integrate_fun(FUN3D((x-mx)*(y-my)));
  double Ixz = myset.integrate_fun(FUN3D((x-mx)*(z-mz)));
  double Iyy = myset.integrate_fun(FUN3D((y-my)*(y-my)));
  double Iyz = myset.integrate_fun(FUN3D((y-my)*(z-mz)));
  double Izz = myset.integrate_fun(FUN3D((z-mz)*(z-mz)));
  NumericVector second({Ixx,Ixy,Ixz,Ixy,Iyy,Iyz,Ixz,Iyz,Izz});
  ret["second"] = NumericMatrix(3,3,second.begin());
  return ret;
}



