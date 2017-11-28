#include "signal.hpp" 
#include "sismograph.hpp"

double pi=2*asin(1.);

signal<double> hamming(size_t n){
  auto N=2*n+1;
  signal<double> v(N,0.,0.001);
  
    for(auto i=0U;i<=n;i++)
      v[i]=0.54-(0.46*cos(pi*double(i)/double(n)));
    for(auto i=0U;i<=n;i++)
      v[n+i]=v[n-i];
  
    return v;
};
using namespace std;
int main() {
  
  auto t=linspace(0.,10.,1000);
  signal<double> y(0.01);
  y.read("../lsst07/fa1-1.7e"s);
  y.write("y.data");
  y.demean( );
  y.write("yorig.data");
  y.energy().write("yy.data");
  y.extract_with_energy_bound(0.05,0.95).write("extract_energy.data");
  return 0;

}
