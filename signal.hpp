#pragma once
#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <complex>
#include <cmath>
#include <cassert>
#include <random>
#include <array>
#include <unsupported/Eigen/FFT>

template<class T >
struct realof{
  using type=T;
};
template<class T>
struct realof<std::complex<T>>{
  using type=  T;
};
////////////////////////////////////////////////////////////////////////////////////////
template<class T>
struct is_complex{
  static constexpr bool value=false;
};
template<class T>
struct is_complex<std::complex<T>>{
  static constexpr bool value=true;
};
////////////////////////////////////////////////////////////////////////////////////////
template<class T>
auto print(T v) ->std::enable_if_t<std::is_arithmetic<T>::value,void> {
  std::cout<<v<<std::endl;
};
template<class T>
auto print(T v) ->std::enable_if_t<is_complex<T>::value,void> {
  std::cout<<v.real()<<'\t'<<v.imag()<<std::endl;
  //  std::cout<<v<<std::endl;
};
template<class T>
auto print(T v) ->std::enable_if<std::is_same<typename T::iterator,typename T::iterator>::value,void > {
  for(auto el:v)
    print(el);
};

/////////////////////////////////////////////////////////////////////////////////////////
enum class method{biased,unbiased};
/////////////////////////////////////////////////////////////////////////////////////////
template<class T=double>
class signal : public std::vector<T> {
  double _dt=0.;
  double _t0=0;
  using Scalar=typename realof<T>::type;
  using Complex=std::complex<Scalar>;
  
public:
  using std::vector<T>::operator=;
  //  using std::vector<T>::vector;
  explicit inline signal(double dt):_dt{dt}{};
  explicit inline signal(size_t sz=0,T defval=static_cast<T>(0.0),double dt=0.0): std::vector<T>(sz,defval),_dt{dt} {};

 
  explicit inline signal(const std::vector<T> lv, const double dt=0.0):std::vector<T>{lv},_dt{dt} {};
  explicit inline signal(std::vector<T>&& rv, const double dt=0.0):std::vector<T>{std::move(rv)},_dt{dt} {};
  
  template<class iter>
  explicit inline signal(iter beg,iter end, const double& dt):std::vector<T>(beg,end),_dt{dt} {};
  
  explicit inline signal(const std::string& filename,double dt):_dt{dt}{read(filename);};  
  inline void read(const std::string& filename) {
    this->clear();
    std::ifstream f(filename.c_str());
    // for(T value;!f.eof();f>>value){
    //   print(value);
    //   this->push_back(value);
    // }
    T value=0.;
    while(!f.eof()) {
      f>>value;
      this->emplace_back(std::move(value));
    }
  };
  inline void write (const std::string& filename) const{
    std::ofstream f(filename.c_str());
    double t=_t0;
    for(const auto& val:*this) {
      f<<t<<'\t'<<val<<'\n';
      // f<<val<<'\n';
      t+=_dt;
    }
  };
  
  inline signal<Complex> fft (void) const {
    double frequency=1./_dt;
    double df=frequency/static_cast<double>(this->size());
    signal<Complex> result(0.,0.,df);
    static   Eigen::FFT<Scalar> FFT;
    FFT.fwd(result,*this);
    result/=(static_cast<double>(this->size())/2.);
    return result;
  };
  inline signal<Scalar>ifft(void) const{
    double time=1./_dt;
    double dtt=time/static_cast<double>(this->size());
    signal<Scalar> result(0.,0.,dtt);
    static Eigen::FFT<Scalar> FFT;
    FFT.inv(result,*this);
    result*=(static_cast<double>(this->size())/2.);
    return result;
  };

  inline signal<Scalar> abs(void) const {
  signal<Scalar> s;s.dt(_dt);
  s.reserve(this->size());
  for(auto el:*this)
    s.emplace_back(std::abs(el) );
  return s;

};

  inline friend signal<Scalar> abs( const signal<T>& v) { return v.abs();};
  inline friend signal<Scalar> abs( signal<T>&& v) {
  signal<Scalar> s{std::move(v)};
  for(auto& el:s)
    el=std::abs(el);
  return s;
};
  ///////////////////////////////////////////////////////////////////////////////////  
  
  inline signal<T> acos(void) const& {
  signal<T> s;s.dt(_dt);
  s.reserve(this->size());   
  for(auto el:*this)
    s.emplace_back(std::acos(el) );
  return s;
};
  inline friend signal<T> acos(const signal<T>& v){ return v.acos();};
  inline signal<T> acos(void) &&{
  for(auto& el:*this)
    el=std::acos(el);
  return *this;
}
  inline friend signal<T> acos( signal<T>&& v) {return std::move(v).acos();};
  //////////////////////////////////////////////////////////////////////////////
  
  inline signal<T> asin(void) const& {
  signal<T> s;s.dt(_dt);
  s.reserve(this->size());   
  for(const auto& el:*this)
    s.emplace_back(std::asin(el) );
  return s;
};
  inline friend signal<T> asin(const signal<T>& v){ return v.asin();};
  inline signal<T> asin(void) &&{
  for(auto& el:*this)
    el=std::asin(el);
  return *this;
}
  inline friend signal<T> asin( signal<T>&& v) {return std::move(v).asin();};
  
  //////////////////////////////////////////////////////////////////////////////////////
  inline signal<T> atan(void) const& {
  signal<T> s;s.dt(_dt);
  s.reserve(this->size());   
  for(auto el:*this)
    s.emplace_back(std::atan(el) );
  return s;
};
  inline friend signal<T> atan(const signal<T>& v){ return v.atan();};
  inline signal<T> atan(void) &&{
  for(auto& el:*this)
    el=std::atan(el);
  return *this;
  };
  inline friend signal<T> atan( signal<T>&& v) { return std::move(v).atan();};
  //////////////////////////////////////////////////////////////////////////////////
  
  inline signal<T> cos(void) const &{
  signal<T> s;s.dt(_dt);
  s.reserve(this->size());   
  for(auto el:*this)
    s.emplace_back(std::cos(el) );
  return s;
};
  inline friend signal<T> cos(const signal<T>& v){ return v.cos();};
  inline signal<T> cos(void) &&{
  for(auto& el:*this)
    el=std::cos(el);
  return *this;
}
  inline friend signal<T> cos( signal<T>&& v) {return std::move(v).cos();};
  //////////////////////////////////////////////////////////////////////////////////////
  
  inline signal<T> cosh(void) const&{
  signal<T> s;s.dt(_dt);
  s.reserve(this->size());   
  for(auto el:*this)
    s.emplace_back(std::cosh(el) );
    return s;
};
  inline friend signal<T> cosh(const signal<T>& v){ return v.cosh();};
  inline signal<T> cosh(void) &&{
  for(auto& el:*this)
    el=std::cosh(el);
  return *this;
}
  inline friend signal<T> cosh( signal<T>&& v) {return std::move(v).cosh();};
  /////////////////////////////////////////////////////////////////////////////////////
  inline signal<T> exp(void) const& {
  signal<T> s;s.dt(_dt);
  s.reserve(this->size());   
  for(auto el:*this)
      s.emplace_back(std::exp(el) );
  return s;
};
  inline friend signal<T> exp(const signal<T>& v){ return v.exp();};
  inline signal<T> exp(void) &&{
  for(auto& el:*this)
    el=std::exp(el);
  return *this;
}
  inline friend signal<T> exp( signal<T>&& v) {  return std::move(v).exp();};
  ///////////////////////////////////////////////////////////////////////////////////////
  inline signal<T> log(void) const& {
  signal<T> s;s.dt(_dt);
  s.reserve(this->size());   
  for(auto el:*this)
    s.emplace_back(std::log(el) );
  return s;
};
  inline friend signal<T> log(const signal<T>& v){ return v.log();};
  inline signal<T> log(void) &&{
  for(auto& el:*this)
    el=std::log(el);
  return *this;
}
  inline friend signal<T> log( signal<T>&& v) { return std::move(v).log();};
  ///////////////////////////////////////////////////////////////////////////  
  inline signal<T> log10(void) const&{
  signal<T> s;s.dt(_dt);
  s.reserve(this->size());   
  for(auto el:*this)
    s.emplace_back(std::log10(el) );
  return s;
};
  inline friend signal<T> log10(const signal<T>& v){ return v.log10();};
  inline signal<T> log10(void) &&{
  for(auto& el:*this)
    el=std::log10(el);
  return *this;
}
  inline friend signal<T> log10( signal<T>&& v) {return std::move(v).log10();};
  //////////////////////////////////////////////////////////////////////////////  
  inline signal<T> sin(void) const&{
  signal<T> s;  s.dt(_dt);
  s.reserve(this->size());   
  for(auto el:*this)
    s.emplace_back(std::sin(el) );
  return s;
};
  inline friend signal<T> sin(const signal<T>& v){ return v.sin();};
  inline signal<T> sin(void) &&{
  for(auto& el:*this)
    el=std::sin(el);
  return *this;
}
  inline friend signal<T> sin( signal<T>&& v) {  return std::move(v).sin();};
  ////////////////////////////////////////////////////////////////////////////////
  inline signal<T> sinh(void) const& {
  signal<T> s;s.dt(_dt);
    s.reserve(this->size());   
    for(auto el:*this)
      s.emplace_back(std::sinh(el) );
    return s;
};
  inline friend signal<T> sinh(const signal<T>& v){ return v.sinh();};
  inline signal<T> sinh(void) &&{
  for(auto& el:*this)
    el=std::sinh(el);
  return *this;
}
  inline friend signal<T> sinh( signal<T>&& v) {  return std::move(v).sinh();};
  //////////////////////////////////////////////////////////////////////////////
  inline signal<T> tan(void) const& {
  signal<T> s;s.dt(_dt);
  s.reserve(this->size());   
  for(auto el:*this)
    s.emplace_back(std::tan(el) );
  return s;
  };
  inline friend signal<T> tan(const signal<T>& v){ return v.tan();};
  inline signal<T> tan(void) &&{
  for(auto& el:*this)
    el=std::tan(el);
  return *this;
}
  inline friend signal<T> tan( signal<T>&& v) {  return std::move(v).tan();};
  /////////////////////////////////////////////////////////////////////////////////
  inline signal<T> tanh(void) const &{
  signal<T> s;s.dt(_dt);
  s.reserve(this->size());   
    for(auto el:*this)
      s.emplace_back(std::tanh(el) );
    return s;
};
  inline friend signal<T> tanh(const signal<T>& v){ return v.tanh();};
  inline signal<T> tanh(void) &&{
  for(auto& el:*this)
    el=std::tanh(el);
  return *this;
}
  inline friend signal<T> tanh( signal<T>&& v) {return std::move(v).tanh();};
  /////////////////////////////////////////////////////////////////////////////////
  inline signal<T> sqrt(void) const& {
  signal<T> s;s.dt(_dt);
  s.reserve(this->size());   
  for(auto el:*this)
    s.emplace_back(std::sqrt(el) );
  return s;
};
  inline friend signal<T> sqrt(const signal<T>& v){ return v.sqrt();};
  inline signal<T> sqrt(void) &&{
  for(auto& el:*this)
    el=std::sqrt(el);
  return *this;
}
  inline friend signal<T> sqrt( signal<T>&& v) {return std::move(v).sqrt();};
  ///////////////////////////////////////////////////////////////////////////////
  inline signal<T> pow(const T& exp)const& {
  signal<T> s;s.dt(_dt);
  s.reserve(this->size());   
  for(auto el:*this)
    s.emplace_back(std::pow(el,exp) );
    return s;
};
  inline friend signal<T> pow(const signal<T>& v,const T& exp) { return v.pow(exp);};
  inline signal<T> pow(const T& exp) &&{
  for(auto& el:*this)
    el=std::pow(el,exp);
  return *this;
  };
  inline friend signal<T> pow( signal<T>&& v,const T& exp) {return std::move(v).pow(exp);};
  //////////////////////////////////////////////////////////////////////////////////////////////////
  inline signal<T>& operator *= (const signal<T>& u) {
  assert(this->size()==u.size());
  for(size_t i=0,sz=this->size();i<sz;++i)
    (*this)[i]*=u[i];
  return *this;
};
  inline signal<T> friend operator * (const signal<T>& v, signal<T>&& u){return u*=v;};
  inline signal<T> friend operator * (signal<T>&& u,const signal<T>& v ){return u*=v;};
  inline signal<T> friend operator * (signal<T>&& u, signal<T>&& v ){return u*=v;};
  inline signal<T> friend operator * (const signal<T>& v, const signal<T>& u){
  signal<T> s(v);
  return s*=u;
}; 
  //////////////////////////////////////////////////////////////////////////////////////////////////
  inline signal<T>& operator *= (const T& u) {
  for(auto& el:*this)
    el*=u;
  return *this;
};
  inline signal<T> friend operator * (signal<T>&& v, const T& u){return v*=u;};
  inline signal<T> friend operator * (const T& u, signal<T>&& v){return v*=u;};
  inline signal<T> friend operator * (const signal<T>& v, const T& u){
  signal<T> s(v);
  return s*=u;
};
  inline signal<T> friend operator * (const T& u,const signal<T>& v){
  signal<T> s(v);
  return s*=u;
};
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  inline signal<T>& operator /= (const T& u) {
  T helper{static_cast<T>(1.)/u};
  for(auto& el:*this)
    el*=helper;
  return *this;
};
  inline signal<T> friend operator / (signal<T>&& v, const T& u){return v/=u;};
  inline signal<T> friend operator / (const T& u, signal<T>&& v){return v/=u;};
  inline signal<T> friend operator / (const signal<T>& v, const T& u){
  signal<T> s(v);
  return s/=u;
};
  inline signal<T> friend operator / (const T& u,const signal<T>& v){
  signal<T> s(v);
  return s/=u;
};
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  inline signal<T>& operator /= (const signal<T>& u) {
  assert(this->size()==u.size());
  for(size_t i=0,sz=this->size();i<sz;++i)
    (*this)[i]/=u[i];
  return *this;
};
  inline signal<T> friend operator / (const signal<T>& v, signal<T>&& u){return u/=v;};
  inline signal<T> friend operator / (signal<T>&& u,const signal<T>& v ){return u/=v;};
  inline signal<T> friend operator / (signal<T>&& u, signal<T>&& v ){return u/=v;};
  inline signal<T> friend operator / (const signal<T>& v, const signal<T>& u){
  signal<T> s(v);
  return s/=u;
};
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  inline signal<T>& operator += (const signal<T>& u) {
  assert(this->size()==u.size());
  for(size_t i=0,sz=this->size();i<sz;++i)
    (*this)[i]+=u[i];
  return *this;
};
  inline signal<T> friend operator + (const signal<T>& v, signal<T>&& u){return u+=v;};
  inline signal<T> friend operator + (signal<T>&& u,const signal<T>& v ){return u+=v;};
  inline signal<T> friend operator + (signal<T>&& u, signal<T>&& v ){return u+=v;};
  inline signal<T> friend operator + (const signal<T>& v, const signal<T>& u){
  signal<T> s(v);
  return s+=u;};
  //////////////////////////////////////////////////////////////////////////////////////////////////////////  
  inline signal<T>& operator -= (const signal<T>& u) {
  assert(this->size()==u.size());
  for(size_t i=0,sz=this->size();i<sz;++i)
    (*this)[i]-=u[i];
  return *this;
};
  inline signal<T> friend operator - (const signal<T>& v, signal<T>&& u){return u-=v;};
  inline signal<T> friend operator - (signal<T>&& u,const signal<T>& v ){return u-=v;};
  inline signal<T> friend operator - (signal<T>&& u, signal<T>&& v ){return u-=v;};
  inline signal<T> friend operator - (const signal<T>& v, const signal<T>& u){
  signal<T> s(v);
  return s-=u;};
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////// 
  inline T sum(void) const {
  auto m= std::accumulate(this->cbegin(),this->cend(),static_cast<T>(0.));
  return m;
};
  inline friend T sum(const signal<T>& v){ return v.sum();};
  template<class Iter>
  inline friend T sum(Iter beg,Iter end){ return std::accumulate(beg,end,static_cast<T>(0.)); };
  
  inline void demean(void) {
  auto m= std::accumulate(this->cbegin(),this->cend(),static_cast<T>(0.));
  m/=static_cast<T>(this->size());
  for(auto& el:*this)
    el-=m;
};
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  inline signal<T> conj(void) const& {
  signal<T> s;s.dt(_dt);
  s.reserve(this->size());   
  for(auto el:*this)
    s.emplace_back(std::conj(el) );
  return s;
  };
  inline friend signal<T> conj(const signal<T>& v){ return v.conj();};
  inline signal<T> conj(void) &&{
  for(auto& el:*this)
    el=std::conj(el);
  return *this;
}
  inline friend signal<T> conj( signal<T>&& v) {  return std::move(v).conj();};
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  inline signal<Scalar> real(void) const {
  signal<Scalar> s;s.dt(_dt);
  s.reserve(this->size());
  for(auto el:*this)
    s.emplace_back(std::real(el) );
  return s;
};
  inline friend signal<Scalar> real( const signal<T>& v) { return v.real();};
  inline friend signal<Scalar> real( signal<T>&& v) {
  signal<Scalar> s{std::move(v)};
  for(auto& el:s)
    el=std::real(el);
  return s;
};
  
  ///////////////////////////////////////////////////////////////////////////////////  
  inline signal<Scalar>imag(void) const {
  signal<Scalar> s;s.dt(_dt);
  s.reserve(this->size());
  for(auto el:*this)
    s.emplace_back(std::imag(el) );
  return s;
};
  inline friend signal<Scalar> imag( const signal<T>& v) { return v.imag();};
  inline friend signal<Scalar> imag( signal<T>&& v) {
  signal<Scalar> s{std::move(v)};
  for(auto& el:s)
    el=std::imag(el);
  return s;
};
  //////////////////////////////////////////////////////////////////////////////////
  
  inline signal<T>bond_pass(double f1,double f2) const {
    assert(f1<f2);
    auto yf=this->fft();
    size_t i1=f1/yf.dt(); //dt in this cas is df
    size_t i2=f2/yf.dt(); //dt in this cas is df
    size_t maxf=yf.size()%2==0?yf.size()/2:(yf.size()/2)+1;
    assert(f2<maxf);
    for(size_t i=0;i<i1;i++)
      yf[i]=0.;
    for(size_t i=i2+1;i<maxf;i++)
      yf[i]=0.;
    return yf.ifft();
  };
  inline signal<T>bond_cut(double f1,double f2 ) const {
    assert(f1<f2);
    auto yf=this->fft();
    size_t i1=f1/yf.dt(); //dt in this cas is df
    size_t i2=f2/yf.dt(); //dt in this cas is df
    size_t maxf=yf.size()%2==0?yf.size()/2:(yf.size()/2)+1;
    assert(f2<maxf);
    for(size_t i=i1;i<=i2;i++)
      yf[i]=0.;
    return yf.ifft();
  };
  inline signal<T>high_pass(double f1) const {
    auto yf=this->fft();
    size_t i1=f1/yf.dt(); //dt in this cas is df
    size_t maxf=yf.size()%2==0?yf.size()/2:(yf.size()/2)+1;
    assert(f1<maxf);
    for(size_t i=0;i<i1;i++)
      yf[i]=0.;
    return yf.ifft();
  };
  inline signal<T>low_pass(double f1) const {
    auto yf=this->fft();
    size_t i1=f1/yf.dt(); //dt in this cas is df
    size_t maxf=yf.size()%2==0?yf.size()/2:(yf.size()/2)+1;
    assert(f1<maxf);
    for(size_t i=i1+1;i<maxf;i++)
      yf[i]=0.;
    return yf.ifft();

  };
  
  template<class Function>
  inline signal<T>filter(const Function& f)const {
    auto yf=this->fft();
    size_t maxf=yf.size()%2==0?yf.size()/2:(yf.size()/2)+1;
    for(size_t i=0;i<maxf;i++)
      yf[i]*=f(i*yf.dt());
    return yf.ifft();
    
  };
  inline signal<Scalar>psd(/*std::vector<double> window*/) {
    //    assert(window.size()%2==1);
   
    return  this->pow(fft().abs(),2.);
  };

  inline signal<T> energy(void) const &{
    signal<T> res(this->size(),0.,this->dt());
    res[0]=std::pow(std::abs((*this)[0]),2.);
    for(size_t i=1;i!=res.size();++i)
      res[i]=res[i-1]+ std::pow(std::abs((*this)[i]),2);
    //    res/=res.back();
    return res;
  };
  inline signal<T> energy(void)&& {
    (*this)[0]= std::pow(std::abs((*this)[0]),2.);
    for(size_t i=1;i!=this->size();++i)
      (*this)[i]=(*this)[i-1]+ std::pow(std::abs((*this)[i]),2.);
    *this/=this->back();
  };
  inline signal<T> extract_with_energy_bound(double f1,double f2) const {
    assert(f1<f2);
    assert(f2<1.);
    double total_energy=0.;
    for(const auto& el:*this)
      total_energy=+std::pow(std::abs(el),2.);
    //    std::cout<<total_energy<<std::endl;
    double temp_energy=0.;
    int i1=-1;
    this->write("secret.data");
    for(size_t i=0;i<this->size();++i) {
      temp_energy+=std::pow((*this)[i],2.);
      if((temp_energy/total_energy)>=f1 ) {
	std::cout<<(total_energy)<<std::endl;
	i1=i;
	break;
      }
    }

    size_t i2=this->size()+1;
    temp_energy=0.;
    for(size_t i=0;i<this->size();++i) {
      temp_energy+=std::pow(std::abs((*this)[i]),2.);
      if((temp_energy/total_energy)>=f2 ) {
	i2=i;
	std::puts("i2");
	break;
      }
    }
    std::cout<<i1<<"\t"<<i2<<std::endl;
    signal<T> res{i2-i1,0.,this->dt()};
    res.t0(this->t0()+i1*this->dt());
    
    std::copy(begin(*this)+i1,end(*this)+i2,
	      std::back_inserter(res) );
    
    return res;
  };
  //inline signal<T>
  
  inline T mean(void) const {
  auto m= std::accumulate(this->cbegin(),this->cend(),static_cast<T>(0.));
  m/=static_cast<T>(this->size());
  return m;
};
  inline friend T mean(const signal<T>& v){ return v.mean();};
  
  inline T variance(void) const {
  assert(this->size()!=0);
  auto m= this->mean();
  size_t n = this->size()-1;
  T s=0.;
  for(auto el:*this)
    s+=(el-m)*(el-m);
  return s/static_cast<T>(n);
};
  inline friend T variance(const signal<T>& v){ return v.variance();};
  
  inline T stddev(void) const {
  return std::sqrt(this->variance());
};
  inline friend T stddev(const signal<T>& v){ return v.stddev();};

  
  inline signal<T> autocor(method meth=method::unbiased) const {
  signal<T> res(this->dt());
  res.reserve(this->size());
  res.t0(0.);
  for(size_t lag=0;lag<this->size();lag++) {
    double s=0.;
    for(size_t i=0;i<this->size()-lag;i++){
      s+=  (*this)[i+lag]* (*this)[i];
    }
    if(meth==method::unbiased)
      s/=static_cast<T>(this->size()-lag);
    else if(meth==method::biased)
      s/=static_cast<T>(this->size());
    res.emplace_back(s);
  }
  assert(this->size()==res.size());
  return res;
};
  inline friend signal<T> autocor(const signal<T>& v,method meth=method::unbiased){ return v.autocor(meth);};

  inline signal<T> xcor(const signal<T>& v,method meth=method::unbiased) const { 
  assert(this->dt()==v.dt());
  assert(this->size()==v.size());
  signal<T> res(this->dt());
  res.reserve(2*this->size()-1);
  res.t0(-this->dt()*static_cast<T>(this->size()-1));
  for(size_t lag=this->size()-1;lag>=1;lag--) {
    double s=0.;
    for(size_t i=0;i<this->size()-lag;i++){
      s+=  (*this)[i]* v[i+lag];
    }
    if(meth==method::unbiased)
      s/=static_cast<T>(this->size()-lag);
    else if(meth==method::biased)
      s/=static_cast<T>(this->size());
    
    res.emplace_back(s);
  }
  for(size_t lag=0;lag<this->size();lag++) {
    double s=0.;
    for(size_t i=0;i<this->size()-lag;i++) {
      s+=  (*this)[i+lag]* v[i];
    }
    if(meth==method::unbiased)
      s/=static_cast<T>(this->size()-lag);
    else if(meth==method::biased)
      s/=static_cast<T>(this->size());
    
    res.emplace_back(s);
  }
  assert( 2*this->size()-1 == res.size()) ;
  return res;
  };
  inline friend signal<T> xcor(const signal<T>& v, const signal<T>& s, method meth=method::unbiased) {  return v.xcor(s,meth); };
  
  inline signal<T> autocov(method meth=method::unbiased) const {
  signal<T> res(this->dt());
  res.reserve(this->size());
  res.t0(0.);
  T m= this->mean();
  for(size_t lag=0;lag<this->size();lag++) {
    double s=0.;
    for(size_t i=0;i<this->size()-lag;i++){
      s+= ( (*this)[i+lag] -m) *( (*this)[i]-m);
    }
    if(meth==method::unbiased)
      s/=static_cast<T>(this->size()-lag);
    else if(meth==method::biased)
      s/=static_cast<T>(this->size());
    res.emplace_back(s);
  }
  assert(this->size()==res.size());
  return res;
  };
  inline friend signal<T> autocov(const signal<T>& v,method meth=method::unbiased) { return v.autocov(meth);};

  inline signal<T> xcov(const signal<T>& v,method meth=method::unbiased) const { 
  assert(this->dt()==v.dt());
  assert(this->size()==v.size());
  signal<T> res(this->dt());
  res.reserve(2*this->size()-1);
  res.t0(-this->dt()*static_cast<T>(this->size()-1));
  T m1=this->mean();
  T m2=v.mean();
  for(size_t lag=this->size()-1;lag>=1;lag--) {
    double s=0.;
    for(size_t i=0;i<this->size()-lag;i++){
      s+=  ((*this)[i]-m1)* (v[i+lag]-m2);
    }
    if(meth==method::unbiased)
      s/=static_cast<T>(this->size()-lag);
    else if(meth==method::biased)
      s/=static_cast<T>(this->size());
    res.emplace_back(s);
  }
  for(size_t lag=0;lag<this->size();lag++) {
    double s=0.;
    for(size_t i=0;i<this->size()-lag;i++){ 
      s+=  ((*this)[i+lag] -m1)* (v[i]-m2);
    }
    if(meth==method::unbiased)
      s/=static_cast<T>(this->size()-lag);
    else if(meth==method::biased)
      s/=static_cast<T>(this->size());
    res.emplace_back(s);
  }
  assert( 2*this->size()-1 == res.size()) ;
  return res;
};
  inline friend signal<T> xcov(const signal<T>& v, const signal<T>& s,method meth=method::unbiased) { return v.xcov(s,meth); };
///////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  inline double t0(void) const {return _t0;};
  inline void t0(double t)  { _t0=t;};
  inline double dt(void) const {return _dt;};
  inline void dt(double dt)  { _dt=dt;};
  inline double fs(void) const {return 1./_dt;};
  inline void fs(double f){_dt=1./f;};
};

template<class T>   //the by defeault dt is (beg-end)/(n-1);
inline  signal<T> linspace (T beg,T end, size_t n) {
  signal<T> v;
  v.reserve(n);
  T diff_n_1=(end-beg)/static_cast<T>(n-1);
  for(size_t i=0;i<n;++i)
    v.emplace_back(beg+(diff_n_1)*static_cast<T>(i));
  v.dt(std::abs(diff_n_1));
  return v;
};

template<class T=double>///note that this function does not gives dt 
inline auto uniform_rand(size_t number,std::array<T,2> v = std::array<double,2>{0.,1.} )  {  
  std::mt19937 eng(time(0));
  std::uniform_real_distribution<T> urd(v[0], v[1]); 
  signal<T> s;
  s.reserve(number);  
  for(size_t i=0;i<number;++i) {
    s.emplace_back(urd(eng));
  }
  
  return s;
};




