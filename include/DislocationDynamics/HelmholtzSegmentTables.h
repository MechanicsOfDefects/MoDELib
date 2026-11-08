/* MoDELib: Helmholtz straight-segment stress, concentration and flux endpoint formulas. */
#ifndef model_HelmholtzSegmentTables_H_
#define model_HelmholtzSegmentTables_H_
#include <DislocationDynamics/HelmholtzSegmentGeometry.h>
#include <DislocationDynamics/StraightDislocationSegmentMoments.h>
#include <array>
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <memory>
#include <mutex>
#include <list>
#include <tuple>
#include <map>
#include <Eigen/Eigenvalues>
#include <numbers>
namespace model { namespace helmholtzSegment {
using V=std::array<double,2>;
constexpr double pi=3.14159265358979323846;
using Rule=model::SegmentGaussRule;
struct Kernel {
 const double h,C,D;
 const bool full;
 std::array<double,32> c{};
 explicit Kernel(double ratio,bool unscreened=false):h(ratio),C(unscreened?2:2*std::pow(h,4)/std::pow(h*h-1,2)),D(unscreened?0:h*h*h/(h*h-1)),full(unscreened){
 if(!std::isfinite(h)||h<4) {
      throw std::invalid_argument("Helmholtz segment screen requires h>=4*ell");
    }
    if(full){c[0]=2;for(int n=1;n<32;++n) {
        c[n]=-c[n-1]/(n+1);
      }
      c[1]=0;return;}
 long double fac=1,H=h,cc=2*H*H*H*H/((H*H-1)*(H*H-1)),dd=H*H*H/(H*H-1);
 for(int n=0;n<32;++n){if(n) {
        fac*=n;
      }
      c[n]=double(cc*(std::pow(-1/H,n+1)-std::pow(-1.L,n+1))/(fac*(n+1))-dd*std::pow(-1/H,n)/fac);}c[1]=0;
 }
 V operator()(double r,double z)const{
  double R=std::hypot(r,z),nx=r/R,A=0,B=0;
  if(R<.75){for(int n=31;n>=3;--n){A=A*R+n*(n-2)*(n-4)*c[n];B=B*R+n*(n-2)*c[n];}}
  else {double f1=0,f2=0,f3=0;double inv=1/R;
   for(int j=0;j<3;++j){double k=j==1?1.:(full?0.:1/h),coef=j==0?C:(j==1?-C:-D),e=coef*std::exp(-k*R);
    if(j<2){f1-=e*inv*(k+inv);f2+=e*inv*(k*k+2*k*inv+2*inv*inv);f3-=e*inv*(k*k*k+3*k*k*inv+6*k*inv*inv+6*inv*inv*inv);}
    else{f1-=k*e;f2+=k*k*e;f3-=k*k*k*e;}}
   if(full) {
        f1+=1;
      }
      B=f2*inv-f1*inv*inv;A=f3-3*B;
  }
  return {A*nx*nx*nx+3*B*nx,B*nx};
 }
 V integral(double r,double ua,double ub,const Rule& q)const{
  V sum{};double half=(ua-ub)/2,mid=(ua+ub)/2;
  for(size_t i=0;i<q.x.size();++i){double t=mid+half*q.x[i];auto v=(*this)(r,r*std::sinh(t));double w=r*std::cosh(t)*half*q.w[i];for(int k=0;k<2;++k) {
        sum[k]+=w*v[k];
      }
    }return sum;
 }
};
// Piecewise tensor Chebyshev polynomial; p=4 is a bicubic interpolant.
// Bicubic coefficients are converted to monomials for nested Horner evaluation.
struct Table {
 bool tail; int p,tiles;double lo=std::log(1e-4),hi,uMax,dx,du;
 std::vector<V>coeff;
 size_t index(int i,int j,int a,int b)const{return ((size_t(i)*tiles+j)*p+a)*p+b;}
 static int checkedOrder(int p){if(p<4||p>16) {
      throw std::invalid_argument("Segment Chebyshev order must be 4..16");
    }
    return p;}
 static int checkedResolution(int n){if(n<4||n>256) {
      throw std::invalid_argument("Segment table resolution must be 4..256");
    }
    return n;}
 Table(int order,int n,const Kernel&kernel,const Rule&q,bool isTail,double maxRho,double maxU):tail(isTail),p(checkedOrder(order)),tiles(checkedResolution(n)),hi(std::log(maxRho)),uMax(maxU),dx((hi-lo)/n),du(maxU/n),coeff(size_t(n)*n*p*p){
  std::vector<double> nodes(p),cosine(p*p);std::vector<V>f(p*p),tmp(p*p);
  for(int a=0;a<p;++a){nodes[a]=std::cos(pi*(a+.5)/p);for(int k=0;k<p;++k) {
        cosine[k*p+a]=std::cos(k*pi*(a+.5)/p);
      }
    }
  for(int i=0;i<n;++i) {
      for(int j=0;j<n;++j){
   for(int a=0;a<p;++a) {
          for(int b=0;b<p;++b){double r=std::exp(lo+(i+.5)*dx+.5*dx*nodes[a]);double u=(j+.5)*du+.5*du*nodes[b];double top=(kernel.full?std::max(u+24.,std::asinh(60/r)):std::asinh(std::sinh(u)+60*kernel.h/r));auto v=tail?kernel.integral(r,top,u,q):kernel.integral(r,u,0,q);for(auto&x:v) {
              x/=r;
            }
            f[a*p+b]=v;}
        }
        for(int a=0;a<p;++a) {
          for(int b=0;b<p;++b){V v{};for(int s=0;s<p;++s) {
              for(int k=0;k<2;++k) {
                v[k]+=f[s*p+b][k]*cosine[a*p+s]*(a?2.:1.)/p;
              }
            }
            tmp[a*p+b]=v;}
        }
        for(int a=0;a<p;++a) {
          for(int b=0;b<p;++b){V v{};for(int s=0;s<p;++s) {
              for(int k=0;k<2;++k) {
                v[k]+=tmp[a*p+s][k]*cosine[b*p+s]*(b?2.:1.)/p;
              }
            }
            coeff[index(i,j,a,b)]=v;}
        }
        if(p==4){ // T0=1, T1=x, T2=2x^2-1, T3=4x^3-3x.
    const double M[4][4]={{1,0,-1,0},{0,1,0,-3},{0,0,2,0},{0,0,0,4}};
    std::array<V,16> out{};
    for(int a=0;a<4;++a) {
            for(int b=0;b<4;++b) {
              for(int s=0;s<4;++s) {
                for(int t=0;t<4;++t) {
                  for(int k=0;k<2;++k) {
                    out[a*4+b][k]+=M[a][s]*M[b][t]*coeff[index(i,j,s,t)][k];
                  }
                }
              }
            }
          }
          for(int a=0;a<4;++a) {
            for(int b=0;b<4;++b) {
              coeff[index(i,j,a,b)]=out[a*4+b];
            }
          }
        }
  }
    }
  }
 // Difference of a patch polynomial, evaluated as dy times its divided
 // difference. This never subtracts nearly equal endpoint values. Across
 // patches sum each polynomial increment, avoiding interpolation jumps.
 V increment(double r,double u0,double delta)const{
  if(!(r>=1e-4&&r<=std::exp(hi)&&u0>=0&&delta>=0&&u0+delta<=uMax+1e-13)) {
      throw std::out_of_range("increment domain");
    }
    V result{};double l=(std::log(r)-lo)/dx;int i=std::min(tiles-1,int(l));double x=2*(l-i)-1;
  for(int j=std::min(tiles-1,int(u0/du));delta>0&&j<tiles;++j){double span=std::min(delta,(j+1)*du-u0);
   if(span<=0) {
        continue;
      }
      double dy=2*span/du,yb=2*(u0/du-j)-1,ya=yb+dy;
   for(int k=0;k<2;++k){double c[16];
    for(int b=0;b<p;++b){if(p==4){double v=0;for(int a=3;a>=0;--a) {
              v=v*x+coeff[index(i,j,a,b)][k];
            }
            c[b]=v;}
     else{double b1=0,b2=0;for(int a=p-1;a>=1;--a){double v=2*x*b1-b2+coeff[index(i,j,a,b)][k];b2=b1;b1=v;}c[b]=x*b1-b2+coeff[index(i,j,0,b)][k];}}
    double v=c[1];
    if(p==4){double dd=1,power=1;for(int n=2;n<p;++n){power*=yb;dd=ya*dd+power;v+=c[n]*dd;}}
    else{double d0=0,d1=1,t0=1,t1=yb;for(int n=2;n<p;++n){double dn=2*ya*d1+2*t1-d0;v+=c[n]*dn;d0=d1;d1=dn;double tn=2*yb*t1-t0;t0=t1;t1=tn;}}
    result[k]+=r*dy*v;
   }
   delta-=span;u0+=span;
  }
  return result;
 }
 V eval(double r,double u)const{
  if(!(r>=1e-4&&r<=std::exp(hi)&&std::abs(u)<=uMax)) {
      throw std::out_of_range("endpoint table domain");
    }
    double l=(std::log(r)-lo)/dx,v=std::abs(u)/du;int i=std::min(tiles-1,int(l)),j=std::min(tiles-1,int(v));double x=2*(l-i)-1,y=2*(v-j)-1;
  V out{};
  for(int k=0;k<2;++k){
   if(p==4){double rows[4];for(int a=0;a<4;++a){double s=0;for(int b=3;b>=0;--b) {
            s=s*y+coeff[index(i,j,a,b)][k];
          }
          rows[a]=s;}double s=0;for(int a=3;a>=0;--a) {
          s=s*x+rows[a];
        }
        out[k]=s;}
   else{double rows[16];for(int a=0;a<p;++a){double b1=0,b2=0;for(int b=p-1;b>=1;--b){double t=2*y*b1-b2+coeff[index(i,j,a,b)][k];b2=b1;b1=t;}rows[a]=y*b1-b2+coeff[index(i,j,a,0)][k];}
    double b1=0,b2=0;for(int a=p-1;a>=1;--a){double t=2*x*b1-b2+rows[a];b2=b1;b1=t;}out[k]=x*b1-b2+rows[0];}
   out[k]*=r*(tail?1:(u<0?-1:(u>0?1:0)));
  }return out;
 }
};
// Pair arithmetic stores f(b) and f(a)-f(b). The difference is propagated
// analytically, rather than recovered by subtracting close function values.
struct Pair {double v,d;Pair(double x=0,double diff=0):v(x),d(diff){} };
inline Pair operator+(Pair a,Pair b){return {a.v+b.v,a.d+b.d};}
inline Pair operator-(Pair a,Pair b){return {a.v-b.v,a.d-b.d};}
inline Pair operator*(Pair a,Pair b){return {a.v*b.v,a.d*b.v+a.v*b.d+a.d*b.d};}
inline Pair operator/(Pair a,Pair b){return {a.v/b.v,(a.d*b.v-a.v*b.d)/(b.v*(b.v+b.d))};}
inline Pair root(Pair a){double r=std::sqrt(a.v);return {r,a.d/(std::sqrt(a.v+a.d)+r)};}
inline Pair exponential(Pair a){double e=std::exp(a.v);return {e,e*std::expm1(a.d)};}
// [Hxx,Hyy,Hxz,Hzz] differences in the segment-local coordinate frame.
inline std::array<double,4> stableHessian(const Kernel&k,double rho,double za,double zb,double dz){
 Pair z(zb,dz),r=root(Pair(rho*rho)+z*z),inv=Pair(1)/r,a,e;
 if(std::max(std::hypot(rho,za),std::hypot(rho,zb))<.75){
  for(int n=31;n>=2;--n) {
      a=a*r+Pair(n*k.c[n]);
    }
    for(int n=31;n>=3;--n) {
      e=e*r+Pair(n*(n-2)*k.c[n]);
    }
    e=e*inv;
 }else{
  Pair i2=inv*inv,i3=i2*inv;
  if(k.full){a=inv;e=Pair(-1)*i3;}
  for(int j=0;j<3;++j){double q=j==1?1.:(k.full?0.:1/k.h),co=j==0?k.C:(j==1?-k.C:-k.D);Pair ex=Pair(co)*exponential(Pair(-q)*r);
   if(j<2){a=a-ex*i2*(Pair(q)+inv);e=e+ex*i3*(Pair(q*q)+Pair(3*q)*inv+Pair(3)*i2);}
   else{a=a-ex*Pair(q)*inv;e=e+ex*i2*(Pair(q*q)+Pair(q)*inv);}
  }
 }
 return {(a+e*Pair(rho*rho)).d,a.d,(e*Pair(rho)*z).d,(a+e*z*z).d};
}
inline std::array<double,4> endpointHessian(const Kernel&k,double rho,double z){
 double R=std::hypot(rho,z),inv=1/R,a=0,d=0;
 if(R<.75){for(int n=31;n>=2;--n){a=a*R+n*k.c[n];d=d*R+n*(n-2)*k.c[n];}}
 else{
  double i2=inv*inv;
  if(k.full){a=inv;d=-inv;}
  for(int j=0;j<3;++j){double q=j==1?1.:(k.full?0.:1/k.h),co=j==0?k.C:(j==1?-k.C:-k.D),e=co*std::exp(-q*R);
   if(j<2){a-=e*i2*(q+inv);d+=e*inv*(q*q+3*q*inv+3*i2);}
   else{a-=e*q*inv;d+=e*(q*q+q*inv);}
  }
 }
 double e=d*inv*inv;return {a+e*rho*rho,a,e*rho*z,a+e*z*z};
}
// Contiguous two-channel contractions. Fixed orders permit compiler unrolling;
// no fast-math or reduced precision is used.
template<int P>inline std::array<V,P> collapse(const Table&t,int i,int j,double x){
 std::array<double,P> tx;tx[0]=1;tx[1]=x;for(int a=2;a<P;++a) {
    tx[a]=2*x*tx[a-1]-tx[a-2];
  }
  std::array<V,P> rows{};
 const V*base=&t.coeff[t.index(i,j,0,0)];
 for(int a=0;a<P;++a) {
    for(int b=0;b<P;++b){rows[b][0]+=tx[a]*base[a*P+b][0];rows[b][1]+=tx[a]*base[a*P+b][1];}
  }
  return rows;
}
template<int P>inline V fastIncrement(const Table&t,double rho,double u,double delta){
 double l=(std::log(rho)-t.lo)/t.dx;int i=std::min(t.tiles-1,int(l));double x=2*(l-i)-1;V out{};
 for(int j=std::min(t.tiles-1,int(u/t.du));delta>0&&j<t.tiles;++j){double span=std::min(delta,(j+1)*t.du-u);
  if(span<=0) {
      continue;
    }
    double dy=2*span/t.du,yb=2*(u/t.du-j)-1,ya=yb+dy;auto rows=collapse<P>(t,i,j,x);
  V v=rows[1];double d0=0,d1=1,t0=1,t1=yb;
  for(int n=2;n<P;++n){double dn=2*ya*d1+2*t1-d0;v[0]+=rows[n][0]*dn;v[1]+=rows[n][1]*dn;d0=d1;d1=dn;double tn=2*yb*t1-t0;t0=t1;t1=tn;}
  out[0]+=rho*dy*v[0];out[1]+=rho*dy*v[1];delta-=span;u+=span;
 }return out;
}
template<int P>inline V fastValue(const Table&t,double rho,double u){
 double l=(std::log(rho)-t.lo)/t.dx,v=std::abs(u)/t.du;int i=std::min(t.tiles-1,int(l)),j=std::min(t.tiles-1,int(v));double x=2*(l-i)-1,y=2*(v-j)-1;auto rows=collapse<P>(t,i,j,x);
 V out=rows[0];double y0=1,y1=y;out[0]+=y*rows[1][0];out[1]+=y*rows[1][1];
 for(int n=2;n<P;++n){double yn=2*y*y1-y0;out[0]+=yn*rows[n][0];out[1]+=yn*rows[n][1];y0=y1;y1=yn;}
 double sign=t.tail?1:(u<0?-1:(u>0?1:0));out[0]*=rho*sign;out[1]*=rho*sign;return out;
}
inline V incrementFast(const Table&t,double r,double u,double delta){
 if(!(r>=1e-4&&r<=std::exp(t.hi)&&u>=0&&delta>=0&&u+delta<=t.uMax+1e-13)) {
    throw std::out_of_range("fast increment domain");
  }
  switch(t.p){case 16:return fastIncrement<16>(t,r,u,delta);case 12:return fastIncrement<12>(t,r,u,delta);case 10:return fastIncrement<10>(t,r,u,delta);case 8:return fastIncrement<8>(t,r,u,delta);default:return t.increment(r,u,delta);}
}
inline V valueFast(const Table&t,double r,double u){
 if(!(r>=1e-4&&r<=std::exp(t.hi)&&std::abs(u)<=t.uMax)) {
    throw std::out_of_range("fast value domain");
  }
  switch(t.p){case 16:return fastValue<16>(t,r,u);case 12:return fastValue<12>(t,r,u);case 10:return fastValue<10>(t,r,u);case 8:return fastValue<8>(t,r,u);default:return t.eval(r,u);}
}

using Vec=Eigen::Vector3d;using Mat=Eigen::Matrix3d;

// Immutable CPU tables survive evaluator/geometry rebuilds. Only lookup and
// construction take the lock; all field queries read shared tables without locks.
struct SharedTables {
 Table lower,tail;
 SharedTables(int order,int resolution,const Kernel& kernel,const Rule& rule,double radius,double limit)
 :lower(order,resolution,kernel,rule,false,radius,limit),
  tail(order,resolution,kernel,rule,true,radius,limit) {}
};
inline std::shared_ptr<const SharedTables> cachedTables(const Kernel& kernel,const Rule& rule,
    int order,int resolution,double maxDistance,double ell,bool& reused)
{
 if(!std::isfinite(ell)||ell<=0||!std::isfinite(maxDistance)||maxDistance<=0) {
  throw std::invalid_argument("Helmholtz table core and extent must be positive and finite");
 }
 Table::checkedOrder(order);Table::checkedResolution(resolution);
 const double radius=std::max(32.,maxDistance/ell);
 const double limit=std::max(12.,std::asinh(maxDistance/(ell*1e-4))+1);
 // Full-kernel tables do not depend on the unused Ewald width, or on mu/nu.
 using Key=std::tuple<bool,double,int,int,double,double>;
 const Key key{kernel.full,kernel.full?0.:kernel.h,order,resolution,radius,limit};
 using Entry=std::pair<Key,std::shared_ptr<const SharedTables>>;
 static std::mutex mutex;
 static std::list<Entry> cache;
 std::lock_guard<std::mutex> lock(mutex);
 for(auto it=cache.begin();it!=cache.end();++it) {
  if(it->first==key) {
   reused=true;
   auto result=it->second;
   cache.splice(cache.begin(),cache,it);
   return result;
  }
 }
 auto result=std::make_shared<const SharedTables>(order,resolution,kernel,rule,radius,limit);
 cache.emplace_front(key,result);
 // Bound retained memory when many widths are explored. Live evaluators retain
 // ownership even if their entry is evicted. The usual two modes stay cached.
 if(cache.size()>4) {
  cache.pop_back();
 }
 reused=false;
 return result;
}

class Evaluator {
 double ell_,mu_,nu_;Kernel kernel;Rule q80{80};
 bool tablesReused_=false;
 std::shared_ptr<const SharedTables> tables_;
 const Table& lower;
 const Table& tail;
 V transverse(double r,double za,double zb,double dz)const{
  if(r==0) {
      return {0,0};
    }
    double ua=std::asinh(za/r),ub=std::asinh(zb/r);
  if(r<1e-4 || r>std::exp(lower.hi) || std::max(std::abs(ua),std::abs(ub))>lower.uMax){
   // Rare out-of-table near-axis queries: integrate bounded transformed panels.
   V out{};int panels=std::max(1,int(std::ceil(std::abs(ua-ub)/2)));
   for(int i=0;i<panels;++i){auto v=kernel.integral(r,ub+(ua-ub)*(i+1)/panels,ub+(ua-ub)*i/panels,q80);for(int c=0;c<2;++c) {
          out[c]+=v[c];
        }
      }return out;
  }
  if(ua*ub>0){double delta=std::asinh(dz*(za+zb)/(za*std::hypot(r,zb)+zb*std::hypot(r,za)));double u0=std::min(std::abs(ua),std::abs(ub));bool useTail=u0>1;
   auto v=incrementFast(useTail?tail:lower,r,u0,std::abs(delta));double sign=(delta<0?-1.:1.)*(useTail?-1.:1.);for(auto&x:v) {
        x*=sign;
      }
      return v;}
  auto a=valueFast(lower,r,ua),b=valueFast(lower,r,ub);return {a[0]-b[0],a[1]-b[1]};
 }
public:
 Evaluator(double ell,double h,double mu,double nu,int order,int resolution,double maxDistance,bool full=false)
 :ell_(ell),mu_(mu),nu_(nu),kernel(h/ell,full),
 tables_(cachedTables(kernel,q80,order,resolution,maxDistance,ell,tablesReused_)),
 lower(tables_->lower),tail(tables_->tail) {}
 bool reusedTables() const { return tablesReused_; }
 Mat stress(const Vec&x,const Segment&g)const{
  Vec chord=(g.p1-g.p0)/ell_;double dz=chord.norm();if(dz==0) {
      return Mat::Zero();
    }
    Vec t=chord/dz;
  Vec r=(x-g.p0)/ell_;double za=r.dot(t),zb=za-dz;Vec perp=r-za*t;perp-=perp.dot(t)*t;double rho=perp.norm();
  if(rho<=64*std::numeric_limits<double>::epsilon()*std::max(1.,r.norm())) {
      rho=0;
    }
    Vec ex=rho>0?Vec(perp/rho):Vec(t.unitOrthogonal()),ey=t.cross(ex);
  auto tr=transverse(rho,za,zb,dz);std::array<double,4> H;
  if(std::abs(dz)<.01*std::min(std::hypot(rho,za),std::hypot(rho,zb))) {
      H=stableHessian(kernel,rho,za,zb,dz);
    } else{auto hessian=[&](double z){if(rho==0&&z==0) {
          return std::array<double,4>{2*kernel.c[2],2*kernel.c[2],0,2*kernel.c[2]};
        }
        return endpointHessian(kernel,rho,z);};auto a=hessian(za),b=hessian(zb);for(int i=0;i<4;++i) {
        H[i]=a[i]-b[i];
      }
    }
  double bx=g.b.dot(ex),by=g.b.dot(ey),bz=g.b.dot(t),gx=tr[0]+tr[1]+H[2],gz=H[0]+H[1]+H[3];
  double f=mu_/(8*pi*ell_),q=mu_/(4*pi*(1-nu_)*ell_);
  double xx=q*by*(tr[0]-gx),yy=q*by*(tr[1]-gx),zz=2*f*gx*by+q*by*(H[2]-gx);
  double xy=-q*bx*tr[1],xz=-f*gz*by+q*by*H[0],yz=f*(gz*bx-gx*bz)-q*bx*H[1];Mat out;
  for(int i=0;i<3;++i) {
      for(int j=i;j<3;++j) {
        out(i,j)=out(j,i)=xx*ex[i]*ex[j]+yy*ey[i]*ey[j]+zz*t[i]*t[j]+xy*(ex[i]*ey[j]+ey[i]*ex[j])+xz*(ex[i]*t[j]+t[i]*ex[j])+yz*(ey[i]*t[j]+t[i]*ey[j]);
      }
    }
    return out;
 }
};
// Each open segment piece is assigned once by its midpoint. Face-coincident
// pieces use floor/clamp ownership, matching the leaf target convention.
inline std::vector<std::vector<Segment>> split(const std::vector<Segment>&segments,const Vec&lo,const Vec&hi,int depth){
 int n=1<<depth;Vec w=(hi-lo)/n;std::vector<std::vector<Segment>> bins(n*n*n);
 for(const auto&s:segments){
  if(!s.p0.allFinite()||!s.p1.allFinite()||!s.b.allFinite()) {
      throw std::invalid_argument("Nonfinite stress segment");
    }
    double tol=64*std::numeric_limits<double>::epsilon()*std::max({1.,lo.cwiseAbs().maxCoeff(),hi.cwiseAbs().maxCoeff()});
  if((s.p0.array()<lo.array()-tol).any()||(s.p1.array()<lo.array()-tol).any()||(s.p0.array()>hi.array()+tol).any()||(s.p1.array()>hi.array()+tol).any()) {
      throw std::out_of_range("Stress segment outside prepared tree");
    }
    Vec d=s.p1-s.p0;if(d.squaredNorm()==0) {
      continue;
    }
    std::vector<double> cuts{0,1};
  for(int c=0;c<3;++c){if(d[c]==0) {
        continue;
      }
      for(int j=1;j<n;++j){double a=(lo[c]+j*w[c]-s.p0[c])/d[c];if(a>0&&a<1) {
          cuts.push_back(a);
        }
      }}
  std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
  for(size_t j=1;j<cuts.size();++j){double a=cuts[j-1],b=cuts[j];Vec mid=s.p0+(a+(b-a)/2)*d;int cell[3];for(int c=0;c<3;++c) {
        cell[c]=std::clamp(int(std::floor((mid[c]-lo[c])/w[c])),0,n-1);
      }
      Vec p0=a==0?s.p0:Vec(s.p0+a*d),p1=b==1?s.p1:Vec(s.p0+b*d);if((p1-p0).squaredNorm()>0) {
        bins[cell[0]+n*(cell[1]+n*cell[2])].push_back({p0,p1,s.b});
      }
    }
 }return bins;
}

    // Primitive differences of the full or screened Helmholtz diffusion Green function.
    // Dimensionless radial distance uses the determinant-one diffusion metric and core length.
    // Each patch stores a polynomial in longitudinal distance, so both its primitive and
    // its first/second-moment primitives are elementary. No line quadrature is used at query time.
    class ConcentrationEvaluator
    {
        struct Table
        {
            int order,resolution;
            double rhoMin=1e-12,rhoMax,logMin,dx,dt;
            std::vector<double> coefficients;
            Table(double ratio,bool screened,int p,int n,double maximum,bool flux)
                :order(p),resolution(n),rhoMin(flux?1e-14:1e-12),rhoMax(std::max(1.0,maximum)),logMin(std::log(rhoMin)),
                 dx((std::log(rhoMax)-logMin)/n),dt(std::asinh(rhoMax/rhoMin)/n),
                 coefficients(std::size_t(n)*n*p*p)
            {
                std::vector<double> nodes(p),cosine(p*p),values(p*p),temporary(p*p);
                for(int a=0;a<p;++a)
                {
                    nodes[a]=std::cos(std::numbers::pi*(a+.5)/p);
                    for(int k=0;k<p;++k)
                    {
                        cosine[k*p+a]=std::cos(k*std::numbers::pi*(a+.5)/p);
                    }
                }
                for(int i=0;i<n;++i)
                {
                    for(int j=0;j<n;++j)
                    {
                        const double zlo=std::sinh(j*dt),zhi=std::sinh((j+1)*dt);
                        for(int a=0;a<p;++a)
                        {
                            const double rho=std::exp(logMin+(i+.5)*dx+.5*dx*nodes[a]);
                            for(int b=0;b<p;++b)
                            {
                                const double z=rho*((zlo+zhi)*.5+(zhi-zlo)*.5*nodes[b]);
                                const double R=std::hypot(rho,z);
                                if(flux)
                                {
                                    // -g'(R)/R; evaluate the cancellation-prone numerator by series.
                                    const auto numerator=[](double z)
                                    {
                                        if(z<.01)
                                        {
                                            return z*z*(.5+z*(-1./3+z*(1./8+z*(-1./30+z*(1./144+z*(-1./840+z/5760))))));
                                        }
                                        return -std::expm1(-z)-z*std::exp(-z);
                                    };
                                    values[a*p+b]=(screened?(numerator(R)-numerator(R/ratio))/(1-1/(ratio*ratio)):numerator(R))/(R*R*R);
                                }
                                else
                                {
                                    values[a*p+b]=screened?
                                        std::exp(-R/ratio)*(-std::expm1(-(1-1/ratio)*R))/((1-1/(ratio*ratio))*R):
                                        -std::expm1(-R)/R;
                                }
                            }
                        }
                        for(int a=0;a<p;++a)
                        {
                            for(int b=0;b<p;++b)
                            {
                                double value=0;
                                for(int k=0;k<p;++k)
                                {
                                    value+=values[k*p+b]*cosine[a*p+k]*(a?2.0:1.0)/p;
                                }
                                temporary[a*p+b]=value;
                            }
                        }
                        for(int a=0;a<p;++a)
                        {
                            for(int b=0;b<p;++b)
                            {
                                double value=0;
                                for(int k=0;k<p;++k)
                                {
                                    value+=temporary[a*p+k]*cosine[b*p+k]*(b?2.0:1.0)/p;
                                }
                                coefficients[index(i,j,a,b)]=value;
                            }
                        }
                    }
                }
            }
            std::size_t index(int i,int j,int a,int b) const
            {
                return ((std::size_t(i)*resolution+j)*order+a)*order+b;
            }
            // Exact polynomial primitive differences on [c-d,c+d]. Taylor shifting
            // avoids subtracting close endpoint values, even for extremely short pieces.
            Eigen::Vector3d primitive(int i,int j,double radial,double c,double d) const
            {
                std::array<double,16> cheb{},power{},t0{},t1{};
                for(int b=0;b<order;++b)
                {
                    double b1=0,b2=0;
                    for(int a=order-1;a>=1;--a)
                    {
                        const double v=2*radial*b1-b2+coefficients[index(i,j,a,b)];b2=b1;b1=v;
                    }
                    cheb[b]=radial*b1-b2+coefficients[index(i,j,0,b)];
                }
                t0[0]=1;t1[1]=1;
                power[0]=cheb[0];power[1]=cheb[1];
                for(int k=2;k<order;++k)
                {
                    std::array<double,16> next{};
                    for(int a=0;a<=k;++a)
                    {
                        next[a]=(a?2*t1[a-1]:0)-t0[a];
                        power[a]+=cheb[k]*next[a];
                    }
                    t0=t1;t1=next;
                }
                std::array<double,16> shifted{};
                // Horner shift p(c+v), without powers or factorials.
                shifted[0]=power[order-1];
                int degree=0;
                for(int k=order-2;k>=0;--k)
                {
                    for(int a=degree+1;a>=1;--a)
                    {
                        shifted[a]=shifted[a]*c+shifted[a-1];
                    }
                    shifted[0]=shifted[0]*c+power[k];++degree;
                }
                double integral=0,centered=0,second=0,dPower=1;
                for(int k=0;k<order;++k)
                {
                    if(k%2==0)
                    {
                        integral+=shifted[k]*dPower/(k+1);
                        second+=shifted[k]*dPower*d*d/(k+3);
                    }
                    else
                    {
                        centered+=shifted[k]*dPower*d/(k+2);
                    }
                    dPower*=d;
                }
                return {2*d*integral,2*d*centered,2*d*second};
            }
            // Return integral g dz and integral (z-midpoint) g dz on a positive interval.
            Eigen::Vector3d positive(double rho,double z,double length) const
            {
                rho=std::max(rho,rhoMin);
                const double radialIndex=(std::log(rho)-logMin)/dx;
                const int i=std::clamp(int(radialIndex),0,resolution-1);
                const double radial=2*(radialIndex-i)-1;
                Eigen::Vector3d result=Eigen::Vector3d::Zero();
                double traversed=0;
                int j=std::clamp(int(std::asinh(z/rho)/dt),0,resolution-1);
                while(traversed<length)
                {
                    const double lo=rho*std::sinh(j*dt),hi=rho*std::sinh((j+1)*dt);
                    const double remaining=length-traversed;
                    const double span=std::min(remaining,std::max(0.0,hi-z));
                    if(span>0)
                    {
                        const double half=(hi-lo)*.5;
                        const double start=(z-lo)/half-1;
                        const double d=span/(2*half);
                        const auto value=primitive(i,j,radial,start+d,d);
                        const double offset=traversed+span*.5-length*.5;
                        result[0]+=half*value[0];
                        result[1]+=half*half*value[1]+offset*half*value[0];
                        result[2]+=half*half*half*value[2]+2*offset*half*half*value[1]+offset*offset*half*value[0];
                        traversed+=span;z+=span;
                    }
                    if(traversed<length)
                    {
                        ++j;
                        if(j==resolution)
                        {
                            throw std::out_of_range("Diffusion segment endpoint exceeds table domain");
                        }
                    }
                }
                return result;
            }
            Eigen::Vector3d integral(double rho,double z,double length) const
            {
                if(rho>rhoMax*(1+1e-12)||std::max(std::abs(z),std::abs(z+length))>rhoMax*(1+1e-12))
                {
                    throw std::out_of_range("Diffusion segment outside prepared table domain");
                }
                if(z>=0)
                {
                    return positive(rho,z,length);
                }
                if(z+length<=0)
                {
                    Eigen::Vector3d value=positive(rho,-(z+length),length);value[1]=-value[1];return value;
                }
                const double leftLength=-z,rightLength=length+z;
                const auto left=positive(rho,0,leftLength),right=positive(rho,0,rightLength);
                return {left[0]+right[0],-left[1]-rightLength*.5*left[0]+right[1]+leftLength*.5*right[0],
                    left[2]+rightLength*left[1]+rightLength*rightLength*.25*left[0]
                    +right[2]+leftLength*right[1]+leftLength*leftLength*.25*right[0]};
            }
        };
        Eigen::Matrix3d metric_;
        double core_,factor_;
        std::shared_ptr<const Table> table_,fluxTable_;
        static std::shared_ptr<const Table> cached(double ratio,bool screened,int order,int resolution,double maximum,bool flux=false)
        {
            using Key=std::tuple<double,bool,int,int,double,bool>;
            static std::mutex mutex;
            static std::list<std::pair<Key,std::shared_ptr<const Table>>> cache;
            static std::map<Key,std::weak_ptr<const Table>> live;
            const Key key(ratio,screened,order,resolution,maximum,flux);
            std::lock_guard<std::mutex> lock(mutex);
            for(auto it=cache.begin();it!=cache.end();++it)
            {
                if(it->first==key)
                {
                    auto value=it->second;cache.splice(cache.begin(),cache,it);return value;
                }
            }
            std::shared_ptr<const Table> value;
            const auto found=live.find(key);
            if(found!=live.end())
            {
                value=found->second.lock();
            }
            if(!value)
            {
                value=std::make_shared<Table>(ratio,screened,order,resolution,maximum,flux);
                live[key]=value;
            }
            if(live.size()>64)
            {
                for(auto it=live.begin();it!=live.end();)
                {
                    if(it->second.expired())
                    {
                        it=live.erase(it);
                    }
                    else
                    {
                        ++it;
                    }
                }
            }
            cache.emplace_front(key,value);
            if(cache.size()>16)
            {
                cache.pop_back();
            }
            return value;
        }
    public:
        ConcentrationEvaluator(const Eigen::Matrix3d& diffusion,double core,double width,
            bool screened,int order,int resolution,double maxPhysicalDistance,bool includeFlux=false):core_(core)
        {
            if(!diffusion.allFinite()||!diffusion.isApprox(diffusion.transpose(),1e-13)||
               !std::isfinite(core)||core<=0||!std::isfinite(width)||width<4*core||
               order<4||order>16||resolution<4||resolution>256||!std::isfinite(maxPhysicalDistance)||maxPhysicalDistance<=0)
            {
                throw std::invalid_argument("Invalid diffusion segment table settings");
            }
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(diffusion);
            if(es.info()!=Eigen::Success||es.eigenvalues().minCoeff()<=0)
            {
                throw std::invalid_argument("Diffusion tensor must be positive definite");
            }
            const double d=std::exp(es.eigenvalues().array().log().mean());
            metric_=d*diffusion.inverse();factor_=1/(4*std::numbers::pi*d);
            const double maximum=maxPhysicalDistance*std::sqrt(d/es.eigenvalues().minCoeff())/core;
            table_=cached(width/core,screened,order,resolution,maximum);
            if(includeFlux)
            {
                fluxTable_=cached(width/core,screened,order,resolution,maximum,true);
            }
        }
        // Returns -D grad(c) for linearly varying source strengths. Number flux
        // additionally divides by omega and the source bias. No query-time quadrature.
        Eigen::Vector3d flux(const Eigen::Vector3d& p0,const Eigen::Vector3d& p1,
            const Eigen::Vector3d& x,double q0,double q1,double begin=0,double end=1) const
        {
            if(!fluxTable_)
            {
                throw std::logic_error("Construct Helmholtz evaluator with flux tables enabled");
            }
            if(!p0.allFinite()||!p1.allFinite()||!x.allFinite()||!std::isfinite(q0)||!std::isfinite(q1)||
               !(begin>=0&&end<=1&&begin<=end))
            {
                throw std::invalid_argument("Invalid Helmholtz flux segment");
            }
            const Eigen::Vector3d chord=p1-p0,r=x-p0;
            const double length=chord.norm(),lm=std::sqrt(chord.dot(metric_*chord));
            if(length==0||begin==end||(q0==0&&q1==0))
            {
                return Eigen::Vector3d::Zero();
            }
            const double center=r.dot(metric_*chord)/(lm*lm),alpha=core_/lm;
            const Eigen::Vector3d perpendicular=r-center*chord;
            const double rho=std::sqrt(std::max(0.,perpendicular.dot(metric_*perpendicular)))/core_;
            const auto moments=fluxTable_->integral(rho,(begin-center)/alpha,(end-begin)/alpha);
            const double mid=(begin+end)*.5,dq=q1-q0,qmid=q0+mid*dq;
            const Eigen::Vector3d rmid=r-mid*chord;
            return length/(4*std::numbers::pi)*alpha/(core_*core_*core_)*
                (qmid*rmid*moments[0]+alpha*(dq*rmid-qmid*chord)*moments[1]-alpha*alpha*dq*chord*moments[2]);
        }
        // Endpoint source shape integrals on parameter interval [begin,end] of p0->p1.
        Eigen::Vector2d concentrationIntegrals(const Eigen::Vector3d& p0,const Eigen::Vector3d& p1,
            const Eigen::Vector3d& x,double begin=0,double end=1) const
        {
            const Eigen::Vector3d chord=p1-p0,r=x-p0;
            const double length=chord.norm(),metricLength=std::sqrt(chord.dot(metric_*chord));
            if(length==0||begin==end)
            {
                return Eigen::Vector2d::Zero();
            }
            if(!p0.allFinite()||!p1.allFinite()||!x.allFinite()||!(begin>=0&&end<=1&&begin<=end))
            {
                throw std::invalid_argument("Invalid diffusion segment endpoints");
            }
            const double center=r.dot(metric_*chord)/(metricLength*metricLength);
            const Eigen::Vector3d perpendicular=r-center*chord;
            const double rho=std::sqrt(std::max(0.0,perpendicular.dot(metric_*perpendicular)))/core_;
            const double z=(begin-center)*metricLength/core_,delta=(end-begin)*metricLength/core_;
            const auto integral=table_->integral(rho,z,delta);
            const double midpoint=(begin+end)*.5;
            const double firstMoment=integral[1]*core_/metricLength;
            return factor_*length/metricLength*Eigen::Vector2d((1-midpoint)*integral[0]-firstMoment,
                                                               midpoint*integral[0]+firstMoment);
        }
    };

}}
#endif
