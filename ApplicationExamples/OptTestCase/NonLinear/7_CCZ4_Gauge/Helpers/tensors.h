#ifndef SMATRIX_H
#define SMATRIX_H

/**
 * This header only classes implement Tensor contractions, Matrix
 * inversions, most importantly to raise indices with the 3-metric which
 * is given as the upper-only representation (6 components).
 * 
 * Currently, we use (again) the basically closed source Pizza Numutils
 * from https://bitbucket.org/whiskydevs/thcsupport/src/311eb8321dbf/PizzaNumUtils/?at=master
 * 
 * But we could easily switch to the open source TensorTemplates
 * https://bitbucket.org/whiskydevs/tensortemplates which were also written
 * in our group here in Frankfurt. I haven't used it here because it has
 * Boost dependencies while Pizza/NumUtils is standalone.
 * 
 * In any way all these operations here are pretty standard.
 *
 * There are a number of shortcomings of this library:
 *   1) It manages the data itself (-> needs copying)
 *   2) It cannot do arbitrary contractions but only the implemented ones.
 * 
 **/



#include <cmath>
#include <cstdlib>
#include <complex>
#include <algorithm>
#include <iostream>

namespace Tensors { /* was namespace Pizza */
	
typedef double pz_real;
typedef std::complex<pz_real> pz_cmplx;

template <class T> inline T sqr(T t) {return t*t;}


typedef pz_real sm_real;
template<int N,bool UP> class sm_tensor2_sym;

//---------------------------------------------------------------------------------------
//  Helper
//---------------------------------------------------------------------------------------
enum zero_literal {ZERO=0}; //we want to write matrix=ZERO but not matrix= 14.0
enum one_literal {ONE=1};

template<int N> class sm_array {
  typedef sm_array<N> me;
  public:
  sm_real v[N];
  sm_array(){}
  sm_real &operator[](int j) {return v[j];}
  const sm_real &operator[](int j) const {return v[j];}

  void assign_prod(const me &a,const sm_real z)
    {for(int i=0;i<N;i++) v[i]=a.v[i]*z;}
  void assign_div(const me &a,sm_real z)
    {assign_prod(a,1.0/z);}
  void assign_sum(const me &a,const me &b)
    {for(int i=0;i<N;i++) v[i]=a.v[i]+b.v[i];}
  void assign_diff(const me &a,const me &b)
    {for(int i=0;i<N;i++) v[i]=a.v[i]-b.v[i];}
  void assign_minus(const me &a)
    {for(int i=0;i<N;i++) v[i]=-a.v[i];}

  void operator+=(const me &a) {assign_sum(*this,a);}
  void operator-=(const me &a) {assign_diff(*this,a);}
  void operator*=(sm_real z) {assign_prod(*this,z);}
  void operator/=(sm_real z) {assign_div(*this,z);}

  void zero() {for(int i=0;i<N;i++) v[i]=0.0;}
};

//---------------------------------------------------------------------------------------
//  N-vector, co/contra variant
//---------------------------------------------------------------------------------------

template<int N,bool UP> class sm_tensor1 {
  typedef sm_tensor1<N,UP> me;
  public:
  sm_array<N> c;
  enum {SIZE=N};
  sm_tensor1(){}
  sm_tensor1(zero_literal z) {c.zero();}
  sm_real &operator()(int j) {return c[j];}
  const sm_real &operator()(int j) const {return c[j];}

  void assign_sum(const me &a,const me &b) {c.assign_sum(a.c,b.c);}
  void assign_diff(const me &a,const me &b) {c.assign_diff(a.c,b.c);}
  void assign_minus(const me &a) {c.assign_minus(a.c);}
  void assign_prod(const me &a,sm_real z) {c.assign_prod(a.c,z);}
  void assign_div(const me &a,sm_real z) {c.assign_div(a.c,z);}

  void operator+=(const me &a) {assign_sum(*this,a);}
  void operator-=(const me &a) {assign_diff(*this,a);}
  void operator*=(sm_real z) {assign_prod(*this,z);}
  void operator/=(sm_real z) {assign_div(*this,z);}

  me operator+(const me &a) const {me e; e.assign_sum(*this,a); return e;}
  me operator-(const me &a) const {me e; e.assign_diff(*this,a); return e;}
  me operator*(sm_real z) const {me e; e.assign_prod(*this,z); return e;}
  me operator/(sm_real z) const {me e; e.assign_div(*this,z); return e;}

  void assign_prod(const sm_tensor2_sym<N,UP> &m,const sm_tensor1<N,!UP> &w);
  sm_real norm() const;
};

template<int N,bool UP>
sm_real sm_tensor1<N,UP>::norm() const
{
  sm_real e=c[0]*c[0];
  for (int i=1;i<N;i++) e+=c[i]*c[i];
  return sqrt(e);
}

template<int N,bool UP>
sm_tensor1<N,UP> operator-(sm_tensor1<N,UP> &v)
{
  sm_tensor1<N,UP> erg;
  erg.assign_minus(v);
  return erg;
}
//---------------------------------------------------------------------------------------
// Contraction v^i w_i  resp. v_i w^i
//---------------------------------------------------------------------------------------

template<int N,bool UP>
inline sm_real operator*(const sm_tensor1<N,UP> &v,const sm_tensor1<N,!UP> &w)
{
  sm_real erg=v(0)*w(0);
  for(int i=1;i<N;i++) erg+=v(i)*w(i);
  return erg;
}

//---------------------------------------------------------------------------------------
// Scalar * vector
//---------------------------------------------------------------------------------------

template<int N,bool UP>
sm_tensor1<N,UP> inline operator*(sm_real z,const sm_tensor1<N,UP> &a) {
  sm_tensor1<N,UP> erg;
  erg.assign_prod(a,z);
  return erg;
}

//---------------------------------------------------------------------------------------
// Contraction m^ij w_j resp. m_ij w^j
//---------------------------------------------------------------------------------------

template<int N,bool UP>
inline void sm_tensor1<N,UP>::assign_prod(const sm_tensor2_sym<N,UP> &m,const sm_tensor1<N,!UP> &w)
{
  for (int i=0; i<N; i++) {
    c[i]=w(0)*m(i,0);
    for (int j=1; j<N; j++) c[i]+=w(j)*m(i,j);
  }
}

template<int N,bool UP>
inline sm_tensor1<N,UP> operator*(const sm_tensor2_sym<N,UP> &m,const sm_tensor1<N,!UP> &w)
{
  sm_tensor1<N,UP> erg;
  erg.assign_prod(m,w);
  return erg;
}

//---------------------------------------------------------------------------------------
// Symmetric tensors m^ij , m_ij
//---------------------------------------------------------------------------------------


/**
 * The indexing tested with Python:
 * > a=lambda i,j: j+(i*(i+1))/2 if j<=i else i+(j*(j+1))/2
 * > from numpy import ndindex
 * > array([a(i,j) for i,j in ndindex(3,3)]).reshape(3,3)
 * > array([[0, 1, 3],
 *          [1, 2, 4],
 *          [3, 4, 5]])
 * In contrast, we typically have the ordering
 *   gxx, gxy, gxz, gyy, gyz, gzz
 **/
template<int N,bool UP> class sm_tensor2_sym {
  typedef sm_tensor2_sym<N,UP> me;
  int index(int i,int j) const {return  j<=i ? j+(i*(i+1))/2 : i+(j*(j+1))/2;}
  public:
  enum {SIZE=(N*(N+1)/2)};
  sm_array<SIZE> c;
  sm_tensor2_sym(){}
  sm_tensor2_sym(zero_literal z) {c.zero();}
  sm_tensor2_sym(one_literal z) {diag(1.0);}

  sm_tensor2_sym(sm_real shadow[SIZE]){
    if(!shadow) return;
    at(0,0) = shadow[0];
    at(0,1) = shadow[1];
    at(0,2) = shadow[2];
    at(1,1) = shadow[3];
    at(1,2) = shadow[4];
    at(2,2) = shadow[5];
  }

  sm_real &at(int i, int j) { return c[index(i,j)];}
  const sm_real &at(int i, int j) const { return c[index(i,j)];}
  sm_real &operator()(int i,int j) {return at(i,j);}
  const sm_real &operator()(int i,int j) const {return at(i,j);}

  void assign_sum(const me &a,const me &b) {c.assign_sum(a.c,b.c);}
  void assign_diff(const me &a,const me &b) {c.assign_diff(a.c,b.c);}
  void assign_minus(const me &a) {c.assign_minus(a.c);}
  void assign_prod(const me &a,sm_real z) {c.assign_prod(a.c,z);}
  void assign_div(const me &a,sm_real z) {c.assign_div(a.c,z);}

  void operator+=(const me &a) {assign_sum(*this,a);}
  void operator-=(const me &a) {assign_diff(*this,a);}
  void operator*=(sm_real z) {assign_prod(*this,z);}
  void operator/=(sm_real z) {assign_div(*this,z);}

  me operator+(const me &a) const {me e; e.assign_sum(*this,a); return e;}
  me operator-(const me &a) const {me e; e.assign_diff(*this,a); return e;}
  me operator*(sm_real z) const {me e; e.assign_prod(*this,z); return e;}
  me operator/(sm_real z) const {me e; e.assign_div(*this,z); return e;}
  me &diag(sm_real d);

  sm_real contract(const sm_tensor1<N,!UP> &vl,const sm_tensor1<N,!UP> &vr) const;
  sm_real contract(const sm_tensor1<N,!UP> &v) const {return contract(v,v);}
  sm_real contract(const sm_tensor2_sym<N,!UP> &m) const;
};

template<int N,bool UP>
sm_tensor2_sym<N,UP> &sm_tensor2_sym<N,UP>::diag(sm_real d) {
  for (int i=0;i<N;i++) {
    (*this)(i,i)=d;
    for (int j=0;j<i;j++) (*this)(i,j)=0.0;
  }
  return *this;
}

template<int N,bool UP>
sm_tensor2_sym<N,UP> operator-(sm_tensor2_sym<N,UP> &v)
{
  sm_tensor2_sym<N,UP> erg;
  erg.assign_minus(v);
  return erg;
}
//---------------------------------------------------------------------------------------
// Scalar * tensor
//---------------------------------------------------------------------------------------

template<int N,bool UP>
inline sm_tensor2_sym<N,UP> operator*(sm_real a,const sm_tensor2_sym<N,UP> &m1){
  sm_tensor2_sym<N,UP> erg;
  erg.assign_prod(m1,a);
  return erg;
}


//---------------------------------------------------------------------------------------
// Contraction v^i w^i m_ij
//---------------------------------------------------------------------------------------

template<int N,bool UP>
inline sm_real sm_tensor2_sym<N,UP>::contract(const sm_tensor1<N,!UP> &v,const sm_tensor1<N,!UP> &w) const
{
  sm_real erg,a;
  int i,j;

  a=v(0)*(*this)(0,0);
  for (j=1;j<N;j++) a+=v(j) * (*this)(0,j);
  erg=w(0)*a;

  for (i=1; i<N;i++) {
    a=v(0)*(*this)(i,0);
    for (j=1;j<N;j++) a+=v(j) * (*this)(i,j);
    erg+=w(i)*a;
  }
  return erg;
}
//---------------------------------------------------------------------------------------
// v_i w_j
//---------------------------------------------------------------------------------------
template<int N,bool UP>
sm_tensor2_sym<N,UP> operator*(const sm_tensor1<N,UP> &v,const sm_tensor1<N,UP> &w)
{
  sm_tensor2_sym<N,UP> erg;
  for (int i=0;i<N;i++) for (int j=0;j<=i;j++) erg(i,j)=v(i)*w(j);
  return erg;
}

//---------------------------------------------------------------------------------------
// Contraction v^ij w_ij (by Sven)
//---------------------------------------------------------------------------------------
template<int N, bool UP>
inline sm_real sm_tensor2_sym<N,UP>::contract(const sm_tensor2_sym<N,!UP> &m) const {
	sm_real erg = 0;
	for(int i=0;i<N;i++) for(int j=0;j<N;j++) erg += (*this)(i,j) * m(i,j);
	return erg;
}


//---------------------------------------------------------------------------------------
// Determinant
//---------------------------------------------------------------------------------------
template<bool UP>
inline sm_real det(const sm_tensor2_sym<3,UP> &lo)
{
  sm_real d =
    -lo(0,2)*lo(0,2)*lo(1,1)
    + 2*lo(0,1)*lo(0,2)*lo(1,2)
    - lo(0,0)*lo(1,2)*lo(1,2)
    - lo(0,1)*lo(0,1)*lo(2,2)
    + lo(0,0)*lo(1,1)*lo(2,2);

  return d;
}


//---------------------------------------------------------------------------------------
// Collection of m^ij & m_ij
//---------------------------------------------------------------------------------------

template<int N> struct sm_tensors2_sym {
  sm_tensor2_sym<N,false> lo;
  sm_tensor2_sym<N,true> up;

  sm_tensors2_sym(){}
  sm_tensors2_sym(sm_real *lobuf, sm_real *upbuf) : lo(lobuf), up(upbuf) {}


  sm_tensor1<N,true> operator*(const sm_tensor1<N,false> &v) const {return up*v;}
  sm_tensor1<N,false> operator*(const sm_tensor1<N,true> &v) const {return lo*v;}
  
  // for metric, actually computing the scalar product
  sm_real contract(const sm_tensor1<N,true> &v,const sm_tensor1<N,true> &w) const
    { return lo.contract(v,w);}
  sm_real contract(const sm_tensor1<N,false> &v,const sm_tensor1<N,false> &w) const
    { return up.contract(v,w);}

  // for metric, ctually computing the norm
  sm_real contract(const sm_tensor1<N,true> &v) const {return lo.contract(v);}
  sm_real contract(const sm_tensor1<N,false> &v) const {return up.contract(v);}
  
  // for mteric, actually computing the traces of m
  sm_real contract(const sm_tensor2_sym<N,true> &m) const { return lo.contract(m); }
  sm_real contract(const sm_tensor2_sym<N,false> &m) const { return up.contract(m); }
};

//---------------------------------------------------------------------------------------
// Metric
//---------------------------------------------------------------------------------------


template<int N> struct sm_metric : public sm_tensors2_sym<N> {
  sm_real vc;

  sm_metric(){}
  sm_metric(sm_real *shadow) : sm_tensors2_sym<N>(shadow,NULL) { complete_from_lower(); }
  sm_real dv() const {return vc;}
  template<bool UP> sm_real norm2(const sm_tensor1<N,UP> &v) const
    {return fabs(contract(v,v));}
  template<bool UP> sm_real norm(const sm_tensor1<N,UP> &v) const
    {return sqrt(norm2(v));}
  template<bool UP> sm_real trace(const sm_tensor2_sym<N,UP> &m) const
    {return this->contract(m);}
  void complete_from_lower();
};

template<>
inline void sm_metric<3>::complete_from_lower()
{
  sm_real d=det(lo);

  up(0,0) = (-lo(1,2)*lo(1,2) + lo(1,1)*lo(2,2) );
  up(0,1) = ( lo(0,2)*lo(1,2) - lo(0,1)*lo(2,2) );
  up(1,1) = (-lo(0,2)*lo(0,2) + lo(0,0)*lo(2,2) );
  up(0,2) = (-lo(0,2)*lo(1,1) + lo(0,1)*lo(1,2) );
  up(1,2) = ( lo(0,1)*lo(0,2) - lo(0,0)*lo(1,2) );
  up(2,2) = (-lo(0,1)*lo(0,1) + lo(0,0)*lo(1,1) );

  up/=d;
  vc=sqrt(fabs(d));
}

//---------------------------------------------------------------------------------------
//Ausgabe-Operatoren
//---------------------------------------------------------------------------------------

template<int N,bool UP>
std::ostream& operator<<(std::ostream &s,const sm_tensor1<N,UP> &v) {
  s << "( " << v(0);
  for (int i=1;i<N;i++) s << ", " <<v(i);
  return s << ")";
}

template<int N,bool UP>
std::ostream& operator<<(std::ostream &s,const sm_tensor2_sym<N,UP> &m) {
  s<<std::endl;
  for (int i=0;i<N;i++) {
    for (int j=0;j<N;j++) s<<"  "<<m(i,j);
    s<<std::endl;
  }
  s<<std::endl;
  return s;
}


//---------------------------------------------------------------------------------------
//Standard 3D variables
//---------------------------------------------------------------------------------------

typedef sm_tensor1<3,true> vec_u;
typedef sm_tensor1<3,false> vec_l;
typedef sm_tensor2_sym<3,true> mats_u;
typedef sm_tensor2_sym<3,false> mats_l;
typedef sm_metric<3> metric;



} // namespace

#endif
