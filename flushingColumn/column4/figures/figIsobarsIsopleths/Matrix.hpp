// Program to test slices and a simple N*M matrix class

// pp 670-674 and 683-684

// No guarantees offered. Constructive comments to bs@research.att.com

#ifndef MATRIX_HPP
#define MATRIX_HPP

#include<iostream>
#include<fstream>
#include<valarray>
#include<algorithm>
#include<numeric>	// for inner_product

#include "Slice_iterator.hpp"

using namespace std;

/**
 * @brief A simple template class to store dense matrix.
 * This class is based on the example in Stroustrup's "" and freely available on the net. \\
 * Basically, the data are store in a sequential valarray in fortran style (and hence by column), then suitable methods allow to access directly such data (via Slice_iter in iterator style) and to extract them in new valarray.\\
 * @note no dimensions checking are performed 
 * @todo: consider using C-style instead of Fortran...
 * @author Andrea Cassioli
 */
template<class Type>
class Matrix {
	
  valarray<Type>* v; /// stores elements by column [00,10,20..01,11,21,,]
  size_t d1;	/// column lenght = # of rows
  size_t d2;	/// row lenght = # of columns


public:
	
  Matrix(){

    d1=0;
    d2=0;
    v=new valarray<Type>(0);
  }

  /**	
     @param[in] x number of rows
     @param[out] y number of columns 
  */	
  Matrix(size_t x, size_t y){
    d1 = x;
    d2 = y;
    v = new valarray<Type>(x*y);
  }

  Matrix(size_t x, size_t y, Type dv){
    d1 = x;
    d2 = y;
    v = new valarray<Type>(dv,x*y);

  }

  //!  Copy constructor avoiding self-copying
  Matrix(const Matrix<Type>& m){
    if(this!=&m){
      d1=m.dim1();
      d2=m.dim2();
      v=new valarray<Type>(m.array());

    }
  }

  //! operator= avoiding self-copying
  Matrix<Type>& operator=(Matrix<Type>& m){
    if(this!=&m){
      d1=m.dim1();
      d2=m.dim2();
      v->resize(d1*d2);
      *v=m.array();
    }
  }

  /// operator= avoiding self-copying
  const Matrix<Type>& operator=(const Matrix<Type>& m){
    if(this!=&m){
      d1=m.dim1();
      d2=m.dim2();
      v->resize(d1*d2);
      *v=m.array();
    }
  }

  //operator= to assign a Type value to all the matriv elements
  Matrix<Type>& operator=(Type d){
    *v=d;
  }

  ~Matrix(){
    delete v;
  }

  //--------------------------
  //dimension handling
		
  size_t size() const { return d1*d2; }
  size_t dim1() const { return d1; }
  size_t dim2() const { return d2; }

  size_t nRows() const {return d1;}
  size_t nCols() const {return d2;}

  ///resizing method: note values are discarded!
  /// @todo we should avoid useless resizing!
  void resize(size_t i, size_t j){
    d1=i;
    d2=j;
    v->resize(d1*d2);
  }
  ///resizing method: note values are discarded!
  /// @todo we should avoid useless resizing!
  void resize(size_t i, size_t j, Type dv){
    d1=i;
    d2=j;
    v->resize(d1*d2, dv);
  }



  //--------------------------
  //Rows access methods
  Slice_iter<Type> row_iter(size_t i){
    return Slice_iter<Type>(v,slice(i,d2,d1));
  }
  Cslice_iter<Type> row_iter(size_t i) const{
    return Cslice_iter<Type>(v,slice(i,d2,d1));
  }
  valarray<Type> row(size_t i){
    valarray<Type> s((*v)[slice(i,d2,d1)]);
    return s;
  }
  const valarray<Type> row(size_t i)const{
    valarray<Type> s((*v)[slice(i,d2,d1)]);
    return s;
  }
		
  //--------------------------
  //column access methods
  Slice_iter<Type> column_iter(size_t i){
    return Slice_iter<Type>(v,slice(i*d1,d1,1));
  }
  Cslice_iter<Type> column_iter(size_t i) const{
    return Cslice_iter<Type>(v,slice(i*d1,d1,1));
  }
  ///a copy of the i-th column is returned
  valarray<Type> column(size_t i){
    valarray<Type> s((*v)[slice(i*d1,d1,1)]);
    /*cout <<"in column("<<i<<")"<<endl;
      cout <<*v;
      cout <<"with d2="<<d2<<" and sum="<<s.sum()<<endl;
      cout <<s;
      cout<<" from matrix"<<endl;
      cout<<*this;*/
    return s;
  }
  ///a const copy of the i-th column is returned
  const valarray<Type> column(size_t i)const{
    valarray<Type> s((*v)[slice(i*d1,d1,1)]);
		
    //cout <<"in column("<<i<<")"<<endl;
    //cout <<"with d2="<<d2<<endl;
    //cout <<s;

    return s;
  }

  //--------------------------
  //element access method
  ///return an iterator-like to the i-th row
  Slice_iter<Type> operator[](size_t i) { return row_iter(i); }	// C-style subscript

  ///return a const  iterator-like to the i-th row
  Cslice_iter<Type> operator[](size_t i) const { return row_iter(i); }

  ///return the (i,j)-th element in C-style subscript
  Type operator()(size_t i, size_t j)const{
    return row_iter(i)[j];
  }
  Type& operator()(size_t i, size_t j){
    return row_iter(i)[j];
  }
	
  //--------------------------
  //whole data access method
  const valarray<Type>& array()const {
    return *v;
  }
  valarray<Type>& array() {
    return *v;
  }

  //operators



  // self applied operator
  ///multiplies all matrix elements  for a given value
  Matrix<Type>& operator*=(Type x){
    *v *= x;
    return *this;
  };

	
  ///adds to all matrix elements a given value
  Matrix<Type>& operator+=(Type x){
    *v += x;
    return *this;
  };

  ///adds (elementwise) to the  matrix a given one
  Matrix<Type>& operator+=(Matrix<Type>& x){
    *v += x.array();
    return *this;
  };

  ///adds (elementwise) to the  matrix a given one
  const Matrix<Type>& operator+=(const Matrix<Type>& x){
    *v += x.array();
    return *this;
  };

	
  ///transpose the matrix
  //assuming 
  void t(Matrix<Type>& a){
    resize(a.nCols(),a.nRows());
    for (int i=0; i<d1; i++)//rows
      for (int j=0; j<d2; j++)//cols
	(*this)[i][j] = a[j][i];
  }	
	
  ///matrix inverter
  void inv(Matrix<Type>& a){
    int i,icol,irow,j,k,l,ll;
    double big,dum,pivinv;
    int n=a.nRows();
    int m=a.nCols();
    //cerr<<"matrix::inv"<<endl; 
    //cerr<<"n="<<d1<<" m="<<d2<<endl;
    //cerr<<"matrix.sise="<<v->size()<<endl;
    (*this)=a;

    valarray<int> indxc(n),indxr(n),ipiv(0,n);

    //cerr<<"gauss first loop"<<endl;
    for (i=0;i<n;i++)
      {
	big=0.0;
	for (j=0;j<n;j++)
	  if (ipiv[j] != 1)
	    for (k=0;k<n;k++)
	      {
		if (ipiv[k] == 0)
		  {
		    if (fabs((*this)[j][k]) >= big)
		      {
			big=fabs((*this)[j][k]);
			irow=j;
			icol=k;
		      }
		  }
	      }

	++(ipiv[icol]);
	if (irow != icol)
	  for (l=0;l<n;l++)
	    swap((*this)[irow][l],(*this)[icol][l]);
	indxr[i]=irow;
	indxc[i]=icol;
	if ((*this)[icol][icol] == 0.0)
	  {
	    cerr<<"GaussJordan::inv Singular Matrix"<<endl;
	    exit(-1);
	  }
	pivinv=1.0/(*this)[icol][icol];
	(*this)[icol][icol]=1.0;
	for (l=0;l<n;l++)
	  (*this)[icol][l] *= pivinv;
	for (ll=0;ll<n;ll++)
	  if (ll != icol)
	    {
	      dum=(*this)[ll][icol];
	      (*this)[ll][icol]=0.0;
	      for (l=0;l<n;l++)
		(*this)[ll][l] -= (*this)[icol][l]*dum;
	    }
      }

    //        	cerr<<"gauss-jordan last loop"<<endl;
    for (l=n-1;l>=0;l--)
      {
	if (indxr[l] != indxc[l])
	  for (k=0;k<n;k++)
	    swap((*this)[k][indxr[l]],(*this)[k][indxc[l]]);
      }

  }
};


//-------------- External operators -----------------

// inner product (practically dot product)
/// a const iterator times a valarray 

/// Cslice_iter<Type> elements times a Type value, returning a 
template<class Type> 
valarray<Type> operator*(Cslice_iter<Type>& v, const Type d)
{
  Cslice_iter<Type> iter=v.begin();
  valarray<Type> res(v.end()-v.begin());
  for (; iter!=v.end(); iter++)
    res[iter-iter.begin()] = (*iter)*d;
  return res;
}

/// Cslice_iter<Type> elements times a Type value
template<class Type> 
valarray<Type> operator*(Slice_iter<Type>& v, Type d)
{
  Slice_iter<Type> iter=v.begin();
  valarray<Type> res(v.end()-v.begin());
  for (; iter!=v.end(); iter++)
    res[iter-iter.begin()] = (*iter)*d;
  return res;
}

template<class Type>
Type mul(const Slice_iter<Type>& v1, const valarray<Type>& v2)
{
  Type res = 0;
  res=inner_product(v1.begin(),v1.end(),&v2[0],Type(0));
  return res;
}


template<class Type>
Type mul(const Cslice_iter<Type>& v1, const valarray<Type>& v2)
{
  Type res = 0;
  for (size_t i = 0; i<v2.size(); i++)
    res+= v1[i]*v2[i];
  return res;
}

template<class Type> 
valarray<Type> operator*(const Matrix<Type>& m, const valarray<Type>& v)
{
  valarray<Type> res(m.nRows());
  for (size_t i = 0; i<m.nRows(); i++)
    res[i] = mul(m.row_iter(i),v);
  return res;
}


template<class Type> 
Matrix<Type> operator*(Matrix<Type>& m, Type v)
{
  Matrix<Type> res(m);
  res*=v;
  return res;
}



template<class Type> 
valarray<Type> operator*(Matrix<Type>& m, valarray<Type>& v)
{
  valarray<Type> res(m.nRows());
  for (size_t i = 0; i<m.nRows(); i++){
    const Slice_iter<Type>& ri=m.row_iter(i);
    res[i] = inner_product(ri.begin(),ri.end(),&v[0],Type(0));	
  }
  return res;
}



// alternative definition of m*v

//valarray<Type> operator*(const Matrix& m, valarray<Type>& v)
template<class Type> 
valarray<Type> mul_mv(const Matrix<Type>& m, valarray<Type>& v)
{

  valarray<Type> res(m.nRows());

  for (size_t i = 0; i<m.Rows(); i++) {
    const Cslice_iter<Type>& ri = m.rowi_iter(i);
    res[i] = inner_product(ri.begin(),ri.end(),&v[0],Type(0));
  }
  return res;
}


template<class Type>
valarray<Type> operator*(valarray<Type>& v, const Matrix<Type>& m)
{

  valarray<Type> res(m.dim1());

  for (size_t i = 0; i<m.nRows(); i++) {
    const Cslice_iter<Type>& ci = m.column(i);
    res[i] = inner_product(ci,ci.end(),&v[0],Type(0));
  }
  return res;
}


template<class Type> 
Matrix<Type> operator+(const Matrix<Type>& m, const Matrix<Type>& v)
{
  Matrix<Type> res(m);
  res+=v;
  return res;
}


template<class Type> 
Matrix<Type> operator*(const Matrix<Type>& m, const Matrix<Type>& v)
{
  Matrix<Type> res(m.nRows(),v.nCols(),0);
  //	cout<<"multiplying m("<<m.nRows()<<";"<<m.nCols()<<") by v("<<v.nRows()<<";"<<v.nCols()<<")"<<endl;
  Type app;
	
  for(size_t i=0;i<m.nRows();i++){
    for (size_t j = 0; j<v.nCols(); j++){
      /*			cout<<"row:"<<endl;
	cout<<m.row_iter(i);
	cout<<"col:"<<endl;
	cout<<v.column_iter(j);*/
      app= inner_product(m.row_iter(i).begin(),m.row_iter(i).end(),v.column_iter(j).begin(), Type(0));
      //			cout<<"res("<<i<<";"<<j<<")"<<"="<<app<<endl;
      res(i,j)=app;
    }
  }
  return res;
}

/// stream operator for matrix
template<class Type>
ostream& operator<<(ostream& os, Matrix<Type>& m)
{
  for(int y=0; y<m.nRows(); y++)
    {
      for(int x=0; x<m.nCols(); x++)
	os<<m[y][x]<<"\t";
      os << "\n";
    }
  return os;
}
/// istream operator for matrix
template<class Type>
std::istream& operator>>(std::istream& os, Matrix<Type>& m)
{
  cout<<"READING A MATRIX"<<endl;
  for(int y=0; y<m.nRows(); y++)
    {
      for(int x=0; x<m.nCols(); x++)
	os>>m[y][x];
                
    }
  return os;
}
/// input stream operator for valarray
template<class Type>
std::istream& operator>>(std::istream& os, valarray<Type>& v)
{
  Type t;
  for (int i = 0; i<v.size(); ++i){
    os >> t;
    v[i]=t;
  }
  return os;
}

/// stream operator for valarray
template<class Type>
ostream& operator<<(ostream& os, const valarray<Type>& v)
{
  os<<endl;
  for (int i = 0; i<v.size(); ++i)
    os << '\t' << v[i]<<endl;
  return os;
}

/// stream operator for Cslice_iter
template<class Type>
ostream& operator<<(ostream& os,const  Cslice_iter<Type>& v)
{
  for (int i = 0; i<v.size(); ++i)
    os << '\t' << v[i];
  os<<endl;
  return os;
}

/// stream operator for Slice_iter
template<class Type>
ostream& operator<<(ostream& os, const Slice_iter<Type>& v)
{
  for (int i = 0; i<v.size(); ++i)
    os << '\t' << v[i];

  os<<endl;
  return os;
}

/**
   a method used to read a LookUp Table from the given string checking dimensions.
   @author Andrea Cassioli
   @return true if operation ends correctly
*/
template<class Type>
bool readTable(ifstream& s, const string &tablename, valarray<Type>& table)
{
  int dim;
  cerr<<" reading table "<< tablename <<endl;
  s.open(tablename.c_str());
  if(!s)
    return false;
  s>>dim;
  cerr<<dim<<" elements need to be read"<<endl;
  if(dim==0)
    {
      s.close();
      return false;
    }
  table.resize(dim);
  for (int i=0; i<table.size(); i++)
    s >> table[i];
  cerr<<"reading complete"<<endl;
  s.close();
  return true;
}

#endif

