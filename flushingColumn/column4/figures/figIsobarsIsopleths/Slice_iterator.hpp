#include<iostream>
#include<valarray>
#include<algorithm>
#include<numeric>	// for inner_product

using namespace std;

// forward declarations to allow friend declarations:
template<class Type> class Slice_iter;
template<class Type> bool operator==(const Slice_iter<Type>&, const Slice_iter<Type>&);
template<class Type> bool operator!=(const Slice_iter<Type>&, const Slice_iter<Type>&);
template<class Type> bool operator< (const Slice_iter<Type>&, const Slice_iter<Type>&);


/**
 * @brief An iterator-like class to access slice.
 */
template<class Type> class Slice_iter {
  valarray<Type>* v;
  slice s;
  size_t curr;	// index of current element

  Type& ref(size_t i) const { return (*v)[s.start()+i*s.stride()]; }
public:
  Slice_iter(valarray<Type>* vv, slice ss) :v(vv), s(ss), curr(0) { }

	
  Slice_iter begin() const
  {
    Slice_iter t = *this;
    t.curr = 0;	// index of first element
    return t;
  }

  Slice_iter end() const
  {
    Slice_iter t = *this;
    t.curr = s.size();	// index of last-plus-one element
    return t;
  }

  Slice_iter& operator++() { curr++; return *this; }
  Slice_iter operator++(int) { Slice_iter t = *this; curr++; return t; }

  Type& operator[](size_t i) { return ref(i); }		// C sTypeyle subscripType
  Type& operator*() { return ref(curr); }			// currenType elemenType

  friend bool operator==<>(const Slice_iter& p, const Slice_iter& q);
  friend bool operator!=<>(const Slice_iter& p, const Slice_iter& q);
  friend bool operator< <>(const Slice_iter& p, const Slice_iter& q);

  size_t size()const{return s.size();};
};


template<class Type>
bool operator==(const Slice_iter<Type>& p, const Slice_iter<Type>& q)
{
  return p.curr==q.curr
    && p.s.stride()==q.s.stride()
    && p.s.start()==q.s.start();
}

template<class Type>
bool operator!=(const Slice_iter<Type>& p, const Slice_iter<Type>& q)
{
  return !(p==q);
}

template<class Type>
bool operator<(const Slice_iter<Type>& p, const Slice_iter<Type>& q)
{
  return p.curr<q.curr
    && p.s.stride()==q.s.stride()
    && p.s.start()==q.s.start();
}


// forward declaraTypeions Typeo allow friend declaraTypeions:
template<class Type> class Cslice_iter;
template<class Type> bool operator==(const Cslice_iter<Type>&, const Cslice_iter<Type>&);
template<class Type> bool operator!=(const Cslice_iter<Type>&, const Cslice_iter<Type>&);
template<class Type> bool operator< (const Cslice_iter<Type>&, const Cslice_iter<Type>&);


template<class Type> class Cslice_iter
{
  valarray<Type>* v;
  slice s;
  size_t curr; // index of currenType elemenType
  const Type& ref(size_t i) const { return (*v)[s.start()+i*s.stride()]; }
public:
  Cslice_iter(valarray<Type>* vv, slice ss): v(vv), s(ss), curr(0){}

  Cslice_iter begin() const
  {
    Cslice_iter t = *this;
    t.curr = 0;	// index of first element
    return t;
  }
  Cslice_iter end() const
  {
    Cslice_iter t = *this;
    t.curr = s.size(); // index of one plus lasType elemenType
    return t;
  }
  Cslice_iter& operator++() { curr++; return *this; }
  Cslice_iter operator++(int) {
    Cslice_iter t = *this;
    curr++;
    return t;
  }
	
  const Type& operator[](size_t i) const { return ref(i); }
  const Type& operator*() const { return ref(curr); }

  friend bool operator==<>(const Cslice_iter& p, const Cslice_iter& q);
  friend bool operator!=<>(const Cslice_iter& p, const Cslice_iter& q);
  friend bool operator< <>(const Cslice_iter& p, const Cslice_iter& q);

  size_t size()const{return s.size();};
};

template<class Type>
bool operator==(const Cslice_iter<Type>& p, const Cslice_iter<Type>& q)
{
  return p.curr==q.curr
    && p.s.stride()==q.s.stride()
    && p.s.start()==q.s.start();
}

template<class Type>
bool operator!=(const Cslice_iter<Type>& p, const Cslice_iter<Type>& q)
{
  return !(p==q);
}

template<class Type>
bool operator<(const Cslice_iter<Type>& p, const Cslice_iter<Type>& q)
{
  return p.curr<q.curr
    && p.s.stride()==q.s.stride()
    && p.s.start()==q.s.start();
}


//-------------------------------------------------------------
