#ifndef DVEC_H_
#define DVEC_H_
#include <vector>
#include <assert.h>
#include <stdio.h>
#include <math.h>



//http://stackoverflow.com/questions/5958572/how-can-i-avoid-stdvector-to-initialize-all-its-elements
// THIS OLD FUNCTION WAS A HACK 
/*
template<typename T>
inline static T create_uninitialized(size_t size, size_t capacity) {
    T v;
#if defined(__GNUC__)
    // Don't say it. I know -_-;
    // Oddly, _M_impl is public in _Vector_base !?
    typedef typename T::value_type     value_type;
    typedef typename T::allocator_type allocator_type;
    typedef std::_Vector_base<value_type, allocator_type> base_type;
    base_type& xb(reinterpret_cast<base_type&>(v));
    value_type* p(new value_type[capacity]);
#if !defined(__EXCEPTIONS)
    size=p?size:0;         // size=0 if p is null
    capacity=p?capacity:0; // capacity=0 if p is null
#endif
    capacity=std::max(size, capacity); // ensure size<=capacity
    xb._M_impl._M_start = p;
    xb._M_impl._M_finish = p+size;
    xb._M_impl._M_end_of_storage = p+capacity;
#else
    // Fallback, for the other compilers
    capacity=std::max(size, capacity);
    v.reserve(capacity);
    v.resize(size);
#endif
    return v;
}
*/





struct Dvec
{
   std::vector<float > data;
   int min,max;
   float minval; // minimum value cache

   inline int init(int min, int max)
   {
      assert(min<max);
   	this->min = min;
   	this->max = max;
      this->data = std::vector<float >(max-min+1,0);
      this->minval = INFINITY;  // by default cache is invalid
//      this->data = create_uninitialized< std::vector<float > >( max-min+1, max-min+1);
      return 0;
   }


   inline Dvec(int min, int max)
   {
      init(min,max);
   }

   inline Dvec(int Numel) 
   {
      init(0,Numel);
   }

   inline Dvec() 
   {
   }

   inline float get_minvalue() {
      if (minval == INFINITY) 
         for(int o=min;o<=max;o++) 
            if (this->operator[](o) < minval) 
               minval=this->operator[](o);
      return minval;

   }


   inline void set(int i, float value) {
      if (i>=this->min && i<=this->max) 
      { 
         int idx = i-this->min;
//         #pragma omp critical
            data[idx]=value;
         minval=INFINITY;  // invalidate  minval cache
      }
   }

   inline void increment(int i, float value) {
      if (i>=this->min && i<=this->max) 
      { 
         int idx = i-this->min;
//         #pragma omp critical
            data[idx]+=value;
         minval=INFINITY;  // invalidate minval cache
      }
   }

   inline void set_nolock(int i, float value) {
      if (i>=this->min && i<=this->max) 
      { 
         int idx = i-this->min;
         data[idx]=value;
         minval=INFINITY;  // invalidate  minval cache
      }
   }

   inline void increment_nolock(int i, float value) {
      if (i>=this->min && i<=this->max) 
      { 
         int idx = i-this->min;
         data[idx]+=value;
         minval=INFINITY;  // invalidate minval cache
      }
   }

   inline float operator[](int i)       { if (i>=this->min && i<=this->max) return this->data[i-this->min]; else return INFINITY;}

};

//int main(){
//   struct Dvec a(10);
//   struct Dvec b(10,20);
//   for(int i=0;i<10;i++){
//      a[0]= i;
//      b[10]= i;
//   }
//   printf("%f\n", a[0]);
//   printf("%f\n", a[13]);
//   printf("%f\n", b[15]);
//}
#endif /* DVEC_H_ */
