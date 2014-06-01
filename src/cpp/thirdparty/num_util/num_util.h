#ifndef NUM_UTIL_H__
#define NUM_UTIL_H__

#if defined(USING_NUMARRAY)

#include <boost/python.hpp>
//#include <Numeric/arrayobject.h>
#include <numarray/arrayobject.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>


namespace num_util{
  //!
  /**
   *A free function that extracts a PyArrayObject from any sequential PyObject.
   *@param x a sequential PyObject wrapped in a Boost/Python 'object'.
   *@return a PyArrayObject wrapped in Boost/Python numeric array.
   */
  boost::python::numeric::array makeNum(boost::python::object x);

  /** 
   *Creates an one-dimensional Numeric array of length n and Numeric type t. 
   * The elements of the array are initialized to zero.
   *@param n an integer representing the length of the array.
   *@param t elements' Numeric type. Default is double.
   *@return a numeric array of size n with elements initialized to zero.
   */
  boost::python::numeric::array makeNum(int n, PyArray_TYPES t);

  /** 
   *Creates a n-dimensional Numeric array with dimensions dimens and Numeric 
   *type t. The elements of the array are initialized to zero.
   *@param dimens a vector of interger specifies the dimensions of the array.
   *@param t elements' Numeric type. Default is double.
   *@return a numeric array of shape dimens with elements initialized to zero.
   */
  boost::python::numeric::array makeNum(std::vector<int> dimens, 
					PyArray_TYPES t);
				      
  /** 
   *Creates an one-dimensional Numeric array from an integer array of size n.
   *@param data an array of integer.
   *@param n an integer indicates the size of the array.
   *@return a Numeric array of size n with elements initialized to data.
   */
  boost::python::numeric::array makeNum(int * data, int n);

  /** 
   *Creates an one-dimensional Numeric array from a float array of size n.
   *@param data an array of float.
   *@param n an integer indicates the size of the array.
   *@return a Numeric array of size n with elements initialized to data.
   */
  boost::python::numeric::array makeNum(float * data, int n);

  /** 
   *Creates an one-dimensional Numeric array from a double array of size n.
   *@param data an array of double.
   *@param n an integer indicates the size of the array.
   *@return a Numeric array of size n with elements initialized to data.
  */
  boost::python::numeric::array makeNum(double * data, int n);

  /** 
   *Creates a n-dimensional Numeric array with dimensions dims and an array
   *of 8-bits.
   *@param data an array of characters.
   *@param dims a vector of interger to represent the dimensions of the array.
   *@return a numeric array of shape dims with elements initialized to data.
  */
  boost::python::numeric::array makeNum(char * data, std::vector<int> dims);

  /** 
   *Creates a n-dimensional Numeric array with dimensions dims and an array
   *of 16-bits.
   *@param data an array of short integers.
   *@param dims a vector of interger to represent the dimensions of the array.
   *@return a numeric array of shape dims with elements initialized to data.
  */
  boost::python::numeric::array makeNum(short * data, std::vector<int> dims);

  /** 
   *Creates a n-dimensional Numeric array with dimensions dims and an array
   *of integer.
   *@param data an array of integer.
   *@param dims a vector of interger to represent the dimensions of the array.
   *@return a numeric array of shape dims with elements initialized to data.
  */
  boost::python::numeric::array makeNum(int * data, std::vector<int> dims);

  /** 
   *Creates a n-dimensional Numeric array with dimensions dims and an array
   *of float number.
   *@param data an array of float.
   *@param dims a vector of interger to represent the dimensions of the array.
   *@return a numeric array of shape dims with elements initialized to data.
   */
  boost::python::numeric::array makeNum(float * data, std::vector<int> dims);

  /** 
   *Creates a n-dimensional Numeric array with dimensions dims and an array
   *of double.
   *@param data an array of double.
   *@param dims a vector of interger to represent the dimensions of the array.
   *@return a numeric array of shape dims with elements initialized to data.
   */
  boost::python::numeric::array makeNum(double * data, std::vector<int> dims);

  /** 
   *Creates a Numeric array from a Numeric array, referencing the data.
   *@param arr a Boost/Python numeric array.
   *@return a numeric array referencing the input array.
   */
  boost::python::numeric::array makeNum(const 
					boost::python::numeric::array& arr);

  /** 
   *A free function that retrieves the Numeric type of a Numeric array.
   *@param arr a Boost/Python numeric array.
   *@return the Numeric type of the array's elements 
   */
  PyArray_TYPES type(boost::python::numeric::array arr);

  /** 
   *Throws an exception if the actual array type is not equal to the expected 
   *type.
   *@param arr a Boost/Python numeric array.
   *@param expected_type an expected Numeric type.
   *@return -----
   */
  void check_type(boost::python::numeric::array arr, 
		  PyArray_TYPES expected_type);

  /** 
   *A free function that retrieves the number of dimensions of a Numeric array.
   *@param arr a Boost/Python numeric array.
   *@return an integer that indicates the rank of an array.
   */
  int rank(boost::python::numeric::array arr);

  /** 
   *Throws an exception if the actual rank is not equal to the expected rank.
   *@param arr a Boost/Python numeric array.
   *@param expected_rank an expected rank of the numeric array.
   *@return -----
   */
  void check_rank(boost::python::numeric::array arr, int expected_rank);
  
  /** 
   *A free function that returns the total size of the array.
   *@param arr a Boost/Python numeric array.
   *@return an integer that indicates the total size of the array.
   */
  int size(boost::python::numeric::array arr);
  
  /** 
   *Throw an exception if the actual total size of the array is not equal to 
   *the expected size.
   *@param arr a Boost/Python numeric array.
   *@param expected_size the expected size of an array.
   *@return -----
   */
  void check_size(boost::python::numeric::array arr, int expected_size);

  /** 
   *Returns the dimensions in a vector.
   *@param arr a Boost/Python numeric array.
   *@return a vector with integer values that indicates the shape of the array.
  */
  std::vector<int> shape(boost::python::numeric::array arr);

  /**
   *Returns the size of a specific dimension.
   *@param arr a Boost/Python numeric array.
   *@param dimnum an integer that identifies the dimension to retrieve.
   *@return the size of the requested dimension.
   */
  int get_dim(boost::python::numeric::array arr, int dimnum);

  /** 
   *Throws an exception if the actual dimensions of the array are not equal to
   *the expected dimensions.
   *@param arr a Boost/Python numeric array.
   *@param expected_dims an integer vector of expected dimension.
   *@return -----
   */
  void check_shape(boost::python::numeric::array arr, 
		   std::vector<int> expected_dims);

  /**
   *Throws an exception if a specific dimension from a Numeric array does not
   *match the expected size.
   *@param arr a Boost/Python numeric array.
   *@param dimnum an integer that specifies which dimension of 'arr' to check.
   *@param dimsize an expected size of the specified dimension.
   *@return -----
  */
  void check_dim(boost::python::numeric::array arr, int dimnum, int dimsize);

  /** 
   *Returns true if the array is contiguous.
   *@param arr a Boost/Python numeric array.
   *@return true if the array is contiguous, false otherwise.
  */
  bool iscontiguous(boost::python::numeric::array arr);

  /** 
   *Throws an exception if the array is not contiguous.
   *@param arr a Boost/Python numeric array.
   *@return -----
  */
  void check_contiguous(boost::python::numeric::array arr);

  /** 
   *Returns a pointer to the data in the array.
   *@param arr a Boost/Python numeric array.
   *@return a char pointer pointing at the first element of the array.
   */
  char* data(boost::python::numeric::array arr);

  /** 
   *Copies data into the array.
   *@param arr a Boost/Python numeric array.
   *@param new_data a char pointer referencing the new data.
   *@return -----
   */
  void copy_data(boost::python::numeric::array arr, char* new_data);
  
  /** 
   *Returns a clone of this array.
   *@param arr a Boost/Python numeric array.
   *@return a replicate of the Boost/Python numeric array.
   */
  boost::python::numeric::array clone(boost::python::numeric::array arr);
  
  /** 
   *Returns a clone of this array with a new type.
   *@param arr a Boost/Python numeric array.
   *@param t PyArray_TYPES of the output array.
   *@return a replicate of 'arr' with type set to 't'.
   */
  boost::python::numeric::array astype(boost::python::numeric::array arr, 
				       PyArray_TYPES t);

  /** 
   *Returns true if the array is a spacesaver array.
   *@param arr a Boost/Python numeric array.
   *@return true if the array is a spacesaver array, false otherwise.
   */
  bool spacesaver(boost::python::numeric::array arr);
  
  /** 
   *Changes the savespace property of the array.
   *@param arr a Boost/Python numeric array.
   *@param set_savespace a boolean value. Default is true.
   *@return -----
   */
  void savespace(boost::python::numeric::array arr, bool set_savespace);

  /** 
   *Returns the reference count of the array.
   *@param arr a Boost/Python numeric array.
   *@return the reference count of the array.
   */
  int refcount(boost::python::numeric::array arr);

  /** 
   *Returns the strides array in a vector of integer.
   *@param arr a Boost/Python numeric array.
   *@return the strides of an array.
   */
  std::vector<int> strides(boost::python::numeric::array arr);

  /** 
   *Throws an exception if the element of a Numeric array is type casted to
   *PyArray_OBJECT.
   *@param newo a Boost/Python object.
   *@return -----
   */
  void check_PyArrayElementType(boost::python::object newo);

  /** 
   *Mapping from a PyArray_TYPE to its corresponding name in string.
   */
  typedef std::map<PyArray_TYPES, std::string> KindStringMap;

  /** 
   *Mapping from a PyArray_TYPE to its corresponding typeID in char.
   */
  typedef std::map<PyArray_TYPES, char> KindCharMap;

  /** 
   *Mapping from a typeID to its corresponding PyArray_TYPE.
   */
  typedef std::map<char, PyArray_TYPES> KindTypeMap;

  /** 
   *Converts a PyArray_TYPE to its name in string.
   *@param t_type a PyArray_TYPES.
   *@return the corresponding name in string.
   */
  std::string type2string(PyArray_TYPES t_type);

  /** 
   *Converts a PyArray_TYPE to its single character typecode.
   *@param t_type a PyArray_TYPES.
   *@return the corresponding typecode in char.
   */
  char type2char(PyArray_TYPES t_type);
  
  /** 
   *Coverts a single character typecode to its PyArray_TYPES.
   *@param e_type a PyArray_TYPES typecode in char.
   *@return its corresponding PyArray_TYPES.
   */
  PyArray_TYPES char2type(char e_type);

  /**
   *Constructs a string which contains a list of elements extracted from the 
   *input vector.
   *@param vec a vector of any type.
   *@return a string that lists the elements from the input vector.
   */
  template <class T>
  inline std::string vector_str(const std::vector<T>& vec);

  /**
   *Throws an exception if the total size computed from a vector of integer
   *does not match with the expected size.
   *@param dims an integer vector of dimensions.
   *@param n an expected size.
   *@return -----
   */
  inline void check_size_match(std::vector<int> dims, int n);

} //  namespace num_util

#endif	//#if defined(USING_NUMARRAY)

#endif
