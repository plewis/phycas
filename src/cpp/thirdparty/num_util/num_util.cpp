#if defined(USING_NUMARRAY)

#if defined(_MSC_VER)
#	pragma warning(disable : 4267) // boost's builtin_converters.hpp casts size_t to int rather than unsigned
#endif

#define PY_ARRAY_UNIQUE_SYMBOL PyArrayHandle
#define NO_IMPORT_ARRAY
#include "num_util.h"

namespace { const char* rcsid = "$Id: num_util.cpp,v 1.2 2005/10/26 22:46:53 plewis Exp $"; }

namespace num_util{

typedef KindStringMap::value_type  KindStringMapEntry;
KindStringMapEntry kindStringMapEntries[] =
  {
    KindStringMapEntry(PyArray_CHAR,   "PyArray_CHAR"),
    KindStringMapEntry(PyArray_UBYTE,  "PyArray_UBYTE"),
    KindStringMapEntry(PyArray_SBYTE,  "PyArray_SBYTE"),
    KindStringMapEntry(PyArray_SHORT,  "PyArray_SHORT"),
    KindStringMapEntry(PyArray_INT,    "PyArray_INT"),
    KindStringMapEntry(PyArray_LONG,   "PyArray_LONG"),
    KindStringMapEntry(PyArray_FLOAT,  "PyArray_FLOAT"),
    KindStringMapEntry(PyArray_DOUBLE, "PyArray_DOUBLE"),
    KindStringMapEntry(PyArray_CFLOAT, "PyArray_CFLOAT"),
    KindStringMapEntry(PyArray_CDOUBLE,"PyArray_CDOUBLE"),
    KindStringMapEntry(PyArray_OBJECT, "PyArray_OBJECT"),
    KindStringMapEntry(PyArray_NTYPES, "PyArray_NTYPES"),
    KindStringMapEntry(PyArray_NOTYPE ,"PyArray_NOTYPE")
  };

typedef KindCharMap::value_type  KindCharMapEntry;
KindCharMapEntry kindCharMapEntries[] =
  {
    KindCharMapEntry(PyArray_CHAR,   'c'),
    KindCharMapEntry(PyArray_UBYTE,  'b'),
    KindCharMapEntry(PyArray_SBYTE,  '1'),
    KindCharMapEntry(PyArray_SHORT,  's'),
    KindCharMapEntry(PyArray_INT,    'i'),
    KindCharMapEntry(PyArray_LONG,   'l'),
    KindCharMapEntry(PyArray_FLOAT,  'f'),
    KindCharMapEntry(PyArray_DOUBLE, 'd'),
    KindCharMapEntry(PyArray_CFLOAT, 'F'),
    KindCharMapEntry(PyArray_CDOUBLE,'D'),
    KindCharMapEntry(PyArray_OBJECT, 'O')
  };
  
typedef KindTypeMap::value_type  KindTypeMapEntry;
KindTypeMapEntry kindTypeMapEntries[] =
  {
    KindTypeMapEntry('c',PyArray_CHAR),
    KindTypeMapEntry('b',PyArray_UBYTE),
    KindTypeMapEntry('1',PyArray_SBYTE),
    KindTypeMapEntry('s',PyArray_SHORT),
    KindTypeMapEntry('i',PyArray_INT),
    KindTypeMapEntry('l',PyArray_LONG),
    KindTypeMapEntry('f',PyArray_FLOAT),
    KindTypeMapEntry('d',PyArray_DOUBLE),
    KindTypeMapEntry('F',PyArray_CFLOAT),
    KindTypeMapEntry('D',PyArray_CDOUBLE),
    KindTypeMapEntry('O',PyArray_OBJECT)
  };

  
int numStringEntries = sizeof(kindStringMapEntries)/sizeof(KindStringMapEntry);
int numCharEntries = sizeof(kindCharMapEntries)/sizeof(KindCharMapEntry);
int numTypeEntries = sizeof(kindTypeMapEntries)/sizeof(KindTypeMapEntry);


using namespace boost::python;
  
static KindStringMap kindstrings(kindStringMapEntries,
                                   kindStringMapEntries + numStringEntries);

static KindCharMap kindchars(kindCharMapEntries,
                                   kindCharMapEntries + numCharEntries);

static KindTypeMap kindtypes(kindTypeMapEntries,
                                   kindTypeMapEntries + numTypeEntries);

//Create a numarray referencing Python sequence object
numeric::array makeNum(object x){
  if (!PySequence_Check(x.ptr())){
    PyErr_SetString(PyExc_ValueError, "expected a sequence");
    throw_error_already_set();
  }
  object obj(handle<>
	     (PyArray_ContiguousFromObject(x.ptr(),PyArray_NOTYPE,0,0)));
  check_PyArrayElementType(obj);
  return extract<numeric::array>(obj); 
}

//Create a one-dimensional Numeric array of length n and Numeric type t
numeric::array makeNum(int n, PyArray_TYPES t=PyArray_DOUBLE){
  object obj(handle<>(PyArray_FromDims(1, &n, t)));
  return extract<numeric::array>(obj);
}
  
//Create a Numeric array with dimensions dimens and Numeric type t
numeric::array makeNum(std::vector<int> dimens, 
		       PyArray_TYPES t=PyArray_DOUBLE){
  object obj(handle<>(PyArray_FromDims(dimens.size(), &dimens[0], t)));
  return extract<numeric::array>(obj);
}

numeric::array makeNum(int * data, int n=0){
  object obj(handle<>(PyArray_FromDims(1, &n, PyArray_INT)));
  char *arr_data = ((PyArrayObject*) obj.ptr())->data; 
  memcpy(arr_data, data, sizeof(int) * n); // copies the input data to 
                                           // PyArrayObject->data
  return extract<numeric::array>(obj);
}

numeric::array makeNum(float * data, int n=0){
  object obj(handle<>(PyArray_FromDims(1, &n, PyArray_FLOAT)));
  char *arr_data = ((PyArrayObject*) obj.ptr())->data; 
  memcpy(arr_data, data, sizeof(float) * n);
  return extract<numeric::array>(obj);
}
  
numeric::array makeNum(double * data, int n=0){
  object obj(handle<>(PyArray_FromDims(1, &n, PyArray_DOUBLE)));
  char *arr_data = ((PyArrayObject*) obj.ptr())->data;
  memcpy(arr_data, data, sizeof(double) * n);
  return extract<numeric::array>(obj);
}

numeric::array makeNum(char * data, std::vector<int> dims){
  int total = 1;
  std::vector<int>::iterator iter = dims.begin();
  while(iter != dims.end()){
    total *= *iter;
    ++iter;
  }  
  object obj(handle<>(PyArray_FromDims(dims.size(),&dims[0], PyArray_SBYTE)));
  char *arr_data = ((PyArrayObject*) obj.ptr())->data;
  memcpy(arr_data, data, sizeof(char) * total);
  return extract<numeric::array>(obj);
}

numeric::array makeNum(short * data, std::vector<int> dims){
  int total = 1;
  std::vector<int>::iterator iter = dims.begin();
  while(iter != dims.end()){
    total *= *iter;
    ++iter;
  }  
  object obj(handle<>(PyArray_FromDims(dims.size(),&dims[0], PyArray_SHORT)));
  char *arr_data = ((PyArrayObject*) obj.ptr())->data;
  memcpy(arr_data, data, sizeof(short) * total);
  return extract<numeric::array>(obj);
}

numeric::array makeNum(int * data, std::vector<int> dims){
  int total = 1;
  std::vector<int>::iterator iter = dims.begin();
  while(iter != dims.end()){
    total *= *iter;
    ++iter;
  }  
  object obj(handle<>(PyArray_FromDims(dims.size(),&dims[0], PyArray_INT)));
  char *arr_data = ((PyArrayObject*) obj.ptr())->data;
  memcpy(arr_data, data, sizeof(int) * total);
  return extract<numeric::array>(obj);
}

numeric::array makeNum(float * data, std::vector<int> dims){
  int total = 1;
  std::vector<int>::iterator iter = dims.begin();
  while(iter != dims.end()){
    total *= *iter;
    ++iter;
  }    
  object obj(handle<>(PyArray_FromDims(dims.size(),&dims[0], PyArray_FLOAT)));
  char *arr_data = ((PyArrayObject*) obj.ptr())->data;
  memcpy(arr_data, data, sizeof(float) * total);
  return extract<numeric::array>(obj);
}

numeric::array makeNum(double * data, std::vector<int> dims){
  int total = 1;
  std::vector<int>::iterator iter = dims.begin();
  while(iter != dims.end()){
    total *= *iter;
    ++iter;
  }  
  object obj(handle<>(PyArray_FromDims(dims.size(),&dims[0], PyArray_DOUBLE)));
  char *arr_data = ((PyArrayObject*) obj.ptr())->data;
  memcpy(arr_data, data, sizeof(double) * total);
  return extract<numeric::array>(obj);
}

numeric::array makeNum(const numeric::array& arr){
  //Returns a reference of arr by calling numeric::array copy constructor.
  //The copy constructor increases arr's reference count.
  return numeric::array(arr);
} 

PyArray_TYPES type(numeric::array arr){
  return char2type(arr.typecode());
}

void check_type(boost::python::numeric::array arr, 
		PyArray_TYPES expected_type){
  PyArray_TYPES actual_type = type(arr);
  if (actual_type != expected_type) {
    std::ostringstream stream;
    stream << "expected Numeric type " << kindstrings[expected_type]
	   << ", found Numeric type " << kindstrings[actual_type] << std::ends;
    PyErr_SetString(PyExc_TypeError, stream.str().c_str());
    throw_error_already_set();
  }
  return;
}

//Return the number of dimensions
int rank(numeric::array arr){
  if(!PyArray_Check(arr.ptr())){
    PyErr_SetString(PyExc_ValueError, "expected a PyArrayObject");
    throw_error_already_set();
  }
  return ((PyArrayObject*) arr.ptr())->nd;
}

void check_rank(boost::python::numeric::array arr, int expected_rank){
  int actual_rank = rank(arr);
  if (actual_rank != expected_rank) {
    std::ostringstream stream;
    stream << "expected rank " << expected_rank 
	   << ", found rank " << actual_rank << std::ends;
    PyErr_SetString(PyExc_RuntimeError, stream.str().c_str());
    throw_error_already_set();
  }
  return;
}

int size(numeric::array arr)
{
  if(!PyArray_Check(arr.ptr())){
    PyErr_SetString(PyExc_ValueError, "expected a PyArrayObject");
    throw_error_already_set();
  }
  return PyArray_Size(arr.ptr());
}

void check_size(boost::python::numeric::array arr, int expected_size){
  int actual_size = size(arr);
  if (actual_size != expected_size) {
    std::ostringstream stream;
    stream << "expected size " << expected_size 
	   << ", found size " << actual_size << std::ends;
    PyErr_SetString(PyExc_RuntimeError, stream.str().c_str());
    throw_error_already_set();
  }
  return;
}

std::vector<int> shape(numeric::array arr){
  std::vector<int> out_dims;
  if(!PyArray_Check(arr.ptr())){
    PyErr_SetString(PyExc_ValueError, "expected a PyArrayObject");
    throw_error_already_set();
  }
  int* dims_ptr = ((PyArrayObject*) arr.ptr())->dimensions;
  int the_rank = rank(arr);
  for (int i = 0; i < the_rank; i++){
    out_dims.push_back(*(dims_ptr + i));
  }
  return out_dims;
}

int get_dim(boost::python::numeric::array arr, int dimnum){
  PHYCAS_ASSERT(dimnum >= 0);
  int the_rank=rank(arr);
  if(the_rank < dimnum){
    std::ostringstream stream;
    stream << "Error: asked for length of dimension ";
    stream << dimnum << " but rank of array is " << the_rank << std::ends;
    PyErr_SetString(PyExc_RuntimeError, stream.str().c_str());       
    throw_error_already_set();
  }
  std::vector<int> actual_dims = shape(arr);
  return actual_dims[dimnum];
}

void check_shape(boost::python::numeric::array arr, std::vector<int> expected_dims){
  std::vector<int> actual_dims = shape(arr);
  if (actual_dims != expected_dims) {
    std::ostringstream stream;
    stream << "expected dimensions " << vector_str(expected_dims)
	   << ", found dimensions " << vector_str(actual_dims) << std::ends;
    PyErr_SetString(PyExc_RuntimeError, stream.str().c_str());
    throw_error_already_set();
  }
  return;
}

void check_dim(boost::python::numeric::array arr, int dimnum, int dimsize){
  std::vector<int> actual_dims = shape(arr);
  if(actual_dims[dimnum] != dimsize){
    std::ostringstream stream;
    stream << "Error: expected dimension number ";
    stream << dimnum << " to be length " << dimsize;
    stream << ", but found length " << actual_dims[dimnum]  << std::ends;
    PyErr_SetString(PyExc_RuntimeError, stream.str().c_str());       
    throw_error_already_set();
  }
  return;
}

bool iscontiguous(numeric::array arr)
{
  //  return arr.iscontiguous();
  return PyArray_ISCONTIGUOUS((PyArrayObject*)arr.ptr());
}

void check_contiguous(numeric::array arr)
{
  if (!iscontiguous(arr)) {
    PyErr_SetString(PyExc_RuntimeError, "expected a contiguous array");
    throw_error_already_set();
  }
  return;
}

char* data(numeric::array arr){
  if(!PyArray_Check(arr.ptr())){
    PyErr_SetString(PyExc_ValueError, "expected a PyArrayObject");
    throw_error_already_set();
  }
  return ((PyArrayObject*) arr.ptr())->data;
}

//Copy data into the array
void copy_data(boost::python::numeric::array arr, char* new_data){
  char* arr_data = data(arr);
  int nbytes = PyArray_NBYTES((PyArrayObject*)arr.ptr());
  for (int i = 0; i < nbytes; i++) {
    arr_data[i] = new_data[i];
  }
  return;
} 

//Return a clone of this array
numeric::array clone(numeric::array arr){
  object obj(handle<>(PyArray_Copy((PyArrayObject*)arr.ptr())));
  return makeNum(obj);
}

  
//Return a clone of this array with a new type
numeric::array astype(boost::python::numeric::array arr, PyArray_TYPES t){
  return (numeric::array) arr.astype(type2char(t));
}

bool spacesaver(numeric::array arr){
  if(!PyArray_Check(arr.ptr())){
    PyErr_SetString(PyExc_ValueError, "expected a PyArrayObject");
    throw_error_already_set();
  }
  return PyArray_ISSPACESAVER(arr.ptr());
}

void savespace(numeric::array arr, bool set_savespace=true){
  if(!PyArray_Check(arr.ptr())){
    PyErr_SetString(PyExc_ValueError, "expected a PyArrayObject");
    throw_error_already_set();
  }
  int flags = ((PyArrayObject*) arr.ptr())->flags;
  if (set_savespace) {
    flags |= SAVESPACE;
  } 
  else {
    flags &= ~SAVESPACE;
  }
  ((PyArrayObject*) arr.ptr())->flags = flags;
  return;
}

std::vector<int> strides(numeric::array arr){
  std::vector<int> out_strides;
  if(!PyArray_Check(arr.ptr())){
    PyErr_SetString(PyExc_ValueError, "expected a PyArrayObject");
    throw_error_already_set();
  }
  int* strides_ptr = ((PyArrayObject*) arr.ptr())->strides;
  int the_rank = rank(arr);
  for (int i = 0; i < the_rank; i++){
    out_strides.push_back(*(strides_ptr + i));
  }
  return out_strides;
}

int refcount(numeric::array arr){
  return (arr.ptr())->ob_refcnt;
}

void check_PyArrayElementType(object newo){
  PyArrayObject* PyArrayPtr = (PyArrayObject*) newo.ptr();
  PyArray_TYPES theType=(PyArray_TYPES) PyArrayPtr->descr->type_num;
  if(theType == PyArray_OBJECT){
      std::ostringstream stream;
      stream << "array elments have been cast to PyArray_OBJECT, "
             << "numhandle can only accept arrays with numerical elements" 
	     << std::ends;
      PyErr_SetString(PyExc_TypeError, stream.str().c_str());
      throw_error_already_set();
  }
  return;
}

std::string type2string(PyArray_TYPES t_type){
  return kindstrings[t_type];
}

char type2char(PyArray_TYPES t_type){
  return kindchars[t_type];
}

PyArray_TYPES char2type(char e_type){
  return kindtypes[e_type];
}

template <class T>
inline std::string vector_str(const std::vector<T>& vec)
{
  std::ostringstream stream;
  stream << "(" << vec[0];

  for(std::size_t i = 1; i < vec.size(); i++){
    stream << ", " << vec[i];
  }
  stream << ")";
  return stream.str();
}

inline void check_size_match(std::vector<int> dims, int n)
{
  int total = 1;
  std::vector<int>::iterator iter = dims.begin();
  
  while(iter != dims.end()){
    total *= *iter;
    ++iter;
  }
  if (total != n) {
    std::ostringstream stream;
    stream << "expected array size " << n
           << ", dimensions give array size " << total << std::ends;
    PyErr_SetString(PyExc_TypeError, stream.str().c_str());
    throw_error_already_set();
  }
  return;
}

} //namespace numarray_bpl

#endif	//#if defined(USING_NUMARRAY)

