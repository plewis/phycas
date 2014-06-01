/* this is a trick to get ipmpar.c compiled as C++
	disutils just uses the suffix to determine which compiler to use
*/
extern "C"
{
#include "thirdparty/dcdflib/src/ipmpar.c"
}
