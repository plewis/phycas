/* this is a trick to get dcdflib.c compiled as C++
	disutils just uses the suffix to determine which compiler to use
*/
extern "C"
{
#include "thirdparty/dcdflib/src/dcdflib.c"
}
