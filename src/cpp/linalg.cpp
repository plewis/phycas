/* this is a trick to get linalg.c compiled as C++
	disutils just uses the suffix to determine which compiler to use
*/
extern "C"
{
#include "linalg.c"
}
