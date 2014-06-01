/*----------------------------------------------------------------------------------------------------------------------
|
|	@file	dls_core.h
|
|	@brief	Dave Swofford's standard definitions.
|
|	This header contains the definitions of several types, preprocessor macros and enum's that are used throughout
|	my programs. 
*/
#ifndef DLS_CORE_H_INCLUDED
#define DLS_CORE_H_INCLUDED

#define PRIVATE static
#define PUBLIC
#define DLS_PRIVATE				/* not static, but should only be called by other "dls" routines */
#define FOREVER ;;				/* makes "for (;;)" more readable */
#define NULLPTR (void *)NULL	/* ensure that "NULL" passed as a stdarg argument will have the right size */

enum { RC_Canceled = -1, RC_OK = 0, RC_Error = 1, RC_User_Abort = 2 };	/* standard function return codes */
enum { Prompt_Ready, Prompt_Pending };									/* prompt status */

#ifndef FALSE
#	define FALSE	0U
#	define TRUE		1U
#endif
#ifndef OFF
#	define OFF		0U
#	define ON		1U
#endif
#ifndef NO
#	define NO		0U
#	define YES		1U
#endif

/*	Macros for global variable definition/declaration:
*/
#ifdef MAIN
#	define EXTERN
#	define INIT(x) = x
#else
#	define EXTERN extern
#	define INIT(x)
#endif

/* MIN and MAX macros:
*/
#ifndef MIN
#	define MIN(x, y)			((x) <= (y)	? (x) : (y))
#endif
#ifndef MAX
#	define MAX(x, y)			((x) >= (y)	? (x) : (y))
#endif

/*	STRUCT_ZERO is a convenience macro used to zero-out a structure at the time it is declared.  It works because
	initializing only the first member of a structure causes all other members to be initialized to zero.  A possible
	drawback is that gcc compiler warnings will be issued if -Wmissing-field-initializers is requested (either
	directly or through -W).  However, I find the convenience of the macro worth giving up that particular warning
	(especially since I use incomplete structure initialization all the time anyway).
*/
#ifdef __cplusplus
#	define STRUCT_ZERO {}
#else
#	define STRUCT_ZERO {0}
#endif

/*	Declare my own "Boolean" types.
|
|	Note: In PAUP* and some of my other programs, sizeof(Bool) MUST be same size as sizeof(int), e.g., for
|	"*(int *)p = x" constructs where the code doesn't know the precise type of the pointer p.  Thus, we don't
|	use enum or C99/nonstandard types such as 'bool' (C++) or Boolean (Macintosh).  We also define a 1-byte
|	Boolean type for situations where we want to minimize the memory footprint of Boolean variables.
*/
typedef unsigned int	Bool;
typedef unsigned char	ShortBool;

/* A few other commonly used structs:
*/
typedef struct fpoint
	{
	float x, y;
	}
	FPoint;

typedef struct intpair
	{
	int i, j;
	}
	IntPair;

typedef struct doublepair
	{
	double a, b;
	}
	DoublePair;
	
#if defined(__GNUC__)
#	define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#endif

/*	Define an UNUSED macro that allows easy suppression of warnings for unused function parameters (e.g., when they
|	must be included in order to satisfy the signature for a callback function).
*/
#if defined(__cplusplus) || defined(__MWERKS__)
#	define UNUSED(x)			/* (omitting name works even for C with Metrowerks) */
#elif defined(__INTEL_COMPILER) || (defined(__GNUC__) && (GCC_VERSION >= 20700))
#	define UNUSED(x) x __attribute__((__unused__))
#else
#	define UNUSED(x) x	/* will have to put up with the warnings, or suppress them using a CC flag */
#endif

/* This macro can be used when it is inconvenient to wrap an unused parameter in the parameter list */
#define SUPPRESS_UNUSED_WARNING(x)	(void)x;

/* Use the WARN_UNUSED_RESULT attribute to generate warnings when a function result is ignored and it's nonsensical
   to do so (e.g., the function has no useful side effects, or it returns allocated memory that needs to be freed). */
#if defined(__INTEL_COMPILER) || (defined(__GNUC__) && (GCC_VERSION >= 30300))
#	define WARN_UNUSED_RESULT __attribute__((__warn_unused_result__))
#else
#	define WARN_UNUSED_RESULT
#endif

#endif	/* DLS_CORE_H_INCLUDED */
