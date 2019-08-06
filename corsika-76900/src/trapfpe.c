#ifdef __GNUC__

# define _GNU_SOURCE 1
# include <features.h>

# if defined(__GLIBC__) && defined(__GLIBC_MINOR__)
#  if ( ( __GLIBC__ == 2 && __GLIBC_MINOR__ >= 2 ) || __GLIBC__ > 2 )

/* Glibc >= 2.2 specific (but hardware-independent?) version */
/* From: http://gcc.gnu.org/ml/gcc-patches/2000-11/msg01099.html */
/* (see also: http://www.fortran-2000.com/ArnaudRecipes/CompilerTricks.html ) */

#   define _GNU_SOURCE 1
#   include <fenv.h>
static void __attribute__ ((constructor)) trapfpe(void)
{
  /* Enable some exceptions. At startup all exceptions are masked. */
  feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}

#  endif
# else

/* x86-hardware-specific version compatible with older library versions */
/* From: http://gcc.gnu.org/ml/gcc-bugs/1999-11/msg00846.html */

# include <fpu_control.h>
static void __attribute__ ((constructor))
trapfpe ()
{
  fpu_control_t cw = 
     _FPU_DEFAULT &
     ~(_FPU_MASK_IM | _FPU_MASK_ZM | _FPU_MASK_OM);
     _FPU_SETCW(cw);
}

# endif

#else

static void trapfpe()
{
}

#endif

#ifdef TEST_TRAPFPE

/* Build with: gcc -DTEST_TRAPFPE trapfpe.c -lm -o trapfpe */
/* Run with:   ./trapfpe
   If catching a divide by zero works you should get the following output:
      Expecting a floating point exception ...
      Floating point exception
   If not, you may get the following output:
      Expecting a floating point exception ...
      No exception thrown: 1.0/0.0 = inf
*/

#include <stdio.h>

int main()
{
   double a, b, c;
   trapfpe();
   a=1.0;
   b=0.0;
   fprintf(stderr,"Expecting a floating point exception ...\n");
   c=a/b;
   fprintf(stderr,"No exception thrown: 1.0/0.0 = %f\n",c);
   return 0;
}

#endif


