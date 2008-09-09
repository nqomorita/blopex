/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* @@@ BLOPEX (version 1.1)                                              */
/* @@@ Copyright 2008 Merico Argentati, Andrew Knyazev,                  */
/* @@@ Ilya Lashuk, Evgueni Ovtchinnikov, and Don McCuan                 */
/* @@@ LGPL Version 3 or above.  See www.gnu.org.                        */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* This code was developed by Merico Argentati, Andrew Knyazev, Ilya Lashuk and Evgueni Ovtchinnikov */

#ifndef HYPRE_MATVEC_FUNCTIONS
#define HYPRE_MATVEC_FUNCTIONS

typedef struct
{
  void*  (*MatvecCreate)  ( void *A, void *x );
  int    (*Matvec)        ( void *matvec_data, double alpha, void *A,
                             void *x, double beta, void *y );
  int    (*MatvecDestroy) ( void *matvec_data );

  void*  (*MatMultiVecCreate)  ( void *A, void *x );
  int    (*MatMultiVec)        ( void *data, double alpha, void *A,
				 void *x, double beta, void *y );
  int    (*MatMultiVecDestroy)  ( void *data );

} HYPRE_MatvecFunctions;

#endif
