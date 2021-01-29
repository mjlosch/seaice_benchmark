C $Header: /u/gcmpack/MITgcm/verification/offline_exf_seaice/code/SEAICE_SIZE.h,v 1.4 2014/06/25 10:54:27 mlosch Exp $
C $Name:  $

#ifdef ALLOW_SEAICE

CBOP
C    !ROUTINE: SEAICE_SIZE.h
C    !INTERFACE:
C #include SEAICE_SIZE.h

C    !DESCRIPTION:
C Contains seaice array-size definition (number of tracers,categories).

C SItrMaxNum :: number of passive tracers to allocate
C nITD       :: number of seaice categories to allocate
CEOP

C-    Maximum Number of categories
      INTEGER nITD
C--
#ifdef SEAICE_ITD
CToM<<<
C nITD defines number of ice thickness categories,
C i.e. size of additional dimension to AREA, HEFF, HSNOW, etc.
C Bitz et al. (2001, JGR) suggest a minimum of nITD = 5
      PARAMETER (nITD = 5)
C>>>ToM
#else
      PARAMETER (nITD = 1)
#endif

C-    Maximum Number of tracers
      INTEGER SItrMaxNum
      PARAMETER(SItrMaxNum = 3 )

#ifdef ALLOW_AUTODIFF_TAMC
      INTEGER iicekey
      INTEGER nEVPstepMax
      PARAMETER ( nEVPstepMax=180 )
      INTEGER NMAX_TICE
      PARAMETER ( NMAX_TICE=10 )
      INTEGER SOLV_MAX_FIXED
      PARAMETER ( SOLV_MAX_FIXED=500 )
      INTEGER MPSEUDOTIMESTEPS
      PARAMETER (MPSEUDOTIMESTEPS=2)
#endif

#endif /* ALLOW_SEAICE */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
