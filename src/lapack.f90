!> \file
!> \author Chris Bradley
!> \brief This module contains the interface descriptions to the LAPACK routines.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand, the University of Oxford, Oxford, United
!> Kingdom and King's College, London, United Kingdom. Portions created
!> by the University of Auckland, the University of Oxford and King's
!> College, London are Copyright (C) 2007-2010 by the University of
!> Auckland, the University of Oxford and King's College, London.
!> All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!>This module contains the interface descriptions to the LAPACK routines.
MODULE LAPACK

  USE KINDS

  IMPLICIT NONE

  INTERFACE

    SUBROUTINE DGESV(N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      USE KINDS
      INTEGER(INTG), INTENT(IN) :: N
      INTEGER(INTG), INTENT(IN) :: NRHS
      INTEGER(INTG), INTENT(IN) :: LDA
      REAL(DP), INTENT(INOUT) :: A(LDA,*)
      INTEGER(INTG), INTENT(OUT) :: IPIV(*)
      INTEGER(INTG), INTENT(IN) :: LDB
      REAL(DP), INTENT(INOUT) :: B(LDB,*)
      INTEGER(INTG), INTENT(OUT) :: INFO
    END SUBROUTINE DGESV

    ! DGESVD - compute the singular value decomposition (SVD) of a
    !  real M-by-N matrix A, optionally computing the left and/or
    !  right singular vectors
    SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
      USE KINDS
      CHARACTER(1) :: JOBU ! Specifies options for computing all or part of the matrix U (options: A,S,O,N)
      CHARACTER(1) :: JOBVT ! Specifies options for computing all or part of the matrix V**T (options: A,S,O,N)
      INTEGER(INTG), INTENT(IN) :: M ! Number of rows in A
      INTEGER(INTG), INTENT(IN) :: N ! Number of columns in A
      REAL(DP), INTENT(INOUT) :: A(LDA,*) ! The matrix to perform the SVD on
      INTEGER(INTG), INTENT(IN) :: LDA ! Leading dimension of A
      REAL(DP), INTENT(OUT) :: S(MIN(M,N)) ! Singular values of A, sorted S(i) >= S(i+1)
      REAL(DP), INTENT(OUT) :: U(LDU,*) ! If JOBU = 'A', U contains the M-by-M orthogonal matrix U
      INTEGER(INTG), INTENT(IN) :: LDU ! Leading dimension of U
      REAL(DP), INTENT(OUT) :: VT(LDVT,N) ! If JOBVT = 'A', VT contains the N-by-N orthogonal matrix V**T
      INTEGER(INTG), INTENT(IN) :: LDVT ! The leading dimension of the array VT
      REAL(DP), INTENT(INOUT) :: WORK(*) ! On exit, if INFO = 0, WORK(1) returns the optimal LWORK
      INTEGER(INTG), INTENT(IN) :: LWORK ! The dimension of the array WORK
      INTEGER(INTG), INTENT(OUT) :: INFO ! 0 if successful exit; < 0 if INFO = -i (the i-th argument had an illegal value); > 0 if DBDSQR did not converge
    END SUBROUTINE DGESVD

    !DSYEV - compute all eigenvalues and, optionally, eigenvectors of a
    !  real symmetric matrix A.
    ! JOBZ is CHARACTER*1
    !     = 'N':  Compute eigenvalues only;
    !     = 'V':  Compute eigenvalues and eigenvectors.
    ! UPLO is CHARACTER*1
    !     = 'U':  Upper triangle of A is stored;
    !     = 'L':  Lower triangle of A is stored.
    ! N is INTEGER
    !      The order of the matrix A.  N >= 0.
    ! A is DOUBLE PRECISION array, dimension (LDA, N)
    !      On entry, the symmetric matrix A.  If UPLO = 'U', the
    !      leading N-by-N upper triangular part of A contains the
    !      upper triangular part of the matrix A.  If UPLO = 'L',
    !      the leading N-by-N lower triangular part of A contains
    !      the lower triangular part of the matrix A.
    !      On exit, if JOBZ = 'V', then if INFO = 0, A contains the
    !      orthonormal eigenvectors of the matrix A.
    !      If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
    !      or the upper triangle (if UPLO='U') of A, including the
    !      diagonal, is destroyed.
    ! LDA is INTEGER
    !      The leading dimension of the array A.  LDA >= max(1,N).
    ! W is DOUBLE PRECISION array, dimension (N)
    !      If INFO = 0, the eigenvalues in ascending order.
    ! WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
    !      On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
    ! LWORK is INTEGER
    !      The length of the array WORK.  LWORK >= max(1,3*N-1).
    !      For optimal efficiency, LWORK >= (NB+2)*N,
    !      where NB is the blocksize for DSYTRD returned by ILAENV.
    !      If LWORK = -1, then a workspace query is assumed; the routine
    !      only calculates the optimal size of the WORK array, returns
    !      this value as the first entry of the WORK array, and no error
    !      message related to LWORK is issued by XERBLA.
    !
    ! INFO is INTEGER
    !      = 0:  successful exit
    !      < 0:  if INFO = -i, the i-th argument had an illegal value
    !      > 0:  if INFO = i, the algorithm failed to converge; i
    !            off-diagonal elements of an intermediate tridiagonal
    !            form did not converge to zero.
    SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
      USE KINDS
      CHARACTER(LEN=1), INTENT(IN) :: JOBZ
      CHARACTER(LEN=1), INTENT(IN) :: UPLO
      INTEGER(INTG), INTENT(IN) :: N
      REAL(DP), INTENT(INOUT) :: A(LDA,N)
      INTEGER(INTG), INTENT(IN) :: LDA
      REAL(DP), INTENT(OUT) :: W(N)
      REAL(DP), INTENT(OUT) :: WORK(LWORK)
      INTEGER(INTG), INTENT(IN) :: LWORK
      INTEGER(INTG), INTENT(OUT) :: INFO
    END SUBROUTINE DSYEV

  END INTERFACE

END MODULE LAPACK
