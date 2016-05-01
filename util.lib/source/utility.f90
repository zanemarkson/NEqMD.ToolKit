!
! vector normalization
!
 subroutine NORMALIZEVECTOR (n, a)

  implicit none

  integer, parameter :: r8  = selected_real_kind (15,307)
  integer, parameter :: i4  = selected_int_kind (9)
  integer, parameter :: i8  = selected_int_kind (18)


  integer ( kind = i4 ), intent (in) :: n
  real    ( kind = r8 ), dimension (n), intent (inout) :: a

  real    ( kind = r8 ) :: summ
  integer ( kind = i4 ) :: i

  summ = 0.0_r8

  do i = 1, n
     summ = summ + a(i)**2
  end do

  !print *, sqrt(summ), summ, 1.0/summ

  a = a / sqrt(summ)

 end subroutine NORMALIZEVECTOR
!
!----------------------------------------------
!
 subroutine MATRIXDIAGONAL(n, a, eigenvalue, eigenvector)

  implicit none

  integer, parameter :: r8  = selected_real_kind (15,307)
  integer, parameter :: i4  = selected_int_kind (9)
  integer, parameter :: i8  = selected_int_kind (18)

  integer ( kind = i4 ), parameter :: nmax = 5000

  integer ( kind = i4 ), intent (in)    :: n
  real    ( kind = r8 ), dimension (n, n), intent (in)    :: a
  real    ( kind = r8 ), dimension (n),    intent (inout) :: eigenvalue
  real    ( kind = r8 ), dimension (n, n), intent (inout) :: eigenvector


  ! local
  real    ( kind = r8 ), dimension ( 3*nmax ) :: work
  real    ( kind = r8 ), dimension (nmax)     :: w
  integer ( kind = i4 ), dimension (nmax)     :: ipiv
  integer ( kind = i4 ) :: lwork, lda, info

  integer ( kind = i4 ) :: i,j

  if (n.gt.nmax) then
     print *, "increase nmax in matrix_matrix_diaginal"
     stop
  end if

  !Write (*,*) "n=", n
  !Write (*, '(3F20.10)') ((a(i,j),j=1,3),i=1,3)
  !stop

  lda = n
  lwork = 3 * n
  eigenvector = a
  call dsyev ('V', 'U', n, eigenvector, lda, W, work, lwork, info)
  if (info.ne.0) then
     print *, " error in diag "
     stop
  end if

  eigenvalue(1:n) = w(1:n)

 end subroutine MATRIXDIAGONAL
!
!------------------------------------------------------------------
!
 subroutine DOTPRODUCT(n, v, u, uv)

  implicit none

  integer, parameter :: r8  = selected_real_kind (15,307)
  integer, parameter :: i4  = selected_int_kind (9)
  integer, parameter :: i8  = selected_int_kind (18)

  integer ( kind = i4 ), intent (in) :: n

  real    ( kind = r8 ), dimension (n), intent (in) :: u
  real    ( kind = r8 ), dimension (n), intent (in) :: v

  real    ( kind = r8 ), intent (out) :: uv

  integer ( kind = i4 ) :: i

  uv = 0.0_r8

  do i = 1, n
     uv = uv + u(i) * v(i)     
  end do

 end subroutine DOTPRODUCT

!
!------------------------------------------------------------------
!
 subroutine  SCHMIDTORTHOGONAL(nvec, nlength, A )

  implicit none

  integer, parameter :: r8  = selected_real_kind (15,307)
  integer, parameter :: i4  = selected_int_kind (9)
  integer, parameter :: i8  = selected_int_kind (18)
  
  integer ( kind = i4 ), intent (in) :: nvec
  integer ( kind = i4 ), intent (in) :: nlength

  real    ( kind = r8 ), dimension (nlength, nvec), intent (inout) :: A

  integer ( kind = i4 ) :: i
  integer ( kind = i4 ) :: j
  integer ( kind = i4 ) :: k
  real    ( kind = r8 ) :: norm
  real    ( kind = r8 ) :: coeff

  real    ( kind = r8 ), dimension (nlength) :: x

  if (nvec .eq. 1) return

  do i = 2, nvec

     x(1:nlength) = A(1:nlength, i)

     do j = 1, i - 1

        call dotproduct(nlength, a(1:nlength,j), x, coeff)
        call dotproduct(nlength, a(1:nlength,j), a(1:nlength,j), norm)

        do k = 1, nlength
           a(k,i) = a(k,i) - coeff * a(k,j) / norm
        end do

     end do

  end do

 end subroutine SCHMIDTORTHOGONAL 














