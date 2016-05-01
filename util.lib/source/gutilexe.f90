!
! Input :
!
!       natom  : integer
!
!       mass   : 1D double precision array mass(natom), UNITS : AMU
!
!       coordxyzinp : 2D double precision array coordxyzinp(3, natom), UNITS : BOHR
!       
!       hessianinp  : 2D double precision array hessiancart(3*natom, 3*natom), UNITS : atomic unuts.
!                     **** NOTE this is NOT mass weighted hessian
!
! Output :
!
!       Dmatrix : 2D double precision array Dmatrix(3*natom, 3*natom)
!
!       Lmatrix : 2D double precision array Lmatrix(3*natom, 3*natom)
!                 **** NOTE only Lmatrix(7:3*natom-6,7:3*natom-6) =\= 0
!                      this is the matrix diagonalize f_int in Eq.7
!                      The eigenvector is columnwise, Lmatrix(7:3*natom-6, k) 
!                      is the k^{th} normal mode
!
!       Dmatrix and Lmatrix are those in Eq.9
!       
!       When you use D*L, don't forget to remormalize it
!
!
! Aug 29, 2014
!  
!


  subroutine GAUSVIB ( natom, mass, coordxyzinp, hessianinp, Dmatrix )
! subroutine GAUSVIB ( natom, mass, coordxyzinp, hessianinp, Dmatrix, Lmatrix)

  implicit none

  integer, parameter :: r8  = selected_real_kind (15,307)
  integer, parameter :: i4  = selected_int_kind (9)
  integer, parameter :: i8  = selected_int_kind (18)

  real    ( kind = r8 ), parameter :: amu     = 1822.88839_r8
  real    ( kind = r8 ), parameter :: autoang = 0.529177249_r8
  real    ( kind = r8 ), parameter :: autoeV  = 27.211384_r8
  real    ( kind = r8 ), parameter :: autocm  = 219474.631_r8

  integer ( kind = i4 ), intent (in) :: natom
  real    ( kind = r8 ), dimension (natom), intent (in) :: mass
  real    ( kind = r8 ), dimension (3, natom), intent (in) :: coordxyzinp
  real    ( kind = r8 ), dimension (3*natom, 3*natom), intent (in) :: hessianinp
  real    ( kind = r8 ), dimension (3*natom, 3*natom), intent (out) :: Dmatrix
  !real    ( kind = r8 ), dimension (3*natom, 3*natom), intent (out) :: Lmatrix

  ! local
  integer ( kind = i4 ) :: ncoord
  integer ( kind = i4 ) :: nmode
  real    ( kind = r8 ), dimension (3, natom) :: coordxyz
  real    ( kind = r8 ), dimension (:),   allocatable :: coord
  real    ( kind = r8 ), dimension (:),   allocatable :: freqcart
  real    ( kind = r8 ), dimension (:),   allocatable :: freqint
  real    ( kind = r8 ), dimension (:),   allocatable :: Inertia
  real    ( kind = r8 ), dimension (:),   allocatable :: masscoord

  real    ( kind = r8 ), dimension (:,:), allocatable :: hessiancart
  real    ( kind = r8 ), dimension (:,:), allocatable :: hessianint
  real    ( kind = r8 ), dimension (:,:), allocatable :: modeAll
  real    ( kind = r8 ), dimension (:,:), allocatable :: modecart
  real    ( kind = r8 ), dimension (:,:), allocatable :: modeint
  real    ( kind = r8 ), dimension (:,:), allocatable :: massmatrix 
  real    ( kind = r8 ), dimension (:,:), allocatable :: momentumInertia
  real    ( kind = r8 ), dimension (:,:), allocatable :: momentumvec

  integer ( kind = i4 ) :: i, j, k, m
  real    ( kind = r8 ) :: summ, sumx, sumy, sumz
  real    ( kind = r8 ) :: PX, PY, PZ

  ncoord = 3 * natom
  nmode  = ncoord - 6

  ! allocate all kinds of arrays
  allocate ( coord(ncoord) )
  allocate ( masscoord(ncoord) )
  allocate ( freqcart(ncoord) )
  allocate ( freqint(nmode) )
  allocate ( Inertia(3) )
  allocate ( hessiancart(ncoord, ncoord))
  allocate ( hessianint(nmode, nmode) )
  allocate ( modecart(ncoord, ncoord) )
  allocate ( modeint(nmode, nmode) )
  allocate ( massmatrix(ncoord,ncoord) )
  allocate ( momentumInertia(3, 3) )
  allocate ( momentumvec(3, 3) )

  hessiancart = hessianinp

  ! mass matrix, in a.u.
  massmatrix = 0.0_r8
  k = 1
  do i = 1, natom
     sumx = 1.0_r8 / sqrt(mass(i)*amu)
     do j = 1, 3
        massmatrix(3*(i-1)+j,3*(i-1)+j) = sumx
        masscoord(k) = sumx
        k = k + 1
     end do
  end do

  ! find the center of mass
  sumx=0.0_r8
  sumy=0.0_r8
  sumz=0.0_r8
  summ=0.0_r8
  do i = 1, natom
     summ = mass(i) + summ
     sumx = sumx + mass(i) * coordxyzinp(1, i)
     sumy = sumy + mass(i) * coordxyzinp(2, i)
     sumz = sumz + mass(i) * coordxyzinp(3, i)
  end do
  
  sumx = sumx / summ
  sumy = sumy / summ
  sumz = sumz / summ

  !Print *, "Center of Mass "
  !Print *, sumx, sumy, sumz
  !Print *, "---------------------"

  ! center of mass cartesian coordinate
  do i = 1, natom
     coordxyz(1, i) = coordxyzinp(1, i) - sumx
     coordxyz(2, i) = coordxyzinp(2, i) - sumy
     coordxyz(3, i) = coordxyzinp(3, i) - sumz
     j = 3 * (i-1) + 1
     coord(j)   = coordxyz(1, i)
     coord(j+1) = coordxyz(2, i)
     coord(j+2) = coordxyz(3, i)
  end do


  ! momentum of inertia, loop over x,y,z (i,j=1,2,3)
  momentumInertia = 0.0_r8
  do i = 1, 3
     j = i + 1
     if (j.eq.4) j=1
     k=j+1
     if (k.eq.4) k=1
     do m = 1, natom
        momentumInertia(i,i) = momentumInertia(i,i) + mass(m) * (coordxyz(j,m)**2 + coordxyz(k,m)**2)
     end do 
  end do

  do i = 1, 2
     do j = i+1, 3
        do m = 1, natom
           momentumInertia(i,j) = momentumInertia(i,j) - mass(m) * coordxyz(i,m) * coordxyz(j,m)
        end do
        momentumInertia(j,i) = momentumInertia(i,j)
     end do
  end do 

  ! diagonal momentum of inertia to obtain the X matrix
  call MATRIXDIAGONAL (3, momentumInertia, Inertia, momentumvec)


  ! prepare the 3 translational eiegenvectors in D matrix
  Dmatrix = 0.0_r8
  do i = 1, 3
     do j = 1, natom
        k = 3 * (j-1) + i
        Dmatrix(k,i) = sqrt(mass(j))
     end do
  end do
  
 
  ! prepare the 3 rotation vector (vector 4~6 in D)
  do m = 1, natom

     summ = sqrt(mass(m))

     PX = coordxyz(1,m) * momentumvec(1,1) + coordxyz(2,m) * momentumvec(1,2) + coordxyz(3,m) * momentumvec(1,3)  
     PY = coordxyz(1,m) * momentumvec(2,1) + coordxyz(2,m) * momentumvec(2,2) + coordxyz(3,m) * momentumvec(2,3)
     PZ = coordxyz(1,m) * momentumvec(3,1) + coordxyz(2,m) * momentumvec(3,2) + coordxyz(3,m) * momentumvec(3,3)

     k = 3 * (m-1) + 1
     Dmatrix(k,   4) = (PY * momentumvec(3,1) - PZ * momentumvec(2,1)) * summ 
     Dmatrix(k+1, 4) = (PY * momentumvec(3,2) - PZ * momentumvec(2,2)) * summ
     Dmatrix(k+2, 4) = (PY * momentumvec(3,3) - PZ * momentumvec(2,3)) * summ

     Dmatrix(k,   5) = (PZ * momentumvec(1,1) - PX * momentumvec(3,1)) * summ
     Dmatrix(k+1, 5) = (PZ * momentumvec(1,2) - PX * momentumvec(3,2)) * summ
     Dmatrix(k+2, 5) = (PZ * momentumvec(1,3) - PX * momentumvec(3,3)) * summ

     Dmatrix(k,   6) = (PX * momentumvec(2,1) - PY * momentumvec(1,1)) * summ
     Dmatrix(k+1, 6) = (PX * momentumvec(2,2) - PY * momentumvec(1,2)) * summ 
     Dmatrix(k+2, 6) = (PX * momentumvec(2,3) - PY * momentumvec(1,3)) * summ

  end do


  ! normalise these 6 translational and rotational vectors
  do i = 1, 6
       call NORMALIZEVECTOR(ncoord, Dmatrix(1:ncoord, i))
  end do 

  ! mass weighted cartesian hessian
  do i = 1, ncoord
     do j = 1, i
        hessiancart(i,j) = hessiancart(i,j) * masscoord(i) * masscoord(j)
        hessiancart(j,i) = hessiancart(i,j)
     end do
  end do
 
  ! diagonal hessian
  call MATRIXDIAGONAL (ncoord, hessiancart, freqcart, modecart)
  
  !print *, " frequency without removing translational and rotational degree of freedom "
  !do i = 7, 20  ! ncoord
  !   print *, sqrt(freqcart(I))*autocm
  !end do

  !j=7
  !Write (*, '(F12.5)') (modecart(i,j),i=1,ncoord)

  !print *, "Eigenvector : "
  !Write (*,'(6F12.5)') ((modecart(i,j),j=7,ncoord),i=1,ncoord)

  !print *, "---------------------------------------------"

  Dmatrix(1:ncoord, 7:ncoord) = modecart(1:ncoord, 7:ncoord)  !???????

  ! Schmidt orthogonalization
  call SCHMIDTORTHOGONAL (ncoord, ncoord, Dmatrix)


  !do i = 1, ncoord-1
  !do j = i+1, ncoord
  !call dotproduct(ncoord, Dmatrix(1:ncoord, i), Dmatrix(1:ncoord, j), summ)
  !print *,i,j, summ
  !end do
  !end do

!  hessiancart = matmul (hessiancart, Dmatrix)
!  hessiancart = matmul (transpose(Dmatrix), hessiancart)

!  hessianint(1:nmode, 1:nmode) = hessiancart(7:ncoord, 7:ncoord) 

  !Write (12,'(E30.16)') ((Dmatrix(i,j),j=1,ncoord),i=1,ncoord)

!  call MATRIX_DIAGONAL (nmode, hessianint, freqint, modeint)

  !print *, " frequency after removing translational and rotational degree of freedom "
  !do i = 1,  nmode
  !   print *, sqrt(freqint(i))*autocm
  !end do
  
  ! create M*D*modeint
!  Lmatrix = 0.0_r8
!  Lmatrix(7:ncoord, 7:ncoord) = modeint(1:nmode, 1:nmode)

  !Write (13,'(E30.16)') ((modecart(i,j),j=1,ncoord),i=1,ncoord)  

  !do i = 1, ncoord
  !   masscoord(i) = masscoord(i) * sqrt(amu)
  !   massmatrix(i,i) = masscoord(i)
  !end do

  !Lmatrix = matmul ( Dmatrix, modecart )
  !Lmatrix = matmul ( massmatrix, Lmatrix )

  !Do i = 1, ncoord
  !   call NORMALIZEVECTOR(ncoord, Lmatrix(1:ncoord,i))    
  !end do  

  !Write (*,'(6F12.5)') ((Lmatrix(i,j),j=7,ncoord),i=1,ncoord)
  deallocate ( coord )
  deallocate ( masscoord )
  deallocate ( freqcart )
  deallocate ( freqint )
  deallocate ( Inertia )
  deallocate ( hessianint )
  deallocate ( modecart )
  deallocate ( modeint )
  deallocate ( massmatrix )
  deallocate ( momentumInertia )
  deallocate ( momentumvec )
  deallocate ( hessiancart )

 end subroutine GAUSVIB

  
