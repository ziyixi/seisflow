!///////////////////////////////////////////////////////////////////////////////
subroutine anchor_point_index(iax,iay,iaz)
!-get the index of anchor points in the GLL points of spectral element
! the topology of the nodes is described in UTILS/chunk_notes_scanned/numbering_convention_27_nodes.tif
! currently anchor points are defined on 3x3x3 hex
!
!-output: iax, iay, iaz
  ! add constants
    integer, parameter :: NGNOD = 27,NGLLX=5,NGLLY=5,NGLLZ=5
    integer, parameter :: dp = kind(1.d0)


  ! index of anchor points
  integer, dimension(NGNOD), intent(out) :: iax,iay,iaz

  ! topology of the control points of the surface element
  integer :: ia
  integer, dimension(NGNOD) :: iaddx, iaddy, iaddz

  ! define topology of the control element
  call hex_nodes(iaddx,iaddy,iaddz)

  do ia = 1, NGNOD

    if (iaddx(ia) == 0) then
      iax(ia) = 1
    else if (iaddx(ia) == 1) then
      iax(ia) = (NGLLX+1)/2
    else if (iaddx(ia) == 2) then
      iax(ia) = NGLLX
    else
      stop 'incorrect value of iaddx'
    endif

    if (iaddy(ia) == 0) then
      iay(ia) = 1
    else if (iaddy(ia) == 1) then
      iay(ia) = (NGLLY+1)/2
    else if (iaddy(ia) == 2) then
      iay(ia) = NGLLY
    else
      stop 'incorrect value of iaddy'
    endif

    if (iaddz(ia) == 0) then
      iaz(ia) = 1
    else if (iaddz(ia) == 1) then
      iaz(ia) = (NGLLZ+1)/2
    else if (iaddz(ia) == 2) then
      iaz(ia) = NGLLZ
    else
      stop 'incorrect value of iaddz'
    endif

  end do ! do ia = 1,NGNOD

end subroutine anchor_point_index


!///////////////////////////////////////////////////////////////////////////////
subroutine hex_nodes(iaddx, iaddy, iaddz)
! index of gll points that define the shape of element
! the topology of the nodes is described in UTILS/chunk_notes_scanned/numbering_convention_27_nodes.tif

  ! add constants
    integer, parameter :: NGNOD = 27,NGLLX=5,NGLLY=5,NGLLZ=5
    integer, parameter :: dp = kind(1.d0)

  integer, dimension(NGNOD), intent(out) :: iaddx,iaddy,iaddz

  if (NGNOD /= 27) then
    stop '[ERROR] hex_nodes: elements should have 27 control nodes'
  endif

! corner nodes
  iaddx(1) = 0; iaddy(1) = 0; iaddz(1) = 0
  iaddx(2) = 2; iaddy(2) = 0; iaddz(2) = 0
  iaddx(3) = 2; iaddy(3) = 2; iaddz(3) = 0
  iaddx(4) = 0; iaddy(4) = 2; iaddz(4) = 0
  iaddx(5) = 0; iaddy(5) = 0; iaddz(5) = 2
  iaddx(6) = 2; iaddy(6) = 0; iaddz(6) = 2
  iaddx(7) = 2; iaddy(7) = 2; iaddz(7) = 2
  iaddx(8) = 0; iaddy(8) = 2; iaddz(8) = 2

! midside nodes (nodes located in the middle of an edge)
  iaddx(9) = 1; iaddy(9) = 0; iaddz(9) = 0
  iaddx(10) = 2; iaddy(10) = 1; iaddz(10) = 0
  iaddx(11) = 1; iaddy(11) = 2; iaddz(11) = 0
  iaddx(12) = 0; iaddy(12) = 1; iaddz(12) = 0
  iaddx(13) = 0; iaddy(13) = 0; iaddz(13) = 1
  iaddx(14) = 2; iaddy(14) = 0; iaddz(14) = 1
  iaddx(15) = 2; iaddy(15) = 2; iaddz(15) = 1
  iaddx(16) = 0; iaddy(16) = 2; iaddz(16) = 1
  iaddx(17) = 1; iaddy(17) = 0; iaddz(17) = 2
  iaddx(18) = 2; iaddy(18) = 1; iaddz(18) = 2
  iaddx(19) = 1; iaddy(19) = 2; iaddz(19) = 2
  iaddx(20) = 0; iaddy(20) = 1; iaddz(20) = 2

! side center nodes (nodes located in the middle of a face)
  iaddx(21) = 1; iaddy(21) = 1; iaddz(21) = 0
  iaddx(22) = 1; iaddy(22) = 0; iaddz(22) = 1
  iaddx(23) = 2; iaddy(23) = 1; iaddz(23) = 1
  iaddx(24) = 1; iaddy(24) = 2; iaddz(24) = 1
  iaddx(25) = 0; iaddy(25) = 1; iaddz(25) = 1
  iaddx(26) = 1; iaddy(26) = 1; iaddz(26) = 2

! center node (barycenter of the eight corners)
  iaddx(27) = 1; iaddy(27) = 1; iaddz(27) = 1

end subroutine hex_nodes


!///////////////////////////////////////////////////////////////////////////////
subroutine xyz2cube_bounded(xyz_anchor, xyz, uvw, misloc, flag_inside)
!-mapping a given point in physical space (xyz) to the 
! reference cube (uvw), 
! and also flag whether the point is inside the cube
! if the point lies outside the element, calculate the bounded (xi,eta,gamma)
! inside or on the surface of the reference unit cube.
!
!-inputs:
! (real) xyz_anchor(3, NGNOD): anchor points of the element
! (real) xyz(3): coordinates of the target point
!
!-outputs:
! (real) uvw(3): local coordinates in reference cube
! (real) misloc: location misfit abs(xyz - XYZ(uvw))
! (logical) flag_inside: flag whether the target point locates inside the element

  ! add constants
    integer, parameter :: NGNOD = 27,NGLLX=5,NGLLY=5,NGLLZ=5
    integer, parameter :: dp = kind(1.d0)

  real(dp), intent(in) :: xyz_anchor(3, NGNOD)
  real(dp), intent(in) :: xyz(3)

  real(dp), intent(out) :: uvw(3)
  real(dp), intent(out) :: misloc 
  logical, intent(out) :: flag_inside

  ! parameters
  ! number of iterations used to locate point inside one element 
  integer, parameter :: niter = 5
  
  ! local variables
  integer :: iter
  real(dp), dimension(3) :: xyzi ! iteratively improved xyz
  real(dp), dimension(3,3) :: DuvwDxyz
  real(dp), dimension(3) :: dxyz, duvw

  real(dp), parameter ::  zero = 0.0_dp, & 
                          one = 1.0_dp, &
                          minus_one = -1.0_dp

  ! DEBUG
  !print *, "[DEBUG] xyz2cube_bounded"
  !print *, "[DEBUG] xyz_anchor(3,NGNOD)=", xyz_anchor

  ! initialize 
  uvw = zero
  flag_inside = .true.

  ! iteratively update local coordinate uvw to approach the target xyz
  do iter = 1, niter

    ! predicted xyzi and Jacobian for the current uvw
    call cube2xyz(xyz_anchor, uvw, xyzi, DuvwDxyz)

    ! compute difference
    dxyz = xyz - xyzi

    ! compute increments
    duvw = matmul(DuvwDxyz, dxyz)

    ! update values
    uvw = uvw + duvw

    ! limit inside the cube
    if (any(uvw < minus_one .or. uvw > one)) then 
      where (uvw < minus_one) uvw = minus_one
      where (uvw > one) uvw = one
      ! set is_inside to false based on the last iteration
      if (iter == niter) then
        flag_inside = .false.
      endif
    endif

  enddo ! do iter_loop = 1,NUM_ITER
  
  ! calculate the predicted position 
  call cube2xyz(xyz_anchor, uvw, xyzi, DuvwDxyz)

  ! residual distance from the target point
  misloc = sqrt(sum((xyz-xyzi)**2))

end subroutine xyz2cube_bounded

subroutine cube2xyz(anchor_xyz, uvw, xyz, DuvwDxyz)
!-map from local coordinate (uvw) to physical position (xyz)
! the shape the element is defined by the anchor points (anchor_xyz)
!
!-input
! anchor_xyz(3): anchor points
! uvw(3): local coordinate 
!
!-output
! xyz(3): map uvw to physical space
! DuvwDxyz(3,3): jacobian matrix

      ! add constants
    integer, parameter :: NGNOD = 27,NGLLX=5,NGLLY=5,NGLLZ=5
    integer, parameter :: dp = kind(1.d0)

  real(dp), intent(in) :: anchor_xyz(3, NGNOD)
  real(dp), intent(in) :: uvw(3)

  real(dp), intent(out) :: xyz(3)
  real(dp), intent(out) :: DuvwDxyz(3, 3)

  ! local variables
  real(dp), dimension(3) :: lag1, lag2, lag3 
  real(dp), dimension(3) :: lag1p, lag2p, lag3p
  real(dp), dimension(NGNOD) :: shape3D 
  real(dp), dimension(NGNOD, 3) :: dershape3D
  real(dp) :: jacobian
  real(dp), dimension(3,3) :: DxyzDuvw

  ! DEBUG
  !print *, "[DEBUG] cube2xyz"
  !print *, "[DEBUG] xyz_anchor(3,NGNOD)=", anchor_xyz
  !stop

  if (NGNOD /= 27) then
    stop "[ERROR] cube2xyz: elements should have 27 control nodes"
  endif

  ! lagrange polynomials of order 3 on [-1,1], with collocation points: -1,0,1 
  lag1 = uvw * (uvw - 1.0_dp) / 2.0_dp
  lag2 = 1.0_dp - uvw**2
  lag3 = uvw * (uvw + 1.0_dp) / 2.0_dp
  
  ! derivative of lagrange polynomials
  lag1p = uvw - 0.5_dp
  lag2p = -2.0_dp * uvw
  lag3p = uvw + 0.5_dp
  
  ! construct the shape function
  shape3D = (/ &
       ! corner center
       lag1(1)*lag1(2)*lag1(3), &
       lag3(1)*lag1(2)*lag1(3), &
       lag3(1)*lag3(2)*lag1(3), & 
       lag1(1)*lag3(2)*lag1(3), &
       lag1(1)*lag1(2)*lag3(3), &
       lag3(1)*lag1(2)*lag3(3), &
       lag3(1)*lag3(2)*lag3(3), &
       lag1(1)*lag3(2)*lag3(3), &
       ! edge center
       lag2(1)*lag1(2)*lag1(3), &
       lag3(1)*lag2(2)*lag1(3), &
       lag2(1)*lag3(2)*lag1(3), &
       lag1(1)*lag2(2)*lag1(3), &
       lag1(1)*lag1(2)*lag2(3), &
       lag3(1)*lag1(2)*lag2(3), &
       lag3(1)*lag3(2)*lag2(3), &
       lag1(1)*lag3(2)*lag2(3), & 
       lag2(1)*lag1(2)*lag3(3), &
       lag3(1)*lag2(2)*lag3(3), &
       lag2(1)*lag3(2)*lag3(3), &
       lag1(1)*lag2(2)*lag3(3), &
       ! face center
       lag2(1)*lag2(2)*lag1(3), &
       lag2(1)*lag1(2)*lag2(3), &
       lag3(1)*lag2(2)*lag2(3), &
       lag2(1)*lag3(2)*lag2(3), &
       lag1(1)*lag2(2)*lag2(3), &
       lag2(1)*lag2(2)*lag3(3), &
       ! body center
       lag2(1)*lag2(2)*lag2(3) /)
                            
  ! derivative of the shape function
  ! corner center
  dershape3D( 1,:) = (/ lag1p(1)*lag1(2)*lag1(3), lag1(1)*lag1p(2)*lag1(3), lag1(1)*lag1(2)*lag1p(3) /)
  dershape3D( 2,:) = (/ lag3p(1)*lag1(2)*lag1(3), lag3(1)*lag1p(2)*lag1(3), lag3(1)*lag1(2)*lag1p(3) /)
  dershape3D( 3,:) = (/ lag3p(1)*lag3(2)*lag1(3), lag3(1)*lag3p(2)*lag1(3), lag3(1)*lag3(2)*lag1p(3) /)
  dershape3D( 4,:) = (/ lag1p(1)*lag3(2)*lag1(3), lag1(1)*lag3p(2)*lag1(3), lag1(1)*lag3(2)*lag1p(3) /)
  dershape3D( 5,:) = (/ lag1p(1)*lag1(2)*lag3(3), lag1(1)*lag1p(2)*lag3(3), lag1(1)*lag1(2)*lag3p(3) /)
  dershape3D( 6,:) = (/ lag3p(1)*lag1(2)*lag3(3), lag3(1)*lag1p(2)*lag3(3), lag3(1)*lag1(2)*lag3p(3) /)
  dershape3D( 7,:) = (/ lag3p(1)*lag3(2)*lag3(3), lag3(1)*lag3p(2)*lag3(3), lag3(1)*lag3(2)*lag3p(3) /)
  dershape3D( 8,:) = (/ lag1p(1)*lag3(2)*lag3(3), lag1(1)*lag3p(2)*lag3(3), lag1(1)*lag3(2)*lag3p(3) /)
  ! edge center
  dershape3D( 9,:) = (/ lag2p(1)*lag1(2)*lag1(3), lag2(1)*lag1p(2)*lag1(3), lag2(1)*lag1(2)*lag1p(3) /)
  dershape3D(10,:) = (/ lag3p(1)*lag2(2)*lag1(3), lag3(1)*lag2p(2)*lag1(3), lag3(1)*lag2(2)*lag1p(3) /)
  dershape3D(11,:) = (/ lag2p(1)*lag3(2)*lag1(3), lag2(1)*lag3p(2)*lag1(3), lag2(1)*lag3(2)*lag1p(3) /)
  dershape3D(12,:) = (/ lag1p(1)*lag2(2)*lag1(3), lag1(1)*lag2p(2)*lag1(3), lag1(1)*lag2(2)*lag1p(3) /)
  dershape3D(13,:) = (/ lag1p(1)*lag1(2)*lag2(3), lag1(1)*lag1p(2)*lag2(3), lag1(1)*lag1(2)*lag2p(3) /)
  dershape3D(14,:) = (/ lag3p(1)*lag1(2)*lag2(3), lag3(1)*lag1p(2)*lag2(3), lag3(1)*lag1(2)*lag2p(3) /)
  dershape3D(15,:) = (/ lag3p(1)*lag3(2)*lag2(3), lag3(1)*lag3p(2)*lag2(3), lag3(1)*lag3(2)*lag2p(3) /)
  dershape3D(16,:) = (/ lag1p(1)*lag3(2)*lag2(3), lag1(1)*lag3p(2)*lag2(3), lag1(1)*lag3(2)*lag2p(3) /)
  dershape3D(17,:) = (/ lag2p(1)*lag1(2)*lag3(3), lag2(1)*lag1p(2)*lag3(3), lag2(1)*lag1(2)*lag3p(3) /)
  dershape3D(18,:) = (/ lag3p(1)*lag2(2)*lag3(3), lag3(1)*lag2p(2)*lag3(3), lag3(1)*lag2(2)*lag3p(3) /)
  dershape3D(19,:) = (/ lag2p(1)*lag3(2)*lag3(3), lag2(1)*lag3p(2)*lag3(3), lag2(1)*lag3(2)*lag3p(3) /)
  dershape3D(20,:) = (/ lag1p(1)*lag2(2)*lag3(3), lag1(1)*lag2p(2)*lag3(3), lag1(1)*lag2(2)*lag3p(3) /)
  ! face center
  dershape3D(21,:) = (/ lag2p(1)*lag2(2)*lag1(3), lag2(1)*lag2p(2)*lag1(3), lag2(1)*lag2(2)*lag1p(3) /)
  dershape3D(22,:) = (/ lag2p(1)*lag1(2)*lag2(3), lag2(1)*lag1p(2)*lag2(3), lag2(1)*lag1(2)*lag2p(3) /)
  dershape3D(23,:) = (/ lag3p(1)*lag2(2)*lag2(3), lag3(1)*lag2p(2)*lag2(3), lag3(1)*lag2(2)*lag2p(3) /)
  dershape3D(24,:) = (/ lag2p(1)*lag3(2)*lag2(3), lag2(1)*lag3p(2)*lag2(3), lag2(1)*lag3(2)*lag2p(3) /)
  dershape3D(25,:) = (/ lag1p(1)*lag2(2)*lag2(3), lag1(1)*lag2p(2)*lag2(3), lag1(1)*lag2(2)*lag2p(3) /)
  dershape3D(26,:) = (/ lag2p(1)*lag2(2)*lag3(3), lag2(1)*lag2p(2)*lag3(3), lag2(1)*lag2(2)*lag3p(3) /)
  ! body center
  dershape3D(27,:) = (/ lag2p(1)*lag2(2)*lag2(3), lag2(1)*lag2p(2)*lag2(3), lag2(1)*lag2(2)*lag2p(3) /)

  ! xyz and Dxyz/Duvw
  xyz = matmul(anchor_xyz, shape3D)
  DxyzDuvw = matmul(anchor_xyz, dershape3D)

  ! adjoint matrix: adj(Dxyz/Duvw)
  DuvwDxyz(1,1) =   DxyzDuvw(2,2)*DxyzDuvw(3,3)-DxyzDuvw(3,2)*DxyzDuvw(2,3)
  DuvwDxyz(2,1) = -(DxyzDuvw(2,1)*DxyzDuvw(3,3)-DxyzDuvw(3,1)*DxyzDuvw(2,3))
  DuvwDxyz(3,1) =   DxyzDuvw(2,1)*DxyzDuvw(3,2)-DxyzDuvw(3,1)*DxyzDuvw(2,2)
 
  DuvwDxyz(1,2) = -(DxyzDuvw(1,2)*DxyzDuvw(3,3)-DxyzDuvw(3,2)*DxyzDuvw(1,3))
  DuvwDxyz(2,2) =   DxyzDuvw(1,1)*DxyzDuvw(3,3)-DxyzDuvw(3,1)*DxyzDuvw(1,3)
  DuvwDxyz(3,2) = -(DxyzDuvw(1,1)*DxyzDuvw(3,2)-DxyzDuvw(3,1)*DxyzDuvw(1,2))
  
  DuvwDxyz(1,3) =   DxyzDuvw(1,2)*DxyzDuvw(2,3)-DxyzDuvw(2,2)*DxyzDuvw(1,3)
  DuvwDxyz(2,3) = -(DxyzDuvw(1,1)*DxyzDuvw(2,3)-DxyzDuvw(2,1)*DxyzDuvw(1,3))
  DuvwDxyz(3,3) =   DxyzDuvw(1,1)*DxyzDuvw(2,2)-DxyzDuvw(2,1)*DxyzDuvw(1,2)

  ! jacobian = det(Dxyz/Duvw)
  jacobian =  DxyzDuvw(1,1)*DuvwDxyz(1,1) &
            + DxyzDuvw(1,2)*DuvwDxyz(2,1) &
            + DxyzDuvw(1,3)*DuvwDxyz(3,1)
  if (jacobian<=0.0_dp) then
    print *, "[ERROR] cube2xyz: 3D Jacobian undefined jacobian=", jacobian
    print *, "anchor_xyz(3,NGNOD)=", anchor_xyz
    print *, "uvw(3)=", uvw
    stop
  endif

  ! DEBUG
  !print *, "[DEBUG] jacobian=", jacobian
  !print *, "[DEBUG] anchor_xyz(3,NGNOD)=", anchor_xyz
  !print *, "[DEBUG] uvw(3)=", uvw
  !stop

  ! inverse matrix: Duvw/Dxyz = inv(Dxyz/Duvw) = adj(DxyzDuvw)/det(DxyzDuvw)  
  DuvwDxyz = DuvwDxyz / jacobian

end subroutine cube2xyz

!///////////////////////////////////////////////////////////////////////////////
subroutine lagrange_poly(x, ngll, xgll, lagrange)
!-get lagrange interpolation coefficients: L_i(x)
!
!-input
! x: coordinates of interpolating point
! ngll: number of colocation points
! xgll(ngll): coordinates of colocation points
!
!-output
! lagrange(ngll): interpolation coeff.
!

      ! add constants
    integer, parameter :: NGNOD = 27,NGLLX=5,NGLLY=5,NGLLZ=5
    integer, parameter :: dp = kind(1.d0)

  real(dp), intent(in) :: x
  integer, intent(in) :: ngll
  real(dp), intent(in) :: xgll(ngll)

  real(dp), intent(out) :: lagrange(ngll)

  ! local variables
  integer :: i
  integer, dimension(ngll) :: ind
  real(dp), dimension(ngll) :: xx, yy

  ! lagrange(ngll) 
  ind = (/(i, i=1,ngll)/)
  xx = x - xgll

  do i = 1, ngll
    yy = xgll(i) - xgll
    lagrange(i) = product(xx/yy, mask=(ind/=i))
  end do

end subroutine lagrange_poly
