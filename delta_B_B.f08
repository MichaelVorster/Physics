program main

  implicit none

  integer, parameter :: dp = kind(1.d0)

! *******************************************
! **                                       **
! **             USER INPUT                **
! **                                       **
! *******************************************

  integer :: time = 14
  character(len=350) :: file_dir = '/home/cronus/vorster/Dust_grains/Kirit_data/' 
  logical :: shock_present = .false.
  character(len=1) :: shock_direction = 'y'
  integer :: avg_shock_pos = 252
! 'cut': number of cells on either side of average shock position that is be excluded
  integer :: cut = 10

! *******************************************
! **                                       **
! **            END USER INPUT             **
! **                                       **
! *******************************************

  integer, save :: nx, ny, nz
  real(dp), dimension(:,:,:), allocatable :: Bx, By, Bz
  real(dp), dimension(:,:,:), allocatable :: B_temp
  real(dp), dimension(:,:,:), allocatable :: B_mag
  integer :: temp_nx, temp_ny, temp_nz
  real(dp) :: B_avg_component(3), delta_B_component(3)
  real(dp) :: component_temp
  real(dp) :: avg_B, delta_B
  real(dp) :: mu_dispersion
  real(dp) :: temp_val

  character(len=350) :: fname
  character(len=10) :: tmp
  integer :: iu = 103

  write (tmp, '(I3,A4)') time, '.BIN'
  fname = trim(file_dir) // '/BB' // adjustl(tmp)
  open(unit=9, file=trim(fname), form='unformatted', status='old',  action='read')
  read(9) nx, ny, nz
  close(9)

  allocate (Bx(nx,ny,nz))
  allocate (By(nx,ny,nz))
  allocate (Bz(nx,ny,nz))
  allocate (B_temp(nx,ny,nz))
  allocate (B_mag(nx,ny,nz))

  ! read turbulence data
  write (tmp, '(I3,A4)') time, '.BIN'
  fname = trim(file_dir) // '/BB' // adjustl(tmp)
  open(unit=iu, file=trim(fname), form='unformatted', status='old',  action='read')
    read(iu) temp_nx, temp_ny, temp_nz
    read(iu) Bx(:,:,:)
    read(iu) By(:,:,:)
    read(iu) Bz(:,:,:)
  close(iu)

  ! calculate average B
  if (shock_present .eqv. .true.) then
    call cut_average(B_avg_component(1), Bx)
    call cut_average(B_avg_component(2), By)
    call cut_average(B_avg_component(3), Bz)
  else
    B_avg_component(1) = sum(Bx)/(nx*ny*nz)
    B_avg_component(2) = sum(By)/(nx*ny*nz)
    B_avg_component(3) = sum(Bz)/(nx*ny*nz)
  endif
  avg_B = sqrt(sum(B_avg_component**2))

  ! calculate magnitude of B at every point
  B_mag = sqrt(Bx**2 + By**2 + Bz**2)

  ! calculate (delta v_para)/v_perp  : Equation C6 in Hoang, Lazarain, Schlickeiser (2012), ApJ 747:54
  temp_val = sum((B_mag - avg_B)**2)/(nx*ny*nz)
  mu_dispersion = temp_val**(1./4.)/avg_B**(1./2.)

  print*, ''
  print('(A, 3ES14.2)'), 'Average B components:', B_avg_component
  print('(A, ES14.2)'), 'Average B:', avg_B

  if (shock_present .eqv. .true.) then
    call cut_delta_component(delta_B_component(1), B_avg_component(1), Bx)
    call cut_delta_component(delta_B_component(2), B_avg_component(2), By)
    call cut_delta_component(delta_B_component(3), B_avg_component(3), Bz)
  else
    delta_B_component(1) = sqrt(sum(Bx**2)/(nx*ny*nz) - B_avg_component(1)**2)
    delta_B_component(2) = sqrt(sum(By**2)/(nx*ny*nz) - B_avg_component(2)**2)
    delta_B_component(3) = sqrt(sum(Bz**2)/(nx*ny*nz) - B_avg_component(3)**2)
  endif
  delta_B = sqrt(sum(delta_B_component**2))

  print*, ''
  print('(A, 3ES14.2)'), 'delta_B components:', delta_B_component
  print('(A, ES14.2)'), 'delta_B:', delta_B

  print*, ''
  print('(A, ES14.2)'), 'delta_B/B', delta_B/avg_B
  print*, ''

  print*, ''
  print('(A, ES14.2)'), 'mu dispersion', mu_dispersion
  print*, ''

  fname = trim(file_dir) // '/delta_B_B.txt'
  open(unit=iu, file=trim(fname), form='formatted', status='replace',  action='write')
    write(iu,'(A, 3ES14.2)') 'Average B components:', B_avg_component
    write(iu,'(A, ES14.2)') 'Average B:', avg_B
    write(iu, '(A)') ''
    write(iu, '(A, 3ES14.2)') 'delta_B components:', delta_B_component
    write(iu,'(A, ES14.2)') 'delta_B:', delta_B
    write(iu, '(A)') ''
    write(iu,'(A, ES14.2)') 'delta_B/B', delta_B/avg_B
  close(iu)
  
  deallocate (Bx)
  deallocate (By)
  deallocate (Bz)
  deallocate (B_temp)
  deallocate (B_mag)


contains


subroutine make_cut_x(component, B)
  real(dp), intent(inout) :: component
  real(dp), intent(inout), dimension(:,:,:) :: B

  component = ( &
    sum(B(1 : avg_shock_pos-cut-1, :, :)) + &
    sum(B(avg_shock_pos+cut+1 : nx, :, :)) &
  )/((nx-(2*cut-1)) * ny * nz)
end subroutine make_cut_x


subroutine make_cut_y(component, B)
  real(dp), intent(inout) :: component
  real(dp), intent(inout), dimension(:,:,:) :: B

  component = ( &
    sum(B(:, 1 : avg_shock_pos-cut-1, :)) + &
    sum(B(:, avg_shock_pos+cut+1 : ny, :)) &
  )/(nx * (ny-(2*cut-1)) * nz)
end subroutine make_cut_y


subroutine make_cut_z(component, B)
  real(dp), intent(inout) :: component
  real(dp), intent(inout), dimension(:,:,:) :: B

  component = ( &
    sum(B(:, :, 1 : avg_shock_pos-cut-1)) + &
    sum(B(:, :, avg_shock_pos+cut+1 : nz)) &
  )/(nx * ny * (nz-(2*cut-1)))
end subroutine make_cut_z


subroutine cut_average(component, B)
  real(dp), intent(inout) :: component
  real(dp), intent(inout), dimension(:,:,:) :: B

  if (shock_direction .eq. 'x') call make_cut_x(component, B)
  if (shock_direction .eq. 'y') call make_cut_y(component, B)
  if (shock_direction .eq. 'z') call make_cut_z(component, B)    
end subroutine cut_average


subroutine cut_delta_component(component, component_avg, B)
  real(dp), intent(inout) :: component
  real(dp), intent(inout) :: component_avg
  real(dp), intent(inout), dimension(:,:,:) :: B

  B_temp = B**2
    if (shock_direction .eq. 'x') call make_cut_x(component_temp, B_temp)
    if (shock_direction .eq. 'y') call make_cut_y(component_temp, B_temp)
    if (shock_direction .eq. 'z') call make_cut_z(component_temp, B_temp)    
  component = sqrt(component_temp - component_avg**2)

end subroutine cut_delta_component

end program main
