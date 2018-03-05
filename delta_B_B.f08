program main

  implicit none

  integer, parameter :: dp = kind(1.d0)

! *******************************************
! **                                       **
! **             USER INPUT                **
! **                                       **
! *******************************************

  integer :: time = 12
  character(len=350) :: file_dir = '/home/cronus/vorster/512_data/bx1.6_cs0.8/initial_turbulence_for_PLUTO/'
  logical :: shock_present = .false.
  character(len=1) :: shock_direction = 'x'
  integer :: avg_shock_pos = 266
! 'cut': number of cells on either side of average shock position that is be excluded
  integer :: cut = 0

! *******************************************
! **                                       **
! **            END USER INPUT             **
! **                                       **
! *******************************************

  integer, save :: nx, ny, nz
  real(dp), dimension(:,:,:), allocatable :: magx, magy, magz
  real(dp), dimension(:,:,:), allocatable :: mag_temp
  integer :: temp_nx, temp_ny, temp_nz
  real(dp) :: B_avg_component(3), delta_B_component(3)
  real(dp) :: component_temp
  real(dp) :: avg_B, delta_B

  character(len=350) :: fname
  character(len=10) :: tmp
  integer :: iu = 103

  write (tmp, '(I3,A4)') time, '.BIN'
  fname = trim(file_dir) // '/BB' // adjustl(tmp)
  open(unit=9, file=trim(fname), form='unformatted', status='old',  action='read')
  read(9) nx, ny, nz
  close(9)

  allocate (magx(nx,ny,nz))
  allocate (magy(nx,ny,nz))
  allocate (magz(nx,ny,nz))
  allocate (mag_temp(nx,ny,nz))

  ! read turbulence data
  write (tmp, '(I3,A4)') time, '.BIN'
  fname = trim(file_dir) // '/BB' // adjustl(tmp)
  open(unit=iu, file=trim(fname), form='unformatted', status='old',  action='read')
    read(iu) temp_nx, temp_ny, temp_nz
    read(iu) magx(:,:,:)
    read(iu) magy(:,:,:)
    read(iu) magz(:,:,:)
  close(iu)

  ! calculate average B
  if (shock_present .eqv. .true.) then
    call cut_average(B_avg_component(1), magx)
    call cut_average(B_avg_component(2), magy)
    call cut_average(B_avg_component(3), magz)
  else
    B_avg_component(1) = sum(magx)/(nx*ny*nz)
    B_avg_component(2) = sum(magy)/(nx*ny*nz)
    B_avg_component(3) = sum(magz)/(nx*ny*nz)
  endif
  avg_B = sqrt(sum(B_avg_component**2))

  print*, ''
  print('(A, 3ES14.2)'), 'Average B components:', B_avg_component
  print('(A, ES14.2)'), 'Average B:', avg_B

  if (shock_present .eqv. .true.) then
    call cut_delta_component(delta_B_component(1), B_avg_component(1), magx)
    call cut_delta_component(delta_B_component(2), B_avg_component(2), magy)
    call cut_delta_component(delta_B_component(3), B_avg_component(3), magz)
  else
    delta_B_component(1) = sqrt(sum(magx**2)/(nx*ny*nz) - B_avg_component(1)**2)
    delta_B_component(2) = sqrt(sum(magy**2)/(nx*ny*nz) - B_avg_component(2)**2)
    delta_B_component(3) = sqrt(sum(magz**2)/(nx*ny*nz) - B_avg_component(3)**2)
  endif
  delta_B = sqrt(sum(delta_B_component**2))

  print*, ''
  print('(A, 3ES14.2)'), 'delta_B components:', delta_B_component
  print('(A, ES14.2)'), 'delta_B:', delta_B

  print*, ''
  print('(A, ES14.2)'), 'delta_B/B', delta_B/avg_B
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
  
  deallocate (magx)
  deallocate (magy)
  deallocate (magz)
  deallocate (mag_temp)


contains


subroutine make_cut_x(component, mag)
  real(dp), intent(inout) :: component
  real(dp), intent(inout), dimension(:,:,:) :: mag

  component = ( &
    sum(mag(1 : avg_shock_pos-cut-1, :, :)) + &
    sum(mag(avg_shock_pos+cut+1 : nx, :, :)) &
  )/((nx-(2*cut-1)) * ny * nz)
end subroutine make_cut_x


subroutine make_cut_y(component, mag)
  real(dp), intent(inout) :: component
  real(dp), intent(inout), dimension(:,:,:) :: mag

  component = ( &
    sum(mag(:, 1 : avg_shock_pos-cut-1, :)) + &
    sum(mag(:, avg_shock_pos+cut+1 : ny, :)) &
  )/(nx * (ny-(2*cut-1)) * nz)
end subroutine make_cut_y


subroutine make_cut_z(component, mag)
  real(dp), intent(inout) :: component
  real(dp), intent(inout), dimension(:,:,:) :: mag

  component = ( &
    sum(mag(:, :, 1 : avg_shock_pos-cut-1)) + &
    sum(mag(:, :, avg_shock_pos+cut+1 : nz)) &
  )/(nx * ny * (nz-(2*cut-1)))
end subroutine make_cut_z


subroutine cut_average(component, mag)
  real(dp), intent(inout) :: component
  real(dp), intent(inout), dimension(:,:,:) :: mag

  if (shock_direction .eq. 'x') call make_cut_x(component, mag)
  if (shock_direction .eq. 'y') call make_cut_y(component, mag)
  if (shock_direction .eq. 'z') call make_cut_z(component, mag)    
end subroutine cut_average


subroutine cut_delta_component(component, component_avg, mag)
  real(dp), intent(inout) :: component
  real(dp), intent(inout) :: component_avg
  real(dp), intent(inout), dimension(:,:,:) :: mag

  mag_temp = mag**2
    if (shock_direction .eq. 'x') call make_cut_x(component_temp, mag_temp)
    if (shock_direction .eq. 'y') call make_cut_y(component_temp, mag_temp)
    if (shock_direction .eq. 'z') call make_cut_z(component_temp, mag_temp)    
  component = sqrt(component_temp - component_avg**2)

end subroutine cut_delta_component

end program main
