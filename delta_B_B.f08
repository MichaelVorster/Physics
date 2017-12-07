program main

  implicit none

  integer, parameter :: dp = kind(1.d0)

  integer :: time = 12
  character(len=350) :: file_dir = '/home/cronus/vorster/Synthetic_turb/combined_turbulence/particle_code_format/'

  integer, save :: nx, ny, nz
  real(dp), dimension(:,:,:), allocatable :: magx, magy, magz
  real(dp), dimension(:), allocatable :: temp_x
  real(dp), dimension(:), allocatable :: temp_y
  real(dp), dimension(:), allocatable :: temp_z
  integer :: temp_nx, temp_ny, temp_nz
  real(dp) :: temp_time, temp_dx, temp_dy, temp_dz
  real(dp) :: B_avg_component(3), delta_B_component(3)
  real(dp) :: avg_B, delta_B

  character(len=350) :: fname
  character(len=10) :: tmp
  integer :: iu = 103

  write (tmp, '(I3,A4)') time, '.BIN'
  fname = trim(file_dir) // 'BB' // adjustl(tmp)
  open(unit=9, file=trim(fname), form='unformatted', status='old',  action='read')
  read(9) nx, ny, nz
  close(9)

  allocate (magx(nx,ny,nz))
  allocate (magy(nx,ny,nz))
  allocate (magz(nx,ny,nz))
  allocate (temp_x(nx))
  allocate (temp_y(ny))
  allocate (temp_z(nz))

  ! read original turbulence data
  write (tmp, '(I3,A4)') time, '.BIN'
  fname = trim(file_dir) // 'BB' // adjustl(tmp)
  open(unit=iu, file=trim(fname), form='unformatted', status='old',  action='read')
  read(iu) temp_nx, temp_ny, temp_nz
  read(iu) magx(:,:,:)
  read(iu) magy(:,:,:)
  read(iu) magz(:,:,:)
  read(iu) temp_time, temp_x, temp_y, temp_z, temp_dx, temp_dy, temp_dz
  close(iu)

  ! calculate average B
  B_avg_component(1) = sum(magx)/(nx*ny*nz)
  B_avg_component(2) = sum(magy)/(nx*ny*nz)
  B_avg_component(3) = sum(magz)/(nx*ny*nz)
  avg_B = sqrt(sum(B_avg_component**2))

  print*, ''
  print*, 'Average B components:', B_avg_component
  print*, 'Average B:', avg_B

  delta_B_component(1) = sqrt(sum(magx**2)/(nx*ny*nz) - (sum(magx)/(nx*ny*nz))**2)
  delta_B_component(2) = sqrt(sum(magy**2)/(nx*ny*nz) - (sum(magy)/(nx*ny*nz))**2)
  delta_B_component(3) = sqrt(sum(magz**2)/(nx*ny*nz) - (sum(magz)/(nx*ny*nz))**2)
  delta_B = sqrt(sum(delta_B_component**2))

  print*, ''
  print*, 'delta_B components:', delta_B_component
  print*, 'delta_B:', delta_B

  print*, ''
  print*, 'delta_B/B', delta_B/avg_B

  deallocate (magx)
  deallocate (magy)
  deallocate (magz)
  deallocate (temp_x)
  deallocate (temp_y)
  deallocate (temp_z)

end program main
