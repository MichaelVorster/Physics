program main

  implicit none

  integer, parameter :: dp = kind(1.d0)
  integer, parameter :: sp = kind(1.0)

! *******************************************
! **                                       **
! **             USER INPUT                **
! **                                       **
! *******************************************

  integer :: time = 12
  character(len=350) :: file_dir = '/home/cronus/vorster/512_data/bx1.6_cs0.8/initial_turbulence_for_PLUTO/test_convert/'

! *******************************************
! **                                       **
! **            END USER INPUT             **
! **                                       **
! *******************************************

  integer :: nx, ny, nz
  real(dp), dimension(:,:,:), allocatable :: vec_x_dp, vec_y_dp, vec_z_dp, scalar_dp
  real(sp), dimension(:,:,:), allocatable :: vec_x_sp, vec_y_sp, vec_z_sp, scalar_sp
  real(dp), dimension(:), allocatable :: x_dp, y_dp, z_dp
  real(sp), dimension(:), allocatable :: x_sp, y_sp, z_sp
  real(dp) :: dx_dp, dy_dp, dz_dp, time_dp
  real(sp) :: dx_sp, dy_sp, dz_sp, time_sp
  integer :: temp_nx, temp_ny, temp_nz
  
  character(len=350) :: fname
  character(len=10) :: tmp
  integer :: iu = 103


  write (tmp, '(I3,A4)') time, '.BIN'
  fname = trim(file_dir) // '/BB' // adjustl(tmp)
  open(unit=iu, file=trim(fname), form='unformatted', status='old',  action='read')
  read(iu) nx, ny, nz
  close(iu)

  allocate (x_dp(nx))
  allocate (y_dp(ny))
  allocate (z_dp(nz))
  allocate (x_sp(nx))
  allocate (y_sp(ny))
  allocate (z_sp(nz))

  call convert_data('BB', 'magnetic field', 'vector')
  call convert_data('DN', 'density', 'scalar')
  call convert_data('VV', 'velocity', 'vector')

  deallocate (x_dp)
  deallocate (y_dp)
  deallocate (z_dp)
  deallocate (x_sp)
  deallocate (y_sp)
  deallocate (z_sp)


  write(*, "(A)") '* Done!'

contains

  subroutine convert_data(field_prefix, name, type)
    character(len=*) :: field_prefix, name, type

    if (type .eq. 'vector') then
      allocate (vec_x_dp(nx,ny,nz))
      allocate (vec_y_dp(nx,ny,nz))
      allocate (vec_z_dp(nx,ny,nz))
    else
      allocate (scalar_dp(nx,ny,nz))
    endif  

  ! read turbulence data
    write(*, "(3A)") '* Converting ', trim(name), ' to single precision...'

    write (tmp, '(I3,A4)') time, '.BIN'
    fname = trim(file_dir) // '/' // trim(field_prefix) // adjustl(tmp)
    open(unit=iu, file=trim(fname), form='unformatted', status='old',  action='read')
      read(iu) temp_nx, temp_ny, temp_nz
      if (type .eq. 'vector') then
        read(iu) vec_x_dp(:,:,:)
        read(iu) vec_y_dp(:,:,:)
        read(iu) vec_z_dp(:,:,:)
      else
        read(iu) scalar_dp(:,:,:)
      endif
      if (name .ne. 'velocity') then
        read(iu) time_dp, x_dp(:), y_dp(:), z_dp(:), dx_dp, dy_dp, dz_dp
      endif
    close(iu)

  ! convert turbulence field
    if (type .eq. 'vector') then
      allocate (vec_x_sp(nx,ny,nz))
      allocate (vec_y_sp(nx,ny,nz))
      allocate (vec_z_sp(nx,ny,nz))

      vec_x_sp = real(vec_x_dp)
      vec_y_sp = real(vec_y_dp)
      vec_z_sp = real(vec_z_dp)

      deallocate (vec_x_dp)
      deallocate (vec_y_dp)
      deallocate (vec_z_dp)
    else
      allocate (scalar_sp(nx,ny,nz))

      scalar_sp = real(scalar_dp)

      deallocate (scalar_dp)
    endif

  ! convert meta data
    x_sp = real(x_dp)
    y_sp = real(y_dp)
    z_sp = real(z_dp)
    
    time_sp = real(time_dp)
    dx_sp = real(dx_dp)
    dy_sp = real(dy_dp)
    dz_sp = real(dz_dp)

  ! write to disk  
    write(*, "(A)") '--> Writing converted data to disk...'

    fname = trim(file_dir) // '/' // trim(field_prefix) // adjustl(tmp)
    open(unit=iu, file=(trim(fname)), form='unformatted', status='replace', action='write')
      write(iu) temp_nx, temp_ny, temp_nz
      if (type .eq. 'vector') then
        write(iu) vec_x_sp(:,:,:)
        write(iu) vec_y_sp(:,:,:)
        write(iu) vec_z_sp(:,:,:)
      else
        write(iu) scalar_sp(:,:,:)
      endif
      write(iu) time_sp, x_sp(:), y_sp(:), z_sp(:), dx_sp, dy_sp, dz_sp
    close(iu)

  ! realease memory for single precision files
    if (type .eq. 'vector') then
      deallocate (vec_x_sp)
      deallocate (vec_y_sp)
      deallocate (vec_z_sp)
    else
      deallocate (scalar_sp)
    endif

  end subroutine convert_data

end program main
