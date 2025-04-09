!> @file mod_wrf_file.f90
!! @brief This module handles WRF file operations, including reading, writing, and defining variables.
!!
!! This module provides utilities for working with WRF NetCDF files, including
!! reading parameters, creating output files, and managing dimensions and variables.
!!
!! @note Ensure that this module is properly included in the relevant parts
!! of the project to access its functionality.

module mod_wrf_file
  use mod_const
  use netcdf
  implicit none

  integer :: n_terms = 7
  integer,parameter :: termV=1,& ! names of omega and height tendency terms
       termT=2,&
       termF=3,&
       termQ=4,&
       termA=5,&
       termVKhi=6,&
       termTKhi=7
  character ( 6 ), dimension ( 7 ), parameter :: &
       omega_term_names &
       = [ 'ome_v ', 'ome_t ', 'ome_f ', 'ome_q ', 'ome_a ', &
       'ome_vs', 'ome_ts' ], &
       htend_term_names &
       = [ 'htv   ', 'htt   ', 'htf   ', 'htq   ', 'hta   ', &
       'htvKhi', 'httKhi' ]
  character ( 8 ), dimension ( 2 ), parameter :: &
       QG_omega_term_names &
       = [ 'ome_v_QG', 'ome_t_QG' ]
  character ( 5 ), parameter :: ome_b_name='ome_b'
  character ( 3 ), parameter :: htend_b_name='htb'
  character ( 9 ), parameter :: ztend_name='ztend_WRF'
  character ( 7 ), parameter :: ome_name='ome_WRF'
  character ( 11 ), dimension ( 4 ), parameter :: rname = &
       [ 'west_east  ', 'south_north', 'vlevs      ', &
       'time       ' ]
  character ( 37 ), dimension ( 2 ), parameter :: QG_omega_long_names = &
       [ 'QG omega due to vorticity advection  ', &
       'QG omega due to temperature advection' ]
  character ( 49 ), dimension ( 7 ), parameter :: omega_long_names = &
       [ 'omega due to vorticity advection                 ', &
       'omega due to temperature advection               ', &
       'omega due to friction                            ', &
       'omega due to diabatic heating                    ', &
       'omega due to ageostrophic vorticity tendency     ', &
       'omega due to vorticity advection by irrot. wind  ',&
       'omega due to temperature advection by irrot. wind' ]
  character ( 59 ), dimension ( 7 ), parameter :: htend_long_names = &
       [ 'height tendency due to vorticity advection                 ', &
       'height tendency due to temperature advection               ', &
       'height tendency due to friction                            ', &
       'height tendency due to diabatic heating                    ', &
       'height tendency due to ageostrophic vorticity tendency     ', &
       'height tendency due to vorticity advection by irrot. wind  ',&
       'height tendency due to temperature advection by irrot. wind' ]
  character ( 49 ), parameter :: ome_b_long_name = &
       'omega due to boundary conditions                 '
  character ( 59 ), parameter :: htend_b_long_name = &
       'height tendency due to boundary conditions                 '

  type wrf_file
     integer :: ncid, dims ( 4 )
     integer :: wrf_cu_phys
     !     real :: dx, dy
     real, dimension ( : ), allocatable :: times, pressure_levels, corpar
     integer, dimension ( : ), allocatable :: xdim,ydim
     integer, dimension ( : ), allocatable :: nlon, nlat, nlev
     real, dimension ( : ), allocatable :: dx, dy, dlev
  end type wrf_file

contains

  !***************************************************************************
  !> @brief Reads simulation parameters from a namelist.
  !!
  !! This subroutine reads simulation parameters from a namelist and populates
  !! the `parameters` type structure.
  !!
  !! @param[out] param Parameters structure to populate.
  !***************************************************************************
  subroutine read_parameters(param)
    type (parameters), intent(out) :: param
    character*140 :: infile, outfile
    character :: mode
    real :: alfa, toler
    integer :: time_1, time_n, ny1, ny2
    logical :: calc_div, calc_omegas, forc, ellipticity_correction, debug

    namelist/PARAMETER_NAMELIST/ infile,outfile,alfa, &
         toler,ny1,ny2,time_1,time_n,&
         mode,calc_omegas,calc_div,debug,forc, &
         ellipticity_correction

    ellipticity_correction=.true.

    read(*,nml=PARAMETER_NAMELIST)
    param % infile = infile
    param % outfile = outfile
    param % alfa = alfa
    param % toler = toler
    param % ny1 = ny1
    param % ny2 = ny2
    param % time_1 = time_1
    param % time_n = time_n
    param % mode = mode
    param % calc_omegas = calc_omegas
    param % calc_div = calc_div
    param % forc = forc
    param % ellipticity_correction = ellipticity_correction

  end subroutine read_parameters

  !***************************************************************************
  !> @brief Creates a new output WRF file.
  !!
  !! This function creates a new NetCDF output file based on the input WRF file,
  !! copying dimensions and defining variables.
  !!
  !! @param[in]  fname      Name of the output file to create.
  !! @param[in]  wrf_infile Input WRF file structure.
  !! @param[in]  mode       Mode of operation ('Q' for QG omega, others for generalized omega).
  !! @param[in]  calc_b     Logical flag to include boundary condition variables.
  !! @param[in]  forc       Logical flag to include forcing variables.
  !! @return     f          Output WRF file structure.
  !***************************************************************************
  function create_out_file ( fname, wrf_infile, mode, calc_b, forc ) result ( f )
    character ( * ),   intent ( in ) :: fname ! file name
    type ( wrf_file ), intent ( in ) :: wrf_infile ! wrf inputfile
    logical,           intent ( in ) :: calc_b, forc
    type ( wrf_file ) :: f
    integer :: dimids(4),i,status,varid,varids(4)
    character :: mode

    print*,"Creating a new output file!"
    ! Create a new netcdf file
    call check( nf90_create ( fname, NF90_CLOBBER, f % ncid ) )

    print*,"Copying dimension information from WRF input file"
    ! Copy dimension information from wrf input file
    f % dims = wrf_infile % dims
    f % pressure_levels = wrf_infile % pressure_levels

    print*,"Creating dimensions and their variables in the new file"
    ! Create dimensions and their variables to the new file
    do i=1,3
       !print*, rname(i), dimids(i), f%dims(i), remove_underscores(rname(i)), f%ncid
       call def_dim(f%ncid, rname(i), dimids(i), f%dims(i), .FALSE.)
    end do
    call def_dim(f%ncid, rname(4), dimids(4), f%dims(4), .TRUE.)

    print*,"Adding axis attributes to the dimensions"
    ! Add axis attributes to dimensions
    do i = 1, 4
       call check ( nf90_inq_varid ( &
            f % ncid, trim ( rname ( i ) ), varid ) )
       varids ( i ) = varid
    end do


    call check( nf90_put_att(f % ncid, varids (1),&
         trim('standard_name'),trim('longitude') ) )
    call check( nf90_put_att(f % ncid, varids (1),&
         trim('units'),trim('degrees_east') ) )
    call check( nf90_put_att(f % ncid, varids (1),&
         trim('axis'),trim('X') ) )
    call check( nf90_put_att(f % ncid, varids (2),&
         trim('standard_name'),trim('latitude') ) )
    call check( nf90_put_att(f % ncid, varids (2),&
         trim('units'),trim('degrees_north') ) )
    call check( nf90_put_att(f % ncid, varids (2),&
         trim('axis'),trim('Y') ) )
    call check( nf90_put_att(f % ncid, varids (3),&
         trim('standard_name'),trim('air_pressure') ) )
    call check( nf90_put_att(f % ncid, varids (3),&
         trim('units'),trim('hPa') ) )
    call check( nf90_put_att(f % ncid, varids (3),&
         trim('positive'),trim('down') ) )
    call check( nf90_put_att(f % ncid, varids (3),&
         trim('axis'),trim('Z') ) )
    call check( nf90_put_att(f % ncid, varids (4),&
         trim('standard_name'),trim('time') ) )
    call check( nf90_put_att(f % ncid, varids (4),&
         trim('units'),trim('hours since 1997-01-01 00:00:00') ) )
    call check( nf90_put_att(f % ncid, varids (4),&
         trim('calendar'),trim('standard') ) )

    print*,"Creating our Omega variables"
    ! Create quasi-geostrophic omega variables
    if (mode.eq.'Q')then
       do i = 1, size ( QG_omega_term_names )
          status = nf90_def_var ( f % ncid, trim ( QG_omega_term_names ( i ) ),&
               NF90_FLOAT, dimids, varid )
          if ( .not. ( status == nf90_enameinuse .or. status == NF90_NOERR ) ) &
               call check ( status )
          call check( nf90_put_att(f % ncid, varid,trim('description'), &
               trim(QG_omega_long_names(i)) ) )
          call check( nf90_put_att(f % ncid, varid,trim('units'),&
               trim('Pa s-1') ) )
       end do

    else

       ! Create generalized omega variables
       do i = 1, size ( omega_term_names )
          status = nf90_def_var ( f % ncid, trim ( omega_term_names ( i ) ), &
               NF90_FLOAT, dimids, varid )
          if ( .not. ( status == nf90_enameinuse .or. status == NF90_NOERR ) ) &
               call check ( status )
          call check( nf90_put_att(f % ncid, varid,trim('description'),&
               trim(omega_long_names(i)) ) )
          call check( nf90_put_att(f % ncid, varid,trim('units'),&
               trim('Pa s-1') ) )
       end do

       if (calc_b) then ! create omega b-variable if wanted
          status = nf90_def_var ( f % ncid, trim ( ome_b_name ), &
               NF90_FLOAT, dimids, varid )
          if ( .not. ( status == nf90_enameinuse .or. status == NF90_NOERR ) ) &
               call check ( status )
          call check( nf90_put_att(f % ncid, varid,trim('description'),&
               trim(ome_b_long_name) ) )
          call check( nf90_put_att(f % ncid, varid,trim('units'),&
               trim('Pa s-1') ) )
       end if

       ! Create ome_WRF variable
       status = nf90_def_var ( f % ncid, trim ( ome_name ), NF90_FLOAT, &
            dimids, varid )
       if ( .not. ( status == nf90_enameinuse .or. status == NF90_NOERR ) ) &
            call check ( status )
       call check( nf90_put_att(f % ncid, varid,trim('description'),&
            trim('omega from WRF') ) )
       call check( nf90_put_att(f % ncid, varid,trim('units'),&
            trim('Pa s-1') ) )

       ! Create height tendency variables
       do i = 1, size ( htend_term_names )
          status = nf90_def_var ( f % ncid, trim ( htend_term_names ( i ) ), &
               NF90_FLOAT, dimids, varid )
          if ( .not. ( status == nf90_enameinuse .or. status == NF90_NOERR ) ) &
               call check ( status )
          call check( nf90_put_att(f % ncid, varid,trim('description'),&
               trim(htend_long_names(i)) ) )
          call check( nf90_put_att(f % ncid, varid,trim('units'),&
               trim('m s-1') ) )
       end do

       if (calc_b) then ! create htend b-variable if wanted
          status = nf90_def_var ( f % ncid, trim ( htend_b_name ), &
               NF90_FLOAT, dimids, varid )
          if ( .not. ( status == nf90_enameinuse .or. status == NF90_NOERR ) ) &
               call check ( status )
          call check( nf90_put_att(f % ncid, varid,trim('description'),&
               trim(htend_b_long_name) ) )
          call check( nf90_put_att(f % ncid, varid,trim('units'),&
               trim('m s-1') ) )
       end if

       ! Create ztend_wrf variable
       status = nf90_def_var ( f % ncid, trim ( ztend_name ), NF90_FLOAT, &
            dimids, varid )
       if ( .not. ( status == nf90_enameinuse .or. status == NF90_NOERR ) ) &
            call check ( status )
       call check( nf90_put_att(f % ncid, varid,trim('description'),&
            trim('height tendency from WRF') ) )
       call check( nf90_put_att(f % ncid, varid,trim('units'),&
            trim('m s-1') ) )
    endif

    ! Create z variable
    status = nf90_def_var ( f % ncid, trim ( 'GHT' ), NF90_FLOAT, &
         dimids, varid )
    if ( .not. ( status == nf90_enameinuse .or. status == NF90_NOERR ) ) &
         call check ( status )
    call check( nf90_put_att(f % ncid, varid,trim('description'),&
         trim('geopotential height') ) )
    call check( nf90_put_att(f % ncid, varid,trim('units'),&
         trim('gpm') ) )

    ! Create psfc variable
    status = nf90_def_var ( f % ncid, trim ( 'PSFC' ), NF90_FLOAT, &
         (/ 1, 2, 4/), varid )
    if ( .not. ( status == nf90_enameinuse .or. status == NF90_NOERR ) ) &
         call check ( status )
    call check( nf90_put_att(f % ncid, varid,trim('description'),&
         trim('sfc pressure') ) )
    call check( nf90_put_att(f % ncid, varid,trim('units'),&
         trim('Pa') ) )

    call define_field3d(f % ncid, 'sigma_elcorr', &
         'Ellipticity correction for sigma', 'x', dimids)

    call define_field3d(f % ncid, 'zeta_elcorr', &
         'Ellipticity correction for zeta', 'x', dimids)

    if (forc) then

       ! Create vadv variable
       status = nf90_def_var ( f % ncid, trim ( 'vadv' ), NF90_FLOAT, &
            dimids, varid )
       if ( .not. ( status == nf90_enameinuse .or. status == NF90_NOERR ) ) &
            call check ( status )
       call check( nf90_put_att(f % ncid, varid,trim('description'),&
            trim('vorticity advection') ) )
       call check( nf90_put_att(f % ncid, varid,trim('units'),&
            trim('1/s^2') ) )
       ! Create tadv variable
       status = nf90_def_var ( f % ncid, trim ( 'tadv' ), NF90_FLOAT, &
            dimids, varid )
       if ( .not. ( status == nf90_enameinuse .or. status == NF90_NOERR ) ) &
            call check ( status )
       call check( nf90_put_att(f % ncid, varid,trim('description'),&
            trim('temperature advection') ) )
       call check( nf90_put_att(f % ncid, varid,trim('units'),&
            trim('K/s') ) )
       ! Create friction variable
       status = nf90_def_var ( f % ncid, trim ( 'fvort' ), NF90_FLOAT, &
            dimids, varid )
       if ( .not. ( status == nf90_enameinuse .or. status == NF90_NOERR ) ) &
            call check ( status )
       call check( nf90_put_att(f % ncid, varid,trim('description'),&
            trim('vorticity tendency due to friction') ) )
       call check( nf90_put_att(f % ncid, varid,trim('units'),&
            trim('1/s^2') ) )
       ! Create diabatic heating variable
       status = nf90_def_var ( f % ncid, trim ( 'diab' ), NF90_FLOAT, &
            dimids, varid )
       if ( .not. ( status == nf90_enameinuse .or. status == NF90_NOERR ) ) &
            call check ( status )
       call check( nf90_put_att(f % ncid, varid,trim('description'),&
            trim('Diabatic heating') ) )
       call check( nf90_put_att(f % ncid, varid,trim('units'),&
            trim('K/s') ) )
       ! Create ageostrophic vorticity tendency variable
       status = nf90_def_var ( f % ncid, trim ( 'ageo' ), NF90_FLOAT, &
            dimids, varid )
       if ( .not. ( status == nf90_enameinuse .or. status == NF90_NOERR ) ) &
            call check ( status )
       call check( nf90_put_att(f % ncid, varid,trim('description'),&
            trim('Ageostrophic vorticity tendency') ) )
       call check( nf90_put_att(f % ncid, varid,trim('units'),&
            trim('1/s^2') ) )
    end if

    ! Stop defining mode
    call check( nf90_enddef ( f % ncid ) )

    ! Copy dimension data to the new file type-variable
    allocate(f % xdim (f%dims(1)))
    allocate(f % ydim (f%dims(2)))

    f % xdim = (/(i, i=1,f%dims(1), 1)/)
    f % ydim = (/(i, i=1,f%dims(2), 1)/)

    print*,"Outputfile created!"

  end function create_out_file

  !***************************************************************************
  !> @brief Defines a 3D field in the NetCDF file.
  !!
  !! This subroutine defines a 3D variable in the NetCDF file with the given
  !! description and units.
  !!
  !! @param[in] ncid    NetCDF file ID.
  !! @param[in] name    Name of the variable.
  !! @param[in] desc    Description of the variable.
  !! @param[in] units   Units of the variable.
  !! @param[in] dimids  Dimension IDs for the variable.
  !***************************************************************************
  subroutine define_field3d(ncid, name, desc, units, dimids)
    integer, intent(in) :: ncid, dimids(4)
    integer :: varid
    character(*), intent(in) :: name, desc, units
    integer :: status

    status = nf90_def_var ( ncid, trim(name), NF90_FLOAT, &
         dimids, varid )
    !status = nf90_def_var ( 1, 'x', NF90_FLOAT, &
    !     [2,2,2,2], 5 )

    if ( .not. ( status == nf90_enameinuse .or. status == NF90_NOERR ) ) &
         call check ( status )
    call check( nf90_put_att(ncid, varid, 'description',&
         trim(desc) ) )
    call check( nf90_put_att(ncid, varid,'units',&
         trim(units) ) )

  end subroutine define_field3d

  !***************************************************************************
  !> @brief Opens an existing WRF input file.
  !!
  !! This function opens an existing WRF NetCDF file and retrieves its dimensions,
  !! attributes, and pressure levels.
  !!
  !! @param[in]  fname Name of the input file to open.
  !! @return     f     WRF file structure with metadata and dimensions populated.
  !***************************************************************************
  function open_wrf_file ( fname ) result ( f )
    character ( * ), intent ( in ) :: fname
    type ( wrf_file ) :: f
    integer :: i, dimid, varid, nres

    print*,"Opening file: ",fname
    call check( nf90_open ( fname, NF90_WRITE, f % ncid ) )

    print*,"Inquiring dimensions from the input file..."
    do i = 1, 4
       call check ( nf90_inq_dimid ( &
            f % ncid, trim ( rname ( i ) ), dimid ) )
       call check ( nf90_inquire_dimension ( &
            f % ncid, dimid, len = f % dims ( i ) ) )
    end do

    ! Number of resolutions in the solving of omega equation
    nres=1+int(log(max(f%dims(1),f%dims(2),f%dims(3))/5.)/log(2.))
    allocate(f % nlon ( nres ), f % nlat ( nres ), f % nlev ( nres ) )
    allocate(f % dx ( nres ), f % dy ( nres ), f % dlev ( nres ) )

    print*,"Getting attributes from the input file..."
    call check ( nf90_get_att ( &
         f % ncid, NF90_GLOBAL, 'DX', f % dx(1) ) )
    call check ( nf90_get_att ( &
         f % ncid, NF90_GLOBAL, 'DY', f % dy(1) ) )
    call check ( nf90_get_att ( &
         f % ncid, NF90_GLOBAL, 'CU_PHYSICS', f % wrf_cu_phys ) )

    print*,"Getting pressure level information from the input file..."
    allocate ( f % pressure_levels ( f % dims ( 3 ) ) )
    call check ( nf90_inq_varid ( f % ncid, 'LEV', varid ) )
    call check ( nf90_get_var ( f % ncid, varid, f % pressure_levels, &
         start = [ 1 ], count = [ size ( f % pressure_levels ) ] ) )

    f % nlon(1) = f % dims(1)
    f % nlat(1) = f % dims(2)
    f % nlev(1) = f % dims(3)
    f % dlev(1) = f % pressure_levels(2) - f % pressure_levels(1)

    ! Number of different resolutions in solving the equation = nres
    ! Choose so that the coarsest grid has at least 5 points
    do i=2,nres
       f % nlon(i) = max(f % nlon(i-1)/2,5)
       f % nlat(i) = max(f % nlat(i-1)/2,5)
       f % nlev(i) = max(f % nlev(i-1)/2,5)
       f % dx(i) = f % dx(1)*real(f % nlon(1))/real(f % nlon(i))
       f % dy(i) = f % dy(1)*real(f % nlat(1))/real(f % nlat(i))
       f % dlev(i) = f % dlev(1)*real(f % nlev(1))/real(f % nlev(i))
    enddo

    print*,"Getting time information from the input file..."
    allocate ( f % times ( f % dims ( 4 ) ) )
    !    f % times = (/ ( i, i = 0, f % dims ( 4 ) - 1, 1 ) /)
    call check ( nf90_inq_varid ( f % ncid, 'XTIME', varid ) )
    call check ( nf90_get_var ( f % ncid, varid, f % times, &
         start = [ 1 ], &
         count = [ size ( f % times ) ] ) )
    f % times = f % times * 60

    print*,"Getting coriolisparameter from the input file..."
    allocate ( f % corpar ( f % dims ( 2 ) ) )
    call check ( nf90_inq_varid ( f % ncid, 'F', varid ) )
    call check ( nf90_get_var ( f % ncid, varid, f % corpar, &
         start = [ 1, 1, 2 ], &
         count = [ 1, size ( f % corpar ), 1 ] ) )

    print*,"Input file opened succesfully!"
  end function open_wrf_file

  !***************************************************************************
  !> @brief Closes a WRF file.
  !!
  !! This subroutine closes an open WRF NetCDF file.
  !!
  !! @param[in] file WRF file structure to close.
  !***************************************************************************
  subroutine close_wrf_file ( file )
    type ( wrf_file ) :: file
    call check ( nf90_close ( file % ncid ) )
  end subroutine close_wrf_file

  !***************************************************************************
  !> @brief Opens an existing WRF output file.
  !!
  !! This function opens an existing WRF NetCDF output file and retrieves its dimensions.
  !!
  !! @param[in]  fname Name of the output file to open.
  !! @return     f     WRF file structure with metadata and dimensions populated.
  !***************************************************************************
  function open_out_file ( fname ) result ( f )
    character ( * ), intent ( in ) :: fname
    type ( wrf_file ) :: f
    integer :: i,dimid

    call check( nf90_open ( fname, NF90_WRITE, f % ncid ) )

    do i = 1, 4
       call check ( nf90_inq_dimid ( &
            f % ncid, trim ( rname ( i ) ), dimid ) )
       !       dimids ( i ) = dimid
       call check ( nf90_inquire_dimension ( &
            f % ncid, dimid, len = f % dims ( i ) ) )
    end do

  end function open_out_file

  !***************************************************************************
  !> @brief Reads a 2D variable from the WRF file.
  !!
  !! This function reads a 2D variable (longitude x latitude) from the WRF file
  !! for a given time step.
  !!
  !! @param[in]  file  WRF file structure.
  !! @param[in]  time  Time index to read data for.
  !! @param[in]  names Array of variable names to read and sum.
  !! @return     real2d 2D array of the variable data.
  !***************************************************************************
  function real2d ( file, time, names )
    type ( wrf_file ), intent ( in ) :: file
    integer, intent ( in ) :: time
    character ( * ), dimension ( : ), intent ( in ) :: names
    real, dimension ( :, : ), allocatable :: real2d
    integer :: i

    real2d = data ( names ( 1 ) )
    do i = 2, size ( names )
       real2d = real2d + data ( names ( i ) )
    end do

  contains

    function data ( name )
      real, dimension ( :, : ), allocatable :: data
      character ( * ) :: name
      integer :: varid

      allocate ( data ( file % dims ( 1 ), file % dims ( 2 ) ) )
      call check ( nf90_inq_varid ( file % ncid, trim ( name ), varid ) )
      call check ( nf90_get_var ( file % ncid, varid, data, &
           start = [ 1, 1, time ], count = [ shape ( data ), 1 ] ) )

    end function data

  end function real2d

  !***************************************************************************
  !> @brief Reads a 3D variable from the WRF file.
  !!
  !! This function reads a 3D variable (longitude x latitude x level) from the WRF file
  !! for a given time step.
  !!
  !! @param[in]  file  WRF file structure.
  !! @param[in]  time  Time index to read data for.
  !! @param[in]  names Array of variable names to read and sum.
  !! @return     real3d 3D array of the variable data.
  !***************************************************************************
  function real3d ( file, time, names )
    type ( wrf_file ), intent ( in ) :: file
    integer, intent ( in ) :: time
    character ( * ), dimension ( : ), intent ( in ) :: names
    real, dimension ( :, :, : ), allocatable :: real3d
    integer :: i

    !    print*,'Reading ',trim(names(1))
    real3d = data ( names ( 1 ) )

    do i = 2, size ( names )
       !       print*,'Reading ',trim(names(i))
       real3d = real3d + data ( names ( i ) )
    end do

  contains

    function data ( name )
      real, dimension ( :, :, : ), allocatable :: data
      character ( * ) :: name
      integer :: varid

      allocate ( data ( &
           file % dims ( 1 ), file % dims ( 2 ), file % dims ( 3 ) ) )

      call check ( nf90_inq_varid ( file % ncid, trim ( name ), varid ) )
      call check ( nf90_get_var ( file % ncid, varid, data, &
           start = [ 1, 1, 1, time ], count = [ shape ( data ), 1 ] ) )
    end function data

  end function real3d

  !***************************************************************************
  !> @brief Checks the status of a NetCDF operation.
  !!
  !! This subroutine checks the status of a NetCDF operation and stops execution
  !! if an error is encountered.
  !!
  !! @param[in] status Status code returned by a NetCDF operation.
  !***************************************************************************
  subroutine check ( status )
    integer, intent ( in ) :: status
    if ( status /= NF90_NOERR ) then
       write(*,*) 'Error in ', trim ( nf90_strerror ( status ) )
       stop 1
    end if
  end subroutine check

  !***************************************************************************
  !> @brief Defines a dimension in the NetCDF file.
  !!
  !! This subroutine defines a dimension in the NetCDF file, optionally as unlimited.
  !!
  !! @param[in] ncid      NetCDF file ID.
  !! @param[in] dim_name  Name of the dimension.
  !! @param[out] dimid    Dimension ID.
  !! @param[in] len       Length of the dimension.
  !! @param[in] unlimit   Logical flag to define the dimension as unlimited.
  !***************************************************************************
  subroutine def_dim(ncid,dim_name,dimid,len,unlimit)
    implicit none
    character(*) :: dim_name
    integer :: ncid,dimid,varid,len
    logical :: unlimit
     print*,"Creating dimension ",trim(dim_name)," with length ",len
    if(unlimit)then
       call check( nf90_def_dim(ncid, trim(dim_name), NF90_UNLIMITED, dimid) )
    else
       call check( nf90_def_dim(ncid, trim(dim_name), len, dimid) )
    endif
    print*,"Passed nf90_def_dim",varid
    call check( nf90_def_var(ncid, trim(dim_name), NF90_REAL, dimid, varid) )

  end subroutine def_dim

  !***************************************************************************
  !> @brief Removes underscores from a string.
  !!
  !! This function removes underscores from a string and replaces them with spaces.
  !!
  !! @param[in]  input_string Input string to process.
  !! @return     output_string Processed string with underscores removed.
  !***************************************************************************
  function remove_underscores(input_string) result(output_string)
    character(*), intent(in) :: input_string
    character(len(input_string)) :: output_string
    integer :: i, j

    j = 1
    do i = 1, len_trim(input_string)
       if (input_string(i:i) /= '_') then
          output_string(j:j) = input_string(i:i)
          j = j + 1
       end if
    end do
    output_string(j:) = ' '  ! Fill the rest with spaces
  end function remove_underscores

end module mod_wrf_file
