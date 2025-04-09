!> @file mod_time_step_loop.f90
!! @brief This module contains the main time-stepping loop and related subroutines.
!!
!! This module implements the main time-stepping loop for the simulation, as well as
!! various utility subroutines for reading, writing, and processing data.
!!
!! @note Ensure that this module is properly included in the relevant parts
!! of the project to access its functionality.

module mod_time_step_loop
  use mod_wrf_file
  use mod_subrs
  use mod_omega
  use mod_const
  use mod_common_subrs
  implicit none

contains

  !***************************************************************************
  !> @brief Main time-stepping loop to calculate each WRF output grid.
  !!
  !! This subroutine implements the main time-stepping loop, reading input data,
  !! solving the QG Omega or Generalized Omega equations, and writing output data.
  !!
  !! @param[in]  wrfin_file Input WRF file structure.
  !! @param[in]  out_file   Output WRF file structure.
  !! @param[in]  param      Simulation parameters.
  !! @param[in]  calc_b     Logical flag to use WRF omegas as boundary conditions.
  !! @param[in]  debug      Logical flag to enable debug output.
  !***************************************************************************
  subroutine time_step_loop ( wrfin_file, out_file, param, calc_b, debug)
    ! This subroutine contains the main time stepping loop. It gets both
    ! input and output files as input arguments.
    type ( wrf_file ), intent(in) :: wrfin_file, out_file
    type ( parameters), intent(in) :: param
    logical,           intent(in) :: debug, calc_b
    real, dimension ( :, :, : ), pointer :: T, u, v, z
    real, dimension ( :, :, :, : ), allocatable :: omegas, hTends, omegas_QG
    real, dimension ( :, :, : ),    allocatable :: dT_dt, du_dt, dv_dt, dz_dt, &
         fx, fy, q, w, mulfact,zeta, &
         zetatend,ukhi,vkhi,sigma,vadv, dudp, dvdp, &
         tadv,fvort,avortt, zeta_elcorr, sigma_elcorr
    real, dimension ( :, : ),       allocatable :: mu_inv, p_sfc
    integer, dimension ( : ),       allocatable :: tdim
    integer :: time, i

    associate ( &
         time_1 => param % time_1, &
         time_n => param % time_n, &
         alfa => param % alfa, &
         toler => param % toler, &
         ny1 => param % ny1, &
         ny2 => param % ny2, &
         mode => param % mode, &
         calc_omegas => param % calc_omegas, &
         calc_div => param % calc_div, &
         forc => param % forc, &
         ellipticity_correction => param % ellipticity_correction, &

         nlon   => wrfin_file % dims ( 1 ), &
         nlat   => wrfin_file % dims ( 2 ), &
         nlev   => wrfin_file % dims ( 3 ), &
         dx     => wrfin_file % dx(1), &
         dy     => wrfin_file % dy(1), &
         dlev     => wrfin_file % dlev(1), &
         p_levs => wrfin_file % pressure_levels, &
         corpar => wrfin_file % corpar )

    if (calc_b) then
       n_terms = n_terms + 1
    end if

    allocate ( tdim ( time_n-time_1+1 ) )
    allocate ( hTends ( nlon, nlat, nlev, n_terms ) )
    allocate ( omegas ( nlon, nlat, nlev, n_terms ) )
    allocate ( omegas_QG ( nlon, nlat, nlev, 3 ) )
    allocate ( ukhi (nlon, nlat, nlev) )
    allocate ( vkhi (nlon, nlat, nlev) )
    allocate ( avortt (nlon, nlat, nlev) )
    allocate ( fvort (nlon, nlat, nlev) )
    allocate ( sigma_elcorr (nlon, nlat, nlev) )
    allocate ( zeta_elcorr (nlon, nlat, nlev) )

    print*,'In time_step_loop'
    call read_T_u_v_z ( wrfin_file, time_1 - 2 )
    call read_T_u_v_z ( wrfin_file, time_1 - 1 )
    i=0
    do time = time_1, time_n
       i=i+1
       print*,'Time step: ', time
       tdim(i)=time
       call read_T_u_v_z ( wrfin_file, time )
       mu_inv = 1.0 / real2d ( wrfin_file, time, [ 'MU ', 'MUB' ]  )
       q      = diabatic_heating ( wrfin_file, time, mu_inv )
       fx     = friction ( wrfin_file, time, 'U', mu_inv )
       fy     = friction ( wrfin_file, time, 'V', mu_inv )
       p_sfc  = real2d ( wrfin_file, time, [ 'PSFC' ]  )
       w      = real3d ( wrfin_file, time, [ 'WW' ]  )
       zeta   = curl_cart ( u, v, dx, dy )
       zetatend = curl_cart ( du_dt, dv_dt, dx, dy )
       mulfact  = calmul(p_sfc, p_levs, nlev)
       sigma = define_sigma(T, p_levs)

       dudp = pder(u,dlev)
       dvdp = pder(v,dlev)

       call calculate_ellipticity_correction (sigma, sigmamin, etamin, zeta, corpar,&
            dudp, dvdp, sigma_elcorr, zeta_elcorr)
       call write3d ( out_file, time-time_1+1, 'zeta_elcorr', zeta_elcorr)
       call write3d ( out_file, time-time_1+1, 'sigma_elcorr', sigma_elcorr)

       !   Calculation of velocity potential
       if(calc_div) call irrotationalWind(u,v,dx,dy,uKhi,vKhi)



       if ( calc_omegas ) then
          call calculate_omegas( wrfin_file, T, u, v, w, z, &
               q, fx, fy, dT_dt, zeta, zetatend, uKhi, vKhi, sigma, &
               mulfact, param, debug, calc_b, &
               omegas, omegas_QG, sigma_elcorr, zeta_elcorr, dudp, dvdp )
       else
          omegas = read_omegas ( out_file, time-time_1+1 )
       end if

       call calculate_tendencies ( omegas, T, u, v, w, z, p_levs, &
            dx, dy, corpar, q, fx, fy, dz_dt, dT_dt, zeta, zetatend, &
            uKhi, vKhi, sigma, mulfact, calc_b, hTends, vadv, tadv, &
            fvort, avortt )

       ! Write data to the output file
       if ( mode .eq. 'Q' ) then
          call write_omegas_QG ( out_file, time-time_1+1, omegas_QG, calc_b )
       else if ( mode .eq. 'G' ) then
          if (calc_omegas) then
             call write_omegas ( out_file, time-time_1+1, calc_b, omegas )
             call write3d ( out_file, time-time_1+1, ome_name, w)
          end if
          call write_height_tendencies(out_file,time-time_1+1,calc_b,hTends)
          call write3d ( out_file, time-time_1+1, ztend_name, dz_dt )
       end if

       call write3d ( out_file, time-time_1+1, 'GHT', z)
       call write2d ( out_file, time-time_1+1, 'PSFC', p_sfc)
       if (forc) then
          call write3d ( out_file, time-time_1+1, 'vadv', vadv)
          call write3d ( out_file, time-time_1+1, 'tadv', tadv)
          call write3d ( out_file, time-time_1+1, 'fvort', fvort)
          call write3d ( out_file, time-time_1+1, 'diab', q)
          call write3d ( out_file, time-time_1+1, 'ageo', avortt)

          !call write3d ( out_file, time-time_1+1, 'mulfact', mulfact)
          call write3d ( out_file, time-time_1+1, 'sigma', sigma)
       end if

    end do
    print*,"Finished time stepping loop."

    ! Write dimension data to the output file
    if ( calc_omegas .or. mode .eq. 'Q' ) then
       call write_dimensions ( out_file, tdim )
    end if
    
  end associate

contains

  !***************************************************************************
  !> @brief Writes dimension data to the output file.
  !!
  !! This subroutine writes the spatial and temporal dimensions to the output file.
  !!
  !! @param[in]  file Output WRF file structure.
  !! @param[in]  tdim Time dimension array.
  !***************************************************************************
  subroutine write_dimensions ( file, tdim )
    ! This subroutine writes dimension data to the output file.
    type ( wrf_file ) :: file
    integer,dimension(:) :: tdim
    integer :: varids(4),varid

    do i = 1, 4
       call check ( nf90_inq_varid ( &
            file % ncid, trim ( rname ( i ) ), varid ) )
       varids ( i ) = varid
    end do

    call check( nf90_put_var(file % ncid, varids(1), file%xdim(:)) )
    call check( nf90_put_var(file % ncid, varids(2), file%ydim(:)) )
    call check( nf90_put_var(file % ncid, varids(3), &
         file%pressure_levels(:)/100.) )
    call check( nf90_put_var(file % ncid, varids(4), tdim(:)) )

  end subroutine write_dimensions

  !***************************************************************************
  !> @brief Reads T, u, v, and z fields from the input file.
  !!
  !! This subroutine reads temperature, zonal wind, meridional wind, and geopotential
  !! height fields from the input file and calculates their tendencies.
  !!
  !! @param[in]  file Input WRF file structure.
  !! @param[in]  time Time index to read data for.
  !***************************************************************************
  subroutine read_T_u_v_z ( file, time )
    ! This subroutine reads T,u,v and Z fields from the input file
    ! and calculate tendencies of them
    type ( wrf_file ) :: file
    integer, intent ( in ) :: time
    real, dimension ( :, :, :, :, : ), allocatable, target, save :: dtvars
    character ( 3 ), dimension ( 4 ), parameter :: &
         dtvar_names = [ 'TT ', 'UU ', 'VV ', 'GHT' ]
    integer, save :: prev = 1, curr = 2, next = 3
    real :: dt2inv
    integer :: i
    print*,'In read_T_u_v_z'

    if ( .not. allocated ( dtvars ) ) then
       associate ( &
            nlon => file % dims ( 1 ), &
            nlat => file % dims ( 2 ), &
            nlev => file % dims ( 3 ) )
       allocate ( dtvars ( nlon, nlat, nlev, 4, 3 ) )
     end associate
  end if
  prev = curr
  curr = next
  next = mod ( next, 3 ) + 1
  T => dtvars ( :, :, :, 1, curr )
  u => dtvars ( :, :, :, 2, curr )
  v => dtvars ( :, :, :, 3, curr )
  z => dtvars ( :, :, :, 4, curr )
  do i = 1, 4
     dtvars ( :, :, :, i, next ) = real3d ( file, time + 1, &
          [ trim ( dtvar_names ( i ) ) ] )
  end do
  if ( time .gt. lbound ( file % times, 1 ) ) then
     dt2inv = 1.0 / ( file % times ( time + 1 ) - file % times ( time - 1 ))
     print*,"dt2inv",dt2inv,file%times(time + 1), file%times(time - 1)
     dT_dt = (dtvars( :, :, :, 1, next) - dtvars( :, :, :, 1,prev)) * dt2inv
     du_dt = (dtvars( :, :, :, 2, next) - dtvars( :, :, :, 2,prev)) * dt2inv
     dv_dt = (dtvars( :, :, :, 3, next) - dtvars( :, :, :, 3,prev)) * dt2inv
     dz_dt = (dtvars( :, :, :, 4, next) - dtvars( :, :, :, 4,prev)) * dt2inv
  end if
end subroutine read_T_u_v_z

  !***************************************************************************
  !> @brief Reads omega terms from the output file.
  !!
  !! This function reads omega terms from the output file for a given time step.
  !!
  !! @param[in]  file Output WRF file structure.
  !! @param[in]  time Time index to read data for.
  !! @return     read_omegas Omega terms (4D array: longitude x latitude x level x term).
  !***************************************************************************
function read_omegas ( file, time )
  type ( wrf_file ) :: file
  integer, intent ( in ) :: time
  real, dimension ( :, :, :, : ), allocatable :: read_omegas
  integer :: i
  associate ( nlon => file % dims ( 1 ), &
       nlat => file % dims ( 2 ), &
       nlev => file % dims ( 3 ), &
       nterms => size ( omega_term_names ) )
  allocate ( read_omegas ( nlon, nlat, nlev, nterms ) )
  do i = 1, nterms
     read_omegas ( :, :, :, i ) &
          = real3d ( file, time, [ trim ( omega_term_names ( i ) ) ] )
  end do
end associate
end function read_omegas

  !***************************************************************************
  !> @brief Computes diabatic heating.
  !!
  !! This function calculates diabatic heating based on input data and pressure levels.
  !!
  !! @param[in]  file   Input WRF file structure.
  !! @param[in]  time   Time index to read data for.
  !! @param[in]  mu_inv Inverse of the dry air mass.
  !! @return     q      Diabatic heating (3D array: longitude x latitude x level).
  !***************************************************************************
function diabatic_heating ( file, time, mu_inv ) result ( q )
type ( wrf_file ), intent ( in ) :: file
integer, intent ( in ) :: time
real, dimension ( :, : ), intent ( in ) :: mu_inv
real, dimension ( :, :, : ), allocatable :: q
integer :: lev
real :: a_lev
allocate ( q ( file % dims ( 1 ), file % dims ( 2 ), &
     file % dims ( 3 ) ) )
associate ( plev => file % pressure_levels )
select case ( file % wrf_cu_phys )
case ( 1 )
   !           q = real3d ( file, time, &
   !                [ 'RTHSHTEN', 'RTHCUTEN', 'RTHRATEN', 'RTHBLTEN' ] )
   q = real3d ( file, time, &
        [ 'RTHCUTEN', 'RTHRATEN', 'RTHBLTEN' ] )
end select
do lev = 1, size ( q, 3 )
   q ( :, :, lev ) = q ( :, :, lev ) * mu_inv
end do
q = q + real3d ( file, time, [ 'H_DIABATIC' ] )
do lev = 1, size ( plev )
   a_lev = ( plev ( lev ) / 100000.0 ) ** ( r / cp )
   q ( :, :, lev ) = q ( :, :, lev ) * a_lev
end do
end associate
end function diabatic_heating

  !***************************************************************************
  !> @brief Loads in friction-related terms from the WRF input file.
  !!
  !! This function retrieves friction terms for the zonal or meridional wind components.
  !!
  !! @param[in]  file      Input WRF file structure.
  !! @param[in]  time      Time index to read data for.
  !! @param[in]  direction Direction ('U' for zonal, 'V' for meridional).
  !! @param[in]  mu_inv    Inverse of the dry air mass.
  !! @return     friction  Friction terms (3D array: longitude x latitude x level).
  !***************************************************************************
function friction ( file, time, direction, mu_inv )
type ( wrf_file ), intent ( in ) :: file
integer, intent ( in ) :: time
character ( 1 ), intent ( in ) :: direction
real, dimension ( :, : ), intent ( in ) :: mu_inv
real, dimension ( :, :, : ), allocatable :: friction
integer :: lev
allocate ( friction ( file % dims ( 1 ), file % dims ( 2 ), &
   file % dims ( 3 ) ) )
select case ( file % wrf_cu_phys )
case ( 1 )
   !         friction = real3d ( file, time, &
   !              [ 'R'//direction//'SHTEN', 'R'//direction//'CUTEN', &
   !              'R'//direction//'BLTEN' ] )
 friction = real3d ( file, time, &
      [ 'R'//direction//'CUTEN', &
      'R'//direction//'BLTEN' ] )
end select
do lev = 1, size ( friction, 3 )
 friction ( :, :, lev ) = friction ( :, :, lev ) * mu_inv
end do
end function friction

  !***************************************************************************
  !> @brief Writes height tendencies to the output file.
  !!
  !! This subroutine writes height tendency terms to the output file.
  !!
  !! @param[in]  file    Output WRF file structure.
  !! @param[in]  time    Time index to write data for.
  !! @param[in]  calc_b  Logical flag to enable additional calculations.
  !! @param[in]  hTends  Height tendencies (4D array: longitude x latitude x level x term).
  !***************************************************************************
subroutine write_height_tendencies ( file, time, calc_b, hTends )
type ( wrf_file ), intent ( in ) :: file
integer, intent ( in ) :: time
logical, intent ( in ) :: calc_b
real, dimension ( :, :, :, : ) :: hTends
integer :: t, varid
associate ( &
   ncid => file % ncid, &
   nlon => file % dims ( 1 ), &
   nlat => file % dims ( 2 ), &
   nlev => file % dims ( 3 ) )
do t = 1, size ( htend_term_names )
 call check ( nf90_inq_varid ( ncid, &
      trim ( htend_term_names ( t ) ), varid ) )
 call check ( nf90_put_var ( ncid, varid, &
      hTends ( :, :, :, t ), &
      [ 1, 1, 1, time ], &
      [ nlon, nlat, nlev, 1 ] ) )
end do

if (calc_b) then
 call check ( nf90_inq_varid ( ncid, &
      trim ( htend_b_name ), varid ) )
 call check ( nf90_put_var ( ncid, varid, &
      hTends ( :, :, :, 8 ), &
      [ 1, 1, 1, time ], &
      [ nlon, nlat, nlev, 1 ] ) )
end if
end associate
end subroutine write_height_tendencies

  !***************************************************************************
  !> @brief Writes omega terms to the output file.
  !!
  !! This subroutine writes omega terms to the output file for a given time step.
  !!
  !! @param[in]  file    Output WRF file structure.
  !! @param[in]  time    Time index to write data for.
  !! @param[in]  calc_b  Logical flag to enable additional calculations.
  !! @param[in]  omegas  Omega terms (4D array: longitude x latitude x level x term).
  !***************************************************************************
subroutine write_omegas ( file, time, calc_b, omegas )
type ( wrf_file ), intent ( in ) :: file
integer, intent ( in ) :: time
logical, intent ( in ) :: calc_b
real, dimension ( :, :, :, : ) :: omegas
integer :: t, varid
associate ( &
 ncid => file % ncid, &
 nlon => file % dims ( 1 ), &
 nlat => file % dims ( 2 ), &
 nlev => file % dims ( 3 ) )
do t = 1, size ( omega_term_names )
call check ( nf90_inq_varid ( ncid, &
    trim ( omega_term_names ( t ) ), varid ) )
call check ( nf90_put_var ( ncid, varid, &
    omegas ( :, :, :, t ), &
    [ 1, 1, 1, time ], &
    [ nlon, nlat, nlev, 1 ] ) )
end do

if (calc_b) then
call check ( nf90_inq_varid ( ncid, &
    trim ( ome_b_name ), varid ) )
call check ( nf90_put_var ( ncid, varid, &
    omegas ( :, :, :, 8 ), &
    [ 1, 1, 1, time ], &
    [ nlon, nlat, nlev, 1 ] ) )
end if

end associate
end subroutine write_omegas

  !***************************************************************************
  !> @brief Writes QG omega terms to the output file.
  !!
  !! This subroutine writes quasi-geostrophic (QG) omega terms to the output file.
  !!
  !! @param[in]  file       Output WRF file structure.
  !! @param[in]  time       Time index to write data for.
  !! @param[in]  omegas_QG  QG omega terms (4D array: longitude x latitude x level x term).
  !! @param[in]  calc_b     Logical flag to enable additional calculations.
  !***************************************************************************
subroutine write_omegas_QG ( file, time, omegas_QG, calc_b )
type ( wrf_file ), intent ( in ) :: file
integer, intent ( in ) :: time
logical, intent(in) :: calc_b
real, dimension ( :, :, :, : ) :: omegas_QG
integer :: t, varid
associate ( &
ncid => file % ncid, &
nlon => file % dims ( 1 ), &
nlat => file % dims ( 2 ), &
nlev => file % dims ( 3 ) )
do t = 1, size ( QG_omega_term_names )
call check ( nf90_inq_varid ( ncid, &
  trim ( QG_omega_term_names ( t ) ), varid ) )
call check ( nf90_put_var ( ncid, varid, &
  omegas_QG ( :, :, :, t ), &
  [ 1, 1, 1, time ], &
  [ nlon, nlat, nlev, 1 ] ) )
end do

if (calc_b) then
call check ( nf90_inq_varid ( ncid, &
  trim ( ome_b_name ), varid ) )
call check ( nf90_put_var ( ncid, varid, &
  omegas ( :, :, :, 8 ), &
  [ 1, 1, 1, time ], &
  [ nlon, nlat, nlev, 1 ] ) )
end if

end associate
end subroutine write_omegas_QG

  !***************************************************************************
  !> @brief Writes a 3D grid to the output file.
  !!
  !! This subroutine writes 3D data (longitude x latitude x level) to the output file.
  !!
  !! @param[in]  file Output WRF file structure.
  !! @param[in]  time Time index to write data for.
  !! @param[in]  name Variable name to write.
  !! @param[in]  data 3D data array to write.
  !***************************************************************************
subroutine write3d ( file, time, name, data )
type ( wrf_file ), intent ( in ) :: file
integer, intent ( in ) :: time
character(*) :: name
real, dimension ( :, :, : ) :: data
integer :: varid

associate ( &
ncid => file % ncid, &
nlon => file % dims ( 1 ), &
nlat => file % dims ( 2 ), &
nlev => file % dims ( 3 ) )

call check ( nf90_inq_varid ( ncid, &
trim ( name ), varid ) )
call check ( nf90_put_var ( ncid, varid, &
data ( :, :, : ), &
[ 1, 1, 1, time ], &
[ nlon, nlat, nlev, 1 ] ) )

end associate

end subroutine write3d

  !***************************************************************************
  !> @brief Writes a 2D grid to the output file.
  !!
  !! This subroutine writes 2D data (longitude x latitude) to the output file.
  !!
  !! @param[in]  file Output WRF file structure.
  !! @param[in]  time Time index to write data for.
  !! @param[in]  name Variable name to write.
  !! @param[in]  data 2D data array to write.
  !***************************************************************************
subroutine write2d ( file, time, name, data )
type ( wrf_file ), intent ( in ) :: file
integer, intent ( in ) :: time
character(*) :: name
real, dimension ( :, : ) :: data
integer :: varid

associate ( &
ncid => file % ncid, &
nlon => file % dims ( 1 ), &
nlat => file % dims ( 2 ) )

call check ( nf90_inq_varid ( ncid, &
trim ( name ), varid ) )
call check ( nf90_put_var ( ncid, varid, &
data, [ 1, 1, time ], &
[ nlon, nlat, 1 ] ) )

end associate

end subroutine write2d
end subroutine time_step_loop

end module mod_time_step_loop
