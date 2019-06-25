MODULE cmorph_v0_x_output_nc

! This module contains the I/O routines related to initializing and
! writing the output for the CMORPH precipitation analyses.

USE netcdf
USE kinds

IMPLICIT NONE

INTEGER (KIND=int_kind), PUBLIC :: &
   ncid, &      ! NetCDF ID
   status       ! NetCDF error status return
          
! Public member functions and subroutines
PUBLIC :: &
   cmorph_create_cdf, &
   cmorph_def_dim, &
   cmorph_def_var, &
   cmorph_write_var, &
   handle_err
   
CONTAINS

!==============================================================================
      SUBROUTINE cmorph_create_cdf(fname_nc)
!==============================================================================
!
! Create the NetCDF dataset and define global attributes
!
!------------------------------------------------------------------------------
      CHARACTER (LEN=*), INTENT(IN) :: &
          fname_nc      ! NetCDF filename
      CHARACTER (LEN=10) :: &
          date
      CHARACTER (LEN=8) :: &
          time
      CHARACTER (LEN=8) :: &
          zone
      INTEGER (KIND=int_kind), DIMENSION(8) :: &
          values
      CHARACTER (LEN=512) :: &
          history_str
      CHARACTER (LEN=255) :: &
          cmd
      CHARACTER (LEN=80) :: &
          prog, host, login, user
!------------------------------------------------------------------------------
      prog = 'cmorph_create_cdf'
!------------------------------------------------------------------------------
! Get current date and time, system host name, command line, and login name
      call date_and_time(VALUES=values)
      write(date,'(i4,''-'',i2.2,''-'',i2.2)') values(1),values(2),values(3)
      write(time,'(i2.2,'':'',i2.2,'':'',i2.2)') values(5),values(6),values(7)
      write(zone,'(''UTC'',i3,''00'')') values(4)/60
      call hostnm(host)
      call get_command(cmd)
      call getlog(login)
      call get_environment_variable("USER", user)
      history_str = date // " " // time // " " // zone // &
                    ": NetCDF generated on host " // &
                    trim(host) // " by user " // trim(user) // &
                    " via command " // trim(cmd)
!------------------------------------------------------------------------------
! Create dataset
      status = nf90_create(fname_nc, NF90_CLOBBER, ncid)
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
!------------------------------------------------------------------------------
! Global attributes
      status = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', 'CF-1.5')
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
      
      status = nf90_put_att(ncid, NF90_GLOBAL, 'title', &
             'CPC Morphing Technique (CMORPH) Global Precipitation Analysis' )
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
      
      status = nf90_put_att(ncid, NF90_GLOBAL, 'version', &
             'CMORPH version 0.x' )
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))

      status = nf90_put_att(ncid, NF90_GLOBAL, 'history', trim(history_str))
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
      
      status = nf90_put_att(ncid, NF90_GLOBAL, 'source', &
             "Precipitation estimates derived from the passive microwave sensors aboard &
              &DMSP 13, 14 & 15 (SSM/I), NOAA-15, 16, 17, & 18 (AMSU-B), and AMSR-E and TMI &
              &aboard NASA Aqua and TRMM")
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
      
      status = nf90_put_att(ncid, NF90_GLOBAL, 'documentation', &
             'http://rda.ucar.edu/datasets/ds502.0')
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))

      status = nf90_put_att(ncid, NF90_GLOBAL, 'references', &
              "http://www.cpc.ncep.noaa.gov/products/janowiak/cmorph.shtml")
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))

      status = nf90_put_att(ncid, NF90_GLOBAL, 'institution', &
               "NetCDF produced at National Center for Atmospheric Research.  &
               &Original binary data produced at NOAA Climate Prediction Center")
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))

      status = nf90_put_att(ncid, NF90_GLOBAL, 'acknowledgments', &
             "The data for this study are from the Research Data Archive (RDA) &
              &which is maintained by the Computational and Information Systems &
              &Laboratory (CISL) at the National Center for Atmospheric Research &
              &(NCAR). NCAR is sponsored by the National Science Foundation (NSF). &
              &The original data are available from the RDA (http://rda.ucar.edu) &
              &in dataset number ds502.0")
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
      
      status = nf90_put_att(ncid, NF90_GLOBAL, 'contact', &
             'Thomas Cram (tcram@ucar.edu)')
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))

!------------------------------------------------------------------------------
! leave define mode
      status = nf90_enddef(ncid)
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))

!==============================================================================
      END SUBROUTINE cmorph_create_cdf
!==============================================================================

!==============================================================================
      SUBROUTINE cmorph_def_dim(lon_len, lat_len, lon, lat, &
                                time_range)
!==============================================================================
!
! Define dimensions of CMORPH data
!
!------------------------------------------------------------------------------
! Input
      INTEGER (KIND=int_kind), INTENT(IN) :: &
          lon_len, &
          lat_len
      REAL (KIND=real_kind), DIMENSION(lon_len), INTENT(IN) :: &
          lon
      REAL (KIND=real_kind), DIMENSION(lat_len), INTENT(IN) :: &
          lat
      INTEGER (KIND=int_kind), DIMENSION(2), INTENT(IN) :: &
          time_range
!------------------------------------------------------------------------------
! Local

! Dimension ID
      INTEGER (KIND=int_kind) :: &
          lon_dim, &   ! Longitude dimension ID
          lat_dim, &   ! Latitude dimension ID
          time_dim     ! Time dimension ID

! Variable ID
      INTEGER (KIND=int_kind) :: &
          lon_id, &
          lat_id, &
          time_id

! Variable shapes
      INTEGER (KIND=int_kind), DIMENSION(1) :: &
          lon_dims, &
          lat_dims, &
          time_dims

! Longitude and latitude ranges
      REAL (KIND=real_kind), DIMENSION(2) :: &
          lon_range, lat_range
          
      CHARACTER (LEN=80) :: &
          prog
!------------------------------------------------------------------------------
      prog = 'cmorph_def_dim'
!------------------------------------------------------------------------------
! Enter define mode
      status = nf90_redef(ncid)
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))

! lon --------------------------------------------------------------------------
        status = nf90_def_dim(ncid, 'lon', lon_len, lon_dim)
        if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))

        lon_dims(1) = lon_dim
        status = nf90_def_var(ncid, 'lon', NF90_REAL, lon_dims, lon_id)
        if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
      
        status = nf90_put_att(ncid, lon_id, 'axis', 'X')
        if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
      
        status = nf90_put_att(ncid, lon_id, 'long_name', 'Longitude')
        if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
      
        status = nf90_put_att(ncid, lon_id, 'standard_name', 'longitude')
        if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
      
        status = nf90_put_att(ncid, lon_id, 'units', 'degree_east')
        if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
      
        lon_range(:) = (/ lon(1), lon(lon_len) /)
        
        status = nf90_put_att(ncid, lon_id, 'actual_range', lon_range)
        if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
    
! lat --------------------------------------------------------------------------
        status = nf90_def_dim(ncid, 'lat', lat_len, lat_dim)
        if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))

        lat_dims(1) = lat_dim
        status = nf90_def_var(ncid, 'lat', NF90_REAL, lat_dims, lat_id)
        if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
      
        status = nf90_put_att(ncid, lat_id, 'axis', 'Y')
        if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
      
        status = nf90_put_att(ncid, lat_id, 'long_name', 'Latitude')
        if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
      
        status = nf90_put_att(ncid, lat_id, 'standard_name', 'latitude')
        if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
      
        status = nf90_put_att(ncid, lat_id, 'units', 'degree_north')
        if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
      
        lat_range(:) = (/ lat(1), lat(lat_len) /)
        
        status = nf90_put_att(ncid, lat_id, 'actual_range', lat_range)
        if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
    
! time ------------------------------------------------------------------------
      status = nf90_def_dim(ncid, 'time', NF90_UNLIMITED, time_dim)
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))

      time_dims(1) = time_dim
      status = nf90_def_var(ncid, 'time', NF90_INT, time_dims, time_id)
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
      
      status = nf90_put_att(ncid, time_id, 'axis', 'T')
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
 
      status = nf90_put_att(ncid, time_id, 'long_name', 'Time')
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
      
      status = nf90_put_att(ncid, time_id, 'standard_name', 'time')
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
      
      status = nf90_put_att(ncid, time_id, 'units', 'seconds since 1970-01-01 00:00 UTC')
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
      
      status = nf90_put_att(ncid, time_id, 'calendar', 'standard')
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
      
      status = nf90_put_att(ncid, time_id, 'actual_range', time_range)
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
    
!------------------------------------------------------------------------------
! leave define mode
      status = nf90_enddef(ncid)
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))

!------------------------------------------------------------------------------
! Enter longitudes
      status = nf90_put_var(ncid, lon_id, lon)
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
! Enter latitudes
      status = nf90_put_var(ncid, lat_id, lat)
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))

!==============================================================================
      END SUBROUTINE cmorph_def_dim
!==============================================================================

!==============================================================================
      SUBROUTINE cmorph_def_var(var_nc, var_long_name, var_units, &
                                var_std_name)
!==============================================================================
!
! Define the CMORPH variables
!
!------------------------------------------------------------------------------

      CHARACTER (LEN=*), INTENT(IN) :: &
          var_nc, &          ! Variable name
          var_long_name, &   ! Variable long name
          var_units          ! Variable units
      CHARACTER (LEN=*), OPTIONAL, INTENT(IN) :: &
          var_std_name       ! Variable standard name

!------------------------------------------------------------------------------
! Local

! Dimension ID
      INTEGER (KIND=int_kind) :: &
          lon_dim, &   ! Longitude dimension ID
          lat_dim, &   ! Latitude dimension ID
          time_dim     ! Time dimension ID

! Variable ID
      INTEGER(KIND=int_kind) :: &
          var_id
! Variable shape
      INTEGER (KIND=int_kind), DIMENSION(3) :: &
          var_dims
! Missing and fill values
      REAL (KIND=real_kind), PARAMETER :: &
          missing = -9999.0, &
          fill_value = -9999.0

      CHARACTER (LEN=80) :: &
          prog
!------------------------------------------------------------------------------
      prog = 'cmorph_def_var'
!------------------------------------------------------------------------------
! Enter define mode
      status = nf90_redef(ncid)
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))

!------------------------------------------------------------------------------
! Get the dimension IDs
      status = nf90_inq_dimid(ncid, 'time', time_dim)
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
      status = nf90_inq_dimid(ncid, 'lat', lat_dim)
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
      status = nf90_inq_dimid(ncid, 'lon', lon_dim)
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
            
      var_dims(1) = lon_dim
      var_dims(2) = lat_dim
      var_dims(3) = time_dim

! define variable      
      status = nf90_def_var(ncid, var_nc, NF90_REAL, var_dims, var_id)
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))

! assign attributes
      status = nf90_put_att(ncid, var_id, 'long_name', var_long_name)
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
      
      if(present(var_std_name)) then
        status = nf90_put_att(ncid, var_id, 'standard_name', var_std_name)
        if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
      endif
        
      status = nf90_put_att(ncid, var_id, 'units', var_units)
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
      
      status = nf90_put_att(ncid, var_id, '_FillValue', fill_value)
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
      
      status = nf90_put_att(ncid, var_id, 'missing_value', missing)
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
      
!------------------------------------------------------------------------------
! leave define mode
      status = nf90_enddef(ncid)
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))

!==============================================================================
      END SUBROUTINE cmorph_def_var
!==============================================================================

!==============================================================================
      SUBROUTINE cmorph_write_var(var_nc, time, values, time_start, val_start)
!==============================================================================

! Input
      CHARACTER (LEN=*), INTENT(IN) :: &
          var_nc
      INTEGER (KIND=int_kind), INTENT(IN) :: &
          time
      REAL (KIND=real_kind), DIMENSION(:,:), INTENT(IN) :: &
          values
      INTEGER (KIND=int_kind), DIMENSION(1), INTENT(IN) :: &
          time_start
      INTEGER (KIND=int_kind), DIMENSION(3), INTENT(IN) :: &
          val_start

! Local
      INTEGER (KIND=int_kind) :: &
          time_id, &
          var_id
          
      CHARACTER (LEN=80) :: &
          prog
!------------------------------------------------------------------------------
      prog = 'cmorph_write_var'

! Inquire time ID and enter time          
      status = nf90_inq_varid(ncid, 'time', time_id)
      if (status /= NF90_NOERR) call handle_err(status, ncid, trim(prog), varid=time_id)

      status = nf90_put_var(ncid, time_id, time, start=time_start)
      if (status /= NF90_NOERR) call handle_err(status, ncid, trim(prog), varid=time_id)

! Inquire variable ID and enter values
      status = nf90_inq_varid(ncid, var_nc, var_id)
      if (status /= NF90_NOERR) call handle_err(status, ncid, trim(prog), varid=var_id)

      status = nf90_put_var(ncid, var_id, values, start=val_start)
      if (status /= NF90_NOERR) call handle_err(status, ncid, trim(prog), varid=var_id)
          
!==============================================================================
      END SUBROUTINE cmorph_write_var
!==============================================================================

!==============================================================================
      SUBROUTINE cdf_close
!==============================================================================

      CHARACTER (LEN=80) :: &
          prog
!------------------------------------------------------------------------------
      prog = 'cdf_close'

      status = nf90_close(ncid)
      if (status /= NF90_NOERR) call handle_err(status,ncid,trim(prog))
      
      RETURN
      
!==============================================================================
      END SUBROUTINE cdf_close
!==============================================================================

!==============================================================================
      SUBROUTINE handle_err(error, ncid, prog, varid)
!==============================================================================

      integer (KIND=int_kind), intent(in) :: &
          error, ncid
      character (LEN=*), intent(in) :: &
          prog
      integer (KIND=int_kind), intent(in), optional :: &
          varid

      if(present(varid)) then
        print *,'error in '// prog // ': ',trim(nf90_strerror(error)), ncid, varid
      else
        print *,'error in '// prog // ': ',trim(nf90_strerror(error)), ncid
      endif
      
      stop 'Stopped'

      RETURN

!==============================================================================
      END SUBROUTINE handle_err
!==============================================================================

END MODULE cmorph_v0_x_output_nc
