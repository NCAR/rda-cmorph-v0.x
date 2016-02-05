!------------------------------------------------------------------
!
!                 CMORPH v0.x 0.25 degree 3-hourly Information
!
! Last Update: 03 December 2014 by Thomas Cram (tcram@ucar.edu)
!
! The 0.25 degree 3-hourly CMORPH can be found on the NCEP ftp server at 
! ftp://ftp.cpc.ncep.noaa.gov/precip/CMORPH_V0.x/RAW/0.25deg-3HLY/
!
!
! Each file is composed of 8 direct access binary ("little_endian") 
! records that are defined as follows:
!
! Record 1: contains the "CMORPH" precipitation estimates for 00 UTC
! Record 2: contains the "CMORPH" precipitation estimates for 03 UTC
! Record 3: contains the "CMORPH" precipitation estimates for 06 UTC
! Record 4: contains the "CMORPH" precipitation estimates for 09 UTC
! Record 5: contains the "CMORPH" precipitation estimates for 12 UTC
! Record 6: contains the "CMORPH" precipitation estimates for 15 UTC
! Record 7: contains the "CMORPH" precipitation estimates for 18 UTC
! Record 8: contains the "CMORPH" precipitation estimates for 21 UTC
!
! All units are "mm/3hr".  Missing data are denoted by values of "-9999."
! Each record contains a 1440 x 480 REAL*4 array of data which is oriented
! from 0.125E EASTward and from -59.875N NORTHward, with a grid increment
! of 0.25 degrees of latitude and longitude.  Thus, the grid locations
! are the centers on a 0.25 degree lat/lon grid box.  Note that these 
! estimates represent spatial averages, so the data are grid-centered, 
! rather than lattice-centered.
!
! For example (1,1) is 0.125E, -59.875N 
!             (2,2) is 0.375E, -59.625N, etc.
!
! Associated GrADS ctl file with this data set is located at:
! ftp://ftp.cpc.ncep.noaa.gov/precip/CMORPH_V0.x/CTL/
! 
! Below is FORTRAN code that will read the data sets:
!
!------------------------------------------------------------------

      PROGRAM cmorph_v0_x_subset
      
      USE kinds
      USE cmorph_v0_x_output_nc
      
      IMPLICIT NONE

      INTEGER (KIND=int_kind), PARAMETER :: &
        IRMAX=1440, &        ! Meridional length of input data array
        JRMAX=480, &         ! Latitudinal length of input data array
        PMAX=1, &
        RSIZE=4
      REAL (KIND=real_kind), PARAMETER :: &
        lon1loc = 0.125, &   ! Latitude location of first grid cell (degrees)
        lon2loc = 359.875, & ! Latitude location of last grid cell (degrees)
        lat1loc = -59.875, & ! Longitude location of first grid cell (degrees)
        lat2loc = 59.875, &  ! Longitude location of last grid cell (degrees)
        resol = 0.25, &      ! Latitude/longitude grid resolution (degrees)
        missing = -9999.     ! Missing value
      INTEGER (KIND=int_kind) :: &
        irec, &              ! Record counter for input data
        outrec, &            ! Record counter for output data
        io, &                ! IOSTAT for configuration file
        starthr, &           ! Start time for temporal subsetting (00Z, 06Z, 12Z, 18Z)
        endhr, &             ! End time for temporal subsetting (00Z, 06Z, 12Z, 18Z)
        time_len, &          ! Number of time records in output
        epoch, &             ! UNIX epoch time stamp
        lon1_ix, lon2_ix, &  ! Indices of longitude subset boundaries
        lat1_ix, lat2_ix, &  ! Indices of latitude subset boundaries
        IRSUB, JRSUB, &      ! Length of lat/lon subset dimensions
        inlength, ihr, &
        i, j, ireg, jreg, &
        timeOUT
      INTEGER (KIND=int_kind), DIMENSION(1) :: &
        time_start
      INTEGER (KIND=int_kind), DIMENSION(3) :: &
        val_start
      INTEGER (KIND=int_kind), DIMENSION(2) :: &
        time_range
      REAL (KIND=real_kind) :: &
        lon1, lon2, lat1, lat2
      REAL (KIND=real_kind), DIMENSION(JRMAX) :: &
        lat_grid             ! Latitude grid centers (degrees)
      REAL (KIND=real_kind), DIMENSION(IRMAX) :: &
        lon_grid             ! Longitude grid centers (degrees)
      REAL (KIND=real_kind), ALLOCATABLE, DIMENSION(:) :: &
        lon_sub, &           ! Longitude subset grid
        lat_sub              ! Latitude subset grid
      REAL*4, DIMENSION(IRMAX, JRMAX, PMAX) :: &
        ssmipropREG          ! Input data array
      REAL (KIND=real_kind), ALLOCATABLE, DIMENSION(:,:) :: &
        cmorphOUT            ! CMORPH precipitation output array
      CHARACTER (LEN=120) :: &
        argv, &              ! Input argument array
        infile, &            ! Input data file (path + file name)
        outfile, &           ! Output data file (path + file name)
        config               ! Configuration file (path + file name)
      CHARACTER (LEN=101) :: &
        crec
      CHARACTER (LEN=20) :: &
        trec
      CHARACTER (LEN=4) :: &
        id
      CHARACTER (LEN=6) :: &
        format
      INTEGER (KIND=int_kind), ALLOCATABLE, DIMENSION(:) :: &
        timestamp
      LOGICAL (KIND=log_kind) :: &
        isfile               ! Logical to check for existence of input/output files
!------------------------------------------------------------------

      irec = 0
      outrec = 0

!------------------------------------------------------------------
!      Input arguments:
!      1. Input data file (path + file name)
!      2. Output data file (path + file name)
!      3. Configuration file (path + file name)

      call getarg(1, argv)
      infile = argv
      call getarg(2, argv)
      outfile = argv
      call getarg(3, argv)
      config = argv

!------------------------------------------------------------------
! Default parameters for subsetting
      starthr = 0
      endhr   = 21
      time_len = 8
      lon1    = lon1loc
      lon2    = lon2loc
      lat1    = lat1loc
      lat2    = lat2loc

!------------------------------------------------------------------
!	open the configuration file and read in parameters
      open (unit=10, file=config, form='formatted')
      do i = 1,10
        read(10, '(A100)', iostat=io) crec
        if (io < 0) exit
        inlength=len_trim(crec)-4
        select case (crec(1:4))
          case ("LATS")
            read (crec,*) id,lat1,lat2
            print *,"Latitude Values: ",lat1,lat2
          case ("LONS")
            read (crec,*) id,lon1,lon2
            print *,"Longitude Values: ",lon1,lon2
          case ("HOUR")
            read (crec,*) id,starthr,endhr
            print *,"Start and end hours: ",starthr,endhr
            time_len = (endhr-starthr)/3 + 1
            ALLOCATE(timestamp(time_len))
          case ("FORM")
            read (crec,*) id,format
            print *,"Output format: ",format
          case ("TIME")
            read (crec,*) id
            print *,id
            do j = 1, time_len
              read(10, '(A20)', iostat=io) trec
              if (io < 0) exit
              read(trec,*) epoch
              timestamp(j) = epoch
              print *,"Unix epoch time: ",j,epoch
            enddo
        end select
      end do
      close(10)

      time_range(:) = (/ timestamp(1), timestamp(time_len) /)

!------------------------------------------------------------------
! Define the CMORPH latitude and longitude grids (0.25 deg resolution)
! and determine the array indices of the lat/lon grid subset boundary

      lon1_ix = 1
      lon2_ix = IRMAX
      lat1_ix = 1
      lat2_ix = JRMAX
      
      do i = 1, IRMAX
        lon_grid(i) = lon1loc + (i-1)*resol
      enddo
      do j = 1, JRMAX
        lat_grid(j) = lat1loc + (j-1)*resol
      enddo

      lon1_ix = MINLOC(lon_grid, DIM = 1, MASK = lon_grid .GE. lon1)
      lon2_ix = MAXLOC(lon_grid, DIM = 1, MASK = lon_grid .LE. lon2)
      lat1_ix = MINLOC(lat_grid, DIM = 1, MASK = lat_grid .GE. lat1)
      lat2_ix = MAXLOC(lat_grid, DIM = 1, MASK = lat_grid .LE. lat2)
      
      IRSUB = lon2_ix - lon1_ix + 1
      JRSUB = lat2_ix - lat1_ix + 1
      
! Allocate memory for lon/lat subgrids and output array
      ALLOCATE (lon_sub(IRSUB))
      ALLOCATE (lat_sub(JRSUB))
      ALLOCATE (cmorphOUT(IRSUB, JRSUB))
      
      lon_sub(:) = lon_grid(lon1_ix:lon2_ix)
      lat_sub(:) = lat_grid(lat1_ix:lat2_ix)
      
!------------------------------------------------------------------
!    Open input regional blended sensor rainfall file for 3
!      hour CPC micro rainfall file

      write(6,'('' Opening input file '',a120)') infile

!    Check if input file exists
      inquire (file=infile, exist=isfile)
      if (isfile == .false.) then
	    write(6,*) 'cmorph_subset: Input file does not exist. Exiting program.'
		call exit()
      endif

!    Open direct access binary input data file (little endian).
      open (unit=15, file=infile, access='direct', status='old', &
            recl=IRMAX*JRMAX*PMAX*RSIZE, form='unformatted', &
            action='read')

!------------------------------------------------------------------
!    Open output data file (binary or NetCDF)
      
      if (format .eq. 'BINARY') then
        open (unit=95, file=outfile, access = 'direct', &
              status='unknown', action='write', &
              form='unformatted', recl=IRSUB*JRSUB*RSIZE)
      else
        call cmorph_create_cdf(outfile)
        call cmorph_def_dim(IRSUB, JRSUB, lon_sub, lat_sub, time_range)
        call cmorph_def_var('cmorph_precip', &
                            'CMORPH precipitation estimate', &
                            'mm (3 hr)^-1')
      endif
      
!------------------------------------------------------------------
!    Loop over 8 3 hourly periods for one day
      do ihr = 0, 21, 3
        irec = irec + 1
        write(6,'(''reading 3-hourly data at hour '',i4)') ihr
        read (15, rec = irec) ssmipropREG

!     Check data
!       write (6,'(/,'' cpc combined MW '')')
!       do jreg = 271, 280
!        write(6,'(15(f12.1))') (ssmipropREG(ireg,jreg,1),ireg=1161,1175)
!       enddo
!       write (6,'(/,'' CMORPH '')')
!       do jreg = 271, 280
!        write(6,'(15(f12.1))') (ssmipropREG(ireg,jreg,2),ireg=1161,1175)
!       enddo

! Skip if outside the requested time period
        if ((ihr .lt. starthr) .or. &
            (ihr .gt. endhr)) cycle

! Write subset data        
        do i = lon1_ix, lon2_ix
        do j = lat1_ix, lat2_ix
          cmorphOUT(i-lon1_ix+1,j-lat1_ix+1) = ssmipropREG(i,j,1)
        enddo
        enddo

        write(6,'(''writing data for hour '',i4)') ihr
        outrec = outrec + 1
        time_start(1) = outrec
        val_start(:) = (/ 1, 1, outrec /)
        timeOUT = timestamp(outrec)

        if (format .eq. 'BINARY') then
          write(95, rec=outrec) cmorphOUT
        else
          call cmorph_write_var('cmorph_precip', timeOUT, cmorphOUT, time_start, val_start)
        endif
        
      enddo  ! end ihr loop
!------------------------------------------------------------------
      
      close (15)
      
      if (format .eq. 'BINARY') close (95)
      if (format .eq. 'NETCDF') call cdf_close()
      
      END PROGRAM cmorph_v0_x_subset
