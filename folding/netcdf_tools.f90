module netcdf_tools

  USE netcdf

  implicit none

  INTRINSIC :: TRIM, ADJUSTL, NINT, SIZE

  interface read_file
    module procedure  read_file_1d
    module procedure  read_file_2d
    module procedure  read_file_3d
    module procedure  read_file_4d
  end interface

  contains

  SUBROUTINE inquire_file(fname,             &
             x_name, y_name, z_name,t_name,  &
             nx, ny, nz, nt)


    IMPLICIT NONE


    ! I/O
    CHARACTER(LEN=*), INTENT(IN)  :: fname   ! filename
    CHARACTER(LEN=*), INTENT(IN)  :: x_name  ! name of dimension in file dimension
    CHARACTER(LEN=*), INTENT(IN)  :: y_name  ! name of dimension in file dimension
    CHARACTER(LEN=*), INTENT(IN)  :: z_name  ! name of dimension in file dimension
    CHARACTER(LEN=*), INTENT(IN)  :: t_name ! name of dimension in file dimension
    INTEGER,          INTENT(OUT) :: nx  ! file dimension lenght
    INTEGER,          INTENT(OUT) :: ny  ! file dimension length
    INTEGER,          INTENT(OUT) :: nz  ! file dimension length
    INTEGER,          INTENT(OUT) :: nt  ! file dimension length

    ! LOCAL
    INTEGER                     :: status
    CHARACTER(LEN=*), PARAMETER :: substr = 'inquire_file'
    INTEGER,SAVE                :: ncid   ! netCDF-ID
    INTEGER                     :: dimid, varid
    LOGICAL                     :: file_exists  ! checking existence of file


    INTEGER, DIMENSION(:), ALLOCATABLE  :: date_file

    CHARACTER(LEN=30) :: name_dim   ! line


    INQUIRE(FILE=TRIM(fname), EXIST=file_exists)
    IF (.not.file_exists) THEN
      WRITE(*,*) 'File ',TRIM(fname),' does NOT exist! skipping....'
      STOP
    ENDIF

    ! OPEN FILE
    CALL NFERR( status, &
         nf90_open(TRIM(fname), NF90_NOWRITE, ncid) &
         ,1)
    ! latitude dimension check
    CALL  NFERR( status, &
         nf90_inq_dimid(ncid, x_name, dimid ) &
         ,2)
    CALL  NFERR( status, &
         nf90_Inquire_Dimension(ncid, dimid, name_dim, nx ) &
         ,3)
    ! longitude dimension check
    CALL  NFERR( status, &
         nf90_inq_dimid(ncid, y_name, dimid ) &
         ,4)
    CALL  NFERR( status, &
         nf90_Inquire_Dimension(ncid, dimid, name_dim, ny ) &
         ,5)
    ! vertical dimension check
    CALL  NFERR( status, &
         nf90_inq_dimid(ncid, z_name, dimid ) &
         ,6)
    CALL  NFERR( status, &
         nf90_Inquire_Dimension(ncid, dimid, name_dim, nz ) &
         ,7)
    ! time dimension check
    CALL  NFERR( status, &
         nf90_inq_dimid(ncid, t_name, dimid ) &
         ,8)
    CALL  NFERR( status, &
         nf90_Inquire_Dimension(ncid, dimid, name_dim, nt ) &
         ,9)


    !CLOSE FILE
    CALL NFERR( status, &
         nf90_close(ncid) &
         ,14)

    ! RETURN
    status = 0

  END SUBROUTINE inquire_file

  ! ------------------------------------------------------------------

  ! ------------------------------------------------------------------------
  SUBROUTINE read_file_1D(fname, varname,  data_file)


    IMPLICIT NONE

    ! I/O
    CHARACTER(LEN=*),          INTENT(IN)  :: fname       ! filename
    CHARACTER(LEN=*),          INTENT(IN)  :: varname         ! variable name
    REAL, DIMENSION(:), INTENT(OUT):: data_file ! INTENT(OUT)

    ! LOCAL
    INTEGER                     :: status
    CHARACTER(LEN=*), PARAMETER :: substr = 'read_file_1d'
    INTEGER,SAVE                :: ncid   ! netCDF-ID
    INTEGER                     :: dimid, varid

    ! OPEN FILE
    CALL NFERR( status, &
         nf90_open(TRIM(fname), NF90_NOWRITE, ncid) &
         ,21)

    CALL  NFERR( status, &
         nf90_inq_varid(ncid, TRIM(varname), varid ) &
         ,22)
    IF (status.ne.0) then
        write(*,*) "variable not found in NetCDF file, skipping"
        STOP
    ENDIF

    CALL  NFERR( status, &
         nf90_get_var(ncid, varid, data_file ) &
         ,23)

    !CLOSE FILE
    CALL NFERR( status, &
         nf90_close(ncid) &
         ,24)

    ! RETURN
    status = 0

  END SUBROUTINE read_file_1D
  ! ------------------------------------------------------------------------
  ! ------------------------------------------------------------------------
  SUBROUTINE read_file_2D(fname, varname,  data_file)


    IMPLICIT NONE

    ! I/O
    CHARACTER(LEN=*),          INTENT(IN)  :: fname       ! filename
    CHARACTER(LEN=*),          INTENT(IN)  :: varname         ! variable name
    REAL, DIMENSION(:,:), INTENT(OUT):: data_file ! INTENT(OUT)

    ! LOCAL
    INTEGER                     :: status
    CHARACTER(LEN=*), PARAMETER :: substr = 'read_file_2d'
    INTEGER,SAVE                :: ncid   ! netCDF-ID
    INTEGER                     :: dimid, varid

    ! OPEN FILE
    CALL NFERR( status, &
         nf90_open(TRIM(fname), NF90_NOWRITE, ncid) &
         ,21)

    CALL  NFERR( status, &
         nf90_inq_varid(ncid, TRIM(varname), varid ) &
         ,22)
    IF (status.ne.0) then
        write(*,*) "variable not found in NetCDF file, skipping"
        STOP
    ENDIF

    CALL  NFERR( status, &
         nf90_get_var(ncid, varid, data_file ) &
         ,23)

    !CLOSE FILE
    CALL NFERR( status, &
         nf90_close(ncid) &
         ,24)

    ! RETURN
    status = 0

  END SUBROUTINE read_file_2D
  ! ------------------------------------------------------------------------
  ! ------------------------------------------------------------------------
  SUBROUTINE read_file_3D(fname, varname,  data_file)


    IMPLICIT NONE

    ! I/O
    CHARACTER(LEN=*),          INTENT(IN)  :: fname       ! filename
    CHARACTER(LEN=*),          INTENT(IN)  :: varname         ! variable name
    REAL, DIMENSION(:,:,:), INTENT(OUT):: data_file ! INTENT(OUT)

    ! LOCAL
    INTEGER                     :: status
    CHARACTER(LEN=*), PARAMETER :: substr = 'read_file_3d'
    INTEGER,SAVE                :: ncid   ! netCDF-ID
    INTEGER                     :: dimid, varid

    ! OPEN FILE
    CALL NFERR( status, &
         nf90_open(TRIM(fname), NF90_NOWRITE, ncid) &
         ,21)

    CALL  NFERR( status, &
         nf90_inq_varid(ncid, TRIM(varname), varid ) &
         ,22)
    IF (status.ne.0) then
        write(*,*) "variable not found in NetCDF file, skipping"
        STOP
    ENDIF

    CALL  NFERR( status, &
         nf90_get_var(ncid, varid, data_file ) &
         ,23)

    !CLOSE FILE
    CALL NFERR( status, &
         nf90_close(ncid) &
         ,24)

    ! RETURN
    status = 0

  END SUBROUTINE read_file_3D
  ! ------------------------------------------------------------------------
  ! ------------------------------------------------------------------------
  SUBROUTINE read_file_4D(fname, varname,  data_file)


    IMPLICIT NONE

    ! I/O
    CHARACTER(LEN=*),          INTENT(IN)  :: fname       ! filename
    CHARACTER(LEN=*),          INTENT(IN)  :: varname         ! variable name
    REAL, DIMENSION(:,:,:,:), INTENT(OUT):: data_file ! INTENT(OUT)

    ! LOCAL
    INTEGER                     :: status
    CHARACTER(LEN=*), PARAMETER :: substr = 'read_file_4d'
    INTEGER,SAVE                :: ncid   ! netCDF-ID
    INTEGER                     :: dimid, varid

    ! OPEN FILE
    CALL NFERR( status, &
         nf90_open(TRIM(fname), NF90_NOWRITE, ncid) &
         ,21)

    CALL  NFERR( status, &
         nf90_inq_varid(ncid, TRIM(varname), varid ) &
         ,22)
    IF (status.ne.0) then
        write(*,*) "variable not found in NetCDF file, skipping"
        STOP
    ENDIF

    CALL  NFERR( status, &
         nf90_get_var(ncid, varid, data_file ) &
         ,23)

    !CLOSE FILE
    CALL NFERR( status, &
         nf90_close(ncid) &
         ,24)

    ! RETURN
    status = 0

  END SUBROUTINE read_file_4D
  ! ------------------------------------------------------------------------

  ! ------------------------------------------------------------------------
  SUBROUTINE read_att(fname, varname, attname, str)


    IMPLICIT NONE

    ! I/O
    CHARACTER(LEN=*),          INTENT(IN)  :: fname      ! filename
    CHARACTER(LEN=*),          INTENT(IN)  :: varname    ! variable name
    CHARACTER(LEN=*),          INTENT(IN)  :: attname    ! attribute name
    CHARACTER(LEN=*),          INTENT(OUT) :: str        ! date string

    ! LOCAL
    INTEGER                     :: status
    CHARACTER(LEN=*), PARAMETER :: substr = 'read_startdate'
    INTEGER,SAVE                :: ncid   ! netCDF-ID
    INTEGER                     :: dimid, varid

    ! OPEN FILE
    CALL NFERR( status, &
         nf90_open(TRIM(fname), NF90_NOWRITE, ncid) &
         ,21)

    CALL  NFERR( status, &
         nf90_inq_varid(ncid, TRIM(varname), varid ) &
         ,22)
    IF (status.ne.0) then
        write(*,*) "variable not found in NetCDF file, skipping"
        STOP
    ENDIF

    CALL  NFERR( status, &
         nf90_get_att(ncid, varid, TRIM(attname), str ) &
         ,23)

    !CLOSE FILE
    CALL NFERR( status, &
         nf90_close(ncid) &
         ,24)

    ! RETURN
    status = 0

  END SUBROUTINE read_att
  ! ------------------------------------------------------------------------
  ! BARTUSEK 20201106 ADDITION: SAME AS ABOVE BUT FOR REAL INSTEAD OF STR ATTRIBUTES IN NETCDF FILES
  ! ------------------------------------------------------------------------
  SUBROUTINE read_att_real(fname, varname, attname, str)


    IMPLICIT NONE

    ! I/O
    CHARACTER(LEN=*),          INTENT(IN)  :: fname      ! filename
    CHARACTER(LEN=*),          INTENT(IN)  :: varname    ! variable name
    CHARACTER(LEN=*),          INTENT(IN)  :: attname    ! attribute name
    REAL,         INTENT(OUT) :: str        ! date string    ! BARTUSEK 20201106 character->real ! old: CHARACTER(LEN=*),          INTENT(OUT) :: str        ! date string

    ! LOCAL
    INTEGER                     :: status
    CHARACTER(LEN=*), PARAMETER :: substr = 'read_startdate'
    INTEGER,SAVE                :: ncid   ! netCDF-ID
    INTEGER                     :: dimid, varid

    ! OPEN FILE
    CALL NFERR( status, &
         nf90_open(TRIM(fname), NF90_NOWRITE, ncid) &
         ,21)

    CALL  NFERR( status, &
         nf90_inq_varid(ncid, TRIM(varname), varid ) &
         ,22)
    IF (status.ne.0) then
        write(*,*) "variable not found in NetCDF file, skipping"
        STOP
    ENDIF

    CALL  NFERR( status, &
         nf90_get_att(ncid, varid, TRIM(attname), str ) &
         ,23)

    !CLOSE FILE
    CALL NFERR( status, &
         nf90_close(ncid) &
         ,24)

    ! RETURN
    status = 0

  END SUBROUTINE read_att_real
  ! ------------------------------------------------------------------------

  ! ------------------------------------------------------------------------
  SUBROUTINE nc_dump(fname, x_data, y_data,z_data,t_data,  &
                     x_units,y_units, z_units,t_units, aps, ak, bk, &
                     label,fold, sfold, mfold, dfold, tp, dp, pmin, pmax)

    USE netcdf

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: fname
    REAL, DIMENSION(:), INTENT(IN) :: x_data, y_data,z_data,t_data
    CHARACTER(LEN=*), INTENT(IN) :: x_units,y_units, z_units,t_units
    REAL, DIMENSION(:,:,:), INTENT(IN) :: aps
    REAL, DIMENSION(:), INTENT(IN) :: ak,bk
    REAL, DIMENSION(:,:,:,:), INTENT(IN) :: label
    REAL, DIMENSION(:,:,:), INTENT(IN) :: fold
    REAL, DIMENSION(:,:,:), INTENT(IN) :: sfold
    REAL, DIMENSION(:,:,:), INTENT(IN) :: mfold
    REAL, DIMENSION(:,:,:), INTENT(IN) :: dfold
    REAL, DIMENSION(:,:,:), INTENT(IN) :: dp
    REAL, DIMENSION(:,:,:), INTENT(IN) :: tp
    REAL, DIMENSION(:,:,:), INTENT(IN) :: pmin
    REAL, DIMENSION(:,:,:), INTENT(IN) :: pmax

    ! LOCAL
    INTEGER        :: status
    CHARACTER(LEN=*), PARAMETER :: substr = 'nc_dump'
    INTEGER :: ncid      ! netCDF-ID
    INTEGER :: nlon,nlat,nlev,ntime
    INTEGER :: dimid_lat, dimid_lon, dimid_lev, dimid_ilev, dimid_time
    INTEGER :: varid_lat, varid_lon, varid_lev, varid_ilev, varid_time
    INTEGER :: varid_aps, varid_hyam, varid_hybm, varid_label
    INTEGER :: varid_fold, varid_sfold, varid_mfold, varid_dfold, varid_tp, varid_dp
    INTEGER :: varid_pmin, varid_pmax

    !
    CHARACTER(LEN=8)         :: date
    CHARACTER(LEN=10)        :: time
    CHARACTER(LEN=5)         :: zone


     nlon  = SIZE(x_data)
     nlat  = SIZE(y_data)
     nlev  = SIZE(z_data)
     ntime = SIZE(t_data)

    ! CREATE NEW FILE
    CALL NFERR(status, &
         nf90_create(TRIM(fname), NF90_CLOBBER, ncid) &
         ,51)

    ! ADD GLOBALE ATTRIBUTES
    ! - VERSION
    CALL NFERR(status, &
         nf90_put_att(ncid, NF90_GLOBAL, 'contact:',     &
         'Andrea Pozzer, MPIC, Mainz') &
         ,52)
    ! - DATE AND TIME
    CALL DATE_AND_TIME(date, time, zone)
    CALL NFERR(status, &
         nf90_put_att(ncid, NF90_GLOBAL, 'date', date) &
         ,53)
    CALL NFERR(status, &
         nf90_put_att(ncid, NF90_GLOBAL, 'time', TRIM(time)//TRIM(zone)) &
         ,54)
    ! DEFINE DIMENSIONS
    CALL NFERR(status, &
         nf90_def_dim(ncid, 'lon', nlon, dimid_lon) &
         ,56)
    CALL NFERR(status, &
         nf90_def_dim(ncid, 'lat', nlat, dimid_lat) &
         ,57)
    CALL NFERR(status, &
         nf90_def_dim(ncid, 'lev', nlev, dimid_lev) &
         ,58)
    CALL NFERR(status, &
         nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid_time) &
         ,59)

    ! DEFINE COORDINATE VARIABLES WITH ATTRIBUTES
    CALL NFERR(status, &
         nf90_def_var(ncid, 'lon', NF90_FLOAT, (/ dimid_lon /), varid_lon) &
         ,60)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_lon, 'long_name', 'longitude') &
         ,61)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_lon, 'units', x_units) &
         ,62)

    CALL NFERR(status, &
         nf90_def_var(ncid, 'lat', NF90_FLOAT, (/ dimid_lat /), varid_lat) &
         ,63)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_lat, 'long_name', 'latitude') &
         ,63)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_lat, 'units', y_units) &
         ,64)

    CALL NFERR(status, &
         nf90_def_var(ncid, 'lev', NF90_FLOAT, (/ dimid_lev /), varid_lev) &
         ,65)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_lev, 'long_name', 'level index') &
         ,66)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_lev, 'units', z_units) &
         ,67)
    CALL NFERR(status, &
         nf90_def_var(ncid, 'time', NF90_FLOAT, (/ dimid_time /), varid_time) &
         ,68)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_time, 'long_name', 'time') &
         ,69)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_time, 'units', t_units) &
         ,70)

    ! DEFINE VARIABLES
    ! - aps
    CALL NFERR(status, &
         nf90_def_var(ncid, 'APS', NF90_FLOAT  &
         , (/ dimid_lon, dimid_lat, dimid_time /), varid_aps) &
         ,78)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_aps, 'long_name' &
         ,'Surface Pressure') &
         ,79)
    ! - ak
    CALL NFERR(status, &
         nf90_def_var(ncid, 'hyam', NF90_FLOAT  &
         , (/ dimid_lev /), varid_hyam) &
         ,78)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_hyam, 'long_name' &
         ,'"hybrid A coefficient at layer midpoints') &
         ,79)
    ! - bk
    CALL NFERR(status, &
         nf90_def_var(ncid, 'hybm', NF90_FLOAT  &
         , (/ dimid_lev /), varid_hybm) &
         ,78)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_hybm, 'long_name' &
         ,'hybrid B coefficient at layer midpoints') &
         ,79)
    ! - label
    CALL NFERR(status, &
         nf90_def_var(ncid, 'label', NF90_FLOAT  &
         , (/ dimid_lon, dimid_lat, dimid_lev, dimid_time /), varid_label) &
         ,78)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_label, 'long_name' &
         ,'TR=1,ST=2,ST.CUTOFF=3,TR.TR.CUTOFF=4,PV.BLOB=5') &
         ,79)
    ! - fold
    CALL NFERR(status, &
         nf90_def_var(ncid, 'fold', NF90_FLOAT  &
         , (/ dimid_lon, dimid_lat, dimid_time /), varid_fold) &
         ,78)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_fold, 'long_name' &
         , 'Tropopause folding: 1=YES, 2=NO') &
         ,79)
    ! - sfold
    CALL NFERR(status, &
         nf90_def_var(ncid, 'sfold', NF90_FLOAT  &
         , (/ dimid_lon, dimid_lat, dimid_time /), varid_sfold) &
         ,78)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_sfold, 'long_name' &
         , 'Shallow tropopause folding (50< >200 hPa) : 1=YES, 2=NO') &
         ,79)
    ! - mfold
    CALL NFERR(status, &
         nf90_def_var(ncid, 'mfold', NF90_FLOAT  &
         , (/ dimid_lon, dimid_lat, dimid_time /), varid_mfold) &
         ,78)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_mfold, 'long_name' &
         , 'Medium Tropopause folding (200< >350 hPa): 1=YES, 2=NO') &
         ,79)
    ! - dfold
    CALL NFERR(status, &
         nf90_def_var(ncid, 'dfold', NF90_FLOAT  &
         , (/ dimid_lon, dimid_lat, dimid_time /), varid_dfold) &
         ,78)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_dfold, 'long_name' &
         , 'Deep Tropopause folding (>=350 hPa): 1=YES, 2=NO') &
         ,79)
    ! - tropopause
    CALL NFERR(status, &
         nf90_def_var(ncid, 'tp', NF90_FLOAT  &
         , (/ dimid_lon, dimid_lat, dimid_time /), varid_tp) &
         ,78)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_tp, 'long_name' &
         , 'Tropopause') &
         ,79)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_tp, 'units' &
         , 'hPa') &
         ,79)
    ! - folding extension
    CALL NFERR(status, &
         nf90_def_var(ncid, 'dp', NF90_FLOAT  &
         , (/ dimid_lon, dimid_lat, dimid_time /), varid_dp) &
         ,78)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_dp, 'long_name' &
         , 'vertical extend of folding') &
         ,79)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_dp, 'units' &
         , 'hPa') &
         ,79)
    CALL NFERR(status, &
         nf90_def_var(ncid, 'pmin', NF90_FLOAT  &
         , (/ dimid_lon, dimid_lat, dimid_time /), varid_pmin) &
         ,78)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_pmin, 'long_name' &
         , 'vertical extend of folding') &
         ,79)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_pmin, 'units' &
         , 'hPa') &
         ,79)
    CALL NFERR(status, &
         nf90_def_var(ncid, 'pmax', NF90_FLOAT  &
         , (/ dimid_lon, dimid_lat, dimid_time /), varid_pmax) &
         ,78)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_pmax, 'long_name' &
         , 'vertical extend of folding') &
         ,79)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_pmax, 'units' &
         , 'hPa') &
         ,79)

    ! SWITCH MODUS
    CALL NFERR(status, &
         nf90_enddef(ncid) &
         ,82)

    ! SAVE COORDINATE VARIBLES
    CALL NFERR(status, &
         nf90_put_var(ncid, varid_lon, x_data)  &
         ,83)
    CALL NFERR(status, &
         nf90_put_var(ncid, varid_lat, y_data)  &
         ,84)
    CALL NFERR(status, &
         nf90_put_var(ncid, varid_lev, z_data)  &
         ,85)
    CALL NFERR(status, &
         nf90_put_var(ncid, varid_time, t_data) &
         ,86)

    CALL NFERR(status, &
         nf90_put_var(ncid, varid_aps, aps) &
         ,89)
    CALL NFERR(status, &
         nf90_put_var(ncid, varid_hyam, ak) &
         ,89)
    CALL NFERR(status, &
         nf90_put_var(ncid, varid_hybm, bk) &
         ,89)
    CALL NFERR(status, &
         nf90_put_var(ncid, varid_label, label) &
         ,89)
    CALL NFERR(status, &
         nf90_put_var(ncid, varid_fold, fold) &
         ,89)
    CALL NFERR(status, &
         nf90_put_var(ncid, varid_sfold, sfold) &
         ,89)
    CALL NFERR(status, &
         nf90_put_var(ncid, varid_mfold, mfold) &
         ,89)
    CALL NFERR(status, &
         nf90_put_var(ncid, varid_dfold, dfold) &
         ,89)
    CALL NFERR(status, &
         nf90_put_var(ncid, varid_tp, tp) &
         ,89)
    CALL NFERR(status, &
         nf90_put_var(ncid, varid_dp, dp) &
         ,89)
    CALL NFERR(status, &
         nf90_put_var(ncid, varid_pmin, pmin) &
         ,89)
    CALL NFERR(status, &
         nf90_put_var(ncid, varid_pmax, pmax) &
         ,89)

    ! CLOSE FILE
    CALL NFERR(status, &
         nf90_close(ncid) &
         ,90)

  END SUBROUTINE nc_dump


  ! ------------------------------------------------------------------
  SUBROUTINE NFERR(status,command,pos)

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    INTEGER,          INTENT(IN) :: command
    INTEGER,          INTENT(IN) :: pos

    status=command
    IF (status /= NF90_NOERR) THEN
       WRITE(*,*) 'netCDF ERROR at position: ', pos
       WRITE(*,*) 'netCDF ERROR status     : ',status
       WRITE(*,*) 'netCDF ERROR            : ',nf90_strerror(status)
    END IF

  END SUBROUTINE NFERR
  ! ------------------------------------------------------------------



end module netcdf_tools
