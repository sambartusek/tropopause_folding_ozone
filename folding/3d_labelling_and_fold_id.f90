PROGRAM clustering

  ! *********************************************************************************
  ! * Interpolate fields to tropopause and assign labels with clustering algorithm  *
  ! * Michael Sprenger / Spring 2006: initial version                               *
  ! * Michael Sprenger / 2012: horizontal propagation of LABEL 2 prohibited in      *
  ! *                          lower troposphere                                    *
  ! * Bojan Skerlak / Spring 2012 / Fall 2012: advanced LABEL 1-5 handling          *
  ! *  Jan 2013:    - change horiz. prop. decision method to height(log(pressure))  *
  ! *               - perform check of label for exotic cases where the pressure    *
  ! *                 is higher than 800 hPa (deep folds)                           *
  ! *               - additional cutoff-special-treatment prohibits identification  *
  ! *                 of strat. air as cutoffs if below horizontal prop. limit      *
  ! *                 (e.g. deep 'funnels' with 'bumps' on the side)                *
  ! *               - additional check in strat. (TH>380) for 'wrong' sign of PV    *
  ! *                 that can lead to errors in calculating TROPO                  *
  ! *               - additional check for PV towers (whole column gets label 2)    *
  ! * July 2013:    - replace fixed tropo_pv,tropo_th by user-defined values        *
  ! * Jan 2014:     - v2.0 Add specific humidity criterion for diabatic PV          *
  ! * Feb 2014:     - v2.1 Correct 'filling' when labels 2/3 or 2/5 meet            *
  ! * Mar 2014:     - v2.2 Filling of 2/3 and 2/5 also propagates downward via      *
  ! *                 'connect' subroutine. Thus, now also dird/diru have to be     *
  ! *                 specified.                                                    *
  ! *               - v2.3 Identify and write out tropopause folds (FOLD) 1/0       *
  ! * Nov 2014:     - remove some commented sections and delete temporary stuff     *
  ! * Dec 2014:     - simply the code to only contain the core algorithm and no     *
  ! *                 input/output section (which depends on data format and        *
  ! *                 model setup)                                                  *
  ! * Andrea Pozzer (MPIC):                                                         *
  ! * Dec 2014:    - Implementation for netcdf I/O, Makefile and namelist reading   *
  ! * Samuel Bartusek:                                                              *
  ! * 2021-2022:   - Adjustments to NetCDF I/O, Tropofold routine, and labelling    *
  ! *                ("Tropopause") routine, for application to ERA5. All changes   *
  ! *                labelled with "BARTUSEK", and labelling routine edits detailed *
  ! *                in Bartusek et al., 2023 Supplementary Information.            *
  ! *********************************************************************************

  !????????????????????????????????????????????????????????????????????????????????
  !????????????????????????????????????????????????????????????????????????????????
  !????????????????????????????????????????????????????????????????????????????????
  !???????????????                                       ??????????????????????????
  !???????????????   1) 3D-LABELLING                     ??????????????????????????
  !???????????????   2) TROPOPAUSE FOLD IDENTIFICATION   ??????????????????????????
  !???????????????                                       ??????????????????????????
  !????????????????????????????????????????????????????????????????????????????????
  !????????????????????????????????????????????????????????????????????????????????
  !MMM?????????????????????????????????????????????????????????????????????????????
  !.....DMMI???????????????????????????????????????????MMMMMMMMMMM?I???????????????
  !.........IMO????????????????????????????????????IMM ...........  8MM8???????????
  !............MM?????????????????????????????????M$       ..............MMMI??????
  !............. MD?????????????????????????????IM.   PMIN .................OMMI??
  !................MI???????????????????????????M . ^      ......................NM
  !.................MO?????????????????????????M... | .............................
  !.................. M???????????????????????M.... | .............................
  !....................M8?????????????????????M.... |     .........................
  !.....................:M????????????????????M ... |  DP .........................
  !.......................MO???????????????????M... |     .........................
  !........................ M??????????????????M~.. | .............................
  !..........................MM????????????????IM.. | .............................
  !............................MI????????????????M. v .............................
  !.............................:MI???????????????M.  .............................
  !...............................DMI??????????????M$ .............................
  !................................ MM??????????????DM ............................
  !.................................. MMI????????????IM ...........................
  !.....................................MMI????????????M ..........................
  !........................................MM??????????IM..........................
  !......................................... $MMI???????M .........................
  !........................................... MMMMMMM7. ..........................
  !.............................................      .............................
  !............................................. PMAX .............................
  !.............................................      .............................
  !................................................................................

   USE mo_f2kcli                    ! command line interface
   USE netcdf_tools

  implicit none

  ! VERSION
  character(len=*), parameter :: VERSION = '1.4'
  ! ... INTERNAL parameter
  integer, parameter :: str_short = 30  ! length of short strings
  integer, parameter :: str_long  = 100 ! length of long strings
  integer, parameter :: str_vlong = 500 ! length of very long strings
  integer, parameter :: iou  = 21       ! I/O unit

  ! FOR COMMAND LINE
  CHARACTER(LEN=256) :: EXE          ! program name
  CHARACTER(LEN=80)  :: CMD          ! argument
  INTEGER            :: NARG         ! number of arguments

  ! VARIABLES
  INTEGER                         :: status       ! status flag
  LOGICAL                         :: l_invert = .FALSE.   ! invert axis
  LOGICAL                         :: l_Pa = .FALSE.   ! convert units Pa->hPa

  ! Set tropopause thresholds
  real :: tropo_pv,tropo_th
  ! diabatically produced PV (Q-threshold), Feb 2014:
  real :: q_th
  ! masking number
  real, parameter :: mdv =-999.0

  ! DATA/VARIABLE
  ! - input data
  real, allocatable, dimension(:) :: x_data, y_data,z_data,t_data
  real, allocatable, dimension(:) :: ak,bk
  real, allocatable, dimension(:,:,:) :: aps
  real, allocatable, dimension(:,:,:,:) :: press
  real, allocatable, dimension(:,:,:,:) :: q_in,pv_in,pt_in
  real, allocatable, dimension(:,:,:,:) :: q,pv,pt
  character(len=str_long) :: x_units,y_units, z_units,t_units
  character(len=str_long) :: aps_units
  ! - final results
  real, allocatable, dimension(:,:,:) :: fold    ! fold yes/no
  real, allocatable, dimension(:,:,:) :: sfold   ! shallow fold yes/no
  real, allocatable, dimension(:,:,:) :: mfold   ! medium fold yes/no
  real, allocatable, dimension(:,:,:) :: dfold   ! deep fold yes/no
  real, allocatable, dimension(:,:,:) :: tp      ! tropopause
  real, allocatable, dimension(:,:,:) :: dp      ! folding depth (Pa)
  real, allocatable, dimension(:,:,:) :: pmin    ! folding depth (Pa)
  real, allocatable, dimension(:,:,:) :: pmax    ! folding depth (Pa)
  real, allocatable, dimension(:,:,:,:) :: label ! folding classification
  real, allocatable, dimension(:,:,:,:) :: label_out

  ! Auxiliary variables
  integer :: stat
  integer :: i,j,k,t
  real :: lev(100)


  !namelist...
  character(len=str_long) file_input
  character(len=str_long) file_output
  character(len=str_short) x_name,y_name,z_name,t_name
  character(len=str_short) aps_name,hyam_name, hybm_name
  character(len=str_short) q_name, pv_name, pt_name

  ! DATA DIMENSION
  integer :: nx,ny,nz,nt

  ! ---------------
  ! BARTUSEK / ADDITION FOR UNPACKING PACKED FILES (USING SCALE_FACTOR AND ADD_OFFSET)
  ! ---------------
  real, allocatable :: aps_sf, ak_sf, bk_sf, q_in_sf, pt_in_sf, pv_in_sf, aps_ao, ak_ao, bk_ao, q_in_ao, pt_in_ao, pv_in_ao
  ! ---------------
  ! BARTUSEK / END
  ! ---------------

  ! -----------------------------------------------------------------
  ! Initialisation paramaters
  ! -----------------------------------------------------------------
  ! Humidity criterion for diabatically produced PV, BOJAN: this is the threshold Q=0.1g/kg
  q_th=0.0001
  tropo_pv = 2.0
  tropo_th = 380.0

  ! -----------------------------------------------------------------
  ! Read command line
  ! -----------------------------------------------------------------
  narg = command_argument_count()    ! number of arguments
  call get_command_argument(0,exe)   ! program name
  !
  if (narg > 1) then
     write(*,*) 'command-line error: too many arguments !'
     call usage(trim(exe))
     stop
  end if
  !
  if (narg == 0) then
     call usage(trim(exe))
     stop
  end if
  !
  call get_command_argument(1,cmd)

  ! -----------------------------------------------------------------
  ! Read namelist file
  ! -----------------------------------------------------------------

  CALL read_nml(status, iou, TRIM(CMD))
  IF (status /= 0) STOP

  ! -----------------------------------------------------------------
  ! Inquire netcdf file
  ! -----------------------------------------------------------------
  CALL inquire_file(TRIM(file_input),     &
             x_name, y_name, z_name,t_name,  &
             nx, ny, nz, nt )
!  IF (status > 0) STOP   ! ERROR


  WRITE(*,*) '==========================================================='
  WRITE(*,*) " FILE DIMENSION: "
  WRITE(*,*) " X  :", nx
  WRITE(*,*) " Y  :", ny
  WRITE(*,*) " Z  :", nz
  WRITE(*,*) " T  :", nt
  WRITE(*,*) '==========================================================='

  ! -----------------------------------------------------------------
  ! allocate data file
  ! -----------------------------------------------------------------
  ALLOCATE(x_data(nx))
  ALLOCATE(y_data(ny))
  ALLOCATE(z_data(nz))
  ALLOCATE(t_data(nt))

  ALLOCATE(aps(nx,ny,nt))
  ALLOCATE(ak(nz))
  ALLOCATE(bk(nz))
  ALLOCATE(press(nx,ny,nz,nt))
  !
  ALLOCATE(q(nx,ny,nz,nt))
  ALLOCATE(pt(nx,ny,nz,nt))
  ALLOCATE(pv(nx,ny,nz,nt))
  ALLOCATE(q_in(nx,ny,nz,nt))
  ALLOCATE(pt_in(nx,ny,nz,nt))
  ALLOCATE(pv_in(nx,ny,nz,nt))

  ALLOCATE(label(nx,ny,nz,nt))
  ALLOCATE(fold(nx,ny,nt))
  ALLOCATE(sfold(nx,ny,nt))
  ALLOCATE(mfold(nx,ny,nt))
  ALLOCATE(dfold(nx,ny,nt))
  ALLOCATE(tp(nx,ny,nt))
  ALLOCATE(dp(nx,ny,nt))
  ALLOCATE(pmin(nx,ny,nt))
  ALLOCATE(pmax(nx,ny,nt))
  ALLOCATE(label_out(nx,ny,nz,nt))

  ! -----------------------------------------------------------------
  ! read netcdf file
  ! -----------------------------------------------------------------
  call read_file(TRIM(file_input), TRIM(x_name), x_data)
  call read_file(TRIM(file_input), TRIM(y_name), y_data)
  call read_file(TRIM(file_input), TRIM(z_name), z_data)
  call read_file(TRIM(file_input), TRIM(t_name), t_data)

  call read_att(TRIM(file_input), TRIM(x_name), "units", x_units)
  call read_att(TRIM(file_input), TRIM(y_name), "units", y_units)
  call read_att(TRIM(file_input), TRIM(z_name), "units", z_units)
  call read_att(TRIM(file_input), TRIM(t_name), "units", t_units)

  ! pressure data
  call read_file(TRIM(file_input), TRIM(aps_name), aps)
  call read_file(TRIM(file_input), TRIM(hyam_name), ak)
  call read_file(TRIM(file_input), TRIM(hybm_name), bk)
  ! presure units
  call read_att(TRIM(file_input), TRIM(aps_name), "units", aps_units)

  ! humidity
  call read_file(TRIM(file_input), TRIM(q_name), q_in)
  ! potential temperature
  call read_file(TRIM(file_input), TRIM(pt_name), pt_in)
  ! potential vorticity
  call read_file(TRIM(file_input), TRIM(pv_name), pv_in)

  ! ---------------
  ! BARTUSEK / ADDITION FOR UNPACKED PACKED FILES (USING SCALE_FACTOR AND ADD_OFFSET)
  ! ---------------
  ALLOCATE(aps_sf)
  ALLOCATE(ak_sf)
  ALLOCATE(bk_sf)
  ALLOCATE(q_in_sf)
  ALLOCATE(pt_in_sf)
  ALLOCATE(pv_in_sf)

  ALLOCATE(aps_ao)
  ALLOCATE(ak_ao)
  ALLOCATE(bk_ao)
  ALLOCATE(q_in_ao)
  ALLOCATE(pt_in_ao)
  ALLOCATE(pv_in_ao)

  ! get scale_factors
  call read_att_real(TRIM(file_input), TRIM(aps_name), "scale_factor", aps_sf)
  call read_att_real(TRIM(file_input), TRIM(hyam_name), "scale_factor", ak_sf)
  call read_att_real(TRIM(file_input), TRIM(hybm_name), "scale_factor", bk_sf)
  call read_att_real(TRIM(file_input), TRIM(q_name), "scale_factor", q_in_sf)
  call read_att_real(TRIM(file_input), TRIM(pt_name), "scale_factor", pt_in_sf)
  call read_att_real(TRIM(file_input), TRIM(pv_name), "scale_factor", pv_in_sf)

  ! get add_offsets
  call read_att_real(TRIM(file_input), TRIM(aps_name), "add_offset", aps_ao)
  call read_att_real(TRIM(file_input), TRIM(hyam_name), "add_offset", ak_ao)
  call read_att_real(TRIM(file_input), TRIM(hybm_name), "add_offset", bk_ao)
  call read_att_real(TRIM(file_input), TRIM(q_name), "add_offset", q_in_ao)
  call read_att_real(TRIM(file_input), TRIM(pt_name), "add_offset", pt_in_ao)
  call read_att_real(TRIM(file_input), TRIM(pv_name), "add_offset", pv_in_ao)

  ! multiply by sf, add ao
  aps = aps*aps_sf + aps_ao
  ak = ak*ak_sf + ak_ao
  ak = ak*100000   !!! BECAUSE BELOW (CREATING PRESS FIELD), PRESS = AK + APS*BK
  bk = bk*bk_sf + bk_ao
  q_in = q_in*q_in_sf + q_in_ao
  pt_in = pt_in*pt_in_sf + pt_in_ao
  pv_in = pv_in*pv_in_sf + pv_in_ao
  ! ---------------
  ! BARTUSEK / END
  ! ---------------

  ! -----------------------------------------------------------------
  ! create pressure field and
  ! invert vertical axis (if needed)
  ! -----------------------------------------------------------------
  ! bottom => high pressure = 1
  ! top => low pressure = nz
  if ((ak(1)+bk(1)*101325).lt.(ak(2)+bk(2)*101325)) then
    l_invert=.TRUE.
    write (*,*) 'Invert vertical axis'
  endif

  if (aps_units == "Pa") then
    write (*,*) 'convert pressure from Pa to hPa'
    l_Pa=.TRUE.
  endif

  if (l_invert) then
     ! to be inverted
     do i=1,nx
        do j=1,ny
          do k=1,nz
            do t=1,nt
            if (l_Pa) then
               press(i,j,nz-k+1,t)= (ak(k)+aps(i,j,t)*bk(k))/100.0 ! Pa -> hPa
            else
               press(i,j,nz-k+1,t)= (ak(k)+aps(i,j,t)*bk(k))
            endif
            q(i,j,nz-k+1,t) = q_in(i,j,k,t)
            pt(i,j,nz-k+1,t)= pt_in(i,j,k,t)
            pv(i,j,nz-k+1,t)= pv_in(i,j,k,t)
            enddo
          enddo
        enddo
     enddo
  ELSE
     ! NOT to be inverted
     do i=1,nx
        do j=1,ny
          do k=1,nz
            do t=1,nt
            if (l_Pa) then
               press(i,j,k,t)= (ak(k)+aps(i,j,t)*bk(k))/100.0 ! Pa -> hPa
            else
               press(i,j,k,t)= ak(k)+aps(i,j,t)*bk(k)
            endif
            enddo
          enddo
        enddo
     enddo
     q(1:nx,1:ny,1:nz,1:nt)  = q_in(1:nx,1:ny,1:nz,1:nt)
     pt(1:nx,1:ny,1:nz,1:nt) = pt_in(1:nx,1:ny,1:nz,1:nt)
     pv(1:nx,1:ny,1:nz,1:nt) = pv_in(1:nx,1:ny,1:nz,1:nt)
  ENDIF

  ! -----------------------------------------------------------------
  ! perform calculation for each timestep
  ! -----------------------------------------------------------------
  ! TIME LOOP
  do t=1,nt


     write(*,*) "time step",t,"/", nt
     ! -------------------------------------------------------------------
     ! Part 1) Run clustering algorithm, get labels and interpolate to tropopause
     ! -------------------------------------------------------------------

     ! ------------------------------
     ! Running clustering algorithm
     ! ------------------------------
     call tropopause (tp(:,:,t),label(:,:,:,t),press(:,:,:,t),pv(:,:,:,t),pt(:,:,:,t),q(:,:,:,t), &
       mdv,y_data(:),x_data(:),nx,ny,nz,tropo_pv,tropo_th,q_th)                ! BARTUSEK / added x_data(:)

     ! -------------------------------------------------------------------
     ! Part 2) Identify tropopause folds from vertical cross-sections
     ! -----------------------------------------------------------------

     ! ------------------------------
     ! Identifying tropopause folds
     ! ------------------------------
     call tropofold (label(:,:,:,t),press(:,:,:,t),pv(:,:,:,t),pt(:,:,:,t),  &
          fold(:,:,t), dp(:,:,t), pmin(:,:,t), pmax(:,:,t), sfold(:,:,t), mfold(:,:,t), dfold(:,:,t),  &
          mdv,nx,ny,nz,tropo_pv,tropo_th)

  enddo
  ! END OF TIME LOOP

  ! -----------------------------------------------------------------
  ! re-invert vertical axis (if needed)
  ! -----------------------------------------------------------------
  if (l_invert) then
     ! to be inverted
     do i=1,nx
        do j=1,ny
          do k=1,nz
            do t=1,nt
            label_out(i,j,nz-k+1,t) = label(i,j,k,t)
            enddo
          enddo
        enddo
     enddo
  ELSE
     ! NOT to be inverted
     label_out(1:nx,1:ny,1:nz,1:nt) = label(1:nx,1:ny,1:nz,1:nt)
  ENDIF

  ! -----------------------------------------------------------------
  ! write netcdf files
  ! -----------------------------------------------------------------
  call nc_dump (TRIM(file_output), x_data, y_data, z_data, t_data,             &
               x_units,y_units, z_units,t_units, aps, ak, bk, label_out, fold, &
               sfold, mfold, dfold, tp, dp, pmin, pmax)


  ! -----------------------------------------------------------------
  ! deallocate data file
  ! -----------------------------------------------------------------
  DEALLOCATE(x_data)
  DEALLOCATE(y_data)
  DEALLOCATE(z_data)
  DEALLOCATE(t_data)

  DEALLOCATE(aps)
  DEALLOCATE(ak)
  DEALLOCATE(bk)
  DEALLOCATE(press)
  !
  DEALLOCATE(q_in)
  DEALLOCATE(pt_in)
  DEALLOCATE(pv_in)
  DEALLOCATE(q)
  DEALLOCATE(pt)
  DEALLOCATE(pv)

  DEALLOCATE(label)
  DEALLOCATE(label_out)
  DEALLOCATE(fold)
  DEALLOCATE(sfold)
  DEALLOCATE(mfold)
  DEALLOCATE(dfold)
  DEALLOCATE(tp)
  DEALLOCATE(dp)
  DEALLOCATE(pmin)
  DEALLOCATE(pmax)

CONTAINS

  ! ------------------------------------------------------------------------
  SUBROUTINE USAGE(EXE)
    CHARACTER (LEN=*) :: EXE
    WRITE(*,*) '--------------------------------------------'
    WRITE(*,*) '3d_labelling_and_fold_id, version: ',VERSION
    WRITE(*,*) 'Authors:  Andrea Pozzer, MPICH  '
    WRITE(*,*) '          Bojan Skerlak, ETH    '
    WRITE(*,*) '          Michael Sprenger, ETH '
    WRITE(*,*) '--------------------------------------------'
    WRITE(*,*) 'Usage: '//TRIM(EXE)//' <namelist-file>'
    WRITE(*,*) 'See README.txt or EMAC.nml for example'
    WRITE(*,*) '--------------------------------------------'
  END SUBROUTINE USAGE
  ! ------------------------------------------------------------------------

 ! ------------------------------------------------------------------------
  SUBROUTINE read_nml(status, iou, fname)

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    INTEGER,          INTENT(IN)  :: iou
    CHARACTER(LEN=*), INTENT(IN)  :: fname

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'read_nml'
    LOGICAL :: lex   ! file exists
    INTEGER :: fstat ! file status

    NAMELIST /CTRL/ file_input,file_output,X_name,Y_name,Z_name,T_name, &
                    APS_name,HYAM_name,HYBM_name,Q_name,PV_name,PT_name

    status = 1 ! ERROR

    WRITE(*,*) '==========================================================='

    ! CHECK IF FILE EXISTS
    INQUIRE(file=TRIM(fname), exist=lex)
    IF (.NOT.lex) THEN
       WRITE(*,*) substr,': FILE DOES NOT EXIST (',TRIM(fname),')'
       status = 1
       RETURN
    END IF

    ! OPEN FILE
    OPEN(iou,file=TRIM(fname))
   ! READ NEMELIST
    WRITE(*,*) 'READING NAMELIST ''CTRL'''//&
         &' FROM '''//TRIM(fname),''' (unit ',iou,') ...'
    !
    READ(iou, NML=CTRL, IOSTAT=fstat)
    !
    IF (fstat /= 0) THEN
       WRITE(*,*) substr,': READ ERROR IN NAMELIST ''CTRL'' (',TRIM(fname),')'
       status = 3  ! READ ERROR IN NAMELIST
       RETURN
    END IF

    WRITE(*,*) ' FILE INPUT  : ', TRIM(file_input)
    WRITE(*,*) ' FILE OUTPUT : ', TRIM(file_output)
    WRITE(*,*) ' X_NAME      : ', TRIM(x_name)
    WRITE(*,*) ' Y_NAME      : ', TRIM(y_name)
    WRITE(*,*) ' Z_NAME      : ', TRIM(z_name)
    WRITE(*,*) ' T_NAME      : ', TRIM(t_name)
    WRITE(*,*) ' APS_NAME    : ', TRIM(aps_name)
    WRITE(*,*) ' HYAM_NAME   : ', TRIM(hyam_name)
    WRITE(*,*) ' HYBM_NAME   : ', TRIM(hybm_name)
    WRITE(*,*) ' Q_NAME      : ', TRIM(q_name)
    WRITE(*,*) ' PV_NAME     : ', TRIM(pv_name)
    WRITE(*,*) ' PT_NAME     : ', TRIM(pt_name)

    ! CLOSE FILE
    CLOSE(iou)

    WRITE(*,*) '==========================================================='

    status = 0

  END SUBROUTINE read_nml
  ! ------------------------------------------------------------------------


end program clustering


! *******************************************************************************
! * Subroutine Section                                                          *
! *******************************************************************************


! --------------------------------------------------------
! Get grid points which are connected to a grid point: BOJAN: this is the 3D-connectedness criterion of the 3D-labelling algorithm
! --------------------------------------------------------

SUBROUTINE connect (outar,label,inx,iny,inz,dirh,dird,diru,nx,ny,nz)

  ! The input array <outar(nx,ny,nz)> is 0/1 labeled. Given a starting point
  ! <inx,iny,inz> where the label is 1, find all points with label 1 which are
  ! connected  and attribute to all these points the label <label>.
  ! The 3D arrays <dir*(nx,ny,nz)> specifie whether an expansion is allowed
  ! dirh=1 where horizontal propagation is allowed,
  ! dird=1 where vertical propagation is allowed downward,
  ! diru=1 where vertical propagation is allowed upward, and
  ! a value of dir*=0 prohibits the respective propagation

  ! Declaration of subroutine parameters
  integer :: nx,ny,nz
  integer :: inx,iny,inz
  integer :: outar(nx,ny,nz)
  integer :: label
  integer :: dirh(nx,ny,nz)
  integer :: diru(nx,ny,nz)
  integer :: dird(nx,ny,nz)

  ! Auxiliary variables
  integer :: il,ir,jb,jf,ku,kd,im,jm,km
  integer :: i,j,k
  integer :: stack
  integer :: indx(nx*ny*nz),indy(nx*ny*nz),indz(nx*ny*nz)

  ! Push the starting elements on the stack
  stack=1
  indx(stack)=inx
  indy(stack)=iny
  indz(stack)=inz
  outar(inx,iny,inz)=label

  ! Define the indices of the neighbouring elements
 100   continue

  il=indx(stack)-1
  if (il.lt.1) il=nx
  ir=indx(stack)+1
  if (ir.gt.nx) ir=1
  jb=indy(stack)-1
  if (jb.lt.1) jb=1
  jf=indy(stack)+1
  if (jf.gt.ny) jf=ny
  ku=indz(stack)+1
  if (ku.gt.nz) ku=nz
  kd=indz(stack)-1
  if (kd.lt.1) kd=1
  im=indx(stack)
  jm=indy(stack)
  km=indz(stack)
  stack=stack-1

  ! Check for index overflow
  if (stack.ge.(nx*ny*nz-nx)) then
     print*,'Stack overflow while clustering...'
     stop
  endif

  ! Mark all connected elements (build up the stack)
  ! Mark level below/above if dird/diru and dirh allowed
  ! BOJAN: U=up, D=down, R=right, L=left, M=mid, F=front, B=back (positions in a 27-element cube of nearest neighbours)

  if ((dirh(im,jm,km).eq.1).and.(dird(im,jm,km).eq.1)) then
     ! below:
     if ( outar(il,jf,kd).eq.1) then
        outar(il,jf,kd)=label
        stack=stack+1
        indx(stack)=il
        indy(stack)=jf
        indz(stack)=kd
     endif
     if ( outar(im,jf,kd).eq.1) then
        outar(im,jf,kd)=label
        stack=stack+1
        indx(stack)=im
        indy(stack)=jf
        indz(stack)=kd
     endif
     if ( outar(ir,jf,kd).eq.1) then
        outar(ir,jf,kd)=label
        stack=stack+1
        indx(stack)=ir
        indy(stack)=jf
        indz(stack)=kd
     endif
     if (outar(il,jm,kd).eq.1) then
        outar(il,jm,kd)=label
        stack=stack+1
        indx(stack)=il
        indy(stack)=jm
        indz(stack)=kd
     endif
     if (outar(ir,jm,kd).eq.1) then
        outar(ir,jm,kd)=label
        stack=stack+1
        indx(stack)=ir
        indy(stack)=jm
        indz(stack)=kd
     endif
     if (outar(il,jb,kd).eq.1) then
        outar(il,jb,kd)=label
        stack=stack+1
        indx(stack)=il
        indy(stack)=jb
        indz(stack)=kd
     endif
     if (outar(im,jb,kd).eq.1) then
        outar(im,jb,kd)=label
        stack=stack+1
        indx(stack)=im
        indy(stack)=jb
        indz(stack)=kd
     endif
     if (outar(ir,jb,kd).eq.1) then
        outar(ir,jb,kd)=label
        stack=stack+1
        indx(stack)=ir
        indy(stack)=jb
        indz(stack)=kd
     endif
  endif
  if ((dirh(im,jm,km).eq.1).and.(diru(im,jm,km).eq.1)) then
     if (outar(il,jf,ku).eq.1) then
        outar(il,jf,ku)=label
        stack=stack+1
        indx(stack)=il
        indy(stack)=jf
        indz(stack)=ku
     endif
     if (outar(im,jf,ku).eq.1) then
        outar(im,jf,ku)=label
        stack=stack+1
        indx(stack)=im
        indy(stack)=jf
        indz(stack)=ku
     endif
     if (outar(ir,jf,ku).eq.1) then
        outar(ir,jf,ku)=label
        stack=stack+1
        indx(stack)=ir
        indy(stack)=jf
        indz(stack)=ku
     endif
     if (outar(il,jm,ku).eq.1) then
        outar(il,jm,ku)=label
        stack=stack+1
        indx(stack)=il
        indy(stack)=jm
        indz(stack)=ku
     endif
     if (outar(ir,jm,ku).eq.1) then
        outar(ir,jm,ku)=label
        stack=stack+1
        indx(stack)=ir
        indy(stack)=jm
        indz(stack)=ku
     endif
     if (outar(il,jb,ku).eq.1) then
        outar(il,jb,ku)=label
        stack=stack+1
        indx(stack)=il
        indy(stack)=jb
        indz(stack)=ku
     endif
     if (outar(im,jb,ku).eq.1) then
        outar(im,jb,ku)=label
        stack=stack+1
        indx(stack)=im
        indy(stack)=jb
        indz(stack)=ku
     endif
     if (outar(ir,jb,ku).eq.1) then
        outar(ir,jb,ku)=label
        stack=stack+1
        indx(stack)=ir
        indy(stack)=jb
        indz(stack)=ku
     endif
  endif
  ! Mark level directly below and above if allowed (dird/u=1)
  if (dird(im,jm,km).eq.1) then
     if (outar(im,jm,kd).eq.1) then
        outar(im,jm,kd)=label
        stack=stack+1
        indx(stack)=im
        indy(stack)=jm
        indz(stack)=kd
     endif
  endif
  if (diru(im,jm,km).eq.1) then
     if (outar(im,jm,ku).eq.1) then
        outar(im,jm,ku)=label
        stack=stack+1
        indx(stack)=im
        indy(stack)=jm
        indz(stack)=ku
     endif
  endif
  ! Mark other points on same level if allowed (dirh=1)
  if (dirh(im,jm,km).eq.1) then
     if (outar(il,jf,km).eq.1) then
        outar(il,jf,km)=label
        stack=stack+1
        indx(stack)=il
        indy(stack)=jf
        indz(stack)=km
     endif
     if (outar(im,jf,km).eq.1) then
        outar(im,jf,km)=label
        stack=stack+1
        indx(stack)=im
        indy(stack)=jf
        indz(stack)=km
     endif
     if (outar(ir,jf,km).eq.1) then
        outar(ir,jf,km)=label
        stack=stack+1
        indx(stack)=ir
        indy(stack)=jf
        indz(stack)=km
     endif
     if (outar(il,jm,km).eq.1) then
        outar(il,jm,km)=label
        stack=stack+1
        indx(stack)=il
        indy(stack)=jm
        indz(stack)=km
     endif
     if (outar(ir,jm,km).eq.1) then
        outar(ir,jm,km)=label
        stack=stack+1
        indx(stack)=ir
        indy(stack)=jm
        indz(stack)=km
     endif
     if (outar(il,jb,km).eq.1) then
        outar(il,jb,km)=label
        stack=stack+1
        indx(stack)=il
        indy(stack)=jb
        indz(stack)=km
     endif
     if (outar(im,jb,km).eq.1) then
        outar(im,jb,km)=label
        stack=stack+1
        indx(stack)=im
        indy(stack)=jb
        indz(stack)=km
     endif
     if (outar(ir,jb,km).eq.1) then
        outar(ir,jb,km)=label
        stack=stack+1
        indx(stack)=ir
        indy(stack)=jb
        indz(stack)=km
     endif
  endif

  if (stack.gt.0) goto 100

  ! Exit point
 700   continue

  return

end subroutine connect

! ---------------------------------------------------------
! Calculate the height of the tropopause: BOJAN: although this subroutine is called tropopause for historical reasons (long story), it actually is the core of the 3D-labelling algorithm :-)
! ---------------------------------------------------------

SUBROUTINE tropopause (f2,f3,p3,pv3,th3,q3,mdv,xlat,xlon,nx,ny,nz,tropo_pv,tropo_th,q_th)    ! BARTUSEK / added xlon

  ! Get the pressure height of the tropopause (f2) and
  ! cluster the atmosphere into labels 1-5 (f3).
  ! Given: 3d field: Pressure <p3(nx,ny,nz)>,
  ! potential vorticity <pv3(nx,ny,nz)>, and potential
  ! temperature <th3(nx,ny,nz)>. The parameters <nx,ny,nz>
  ! characterize the grid. The missing data value is <mdv>.
  ! Bojan July 2013: tropo_pv and tropo_th have to be specified!

  implicit none

  ! Declaration of subroutine parameters
  integer,intent(in) :: nx,ny,nz
  real,intent(out) :: f3(nx,ny,nz)
  real,intent(out) :: f2(nx,ny)
  real,intent(in) :: p3(nx,ny,nz)
  real,intent(in) :: q3(nx,ny,nz)
  real,intent(inout) :: pv3(nx,ny,nz)
  real,intent(in) :: th3(nx,ny,nz)
  real,intent(in) :: xlat(ny)
  real,intent(in) :: xlon(nx)              ! BARTUSEK / added
  integer :: dirsh(nx,ny,nz),dirsd(nx,ny,nz),dirsu(nx,ny,nz)
  integer :: dirth(nx,ny,nz),dirtd(nx,ny,nz),dirtu(nx,ny,nz)
  integer :: dirdh(nx,ny,nz),dirdd(nx,ny,nz),dirdu(nx,ny,nz)
  integer :: dirfh(nx,ny,nz),dirfd(nx,ny,nz),dirfu(nx,ny,nz)
  real :: mdv

  ! Set tropopause thresholds
  real :: tropo_pv,tropo_th
  ! diabatically produced PV (Q-threshold), Feb 2014:
  real :: q_th

  ! Set permission for expansion of stratospheric label. A horizontal
  ! expansion is forbidden if less than <forbid_h> of the grid column
  ! below a grid point is tropospheric
  real :: forbid_h
  parameter (forbid_h=0.5)

  ! Internal variables
  integer :: i,j,k,l
  real :: pv2(nx,ny)
  real :: th2(nx,ny)
  real :: pvn(nx,ny),pvs(nx,ny)
  real :: lat
  integer :: st(nx,ny,nz),tr(nx,ny,nz),de(nx,ny,nz),fi(nx,ny,nz)
  real :: sign
  real :: lev(nx,ny)
  integer :: tschk,kabove,kbelow,vertpvchk
  real :: below, total

  ! ----- STEP 1a ----- BOJAN: these steps should be the same as described in the appendix of my PhD thesis

  ! Mark all points above 50 hPa as stratospheric
  do i=1,nx
     do j=1,ny
        do k=1,nz
           if (xlat(j).ge.0.) then
              sign=1.
           else
              sign=-1.
           endif
           if (p3(i,j,k).lt.50.) then
              pv3(i,j,k)=sign*1.5*tropo_pv
           endif
        enddo
     enddo
  enddo

  ! Mark points with PV > PV threshold (e.g. 2 pvu) OR TH > TH threshold (e.g. 380 K)
  do i=1,nx
     do j=1,ny
        do k=1,nz
           if ((abs(pv3(i,j,k)).ge.tropo_pv).or.(th3(i,j,k).ge.tropo_th)) then
              st(i,j,k)=1 ! 'stratosphere'
              tr(i,j,k)=0 ! 'tropopsphere'
              de(i,j,k)=1 ! additional array for de-treatment
           else
              st(i,j,k)=0
              tr(i,j,k)=1
              de(i,j,k)=0
           endif
        enddo
     enddo
  enddo

  ! ----- STEP 1b -----

  ! Set the expansion permissions for the individual labels (1=allowed, 0=forbidden)
  ! h=horizontal, v=vertical
  do i=1,nx
     do j=1,ny
        do k=1,nz

           ! Set permissions for stratospheric (st) label: always vertical, never horizontal
           dirsh(i,j,k) = 0
           dirsd(i,j,k) = 1
           dirsu(i,j,k) = 1

           ! Set permissions for fill up thing label (fi): always horizontal, only downward, never upward
           dirfh(i,j,k) = 1
           dirfd(i,j,k) = 1
           dirfu(i,j,k) = 0

           ! Set permissions for tropospheric label (can always propagate)
           dirth(i,j,k) = 1
           dirtd(i,j,k) = 1
           dirtu(i,j,k) = 1

           ! Set permissions for deep fold (de) label: always horizontal, vertical only if air is dry
           dirdh(i,j,k) = 1
           ! Allow vertical propagation of de label only for dry air (Q < Q-threshold)
           if ( q3(i,j,k).le.q_th ) then
              dirdd(i,j,k) = 1
              dirdu(i,j,k) = 1
           else
              dirdd(i,j,k) = 0
              dirdu(i,j,k) = 0
           endif
           !endif

        enddo
     enddo
  enddo

  ! ----- STEP 2a -----

  do i=1,nx
     do j=1,ny
       ! Determine stratosphere by connectivity criterion: for all points (x/y) at highest
       ! model level with st=1, find all connected points with st=1 => they get label st=2
       ! this label can only propagate vertically and it is important that this happens before
       ! label st=5 can propagate from the bottom and overwrite st=1 with st=5
        if ((st(i,j,nz).eq.1).and.(p3(i,j,nz).lt.500.)) then
           call connect (st,2,i,j,nz,dirsh,dirsd,dirsu,nx,ny,nz)
        endif
     enddo
  enddo

  do i=1,nx
     do j=1,ny

  ! ----- STEP 2b -----

        ! for de array, do the same thing. The big difference is that the horizontal propagation is
        ! always allowed. In the vertical direction, it is only allowed if RH<50%, ensuring that
        ! very deep tropopause folds, which are not connected by st are reached by de.

        ! This is needed because for some cases, stratospheric air (e.g. at the boundary of a deep-reaching
        ! streamer) is not identified as such (the label st=2 is not allowed to propagate
        ! horizontally). The de=2 label can always propagate horizontally and is in the end
        ! used to distinguish 'real' cutoffs (with st=1 and de=1) from
        ! deep reaching stratospheric air (with st=1 and de=2)
        if ((de(i,j,nz).eq.1).and.(p3(i,j,nz).lt.500.)) then
           call connect (de,2,i,j,nz,dirdh,dirdd,dirdu,nx,ny,nz)
        endif

  ! ----- STEP 2c -----

        ! Determine troposphere analogeously: for all points (x/y) at lowest model
        ! level with tr=1 => all connected points get label tr=2
        if ((tr(i,j,1).eq.1).and.(p3(i,j,1).gt.300.)) then
           call connect (tr,2,i,j,1,dirth,dirtd,dirtu,nx,ny,nz)
        endif

  ! ----- STEP 2d -----

        ! Determine surface-bound PV blobs (|PV|>2): for all points (x/y) at lowest
        ! model level with st=1 => all connected points get label st=5
        ! This label can always propagate! (use dirth,dirtv)
        if ((st(i,j,1).eq.1).and.(p3(i,j,1).gt.300.)) then
           call connect (st,5,i,j,1,dirth,dirtd,dirtu,nx,ny,nz)
        endif

     enddo
  enddo

  ! ----- STEP 2e -----

  ! Remove tropospheric blobs in surface-bound PV blobs (tr=1 surrounded by st=5)
  do i=1,nx
     do j=1,ny
        do k=1,nz

           if (tr(i,j,k).eq.1) then
              tschk=0
              do l=k,nz
              if (tschk.eq.0) then
                 ! Check if air above is stratosphere (st=2) or top of atmosphere (nz)
                 if ( (st(i,j,l).eq.2) .or. (l.eq.nz) ) then
                    tschk=2
                 ! Check if air above is surface-bound PV blob (st=5)
                 elseif (st(i,j,l).eq.5) then
                    tschk=5
                 endif
              endif
              enddo
              ! Mark this blob as tr=5
              if (tschk.eq.5) then
                 call connect (tr,5,i,j,k,dirth,dirtd,dirtu,nx,ny,nz)
              endif
           endif

        enddo
     enddo
  enddo


  ! ----- STEP 3 -----

  ! --------------------------------------------------------------------------------------------------
  ! Set the stratospheric and troposhperic labels (commented out version is 1=strat/0=everything else)
  ! --------------------------------------------------------------------------------------------------
  do i=1,nx
     do j=1,ny
        do k=1,nz
           f3(i,j,k)=0.

           if (th3(i,j,k).gt.tropo_th) then
              f3(i,j,k) = 2.                     ! stratosphere TH (2)
           !   f3(i,j,k) = 1.                     ! stratosphere TH (1)
           else if (st(i,j,k).eq.2) then
              f3(i,j,k) = 2.                     ! stratosphere PV (2)
           !   f3(i,j,k) = 1.                     ! stratosphere PV (1)
           else if (tr(i,j,k).eq.2) then
              f3(i,j,k) = 1.                     ! troposphere (1)
           !   f3(i,j,k) = 0.                     ! troposphere (0)
           else if ((st(i,j,k).eq.1).and.(de(i,j,k).eq.1)) then
              f3(i,j,k) = 3.                     ! strat cutoff (3)
           !   f3(i,j,k) = 0.                     ! strat cutoff (0)
           else if ((st(i,j,k).eq.1).and.(de(i,j,k).eq.2)) then
              f3(i,j,k) = 2.                     ! deep reaching stratospheric air (2)
           !   f3(i,j,k) = 1.                     ! deep reaching stratospheric air (1)
           else if (tr(i,j,k).eq.1) then
              f3(i,j,k) = 4.                     ! trop cutoff (4)
           !   f3(i,j,k) = 0.                     ! trop cutoff (0)
           else if (st(i,j,k).eq.5) then
              f3(i,j,k) = 5.                     ! surface-bound PV blob (5)
           !   f3(i,j,k) = 0.                     ! surface-bound PV blob (0)
           else if (tr(i,j,k).eq.5) then
              f3(i,j,k) = 1. ! tropospheric air within surface-bound PV blob (1)
           !   f3(i,j,k) = 0. ! tropospheric air within surface-bound PV blob (0)
           endif

        enddo
     enddo
  enddo

  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------
  ! ----- STEP 6 -----
  do i=1,nx
     do j=1,ny
        do k=1,nz
           fi(i,j,k)=0
           if (( f3(i,j,k).eq.2 ) .or. ( f3(i,j,k).eq.3 ).or.( f3(i,j,k).eq.5 )) then
              fi(i,j,k)=1 ! mark point as possible horizontal fill-up point
           endif
         end do
      end do
   end do

   ! BARTUSEK / SKIP FOLLOWING "AESTHETIC" PART; PRODUCES SOME MIS-LABELLING
   goto 777    ! BARTUSEK / added
   ! Horizontally fill up 2/3 and 2/5 patches (only aestetically important ;-))
   do i=1,nx
      do j=1,ny
         do k=1,nz
            if ( f3(i,j,k).eq.3 ) then ! label 3
               call connect (fi,3,i,j,k,dirfh,dirfd,dirfu,nx,ny,nz)
            endif

            if ( f3(i,j,k).eq.5 ) then ! label 3
               call connect (fi,5,i,j,k,dirfh,dirfd,dirfu,nx,ny,nz)
            endif
         enddo
      enddo
   enddo
   do i=1,nx
      do j=1,ny
         do k=1,nz
            if ( fi(i,j,k).eq.3  ) then ! label 3 that was horizontally filled up
               f3(i,j,k)=3.
               do l=1,k-1 ! vertically fill up points below if label is 2
                  if ( f3(i,j,l).eq.2 ) then
                     f3(i,j,l)=3.
                  end if
               end do
            endif
            if ( fi(i,j,k).eq.5  ) then ! label 5 that was horizontally filled up
               f3(i,j,k)=5.
            endif
         enddo
      enddo
   enddo
   777 continue    ! BARTUSEK / added

  ! CHECK that all points obtained a label!
  do i=1,nx
     do j=1,ny
        do k=1,nz
           if ( f3(i,j,k).eq.0 ) then           ! no label assigned
              print*,'Error: no label assigned'
              print*,'i,j,k = ',i,j,k
              print*,'p3(i,j,k) = ',p3(i,j,k)
              print*,'st(i,j,k) = ',st(i,j,k)
              print*,'tr(i,j,k) = ',tr(i,j,k)
              print*,'de(i,j,k) = ',de(i,j,k)
              stop
           endif
        enddo
     enddo
  enddo

  ! --------------------------------------------------------------------------------------------------
  ! DONE with assigning labels, now interpolate field to tropopause pressure
  ! --------------------------------------------------------------------------------------------------
   ! BOJAN: here, we calculate the height of the tropopause (folds are ignored, see PhD thesis)
  ! Remove stratospheric blobs in troposphere and tropospheric blobs in stratosphere
  ! Brute force: Change PV values above and below tropo threshold
  do i=1,nx
     do j=1,ny
        do k=1,nz
           if (xlat(j).ge.0.) then
              sign=1.
           else
              sign=-1.
           endif

           ! Set PV of tropo. cutoff in stratosphere
           if ( (st(i,j,k).eq.0).and.(tr(i,j,k).eq.1) ) then
              pv3(i,j,k)=sign*1.5*tropo_pv
           ! Set PV of exotic cases where PV < -2 pvu (PV > +2 pvu) in NH (SH) => would otherwise be detected as TP
           else if ( ((sign .gt. 0.).and.(pv3(i,j,k).lt.-tropo_pv).or.(sign .lt. 0.).and.(pv3(i,j,k).gt.tropo_pv)).and.(th3(i,j,k).gt.tropo_pv) ) then
              pv3(i,j,k)=sign*1.5*tropo_pv
           ! Set PV of strat. cutoff in trop. or surface bound PV (only below 380K)
           !else if ( (((st(i,j,k).eq.1).and.(de(i,j,k).eq.1)).or.(st(i,j,k).eq.5)).and.(th3(i,j,k).le.tropo_th) ) then
!           else if ( (f3(i,j,k).eq.3).or.(f3(i,j,k).eq.5) ) then        						! BARTUSEK / commented out -> causes folds affected by Sprenger/Bartusek edit below to fail to be identified (re-implemented in edit featuring "0.5*tropo_pv" below)
!              pv3(i,j,k)=sign*0.5*tropo_pv														! BARTUSEK / commented out
!           endif																				! BARTUSEK / commented out
!           else if ( (f3(i,j,k).eq.3).or.(f3(i,j,k).eq.5).and.( q3(i,j,k).ge.q_th ) ) then  	! BARTUSEK / commented out     
!              pv3(i,j,k)=sign*0.5*tropo_pv														! BARTUSEK / commented out
           endif
        enddo
     enddo
  enddo

  ! Calculate the height of the dynamical tropopause (NH:2pvu, SH:-2pvu)
  do i=1,nx
     do j=1,ny
        lev(i,j)=tropo_pv
    enddo
  enddo
  call thipo_lev(p3,pv3,lev,pvn,nx,ny,nz,mdv,0) ! for NH
  do i=1,nx
     do j=1,ny
        do k=1,nz
           pv3(i,j,k)=-pv3(i,j,k) ! change sign of PV
        enddo
     enddo
  enddo
  call thipo_lev(p3,pv3,lev,pvs,nx,ny,nz,mdv,0) ! for SH
  do i=1,nx
     do j=1,ny
        do k=1,nz
           pv3(i,j,k)=-pv3(i,j,k) ! change sign of PV back
        enddo
     enddo
  enddo
  do j=1,ny
     do i=1,nx
        if (xlat(j).ge.0.) then
           pv2(i,j)=pvn(i,j) ! NH
        else
           pv2(i,j)=pvs(i,j) ! SH
        endif
     enddo
  enddo

  ! Special treatment for PV towers (in whole column, |PV|>2)
  do i=1,nx
     do j=1,ny
        if (pv2(i,j).eq.mdv) then
           vertpvchk=0
           do k=1,nz
              if (abs(pv3(i,j,k)).ge.tropo_pv) vertpvchk=vertpvchk+1
           enddo
           if(vertpvchk.eq.nz) then
              !print*,'PV tower! Tropopause reaches ground: ',p3(i,j,1)
              pv2(i,j)=p3(i,j,1)
           endif
        endif
     enddo
  enddo

  ! Calculate the height of the 380 K surface
  do i=1,nx
     do j=1,ny
        lev(i,j)=tropo_th
    enddo
  enddo
  call thipo_lev(p3,th3,lev,th2,nx,ny,nz,mdv,0)

  ! Set unreasonable values to mdv
  do i=1,nx
     do j=1,ny
        if ((th2(i,j).lt.40.).or.(th2(i,j).gt.800.)) then ! Out of 'regular' limits => check label
           if (th2(i,j).ne.mdv) then
              !print*,'Warning: p(TP) (TH) out of limits:',th2(i,j)
              !print*,'i,j = ',i,j
              !print*,'LABEL profile:',f3(i,j,:)
              do k=1,nz
                 if (p3(i,j,k).ge.th2(i,j)) then
                    kbelow=k
                 endif
                 if (p3(i,j,nz+1-k).le.th2(i,j)) then
                    kabove=nz+1-k
                 endif
              enddo
              if ((f3(i,j,kbelow).eq.1).and.(f3(i,j,kabove).eq.2)) then
                 !print*,'Point is okay!'
              else
                 !print*,'Point is not okay => set to mdv'
                 th2(i,j)=mdv
              endif
           endif
        endif
        if ((pv2(i,j).lt.40.).or.(pv2(i,j).gt.800.)) then ! Out of 'regular' limits => check label
           if (pv2(i,j).ne.mdv) then
              !print*,'Warning: p(TP) (PV) out of limits:',pv2(i,j)
              !print*,'i,j = ',i,j
              !print*,'LABEL profile:',f3(i,j,:)
              do k=1,nz
                 if (p3(i,j,k).ge.pv2(i,j)) then
                    kbelow=k
                 endif
                 if (p3(i,j,nz+1-k).le.pv2(i,j)) then
                    kabove=nz+1-k
                 endif
              enddo
              if ((f3(i,j,kbelow).eq.1).and.(f3(i,j,kabove).eq.2)) then
                 !print*,'Point is okay!'
              elseif ((kabove.eq.1).and.(kbelow.eq.1)) then
                 !print*,'PV tower! Tropopause reaches ground: ',pv2(i,j)
              else
                 !print*,'Point is not okay => set to mdv'
                 pv2(i,j)=mdv
              endif
           endif

        endif
     enddo
  enddo

  ! --------------------------------------------------------------------------------------------------
  ! Set the resulting tropopause height (maximum pressure: 380 K, 2 PVU)
  ! --------------------------------------------------------------------------------------------------
  do i=1,nx
     do j=1,ny
        if ((th2(i,j).ne.mdv).and.(pv2(i,j).eq.mdv)) then
           f2(i,j)=th2(i,j)
        else if ((th2(i,j).eq.mdv).and.(pv2(i,j).ne.mdv)) then
           f2(i,j)=pv2(i,j)
        else if (th2(i,j).gt.pv2(i,j)) then
           f2(i,j)=th2(i,j)
        else
           f2(i,j)=pv2(i,j)
        endif
     enddo
  enddo
  ! --------------------------------------------------------------------------------------------------
  ! DONE
  ! --------------------------------------------------------------------------------------------------

  !   ! -------------------------------------------------------------------------------------------------
  !   ! Handle label 5 that propagates up to the stratosphere                  SPRENGER ADDITION 12/20/20    ! BARTUSEK / commented out -> implemented within BARTUSEK addition just below
  !   ! Compare distance to tropopause to distance from surface: the smaller wins
  !   ! -------------------------------------------------------------------------------------------------
  !   goto 2008   !!!! Skipping to see if tibetanedit is good enough to replace
  !   do i=1,nx
  !    do j=1,ny
  !      do k=1,nz
  !        if ( f3(i,j,k).eq.5 ) then
  !           if ( (p3(i,j,k)-f2(i,j)).lt.(p3(i,j,1)-p3(i,j,k)) ) then
  !               f3(i,j,k) = 2.
  !           endif
  !        endif
  !      enddo
  !      enddo
  !   enddo
  !   2008 continue
  !   ! --------------------------------------------------------------------------------------------------
  !   ! DONE
  !   ! --------------------------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------------------------
  ! BARTUSEK / ADDITION; Handle label 5 that propagates up too high, and disallow 2 within q threshold      
  ! -------------------------------------------------------------------------------------------------
  do i=1,nx
   do j=1,ny
     do k=1,nz
       !!!!! First change 2s that are outside of q_th -->5
       if ( ( f3(i,j,k).eq.2 ).and.( q3(i,j,k).ge.q_th ) ) then
          f3(i,j,k) = 5.
       endif
       !!!!! Change 5s/3s that are (a) within q_th, and (b) are connected to strat, and (c) are closer (in pressure) to tropopause than ground (except over Tibetan Plateau; see lat/lon ranges) -->2
       if ( (( f3(i,j,k).eq.5 ).or.( f3(i,j,k).eq.3 )).and.( de(i,j,k).eq.2 ).and.( q3(i,j,k).le.q_th ) ) then
           if ( ( xlat(j).ge.25 ).and.( xlat(j).le.40 ).and.( xlon(i).ge.75 ).and.( xlon(i).le.110 ) ) then
              f3(i,j,k) = 2.
           else
              if ( ( (p3(i,j,k)-f2(i,j)).lt.(p3(i,j,1)-p3(i,j,k)) ).and.( (p3(i,j,1)-f2(i,j)).ge.200 ) ) then
                 f3(i,j,k) = 2.
              endif
           endif
       endif
     enddo
   enddo
  enddo
  ! --------------------------------------------------------------------------------------------------
  ! BARTUSEK / END
  ! --------------------------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------------------------
  ! BARTUSEK / ADDITION; Reset PV in label 3/5 to 0.5*tropo_pv (as previously done above but commented out)
  ! -------------------------------------------------------------------------------------------------
  do i=1,nx
     do j=1,ny
        do k=1,nz

           if ( (f3(i,j,k).eq.3).or.(f3(i,j,k).eq.5) ) then
              pv3(i,j,k)=sign*0.5*tropo_pv
           endif

        enddo
     enddo
  enddo
  ! --------------------------------------------------------------------------------------------------
  ! BARTUSEK / END
  ! --------------------------------------------------------------------------------------------------

end subroutine tropopause

SUBROUTINE tropofold (label3,p3,pv3,th3,fold,dp, pmin, pmax, sfold,mfold,dfold,mdv,nx,ny,nz,tropo_pv,tropo_th)

  ! -------------------------------------------------------------------------------
  ! Part 2) Identify tropopause folds (multiple crossings of the tropopause in one column)
  ! -------------------------------------------------------------------------------

  ! Given the label field, determine whether a tropopause fold is present.
  ! The parameters <nx,ny,nx> characterize the grid. The missing data value is <mdv>.

  implicit none

  ! Declaration of subroutine parameters
  integer :: nx,ny,nz
  real :: mdv
  real :: label3(nx,ny,nz)
  real :: p3(nx,ny,nz)
  real :: pv3(nx,ny,nz)
  real :: th3(nx,ny,nz)
  real :: pre0,pre1,pre2,ok1,ok2
  real :: tropo_pv,tropo_th
  real :: fold(nx,ny)
  real :: dp(nx,ny)
  real :: pmin(nx,ny)
  real :: pmax(nx,ny)
  real :: sfold(nx,ny)
  real :: mfold(nx,ny)
  real :: dfold(nx,ny)

  ! Aux integers
  integer :: i,j,k,l,k0,k1,k2,ncross

  ! Adjust fields (integer labels, absolute value of PV)
  do i=1,nx
     do j=1,ny
        do k=1,nz
           pv3(i,j,k) = abs(pv3(i,j,k))
        enddo
     enddo
  enddo

  ! Definition of labels
  !     1 : troposphere
  !     2 : stratosphere
  !     3 : stratosperic cutoff
  !     4 : tropospheric cutoff
  !     5 : surface-bound PV blobs

  do i=1,nx
  do j=1,ny
     k0 = 0
     k1 = 0
     k2 = 0
     ncross = 0

     do k=nz-1,1,-1

        ! strat-trop transition
        if ( (label3(i,j,k+1).eq.2).and.(label3(i,j,k).eq.1) ) then

           ncross = ncross + 1
           if ( k2.eq.0 ) then
              k2 = k
           elseif ( k1.ne.0 ) then
              k0 = k
           endif
        endif
        ! special case for folds that are affected by the q-criterion (q<0.0001 only)				! BARTUSEK / commented out -> implemented within edits below
        ! example (top to bottom): 222211122231111 need to recognize the 2-3-1 transition			! BARTUSEK / commented out 
        ! transitions like 222222223111 should not be recognized, thus only if k2.ne.0				! BARTUSEK / commented out 
!        if ( (label3(i,j,k+1).eq.2).and.(label3(i,j,k).eq.3) ) then								! BARTUSEK / commented out 
!																									! BARTUSEK / commented out 
!           ncross = ncross + 1																		! BARTUSEK / commented out 
!           if (( k2.ne.0 ) .and. ( k1.ne.0 )) then													! BARTUSEK / commented out 
!              k0 = k																				! BARTUSEK / commented out 
!           endif																					! BARTUSEK / commented out 
!        endif																						! BARTUSEK / commented out 

        ! --------------
        ! BARTUSEK / ADDITION; RECOGNIZE CASES EXPANDING ON "SPECIAL CASE" ABOVE, BUT EITHER (a) WITH 5/3 BETWEEN 2 AND 1 OR (b) WITH 5/3 BETWEEN 1 AND 2
        ! --------------
        ! Example (top to bottom): 22225111122251111 need to recognize the strat-trop 2-5-1 transition(s)    ! Note, this won't create additional spurious crossings; 2-to-5 is the strat-to-trop, and the 5-to-1 won't be caught by anything here
        if ( (label3(i,j,k+1).eq.2).and.(label3(i,j,k).eq.5) ) then
           ncross = ncross + 1
           if ( k2.eq.0 ) then   !!! Set upper strat-trop transition
              k2 = k
           elseif (( k2.ne.0 ) .and. ( k1.ne.0 )) then
              k0 = k    !!! Set lower strat-trop transition
           endif
        endif

        ! Example (top to bottom): 22223111122231111 need to recognize the strat-trop 2-3-1 transition(s)    ! Note, this won't create additional spurious crossings; 2-to-3 is the strat-to-trop, and the 3-to-1 won't be caught by anything here
        if ( (label3(i,j,k+1).eq.2).and.(label3(i,j,k).eq.3) ) then
           ncross = ncross + 1
           if ( k2.eq.0 ) then   !!! Set upper strat-trop transition
              k2 = k
           elseif (( k2.ne.0 ) .and. ( k1.ne.0 )) then
              k0 = k    !!! Set lower strat-trop transition
           endif
        endif

        ! Recognize trop-strat (3-2) transition for (top to bottom) …221323211…, so need k2.ne.0
        if ( (label3(i,j,k+1).eq.3).and.(label3(i,j,k).eq.2) ) then
           ncross = ncross + 1
           if ( ( k2.ne.0 ).and.( k1.eq.0 ) ) then
              k1 = k
           endif
        endif

        ! Recognize trop-strat (5-2) transition for (top to bottom) …221525211…, so need k2.ne.0
        if ( (label3(i,j,k+1).eq.5).and.(label3(i,j,k).eq.2) ) then
           ncross = ncross + 1
           if ( ( k2.ne.0 ).and.( k1.eq.0 ) ) then
              k1 = k
           endif
        endif
        ! --------------
        ! BARTUSEK / END
        ! --------------

        ! trop-strat transition
        if ( (label3(i,j,k+1).eq.1).and.(label3(i,j,k).eq.2) ) then

           ncross = ncross + 1
           if ( ( k2.ne.0 ).and.( k1.eq.0 ) ) then
              k1 = k
           endif
        endif
     enddo

     ! Skip further steps if we have a single tropopause
     if ( ncross.le.2 ) goto 100

     ! --------------
     ! BARTUSEK / ADDITION; Troubleshooting fix
     ! --------------     
     ! If mismatch between ncross and # of non-zero k's (can extremely rarely happen with non-orthodox labelling), no fold identified
     if ( k0.eq.0 ) goto 100
     if ( k1.eq.0 ) goto 100
     if ( k2.eq.0 ) goto 100
     ! --------------
     ! BARTUSEK / END
     ! --------------

     ! Get the exact pressures at the crossings
     pre0 = 0.
     pre1 = 0.
     pre2 = 0.
     ! lowest point (0) => pre0
     ok1 = ( pv3(i,j,k0+1)-tropo_pv ) * ( pv3(i,j,k0)-tropo_pv )
     ok2 = ( th3(i,j,k0+1)-tropo_th ) * ( th3(i,j,k0)-tropo_th )

     if ( ok1.le.0. ) then
        pre0 = p3(i,j,k0) + ( p3(i,j,k0+1) - p3(i,j,k0) ) * &
              	      ( tropo_pv - pv3(i,j,k0) )/( pv3(i,j,k0+1) - pv3(i,j,k0) )
     elseif ( ok2.le.0. ) then
        pre0 = p3(i,j,k0) + ( p3(i,j,k0+1) - p3(i,j,k0) ) * &
            	      ( tropo_th - th3(i,j,k0) )/( th3(i,j,k0+1) - th3(i,j,k0) )
     endif

     ! middle point (1) => pre1
     if ( k1.eq.0 ) then
        print *, label3(i,j,:)
        print "(i6.3, i6.3, i6.3)", k0, k1, k2
        !print "(f15.7, f15.7)", xlat(j), xlon(i)
     endif
     ok1 = ( pv3(i,j,k1+1)-tropo_pv ) * ( pv3(i,j,k1)-tropo_pv )
     ok2 = ( th3(i,j,k1+1)-tropo_th ) * ( th3(i,j,k1)-tropo_th )
     if ( ok1.le.0. ) then
        pre1 = p3(i,j,k1) + ( p3(i,j,k1+1) - p3(i,j,k1) ) * &
              	      ( tropo_pv - pv3(i,j,k1) )/( pv3(i,j,k1+1) - pv3(i,j,k1) )
     elseif ( ok2.le.0. ) then
        pre1 = p3(i,j,k1) + ( p3(i,j,k1+1) - p3(i,j,k1) ) * &
            	      ( tropo_th - th3(i,j,k1) )/( th3(i,j,k1+1) - th3(i,j,k1) )
     endif

     ! upper point (2) => pre2
     ok1 = ( pv3(i,j,k2+1)-tropo_pv ) * ( pv3(i,j,k2)-tropo_pv )
     ok2 = ( th3(i,j,k2+1)-tropo_th ) * ( th3(i,j,k2)-tropo_th )
     if ( ok1.le.0. ) then
        pre2 = p3(i,j,k2) + ( p3(i,j,k2+1) - p3(i,j,k2) ) * &
              	      ( tropo_pv - pv3(i,j,k2) )/( pv3(i,j,k2+1) - pv3(i,j,k2) )
     elseif ( ok2.le.0. ) then
        pre2 = p3(i,j,k2) + ( p3(i,j,k2+1) - p3(i,j,k2) ) * &
            	      ( tropo_th - th3(i,j,k2) )/( th3(i,j,k2+1) - th3(i,j,k2) )
     endif
     
     ! Decide whether all pressure values are ok
     if ( ( pre0.lt.p3(i,j,k0+1) ).or. &
          ( pre0.gt.p3(i,j,k0  ) ).or. &
          ( pre1.lt.p3(i,j,k1+1) ).or. &
          ( pre1.gt.p3(i,j,k1  ) ).or. &
          ( pre2.lt.p3(i,j,k2+1) ).or. &
          ( pre2.gt.p3(i,j,k2  ) ) ) then
          print "(f6.3,f6.3,f6.3)",pre0,pre1,pre2
          goto 100
     endif

     ! Everything is fine: remember the fold
     fold(i,j) = 1.
     dp(i,j) = pre1 - pre2
     pmin(i,j) = pre2
     pmax(i,j) = pre0

     ! Exit point for loop
100  continue

  enddo ! y
  enddo ! x

  where((fold .gt. 0.0) .and. (dp .ge. 50.0) .and. (dp .lt. 200.0))
     sfold=1.0
  elsewhere((fold .gt. 0.0) .and. (dp .ge. 200.0) .and. (dp .lt. 350.0))
     mfold=1.0
  elsewhere((fold .gt. 0.0) .and. (dp .ge. 350.0))
     dfold=1.0
  end where

end subroutine tropofold

! -----------------------------------------------------------
! Interpolation onto theta surface (top -> down): BOJAN: not sure you have all the interpolation subroutines used here, contact us or feel free to use any other interpolation routine that does the job :-)
! -----------------------------------------------------------
subroutine thipo_lev(var3d,th3d,lev,var,nx,ny,nz,mdv,mode)

  ! Interpolates the 3d variable var3d on the isentropic surface
  ! defined by lev. The interpolated field is returned as var.
  ! th3d denotes the 3d theta array.
  ! mode determines the way of vertical interpolation:
  ! mode=0 is for linear interpolation
  ! mode=1 is for logarithmic interpolation

  integer :: nx,ny,nz,mode
  real :: lev(nx,ny),mdv
  real :: var3d(nx,ny,nz),th3d(nx,ny,nz),var(nx,ny)

  integer :: i,j,k
  real :: kind
  real :: int3dm

  do i=1,nx
     do j=1,ny

        if (lev(i,j).ne.mdv) then

           kind=0.
           do k=nz-1,1,-1
              if ((th3d(i,j,k).le.lev(i,j)).and.(th3d(i,j,k+1).ge.lev(i,j))) then
                 kind=float(k)+(th3d(i,j,k)-lev(i,j))/(th3d(i,j,k)-th3d(i,j,k+1))
                 goto 100
              endif
           enddo
 100       continue

           if (kind.eq.0) then
              var(i,j)=mdv
           else
              var(i,j)=int3dm(var3d,nx,ny,nz,float(i),float(j),kind,mdv)
           endif

        else

           var(i,j)=mdv

        endif

     enddo
  enddo

  return
end subroutine thipo_lev

! -----------------------------------------------------------
! Interpolation onto theta surface (bottom -> up)
! -----------------------------------------------------------

subroutine pipo_lev(var3d,p3d,lev,var,nx,ny,nz,mdv,mode)

  ! Interpolates the 3d variable var3d on the pressure surface
  ! defined by lev. The interpolated field is returned as var.
  ! p3d denotes the 3d pressure array.
  ! mode determines the way of vertical interpolation:
  ! mode=0 is for linear interpolation
  ! mode=1 is for logarithmic interpolation

  integer :: nx,ny,nz,mode
  real :: lev(nx,ny),mdv
  real :: var3d(nx,ny,nz),p3d(nx,ny,nz),var(nx,ny)

  integer :: i,j,k
  real :: kind
  real :: int3dm

  do i=1,nx
     do j=1,ny

        if (lev(i,j).ne.mdv) then

           kind=0.
           do k=1,nz-1
              if ((p3d(i,j,k).ge.lev(i,j)).and.(p3d(i,j,k+1).le.lev(i,j))) then
                 kind=float(k)+(p3d(i,j,k)-lev(i,j))/(p3d(i,j,k)-p3d(i,j,k+1))
                 goto 100
              endif
           enddo
 100       continue

           if (kind.eq.0.) then
              var(i,j)=mdv
           else
              var(i,j)=int3dm(var3d,nx,ny,nz,float(i),float(j),kind,mdv)
           endif

        else

           var(i,j)=mdv

        endif

     enddo
  enddo

  return
end subroutine pipo_lev


! ------------------------------------------------------------------
! Auxiliary routines
! ------------------------------------------------------------------
