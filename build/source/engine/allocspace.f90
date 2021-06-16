! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2020 NCAR/RAL; University of Saskatchewan; University of Washington
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module allocspace_module

! data types
USE nrtype

! provide access to the derived types to define the data structures
USE data_types,only:&
                    ! final data vectors
                    dlength,             & ! var%dat
                    ilength,             & ! var%dat
                    ! no spatial dimension
                    var_i,               & ! x%var(:)            (i4b)
                    var_i8,              & ! x%var(:)            integer(8)
                    var_d,               & ! x%var(:)            (dp)
                    var_flagVec,         & ! x%var(:)%dat        (logical)
                    var_ilength,         & ! x%var(:)%dat        (i4b)
                    var_dlength,         & ! x%var(:)%dat        (dp)
                    ! gru dimension
                    gru_int,             & ! x%gru(:)%var(:)     (i4b)
                    gru_int8,            & ! x%gru(:)%var(:)     integer(8)
                    gru_double,          & ! x%gru(:)%var(:)     (dp)
                    gru_intVec,          & ! x%gru(:)%var(:)%dat (i4b)
                    gru_doubleVec,       & ! x%gru(:)%var(:)%dat (dp)
                    ! gru+hru dimension
                    gru_hru_int,         & ! x%gru(:)%hru(:)%var(:)     (i4b)
                    gru_hru_int8,        & ! x%gru(:)%hru(:)%var(:)     integer(8)
                    gru_hru_double,      & ! x%gru(:)%hru(:)%var(:)     (dp)
                    gru_hru_intVec,      & ! x%gru(:)%hru(:)%var(:)%dat (i4b)
                    gru_hru_doubleVec,   & ! x%gru(:)%hru(:)%var(:)%dat (dp)
                    ! gru+hru+z dimension
                    gru_hru_z_vLookup      ! x%gru(:)%hru(:)%z(:)%var(:)%lookup (dp)

! metadata structure
USE data_types,only:var_info               ! data type for metadata

! access missing values
USE globalData,only:integerMissing         ! missing integer
USE globalData,only:realMissing            ! missing double precision number

USE globalData,only: nTimeDelay            ! number of timesteps in the time delay histogram
USE globalData,only: nBand                 ! number of spectral bands

! access variable types
USE var_lookup,only:iLookVarType           ! look up structure for variable typed
USE var_lookup,only:maxvarFreq             ! allocation dimension (output frequency)

! privacy
implicit none
private
public::allocGlobal
public::allocLocal
public::resizeData

! -----------------------------------------------------------------------------------------------------------------------------------
contains

 ! ************************************************************************************************
 ! public subroutine allocGlobal: allocate space for global data structures
 ! ************************************************************************************************
 subroutine allocGlobal(metaStruct,dataStruct,err,message)
 ! NOTE: safety -- ensure only used in allocGlobal
 USE globalData,only: gru_struc     ! gru-hru mapping structures
 implicit none
 ! input
 type(var_info),intent(in)       :: metaStruct(:)  ! metadata structure
 ! output
 class(*),intent(out)            :: dataStruct     ! data structure
 integer(i4b),intent(out)        :: err            ! error code
 character(*),intent(out)        :: message        ! error message
 ! local variables
 logical(lgt)                    :: check          ! .true. if structure is already allocated
 integer(i4b)                    :: iHRU           ! loop index through HRUs
 integer(i4b)                    :: iGRU           ! loop index through GRUs
 integer(i4b)                    :: nGRU           ! number of GRUs
 logical(lgt)                    :: spatial        ! spatial flag
 character(len=256)              :: cmessage       ! error message of the downwind routine
 ! initialize error control
 err=0; message='allocGlobal/'

 ! initialize allocation check
 check=.false.

 ! get the number of GRUs
 nGRU = size(gru_struc)

 ! * allocate GRU dimension
 select type(dataStruct)
  ! gru dimension only
  class is (gru_int);           if(allocated(dataStruct%gru))then; check=.true.; else; allocate(dataStruct%gru(nGRU),stat=err); end if
  class is (gru_int8);          if(allocated(dataStruct%gru))then; check=.true.; else; allocate(dataStruct%gru(nGRU),stat=err); end if
  class is (gru_intVec);        if(allocated(dataStruct%gru))then; check=.true.; else; allocate(dataStruct%gru(nGRU),stat=err); end if
  class is (gru_double);        if(allocated(dataStruct%gru))then; check=.true.; else; allocate(dataStruct%gru(nGRU),stat=err); end if
  class is (gru_doubleVec);     if(allocated(dataStruct%gru))then; check=.true.; else; allocate(dataStruct%gru(nGRU),stat=err); end if
  ! gru+hru dimensions
  class is (gru_hru_int);       if(allocated(dataStruct%gru))then; check=.true.; else; allocate(dataStruct%gru(nGRU),stat=err); end if
  class is (gru_hru_int8);      if(allocated(dataStruct%gru))then; check=.true.; else; allocate(dataStruct%gru(nGRU),stat=err); end if
  class is (gru_hru_intVec);    if(allocated(dataStruct%gru))then; check=.true.; else; allocate(dataStruct%gru(nGRU),stat=err); end if
  class is (gru_hru_double);    if(allocated(dataStruct%gru))then; check=.true.; else; allocate(dataStruct%gru(nGRU),stat=err); end if
  class is (gru_hru_doubleVec); if(allocated(dataStruct%gru))then; check=.true.; else; allocate(dataStruct%gru(nGRU),stat=err); end if
 ! gru+hru+z dimensions
  class is (gru_hru_z_vLookup); if(allocated(dataStruct%gru))then; check=.true.; else; allocate(dataStruct%gru(nGRU),stat=err); end if
 end select

 ! check errors
 if(check) then; err=20; message=trim(message)//'GRU structure was unexpectedly allocated already'; return; end if
 if(err/=0)then; err=20; message=trim(message)//'problem allocating GRU dimension'; return; end if

 ! * allocate HRU dimension
 do iGRU=1,nGRU
  ! allocate the HRU dimension
  select type(dataStruct)
   class is (gru_hru_int);       if(allocated(dataStruct%gru(iGRU)%hru))then; check=.true.; else; allocate(dataStruct%gru(iGRU)%hru(gru_struc(iGRU)%hruCount),stat=err); end if
   class is (gru_hru_int8);      if(allocated(dataStruct%gru(iGRU)%hru))then; check=.true.; else; allocate(dataStruct%gru(iGRU)%hru(gru_struc(iGRU)%hruCount),stat=err); end if
   class is (gru_hru_intVec);    if(allocated(dataStruct%gru(iGRU)%hru))then; check=.true.; else; allocate(dataStruct%gru(iGRU)%hru(gru_struc(iGRU)%hruCount),stat=err); end if
   class is (gru_hru_double);    if(allocated(dataStruct%gru(iGRU)%hru))then; check=.true.; else; allocate(dataStruct%gru(iGRU)%hru(gru_struc(iGRU)%hruCount),stat=err); end if
   class is (gru_hru_doubleVec); if(allocated(dataStruct%gru(iGRU)%hru))then; check=.true.; else; allocate(dataStruct%gru(iGRU)%hru(gru_struc(iGRU)%hruCount),stat=err); end if
 class is (gru_hru_z_vLookup); if(allocated(dataStruct%gru(iGRU)%hru))then; check=.true.; else; allocate(dataStruct%gru(iGRU)%hru(gru_struc(iGRU)%hruCount),stat=err); end if
   class default  ! do nothing: It is acceptable to not be any of these specified cases
  end select
  ! check errors
  if(check) then; err=20; message=trim(message)//'HRU structure was unexpectedly allocated already'; return; end if
  if(err/=0)then; err=20; message=trim(message)//'problem allocating HRU dimension'; return; end if
 end do

 ! * allocate local data structures where there is a spatial dimension
 gruLoop: do iGRU=1,nGRU

  ! initialize the spatial flag
  spatial=.false.

  ! loop through HRUs
  hruLoop: do iHRU=1,gru_struc(iGRU)%hruCount

   ! get the number of snow and soil layers
   associate(&
   nSnow => gru_struc(iGRU)%hruInfo(iHRU)%nSnow, & ! number of snow layers for each HRU
   nSoil => gru_struc(iGRU)%hruInfo(iHRU)%nSoil  ) ! number of soil layers for each HRU

   ! allocate space for structures WITH an HRU dimension
   select type(dataStruct)
    class is (gru_hru_int);       call allocLocal(metaStruct,dataStruct%gru(iGRU)%hru(iHRU),nSnow,nSoil,err,cmessage); spatial=.true.
    class is (gru_hru_int8);      call allocLocal(metaStruct,dataStruct%gru(iGRU)%hru(iHRU),nSnow,nSoil,err,cmessage); spatial=.true.
    class is (gru_hru_intVec);    call allocLocal(metaStruct,dataStruct%gru(iGRU)%hru(iHRU),nSnow,nSoil,err,cmessage); spatial=.true.
    class is (gru_hru_double);    call allocLocal(metaStruct,dataStruct%gru(iGRU)%hru(iHRU),nSnow,nSoil,err,cmessage); spatial=.true.
    class is (gru_hru_doubleVec); call allocLocal(metaStruct,dataStruct%gru(iGRU)%hru(iHRU),nSnow,nSoil,err,cmessage); spatial=.true.
    class is (gru_hru_z_vLookup); spatial=.true. ! (special case, allocate space separately later)
    class default; exit hruLoop
   end select

   ! error check
   if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

   ! end association to info in data structures
   end associate

  end do hruLoop ! loop through HRUs

  ! allocate space for structures *WITHOUT* an HRU dimension
  select type(dataStruct)
   class is (gru_double);    call allocLocal(metaStruct,dataStruct%gru(iGRU),nSnow=0,nSoil=0,err=err,message=cmessage); spatial=.true.
   class is (gru_doubleVec); call allocLocal(metaStruct,dataStruct%gru(iGRU),nSnow=0,nSoil=0,err=err,message=cmessage); spatial=.true.
   class default
    if(.not.spatial) exit gruLoop  ! no need to allocate spatial dimensions if none exist for a given variable
    cycle gruLoop  ! can have an HRU dimension if we get to here
  end select

  ! error check
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

 end do gruLoop ! loop through GRUs

 ! * allocate local data structures where there is no spatial dimension
 select type(dataStruct)
  class is (var_i);         call allocLocal(metaStruct,dataStruct,err=err,message=cmessage)
  class is (var_i8);        call allocLocal(metaStruct,dataStruct,err=err,message=cmessage)
  class is (var_d);         call allocLocal(metaStruct,dataStruct,err=err,message=cmessage)
  class is (var_ilength);   call allocLocal(metaStruct,dataStruct,err=err,message=cmessage)
  class is (var_dlength);   call allocLocal(metaStruct,dataStruct,err=err,message=cmessage)
  ! check identified the data type
  class default; if(.not.spatial)then; err=20; message=trim(message)//'unable to identify derived data type'; return; end if
 end select

 ! error check
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

 end subroutine allocGlobal

 ! ************************************************************************************************
 ! public subroutine allocLocal: allocate space for local data structures
 ! ************************************************************************************************
 subroutine allocLocal(metaStruct,dataStruct,nSnow,nSoil,err,message)
 implicit none
 ! input-output
 type(var_info),intent(in)        :: metaStruct(:)  ! metadata structure
 class(*),intent(inout)           :: dataStruct     ! data structure
 ! optional input
 integer(i4b),intent(in),optional :: nSnow          ! number of snow layers
 integer(i4b),intent(in),optional :: nSoil          ! number of soil layers
 ! output
 integer(i4b),intent(out)         :: err            ! error code
 character(*),intent(out)         :: message        ! error message
 ! local
 logical(lgt)                     :: check          ! .true. if the variables are allocated
 integer(i4b)                     :: nVars          ! number of variables in the metadata structure
 integer(i4b)                     :: nLayers        ! total number of layers
 character(len=256)               :: cmessage       ! error message of the downwind routine
 ! initialize error control
 err=0; message='allocLocal/'

 ! get the number of variables in the metadata structure
 nVars = size(metaStruct)

 ! check if nSnow and nSoil are present
 if(present(nSnow) .or. present(nSoil))then
  ! check both are present
  if(.not.present(nSoil))then; err=20; message=trim(message)//'expect nSoil to be present when nSnow is present'; return; end if
  if(.not.present(nSnow))then; err=20; message=trim(message)//'expect nSnow to be present when nSoil is present'; return; end if
  nLayers = nSnow+nSoil

 ! It is possible that nSnow and nSoil are actually needed here, so we return an error if the optional arguments are missing when needed
 else
  select type(dataStruct)
   class is (var_flagVec); err=20
   class is (var_ilength); err=20
   class is (var_dlength); err=20
  end select
  if(err/=0)then; message=trim(message)//'expect nSnow and nSoil to be present for variable-length data structures'; return; end if
 end if

 ! initialize allocation check
 check=.false.

 ! allocate the dimension for model variables
 select type(dataStruct)
  class is (var_i);       if(allocated(dataStruct%var))then; check=.true.; else; allocate(dataStruct%var(nVars),stat=err); end if; return
  class is (var_i8);      if(allocated(dataStruct%var))then; check=.true.; else; allocate(dataStruct%var(nVars),stat=err); end if; return
  class is (var_d);       if(allocated(dataStruct%var))then; check=.true.; else; allocate(dataStruct%var(nVars),stat=err); end if; return
  class is (var_flagVec); if(allocated(dataStruct%var))then; check=.true.; else; allocate(dataStruct%var(nVars),stat=err); end if
  class is (var_ilength); if(allocated(dataStruct%var))then; check=.true.; else; allocate(dataStruct%var(nVars),stat=err); end if
  class is (var_dlength); if(allocated(dataStruct%var))then; check=.true.; else; allocate(dataStruct%var(nVars),stat=err); end if
  class default; err=20; message=trim(message)//'unable to identify derived data type for the variable dimension'; return
 end select
 ! check errors
 if(check) then; err=20; message=trim(message)//'structure was unexpectedly allocated already'; return; end if
 if(err/=0)then; err=20; message=trim(message)//'problem allocating'; return; end if

 ! allocate the dimension for model data
 select type(dataStruct)
  class is (var_flagVec); call allocateDat_flag(metaStruct,nSnow,nSoil,nLayers,dataStruct,err,cmessage)
  class is (var_ilength); call allocateDat_int( metaStruct,nSnow,nSoil,nLayers,dataStruct,err,cmessage)
  class is (var_dlength); call allocateDat_rkind(  metaStruct,nSnow,nSoil,nLayers,dataStruct,err,cmessage)
  class default; err=20; message=trim(message)//'unable to identify derived data type for the data dimension'; return
 end select

 ! check errors
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

 end subroutine allocLocal

 ! ************************************************************************************************
 ! public subroutine resizeData: resize data structure
 ! ************************************************************************************************
 subroutine resizeData(metaStruct,dataStructOrig,dataStructNew,copy,err,message)
 implicit none
 ! input
 type(var_info),intent(in)             :: metaStruct(:)  ! metadata structure
 class(*)      ,intent(in)             :: dataStructOrig ! original data structure
 ! output
 class(*)      ,intent(inout)          :: dataStructNew  ! new data structure
 ! control
 logical(lgt)  ,intent(in)   ,optional :: copy           ! flag to copy data
 integer(i4b)  ,intent(out)            :: err            ! error code
 character(*)  ,intent(out)            :: message        ! error message
 ! local
 integer(i4b)                          :: iVar           ! number of variables in the structure
 integer(i4b)                          :: nVars          ! number of variables in the structure
 logical(lgt)                          :: isCopy         ! flag to copy data (handles absence of optional argument)
 character(len=256)                    :: cmessage       ! error message of the downwind routine
 ! initialize error control
 err=0; message='resizeData/'

 ! get the copy flag
 if(present(copy))then
  isCopy = copy
 else
  isCopy = .false.
 endif

 ! get the number of variables in the data structure
 nVars = size(metaStruct)

 ! check that the input data structure is allocated
 select type(dataStructOrig)
  class is (var_ilength); err=merge(0, 20, allocated(dataStructOrig%var))
  class is (var_dlength); err=merge(0, 20, allocated(dataStructOrig%var))
  class default; err=20; message=trim(message)//'unable to identify type of data structure'; return
 end select
 if(err/=0)then; message=trim(message)//'input data structure dataStructOrig%var'; return; end if

 ! allocate the dimension for model variables
 select type(dataStructNew)
  class is (var_ilength); if(.not.allocated(dataStructNew%var)) allocate(dataStructNew%var(nVars),stat=err)
  class is (var_dlength); if(.not.allocated(dataStructNew%var)) allocate(dataStructNew%var(nVars),stat=err)
  class default; err=20; message=trim(message)//'unable to identify derived data type for the variable dimension'; return
 end select
 if(err/=0)then; message=trim(message)//'problem allocating space for dataStructNew%var'; return; end if

 ! loop through variables
 do iVar=1,nVars

  ! resize and copy data structures
  select type(dataStructOrig)

   ! double precision
   class is (var_dlength)
    select type(dataStructNew)
     class is (var_dlength); call copyStruct_rkind( dataStructOrig%var(iVar),dataStructNew%var(iVar),isCopy,err,cmessage)
     class default; err=20; message=trim(message)//'mismatch data structure for variable'//trim(metaStruct(iVar)%varname); return
    end select

   ! integer
   class is (var_ilength)
    select type(dataStructNew)
     class is (var_ilength); call copyStruct_i4b(dataStructOrig%var(iVar),dataStructNew%var(iVar),isCopy,err,cmessage)
     class default; err=20; message=trim(message)//'mismatch data structure for variable'//trim(metaStruct(iVar)%varname); return
    end select

   ! check
   class default; err=20; message=trim(message)//'unable to identify type of data structure'; return
  end select
  if(err/=0)then; message=trim(message)//trim(cmessage)//' ('//trim(metaStruct(iVar)%varname)//')'; return; end if

 end do  ! looping through variables in the data structure

 end subroutine resizeData

 ! ************************************************************************************************
 ! private subroutine copyStruct_rkind: copy a given data structure
 ! ************************************************************************************************
 subroutine copyStruct_rkind(varOrig,varNew,copy,err,message)
 ! dummy variables
 type(dlength),intent(in)    :: varOrig        ! original data structure
 type(dlength),intent(inout) :: varNew         ! new data structure
 logical(lgt) ,intent(in)    :: copy           ! flag to copy data
 integer(i4b) ,intent(out)   :: err            ! error code
 character(*) ,intent(out)   :: message        ! error message
 ! local
 logical(lgt)                :: allocatedOrig  ! flag to denote if a given variable in the original data structure is allocated
 logical(lgt)                :: allocatedNew   ! flag to denote if a given variable in the new data structure is allocated
 integer(i4b)                :: lowerBoundOrig ! lower bound of a given variable in the original data structure
 integer(i4b)                :: upperBoundOrig ! upper bound of a given variable in the original data structure
 integer(i4b)                :: lowerBoundNew  ! lower bound of a given variable in the new data structure
 integer(i4b)                :: upperBoundNew  ! upper bound of a given variable in the new data structure
 ! initialize error control
 err=0; message='copyStruct_rkind/'

 ! get the information from the data structures
 call getVarInfo(varOrig,allocatedOrig,lowerBoundOrig,upperBoundOrig)
 call getVarInfo(varNew, allocatedNew, lowerBoundNew, upperBoundNew)

 ! check that the variable of the original data structure is allocated
 if(.not.allocatedOrig)then
  message=trim(message)//'variable in the original data structure is not allocated'
  err=20; return
 endif

 ! re-size data structure if necessary
 if(lowerBoundOrig/=lowerBoundNew .or. upperBoundOrig/=upperBoundNew .or. .not.allocatedNew)then

  ! deallocate space (if necessary)
  if(allocatedNew) deallocate(varNew%dat)

  ! allocate space
  allocate(varNew%dat(lowerBoundOrig:upperBoundOrig), stat=err)
  if(err/=0)then; message=trim(message)//'problem allocating'; return; endif

 endif  ! if need to resize

 ! copy the data structure
 if(copy)then
  varNew%dat(:) = varOrig%dat(:)

 ! initialize the data structure to missing values
 else
  varNew%dat(:) = realMissing
 endif

 ! internal routines
 contains

  ! internal subroutine getVarInfo: get information from a given data structure
  subroutine getVarInfo(var,isAllocated,lowerBound,upperBound)
  ! input
  type(dlength),intent(in)         :: var            ! data vector for a given variable
  ! output
  logical(lgt),intent(out)         :: isAllocated    ! flag to denote if the data vector is allocated
  integer(i4b),intent(out)         :: lowerBound     ! lower bound
  integer(i4b),intent(out)         :: upperBound     ! upper bound
  ! local
  integer(i4b),dimension(1)        :: lowerBoundVec  ! lower bound vector
  integer(i4b),dimension(1)        :: upperBoundVec  ! upper bound vector
  ! initialize error control
  err=0; message='getVarInfo/'

  ! check that the input data structure is allocated
  isAllocated = allocated(var%dat)

  ! if allocated then get the bounds
  ! NOTE: also convert vector to scalar
  if(isAllocated)then
   lowerBoundVec=lbound(var%dat); lowerBound=lowerBoundVec(1)
   upperBoundVec=ubound(var%dat); upperBound=upperBoundVec(1)

  ! if not allocated then return zero bounds
  else
   lowerBound=0
   upperBound=0
  endif ! (check allocation)

  end subroutine getVarInfo

 end subroutine copyStruct_rkind

 ! ************************************************************************************************
 ! private subroutine copyStruct_i4b: copy a given data structure
 ! ************************************************************************************************
 subroutine copyStruct_i4b(varOrig,varNew,copy,err,message)
 ! dummy variables
 type(ilength),intent(in)    :: varOrig        ! original data structure
 type(ilength),intent(inout) :: varNew         ! new data structure
 logical(lgt) ,intent(in)    :: copy           ! flag to copy data
 integer(i4b) ,intent(out)   :: err            ! error code
 character(*) ,intent(out)   :: message        ! error message
 ! local
 logical(lgt)                :: allocatedOrig  ! flag to denote if a given variable in the original data structure is allocated
 logical(lgt)                :: allocatedNew   ! flag to denote if a given variable in the new data structure is allocated
 integer(i4b)                :: lowerBoundOrig ! lower bound of a given variable in the original data structure
 integer(i4b)                :: upperBoundOrig ! upper bound of a given variable in the original data structure
 integer(i4b)                :: lowerBoundNew  ! lower bound of a given variable in the new data structure
 integer(i4b)                :: upperBoundNew  ! upper bound of a given variable in the new data structure
 ! initialize error control
 err=0; message='copyStruct_i4b/'

 ! get the information from the data structures
 call getVarInfo(varOrig,allocatedOrig,lowerBoundOrig,upperBoundOrig)
 call getVarInfo(varNew, allocatedNew, lowerBoundNew, upperBoundNew)

 ! check that the variable of the original data structure is allocated
 if(.not.allocatedOrig)then
  message=trim(message)//'variable in the original data structure is not allocated'
  err=20; return
 endif

 ! re-size data structure if necessary
 if(lowerBoundOrig/=lowerBoundNew .or. upperBoundOrig/=upperBoundNew .or. .not.allocatedNew)then

  ! deallocate space (if necessary)
  if(allocatedNew) deallocate(varNew%dat)

  ! allocate space
  allocate(varNew%dat(lowerBoundOrig:upperBoundOrig), stat=err)
  if(err/=0)then; message=trim(message)//'problem allocating'; return; endif

 endif  ! if need to resize

 ! copy the data structure
 if(copy)then
  varNew%dat(:) = varOrig%dat(:)

 ! initialize the data structure to missing values
 else
  varNew%dat(:) = integerMissing
 endif

 ! internal routines
 contains

  ! internal subroutine getVarInfo: get information from a given data structure
  subroutine getVarInfo(var,isAllocated,lowerBound,upperBound)
  ! input
  type(ilength),intent(in)         :: var            ! data vector for a given variable
  ! output
  logical(lgt),intent(out)         :: isAllocated    ! flag to denote if the data vector is allocated
  integer(i4b),intent(out)         :: lowerBound     ! lower bound
  integer(i4b),intent(out)         :: upperBound     ! upper bound
  ! local
  integer(i4b),dimension(1)        :: lowerBoundVec  ! lower bound vector
  integer(i4b),dimension(1)        :: upperBoundVec  ! upper bound vector
  ! initialize error control
  err=0; message='getVarInfo/'

  ! check that the input data structure is allocated
  isAllocated = allocated(var%dat)

  ! if allocated then get the bounds
  ! NOTE: also convert vector to scalar
  if(isAllocated)then
   lowerBoundVec=lbound(var%dat); lowerBound=lowerBoundVec(1)
   upperBoundVec=ubound(var%dat); upperBound=upperBoundVec(1)

  ! if not allocated then return zero bounds
  else
   lowerBound=0
   upperBound=0
  endif ! (check allocation)

  end subroutine getVarInfo

 end subroutine copyStruct_i4b


 ! ************************************************************************************************
 ! private subroutine allocateDat_rkind: initialize data dimension of the data structures
 ! ************************************************************************************************
 subroutine allocateDat_rkind(metadata,nSnow,nSoil,nLayers, & ! input
                           varData,err,message)            ! output
 ! access subroutines
 USE get_ixName_module,only:get_varTypeName       ! to access type strings for error messages

 implicit none
 ! input variables
 type(var_info),intent(in)         :: metadata(:) ! metadata structure
 integer(i4b),intent(in)           :: nSnow       ! number of snow layers
 integer(i4b),intent(in)           :: nSoil       ! number of soil layers
 integer(i4b),intent(in)           :: nLayers     ! total number of soil layers in the snow+soil domian (nSnow+nSoil)
 ! output variables
 type(var_dlength),intent(inout)   :: varData     ! model variables for a local HRU
 integer(i4b),intent(out)          :: err         ! error code
 character(*),intent(out)          :: message     ! error message
 ! local variables
 integer(i4b)                      :: iVar        ! variable index
 integer(i4b)                      :: nVars       ! number of variables in the metadata structure

! initialize error control
 err=0; message='allocateDat_rkind/'

 ! get the number of variables in the metadata structure
 nVars = size(metadata)

 ! loop through variables in the data structure
 do iVar=1,nVars

  ! check allocated
  if(allocated(varData%var(iVar)%dat))then
   message=trim(message)//'variable '//trim(metadata(iVar)%varname)//' is unexpectedly allocated'
   err=20; return

  ! allocate structures
  ! NOTE: maxvarFreq is the number of possible output frequencies
  !        -- however, this vector must store two values for the variance calculation, thus the *2 in this allocate
  !            (need enough space in the event that variance is the desired statistic for all output frequencies)
  else
   select case(metadata(iVar)%vartype)
    case(iLookVarType%scalarv); allocate(varData%var(iVar)%dat(1),stat=err)
    case(iLookVarType%wLength); allocate(varData%var(iVar)%dat(nBand),stat=err)
    case(iLookVarType%midSnow); allocate(varData%var(iVar)%dat(nSnow),stat=err)
    case(iLookVarType%midSoil); allocate(varData%var(iVar)%dat(nSoil),stat=err)
    case(iLookVarType%midToto); allocate(varData%var(iVar)%dat(nLayers),stat=err)
    case(iLookVarType%ifcSnow); allocate(varData%var(iVar)%dat(0:nSnow),stat=err)
    case(iLookVarType%ifcSoil); allocate(varData%var(iVar)%dat(0:nSoil),stat=err)
    case(iLookVarType%ifcToto); allocate(varData%var(iVar)%dat(0:nLayers),stat=err)
    case(iLookVarType%parSoil); allocate(varData%var(iVar)%dat(nSoil),stat=err)
    case(iLookVarType%routing); allocate(varData%var(iVar)%dat(nTimeDelay),stat=err)
    case(iLookVarType%outstat); allocate(varData%var(iVar)%dat(maxvarfreq*2),stat=err)
    case(iLookVarType%unknown); allocate(varData%var(iVar)%dat(0),stat=err)  ! unknown = special (and valid) case that is allocated later (initialize with zero-length vector)
    case default
     err=40; message=trim(message)//"1. unknownVariableType[name='"//trim(metadata(iVar)%varname)//"'; type='"//trim(get_varTypeName(metadata(iVar)%vartype))//"']"
     return
   end select
   ! check error
   if(err/=0)then; err=20; message=trim(message)//'problem allocating variable '//trim(metadata(iVar)%varname); return; end if
   ! set to missing
   varData%var(iVar)%dat(:) = realMissing
  end if  ! if not allocated

 end do  ! looping through variables

 end subroutine allocateDat_rkind

 ! ************************************************************************************************
 ! private subroutine allocateDat_int: initialize data dimension of the data structures
 ! ************************************************************************************************
 subroutine allocateDat_int(metadata,nSnow,nSoil,nLayers, & ! input
                            varData,err,message)            ! output
 USE get_ixName_module,only:get_varTypeName       ! to access type strings for error messages
 implicit none
 ! input variables
 type(var_info),intent(in)         :: metadata(:) ! metadata structure
 integer(i4b),intent(in)           :: nSnow       ! number of snow layers
 integer(i4b),intent(in)           :: nSoil       ! number of soil layers
 integer(i4b),intent(in)           :: nLayers     ! total number of soil layers in the snow+soil domian (nSnow+nSoil)
 ! output variables
 type(var_ilength),intent(inout)   :: varData     ! model variables for a local HRU
 integer(i4b),intent(out)          :: err         ! error code
 character(*),intent(out)          :: message     ! error message
 ! local variables
 integer(i4b)                      :: iVar        ! variable index
 integer(i4b)                      :: nVars       ! number of variables in the metadata structure

! initialize error control
 err=0; message='allocateDat_int/'

 ! get the number of variables in the metadata structure
 nVars = size(metadata)

! loop through variables in the data structure
 do iVar=1,nVars

  ! check allocated
  if(allocated(varData%var(iVar)%dat))then
   message=trim(message)//'variable '//trim(metadata(iVar)%varname)//' is unexpectedly allocated'
   err=20; return

  ! allocate structures
  ! NOTE: maxvarFreq is the number of possible output frequencies
  !        -- however, this vector must store two values for the variance calculation, thus the *2 in this allocate
  !            (need enough space in the event that variance is the desired statistic for all output frequencies)
  else
   select case(metadata(iVar)%vartype)
    case(iLookVarType%scalarv); allocate(varData%var(iVar)%dat(1),stat=err)
    case(iLookVarType%wLength); allocate(varData%var(iVar)%dat(nBand),stat=err)
    case(iLookVarType%midSnow); allocate(varData%var(iVar)%dat(nSnow),stat=err)
    case(iLookVarType%midSoil); allocate(varData%var(iVar)%dat(nSoil),stat=err)
    case(iLookVarType%midToto); allocate(varData%var(iVar)%dat(nLayers),stat=err)
    case(iLookVarType%ifcSnow); allocate(varData%var(iVar)%dat(0:nSnow),stat=err)
    case(iLookVarType%ifcSoil); allocate(varData%var(iVar)%dat(0:nSoil),stat=err)
    case(iLookVarType%ifcToto); allocate(varData%var(iVar)%dat(0:nLayers),stat=err)
    case(iLookVarType%routing); allocate(varData%var(iVar)%dat(nTimeDelay),stat=err)
    case(iLookVarType%outstat); allocate(varData%var(iVar)%dat(maxvarFreq*2),stat=err)
    case(iLookVarType%unknown); allocate(varData%var(iVar)%dat(0),stat=err)  ! unknown=special (and valid) case that is allocated later (initialize with zero-length vector)
    case default; err=40; message=trim(message)//"unknownVariableType[name='"//trim(metadata(iVar)%varname)//"'; type='"//trim(get_varTypeName(metadata(iVar)%vartype))//"']"; return
   end select
   ! check error
   if(err/=0)then; err=20; message=trim(message)//'problem allocating variable '//trim(metadata(iVar)%varname); return; end if
   ! set to missing
   varData%var(iVar)%dat(:) = integerMissing
  end if  ! if not allocated

 end do  ! looping through variables

 end subroutine allocateDat_int

 ! ************************************************************************************************
 ! private subroutine allocateDat_flag: initialize data dimension of the data structures
 ! ************************************************************************************************
 subroutine allocateDat_flag(metadata,nSnow,nSoil,nLayers, & ! input
                             varData,err,message)            ! output
 USE get_ixName_module,only:get_varTypeName       ! to access type strings for error messages
 implicit none
 ! input variables
 type(var_info),intent(in)         :: metadata(:) ! metadata structure
 integer(i4b),intent(in)           :: nSnow       ! number of snow layers
 integer(i4b),intent(in)           :: nSoil       ! number of soil layers
 integer(i4b),intent(in)           :: nLayers     ! total number of soil layers in the snow+soil domian (nSnow+nSoil)
 ! output variables
 type(var_flagVec),intent(inout)   :: varData     ! model variables for a local HRU
 integer(i4b),intent(out)          :: err         ! error code
 character(*),intent(out)          :: message     ! error message
 ! local variables
 integer(i4b)                      :: iVar        ! variable index
 integer(i4b)                      :: nVars       ! number of variables in the metadata structure

! initialize error control
 err=0; message='allocateDat_flag/'

 ! get the number of variables in the metadata structure
 nVars = size(metadata)

! loop through variables in the data structure
 do iVar=1,nVars

  ! check allocated
  if(allocated(varData%var(iVar)%dat))then
   message=trim(message)//'variable '//trim(metadata(iVar)%varname)//' is unexpectedly allocated'
   err=20; return

  ! allocate structures
  ! NOTE: maxvarFreq is the number of possible output frequencies
  !        -- however, this vector must store two values for the variance calculation, thus the *2 in this allocate
  !            (need enough space in the event that variance is the desired statistic for all output frequencies)
  else
   select case(metadata(iVar)%vartype)
    case(iLookVarType%scalarv); allocate(varData%var(iVar)%dat(1),stat=err)
    case(iLookVarType%wLength); allocate(varData%var(iVar)%dat(nBand),stat=err)
    case(iLookVarType%midSnow); allocate(varData%var(iVar)%dat(nSnow),stat=err)
    case(iLookVarType%midSoil); allocate(varData%var(iVar)%dat(nSoil),stat=err)
    case(iLookVarType%midToto); allocate(varData%var(iVar)%dat(nLayers),stat=err)
    case(iLookVarType%ifcSnow); allocate(varData%var(iVar)%dat(0:nSnow),stat=err)
    case(iLookVarType%ifcSoil); allocate(varData%var(iVar)%dat(0:nSoil),stat=err)
    case(iLookVarType%ifcToto); allocate(varData%var(iVar)%dat(0:nLayers),stat=err)
    case(iLookVarType%routing); allocate(varData%var(iVar)%dat(nTimeDelay),stat=err)
    case(iLookVarType%outstat); allocate(varData%var(iVar)%dat(maxvarFreq*2),stat=err)
    case(iLookVarType%unknown); allocate(varData%var(iVar)%dat(0),stat=err)  ! unknown=special (and valid) case that is allocated later (initialize with zero-length vector)
    case default; err=40; message=trim(message)//"unknownVariableType[name='"//trim(metadata(iVar)%varname)//"'; type='"//trim(get_varTypeName(metadata(iVar)%vartype))//"']"; return
   end select
   ! check error
   if(err/=0)then; err=20; message=trim(message)//'problem allocating variable '//trim(metadata(iVar)%varname); return; end if
   ! set to false
   varData%var(iVar)%dat(:) = .false.
  end if  ! if not allocated

 end do  ! looping through variables

 end subroutine allocateDat_flag

end module allocspace_module
