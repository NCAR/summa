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

MODULE data_types
 ! used to define model data structures
 USE nrtype, integerMissing=>nr_integerMissing
 USE var_lookup,only:maxvarFreq
 USE var_lookup,only:maxvarStat
 implicit none
 ! constants necessary for variable defs
 private

 ! ***********************************************************************************************************
 ! Define the model decisions
 ! ***********************************************************************************************************
 ! the model decision structure
 type,public  :: model_options
  character(len=64)                      :: cOption   = 'notPopulatedYet'
  character(len=64)                      :: cDecision = 'notPopulatedYet'
  integer(i4b)                           :: iDecision = integerMissing
 end type model_options

 ! ***********************************************************************************************************
 ! Define metadata for model forcing datafile
 ! ***********************************************************************************************************
 ! define a derived type for the data in the file
 type,public  :: file_info
  character(len=256)                     :: filenmData='notPopulatedYet' ! name of data file
  integer(i4b)                           :: nVars                    ! number of variables in the file
  integer(i4b)                           :: nTimeSteps               ! number of variables in the file
  integer(i4b),allocatable               :: var_ix(:)                ! index of each forcing data variable in the data structure
  integer(i4b),allocatable               :: data_id(:)               ! netcdf variable id for each forcing data variable
  character(len=256),allocatable         :: varName(:)               ! netcdf variable name for each forcing data variable
  real(rkind)                               :: firstJulDay              ! first julian day in forcing file
  real(rkind)                               :: convTime2Days            ! factor to convert time to days
 end type file_info

 ! ***********************************************************************************************************
 ! Define metadata on model parameters
 ! ***********************************************************************************************************
 ! define a data type to store model parameter information
 type,public  :: par_info
  real(rkind)                               :: default_val              ! default parameter value
  real(rkind)                               :: lower_limit              ! lower bound
  real(rkind)                               :: upper_limit              ! upper bound
 endtype par_info

 ! ***********************************************************************************************************
 ! Define variable metadata
 ! ***********************************************************************************************************
 ! define derived type for model variables, including name, description, and units
 type,public :: var_info
  character(len=64)                      :: varname   = 'empty'         ! variable name
  character(len=128)                     :: vardesc   = 'empty'         ! variable description
  character(len=64)                      :: varunit   = 'empty'         ! variable units
  integer(i4b)                           :: vartype   = integerMissing  ! variable type
  integer(i4b),dimension(maxvarFreq)     :: ncVarID   = integerMissing  ! netcdf variable id (missing if frequency is not desired)
  integer(i4b),dimension(maxvarFreq)     :: statIndex = integerMissing  ! index of desired statistic for temporal aggregation
  logical(lgt)                           :: varDesire = .false.         ! flag to denote if the variable is desired for model output
 endtype var_info

 ! define extended data type (include indices to map onto parent data type)
 type,extends(var_info),public :: extended_info
  integer(i4b)                           :: ixParent         ! index in the parent data structure
 endtype extended_info

 ! define extended data type (includes named variables for the states affected by each flux)
 type,extends(var_info),public :: flux2state
  integer(i4b)                           :: state1           ! named variable of the 1st state affected by the flux
  integer(i4b)                           :: state2           ! named variable of the 2nd state affected by the flux
 endtype flux2state

 ! ***********************************************************************************************************
 ! Define summary of data structures
 ! ***********************************************************************************************************
 ! data structure information
 type,public :: struct_info
  character(len=32)                      :: structName  ! name of the data structure
  character(len=32)                      :: lookName    ! name of the look-up variables
  integer(i4b)                           :: nVar        ! number of variables in each data structure
 end type struct_info

 ! ***********************************************************************************************************
 ! Define data types to map between GRUs and HRUs
 ! ***********************************************************************************************************

 ! hru info data structure
 type, public :: hru_info
  integer(i4b)                      :: hru_nc                   ! index of the hru in the netcdf file
  integer(i4b)                      :: hru_ix                   ! index of the hru in the run domain
  integer(8)                        :: hru_id                   ! id (non-sequential number) of the hru
  integer(i4b)                      :: nSnow                    ! number of snow layers
  integer(i4b)                      :: nSoil                    ! number of soil layers
 endtype hru_info

 ! define mapping from GRUs to the HRUs
 type, public :: gru2hru_map
  integer(8)                        :: gru_id                   ! id of the gru
  integer(i4b)                      :: hruCount                 ! total number of hrus in the gru
  type(hru_info), allocatable       :: hruInfo(:)               ! basic information of HRUs within the gru
  integer(i4b)                      :: gru_nc                   ! index of gru in the netcdf file
 endtype gru2hru_map

 ! define the mapping from the HRUs to the GRUs
 type, public :: hru2gru_map
  integer(i4b)                      :: gru_ix                   ! index of gru which the hru belongs to
  integer(i4b)                      :: localHRU_ix              ! index of a hru within a gru (start from 1 per gru)
 endtype hru2gru_map

 ! ***********************************************************************************************************
 ! Define hierarchal derived data types
 ! ***********************************************************************************************************
  ! define derived types to hold look-up tables for each soil layer
 ! ** double precision type
 type, public :: dLookup
  real(dp),allocatable                :: lookup(:)   ! lookup(:)
 endtype dLookup
 ! ** double precision type for a variable number of soil layers; variable length
 type, public :: vLookup
  type(dLookup),allocatable           :: var(:)      ! var(:)%lookup(:)
 endtype vLookup
 ! ** double precision type for a variable number of soil layers
 type, public :: zLookup
  type(vLookup),allocatable           :: z(:)        ! z(:)%var(:)%lookup(:)
 endtype zLookup
 ! ** double precision type for a variable number of soil layers
 type, public :: hru_z_vLookup
  type(zLookup),allocatable           :: hru(:)      ! hru(:)%z(:)%var(:)%lookup(:)
 endtype hru_z_vLookup
 ! ** double precision type for a variable number of soil layers
 type, public :: gru_hru_z_vLookup
  type(hru_z_vLookup),allocatable     :: gru(:)      ! gru(:)%hru(:)%z(:)%var(:)%lookup(:)
 endtype gru_hru_z_vLookup
 ! define derived types to hold multivariate data for a single variable (different variables have different length)
 ! NOTE: use derived types here to facilitate adding the "variable" dimension
 ! ** double precision type
 type, public :: dlength
  real(rkind),allocatable                :: dat(:)    ! dat(:)
 endtype dlength
 ! ** integer type (4 byte)
 type, public :: ilength
  integer(i4b),allocatable            :: dat(:)    ! dat(:)
 endtype ilength
 ! ** integer type (8 byte)
 type, public :: i8length
  integer(8),allocatable              :: dat(:)    ! dat(:)
 endtype i8length
 ! ** logical type
 type, public :: flagVec
  logical(lgt),allocatable            :: dat(:)    ! dat(:)
 endtype flagVec

 ! define derived types to hold data for multiple variables
 ! NOTE: use derived types here to facilitate adding extra dimensions (e.g., spatial)

 ! ** double precision type of variable length
 type, public :: var_dlength
  type(dlength),allocatable           :: var(:)    ! var(:)%dat
 endtype var_dlength
 ! ** integer type of variable length (4 byte)
 type, public :: var_ilength
  type(ilength),allocatable           :: var(:)    ! var(:)%dat
 endtype var_ilength
 ! ** integer type of variable length (8 byte)
 type, public :: var_i8length
  type(i8length),allocatable          :: var(:)    ! var(:)%dat
 endtype var_i8length
 ! ** logical type of variable length
 type, public :: var_flagVec
  type(flagVec),allocatable           :: var(:)    ! var(:)%dat
 endtype var_flagVec

 ! ** double precision type of fixed length
 type, public :: var_d
  real(rkind),allocatable                :: var(:)    ! var(:)
 endtype var_d
 ! ** integer type of fixed length (4 byte)
 type, public :: var_i
  integer(i4b),allocatable            :: var(:)    ! var(:)
 endtype var_i
 ! ** integer type of fixed length (8 byte)
 type, public :: var_i8
  integer(8),allocatable              :: var(:)    ! var(:)
 endtype var_i8

 ! ** double precision type of fixed length
 type, public :: hru_d
  real(rkind),allocatable                :: hru(:)    ! hru(:)
 endtype hru_d
 ! ** integer type of fixed length (4 byte)
 type, public :: hru_i
  integer(i4b),allocatable            :: hru(:)    ! hru(:)
 endtype hru_i
 ! ** integer type of fixed length (8 byte)
 type, public :: hru_i8
  integer(8),allocatable              :: hru(:)    ! hru(:)
 endtype hru_i8

 ! define derived types to hold JUST the HRU dimension
 ! ** double precision type of variable length
 type, public :: hru_doubleVec
  type(var_dlength),allocatable      :: hru(:)     ! hru(:)%var(:)%dat
 endtype hru_doubleVec
 ! ** integer type of variable length (4 byte)
 type, public :: hru_intVec
  type(var_ilength),allocatable      :: hru(:)     ! hru(:)%var(:)%dat
 endtype hru_intVec
 ! ** integer type of variable length (8 byte)
 type, public :: hru_int8Vec
  type(var_i8length),allocatable     :: hru(:)     ! hru(:)%var(:)%dat
 endtype hru_int8Vec
 ! ** double precision type of fixed length
 type, public :: hru_double
  type(var_d),allocatable            :: hru(:)     ! hru(:)%var(:)
 endtype hru_double
 ! ** integer type of fixed length (4 byte)
 type, public :: hru_int
  type(var_i),allocatable            :: hru(:)     ! hru(:)%var(:)
 endtype hru_int
 ! ** integer type of fixed length (8 byte)
 type, public :: hru_int8
  type(var_i8),allocatable           :: hru(:)     ! hru(:)%var(:)
 endtype hru_int8

 ! define derived types to hold JUST the HRU dimension
 ! ** double precision type of variable length
 type, public :: gru_doubleVec
  type(var_dlength),allocatable      :: gru(:)     ! gru(:)%var(:)%dat
 endtype gru_doubleVec
 ! ** integer type of variable length (4 byte)
 type, public :: gru_intVec
  type(var_ilength),allocatable      :: gru(:)     ! gru(:)%var(:)%dat
 endtype gru_intVec
 ! ** integer type of variable length (8 byte)
 type, public :: gru_int8Vec
  type(var_i8length),allocatable     :: gru(:)     ! gru(:)%var(:)%dat
 endtype gru_int8Vec
 ! ** double precision type of fixed length
 type, public :: gru_double
  type(var_d),allocatable            :: gru(:)     ! gru(:)%var(:)
 endtype gru_double
 ! ** integer type of variable length (4 byte)
 type, public :: gru_int
  type(var_i),allocatable            :: gru(:)     ! gru(:)%var(:)
 endtype gru_int
 ! ** integer type of variable length (8 byte)
 type, public :: gru_int8
  type(var_i8),allocatable           :: gru(:)     ! gru(:)%var(:)
 endtype gru_int8

 ! define derived types to hold BOTH the GRU and HRU dimension
 ! ** double precision type of variable length
 type, public :: gru_hru_doubleVec
  type(hru_doubleVec),allocatable    :: gru(:)     ! gru(:)%hru(:)%var(:)%dat
 endtype gru_hru_doubleVec
 ! ** integer type of variable length (4 byte)
 type, public :: gru_hru_intVec
  type(hru_intVec),allocatable       :: gru(:)     ! gru(:)%hru(:)%var(:)%dat
 endtype gru_hru_intVec
 ! ** integer type of variable length (8 byte)
 type, public :: gru_hru_int8Vec
  type(hru_int8Vec),allocatable      :: gru(:)     ! gru(:)%hru(:)%var(:)%dat
 endtype gru_hru_int8Vec
 ! ** double precision type of fixed length
 type, public :: gru_hru_double
  type(hru_double),allocatable       :: gru(:)     ! gru(:)%hru(:)%var(:)
 endtype gru_hru_double
 ! ** integer type of variable length (4 byte)
 type, public :: gru_hru_int
  type(hru_int),allocatable          :: gru(:)     ! gru(:)%hru(:)%var(:)
 endtype gru_hru_int
 ! ** integer type of variable length (8 byte)
 type, public :: gru_hru_int8
  type(hru_int8),allocatable         :: gru(:)     ! gru(:)%hru(:)%var(:)
 endtype gru_hru_int8
 ! ** double precision type of fixed length
 type, public :: gru_d
  type(hru_d),allocatable            :: gru(:)    ! gru(:)%hru(:)
 endtype gru_d
 ! ** integer type of fixed length
 type, public :: gru_i
  type(hru_i),allocatable            :: gru(:)    ! gru(:)%hru(:)
 endtype gru_i

END MODULE data_types

