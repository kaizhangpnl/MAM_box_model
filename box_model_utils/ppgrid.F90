! ppgrid.F90
!    This F90 module file is a special version of the equivalent ACME (and CAM5) module.
!    It provides the functionality needed by the cambox offline code
!    that is used for development and testing of the modal aerosol module (MAM),
!    but (in most cases) not all the functionality of the equivalent ACME module.
!    Also, it may have been taken from a version of CAM5 that was older
!    than ACME-V0 (i.e., pre 2014).

      module ppgrid

      implicit none

      public

      integer, parameter :: pcols = PCOLS
      integer, parameter :: pver  = PVER
      integer, parameter :: psubcols = 1

      integer, parameter :: begchunk = 1
      integer, parameter :: endchunk = 2

      end module ppgrid
