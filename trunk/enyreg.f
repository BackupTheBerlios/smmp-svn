! This file contains enyreg
!
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
! *******************************
      real*8 function enyreg(nml)

! ----------------------------------------------------
!
! PURPOSE: sum( ( R_i - R^ref_j )**2 )
!
!    with: R_i     - atom position i in SMMP structure
!          R^ref_j - corresponding atom j in PDB str.
!
! CALLS: none
!
! ----------------------------------------------------

      include 'INCL.H'
      include 'INCP.H'


      double precision eny

      integer i, nml, j

      eny = 0.d0

      do i=iatrs1(irsml1(nml)),iatrs2(irsml2(nml))

        j=ixatp(i)
        if (j.gt.0) then  ! corresp. atom in ref. structure

          eny = eny + (xat(i)-xatp(j))**2+(yat(i)-yatp(j))**2+
     &                (zat(i)-zatp(j))**2
        endif
      enddo   ! atoms

      enyreg = eny

      return
      end

