! **************************************************************
!
! This file contains the subroutines: contacts,c_alfa,c_cont
!
! Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
!                      Shura Hayryan, Chin-Ku 
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
!
! **************************************************************
      
       subroutine contacts(ncn,nham2,dham)

! ..............................................................
!
! Calculates number of contacts in given conformation, number of
! contacts which are the same in given and reference conformation,
! and the hamming distance between given conformation and the 
! reference conformation.
!
! CALLS: c_cont
! ..............................................................

      include 'INCL.H'

!f2py integer, intent(out) :: ncn, nham2 
!f2py double precistion, intent(out) :: dham

      call c_cont(1,ncode)

      ncn=0
      nham=0
      nham2=0
      nresi=irsml2(1)-irsml1(1)+1
      do i=1,nresi
        do j=nresi,i+3,-1

        if (ijcont(i,j).eq.1) then
          ncn=ncn+1
          if (iref(i,j).eq.1) nham2 = nham2+1
        end if

        nham = nham + abs(ijcont(i,j)-iref(i,j))
       end do
      end do

      if (ncn.ne.0.and.nci.ne.0) then
        dham = float(nham)/float(ncn)/float(nci)
      else
        dham = 1.0
      end if

      return
      end 


! ********************************* 
      subroutine c_alfa(nmol,ncode)

! ......................................................
!    Calculates the indices of C-alpha atoms and 
!    stores in the array ind_alf(mxrs)
!                        
!    Usage: call c_alfa(nmol,ncode)
!
!           nmol - index of the molecule
!           ncode ---> not in use in the current version
!
!    OUTPUT:  ind_alf(mxrs)
!
! CALLS: none
! ......................................................

      include 'INCL.H'
              
      do n_res=irsml1(nmol),irsml2(nmol) ! Over res. 
        do ia=iatrs1(n_res),iatrs2(n_res) ! Over the atoms of res. 

!     Check for C_alpha atoms

         if (nmat(ia)(1:2).eq.'ca') then
           ind_alf(n_res)=ia
         endif

       enddo ! Over the atoms of res. 
      enddo ! Over the res. 
 
      return
      end

! **********************************
      subroutine c_cont (nmol,ncode)

!..............................................................
!  Calculates the matrix of contacts between aminoacid residues 
!  of the molecule "nmol" according to  L.Mirny and E.Domany, 
!  PROTEINS:Structure, Function, and Genetics 26:391-410 (1996)
!                
!  Two residues are in contact if their C_alpha atoms are
!  closer than 8.5 Angstrem
!
!  Usage: call c_cont(nmol,ncode)
!
!       Where nmol is the index of the molecule (always 1, in the 
!       current version of SMM)
!       ncode ---> not in use in the current version
!
!  IMPORTANT: Before the first call of this subroutine  "c_alfa"
!          must be called to calculate the inices of C_alpha atoms.
!          (ONLY ONCE)
!
!   OUTPUT: The output of this routine is the contact matrix
!          ijcont(mxrs,mxrs) 
!
!              ijcont(i,j)=0---> residues i and j are not in contact
!              ijcont(i,j)=1---> ---------''----- are in contact
!              ijcont(i,j)=2---> residues i and j are adjacent
!
!    NOTE:  Adjacent residues are always in contact (and therefore not
!           counted)
!
!         Here "mxrs" is the maximum number of residues for SMM
!         Obviously, this subroutine calculates only NxN part
!         of that matrix, N -is the number of res. in "nmol"
! 
! CALLS:  none
!..............................................................

       include 'INCL.H'

       rcut=8.5   ! Domany
              
              
       do nr_i=irsml1(nmol),irsml2(nmol) ! Over res. i

          ijcont(nr_i,nr_i)=2
          if(nr_i+1.le.irsml2(nmol)) then
              ijcont(nr_i,nr_i+1)=2
              ijcont(nr_i+1,nr_i)=2
              if(nr_i+2.le.irsml2(nmol)) then
                 ijcont(nr_i,nr_i+2)=2
                 ijcont(nr_i+2,nr_i)=2
               end if
           end if

           do nr_j=nr_i+3,irsml2(nmol) ! Over res. j 

!             write(*,'(2i3)'),nr_i,nr_j

              ic=0

              ialf=ind_alf(nr_i)
              jalf=ind_alf(nr_j)

              rij2=(xat(ialf)-xat(jalf))**2
     &              +(yat(ialf)-yat(jalf))**2
     &                   + (zat(ialf)-zat(jalf))**2
              if(sqrt(rij2).lt.rcut) ic=1

!             write(*,'(2i3)'),nr_i,nr_j

              ijcont(nr_i,nr_j)=ic
              ijcont(nr_j,nr_i)=ic ! The matrix is symmetrical

           end do ! Over res. j
       end do ! Over res. i

       return
       end

