! **************************************************************
!
! This file contains the subroutines: mklist,quench
!
! Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
!                      Shura Hayryan, Chin-Ku 
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
!
! **************************************************************

      subroutine mklist(nml)

! ......................................................................
! PURPOSE: Compile interaction lists ('1-4' according to Scheraga)
!
! CALLS: quench
! ......................................................................
! TODO: Calculate van-der-Waals regions over all molecules.
      include 'INCL.H'
      integer nml

      integer mxh, mx2, l1st1, l1st2, l2nd1, l2nd2, l1i, l2i
      parameter (mxh=50,         ! max. # of atom regions
     &           mx2=50)

      logical ovlp,quench

      dimension l1st1(mxh),l1st2(mxh),l2nd1(mxh),l2nd2(mxh)
     &         ,l1i(mxbd),l2i(mx2)

      integer ia, i, i1s, i2s, ibd, ib, ifivr, ifiat, ilaat, im, ilavr,
     &        iob, io, ioiob, iow, it, is, iv, jbd, j, n1i, n14, n1st,
     &        n2nd, n2i, ntlms, nvw

! _______________________ indices of 1st vdw-region/14-partner for 'nml'
      if (nml.eq.1) then
        ivwml1(1)=1
        i14ml1(1)=1
      else
        ivwml1(nml)=ivwml1(nml-1)+nvwml(nml-1)
        i14ml1(nml)=i14ml1(nml-1)+n14ml(nml-1)
      endif

      ntlms=nmsml(nml)
      if (ntlms.eq.0) then
        write (*,'(a,i4)')
     &           ' mklist> No mov. sets defined in molecule #',nml
        nvwml(nml)=0
        n14ml(nml)=0
        return
      endif

      nvw=ivwml1(nml)-1     ! # of vdw-regions we have so far
      n14=i14ml1(nml)-1     ! # of 14-partners      -"-
! First atom in molecule
      ifiat=iatrs1(irsml1(nml))
! Last atom in molecule
      ilaat=iatrs2(irsml2(nml))
! First variable in molecule
      ifivr=ivrml1(nml)
! Last variable in molecule
      ilavr=ifivr+nvrml(nml)-1
! ____________________________ initialize: 1st vdw-region & 14-partner per atom 
      do i=ifiat,ilaat
        ivwat1(i)=0         !!!  for some atoms ...
        i14at1(i)=0         !!!  ... remains = 0
      enddo
      n1st=1                ! initialize 1ST list of interact. partners:
      l1st1(1)=ifiat        !   one region including ALL atoms
      l1st2(1)=ilaat        !   of molecule 'nml'

      i1s=imsml1(nml)+ntlms ! 1st mov.set of molecule 'nml+1'

      do io=ilavr,ifivr,-1  ! ====== from last -> first variable in 'nml'
        iv=iorvr(io)        ! ====== according to 'descendent' order
        it=ityvr(iv)  ! type of var.
        i2s=i1s-1
        i1s=imsvr1(iv)
        if ((i2s-i1s+1).gt.0) then

! ____________ exclude mov.sets of var. 'iv' from 1ST list of interact.partn.
          do is=i1s,i2s
            ovlp=quench(latms1(is),latms2(is),n1st,mxh,l1st1,l1st2)
          enddo
! _______________________________ intitialize 2ND list with current 1ST list
          do i=1,n1st
            l2nd1(i)=l1st1(i)
            l2nd2(i)=l1st2(i)
          enddo
          n2nd=n1st
! _________________________________ exclude 'ib' of var. 'iv' from 2ND list
          ib=iowat(iatvr(iv))
          ovlp=quench(ib,ib,n2nd,mxh,l2nd1,l2nd2)

          ovlp=.false.
          iob=iowat(ib)
          n2i=0
          if (iob.gt.0) then                          ! 'iob' exists
            ovlp=quench(iob,iob,n2nd,mxh,l2nd1,l2nd2) ! & in 2ND list

! _____ atoms branching from 'iob': into GENERAL list of 1-4 partners
            do i=1,nbdat(iob)
              ibd=ibdat(i,iob)
              if (ibd.ne.ib.and.iowat(ibd).eq.iob.and.
     &            quench(ibd,ibd,n2nd,mxh,l2nd1,l2nd2) ) then
                n2i=n2i+1
                if (n2i.gt.mx2) then
                  write (*,'(a,i3,2a)')  ' mklist> Molecule # ',nml,
     &                ': too many atoms bound to ',nmat(iob)
                  stop
                endif
                l2i(n2i)=ibd
              endif
            enddo  ! ... branches of 'iob'
! ____________________________ check for further '1-4' partners
!                              connected to branches 'l2i'
            do i=1,n2i
              ia=l2i(i)
              im=ixmsat(ia)
              if (im.gt.0) then
                do j=latms1(im),latms2(im)
                  if (ia.ne.j.and.
     &                quench(j,j,n2nd,mxh,l2nd1,l2nd2) ) then
                    n2i=n2i+1
                    if (n2i.gt.mx2) then
                      write (*,'(a,i3,a)')  ' mklist> Molecule # '
     &                         ,nml,': too many atoms in list L2I'
                      stop
                    endif
                    l2i(n2i)=j
                  endif
                enddo
              endif
            enddo

! ____ If 'iow(iob)' exists and in 2ND list: into GENERAL list of 1-4 partners
            ioiob=iowat(iob)    !  existence of iow( iow(base) )
            if (ioiob.gt.0) then
              if( quench(ioiob,ioiob,n2nd,mxh,l2nd1,l2nd2) ) then
                n2i=n2i+1
                if (n2i.gt.mx2) then
                  write (*,'(a,i3,2a)')  ' mklist> Molecule # '
     &             ,nml,': too many atoms bound to ',nmat(iob)
                  stop
                endif
                l2i(n2i)=ioiob
              endif
            else
              ioiob=-10
            endif

          else
            iob=-10
            ioiob=-10
          endif

! ______ Atoms bound to 'ib' & in 2ND list(=are NOT in m.s of 'iv'):
!                    exclude from 2ND list & put in list 'l1i'
          n1i=0
          do i=1,nbdat(ib)
            ibd=ibdat(i,ib)
            if (iowat(ibd).eq.ib.and.
     &        quench(ibd,ibd,n2nd,mxh,l2nd1,l2nd2) ) then
              n1i=n1i+1
              if (n1i.gt.mxbd) then
                write (*,'(a,i3,2a)')  ' mklist> Molecule # ',nml,
     &             ': too many atoms bound to ',nmat(ib)
                stop
              endif
              l1i(n1i)=ibd
! _______ add atoms branching from 'l1i'-atoms to GENERAL list 1-4 partners
              do j=1,nbdat(ibd)
                jbd=ibdat(j,ibd)
                if (iowat(jbd).eq.ibd.and.
     &              quench(jbd,jbd,n2nd,mxh,l2nd1,l2nd2) ) then
                  n2i=n2i+1
                  if (n2i.gt.mx2) then
                    write (*,'(a,i3,2a)')  ' mklist> Molecule # ',nml,
     &              ': too many atoms bound to branches of ',nmat(ib)
                    stop
                  endif
                  l2i(n2i)=jbd
                endif
              enddo
            endif
          enddo
! _____________________________ check for further '1-4' partners
!                               belonging to moving set of base 'ib'
          im=ixmsat(ib)
          if (im.gt.0) then
            do i=latms1(im),latms2(im)
              if (quench(i,i,n2nd,mxh,l2nd1,l2nd2) ) then
                n2i=n2i+1
                if (n2i.gt.mx2) then
                  write (*,'(a,i3,a)')  ' mklist> Molecule # ',nml,
     &            ': too many atoms n list L2I '
                  stop
                endif
                l2i(n2i)=i
              endif
            enddo
          endif

          do is=i1s,i2s
            do i=latms1(is),latms2(is)  ! ============= atoms in m.s of 'iv'
! ________________________________________ Current 2ND list -> VdW-interact.
              if ((nvw+n2nd).gt.mxvw) then
                write (*,'(a,i4,a,i5)') ' mklist> Molecule # ',nml,
     &            ': Number of vdw-domains > ',mxvw
                stop
              endif
              ivwat1(i)=nvw+1      ! first and last vdW-domain ..
              ivwat2(i)=nvw+n2nd   ! .. per atom
              do j=1,n2nd
                nvw=nvw+1
                ixatvw(nvw)=i
                lvwat1(nvw)=l2nd1(j)
                lvwat2(nvw)=l2nd2(j)
              enddo  ! ... vdW-domains
! _________________________________________ General list of 1-4 partners
              if ((n14+n2i).gt.mx14) goto 1
              i14at1(i)=n14+1
              do j=1,n2i
                n14=n14+1
                ixat14(n14)=i
                l14at(n14)=l2i(j)
              enddo
! __________________________________ Special cases of 1-4 interactions
!                                    (list l1i, atoms iob,ib)
              iow=iowat(i)
              if (iow.ne.ib) then
                if ((n14+n1i).gt.mx14) goto 1
                do j=1,n1i  ! _____ branches of 'ib' NOT in m.s
                  n14=n14+1
                  ixat14(n14)=i
                  l14at(n14)=l1i(j)
                enddo
                if (ovlp.and.(it.eq.1.or.it.eq.2)) then ! _____ iob:
                  n14=n14+1                             ! b.lengths/angles
                  if (n14.gt.mx14) goto 1
                  ixat14(n14)=i
                  l14at(n14)=iob
                endif
                if (iowat(iow).ne.ib.and.it.eq.1) then ! ___ ib:
                  n14=n14+1                                 ! b.length
                  if (n14.gt.mx14) goto 1
                  ixat14(n14)=i
                  l14at(n14)=ib
                endif
              endif  ! ... spec. case
              i14at2(i)=n14
            enddo  ! ... atoms for moving set 'is'
          enddo  ! ... m.s for var. 'iv'

        endif  ! if there are moving sets 
      enddo  ! ... variables

      nvwml(nml)=nvw-ivwml1(nml)+1
      n14ml(nml)=n14-i14ml1(nml)+1
! _________________________________ some cleaning up
      do i=ifiat,ilaat
        if (ivwat1(i).le.0) then
          ivwat1(i)=1
          ivwat2(i)=0
        endif
        if (i14at1(i).le.0) then
          i14at1(i)=1
          i14at2(i)=0
        endif
      enddo

! ____________________________________________ Summary
!      do i=ifiat,ilaat
!        write (*,'(3a,i5,a)') ' ######## atom ',nmat(i),'(',i,')'
!        iv1=ivwat1(i)
!        iv2=ivwat2(i)
!        if (iv1.le.iv2) then
!          write(*,'(a)') ' ---> vdW :'
!          do j=iv1,iv2
!            write (*,'(i5,a,i5)') lvwat1(j),'-',lvwat2(j)
!          enddo
!        endif
!        i41=i14at1(i)
!        i42=i14at2(i)
!        if (i41.le.i42) then
!          write(*,'(a)') ' ---> 1-4 :'
!          write(*,'(10i5)') (l14at(j),j=i41,i42)
!        endif
!      enddo

      return

    1 write (*,'(a,i4,a,i5)') ' mklist> Molecule # ',nml,
     &                     ': Number of 1-4 interactions > ',mx14
      stop
      end
! *********************************************
      logical function quench(i1,i2,n,mx,l1,l2)

! ....................................................
! PURPOSE:  Correct size/number (n) of index ranges
!           given by lists 'l1' & 'l2' in order to
!           EXCLUDE overlaps with range 'i1-i2'
!
!           quench = true, if any overlap was obtained
!
! CALLS: none
!
! ....................................................

      implicit none
      integer mx, j, n, j1, l1, j2, l2, i1, i2, ja, k

      dimension l1(mx),l2(mx)

      quench=.false.  ! initialize

      j=1
      do while (j.le.n)  ! while there are sets
        j1=l1(j)
        j2=l2(j)

        if (i1.le.j2.and.i2.ge.j1) then  ! Overlap
          quench=.true.

          ja=0
          if (i1.gt.j1) then
            ja=1
            l2(j)=i1-1
          endif
          if (i2.lt.j2) then
            if (ja.gt.0) then  ! +1 set
              n=n+1
              if (n.gt.mx) then
                write (*,'(a)') ' quench> too many sets'
                stop
              endif
              do k=n,j+2,-1  ! shift sets
                l1(k)=l1(k-1)
                l2(k)=l2(k-1)
              enddo
              l2(j+1)=j2
            endif
            l1(j+ja)=i2+1
            ja=ja+1
          endif

          if (ja.eq.0) then  ! -1 set
            n=n-1
            do k=j,n
              l1(k)=l1(k+1)
              l2(k)=l2(k+1)
            enddo
          else
            j=j+ja
          endif

        else  ! No overlap
          j=j+1
        endif 

      enddo  ! ... sets

      return
      end

