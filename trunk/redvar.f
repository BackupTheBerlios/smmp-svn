!**************************************************************
!
! This file contains the subroutines: redvar
!
! Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
!                      Shura Hayryan, Chin-Ku 
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
!
! **************************************************************


      subroutine redvar

! ...................................................................
!
! PURPOSE: Read global parameters for molecules from lines
!
!          +--------------------------------------------------+
!          |@ molecule no. : six floats separated by commas   |
!          +--------------------------------------------------+
!
!          NB: 1) if omit field with molecule no. assume: nml=1
!              2) last 3 float are angles in deg.
!
!          Read and interpret file to SET and FIX internal variables
!          by commands:
!
!          +-----------------------------------------+
!          |  molecule : residue : variable : value  |
!          +-----------------------------------------+
!
!        * Lines containing '&' assign FIXED variable(s), they will
!          not be varied during subsequent minimization etc.
!
!        * Empty LINES or lines containing '#' are ignored
!        * Several commands on same line must be separated by ';'
!        * Empty COMMANDS, i.e. ' : : ' are ignored
!        * All spaces are not significant and are therefore ignored
!
!        * A command consists of up to 4 (maxfld) fields, separated
!          by ':'
!
!        - last field     : value for VARIABLE (REAL)
!                           ! should never be empty
!        - 1st before last: name(s) (CHAR) or index(ices) of VARIABLE(S)
!        - 2nd before last: name(s) or index(ices) of RESIDUE(S)
!        - 3rd before last: name or number(ices) of MOLECULE(S)
!
!        * molecules, residues, variables can be identified, either by,
!          INDICES (zones 'n1-n2' possible) or NAMES
!
!        * several identifiers in a field can be separated by ','
!
!        * INDICES: for residues  - refer to numbering within molecule
!                 : for variables - refer to numbering within residue
!        * ZONES:   '-n2' indicates '1-n2'
!                   'n1-' indicates 'n1-(all)'
!        * NAMES or their ends can be indicated by wild-card '*'
!                are case-sensitive
!
!        Example:  phi:-65; psi:-45  >set all phi=-65, all psi=-45
!                  om*: 180 &  >set all omg, omt ... to 180 & fix them
!                  5 : x* : -60  >set all xi-angles of residue 5 to 60
!
! CALLS: setvar,extstr,iendst,ibegst,iopfil,iredin,iredrl
! ......................................................................


      include 'INCL.H'

!     functions
      integer iopfil, iendst, iredin, iredrl, ibegst
! maxfld: max. # of fields in one command
! maxide: max. # of identifiers in a field
! maxcmd: max. # of commands to be interpreted
! ilrg:   a large integer

      integer maxfld, maxide, maxcmd, ilrg, ifdend, icb, i, ib, ibz
      integer ife, ifb, ide, id, ice, ie, ieh, ieh1, iez, ihz, ifx, ifld
      integer ile, ii, ihy, ilb, iml, it, in ,io, inum, kbz, iv, ity, jb
      integer j, k, l, kez, kk, kv, lez, lbz,ll, ll1, ll2, n, nfi, ncmd
      integer nfld, nide, nml, ntlvr
      double precision vlvrx, rn, vr, val
      parameter (maxfld=4,
     &           maxide=30,
     &           maxcmd=5000,
     &           ilrg=1000000)

      character spcm,spfd,spcc,sphy,cmt,wdc,sfix,blnk, sglp,
     &          line*132,lincmd*132,linfld(maxfld)*132,linide*132,
     &          linh*132,strg(6)*17
      dimension ifdend(maxfld),vlvrx(mxvr),rn(6)
      logical fix,did,exa,forml(mxml),forrs(mxrs),forvr(mxvr),
     &        stvr(mxvr)
      data spcm/';'/,spfd/':'/,spcc/','/,sphy/'-'/,cmt/'#'/,wdc/'*'/,
     &     sfix/'&'/,blnk/' '/, sglp/'@'/


! ___________________________________ Checks
      ntlvr=ivrml1(ntlml)+nvrml(ntlml)-1
      if (ntlvr.eq.0) then
        write(logString,*)'redvar> No variables defined in molecule(s)'
        return
      endif
! ___________________________________ Initialize

      io=iopfil(lunvar,varfil,'old','formatted')
      if (io.eq.0) then
        write (logString, '(a,/,a,i3,2a)') 
     &    ' redvar> ERROR opening file to set variables:',
     &    ' LUN=',lunvar,' FILE=',varfil(1:iendst(varfil))
        stop
      elseif (io.eq.-1) then
        return
      endif
! ___________________________________ Initialization
      do i=1,ntlml
        forml(i)=.true.
        do j=irsml1(i),irsml2(i)
          forrs(j)=.true.
        enddo
      enddo

      do i=1,ntlvr
        fxvr(i)=.false.
        forvr(i)=.true.
        stvr(i)=.false.

        it=ityvr(i)  ! var. type
        if (it.eq.3) then      ! torsion
          vr=toat(iatvr(i))
        elseif (it.eq.2) then  ! b.angle
          vr=baat(iatvr(i))
        elseif (it.eq.1) then  ! b.length
          vr=blat(iatvr(i))
        else
          write(logString,*) 'redvar>  unknown variable type: ',it,' !'
          stop
        endif
        vlvrx(i)=vr
      enddo

      ncmd=0

    1 read (lunvar,'(a)',end=2) line
      ile=iendst(line)
! _________________________________ ! ignore empty and commentary lines
      if (ile.gt.0.and.index(line(1:ile),cmt).le.0) then

! _________________________________________ Global variables
        ilb = index(line(1:ile),sglp)+1
        if (ilb.ge.2) then

          if (index(line(ilb:ile),spfd).gt.0) then  ! field with mol.#

            call extstr(spfd,ilb,ile,line,lincmd,l)

            if (iredin(lincmd,nml).le.0.or.
     &          nml.le.0.or.nml.gt.ntlml) then
              write (logString, *) 
     &                 'redvar> ','Incorrect molecule number >',
     &                    lincmd(1:l),'<  Must be in range [1,',
     &                    ntlml,'] !'
              close(lunvar)
              stop
            endif

          else
            nml = 1  ! assume mol. #1
          endif

          l=ile-ilb+1
          if (l.le.0) goto 105
          lincmd=blnk
          lincmd(1:l)=line(ilb:ile)

          k = 1
          do i = 1,5  ! try to read 5 parameters
            call extstr(spcc,k,l,lincmd,linh,n)
            if (k.gt.l.or.iredrl(linh,rn(i)).le.0) goto 105
          enddo

          n=l-k+1  ! try 6th parameter
          if (n.le.0) goto 105
          linh=blnk
          linh(1:n)=lincmd(k:l)
          if (iredrl(linh,rn(6)).le.0) goto 105

! ---------------------------------------- check global angles
          if (  abs(rn(4)).gt.(1.8d2+1d-6)
     &     .or. abs(rn(5)).gt.(9d1+1d-6) 
     &     .or. abs(rn(6)).gt.(1.8d2+1d-6)
     &    ) goto 106
           
          do i = 1,3
            gbpr(i,nml) = rn(i)
          enddo
          do i = 4,6
            gbpr(i,nml) = rn(i)*cdr
          enddo

          goto 1

        endif  ! global vars


        ilb=1

        do while (ilb.le.ile)   ! ________________________ Commands
          call extstr(spcm,ilb,ile,line,lincmd,ice)

          if (ice.gt.0) then   ! ignore empty commands
            ncmd=ncmd+1
            if (ncmd.gt.maxcmd) goto 101

            ifx=index(lincmd(1:ice),sfix)
            if (ifx.gt.0) then      ! check for commands to fix variables
              fix=.true.
              lincmd(ifx:ifx)=blnk
              if (ifx.eq.ice) then
                ice=ice-1
                if (ice.eq.0) then    ! fix all
                  ice=1
                  lincmd(1:1)=wdc
                endif
              endif
            else
              fix=.false.
            endif

! _________________________________________ Extract Command Fields
            nfld=0
            icb=1
            do while (icb.le.ice)
              nfld=nfld+1
              if (nfld.gt.maxfld) goto 100
              call extstr(spfd,icb,ice,lincmd,linfld(nfld),ifdend(nfld))

              if (ifdend(nfld).le.0) then   ! empty field means 'all'
                linfld(nfld)(1:1)=wdc
                ifdend(nfld)=1
              endif

            enddo
! _______________________________ Interpret Command Fields (except last)
            do i=1,nfld-1
              ii=i
              ifld=nfld-i

              if (ifld.eq.3) then        !  Initialize Molecules
                do j=1,ntlml
                  forml(j)=.false.
                enddo
              elseif (ifld.eq.2) then    !  Initialize Residues
                do j=1,ntlml
                  do k=irsml1(j),irsml2(j)
                    forrs(k)=.false.
                  enddo
                enddo
              elseif (ifld.eq.1) then    !  Initialize Variables
                do j=1,ntlvr
                  forvr(j)=.false.
                enddo
              endif
! __________________________________ Identifiers in field
              nide=0
              ifb=1
              ife=ifdend(i)
              do while (ifb.le.ife)
                nide=nide+1
                if (nide.gt.maxide) goto 103
                call extstr(spcc,ifb,ife,linfld(ii),linide,ide)
                if (ide.le.0.or.linide(1:1).eq.wdc) then   ! ... All
                  if (ifld.eq.3) then        !  Mol.
                    do j=1,ntlml
                      forml(j)=.true.
                    enddo
                  elseif (ifld.eq.2) then    !  Res.
                    do j=1,ntlml
                      if (forml(j)) then
                        do k=irsml1(j),irsml2(j)
                          forrs(k)=.true.
                        enddo
                      endif
                    enddo
                  elseif (ifld.eq.1) then    !  Var.
                    do j=1,ntlml
                      if (forml(j)) then
                        do k=irsml1(j),irsml2(j)
                          if (forrs(k)) then
                            ll=ivrrs1(k)
                            do l=ll,ll+nvrrs(k)-1
                              forvr(l)=.true.
                            enddo
                          endif
                        enddo
                      endif
                    enddo
                  endif

                else  ! ...................... Identifier .ne. wdc

                  ihy=index(linide(1:ide),sphy)  ! ? zone of numbers

                  if (ihy.le.0) then        ! _____ No zone
                    if (iredin(linide,inum).gt.0) then  ! ... number
                      if (ifld.eq.3) then        !  Mol.

! ################### impossible # (inum) of molecule

                        if (inum.le.0.or.inum.gt.ntlml) then
                          write (logString, *) ' # 1: ',inum
                          goto 104
                        endif

                        forml(inum)=.true.
                      elseif (ifld.eq.2) then    !  Res.
                        do j=1,ntlml
                          if (forml(j)) then
                            nfi=irsml1(j)
                            k=inum+nfi-1

! ################### impossible # of residue (inum) in molecule

                            if (k.lt.nfi.or.k.gt.irsml2(j)) then
                              write (logString, *) ' # 2: ',inum
                              goto 104
                            endif 

                            forrs(k)=.true.
                          endif
                        enddo
                      elseif (ifld.eq.1) then    !  Var.
                        do j=1,ntlml
                          if (forml(j)) then
                            do k=irsml1(j),irsml2(j)
                              if (forrs(k)) then
                                nfi=ivrrs1(k)
                                l=inum+nfi-1

! ################### impossible # of variable (inum) in residue

                                if (l.lt.nfi.or.
     &                              l.gt.nfi+nvrrs(k)-1) then
                                  write (logString, *) ' # 3: ',inum
                                  goto 104
                                endif

                                forvr(l)=.true.
                              endif
                            enddo
                          endif
                        enddo
                      endif

                    else                                     ! ... Name
                      if (linide(ide:ide).eq.wdc) then
                        id=ide-1
                        exa=.false.
                      else                 ! exact match of names
                        id=ide
                        exa=.true.
                      endif

                      if (ifld.eq.3) then        !  Mol.
                        do j=1,ntlml
                          ib=ibegst(nmml(j))
                          if (ib.gt.0) then
                            linh=blnk
                            ieh=iendst(nmml(j))
                            ieh1=ieh-ib+1
                            linh(1:ieh1)=nmml(j)(ib:ieh)
                            if (((exa.and.ieh1.eq.id).or.
     &                         (.not.exa.and.ieh1.ge.id)).and.
     &                         linh(1:id).eq.linide(1:id))
     &                         forml(j)=.true.
                          endif
                        enddo        
                      elseif (ifld.eq.2) then    !  Res.
                        do j=1,ntlml
                          if (forml(j)) then
                            do k=irsml1(j),irsml2(j)
                              ib=ibegst(seq(k))
                              if (ib.gt.0) then
                                linh=blnk
                                ieh=iendst(seq(k))
                                ieh1=ieh-ib+1
                                linh(1:ieh1)=seq(k)(ib:ieh)
                                if (((exa.and.ieh1.eq.id).or.
     &                             (.not.exa.and.ieh1.ge.id))
     &                          .and.linh(1:id).eq.linide(1:id))
     &                            forrs(k)=.true.
                              endif
                            enddo
                          endif
                        enddo
                      elseif (ifld.eq.1) then    !  Var.
                        do j=1,ntlml
                          if (forml(j)) then
                            do k=irsml1(j),irsml2(j)
                              if (forrs(k)) then
                                ll=ivrrs1(k)
                                do l=ll,ll+nvrrs(k)-1
                                  ib=ibegst(nmvr(l))
                                  if (ib.gt.0) then
                                    linh=blnk
                                    ieh=iendst(nmvr(l))
                                    ieh1=ieh-ib+1
                                    linh(1:ieh1)=nmvr(l)(ib:ieh)
                                    if (((exa.and.ieh1.eq.id)
     &                         .or.(.not.exa.and.ieh1.ge.id))
     &                          .and.linh(1:id).eq.linide(1:id))
     &                              forvr(l)=.true.
                                  endif
                                enddo
                              endif
                            enddo
                          endif
                        enddo
                      endif

                    endif

                  else                                       ! ___ Zone

! ################### impossible zone '-' (without integer)

                    if (ide.eq.1.and.ihy.eq.ide) then
                      write (logString, *) ' # 4: ',ide
                      goto 104
                    endif

                    if (ihy.eq.1) then
                      ibz=1
                    else
                      linh=blnk
                      linh=linide(1:ihy-1)

! ################### impossible (to read) integer before '-'

                      if (iredin(linh,ibz).le.0.or.ibz.le.0)
     &                  then
                        write (logString, *) ' # 5 '
                        goto 104
                      endif

                    endif
                    if (ihy.eq.ide) then
                      iez=ilrg
                    else
                      linh=blnk
                      linh=linide(ihy+1:ide)

! ################### impossible (to read) integer after '-'

                      if (iredin(linh,iez).le.0.or.iez.le.0.or.
     &                  iez.lt.ibz) then
                        write (logString, *) ' # 6 '
                        goto 104
                      endif

                    endif

                    if (ifld.eq.3) then        !  Mol.
                      if (iez.gt.ntlml) iez=ntlml
                      do j=ibz,iez
                        forml(j)=.true.
                      enddo
                    elseif (ifld.eq.2) then    !  Res.
                      do j=1,ntlml
                        if (forml(j)) then
                          kbz=irsml1(j)+ibz-1
                          kez=irsml1(j)+iez-1
                          if (kez.gt.irsml2(j)) then
                            kk=irsml2(j)
                          else
                            kk=kez
                          endif
                          do k=kbz,kk
                            forrs(k)=.true.
                          enddo
                        endif
                      enddo
                    elseif (ifld.eq.1) then    !  Var.
                      do j=1,ntlml
                        if (forml(j)) then
                          do k=irsml1(j),irsml2(j)
                            kv=nvrrs(k)
                            if (forrs(k).and.kv.gt.0) then
                              ll=ivrrs1(k)
                              lbz=ll+ibz-1
                              if (iez.gt.kv) then
                                lez=ll+kv-1
                              else
                                lez=ll+iez-1
                              endif
                              do l=lbz,lez
                                forvr(l)=.true.
                              enddo
                            endif
                          enddo
                        endif
                      enddo
                    endif

                  endif
                endif

              enddo  ! ... identifiers
            enddo  ! ... Fields (excl. value)

! _____________________________________________________ Execute Command

            if (iredrl(linfld(nfld),val).gt.izero) then  ! Read Value
              did=.false.
              do i=1,ntlvr
                if (forvr(i)) then
                  did=.true.
                  vlvrx(i)=val

                  fxvr(i)=fix

                  stvr(i)=.true.
                endif
              enddo
              if (.not.did) write (logString, '(3a)')
     &          ' redvar> No variables affected by command >',
     &          lincmd(1:ice),'<'
            else

              ll1=ibegst(linfld(nfld))
              ll2=iendst(linfld(nfld))
              write (logString, *) 'll1,ll2, linfld(nfld): ',ll1,ll2,
     &           '>',linfld(nfld)(ll1:ll2),'<'

              goto 102
            endif

          endif
        enddo  ! ... Commands at one line
      endif
      goto 1

    2 close(lunvar)
! __________________________ Summary
      iv=0
      do i=1,ntlml

        ie=iendst(nmml(i))

        do j =1,6
          if (gbpr(j,i).ne.zero) then

             do k = 1,3
                write(strg(k),'(f17.6)') gbpr(k,i)
             enddo
             do k = 4,6
                write(strg(k),'(f17.6)') (gbpr(k,i)*crd)
             enddo

             write (logString, '(3a,/,1x,5(a,2x),a)') ' redvar> ',
     &                               nmml(i)(1:ie),
     &                                    ' with global parameters:',
     &                              (strg(k)(ibegst(strg(k)):),k=1,6)
             call setvar(i,vlvrx)
             goto 3
          endif
        enddo

    3   if (nvrml(i).gt.0) then
          iml=i
          did=.false.
          in=0
          jb=irsml1(i)-1
          do j=irsml1(i),irsml2(i)
            kk=ivrrs1(j)
            do k=kk,kk+nvrrs(j)-1
              iv=iv+1
              if (stvr(iv)) then
                did=.true.
                if (fxvr(iv)) then
                  write(logString,'(3a,i4,1x,4a,f10.3,a)') ' redvar> ',
     &                nmml(i)(1:ie),': residue ',j-jb,seq(j),
     &                ': ',nmvr(iv),' set ',vlvrx(iv),'   Fixed'
                else
                  write (logString, '(3a,i4,1x,4a,f10.3)') ' redvar> ',
     &                nmml(i)(1:ie),': residue ',j-jb,seq(j),
     &                ': ',nmvr(iv),' set ',vlvrx(iv)
                endif
                ity=ityvr(iv)
                if (ity.eq.3.or.ity.eq.2)
     &            vlvrx(iv)=vlvrx(iv)*cdr             ! angles
                  
              else
                in=in+1
              endif
            enddo
          enddo
          if (did) then
            if (in.gt.0) write (logString, '(3a,i5,a)')
     &        ' redvar> Molecule ',nmml(i)(1:ie),': ',in,
     &        ' variable(s) remain unchanged'
            call setvar(iml,vlvrx)
          else
            write (logString, '(3a)') ' redvar> Molecule ',
     &        nmml(i)(1:ie),': No internal variables changed'
          endif
        endif
      enddo

      return
! ____________________________________________________________ Errors
  100 write (logString, '(3a)') ' redvar> Cannot interpret command >',
     &                 lincmd(1:ice),'<'
      close(lunvar)
      stop
  101 write (logString, '(a,i5,a)') ' redvar> Command number ',ncmd,
     &       ' reached'
      close(lunvar)
      stop
  102 write (logString, '(3a)') ' redvar> Cannot read value from >',
     &                   lincmd(1:ice),'<'
      close(lunvar)
      stop
  103 write (logString, '(a,i3,3a)') ' redvar> Cannot read >',maxide,
     &         ' identifiers from >',linfld(ii)(1:ife),'<'
      close(lunvar)
      stop
  104 write (logString, '(5a)') ' redvar> Error in identifier >',
     &            linide(1:ide),'< of command >',lincmd(1:ice),'<'
      close(lunvar)
      stop
  105 write (logString, '(a,/,a,/,2a,/)') 
     &                           ' redvar> line with global paramters:',
     &                           line(1:ile),' must contain 6 floating',
     &                           ' point numbers separated by commas !'
      close(lunvar)
      stop

  106 write (logString, '(a,/,a,/,2a,/)') 
     &                           ' redvar> line with global paramters:',
     &                           line(1:ile),' angles must be inside ',
     &'ranges [-180,180], [-90,90], and [-180,180] Deg., respectively !'
      close(lunvar)
      stop

      end
