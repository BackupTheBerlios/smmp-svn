!**************************************************************
!
! This file contains the subroutines: redseq
!
! Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
!                      Shura Hayryan, Chin-Ku 
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
!
! **************************************************************


      subroutine redseq

! ............................................................
! PURPOSE: read 'lunseq' 'seqfil', extract names of molecules,
!          sequences
!
! Molecules are separated by lines containing char. '#',
!           a name for the molecule may follow '#' on this line
! Residue   names can be of 1-4 characters to be separated by ' '
!
! Returns: ntlml,nmml,irsml1,irsml2,seq
!
! CALLS:   ibegst,iendst,iopfil,tolost
! ............................................................

      include 'INCL.H'

      ! Function definitions
      integer iofil, iendst, ibegst, iopfil

      character blnk,res*4,line*132,hlin*132
      integer i, i1, i2, ifirs, ie, id, ib, ic, ii, l, lg, nln, nrs
      data blnk/' '/

      if (iopfil(lunseq,seqfil,'old','formatted').le.izero) then
        write (logString, '(a,/,a,i3,2a)') 
     &    ' redseq> ERROR opening sequence file:',
     &      ' LUN=',lunseq,' FILE=',seqfil(1:iendst(seqfil))
        stop
      endif

!      ntlml=0
      if (ntlml.gt.0) then
        nrs = irsml2(ntlml)
      else
        nrs = 0
      endif
      nln=0
    1 line=blnk
      nln=nln+1
      read (lunseq,'(a)',err=4,end=3) line
      lg=iendst(line)

      if (lg.gt.0) then  ! line not empty

        ic = index(line(1:lg),'#')

        if (ic.gt.0) then  ! found '#'

! ____________________________________ new molecule

          if (ntlml.gt.0) then   ! check previous molecule

            irsml2(ntlml) = nrs

            if ((nrs-irsml1(ntlml)+1).eq.0) then
              write (logString, '(2a)') ' redseq> IGNORE molecule: ',
     &                        nmml(ntlml)(1:iendst(nmml(ntlml)))
              ntlml=ntlml-1
            endif
          endif
          ntlml=ntlml+1
          if (ntlml.gt.mxml) then
            write (logString, '(a,i4,2a)') 
     &                          ' redseq> NUMBER of molecules > '
     &                          ,mxml,' in ',seqfil(1:iendst(seqfil))
            close(lunseq)
            stop
          endif
          nmml(ntlml)=blnk
          irsml1(ntlml)=nrs+1
          ic=ic+1

          if (ic.le.lg) then
! ___________________________________ extract name of molecule

            hlin=blnk
            hlin(1:)=line(ic:lg)
            i1=ibegst(hlin)
            if (i1.gt.0) then
              i2=iendst(hlin)
              l=i2-i1+1
              if (l.gt.80) i2=i2-l+80
              nmml(ntlml)(1:)=hlin(i1:i2)
            else
              write(nmml(ntlml)(1:4),'(i4)') ntlml
            endif
          else
            write(nmml(ntlml)(1:4),'(i4)') ntlml
          endif

        else  ! no '#'

! _________________________________________ sequence

          ib=ibegst(line)
          if (ib.eq.0) goto 1   ! empty line

          if (ntlml.eq.0) then
            ntlml=1
            nmml(1)(1:)='   1'
          endif

          ie=iendst(line)
! ___________________________________ extract names of residues
    2     id=index(line(ib:ie),blnk)-1   ! find next separator
          if (id.gt.0) then
            ii=ib+id-1
          else
            ii=ie
            id=ie-ib+1
          endif

          if (id.gt.4) then
            write (logString, '(4a)') ' redseq> INVALID residue NAME >',
     &                       line(ib:ii),'< in ',
     &      seqfil(1:iendst(seqfil))
            close(lunseq)
            stop
          else

            nrs=nrs+1
            if (nrs.gt.mxrs) then
              write (logString, '(a,i4,2a)') 
     &                       ' redseq> NUMBER of residues > '
     &                       ,mxrs,' in ',seqfil(1:iendst(seqfil))
              close(lunseq)
              stop
            endif

            seq(nrs)=blnk
            seq(nrs)(1:)=line(ib:ii)

          endif

          ii=ii+1
          do i=ii,ie
            if (line(i:i).ne.blnk) then
              ib=i
              goto 2
            endif
          enddo

        endif

      endif

      goto 1

    3 close(lunseq)
! ___________________________________ output

      if (nrs.eq.0) then
        write (logString, '(2a)') ' redseq> no residues found in ',
     &                   seqfil(1:iendst(seqfil))
        stop
      else

        do i=1,ntlml

          ifirs=irsml1(i)

          if (.not.flex) then
            res=seq(ifirs)
            call tolost(res)
            if (res(1:3).eq.'pro') seq(irsml1(i))='pron'  ! only ECEPP/3
          endif

          if (i.eq.ntlml) then   ! Check last molecule
            if ((nrs-ifirs+1).eq.0) then
              write (logString, '(2a)') ' redseq> IGNORE molecule '
     &                        ,nmml(ntlml)(1:iendst(nmml(ntlml)))
              ntlml=ntlml-1
              if (ntlml.eq.0) then
                write (logString, '(2a)') 
     &             ' redseq> no residues found in ',
     &             seqfil(1:iendst(seqfil))
                stop
              endif
              return
            endif
            irsml2(i)=nrs
          endif

!          write (logString, '(/,a,i4,2a)') ' redseq> ',irsml2(i)-irsml1(i)+1,
!     &           ' residue(s) in molecule: ',
!     &           nmml(i)(1:iendst(nmml(i)))
!          write (logString, '(15(1x,a))') (seq(j),j=irsml1(i),irsml2(i))

        enddo

      endif
      return
! _______________________________________________ error

    4 write (logString, '(a,i4,2a)') ' redseq> ERROR reading line No. ',nln,
     & ' in ',seqfil(1:iendst(seqfil))
      close(lunseq)
      stop

      end

