!**************************************************************
!
! This file contains the subroutines: extstr,ibegst,iendst,
!                                     iredin,iredrl,iopfil,
!                                     tolost,toupst
!
! Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
!                      Shura Hayryan, Chin-Ku 
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
!
! **************************************************************


      subroutine extstr(spr,ib,ie,str,strn,l)

! ..........................................................
! PURPOSE:  Extract substring preceeding separator 'spr'
!           from 'str' searching from position 'ib' up to
!           position 'ie' and put it into 'strn(1:l)'.
!           'ib' is shifted to position following 'spr' or
!           to 'ie+1', if 'spr' is not found
!
!          ! 'spr' should not be blank
!
! CALLS: ibegst,iendst
! ..........................................................
      integer ib, ie, l
      integer ln, is, ibegst, j, iendst

      character spr,blnk,str*(*),strn*(*)
      character*255 logString
      
      data blnk/' '/
      
      integer i, ic, ish, ii

      if (spr.eq.blnk) then
        write (logString, *) ' extstr> Separator should not be blank'
        stop
      endif

      l=0
      ln=len(strn)
      strn=blnk
      is=index(str(ib:ie),spr)  ! position of spr

      if (is.lt.1) then          ! _________ no separator

        l=ie-ib+1
        if (ln.lt.l) goto 1
        strn(1:l)=str(ib:ie)
        ib=ie+1
      elseif (is.eq.1) then      ! _________ empty substring
        ib=ib+1
        return
      else                       ! _________ found separator
        l=is-1
        if (ln.lt.l) goto 1
        strn(1:l)=str(ib:ib+l-1)
        ib=ib+is
      endif

      i=ibegst(strn)

      if (i.lt.1) then   ! empty substring

        l=0
        strn=blnk
! ____________________________ make string in 'strn' left justified
      elseif (i.gt.1) then
        j=iendst(strn)
        l=j-i+1
        strn(1:l)=strn(i:j)
        strn(l+1:ln)=blnk
      else
        l=iendst(strn)
      endif

      return
! ______________________________________________________________ Error
    1 write (logString, '(a)') 
     &   ' extstr> Substring to be extracted is too long !'
      stop

      end
! **********************************
      integer*4 function ibegst(str)

! .............................................................
! PURPOSE: returns position of 1st non-blank character in 'str'
!
! CALLS: none
!
! .............................................................

      implicit none
      integer i
      character blnk,str*(*)
      data blnk/' '/

      do i=1,len(str)
        if (str(i:i).ne.blnk) then
          ibegst=i
          return
        endif
      enddo

      ibegst=0

      return
      end
! **********************************
      integer*4 function iendst(str)

! ..............................................................
! PURPOSE: returns position of last non-blank character in 'str'
!
! CALLS: none
!
! ..............................................................

      implicit none
      integer i
      character blnk,str*(*)
      data blnk/' '/

      do i=len(str),1,-1
        if (str(i:i).ne.blnk) then
          iendst=i
          return
        endif
      enddo

      iendst=0

      return
      end
! **************************************
      integer*4 function iredin(line,in)

! ..........................................
! PURPOSE: Read integer*4 value 'in' from 'line'
!          with format 'i9'
!
!          iredin=0 : error status
!          iredin=1 : success
!
! CALLS: ibegst,iendst
! ..........................................

      implicit none
      integer mxd, ib, ibegst, ie, iendst, il, i0, i9, i, ii, in
      parameter (mxd=9)                 ! max. # of digits

      character blnk,value*(mxd),line*(*)
      data blnk/' '/

      iredin=0
      ib=ibegst(line)
      if (ib.gt.0) then
        ie=iendst(line)
        il=ie-ib
        if (il.lt.mxd) then
          i0=ichar('0')
          i9=ichar('9')
          do i=ib,ie
            ii=ichar(line(i:i))
            if (ii.lt.i0.or.ii.gt.i9) goto 1
          enddo
          value=blnk
          value(mxd-il:mxd)=line(ib:ie)
          read(value,'(i9)',err=1) in
          iredin=1
        endif
      endif
    1 return
      end
! *************************************
      integer*4 function iredrl(line,r)

! ..........................................
! PURPOSE: Read real*8 value 'r' from 'line'
!          with format 'd17.6'
!
!          iredrl=0 : error status
!          iredrl=1 : success
!
! CALLS: ibegst,iendst
! ..........................................

      implicit none
      integer mxd, mxap, mxip, ib, ibegst, ie, iendst, il, ip, ibp
      parameter (mxd =17,   ! max. # of digits
     &           mxap= 6,   ! max. # of digits after period
     &           mxip=mxd-mxap)
 
      real*8 r
      character per,blnk,value*(mxd),line*(*)
      data per/'.'/,blnk/' '/

      iredrl=0

      ib=ibegst(line)
      if (ib.gt.0) then
        ie=iendst(line)
        if (index(line(ib:ie),',').gt.0) return
        il=ie-ib+1
        ip=index(line,per)
        value=blnk
        if (ip.gt.0) then       !  found period
          ibp=ip-ib
          if (il.le.mxd.and.ibp.lt.mxip.and.ie-ip.le.mxap) then
            value(mxip-ibp:)=line(ib:ie)
            read (value,'(d17.6)',err=1) r
            iredrl=1
          endif
        else                    !  no period
          if (il.lt.mxip) then
            value(mxip-il:)=line(ib:ie)//per
            read (value,'(d17.6)',err=1) r
            iredrl=1
          endif
        endif
      endif

    1 return
      end
! **************************
      subroutine tolost(str)

! ..........................................
!  PURPOSE:  converts 'string' to lower-case
!  INPUT:    str - string to be converted
!  CALLS:    ibegst,iendst
! ..........................................

      include 'INCL.H'

      character*(*) str
      
      integer ibegst, iendst
      
      integer i, ic, ii, ish
      ii=ibegst(str)
      if (ii.gt.0) then
        ish=idupa-idloa
        do i=ii,iendst(str)
          ic=ichar(str(i:i))
          if (ic.ge.idupa.and.ic.le.idupz) str(i:i)=char(ic-ish) 
        enddo
      endif

      return
      end
! **************************
      subroutine toupst(str)

! ..........................................
!  PURPOSE:  converts 'string' to upper-case
!  INPUT:    str - string to be converted
!  CALLS:    ibegst,iendst
! ..........................................

      include 'INCL.H'

      character str*(*)
      
      integer iendst, ibegst
      
      integer i, ii, ic, ish

      ii=ibegst(str)
      if (ii.gt.0) then
        ish=idupa-idloa
        do i=ii,iendst(str)
          ic=ichar(str(i:i))
          if (ic.ge.idloa.and.ic.le.idloz) str(i:i)=char(ic+ish) 
        enddo
      endif

      return
      end
! *****************************************************
      integer*4 function iopfil(lun,filnam,stat,format)

! ........................................................
! PURPOSE: open 'lun' with 'filnam' 'stat' 'format'
!
!          returns: 1 = file successful opened
!                   0 = error during open of existing file
!                  -1 = file does not exist
!
! CALLS: ibegst
! ........................................................

      integer lun
      integer i, ibegst, j, k
      logical exs
      character*(*) filnam,stat,format
              
      iopfil=0
      
      if (lun.gt.0.and.lun.lt.100) then
        i=ibegst(filnam)
        if (i.gt.0) then
          inquire(file=filnam(i:),exist=exs)
          if (exs) then
            j=ibegst(stat)
            k=ibegst(format)
            if (j.gt.0.and.k.gt.0) then
              open(lun,file=filnam(i:),status=stat(j:),
     &             form=format(k:),err=1)
              iopfil=1
            endif
          else
            iopfil=-1
          endif
        endif
      endif

    1 return
      end

