! **************************************************************
!
! This file contains the subroutines: mincjg
!
! Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
!                      Shura Hayryan, Chin-Ku
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
!
! ********************************************************************
      subroutine mincjg(n,mxn,x,f,g,acur,d,xa,ga,dt,yt,gt,maxfun,nfun)

! ....................................................................
!
!  Conjugate Gradient Minimizer
!
!  INPUT:   X,F,G - variables, value of FUNC, gradient at START/
!           ACUR - convergence is assumed if ACUR > SUM ( G(I)**2 )
!           MAXFUN - maximum overall number of function calls
!
!  OUTPUT:  X,F,G - variables, value of FUNC, gradient at MINIMUM
!           NFUN  - overall number of function calls used
!
!  ARRAYS:  D,XA,GA,YT,DT,GT - dimension N
!
!  CALLS:   MOVE - calculate function & its gradients for current X
!
!  PARAMETERS:  AMF    - rough estimate of first reduction in F, used
!                        to guess initial step of 1st line search
!               MXFCON - see 'ier=4'
!               MAXLIN -
!
!  DIAGNOSTICS (ier)
!
!           = 0: minimization completed successfully
!           = 1: number of steps reached MAXFUN
!           = 2: line search was abandoned
!           = 3: search direction is uphill
!           = 4: two consecutive line searches failed to reduce F
! ....................................................................

      double precision amf, tol, eps, gsqrd, fmin, f, t, g, ga, d, xa, x
      double precision gsq2, gnew, dfpr, stmin, fuit, gdit, gt, gdmi
      double precision sbound, stepch, step, wo, sum, fch, acur, ddspln
      double precision gspln, beta, gama, yt, dt, gamden, di

      integer mxfcon, maxlin, mxn, ier, iter, iterfm, iterrs, nfun
      integer nfopt, i, n, nfbeg, iretry, maxfun

      parameter (AMF = 10.d0,
     &           MXFCON = 2,
     &           MAXLIN = 5,
     &           TOL = 1.d-7,  ! controls 'stepch'
     &           EPS = .7d0)

      dimension x(mxn),g(mxn),
     &          d(mxn),xa(mxn),ga(mxn),dt(mxn),yt(mxn),gt(mxn)

      character(255) logString
      
      ier = 0
      iter = 0
      iterfm = 0
      iterrs = 0
      nfun = 0
      nfopt = nfun

      gsqrd = 0.d0
      fmin = f

      do i=1,n
        t = g(i)
        ga(i) = t
        d(i) = -t
        gsqrd = gsqrd + t**2
        xa(i) = x(i)
      enddo

      gsq2 = .2d0 * gsqrd
      gnew = -gsqrd
      dfpr = AMF
      stmin = AMF / gsqrd

    1 iter = iter + 1

      fuit = f
      gdit = 0.d0

      do i=1,n
        t = g(i)
        gt(i) = t
        gdit = gdit + d(i) * t
      enddo

      if ( gdit .ge. 0.d0 ) then
        write (logString, *) ' mincjg>  search direction is uphill'
        ier = 3
        goto 6
      endif

      gdmi = gdit
      sbound = -1.d0
      nfbeg = nfun
      iretry = -1

      stepch = min( stmin, abs ( (dfpr/gdit) ) )
      stmin = 0.d0

    2 step = stmin + stepch
      wo = 0.d0

      do i=1,n
        t = stepch * d(i)
        wo = max( wo, abs( t ) )
        x(i) = xa(i) + t
      enddo

      if ( wo .gt. TOL )  then

        nfun = nfun + 1
        call move(nfun,n,f,x,g)

        gnew = 0.d0
        sum = 0.d0

        do i=1,n
          t = g(i)
          gnew = gnew + d(i) * t
          sum = sum + t**2
        enddo

        fch = f - fmin

        if ( fch .le. 0.d0 ) then

          if ( fch .lt. 0.d0 .or. (gnew/gdmi) .ge. -1.d0 ) then

            fmin = f
            gsqrd = sum
            gsq2 = .2d0 * gsqrd
            nfopt = nfun

            do i=1,n
              xa(i) = x(i)
              ga(i) = g(i)
            enddo

          endif

          if ( sum .le. ACUR ) return  ! normal end

        endif

        if ( nfun .eq. MAXFUN ) then
          ier = 1
          return
        endif

      else  !   stepch is effectively zero

        if ( nfun .gt. (nfbeg + 1) .or.
     &       abs(gdmi/gdit) .gt. EPS ) then

          ier=2
          write (logString, *) 
     &       ' mincjg>  too small step in search direction'
        endif

        goto 6

      endif

      wo = (fch + fch) / stepch - gnew - gdmi
      ddspln = (gnew - gdmi) / stepch

      if ( nfun .le. nfopt ) then

        if ( gdmi * gnew .le. 0.d0 )  sbound = stmin

        stmin = step
        gdmi = gnew
        stepch = -stepch

      else

        sbound = step

      endif

      if ( fch .ne. 0.d0 )  ddspln = ddspln + (wo + wo) / stepch

      if ( gdmi .eq. 0.d0 )  goto 6

      if ( nfun .le. (nfbeg + 1) )  goto 4

      if ( abs(gdmi/gdit) .le. EPS )  goto 6

    3 if ( nfun .ge. (nfopt + MAXLIN) ) then

        ier = 2
        goto 6

      endif

    4 if ( sbound .lt. -.5d0 ) then
        stepch = 9.d0 * stmin
      else
        stepch = .5d0 * ( sbound - stmin )
      endif

      gspln = gdmi + stepch * ddspln

      if ( (gdmi * gspln) .lt. 0.d0 )  stepch = stepch * gdmi /
     &                                          (gdmi - gspln)

      goto 2

    5 sum = 0.d0
      do i=1,n
        sum = sum + g(i) * gt(i)
      enddo

      beta = (gsqrd - sum) / (gdmi - gdit)

      if ( abs(beta * gdmi) .gt. gsq2) then
        iretry = iretry + 1
        if (iretry .le. 0) goto 3
      endif

      if ( f .lt. fuit )  iterfm = iter

      if ( iter .ge. (iterfm + MXFCON) ) then
        ier = 4
        write (logString, *) 
     &    ' mincjg>  line search failed to reduce function'
        return
      endif

      dfpr = stmin * gdit

      if ( iretry .gt. 0 ) then

        do i=1,n
          d(i) = -g(i)
        enddo

        iterrs = 0

      else

        if (iterrs .ne. 0 .and. (iter - iterrs) .lt. (n-1) .and.
     &      abs(sum) .lt. gsq2 ) then

          gama = 0.d0
          sum = 0.d0

          do i=1,n
            t = g(i)
            gama = gama + t * yt(i)
            sum = sum + t * dt(i)
          enddo

          gama = gama / gamden

          if ( abs( (beta * gdmi + gama * sum) ) .lt. gsq2) then

            do i=1,n
              d(i) = -g(i) + beta * d(i) + gama * dt(i)
            enddo

            goto 1

          endif

        endif

        gamden = gdmi - gdit

        do i=1,n
          t = g(i)
          di = d(i)
          dt(i) = di
          yt(i) = t - gt(i)
          d(i) = -t + beta * di
        enddo

        iterrs = iter

      endif

      goto 1

    6 if ( nfun .ne. nfopt ) then

        f = fmin

        do i=1,n
          x(i) = xa(i)
          g(i) = ga(i)
        enddo

      endif

      if ( ier .eq. 0 )  goto 5

      nfun = nfun + 1
      call move(nfun,n,f,x,g)

      write (logString, *) ' mincjg> ier = ',ier

      return
      end

