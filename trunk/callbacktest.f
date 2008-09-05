      subroutine metropolis(eol,enw,dummy)
!f2py real*8 intent(in,out) eol
!f2py real*8 intent(in,out) enw
        external dummy
        delta =  dummy(enw) - dummy(eol)
        write (*,*) delta
      end