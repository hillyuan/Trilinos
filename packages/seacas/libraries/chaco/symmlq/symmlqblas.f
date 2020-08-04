C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     symmlqblas  fortran
*
*     daxpy    dcopy    ddot     dnrm2
*
*** from netlib, Thu May 16 21:00:13 EDT 1991 ***
*** Declarations of the form dx(1) changed to dx(*)
*
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine daxpy(n,da,dx,incx,dy,incy)

c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.

      double precision dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n

      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20

c        code for unequal increments or equal increments
c          not equal to 1

      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return

c        code for both increments equal to 1

c        clean-up loop

   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end

*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine  dcopy(n,dx,incx,dy,incy)

c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.

      double precision dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n

      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20

c        code for unequal increments or equal increments
c          not equal to 1

      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return

c        code for both increments equal to 1

c        clean-up loop

   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return
      end

*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      double precision function ddot(n,dx,incx,dy,incy)

c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.

      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n

      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20

c        code for unequal increments or equal increments
c          not equal to 1

      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return

c        code for both increments equal to 1

c        clean-up loop

   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end

*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      double precision function dnrm2 ( n, dx, incx)
      integer          next
      double precision   dx(*), cutlo, cuthi, hitest, sum, xmax,zero,one
      data   zero, one /0.0d0, 1.0d0/

c     euclidean norm of the n-vector stored in dx() with storage
c     increment incx .
c     if    n .le. 0 return with result = 0.
c     if n .ge. 1 then incx must be .ge. 1

c           c.l.lawson, 1978 jan 08

c     four phase method     using two built-in constants that are
c     hopefully applicable to all machines.
c         cutlo = maximum of  dsqrt(u/eps)  over all known machines.
c         cuthi = minimum of  dsqrt(v)      over all known machines.
c     where
c         eps = smallest no. such that eps + 1. .gt. 1.
c         u   = smallest positive no.   (underflow limit)
c         v   = largest  no.            (overflow  limit)

c     brief outline of algorithm..

c     phase 1    scans zero components.
c     move to phase 2 when a component is nonzero and .le. cutlo
c     move to phase 3 when a component is .gt. cutlo
c     move to phase 4 when a component is .ge. cuthi/m
c     where m = n for x() real and m = 2*n for complex.

c     values for cutlo and cuthi..
c     from the environmental parameters listed in the imsl converter
c     document the limiting values are as follows..
c     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
c                   univac and dec at 2**(-103)
c                   thus cutlo = 2**(-51) = 4.44089e-16
c     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
c                   thus cuthi = 2**(63.5) = 1.30438e19
c     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
c                   thus cutlo = 2**(-33.5) = 8.23181d-11
c     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19
c     data cutlo, cuthi / 8.232d-11,  1.304d19 /
c     data cutlo, cuthi / 4.441e-16,  1.304e19 /
      data cutlo, cuthi / 8.232d-11,  1.304d19 /

      if(n .gt. 0) go to 10
         dnrm2  = zero
         go to 300

   10 assign 30 to next
      sum = zero
      nn = n * incx
c                                                 begin main loop
      i = 1
   20    go to next,(30, 50, 70, 110)
   30 if( dabs(dx(i)) .gt. cutlo) go to 85
      assign 50 to next
      xmax = zero

c                        phase 1.  sum is zero

   50 if( dx(i) .eq. zero) go to 200
      if( dabs(dx(i)) .gt. cutlo) go to 85

c                                prepare for phase 2.
      assign 70 to next
      go to 105

c                                prepare for phase 4.

  100 i = j
      assign 110 to next
      sum = (sum / dx(i)) / dx(i)
  105 xmax = dabs(dx(i))
      go to 115

c                   phase 2.  sum is small.
c                             scale to avoid destructive underflow.

   70 if( dabs(dx(i)) .gt. cutlo ) go to 75

c                     common code for phases 2 and 4.
c                     in phase 4 sum is large.  scale to avoid overflow.

  110 if( dabs(dx(i)) .le. xmax ) go to 115
         sum = one + sum * (xmax / dx(i))**2
         xmax = dabs(dx(i))
         go to 200

  115 sum = sum + (dx(i)/xmax)**2
      go to 200

c                  prepare for phase 3.

   75 sum = (sum * xmax) * xmax

c     for real or d.p. set hitest = cuthi/n
c     for complex      set hitest = cuthi/(2*n)

   85 hitest = cuthi/float( n )

c                   phase 3.  sum is mid-range.  no scaling.

      do 95 j =i,nn,incx
      if(dabs(dx(j)) .ge. hitest) go to 100
   95    sum = sum + dx(j)**2
      dnrm2 = dsqrt( sum )
      go to 300

  200 continue
      i = i + incx
      if ( i .le. nn ) go to 20

c              end of main loop.

c              compute square root and adjust for scaling.

      dnrm2 = xmax * dsqrt(sum)
  300 continue
      return
      end
