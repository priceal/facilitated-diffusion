!---------------------------------------------------------------!
! subtroutines for random walk processes in 1D
!---------------------------------------------------------------!

!---------------------------------------------------------------!
! Block data for the fpd routines. These values are coefficients
! and exponent factors for the expansions of the fpds.
!---------------------------------------------------------------!
block data fpd

 real(8) :: coefficient1d(1:10), exponent1d(1:10),coefficient2d(1:30), exponent2d(1:30)
 
 common /block_fpd/ coefficient1d, exponent1d, coefficient2d, exponent2d

 data coefficient1d/6.283185307179586, -18.84955592153876, 31.41592653589793, -43.982297150257104, &
  56.548667764616276, -69.11503837897544, 81.68140899333463, -94.24777960769379, &
  106.81415022205297, -119.38052083641213/
  
 data exponent1d/-4.934802200544679, -44.41321980490211, -123.37005501361698, -241.80530782668927, &
  -399.718978244119, -597.1110662659062, -833.9815718920507, -1110.3304951225527, &
  -1426.1578359574123, -1781.463594396629/

 data coefficient2d/4.6322588,-16.22288877,31.87936853,-50.72504085,72.28843166,&
 -96.26154154,122.42245783,-150.60134813,180.66276676,-212.49539568,246.00562378,&
 -281.1132937,317.748755,-355.85075027,395.36485256,-436.24228262,478.4389941,&
 -521.91495476,566.63357022,-612.56121906,659.66687049,-707.92176748,757.29916201,&
 -807.77409173,859.32318927,-911.92451926,965.55743823,-1020.20247172,1075.84120925,-1132.45620899/
  
 data exponent2d/-2.891592987,-15.235631170,-37.443503370,-69.520142225,-111.466151831,&
  -163.281676504,-224.966764262,-296.521434810,-377.945697336,-469.239556786,-570.40301543,&
  -681.43607529,-802.33873718,-933.11100189,-1073.75286977,-1224.26434126, -1384.64541603,&
  -1554.89609516, -1735.01637786, -1925.00626441,-2124.86575531,-2334.5948504 ,-2554.1935494,&
  -2783.66185288, -3022.99976089,-3272.20727287,-3531.28438956,-3800.23110965,-4079.0474352,-4367.73336409/
  
end block data fpd


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! Walk on green's functions methods
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------!
! Solves 1D first passage problem for walker released at point xo
! in unit interval using walk on line algorithm.
! returns x, the final position which is either
! 0 or 1 and returns the exit time. Unit of time is the time it 
! takes to diffuse <Dx^2> = 1.  t = L^2/(2D)
!
! xo = the initial position 0 < xo < 1 (input)
! x = final position, either 0 or 1    (output)
! time = the exit time                 (output)
!---------------------------------------------------------------
subroutine wol(xo, x, time)

 implicit none
 
 real(8), intent(in) :: xo
 real(8), intent(out) :: x, time
 
 real(8)    :: rndm, displacement, delt
 
 x = xo
 time = 0.0d0
 do while ( x > 0.0d0 .AND. x < 1.0d0 )
    
    displacement = min(x, 1.0d0-x)   
    call vipull(rndm)
    x = x + (1.0d0-2.0d0*anint(rndm))*displacement   
    call randfpd1d(delt)
    time = time + delt*4.0d0*displacement*displacement
!    write (*,*) delt
    
 end do
     
 return
 
end subroutine wol

!------------------------------------------------------------------!
! Solves 2D first passage problem for walker released at point ro in
! annular region using walk on circle algorithm.  capture radius = 1
! escape radius is given as argument. The penetration depth is
! fixed in the code--but is a parameter in the algorithm. The unit
! of time is set by 2Dt=Rcap^2. Since Rcap = 1 when t = 1, 
! 2D = 1.
!
! ro = the initial position 1 < ro < Rmax (input)
! Rmax = the escape radius (input)
! x,y = final position   (output)
! time = the exit time                 (output)
!
! a = inner capture radius
! d = outer capture radius
!---------------------------------------------------------------
subroutine woc(ro, Rmax, x, y, time)

 implicit none
 
 real(8), intent(in) :: ro, Rmax
 real(8), intent(out) :: x, y, time
 
 real(8)    :: r,rnd1,rnd2,mag,a,d,displacement, delt
 
 ! penetration depths
 a = 0.99d0
 d = 1.01d0*Rmax
 
 ! initialize coordinates
 x = ro
 y = 0.0d0
 r = ro
 time = 0.0d0
! write (*,*) a, ro, Rmax, d
 do while ( r > 1 .AND. r < Rmax )
    
!    write (*,*) x, y, time
    
    displacement = min(d-r, r-a)   
!    write (*,*) displacement, d-r, r-a
    ! use trick to get (delx,dely) uniformly distributed on unit circle
    call boxmuller(rnd1,rnd2) !gaussian variates
    mag = dsqrt(rnd1*rnd1+rnd2*rnd2)/displacement
    x = x + rnd1/mag
    y = y + rnd2/mag
    r = dsqrt(x*x+y*y)
    
    ! grab the time step---must scale by displacement^2 
    call randfpd2d(delt)
    time = time + delt*displacement*displacement
    
 end do
     
 return
 
end subroutine woc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! facilitated diffusion methods
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------!
! Returns a hop. Rcap = capture radius. The unit
! of time is set by 2Dt=Rcap^2. I.e. t = Rcap^2/(2D)
! The unit of distance = Rcap
!
! ro = the initial position 1 < ro < Rmax (input)
! Rmax = the escape radius                (input)
! x,y,z = final position                  (output)
! time = the exit time                    (output)
!
!---------------------------------------------------------------
subroutine hop(ro, Rmax, x, y, z, time)

 implicit none
 
 real(8), intent(in) :: ro, Rmax
 real(8), intent(out) :: x, y, z, time
 
 real(8)    :: rnd1,rnd2
 
 call woc(ro, Rmax, x, y, time)
 call boxmuller(rnd1,rnd2)
 z = dsqrt(time)*rnd1
     
 return
 
end subroutine hop

!------------------------------------------------------------------!
! Returns a slide in infinite domain. The unit of time is set by 
! t = 1/koff, where koff is the off rate. So unit of time = mean 
! residence time.  the unit of distance d is the RMS displacement 
! during one residence time, set by 2D/koff = d^2. Therefore unit 
! distance = sqrt(2D/koff). initial x = 0.
!
! x = final position                (output)
! t = the exit time                 (output)
!
!---------------------------------------------------------------
subroutine slide(x, t)

 implicit none
 
 real(8), intent(out) :: x, t
 
 real(8)    :: rnd1,rnd2
 
 call randexp(t)
 call boxmuller(rnd1,rnd2)
 x = dsqrt(t)*rnd1
     
 return
 
end subroutine slide

!------------------------------------------------------------------!
! Returns a slide in a finite interval with reflecting BCs.
! The unit of time is set by time to diffuse the interval, 2Dt = L^2
! therefore unit of time = L^2/(2D). 
! the unit of distance is the interval length.
!
! zo = initial position             (input)
! tau = mean residence time          (input)
! z = final position                (output)
! t = the exit time                 (output)
!
!---------------------------------------------------------------
subroutine slide_ref(zo, tau, z, t)

 implicit none
 
 real(8), intent(in) :: zo, tau
 real(8), intent(out) :: z, t
 
 call randexp(t) 
 t = tau*t  ! rescale to units of diffusion time
 call randpdf1d(zo, t, z)
     
 return
 
end subroutine slide_ref









!------------------------------------------------------------------!
! Returns a slide in a finite interval with reflecting BCs.
! The unit of time is set by time to diffuse the interval, 2Dt = L^2
! therefore unit of time = L^2/(2D). 
! the unit of distance is the interval length.
!
! zo = initial position             (input)
! tau = mean residence time          (input)
! z = final position                (output)
! t = the exit time                 (output)
!
!---------------------------------------------------------------
!subroutine slide_abs(zo, tau, z, t)

! implicit none
 
! real(8), intent(in) :: zo, tau
! real(8), intent(out) :: z, t
 
! call randexp(toff) 
! toff = tau*toff  ! rescale to units of diffusion time

! xo = 0.5d0*zo ! rescale distance to L=2
! call wol(xo,x,tabs)
! tabs = 0.25d0*tabs  ! rescale time back to L=1
 
!# if (toff < tau) then
!#   need to create a random # generator for this problem!
!# else
!   z= 0.0d0
     
! return
 
!end subroutine slide_abs





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! methods for calculating various PDFs and FPDs
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!----------------------------------------------------------------!
! this returns the probability density for a particle diffusing
! in a unit interval with reflecting BCs. unit time is
! time it takes to diffuse <Dx^2> = 1 .
!
! z = the final position, 0 < z < 1    (input)
! zo = the initial position, 0 < zo < 1 (input)
! tau = the time interval   (input)
! answer = the probability density (output)
!----------------------------------------------------------------!
subroutine pdf1d(z, zo, tau, answer)

 implicit none

 real(8), intent(in) :: z, zo, tau
 real(8), intent(out) :: answer

 integer :: N, m
 real(8) :: prefactor

 N = 80 ! number of terms in series representation
 answer = 1.0d0
 do m = 1, N
   prefactor = 2.0d0*dcos(m*3.141592653589793d0*zo)
   prefactor = prefactor*dexp(-4.934802200544679d0*m*m*tau)
   answer = answer + prefactor*dcos(m*3.141592653589793d0*z)
 end do

 return
 
end subroutine pdf1d

!----------------------------------------------------------------!
! this returns the probability density for first passage time
! to the limits of a 1D interval of unit length for a random walker released
! at the center of the interval. The unit of time is the time it
! takes the walker to have an RMS displacement equal to the 
! length of the interval (=1).
!
! tau = the exit time                        (input)
! answer = the exit time probability density (output)
!----------------------------------------------------------------!
subroutine fpd1d(tau, answer)

 implicit none

 real(8), intent(in) :: tau
 real(8), intent(out) :: answer
 
 ! the common block data used for expansion
 real(8) :: coefficient1d(1:10), exponent1d(1:10),coefficient2d(1:30), exponent2d(1:30)
 common /block_fpd/ coefficient1d, exponent1d, coefficient2d, exponent2d
 
 integer :: N, i

 N = 10 ! number of terms in series representation (see block data)
 ! must deal with non-uniform convergence at t = 0. below tauMIN
 ! set answer = 0. tauMIN depends on number of terms in series!
 answer = 0.0d0
 if (tau > 0.0072014d0) then
   do i = 1, N
     answer = answer + coefficient1d(i)*dexp(exponent1d(i)*tau)
   end do
 end if

 return
 
end subroutine fpd1d


!----------------------------------------------------------------!
! this returns the probability density for first passage time
! to the limits of a circular region of radius 1.0 for a particle
! release at the center. The unit of time is the time it
! takes the walker to have an RMS displacement equal to the 
! radius of the circle (=1).
!
! tau = the exit time                        (input)
! answer = the exit time probability density (output)
!----------------------------------------------------------------!
subroutine fpd2d(tau, answer)

 implicit none

 real(8), intent(in) :: tau
 real(8), intent(out) :: answer
 
 ! the common block data used for expansion
 real(8) :: coefficient1d(1:10), exponent1d(1:10),coefficient2d(1:30), exponent2d(1:30)
 common /block_fpd/ coefficient1d, exponent1d, coefficient2d, exponent2d
 
 integer :: N, i

 N = 30 ! number of terms in series representation (see block data)
 ! must deal with non-uniform convergence at t = 0. below tauMIN
 ! set answer = 0. tauMIN depends on number of terms in series!
 answer = 0.0d0
 if (tau > 0.02284) then
   do i = 1, N
     answer = answer + coefficient2d(i)*dexp(exponent2d(i)*tau)
!   write (*,*) answer, coefficient2d(i), exponent2d(i)
   end do
 end if

 return
 
end subroutine fpd2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! methods for generating random #'s from various PDFs and FPDs
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!----------------------------------------------------------------!
! this returns a random sample from the pdf for the 1D interval with
! reflecting BCs. uses rejection method---contains a "bug" that
! invalidates result near boundaries for longer times
!
! zo = the initial position (input)
! tau = the time interval   (input)
! z = the final position (output)
!# there is a problem whenever the slope at z=0 or 1
!# is 0, then the randomly generated value does not
!# seem to follow distribution well. could be problem with
!non-uniform convergence... but

!----------------------------------------------------------------!
subroutine randpdf1d_old(zo, t, z)

 implicit none

 real(8), intent(in) :: zo, t
 real(8), intent(out) :: z
 
 real(8) :: y, p, ymax, rndm

 ! first generate a uniform variant from interval [0,1] and call
 ! pd1d to get max probability value = ymax
 call pdf1d(zo, zo, t, ymax) ! this returns the maximum value of pd1d
 call vipull(rndm)
 y = ymax * rndm
 call vipull(z)
 call pdf1d(z, zo, t, p)
 
 ! repeat until exit condition is reached
 do while (y > p)
 ! below does not actually work, as when slope = 0 at boundaries,
 ! the max value occurs at boundaries.
   call pdf1d(zo, zo, t, ymax) ! this returns the maximum value of pd1d
   call vipull(rndm)
   y = ymax * rndm
   call vipull(z)
   call pdf1d(z, zo, t, p)
   
 end do
 
 return
 
end subroutine randpdf1d_old

!----------------------------------------------------------------!
! this returns a random sample from the pdf for the 1D interval.
! reflecting BCs. Uses gaussian PDf and reflection of points outside
! unit interval---fast and reliable! unit of distance = interval length
! unit of time = time to diffuse interval = L^2/(2D)
!
! zo = the initial position, 0 < zo < 1 (input)
! t = the time interval   (input)
! z = the final position, 0 < z < 1 (output)
!----------------------------------------------------------------!
subroutine randpdf1d(zo, t, z)

 implicit none

 real(8), intent(in) :: zo, t
 real(8), intent(out) :: z
 
 real(8) :: rnd1, rnd2, half, delx

 ! generate a gaussian random #
 call boxmuller(rnd1,rnd2)
 delx = rnd1*dsqrt(t)
 z = zo + delx
 
 ! now fold over position into unit interval
 half = 0.5d0*z
! write(*,*) rnd1, half, floor(half)
 z = 2.0d0*(half-floor(half))
 if (z>1.0) then
   z = 2.0d0-z 
 end if
 
 return
 
end subroutine randpdf1d

!----------------------------------------------------------------!
! this returns a random sample from the fpd for the 1D exit problem
! particle released at center (x=0.5)
!
! time = returned sample (output)
!----------------------------------------------------------------!
subroutine randfpd1d(time)

 implicit none

 real(8), intent(out) :: time
 
 real(8) :: y, f, ymax, rndm, t

 ! first generate a random exponentially distributed number and calculate
 ! the max value which is greater than the fpd value at that point
 call randexp(t)
 ymax = 6.284d0 * dexp(-1.0d0*t) ! this is the upper bounding function
! write (*,*) t
 
 ! now generate a potential return time and y value, 0 < y < ymax
 call vipull(rndm)
 y = ymax * rndm
 time = t / 4.934802200544679 ! this uses exponent1d(1) value
 call fpd1d(time,f)
 
 ! repeat until exit condition is reached
 do while (y > f)
   call randexp(t)
   ymax = 6.284d0 * dexp(-1.0d0*t)
   call vipull(rndm)
   y = ymax * rndm
   time = t / 4.934802200544679 
   call fpd1d(time,f)
 end do
 
 return
 
 end subroutine randfpd1d

!------------------------------------------------------
! returns a random number picked from an exponential
! distribution with unit mean.
!------------------------------------------------------
subroutine randexp(answer)

 implicit none
 
 real(8), intent(out) :: answer
 
 real(8) :: rndm

 call vipull(rndm)
 answer = -dlog(rndm)
! write (*,*) rndm, answer
 
 return

end subroutine


!----------------------------------------------------------------!
! this returns a random sample from the fpd for the 2D exit problem
! particle released at center of circle with R=1. unit of time
! is time it takes to diffuse <Dx^2> = 1
! uses rejection method with a double exponential envelope fn.
!
! time = returned sample (output)
!----------------------------------------------------------------!
subroutine randfpd2d(t)

 implicit none

 real(8), intent(out) :: t
 
 real(8) :: y, f, ymax, rndm, t1, t2

 ! first generate a random distributed number from the envelope fn
 ! ymax(t) = 3.2*(exp(-t/0.346)-exp(-t/0.15))/(0.346+0.15)
 call randexp(t1)
 call randexp(t2)
 t = 0.346d0*t1 + 0.15d0*t2
 
 ! now determine envelope value and the fpd value at that time
 ymax = 6.451612903d0*( dexp(t/(-0.346d0))-dexp(t/(-0.15d0)) )
 call fpd2d(t,f)
 
 ! now generate a random y value, 0 < y < ymax
 call vipull(rndm)
 y = ymax * rndm
 
 ! repeat until exit condition is reached
 do while (y > f)
 
   ! first generate a random distributed number from the envelope fn
   ! ymax(t) = 3.2*(exp(-t/0.346)-exp(-t/0.15))/(0.346+0.15)
   call randexp(t1)
   call randexp(t2)
   t = 0.346d0*t1 + 0.15d0*t2
 
   ! now determine envelope value and the fpd value at that time
   ymax = 6.451612903d0*( dexp(t/(-0.346d0))-dexp(t/(-0.15d0)) )
   call fpd2d(t,f)
 
   ! now generate a random y value, 0 < y < ymax
   call vipull(rndm)
   y = ymax * rndm
 
 end do
 
 return
 
 end subroutine randfpd2d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! core methods
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-----------------------------!

subroutine random(rand)

 implicit none

 real(8), intent(out) :: rand

 real(8) :: vi(0:31), yi
 integer(8) :: randint
 common/block_random/vi,yi,randint

 real(8) :: dmul

 dmul = 0.5d0/dble(2*(2_8**62-1)+1)
 randint = 2862933555777941757_8*randint+1013904243
 rand = 0.5d0+dble(randint)*dmul
 
 return

end subroutine random

!-----------------------------!

subroutine randseed(nseed)

 implicit none

 integer(8), intent(in) :: nseed
 
 real(8) :: vi(0:31), yi
 integer(8) :: randint
 common/block_random/vi,yi,randint

 randint = nseed

 return
 
end subroutine randseed

!-----------------------------!

subroutine viseed()

 implicit none

 real(8) :: vi(0:31), yi
 integer(8) :: randint
 common/block_random/vi,yi,randint

 real(8) :: rnum
 integer :: ii
 
 do ii = 0, 31
   call random(rnum)
   vi(ii) = rnum
 end do
 call random(rnum)
 yi = rnum
 
 return

end subroutine viseed


!-----------------------------!

subroutine vipull(vivalue)

 implicit none

 real(8), intent(out) :: vivalue
 
 real(8) :: vi(0:31), yi
 integer(8) :: randint
 common/block_random/vi,yi,randint

 real(8) :: randreal
 integer :: ii
 
 ii = yi*32
 vivalue = vi(ii)
 yi = vi(ii)
 call random(randreal); vi(ii) = randreal
 
 return

end subroutine vipull


!-------------------------------!
subroutine boxmuller(rand1,rand2)
 
 implicit none

 real(8), intent(out) :: rand1, rand2
 
 real(8) :: junk1, junk2

 call vipull(junk1)
 call vipull(junk2)
 rand1 = dsqrt(-2.0d0*dlog(junk1))*dcos(6.28318530717958656*junk2)
 rand2 = dsqrt(-2.0d0*dlog(junk1))*dsin(6.28318530717958656*junk2)
 
 return

end subroutine boxmuller
 

