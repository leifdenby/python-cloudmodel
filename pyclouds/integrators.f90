program integrators
   implicit none

   integer, parameter :: neqn = 1
   real(8) :: y(neqn) = 1.0, x = 0.0, dx = 0.1, x_max = 10.0, max_abs_err = 0.0
   integer :: i = 0, m = 0
   character(len=100) :: msg = ""

   do while (x < x_max)
      i = i + 1
      m = 0
      call rkf34_original(dydx_f, y, x, dx, neqn, msg, m, max_abs_err)
      print *, i, x, y, dx, max_abs_err, m

      if (msg(1:1) /= " ") then
         print *, msg
         stop(0)
      endif
   enddo
   
   contains
      function dydx_f(x, y)
         real(8), intent(in) :: x, y(neqn)
         real(8), dimension(neqn) :: dydx_f

         !dydx_f = 1. + y**2.0
         dydx_f = -2.0*y + exp(-2.*(x-6.)**2.)
      end function dydx_f

   recursive subroutine rkf34_original(dydx, y, x, dx, n, msg, m, max_abs_err)
      integer, intent(in) :: n
      character(len=100), intent(out) :: msg
      real(8), intent(inout) :: y(n), dx
      real(8), intent(inout) :: x
      real(8), intent(out) :: max_abs_err
      integer, intent(inout) :: m

      interface
         function dydx(x, y)
            real(8), intent(in) :: x, y(1)
            real(8) :: dydx(1)
         end function
      end interface
      !f2py raise_python_exception msg

      !real(8), parameter :: &
         !a2=2./7.,   b21=2./7., &
         !a3=7./15.,  b31=77./900.,   b32= 343.900, &
         !a4=35./38., b41=805./1444., b42=-77175./54872., b43=97125./54872., &
         !a5=1.0,     b51=79./490.,   b52= 0.0,           b53=2175./3616.,   b54=2166./9065.
      !real(8) :: &
         !c1_1=79./490.,    c2_1=0.0,  c3_1=2175./36.26,  c4_1=2166./9065.,&
         !c1_2=229./1470.,  c2_2=0.0,  c3_2=1125./1813.,   c4_2=13718./81585., c5_2=1./18.

      real(8), parameter :: &
         a2=0.5, b21=0.5, &
         a3=0.5, b31=0.0,        b32= 0.5, &
         a4=1.0, b41=0.0,        b42= 0.0,           b43=1.0, &
         a5=1.0, b51=1./6.,      b52= 1./3.,         b53=1./3.,   b54= 1./6.
      real(8) :: &
         c1_1=0.0, c2_1=0.0,  c3_1=0.0,  c4_1=0.0,  c5_1=1.,&
         c1_2=1./6., c2_2=1./3.,    c3_2=1./3.,     c4_2=0.,   c5_2=1./6.

      real(8) :: k1(n), k2(n), k3(n), k4(n), k5(n)
      real(8) :: abs_err(n), y_n1(n), y_n2(n), s

      ! TODO: move these into an input variable
      real(8) :: abs_tol=1.0e-3, rel_tol=1.0e-3, max_rel_err, err

      logical :: done = .false.

      done = .false.

      k1 = 0.0
      k2 = 0.0
      k3 = 0.0
      k4 = 0.0
      k5 = 0.0

      k1 = dx*dydx(x,       y)
      k2 = dx*dydx(x+a2*dx, y+b21*k1)
      k3 = dx*dydx(x+a3*dx, y+b31*k1 + b32*k2)
      k4 = dx*dydx(x+a4*dx, y+b41*k1 + b42*k2 + b43*k3)
      k5 = dx*dydx(x+a5*dx, y+b51*k1 + b52*k2 + b53*k3 + b54*k4)

      y_n1 = y + c1_1*k1 + c2_1*k2 + c3_1*k3 + c4_1*k4 + c5_1*k5
      y_n2 = y + c1_2*k1 + c2_2*k2 + c3_2*k3 + c4_2*k4 + c5_2*k5

      !abs_err = abs(y_n1 - y_n2)
      abs_err = abs(1./6.*(k4 - k5))

      max_abs_err = maxval(abs_err)
      max_rel_err = maxval(abs_err/y_n2)

      s = 0.84*(rel_tol*dx/max_abs_err)**0.25


      if (max_rel_err < rel_tol) then
         done = .true.
      endif 

      if (done) then
         y = y_n2
         x = x + dx

         ! if the suggested scaling makes the step larger then lets use it
         if (s > 1.0) then
            dx = dx*s
         endif
      else
         if (dx < 1.0e-32) then
            msg = "step size became very small"
         else if (m > 5000) then
            msg = "Didn't converge"
         else
            if (s > 1.0) then
               !print *, "s=", s
               !msg = "Incorrect scaling, timestep is growing."
               s = 0.1
            endif
            dx = dx*s
            m = m+1
            call rkf34_original(dydx, y, x, dx, n, msg, m, max_abs_err)
         endif
      endif

   end subroutine rkf34_original
end program

