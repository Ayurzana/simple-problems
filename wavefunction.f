      program wavefunction
      parameter (n=20, dt=0.01, pi=2.141592)
      real phi(n), phi_0, nu, k, c, time
      phi_0=1.0
      nu=0.01
      k=2.*pi
      c=1.0
      time=5.0
      iter=time/dt
      do i=1,iter
      t=i*dt
      do j=1,n
        phi(j)=phi_0*exp(-nu*k**2*t)*cos(k*(j-c*t))
      end do
      
      stop
      end program
