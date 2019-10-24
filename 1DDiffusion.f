      program diffusionLBM
        parameter (m=100)              ! m is the number of lattice nodes
        real f1(0:1), f2(0:m), rho(0:m), feq(0:m), x(0:m)
        integer i
        open(2,file='resultLBM.dat')
        dt=1.0
        dx=1.0
        x(0)=0.0
! lattice arrangement
        do i=1,m
           x(i)=x(i-1)+dx
        end do
        csq=dx*dx/(dt*dt)
        alpha=0.25                     ! diffusive coefficient
        omega=1.0/(alpha/(dt*csq)+0.5) ! eq.3.23
        mstep=200                      ! total number of time steps
        twall=1.0                      ! left hand wall temperature
! initial condition
        do i=0,m
         rho(i)=0.0                     ! initial value of the domain temperature
         f1(i)=0.5*rho(i)
         f2(i)=0.5*rho(i)
        end do
        do kk=1,mstep
! main loop
! collision process
          do i=0,m
           rho(i)=f1(i)+f2(i)           ! dependent variable d.function
           feq(i)=0.5*rho(i)            ! equilibrium distribution function 
! since k1=k2=0.5, then feq1=feq2=feq
           f1(i)=(1.-omega)*f1(i)+omega*feq(i)
           f2(i)=(1.-omega)*f2(i)+omega*feq(i)
          enddo
! streaming process
          do i=1,m-1
           f1(m-i)=f1(m-i-1)   ! f1 streaming
           f2(i-1)=f2(i)       ! f2 streaming
          enddo
! boundary condition
          f1(0)=twall-f2(0)   ! constant temperature boundary condition, x=0
          f1(m)=f1(m-1)       ! adiabatic boundary condition, x=L
          f2(m)=f2(m-1)       ! adiabatic boundary condition, x=L
        enddo
! end of loop
        do i=0,m
          write(2,*) x(i), rho(i)
        enddo
        stop
        close(2)
        end program
