      program wavefunction
      parameter (n=20, dt=0.01, pi=3.141592)
      real phi(n)
      aphi_0=1.0
      anu=0.01
      ak=2.*pi ! энэ хэмжигдэхүүнүүдийг өөрчилнө. 
      c=1.0
      dx=1./n ! торын алхам
      atime=5.0
      iter=atime/dt

      open(10,file='data.dat')
      do i=0,iter
        a=dt*i
        do j=0,n
        ax=j*dx
        phi(j)=aphi_0*exp(-anu*ak**2*a)*cos(ak*(ax-c*a))
      if(mod(i,100).eq.0) write(10,*) j*dx, phi(j)
        end do
      end do  
      stop
      end program  
