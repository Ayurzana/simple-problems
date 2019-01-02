! Шүүрэлтийн тооцоо хийх фортран програм
! Кодолсон Б.Аюурзана 2018 он
! Лекцийн бодлого 1
! Ус үл нэвтрүүлэх буурьтай, нэгэн төрлийн шороон боомт
      program Shuurelt
      parameter (m=1000000, alh=0.0001)
      real Hb, Ht, bb, m1, m2, deltl, btaw, ntaw
      real kq, lsh, vgar, vzuw
      real h1, h(m), l(m), ho
      real qsh
	  
      print*, "Boomtiin undur, Hb [m]"
      read(*,*) Hb
      print*, "Boomtiin turelt, Ht [m]"
      read(*,*) Ht
      print*, "Boomtiin hyriin urgun, bb [m]"
      read(*,*) bb
      print*, "Boomtiin deed dood naluu, m1, m2"
      read(*,*) m1, m2
      print*, "Tawtsangiin urgun, too- bt, nt"
      read(*,*) btaw, ntaw
      print*, "Ikh biye shuureltiin coef, kq [m/s]"
      read(*,*) kq
      print*, "Dood hashits dahi usnii gun, ho [m]"
      read(*,*) ho

      open(33,file='xy.dat')
      open(32,file='damsize.dat')

! Бодолт
      ! Босоо шугамын байрлал
      deltl=(m1/(2.*m1+1.))*Ht
      print*, "Tenhlegiin bairlal d_L", deltl, "[m]"
      ! тооцооны хэсгийн урт
      lsh=deltl+m1*(Hb-Ht)+bb+btaw*ntaw+Hb*m2
      print*, "Tootsoonii hesgiin urt", lsh, "[m]"
      h1=lsh/m2-sqrt((lsh/m2)**2-(Ht-ho)**2)
      print*, "Shuureltiin usnii garah undur, h1", h1
      ! шүүрэлтийн усны зарцуулга
      qsh=kq*(Ht**2-h1**2)/(2.*(lsh-m2*h1))
	  
      print*, "Huwiin zartsuulga q=", qsh, "[m3/s.m]"
	  
      l=0.
      h(1)=0.
	  do i=1,m
         h(i)=sqrt(Ht**2-2.*qsh*l(i)/kq)
             write(33,*) l(i)+ht*m1-deltl, h(i)
		 if(l(i).ge.lsh-m2*h1) then
             i_gar=i
             go to 22
         else
             l(i+1)=l(i)+100.*alh
         end if
      end do
!   22 vgar=kq*(h(i_gar-10)-h(i_gar-20))/(l(i_gar-20)-l(i_gar-10))
   22 vgar=
     &  kq*(h(i_gar)-h(i_gar-i_gar/100))/(l(i_gar-i_gar/100)-l(i_gar))
         print*, "Iteration number", i_gar
      vzuw=sqrt(kq/15./m2)

         if(vgar.gt.vzuw) then
      print*, "Siiregjilt uusne. Garaltiin hurd=", vgar, "[m/day]"
      print*, "              Zuwshuurugduh hurd=", vzuw, "[m/day]"
         else
      print*, "Siiregjilt uusehgui. Garaltiin hurd=", vgar, "[m/day]"
      print*, "                 Zuwshuurugduh hurd=", vzuw, "[m/day]"
         end if
       write(33,*) lsh+ht*m1-deltl, ho
! Геометр бичилт
      write(32,*) ht*m1, ht
      write(32,*) 0, ht
      write(32,*) 0, 0
      write(32,*) hb*m1, hb
      write(32,*) hb*m1+bb, hb
      write(32,*) hb*m1+bb+hb*m2, 0	   
	   
       call system("gnuplot -plot xy.plt")
      close (33)

      close (32)
      end program Shuurelt
