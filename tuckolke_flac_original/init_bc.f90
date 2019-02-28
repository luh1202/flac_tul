!
!    Variable Boundary conditions (Initialization)
!
       subroutine init_bc 

       include 'precision.inc' 
       include 'arrays.inc' 
       include 'params.inc'
!
!      ---- bc(nh*2)-for each nodal degree-of-freedom this is assigned zero
!      value,unless a boundary condition is applied,in which case it is
!      diven the value of the boundary condition.
!
!      ---- ncod(nh*2)-the entries of this array are coded zero for a nodal
!      degree-of-freedom with no applied boundary conditions ,and unity
!      for an applied boundary condition.
! 
       bc = 0.
       ncod = 0
       bcstress = 0.
       nopbou = 0  
       ncodbou = 0
       nebou = 0
!---------Number of side (nofside) -------------------------------
!              4
!       !----------------!
!     1 !                ! 3   nz
!       !                !
!       !----------------!
!              2
!              nx
!-----------------------------------------------------------------

!  nofbc        number of boundary conditions         
!  nofside      1-left,2-bottom,3-right,4-top 
!  nbc1,nbc2    range of application of boundary condition 
!  nbc          type of boundary condition:  0-no,1-vel,2-stress,3-periodic 
!  0  - no condition 
!  10 - velx           01 - vely 
!  20 - normal stress  02 - shear stress 
!  30 - vely           40 - shear stress out of the plane
! ---------------------------------------------------------------- 
! function of boundary conditions: 
!   a + bx + cx**2 + d cos (2pi ex) + f sin (2pi gx) 
! where x is undimensional. i.e. x = (x - x(na1)) / (x(na2) - x(na1)) 
! ----------------------------------------------------------------
! nnop is numeration pointer  only for boundary arrays
! as nebou, bcstress, nopbou, ncodbou 

       nnop = 0 
       do 1 i = 1, nofbc
         if(nbc(i).eq.nbc(i-1).and.nofside(i).eq.nofside(i-1)) then  
!           ndbc1 = ndbc				!G.Ito removed 3/17/06
!           ndbc = (ndbc)+nbc2(i) - nbc1(i) + 1
	    if (nbc(i).gt.10.and.nofside(i).ne.1) then
       	      write(*,*) 'Mixed b.c. only implemented for velocity on left side:  G. Ito 3/17/06'
	      stop
	    endif
         else
           ndbc = nbc2(i) - nbc1(i) +1
         endif
!
! Number of bc for stresses
!
         if (nbc(i) .eq. 2 .or. nbc(i) .eq. 20.or.nbc(i).eq.40 ) then
           ndbc = ndbc - 1 
           if(nbc2(i).lt.nz) ndbc = ndbc +1
         endif

         if ( ndbc .le. 0 ) write (6,*)'***NOT correct range of BC**' 
          
!
!-------------------------------------------------------------
!       left   (only for left node)
!-------------------------------------------------------------
!
         if ( nofside(i) .eq. 1 ) then

                 x1 = cord ( nbc1(i),1, 2)
                 x2 = cord ( nbc2(i),1, 2)

               do 2 n = 1,ndbc 
                  k = i
!   if(nbc(i).eq.nbc(i-1).and.nofside(i).eq.nofside(i-1).and.n.le.ndbc1)then !G.Ito 3/17/06
!                  k=i-1
!                  endif
                  numbp  = n + nbc1(i)-1
                  numbp1 = numbp+1 
                  x  = (cord (numbp,1,2)  - x1)/(x2-x1)   
                  xn = (cord (numbp1,1,2) - x1)/(x2-x1)   
                  xa = 0.5 * (x+xn)
           if (nbc(i).eq.1.or.nbc(i).eq.10.or.nbc(i).eq.30)  & 
          call velbc (k,numbp,x)  

           if (nbc(i).eq.2.or.nbc(i).eq.20.or.nbc(i).eq.40) then
                   nnop = nnop + 1
           call stressbc (i,nnop,numbp,numbp1,xa)  
                endif

 2             continue
         endif

!-------------------------------------------------------------
!        bottom
!-------------------------------------------------------------
           if ( nofside(i) .eq. 2 ) then

                 x1 = cord ( nz, nbc1(i), 1)
                 x2 = cord (nz, nbc2(i), 1)

               do 3 n = 1,ndbc 
                  k = i
   if(nbc(i).eq.nbc(i-1).and.nofside(i).eq.nofside(i-1).and.n.le.ndbc1)then 
                  k=i-1
                  endif
                  numbp  = n   
                  numbp1 = n + 1 
                  x = (cord (nz,numbp,1) - x1)/(x2-x1)   
                  xn = (cord (nz,numbp1,1) - x1)/(x2-x1)   
                  xa = 0.5 * (x+xn)

           if (nbc(i).eq.1.or.nbc(i).eq.10.or.nbc(i).eq.30)  & 
           call velbc (k,numbp,x)  
           if (nbc(i).eq.2.or.nbc(i).eq.20.or.nbc(i).eq.40) then
                   nnop = nnop + 1
           call stressbc (i,nnop,numbp,numbp1,xa)   
                endif

 3             continue
           endif

!-------------------------------------------------------------
!        right  (only for the right node)
!-------------------------------------------------------------

           if ( nofside(i) .eq. 3 ) then

                 x1 = cord ( nbc1(i),nx, 2)
                 x2 = cord ( nbc2(i),nx, 2)

               do 4 n = 1,ndbc 
                  k = i
   if(nbc(i).eq.nbc(i-1).and.nofside(i).eq.nofside(i-1).and.n.le.ndbc1)then 
                  k=i-1
                  endif

                  numbp  =  n  
                  numbp1 = n + 1 

                  x  = (cord (numbp,nx,2) - x1)/(x2-x1)   
                  xn = (cord (numbp1,nx,2) - x1)/(x2-x1)   
                  xa = 0.5 * (x+xn)

           if (nbc(i).eq.1.or.nbc(i).eq.10.or.nbc(i).eq.30)  & 
          call velbc (k,numbp,x)  

           if (nbc(i).eq.2.or.nbc(i).eq. 20.or.nbc(i).eq.40) then
                   nnop = nnop + 1
           call stressbc (i,nnop,numbp,numbp1,xa)      
                endif

 4             continue

           endif
!-------------------------------------------------------------
!        top
!-------------------------------------------------------------
           if ( nofside(i) .eq. 4 ) then


                 x1 = cord ( 1,nbc1(i), 1)
                 x2 = cord ( 1,nbc2(i), 1)

               do 5 n = 1,ndbc 
                  k = i
   if(nbc(i).eq.nbc(i-1).and.nofside(i).eq.nofside(i-1).and.n.le.ndbc1)then 
                  k=i-1
                  endif

                  numbp  = n 
                  numbp1 = n + 1 

                  x = (cord (1,numbp,1) - x1)/(x2-x1)   
                  xn = (cord (1,numbp1,1) - x1)/(x2-x1)   
                  xa = 0.5 * (x+xn)


           if (nbc(i).eq.1.or.nbc(i).eq.10.or.nbc(i).eq.30)   & 
           call velbc (k,numbp,x)  

           if (nbc(i).eq.2.or.nbc(i).eq.20.or.nbc(i).eq.40) then
                   nnop = nnop + 1
           call stressbc (i,nnop,numbp,numbp1,xa)       
                endif

 5             continue

           endif


 1     continue
!
! Maximum node with applied stress
!
                   nopbmax = nnop
!
! Maximum applied velocity and/or 
! Convert constant strain INTERNAL B.C.
! to real velocities 
!
       call vbcal

       return
       end

!----------------------------------------------------------------
       subroutine stressbc (i,n,numbp,numbp1,x)      
      
       include 'precision.inc'
       include 'arrays.inc'
       include 'params.inc'
  
       pi2 = 2. * 3.14159
 
       fun =  bca(i) + bcb(i)*x + bcc(i)*x*x  & 
     + (bcd(i)*cos (pi2*bce(i)*x) + bcf(i)*sin (pi2*bcg(i)*x))  &
      *exp(-((x-bci(i))*bch(i))**2 )
       if (nofside(i).eq.1) then
       nopbou (n,1) = 1 
       nopbou (n,2) = 1 
       nopbou (n,3) = numbp
       nopbou (n,4) = numbp1
       endif 
       if (nofside(i).eq.2) then
       nopbou (n,1) = numbp 
       nopbou (n,2) = numbp1 
       nopbou (n,3) = nz 
       nopbou (n,4) = nz 
       endif 
       if (nofside(i).eq.3) then
       nopbou (n,1) = nx 
       nopbou (n,2) = nx 
       nopbou (n,3) = numbp 
       nopbou (n,4) = numbp1 
       endif 
       if (nofside(i).eq.4) then
       nopbou (n,1) = numbp 
       nopbou (n,2) = numbp1 
       nopbou (n,3) = 1 
       nopbou (n,4) = 1 
       endif 
! - normal component
             if (nbc(i) .eq. 20 ) then
       ncodbou   (n,1) = 1
       bcstress  (n,1) = fun
             endif
 
! - shear component
             if (nbc(i) .eq. 2 ) then
       ncodbou  (n,2) = 1
       bcstress (n,2) = fun
             endif

! - shear component out of plane
!            if (nbc(i) .eq. 40 ) then
!      ncodbou  (n,3) = 1
!      bcstress (n,3) = fun
!             endif

           return
           end

!----------------------------------------------------------------
       subroutine velbc (i,numbp,x)  

       include 'precision.inc' 
       include 'arrays.inc'
       include 'params.inc'
                                  
       pi2 = 2. * 3.14159 
 
       fun =  bca(i) + bcb(i)*x + bcc(i)*x*x   & 
     + (bcd(i)*cos (pi2*bce(i)*x) + bcf(i)*sin (pi2*bcg(i)*x))    &
     *exp(-((x-bci(i))*bch(i))**2)

       if (nofside(i).eq.1) then
       ii1 = 1 
       jj1 = numbp
       endif 
       if (nofside(i).eq.2) then
       ii1 = numbp 
       jj1 = nz 
       endif 
       if (nofside(i).eq.3) then
       ii1 = nx 
       jj1 = numbp 
       endif 
       if (nofside(i).eq.4) then
       ii1 = numbp 
       jj1 = 1 
       endif 
! - x component 
             if (nbc(i) .eq. 10 ) then  
       ncod(jj1,ii1,1) = 1
!       if (abs(bc(jj1,ii1,1)).gt.0.) then  !G.Ito removed 3/17/06
!        fun = bc(jj1,ii1,1)
!       endif  
       bc(jj1,ii1,1) = fun  
             endif 
 
! - z component
             if (nbc(i) .eq. 1 ) then
       ncod(jj1,ii1,2) = 1
       bc  (jj1,ii1,2) = fun
!       write(*,*) ncod(jj1,ii1,2),bc(jj1,ii1,2)
             endif

! - y component 

!            if (nbc(i) .eq. 30 ) then
!       ncod(numbp,3) = 1
!       bc  (numbp,3) = fun
!        write(*,*) numbp,ncod(numbp,3),bc(numbp,3)
!             endif

           return
           end

!
! Find Maximum velocity applied to the boundaries
!
       subroutine vbcal 

       include 'precision.inc'
       include 'arrays.inc'
       include 'params.inc'

       vbc = 0.

       do 1 i = 1,nx
       do 1 j = 1,nz 
       do 1 k = 1,2
       if (ncod (j,i,k) .eq.1.or.ncod(j,i,k).eq.10.or.ncod(j,i,k).eq.30)  &
       vbc = max (vbc, abs(bc(j,i,k)) ) 
 1     continue
       
       if(vbc.eq.0.) vbc=3.e-10

       open(13,file = 'vbc.s')
       write(13,*) vbc
       close(13)
       if (ny_inject.eq.3) then
           vbc = vbc*(1-rate_inject)   
       endif
       return 
       end








