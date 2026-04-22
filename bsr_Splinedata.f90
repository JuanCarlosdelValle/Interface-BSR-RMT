!======================================================================
!     UTILITY       BSR _SPLINEDATA                       version 1
!
!               C O P Y R I G H T -- 2024
!
!     Written by:   Juan C. del Valle
!                   email: jcdvaller@gmail.com 
!
!======================================================================
!   based on knot.dat, it creates the Splinedata file
!======================================================================
!
!      RMATRIXII: nc,k,nb,ndim,nt
!
!
!     INPUT FILES:
!
!      knot.dat      - B-spline information     
!   
!    In brackets the name of the variable in RMATRIXII is displayed. 
!
!        nv-1(nc) : the number of intervals.
!           ks(k) : order of splines.
!          ns(nb) : total number of splines. 
!      ns-1(ndim) : maximum principal quantum number of the orbitals. According to Hugo: max number of functions per l.
!      ns+ks (nt) : number of knots/the length of the knot set.
!
! 
!     OUTPUT FILES:
!
!     Splinedata     - unformatted file
!----------------------------------------------------------------------
      Use spline_param; Use spline_atomic; Use spline_grid; Use spline_orbitals

      Integer :: nc, ndim, nt
  
!----------------------------------------------------------------------
! ... sets up grid points and initializes according to the file
! ... "knot.dat" :

       Call define_grid(z);       
        ndim = ns - 1  
        nt = ns + ks    
        nc = nv + 1  !WARNING: the " + 1" is inserted to follow RMATRIXII convention.      
 
!----------------------------------------------------------------------
       
       Open(unit=97,file='Splinedata',form='unformatted')
       Open(unit=997,file='Splinedata-fmt',form='formatted')

       Write(97) nc, ks, ns, ndim, nt
       Write(97) (t(i),i=1,nt)        

       Write(997,9997) nc, ks, ns, ndim, nt   
9997 Format(10i5)

       Write(997,9998) (t(i),i=1,nt)      
9998 Format(1p10e11.3)
       
       close(unit=97)
       close(unit=997) 
      
      End   !  utility BSR _SPLINEDATA


