!======================================================================
!     UTILITY       BSR _SPLINEWAVES                       version 1
!
!               C O P Y R I G H T -- 2024
!
!     Written by:   Juan C. del Valle
!                   email: jcdvaller@gmail.com 
!
!======================================================================
!   based on rsol.nnn files, it creates the Splinewaves file
!======================================================================
!
!     INPUT FILES:
!
!      rsol.nnn      - file  with coefficients of the orbitals in the  B-Spline representation     
!                      (the first coefficient is removed here to fit the RMII style)
!
!    In brackets the name of the variable in RMATRIXII is displayed. 
!
!     OUTPUT FILES:
!
!     Splinewaves     - unformatted file
!-------------------------------------------------------------------------
!     WARNING: CURRENTLY THIS UTILITY DOES NOT WORK IF PERTURBERS ARE USED
!----------------------------------------===------------------------------
!  USAGE: bsr_Splinewaves nlsp=Number_of_Partial_Waves
!-------------------------------------------------------------------------
      Implicit real(8) (A-H,O-Z)

      Character(80) :: AF_r = 'rsol.nnn',  AF_all = 'Splinewaves' 
      Character(80) :: AF, AF1,AF2
      Integer, allocatable :: nchn(:), nstk(:)
      Real(8), allocatable :: v(:)
      Real(8), allocatable :: vred(:) ! Reduced WF in B-Spline representation, where the first coefficient/spline is removed.
      Integer, allocatable :: skpL(:)  ! List with the positions of the first coefficients.

      iarg = COMMAND_ARGUMENT_COUNT()
      if(iarg.gt.0) Call GET_COMMAND_ARGUMENT(1,AF)

      if(AF.eq.'?') then
        write(*,'(/a)') 'collect_rsol merges the set of rsol.nnn file in one Splinewaves'
        write(*,'(/a)') '     rsol.nnn  -->  Splinewaves'
        write(*,'(/a)') 'Call as:  collect_rsol  nlsp=..,  where nlsp the number of partial waves'
        write(*,'(/a)') 
        Stop 
      end if        
                                                        
      Call Read_iarg('nlsp',nlsp)
      if(nlsp.le.0) Stop 'nlsp <= 0'

      iout = 2
      open(iout,file=AF_all,form='UNFORMATTED')

      Allocate(nchn(nlsp), nstk(nlsp))

      Do ilsp = 1,nlsp

       write(AF,'(a,i3.3)') 'rsol.',ilsp
       Call Check_file(AF)
       inp=1; open(inp,file=AF,form='UNFORMATTED')

       read(inp)  nhm,khm,kch,kcp,ns 
       read(inp)  (e,i=1,khm)

       Allocate(v(kch*ns))

!       INITIALIZING THE REDUCED ARRAY
         if (allocated(vred)) then
            deallocate(vred)
         end if
         allocate(vred(kch*(ns-1)))

!       INITIALIZING LIST WITH POSITIONS 
         if (allocated(skpL)) then
            deallocate(skpL)
         end if
         allocate(skpL(kch))

       do ich =1, kch
          skpL(ich) = (ich-1)*ns+1    
       end do

       Do  j = 1,khm 
           read(inp) v

!          INSERTION TO REMOVE THE FIRST SPLINE

           ishift = 0
             
           Do k = 1, kch*(ns-1)
              if (any(skpL == k)) then 
                  ishift = ishift + 1 
              end if   
                  vred(k) = v(k+ishift) 
           End do
!
       write(iout) (vred(i),i=1,kch*(ns-1))
       End do 
       Deallocate(v)
       Deallocate(vred)

       nchn(ilsp) = kch
       nstk(ilsp) = khm

      End do

      ipri = 3; open(ipri,file='bsr_Splinewaves.log')

      write(ipri,'(a)') 'nlsp (inast), number of partial waves:'
      write(ipri,'(i5)') nlsp

      write(ipri,'(a)') 'bspl_ndim:'
      write(ipri,'(i5)') ns

      nchmx = maxval(nchn)
      write(ipri,'(a)') 'nchmx, maximum number of channels:'
      write(ipri,'(i5)') nchmx

      nstmx = maxval(nstk)
      write(ipri,'(a)') 'nstmx, maximum number of RM solutions:'
      write(ipri,'(i5)') nstmx

      write(ipri,'(a)') 'nchn(1:nlsp) - number of channels for each partial wave:'
      write(ipri,'(10i10)') nchn

      write(ipri,'(a)') 'nstk(1:nlsp) - number of RM solutions for each partial wave:'
      write(ipri,'(10i10)') nstk


      End ! program

      

