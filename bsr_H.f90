!======================================================================
!     UTILITY       BSR_H                       version 1
!
!               C O P Y R I G H T -- 2024
!
!     Written by:   Juan C. del Valle
!                   email: jcdvaller@gmail.com 
!
!======================================================================
!   Based on h.nnn, it creates the H file in the same way RMII does
!   IN RMII, the subroutines write1, write2, writcf, and write3 
!   (see source/ham/farm_intfce.f90) are in charge of writting the H-file
!
!   RMII input: farm =.true.  
!
!======================================================================
!
!     INPUT FILES:
!
! h.nnn - Diagonalized Hamiltonian matrix data for the partial wave nnn
!
!     OUTPUT FILES:
!
!    H - unformatted file
!----------------------------------------------------------------------

      Implicit real(8) (A-H,O-Z)
 
      Real(8),allocatable :: ENAT(:), ENAT1(:), VALUE(:)
      Integer,allocatable :: LAT(:), LAT1(:), ISAT(:), ISAT1(:), &
                             NCONAT(:), L2P(:), IPAT(:), IPAT1(:), jkch(:)
      Real(8),allocatable :: COEFF(:,:), WMAT(:,:), CF(:,:,:)  

      Integer :: ih1    = 1
      Integer :: ih2    = 999
      Integer :: MRANG2 = 201 !WHAT IS THIS?
      Integer :: folder = 0
      Integer :: LRANG2 = 0 
      Character(80) :: AF
      Character(len=100) :: errmsg
! ... Check the arguments:

      Call get_command_argument(1,AF)  
      if(AF.eq.'?') then    !  help section 
       write(*,*)
       write(*,*) 'bsr_H merges the different h.nnn files: '
       write(*,*)
       write(*,*) 'h.001 + h.002 + ... + h.nnn --->  H'
       write(*,*)
       write(*,*) 'Call as:   bsr_H  [ih1=.. ih2=.. folder=..]'
       write(*,*)
       write(*,*) 'with default values for partial h.nnn: ih1=1, ih2=999'
       write(*,*)
       write(*,*) 'if folder <> 0, we supposed paths:  nnn/h.nnn'
       write(*,*)
       Stop ' '
      end if

! ... files:

       in=1 !input
       iout=2; Open(iout,file='H',form='unformatted') !output
  

      Call Read_iarg('ih1',ih1)
      if(ih1.le.0) Stop 'ih1 <= 0' 

      Call Read_iarg('ih2',ih2)
      if(ih2.le.0) Stop 'ih2 <= 0'

!----------------------------------------------------------------
! ... Finding LRANG2 in h.nnn files  
           
       Do ih = ih1,ih2
           write(AF,'(a,i3.3)') 'h.',ih
           if(Icheck_file(AF).eq.0) then 
               Write(errmsg, '(a,i3.3,a)') 'file h.', ih, ' not found.'
               Stop Trim(errmsg)
           end if 
           Open(in,file=AF,form='UNFORMATTED',status='OLD')
           read(in) NELC, NZ, LRANG2DUMMY, LAMAX, NAST, RA, BSTO
           
           If ( LRANG2 .lt. LRANG2DUMMY ) then
                LRANG2 = LRANG2DUMMY
           End if
          

        end do
!---------------------------------------------------------------


! ... Cycle over different h.nnn -files:

      istart= 1
Do ih = ih1,ih2

       write(AF,'(a,i3.3)') 'h.',ih
       if(folder.ne.0) write(AF,'(i3.3,a,i3.3)') ih,'/h.',ih
       if(Icheck_file(AF).eq.0) Cycle
       Open(in,file=AF,form='UNFORMATTED',status='OLD')

!----------------------------------------------------------------------
!                                                   target information:

       read(in) NELC, NZ, LRANG2_CURRENT, LAMAX, NAST, RA, BSTO

      if(istart.eq.1)  then ! Main if:  To guarantee that the first line is written once

        Allocate (COEFF(3,MRANG2), ENAT(NAST),LAT(NAST),ISAT(NAST), &
                  ENAT1(NAST), LAT1(NAST), ISAT1(NAST),NCONAT(NAST), &
                  IPAT(NAST), IPAT1(NAST) )

        read(in) ENAT
        read(in) LAT                                                                     
        read(in) ISAT,IPAT

        COEFF = 0.d0
        !LRANG2 = min(MRANG2,MaxLRANG2)
        read(in) ((COEFF(K,L),K=1,3),L=1,LRANG2_CURRENT) ! dummy Buttle corrections, they are not used (this line is for reading purposes )

        NELC1=NELC; NZ1=NZ; NAST1=NAST; ENAT1=ENAT; LAT1=LAT; ISAT1=ISAT; IPAT1=IPAT 
        LAMAX1=LAMAX; RA1=RA; BSTO1=BSTO

        write(iout) NELC, NZ, LRANG2, LAMAX, NAST, RA, BSTO ! General information 
 !       print*,  NELC, NZ, LRANG2, LAMAX, NAST, RA, BSTO 
        
        write(iout) (ENAT(i),i=1,NAST) ! target state energies     
!        print*, (ENAT(i),i=1,NAST)        

        write(iout) LAT ! L of atomic target states
!        print*, LAT        
 
        write(iout) ISAT !2S+1 of atomic target states
 !       print*, ISAT


        write(iout) (0.d0,0.d0,0.d0,i=1,LRANG2)! Buttle corrections (not Used in BSR), (definition as in subroutine H_OUT) 
!        print*, (0.d0,0.d0,0.d0,i=1,LRANG2)        

        istart = 0 ! prevents writing again the first 5 lines 
         
       else ! Checking that general information in h files is the same

        if(NAST1.ne.NAST) then
          write(*,*) ' NAST <> NAST1 for  ',AF; Stop
        end if
        if(NZ1.ne.NZ) then
          write(*,*) ' NZ <> NZ1 for  ',AF; Stop
        end if
        if(LAMAX1.ne.LAMAX) then
          write(*,*) ' LAMAX <> LAMAX1 for  ',AF; Stop
        end if
        if(RA1.ne.RA) then
          write(*,*) ' RA <> RA1 for  ',AF; Stop
        end if
        if(BSTO1.ne.BSTO) then
          write(*,*) ' BSTO <> BSTO1 for  ',AF; Stop
        end if     
    
        read(in) ENAT
        read(in) LAT
        read(in) ISAT,IPAT
        LM = min(MRANG2,LRANG2)
        read(in) ((COEFF(K,L),K=1,3),L=1,LRANGE2_CURRENT)

        Do N=1,NAST
         if(ENAT1(N).ne.ENAT(N)) then
          write(*,*) ' ENAT <> ENAT1 for ',AF; Stop
         end if
         if(LAT1(N).ne.LAT(N)) then
          write(*,*) ' LAT <> LAT1 for ',AF; Stop
         end if
         if(ISAT1(N).ne.ISAT(N)) then
          write(*,*) ' ISAT <> ISAT1 for ',AF; Stop
         end if
         if(IPAT1(N).ne.IPAT(N)) then
          write(*,*) ' IPAT <> IPAT1 for ',AF; Stop
         end if
        End do

 end if ! end of main if 

!----------------------------------------------------------------------
!                                               scattering information:

1 read(in) LRGL, NSPN, NPTY, NCHAN, MNP2, MORE0 

MORE=2; if(ih.eq.ih2.and.MORE0.eq.0) MORE=0 ! MORE = 2 in RMII, see ham.f90. 
!MORE=1 (original set-up)

Allocate (L2P(nchan), CF(nchan,nchan,lamax), VALUE(MNP2),&
                WMAT(NCHAN,MNP2),jkch(nchan))

      read(in) (NCONAT(N), N=1,NAST)
      read(in) (L2P(I), I=1,NCHAN),(jkch(i),i=1,nchan)
      read(in) (((CF(I,J,K), I=1,NCHAN), J=1,NCHAN), K=1,LAMAX)
      read(in) (VALUE(K),K=1,MNP2)
      read(in) ((WMAT(I,K),I=1,NCHAN), K=1,MNP2)

      write(iout) LRGL, NSPN, NPTY, NCHAN, MNP2, MORE ! Information about the partial wave  !MORE ISPOTENTIAL SOURCE OF ERROR
      print*, LRGL, NSPN, NPTY, NCHAN, MNP2, MORE 

      write(iout) (NCONAT(N), N=1,NAST)
 !     print*, (NCONAT(N), N=1,NAST)

      write(iout) (L2P(I), I=1,NCHAN) ! angular momentum
 !     print*, (L2P(I), I=1,NCHAN) 

      write(iout)  CF(1:NCHAN,1:NCHAN,1:LAMAX)! asymptotic coefficients  
!      print*, CF(1:NCHAN,1:NCHAN,1:LAMAX)
!       do lam = 1, 4
!        write(6,'(a,i5)') 'lam = ', lam 
!!        do i=1,nchan
 !        write(6,'(a,i5)') 'nchan = ', nchan 
 !        write(6,1091) CF(i,1:NCHAN,lam)
 !        write(6,*) 
 !       enddo
 !       write(6,*)
 !      enddo
!1091   format(1p100e11.3)       

!      write(iout) (VALUE(K),K=1,MNP2)! eigenvalues in descending order (RMII uses ascending) !POTENTIAL SOURCE OF ERROR
!      print*, (VALUE(K),K=1,MNP2)!
       
       write(iout) (VALUE(MNP2-K+1),K=1,MNP2) ! eigenvalues in ascending order
!       print*,  (VALUE(MNP2-K+1),K=1,MNP2) 


!      write(iout) ((WMAT(I,K),I=1,NCHAN), K=1,MNP2) ! surface amplitudes  (is this the correct order?)   !POTENTIAL SOURCE OF ERROR
!      print*, ((WMAT(I,K),I=1,NCHAN), K=1,MNP2) 

       write(iout)  ((WMAT(I,MNP2-K+1),I=1,NCHAN), K=1,MNP2) ! surface amplitudes
!       print*,  ((WMAT(I,MNP2-K+1),I=1,NCHAN), K=1,MNP2)

      Deallocate(L2P, CF, VALUE, WMAT, jkch) 

      istart = 0
      if(MORE0.ne.0) go to 1
      End do  ! over ih

      End  ! sum_hh
