!======================================================================
!     bsr_d - utility for calculation of dipole matrix elements
!              between R-matrix states saved in rsol.nnn files
!           - it also calculates the dipole matrix elements between
!             the target states 
!           - USE ONLY FOR LS COUPLING (Nov. 2024) 
!======================================================================
!     INPUT:    target, cfg.nnn, rsol.nnn, H (coming from bsr_H)
!
!     OUTPUT:   d, d00 (unformatted files)
!
!     SYSTEM CALLS:   MULT3,  BSR_DMAT3
!---------------------------------------------------------------------
!     this utilty is used as preparation for RMT calculations
!---------------------------------------------------------------------

      Use target
      Use channels

      Implicit real(8) (A-H,O-Z)

      Character(1)  :: blank = ' '
      Character(80) :: AS,BS,BI,BJ,AF,AFD,AFH
      Real(8), allocatable :: DL(:,:), DV(:,:) ! Dipole matrix elements for rsol.nnn 
      Real(8), allocatable :: TDL(:), TDV(:) ! Dipole matrix elements for target
      Real(8), allocatable :: CF(:,:,:) ! long-range coefficients

! ... short instructions:
      
       Call get_command_argument(1,AF)  
       if(AF.eq.'?') then
       write(*,'(/a)') 'bsr_d - utility for calculation of dipole matrix elements '   
       write(*,'(a)')  'between R-matrix states saved in rsol.nnn files'
       write(*,'(a)')  'In addition, it calculates the dipole elements between target states'
       write(*,'(/a)') 'INPUT:    target, cfg.nnn, rsol.nnn, H (coming from bsr_H)'
       write(*,'(/a)') 'OUTPUT:   dd_nnn_mmm, d00,  where nnn,mmm - indexes of partial waves'
       write(*,'(/a)') 'Call as:  bsr_dd  [nlsp=..]'
       write(*,'(/a)') 'where you may indicate the maximum partial waves "nlsp"'
       write(*,*)
       Stop 
      end if



      
!      klsp=0;
      
!      Call Read_iarg('nlsp',klsp)
!      if (klsp.le.0) stop ' nlsp<=0'
      nut=1;
      Open(nut,file='target',status='OLD')
      Call R_target (nut) ! reads target information, e.g., orbital angular momentum, spin, and parity
                          ! ltarg(i),istarg(i),iptarg(i) (i =1,ntarg) 
 

      Call R_channels(nut) ! channel(s) information, nlsp, ... 
      np=ntarg; Call Read_ipar(nut,'np',np) ! number of target states 
      ni=ntarg; Call Read_ipar(nut,'ni',ni) ! number of target states
  
      allocate(TDL(ntarg*ntarg),TDV(ntarg*ntarg)); TDL=0.d0; TDV=0.d0 
 
 
      kpol = 1; icount = 0

      ipri = 3; open(ipri,file='bsr_d.log')
      
      write(ipri,'(a)') 'B S R _ D'
      write(ipri,'(a)') '*********' 
      write(ipri,*)
      write(ipri,'(a)') repeat('=', 72)
      write(ipri,'(a)') 'Target'
      write(ipri,'(3(A6,1X))') 'NELC', 'NZ', 'NTARG'
      write(ipri,'(3(I6,1X))')  NELC, NZ, ntarg
      write(ipri,*)
      write(ipri,'(a)') 'Transitions according to s. rules'
      write(ipri,*)
      write(ipri,'(/2a5,5x,a/)') 'idx1', 'idx2',  'states involved'


!----------------------------------------------------------------------
!                      TARGET DIPOLE MATRIX ELEMENTS
!----------------------------------------------------------------------


! ... selection rules N-electron system following Oleg's fashion

        icount = 0        

        do i = 1, ntarg-1; do j = i+1, ntarg
           if(iptarg(i).eq.iptarg(j)) Cycle
           if(istarg(i).ne.istarg(j)) Cycle
           if(ITRA(ltarg(i),kpol,ltarg(j)).eq.0) Cycle ! ITRA(i1,i2,i3) checks triangle relations for momentums i1,i2,i3
           if(ltarg(i)+ltarg(j).lt.kpol) Cycle          
           icount = icount + 1
           write(ipri,'(2i5,5x,3a)') i, j, trim(BFT(i)), ' <-> ', trim(BFT(j))
         end do; end do


!----------------------------------------------------------------------




! ... read  assymtotic coefficients from H.DAT file

        in=2 ! unit for H file

        Open(in,file="H",status='OLD',form='UNFORMATTED')

        read(in) NELC, NZ, LRANG2, LAMAX, NAST, RA, BSTO ! general information
        read(in) (E,N=1,NAST) ! target energies, not used from H file
        read(in) (L,N=1,NAST) ! target orbital angular momentum, not used from H file
        read(in) (I,N=1,NAST) ! target spin, not used from H file
        read(in) ((C,K=1,3),L=1,LRANG2) ! Buttle corrections

      write(ipri,*)
      write(ipri,'(a,i0)') 'icount, total number of allowed transitions: ',icount
      write(ipri,'(a)') repeat('=', 72)



! ... cycle over partial waves

icount=0;
 1 Continue

      icount=icount+1 ! counting the number of partial waves 
     
      read(in,end=2) IL2, IS2, IP2, NCHAN, NHM, MORE
      IP2 = (-1)**IP2

      klsp = 0
       
      Do i=1,nlsp
       if(IL2.ne.Lpar(i)) Cycle 
       if(IS2.ne.ISpar(i)) Cycle
       if(IP2.ne.ipar(i)) Cycle
       klsp=i; Exit  ! Identifying number of partial wave in the H file
      End do
     
!     print*,klsp  
      
      if(klsp.eq.0) Stop 'unknown partial wave in H'
      if(nch(klsp).ne.nchan) Stop 'nchan <> nch' ! checking consistency between the number of channels: H file vs target

       read(in) (NCONAT, N=1,NAST)
       read(in) (L2P, I=1,NCHAN)
   
      if(allocated(CF)) Deallocate(CF); Allocate (CF(nchan,nchan,lamax)) ! array for long-range coefficients 
      read(in) (((CF(I,J,K), I=1,nchan), J=1,nchan), K=1,LAMAX) ! read coefficients
      read (in) (E,i=1,NHM)      ! skip eigenvalues, not needed
      read (in) ((W,i=1,nchan),j=1,NHM)   ! skip surface amplitudes, not needed

! ... calculation of the target dipole matrix elements

     eps = 1.D-6 ! tolerance 
       

!     print*, "Calling d_values for klsp = ", klsp 
      Call d_values(nut,klsp,np,ni,ipri,eps,lamax,nchan,CF,TDL,TDV)
!     print*,"end of call "     

      if(MORE.gt.0) go to 1 


 2 Continue

!      print*, icount

      if(nlsp.ne.icount) Stop 'did you miss some partial waves?' ! checking consistency

!   write(666,6666) MDLT(1:ntarg,1:ntarg)
!6666 FORMAT(10E15.5)

      write(666,1112) TDL
      write(666,1112) TDV

!stop 'here after reading'






!----------------------------------------------------------------------
!                          DIPOLE MATRIX ELEMENTS
!----------------------------------------------------------------------

      
      icount = 0 ! restarting counter

      write(ipri,*)
      write(ipri,'(a)') repeat('=', 72) 
      write(ipri,'(a)') 'N+1 electron system'
      write(ipri,*)
      write(ipri,'(a,i0)') 'nlsp (inast), number of partial waves: ', nlsp
      write(ipri,'(a)') repeat('=', 72)
      write(ipri,*)
      write(ipri,'(a/)') 'Transitions (L 2S+1 P) <-> (L 2S+1 P):'
      
     Do i=1,klsp ! LOOP OVER PARTIAL WAVE

       if(ipar(i).ne.1) Cycle

       write(BI,'(a,i3.3)') 'rsol.' ,i 
       if(Icheck_file(BI).eq.0) Cycle


       write(BI,'(a,i3.3)') 'cfg.' ,i 
   

 
       Do j=1,klsp   ! LOOP OVER PARTIAL WAVE

        write(BJ,'(a,i3.3)') 'rsol.' ,j 
        if(Icheck_file(BJ).eq.0) Cycle
        

        write(BJ,'(a,i3.3)') 'cfg.' ,j 
        
!--------------------------------------------------------------------------
!                       Selection Rules  N+1-electron system 
!--------------------------------------------------------------------------
             
        if(ipar(i).eq.ipar(j)) Cycle
        if(ispar(i).ne.ispar(j)) Cycle
        if(ITRA(lpar(i),kpol,lpar(j)).eq.0) Cycle
        if(lpar(i)+lpar(j).lt.kpol) Cycle   
      
        write(ipri, '(a,1x,a,1x,a)') TRIM(BI), '<->', TRIM(BJ) 
        write(ipri, '(a,1x,a,1x,a)') trim(AFP(i)), '<-> ',trim(AFP(j)) 
        write(ipri,'(i0,1x,i0,1x,i0,a,i0,1x,i0,1x,i0)') lpar(i),ispar(i),ipar(i),' <-> ',lpar(j),ispar(j),ipar(j)
        write(ipri,*) 

!--------------------------------------------------------------------------
!                        Angular Integrations: CALL MULT
!--------------------------------------------------------------------------

        write(AS,'(3(a,1x),a,2i3.3)') '../../BIN/mult3',trim(BJ),trim(BI),'E1 mult_bnk_',i,j
!       write(*,*) trim(AS)
        
        Call System(AS)

        write(AS,'(a,a,2i3.3,a)') 'cp ','mult_bnk_',i,j,' mult_bnk'
        Call System(AS)
   
!-------------------------------------------------------------------------


        icount = icount + 1     !  index for transition


!--------------------------------------------------------------------------
!                    Calculation of Dipole Matrix Elements      
!--------------------------------------------------------------------------

        write(AF,'(a,i3.3)') 'D',icount
        write(AS,'(a,a,a,a,a,a,a)')  &
        '../../BIN/bsr_dmat3 ',trim(BJ),' ',trim(BI),' b d  AF_dd=',trim(AF)
!       write(*,*) trim(AS)
        Call System(AS)



      write(81,*)  icount, ispar(j),lpar(j),(1-ipar(j))/2, &
                             ispar(i),lpar(i),(1-ipar(i))/2
      End do;
       
      End do

!--------------------------------------------------------------------------
!  Getting relevant parameters (noterm, kstate1, kstate2) for allocation    
!--------------------------------------------------------------------------    
      
       idout = 5;  open(idout,file='d',form='UNFORMATTED',STATUS='unknown')    
      
        write(idout) TDL ! target dipoles in lenght form
        write(idout) TDV  ! target dipoles in velocity form   

       Do ic = 1, icount
          write(AFD,'(a,i3.3)') 'D' ,ic
         idinp=12;  open(idinp,file=AFD,form='UNFORMATTED', STATUS='OLD')
          read(idinp) noterm,nstate2,kch2,JLT2,JPT2,JST2,nstate1,kch1,JLT1
!          print*, noterm(ic),nstate2(ic),nstate1(ic)
        
          Allocate(DL(nstate2,nstate1),DV(nstate2,nstate1))
           
      n2 = 0

      Do jk = 1, 1 + (nstate1-1)/noterm
         n1=n2+1; n2 = min(n2+noterm,nstate1)
         m2 = 0
         Do ik = 1, 1 + (nstate2-1)/noterm
            m1=m2+1; m2 = min(m2+noterm,nstate2)
      read(idinp) DL(m1:m2,n1:n2), DV(m1:m2, n1:n2)   
      End do
      End do
      

!--------------------------------------------------------------------------------------
!                                  d file ala RMATRIX II
!--------------------------------------------------------------------------------------

      ! FORMATTED FILES

!      write(666,*) noterm, nstate2, kch2, nstate1, kch1
!      write(666,1112) DL
!      write(666,1112) DV  
       
      write(idout) noterm, nstate2, kch2, nstate1, kch1
      write(idout) DL ! RMII: dip1
      write(idout) DV ! RMII: dip2

 
1112 FORMAT(10F8.3)


       Deallocate(DL)
       Deallocate(DV)
       close(idinp)

      End Do
      
      close(idout)
      
      write(ipri,*) 
      write(ipri,'(a)') repeat('=', 72)
      write(ipri,'(a,i0)') 'icount, total  number of allowed transitions: ',icount
      write(ipri,'(a)') repeat('=', 72)   




!--------------------------------------------------------------------------------------
!                                  D00 file
!--------------------------------------------------------------------------------------

      open(82,file='d00',form='UNFORMATTED')
      open(123,file='D00',form='FORMATTED')

      write(82) icount
      write(123,*) icount
      rewind(81)
    
      Do i = 1,icount
        read(81,*) i0,i1,i2,i3,i4,i5,i6
        write(82)     i1,i2,i3,i4,i5,i6
        write(123,*)     i1,i2,i3,i4,i5,i6
      End do
      
      Close(82)
      Close(81)
     
       Close(81,status='DELETE') 

      



      close(ipri)
      End ! program bsr_dd




!==========================================================================
      Subroutine d_values(nut,klsp,np,ni,ipri,eps_acf,km,nchan,ACF,TDL,TDV)
!==========================================================================
!     define the dipole matrix elements  between target states based on the
!     given asimptotic coefficients ACF for k=1 
!-------------------------------------------------------------------------
      Use zconst, only: c_au, time_au
      Use target
      Use channel

      Implicit none
      Real(8) :: ACF(nchan,nchan,km)
      Real(8) :: AK(ntarg,ntarg),  de(ntarg,ntarg)
      Real(8), intent(out) ::  TDL(ntarg), TDV(ntarg)
      Real(8) :: DLT(ntarg,ntarg), DVT(ntarg,ntarg)
      Real(8) :: RDME(ntarg,ntarg) 
      Integer :: IP(ntarg,ntarg)
      Real(8), external :: Z_6jj, Reduce_factor
      Real(8) :: S,SS, g1,g2, eps_acf,a,f
      Integer :: i,j, i1,i2,nt, ipri, km, nchan, nut, klsp, k,ij, np,ni
      Character(100) :: line
      
      Call R_channel(nut,klsp)

       AK=0.d0; IP=0; DLT=0.d0; DVT=0.d0
       Do i=1,nch-1; Do j=i+1,nch
          S=ACF(i,j,1)/2.d0;
          i1=iptar(i); i2=iptar(j)
           if(abs(S).lt.eps_acf) Cycle
           if(i1.gt.ni.and.i2.gt.np) Cycle

       SS = Reduce_factor(i,j,1)

       if(abs(SS).lt.eps_acf) Cycle
       S = S / SS
       DLT(i1,i2) = DLT(i1,i2) + S
 
       de(i2,i1)=Etarg(i2)-Etarg(i1)

!    if(istarg(i1).ne.0) then
!        S = S*S*istarg(i1)
!        g1 = (2*ltarg(i1)+1) * istarg(i1)
!        g2 = (2*ltarg(i2)+1) * istarg(i2)
!       else
!        S = S*S
!        g1 = jtarg(i1)
!        g2 = jtarg(i2)
!       end if

!       f = 2.d0/3.d0*de(i2,i1)*S /g1
!       a = 4.d0/3.d0*de(i2,i1)**3*S/c_au**3/time_au /g2

       DLT(i2,i1) = DLT(i2,i1) + S

       IP(i1,i2) = IP(i1,i2) + 1  ! this is a factor that corrects overcounting, Oleg's heritage
      End do; End do


      nt = 0
      Do i1=1,ntarg; Do i2=i1,ntarg
       if(IP(i1,i2).eq.0) Cycle
       DLT(i1,i2) = DLT(i1,i2) / IP(i1,i2)
       DLT(i2,i1) = DLT(i2,i1) / IP(i1,i2)
!       AK(i1,i2) = AK(i1,i2) / IP(i1,i2)
       nt = nt + 1
       DVT(i1,i2) = DLT(i1,i2)*de(i2,i1)
      End do; End do


! ... print results in log file:

      if(nt.gt.0) &
      write(ipri,'(/a,i0)') 'Target dipole elements found in partial wave no. ', klsp 
      write(ipri,'(/2a5,2a15,5x,a/)') 'idx1', 'idx2', 'length form', 'velocity form', 'states involved'

      Do i=1,ntarg-1; Do j=i+1,ntarg
       if(DLT(i,j).eq.0) cycle  
      
!      if(AK(i,j).eq.0.d0) then 
!      print*, 'I am not considering',  trim(BFT(i)), ' <-> ', trim(BFT(j))  
!      Cycle
!      end if 

        write(ipri,'(2i5,E15.5,E15.5 ,5x,3a)') &
        i,j, DLT(i,j),DVT(i,j), trim(BFT(i)), ' <-> ', trim(BFT(j))

        TDL((j-1)*ntarg+i)=DLT(i,j);  TDV((j-1)*ntarg+i)=DVT(i,j)      
        TDL((i-1)*ntarg+j)=DLT(i,j);  TDV((i-1)*ntarg+j)=-DVT(i,j) 
     
      End do; End do

        write(ipri,*)
        write(ipri,*)'------------------------------------------------------------------------'

      End  Subroutine d_values


!======================================================================
      Real(8) Function Reduce_factor(ich,jch,k)
!======================================================================
!     define factor connecting reduced dipole matrix element with the
!     asymptotic coefficient for multipole index k
!----------------------------------------------------------------------
      Use target
      Use channel

      Implicit none
      Integer, intent(in) :: ich,jch,k
      Integer :: ll1,ll2,jj1,jj2,it,jt,L1,L2,S1,S2,J1,J2,LT,JJ, ip,jp,kz
      Real(8) :: S,SS, zero = 0.d0
      Real(8), external :: ZCLKL, Z_6jj, Z_6j

      Reduce_factor = zero

      ll1 = lch(ich);  ll2 = lch(jch)
      S = ZCLKL(ll1,k,ll2)
      if(S.eq.zero) Return

      jj1 = jkch(ich); jj2 = jkch(jch)
      it = iptar(ich); jt = iptar(jch)
      L1 = ltarg(it);  L2 = ltarg(jt)
      S1 = istarg(it); S2 = istarg(jt)
      J1 = jtarg(it);  J2 = jtarg(jt)
      ip = iptarg(it); jp = iptarg(jt)
      kz = 0; if(ip.ne.jp) kz=1;  kz = kz - k
      if(MOD(kz,2).ne.0) then; Reduce_factor = zero; Return; end if

      LT = lpar; JJ =jpar
      if(coupling.eq.'LS') then

       if(S1.ne.S2) then; Reduce_factor = zero; Return; end if
       S = S * Z_6jj(L1,ll1,LT,ll2,L2,k)
       kz = L2+ll1+LT
       S = S * (-1) ** kz

      elseif(coupling.eq.'JK') then

       if(jj1.ne.jj2) Return
       S = S * Z_6j(ll2+ll2+1,ll1+ll1+1,k+k+1,J1,J2,jj1)
       kz = J2+JJ+JJ-ll1-ll1-jj1+1; kz=kz/2
       S = S * (-1) ** kz

      elseif(coupling.eq.'JJ') then

       SS = jj1*jj2;   S = S * sqrt(SS)
       S = S * Z_6j(ll1+ll1+1,2,jj1,jj2,k+k+1,ll2+ll2+1)
       S = S * Z_6j(J1,jj1,JJ,jj2,J2,k+k+1)
       kz = J2+jj1+JJ+ll1+ll1+1+jj2+k+k; kz = kz/2
       S = S * (-1) ** kz

      end if

      Reduce_factor = S

      End Function Reduce_factor




