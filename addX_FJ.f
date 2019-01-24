!-------- Calculate cross-sections to use in the Fast-J photolysis scheme --------
!
! Louis Marelle, 2018/10/24, based on mprather's routine for Fast-JX
!
!----addX_FJX.f  -- with user supplied subroutine that supplies X-section x Q-yield
!----     generate fast-JX 18-bin X-sections
!---------revised and updated for (mprather,2/2013)
      implicit none
      ! Number of interpolation wavelengths in wavel-bins.dat (can be
      ! larger than NB in wavel-bins.dat)
      integer, parameter :: NB_ = 100
      ! Number of solar wavelength in solar-p05nm-UCI.dat
      integer, parameter :: NS_ = 40000
      ! Number of wavelengths in XO3-p05nm-UCI.dat
      integer, parameter :: NZ_ = 13550
      ! Number of Fast-J wavelengths (18 for Fast-JX, but the 7 longer
      ! wavelengths are Fast-J wavelengths)
      integer, parameter :: NJ_ = 18

      real*8 SRB(15,NJ_)
      real*8, dimension(NB_+1) :: WBIN
      real*8, dimension(NB_) :: FBIN, ABIN
      real*8, dimension(NJ_) :: FFBIN,AABIN
      integer IJX(NB_), ITT
      integer NB,I,J,J1,J2,K,K1,K2
      integer INIT
      real*8 W(NS_),F(NS_)
      integer IBINJ(NS_)
      real*8 WZ(NZ_),X(NZ_,3)
      real*8 W1,W2, TT,XP,XM, WW,XNEW
      character*6 TITLNEW
      character*4 TITLET

      open (1, file='wavel-bins.dat', status='OLD')
        SRB(:,:) = 0.d0
        read(1,'(i5)') NB
        if (NB .gt. NB_) stop
        read(1,'(5x,f8.3)') (WBIN(I), I=1,NB+1)
        read(1,*)
        read(1,*)
        read(1,'(2x,15f5.1)') ((SRB(I,J),I=1,15),J=1,8)
        read(1,*)
        read(1,'(5x,i5)') (IJX(I),I=16,NB)
      close (1)

      open (2, file='solar-p05nm-UCI.dat', status='OLD')
        read(2,*)
        read(2,*)
        read(2,'(f10.4,e10.3)') (W(J),F(J), J=1,NS_)
      close (2)

      open (3, file='XO3-p05nm-UCI.dat', status='OLD')
        read(3,*)
        read(3,*)
        read(3,'(f10.4,3e10.3)') (WZ(J),X(J,1),X(J,2),X(J,3), J=1,NZ_)
      close (3)

!!!!!!!!!!!!!!!!!!!!initialization call to user subroutine!!!!!!!!!!!!!!
        INIT = 0
      call XBRONO2 (WW,TT,XP,XM, X,INIT,TITLNEW)



!---synchronize with the BRONO2 cross sections (whether done or not)
!---will loop K=1:NZ_  (J = K - 1 + K1), wavel = WZ(K)
        do J=1,NS_
          if (WZ(1) .eq. W(J)) goto 10
        enddo
          write(6,*) ' cannot synch the solar & Xsections'
          stop
   10   K1 = J
        K2 = min( NS_, NZ_+K1-1)
!          write(6,'(a,2f10.3,i10)') ' synch:',WZ(1),W(K1),K1
!          write(6,'(a,2f10.3,i10)') ' synch:',WZ(NZ_),W(K2),K2


!---now assign bin #(I=1:77) to each p05nm microbin J (1:40000)
        IBINJ(:) = 0
      do I=1,NB
         W1 = WBIN(I)
         W2 = WBIN(I+1)
        do J=1,NS_
          if (W(J) .gt. W1) goto 11
        enddo
          J = NS_ + 1
   11     J1 = J
       do J=J1,NS_
          if (W(J) .gt. W2) goto 12
        enddo
          J = NS_ + 1
   12     J2 = J-1
        do J=J1,J2
          IBINJ(J) = I
        enddo
!          write(6,'(i5,2f9.3,2i9,2f9.3)') I, W1,W2, J1,J2,W(J1),W(J2)
      enddo


!>>>>be aware that this binning does not interpolate and is OK for large bins
!          it has 7% error in the very short wavel S-R bins of pratmo.
!>>>>it should be fine for weighting cross sections!

!  Total Q-yld from Stern-Volmer converted (as in acetone) to 3 different trop levels:
!         alt=   0 km      5 km      13 km
!           M=  2.46E+19   1.50E+19   5.8E+18
!           T=     295K      272K       220K
!           P =    999 hPa   566 hPa    177 hPa
!   Only apply a single overall quantum yield PHI
!        PHI = exp[ -0.055 * (w-308) ] /( 5.5 + 9.2e-19*[M])
!   Thus there are 3 tables for MeVK, each designated by T(K)


!---this looping is set for Temperature only
        XP = 999.d0
        XM = 2.46d19

!!!!!!!!!!!!!!!!!!!major temperature-density loop !!!!!!!!!!!!!!!!!!!!!!
      do ITT =2,2

        if (ITT.eq.1) then
          TITLET = ' 220'
          TT = 220.d0
        else if (ITT.eq.2) then
          TITLET = ' 230'
          TT = 230.d0
        else
          stop
        endif
!---now ready to do any flux-weighted means over the 77-pratmo bins
         FBIN(:) = 0.d0
         ABIN(:) = 0.d0

!!!!!!!!!!!!!!!!!!!primary high-resolution wavelength loop!!!!!!!!!!!!!!
       do J = K1,K2
        K = J - K1 + 1
        I = IBINJ(J)
        if (I .gt. 0) then

        call XBRONO2 (W(J), TT,XP,XM, XNEW, INIT, TITLNEW)

          FBIN(I) = FBIN(I) + F(J)
          ABIN(I) = ABIN(I) + F(J)*XNEW
        endif
       enddo
       do I=1,NB
        if (FBIN(I) .gt. 0.d0) ABIN(I) = ABIN(I)/FBIN(I)
       enddo
!---write out UCI std 77-bin data
!       write(6,'(a6,a4/(1p,8e10.3))') TITLNEW,TITLET, ABIN

!!!!!!!!!!!!!!!!!!!secondary sum 77-bin pratmo ==> 18-bin fast-JX!!!!!!!
!---combine fast-JX bins:
!---    non-SR bands (16:NB) are assigned a single JX bin
!---    SR bands are split (by Opacity Distrib Fn) into a range of JX bins
        FFBIN(:) = 0.d0
        AABIN(:) = 0.d0
      do I=16,NB
        J = IJX(I)
        FFBIN(J) = FFBIN(J) + FBIN(I)
        AABIN(J) = AABIN(J) + FBIN(I)*ABIN(I)
      enddo
      do I=1,15
        do J=1,NJ_
          FFBIN(J) = FFBIN(J) + FBIN(I)*SRB(I,J)
          AABIN(J) = AABIN(J) + FBIN(I)*ABIN(I)*SRB(I,J)
        enddo
      enddo
      do J=1,NJ_
        if (FFBIN(J) .gt. 0.d0) AABIN(J) = AABIN(J)/FFBIN(J)
      enddo


!!!!!!!!!!!!!!!!!!!save UCI fast-JX data bins!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(6,'(a6,a4,1p,7e10.3/10x)')
     &     TITLNEW, TITLET, AABIN(NJ_-6:NJ_)

      enddo
      stop
      end


!-------------sample subroutine for fast-JX Xsection generation---------
!-----------------------------------------------------------------------
      subroutine XBRONO2 (WW,TT,PP,MM, XXWT,INIT, TITLNEW)
!-----------------------------------------------------------------------
!   WW = wavelength (nm)
!   TT = temerature (K) for interpolation
!   PP = pressure (hPa) can be used in Stern-Volmer formula if need be
!   MM = air density (#/cm3), ditto
!   XXWT = cross section (cm2) as a function of WW and TT (and PP, MM)
!   INIT = initialization:
!     if INIT.eq.0 .then reads in any tables and sets Xsect name to TITLNEW

      implicit none
      real*8, intent(in) :: WW, TT, PP, MM
      integer, intent(inout) :: INIT
      real*8, intent(out) :: XXWT
      character*6, intent(out) :: TITLNEW
      character*80 FTBL,TABLE,FORMW
      real*8 W(999), XW(999),WWL,TTL,XBFACT,XTFACT,XXW,FW
      real*8 TFACT,T3(3),XT3(999,3)
      integer NW,NB,N,I,IW,IT

      if(INIT .eq. 0) then
        ! Read cross section file
        TITLNEW = 'BRONO2 '
        FTBL = 'XBRONO2_230K_IUPACtbl.dat'
        open (3, file=FTBL, status='OLD')
          ! Read header and write it
          read(3,'(a)') TABLE
          write(6,'(2a/a)') ' openfile=',FTBL, TABLE
          ! Read number of lines/wavelengths and format
          read(3,'(i4,1x,a)') NW,FORMW
        ! Read wavelength and cross section
        do N=1,NW
         read(3,FORMW) W(N),XW(N)
        enddo
        close(3)
        INIT = 1
      else
        TTL = min(T3(3), max(T3(1), TT))
        if (TT .gt. T3(2)) then
          IT = 2
        else
          IT = 1
        endif
        TFACT = (TTL - T3(IT))/(T3(IT+1) - T3(IT))
!---interpolate X-section vs. Wavelength, T=298
        IW = 1
        do I = 2,NW-1
          if (WW .gt. W(I)) IW = I
        enddo
        FW = (WW - W(IW))/(W(IW+1) - W(IW))
        FW = min(1.d0, max(0.d0, FW))
        XXW = XW(IW) + FW*(XW(IW+1)-XW(IW))
        XXWT = XXW  * 1.d-20
        ! If the wavelength is outside the wavelength range force to 0
        if (WW .gt. W(NW)) XXWT = 0.0
        ! write(*, '(4e13.4,1x,2i4)') WW,XXWT,W(IW),W(IW+1),IW,NW
      endif
      return
      end
