C **********************************************************************
C "EL2000C.F" - PROGRAM COMPUTING AN ORBIT OF A METEOR, WHICH GEOCENTRIC
C               RADIANT AND VELOCITY ARE KNOWN
C               VERSION: January 20, 2009
C             - TO EQUINOX 2000.0
C **********************************************************************
C
C MAIN PROGRAM
C
C   INPUT: all data on meteors in the standardized format
C
C   OUTPUT:
C      EA(1) - PERIHELION DISTANCE [AU]
C      EA(2) - ECCENTRICITY
C      EA(3) - ARGUMENT OF PERIHELION [DEGREES]
C      EA(4) - LONGITUDE OF NODE [DEGREES]
C      EA(5) - INCLINATION TO ECLIPTIC [DEGREES]
C      EA(6) - (empty)
C
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      IMPLICIT REAL*8(A-H,K-Z)
      IMPLICIT INTEGER(I,J)
      CHARACTER*4 MON
      character*24 filin,fil15,fil16
      character*72 ret72
      DIMENSION EA(6),MON(12),EB(6),R(6),WEL(7),iau(5000),isol(5000)
      dimension sollam(5000),al(5000),dl(5000),vg(5000)
      dimension q(5000),e(5000),som(5000),gom(5000),sk(5000),jnn(5000)

      PI=4.0D0*DATAN(1.0D0)
      PI180=PI/0.180D3
      GAUSS=1.7202098950D-2
      G2=GAUSS*GAUSS
      GMS=0.2959122082855911E-03
      GMB=0.8997011346712499E-09
      AU=0.1495978706910000E+09  
      AUKM=0.149597870D9
      D24=8.6400D4
      VELU=D24/AUKM
C      EPS=2.344579D1*PI180 ... 1950.0
      EPS=2.3439291D1*PI180
      CEPS=DCOS(EPS)
      SEPS=DSIN(EPS)
C
      MON(1)='JAN.'
      MON(2)='FEB.'
      MON(3)='MAR.'
      MON(4)='APR.'
      MON(5)='MAY '
      MON(6)='JUNE'
      MON(7)='JULY'
      MON(8)='AUG.'
      MON(9)='SEP.'
      MON(10)='OCT.'
      MON(11)='NOV.'
      MON(12)='DEC.'
C
c Tolerance limits:
      open(unit=9,file='inparams.ele',access='sequential')
        read(9,*) ret72
        read(9,*) filin
        read(9,*) ret72
        read(9,*) fil15
        read(9,*) ret72
        read(9,*) fil16
        read(9,*) ret72
        read(9,*) DELQ
        read(9,*) ret72
        read(9,*) DELE
        read(9,*) ret72
        read(9,*) DELSO
        read(9,*) ret72
        read(9,*) DELGO
        read(9,*) ret72
        read(9,*) DELSK
      close(unit=9)

      open(unit=10,file=filin,access='sequential')
      j=0
  10  continue
      j=j+1
      read(10,*,end=20) iau(j),isol(j),sollam(j),al(j),dl(j),vg(j),q(j),
     *e(j),som(j),gom(j),sk(j),jnn(j)
      goto 10
  20  continue
      jall=j-1
      close(unit=10)

      OPEN(UNIT=15,FILE=fil15,ACCESS='SEQUENTIAL')
      OPEN(UNIT=16,FILE=fil16,ACCESS='SEQUENTIAL')
      WRITE(15,500)
      WRITE(15,510)
C
      DO 100 i=1,jall
        zemlam=sollam(i)+1.80d2
        if(zemlam.ge.3.60d2) zemlam=zemlam-3.60d2
        ttv=(zemlam*pi180-1.75347031435d0)/6.28307584918d3
        t20=ttv*1.0d1
        call EARTHrv(t20,X,Y,Z,VXE,VYE,VZE,WEL)
c        ET=3.65250d5*ttv+2.4515450d6
c        J=3
c        CALL  Pleph ( ET, J, NCTR, R )
C
cC The rectangular ecliptical coordinates of the Earth/meteor:
c        X=R(1)
c        Y=R(2)*CEPS+R(3)*SEPS
c        Z=-R(2)*SEPS+R(3)*CEPS
        RE=DSQRT(X*X+Y*Y+Z*Z)
cC The rectangular ecliptical components of the Earth's velocity vector:
c        VXE=R(4)
c        VYE=R(5)*CEPS+R(6)*SEPS
c        VZE=-R(5)*SEPS+R(6)*CEPS
cC
        VXR=VXE
        VYR=VYE*CEPS-VZE*SEPS
        VZR=VYE*SEPS+VZE*CEPS
        SDZ=VZR
        CDZ=DSQRT(1.0D0-SDZ*SDZ)
        ALZ=DATAN(VYR/VXR)
        ALZ=DABS(ALZ)
        IF(VXR.LT.0.0D0.AND.VYR.GE.0.0D0) ALZ=PI-ALZ
        IF(VXR.LT.0.0D0.AND.VYR.LT.0.0D0) ALZ=PI+ALZ
        IF(VXR.GE.0.0D0.AND.VYR.LT.0.0D0) ALZ=2.0D0*PI-ALZ
C
        YMIN=9.999D99
        YMAX=0.0D0
C
        CALA=DCOS(al(i)*PI180)
        SALA=DSIN(al(i)*PI180)
        CDLA=DCOS(dl(i)*PI180)
        SDLA=DSIN(dl(i)*PI180)
        VXG=vg(i)*VELU*CALA*CDLA
        VYR=vg(i)*VELU*SALA*CDLA
        VZR=vg(i)*VELU*SDLA
        VYG=VYR*CEPS+VZR*SEPS
        VZG=-VYR*SEPS+VZR*CEPS
        VX=-(VXG-VXE)
        VY=-(VYG-VYE)
        VZ=-(VZG-VZE)
C
        VVH=DSQRT(VX*VX+VY*VY+VZ*VZ)/VELU
C
        CALL ELEMS(t20,X,Y,Z,VX,VY,VZ,EB)
C
CC      CGAM=SDZ*SDLA+CDZ*CDLA*DCOS(ALZ-YY(I,11)*PI180)
CC      VE2L=8.8804D2*(2.0D0/RE-1.0D0)
CC      VEL=DSQRT(VE2L)
CC      CGAM=(VE2L+YY(I,14)*YY(I,14)-YY(I,15)*YY(I,15))/2.0D0/VEL/YY(I,14)
CC      SGAM=DSQRT(1.0D0-CGAM*CGAM)
CC      GAMMA=DATAN(SGAM/CGAM)
CC      IF(GAMMA.LT.0.0D0) GAMMA=PI+GAMMA
CC      GAMMA=GAMMA/PI180
C
 180  CONTINUE
        QW=EB(1)
        AW1=(1.0D0-EB(2))/EB(1)
        AW=1.0D0/AW1
        QQW=AW*(1.0D0+EB(2))
        EW=EB(2)
        SKW=EB(5)
        SOW=EB(3)
        GOW=EB(4)
        PIW=EB(3)+EB(4)
        IF(PIW.GE.3.60D2) PIW=PIW-3.60D2
CC      GAMW=GAMMA
        VINFW=DSQRT(vg(i)*vg(i)+1.249924D2)
        XX=VINFW
        IF(XX.LT.YMIN) YMIN=XX
        IF(XX.GT.YMAX) YMAX=XX
C      QS=QW
C      AS=AW
C      ES=EW
C      SKS=SKW
C      SOS=SOW
C      GOS=GOW
C      PIS=PIW
C      LAMS=GAMW
C      MNS=MNO
C      VVHS=VVH
C      VIS=VINFW
C  22  CONTINUE
C  44  CONTINUE
C  66  CONTINUE
      WRITE(15,610) iau(i),isol(i),q(i),e(i),som(i),gom(i),sk(i)
      WRITE(15,650) QW,EW,SOW,GOW,SKW
C
CC      IF(FQ(I,6).EQ.' 1') THEN
CC      IF(FE(I,6).EQ.' 1') THEN
CC      DELX=DYY(I,6)
CC      ELSE
CC      DELX=1.0D-2*YY(I,6)
CC      END IF
CC      IF(DABS(YY(I,6)-GAMW).GT.DELX) WRITE(16,*) CID(I),' [',ANO(I),']: 
CC     *  LS_obs,   LS_calc = ',YY(I,6),GAMW
CCCC      ELSE
CCCC      WRITE(16,*) CID(I),' [',ANO(I),']:  impossible checking of q; obse
CCCC     *rved value is missing .'
CC      END IF
      IF(DABS(q(i)-QW).GT.DELQ) WRITE(16,*) iau(i),isol(i),'  q_obs,   q
     *_calc = ',q(i),QW
      IF(DABS(e(i)-EW).GT.DELE) WRITE(16,*) iau(i),isol(i),'  e_obs,   e
     *_calc = ',e(i),EW
      IF(DABS(sk(i)-SKW).GT.DELSK) WRITE(16,*) iau(i),isol(i),'  i_obs, 
     *  i_calc = ',sk(i),SKW
      IF(DABS(som(i)-SOW).GT.DELSO) WRITE(16,*) iau(i),isol(i),'  arg_ob
     *s, arg_calc = ',som(i),SOW
      IF(DABS(gom(i)-GOW).GT.DELGO) WRITE(16,*) iau(i),isol(i),'  nod_ob
     *s, nod_calc = ',gom(i),GOW
 100  CONTINUE
      CLOSE(UNIT=10)
      CLOSE(UNIT=15)
      CLOSE(UNIT=16)
C
      STOP
C
 500  FORMAT(' IAU sol.  q        e       arg.     node       i')
 510  FORMAT('----------------------------------------------------')
 600  FORMAT(A8,1X,A6,2X,I4,I3,2F9.5,F8.4,F9.5,2F9.3,F11.5)
 610  FORMAT(i4,i3,2f9.5,3f9.3)
 620  FORMAT(A8,1X,A6,2X,I4,I3,2F9.5,8X,F9.5,2F9.3,F11.5)
 630  FORMAT(A8,1X,A6,2X,I4,I3,F9.5,9X,F8.4,F9.5,2F9.3,F11.5)
 632  FORMAT(A8,1X,A6,18X,F9.5,F8.4,F9.5,2F9.3,F11.5,/)
 650  FORMAT(7x,2f9.5,3f9.3,/)
C
      END 
C ---------------------------------------------------------------------
      SUBROUTINE ELEMS(T20,X,Y,Z,VX,VY,VZ,EA)
C ---------------------------------------------------------------------
C      - CALCULATION OF ORBIT OF A BODY, WHICH POSITION IN A GIVEN TIME
C        IS CHARACTERIZED WITH RADIUS AND VELOCITY VECTORS
C
C   INPUT QUANTITIES:
C        T20 - TIME, WHEN BODY'S POSITION IS KNOWN [JULIAN CENTURIES
C              FROM 2000.0; JD = 2451545.0]
C        X, Y, Z - RECTANGULAR ECLIPTIC HELIOCENTRIC COORDINATES OF
C                  THE BODY IN TIME T20 [AU]
C        VX, VY, VZ - COMPONENTS OF BODY'S VELOCITY VECTOR IN THE EC-
C                     LIPTIC HELIOCENTRIC COORDINATE SYSTEM IN TIME T20
C                     [AU/DAY]
C
C   OUTPUT QUAMTITIES (ELEMENTS OF BODY'S ORBIT):
C        EA(1) - PERIHELION DISTANCE [AU]
C        EA(2) - ECCENTRICITY
C        EA(3) - ARGUMENT OF PERIHELION [DEGREES]
C        EA(4) - LONGITUDE OF NODE [DEGREES]
C        EA(5) - INCLINATION TO ECLIPTIC [DEGREES]
C        EA(6) - TIME OF PERIHELION PASSAGE [JULIAN CENTURIES FROM
C                2000.0; JD = 2451545.0]
C
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION EA(6)
C
      PI=4.0D0*DATAN(1.0D0)
      PI05=PI/2.0D0
      PI2=2.0D0*PI
      PI180=PI/0.180D3
      GAUSS=0.17202098950D-1
      G2=GAUSS*GAUSS
C
      R=DSQRT(X*X+Y*Y+Z*Z)
      V=DSQRT(VX*VX+VY*VY+VZ*VZ)
      CVFI=-(X*VX+Y*VY+Z*VZ)/R/V
      SVFI2=1.0D0-CVFI*CVFI
      YMEN=2.0D0/R-V*V/G2
      IF(DABS(YMEN).GT.1.0D-12) GOTO 10
      EA(1)=R*SVFI2
      EA(2)=1.0D0
      GOTO 20
  10  CONTINUE
      A=1.0D0/YMEN
      EA(2)=DSQRT(1.0D0-R*(2.0D0-R/A)*SVFI2/A)
      EA(1)=A*(1.0D0-EA(2))
C
  20  CONTINUE
      IF(EA(2).LT.1.0D-12) GOTO 40
      CF0=(EA(1)*(1.0D0+EA(2))/R-1.0D0)/EA(2)
      IF(DABS(CF0).GT.1.000001D0) WRITE(5,*) 'ERROR IN ELEMS; CF0'
      IF(DABS(CF0).GT.1.000001D0) STOP
      IF(CF0.GT.1.0D0) CF0=0.999999999D0
      IF(CF0.LT.-1.0D0) CF0=-0.999999999D0
      SF0=DSQRT(1.0D0-CF0*CF0)
      IF(DABS(CF0).LT.1.0D-12) GOTO 30
      F0=DATAN(SF0/CF0)
      IF(F0.LT.0.0D0) F0=PI+F0
      F=F0
      IF(CVFI.GT.0.0D0) F=PI2-F0
      GOTO 40
  30  F=PI/2.0D0
C
  40  CONTINUE
      BX=Y*VZ-Z*VY
      BY=Z*VX-X*VZ
      BZ=X*VY-Y*VX
      CI=BZ/DSQRT(BX*BX+BY*BY+BZ*BZ)
      IF(DABS(CI).GT.1.000001D0) WRITE(5,*) ' ERROR IN SUBROUTINE ELEMS 
     *(CALCULATION OF "CI")'
      IF(DABS(CI).GT.1.000001D0) STOP
      IF(CI.GT.1.0D0) CI=0.999999999D0
      IF(CI.LT.-1.0D0) CI=-0.999999999D0
      SI=DSQRT(1.0D0-CI*CI)
      IF(DABS(CI).LT.1.0D-12) GOTO 42
      YI=DATAN(SI/CI)
      IF(YI.LT.0.0D0) YI=PI+YI
      EA(5)=YI/PI180
      GOTO 45
  42  EA(5)=9.0D1
C
  45  CONTINUE
      CGO=-BY/DSQRT(BX*BX+BY*BY)
      SGO=BX/DSQRT(BX*BX+BY*BY)
      IF(DABS(CGO).LT.1.0D-12) GOTO 47
      YGO=DATAN(SGO/CGO)
      YGO=DABS(YGO)
      GOTO 50
  47  YGO=PI/2.0D0
  50  CONTINUE
      IF(CGO.GE.0.0D0.AND.SGO.LT.0.0D0) YGO=PI2-YGO
      IF(CGO.LT.0.0D0.AND.SGO.GE.0.0D0) YGO=PI-YGO
      IF(CGO.LT.0.0D0.AND.SGO.LT.0.0D0) YGO=PI+YGO
      EA(4)=YGO/PI180
C
      IF(EA(2).LT.1.0D-12) GOTO 60
      CSOF=(X*CGO+Y*SGO)/R
      IF(DABS(CI).LT.1.0D-12) GOTO 52
      SSOF=(Y*CGO-X*SGO)/R/CI
      GOTO 54
  52  SSOF=Z/R/SI
  54  CONTINUE
      IF(DABS(CSOF).LT.1.0D-12) GOTO 56
      SOF=DATAN(SSOF/CSOF)
      SOF=DABS(SOF)
      GOTO 58
  56  SOF=PI/2.0D0
  58  CONTINUE
      IF(CSOF.GE.0.0D0.AND.SSOF.LT.0.0D0) SOF=PI2-SOF
      IF(CSOF.LT.0.0D0.AND.SSOF.GE.0.0D0) SOF=PI-SOF
      IF(CSOF.LT.0.0D0.AND.SSOF.LT.0.0D0) SOF=PI+SOF
      YSO=SOF-F
      IF(YSO.LT.0.0D0) YSO=YSO+PI2
      EA(3)=YSO/PI180
      GOTO 70
  60  CONTINUE
      EA(3)=0.0D0
C
  70  CONTINUE
C
      RETURN
C
      END
C ----------------------------------------------------------------------
      SUBROUTINE EARTHrv(T20,X,Y,Z,VXZ,VYZ,VZZ,WEL)
C ----------------------------------------------------------------------
C     - CALCULATION OF QUANTITIES CHARACTERIZING THE POSITION OF THE
C       EARTH AND ITS ORBIT AT GIVEN TIME (RECTANGULAR COORDINATES OF
C       THE EARTH WERE CALCULATED BY: BRETAGNON P.: 1982, ASTRON.
C       ASTROPHYS. 114, 278.)
C
C INPUT:
C    T20 - THE TIME [JULIAN CENTURIES FROM JD = 2451545.0]
C
C OUTPUT = THE VALUES VALID IN TIME "T20" (EQUINOX 2000.0):
C    X, Y, Z - HELIOCENTRIC RECTANGULAR COORDINATES REFERRED TO THE
C              ECLIPTIC [AU]
C    WEL(1) - SEMI-MAJOR AXIS [AU]
C    WEL(2) - ECCENTRICITY
C    WEL(3) - ECLIPTIC MEAN LONGITUDE [RADIANS]
C    WEL(4) - LONGITUDE OF PERIHELION [RADIANS]
C    WEL(5) - LONGITUDE OF ASCENDING NODE [RADIANS]
C    WEL(6) - INCLINATION TO THE ECLIPTIC [RADIANS]
C    WEL(7) - TRUE ANOMALY [RADIANS]
C
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      IMPLICIT REAL*8 (A-H,K-Z)
      DIMENSION WEL(7)
C
      PI=4.0D0*DATAN(1.0D0)
      PI2=2.0D0*PI
      GAUSS=1.7202098950D-2
C
      WEL(1)=1.00000101778D0
C
      WEL(3)=1.75347031435D0+6.28307584918D2*T20
      WEL(3)=WEL(3)-9.9189D-8*T20*T20+7.3D-13*T20*T20*T20
   3  CONTINUE
      IF(WEL(3).GE.0.0D0) GOTO 5
      WEL(3)=WEL(3)+PI2
      GOTO 3
   5  CONTINUE
      IF(WEL(3).LT.PI2) GOTO 7
      WEL(3)=WEL(3)-PI2
      GOTO 5
C
   7  WK=-3.7408165D-3-8.2266699D-5*T20
      WK=WK+2.748939D-7*T20*T20+1.04217D-9*T20*T20*T20
      WH=1.628447663D-2-6.2030259D-5*T20
      WH=WH-3.353888D-7*T20*T20+7.1185D-10*T20*T20*T20
      WQ=-1.13469002D-4*T20+1.237314D-7*T20*T20+1.2705D-9*T20*T20*T20
      WP=1.0180391D-5*T20+4.701998D-7*T20*T20-5.3829D-10*T20*T20*T20
C 
      WEL(2)=DSQRT(WK*WK+WH*WH)
C
      IF(DABS(WK).LT.1.0D-12) GOTO 210
      WEL(4)=DATAN(WH/WK)
      WEL(4)=DABS(WEL(4))
      IF(WH.GE.0.0D0.AND.WK.LT.0.0D0) WEL(4)=PI-WEL(4)
      IF(WH.LT.0.0D0.AND.WK.LT.0.0D0) WEL(4)=PI+WEL(4)
      IF(WH.LT.0.0D0.AND.WK.GE.0.0D0) WEL(4)=PI2-WEL(4)
C
  10  CONTINUE
      IF(DABS(WQ).LT.1.0D-12) GOTO 220
      WEL(5)=DATAN(WP/WQ)
      WEL(5)=DABS(WEL(5))
      IF(WP.GE.0.0D0.AND.WQ.LT.0.0D0) WEL(5)=PI-WEL(5)
      IF(WP.LT.0.0D0.AND.WQ.LT.0.0D0) WEL(5)=PI+WEL(5)
      IF(WP.LT.0.0D0.AND.WQ.GE.0.0D0) WEL(5)=PI2-WEL(5)
  20  WG=DSQRT(WP*WP+WQ*WQ)
      POM=1.0D0-WG*WG
      IF(DABS(POM).LT.1.0D-12) GOTO 230
      POM=WG/DSQRT(POM)
      WEL(6)=2.0D0*DATAN(POM)
      IF(WEL(6).LT.0.0D0) WEL(6)=PI+WEL(6)
C
  30  CONTINUE
      WM=WEL(3)-WEL(4)
      WEL(7)=WM+(2.0D0-WEL(2)*WEL(2)/4.0D0)*WEL(2)*DSIN(WM)
      WEL(7)=WEL(7)+5.0D0*WEL(2)*WEL(2)/4.0D0*DSIN(2.0D0*WM)
      WEL(7)=WEL(7)+1.3D1*WEL(2)*WEL(2)*WEL(2)/1.2D1*DSIN(3.0D0*WM)
  40  CONTINUE
      IF(WEL(7).GE.0.0D0) GOTO 50
      WEL(7)=WEL(7)+PI2
      GOTO 40
  50  CONTINUE
      IF(WEL(7).LT.PI2) GOTO 60
      WEL(7)=WEL(7)-PI2
      GOTO 50
  60  RZ=WEL(1)*(1.0D0-WEL(2)*WEL(2))/(1.0D0+WEL(2)*DCOS(WEL(7)))
      W2=WEL(7)+WEL(4)-WEL(5)
      C2O=DCOS(W2)
      S2O=DSIN(W2)
      CVO=DCOS(WEL(5))
      SVO=DSIN(WEL(5))
      CSK=DCOS(WEL(6))
      SSK=DSIN(WEL(6))
      X=RZ*(C2O*CVO-S2O*CSK*SVO)
      Y=RZ*(C2O*SVO+S2O*CSK*CVO)
      Z=RZ*S2O*SSK
      VELZ=GAUSS*DSQRT(2.0D0/RZ-1.0D0/WEL(1))

      MEN2=1.0D0+WEL(2)*WEL(2)+2.0D0*WEL(2)*DCOS(wel(7))
      MEN=DSQRT(MEN2)
      CVFI=-WEL(2)*DSIN(WEL(7))/MEN
      SVFI=(1.0D0+WEL(2)*DCOS(WEL(7)))/MEN
      VFI=DATAN(SVFI/CVFI)
      VFI=DABS(VFI)
      if(CVFI.LT.0.0D0.AND.SVFI.GE.0.0D0) VFI=PI-VFI
      if(CVFI.LT.0.0D0.AND.SVFI.LT.0.0D0) VFI=PI+VFI
      if(CVFI.GE.0.0D0.AND.SVFI.LT.0.0D0) VFI=PI2-VFI
      SOM=WEL(4)-WEL(5)
      C3Q=DCOS(SOM+WEL(7)-VFI)
      S3Q=DSIN(SOM+WEL(7)-VFI)
      CGO=DCOS(WEL(5))
      SGO=DSIN(WEL(5))
      CSK=DCOS(WEL(6))
      SSK=DSIN(WEL(6))
      VXZ=-VELZ*(C3Q*CGO-S3Q*CSK*SGO)
      VYZ=-VELZ*(C3Q*SGO+S3Q*CSK*CGO)
      VZZ=-VELZ*S3Q*SSK
      RETURN
C
 210  WEL(4)=PI/2.0D0
      IF(WH.LT.0.0D0) WEL(4)=1.5D0*PI
      GOTO 10
 220  WEL(5)=PI/2.0D0
      IF(WP.LT.0.0D0) WEL(5)=1.5D0*PI
      GOTO 20
 230  WEL(6)=PI
      GOTO 30
C
      END
