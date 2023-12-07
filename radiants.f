C **********************************************************************
C "METHODS.FOR" - PROGRAM COMPUTING THEORETICAL METEOR SHOWER RADIANTS
C                 (VERSION 2)
C               - TO EQUINOX 2000.0
C **********************************************************************
C
C STATUS OF USING:
C ====================
C
C   1. THE USE OF THE PROGRAM AS A SUBJECT OF COMMERCE IS NOT ALLOWED.
C   2. FOR OTHER PURPOSES, THE AUTHORS PERMIT TO COPY AND USE THE
C      PROGRAM OR ITS PART(S) FREELY;
C   3. EACH USING OF THE PROGRAM OR ANY OF ITS PARTS (SUBROUTINES)
C      SHOULD BE QUOTED IN THE OFFICIAL PRESENTATION; THE PROGRAM WAS
C      INTRODUCED IN THE "ASTRONOMY AND ASTROPHYSICS, LETTERS", PLEASE
C      SEE: NESLUSAN L., SVOREN J., PORUBCAN V.: ASTRON. ASTROPHYS., 
C      SUBMITTED.
C
C **********************************************************************
C
C BRIEF DESCRIPTION:
C ==================
C
C      THE PROGRAM COMPUTES A THEORETICAL RADIANT OF A METEOR SHOWER,
C      IF THE ORBIT OF AN ASSUMED PARENT BODY IS WELL-KNOWN
C       
C      THE PROGRAM INCLUDES SIX DIFFERENT METHODS:
C      - (Q) Q-ADJUSTMENT (HASEGAWA; 1990)
C      - (B) VARIATION OF PERIHELION DISTANCE AND EXCENTRICITY (SVOREN
C        ET AL; 1993)
C      - (W) ROTATION OF THE LINE OF APSIDES (STEEL AND BAGGALEY; 1985)
C      - (A) ROTATION AROUND THE LINE OF APSIDES (SVOREN ET AL; 1993)
C      - (H) OMEGA-ADJUSTMENT (HASEGAWA; 1990)
C      - (P) PORTER'S METHOD (PORTER; 1952)
C     (FOR MORE DETAILED DESCRIPTION AND REFERENCES, PLEASE, SEE IN THE
C     ABOVE MENTIONED PAPERS. REFERENCES ARE ALSO GIVEN IN BRIEF INDI-
C     VIDUAL DESCRIPTIONS OF SUBROUTINES BELOW.)
C
C ---------------------------------------------------------------------
C MAIN PROGRAM
C
C   INPUT:
C   ORBITAL ELEMENTS OF THE ORBIT OF THE PARENT BODY
C      EA(1) - PERIHELION DISTANCE [AU]
C      EA(2) - ECCENTRICITY
C      EA(3) - ARGUMENT OF PERIHELION [DEGREES]
C      EA(4) - LONGITUDE OF NODE [DEGREES]
C      EA(5) - INCLINATION TO ECLIPTIC [DEGREES]
C      YEAR, MONTH, DAY - TIME OF PERIHELION PASSAGE OF THE PARENT BODY
C      (CONVERTED TO: EA(6) - THE TIME IN JULIAN CENTURIES FROM
C                             JD = 2433282.423)
C
C   OUTPUT:
C   THE QUANTITIES DESCRIBING THEORETICAL RADIANTS AT BOTH ARCS (PRE-
C   PERIHELION AND POST-PERIHELION) OF THE ORBIT OF THE PARENT BODY;
C   THE OUTPUT QUANTITIES AT THE POST-PERIHELION ARC ARE:
C      R1(1) - SOLAR LONGITUDE OF THE MAXIMUM OF POTENTIAL SHOWER 
C              [DEGRRES]
C      R1(2) - RIGHT ASCENSION OF THE PREDICTED RADIANT [DEGREES]
C      R1(3) - DECLINATION OF THE RADIANT [DEGREES]
C      R1(4) - PREDICTED GEOCENTRIC VELOCITY OF METEORS [KM/S]
C      R1(5) - PREDICTED HELIOCENTRIC VELOCITY OF METEORS [KM/S]
C      R1(6) - PREDICTED TIME OF THE SHOWER MAXIMUM [JULIAN CENTURIES
C              FROM 1950.0]
C      WM1(1) - MOMENT OF MINIMUM DISTANCE BETWEEN THE EARTH AND THE
C               ORBIT OF THE PARENT BODY [JULIAN CENTURIES FROM 1950.0]
C      WM1(2) - MINIMUM DISTANCE BETWEEN POST-PERIHELION ARC OF THE
C               PARENT BODY ORBIT AND EARTH'S ORBIT [AU]
C      WM1(3) - DISTANCE OF THE PARENT BODY FROM THE EARTH AT THE MOMENT
C               OF MAXIMUM [AU]
C      WM1(4) - TIME INTERVAL BETWEEN THE PASSAGES OF THE PARENT BODY
C               AND EARTH ACROSS THE NEAREST POINTS OF THE POST-PERI-
C               HELION ARC OF THE PARENT BODY ORBIT AND EARTH'S ORBIT
C               [DAYS]
C      WM1(5) - TRUE ANOMALY OF THE PARENT BODY AT THE MOMENT OF ITS
C               NEAREST POST-PERIHELION APPROACH TO THE EARTH'S ORBIT
C               [DEGREES]
C      WM1(6) - TRUE ANOMALY OF THE EARTH AT THE MOMENT OF ITS NEAREST
C               APPROACH TO THE POST-PERIHELION ARC OF THE PARENT BODY
C               ORBIT [DEGREES]
C   ORBITAL ELEMENTS OF MODIFIED ORBIT CROSSING THE EARTH'S ORBIT
C      EB1(1) - PERIHELION DISTANCE [AU]
C      EB1(2) - ECCENTRICITY
C      EB1(3) - ARGUMENT OF PERIHELION [DEGREES]
C      EB1(4) - LONGITUDE OF NODE [DEGREES]
C      EB1(5) - INCLINATION TO ECLIPTIC [DEGREES]
C      DD1 - D-CRITERION CHARACTERIZING THE DIFFERENCE BETWEEN THE ORBIT
C            OF THE PARENT BODY (INPUT ORBIT) AND MODIFIED ORBIT
C            CROSSING THE EARTH'S ORBIT (OUTPUT ORBIT)
C   THE ANALOGOUS QUANTITIES CONCERNING THE PRE-PERIHELION ARC ARE
C   STORED IN THE FIELDS OF VALUES "R2(I)", "WM2(I)", AND "EB2(I)",
C   WHERE I=1, 2,..., 6 (5, RESPECTIVELY), AND IN VALUE "DD2".
C
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER  I,IM,IYR,J
C
      CHARACTER*1 MET
      character*4 MON
      character*24 filin,fil15,fil16,fil17
      character*72 ret72
      DIMENSION EA(6),MON(12),EB(6),MET(6)
      DIMENSION EB1(6),EB2(6),R1(6),R2(6),RY1(6,6),RY2(6,6)
      DIMENSION EY1(6,6),EY2(6,6),DDY1(6),DDY2(6)
      DIMENSION WM1(6),WM2(6),iau(5000),isol(5000)
      dimension sollam(5000),al(5000),dl(5000),vg(5000)
      dimension q(5000),e(5000),som(5000),gom(5000),sk(5000),jnn(5000)
C
      pi=4.0d0*datan(1.0d0)
      pi180=pi/1.80d2

c Tolerance limits:
      open(unit=9,file='inparams.rad',access='sequential')
        read(9,*) ret72
        read(9,*) filin
        read(9,*) ret72
        read(9,*) fil15
        read(9,*) ret72
        read(9,*) fil16
        read(9,*) ret72
        read(9,*) DELLS
        read(9,*) ret72
        read(9,*) DELAL
        read(9,*) ret72
        read(9,*) DELDL
        read(9,*) ret72
        read(9,*) DELVG
        read(9,*) ret72
        read(9,*) iyn
      close(unit=9)

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
      MET(1)='Q'
      MET(2)='B'
      MET(3)='W'
      MET(4)='A'
      MET(5)='H'
      MET(6)='P'

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
      if(iyn.eq.1) OPEN(UNIT=17,FILE='debug.d',ACCESS='SEQUENTIAL')
      WRITE(15,500)
      WRITE(15,510)

      IROKI=2000
      DEN=1.0D0
      IMES=1
      CALL JULIAN(IROKI,IMES,DEN,T19,TMIN)
      IROKI=IROKI+1
      CALL JULIAN(IROKI,IMES,DEN,T19,TMAX)
      IROKI=IROKI-1

      DO 100 i=1,jall
      WRITE(*,*) 'METEOR SHOWER - IAU No., sol. = ',iau(i),isol(i)
      if(iyn.eq.1) WRITE(17,700) iau(i),isol(i)
      EA(1)=q(i)
      EA(2)=e(i)
      EA(3)=som(i)
      EA(4)=gom(i)
      EA(5)=sk(i)
C
c      IROK=YY(I,3)
c      IMES=YY(I,4)
c      DEN=YY(I,5)
cClll      DEN=17.0
c      CALL JULIAN(IROK,IMES,DEN,T19,ETT)
c      EA(6)=ETT

      CALL APPROACH(TMIN,TMAX,EA,WM1,WM2)
      if(iyn.eq.1) WRITE(17,*) 'MOID = ',WM1(2),WM2(2)

      zemlam=sollam(i)+1.80d2
      if(zemlam.ge.3.60d2) zemlam=zemlam-3.60d2
      ttv=(zemlam*pi180-1.75347031435d0)/6.28307584918d3
      ETT=3.65250d5*ttv+2.4515450d6

 210  CONTINUE
      CALL QMETH(ETT,EA,R1,EB1,DD1,R2,EB2,DD2)
Ct      WRITE(*,*) 'R1(6) = ',R1(6)
Ct      WRITE(*,*) 'R1(1), GO = ',R1(1),EA(4)
Ct      WRITE(*,*) 'Von z QMETH.'
      if(iyn.eq.1) then
      WRITE(17,*) '     Q-method:'
      WRITE(17,*) 'q = ',EA(1),EB1(1),EB2(1)
      WRITE(17,*) 'e = ',EA(2),EB1(2),EB2(2)
      WRITE(17,*) 'omega = ',EA(3),EB1(3),EB2(3)
      WRITE(17,*) 'Omega = ',EA(4),EB1(4),EB2(4)
      WRITE(17,*) 'i = ',EA(5),EB1(5),EB2(5)
      WRITE(17,*) 'D_1, D_2 = ',DD1,DD2
      WRITE(17,*) 'lambda_sun = ',sollam(i),R1(1),R2(1)
      WRITE(17,*) 'alpha = ',al(i),R1(2),R2(2)
      WRITE(17,*) 'delta = ',dl(i),R1(3),R2(3)
      WRITE(17,*) 'Vg = ',vg(i),R1(4),R2(4)
      WRITE(17,*) 'Vh(calc.) = ',R1(5),R2(5)
      endif
      DO 215 J=1,6
      RY1(1,J)=R1(J)
      RY2(1,J)=R2(J)
      EY1(1,J)=EB1(J)
      EY2(1,J)=EB2(J)
 215  CONTINUE
      DDY1(1)=DD1
      DDY2(1)=DD2
Cmmm      J1=1
Cmmm      J2=1
Cmmm      GOTO 283
C
 220  CONTINUE
      CALL BMETH(ETT,EA,R1,EB1,DD1,R2,EB2,DD2)
Ct      WRITE(*,*) 'Von z BMETH.'
      if(iyn.eq.1) then
      WRITE(17,*) '     B-method:'
      WRITE(17,*) 'q = ',EA(1),EB1(1),EB2(1)
      WRITE(17,*) 'e = ',EA(2),EB1(2),EB2(2)
      WRITE(17,*) 'omega = ',EA(3),EB1(3),EB2(3)
      WRITE(17,*) 'Omega = ',EA(4),EB1(4),EB2(4)
      WRITE(17,*) 'i = ',EA(5),EB1(5),EB2(5)
      WRITE(17,*) 'D_1, D_2 = ',DD1,DD2
      WRITE(17,*) 'lambda_sun = ',sollam(i),R1(1),R2(1)
      WRITE(17,*) 'alpha = ',al(i),R1(2),R2(2)
      WRITE(17,*) 'delta = ',dl(i),R1(3),R2(3)
      WRITE(17,*) 'Vg = ',vg(i),R1(4),R2(4)
      WRITE(17,*) 'Vh(calc.) = ',R1(5),R2(5)
      endif
      DO 225 J=1,6
      RY1(2,J)=R1(J)
      RY2(2,J)=R2(J)
      EY1(2,J)=EB1(J)
      EY2(2,J)=EB2(J)
 225  CONTINUE
      DDY1(2)=DD1
      DDY2(2)=DD2
C
 230  CONTINUE
      CALL WMETH(TMIN,TMAX,EA,R1,EB1,DD1,R2,EB2,DD2)
Ct      WRITE(*,*) 'Von z WMETH.'
      if(iyn.eq.1) then
      WRITE(17,*) '     W-method:'
      WRITE(17,*) 'q = ',EA(1),EB1(1),EB2(1)
      WRITE(17,*) 'e = ',EA(2),EB1(2),EB2(2)
      WRITE(17,*) 'omega = ',EA(3),EB1(3),EB2(3)
      WRITE(17,*) 'Omega = ',EA(4),EB1(4),EB2(4)
      WRITE(17,*) 'i = ',EA(5),EB1(5),EB2(5)
      WRITE(17,*) 'D_1, D_2 = ',DD1,DD2
      WRITE(17,*) 'lambda_sun = ',sollam(i),R1(1),R2(1)
      WRITE(17,*) 'alpha = ',al(i),R1(2),R2(2)
      WRITE(17,*) 'delta = ',dl(i),R1(3),R2(3)
      WRITE(17,*) 'Vg = ',vg(i),R1(4),R2(4)
      WRITE(17,*) 'Vh(calc.) = ',R1(5),R2(5)
      endif
      DO 235 J=1,6
      RY1(3,J)=R1(J)
      RY2(3,J)=R2(J)
      EY1(3,J)=EB1(J)
      EY2(3,J)=EB2(J)
 235  CONTINUE
      DDY1(3)=DD1
      DDY2(3)=DD2
C
 240  CONTINUE
      CALL AMETH(TMIN,TMAX,EA,R1,EB1,DD1,R2,EB2,DD2)
Ct      WRITE(*,*) 'Von z AMETH.'
      if(iyn.eq.1) then
      WRITE(17,*) '     A-method:'
      WRITE(17,*) 'q = ',EA(1),EB1(1),EB2(1)
      WRITE(17,*) 'e = ',EA(2),EB1(2),EB2(2)
      WRITE(17,*) 'omega = ',EA(3),EB1(3),EB2(3)
      WRITE(17,*) 'Omega = ',EA(4),EB1(4),EB2(4)
      WRITE(17,*) 'i = ',EA(5),EB1(5),EB2(5)
      WRITE(17,*) 'D_1, D_2 = ',DD1,DD2
      WRITE(17,*) 'lambda_sun = ',sollam(i),R1(1),R2(1)
      WRITE(17,*) 'alpha = ',al(i),R1(2),R2(2)
      WRITE(17,*) 'delta = ',dl(i),R1(3),R2(3)
      WRITE(17,*) 'Vg = ',vg(i),R1(4),R2(4)
      WRITE(17,*) 'Vh(calc.) = ',R1(5),R2(5)
      endif
      DO 245 J=1,6
      RY1(4,J)=R1(J)
      RY2(4,J)=R2(J)
      EY1(4,J)=EB1(J)
      EY2(4,J)=EB2(J)
 245  CONTINUE
      DDY1(4)=DD1
      DDY2(4)=DD2
C
 250  CONTINUE
      CALL HMETH(TMIN,TMAX,EA,R1,EB1,DD1,R2,EB2,DD2)
Ct      WRITE(*,*) 'Von z HMETH.'
      if(iyn.eq.1) then
      WRITE(17,*) '     H-method:'
      WRITE(17,*) 'q = ',EA(1),EB1(1),EB2(1)
      WRITE(17,*) 'e = ',EA(2),EB1(2),EB2(2)
      WRITE(17,*) 'omega = ',EA(3),EB1(3),EB2(3)
      WRITE(17,*) 'Omega = ',EA(4),EB1(4),EB2(4)
      WRITE(17,*) 'i = ',EA(5),EB1(5),EB2(5)
      WRITE(17,*) 'D_1, D_2 = ',DD1,DD2
      WRITE(17,*) 'lambda_sun = ',sollam(i),R1(1),R2(1)
      WRITE(17,*) 'alpha = ',al(i),R1(2),R2(2)
      WRITE(17,*) 'delta = ',dl(i),R1(3),R2(3)
      WRITE(17,*) 'Vg = ',vg(i),R1(4),R2(4)
      WRITE(17,*) 'Vh(calc.) = ',R1(5),R2(5)
      endif
      DO 255 J=1,6
      RY1(5,J)=R1(J)
      RY2(5,J)=R2(J)
      EY1(5,J)=EB1(J)
      EY2(5,J)=EB2(J)
 255  CONTINUE
      DDY1(5)=DD1
      DDY2(5)=DD2
C
 260  CONTINUE
      CALL PMETH(WM1,WM2,EA,R1,EB1,DD1,R2,EB2,DD2)
Ct      WRITE(*,*) 'Von z PMETH.'
      if(iyn.eq.1) then
      WRITE(17,*) '     P-method:'
      WRITE(17,*) 'q = ',EA(1),EB1(1),EB2(1)
      WRITE(17,*) 'e = ',EA(2),EB1(2),EB2(2)
      WRITE(17,*) 'omega = ',EA(3),EB1(3),EB2(3)
      WRITE(17,*) 'Omega = ',EA(4),EB1(4),EB2(4)
      WRITE(17,*) 'i = ',EA(5),EB1(5),EB2(5)
      WRITE(17,*) 'D_1, D_2 = ',DD1,DD2
      WRITE(17,*) 'lambda_sun = ',sollam(i),R1(1),R2(1)
      WRITE(17,*) 'alpha = ',al(i),R1(2),R2(2)
      WRITE(17,*) 'delta = ',dl(i),R1(3),R2(3)
      WRITE(17,*) 'Vg = ',vg(i),R1(4),R2(4)
      WRITE(17,*) 'Vh(calc.) = ',R1(5),R2(5)
      WRITE(17,*) ' '
      endif
      DO 265 J=1,6
      RY1(6,J)=R1(J)
      RY2(6,J)=R2(J)
      EY1(6,J)=EB1(J)
      EY2(6,J)=EB2(J)
 265  CONTINUE
      DDY1(6)=DD1
      DDY2(6)=DD2
C
      DD1=DDY1(1)
      J1=1
      DO 270 J=2,6
      IF(DDY1(J).LT.DD1) J1=J
      IF(DDY1(J).LT.DD1) DD1=DDY1(J)
 270  CONTINUE
      DD2=DDY2(1)
      J2=1
      DO 280 J=2,6
      IF(DDY2(J).LT.DD2) J2=J
      IF(DDY2(J).LT.DD2) DD2=DDY2(J)
 280  CONTINUE
C
 283  CONTINUE
      IF(DD1.LT.9.9999D1.OR.DD2.LT.9.9999D1) GOTO 80
      YQ=0.0D0
      YE=0.0D0
      YSO=0.0D0
      YGO=0.0D0
      YI=0.0D0
      YLAM=0.0D0
      ALB=0.0D0
      DLB=0.0D0
      VGB=0.0D0
      VHB=0.0D0
      GOTO 180
C
  80  CONTINUE
C ???      IF(DABS(RY2(J2,6)-EA(6)).LT.DABS(RY1(J1,6)-EA(6))) GOTO 135
      IF(DD2.LT.DD1) GOTO 135
      TJD=RY1(J1,6)
      CALL INVJUL(TJD,IRK,IMS,DNI)
Ct      WRITE(*,*) 'TJD, IRK, IMS, DNI = ',TJD,IRK,IMS,DNI
      YQ=EY1(J1,1)
      YE=EY1(J1,2)
      YSO=EY1(J1,3)
      YGO=EY1(J1,4)
      YI=EY1(J1,5)
      YLAMB=RY1(J1,1)
      ALB=RY1(J1,2)
      DLB=RY1(J1,3)
      VGB=RY1(J1,4)
      VHB=RY1(J1,5)
      JMET=J1
      GOTO 180
 135  CONTINUE      
      TJD=RY2(J2,6)
      CALL INVJUL(TJD,IRK,IMS,DNI)
      YQ=EY2(J2,1)
      YE=EY2(J2,2)
      YSO=EY2(J2,3)
      YGO=EY2(J2,4)
      YI=EY2(J2,5)
      YLAMB=RY2(J2,1)
      ALB=RY2(J2,2)
      DLB=RY2(J2,3)
      VGB=RY2(J2,4)
      VHB=RY2(J2,5)
      JMET=J2
C
 180  CONTINUE
      YLAMA=gom(i)
      IF(DABS(YLAMA-YLAMB).GT.1.80D2) YLAMA=YLAMA+1.80D2
      IF(YLAMA.GE.3.60D2) YLAMA=YLAMA-3.60D2
      WRITE(15,600) iau(i),isol(i),YLAMA,al(i),dl(i),vg(i)
      WRITE(15,650) YLAMB,ALB,DLB,VGB,MET(JMET)
C
      lamx=sollam(i)
      if(lamx.lt.9.0d1.and.YLAMB.gt.2.70d2) lamx=lamx+3.60d2
      if(lamx.gt.2.70d2.and.YLAMB.lt.9.0d1) lamx=lamx-3.60d2
      IF(DABS(YLAMB-lamx).GT.DELLS) THEN
      WRITE(16,*) iau(i),isol(i),'   LS_obs,   LS_calc = ',sollam(i),YLA
     *MB
      END IF
      alx=al(i)
      if(alx.lt.9.0d1.and.ALB.gt.2.70d2) alx=alx+3.60d2
      if(alx.gt.2.70d2.and.ALB.lt.9.0d1) alx=alx-3.60d2
      IF(DABS(ALB-alx).GT.DELAL) THEN
      WRITE(16,*) iau(i),isol(i),'   RA_obs,   RA_calc = ',al(i),ALB
      END IF
      IF(DABS(DLB-dl(i)).GT.DELDL) THEN
      WRITE(16,*) iau(i),isol(i),'  DEC_obs,  DEC_calc = ',dl(i),DLB
      END IF
      IF(DABS(VGB-vg(i)).GT.DELVG) THEN
      WRITE(16,*) iau(i),isol(i),'   Vg_obs,   Vg_calc = ',vg(i),VGB
      END IF
c      IF(DABS(VHB-vh(i)).GT.DELVG) THEN
c      WRITE(16,*) iau(i),isol(i),'   Vh_obs,   Vh_calc = ',vh(i),VHB
c      END IF
CC      CLOSE(UNIT=15)
C
  95  CONTINUE
 100  CONTINUE
      CLOSE(UNIT=15)
      CLOSE(UNIT=16)
      if(iyn.eq.1) CLOSE(UNIT=17)
C
      STOP
C
 500  FORMAT(' IAU sol.  LS       RA      DEC      Vg   method')
 510  FORMAT('------------------------------------------------')
 600  FORMAT(i4,i3,3f9.3,f8.3)
 610  FORMAT(A8,1X,A6,2X,I4,I3,F9.5,3F9.3,F8.3)
 620  FORMAT(A8,1X,A6,2X,I4,I3,F9.5,9X,2F9.3,2F8.3)
 650  FORMAT(7x,3f9.3,f8.3,3x,a1,/)
 700  format(i4,i3)
      END 
C ---------------------------------------------------------------------
      SUBROUTINE QMETH(ETT,EA,R1,EB1,DD1,R2,EB2,DD2)
C ---------------------------------------------------------------------
C      - (Q) Q-ADJUSTMENT
C     (HASEGAWA I.: 1990, PUBL. ASTRON. SOC. JAPAN 42, 175.)
C
C     "TMIN", "TMAX" - BEGINNING AND END OF YEAR OF INVESTIGATION OF THE
C     POTENTIAL SHOWER GIVEN IN JULIAN CENTURIES FROM 1950.0. NEXT
C     INPUT (FIELD OF VALUES "EA(I)", I = 1, 2,..., 6) AND OUTPUT (FIEL-
C     DS OF VALUES "R1(I)", "R2(I)", "EB1(I)", AND "EB2(I)", I = 1, 2,
C     3,..., 6 OR 5, RESPECTIVELY), PLEASE, SEE THE MAIN PROGRAM. AS
C     THE OUTPUT, THERE ARE MOREOVER THE VALUES "DD1" AND "DD2" OF D-
C     DISCRIMINANT CORRESPONDING TO THE POST- AND PRE-PERIHELION ARC
C     OF THE PARENT BODY ORBIT, RESPECTIVELY.
C
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER  I,IM,IYR,J
      DIMENSION EA(6),EB1(5),EB2(5),R1(6),R2(6),WEL(7),EZ(6)
C
      PI=4.0D0*DATAN(1.0D0)
      PI2=2.0D0*PI
      PI180=PI/0.18D3
      EPS=2.3439291D1
      CEPS=DCOS(PI180*EPS)
      SEPS=DSIN(PI180*EPS)
      NTARG=3
      NCENT=11
C
      DO 50 J=2,5
      EB1(J)=EA(J)
      EB2(J)=EA(J)
 50   CONTINUE
      DT0=(TMAX-TMIN)/3.6525D2
      DT05=DT0/2.0D0
      DD1=0.999D99
      DD2=0.999D99
C
      IF(EA(3).GT.0.18D3) GOTO 100
      FB=0.18D3-EA(3)
      IF(FB.GT.1.799999999D2.AND.EA(2).GT.0.9999999D0) GOTO 200
      VO0=EA(4)+0.18D3
      IF(VO0.GE.0.36D3) VO0=VO0-0.36D3
      GOTO 110
 100  FB=0.36D3-EA(3)
      VO0=EA(4)
 110  VO0=VO0*PI180
      DWM=2.0D0*PI2
      TA=TMIN
      TB=TMAX
      DT=DT0
C
 120  CONTINUE
C>>>      DO 130 ETT=TA,TB,DT
c      CALL PLEPH(ETT,NTARG,NCENT,RRD)
      T20=(ETT-2.451545D6)/3.6525D4
      CALL EARTHrv(T20,XZ,YZ,ZZ,VXZ,VYZ,VZZ,WEL)
c      XZ=RRD(1)
c      YZ=RRD(2)*CEPS+RRD(3)*SEPS
C      ZZ=-RRD(2)*SEPS+RRD(3)*CEPS
C      VXZ=RRD(4)
C      VYZ=RRD(5)*CEPS+RRD(6)*SEPS
C      VZZ=-RRD(5)*SEPS+RRD(6)*CEPS
C      CALL ELEMS(ETT,XZ,YZ,VZ,VXZ,VYZ,VZZ,EZ)
      WL1=WEL(7)+WEL(4)
      IF(WL1.GE.PI2) WL1=WL1-PI2
      WL1=(WL1+PI)/PI180
CC      WRITE(*,*) 'Old value = ',WL1
CCC      T15=(ETT-2.433282423D6)/3.6525D4
CCC 121  CONTINUE
CCC      IF(WEL(3).LT.PI2) GOTO 123
CCC      WEL(3)=WEL(3)-PI2
CCC      GOTO 121
CCC 123  CONTINUE
C
      IF(DABS(XZ).LT.1.0D-12) THEN
      WL=PI/2.0D0
      IF(YZ.LT.0.0D0) WL=WL+PI
      ELSE
      WL=DATAN(YZ/XZ)
      WL=DABS(WL)
      IF(XZ.LT.0.0D0.AND.YZ.GE.0.0D0) WL=PI-WL
      IF(XZ.LT.0.0D0.AND.YZ.LT.0.0D0) WL=PI+WL
      IF(XZ.GE.0.0D0.AND.YZ.LT.0.0D0) WL=PI2-WL
      END IF
Ct      WLOUT=(WL+PI)/PI180
Ct      IF(WLOUT.GE.3.60D2) WLOUT=WLOUT-3.60D2
Ct      WRITE(*,*) 'JPL eph., lam = ',WLOUT
CCC      WLNEW=EZ(6)+EZ(3)+EZ(4)
CCC      IF(WLNEW.GE.PI2) WLNEW=WLNEW-PI2
CCC      IF(WLNEW.GE.PI2) WLNEW=WLNEW-PI2
CCC      WLNEW=(WLNEW+PI)/PI180
CCC      WRITE(*,*) 'New value = ',WLNEW
CCC      WRITE(*,*) ' '
CCC      STOP
C>>>     IF(DABS(WL-VO0).GT.DWM) GOTO 130
Cx      EZ(1)=WEL(1)
Cx      EZ(2)=WEL(2)
Cx      EZ(3)=WEL(4)-WEL(5)
Cx      EZ(4)=WEL(5)
Cx      EZ(5)=WEL(6)
Cx      FZM=WEL(7)
Cx      EZ(1)=EZ(1)/(1.0D0-EZ(2))
Cx      FZM=EZ(6)
c      RZM=DSQRT(RRD(1)*RRD(1)+RRD(2)*RRD(2)+RRD(3)*RRD(3))
      RZM=DSQRT(XZ*XZ+YZ*YZ+ZZ*ZZ)
      R1(1)=WL/PI180
      R1(6)=ETT
      DWM=DABS(WL-VO0)
C>>> 130  CONTINUE
C>>>      IF(DT.LT.DT05) GOTO 150
C>>>      TA=R1(6)-DT
C>>>      TB=R1(6)+DT
C>>>      DT=0.1D-1*DT
C>>>      GOTO 120
C
 150  QWV=FB*PI180
Ct      WRITE(*,*) 'R1(6) = ',R1(6)
Ct      WLOUT=R1(1)+1.80D2
Ct      IF(WLOUT.GE.3.60D2) WLOUT=WLOUT-3.60D2
Ct      WRITE(*,*) 'JPL eph., lam = ',WLOUT
Cx????      IF(EZ(3).LT.0.0D0) EZ(3)=EZ(3)+PI2
Cx      RZM=EZ(1)*(1.0D0-EZ(2)*EZ(2))/(1.0D0+EZ(2)*DCOS(FZM))
      EB1(1)=RZM*(1.0D0+EA(2)*DCOS(QWV))/(1.0D0+EA(2))
      IF(EB1(1).LT.0.0D0) GOTO 200
      DD1=DABS(EA(1)-EB1(1))
      R1(1)=R1(1)-0.18D3
      IF(R1(1).LT.0.0D0) R1(1)=R1(1)+0.36D3
      VYR=VYZ*CEPS-VZZ*SEPS
      VZR=VYZ*SEPS+VZZ*CEPS
      CALL RADIANT(EB1,FB,VXZ,VYR,VZR,R1)
C
 200  CONTINUE
      IF(EA(3).GT.0.18D3) GOTO 205
      FB=0.36D3-EA(3)
      IF(FB.LT.1.800000001D2.AND.EA(2).GT.0.9999999D0) GOTO 300
      VO0=EA(4)
      GOTO 210
 205  FB=0.54D3-EA(3)
      VO0=EA(4)+0.18D3
      IF(VO0.GE.0.36D3) VO0=VO0-0.36D3
 210  VO0=VO0*PI180
      DWM=2.0D0*PI2
      TA=TMIN
      TB=TMAX
      DT=DT0
C
 220  CONTINUE
C>>>      DO 230 ETT=TA,TB,DT
c      CALL PLEPH(ETT,NTARG,NCENT,RRD)
      T20=(ETT-2.451545D6)/3.6525D4
      CALL EARTHrv(T20,XZ,YZ,ZZ,VXZ,VYZ,VZZ,WEL)
      IF(DABS(XZ).LT.1.0D-12) THEN
      WL=PI/2.0D0
      IF(YZ.LT.0.0D0) WL=WL+PI
      ELSE
      WL=DATAN(YZ/XZ)
      WL=DABS(WL)
      IF(XZ.LT.0.0D0.AND.YZ.GE.0.0D0) WL=PI-WL
      IF(XZ.LT.0.0D0.AND.YZ.LT.0.0D0) WL=PI+WL
      IF(XZ.GE.0.0D0.AND.YZ.LT.0.0D0) WL=PI2-WL
      END IF
C      WL=WEL(7)+WEL(4)
C      IF(WL.GE.PI2) WL=WL-PI2
C>>>      IF(DABS(WL-VO0).GT.DWM) GOTO 230
Cx      EZ(1)=WEL(1)
Cx      EZ(2)=WEL(2)
CX      EZ(3)=WEL(4)-WEL(5)
Cx      EZ(4)=WEL(5)
Cx      EZ(5)=WEL(6)
Cx      FZM=WEL(7)
Cx      EZ(1)=EZ(1)/(1.0D0-EZ(2))
Cx      FZM=EZ(6)
      RZM=DSQRT(XZ*XZ+YZ*YZ+ZZ*ZZ)
      R2(1)=WL/PI180
      R2(6)=ETT
      DWM=DABS(WL-VO0)
C>>> 230  CONTINUE
C>>>      IF(DT.LT.DT05) GOTO 250
C>>>      TA=R2(6)-DT
C>>>      TB=R2(6)+DT
C>>>      DT=0.1D-1*DT
C>>>      GOTO 220
C
 250  QWV=FB*PI180
Cx ????      IF(EZ(3).LT.0.0D0) EZ(3)=EZ(3)+PI2
Cx      RZM=EZ(1)*(1.0D0-EZ(2)*EZ(2))/(1.0D0+EZ(2)*DCOS(FZM))
      EB2(1)=RZM*(1.0D0+EA(2)*DCOS(QWV))/(1.0D0+EA(2))
      IF(EB2(1).LT.0.0D0) GOTO 300
      DD2=DABS(EA(1)-EB2(1))
      R2(1)=R2(1)-0.18D3
      IF(R2(1).LT.0.0D0) R2(1)=R2(1)+0.36D3
      VYR=VYZ*CEPS-VZZ*SEPS
      VZR=VYZ*SEPS+VZZ*CEPS
      CALL RADIANT(EB2,FB,VXZ,VYR,VZR,R2)
C
 300  CONTINUE
      RETURN
C
      END
C ---------------------------------------------------------------------
      SUBROUTINE BMETH(ETT,EA,R1,EB1,DD1,R2,EB2,DD2)
C ---------------------------------------------------------------------
C      - (B) VARIATION OF PERIHELION DISTANCE AND EXCENTRICITY
C     (SVOREN J., NESLUSAN L., PORUBCAN V.: 1993, CONTRIB. ASTRON. OBS.
C     SKALNATE PLESO 23, 23.)
C
C     "TMIN", "TMAX" - BEGINNING AND END OF YEAR OF INVESTIGATION OF THE
C     POTENTIAL SHOWER GIVEN IN JULIAN CENTURIES FROM 1950.0. NEXT
C     INPUT (FIELD OF VALUES "EA(I)", I = 1, 2,..., 6) AND OUTPUT (FIEL-
C     DS OF VALUES "R1(I)", "R2(I)", "EB1(I)", AND "EB2(I)", I = 1, 2,
C     3,..., 6 OR 5, RESPECTIVELY), PLEASE, SEE THE MAIN PROGRAM. AS
C     THE OUTPUT, THERE ARE MOREOVER THE VALUES "DD1" AND "DD2" OF D-
C     DISCRIMINANT CORRESPONDING TO THE POST- AND PRE-PERIHELION ARC
C     OF THE PARENT BODY ORBIT, RESPECTIVELY.
C
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER  I,IM,IYR,J
      DIMENSION EA(6),EB1(5),EB2(5),R1(6),R2(6),WEL(7)
Cx      ,EZ(6)
C
      PI=4.0D0*DATAN(1.0D0)
      PI2=2.0D0*PI
      PI180=PI/0.18D3
      EPS=2.3439291D1
      CEPS=DCOS(PI180*EPS)
      SEPS=DSIN(PI180*EPS)
      NTARG=3
      NCENT=11
C
      DO 50 J=3,5
      EB1(J)=EA(J)
      EB2(J)=EA(J)
 50   CONTINUE
      DT0=(TMAX-TMIN)/3.6525D2
      DT05=DT0/2.0D0
      DD1=0.999D99
      DD2=0.999D99
C
      IF(EA(3).GT.0.18D3) GOTO 100
      FB=0.18D3-EA(3)
      IF(FB.GT.1.799999999D2.AND.EA(2).GT.0.9999999D0) GOTO 200
      VO0=EA(4)+0.18D3
      IF(VO0.GE.0.36D3) VO0=VO0-0.36D3
      GOTO 110
 100  FB=0.36D3-EA(3)
      VO0=EA(4)
 110  VO0=VO0*PI180
      DWM=2.0D0*PI2
      TA=TMIN
      TB=TMAX
      DT=DT0
C
 120  CONTINUE
C>>>      DO 130 ETT=TA,TB,DT
c      CALL PLEPH(ETT,NTARG,NCENT,RRD)
      T20=(ETT-2.451545D6)/3.6525D4
      CALL EARTHrv(T20,XZ,YZ,ZZ,VXZ,VYZ,VZZ,WEL)
      IF(DABS(XZ).LT.1.0D-12) THEN
      WL=PI/2.0D0
      IF(YZ.LT.0.0D0) WL=WL+PI
      ELSE
      WL=DATAN(YZ/XZ)
      WL=DABS(WL)
      IF(XZ.LT.0.0D0.AND.YZ.GE.0.0D0) WL=PI-WL
      IF(XZ.LT.0.0D0.AND.YZ.LT.0.0D0) WL=PI+WL
      IF(XZ.GE.0.0D0.AND.YZ.LT.0.0D0) WL=PI2-WL
      END IF
C      WL=WEL(7)+WEL(4)
C      IF(WL.GE.PI2) WL=WL-PI2
C>>>      IF(DABS(WL-VO0).GT.DWM) GOTO 130
Cx      EZ(1)=WEL(1)
Cx      EZ(2)=WEL(2)
Cx      EZ(3)=WEL(4)-WEL(5)
Cx      EZ(4)=WEL(5)
Cx      EZ(5)=WEL(6)
Cx      FZM=WEL(7)
Cx      EZ(1)=EZ(1)/(1.0D0-EZ(2))
Cx      FZM=EZ(6)
      RZM=DSQRT(XZ*XZ+YZ*YZ+ZZ*ZZ)
      VYR=VYZ*CEPS-VZZ*SEPS
      VZR=VYZ*SEPS+VZZ*CEPS
      R1(1)=WL/PI180
      R1(6)=ETT
      DWM=DABS(WL-VO0)
C>>> 130  CONTINUE
C>>>      IF(DT.LT.DT05) GOTO 150
C>>>      TA=R1(6)-DT
C>>>      TB=R1(6)+DT
C>>>      DT=0.1D-1*DT
C>>>      GOTO 120
C
 150  DCMIN=0.999D99
      CF=DCOS(FB*PI180)
Cx ???      IF(EZ(3).LT.0.0D0) EZ(3)=EZ(3)+PI2
Cx      RZM=EZ(1)*(1.0D0-EZ(2)*EZ(2))/(1.0D0+EZ(2)*DCOS(FZM))
C >      DO 160 EC=0.0D0,5.0D0,0.1D0
      DO 160 IEC=0,50
      EC=IEC/1.0D1
      DE=EC-EA(2)
      QX=RZM*(1.0D0+EC*CF)/(1.0D0+EC)
      DQ=EA(1)-QX
      DDX=DSQRT(DQ*DQ+DE*DE)
      IF(DDX.GE.DCMIN) GOTO 160
      DCMIN=DDX
      EB1(2)=EC
      EB1(1)=QX
 160  CONTINUE
C
      ED=EB1(2)-0.2D0
      IF(ED.LT.0.0D0) ED=0.0D0
      EU=EB1(2)+0.2D0
C >      DO 170 EC=ED,EU,0.1D-3
      IEU=(EU+1.0D-11-ED)/0.1D-3
      DO 170 IEC=0,IEU
      EC=ED+IEC*0.1D-3
      DE=EC-EA(2)
      QX=RZM*(1.0D0+EC*CF)/(1.0D0+EC)
      DQ=EA(1)-QX
      DDX=DSQRT(DQ*DQ+DE*DE)
      IF(DDX.GE.DCMIN) GOTO 170
      DCMIN=DDX
      EB1(2)=EC
      EB1(1)=QX
 170  CONTINUE
      R1(1)=R1(1)-0.18D3
      IF(R1(1).LT.0.0D0) R1(1)=R1(1)+0.36D3
      CALL RADIANT(EB1,FB,VXZ,VYR,VZR,R1)
      DD1=DCMIN
C
 200  CONTINUE
      IF(EA(3).GT.0.18D3) GOTO 205
      FB=0.36D3-EA(3)
      IF(FB.LT.1.800000001D2.AND.EA(2).GT.0.9999999D0) GOTO 300
      VO0=EA(4)
      GOTO 210
 205  FB=0.54D3-EA(3)
      VO0=EA(4)+0.18D3
      IF(VO0.GE.0.36D3) VO0=VO0-0.36D3
 210  VO0=VO0*PI180
      DWM=2.0D0*PI2
      TA=TMIN
      TB=TMAX
      DT=DT0
C
 220  CONTINUE
C>>>      DO 230 ETT=TA,TB,DT
c      CALL PLEPH(ETT,NTARG,NCENT,RRD)
      T20=(ETT-2.451545D6)/3.6525D4
      CALL EARTHrv(T20,XZ,YZ,ZZ,VXZ,VYZ,VZZ,WEL)
      IF(DABS(XZ).LT.1.0D-12) THEN
      WL=PI/2.0D0
      IF(YZ.LT.0.0D0) WL=WL+PI
      ELSE
      WL=DATAN(YZ/XZ)
      WL=DABS(WL)
      IF(XZ.LT.0.0D0.AND.YZ.GE.0.0D0) WL=PI-WL
      IF(XZ.LT.0.0D0.AND.YZ.LT.0.0D0) WL=PI+WL
      IF(XZ.GE.0.0D0.AND.YZ.LT.0.0D0) WL=PI2-WL
      END IF
C      WL=WEL(7)+WEL(4)
C      IF(WL.GE.PI2) WL=WL-PI2
C>>>      IF(DABS(WL-VO0).GT.DWM) GOTO 230
Cx      EZ(1)=WEL(1)
Cx      EZ(2)=WEL(2)
Cx      EZ(3)=WEL(4)-WEL(5)
Cx      EZ(4)=WEL(5)
Cx      EZ(5)=WEL(6)
Cx      FZM=WEL(7)
Cx      EZ(1)=EZ(1)/(1.0D0-EZ(2))
Cx      FZM=EZ(6)
      RZM=DSQRT(XZ*XZ+YZ*YZ+ZZ*ZZ)
      VYR=VYZ*CEPS-VZZ*SEPS
      VZR=VYZ*SEPS+VZZ*CEPS
      R2(1)=WL/PI180
      R2(6)=ETT
      DWM=DABS(WL-VO0)
C>>> 230  CONTINUE
C>>>      IF(DT.LT.DT05) GOTO 250
C>>>      TA=R2(6)-DT
C>>>      TB=R2(6)+DT
C>>>      DT=0.1D-1*DT
C>>>      GOTO 220
C
 250  DCMIN=0.999D99
      CF=DCOS(FB*PI180)
Cx ???      IF(EZ(3).LT.0.0D0) EZ(3)=EZ(3)+PI2
Cx      RZM=EZ(1)*(1.0D0-EZ(2)*EZ(2))/(1.0D0+EZ(2)*DCOS(FZM))
C >      DO 260 EC=0.0D0,5.0D0,0.1D0
      DO 260 IEC=0,50
      EC=IEC/1.0D1
      DE=EC-EA(2)
      QX=RZM*(1.0D0+EC*CF)/(1.0D0+EC)
      DQ=EA(1)-QX
      DDX=DSQRT(DQ*DQ+DE*DE)
      IF(DDX.GE.DCMIN) GOTO 260
      DCMIN=DDX
      EB2(2)=EC
      EB2(1)=QX
 260  CONTINUE
C
      ED=EB2(2)-0.2D0
      IF(ED.LT.0.0D0) ED=0.0D0
      EU=EB2(2)+0.2D0
C >      DO 270 EC=ED,EU,0.1D-3
      IEU=(EU+1.0D-11-ED)/0.1D-3
      DO 270 IEC=0,IEU
      EC=ED+IEC*0.1D-3
      DE=EC-EA(2)
      QX=RZM*(1.0D0+EC*CF)/(1.0D0+EC)
      DQ=EA(1)-QX
      DDX=DSQRT(DQ*DQ+DE*DE)
      IF(DDX.GE.DCMIN) GOTO 270
      DCMIN=DDX
      EB2(2)=EC
      EB2(1)=QX
 270  CONTINUE
      R2(1)=R2(1)-0.18D3
      IF(R2(1).LT.0.0D0) R2(1)=R2(1)+0.36D3
      CALL RADIANT(EB2,FB,VXZ,VYR,VZR,R2)
      DD2=DCMIN
C
 300  CONTINUE
      RETURN
C
      END
C ---------------------------------------------------------------------
      SUBROUTINE WMETH(TMIN,TMAX,EA,R1,EB1,DD1,R2,EB2,DD2)
C ---------------------------------------------------------------------
C      - (W) ROTATION OF THE LINE OF APSIDES
C     (STEEL D.I., BAGGALEY J.W.: 1985, MON. NOT. R. ASTRON. SOC. 212,
C     817.)
C
C     "TMIN", "TMAX" - BEGINNING AND END OF YEAR OF INVESTIGATION OF THE
C     POTENTIAL SHOWER GIVEN IN JULIAN CENTURIES FROM 1950.0. NEXT
C     INPUT (FIELD OF VALUES "EA(I)", I = 1, 2,..., 6) AND OUTPUT (FIEL-
C     DS OF VALUES "R1(I)", "R2(I)", "EB1(I)", AND "EB2(I)", I = 1, 2,
C     3,..., 6 OR 5, RESPECTIVELY), PLEASE, SEE THE MAIN PROGRAM. AS
C     THE OUTPUT, THERE ARE MOREOVER THE VALUES "DD1" AND "DD2" OF D-
C     DISCRIMINANT CORRESPONDING TO THE POST- AND PRE-PERIHELION ARC
C     OF THE PARENT BODY ORBIT, RESPECTIVELY.
C
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER  I,IM,IYR,J
      DIMENSION EA(6),EB1(5),EB2(5),R1(6),R2(6),WEL(7)
Cx      ,EZ(6)
C
      PI=4.0D0*DATAN(1.0D0)
      PI180=PI/0.18D3
      PI05=PI/2.0D0
      PI2=2.0D0*PI
      EPS=2.3439291D1
      CEPS=DCOS(PI180*EPS)
      SEPS=DSIN(PI180*EPS)
      NTARG=3
      NCENT=11
C
      DO 10 J=1,5
      EB1(J)=EA(J)
      EB2(J)=EA(J)
 10   CONTINUE      
      IF(EA(2).GT.0.999D0) GOTO 40
      IF(EA(1).GT.1.1D0) GOTO 20
      IF(EA(1)*(1.0D0+EA(2))/(1.0D0-EA(2)).LT.0.9D0) GOTO 20
      IF(EA(2).LT.1.0D-12) GOTO 20
      GOTO 40
 20   DD1=0.999D99
      DD2=0.999D99
      RETURN
C
 40   DT0=(TMAX-TMIN)/3.6525D2
      DT05=DT0/2.0D0
C
      IF(EA(3).GT.0.18D3) GOTO 100
      VO0=EA(4)+0.18D3
      IF(VO0.GE.0.36D3) VO0=VO0-0.36D3
      GOTO 110
 100  VO0=EA(4)
 110  VO0=VO0*PI180
      DWM=2.0D0*PI2
      TA=TMIN
      TB=TMAX
      DT=DT0
C
 120  CONTINUE
C >      DO 130 ETT=TA,TB,DT
      IETTU=(TB+1.0D-11-TA)/DT
      DO 130 IETT=0,IETTU
      ETT=TA+IETT*DT
c      CALL PLEPH(ETT,NTARG,NCENT,RRD)
      T20=(ETT-2.451545D6)/3.6525D4
      CALL EARTHrv(T20,XZ,YZ,ZZ,VXZ,VYZ,VZZ,WEL)
      IF(DABS(XZ).LT.1.0D-12) THEN
      WL=PI/2.0D0
      IF(YZ.LT.0.0D0) WL=WL+PI
      ELSE
      WL=DATAN(YZ/XZ)
      WL=DABS(WL)
      IF(XZ.LT.0.0D0.AND.YZ.GE.0.0D0) WL=PI-WL
      IF(XZ.LT.0.0D0.AND.YZ.LT.0.0D0) WL=PI+WL
      IF(XZ.GE.0.0D0.AND.YZ.LT.0.0D0) WL=PI2-WL
      END IF
C      WL=WEL(7)+WEL(4)
C      IF(WL.GE.PI2) WL=WL-PI2
      IF(DABS(WL-VO0).GT.DWM) GOTO 130
Cx      EZ(1)=WEL(1)
Cx      EZ(2)=WEL(2)
Cx      EZ(3)=WEL(4)-WEL(5)
Cx      EZ(4)=WEL(5)
Cx      EZ(5)=WEL(6)
Cx      FZM=WEL(7)
Cx      EZ(1)=EZ(1)/(1.0D0-EZ(2))
Cx      FZM=EZ(6)
      RZM=DSQRT(XZ*XZ+YZ*YZ+ZZ*ZZ)
      VYR=VYZ*CEPS-VZZ*SEPS
      VZR=VYZ*SEPS+VZZ*CEPS
      R1(1)=WL/PI180
      R1(6)=ETT
      DWM=DABS(WL-VO0)
 130  CONTINUE
      IF(DT.LT.DT05) GOTO 150
      TA=R1(6)-DT
      TB=R1(6)+DT
      DT=0.1D-1*DT
      GOTO 120
C
 150  CONTINUE
Cx ????      IF(EZ(3).LT.0.0D0) EZ(3)=EZ(3)+PI2
Cx      RZM=EZ(1)*(1.0D0-EZ(2)*EZ(2))/(1.0D0+EZ(2)*DCOS(FZM))
      CFB=(EA(1)*(1.0D0+EA(2))/RZM-1.0D0)/EA(2)
      DD1=0.999D99
      IF(DABS(CFB).GT.1.000001D0) GOTO 200
      IF(CFB.LT.-1.0D0) CFB=-1.0D0
      IF(CFB.GT.1.0D0) CFB=1.0D0
      IF(DABS(CFB).LT.1.0D-12) GOTO 160
      SFB=DSQRT(1.0D0-CFB*CFB)
      FB=DATAN(SFB/CFB)
      IF(FB.LT.0.0D0) FB=FB+PI
      GOTO 170
 160  FB=PI05
 170  FB=FB/PI180
      EB1(3)=0.18D3-FB
      IF(EA(3).GT.0.18D3) EB1(3)=0.36D3-FB
      CALL DCRIT(EA,EB1,DD1)
      R1(1)=R1(1)-0.18D3
      IF(R1(1).LT.0.0D0) R1(1)=R1(1)+0.36D3
      CALL RADIANT(EB1,FB,VXZ,VYR,VZR,R1)
C
 200  CONTINUE
      IF(EA(3).GT.0.18D3) GOTO 205
      VO0=EA(4)
      GOTO 210
 205  VO0=EA(4)+0.18D3
      IF(VO0.GE.0.36D3) VO0=VO0-0.36D3
 210  VO0=VO0*PI180
      DWM=2.0D0*PI2
      TA=TMIN
      TB=TMAX
      DT=DT0
C
 220  CONTINUE
C >      DO 230 ETT=TA,TB,DT
      IETTU=(TB+1.0D-11-TA)/DT
      DO 230 IETT=0,IETTU
      ETT=TA+IETT*DT
c      CALL PLEPH(ETT,NTARG,NCENT,RRD)
      T20=(ETT-2.451545D6)/3.6525D4
      CALL EARTHrv(T20,XZ,YZ,ZZ,VXZ,VYZ,VZZ,WEL)
      IF(DABS(XZ).LT.1.0D-12) THEN
      WL=PI/2.0D0
      IF(YZ.LT.0.0D0) WL=WL+PI
      ELSE
      WL=DATAN(YZ/XZ)
      WL=DABS(WL)
      IF(XZ.LT.0.0D0.AND.YZ.GE.0.0D0) WL=PI-WL
      IF(XZ.LT.0.0D0.AND.YZ.LT.0.0D0) WL=PI+WL
      IF(XZ.GE.0.0D0.AND.YZ.LT.0.0D0) WL=PI2-WL
      END IF
C      WL=WEL(7)+WEL(4)
C      IF(WL.GE.PI2) WL=WL-PI2
      IF(DABS(WL-VO0).GT.DWM) GOTO 230
Cx      EZ(1)=WEL(1)
Cx      EZ(2)=WEL(2)
Cx      EZ(3)=WEL(4)-WEL(5)
Cx      EZ(4)=WEL(5)
Cx      EZ(5)=WEL(6)
Cx      FZM=WEL(7)
Cx      EZ(1)=EZ(1)/(1.0D0-EZ(2))
Cx      FZM=EZ(6)
      RZM=DSQRT(XZ*XZ+YZ*YZ+ZZ*ZZ)
      VYR=VYZ*CEPS-VZZ*SEPS
      VZR=VYZ*SEPS+VZZ*CEPS
      R2(1)=WL/PI180
      R2(6)=ETT
      DWM=DABS(WL-VO0)
 230  CONTINUE
      IF(DT.LT.DT05) GOTO 250
      TA=R2(6)-DT
      TB=R2(6)+DT
      DT=0.1D-1*DT
      GOTO 220
C
 250  CONTINUE
Cx ????      IF(EZ(3).LT.0.0D0) EZ(3)=EZ(3)+PI2
Cx      RZM=EZ(1)*(1.0D0-EZ(2)*EZ(2))/(1.0D0+EZ(2)*DCOS(FZM))
      CFB=(EA(1)*(1.0D0+EA(2))/RZM-1.0D0)/EA(2)
      DD2=0.999D99
      IF(DABS(CFB).GT.1.000001D0) GOTO 300
      IF(CFB.LT.-1.0D0) CFB=-1.0D0
      IF(CFB.GT.1.0D0) CFB=1.0D0
      IF(DABS(CFB).LT.1.0D-12) GOTO 260
      SFB=DSQRT(1.0D0-CFB*CFB)
      FB=DATAN(SFB/CFB)
      IF(FB.LT.0.0D0) FB=FB+PI
      GOTO 270
 260  FB=PI05
 270  FB=0.36D3-FB/PI180
      EB2(3)=0.36D3-FB
      IF(EA(3).GT.0.18D3) EB2(3)=0.54D3-FB
      CALL DCRIT(EA,EB2,DD2)
      R2(1)=R2(1)-0.18D3
      IF(R2(1).LT.0.0D0) R2(1)=R2(1)+0.36D3
      CALL RADIANT(EB2,FB,VXZ,VYR,VZR,R2)
C
 300  CONTINUE
      RETURN
C
      END
C ---------------------------------------------------------------------
      SUBROUTINE AMETH(TMIN,TMAX,EA,R1,EB1,DD1,R2,EB2,DD2)
C ---------------------------------------------------------------------
C      - (A) ROTATION AROUND THE LINE OF APSIDES
C     (SVOREN J., NESLUSAN L., PORUBCAN V.: 1993, CONTRIB. ASTRON. OBS.
C     SKALNATE PLESO 23, 23.)
C
C     "TMIN", "TMAX" - BEGINNING AND END OF YEAR OF INVESTIGATION OF THE
C     POTENTIAL SHOWER GIVEN IN JULIAN CENTURIES FROM 1950.0. NEXT
C     INPUT (FIELD OF VALUES "EA(I)", I = 1, 2,..., 6) AND OUTPUT (FIEL-
C     DS OF VALUES "R1(I)", "R2(I)", "EB1(I)", AND "EB2(I)", I = 1, 2,
C     3,..., 6 OR 5, RESPECTIVELY), PLEASE, SEE THE MAIN PROGRAM. AS
C     THE OUTPUT, THERE ARE MOREOVER THE VALUES "DD1" AND "DD2" OF D-
C     DISCRIMINANT CORRESPONDING TO THE POST- AND PRE-PERIHELION ARC
C     OF THE PARENT BODY ORBIT, RESPECTIVELY.
C
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER  I,IM,IYR,J
      DIMENSION EA(6),EB1(5),EB2(5),R1(6),R2(6),WEL(7)
Cx      ,EZ(6)
C
      PI=4.0D0*DATAN(1.0D0)
      PI05=PI/2.0D0
      PI2=2.0D0*PI
      PI180=PI/0.18D3
      EPS=2.3439291D1
      CEPS=DCOS(PI180*EPS)
      SEPS=DSIN(PI180*EPS)
      NTARG=3
      NCENT=11
C
      DD1=0.999D99
      DD2=0.999D99
      IF(EA(2).LT.1.0D-12) GOTO 5
      IF(EA(2).GT.0.999D0) GOTO 10
      IF(EA(1).LT.1.1D0) GOTO 10
      IF(EA(1)*(1.0D0+EA(2))/(1.0D0-EA(2)).GT.0.9D0) GOTO 10
   5  CONTINUE
      RETURN
C
  10  EB1(1)=EA(1)
      EB1(2)=EA(2)
      EB2(1)=EA(1)
      EB2(2)=EA(2)
      CSO=DCOS(PI180*EA(3))
      SSO=DSIN(PI180*EA(3))
      CHO=DCOS(PI180*EA(4))
      SHO=DSIN(PI180*EA(4))
      CSK=DCOS(PI180*EA(5))
      SSK=DSIN(PI180*EA(5))
      XQ=CSO*CHO-SSO*CSK*SHO
      YQ=CSO*SHO+SSO*CSK*CHO
      ZQ=SSO*SSK
      RQ=DSQRT(XQ*XQ+YQ*YQ+ZQ*ZQ)
      IF(DABS(RQ-1.0D0).GT.1.0D-12) WRITE(*,11)
 11   FORMAT(/,'  ERROR IN SUBROUTINE `AMETH` (RQ.NE.1).')
      IF(DABS(RQ-1.0D0).GT.1.0D-12) STOP
C
      DT0=(TMAX-TMIN)/3.6525D2
      DT05=DT0/2.0D0
C
C >      DO 13 ETT=TMIN,TMAX,DT05
      IETTU=(TMAX+1.0D-11-TMIN)/DT05
      DO 13 IETT=0,IETTUD
      ETT=TMIN+IETT*DT05
c      CALL PLEPH(ETT,NTARG,NCENT,RRD)
      T20=(ETT-2.451545D6)/3.6525D4
      CALL EARTHrv(T20,XZ,YZ,ZZ,VXZ,VYZ,VZZ,WEL)
Cx      XZ=RRD(1)
Cx      YZ=RRD(2)*CEPS+RRD(3)*SEPS
Cx      ZZ=-RRD(2)*SEPS+RRD(3)*CEPS
Cx      VXZ=RRD(4)
Cx      VYZ=RRD(5)*CEPS+RRD(6)*SEPS
Cx      VZZ=-RRD(5)*SEPS+RRD(6)*CEPS
Cx      CALL ELEMS(ETT,XZ,YZ,VZ,VXZ,VYZ,VZZ,EZ)
Cx      EZ(1)=EZ(1)/(1.0D0-EZ(2))
Cx      FZM=EZ(6)
Cx      RZ0=WEL(1)*(1.0D0-WEL(2)*WEL(2))/(1.0D0+WEL(2)*DCOS(WEL(7)))
      RZ0=DSQRT(XZ*XZ+YZ*YZ+ZZ*ZZ)
      CF0=(EA(1)*(1.0D0+EA(2))/RZ0-1.0D0)/EA(2)
      IF(DABS(CF0).LE.1.0D0) GOTO 15
 13   CONTINUE
      GOTO 5
C
 15   CONTINUE
      IF(DABS(CF0).LT.1.0D-9) GOTO 17
      IF(DABS(CF0).GT.1.000000001D0) GOTO 200
      IF(CF0.GT.1.0D0) CF0=1.0D0
      IF(CF0.LT.-1.0D0) CF0=-1.0D0
      SF0=DSQRT(1.0D0-CF0*CF0)
      F0=DATAN(SF0/CF0)
      IF(F0.LT.0.0D0) F0=F0+PI
      GOTO 20
  17  F0=PI05
  20  EB1(3)=PI-F0
      IF(EA(3).GT.0.18D3) EB1(3)=PI2-F0
      CSOB=DCOS(EB1(3))
      SSOB=DSIN(EB1(3))
      EB1(3)=EB1(3)/PI180
      QKOB=SSOB*SSOB+CSOB*CSOB
      IF(DABS(QKOB-1.0D0).GT.1.0D-12) WRITE(*,30)
  30  FORMAT(/,'  ERROR IN SUBROUTINE `AMETH` (QKOB.NE.1).')
      IF(DABS(QKOB-1.0D0).GT.1.0D-12) STOP
      IF(DABS(SSOB).LT.1.0D-9.AND.DABS(SSO*SSK).GT.1.0D-9) GOTO 200
      IF(DABS(SSOB).LT.1.0D-9) GOTO 40
      SIB=SSO*SSK/SSOB
      IF(DABS(SIB).GT.1.0000000001D0) GOTO 200
      GOTO 60
  40  SIB=SSK
  60  CONTINUE
      IF(SIB.GT.1.0D0) SIB=1.0D0
      IF(SIB.LT.-1.0D0) SIB=-1.0D0
      CIB=DSQRT(1.0D0-SIB*SIB)
      JIB=1
  65  CONTINUE
      IF(DABS(CIB).LT.1.0D-12) EB1(5)=PI05
      IF(DABS(CIB).LT.1.0D-12) GOTO 70
      EB1(5)=DATAN(SIB/CIB)
      IF(CIB.LT.0.0D0) EB1(5)=PI+EB1(5)
  70  CONTINUE
      EB1(5)=EB1(5)/PI180
      WEVO=CSOB*CSOB+SSOB*SSOB*CIB*CIB
C      IF(DABS(CSOB).LT.1.0D-12.AND.DABS(CIB).LT.1.0D-12) DD1=0.999D99
      IF(DABS(CSOB).LT.1.0D-12.AND.DABS(CIB).LT.1.0D-12) GOTO 200
      QQS=(YQ*CSOB-XQ*SSOB*CIB)/WEVO
      QQC=(XQ*CSOB+YQ*SSOB*CIB)/WEVO
      QKONTR=QQS*QQS+QQC*QQC
      IF(DABS(QKONTR-1.0D0).GT.1.0D-12) WRITE(*,80)
  80  FORMAT(/,'  ERROR IN SUBROUTINE `AMETH` ((SIN**2(EB(4)) + COS**(EB
     *(4))).GT.1).')
      IF(DABS(QKONTR-1.0D0).GT.1.0D-12) STOP     
      IF(DABS(QQC).LT.1.0D-12) EB1(4)=PI05
      IF(DABS(QQC).LT.1.0D-12) GOTO 90
      EB1(4)=DATAN(QQS/QQC)
      EB1(4)=DABS(EB1(4))
      IF(QQC.GE.0.0D0.AND.QQS.LT.0.0D0) EB1(4)=PI2-EB1(4)
      IF(QQC.LT.0.0D0.AND.QQS.GE.0.0D0) EB1(4)=PI-EB1(4)
      IF(QQC.LT.0.0D0.AND.QQS.LT.0.0D0) EB1(4)=PI+EB1(4)
  90  CONTINUE
      EB1(4)=EB1(4)/PI180
      CALL DCRIT(EA,EB1,DD1)
      IF(JIB.EQ.2) GOTO 110
      SOBO=EB1(3)
      HOBO=EB1(4)
      SKBO=EB1(5)
      DDO=DD1
 100  JIB=2
      CIB=-CIB
      GOTO 65
 110  CONTINUE
      TA=TMIN
      TB=TMAX
      DT=DT0
      IF(DD1.LT.DDO) GOTO 120
      EB1(3)=SOBO
      EB1(4)=HOBO
      EB1(5)=SKBO
      DD1=DDO
C
 120  DWM=2.0D0*PI2
      VO0=EB1(4)*PI180+PI
      IF(EB1(3).GT.0.18D3) VO0=EB1(4)*PI180
      IF(VO0.GE.PI2) VO0=VO0-PI2
C >      DO 130 ETT=TA,TB,DT
      IETTU=(TB+1.0D-11-TA)/DT
      DO 130 IETT=0,IETTU
      ETT=TA+IETT*DT
c      CALL PLEPH(ETT,NTARG,NCENT,RRD)
      T20=(ETT-2.451545D6)/3.6525D4
      CALL EARTHrv(T20,XZ,YZ,ZZ,VXZ,VYZ,VZZ,WEL)
      IF(DABS(XZ).LT.1.0D-12) THEN
      WL=PI/2.0D0
      IF(YZ.LT.0.0D0) WL=WL+PI
      ELSE
      WL=DATAN(YZ/XZ)
      WL=DABS(WL)
      IF(XZ.LT.0.0D0.AND.YZ.GE.0.0D0) WL=PI-WL
      IF(XZ.LT.0.0D0.AND.YZ.LT.0.0D0) WL=PI+WL
      IF(XZ.GE.0.0D0.AND.YZ.LT.0.0D0) WL=PI2-WL
      END IF
C      WL=WEL(7)+WEL(4)
C      IF(WL.GE.PI2) WL=WL-PI2
      IF(DABS(WL-VO0).GT.DWM) GOTO 130
Cx      EZ(1)=WEL(1)
Cx      EZ(2)=WEL(2)
Cx      EZ(3)=WEL(4)-WEL(5)
Cx      EZ(4)=WEL(5)
Cx      EZ(5)=WEL(6)
Cx      FZM=WEL(7)
Cx      EZ(1)=EZ(1)/(1.0D0-EZ(2))
Cx      FZM=EZ(6)
      RZ=DSQRT(XZ*XZ+YZ*YZ+ZZ*ZZ)
      VYR=VYZ*CEPS-VZZ*SEPS
      VZR=VYZ*SEPS+VZZ*CEPS
      R1(1)=WL/PI180
      R1(6)=ETT
      DWM=DABS(WL-VO0)
 130  CONTINUE
      IF(DT.LT.DT05) GOTO 150
      TA=R1(6)-DT
      TB=R1(6)+DT
      DT=0.1D-1*DT
      GOTO 120
C
 150  CONTINUE
Cx ????      IF(EZ(3).LT.0.0D0) EZ(3)=EZ(3)+PI2
Cx      RZ=EZ(1)*(1.0D0-EZ(2)*EZ(2))/(1.0D0+EZ(2)*DCOS(FZM))      
      IF(DABS(RZ-RZ0).LT.1.0D-4) GOTO 160
      RZ0=RZ
      CF0=(EA(1)*(1.0D0+EA(2))/RZ0-1.0D0)/EA(2)
      IF(DABS(CF0).LE.1.0D0) GOTO 15
      GOTO 200
 160  CONTINUE
      IF(DABS(CF0).LT.1.0D-9) GOTO 170
      IF(DABS(CF0).GT.1.000000001D0) GOTO 200
      IF(CF0.GT.1.0D0) CF0=1.0D0
      IF(CF0.LT.-1.0D0) CF0=-1.0D0
      SF0=DSQRT(1.0D0-CF0*CF0)
      F0=DATAN(SF0/CF0)
      IF(F0.LT.0.0D0) F0=F0+PI
      GOTO 180
 170  F0=PI05
 180  EB1(3)=(PI-F0)/PI180
      IF(EA(3).GT.0.18D3) EB1(3)=(PI2-F0)/PI180
      FB=F0/PI180
      R1(1)=R1(1)-0.18D3
      IF(R1(1).LT.0.0D0) R1(1)=R1(1)+0.36D3
      CALL RADIANT(EB1,FB,VXZ,VYR,VZR,R1)
C
C   2-ND ARC
C
 200  CONTINUE
      IF(DABS(CF0).LT.1.0D-9) GOTO 217
      IF(DABS(CF0).GT.1.000000001D0) GOTO 400
      IF(CF0.GT.1.0D0) CF0=1.0D0
      IF(CF0.LT.-1.0D0) CF0=-1.0D0
      SF0=DSQRT(1.0D0-CF0*CF0)
      F0=DATAN(SF0/CF0)
      IF(F0.LT.0.0D0) F0=F0+PI
      GOTO 220
 217  F0=PI05
 220  F0=PI2-F0
      EB2(3)=PI2-F0
      IF(EA(3).GT.0.18D3) EB2(3)=PI2+PI-F0
      CSOB=DCOS(EB2(3))
      SSOB=DSIN(EB2(3))
      EB2(3)=EB2(3)/PI180
      QKOB=SSOB*SSOB+CSOB*CSOB
      IF(DABS(QKOB-1.0D0).GT.1.0D-12) WRITE(*,230)
 230  FORMAT(/,'  ERROR IN SUBROUTINE `AMETH` (QKOB.NE.1).')
      IF(DABS(QKOB-1.0D0).GT.1.0D-12) STOP
      IF(DABS(SSOB).LT.1.0D-12.AND.DABS(SSO).GT.1.0D-9) GOTO 400
      IF(DABS(SSOB).LT.1.0D-12) GOTO 245
      SIB=SSO*SSK/SSOB
      IF(DABS(SIB).GT.1.0000000001D0) GOTO 400
      GOTO 260
 245  SIB=SSK
 260  CONTINUE
      IF(SIB.GT.1.0D0) SIB=1.0D0
      IF(SIB.LT.-1.0D0) SIB=-1.0D0
      CIB=DSQRT(1.0D0-SIB*SIB)
      JIB=1
 265  CONTINUE
      IF(DABS(CIB).LT.1.0D-12) EB2(5)=PI05
      IF(DABS(CIB).LT.1.0D-12) GOTO 270
      EB2(5)=DATAN(SIB/CIB)
      IF(CIB.LT.0.0D0) EB2(5)=PI+EB2(5)
 270  CONTINUE
      EB2(5)=EB2(5)/PI180
      WEVO=CSOB*CSOB+SSOB*SSOB*CIB*CIB
C      IF(DABS(CSOB).LT.1.0D-12.AND.DABS(CIB).LT.1.0D-12) DD2=0.999D99
      IF(DABS(CSOB).LT.1.0D-12.AND.DABS(CIB).LT.1.0D-12) GOTO 400
      QQS=(YQ*CSOB-XQ*SSOB*CIB)/WEVO
      QQC=(XQ*CSOB+YQ*SSOB*CIB)/WEVO
      QKONTR=QQS*QQS+QQC*QQC
      IF(DABS(QKONTR-1.0D0).GT.1.0D-12) WRITE(*,280)
 280  FORMAT(/,'  ERROR IN SUBROUTINE `AMETH` ((SIN**2(EB(4)) + COS**(EB
     *(4))).GT.1).')
      IF(DABS(QKONTR-1.0D0).GT.1.0D-12) STOP     
      IF(DABS(QQC).LT.1.0D-12) EB2(4)=PI05
      IF(DABS(QQC).LT.1.0D-12) GOTO 290
      EB2(4)=DATAN(QQS/QQC)
      EB2(4)=DABS(EB2(4))
      IF(QQC.GE.0.0D0.AND.QQS.LT.0.0D0) EB2(4)=PI2-EB2(4)
      IF(QQC.LT.0.0D0.AND.QQS.GE.0.0D0) EB2(4)=PI-EB2(4)
      IF(QQC.LT.0.0D0.AND.QQS.LT.0.0D0) EB2(4)=PI+EB2(4)
 290  CONTINUE
      EB2(4)=EB2(4)/PI180
      CALL DCRIT(EA,EB2,DD2)
      IF(JIB.EQ.2) GOTO 310
      SOBO=EB2(3)
      HOBO=EB2(4)
      SKBO=EB2(5)
      DDO=DD2
 300  JIB=2
      CIB=-CIB
      GOTO 265
 310  CONTINUE
      TA=TMIN
      TB=TMAX
      DT=DT0
      IF(DD2.LT.DDO) GOTO 320
      EB2(3)=SOBO
      EB2(4)=HOBO
      EB2(5)=SKBO
      DD2=DDO
C
 320  DWM=2.0D0*PI2
      VO0=EB2(4)*PI180
      IF(EB2(3).GT.0.18D3) VO0=EB2(4)*PI180+PI
      IF(VO0.GE.PI2) VO0=VO0-PI2
C >      DO 330 ETT=TA,TB,DT
      IETTU=(TB+1.0D-11-TA)/DT
      DO 330 IETT=0,IETTU
      ETT=TA+IETT*DT
c      CALL PLEPH(ETT,NTARG,NCENT,RRD)
      T20=(ETT-2.451545D6)/3.6525D4
      CALL EARTHrv(T20,XZ,YZ,ZZ,VXZ,VYZ,VZZ,WEL)
      IF(DABS(XZ).LT.1.0D-12) THEN
      WL=PI/2.0D0
      IF(YZ.LT.0.0D0) WL=WL+PI
      ELSE
      WL=DATAN(YZ/XZ)
      WL=DABS(WL)
      IF(XZ.LT.0.0D0.AND.YZ.GE.0.0D0) WL=PI-WL
      IF(XZ.LT.0.0D0.AND.YZ.LT.0.0D0) WL=PI+WL
      IF(XZ.GE.0.0D0.AND.YZ.LT.0.0D0) WL=PI2-WL
      END IF
C      WL=WEL(7)+WEL(4)
C      IF(WL.GE.PI2) WL=WL-PI2
      IF(DABS(WL-VO0).GT.DWM) GOTO 330
Cx      EZ(1)=WEL(1)
Cx      EZ(2)=WEL(2)
Cx      EZ(3)=WEL(4)-WEL(5)
Cx      EZ(4)=WEL(5)
Cx      EZ(5)=WEL(6)
Cx      FZM=WEL(7)
Cx      EZ(1)=EZ(1)/(1.0D0-EZ(2))
Cx      FZM=EZ(6)
      RZ=DSQRT(XZ*XZ+YZ*YZ+ZZ*ZZ)
      VYR=VYZ*CEPS-VZZ*SEPS
      VZR=VYZ*SEPS+VZZ*CEPS
      R2(1)=WL/PI180
      R2(6)=ETT
      DWM=DABS(WL-VO0)
 330  CONTINUE
      IF(DT.LT.DT05) GOTO 350
      TA=R2(6)-DT
      TB=R2(6)+DT
      DT=0.1D-1*DT
      GOTO 320
C
 350  CONTINUE
Cx ????      IF(EZ(3).LT.0.0D0) EZ(3)=EZ(3)+PI2
Cx      RZ=EZ(1)*(1.0D0-EZ(2)*EZ(2))/(1.0D0+EZ(2)*DCOS(FZM))      
      IF(DABS(RZ-RZ0).LT.1.0D-4) GOTO 360
      RZ0=RZ
      CF0=(EA(1)*(1.0D0+EA(2))/RZ0-1.0D0)/EA(2)
      IF(DABS(CF0).LE.1.0D0) GOTO 200
      GOTO 400
 360  CONTINUE
      IF(DABS(CF0).LT.1.0D-9) GOTO 370
      IF(DABS(CF0).GT.1.000000001D0) GOTO 400
      IF(CF0.GT.1.0D0) CF0=1.0D0
      IF(CF0.LT.-1.0D0) CF0=-1.0D0
      SF0=DSQRT(1.0D0-CF0*CF0)
      F0=DATAN(SF0/CF0)
      IF(F0.LT.0.0D0) F0=F0+PI
      GOTO 380
 370  F0=PI05
 380  F0=PI2-F0
      EB2(3)=(PI2-F0)/PI180
      IF(EA(3).GT.0.18D3) EB2(3)=(PI2+PI-F0)/PI180
      FB=F0/PI180
      R2(1)=R2(1)-0.18D3
      IF(R2(1).LT.0.0D0) R2(1)=R2(1)+0.36D3
      CALL RADIANT(EB2,FB,VXZ,VYR,VZR,R2)
C
 400  CONTINUE
      RETURN
C
      END
C ---------------------------------------------------------------------
      SUBROUTINE HMETH(TMIN,TMAX,EA,R1,EB1,DD1,R2,EB2,DD2)
C ---------------------------------------------------------------------
C      - (H) OMEGA-ADJUSTMENT
C     (HASEGAWA I.: 1990, PUBL. ASTRON. SOC. JAPAN 42, 175.)
C
C     "TMIN", "TMAX" - BEGINNING AND END OF YEAR OF INVESTIGATION OF THE
C     POTENTIAL SHOWER GIVEN IN JULIAN CENTURIES FROM 1950.0. NEXT
C     INPUT (FIELD OF VALUES "EA(I)", I = 1, 2,..., 6) AND OUTPUT (FIEL-
C     DS OF VALUES "R1(I)", "R2(I)", "EB1(I)", AND "EB2(I)", I = 1, 2,
C     3,..., 6 OR 5, RESPECTIVELY), PLEASE, SEE THE MAIN PROGRAM. AS
C     THE OUTPUT, THERE ARE MOREOVER THE VALUES "DD1" AND "DD2" OF D-
C     DISCRIMINANT CORRESPONDING TO THE POST- AND PRE-PERIHELION ARC
C     OF THE PARENT BODY ORBIT, RESPECTIVELY.
C
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER  I,IM,IYR,J
      DIMENSION EA(6),EB1(5),EB2(5),R1(6),R2(6),WEL(7)
Cx      ,EZ(6)
C
      PI=4.0D0*DATAN(1.0D0)
      PI180=PI/0.18D3
      PI05=PI/2.0D0
      PI2=2.0D0*PI
      EPS=2.3439291D1
      CEPS=DCOS(PI180*EPS)
      SEPS=DSIN(PI180*EPS)
      NTARG=3
      NCENT=11
C
      GOTO 10
   5  DD1=0.999D99
      DD2=0.999D99
      RETURN
  10  EB1(1)=EA(1)
      EB1(2)=EA(2)
      EB2(1)=EA(1)
      EB2(2)=EA(2)
      CSK=DCOS(PI180*EA(5))
      SSK=DSIN(PI180*EA(5))
C
      DT0=(TMAX-TMIN)/3.6525D2
      DT05=DT0/2.0D0
C
C >      DO 13 ETT=TMIN,TMAX,DT05
      IETTU=(TMAX+1.0D-11-TMIN)/DT05
      DO 13 IETT=0,IETTU
      ETT=TMIN+IETT*DT05
c      CALL PLEPH(ETT,NTARG,NCENT,RRD)
      T20=(ETT-2.451545D6)/3.6525D4
      CALL EARTHrv(T20,XZ,YZ,ZZ,VXZ,VYZ,VZZ,WEL)
      RZ0=DSQRT(XZ*XZ+YZ*YZ+ZZ*ZZ)
      CF0=(EA(1)*(1.0D0+EA(2))/RZ0-1.0D0)/EA(2)
      IF(DABS(CF0).LE.1.0D0) GOTO 15
 13   CONTINUE
      GOTO 5
C
 15   CONTINUE
      IF(DABS(CF0).LT.1.0D-9) GOTO 17
      SF0=DSQRT(1.0D0-CF0*CF0)
      F0=DATAN(SF0/CF0)
      IF(F0.LT.0.0D0) F0=F0+PI
      GOTO 20
  17  F0=PI05
  20  U=EA(3)*PI180+F0
      IF(U.GE.PI2) U=U-PI2
      CU=DCOS(U)
      SU=DSIN(U)
      IF(DABS(CU).LT.1.0D-12) HLVO=PI05
      IF(DABS(CU).LT.1.0D-12) GOTO 30
      HLVO=DATAN(SU*CSK/CU)
      HLVO=DABS(HLVO)
  30  CONTINUE
      IF(SU*CSK.GE.0.0D0.AND.CU.LT.0.0D0) HLVO=PI-HLVO
      IF(SU*CSK.LT.0.0D0.AND.CU.GE.0.0D0) HLVO=PI2-HLVO
      IF(SU*CSK.LT.0.0D0.AND.CU.LT.0.0D0) HLVO=PI+HLVO
      HL=HLVO+EA(4)*PI180
      IF(HL.GE.PI2) HL=HL-PI2
C
      TA=TMIN
      TB=TMAX
      DT=DT0
      DWM=2.0D0*PI2
C
  40  CONTINUE
C >      DO 50 ETT=TA,TB,DT
      IETTU=(TB+1.0D-11-TA)/DT
      DO 50 IETT=0,IETTU
      ETT=TA+IETT*DT
c      CALL PLEPH(ETT,NTARG,NCENT,RRD)
      T20=(ETT-2.451545D6)/3.6525D4
      CALL EARTHrv(T20,XZ,YZ,ZZ,VXZ,VYZ,VZZ,WEL)
      IF(DABS(XZ).LT.1.0D-12) THEN
      WL=PI/2.0D0
      IF(YZ.LT.0.0D0) WL=WL+PI
      ELSE
      WL=DATAN(YZ/XZ)
      WL=DABS(WL)
      IF(XZ.LT.0.0D0.AND.YZ.GE.0.0D0) WL=PI-WL
      IF(XZ.LT.0.0D0.AND.YZ.LT.0.0D0) WL=PI+WL
      IF(XZ.GE.0.0D0.AND.YZ.LT.0.0D0) WL=PI2-WL
      END IF
C      WL=WEL(7)+WEL(4)
C      IF(WL.GE.PI2) WL=WL-PI2
      IF(DABS(WL-HL).GT.DWM) GOTO 50
Cx      EZ(1)=WEL(1)
Cx      EZ(2)=WEL(2)
Cx      EZ(3)=WEL(4)-WEL(5)
Cx      EZ(4)=WEL(5)
Cx      EZ(5)=WEL(6)
Cx      FZM=WEL(7)
Cx      EZ(1)=EZ(1)/(1.0D0-EZ(2))
Cx      FZM=EZ(6)
      RZ=DSQRT(XZ*XZ+YZ*YZ+ZZ*ZZ)
      VYR=VYZ*CEPS-VZZ*SEPS
      VZR=VYZ*SEPS+VZZ*CEPS
      R1(1)=WL/PI180
      R1(6)=ETT
      DWM=DABS(WL-HL)
  50  CONTINUE
      IF(DT.LT.DT05) GOTO 60
      TA=R1(6)-DT
      TB=R1(6)+DT
      DT=0.1D-1*DT
      GOTO 40
C
  60  CONTINUE
Cx ????      IF(EZ(3).LT.0.0D0) EZ(3)=EZ(3)+PI2
Cx      RZ=EZ(1)*(1.0D0-EZ(2)*EZ(2))/(1.0D0+EZ(2)*DCOS(FZM))
      CF0=(EA(1)*(1.0D0+EA(2))/RZ-1.0D0)/EA(2)
      IF(DABS(CF0).GT.1.000000001D0) GOTO 5
      IF(DABS(RZ-RZ0).LT.1.0D-4) GOTO 65
      RZ0=RZ
      GOTO 15
  65  CONTINUE
      IF(DABS(CF0).LT.1.0D-9) GOTO 67
      IF(CF0.GT.1.0D0) CF0=1.0D0
      IF(CF0.LT.-1.0D0) CF0=-1.0D0
      SF0=DSQRT(1.0D0-CF0*CF0)
      F0=DATAN(SF0/CF0)
      IF(F0.LT.0.0D0) F0=F0+PI
      GOTO 70
  67  F0=PI05
C
  70  CONTINUE
      EB1(3)=PI2-F0
      EB1(4)=HL
      U=U/PI180
      IF(U.LE.9.0D1.OR.U.GT.2.7D2) GOTO 80
      EB1(3)=EB1(3)-PI
      IF(EB1(3).LT.0.0D0) EB1(3)=EB1(3)+PI2
      EB1(4)=EB1(4)+PI
      IF(EB1(4).GE.PI2) EB1(4)=EB1(4)-PI2
  80  EB1(3)=EB1(3)/PI180
      EB1(4)=EB1(4)/PI180
C
      SIB=SSK*DCOS(HLVO)
      IF(DABS(SIB).GT.1.0000000001D0) GOTO 90
      GOTO 100
  90  CONTINUE
      WRITE(*,95)
  95  FORMAT(/,'  ERROR IN SUBROUTINE `HMETH` (ABS(SIN(EB(5))).GT.1).')
      STOP
 100  CONTINUE
      IF(SIB.GT.1.0D0) SIB=1.0D0
      IF(SIB.LT.-1.0D0) SIB=-1.0D0
      CIB=DSQRT(1.0D0-SIB*SIB)
      IF(DABS(CIB).LT.1.0D-12) EB1(5)=PI05
      IF(DABS(CIB).LT.1.0D-12) GOTO 110
      EB1(5)=DATAN(SIB/CIB)
      IF(EB1(5).LT.0.0D0) EB1(5)=PI+EB1(5)
 110  EB1(5)=EB1(5)/PI180
      IF(EA(5).GT.9.0D1.AND.EB1(5).LT.9.0D1) EB1(5)=1.8D2-EB1(5)
      IF(EA(5).LT.9.0D1.AND.EB1(5).GT.9.0D1) EB1(5)=1.8D2-EB1(5)
C
      CALL DCRIT(EA,EB1,DD1)
      FB=F0/PI180
      R1(1)=R1(1)-0.18D3
      IF(R1(1).LT.0.0D0) R1(1)=R1(1)+0.36D3
      CALL RADIANT(EB1,FB,VXZ,VYR,VZR,R1)
C
C   2-ND ARC
C
 200  CONTINUE
      IF(DABS(CF0).LT.1.0D-9) GOTO 217
      SF0=DSQRT(1.0D0-CF0*CF0)
      F0=DATAN(SF0/CF0)
      IF(F0.LT.0.0D0) F0=F0+PI
      F0=PI2-F0
      GOTO 220
 217  F0=PI2-PI05
 220  U=EA(3)*PI180+F0
      IF(U.GE.PI2) U=U-PI2
      CU=DCOS(U)
      SU=DSIN(U)
      IF(DABS(CU).LT.1.0D-12) HLVO=PI05
      IF(DABS(CU).LT.1.0D-12) GOTO 230
      HLVO=DATAN(SU*CSK/CU)
      HLVO=DABS(HLVO)
 230  CONTINUE
      IF(SU*CSK.GE.0.0D0.AND.CU.LT.0.0D0) HLVO=PI-HLVO
      IF(SU*CSK.LT.0.0D0.AND.CU.GE.0.0D0) HLVO=PI2-HLVO
      IF(SU*CSK.LT.0.0D0.AND.CU.LT.0.0D0) HLVO=PI+HLVO
      HL=HLVO+EA(4)*PI180
      IF(HL.GE.PI2) HL=HL-PI2
C
      TA=TMIN
      TB=TMAX
      DT=DT0
      DWM=2.0D0*PI2
C
 240  CONTINUE
C >      DO 250 ETT=TA,TB,DT
      IETTU=(TB+1.0D-11-TA)/DT
      DO 250 IETT=0,IETTU
      ETT=TA+IETT*DT
c      CALL PLEPH(ETT,NTARG,NCENT,RRD)
      T20=(ETT-2.451545D6)/3.6525D4
      CALL EARTHrv(T20,XZ,YZ,ZZ,VXZ,VYZ,VZZ,WEL)
      IF(DABS(XZ).LT.1.0D-12) THEN
      WL=PI/2.0D0
      IF(YZ.LT.0.0D0) WL=WL+PI
      ELSE
      WL=DATAN(YZ/XZ)
      WL=DABS(WL)
      IF(XZ.LT.0.0D0.AND.YZ.GE.0.0D0) WL=PI-WL
      IF(XZ.LT.0.0D0.AND.YZ.LT.0.0D0) WL=PI+WL
      IF(XZ.GE.0.0D0.AND.YZ.LT.0.0D0) WL=PI2-WL
      END IF
C      WL=WEL(7)+WEL(4)
C      IF(WL.GE.PI2) WL=WL-PI2
      IF(DABS(WL-HL).GT.DWM) GOTO 250
Cx      EZ(1)=WEL(1)
Cx      EZ(2)=WEL(2)
Cx      EZ(3)=WEL(4)-WEL(5)
Cx      EZ(4)=WEL(5)
Cx      EZ(5)=WEL(6)
Cx      FZM=WEL(7)
Cx      EZ(1)=EZ(1)/(1.0D0-EZ(2))
Cx      FZM=EZ(6)
      RZ=DSQRT(XZ*XZ+YZ*YZ+ZZ*ZZ)
      VYR=VYZ*CEPS-VZZ*SEPS
      VZR=VYZ*SEPS+VZZ*CEPS
      R2(1)=WL/PI180
      R2(6)=ETT
      DWM=DABS(WL-HL)
 250  CONTINUE
      IF(DT.LT.DT05) GOTO 260
      TA=R2(6)-DT
      TB=R2(6)+DT
      DT=0.1D-1*DT
      GOTO 240
C
 260  CONTINUE
Cx ?????      IF(EZ(3).LT.0.0D0) EZ(3)=EZ(3)+PI2
Cx      RZ=EZ(1)*(1.0D0-EZ(2)*EZ(2))/(1.0D0+EZ(2)*DCOS(FZM))
      CF0=(EA(1)*(1.0D0+EA(2))/RZ-1.0D0)/EA(2)
      IF(DABS(CF0).GT.1.000000001D0) GOTO 5
      IF(DABS(RZ-RZ0).LT.1.0D-4) GOTO 265
      RZ0=RZ
      GOTO 200
 265  CONTINUE
      IF(DABS(CF0).LT.1.0D-9) GOTO 267
      IF(CF0.GT.1.0D0) CF0=1.0D0
      IF(CF0.LT.-1.0D0) CF0=-1.0D0
      SF0=DSQRT(1.0D0-CF0*CF0)
      F0=DATAN(SF0/CF0)
      IF(F0.LT.0.0D0) F0=F0+PI
      F0=PI2-F0
      GOTO 270
 267  F0=PI2-PI05
C
 270  CONTINUE
      EB2(3)=PI2-F0
      EB2(4)=HL
      U=U/PI180
      IF(U.LE.9.0D1.OR.U.GT.2.7D2) GOTO 280
      EB2(3)=EB2(3)-PI
      IF(EB2(3).LT.0.0D0) EB2(3)=EB2(3)+PI2
      EB2(4)=EB2(4)+PI
      IF(EB2(4).GE.PI2) EB2(4)=EB2(4)-PI2
 280  EB2(3)=EB2(3)/PI180
      EB2(4)=EB2(4)/PI180
C
      SIB=SSK*DCOS(HLVO)
      IF(DABS(SIB).GT.1.0000000001D0) GOTO 290
      GOTO 300
 290  CONTINUE
      WRITE(*,95)
      STOP
 300  CONTINUE
      IF(SIB.GT.1.0D0) SIB=1.0D0
      IF(SIB.LT.-1.0D0) SIB=-1.0D0
      CIB=DSQRT(1.0D0-SIB*SIB)
      IF(DABS(CIB).LT.1.0D-12) EB2(5)=PI05
      IF(DABS(CIB).LT.1.0D-12) GOTO 310
      EB2(5)=DATAN(SIB/CIB)
      IF(EB2(5).LT.0.0D0) EB2(5)=PI+EB2(5)
 310  EB2(5)=EB2(5)/PI180
      IF(EA(5).GT.9.0D1.AND.EB2(5).LT.9.0D1) EB2(5)=1.8D2-EB2(5)
      IF(EA(5).LT.9.0D1.AND.EB2(5).GT.9.0D1) EB2(5)=1.8D2-EB2(5)
C
      CALL DCRIT(EA,EB2,DD2)
      FB=F0/PI180
      R2(1)=R2(1)-0.18D3
      IF(R2(1).LT.0.0D0) R2(1)=R2(1)+0.36D3
      CALL RADIANT(EB2,FB,VXZ,VYR,VZR,R2)
C
 400  CONTINUE
      RETURN
C
      END
C ---------------------------------------------------------------------
      SUBROUTINE PMETH(WM1,WM2,EA,R1,EB1,DD1,R2,EB2,DD2)
C ---------------------------------------------------------------------
C      - (P) PORTER'S METHOD
C     (PORTER J.G.: 1952, COMETS AND METEOR STREAMS, CHAPMAN AND HALL
C     LTD., LONDON.)
C
C     INPUT (FIELD OF VALUES "EA(I)", I = 1, 2,..., 6) AND OUTPUT (FIEL-
C     DS OF VALUES "R1(I)", "R2(I)", "WM1(I)", WM2(I)", "EB1(I)", AND
C     "EB2(I)", I = 1, 2, 3,..., 6 OR 5, RESPECTIVELY), PLEASE, SEE THE
C     MAIN PROGRAM. AS THE OUTPUT, THERE ARE MOREOVER THE VALUES "DD1"
C     AND "DD2" OF D-DISCRIMINANT CORRESPONDING TO THE POST- AND PRE-
C     PERIHELION ARC OF THE PARENT BODY ORBIT, RESPECTIVELY.
C
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER  I,IM,IYR,J
      DIMENSION EA(6),EB1(5),EB2(5),R1(6),R2(6),WM1(6),WM2(6),WEL(7)
Cx      DIMENSION EZ(6)
C
      PI=4.0D0*DATAN(1.0D0)
      PI05=PI/2.0D0
      PI2=2.0D0*PI
      PI180=PI/0.18D3
      GAUSS=.17202098950D-1
      GS2=GAUSS*GAUSS
      EPS=2.3439291D1
      CEPS=DCOS(PI180*EPS)
      SEPS=DSIN(PI180*EPS)
      NTARG=3
      NCENT=11
C
      CHO=DCOS(EA(4)*PI180)
      SHO=DSIN(EA(4)*PI180)
      CSK=DCOS(EA(5)*PI180)
      SSK=DSIN(EA(5)*PI180)
C
      DD1=0.999D99
      IF(WM1(2).GT.1.0D2.OR.WM1(5).LT.0.0D0) GOTO 200
      CVK=DCOS(WM1(5)*PI180)
      SVK=DSIN(WM1(5)*PI180)
      WVFI=DSQRT(1.0D0+EA(2)*EA(2)+2.0D0*EA(2)*CVK)
      CHFI=-EA(2)*SVK/WVFI
      SHFI=(1.0D0+EA(2)*CVK)/WVFI
      QKONTR=SHFI*SHFI+CHFI*CHFI
      IF(DABS(QKONTR-1.0D0).GT.1.0D-12) WRITE(*,5)
   5  FORMAT(/,'  ERROR IN SUBROUTINE `PMETH` ((SIN**2(HFI) + COS**2(VFI
     *).NE.1).')
      IF(DABS(QKONTR-1.0D0).GT.1.0D-12) STOP     
      IF(SHFI.LT.-1.0D-12) WRITE(*,7)
   7  FORMAT(/,'  ERROR IN SUBROUTINE `PMETH` (ANGLE `HFI` > PI).')
      IF(SHFI.LT.-1.0D-12) STOP     
      IF(DABS(CHFI).LT.1.0D-12) HFI=PI05
      IF(DABS(CHFI).LT.1.0D-12) GOTO 10
      HFI=DATAN(SHFI/CHFI)
      IF(HFI.LT.0.0D0) HFI=PI+HFI
  10  TMN1=WM1(1)
c      CALL PLEPH(TMN1,NTARG,NCENT,RRD)
      T20=(TMN1-2.451545D6)/3.6525D4
      CALL EARTHrv(T20,XZ,YZ,ZZ,VXZ,VYZ,VZZ,WEL)
      IF(DABS(XZ).LT.1.0D-12) THEN
      WL=PI/2.0D0
      IF(YZ.LT.0.0D0) WL=WL+PI
      ELSE
      WL=DATAN(YZ/XZ)
      WL=DABS(WL)
      IF(XZ.LT.0.0D0.AND.YZ.GE.0.0D0) WL=PI-WL
      IF(XZ.LT.0.0D0.AND.YZ.LT.0.0D0) WL=PI+WL
      IF(XZ.GE.0.0D0.AND.YZ.LT.0.0D0) WL=PI2-WL
      END IF
C      WL=WEL(7)+WEL(4)
C      IF(WL.GE.PI2) WL=WL-PI2
Cx      EZ(1)=WEL(1)
Cx      EZ(2)=WEL(2)
Cx      EZ(3)=WEL(4)-WEL(5)
Cx      IF(EZ(3).LT.0.0D0) EZ(3)=EZ(3)+PI2
Cx      EZ(4)=WEL(5)
Cx      EZ(5)=WEL(6)
Cx      FZM=WEL(7)
Cx      EZ(1)=EZ(1)/(1.0D0-EZ(2))
Cx      FZM=EZ(6)
      RZ=DSQRT(XZ*XZ+YZ*YZ+ZZ*ZZ)
      VYR=VYZ*CEPS-VZZ*SEPS
      VZR=VYZ*SEPS+VZZ*CEPS
      R1(1)=WL/PI180
      R1(6)=TMN1
Cx      RZ=DSQRT(X*X+Y*Y+Z*Z)
      IF(EA(1).GT.RZ) GOTO 20
      IF(EA(2).GE.0.9999999D0) GOTO 30
      QQK=EA(1)*(1.0D0+EA(2))/(1.0D0-EA(2))
      IF(QQK.GE.RZ) GOTO 30
      AV=GAUSS*DSQRT(2.0D0/QQK-(1.0D0-EA(2))/EA(1))
      GOTO 40
  20  CONTINUE
      AV=GAUSS*DSQRT((1.0D0+EA(2))/EA(1))
      GOTO 40
  30  CONTINUE
      AV=GAUSS*DSQRT(2.0D0/RZ-(1.0D0-EA(2))/EA(1))
  40  AR3S=(WM1(5)+EA(3))*PI180-HFI
      C3S=DCOS(AR3S)
      S3S=DSIN(AR3S)
      VX=-AV*(C3S*CHO-S3S*CSK*SHO)
      VY=-AV*(C3S*SHO+S3S*CSK*CHO)
      VZ=-AV*S3S*SSK
C
      WMEN=(XZ*VY-YZ*VX)*(XZ*VY-YZ*VX)+(YZ*VZ-ZZ*VY)*(YZ*VZ-ZZ*VY)
      WMEN=DSQRT(WMEN+(ZZ*VX-XZ*VZ)*(ZZ*VX-XZ*VZ))
      CIB=(XZ*VY-YZ*VX)/WMEN
      SIB=DSQRT(1.0D0-CIB*CIB)
      IF(DABS(CIB).LT.1.0D-12) GOTO 50
      EB1(5)=DATAN(SIB/CIB)
      IF(EB1(5).LT.0.0D0) EB1(5)=PI+EB1(5)
      GOTO 60
  50  EB1(5)=PI05
  60  WMEN=(YZ*VZ-ZZ*VY)*(YZ*VZ-ZZ*VY)+(ZZ*VX-XZ*VZ)*(ZZ*VX-XZ*VZ)
      WMEN=DSQRT(WMEN)
      CVO=(XZ*VZ-ZZ*VX)/WMEN
      SVO=(YZ*VZ-ZZ*VY)/WMEN
      IF(DABS(CVO).LT.1.0D-12) GOTO 70
      EB1(4)=DATAN(SVO/CVO)
      EB1(4)=DABS(EB1(4))
      GOTO 80
  70  EB1(4)=PI05
  80  CONTINUE
      IF(CVO.GE.0.0D0.AND.SVO.LT.0.0D0) EB1(4)=PI2-EB1(4)
      IF(CVO.LT.0.0D0.AND.SVO.GE.0.0D0) EB1(4)=PI-EB1(4)
      IF(CVO.LT.0.0D0.AND.SVO.LT.0.0D0) EB1(4)=PI+EB1(4)
C
      AKI=2.0D0/RZ-AV*AV/GS2
      EB1(2)=DSQRT(1.0D0-RZ*AKI*(2.0D0-RZ*AKI)*SHFI*SHFI)
      IF(DABS(EB1(2)-1.0D0).LT.1.0D-9) GOTO 90
      EB1(1)=(1.0D0-EB1(2))/AKI
      GOTO 100
  90  EB1(1)=RZ*SHFI*SHFI
 100  CF=(EB1(1)*(1.0D0+EB1(2))/RZ-1.0D0)/EB1(2)
      IF(CF.GE.0.999999999D0) GOTO 115
      IF(CF.LE.-0.999999999D0) GOTO 117
      SF=DSQRT(1.0D0-CF*CF)
      IF(DABS(CF).LT.1.0D-12) GOTO 110
      F=DATAN(SF/CF)
      IF(F.LT.0.0D0) F=PI+F
      GOTO 120
 110  F=PI05
      GOTO 120
 115  F=0.0D0
      GOTO 120
 117  F=PI
 120  CONTINUE
      IF(HFI.LT.PI05) F=PI2-F
C
      C2S=(XZ*CVO+YZ*SVO)/RZ
      IF(DABS(CIB).LT.1.0D-9) GOTO 130
      S2S=(YZ*CVO-XZ*SVO)/RZ/CIB
      GOTO 140
 130  S2S=ZZ/RZ/SIB
 140  CONTINUE
      IF(DABS(C2S).LT.1.0D-12) GOTO 150
      U2S=DATAN(S2S/C2S)
      U2S=DABS(U2S)
      GOTO 160
 150  U2S=PI05
 160  CONTINUE
      IF(C2S.GE.0.0D0.AND.S2S.LT.0.0D0) U2S=PI2-U2S
      IF(C2S.LT.0.0D0.AND.S2S.GE.0.0D0) U2S=PI-U2S
      IF(C2S.LT.0.0D0.AND.S2S.LT.0.0D0) U2S=PI+U2S
      EB1(3)=U2S-F
      IF(EB1(3).LT.0.0D0) EB1(3)=EB1(3)+PI2
C
      EB1(3)=EB1(3)/PI180
      EB1(4)=EB1(4)/PI180
      EB1(5)=EB1(5)/PI180
      CALL DCRIT(EA,EB1,DD1)
      FB=F/PI180
      R1(1)=R1(1)-0.18D3
      IF(R1(1).LT.0.0D0) R1(1)=R1(1)+0.36D3
      CALL RADIANT(EB1,FB,VXZ,VYR,VZR,R1)
C
C   2-ND ARC
C
 200  CONTINUE
      DD2=0.999D99
      IF(WM2(2).GT.1.0D2.OR.WM2(5).GT.0.36D3) RETURN
      CVK=DCOS(WM2(5)*PI180)
      SVK=DSIN(WM2(5)*PI180)
      WVFI=DSQRT(1.0D0+EA(2)*EA(2)+2.0D0*EA(2)*CVK)
      CHFI=-EA(2)*SVK/WVFI
      SHFI=(1.0D0+EA(2)*CVK)/WVFI
      QKONTR=SHFI*SHFI+CHFI*CHFI
      IF(DABS(QKONTR-1.0D0).GT.1.0D-12) WRITE(*,5)
      IF(DABS(QKONTR-1.0D0).GT.1.0D-12) STOP     
      IF(SHFI.LT.-1.0D-12) WRITE(*,7)
      IF(SHFI.LT.-1.0D-12) STOP     
      IF(DABS(CHFI).LT.1.0D-12) HFI=PI05
      IF(DABS(CHFI).LT.1.0D-12) GOTO 210
      HFI=DATAN(SHFI/CHFI)
      IF(HFI.LT.0.0D0) HFI=PI+HFI
 210  TMN2=WM2(1)
c      CALL PLEPH(TMN2,NTARG,NCENT,RRD)
      T20=(TMN2-2.451545D6)/3.6525D4
      CALL EARTHrv(T20,XZ,YZ,ZZ,VXZ,VYZ,VZZ,WEL)
      IF(DABS(XZ).LT.1.0D-12) THEN
      WL=PI/2.0D0
      IF(YZ.LT.0.0D0) WL=WL+PI
      ELSE
      WL=DATAN(YZ/XZ)
      WL=DABS(WL)
      IF(XZ.LT.0.0D0.AND.YZ.GE.0.0D0) WL=PI-WL
      IF(XZ.LT.0.0D0.AND.YZ.LT.0.0D0) WL=PI+WL
      IF(XZ.GE.0.0D0.AND.YZ.LT.0.0D0) WL=PI2-WL
      END IF
C      WL=WEL(7)+WEL(4)
C      IF(WL.GE.PI2) WL=WL-PI2
Cx      EZ(1)=WEL(1)
Cx      EZ(2)=WEL(2)
Cx      EZ(3)=WEL(4)-WEL(5)
Cx      IF(EZ(3).LT.0.0D0) EZ(3)=EZ(3)+PI2
Cx      EZ(4)=WEL(5)
Cx      EZ(5)=WEL(6)
Cx      FZM=WEL(7)
Cx      EZ(1)=EZ(1)/(1.0D0-EZ(2))
Cx      FZM=EZ(6)
      VYR=VYZ*CEPS-VZZ*SEPS
      VZR=VYZ*SEPS+VZZ*CEPS
      R2(1)=WL/PI180
      R2(6)=TMN2
      RZ=DSQRT(XZ*XZ+YZ*YZ+ZZ*ZZ)
      IF(EA(1).GT.RZ) GOTO 220
      IF(EA(2).GE.0.9999999D0) GOTO 230
      QQK=EA(1)*(1.0D0+EA(2))/(1.0D0-EA(2))
      IF(QQK.GE.RZ) GOTO 230
      AV=GAUSS*DSQRT(2.0D0/QQK-(1.0D0-EA(2))/EA(1))
      GOTO 240
 220  CONTINUE
      AV=GAUSS*DSQRT((1.0D0+EA(2))/EA(1))
      GOTO 240
 230  CONTINUE
      AV=GAUSS*DSQRT(2.0D0/RZ-(1.0D0-EA(2))/EA(1))
 240  AR3S=(WM2(5)+EA(3))*PI180-HFI
      C3S=DCOS(AR3S)
      S3S=DSIN(AR3S)
      VX=-AV*(C3S*CHO-S3S*CSK*SHO)
      VY=-AV*(C3S*SHO+S3S*CSK*CHO)
      VZ=-AV*S3S*SSK
C
      WMEN=(XZ*VY-YZ*VX)*(XZ*VY-YZ*VX)+(YZ*VZ-ZZ*VY)*(YZ*VZ-ZZ*VY)
      WMEN=DSQRT(WMEN+(ZZ*VX-XZ*VZ)*(ZZ*VX-XZ*VZ))
      CIB=(XZ*VY-YZ*VX)/WMEN
      SIB=DSQRT(1.0D0-CIB*CIB)
      IF(DABS(CIB).LT.1.0D-12) GOTO 250
      EB2(5)=DATAN(SIB/CIB)
      IF(EB2(5).LT.0.0D0) EB2(5)=PI+EB2(5)
      GOTO 260
 250  EB2(5)=PI05
 260  WMEN=(YZ*VZ-ZZ*VY)*(YZ*VZ-ZZ*VY)+(ZZ*VX-XZ*VZ)*(ZZ*VX-XZ*VZ)
      WMEN=DSQRT(WMEN)
      CVO=(XZ*VZ-ZZ*VX)/WMEN
      SVO=(YZ*VZ-ZZ*VY)/WMEN
      IF(DABS(CVO).LT.1.0D-12) GOTO 270
      EB2(4)=DATAN(SVO/CVO)
      EB2(4)=DABS(EB2(4))
      GOTO 280
 270  EB2(4)=PI05
 280  CONTINUE
      IF(CVO.GE.0.0D0.AND.SVO.LT.0.0D0) EB2(4)=PI2-EB2(4)
      IF(CVO.LT.0.0D0.AND.SVO.GE.0.0D0) EB2(4)=PI-EB2(4)
      IF(CVO.LT.0.0D0.AND.SVO.LT.0.0D0) EB2(4)=PI+EB2(4)
C
      AKI=2.0D0/RZ-AV*AV/GS2
      EB2(2)=DSQRT(1.0D0-RZ*AKI*(2.0D0-RZ*AKI)*SHFI*SHFI)
      IF(DABS(EB2(2)-1.0D0).LT.1.0D-9) GOTO 290
      EB2(1)=(1.0D0-EB2(2))/AKI
      GOTO 300
 290  EB2(1)=RZ*SHFI*SHFI
 300  CF=(EB2(1)*(1.0D0+EB2(2))/RZ-1.0D0)/EB2(2)
      IF(CF.GE.0.999999999D0) GOTO 315
      IF(CF.LE.-0.999999999D0) GOTO 317
      SF=DSQRT(1.0D0-CF*CF)
      IF(DABS(CF).LT.1.0D-12) GOTO 310
      F=DATAN(SF/CF)
      IF(F.LT.0.0D0) F=PI+F
      GOTO 320
 310  F=PI05
      GOTO 320
 315  F=0.0D0
      GOTO 320
 317  F=PI
 320  CONTINUE
      IF(HFI.LT.PI05) F=PI2-F
C
      C2S=(XZ*CVO+YZ*SVO)/RZ
      IF(DABS(CIB).LT.1.0D-9) GOTO 330
      S2S=(YZ*CVO-XZ*SVO)/RZ/CIB
      GOTO 340
 330  S2S=ZZ/RZ/SIB
 340  CONTINUE
      IF(DABS(C2S).LT.1.0D-12) GOTO 350
      U2S=DATAN(S2S/C2S)
      U2S=DABS(U2S)
      GOTO 360
 350  U2S=PI05
 360  CONTINUE
      IF(C2S.GE.0.0D0.AND.S2S.LT.0.0D0) U2S=PI2-U2S
      IF(C2S.LT.0.0D0.AND.S2S.GE.0.0D0) U2S=PI-U2S
      IF(C2S.LT.0.0D0.AND.S2S.LT.0.0D0) U2S=PI+U2S
      EB2(3)=U2S-F
      IF(EB2(3).LT.0.0D0) EB2(3)=EB2(3)+PI2
C
      EB2(3)=EB2(3)/PI180
      EB2(4)=EB2(4)/PI180
      EB2(5)=EB2(5)/PI180
      CALL DCRIT(EA,EB2,DD2)
      FB=F/PI180
      R2(1)=R2(1)-0.18D3
      IF(R2(1).LT.0.0D0) R2(1)=R2(1)+0.36D3
      CALL RADIANT(EB2,FB,VXZ,VYR,VZR,R2)
C
      RETURN
C
      END
C ---------------------------------------------------------------------
      SUBROUTINE RADIANT(EB,FB,VXZ,VYZ,VZZ,RR)
C ---------------------------------------------------------------------
C      - CALCULATION OF METEOR SHOWER RADIANT, IF THE FOLLOWING INPUT
C     QUANTITIES OF THE ORBIT CROSSING THE EARTH'S ORBIT ARE WELL-KNOWN:
C
C        EB(1) - PERIHELION DISTANCE [AU]
C        EB(2) - ECCENTRICITY
C        EB(3) - ARGUMENT OF PERIHELION [DEGREES]
C        EB(4) - LONGITUDE OF NODE [DEGREES]
C        EB(5) - INCLINATION TO ECLIPTIC [DEGREES]
C        FB - TRUE ANOMALY OF THE POINT OF INTERSECTION OF CALCULATED
C            SHOWER AND EARTH ORBITS [DEGREES]
C        VXZ, VYZ, VZZ - HELIOCENTIC EQUATORIAL RECTANULAR COMPONENTS
C           IF THE EARTH [AU/DAY]
C
C     OUTPUT:
C     THE QUANTITIES DESCRIBING THE RADIANT OF POTENTIAL METEOR SHOWER
C        RR(2) - RIGHT ASCENSION OF THE PREDICTED RADIANT [DEGREES]
C        RR(3) - DECLINATION OF THE RADIANT [DEGREES]
C        RR(4) - PREDICTED GEOCENTRIC VELOCITY OF METEORS [KM/S]
C        RR(5) - PREDICTED HELIOCENTRIC VELOCITY OF METEORS [KM/S]
C
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION RR(6),EB(5)
Cx      ,EZ(6)
C
      PI=4.0D0*DATAN(1.0D0)
      PI05=PI/2.0D0
      PI2=2.0D0*PI
      PI180=PI/0.18D3
      GAUSS=0.17202098950D-1
      AUKM=0.149597870D9
      D24=8.6400D4
C      EPS=2.344579D1......1950.0
      EPS=2.3439291D1
C
      CEPS=DCOS(PI180*EPS)
      SEPS=DSIN(PI180*EPS)
      CHO=DCOS(PI180*EB(4))
      SHO=DSIN(PI180*EB(4))
      CSK=DCOS(PI180*EB(5))
      SSK=DSIN(PI180*EB(5))
      RK=EB(1)*(1.0D0+EB(2))/(1.0D0+EB(2)*DCOS(FB*PI180))
Cx      CHZ=DCOS(EZ(4))
Cx      SHZ=DSIN(EZ(4))
Cx      CSZ=DCOS(EZ(5))
Cx      SSZ=DSIN(EZ(5))
Cx      RZ=EZ(1)*(1.0D0-EZ(2)*EZ(2))/(1.0D0+EZ(2)*DCOS(FZ))
C    
      VHR=PI180*FB
      CVH=DCOS(VHR)
      SVH=DSIN(VHR)
      WHFI=DSQRT(1.0D0+EB(2)*EB(2)+2.0D0*EB(2)*CVH)
      CHFI=-EB(2)*SVH/WHFI
      SHFI=(1.0D0+EB(2)*CVH)/WHFI
      QKONTR=SHFI*SHFI+CHFI*CHFI
      IF(DABS(QKONTR-1.0D0).GT.1.0D-9) THEN
      WRITE(*,100)
 100  FORMAT(/,'  ERROR IN SUBROUTINE `RADIANT` (ABS(SIN**2(HFI) + COS**
     *2(HFI)).GT.1)')
      WRITE(*,*) 'QKONTR =',QKONTR
      ENDIF
      IF(DABS(QKONTR-1.0D0).GT.1.0D-9) STOP
      IF(DABS(CHFI).LT.1.0D-9) HFI=PI05
      IF(DABS(CHFI).LT.1.0D-9) GOTO 120
      HFI=DATAN(SHFI/CHFI)
      HFI=DABS(HFI)
      IF(CHFI.GE.0.0D0.AND.SHFI.LT.0.0D0) HFI=PI2-HFI
      IF(CHFI.LT.0.0D0.AND.SHFI.GE.0.0D0) HFI=PI-HFI
      IF(CHFI.LT.0.0D0.AND.SHFI.LT.0.0D0) HFI=PI+HFI
 120  CONTINUE
      CFSO=DCOS(PI180*(FB+EB(3)))
      SFSO=DSIN(PI180*(FB+EB(3)))
      XK=RK*(CFSO*CHO-SFSO*CSK*SHO)
      YK=RK*(CFSO*SHO+SFSO*CSK*CHO)
      ZK=RK*SFSO*SSK
      C3=DCOS(PI180*(FB+EB(3))-HFI)
      S3=DSIN(PI180*(FB+EB(3))-HFI)
      AV=GAUSS*DSQRT(2.0D0/RK-(1.0D0-EB(2))/EB(1))
      VXE=-AV*(C3*CHO-S3*CSK*SHO)
      VYE=-AV*(C3*SHO+S3*CSK*CHO)
      VZE=-AV*S3*SSK
      VE=DSQRT(VXE*VXE+VYE*VYE+VZE*VZE)
      IF(DABS(AV-VE).GT.1.0D-9) WRITE(*,140)
 140  FORMAT(/,'  ERROR IN SUBROUTINE `RADIANT` (AV.NE.VE).')
      IF(DABS(AV-VE).GT.1.0D-9) STOP
      VXK=VXE
      VYK=VYE*CEPS-VZE*SEPS
      VZK=VYE*SEPS+VZE*CEPS
      VH=DSQRT(VXK*VXK+VYK*VYK+VZK*VZK)
      IF(DABS(AV-VH).GT.1.0D-9) WRITE(*,160)
 160  FORMAT(/,'  ERROR IN SUBROUTINE `RADIANT` (AV.NE.VH).')
      IF(DABS(AV-VH).GT.1.0D-9) STOP
C
Cx      CVH=DCOS(FZ)
Cx      SVH=DSIN(FZ)
Cx      WHFI=DSQRT(1.0D0+EZ(2)*EZ(2)+2.0D0*EZ(2)*CVH)
Cx      CHFI=-EZ(2)*SVH/WHFI
Cx      SHFI=(1.0D0+EZ(2)*CVH)/WHFI
Cx      QKONTR=SHFI*SHFI+CHFI*CHFI
Cx      IF(DABS(QKONTR-1.0D0).GT.1.0D-12) WRITE(*,200)
Cx 200  FORMAT(/,'  ERROR IN SUBROUTINE `RADIANT` (ABS(SIN**2(HFI) + COS**
Cx     *2(HFI)).GT.1)')
Cx      IF(DABS(QKONTR-1.0D0).GT.1.0D-12) STOP
Cx      IF(DABS(CHFI).LT.1.0D-12) HFI=PI05
Cx      IF(DABS(CHFI).LT.1.0D-12) GOTO 220
Cx      HFI=DATAN(SHFI/CHFI)
Cx      HFI=DABS(HFI)
Cx      IF(CHFI.GE.0.0D0.AND.SHFI.LT.0.0D0) HFI=PI2-HFI
Cx      IF(CHFI.LT.0.0D0.AND.SHFI.GE.0.0D0) HFI=PI-HFI
Cx      IF(CHFI.LT.0.0D0.AND.SHFI.LT.0.0D0) HFI=PI+HFI
Cx 220  CONTINUE
Cx      CFSO=DCOS(FZ+EZ(3))
Cx      SFSO=DSIN(FZ+EZ(3))
Cx      XZ=RZ*(CFSO*CHZ-SFSO*CSZ*SHZ)
Cx      YZ=RZ*(CFSO*SHZ+SFSO*CSZ*CHZ)
Cx      ZZ=RZ*SFSO*SSZ
Cx      C3=DCOS(FZ+EZ(3)-HFI)
Cx      S3=DSIN(FZ+EZ(3)-HFI)
Cx      AVZ=GAUSS*DSQRT(2.0D0/RZ-1.0D0/EZ(1))
Cx      VXE=-AVZ*(C3*CHZ-S3*CSZ*SHZ)
Cx      VYE=-AVZ*(C3*SHZ+S3*CSZ*CHZ)
Cx      VZE=-AVZ*S3*SSZ
Cx      VEZ=DSQRT(VXE*VXE+VYE*VYE+VZE*VZE)
Cx      IF(DABS(AVZ-VEZ).GT.1.0D-12) WRITE(*,240)
Cx 240  FORMAT(/,'  ERROR IN SUBROUTINE `RADIANT` (AVZ.NE.VEZ).')
Cx      IF(DABS(AVZ-VEZ).GT.1.0D-12) STOP
Cx      VXZ=VXE
Cx      VYZ=VYE*CEPS-VZE*SEPS
Cx      VZZ=VYE*SEPS+VZE*CEPS
Cx      VHZ=DSQRT(VXZ*VXZ+VYZ*VYZ+VZZ*VZZ)
Cx      IF(DABS(AVZ-VHZ).GT.1.0D-12) WRITE(*,260)
Cx 260  FORMAT(/,'  ERROR IN SUBROUTINE `RADIANT` (AVZ.NE.VHZ).')
Cx      IF(DABS(AVZ-VHZ).GT.1.0D-12) STOP
C
      VX=-(VXK-VXZ)
      VY=-(VYK-VYZ)
      VZ=-(VZK-VZZ)
      VG=DSQRT(VX*VX+VY*VY+VZ*VZ)
      IF(DABS(VZ).LT.1.0D-9) D=0.0D0
      IF(DABS(VZ).LT.1.0D-9) GOTO 310
      HELVEC=DSQRT(1.0D0-VZ*VZ/VG/VG)
      D=DATAN(VZ/HELVEC/VG)
 310  CONTINUE
      IF(DABS(D-PI05).LT.1.0D-9.OR.DABS(D+PI05).LT.1.0D-9) A=0.0D0
      IF(DABS(D-PI05).LT.1.0D-9.OR.DABS(D+PI05).LT.1.0D-9) GOTO 340
      CDVEC=DCOS(D)
      QQS=VY/VG/CDVEC
      QQC=VX/VG/CDVEC
      QKONTR=QQS*QQS+QQC*QQC
      IF(DABS(QKONTR-1.0D0).GT.1.0D-9) WRITE(*,320)
 320  FORMAT(/,'  ERROR IN SUBROUTINE `RADIANT` (SIN**2(A) + COS**2(A).G
     *T.1).')
      IF(DABS(QKONTR-1.0D0).GT.1.0D-9) STOP
      IF(DABS(QQC).LT.1.0D-9) A=PI05
      IF(DABS(QQC).LT.1.0D-9) GOTO 340
      A=DATAN(QQS/QQC)
      A=DABS(A)
      IF(QQC.GE.0.0D0.AND.QQS.LT.0.0D0) A=PI2-A
      IF(QQC.LT.0.0D0.AND.QQS.GE.0.0D0) A=PI-A
      IF(QQC.LT.0.0D0.AND.QQS.LT.0.0D0) A=PI+A
 340  CONTINUE
      RR(2)=A/PI180
      RR(3)=D/PI180
      RR(4)=VG*AUKM/D24
      RR(5)=AV*AUKM/D24
      RETURN
C
      END
C ---------------------------------------------------------------------
      SUBROUTINE DCRIT(EA,EB,DD)
C ---------------------------------------------------------------------
C      - COMPUTATION OF D-DISCRIMINANT OF THE SOUTHWORTH-HAWKINS' D-CRI-
C     TERION EVALUATING SIMILARITY OF TWO ORBITS (SOUTHWORTH R.B.,
C     HAWKINS G.S.: 1963, SMITHSON. CONTR. ASTROPHYS. 7, 261.)
C
C INPUT:
C     EA(1), EB(1) - PERIHELION DISTANCES OF THE ORBITS [AU]
C     EA(2), EB(2) - ECCENTRICITIES
C     EA(3), EB(3) - ARGUMENTS OF PERIHELIA [DEGREES]
C     EA(4), EB(4) - LONGITUDES OF NODES [DEGREES]
C     EA(5), EB(5) - INCLINATIONS TO ECLIPTIC [DEGREES]
C
C OUTPUT:
C     DD - RESULTANT VALUE OF D-DISCRIMINANT
C
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION EA(6),EB(5)
C
      PI=4.0D0*DATAN(1.0D0)
      PI180=PI/0.18D3
C
      DQ=EA(1)-EB(1)
      DE=EA(2)-EB(2)
      DSO=(EA(3)-EB(3))*PI180
      DHO=(EA(4)-EB(4))*PI180
      DSK=(EA(5)-EB(5))*PI180
      SSKA=DSIN(EA(5)*PI180)
      SSKB=DSIN(EB(5)*PI180)
      POMSK=DSIN(DSK/2.0D0)
      POMHO=DSIN(DHO/2.0D0)
      DAI=4.0D0*POMSK*POMSK+SSKA*SSKB*4.0D0*POMHO*POMHO
      IF(DAI.GE.3.9999999D0.AND.DAI.LT.4.0000001D0) DAI=3.9999999D0
      IF(DAI.LT.4.0D0) GOTO 700
      WRITE(*,600)
 600  FORMAT(/,'  ERROR IN SUBROUTINE `DCRIT` (QUANTITY DAI.GE.4).')
      WRITE(*,*) DAI
      STOP
 700  SECSK=1.0D0/SQRT(1.0D0-DAI/4.0D0)
      ARGS=DCOS((EA(5)+EB(5))*PI180/2.0D0)*POMHO*SECSK
      JSGAR=-1
      IF(ARGS.GE.0.0D0) JSGAR=1
      POM=DABS(DHO)
      JSIGN=2
      IF(POM.GT.PI) JSIGN=-2
      POM=DABS(ARGS)
      IF(POM.LT.1.0000001D0) GOTO 720
      WRITE(*,620)
 620  FORMAT(/,'  ERROR IN SUBROUTINE `DCRIT` (ABS(ARGS).GT.1).')
      STOP
 720  CONTINUE
      IF(POM.LT.1) GOTO 740
      DBP=DSO+JSIGN*JSGAR*PI/2.0D0
      GOTO 760
 740  ODM=ARGS/DSQRT(1.0D0-ARGS*ARGS)
      DBP=DSO+JSIGN*DATAN(ODM)
 760  POM=DSIN(DBP/2.0D0)
      DD=DSQRT(DQ*DQ+DE*DE+DAI+(EA(2)+EB(2))*(EA(2)+EB(2))*POM*POM)
      RETURN
C
      END
C ---------------------------------------------------------------------
      SUBROUTINE APPROACH(TMIN,TMAX,EA,WM1,WM2)
C ---------------------------------------------------------------------
C      - CALCULATION OF SOME QUANTITIES CHARACTERIZING THE MEAN ORBIT OF
C        THE METEOR SHOWER
C        
C INPUT:
C    "TMIN", "TMAX" - BEGINNING AND END OF YEAR OF INVESTIGATION OF THE
C       POTENTIAL SHOWER GIVEN IN JULIAN CENTURIES FROM 1950.0.
C    ORBITAL ELEMENTS OF THE ORBIT OF THE PARENT BODY
C       EA(1) - PERIHELION DISTANCE [AU]
C       EA(2) - ECCENTRICITY
C       EA(3) - ARGUMENT OF PERIHELION [DEGREES]
C       EA(4) - LONGITUDE OF NODE [DEGREES]
C       EA(5) - INCLINATION TO ECLIPTIC [DEGREES]
C       YEAR, MONTH, DAY - TIME OF PERIHELION PASSAGE OF THE PARENT BODY
C       (CONVERTED TO: EA(6) - THE TIME IN JULIAN CENTURIES FROM
C                              JD = 2433282.423)
C
C OUTPUT:
C   SOME QUANTITIES DESCRIBING THEORETICAL RADIANTS AT BOTH ARCS (PRE-
C   PERIHELION AND POST-PERIHELION) OF THE ORBIT OF THE PARENT BODY;
C   THE OUTPUT QUANTITIES AT THE POST-PERIHELION ARC ARE:
C      WM1(1) - MOMENT OF MINIMUM DISTANCE BETWEEN THE EARTH AND THE
C               ORBIT OF THE PARENT BODY [JULIAN CENTURIES FROM 1950.0]
C      WM1(2) - MINIMUM DISTANCE BETWEEN POST-PERIHELION ARC OF THE
C               PARENT BODY ORBIT AND EARTH'S ORBIT [AU]
C      WM1(3) - DISTANCE OF THE PARENT BODY FROM THE EARTH AT THE MOMENT
C               OF MAXIMUM [AU]
C      WM1(4) - TIME INTERVAL BETWEEN THE PASSAGES OF THE PARENT BODY
C               AND EARTH ACROSS THE NEAREST POINTS OF POST-PERIHELION
C               ARC OF THE PARENT BODY ORBIT AND EARTH'S ORBIT [DAYS]
C      WM1(5) - TRUE ANOMALY OF THE PARENT BODY IN THE MOMENT OF ITS
C               NEAREST POST-PERIHELION APPROACH TO THE EARTH'S ORBIT
C               [DEGREES]
C      WM1(6) - TRUE ANOMALY OF THE EARTH IN THE MOMENT OF ITS NEAREST
C               APPROACH TO THE POST-PERIHELION ARC OF THE PARENT BODY
C               ORBIT [DEGREES]
C   THE OUTPUT QUANTITIES AT THE PRE-PERIHELION ARC ARE ANALOGOUSLY
C   STORED IN FIELD "WM2(I)", WHERE "I = 1, 2, 3,..., 6".
C
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I,J)
      DIMENSION EA(6),WM1(6),WM2(6),WEL(7)
Cx      ,EZ(6)
C
      PI=4.0D0*DATAN(1.0D0)
      PI180=PI/0.18D3
      GAUSS=.17202098950D-1
      EPS=2.3439291D1*PI180
      CEPS=DCOS(EPS)
      SEPS=DSIN(EPS)
      NTARG=3
      NCENT=11
C
      TMN1=(TMIN+TMAX)/2.0D0
      TMN2=TMN1
      DT0=(TMAX-TMIN)/3.6525D2
C
      CHO=DCOS(PI180*EA(4))
      SHO=DSIN(PI180*EA(4))
      CSK=DCOS(PI180*EA(5))
      SSK=DSIN(PI180*EA(5))
      PRM=EA(1)*(1.0D0+EA(2))
      IF(EA(2).GT.0.99999999D0) GOTO 30
      A=EA(1)/(1.0D0-EA(2))
      PP05=PI*A*DSQRT(A)/GAUSS
      PP=2.0D0*PP05
C
C  === 1-ST ARC ===
C
  30  VKD=0.0D0
      VKU=0.179D3
      DVK=1.0D0
      TA=TMIN
      TB=TMAX
      DT=DT0
C
  50  WM1(2)=1.111D3
C >      DO 200 VK=VKD,VKU,DVK
      IVKU=(VKU+1.0D-11-VKD)/DVK
      DO 200 IVK=0,IVKU
      VK=VKD+IVK*DVK
      WEN=1.0D0+EA(2)*DCOS(PI180*VK)
      RK=1.0D6
      IF(DABS(WEN).LT.1.0D-12) GOTO 100
      RK=PRM/WEN
 100  CVS=DCOS(PI180*(VK+EA(3)))
      SVS=DSIN(PI180*(VK+EA(3)))
      XK=RK*(CVS*CHO-SVS*CSK*SHO)
      YK=RK*(CVS*SHO+SVS*CSK*CHO)
      ZK=RK*SVS*SSK
C
C >      DO 140 ETT=TA,TB,DT
      IETTU=(TB+1.0D-11-TA)/DT
      DO 140 IETT=0,IETTU
      ETT=TA+IETT*DT
c      CALL PLEPH(ETT,NTARG,NCENT,RRD)
      T20=(ETT-2.451545D6)/3.6525D4
      CALL EARTHrv(T20,XZ,YZ,ZZ,VXZ,VYZ,VZZ,WEL)
      DR=DSQRT((XK-XZ)*(XK-XZ)+(YK-YZ)*(YK-YZ)+(ZK-ZZ)*(ZK-ZZ))
      IF(DR.GT.1.0D2) GOTO 240
      IF(DR.GE.WM1(2)) GOTO 140
      WM1(2)=DR
      WM1(5)=VK
Cx      WM1(6)=WEL(7)
C-      WM1(6)=EZ(6)
      TMN1=ETT
 140  CONTINUE
 200  CONTINUE
C      BOTH LABELS 200 AND 240 ARE NECESSARY!
 240  CONTINUE
      IF(DVK.GT.0.5D0) GOTO 245
      IF(DABS(WM1(2)-RMOLD).LT.1.0D-5) GOTO 260
 245  VKD=WM1(5)-1.0D0*DVK
      VKU=WM1(5)+1.0D0*DVK
      TA=TMN1-1.0D0*DT
      TB=TMN1+1.0D0*DT
      DVK=0.1D-1*DVK
      DT=0.1D-1*DT
      RMOLD=WM1(2)
      GOTO 50
C
 260  CONTINUE
C
      IF(WM1(2).GT.1.0D2) GOTO 530
Cx      DTM=3.6525D4*(TMN1-EA(6))
C?      DTM=TMN1-EA(6)
C?      READ(*,*) JXXX
C?      EQX=EA(1)
C?      EEX=EA(2)
C?      CALL KEPLER(EQX,EEX,DTM,QWV)
C?      WENK=1.0D0+EA(2)*DCOS(QWV)
C?      RK=PRM/WENK
C?      CVS=DCOS(QWV+PI180*EA(3))
C?      SVS=DSIN(QWV+PI180*EA(3))
C?      XK=RK*(CVS*CHO-SVS*CSK*SHO)
C?      YK=RK*(CVS*SHO+SVS*CSK*CHO)
C?      ZK=RK*SVS*SSK
C?      CALL PLEPH(TMN1,NTARG,NCENT,RRD)
C      CALL ECLXYZ(TMN1,XZ,YZ,ZZ,WEL)
C?      XZ=RRD(1)
C?     YZ=RRD(2)*CEPS+RRD(3)*SEPS
C?      ZZ=-RRD(2)*SEPS+RRD(3)*CEPS
C      VXZ=RRD(4)
C      VYZ=RRD(5)*CEPS+RRD(6)*SEPS
C      VZZ=-RRD(5)*SEPS+RRD(6)*CEPS
C      CALL ELEMS(ETT,XZ,YZ,VZ,VXZ,VYZ,VZZ,EZ)
C?      DR=DSQRT((XK-XZ)*(XK-XZ)+(YK-YZ)*(YK-YZ)+(ZK-ZZ)*(ZK-ZZ))
C?      WM1(3)=DR
C
C???      WM1(4)=DTM
C?      WM1(4)=TMN1
C?      QWV=WM1(5)*PI180
C?      CALL TIMEANOM(EQX,EEX,QWV,DTT)
C?      WM1(4)=WM1(4)-DTT
C?      IF(EA(2).GT.0.99999999D0) GOTO 270
C? 263  CONTINUE
C?      IF(WM1(4).LE.PP05) GOTO 265
C?      WM1(4)=WM1(4)-PP
C?      GOTO 263
C? 265  CONTINUE
C?      IF(WM1(4).GE.-PP05) GOTO 270
C?      WM1(4)=WM1(4)+PP
C?      GOTO 265
C
 270  WM1(1)=TMN1
Ct      WRITE(*,*) 'TMN1, EA(6) = ',TMN1,EA(6)
C
C  === 2-ND ARC ===
C
 530  VKD=0.359D3
      VKU=0.181D3
      DVK=-1.0D0
      TA=TMIN
      TB=TMAX
      DT=DT0
C
 550  WM2(2)=1.111D3
C >      DO 700 VK=VKD,VKU,DVK
      IVKU=(VKU+1.0D-11-VKD)/DVK
      DO 700 IVK=0,IVKU
      VK=VKD+IVK*DVK
      WEN=1.0D0+EA(2)*DCOS(PI180*VK)
      RK=1.0D6
      IF(DABS(WEN).LT.1.0D-12) GOTO 600
      RK=PRM/WEN
 600  CVS=DCOS(PI180*(VK+EA(3)))
      SVS=DSIN(PI180*(VK+EA(3)))
      XK=RK*(CVS*CHO-SVS*CSK*SHO)
      YK=RK*(CVS*SHO+SVS*CSK*CHO)
      ZK=RK*SVS*SSK
C
C >      DO 640 ETT=TA,TB,DT
      IETTU=(TB+1.0D-11-TA)/DT
      DO 640 IETT=0,IETTU
      ETT=TA+IETT*DT
c      CALL PLEPH(ETT,NTARG,NCENT,RRD)
Ct      WRITE(*,*) 'APPROACH-3: Von z PLEPH.'
      T20=(ETT-2.451545D6)/3.6525D4
      CALL EARTHrv(T20,XZ,YZ,ZZ,VXZ,VYZ,VZZ,WEL)
      DR=DSQRT((XK-XZ)*(XK-XZ)+(YK-YZ)*(YK-YZ)+(ZK-ZZ)*(ZK-ZZ))
      IF(DR.GT.1.0D2) GOTO 740
      IF(DR.GE.WM2(2)) GOTO 640
      WM2(2)=DR
      WM2(5)=VK
Cx      WM2(6)=WEL(7)
C-      WM2(6)=EZ(6)
      TMN2=ETT
 640  CONTINUE
 700  CONTINUE
C      BOTH LABELS 700 AND 740 ARE NECESSARY!
 740  CONTINUE
      IF(DVK.GT.0.5D0) GOTO 745
      IF(DABS(WM2(2)-RMOLD).LT.1.0D-5) GOTO 760
 745  VKD=WM2(5)-1.0D0*DVK
      VKU=WM2(5)+1.0D0*DVK
      TA=TMN2-1.0D0*DT
      TB=TMN2+1.0D0*DT
      DVK=0.1D-1*DVK
      DT=0.1D-1*DT
      RMOLD=WM2(2)
      GOTO 550
C
 760  CONTINUE
C
      IF(WM2(2).GT.1.0D2) GOTO 800
Cx      DTM=3.6525D4*(TMN2-EA(6))
C?      DTM=TMN2-EA(6)
C?      EQX=EA(1)
C?      EEX=EA(2)
C?      CALL KEPLER(EQX,EEX,DTM,QWV)
C?      WENK=1.0D0+EA(2)*DCOS(QWV)
C?      RK=PRM/WENK
C?      CVS=DCOS(QWV+PI180*EA(3))
C?      SVS=DSIN(QWV+PI180*EA(3))
C?      XK=RK*(CVS*CHO-SVS*CSK*SHO)
C?      YK=RK*(CVS*SHO+SVS*CSK*CHO)
C?      ZK=RK*SVS*SSK
C?      CALL PLEPH(TMN2,NTARG,NCENT,RRD)
C      CALL ECLXYZ(TMN2,XZ,YZ,ZZ,WEL)
C?      XZ=RRD(1)
C?      YZ=RRD(2)*CEPS+RRD(3)*SEPS
C?      ZZ=-RRD(2)*SEPS+RRD(3)*CEPS
C      VXZ=RRD(4)
C      VYZ=RRD(5)*CEPS+RRD(6)*SEPS
C      VZZ=-RRD(5)*SEPS+RRD(6)*CEPS
C      CALL ELEMS(ETT,XZ,YZ,VZ,VXZ,VYZ,VZZ,EZ)
C?      DR=DSQRT((XK-XZ)*(XK-XZ)+(YK-YZ)*(YK-YZ)+(ZK-ZZ)*(ZK-ZZ))
C?      WM2(3)=DR
C
C???      WM2(4)=DTM
C?      WM2(4)=TMN2
C?      QWV=WM2(5)*PI180
C?      CALL TIMEANOM(EQX,EEX,QWV,DTT)
C?      WM2(4)=WM2(4)-DTT
C?      IF(EA(2).GT.0.99999999D0) GOTO 770
C? 763  CONTINUE
C?      IF(WM2(4).LE.PP05) GOTO 765
C?      WM2(4)=WM2(4)-PP
C?      GOTO 763
C? 765  CONTINUE
C?      IF(WM2(4).GE.-PP05) GOTO 770
C?      WM2(4)=WM2(4)+PP
C?      GOTO 765
C
 770  WM2(1)=TMN2
 800  CONTINUE
      RETURN
C
      END
C --------------------------------------------------------------------
      SUBROUTINE JULIAN(IROK,IMES,DEN,T20,ETT)
C --------------------------------------------------------------------
C     - CALCULATION OF JULIAN DATE FROM CIVIL DATE
C
C INPUT:
C    IROK, IMES, DEN - YEAR, MONTH, AND DAY OF THE DATE
C
C OUTPUT:
C    T20 - GIVEN TIME ACCOUNTED FROM THE BEGINNING OF EPOCH 1900.0
C          (JD = 2415020.0) AND EXPRESSED IN JULIAN CENTURIES
C    ETT - JDT OF THE GIVEN DATE
C
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      IMPLICIT REAL*8 (A-H,K-Z)
      IMPLICIT INTEGER (I,J)
      DIMENSION JMESI(13)
      DATA JMESI/0,31,28,31,30,31,30,31,31,30,31,30,31/
C
      IF(IROK.GT.1582) GOTO 101
      IF(IROK.EQ.1582) GOTO 102
      IF(IROK.GT.0.AND.IROK.LT.1582) GOTO 103
      IF(IROK.GE.-4713.AND.IROK.LT.0) GOTO 104
C
      ETT=(IROK+4713)*365
      ETT=ETT-INT((-IROK-4712)/4)
  100 T0JD=-0.5D0
      GOTO 210
  101 ETT=(IROK-1583)*365
      ETT=ETT+INT((IROK-1581)/4)
      ETT=ETT-INT((IROK-1501)/100)+INT((IROK-1201)/400)
      T0JD=2.2992385D6
      GOTO 210
  102 ETT=0.0D0
      IF(IMES.EQ.10.AND.DEN.GT.1.49999999D1) ETT=-1.0D1
      IF(IMES.GT.10) ETT=-1.0D1
      T0JD=2.2988835D6
      GOTO 210
  103 ETT=(IROK-1)*365
      ETT=ETT+INT((IROK-1)/4)
      T0JD=1.7214235D6
      GOTO 210
  104 ETT=(IROK+4713)*365
      ETT=ETT+INT((IROK+4715)/4)
      T0JD=-0.5D0
  210 YY=IROK/1.00D2-INT(IROK/100)
      ZZ=IROK/4.00D2-INT(IROK/400)
      IF(IROK.GT.1582.AND.YY.EQ.0.0D0.AND.ZZ.GT.0.0D0) GOTO 220
      YY=IROK/4.0D0-INT(IROK/4)
      IF(YY.EQ.0.0D0) JMESI(3)=29
      IF(IROK.EQ.-1) JMESI(3)=29
  220 CONTINUE
      DO 225 JJM=1,IMES-1
      ETT=ETT+JMESI(JJM+1)
  225 CONTINUE
      ETT=ETT+DEN-1.0D0+T0JD
      JMESI(3)=28
C
      T20=(ETT-2.415020D6)/3.6525D4
CC      ETT=(ETT-2.433282423D6)/3.6525D4
      RETURN
C
      END
C ---------------------------------------------------------------------
      SUBROUTINE INVJUL(TJD,IRK,IMS,DNI)
C ---------------------------------------------------------------------
C     - CALCULATION OF CIVIL DATE FROM JULIAN DATE (INVERSION OPERATION
C       TO THAT EXECUTED BY SUBROUTINE "JULIAN")
C
C INPUT:
C    TJD - TIME FROM THE BEGINNING OF EPOCH 2000.0 (JD = 2433282.423)
C          IN JULIAN CENTURIES (IT CAN BE CALCULATED FROM JULIAN DATE
C          "JD" IN DAYS AS: TJD=(JD-2433282.423)/36525)
C
C OUTPUT:
C    IRK, IMS, DNI - YEAR, MONTH, AND DAY OF THE CORRESPONDING CIVIL
C                    DATE
C
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      IMPLICIT REAL*8 (D,T)
      IMPLICIT INTEGER (I)
      DIMENSION IMESI(13)
      DATA IMESI/0,31,28,31,30,31,30,31,31,30,31,30,31/
C
CC      DNI=TJD*3.6525D4+2.433282423D6+1.0D0
      DNI=TJD+1.0D0
      IF(DNI.LT.0.5D0) GOTO 16
      IF(DNI.GE.0.5D0.AND.DNI.LT.1.7214245D6) GOTO 15
      IF(DNI.GE.1.7214245D6.AND.DNI.LT.2.2991615D6) GOTO 14
      IF(DNI.GE.2.2991615D6.AND.DNI.LT.2.2992395D6) GOTO 13
      IF(DNI.GE.2.2992395D6.AND.DNI.LT.2.4515455D6) GOTO 12
C
C  FROM 2000-JAN-1.0
  11  CONTINUE
      DNI=DNI-2.4515445D6
      IRK=INT(DNI/3.65242198D2)
      DN0=365*IRK+INT((IRK+3)/4.0D0)-INT((IRK-1)/1.0D2)+INT((IRK-1)/4.0D
     *2)
      DNI=DNI-DN0
      IF(DNI.LT.1.0D0) GOTO 30
      IRK=IRK+2000
      GOTO 20
  30  IRK=IRK-1
      DNI=DNI+3.65D2
      DH4=IRK/4.0D0-INT(IRK/4.0D0)
      DH100=IRK/1.0D2-INT(IRK/1.0D2)
      IF(DABS(DH4).LT.1.0D-6.AND.DH100.GT.0.0D0) DNI=DNI+1.0D0
      DH400=IRK/4.0D2-INT(IRK/4.0D2)
      IF(DABS(DH400).LT.1.0D-6) DNI=DNI+1.0D0
      IRK=IRK+2000
      GOTO 20
C
C   FROM 1583-JAN-1.0 TO 2000-JAN-1.0
  12  CONTINUE 
      DNI=DNI-2.2992385D6
      IRK=INT(DNI/3.65242198D2)
      DN0=365*IRK+INT((IRK+2)/4.0D0)-INT((IRK+82)/1.0D2)+INT((IRK+382)/4
     *.0D2)
      DNI=DNI-DN0
      IF(DNI.LT.1.0D0) GOTO 35
      IRK=IRK+1583
      GOTO 20
  35  IRK=IRK-1
      DH4=(IRK+3)/4.0D0-INT((IRK+3)/4.0D0)
      DH100=(IRK+83)/1.0D2-INT((IRK+83)/1.0D2)
      DNI=DNI+3.65D2
      IF(DABS(DH4).LT.1.0D-6.AND.DH100.GT.0.0D0) DNI=DNI+1.0D0
      DH400=(IRK+383)/4.0D2-INT((IRK+383)/4.0D2)
      IF(DABS(DH400).LT.1.0D-6) DNI=DNI+1.0D0
      IRK=IRK+1583
      GOTO 20
C
C   FROM 1582-OCT-15.0 TO 1583-JAN-1.0
  13  CONTINUE 
      IRK=1582
      DNI=DNI-2.2991605D6+1.4D1
      IF(DNI.GE.3.2D1) GOTO 40
      IMS=10
      GOTO 26
  40  CONTINUE
      IF(DNI.GE.6.2D1) GOTO 45
      DNI=DNI-3.1D1
      IMS=11
      GOTO 26
  45  DNI=DNI-6.1D1
      IMS=12
      GOTO 26
C
C   FROM 1-JAN-1.0 TO 1582-OCT-15.0
  14  CONTINUE
      DNI=DNI-1.7214235D6
      IRK=INT(DNI/3.6525D2)+1
      DN0=3.65D2*(IRK-1)+INT((IRK-1)/4.0D0)
      DNI=DNI-DN0
      IF(DNI.GE.1.0D0) GOTO 20
      IRK=IRK-1
      DH4=IRK/4.0D0-INT(IRK/4.0D0)
      DNI=DNI+3.65D2
      IF(DABS(DH4).LT.1.0D-6) DNI=DNI+1.0D0
      GOTO 20
C
C   FROM -4713-JAN-1.0 TO 1-JAN-1.0
  15  CONTINUE 
      DNI=DNI+0.5D0
      IRK=INT(DNI/3.6525D2)
      DN0=365*IRK+INT((IRK+2)/4.0D0)
      IF(DNI.GT.1.721118D6) DNI=DNI-1.0D0
      DNI=DNI-DN0
      IF(DNI.LT.1.0D0) GOTO 50
      IRK=IRK-4713
      GOTO 20
  50  IRK=IRK-1
      DH4=(IRK+3)/4.0D0-INT((IRK+3)/4.0D0)
      DNI=DNI+3.65D2
      IF(DABS(DH4).LT.1.0D-6) DNI=DNI+1.0D0
      IRK=IRK-4713
      GOTO 20
C
C   TO -4713-JAN-1.0
  16  CONTINUE
      DNI=DNI+1.9307115D6
      IRK=INT(DNI/3.6525D2)
      DN0=365*IRK+INT(IRK/4.0D0)
      DNI=DNI-DN0
      IF(DNI.LT.1.0D0) GOTO 55
      IRK=IRK-9999
      GOTO 20
  55  IRK=IRK-1
      DH4=(IRK+1)/4.0D0-INT((IRK+1)/4.0D0)
      DNI=DNI+3.65D2
      IF(DABS(DH4).LT.1.0D-6) DNI=DNI+1.0D0
      IRK=IRK-9999
C
  20  CONTINUE
      DH4=IRK/4.0D0-INT(IRK/4.0D0)
      DH100=IRK/1.0D2-INT(IRK/1.0D2)
      DH400=IRK/4.0D2-INT(IRK/4.0D2)
      IF(DABS(DH4).LT.1.0D-6) IMESI(3)=29
      IF(IRK.GT.1582.AND.DABS(DH100).LT.1.0D-6) IMESI(3)=28
      IF(IRK.GT.1582.AND.DABS(DH400).LT.1.0D-6) IMESI(3)=29
      IF(IRK.EQ.0) IRK=-1
      IMS=1
  23  CONTINUE
      IF(IMS.NE.13) GOTO 60
      IMS=12
      DNI=DNI-1.0D0
  60  CONTINUE
      IF(IMESI(IMS+1)+1.GT.DNI) GOTO 26
      DNI=DNI-IMESI(IMS+1)
      IMS=IMS+1
      GOTO 23
  26  CONTINUE
      IMESI(3)=28
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
