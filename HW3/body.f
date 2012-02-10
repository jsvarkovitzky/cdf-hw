      Program main
      parameter (idm=500)     ! maximum number of grid points on airfoil
      real x(idm),y(idm)      ! x,y coordinates of grid points on airfoil.
                              ! airfoil is horizontally oriented with
                              ! leading edge on left, trailing edge on right
                              ! grid points start from lower trailing edge
                              ! and loops from bottom, through leading edge
                              ! to top trailing edge.
      integer nnode           ! number of grid points on airfoil
      integer naca            ! NACA 4-digit serier number, say, 0012, etc.
c
c     input nnode to subroutine body, the resulting airfoil coordinates wil
c     be put in the arrays x and y. Leading edge is at origin. Chord length
c     is normalized to be 1.
      
      nnode = 129
      call body(nnode,x,y)
      open(UNIT=10,FILE="body.dat")
      do i=1,nnode
        write(6,601) i, x(i),y(i)
        write(10,601) i, x(i),y(i)
      enddo
      stop
  601 format(1x/I3,2x,2E16.8)
      end
C
      SUBROUTINE BODY(nnode,x,y)
      integer nnode,npanel
      integer naca
      real x(nnode),y(nnode)
C
C     READ IN OR CALCULATE GRID POINTS on an airfoil
C
      npanel = nnode -1
      pi     = 4.*atan(1.)
      write(*,*)
      write(*,*) '       Choice of Airfoils '
      write(*,*) '1 --- NACA 4-digit Airfoil'
      write(*,*) '2 --- Your Camber and Thickness Definition'
      write(*,*) '3 --- Read in Data points from Airfoil.dat'
      write(*,*) ' Enter Your Choice:'
      read(*,*) IOPT
      IF(IOPT .EQ. 3) THEN
        open(13,file='airfoil.dat',form='formatted')
        read(13,*) NNODE
        DO 5 N=1,NNODE
        read(13,*) X(N),Y(N)
    5   CONTINUE
        GO TO 100
      ELSE IF(IOPT .EQ. 1) THEN
        write(*,*) 'ENTER THE 4-digit'
        read(*,*) NACA
        IF( NACA .GE. 10000) THEN
          write(*,*) 'More than 4-digit, STOP'
          stop
        ENDIF
        N1   = NACA/1000
        N2   = (NACA -N1*1000)/100
        N34  = NACA -N1*1000 -N2*100
        CMAX = FLOAT(N1)/100.
        PMAX = FLOAT(N2)/10.
        TMAX = FLOAT(N34)/100.
      ELSE 
        write(*,*) 'Not implemented yet, stop'
          stop
      ENDIF
      NL = NPANEL/2 + 1
      NLP = NL + 1
      DTH  = 2.*PI/FLOAT(NPANEL)
      DO 10 I = 1,NL
      XX     = 0.5 +.5*COS(DTH*(I-1))
      CALL CAMBER(CMAX,PMAX,XX,YCAM,CTH,STH)
      CALL  THICK(TMAX,XX,YT)
      XL     = XX   +YT*STH
      YL     = YCAM -YT*CTH
      XU     = XX   -YT*STH
      YU     = YCAM +YT*CTH
      X(I)   = XL
      Y(I)   = YL
      IU     = NNODE -I +1
      X(IU)  = XU
      Y(IU)  = YU
   10 CONTINUE
  100 CONTINUE
      X(NNODE) = .5*(X(1) +X(NNODE))
      Y(NNODE) = .5*(Y(1) +Y(NNODE))
      X(1)     = X(NNODE)
      Y(1)     = Y(NNODE)
      XMAX     = X(1)
      XMIN     = X(1)
      DO N=1,NNODE
        XMIN     = MIN(XMIN,X(N))
        XMAX     = MAX(XMAX,X(N))
      ENDDO
      CHORD      = XMAX -XMIN
      DO I=1,NNODE
        X(I)     = (X(I) -XMIN)/CHORD
        Y(I)     = (Y(I) -Y(1))/CHORD
      ENDDO
      XM         = .25
      YM         = Y(1)
      RETURN
      END
C
      SUBROUTINE THICK(TMAX,X,YT)
C
C     Calculate thickness for NACA 4-digit Airfoil
C
C     TMAX --- Maximum Thickness % chord
C
      X0    =  1.008930*X
      YT    =  (.29690*SQRT(X0)-.12600*X0-.35160*X0**2)
     .        +(.28430*X0**3-.10150*X0**4)
      YT    =  (TMAX/.20)*YT
      RETURN
      END
C
      SUBROUTINE CAMBER(CMAX,P,X,YC,CTH,STH)
C
C     Calculate Camber for NACA 4-digit Airfoil
C     CMAX --- Maximum Camber % chord
C     P    --- Location of Maximum Camber % chord
C     CTH  --- COS(TH)
C     STH  --- SIN(TH)
C
      IF(X.LT.P) THEN
        YC   = (CMAX/P**2)*(2.*P*X -X**2)
        YCP  = (CMAX/P**2)*2.*(P -X)
      ELSE
        YC   = (CMAX/(1. -P)**2)*((1. -2.*P) +2.*P*X -X**2)
        YCP  = (CMAX/(1. -P)**2)*2.*(P -X)
      ENDIF
      TH   = ATAN(YCP)
      CTH  = COS(TH)
      STH  = SIN(TH)
      RETURN
      END
