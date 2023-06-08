PROGRAM Convective
Implicit none

INTEGER, parameter:: IO = 12 ! input-output unit      
INTEGER NX,NT,I,J,ID,m, ISCHEME, ITASK, Burgers
REAL,ALLOCATABLE :: U(:),UN(:),X(:),UN1(:)
REAL L,h,CFL,dt,t,Time
REAL C, C0, C1, pi, k,G,FE, beta

WRITE(*,*) 'Read input file' 
OPEN(IO,FILE='Input.txt')
READ(IO,*) L
READ(IO,*) m
READ(IO,*) C
READ(IO,*) C0, C1
READ(IO,*) NX
READ(IO,*) Time
READ(IO,*) CFL
READ(IO,*) ISCHEME
READ(IO,*) ITASK
READ(IO,*) Burgers
CLOSE(IO)      

ALLOCATE(U(0:NX+1),UN(0:NX+1),X(0:NX+1), UN1(0:NX+1))

pi=3.14159265359d0
k = m*pi/L
h = 1./REAL(NX-1) 
dt = h*CFL/C  
NT = Time/dt

beta=k*h
G=sqrt((1-CFL*(1-cos(beta)))**2+(CFL*sin(beta))**2)
FE=-atan((CFL*sin(beta))/(1-CFL*(1-cos(beta))))
print*, 'G=', G, 'FE=', FE

WRITE(*,*) 'L=',L, 'h=', h, 'NX=', NX
WRITE(*,*) 'CFL=', CFL, 'dt=', dt, 'Time=', Time, 'NT=', NT

X(0)=-h
DO I=1, NX+1
      X(I)=X(I-1)+h
END DO

U(:)=0.0
UN(:)=0.0


CALL InitValue(U, X, k, NX, ITASK)
CALL BoundValue(U, NX)     

if (Burgers==0) then
!-------------------------  Solve equation ------------------
if (ISCHEME==1) then
      t=0.0d0
      DO J=1,NT
            DO I=1,NX
                  UN(I)=U(I)-dt/2/h*((C+abs(C))*(U(i)-U(i-1))+(C-abs(C))*(U(i+1)-U(i)))
            END DO
            UN(1)=UN(NX)
            U=UN
            t=t+dt
            CALL BoundValue(U, NX)
      end do

elseif (ISCHEME==8) then

      DO I=1,NX
            UN(I)=U(I)-dt/2/h*((C+abs(C))*(U(i)-U(i-1))+(C-abs(C))*(U(i+1)-U(i)))
      END DO
      UN(1)=UN(NX)
      CALL BoundValue(UN, NX)
      t=dt
      do j=2,NT
            do i=1,NX
                  UN1(i)=UN(i)-UN(i-1)+U(i-1)-2*dt/h*C*(UN(i)-UN(i-1))
            end do
            U=UN
            UN=UN1
            CALL BoundValue(UN, NX)
            t=t+dt
      end do 
end if
elseif(Burgers==1) then
      C=1
      if (ISCHEME==1) then
            t=0.0d0
            DO J=1,NT
                  DO I=1,NX
                        UN(I)=U(I)-dt/2/h*((C+abs(C))*(U(i)**2/2-U(i-1)**2/2)+(C-abs(C))*(U(i+1)**2/2-U(i)**2/2))
                  END DO
                  UN(1)=UN(NX)
                  U=UN
                  t=t+dt
                  CALL BoundValue(U, NX)
            end do
      
      
      elseif (ISCHEME==8) then
      
            DO I=1,NX
                  UN(I)=U(I)-dt/2/h*((C+abs(C))*(U(i)**2/2-U(i-1)**2/2)+(C-abs(C))*(U(i+1)**2/2-U(i)**2/2))
            END DO
            UN(1)=UN(NX)
            CALL BoundValue(UN, NX)
            t=dt
            do j=2,NT
                  do i=1,NX
                        UN1(i)=UN(i)-UN(i-1)+U(i-1)-2*dt/h*C*(UN(i)**2/2-UN(i-1)**2/2)
                  end do
                  U=UN
                  UN=UN1
                  CALL BoundValue(UN, NX)
                  t=t+dt
            end do 
      end if
end if

OPEN(IO,FILE='Res.dat')
      call Output(X, U, NX, IO)
CLOSE(IO)


do i=1,NX
      UN(i)= abs(G)**NT*SIN(k*X(I)+FE*NT)
      UN1(i) = sin(k*X(I)-k*NT*dt*C)
end do

OPEN(IO,FILE='Reschisl.dat')
call Output(X, U, NX, IO)
CLOSE(IO) 

OPEN(IO,FILE='Reschislan.dat')
call Output(X, UN, NX, IO)
CLOSE(IO)

OPEN(IO,FILE='Resan.dat')
call Output(X, UN1, NX, IO)
CLOSE(IO)

OPEN(IO, FILE='DissipativeForce.dat')
beta=0
CFL=0
do i=1,101
      CFL=(i-1)*0.01
      do j=1,101
      beta=(j-1)*0.03
      G=sqrt((1-CFL*(1-cos(beta)))**2+(CFL*sin(beta))**2)
      write(IO,*) CFL, beta, G
      end do
end do
CLOSE(IO)

END PROGRAM

!----------------------- Set Initial Value -----------------           
subroutine InitValue(U, X, k, N, ITASK)
IMPLICIT NONE
INTEGER :: i, N, ITASK
REAL, DIMENSION(0:N+1) :: U, X
REAL :: k

if (ITASK==1) then
      do i=1,N
            U(I)=SIN(k*X(I))
      END DO
elseif (ITASK==9) then
      do i=1,N
            if (X(i)<0.4) then
                  U(i)=1.0
            elseif ((X(i)>=0.4).and.(X(i)<=0.8)) then 
                  U(i) = 5*X(I)-3
            else 
                  U(i) = 1.0
            end if

      end do
end if

end subroutine

!----------------------- Set Boundary Condition ------------                  
subroutine BoundValue(U, N)
IMPLICIT NONE
INTEGER :: N
REAL, DIMENSION(0:N+1) :: U
U(0)=U(N-1)
U(N+1)=U(2)
U(1)=U(N)
end subroutine

!----------------------- Output results ------------                  
subroutine Output(X, U, NX, IO)
IMPLICIT NONE
integer IO, NX, I
real, dimension(0:NX+1) :: X
real, dimension(0:NX+1) :: U 
DO i=1,NX
      WRITE(IO, *) X(i), U(i)
END DO
end subroutine