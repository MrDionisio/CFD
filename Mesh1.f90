program Mesh1
    implicit none
    integer i, Ny, Nt, IScheme, IO, k
    real, ALLOCATABLE :: UN(:), UN_1(:), Y(:)
    real :: h, V, A, dt, t, Time, dy

    allocate(UN(0:Ny+1), UN_1(0:Ny+1), Y(Ny))

    open(1,file='input.txt')
    read(1,*) h
    read(1,*) V
    read(1,*) Time 
    read(1,*) A 
    read(1,*) Ny
    read(1,*) IScheme
    close(1)


    print*, 'h', h, 'V', V, 'Time', Time,'A', A, 'N', Ny, IScheme

    dy=h/Ny 
    dt=dy**2/2/V
    Nt=INT(Time/dt) 

    print*, 'dy', dy, 'dt', dt, 'Nt' ,Nt

    Y(1)=0.0
    DO I=2, Ny-1
        Y(I)=Y(I-1)+dy
    END DO
    Y(Ny)=h

    UN=0.0
    
    
    t=0.0
    UN_1=UN 
    do k=1,Nt
        call MethodScheme(Ischeme, UN, UN_1, Ny, V, dt, dy, A)
        call Boundary(UN, Ny) 
        UN=UN_1
        t=t+dt
    end do

    do k=1,Ny
        UN_1(k)=A/2/V*((h/2)**2-(Y(k)-h/2)**2)
    end do

    IO=1
    open(io, file='Res.dat')
    call Output(IO, Y, UN, Ny)
    close(IO)    
    print*, 'Time', t, Nt
end program Mesh1



subroutine Boundary(U, Ny)
    integer Ny
    real, dimension(Ny) :: U
    U(1)=0
    U(Ny)=0 
end subroutine

subroutine MethodScheme(IScheme, U, U1, N, a, dt, dy, P)
    integer i, IScheme, N
    real a, dt, dy, P, A1, A2, A3
    real, dimension(0:N+1) :: U, U1
    if (IScheme==1) then
        do i=2,N-1
            U1(I)=U(i)+dt*(a/dy**2*(U(i-1)-2*U(i)+U(i+1))+P)
        end do
    end if
    if (IScheme==17) then
        do i=2,N-1
            A1 = a/dy**2*(U(i-1)-2*U(i)+U(i+1))
            A2 = (-a**2*dt/2+a*dy**2/12)
            A3 = (U(i+2)-4*U(i+1)+6*U(i)-4*U(i-1)+U(i-2))/dy**4
            U1(I)=U(i)+dt*(A1-A2*A3+P) 
        end do
        U1(0) = 2*U1(1)-U1(2)
        U1(N+1) = 2*U1(N)-U1(N-1)        
    end if
end subroutine



SUBROUTINE Output(IO, X, U, NX)
    IMPLICIT NONE
    integer IO, NX, I
    real, dimension(0:NX+1) :: X, U 
    DO i=1,NX
          WRITE(IO, *) X(i-1), U(i)
    END DO
END SUBROUTINE