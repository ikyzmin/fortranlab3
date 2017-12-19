
module svd
contains
	subroutine print_vec(vec)
		double precision, allocatable,dimension(:) :: vec
		n = size(vec)
		do i=1,n
			print *,vec(i)
		enddo
	end subroutine print_vec

	subroutine print_matr(matr,M,N,stream)
		double precision, allocatable,dimension(:,:) :: matr
		integer :: m,n,stream
	
		do i=1,M
			write(stream,*),(matr(i,j), j=1,N)
		enddo
	end subroutine print_matr

	subroutine diag(SIGMA,DIAG_SIGMA)
	double precision, allocatable,dimension(:) :: SIGMA
	double precision , allocatable,dimension(:,:) :: DIAG_SIGMA
	m = size(SIGMA)
	do i=1,m
	DIAG_SIGMA(i,i) = SIGMA(i)
	enddo
	end subroutine 

	subroutine diag_solve(DIAG_SIGMA)
	double precision, allocatable,dimension(:) :: SIGMA
	double precision , allocatable,dimension(:,:) :: U,DIAG_SIGMA
	m = size(DIAG_SIGMA,2)
	do i=1,m
	DIAG_SIGMA(i,i) = 1.0/DIAG_SIGMA(i,i)
	enddo
	end subroutine 
	
	function mean(X)
	double precision , allocatable,dimension(:,:) :: X
	double precision,allocatable,dimension(:) :: mean
	double precision :: tempMean
	m = size(X,1)
	n = size(X,2)
	allocate(mean(n))
	do i=1,n
	tempMean = 0
	do j=1,m
	tempMean = tempMean + X(j,i)
	enddo
	mean(i) = tempMean/m
	enddo
	end function mean

	function rep(X,ntimes)
	double precision, allocatable,dimension(:) :: X
	double precision , allocatable,dimension(:,:) ::rep
	integer :: ntimes
	n = size(X)
	allocate(rep(ntimes,n))
	do i=1,ntimes
	rep(i,:) =X 
	enddo
	end function rep

end module svd

 SUBROUTINE PRINT_MATRIX( DESC, M,N, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M,N, LDA
      double precision             A( LDA, * )

      INTEGER          I, J

      WRITE(*,*)
      WRITE(*,*) DESC
      DO I = 1, M
         print *,( A( I, J ), J = 1, N )
      END DO

      RETURN
      END

program lapack_lab
	use svd
	double precision,allocatable,dimension(:) :: X1,X2,SIGMA,ALPHA, BETA, mX
	double precision,allocatable,dimension(:,:) :: S,U,V,X,VT,NEW_S,INIT,DIAG_SIGMA, Q, B,VGSVD,A,CF,dY,Y1,U1,VT1
	integer :: N,M,P
	integer :: LDA,LWORK,LDTV,INFO, LDB, LDQ, LDV, K, L
	CHARACTER :: JOBU, JOBVT, JOBV, JOBQ
	DOUBLE PRECISION  :: wdum(1)
	!REAL :: WORK( 1000 )
	double precision,ALLOCATABLE,dimension(:) :: WORK
 	integer, Allocatable :: iwork(:)
	open (1, FILE='u.mat')
	open (2,FILE="a.mat")
	open (3,FILE = "bef.mat")
	open (4,FILE = "svd.mat")
	N =2
	M = 500
	JOBU='A'
      	JOBVT='A'
	p = m-1
	LDB = p
    	LDV = p
    	LDQ = n
      	LDA=M
      	LDU=M
      	LDVT=N
	allocate(X1(M),X2(M),A(M,N),CF(N,N),INIT(M,N),S(M,N),NEW_S(M,N),SIGMA(min(m,n)),U(M,M),V(LDVT,N),VT(LDVT,N),DIAG_SIGMA(M,N), B(LDB, n), BETA(n),Q(LDQ,n),VGSVD(LDV,p), ALPHA(n),iwork(n),dY(M,M),Y1(M,M),U1(M,M),VT1(M,M))
 	allocate(mX(2))
	call random_number(X1)
	call random_number(X2)
	call random_number(CF)
	print *,"X1 = "
	!call print_vec(X1)
	print *,"X2 = "
	!call print_vec(X2)
	S(:,1) = X1
	S(:,2) = X2
	!S(1,1)=1
	!S(2,1)=0
	!S(3,1)=4
	!S(1,2)=-2
	!S(2,2)=1
	!S(3,2)=4
	call print_matr(S,M,N,3)
	A = matmul(S,CF)
	print *,"coefficents = "
	!call print_matr(CF,N,N)

	print *,"S = "
	!call print_matr(S,M,N)
	print *,"S*COEF = "
	!call print_matr(A,M,N)
	
	!LWORK=MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N))
	!ALLOCATE(work(lwork))
	mX = mean(A)
	A = A-spread(mX,1,m);

	INIT = A
	call print_matr(A,M,N,2)
	LWORK=MAX(1,3*MIN(m,n)+MAX(m,n),5*MIN(m,n))
        ALLOCATE(work(lwork))
  do i=1,m-1
        B(i, 1) = -A(i, 1)+A(i+1, 1)
        B(i, 2) = -A(i, 2)+A(i+1, 2)
    end do
print *,"A = "
     ! CALL PRINT_MATR(A,M,N )
print *,"B = "
      !CALL PRINT_MATR(B,LDB,N)

CALL DGESVD(JOBU, JOBVT, M, N, A, LDA, SIGMA, U, LDU, VT, LDVT,WORK, LWORK,INFO)
	call diag(SIGMA,DIAG_SIGMA)
	call diag_solve(DIAG_SIGMA)
Y1 = matmul(INIT,matmul(VT,transpose(DIAG_SIGMA)))
print *,"B = "
dY = Y1
LWORK=MAX(1,3*m+m,5*m)
deallocate(work)
allocate(work(lwork))
CALL DGESVD(JOBU, JOBVT, M, M, dY, LDA, SIGMA, U1, M, VT1, M,WORK, LWORK,INFO)
Y1 = matmul(Y1,VT1);
print *,"B = "
call print_matr(Y1,M,M,4)


call DGGSVD3('U','V', 'Q', m, n, p, K, L, INIT, LDA, B, LDB, ALPHA, BETA, U, LDU, VGSVD, LDV, Q, LDQ,work, lwork, iwork, INFO)
 	IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm computing GSVD failed to converge.'
         STOP
      END IF
	V = transpose(VT)

print *,"AFTER GSVD A = "
      !CALL PRINT_MATR(A,M,N )
print *,"AFTER GSVD B = "
     ! CALL PRINT_MATR(B,LDB,N)
      call print_matr(U,M,M,1)
      close(1)
      close(2)
      close(3)
      close(4)
      !CALL PRINT_MATRIX( 'U= ', M, M, U, M )

      !CALL PRINT_MATRIX( 'V= ', LDV, P, VGSVD, LDV )

      !CALL PRINT_MATRIX( 'Q= ', LDQ, N, Q, LDQ )


	deallocate(X1,X2,A,S)
	
end program lapack_lab
