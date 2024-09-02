! Version: 2023-04-05
! program introgression.f90

! Uses periodic boundary conditions
! Reads file with barriers
! Number of individuals grows until N(t) >= N
! sex separation in vector s()
! g(i,j) = nuclear gene j of individual i
! gm(i,j) = mitochondrial gene j o individual i

! let R increases up to R + nrmax to find surrugate mother
! if no new mother is available within R+nrmax, individual dies.

! Initial mitochondrial genome is a mix of 0's and 1's 

! Marcus A.M. de Aguiar - Debora Princepe

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! module defining global variables
MODULE globals
INTEGER, SAVE :: iseed(33)
INTEGER(1), ALLOCATABLE, SAVE :: g(:,:),gm(:,:),s(:),fs(:,:)
INTEGER(1), ALLOCATABLE, SAVE :: gm0(:,:)
INTEGER, SAVE :: iibound
INTEGER, SAVE :: ndim,nf1,nf2,nb,nbm,mnpm,ineighborg,ineighbor,ineighbort,nrmax,iitime,ntime
INTEGER, SAVE :: igt,igt_old,igtfull,igtfull_old,igtm
INTEGER, SAVE :: inisbit,nc,nct,nco,ifs,ineighborden
INTEGER, SAVE :: iini,ifin,jini,jfin
INTEGER, ALLOCATABLE, SAVE :: x(:),y(:),neig(:),neigsp(:),neigspt(:),ispecies(:,:),ispeciesm(:,:)
INTEGER, ALLOCATABLE, SAVE :: worg(:,:,:),nworg(:,:)
INTEGER, ALLOCATABLE, SAVE :: ndq(:),idqx(:,:),idqy(:,:),ispv(:)
INTEGER, ALLOCATABLE, SAVE :: t(:,:),tp(:,:),p1(:),p2(:),ispidx(:),ispidxold(:)
INTEGER, ALLOCATABLE, SAVE :: vm(:),vm0(:)
INTEGER, SAVE :: max_species,fracrand, fullcount
INTEGER, ALLOCATABLE, SAVE :: ext_time_old(:),ext_time(:),t_sp(:,:),t_sp_old(:,:)
INTEGER, ALLOCATABLE, SAVE :: size_ext(:),size_ext_old(:),ispv_old(:)
INTEGER, ALLOCATABLE, SAVE :: sister_sp(:),sister_sp_old(:)
REAL, SAVE :: rg,qmat,aux,rho0,radius,mut,diff,mutmit,fmax,fdelta,width,inifrac
REAL, ALLOCATABLE, SAVE :: fitness(:),distmitonuc(:)

character(200),save     :: folder,path,filename
integer, save           :: counter
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module that reads and writes the populations
MODULE readwrite
INTEGER, SAVE :: iread
CHARACTER*25, SAVE :: name
CHARACTER*30, SAVE :: name1
CONTAINS
	
	SUBROUTINE FORBIDEN
	USE globals
		fs = 0  
		OPEN(UNIT=7,FILE='forbiden-sites.in',STATUS='OLD')
		ifs = 0
		DO 
			READ(7,*) i,j
			IF(i == -1) EXIT			
			i=i+1
			j=nf2-j
			fs(i,j) = 1
			ifs = ifs + 1  !number of forbiden sites
		END DO
		CLOSE(7)
	END SUBROUTINE FORBIDEN
    
	SUBROUTINE READPOP
	USE globals
	! initialize population
    
    IF(inisbit == 0) THEN ! random distribution
        nct = nco
        DO i=1,nct   
            ictr = 0
            DO WHILE (ictr == 0)
                CALL RANDOM_NUMBER(aux)
                ii = int(aux*nf1)+1
                CALL RANDOM_NUMBER(aux)
                jj = int(aux*nf2)+1					
                IF (fs(ii,jj) == 1) CYCLE
                x(i) = ii
                y(i) = jj
                ictr = 1					
            END DO
            CALL RANDOM_NUMBER(aux)
            s(i) = 1
            IF(aux < 0.5) s(i) = 0
            IF(fracrand == 1) THEN   ! 1's are randomly positioned
                DO jj=1,nbm
                    CALL RANDOM_NUMBER(aux)
                    IF(aux < inifrac) gm(i,jj) = 1
                END DO
            ELSE                     ! first loci are set to 1
                nones = int(inifrac*nbm)
                DO jj=1,nones
                    gm(i,jj) = 1
                END DO
            END IF
        END DO
    ELSE ! localized distribution
        nct = (ifin-iini)*(jfin-jini)*rho0
        i = 0			
        DO WHILE (i < nct)
            CALL RANDOM_NUMBER(aux)
            ii = iini + int(aux*(ifin-iini+1))
            CALL RANDOM_NUMBER(aux)
            jj = jini + int(aux*(jfin-jini+1))				
            IF (fs(ii,jj) == 1) CYCLE
            i = i + 1
            x(i) = ii
            y(i) = jj
            CALL RANDOM_NUMBER(aux)
            s(i) = 1
            IF(aux < 0.5) s(i) = 0
        END DO
    END IF

	END SUBROUTINE READPOP

	SUBROUTINE WRITEPOP
	USE GLOBALS
	name='new.dat'
	name1 = 'pop-'//name

	write(filename,trim('(a,a)'))trim(path),trim(name1)
	OPEN(UNIT=10,FILE=filename,STATUS='UNKNOWN',POSITION='REWIND')
	WRITE(10,*) iitime,nct,nf1,nf2
	WRITE(10,*) mut,mutmit,diff
	WRITE(10,*) radius,rg,nb,nbm,mnpm
	DO i=1,nct
		WRITE (10,901) x(i),y(i),(g(i,j),j=1,nb)
	END DO
    DO i=1,nct
        WRITE (10,902) (gm(i,j),j=1,nbm)
    END DO
    WRITE (10,903) (s(i),i=1,nct)

	close(10)

    901 FORMAT(i4,1x,i4,1x,200000i1)
    902 FORMAT(200000i1)
    903 FORMAT(200000(i1,1x))
	END SUBROUTINE WRITEPOP


    SUBROUTINE WRITEFITNESS(current_time)
    USE GLOBALS
    integer   :: current_time
        
    write(filename,trim('(a,"fitness_",I4.4".dat")'))trim(path),current_time
    OPEN(unit=5,file=filename,status='unknown')
    do ii=1,nct
        write(5,902) fitness(ii),distmitonuc(ii),ispidx(ii)
    end do
    CLOSE(5)
902 FORMAT(2(2x,f7.5),2x,i5)
    END SUBROUTINE WRITEFITNESS

END MODULE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module random seed
MODULE randomseed
IMPLICIT NONE
INTEGER :: iseed(12)
CONTAINS

    SUBROUTINE initialize
    USE globals, ONLY: aux,iseed
    OPEN(UNIT=50,FILE='seed.in',STATUS='OLD')  
	   READ(50,*) iseed
    CLOSE(50)
    CALL RANDOM_SEED(put=iseed)
    CALL RANDOM_NUMBER(aux)
    END SUBROUTINE initialize
    
    SUBROUTINE SAVESEED
    USE globals, ONLY: aux,iseed,filename,path
    
    open(unit=50,file='seed.in',status='old')
	read(50,*) iseed
	write(filename,trim('(a,"seed.dat")'))trim(path)
	open(10,file=filename)
	write(10,*) iseed
	close(10)
	close(50)
	
	CALL RANDOM_SEED(get=iseed)
    OPEN(UNIT=50,FILE='seed.in',STATUS='OLD', POSITION='REWIND')
        WRITE (50,*) iseed
    close(50)
    
    END SUBROUTINE SAVESEED  
END MODULE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! begin main program
PROGRAM introgression
USE globals
USE randomseed
USE readwrite

INTEGER(1), ALLOCATABLE :: gl(:,:),gml(:,:),sl(:)
INTEGER, ALLOCATABLE :: indxorg(:),xl(:),yl(:)
INTEGER xj(32),yj(32),ichoose(1)
INTEGER iix,iiy,ix1,iy1,idensity
REAL, ALLOCATABLE :: neighfit(:)

read(*,*) counter
write(folder, '( "run_",I2.2)') counter
call system('mkdir ' // trim(folder))
path = trim('./')//trim(folder)//trim('/')

! Read input data
! Note: nf1 >= nf2
OPEN(UNIT=7,FILE='input.in',STATUS='OLD',POSITION='REWIND')
READ(7,*) ntime,nco,nf1,nf2
READ(7,*) mut,mutmit,diff,qmat
READ(7,*) radius,rg,nb,nbm,mnpm,nrmax
READ(7,*) iread,inisbit,njump,nmax
READ(7,*) iini,ifin,jini,jfin
READ(7,*) width,inifrac,fracrand
CLOSE(7)

! diffusion table
CALL JUMPTABLE(xj,yj)

! initialize random number generator
CALL initialize

write(filename,trim('(a,"number0.dat")'))trim(path)
OPEN(UNIT=15,FILE=filename,STATUS='unknown')

! open 'count mutations' files
write(filename,trim('(a,"count-mut.dat")'))trim(path)
OPEN(unit=21,file=filename, STATUS='unknown')


! model parameters and initializations
idensity = 15   ! maximum number of individuals per site before  mating is idensity/3
rho0 = float(nco)/float(nf1*nf2-ifs)
nc = 1.2*nco  !allocate vectors with a margin for nct > nco

ALLOCATE (fs(nf1,nf2))
CALL FORBIDEN

ALLOCATE (g(nc,nb),gl(nc,nb))
ALLOCATE (gm(nc,nbm),gml(nc,nbm),s(nc),sl(nc))
ALLOCATE (x(nc),y(nc),xl(nc),yl(nc))
ALLOCATE (neig(nc),neigsp(nc),neigspt(nc))
ALLOCATE (ispv(nc),ispecies(nc,nc))
ALLOCATE (worg(nf1,nf2,idensity),nworg(nf1,nf2),indxorg(nc))
ALLOCATE (t(nc,nc),tp(nc,nc),p1(nc),p2(nc),ispidx(nc),ispidxold(nc))
ALLOCATE (fitness(nc),distmitonuc(nc))
ALLOCATE (neighfit(nc))

max_species = 15000

ALLOCATE (ext_time(max_species),ext_time_old(max_species))
ALLOCATE (t_sp(max_species,max_species),t_sp_old(max_species,max_species))
ALLOCATE (ispv_old(nc),size_ext(max_species),size_ext_old(max_species))
ALLOCATE (sister_sp(max_species),sister_sp_old(max_species))

g = 0  ! identical genomes
gl = 0
gm = 0 ! identical mitochondria
gml = 0
x = 0
y = 0
xl = 0
yl = 0
t = 0  ! initial times to common ancestor
tp = 0 ! aux variables
p1 = 0 ! first parent
p2 = 0 ! second parent
ispidx = 1
igt_old = 1
size_ext = 0  ! population sizes at time of extinction
sister_sp = 0 ! sister species for extinct species
sister_sp_old = 0

! initialize populations
CALL READPOP

! save initial mitochondrial genomes for substitution rate
ALLOCATE (gm0(nc,nbm),vm(nc),vm0(nc))
gm0 = gm
vm = 0
DO i=1,nct
    vm0(i) = i
END DO

WRITE(6,*) 'total area =',nf1*nf2
WRITE(6,*) 'blocked area =',ifs
WRITE(6,*) 'available area =',nf1*nf2-ifs
WRITE(6,*) 'average density =', rho0
WRITE(6,*) 'average number of individuals in S =',3.1416*rho0*radius**2
WRITE(6,*) 'inital number of individuals =',nct
WRITE(6,*) 'total integration time =',ntime
WRITE(6,*) 'width of selection =',width

! Place the organisms by location in the poster space
nworg =0
worg = 0
indxorg = 0
DO i=1,nct     
	nworg(x(i),y(i)) = nworg(x(i),y(i)) + 1
	worg(x(i),y(i),nworg(x(i),y(i))) = i
	indxorg(i) = nworg(x(i),y(i))
END DO

! Set up template for the locations that are part of the neighborhood for a radius
ndim = 4*(radius+nrmax)**2
ALLOCATE (idqx(nrmax,ndim),idqy(nrmax,ndim),ndq(nrmax))
idqx = 0
idqy = 0
ndq = 0
DO irad = 1,nrmax 
	idq = 0
	rs = radius + irad - 1
	rs2 = rs**2
	nj = int(rs)+1
	DO itx=-nj,nj               
		DO ity=-nj,nj      
			IF(itx*itx+ity*ity <= rs2) THEN
				idq = idq + 1
				idqx(irad,idq) = itx
				idqy(irad,idq) = ity
			END IF
		END DO
	END DO
	ndq(irad) = idq
END DO

write(6,*) 'initial pop=',nct
write(6,*)
write(6,*) 'partial time   total time   total pop.  species   ind./species' 
write(6,*)

!
CALL POPFITNESS
write(6,*) 'initial average fitness =',sum(fitness)/float(nct)

!
! Time evolution: mating, mutation and diffusion
!
ntest = nco
IF (nmax /= 0) ntest = nmax
DO j=1,ntime
	!Mating
    knext = 0   ! count individuals of the next generation
    CALL POPFITNESS
	looppop: DO k=1,nct
		CALL FINDNEIG(k,1)
        nlocal = rho0*iibound*0.7
        ! check if extra offspring is produced
        iextra = 0
        IF(nct < ntest .AND. ineighbort < nlocal) THEN
            nthere = nworg(x(k),y(k))
            IF(ineighborg >= mnpm .AND. nthere < idensity/3) THEN
                iextra = 1
            END IF
        END IF

        IF(iextra == 1) THEN
            knext = knext + 1
            !
            ! choose mate according to fitness
            neighfit = 0.0
            neighfit(1) = fitness(neig(1))
            DO i=2,ineighborg
                neighfit(i) = neighfit(i-1)+fitness(neig(i))
            END DO
            CALL RANDOM_NUMBER(aux)
	        IF(neighfit(ineighborg) /= 0.0) THEN
		        neighfit = neighfit/neighfit(ineighborg)
		        ichoose = minloc(neighfit-aux,mask=(neighfit - aux > 0))
	        ELSE
		        ichoose = int(aux*ineighborg) + 1
	        END IF
            kmate = neig(ichoose(1))
            ! end choosing mate
            !
            DO kc=1,nb
                CALL RANDOM_NUMBER(aux)
                IF(aux < 0.5) THEN
                    gl(knext,kc) = g(k,kc)
                ELSE
                    gl(knext,kc) = g(kmate,kc)
                END IF
                CALL RANDOM_NUMBER(aux)
                IF(aux < mut) gl(knext,kc) = 1-gl(knext,kc)   ! Mutation
            END DO
            ifemale = k
            imale = kmate
            IF(s(k) == 1) THEN
                ifemale = kmate
                imale = k
            END IF
            gml(knext,1:nbm) = gm(ifemale,1:nbm)
            DO kc=1,nbm
                CALL RANDOM_NUMBER(aux)
                IF(aux < mutmit) gml(knext,kc) = 1-gml(knext,kc)   ! Mutation
            END DO
            xl(knext) = x(k)
            yl(knext) = y(k)
            sl(knext) = 1
            vm(knext) = vm0(ifemale)
            CALL RANDOM_NUMBER(aux)
            IF(aux < 0.5) sl(knext) = 0
            ! save parents
            p1(knext) = ifemale
            p2(knext) = imale

            ! regular offspring - same pair k-kmate
            knext = knext + 1
            DO kc=1,nb
                CALL RANDOM_NUMBER(aux)
                IF(aux < 0.5) THEN
                    gl(knext,kc) = g(k,kc)
                ELSE
                    gl(knext,kc) = g(kmate,kc)
                END IF
                CALL RANDOM_NUMBER(aux)
                IF(aux < mut) gl(knext,kc) = 1-gl(knext,kc)   ! Mutation
            END DO
            gml(knext,1:nbm) = gm(ifemale,1:nbm)
            DO kc=1,nbm
                CALL RANDOM_NUMBER(aux)
                IF(aux < mutmit) gml(knext,kc) = 1-gml(knext,kc)   ! Mutation
            END DO
            xl(knext) = x(k)
            yl(knext) = y(k)
            sl(knext) = 1
            vm(knext) = vm0(ifemale)
            CALL RANDOM_NUMBER(aux)
            IF(aux < 0.5) sl(knext) = 0
            ! save parents
            p1(knext) = ifemale
            p2(knext) = imale
		ELSE
			! only regular offspring
			irad = 1
			DO WHILE(ineighbor < 2)
				irad = irad + 1
				IF (irad == nrmax + 1) THEN
                    CYCLE looppop
				END IF
				CALL FINDNEIG(k,irad)
			END DO
			kmother = k 
			irad = 1
			CALL FINDMATE(k,kmother,kmate,irad)
			IF(ineighborg < 2) THEN
				CYCLE looppop
			ELSE
                knext = knext + 1
                DO kc=1,nb
                    CALL RANDOM_NUMBER(aux)
                    IF(aux < 0.5) THEN
                        gl(knext,kc) = g(kmother,kc)
                    ELSE
                        gl(knext,kc) = g(kmate,kc)
                    END IF
                    CALL RANDOM_NUMBER(aux)
                    IF(aux < mut) gl(knext,kc) = 1-gl(knext,kc)   ! Mutation
                END DO
                ifemale = kmother
                imale = kmate
                IF(s(kmother) == 1) THEN
                    ifemale = kmate
                    imale = kmother
                END IF

                gml(knext,1:nbm) = gm(ifemale,1:nbm)
                DO kc=1,nbm
                    CALL RANDOM_NUMBER(aux)
                    IF(aux < mutmit) gml(knext,kc) = 1-gml(knext,kc)   ! Mutation
                END DO
                xl(knext) = x(k)
                yl(knext) = y(k)
                sl(knext) = 1
                vm(knext) = vm0(ifemale)
                CALL RANDOM_NUMBER(aux)
                IF(aux < 0.5) sl(knext) = 0
                ! save parents
                p1(knext) = ifemale
                p2(knext) = imale
			END IF
		END IF
	END DO looppop

	!update number of individuals, genomes and sex
    nct = knext
	g = gl
	gl = 0
    gm = gml
    gml = 0
    s = sl
    sl = 0
    vm0 = vm

    ! update time matrix
    do ii=1,nct-1
        do jj=ii+1,nct
            it11 = t(p1(ii),p1(jj))
            it12 = t(p1(ii),p2(jj))
            it21 = t(p2(ii),p1(jj))
            it22 = t(p2(ii),p2(jj))
            itlm = min(it11,it12,it21,it22)
            tp(ii,jj) = itlm + 1
            tp(jj,ii) = itlm + 1
        end do
    end do
    t = tp
    tp = 0

	!Diffusion
	IF(diff /= 0.0) THEN
		DO i=1,nct 
			CALL RANDOM_NUMBER(aux)
			IF(aux < diff) THEN
				ibound = 0	
				DO WHILE(ibound == 0)
					CALL RANDOM_NUMBER(aux)
					jjump = INT(njump*aux)+1
					iix = xl(i)+xj(jjump)
					iiy = yl(i)+yj(jjump)
					IF (iix > nf1) iix = iix - nf1
					IF (iix < 1) iix = nf1 + iix
					IF (iiy > nf2) iiy = iiy - nf2
					IF (iiy < 1) iiy = nf2 + iiy
					IF (fs(iix,iiy) == 1) CYCLE
					xl(i) = iix
					yl(i) = iiy
					ibound = 1
				END DO
			END IF	
		END DO
	END IF

    x = xl
    y = yl


    ! Place the organisms by location in the poster space
    nworg =0
    worg = 0
    indxorg = 0
    DO i=1,nct
        nworg(x(i),y(i)) = nworg(x(i),y(i)) + 1
        worg(x(i),y(i),nworg(x(i),y(i))) = i
        indxorg(i) = nworg(x(i),y(i))
    END DO

! calculate species every generation starting at 2
	IF (j >= 2) THEN
        ispidxold = ispidx
        igtfull_old = igtfull
        ext_time_old = ext_time
        t_sp_old = t_sp
        ispv_old = ispv
        size_ext_old = size_ext
        sister_sp_old = sister_sp
        CALL FINDSPECIES
        CALL POPFITNESS
        write(15,*) j+iitime, igt, sum(distmitonuc)/float(nct), sum(fitness)/float(nct)
        write(6,*) j,j+iitime,nct,igt,(ispv(k),k=1,igt)

        CALL COUNTMUT
        write(6,*)
	END IF

END DO  ! end loop in time
CLOSE(15)
CLOSE(21)


CALL FINDSPECIES
write(6,*)
write(6,*) 'species'
write(6,100) igt
write(6,100) (ispv(k),k=1,igt)

CALL SUBSTITUTION_RATE

100 FORMAT(50(i4,2x))


iptime = ntime
iitime = iitime + iptime
CALL WRITEPOP
CALL SAVESEED

END PROGRAM introgression


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE JUMPTABLE(xj,yj)
integer xj(32),yj(32)
xj(1)=0;    yj(1)=1;
xj(2)=1;    yj(2)=0;
xj(3)=0;    yj(3)=-1;
xj(4)=-1;   yj(4)=0;  !first 4 neighbors, dist-max = 1

xj(5)=1;    yj(5)=1;
xj(6)=1;    yj(6)=-1;
xj(7)=-1;   yj(7)=-1;
xj(8)=-1;   yj(8)=1;
xj(9)=0;    yj(9)=2;
xj(10)=2;   yj(10)=0;
xj(11)=0;   yj(11)=-2;
xj(12)=-2;  yj(12)=0;
xj(13)=1;   yj(13)=2;
xj(14)=2;   yj(14)=1;
xj(15)=2;   yj(15)=-1;
xj(16)=1;   yj(16)=-2;
xj(17)=-1;  yj(17)=-2;
xj(18)=-2;  yj(18)=-1;
xj(19)=-2;  yj(19)=1;
xj(20)=-1;  yj(20)=2; !first 20 neighbors, dist-max = 2

xj(21)=0;   yj(21)=3;
xj(22)=3;   yj(22)=0;
xj(23)=0;   yj(23)=-3;
xj(24)=-3;  yj(24)=0;
xj(25)=1;   yj(25)=3;
xj(26)=3;   yj(26)=1;
xj(27)=3;   yj(27)=-1;
xj(28)=1;   yj(28)=-3;
xj(29)=-1;  yj(29)=-3;
xj(30)=-3;  yj(30)=-1;
xj(31)=-3;  yj(31)=1;
xj(32)=-1;  yj(32)=3; !first 32 neighbors, dist-max = 3

END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For each individual k, find a kmother and its mate      !
! Mother is usually k, but not always. 
! At the beginning NEIG was called with irad = 1.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDMATE(k,kmother,kmate,irad)
USE globals
INTEGER, INTENT(IN) :: k
INTEGER, INTENT(OUT) :: kmate
INTEGER, INTENT(INOUT) :: kmother
INTEGER ichoose(1),imother,ineighbork
REAL, ALLOCATABLE :: neighfit(:)

ALLOCATE (neighfit(nc))
neighfit = 0.0

! choose a new mother if aux < fitness or if potential mates < P
! the choice may involve increasing irad until new mother has at least 2 potential mates
!
CALL RANDOM_NUMBER(aux)
qnorm = qmat
IF(fdelta /= 0.0) qnorm = qmat*(1+2*(fdelta-fitness(k)))

IF(aux < qnorm .OR. ineighborg < mnpm) THEN
    ! CHOOSE MOTHER ACCORDING TO FITNESS
    ! compute cumulative normalized fitness vector
    neighfit(1) = fitness(neigspt(1))
    DO i=2,ineighbort
        neighfit(i) = neighfit(i-1)+fitness(neigspt(i))
    END DO
    neighfit = neighfit/neighfit(ineighbort)
    CALL RANDOM_NUMBER(aux)
    ichoose = minloc(neighfit-aux,mask=(neighfit - aux > 0))
    kmother = neigspt(ichoose(1))  ! choose a new mother
    CALL FINDNEIG(kmother,irad) ! find mother's neighbors
    icc = 0
    DO WHILE (ineighborg < 2 .AND. irad < nrmax)  ! choose a neighbor mother with ineighborg > 1
        IF(icc > 1) THEN
            irad = irad + 1  ! increase search radius
            icc = 0
            CALL FINDNEIG(k,irad)
            ! compute cumulative normalized fitness vector
            IF(ineighbort /= 0) THEN
                neighfit(1) = fitness(neigspt(1))
                DO i=2,ineighbort
                    neighfit(i) = neighfit(i-1)+fitness(neigspt(i))
                END DO
		        fittest = neighfit(ineighbort)
		        IF(fittest /= 0.0) THEN
                    neighfit = neighfit/neighfit(ineighbort)
		        END IF
            END IF
        END IF
        CALL RANDOM_NUMBER(aux)  !try a different neighbor within the same irad if icc=0
	    IF(fittest /= 0.0) THEN
        	ichoose = minloc(neighfit-aux,mask=(neighfit - aux > 0))
	    ELSE
		    ichoose = int(aux*ineighbort) + 1
	    END IF
        kmother = neigspt(ichoose(1))
        ! if there are no spatial neighbors kmother=0 and need to increase neighborhood size
        IF(kmother /= 0) THEN
            CALL FINDNEIG(kmother,irad) ! find mother's neighbors
        END IF
        icc = icc + 1
    END DO
    IF(irad > nrmax) WRITE(6,*) 'picking last individual in search'
END IF


! choose a mate for the mother
!
IF (ineighborg >= 2) THEN
    ! choose mate according to fitness
    neighfit = 0.0
    neighfit(1) = fitness(neig(1))
    DO i=2,ineighborg
        neighfit(i) = neighfit(i-1)+fitness(neig(i))
    END DO
    CALL RANDOM_NUMBER(aux)
    IF(neighfit(ineighborg) /= 0.0) THEN
	    neighfit = neighfit/neighfit(ineighborg)
	    ichoose = minloc(neighfit-aux,mask=(neighfit - aux > 0))
    ELSE
	    ichoose = int(ineighborg*aux) + 1
    END IF
    kmate = neig(ichoose(1))
    ! end choosing mate
END IF

END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Find spatial and genetic neighbors of kmother !
!   only neighbors of oposite sex are selected   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDNEIG(kmother,irad)
USE globals
INTEGER, INTENT(IN) :: kmother

ineighbor = 0
ineighborg = 0
ineighbort = 0
iibound = 0
neig = 0    ! genetic neighbors
neigsp = 0  ! spatial neighbors of oposite sex
neigspt = 0 ! total spatial neighbors
! the focal individual kmother is not included as a neighbor of itself
ix = x(kmother)
iy = y(kmother)
loop1: DO isite = 1,ndq(irad) 
	ix1 = ix+idqx(irad,isite)
	iy1 = iy+idqy(irad,isite)
	IF (ix1 > nf1) ix1 = ix1 - nf1
	IF (ix1 < 1) ix1 = nf1 + ix1
	IF (iy1 > nf2) iy1 = iy1 - nf2
	IF (iy1 < 1) iy1 = nf2 + iy1
    IF (fs(ix1,iy1) == 1) CYCLE loop1
	iibound = iibound + 1 ! number of actual sites in the neighborhood
	loop2: DO iworg = 1,nworg(ix1,iy1)
        korg = worg(ix1,iy1,iworg)
        if(korg /= kmother) then
            ineighbort = ineighbort + 1
            neigspt(ineighbort) = korg
            if(s(korg) /= s(kmother)) then
                ineighbor = ineighbor + 1
                neigsp(ineighbor) = korg
                dista = 0
                DO l=1,nb
                    dista = dista + ABS(g(kmother,l)-g(neigsp(ineighbor),l))
                    IF(dista > rg) CYCLE loop2
                END DO
                ineighborg = ineighborg + 1
                neig(ineighborg) = neigsp(ineighbor)
            end if
        end if
	END DO loop2
END DO loop1

END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find species using nuclear chromosome         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDSPECIES
USE globals, ONLY: nct,nb,rg,g,igt,ispecies,ispv,ispidx
INTEGER, ALLOCATABLE :: species(:),auxy1(:),auxy2(:)
INTEGER iix,iiy

ALLOCATE (species(nct),auxy1(nct),auxy2(nct))

itot = 0  ! count total population in groups
igt = 0   ! count number of groups
i2 = nct
DO i=1,i2  ! initialize aux2 -- contains all individuals not yet classified
	auxy2(i) = i
END DO 
	
DO WHILE (itot < nct)

	icr = auxy2(1)   !take first individual and find its species
	isp = 0
	ispold = 1
	i1 = 0
	auxy1 = 0
	loop1: DO i=1,i2
		ii = auxy2(i)
		dista = 0
		DO l=1,nb
			IF(g(icr,l) /= g(ii,l)) dista = dista + 1  
			IF(dista > rg) THEN
				i1 = i1 + 1
				auxy1(i1) = ii      !put creatures with dist > rg into aux1
				CYCLE loop1
			END IF
		END DO
		isp = isp + 1
		species(isp) = ii   !collect individuals with dist <= rg from icr
	END DO loop1

    !check if individuals in aux1 have to be included; put the rest in aux2
	itest = 1
	DO WHILE(itest /= 0)
		i2 = 0
		auxy2 = 0
		itest = 0
		isp0 = isp
		IF(i1 /= 0) THEN
			loop2:	DO i=1,i1
			DO ji=ispold+1,isp0  
				dista = 0
				DO l=1,nb
					IF(g(auxy1(i),l) /= g(species(ji),l)) dista = dista + 1  !HERE
					IF(dista > rg) EXIT
				END DO
				IF(dista <= rg) THEN
					isp = isp + 1 
					species(isp) = auxy1(i)   ! colect the aux1 individual
					itest = 1                 ! indicates that the process has to be repeated
					CYCLE loop2
				END IF
			END DO
			i2 = i2 + 1
			auxy2(i2) = auxy1(i)  ! put individual in aux2
			END DO loop2
		END IF
		auxy1 = auxy2   ! aux1 contains the creatures not in the species
		i1 = i2
		ispold = isp0
	END DO

	itot = itot + isp    !total number of individuals classified into species
	igt = igt + 1        !number of species

	! save species info
	DO i=1,isp
		ispecies(igt,i) = species(i)
        ispidx(species(i)) = igt
	END DO
	ispv(igt) = isp          ! number of individuals in species

END DO

END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the fitness of kmother's neighbors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE POPFITNESS
USE globals

fitness = 0.0
anbm = float(nbm)

w2 = width**2
DO k=1,nct
    dista = 0
    DO l=1,nbm
        dista = dista + ABS(g(k,l)-gm(k,l))
    END DO
    dista = dista/anbm
    distmitonuc(k) = dista
    if (width==0) then
        fitness(k) = 1.0
    else
        fitness(k) = exp(-0.5*dista**2/w2)
    end if
END DO


fmax = maxval(fitness, mask=fitness > 0)
fmin = minval(fitness, mask=fitness > 0)
fdelta = (fmax+fmin)/2

END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Count mutations at every generation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE COUNTMUT
USE globals

INTEGER, ALLOCATABLE :: mutmit1(:),mutmit2(:),mutnuc1(:),mutnuc2(:),mutnucneutral(:)
INTEGER :: sub01,sub10

ALLOCATE (mutmit1(nct),mutmit2(nct),mutnuc1(nct),mutnuc2(nct),mutnucneutral(nct))

mutnuc1 = 0; mutnuc2 = 0
mutmit1 = 0; mutmit2 = 0
mutnucneutral = 0

nones = int(inifrac*nbm)
DO j=1,nct
    DO k=1,nones
        mutnuc1(j) = mutnuc1(j) + g(j,k)
    END DO
    DO k=nones+1,nbm
        mutnuc2(j) = mutnuc2(j) + g(j,k)
    END DO
    DO k=nbm+1,nb
        mutnucneutral(j) = mutnucneutral(j) + g(j,k)
    END DO
    DO k=1,nones
        mutmit1(j) = mutmit1(j) + gm(j,k)
    END DO
    DO k=nones+1,nbm
        mutmit2(j) = mutmit2(j) + gm(j,k)
    END DO
END DO

sub10 = 0.0; sub01 = 0.0
DO i=1,nct
    DO j=1,nbm
        gm0ij = gm0(vm(i),j)
        IF(gm(i,j) /= gm0ij) THEN
		IF(gm0ij == 0) THEN
			sub01 = sub01 + 1.0
		ELSE
			sub10 = sub10 + 1.0
		END IF
	END IF
    END DO
END DO

write(21,*) sum(mutmit1), sum(mutmit2), sum(mutnuc1), sum(mutnuc2), sum(mutnucneutral), sub10, sub01

END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the average substitution rate per locus per unit time  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SUBSTITUTION_RATE
USE globals
subrate10 = 0.0
subrate01 = 0.0
DO i=1,nct
    DO j=1,nbm
        gm0ij = gm0(vm(i),j)
        IF(gm(i,j) /= gm0ij) THEN
		IF(gm0ij == 0) THEN
			subrate01 = subrate01 + 1.0
		ELSE
			subrate10 = subrate10 + 1.0
		END IF
	END IF
    END DO
END DO
subrate = subrate01 + subrate10

subrate01 = subrate01/((1.0-inifrac)*float(nct*nbm*ntime))
subrate10 = subrate10/(inifrac*float(nct*nbm*ntime))
subrate = subrate/float(nct*nbm*ntime)



write(filename,trim('(a,"survivals.dat")'))trim(path)
OPEN(unit=35,file=filename)
DO i=1,nct
    WRITE(35,*) vm(i)
END DO
END

