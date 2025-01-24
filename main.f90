MODULE functions
IMPLICIT NONE
! All functions needed to solve the ionization state linear system that we consider
CONTAINS

REAL*8 FUNCTION f_alpha_h_p(T)
IMPLICIT NONE
REAL*8:: T

f_alpha_h_p = 8.4e-11 * T**(-0.5d0) * (T/1000.d0)**(-0.2d0) * 1/(1+(T/1e6)**0.7d0)

END FUNCTION f_alpha_h_p

REAL*8 FUNCTION f_alpha_he_p(T)
IMPLICIT NONE 
REAL*8:: T

f_alpha_he_p = 1.5e-10 * T**(-0.6353d0)

END FUNCTION f_alpha_he_p

REAL*8 FUNCTION f_alpha_d(T)
IMPLICIT NONE
REAL*8:: T

f_alpha_d = 1.9e-3 * T**(-1.5d0)* EXP(-470000.d0/T) * (1+0.3*EXP(-94000.d0/T))

END FUNCTION f_alpha_d

REAL*8 FUNCTION f_alpha_he_pp(T)
IMPLICIT NONE
REAL*8:: T

f_alpha_he_pp = 3.36e-10 * T**(-0.5) * (T/1000.d0)**(-0.2) * 1/(1+(T/1e6)**0.7)

END FUNCTION f_alpha_he_pp

REAL*8 FUNCTION f_gamma_h0(T)
IMPLICIT NONE
REAL*8:: T

f_gamma_h0 = 5.85e-11 * SQRT(T) * EXP(-157809.1d0/T) * 1/(1+SQRT(T/100000.d0))

END FUNCTION f_gamma_h0

REAL*8 FUNCTION f_gamma_he0(T)
IMPLICIT NONE
REAL*8:: T

f_gamma_he0 = 2.38e-11 * SQRT(T) * EXP(-285335.4d0/T) * 1/(1+SQRT(T/100000.d0))

END FUNCTION f_gamma_he0

REAL*8 FUNCTION f_gamma_he_p(T)
IMPLICIT NONE
REAL*8:: T
f_gamma_he_p = 5.68e-12 * SQRT(T) * EXP(-631515.d0/T) * 1/(1+SQRT(T/100000.d0))
END FUNCTION f_gamma_he_p

END MODULE functions

MODULE ionization
USE functions
IMPLICIT NONE
INTEGER:: j,k
REAL*8, PUBLIC:: states(6), X, Y, yhe, set(2), n_H 
REAL*8, PARAMETER, PUBLIC:: mp = 1.67e-24, gammaa = 5.d0/3.d0, k_boltz = 1.380649e-16
CHARACTER(4), PUBLIC:: x_name
CHARACTER(5), PUBLIC:: h_name

CONTAINS 

FUNCTION set_state()
! This function sets X and n_H from user input
IMPLICIT NONE
REAL*8:: set_state(2)
X = 0.d0
PRINT*, "Input fraction of Hydrogen mass X: "
READ(*,*) X 

PRINT*, "Input total numeric hydrogen density n_H: " 
READ(*,*) n_H

set_state(1) = X
set_state(2) = n_H
WRITE(x_name,'(f4.2)')X
WRITE(h_name,'(f5.2)')n_H
Y = 1 - X 
yhe = Y / (4*(1-Y))

END FUNCTION

FUNCTION ionization_state(T)
! Incorporated the linear system resolution into a unique function of T, so I can call n_x(T)
IMPLICIT NONE
REAL*8, DIMENSION(6):: ionization_state
REAL*8:: incognite(6, 6), t_noti(6)
REAL*8:: alpha_h_p, alpha_he_p, alpha_d, alpha_he_pp, gamma_h0, gamma_he0, gamma_he_p
REAL*8, INTENT(in):: T


IF (X==0) THEN
set = set_state()
X = set(1)
n_H = set(2)
END IF

! Ionization parameters as function of temperature:
alpha_h_p = f_alpha_h_p(T)
alpha_he_p = f_alpha_he_p(T)
alpha_he_pp = f_alpha_he_pp(T)
alpha_d= f_alpha_d(T)
gamma_h0 = f_gamma_h0(T)
gamma_he0= f_gamma_he0(T)
gamma_he_p = f_gamma_he_p(T)


DO k = 1, 6
  t_noti(k) = 0.d0
  DO j = 1, 6
  incognite(k,j) = 0.d0 
  END DO
ENDDO

! Defining the matrices to solve
incognite(1,1) = gamma_h0
incognite(1,2) = -alpha_h_p
incognite(2,3) = gamma_he0
incognite(2,4) = -alpha_he_p - alpha_d
incognite(3,5) = alpha_he_pp
incognite(3,4) = -gamma_he_p
incognite(4,1) = 1.d0
incognite(4,2) = 1.d0
incognite(5,2) = 1.d0
incognite(5,4) = 1.d0
incognite(5,5) = 2.d0
incognite(5,6) = -1.d0
incognite(6,3) = 1.d0
incognite(6,4) = 1.d0
incognite(6,5) = 1.d0
t_noti(4) = 1
t_noti(6) = Y / 4*(1-Y)

CALL solve_linear_system(incognite, t_noti, states, 6)

ionization_state = states

END FUNCTION ionization_state


! Functions that define the ionization state at any T
REAL*8 FUNCTION n_H0(T)
IMPLICIT NONE
REAL*8, INTENT(in):: T
REAL*8, DIMENSION(6):: s

s = ionization_state(T)
n_H0 = s(1) * n_H

END FUNCTION n_H0

REAL*8 FUNCTION n_H_p(T)
IMPLICIT NONE
REAL*8, INTENT(in):: T
REAL*8, DIMENSION(6):: s

s = ionization_state(T)
n_H_p = s(2) * n_H

END FUNCTION n_H_p

REAL*8 FUNCTION n_He0(T)
IMPLICIT NONE
REAL*8, INTENT(in):: T
REAL*8, DIMENSION(6):: s

s = ionization_state(T)
n_He0 = s(3) * n_H

END FUNCTION n_He0

REAL*8 FUNCTION n_He_p(T)
IMPLICIT NONE
REAL*8, INTENT(in):: T
REAL*8, DIMENSION(6):: s

s = ionization_state(T)
n_He_p = s(4) * n_H

END FUNCTION n_He_p

REAL*8 FUNCTION n_He_pp(T)
IMPLICIT NONE
REAL*8, INTENT(in):: T
REAL*8, DIMENSION(6):: s

s = ionization_state(T)
n_He_pp = s(5) * n_H

END FUNCTION n_He_pp

REAL*8 FUNCTION n_e(T)
IMPLICIT NONE
REAL*8, INTENT(in):: T
REAL*8, DIMENSION(6):: s

s = ionization_state(T)
n_e = s(6) * n_H

END FUNCTION n_e



! Functions that define the cooling curve at any T
REAL*8 FUNCTION coll_ex_hydrogen(T)
IMPLICIT NONE
REAL*8:: T

coll_ex_hydrogen = 7.50e-19 * EXP(-118348.d0/T) * (1/(1 + SQRT(T/100000))) * n_e(T) * n_H0(T)

END FUNCTION coll_ex_hydrogen

REAL*8 FUNCTION coll_ex_helium(T)
IMPLICIT NONE
REAL*8:: T
coll_ex_helium = 5.54e-17 * T**(-0.397d0) * EXP(-473638.d0/T) * (1/(1 + SQRT(T/100000))) * n_e(T) * n_He_p(T)
END FUNCTION coll_ex_helium

REAL*8 FUNCTION coll_ion_hydrogen(T)
IMPLICIT NONE
REAL*8:: T
coll_ion_hydrogen = 1.27e-21 * SQRT(T) * EXP(-157809.1d0/T) * (1/(1 + SQRT(T/100000))) * n_e(T) * n_H0(T)
END FUNCTION coll_ion_hydrogen

REAL*8 FUNCTION coll_ion_helium(T)
IMPLICIT NONE
REAL*8:: T
coll_ion_helium = 9.38e-22 * SQRT(T) * EXP(-285335.4d0/T) * (1/(1 + SQRT(T/100000))) * n_e(T) * n_He0(T)
END FUNCTION coll_ion_helium

REAL*8 FUNCTION coll_ion_helium_p(T)
IMPLICIT NONE
REAL*8:: T
coll_ion_helium_p = 4.95e-22 * SQRT(T) * EXP(-631515.d0/T) * (1/(1 + SQRT(T/100000))) * n_e(T) * n_He_p(T)
END FUNCTION coll_ion_helium_p

REAL*8 FUNCTION recombination_hydrogen_p(T)
IMPLICIT NONE
REAL*8:: T
recombination_hydrogen_p = 8.70e-27 * SQRT(T) * ((T/1000)**(-0.2d0)) * 1/(1 + (T/1000000)**0.7d0) * n_e(T) * n_H_p(T)
END FUNCTION recombination_hydrogen_p

REAL*8 FUNCTION recombination_helium_p(T)
IMPLICIT NONE
REAL*8:: T
recombination_helium_p = 1.55e-26 * T**(0.3647d0) * n_e(T) * n_He_p(T)
END FUNCTION recombination_helium_p

REAL*8 FUNCTION recombination_helium_pp(T)
IMPLICIT NONE
REAL*8:: T

recombination_helium_pp = 3.48e-26 * SQRT(T) * ((T/1000)**-0.2) * (1/(1+(T/1000000)**0.7)) * n_e(T) * n_He_pp(T)

END FUNCTION recombination_helium_pp

REAL*8 FUNCTION dielectric_recomb(T)
iMPLICIT NONE
REAL*8:: T
dielectric_recomb = 1.24e-13 * T**-1.5d0 * EXP(-470000.d0/T) * (1+0.3d0*EXP(-94000.d0/T)) * n_e(T) * n_He_p(T)
END FUNCTION dielectric_recomb

REAL*8 FUNCTION free_free(T)
IMPLICIT NONE
REAL*8:: T
free_free = 1.42e-27 * SQRT(T) * (1.1d0 +0.34d0*EXP((-(5.5d0-LOG(T))**2) / 3.d0))&
 * n_e(T) * (n_H_p(T) + n_He_p(T) + 4*n_He_pp(T))
END FUNCTION free_free

REAL*8 FUNCTION cooling_func(T)
IMPLICIT NONE
REAL*8:: T
cooling_func = (coll_ex_hydrogen(T) + coll_ex_helium(T) + coll_ion_hydrogen(T) + coll_ion_helium(T)+&
coll_ion_helium_p(T)+recombination_hydrogen_p(T)+recombination_helium_p(T)+recombination_helium_pp(T)+&
dielectric_recomb(T)+free_free(T)) / n_H**2.d0
END FUNCTION cooling_func



! Functions that define the internal energy at any T
REAL*8 FUNCTION mu(T)

IMPLICIT NONE
REAL*8:: T

mu = (1+4*yhe) / (1+yhe+n_e(T)/n_H)

END FUNCTION

REAL*8 FUNCTION internal_energy(T)
IMPLICIT NONE
REAL*8:: T

internal_energy = (1/(gammaa - 1)) * ((k_boltz*T)/(mp*mu(T)))

END FUNCTION internal_energy

REAL*8 FUNCTION cooling_time(T)
IMPLICIT NONE
REAL*8:: T, rho

rho = n_H * mp / X
cooling_time = (internal_energy(T) * rho) / (cooling_func(T)*n_H**2.d0)

END FUNCTION cooling_time

REAL*8 FUNCTION u_t(T, u)
IMPLICIT NONE
REAL*8:: T, u

u_t = T - ((gammaa - 1) * (mp * mu(T) * u) / k_boltz)

END FUNCTION u_t

! This temperature functions returns T given a value of u + an approximate guess. Secant method
REAL*8 FUNCTION temperature(t0, u0)
  IMPLICIT NONE
  REAL*8:: u0, xn, t0
  REAL*8:: x0,x1,f0,f1
  REAL*8, PARAMETER:: tol=1.E-5
  INTEGER:: i

  i=0
  x0=t0              ! Starting guesses
  x1 = t0*1.1d0      
  2 f0=u_t(x0, u0)
  f1=u_t(x1, u0)
  xn = x0 - (f0*(x1-x0))/(f1-f0)    ! New point to iter from 
  i=i+1
  x0=x1
  x1=xn
  IF(i > 200)THEN
    PRINT*, "no res"
  ENDIF
  IF(ABS((x1-x0)/x0) .GT. tol) GO TO 2 ! Iterate this until the two guesses are around the zero with tol precision
    temperature = xn
  RETURN
  STOP
END FUNCTION temperature

END MODULE ionization

MODULE RHS
! Right hand side of the equation du/dt
USE ionization
  IMPLICIT NONE
  
  CONTAINS 
  
  SUBROUTINE dudt(guess, u, L)
  IMPLICIT NONE
  REAL*8, INTENT(in):: u
  REAL*8, INTENT(inout):: guess
  REAL*8, INTENT(out):: L
  
  L = - cooling_func(temperature(guess, u)) / cooling_func(1000000.d0) ! change if you change starting temperature
  
  END SUBROUTINE dudt
END MODULE RHS 

MODULE ODE
USE RHS
IMPLICIT NONE

CONTAINS 

SUBROUTINE RK2(guess, u_prev, u_next)
  IMPLICIT NONE
  REAL*8, INTENT(in):: u_prev
  REAL*8, INTENT(inout):: guess
  REAL*8, INTENT(out):: u_next
  REAL*8:: k1, k2, h
  
  h = 0.05d0 * cooling_time(temperature(guess, u_prev)) * n_H ! frazione del tcool attuale
  
  CALL dudt(guess, u_prev, k1)
  
  u_next = u_prev + h * k1  
  
  CALL dudt(guess, u_next, k2)
  
  u_next = u_prev + 0.5 * h * (k1+ k2)
  

END SUBROUTINE RK2

END MODULE ODE

PROGRAM cooling_exam
USE functions
USE ionization
USE ODE

IMPLICIT NONE
INTEGER, PARAMETER:: n = 200
INTEGER::i, iter
REAL*8, DIMENSION(0:4500):: energy_over_time
REAL*8, ALLOCATABLE:: timestep2(:), temperature_over_time(:), energy(:)
REAL*8:: ts, ionization_states(n,6)
REAL*8:: T(n), guess, time_0, u_0


! Uniformly spaced temperature points
ts = 4.d0 / 200.d0
DO i = 1, n
  T(i) = 10**(4 + ts * i)
ENDDO

!----------------------------------------------------------------------------------------------------------------------


! Step one: saving ionization state to file; display start and end on screen. Plot with python "graph_ionization.py"
DO i = 1,n
ionization_states(i,:) =  ionization_state(T(i))  ! Everything below this happens at the selected X
ENDDO

PRINT*, "Displaying ionization state at interval bounds 10^8 and 10^4"
WRITE(*,*)
PRINT*, "Starting ionization state at ", T(1),"K: "
WRITE(*,*) "Neutral Hydrogen: ", n_H0(T(1)) 
WRITE(*,*) "Ionized Hydrogen: ", n_H_p(T(1))
WRITE(*,*) "Neutral Helium: ", n_He0(T(1))
WRITE(*,*) "Once Ionized Helium: ", n_He_p(T(1))
WRITE(*,*) "Doubly Ionized Helium: ", n_He_pp(T(1))
WRITE(*,*) "Free Electrons in the gas: ", n_e(T(1))
WRITE(*,*)
PRINT*, "End ionization state at ", T(n),"K: "
WRITE(*,*) "Neutral Hydrogen: ", n_H0(T(n)) 
WRITE(*,*) "Ionized Hydrogen: ", n_H_p(T(n))
WRITE(*,*) "Neutral Helium: ", n_He0(T(n))
WRITE(*,*) "Once Ionized Helium: ", n_He_p(T(n))
WRITE(*,*) "Doubly Ionized Helium: ", n_He_pp(T(n))
WRITE(*,*) "Free Electrons in the gas: ", n_e(T(n))
WRITE(*,*)

CALL save_results("ionization_X"//x_name//"_nH"//h_name//".txt", 6, n, T, ionization_states)

PRINT*, "Complete ionization status saved in ionization_X"//x_name//"_nH"//h_name//".txt"
PRINT*, "Execute graph_ionization.py to plot all states..."
WRITE(*,*)
!----------------------------------------------------------------------------------------------------------------------



!----------------------------------------------------------------------------------------------------------------------
! Step two: saving cooling function to file. Plot with python "graph_cooling.py"
OPEN(33, file="cooling_X"//x_name//"_nH"//h_name//".txt", status="unknown")

DO i = 1, n
WRITE(33,'(12(1pe20.10))') LOG10(T(i)), LOG10(coll_ex_hydrogen(T(i))),LOG10(coll_ex_helium(T(i))), &
LOG10(coll_ion_hydrogen(T(i))), LOG10(coll_ion_helium(T(i))), LOG10(coll_ion_helium_p(T(i))),&
LOG10(recombination_hydrogen_p(T(i))),LOG10(recombination_helium_p(T(i))),&
LOG10(recombination_helium_pp(T(i))),LOG10(dielectric_recomb(T(i))),LOG10(free_free(T(i))),LOG10(cooling_func(T(i)))
ENDDO

PRINT*, "Complete cooling function saved in ", "cooling_X"//x_name//"_nH"//h_name//".txt"
PRINT*, "Execute graph_cooling.py to plot all cooling curves..."
WRITE(*,*) "---------------------------------------------------------------------"
!----------------------------------------------------------------------------------------------------------------------



! Step three: secant method for T(u) in ionization module


!----------------------------------------------------------------------------------------------------------------------
! Step four: integrating with variable timestep and saving temperature evolution for timestep

PRINT*, "Starting conditions:"
WRITE(*,*)
PRINT*, "Fraction of hydrogen in the gas X: ", X
PRINT*, "Hydrogen density nH: ", n_H
PRINT*, "Initial temperature T: ", 1000000
WRITE(*,*)
guess = 1000000.d0
time_0 = cooling_time(1000000.d0)
u_0 = internal_energy(1000000.d0)
energy_over_time(0) = u_0
iter = 0


DO WHILE (temperature(guess, energy_over_time(iter)) > 0.d0)  
  iter = iter + 1
  CALL RK2(guess, energy_over_time(iter-1), energy_over_time(iter)) 
END DO


ALLOCATE(energy(0:iter), temperature_over_time(0:iter), timestep2(0:iter))

timestep2(0) = 0.d0

DO i = 1, iter
  timestep2(i) = timestep2(i-1) + time_0/iter
END DO

DO i = 0, iter
  energy(i) = MAX(energy_over_time(i), internal_energy(10000.d0)) 
  temperature_over_time(i) = temperature(guess, energy(i))
END DO
!----------------------------------------------------------------------------------------------------------------------



!----------------------------------------------------------------------------------------------------------------------
! Displaying temperature evolution on screen and saving it on file. 

WRITE(*,*) "      Time (t0)", "               Temperature (K)"
WRITE(*,*) timestep2(1)/time_0, temperature_over_time(1)
WRITE(*,*) timestep2(2)/time_0, temperature_over_time(2)
WRITE(*,*) timestep2(3)/time_0, temperature_over_time(3)
WRITE(*,*) "."
WRITE(*,*) "."
WRITE(*,*) "."
WRITE(*,*) timestep2(iter-2)/time_0, temperature_over_time(iter-2)
WRITE(*,*) timestep2(iter-1)/time_0, temperature_over_time(iter-1)
WRITE(*,*) timestep2(iter)/time_0, temperature_over_time(iter)
WRITE(*,*)


! File with time and temperature
OPEN(77, file="evolution_temp_X"//x_name//"_nH"//h_name//".txt", status="unknown")

DO i = 0,iter
WRITE(77,'(3(1pe20.10))') timestep2(i)/time_0, temperature_over_time(i)
ENDDO
CLOSE(77)

! Get the first instance of temperature at ten thousand K and see what the time was
DO i = 1, iter
IF ( ABS(temperature_over_time(i) - 10000.d0) < 0.1d0) THEN
GO TO 6
ENDIF
ENDDO

6 PRINT*, "The simulated gas went from a 1,000,000 K to 10,000 K in ", &
timestep2(i) / 24.d0 / 365.d0 /60.d0/60.d0/1000000.d0, " Myr"


PRINT*, "Analitically calculated cooling time for 1,000,000 K to 10,700 K: ", (time_0-cooling_time(10700.d0))&
/ 24.d0 / 365.d0 / 60.d0 / 60.d0 / 1000000.d0, " Myr"

WRITE(*,*)
PRINT*, "Complete temperature evolution saved in evolution_temp_X"//x_name//"_nH"//h_name//".txt"

END PROGRAM

SUBROUTINE solve_linear_system(a,b,x,n)
!===========================================================
! Gauss elim w scaling and pivoting
! Input A(n,n)|b(n) output x(n). Input not preserved
!===========================================================
IMPLICIT NONE
INTEGER, INTENT(in):: n
REAL*8, INTENT(inout):: a(n,n), b(n)
REAL*8, INTENT(out):: x(n)
REAL*8:: s(n), c, pivot, cache
INTEGER:: i, j, k, l
! Forward elimination
DO k=1, n-1
  DO i=k,n 
    s(i) = 0.0
    DO j=k,n 
      s(i) = MAX(s(i),ABS(a(i,j)))  ! get max element of row
    ENDDO
   ENDDO
! Pivoting
    pivot = abs(a(k,k)/s(k))
    l=k
    DO j=k+1,n
      IF(ABS(a(j,k)/s(j)) > pivot) THEN
      pivot = ABS(a(j,k)/s(j))
      l=j
      ENDIF
    ENDDO
! Check if the system has a sigular matrix
IF (pivot == 0.0) THEN
WRITE(*,*) "Determinant is 0"
RETURN
ENDIF

IF (l /= k) THEN
  DO j=k,n
    cache = a(k,j)
    a(k,j) = a(l,j)
    a(l,j) = cache
  ENDDO
  cache = b(k)
  b(k) = b(l)
  b(l) = cache
ENDIF

DO i=k+1,n
  c=a(i,k)/a(k,k)
  a(i,k) = 0.0
  b(i)=b(i)- c*b(k)
  DO j=k+1,n
    a(i,j) = a(i,j)-c*a(k,j)
  ENDDO
ENDDO

ENDDO
! Backwards Sub
x(n) = b(n)/a(n,n)
DO i=n-1,1,-1
  c=0.0
  DO j=i+1,n
    c= c + a(i,j)*x(j)
  ENDDO
  x(i) = (b(i)- c)/a(i,i)
ENDDO
END SUBROUTINE solve_linear_system

SUBROUTINE save_results(filename, neq, npoints, x, y)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: npoints, neq
    REAL*8, DIMENSION(npoints), INTENT(IN) :: x
    REAL*8, DIMENSION(npoints, neq), INTENT(IN) :: y

    INTEGER :: i, j
    CHARACTER(LEN=25) :: forstring

    WRITE(forstring, '(a,i3,a)') '(', neq+1, '(1pe20.10))'

    OPEN(11, FILE=filename)

    DO i=1, npoints
      WRITE(11, forstring) x(i), (y(i,j), j=1, neq)
    END DO
    
    CLOSE(11)
END SUBROUTINE save_results


