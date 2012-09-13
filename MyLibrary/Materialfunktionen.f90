! *********************************************************************************
! *  XTwoodSim, A Solver for the Finite Element Software Elmer
! *  
! *  Copyright 11th September 2012 - , Uwe Jaschke
! *  
! *  This file is part of XTwoodSim.
! *  
! *      XTwoodSim is free software: you can redistribute it and/or modify
! *      it under the terms of the GNU General Public License as published by
! *      the Free Software Foundation, either version 3 of the License, or
! *      (at your option) any later version.
! *  
! *      XTwoodSim is distributed in the hope that it will be useful,
! *      but WITHOUT ANY WARRANTY; without even the implied warranty of
! *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *      GNU General Public License for more details.
! *  
! *      You should have received a copy of the GNU General Public License
! *      along with XTwoodSim.  If not, see <http://www.gnu.org/licenses/>.
! *  
! *      Diese Datei ist Teil von XTwoodSim.
! *  
! *      XTwoodSim ist Freie Software: Sie können es unter den Bedingungen
! *      der GNU General Public License, wie von der Free Software Foundation,
! *      Version 3 der Lizenz oder (nach Ihrer Option) jeder späteren
! *      veröffentlichten Version, weiterverbreiten und/oder modifizieren.
! *  
! *      XTwoodSim wird in der Hoffnung, dass es nützlich sein wird, aber
! *      OHNE JEDE GEWÄHELEISTUNG, bereitgestellt; sogar ohne die implizite
! *      Gewährleistung der MARKTFÄHIGKEIT oder EIGNUNG FÜR EINEN BESTIMMTEN ZWECK.
! *      Siehe die GNU General Public License für weitere Details.
! *  
! *      Sie sollten eine Kopie der GNU General Public License zusammen mit diesem
! *      Programm erhalten haben. Wenn nicht, siehe <http://www.gnu.org/licenses/>.
! *********************************************************************************
    ! DOCUMENTATION
    !> gibt die relative Feuchtigkeit passend zur Temperatur 
    !! und dem Feuchtegehalt nach Zurwitz et. al (Avramidis1989)
    !! zurück
    !!
    !> Zurwitz definiert folgenden Zusammenhang
    !! \f$ MC = \left[  -T \frac{ln(1-h)}{c_2 \left(1-\frac{T}{T_c}  \right)^{c_1}}\right]^{\frac{1}{c_3 T^{c_4}}}\f$
    !! \f$MC\f$ und die Temperatur sind für jeden Knoten bekannt und es folgt
    !! nach umstellen:
    !! \f$RH = 1 - \exp{\left(\frac{1}{-T} M^{c_3 T^{c_4}}c_2\left(1-\frac{T}{T_c}\right)^{c_1}\right)} \f$
FUNCTION RH_Zurwitz(model, n, mc) RESUlT(RH)
    USE DefUtils
    IMPLICIT NONE

    Type(model_t) :: model
    INTEGER :: n !< Knotennummer
    REAL(KIND=dp), INTENT(in) :: mc
    REAL(KIND=dp) :: RH !< relative humidity

    REAL(KIND=dp) :: c1 = -6.46,   &
        &            c2 = 0.13,    &
        &            c3 = 0.11d03, &
        &            c4 = -0.75,   &
        &            Tc = 647.1 

    TYPE(Variable_t), POINTER :: TempVar
    INTEGER, POINTER :: TempPerm(:)
    REAL(KIND=dp), POINTER :: Temperature(:)
    REAL(KIND=dp) :: localTemp

    TempVar => VariableGet(Model%Solver%Mesh%Variables, 'Temperature' )
    IF ( ASSOCIATED( TempVar) ) THEN
    TempPerm => TempVar % Perm
    Temperature => TempVar % Values
    !!!! stop if temperature field has not been found !!!!
    ELSE
    CALL Fatal('MyOwnSolver', 'No variable Temperature found')
    ENDIF
    
    localTemp = Temperature(TempPerm(n))

    !localTemp = localTemp - 273.15 !Konversion von Kelvin zu Grad

    ! T in Kelvin, mc in %
    RH = 100*(1-EXP(-1/localTemp * mc**(c3*localTemp**c4)*c2*(1-localTemp/Tc)**c1))
    !asprint *, 'At node ', n, 'we have a Temperature of ', localTemp, ' and a mc&
    !    & of', mc, 'which results in a RH of ', RH

    RETURN
END FUNCTION RH_Zurwitz

    ! DOCUMENTATION
    !> gibt \f$\frac{\delta\omega}{\delta H}\f$  passend zur Temperatur 
    !! und dem Feuchtegehalt nach Zurwitz et. al (Avramidis1989)
    !! zurück
    !!
    !> Zurwitz definiert folgenden Zusammenhang
    !! \f$ MC = \left[  -T \frac{ln(1-h)}{c_2 \left(1-\frac{T}{T_c}  \right)^{c_1}}\right]^{\frac{1}{c_3 T^{c_4}}}\f$
    !! \f$MC\f$ und die Temperatur sind für jeden Knoten bekannt und es folgt
    !! nach umstellen:
    !! \f$ \frac{\delta\omega}{\delta H}\f$= \frac{1}{c_3
    !! T^{c_4}}\left(\frac{T}{\left(1-h\right)c_2\left(1-\frac{T}{T_c}^{c_1}\right)^{\frac{1-c_3
    !! T^{c_4}}{c_3 T^{c_4}}\f$
FUNCTION dwdH_Zurwitz(model, n, mc) RESUlT(dwdH)
    USE DefUtils
    !USE myMaterialfunctions
    IMPLICIT NONE

    Type(model_t) :: model
    INTEGER :: n !< Knotennummer
    REAL(KIND=dp), INTENT(in) :: mc
    REAL(KIND=dp) :: dwdH !< relative humidity

    REAL(KIND=dp) :: RH_Zurwitz
    
    REAL(KIND=dp) :: c1 = -6.46, &
        &            c2 = 0.13, &
        &            c3 = 0.11d03, &
        &            c4 = -0.75, &
        &            Tc = 647.1 

    TYPE(Variable_t), POINTER :: TempVar
    INTEGER, POINTER :: TempPerm(:)
    REAL(KIND=dp), POINTER :: Temperature(:)
    REAL(KIND=dp) :: localTemp

    TempVar => VariableGet(Model%Solver%Mesh%Variables, 'Temperature' )
    IF ( ASSOCIATED( TempVar) ) THEN
    TempPerm => TempVar % Perm
    Temperature => TempVar % Values
    !!!! stop if temperature field has not been found !!!!
    ELSE
    CALL Fatal('MyOwnSolver', 'No variable Temperature found')
    ENDIF
    
    localTemp = Temperature(TempPerm(n))

    ! T in Kelvin, mc in %
    RH_Zurwitz = 100*(1-EXP(-1/localTemp * mc**(c3*localTemp**c4)*c2*(1-localTemp/Tc)**c1))

    dwdH =  1/(c3*localTemp**c4)* &
        & (localTemp/( &
        &  (1-RH_Zurwitz)*& 
        &   c2*(1-localTemp/Tc)**c1)&
        & )**((1-c3*localTemp**c4)/(c3*localTemp**c4))   

    RETURN
END FUNCTION dwdH_Zurwitz

    ! berechnet EMC für gegebene Luftfeuchte und Temp 
    ! nach Zurwitz (Avramidis1989)
    ! emc = [-T ln(1-rh)/(c_2(1-T/T_c)^{c_1})]^{1/(c_3 T^{c_4})
FUNCTION EMC_Zurwitz(model, n, time) RESULT(emc)
    USE DefUtils
    IMPLICIT NONE

    Type(model_t) :: model
    INTEGER :: n !< Knotennummer
    REAL(KIND=dp), INTENT(in) ::time
    REAL(KIND=dp) :: emc, rh

    REAL(KIND=dp),PARAMETER :: c1 = -6.46, &
        &            c2 = 0.13, &
        &            c3 = 0.11d03, &
        &            c4 = -0.75, &
        &            Tc = 647.1 

    TYPE(Variable_t), POINTER :: TempVar
    INTEGER, POINTER :: TempPerm(:)
    REAL(KIND=dp), POINTER :: Temperature(:)
    REAL(KIND=dp) :: localTemp

    REAL(KIND=dp) :: x1,x2,y1,y2

    IF (time.le.86400) THEN
        rh = 0.9
    ELSE IF (time.le.345600) THEN
        rh = 0.8
    ELSE IF (time.le.432000)THEN
        x1=345600
        x2=432000
        y1=0.8
        y2=0.6
        rh = (y1-y2)/(x1-x2)*time + (x1*y2-x2*y1)/(x1-x2)
    ELSE IF (time.le.1209600)THEN
        x1=432000
        x2=1209600
        y1=0.6
        y2=0.4
        rh = (y1-y2)/(x1-x2)*time + (x1*y2-x2*y1)/(x1-x2)
    ELSE IF (time.gt.1209600) THEN
        rh=0.55
    ENDIF

    TempVar => VariableGet(Model%Solver%Mesh%Variables, 'Temperature' )
    IF ( ASSOCIATED( TempVar) ) THEN
    TempPerm => TempVar % Perm
    Temperature => TempVar % Values
    !!!! stop if temperature field has not been found !!!!
    ELSE
    CALL Fatal('MyOwnSolver', 'No variable Temperature found')
    ENDIF

    localTemp = Temperature(TempPerm(n))

    emc =0.01*(-localTemp*LOG(1-rh)/(c2*(1-localTemp/Tc)**c1))**(1/(c3*localTemp**c4))

    RETURN
END FUNCTION

    ! DOCUMENTATION
    !> Specific Heat capacity according to Olek2003 in \f$J/(Kg K)\f$
    !!
    !> \f$c = \frac{0.0022}{1+0.01*M}*T^2 + \frac{3.32*0.01*M+2.95}{1+0.01M}T+
    !! \frac{4057*0.01M+526}{1+0.01M}\f$
    !! \f$M\f$ in %
    !! \f$T\f$ in Kelvin
FUNCTION cDeliiski(model,n,mc) RESULT(c)
    USE DefUtils
    
    TYPE(model_t), INTENT(in) :: model
    INTEGER, INTENT(in) :: n
    REAL(KIND=dp), INTENT(in) :: mc

    REAL(KIND=dp) :: c

    TYPE(Variable_t), POINTER :: TempVar
    INTEGER, POINTER :: TempPerm(:)
    REAL(KIND=dp), POINTER :: Temperature(:)
    REAL(KIND=dp) :: T
    
    TempVar => VariableGet(Model%Solver%Mesh%Variables, 'Temperature' )
    IF ( ASSOCIATED( TempVar) ) THEN
    TempPerm => TempVar % Perm
    Temperature => TempVar % Values
    !!!! stop if temperature field has not been found !!!!
    ELSE
    CALL Fatal('MyOwnSolver', 'No variable Temperature found')
    ENDIF
    
    T = Temperature(TempPerm(n))

    ! T in Kelvin, mc in %
    c=((0.0022)/(1+mc))*T**2 + ((3.32*mc+2.95)/(1+mc))*T + (4057*mc+526)/(1+mc) 
    RETURN
END FUNCTION cDeliiski

    !> Activierungsenergie gebundenen Wassers als Funktion
    !! des Feuchtegehalt
    !!
    !> In Eriksson2006 wir für die Aktivierungsenergie
    !! gebundenen Wassers \f$E_b = 500-290*MC \left[\frac{J}{mol}\right]
FUNCTION Eb(model,n,mc) RESUlT(e)
    USE DefUtils
    
    TYPE(model_t), INTENT(in) :: model
    INTEGER, INTENT(in) :: n
    REAL(KIND=dp), INTENT(in) :: mc

    REAL(KIND=dp) :: e

    e = 500-290*mc
    RETURN
END FUNCTION Eb

    ! DOCUMENTATION
    !> Heat conductivity tensor of european beech wood according to Olek2003 in \f$W/(m K)\f$
    !!
    !> \f$ k_T = 0.19933+0.18888*10^{-3}*(T-293.15)\f$
    !! k_R = 0.19958+0.33211*10^{-3}(T-293.15) \f$
    !! k_L = 0.29937+0.70147*10^{-3}(T-293.15) \f$
    !! \f$ diag(D_T) = (k_T, k_R, k_L)\f$ 
SUBROUTINE D_T_Beech_IHTP(model,n,T,D_T) 
    USE DefUtils

    TYPE(model_t) :: model
    INTEGER :: n
    REAL(KIND=dp) :: T

    REAL(KIND=dp), DIMENSION(:,:), POINTER :: D_T 

    REAL(KIND=dp) :: k_T, k_R, k_L

    T = T -273.15

    D_T = 0.0d0
    !T in Celsius
    k_T = 0.19933+0.18888*10**(-3)*(T) 
    k_R = 0.19958+0.33211*10**(-3)*(T)
    k_L = 0.29937+0.70147*10**(-3)*(T) 

    D_T(1,1) = k_R
    D_T(2,2) = k_T
    D_T(3,3) = k_L
END SUBROUTINE D_T_Beech_IHTP

    ! DOCUMENTATION
    !> Heat conductivity tensor of scots pine wood according to Olek2003 in \f$W/(m K)\f$
    !!
    !> \f$ k_T = 0.1989+0.8313*10^{-4}*(T-293.15)\f$
    !! k_R = 0.1990+0.8393*10^{-4}(T-293.15) \f$
    !! k_L = 0.2991+0.6184*10^{-4}(T-293.15) \f$
    !! \f$ diag(D_T) = (k_T, k_R, k_L)\f$ 
SUBROUTINE D_T_Pine_IHTP(model,n,T,D_T) 
    USE DefUtils

    TYPE(model_t) :: model
    INTEGER :: n
    REAL(KIND=dp) :: T

    REAL(KIND=dp), DIMENSION(:,:), POINTER :: D_T 

    REAL(KIND=dp) :: k_T, k_R, k_L

    T = T - 273.15
    D_T = 0.0d0
    !T in Celsius
    k_T = 0.1989+0.8313*10 **(-4)*(T)
    k_R = 0.1990+0.8393*10 **(-4)*(T) 
    k_L = 0.2991+0.6184*10 **(-4)*(T) 

    D_T(1,1) = k_R
    D_T(2,2) = k_T
    D_T(3,3) = k_L
END SUBROUTINE D_T_Pine_IHTP

    ! DOCUMENTATION
    !> Diffusion Coefficient matrix according to Eriksson2006
    !!
    !> \F$D_{\omega} = 2*10^{-9} m^2/s * e^{0.0641+0.04867*\omega}\f$
SUBROUTINE D_w_Eriksson(model,n,mc,Diffusivity) 
    USE DefUtils

    TYPE(model_t) :: model
    INTEGER :: n
    REAL(KIND=dp) :: mc

    REAL(KIND=dp), DIMENSION(:,:),POINTER :: Diffusivity

    REAL(KIND=dp) :: k, d_w0

    TYPE(Variable_t), POINTER :: TempVar
    INTEGER, POINTER :: TempPerm(:)
    REAL(KIND=dp), POINTER :: Temperature(:)
    REAL(KIND=dp) :: T
    
    TempVar => VariableGet(Model%Solver%Mesh%Variables, 'Temperature' )
    IF ( ASSOCIATED( TempVar) ) THEN
    TempPerm => TempVar % Perm
    Temperature => TempVar % Values
    !!!! stop if temperature field has not been found !!!!
    ELSE
    CALL Fatal('MyOwnSolver', 'No variable Temperature found')
    ENDIF
    
    T = Temperature(TempPerm(n))

    d_w0 = 2d-9

    Diffusivity=0.0d0
    
    ! hier mc korriegieren, da aus Eriksson2006 interpoliert wurde
    ! dort wird die alte skallierung verwendet
    ! mc in % * 100
    mc = mc*100
    k = EXP(-30.71 + 5.32*mc/450 + 0.0266 *T)
    !k= d_w0

    Diffusivity(1,1) = 3*k !radial
    Diffusivity(2,2) = 1*k !tangential
    Diffusivity(3,3) = 10*k !longitudinal

    !print *, Diffusivity(1,1), Diffusivity(2,2), Diffusivity(3,3), d_w0
END SUBROUTINE D_w_Eriksson
