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
!> \addtogroup Eriksson2006DGL_Class
      !> \ingroup Eriksson2006DGLImplentierung
      !! \{
      !> \brief Numerische Implementierung des gekoppelten Feldproblems
      !! nach Eriksson2006DGL 
      !!
      !> Berechnug der Elementmatrizen nach der Galerkin Methode für das 
      !! gekoppelte Temperatur-Feuchtigkeitsfeldproblem nach Eriksson2006. 
      !! Das Feldproblem wird nach einem Staggerd Iteration Schema gelöst,
      !! wobei die beiden Felder getrennt gelöst werden
MODULE Eriksson2006DGL_Class
      USE DefUtils
      USE LokaleGekoppelteMatrizenClass
      USE RandbedingungClass
      USE MaterialSammlung_Class
      USE geometrischeEigenschaften_Class
      USE MyStdModules
            
      IMPLICIT NONE
      
      PRIVATE
      Public::neu, berechne
      TYPE, Public :: Eriksson2006DGL_t
          PRIVATE
          TYPE(LokaleGekoppelteMatrizen_t), POINTER :: LokaleGekoppelteMatrizen
          TYPE(MaterialSammlung_t), POINTER :: materialWerte => NULL()
          TYPE(geometrischeEigenschaften_t), POINTER :: geometrie => NULL()
          TYPE(Randbedingung_t), POINTER :: temperaturRandbedingung=>NULL()
          TYPE(Randbedingung_t), POINTER :: feuchtigkeitRandbedingung=>NULL()
      END TYPE Eriksson2006DGL_t

      REAL(KIND=dp), DIMENSION(:,:), POINTER, private :: Lsg

      interface neu 
          MODULE PROCEDURE initialisiereDGLGebietselement
          MODULE PROCEDURE initialisiereDGLRandelement
      END interface
      interface berechne
          MODULE PROCEDURE berechneAlleMatrizen
      END interface
      interface gibLokalenWert
          MODULE PROCEDURE gibLokaleMaterialmatrix
          MODULE PROCEDURE gibLokalenMaterialwert
      END interface
  CONTAINS
      SUBROUTINE initialisiereDGLGebietselement(this,LokaleGekoppelteMatrizen,materialWerte, &
          & feuchtigkeitRandbedingung, temperaturRandbedingung, geometrie)
          TYPE(Eriksson2006DGL_t), INTENT(inout) :: this
          TYPE(LokaleGekoppelteMatrizen_t), INTENT(in), TARGET :: LokaleGekoppelteMatrizen
          TYPE(MaterialSammlung_t), INTENT(in), TARGET :: materialWerte
          TYPE(Randbedingung_t), INTENT(in), TARGET :: feuchtigkeitRandbedingung
          TYPE(Randbedingung_t), INTENT(in), TARGET :: temperaturRandbedingung
          TYPE(geometrischeEigenschaften_t), INTENT(in), TARGET :: geometrie

          this%LokaleGekoppelteMatrizen => LokaleGekoppelteMatrizen

          ! Nötig um auf VolumenRandbedingung zu prüfen
          this%feuchtigkeitRandbedingung => feuchtigkeitRandbedingung
          this%temperaturRandbedingung => temperaturRandbedingung

          this%geometrie => geometrie
          this%materialWerte => materialWerte

          Write(MESSAGE, *) 'density',this%materialWerte%dichte
          CALL INFO('Eriksson2006DGL_Class:initialisiereDGLGebietselement', MESSAGE, level=9)
          Write(MESSAGE, *) &
              & 'spezifischeWaerme',this%materialWerte%spezifischeWaerme
          CALL INFO('Eriksson2006DGL_Class:initialisiereDGLGebietselement', MESSAGE, level=9)
          Write(MESSAGE, *) 'Diffusivity',&
              & this%materialWerte%diffusionskoeffizientTensor(:,:,1)
          CALL INFO('Eriksson2006DGL_Class:initialisiereDGLGebietselement', MESSAGE, level=9)
          Write(MESSAGE, *) 'Conductivity',&
              & this%materialWerte%konduktivitaetTensor(:,:,1)
          CALL INFO('Eriksson2006DGL_Class:initialisiereDGLGebietselement', MESSAGE, level=9)
          Write(MESSAGE, *) 'Eb',this%materialWerte%E_b
          CALL INFO('Eriksson2006DGL_Class:initialisiereDGLGebietselement', MESSAGE, level=9)
          Write(MESSAGE, *) 'R',this%materialWerte%R
          CALL INFO('Eriksson2006DGL_Class:initialisiereDGLGebietselement', MESSAGE, level=9)
          Write(MESSAGE, *) 'RH',this%materialWerte%relativeHumidity
          CALL INFO('Eriksson2006DGL_Class:initialisiereDGLGebietselement', MESSAGE, level=9)
          Write(MESSAGE, *) 'dwdH',this%materialWerte%dwdH
          CALL INFO('Eriksson2006DGL_Class:initialisiereDGLGebietselement', MESSAGE, level=9)
          !> \todo Gleichungen für dieses konkrete element aufbereiten
          CALL INFO('Eriksson2006DGL_Class', 'initialisiereDGLGebietselement', &
              & level=5)
      END SUBROUTINE initialisiereDGLGebietselement

      SUBROUTINE initialisiereDGLRandelement(this, LokaleGekoppelteMatrizen,&
          & feuchtigkeitRandbedingung, temperaturRandbedingung, geometrie)
          TYPE(Eriksson2006DGL_t), INTENT(inout) :: this
          TYPE(LokaleGekoppelteMatrizen_t), INTENT(in), TARGET :: LokaleGekoppelteMatrizen
          TYPE(Randbedingung_t), INTENT(in), TARGET :: feuchtigkeitRandbedingung
          TYPE(Randbedingung_t), INTENT(in), TARGET :: temperaturRandbedingung
          TYPE(geometrischeEigenschaften_t), INTENT(in), TARGET :: geometrie

          this%LokaleGekoppelteMatrizen => LokaleGekoppelteMatrizen
          this%feuchtigkeitRandbedingung => feuchtigkeitRandbedingung
          this%temperaturRandbedingung => temperaturRandbedingung
          this%geometrie => geometrie

          CALL INFO('Eriksson2006DGL_Class', 'initialisiereDGLRandelement',&
              & level=5)
      END SUBROUTINE initialisiereDGLRandelement

      SUBROUTINE berechneAlleMatrizen(this, vorherigeLsg, istGebietselement)
          TYPE(Eriksson2006DGL_t), INTENT(inout) :: this
          REAL(KIND=dp), DIMENSION(:,:), TARGET, INTENT(in) :: vorherigeLsg
          LOGICAL, INTENT(in) :: istGebietselement
          LOGICAL :: stat
          INTEGER :: t !< Laufvariable

          Lsg => vorherigeLsg

              !> Berechne unterschiedliche Matrizen für Gebietselemente
              !! und Randgebietselemente
          IF(istGebietselement) THEN
                  !> numerische Integration
                  !! Im folgenden eine mathematische Beschreibung
                  !! am Beispiel einer Submatrix des gekoppelten Systems.
                  !! \f$ K_{aa} = \int_{Omega} B^T D_{\omega} B d\Omega\f$
                  !! Diese würde wie folgt numerisch integriert werden
                  !! \f$ K_{aa} \approx \sum_{t=1}^{N_{IP}} \left( \sqrt{ds^2}
                  !! \sqrt{det\left( J^T*J \right)} B^T D_{\omega} B \right)\mid_{IP}\f$
                  !! mit 
                  !!  - \f$J\f$ Jakobimatrix
                  !!  - \f$\sqrt{ds^2}\f$ Skalierung auf Modellkoordinatensystem
              DO t=1, this%geometrie%IP%n
                  WRITE(MESSAGE,*)'Integrationspunkt ', t, 'of ', this%geometrie%IP%n
                  CALL INFO('Eriksson2006DGL_Class',MESSAGE, level=8)
                  stat = Elementinfo(this%geometrie%element, &
                      &              this%geometrie%Nodes, &
                      &              this%geometrie%IP%U(t), &
                      &              this%geometrie%IP%V(t), &
                      &              this%geometrie%IP%W(T), &
                      &              this%geometrie%detJ, &
                      &              this%geometrie%N(1,:),&
                      &              this%geometrie%B,&
                      &              this%geometrie%dBdx, &
                      &              .false. )
                  !WRITE(MESSAGE, *) 'N', this%geometrie%N
                  !CALL INFO('Eriksson2006DGL_Class', MESSAGE, level=9)
                  ! aktualisiert Kaa Submatrix
                  this%LokaleGekoppelteMatrizen%Kaa = &
                      & this%LokaleGekoppelteMatrizen%Kaa + &
                      & K_omega(this,t)

                  ! aktualisiert Kab Submatrix
                  this%LokaleGekoppelteMatrizen%Kab = &
                      & this%LokaleGekoppelteMatrizen%Kab + &
                      & K_omega_T(this,t)
                  !Call printArray(this%LokaleGekoppelteMatrizen%Kab, 'Kab')

                  ! aktualisiert Kbb Submatrix
                  this%LokaleGekoppelteMatrizen%Kbb = &
                      & this%LokaleGekoppelteMatrizen%Kbb + &
                      & K_T(this,t)
                  ! aktualisiere Caa Submatrix
                  this%LokaleGekoppelteMatrizen%Caa = &
                      & this%LokaleGekoppelteMatrizen%Caa + &
                      & C_aa(this,t)
                  ! aktualisiere Cbb Submatrix
                  this%LokaleGekoppelteMatrizen%Cbb = &
                      & this%LokaleGekoppelteMatrizen%Cbb + &
                      & C_bb(this,t)

                  ! aktualisiere Cba Submatrix
                  this%LokaleGekoppelteMatrizen%Cba = &
                      & this%LokaleGekoppelteMatrizen%Cba + &
                      & C_ba(this,t)

                  ! aktualisiere fa Untervektor
                  IF(istVolumen(this%feuchtigkeitRandbedingung))THEN
                      CALL INFO('Eriksson2006DGL_Class', &
                          & 'Feuchtigkeitsvolumenlast', level=9)
                      this%LokaleGekoppelteMatrizen%fa = &
                          & this%LokaleGekoppelteMatrizen%fa + &
                          & f_omega_l(this,t)
                  ENDIF 
                  IF(istVolumen(this%temperaturRandbedingung)) THEN
                      CALL INFO('Eriksson2006DGL_Class', &
                          & 'Temperaturvolumenlast', level=9)
                      this%LokaleGekoppelteMatrizen%fb = &
                          & this%LokaleGekoppelteMatrizen%fb + &
                          & f_T_l(this,t)
                  ENDIF
                  !Call  printArray(this%LokaleGekoppelteMatrizen%Kaa, 'Kaa')
              END DO 
          Else IF(.not.istGebietselement) THEN 
              Do t=1,this%geometrie%IP%n
                  WRITE(MESSAGE,*)'Integrationspunkt ', t, 'of ', this%geometrie%IP%n
                  CALL INFO('Eriksson2006DGL_Class',MESSAGE, level=8)
                  stat = Elementinfo(&
                      &this%geometrie%element, &
                      &this%geometrie%Nodes, &
                      &this%geometrie%IP%U(t), &
                      &this%geometrie%IP%V(t), &
                      &this%geometrie%IP%W(T), &
                      &this%geometrie%detJ, &
                      &this%geometrie%N(1,:),&
                      &this%geometrie%B,&
                      &this%geometrie%dBdx, &
                      &.false. )
                  !aktualisiere Kaa
                  IF(istCauchy(this%feuchtigkeitRandbedingung)) THEN
                      CALL INFO('Eriksson2006DGL_Class', &
                          & 'Cauchy Feucht', level=9)
                      this%LokaleGekoppelteMatrizen%Kaa = &
                          & this%LokaleGekoppelteMatrizen%Kaa + &
                          & K_omega_c(this,t)
                      this%LokaleGekoppelteMatrizen%fa = &
                          & this%LokaleGekoppelteMatrizen%fa +&
                          & f_omega_c(this,t)
                  ENDIF
                  !aktualisiere Kbb
                  IF(istCauchy(this%temperaturRandbedingung)) THEN
                      CALL INFO('Eriksson2006DGL_Class', &
                          & 'Cauchy Temp', level=9)
                      this%LokaleGekoppelteMatrizen%Kbb = &
                          & this%LokaleGekoppelteMatrizen%Kbb + &
                          & K_T_c(this,t)
                      this%LokaleGekoppelteMatrizen%fb = &
                          & this%LokaleGekoppelteMatrizen%fb +&
                          & f_T_c(this,t)
                  ENDIF
                  IF(istNeumann(this%feuchtigkeitRandbedingung)) THEN
                      CALL INFO('Eriksson2006DGL_Class', &
                          & 'Neumann Feucht', level=9)
                      this%LokaleGekoppelteMatrizen%fa = &
                          & this%LokaleGekoppelteMatrizen%fa - &
                          & f_omega_N(this,t)
                  ENDIF
                  IF(istNeumann(this%temperaturRandbedingung)) THEN
                      CALL INFO('Eriksson2006DGL_Class', &
                          & 'Neumann Temp', level=9)
                      this%LokaleGekoppelteMatrizen%fb = &
                          & this%LokaleGekoppelteMatrizen%fb - &
                          & f_T_N(this,t)
                  ENDIF
                  !CALL printArray(this%LokaleGekoppelteMatrizen%fb, 'fb')
              END DO
          Else
              CALL FATAL('Eriksson2006DGL_Class', 'berechneAlleMatrizen kann &
                  & nicht bestimmen, ob die Matrizen für ein Gebiets- oder &
                  & Randgebietselemente berechnet werden sollen')
          ENDIF
      END SUBROUTINE berechneAlleMatrizen

          ! DOCUMENTATION
          !> gib den Lokalen Wert von den Knotenwerten, gewichtet nach
          !! den Basisfunktionen zurück
      FUNCTION gibLokaleMaterialmatrix(Knotenwerte, basisfunktionen)
          REAL(KIND=dp), DIMENSION(:,:,:), INTENT(in) :: Knotenwerte
          REAL(KIND=dp), DIMENSION(:), INTENT(in) :: basisfunktionen
          REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE:: gibLokaleMaterialmatrix

          INTEGER :: n,a,b,i
          REAL(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: dummy
          
          IF(SIZE(Knotenwerte,3).eq.SIZE(basisfunktionen)) THEN
              n = SIZE(Knotenwerte,3)
              a = SIZE(Knotenwerte,1)
              b = SIZE(Knotenwerte,2)
              ALLOCATE(gibLokaleMaterialmatrix(a,b), &
                  &    dummy(a,b,n))

              DO i=1, n
                  dummy(:,:,i) = Knotenwerte(:,:,i)*basisfunktionen(i)
              END DO
              gibLokaleMaterialmatrix = SUM(dummy,3)
          ELSE
              CALL FATAL('Eriksson2006DGL_Class','gibLokalenMaterialwert:&
                  & unterschiedlich lange Vektoren')
              ALLOCATE(gibLokaleMaterialmatrix(1,1))
              gibLokaleMaterialmatrix = 0.0d0
          ENDIF
          RETURN
      END FUNCTION gibLokaleMaterialmatrix
      
          ! DOCUMENTATION
          !> gib den Lokalen Wert von den Knotenwerten, gewichtet nach
          !! den Basisfunktionen zurück
      FUNCTION gibLokalenMaterialwert(Knotenwerte, basisfunktionen)
          REAL(KIND=dp), DIMENSION(:), INTENT(in) :: Knotenwerte
          REAL(KIND=dp), DIMENSION(:), INTENT(in) :: basisfunktionen
          REAL(KIND=dp):: gibLokalenMaterialwert

          INTEGER :: n
          IF(SIZE(Knotenwerte).eq.SIZE(basisfunktionen)) THEN
              n = SIZE(Knotenwerte)
              gibLokalenMaterialwert=SUM(Knotenwerte(1:n)*basisfunktionen(1:n))
          ELSE
              CALL FATAL('Eriksson2006DGL_Class','gibLokalenMaterialwert:&
                  & unterschiedlich lange Vektoren')
              gibLokalenMaterialwert = 0.0d0
          ENDIF
          RETURN
      END FUNCTION gibLokalenMaterialwert

      !> \f$K_{\omega} = \int_{\Omega} B^T D_{\omega} B d\Omega\f$
      FUNCTION K_omega(this,t)
          TYPE(Eriksson2006DGL_t), INTENT(IN) :: this
          INTEGER, INTENT(in) :: t

          REAL(KIND=dp),& !:: K_omega
              & DIMENSION(&
              &    SIZE(this%LokaleGekoppelteMatrizen%Kaa,1),&
              &    SIZE(this%LokaleGekoppelteMatrizen%Kaa,2)):: K_omega

          If(ASSOCIATED(this%materialWerte%diffusionskoeffizientTensor)) THEN
              K_omega = MATMUL(MATMUL(this%geometrie%B,&
                  & gibLokalenWert(this%materialWerte%diffusionskoeffizientTensor, &
                  &                this%geometrie%N(1,:)) ), &
                  & Transpose(this%geometrie%B))*&
                  & this%geometrie%detJ*this%geometrie%IP%s(t)
          ELSE IF(ASSOCIATED(this%materialWerte%diffusionskoeffizient))THEN
              K_omega = MATMUL(this%geometrie%B,&
                  & Transpose(this%geometrie%B))*&
                  & gibLokalenWert(this%materialWerte%diffusionskoeffizient, &
                  &                this%geometrie%N(1,:)) * &
                  & this%geometrie%detJ*this%geometrie%IP%s(t)
          Else
              CALL WARN('Eriksson2006DGL_Class','Für die Berechnung der &
                  & Matrix K_Omega fehlt der nötige Materialparameter &
                  & Diffusivity (diffusionskoeffizientTensor oder diffusionskoeffizient&
                  & ). Dieser wird auf 1 gesetzt')
              K_omega = this%geometrie%B*Transpose(this%geometrie%B)*&
                  & this%geometrie%detJ*this%geometrie%IP%s(t)
          ENDIF
          !CALL printArray(MATMUL(MATMUL(this%geometrie%B,&
          !        & gibLokalenWert(this%materialWerte%diffusionskoeffizientTensor, &
          !        &                this%geometrie%N(1,:)) ), &
          !        & Transpose(this%geometrie%B)), 'B.T*D*B')
          !print *, '(detjS, t, s(t)) ',this%geometrie%detJ*this%geometrie%IP%s(t), t, this%geometrie%IP%s(t)
          RETURN
      END FUNCTION K_omega

      !> \f$K_{T} = \int_{\Omega} B^T D_{T} B d\Omega\f$
      FUNCTION K_T(this,t)
          TYPE(Eriksson2006DGL_t), INTENT(IN) :: this
          INTEGER, INTENT(in) :: t

          REAL(KIND=dp),& ! :: K_T
              & DIMENSION(&
              &    SIZE(this%LokaleGekoppelteMatrizen%Kbb,1),&
              &    SIZE(this%LokaleGekoppelteMatrizen%Kbb,2)):: K_T
          
          If(ASSOCIATED(this%materialWerte%KonduktivitaetTensor)) THEN
              K_T = MATMUL(MATMUL(this%geometrie%B,&
                  & gibLokalenWert(this%materialWerte%KonduktivitaetTensor, &
                  &                this%geometrie%N(1,:)) ), &
                  & Transpose(this%geometrie%B))*&
                  & this%geometrie%detJ*this%geometrie%IP%s(t)
          ELSE IF(ASSOCIATED(this%materialWerte%Konduktivitaet))THEN
              K_T = MATMUL(this%geometrie%B,&
                  & Transpose(this%geometrie%B))*&
                  & gibLokalenWert(this%materialWerte%Konduktivitaet, &
                  &                this%geometrie%N(1,:)) * &
                  & this%geometrie%detJ*this%geometrie%IP%s(t)
          Else
              CALL WARN('Eriksson2006DGL_Class','Für die Berechnung der &
                  & Matrix K_T fehlt der nötige Materialparameter &
                  & Diffusivity (KonduktivitaetTensor oder Konduktivitaet&
                  & ). Dieser wird auf 1 gesetzt')
              K_T = this%geometrie%B*Transpose(this%geometrie%B)*&
                  & this%geometrie%detJ*this%geometrie%IP%s(t)
          ENDIF
          RETURN
      END FUNCTION K_T
      
      FUNCTION K_omega_T(this,t)
          TYPE(Eriksson2006DGL_t), INTENT(IN) :: this
          INTEGER, INTENT(in) :: t
          REAL(KIND=dp), & ! :: K_omega_T
              & DIMENSION(&
              &    SIZE(this%LokaleGekoppelteMatrizen%Kab,1),&
              &    SIZE(this%LokaleGekoppelteMatrizen%kab,2)):: K_omega_T

          If(ASSOCIATED(this%materialWerte%diffusionskoeffizientTensor)) THEN
              K_omega_T = MATMUL(MATMUL(this%geometrie%B,&
                  & gibLokalenWert(this%materialWerte%diffusionskoeffizientTensor, &
                  &                this%geometrie%N(1,:)) ), &
                  & Transpose(this%geometrie%B))*K_omega_T_Sorreteinfluss(this,t)*&
                  & this%geometrie%detJ*this%geometrie%IP%s(t)
          ELSE IF(ASSOCIATED(this%materialWerte%diffusionskoeffizient))THEN
              K_omega_T = MATMUL(this%geometrie%B,&
                  & Transpose(this%geometrie%B))*K_omega_T_Sorreteinfluss(this,t)*&
                  & gibLokalenWert(this%materialWerte%diffusionskoeffizient, &
                  &                this%geometrie%N(1,:)) * &
                  & this%geometrie%detJ*this%geometrie%IP%s(t)
          Else
              CALL WARN('Eriksson2006DGL_Class','Für die Berechnung der &
                  & Matrix K_Omega fehlt der nötige Materialparameter &
                  & Diffusivity (diffusionskoeffizientTensor oder diffusionskoeffizient&
                  & ). Dieser wird auf 1 gesetzt')
              K_omega_T = this%geometrie%B*Transpose(this%geometrie%B)*&
                  & K_omega_T_Sorreteinfluss(this,t)*&
                  & this%geometrie%detJ*this%geometrie%IP%s(t)
          ENDIF
          !CALL printArray(K_omega_T, 'K_omega_T')
          RETURN
      END FUNCTION K_omega_T

          !> Mathematische Beschreibung des Einflusses des Sorret-Effektes
          !! auf die Submatrix \f$K_{ab}\f$
          !! 
          !> \f$k_{\omega T} = \frac{H}{R T}\frac{\delta\omega}{\delta
          !! H}\frac{E_b}{T}\f$
      FUNCTION K_omega_T_Sorreteinfluss(this,t)
          TYPE(Eriksson2006DGL_t), INTENT(IN) :: this
          INTEGER, INTENT(in) :: t
          REAL(KIND=dp) :: K_omega_T_Sorreteinfluss
          K_omega_T_Sorreteinfluss = &
              & gibLokalenWert(&
              &   this%materialWerte%relativeHumidity / &
              &   (this%materialWerte%R*(Lsg(2,:)-272.15))* &
              &   this%materialWerte%dwdH * &
              &   this%materialWerte%E_b/(Lsg(2,:)-272.15), &
              &   this%geometrie%N(1,:))

          RETURN 
      END FUNCTION K_omega_T_Sorreteinfluss

      FUNCTION C_aa(this,t)
          TYPE(Eriksson2006DGL_t), INTENT(IN) :: this
          INTEGER, INTENT(in) :: t
          REAL(KIND=dp), & ! :: C_aa
              & DIMENSION(&
              &    SIZE(this%LokaleGekoppelteMatrizen%Caa,1),&
              &    SIZE(this%LokaleGekoppelteMatrizen%Caa,2)):: C_aa

          C_aa = MATMUL(Transpose(this%geometrie%N),this%geometrie%N)*&
              &  this%geometrie%detJ*this%geometrie%IP%s(t)

          RETURN
      END FUNCTION C_aa
      
      FUNCTION C_bb(this,t)
          TYPE(Eriksson2006DGL_t), INTENT(IN) :: this
          INTEGER, INTENT(in) :: t
          REAL(KIND=dp), & ! :: C_bb
              & DIMENSION(&
              &     SIZE(this%LokaleGekoppelteMatrizen%Cbb,1),&
              &     SIZE(this%LokaleGekoppelteMatrizen%cbb,2)):: C_bb

          C_bb = MATMUL(&
              &     Transpose(&
              &         this%geometrie%N),&
              &     this%geometrie%N &
              &  )*&
              &  gibLokalenWert(&
              &      this%materialWerte%spezifischeWaerme, &
              &      this%geometrie%N(1,:)&
              &  )*&
              &  gibLokalenWert(&
              &      this%materialWerte%dichte, &
              &      this%geometrie%N(1,:)&
              &  )* &
              &  this%geometrie%detJ*this%geometrie%IP%s(t)

          RETURN
      END FUNCTION C_bb

      FUNCTION C_ba(this,t)
          TYPE(Eriksson2006DGL_t), INTENT(IN) :: this
          INTEGER, INTENT(in) :: t
          REAL(KIND=dp),& ! :: C_ba
              & DIMENSION(&
              &     SIZE(this%LokaleGekoppelteMatrizen%Cba,1),&
              &     SIZE(this%LokaleGekoppelteMatrizen%Cba,2)):: C_ba

          C_ba = MATMUL(Transpose(this%geometrie%N),this%geometrie%N)*&
              & gibLokalenWert(this%materialWerte%E_b, &
              &     this%geometrie%N(1,:))/(-0.018)*&
              &  this%geometrie%detJ*this%geometrie%IP%s(t)

          RETURN
      END FUNCTION C_ba

          !> Integrant für Kraftuntervektor der Volumenlast
          !! des Feuchtigkeitsfeldproblem
      FUNCTION f_omega_l(this,t)
          TYPE(Eriksson2006DGL_t), INTENT(in) :: this
          INTEGER, INTENT(in) :: t

          REAL(KIND=dp), DIMENSION(& ! :: f_omega_l
              &SIZE(this%geometrie%N)):: f_omega_l

          f_omega_l = &
              &       this%geometrie%N(1,:)*&
              &       gibLokalenWert(&
              &         gibRandbedingungVolumen(this%feuchtigkeitRandbedingung), &
              &         this%geometrie%N(1,:)&
              &       )* &
              &       this%geometrie%detJ * this%geometrie%IP%s(t)
          RETURN
      END FUNCTION f_omega_l

          !> Integrant für Kraftuntervektor der Volumenlast
          !! des Temperaturfeldproblems
      FUNCTION f_T_l(this,t)
          TYPE(Eriksson2006DGL_t), INTENT(in) :: this
          INTEGER, INTENT(in) :: t

          REAL(KIND=dp),& ! :: f_T_l
              & DIMENSION(SIZE(this%geometrie%N)):: f_T_l

          f_T_l = this%geometrie%N(1,:) * gibLokalenWert(&
              & gibRandbedingungVolumen(this%temperaturRandbedingung), &
              & this%geometrie%N(1,:)) * &
              & this%geometrie%detJ * this%geometrie%IP%s(t)
          RETURN
      END FUNCTION f_T_l

          !> Integrant der Kraftsubmatrix für cauchy-
          !! Randbedingung des Feuchtigkeitsfeldproblems
      FUNCTION K_omega_c(this,t)
          TYPE(Eriksson2006DGL_t), INTENT(in) :: this
          INTEGER, INTENT(in) :: t

          REAL(KIND=dp), & ! :: K_omega_c
              & DIMENSION(&
              &   SIZE(this%LokaleGekoppelteMatrizen%Kaa,1),&
              &   SIZE(this%LokaleGekoppelteMatrizen%Kaa,2)) :: K_omega_c

          REAL(KIND=dp), DIMENSION(:), POINTER :: WertCauchy

          WertCauchy => gibRandbedingungCauchy(this%feuchtigkeitRandbedingung,2)

          !K_omega_c = MATMUL(Transpose(dummy),dummy)*&
          K_omega_c = MATMUL(Transpose(this%geometrie%N),this%geometrie%N)*&
              & gibLokalenWert(WertCauchy, & ! alpha
              &     this%geometrie%N(1,:))*&
              &  this%geometrie%detJ*this%geometrie%IP%s(t)
          RETURN          
      END FUNCTION K_omega_c
      
          !> Integrant der Kraftsubmatrix für cauchy-
          !! Randbedingung des Feuchtigkeitsfeldproblems
      FUNCTION K_T_c(this,t)
          TYPE(Eriksson2006DGL_t), INTENT(in) :: this
          INTEGER, INTENT(in) :: t

          REAL(KIND=dp), & ! :: K_T_c
              & DIMENSION(&
              &   SIZE(this%LokaleGekoppelteMatrizen%Kbb,1),&
              &   SIZE(this%LokaleGekoppelteMatrizen%Kbb,2)) :: K_T_c

          REAL(KIND=dp), DIMENSION(this%geometrie%anzahlKnoten):: WertCauchy

          WertCauchy = gibRandbedingungCauchy(this%temperaturRandbedingung,2)

          IF(this%geometrie%anzahlKnoten.eq.1)THEN
              K_T_c(1,1) = &
                  & WertCauchy(1)* &
                  & this%geometrie%detJ*this%geometrie%IP%s(t)
          ELSE
              K_T_c = MATMUL(Transpose(this%geometrie%N),this%geometrie%N)*&
                  & gibLokalenWert(WertCauchy, & ! beta
                  &     this%geometrie%N(1,:))*&
                  &  this%geometrie%detJ*this%geometrie%IP%s(t)
          ENDIF
          RETURN          
      END FUNCTION K_T_c

      FUNCTION f_omega_c(this,t) RESULT(res)
          TYPE(Eriksson2006DGL_t), INTENT(in) :: this
          INTEGER, INTENT(in) :: t

          REAL(KIND=dp), DIMENSION(this%geometrie%anzahlKnoten) :: res

          REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: WertCauchy
          INTEGER, POINTER :: n

          n => this%geometrie%anzahlKnoten
          ALLOCATE(WertCauchy(2,n))

          WertCauchy(1,1:n) = gibRandbedingungCauchy(&
              &                 this%feuchtigkeitRandbedingung,1)
          WertCauchy(2,1:n) = gibRandbedingungCauchy(&
              &                 this%feuchtigkeitRandbedingung,2)

          res = this%geometrie%N(1,1:n)*&
              &   gibLokalenWert(&
              &     WertCauchy(1,1:n),&
              &     this%geometrie%N(1,1:n))* &
              &   gibLokalenWert(&
              &     WertCauchy(2,1:n),&
              &     this%geometrie%N(1,1:n))* &
              &   this%geometrie%detJ*this%geometrie%IP%s(t)

          DEALLOCATE(WertCauchy)
          RETURN
      END FUNCTION f_omega_c

      FUNCTION f_T_c(this,t) RESULT(res)
          TYPE(Eriksson2006DGL_t), INTENT(in) :: this
          INTEGER, INTENT(in) :: t

          REAL(KIND=dp), DIMENSION(this%geometrie%anzahlKnoten) :: res
          REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE:: WertCauchy

          INTEGER, POINTER :: n

          n => this%geometrie%anzahlKnoten
          ALLOCATE(WertCauchy(2,n))
          WertCauchy(1,1:n) = gibRandbedingungCauchy(&
              &              this%temperaturRandbedingung,1)
          WertCauchy(2,1:n) = gibRandbedingungCauchy(&
              &              this%temperaturRandbedingung,2)

          res(1:n) = this%geometrie%N(1,1:n)*&
          & gibLokalenWert(WertCauchy(1,1:n),this%geometrie%N(1,1:n))*&
          & gibLokalenWert(WertCauchy(2,1:n),this%geometrie%N(1,1:n))*&
          & this%geometrie%detJ*this%geometrie%IP%s(t)
          DEALLOCATE(WertCauchy)

          RETURN
      END FUNCTION f_T_c

      FUNCTION f_omega_N(this,t) RESULT(res)
          TYPE(Eriksson2006DGL_t), INTENT(in) :: this
          INTEGER, INTENT(in) :: t

          REAL(KIND=dp), DIMENSION(this%geometrie%anzahlKnoten) :: res

          REAL(KIND=dp), DIMENSION(:), POINTER :: WertNeumann
          INTEGER, POINTER :: n

          n => this%geometrie%anzahlKnoten

          WertNeumann => gibRandbedingungNeumann(this%feuchtigkeitRandbedingung)

          WRITE(MESSAGE,*) 'WertNeumann Moist=', WertNeumann
          CALL INFO('Eriksson2006DGL_Class', MESSAGE, level=9)

          res(1:n) = this%geometrie%N(1,1:n) * &
               &      gibLokalenWert(WertNeumann, this%geometrie%N(1,1:n))*&
               &      this%geometrie%detJ*this%geometrie%IP%s(t)
          !res(1:n) = this%geometrie%N(1,1:n) * &
          !     &     WertNeumann(1:n)*&
          !     &     this%geometrie%detJ*this%geometrie%IP%s(t)
          RETURN
      END FUNCTION f_omega_N
      
      FUNCTION f_T_N(this,t) RESULT(res)
          TYPE(Eriksson2006DGL_t), INTENT(in) :: this
          INTEGER, INTENT(in) :: t

          REAL(KIND=dp), DIMENSION(this%geometrie%anzahlKnoten) :: res

          REAL(KIND=dp), DIMENSION(:), POINTER :: WertNeumann
          INTEGER, POINTER :: n

          n => this%geometrie%anzahlKnoten

          WertNeumann => gibRandbedingungNeumann(this%temperaturRandbedingung)

          WRITE(MESSAGE,*) 'WertNeumann Heat=', WertNeumann
          CALL INFO('Eriksson2006DGL_Class', MESSAGE, level=9)

          res = this%geometrie%N(1,1:n) * &
               &      gibLokalenWert(WertNeumann, this%geometrie%N(1,1:n))*&
               &      this%geometrie%detJ*this%geometrie%IP%s(t)
          !res = this%geometrie%N(1,1:n) * &
          !     &      WertNeumann(1:n)*&
          !     &      this%geometrie%detJ*this%geometrie%IP%s(t)
          RETURN
      END FUNCTION f_T_N
END MODULE Eriksson2006DGL_Class

