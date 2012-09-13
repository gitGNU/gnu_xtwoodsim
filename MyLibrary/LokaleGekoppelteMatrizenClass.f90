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
    !> \defgroup LokaleGekoppelteMatrizenClass LokaleGekoppelteMatrizenClass
    !> \ingroup Elementeigenschaften
    !! \{
    !! \brief Verwaltung gekoppelte Elementmatrizen der FEM
    !!
    !> Diese Elementmatrizen werden in einem Type gesammelt
    !! und geschlossen behandelt. Diese Klasse erlaubt die Modifikation der
    !! Submatrizen. 
    !> \todo 
    !!  * überlege, welche Anordnung der Elementen
    !!    nach außen und innen sinnvoll ist, denn laut!!    Elmer Forum ist
    !!    für direkte Solver eine Blockanordnung der
    !!    Submatrizen ungeeignet, Iterative Solver
    !!    kommen damit wohl klar
MODULE LokaleGekoppelteMatrizenClass
      USE DefUtils !<Module profided by ELMER 
      USE LokaleMatrizenClass !< bereitgestellt in GeneralSolverPackage
      implicit none

          ! DOCUMENTATION
          !> \struct LokaleGekoppelteMatrizen_t
          !! Type Struktur zur Sammlung der Systemmatrizen
          !! und deren Submatrizen
          !! 
          !> vereint die MASS, STIFF and FORCE Arrays, sowie Pointer auf
          !! die Submatrizen
          !! Kaa, Kab, Kba, Kbb, Caa, Cab, Cba, Cbb, fa und fb
          !! in einem TYPE um Parameterlisten übersichtlicher zu
          !! gestalten und Functionen gleichzeitig auf alle Arrays 
      TYPE, PUBLIC:: LokaleGekoppelteMatrizen_t
          !PRIVATE
          TYPE(LokaleMatrizen_t), private :: LokaleMatrizen
          REAL( KIND=dp ), POINTER, DIMENSION(:,:):: Kaa, Kab,&
               & Kba, Kbb, Caa, Cab, Cba, Cbb
          REAL( KIND=dp ), POINTER, DIMENSION(:)  :: fa, fb
      END TYPE

          !> alles privat, Kommunikation zwischen Modulen 
          !! erfolgt über Interfaces
      PRIVATE:: addiere, setzeMatrizen, setzeAufWert, &
          & gibMASS, gibSTIFF, gibFORCE, initMatrizen, nullifiziereMatrizen
          !> Interfaces als Schnittstelle nach Ausen
          !! ermöglichen das sicher Überladen von
          !! Funktionen und Subroutinen
          !> \{
      INTERFACE operator(+)
          MODULE PROCEDURE addiere
      end interface
      INTERFACE neu
          MODULE procedure setzeMatrizen
          MODULE procedure initialisiereMatrizen
      end INTERFACE
      INTERFACE setze
          MODULE procedure setzeAufWert
      end INTERFACE
      INTERFACE nullifiziere
          MODULE PROCEDURE nullifiziereMatrizen
      END INTERFACE
      INTERFACE MASS
          MODULE PROCEDURE gibMASS
      END INTERFACE
      INTERFACE FORCE
          MODULE PROCEDURE gibFORCE
      END INTERFACE
      INTERFACE STIFF
          MODULE PROCEDURE gibSTIFF
      END INTERFACE
      Interface schuettleMatrizen
          MODULE PROCEDURE schuettleAlleMatrizen
      END INTERFACE

          !> \}
CONTAINS
          ! DOCUMENTATION
          !> gibSTIFF gibt die Steifigkeitsmatrix zurück
          !> \param this
      FUNCTION gibSTIFF(this)
          Type(LokaleGekoppelteMatrizen_t), INTENT(in):: this
          REAL(KIND=dp), POINTER :: gibSTIFF(:,:)
              gibSTIFF => STIFF(this % LokaleMatrizen)
          RETURN
      END FUNCTION gibSTIFF
      
          ! DOCUMENTATION
          !> gibMASS gibt die Massematrix zurück
      FUNCTION gibMASS(this)
          Type(LokaleGekoppelteMatrizen_t), INTENT(in):: this
          REAL(KIND=dp),  POINTER:: gibMASS(:,:)

              gibMASS => Mass(this % LokaleMatrizen)
          RETURN
      END FUNCTION gibMASS

          ! DOCUMENTATION
          !> gibFORCE gibt den Kraftvektor zurück
      FUNCTION gibFORCE(this)
          Type(LokaleGekoppelteMatrizen_t), INTENT(in):: this
          REAL(KIND=dp), POINTER:: gibFORCE(:)
              gibFORCE => Force(this % LokaleMatrizen)
          RETURN
      END FUNCTION gibFORCE
      
          ! DOCUMENTATION
          !> \brief setzeAufWert setzt die gesammten Matrizen auf Werte
          !! 
          !> FORCE, STIFF und MASS werden geschlossen auf Werte gesetzt.
          !! Das ist z.B. nützlich, um die Matrizen schnell auf 0 zu
          !! setzen
          !> \param this Sammlung der Matrizen 
          !> \param wertFORCE 0.0d0, falls nicht vorhanden
      SUBROUTINE setzeAufWert(this,wertFORCE, wertSTIFF, wertMASS)
          Type(LokaleGekoppelteMatrizen_t), INTENT(inout) :: this
          REAL(KIND=dp), DIMENSION(:), OPTIONAL:: wertFORCE
          REAL(KIND=dp), DIMENSION(:,:), OPTIONAL:: wertSTIFF!< alter Wert falls nicht vorhanden
          REAL(KIND=dp), DIMENSION(:,:), OPTIONAL:: wertMASS !< alter Wert falls nicht vorhanden
          
          IF(.not.PRESENT(wertForce)) wertForce = FORCE(this)
          IF(.not.PRESENT(wertSTIFF)) wertSTIFF = STIFF(this)
          IF(.not.Present(wertMass)) wertMass = Mass(this)
          CALL setze(this % LokaleMatrizen, wertForce, &
              & wertSTIFF, wertMass )
      END SUBROUTINE setzeAufWert

      SUBROUTINE nullifiziereMatrizen(this)
          TYPE(LokaleGekoppelteMatrizen_t), INTENT(inout) :: this

          CALL nullifiziere(this%LokaleMatrizen)
      END SUBROUTINE nullifiziereMatrizen
      
          ! DOCUMENTATION
          !> setzeMatrizen setzt die Systemmatrizen 

          !> und prüft die Matrizengrößen der Eingangsmatrizen
          !! , ruft initMatrizen auf und ordnet die Matrizen in den 
          !! LokaleMatrizen_t ein
          !> \todo füge weitere Routine für skalaren Input ein
      SUBROUTINE setzeMatrizen( this, aFORCE, aSTIFF, aMASS)
          Type(LokaleGekoppelteMatrizen_t), INTENT(inout):: this
          REAL(KIND=dp), INTENT(in) :: aSTIFF(:,:), aFORCE(:)
          REAL(KIND=dp), INTENT(in):: aMASS(:,:) !Make OPTIONAL

          CALL neu(this%LokaleMatrizen, aFORCE, aSTIFF, aMASS)
          CALL initMatrizen(this)
      END SUBROUTINE setzeMatrizen

      SUBROUTINE initialisiereMatrizen(this, dof)
          TYPE(LokaleGekoppelteMatrizen_t), INTENT(inout) :: this
          INTEGER, INTENT(in) :: dof
          
          REAL(KIND=dp), DIMENSION(:,:), POINTER :: matrix
          REAL(KIND=dp), DIMENSION(:), POINTER :: vektor
          INTEGER :: a,b

          CALL neu(this%LokaleMatrizen,dof)

          matrix => STIFF(this)
          a = SIZE(STIFF(this),1)
          b = SIZE(STIFF(this),2)
          WRITE(MESSAGE,*) ' LokaleGekoppelteMatrizen der Dimension ',&
              & a,'x',b,' initialisiert'
          CALL INFO('LokaleGekoppelteMatrizenClass', MESSAGE, level=5)

          matrix => STIFF(this)
          this%Kaa => matrix(1:a/2, 1:b/2)
          this%Kba => matrix(a/2+1:a, 1:b/2)
          this%Kab => matrix(1:a/2, b/2+1:b)
          this%Kbb => matrix(a/2+1:a, b/2+1:b)
          NULLIFY(matrix)

          matrix => MASS(this)
          this%Caa => matrix(1:a/2, 1:b/2)
          this%Cba => matrix(a/2+1:a, 1:b/2)
          this%Cab => matrix(1:a/2, b/2+1:b)
          this%Cbb => matrix(a/2+1:a, b/2+1:b)
          NULLIFY(matrix)

          vektor => FORCE(this)
          this%fa => vektor(1:a/2)
          this%fb => vektor(a/2+1:a)
          NULLIFY(vektor)
      END SUBROUTINE initialisiereMatrizen

          ! DOCUMENTATION
          !> initMatrizen fügt den Matrizen Pointer auf deren
          !! Submatrizen hinzu
          !!
          !> Die Matrizen sind Blockweise geteilt, sodass das 
          !! lineare Gleichungssystem die Form 
          !! \f[ 
          !!      \left[\begin{array}{c}
          !!            f_a \\
          !!            f_b
          !!            \end{array}\right] =
          !!      \left[\begin{array}{cc}
          !!             K_{aa}&K_{ab}\\
          !!             K_{ba}&K_{bb}
          !!            \end{array}\right]
          !!      \left[\begin{array}{c}
          !!            a \\
          !!            b
          !!            \end{array}\right] +
          !!      \left[\begin{array}{cc} 
          !!             c_{aa} & c_{ab} \\
          !!             c_{ba} & c_{bb}
          !!            \end{array}\right] 
          !!      \left[\begin{array}{c}
          !!            \dot{a} \\
          !!            \dot{b}
          !!            \end{array\right]}
          !! \f]
          !! hat. Pointer auf die Submatrizen sind im LokaleGekoppelteMatrizen_t
          !! gespeichert
      SUBROUTINE initMatrizen(this)
          TYPE( LokaleGekoppelteMatrizen_t ), INTENT(inout):: this
          REAL(KIND=dp), DIMENSION(:,:), POINTER :: matrix
          REAL(KIND=dp), DIMENSION(:), POINTER :: vektor
          INTEGER :: a,b

          print *, 'asd'
          matrix => STIFF(this)
          a = SIZE(STIFF(this),1)
          b = SIZE(STIFF(this),2)

          matrix => STIFF(this)
          this%Kaa => matrix(1:a/2, 1:b/2)
          this%Kab => matrix(a/2+1:a, 1:b/2)
          this%Kba => matrix(1:a/2, b/2+1:b)
          this%Kbb => matrix(a/2+1:a, b/2+1:b)

          matrix => MASS(this)
          this%Caa => matrix(1:a/2, 1:b/2)
          this%Cab => matrix(a/2+1:a, 1:b/2)
          this%Cba => matrix(1:a/2, b/2+1:b)
          this%Cbb => matrix(a/2+1:a, b/2+1:b)
          NULLIFY(matrix)

          vektor => FORCE(this)
          this%fa => vektor(1:a/2)
          this%fb => vektor(a/2+1:a)
          NULLIFY(vektor)
      END SUBROUTINE initMatrizen
       
           ! DOCUMENTATION
           !> \brief addiert zwei LokaleMatrizen
           !> indem die einzelenen Arrays (FORCE, STIFF, MASS) 
           !! einzeln und elementweise addiert werden
       TYPE(LokaleGekoppelteMatrizen_t) FUNCTION addiere(this,LM2) 
           TYPE(LokaleGekoppelteMatrizen_t), INTENT(in) :: this, LM2

           addiere%LokaleMatrizen = this%LokaleMatrizen + LM2%LokaleMatrizen
           Return 
       END FUNCTION addiere

       SUBROUTINE schuettleAlleMatrizen(this)
           TYPE(LokaleGekoppelteMatrizen_t), INTENT(inout) :: this

           CALL schuettleMatrizen(this%LokaleMatrizen)

       END SUBROUTINE schuettleAlleMatrizen

END MODULE LokaleGekoppelteMatrizenClass
!> \}
