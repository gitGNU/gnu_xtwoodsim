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
    !> \defgroup LokaleMatrizenClass LokaleMatrizenClass
    !> \ingroup Elementeigenschaften
    !! \{
    !! \brief Verwaltung Elementmatrizen der FEM
    !!
    !> Diese Elementmatrizen werden in einem Type gesammelt
    !! und geschlossen behandelt. Dazu stellt diese KLasse diverse 
    !! Funktionen bereit um Matrizen zu manipulieren.
MODULE LokaleMatrizenClass
      USE DefUtils !<Module profided by ELMER 
      implicit none

          !> \struct LokaleMatrizen_t
          !! Type Struktur zur Sammlung der Systemmatrizen
          !!
          !> vereint die MASS, STIFF and FORCE Arrays
          !! in einem TYPE um Parameterlisten übersichtlicher zu
          !! gestalten und Functionen gleichzeitig auf alle Arrays 
          !! anzuwenden. Es werden allokierbare POINTER als Behelflösung
          !! verwendet, da allokierbare Arrays innerhalb von TYPE Definitionen
          !! nicht erlaubt sind.
          !! Die Matrizen mit Shake* geben den für elmer umsortierten
          !! Matrizenelementanordnungen wieder. 
      TYPE, PUBLIC :: LokaleMatrizen_t
          private
          REAL( KIND=dp ), private, POINTER, DIMENSION(:,:):: MASS, STIFF
          REAL( KIND=dp ), private, POINTER, DIMENSION(:,:):: ShakedMASS,&
               & ShakedSTIFF
          REAL( KIND=dp ), private, POINTER, DIMENSION(:)  :: FORCE
          REAL( KIND=dp ), private, POINTER, DIMENSION(:)  :: ShakedFORCE
          LOGICAL, private:: AllocationsDone=.FALSE.
      END TYPE

      private :: initMatrizen, setzeMatrizen, addiere, &
              & gibSTIFF, schuettleAlleMatrizen, &
              & gibMASS, gibFORCE, setzeAufWert, nullifiziereMatrizen
    
      interface operator(+)
          MODULE PROCEDURE addiere
      end interface
      interface neu
          MODULE procedure setzeMatrizen
          MODULE procedure initMatrizen
          MODULE procedure initMatrizenKurz
      end INTERFACE
      interface setze
          MODULE procedure setzeAufWert
      end INTERFACE
      INTERFACE nullifiziere
          MODULE PROCEDURE nullifiziereMatrizen
      END INTERFACE
      INTERFACE MASS
          MODULE procedure gibMASS
      END INTERFACE
      INTERFACE FORCE
          MODULE procedure gibFORCE
      END INTERFACE
      INTERFACE STIFF
          MODULE procedure gibSTIFF
      END INTERFACE
      Interface schuettleMatrizen
          MODULE PROCEDURE schuettleAlleMatrizen
      END INTERFACE
CONTAINS
          ! DOCUMENTATION
          !> gibSTIFF gibt die Steifigkeitsmatrix zurück
      FUNCTION gibSTIFF(this)
          !> LokaleMatrizen_t:: this
          !! Returns REAL(KIND=dp), DIMENSION(SHAPE(STIFF)
          TYPE(LokaleMatrizen_t), INTENT(in):: this
          REAL(KIND=dp), POINTER :: gibSTIFF(:,:)
          IF( this%AllocationsDone ) THEN
              gibSTIFF => this % STIFF
          ELSE
              print *, 'Allokation der Lokalen Matrizen nicht &&
                      && erfolgreich durchgeführt'
          ENDIF
          RETURN
      END FUNCTION gibSTIFF
      
          ! DOCUMENTATION
          !> gibMASS gibt die Massematrix zurück
      FUNCTION gibMASS(this)
          TYPE(LokaleMatrizen_t), INTENT(in):: this
          REAL(KIND=dp),  POINTER:: gibMASS(:,:)

          IF( this%AllocationsDone ) THEN
              gibMASS => this % MASS
          ELSE
              print *, 'Allokation der Lokalen Matrizen nicht &&
                      && erfolgreich durchgeführt'
          ENDIF
          RETURN
      END FUNCTION gibMASS

          ! DOCUMENTATION
          !> gibFORCE gibt den Kraftvektor zurück
      FUNCTION gibFORCE(this)
          TYPE(LokaleMatrizen_t), INTENT(in):: this
          REAL(KIND=dp), POINTER:: gibFORCE(:)
          IF( this%AllocationsDone ) THEN
              gibFORCE => this % FORCE
          ELSE
              print *, 'Allokation der Lokalen Matrizen nicht &&
                      && erfolgreich durchgeführt'
          ENDIF
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
          TYPE(LokaleMatrizen_t), INTENT(inout) :: this
          REAL(KIND=dp), DIMENSION(:), OPTIONAL:: wertFORCE!< alter Wert falls nicht vorhanden
          REAL(KIND=dp), DIMENSION(:,:), OPTIONAL:: wertSTIFF!< alter Wert falls nicht vorhanden
          REAL(KIND=dp), DIMENSION(:,:), OPTIONAL:: wertMASS !< alter Wert falls nicht vorhanden

          IF (PRESENT(wertForce).and.PRESENT(wertSTIFF).and.PRESENT(wertMASS)) THEN
              print *, 'setzeAufWert TRUE START'
              this%Force = wertForce
              this%STIFF = wertSTIFF
              this%MASS  = wertMass
              print *, 'setzeAufWert TRUE END'
          ELSE
              print *, 'setzeAufWert FALSE START'
              IF (.not.PRESENT(wertForce)) wertFORCE=FORCE(this)
              IF (.not.PRESENT(wertSTIFF)) wertSTIFF=STIFF(this)
              IF (.not.PRESENT(wertMass)) wertMASS=MASS(this)
              this%Force = wertForce
              this%STIFF = wertSTIFF
              this%MASS  = wertMass
              print *, 'setzeAufWert FALSE END'
          ENDIF
      END SUBROUTINE setzeAufWert
      
          ! DOCUMENTATION
          !> setzeMatrizen setzt die Systemmatrizen 

          !> und prüft die Matrizengrößen der Eingangsmatrizen
          !! , ruft initMatrizen auf und ordnet die Matrizen in den 
          !! LokaleMatrizen_t ein
          !> \todo mache diese Funktion unabhängig von Eingangsmatrizen
          !! übergebe stattdessen die Anzahl der Freiheitsgrade
          !! pro Element
      SUBROUTINE setzeMatrizen( this, aFORCE, aSTIFF, aMASS)
          TYPE(LokaleMatrizen_t), INTENT(inout):: this
          REAL(KIND=dp), INTENT(in) :: aSTIFF(:,:), aFORCE(:)
          REAL(KIND=dp), INTENT(in):: aMASS(:,:) !Make OPTIONAL

          IF (.not.this%AllocationsDone) THEN
              IF (SIZE(aSTIFF,1)==SIZE(aMASS,1) .AND. &
                  SIZE(aSTIFF,2)==SIZE(aMASS,2) .AND. &
                  SIZE(aSTIFF,1)==SIZE(aSTIFF,2) .AND. &
                  SIZE(aSTIFF,1)==SIZE(aFORCE) ) THEN
                  CALL initMatrizen(this,SHAPE(aMASS))
              ELSE
                  WRITE( MESSAGE, *) &
                  & 'Matrizengrößen fehlerhaft! &
                  & Das Gleichungssystem kann nicht gelöst &
                  & werden.', &
                  & 'SIZE(FORCE)=', SHAPE(aFORCE), &
                  & 'SIZE(MASS)=',  SHAPE(aMASS), &
                  & 'SIZE(STIFF)=', SHAPE(aSTIFF)
                  CALL Fatal('setzeMatrizen', MESSAGE)
              ENDIF
          ENDIF
          this%STIFF = aSTIFF
          this%MASS = aMASS
          this%FORCE = aFORCE
      END SUBROUTINE setzeMatrizen

      SUBROUTINE nullifiziereMatrizen(this)
          TYPE(LokaleMatrizen_t), INTENT(inout) :: this

          this%STIFF = 0.0d0
          this%MASS  = 0.0d0
          this%FORCE = 0.0d0
      END SUBROUTINE nullifiziereMatrizen

      SUBROUTINE initMatrizenKurz(this, dof)
          TYPE( LokaleMatrizen_t ), INTENT(inout) :: this
          INTEGER, INTENT(in) :: dof
          CALL initMatrizen(this, (/ dof, dof /) )
      END SUBROUTINE initMatrizenKurz
          
          ! DOCUMENTATION
          !> initMatrizen allokiert Speicherplatz für Systemmatrizen
          !!
          !> Um Types mit dynamische gro0en Elementen+
          !! zu erzeugen ist es nötig POINTER zu verwenden,
          !! allokierbare Arrays innerhalb von Type Definitionen
          !! nicht erlaubt sind. Die POINTER bekommenim folgenden
          !! ähnlich allokierbarer Arrays einen Speicherbereich
          !! ihrer Größe entsprechend zugewiesen.
      SUBROUTINE initMatrizen(this,maxElementMatrixSize)
          TYPE( LokaleMatrizen_t ), INTENT(inout):: this
          INTEGER, DIMENSION(2), INTENT(in):: maxElementMatrixSize
          INTEGER:: a,b
          INTEGER:: error

          a = maxElementMatrixSize(1)
          b = maxElementMatrixSize(2)

          IF (this%AllocationsDone) THEN
              DEALLOCATE(this%MASS, this%STIFF, this%FORCE)
              NULLIFY(this%MASS, this%STIFF, this%FORCE)
          ENDIF
          ALLOCATE(  this%MASS(a,b), &
                   this%STIFF(a,b), &
                   this%FORCE(a), stat=error )
          IF (error.ne.0) THEN
              CALL FATAL('initMatrizen','Allokation der Lokalen Matrizen &&
                      && fehlgeschlagen')
          ELSE
              this%AllocationsDone=.TRUE.
          ENDIF
      END SUBROUTINE initMatrizen

      SUBROUTINE schuettleAlleMatrizen(this)
          Type(LokaleMatrizen_t), INTENT(inout) :: this

          INTEGER:: a,b
          INTEGER :: i,j,p,q !Zählvariablen
          REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: M,S
          REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: F

          ALLOCATE(M(SIZE(this%MASS,1),SIZE(this%MASS,2)), &
              & S(SIZE(this%STIFF,1),SIZE(this%STIFF,2)), &
              & F(SIZE(this%FORCE)))
          
          a = SIZE(this%STIFF,1)
          b = SIZE(this%STIFF,2)

          M = this%MASS
          S = this%STIFF
          F = this%FORCE

          DO j=1, b/2
              q=2*j-1
              DO i=1, a/2
                   p = 2*i-1
                   this%MASS(p,q) = M(i,j)
                   this%MASS(p,q+1) = M(i,b/2+j)
                   this%MASS(p+1,q+1) = M(a/2+i,b/2+j)
                   this%MASS(p+1,q) = M(a/2+i,j)

                   this%STIFF(p,q) = S(i,j)
                   this%STIFF(p,q+1) = S(i,b/2+j)
                   this%STIFF(p+1,q+1) = S(a/2+i,b/2+j)
                   this%STIFF(p+1,q) = S(a/2+i,j)
               END DO
               this%FORCE(q)= F(j)
               this%FORCE(q+1) = F(a/2+j)
           END DO
       END SUBROUTINE schuettleAlleMatrizen

       
           !> \brief addiert zwei LokaleMatrizen
           !> indem die einzelenen Arrays (FORCE, STIFF, MASS) 
           !! einzeln und elementweise addiert werden
      TYPE(LokaleMatrizen_t) FUNCTION addiere(this,LM2) 
           TYPE(LokaleMatrizen_t), INTENT(in) :: this, LM2
           CALL setzeMatrizen(addiere, &
                             & this%FORCE + LM2%FORCE, &
                             & this%STIFF + LM2%STIFF, &
                             & this%MASS + LM2%MASS )
           Return 
      END FUNCTION addiere

          !> \}
END MODULE LokaleMatrizenClass
!> \}

    !Materialliste durch dynamisch verlinkte Pointerliste
        !> \defgroup MaterialClass MaterialClass
        !> \ingroup Elementeigenschaften 
              !!\{
        !> \class MaterialClass MaterialClass.f90 "MaterialClass.f90"
              !! \brief Struktur für einen Materialparameter
    !MODULE MaterialClass
    !      USE DefUtils
    !      
    !          ! DOCUMENTATION
    !          !> \struct MaterialClass::material_t
    !          !! Sammlung der Materialwerte für ein spezifisches Material
    !      TYPE materialParameter_t
    !          CHARACTER(len=80):: parameterName !< siehe SOLVER.KEYWORDS
    !          REAL(KIND=dp), POINTER, DIMENSION(:) :: wertSkalar=>NULL()
    !          REAL(KIND=dp), POINTER, DIMENSION(:,:) :: wertVektor=>NULL()
    !          TYPE(Materialparameter_t), POINTER :: &
    !              & naechsterMaterialparameter=>NULL()
    !      END TYPE
    !
    !      INTERFACE neu
    !          MODULE procedure initilisiereMaterialParameter
    !      END INTERFACE
    !      INTERFACE wert 
    !          MODULE Procedure setzeSkalar
    !          MODULE Procedure setzeVektor
    !      END INTERFACE wert
    !  CONTAINS
    !      SUBROUTINE initilisiereMaterialParameter(this, elmerName)
    !          TYPE(Materialparameter_t), INTENT(inout) :: this
    !          CHARACTER(*) :: elmerName
    !          this%parameterName = elmerName
    !      END SUBROUTINE initilisiereMaterialParameter
    !
    !      SUBROUTINE setzeSkalar(this, skalar)
    !          TYPE(Materialparameter_t), INTENT(inout) :: this
    !          REAL(KIND=dp), DIMENSION(:), INTENT(in) :: skalar
    !      END SUBROUTINE setzeSkalar
    !      
    !      SUBROUTINE setzeVektor(this, vektor)
    !          TYPE(Materialparameter_t), INTENT(inout) :: this
    !          REAL(KIND=dp), DIMENSION(:,:), INTENT(in) :: Vektor
    !      END SUBROUTINE setzeVektor
    !
    !      SUBROUTINE setzePointerAufNaechstenParameter(this, &
    !          & naechsterMaterialparameter)
    !          TYPE(Materialparameter_t), INTENT(inout) :: this
    !          TYPE(Materialparameter_t), POINTER :: &
    !              & naechsterMaterialparameter
    !          
    !          this%naechsterMaterialparameter => naechsterMaterialparameter
    !
    !      END SUBROUTINE setzePointerAufNaechstenParameter
    !END MODULE

        !> \class MaterialSammlungClass MaterialSammlungClass "GeneralSolverPackage.f90"
        !> \brief Sammelt die Materialeigenschaften in einer Baumstruktur
        !!
        !> materialSammlung_t ist die Wurzel einer Baumstruktur, welche die 
        !! Materialwerte für skalare, isotrope und anisotrope Parameter in
        !! linked lists ablegt. Ein beliebiger Wert kann aus den Linked lists
        !! geholt oder weitere hinzugefügt werden. Durch die Struktur muss
        !! die Anzahl der Parameter nicht vor der Programmausführung bekannt 
        !! sein
    !MODULE MaterialSammlungClass
    !      USE DefUtils
    !      USE MaterialClass
    !
    !          !> \struct ElementClass::material_t
    !          !! Sammlung der für die Simulation nötigen Materialkennwerte
    !          !!
    !          !> Die Materialwerte werden in einer Baum- bzw Pointerliste abgelegt
    !      TYPE materialSammlung_t
    !          PRIVATE
    !          TYPE(Materialparameter_t), POINTER, DIMENSION(:):: skalare
    !          TYPE(Materialparameter_t), POINTER, DIMENSION(:):: isotrop
    !          TYPE(Materialparameter_t), POINTER, DIMENSION(:):: anisotrop
    !          INTEGER :: elementDoF
    !      END TYPE materialSammlung_t
    !
    !      PRIVATE:: neueMaterialsammlung
    !    !     INTERFACE neu 
    !    !         MODULE Procedure neueMaterialsammlung
    !    !     END INTERFACE neu
    !    ! CONTAINS
    !    !     SUBROUTINE fuegeSkalarhinzu(this,parametername,wertSkalar)
    !    !         TYPE(materialSammlung_t), INTENT(inout) :: this
    !    !         CHARACTER(*), INTENT(in) :: parametername
    !    !         REAL(KIND=dp), POINTER, DIMENSION(:) :: wertSkalar
    !
    !    !         TYPE(materialParameter_t), save :: neuerParameter
    !
    !    !         !Ein neuer Parameter landet immer vorne in der Liste und 
    !    !         !wird nach skalar "hineingeschoben"
    !    !         CALL neu(neuerParameter, parametername) 
    !    !         CALL wert(neuerParameter, wertSkalar)
    !    !         CALL setzePointerAufNaechstenParameter(neuerParameter, &
    !    !             & this%skalare%naechsterMaterialparameter)
    !    !         this%skalare%naechsterMaterialparameter => neuerParameter
    !    !     END SUBROUTINE fuegeSkalarhinzu
    !
    !    !     SUBROUTINE fuegeIsotrophinzu(this,parametername,wertIsotrop)
    !    !         TYPE(materialSammlung_t), INTENT(inout) :: this
    !    !         CHARACTER(*), INTENT(in) :: parametername
    !    !         REAL(KIND=dp), POINTER, DIMENSION(:) :: wertIsotrop
    !    !     END SUBROUTINE fuegeIsotrophinzu
    !
    !    !     SUBROUTINE fuegeAnisotrophinzu(this,parametername,wertAnisotrop)
    !    !         TYPE(materialSammlung_t), INTENT(inout) :: this
    !    !         CHARACTER(*), INTENT(in) :: parametername
    !    !         REAL(KIND=dp), POINTER, DIMENSION(:) :: wertAnisotrop
    !    !     END SUBROUTINE fuegeAnisotrophinzu
    !
    !    !     
    !    !     SUBROUTINE neueMaterialsammlung(this, elementDoF)
    !    !         TYPE(materialSammlung_t), INTENT(inout) :: this
    !    !         INTEGER, INTENT(in) :: elementDoF
    !
    !    !         this%elementDoF = elementDoF
    !    !         ! Dummyargumente um einen Zugriff zu haben
    !
    !    !     END SUBROUTINE neueMaterialsammlung 
    !END MODULE MaterialSammlungClass

    !> \class NeumannRandbedingungClass RandbedingungClass.f90
    !! "RandbedingungClass.f90"
      !> \brief Einbinden von Neuman Randbedingungen
      !!
      !> Die Klasse wird dynamisch eingebunden, falls für das 
      !! Element Neuman-Randbedingungen definiert sind
      !!
      !> \todo
      !!  - Derived Type ausbauen
      !!  - Proceduren hinzufügen
MODULE NeumannRandbedingungClass
      USE DefUtils

      implicit none
          !> \struct NeumannRandbedingung::NeumannRandbedingung_t
          !! Beinhaltet den Wert der Neuman-Randbedingungen 
          !!
          !> \param Wert speichert den Wert der Neuman-Randbedingung 
          !!        für ein Element
      TYPE, PUBLIC :: NeumannRandbedingung_t
          PRIVATE
          Real(Kind=dp), DIMENSION(:), POINTER  :: Wert
      END TYPE
      PRIVATE
      PUBLIC:: setze, gibWert
      INTERFACE setze 
          MODULE PROCEDURE setzeWert
      END INTERFACE
      INTERFACE gibWert
          MODULE PROCEDURE gibWertNeumann
      END INTERFACE
  CONTAINS
          !> Setzt den Wert für eine Volumenlast für die 
          !! Knoten eines Elements
      SUBROUTINE setzeWert(this, Wert)
          TYPE(NeumannRandbedingung_t), INTENT(inout) :: this
          REAL(kind=dp), DIMENSION(:), TARGET :: Wert
          this%Wert => Wert
          WRITE(MESSAGE,*) 'Setze Neumann Randbedingung. wert = ', &
          & Wert, this%Wert 
          CALL INFO('NeumannRandbedingungClass', MESSAGE, level=9)
      END SUBROUTINE setzeWert
          !> Gibt einen Pointer auf ein Array mit den Knotenwerten
          !! für eine Volumenlast zurück
      FUNCTION gibWertNeumann(this) Result(wert)
          TYPE(NeumannRandbedingung_t), INTENT(in) :: this
          REAL(kind=dp), DIMENSION(:), POINTER :: Wert
          Wert => this%Wert
          Return
      END FUNCTION gibWertNeumann
END MODULE NeumannRandbedingungClass

      !> \class CauchyRandbedingungClass RandbedingungClass.f90
      !! "RandbedingungClass.f90"
      !> \brief Einbinden von Cauchy- Randbedingungen
      !!
      !> Diese Klasse wird dynamisch eingebunden, falls für 
      !! das Element Cauchy- Randbedingungen definiert sind
MODULE CauchyRandbedingungClass
      USE DefUtils
      USE MyStdModules

      implicit none
          !> \struct CauchyRandbedingungClass::CauchyRandbedingung_t
          !! Beinhaltet Werte der Cauchy Randbedingung
      TYPE, PUBLIC:: CauchyRandbedingung_t
          PRIVATE
          REAL(KIND=dp), DIMENSION(:), POINTER  :: bezugswert=>NULL()
          REAL(KIND=dp), DIMENSION(:), POINTER  :: skalierfaktor=>NULL()
      END TYPE
      PRIVATE
      PUBLIC:: setze, gibWert
      INTERFACE setze 
          MODULE PROCEDURE setzeWert
      END INTERFACE
      INTERFACE gibWert
          MODULE PROCEDURE gibWertCauchy
      END INTERFACE
  CONTAINS
          !> Setzt den Wert für eine Volumenlast für die 
          !! Knoten eines Elements
      SUBROUTINE setzeWert(this, WertRef, S)
          TYPE(CauchyRandbedingung_t), INTENT(inout) :: this
          REAL(kind=dp), DIMENSION(:), TARGET :: WertRef, S
          this%bezugswert => WertRef
          this%skalierfaktor => S
          WRITE(MESSAGE,*) 'Setze Cauchy Randbedingung. bezugswert = ', &
          & WertRef, this%bezugswert, ' skalierfaktor = ', s, this%skalierfaktor
          CALL INFO('CauchyRandbedingungClass', MESSAGE, level=9)
      END SUBROUTINE setzeWert
          !> Gibt einen Pointer auf ein Array mit den Knotenwerten
          !! für eine Volumenlast zurück
      FUNCTION gibWertCauchy(this,i) Result(wert)
          TYPE(CauchyRandbedingung_t), INTENT(in) :: this
          INTEGER, INTENT(in) :: i
          REAL(kind=dp), POINTER, DIMENSION(:):: Wert

          IF(i==1) THEN 
              Wert => this%bezugswert
          ELSE IF(i==2) THEN
              Wert => this%skalierfaktor
          ELSE
              Wert => NULL()
              CALL WARN('CauchyRandbedingungClass', 'gibWertCauchy(this,i)&
                  & i darf nur 1 für den bezugswert oder 2 für den &
                  & skalierfaktor sein. Alle anderen Werte enden in diesem &
                  & Fehler und Wert=>NULL()!')
          ENDIF

          Return
      END FUNCTION gibWertCauchy
END MODULE CauchyRandbedingungClass

      !> \class DirichletRandbedingung RandbedingungClass 
      !! "RandbedingungClass.f90"
      !> \brief Einbinden von Dirichlet- Randbedingungen
      !!
      !> Diese Klasse wird dynamisch eingebunden, falls für 
      !! das Element Dirichlet- Randbedingungen definiert sind
      !!
      !> \todo
      !!  - Derived Type ausbauen
      !!  - Proceduren hinzufügen
MODULE DirichletRandbedingungClass
      USE DefUtils

      implicit none
          !> \struct DirichletRandbedingung::DirichletRandbedingung_t
          !! Beinhaltet den Wert der DirichletRandbedingungen 
          !!
          !> \param Wert speichert den Wert der DirichletRandbedingung 
          !!        für ein Element
      TYPE, PUBLIC:: DirichletRandbedingung_t
          PRIVATE
          Real(Kind=dp), DIMENSION(:), POINTER  :: Wert
      END TYPE
      PRIVATE
      PUBLIC:: setze, gibWert
      INTERFACE setze 
          MODULE PROCEDURE setzeWert
      END INTERFACE
      INTERFACE gibWert
          MODULE PROCEDURE gibWertDirichlet
      END INTERFACE
  CONTAINS
          !> Setzt den Wert für eine Volumenlast für die 
          !! Knoten eines Elements
      SUBROUTINE setzeWert(this, Wert)
          TYPE(DirichletRandbedingung_t), INTENT(inout) :: this
          REAL(kind=dp), DIMENSION(:), TARGET :: Wert
          this%Wert => Wert
      END SUBROUTINE setzeWert
          !> Gibt einen Pointer auf ein Array mit den Knotenwerten
          !! für eine Volumenlast zurück
      FUNCTION gibWertDirichlet(this) Result(wert)
          TYPE(DirichletRandbedingung_t), INTENT(in) :: this
          REAL(kind=dp), DIMENSION(:), POINTER :: Wert
          Wert => this%Wert
          Return
      END FUNCTION gibWertDirichlet
END MODULE DirichletRandbedingungClass

      !> \class VolumenRandbedingungClass RandbedingungClass.f90 
      !!        "RandbedingungClass.f90"
      !> \brief Einbinden von Volumenlasten
MODULE VolumenRandbedingungClass
      USE DefUtils

      implicit none

          !> \struct VolumenRandbedingung::VolumenRandbedingung_t
          !! Beinhaltete Werte der Volumenquellen
      TYPE, PUBLIC:: VolumenRandbedingung_t
          PRIVATE
          Real(Kind=dp), DIMENSION(:), POINTER  :: Wert
      END TYPE
      PRIVATE
      PUBLIC:: setze, gibWert
      INTERFACE setze 
          MODULE PROCEDURE setzeWert
      END INTERFACE
      INTERFACE gibWert
          MODULE PROCEDURE gibWertVolumen
      END INTERFACE
  CONTAINS
          !> Setzt den Wert für eine Volumenlast für die 
          !! Knoten eines Elements
      SUBROUTINE setzeWert(this, Wert)
          TYPE(VolumenRandbedingung_t), INTENT(inout) :: this
          REAL(kind=dp), DIMENSION(:), TARGET :: Wert
          this%Wert => Wert
      END SUBROUTINE setzeWert
          !> Gibt einen Pointer auf ein Array mit den Knotenwerten
          !! für eine Volumenlast zurück
      FUNCTION gibWertVolumen(this) Result(wert)
          TYPE(VolumenRandbedingung_t), INTENT(in) :: this
          REAL(kind=dp), DIMENSION(:), POINTER :: Wert
          Wert => this%Wert
          Return
      END FUNCTION gibWertVolumen
END MODULE VolumenRandbedingungClass

     !> \class RandbedingungClass RandbedingungClass.f90
     !! "RandbedingungClass.f90"
     !> \brief Beschreibung der Randbedingungen
     !!
     !>  
     !! Die Klasse RandbedingungClass benutzt Runtime-Polymorphismus
     !! um mit einer Klasse verschiedene Randbedingungen zu nutzen. 
     !! Dazu wird ein Derived Type erzeugt, der Pointer auf alle 
     !! möglichen speziellen Randbedingungsklassen enthält. In dem 
     !! Modul werden auch alle in diesen Klassen vorhandenen class
     !! member functions (Prozeduren) überschrieben und beinhalten 
     !! nun Entscheidungsalgorithmen, welche Klassen aufgerufen werden.
     !! Der Vorteil von Runtime-Polymorphismus ist, dass die Klassen 
     !! , die Randbedingungen nutzen wollen, wie BauteilElementClass 
     !! oder Eriksson2006DGL, nicht wissen müssen, welche spezielle 
     !! Randbedingungen behandelt werden muss. Damit werden
     !! Aufgaben ausgelagert und besser strukturierbar. 
 MODULE RandbedingungClass
       USE DefUtils
       USE NeumannRandbedingungClass
       USE CauchyRandbedingungClass
       USE DirichletRandbedingungClass
       USE VolumenRandbedingungClass
       USE MyStdModules
       
       implicit none
 
       PRIVATE
       PUBLIC:: istVolumen, istNeumann, istCauchy, istDirichlet, &
           & gibRandbedingungDirichlet, gibRandbedingungNeumann, &
           & gibRandbedingungCauchy, gibRandbedingungVolumen, &
           & neu, setze
       TYPE, PUBLIC :: Randbedingung_t
           PRIVATE
           TYPE(VolumenRandbedingung_t), POINTER :: VolumenRandbedingung=>NULL()
           TYPE(NeumannRandbedingung_t), POINTER :: NeumannRandbedingung=>NULL()
           TYPE(CauchyRandbedingung_t), POINTER :: CauchyRandbedingung=>NULL()
           TYPE(DirichletRandbedingung_t), POINTER :: DirichletRandbedingung=>NULL()
       END TYPE
       INTERFACE neu 
           MODULE PROCEDURE initialisiereRandbedingungDirichlet
           MODULE PROCEDURE initialisiereRandbedingungNeumann
           MODULE PROCEDURE initialisiereRandbedingungCauchy
           MODULE PROCEDURE initialisiereRandbedingungVolumen
       END INTERFACE neu
       INTERFACE setze 
           MODULE PROCEDURE setzeRandbedingungSkalar
           MODULE PROCEDURE setzeRandbedingungCauchy
       END INTERFACE setze
   CONTAINS
       SUBROUTINE pruefeAssoziation(this)
           TYPE(Randbedingung_t), INTENT(in) :: this
           IF(ASSOCIATED(this%CauchyRandbedingung).or. &
                & ASSOCIATED(this%VolumenRandbedingung).or. &
                & ASSOCIATED(this%NeumannRandbedingung).or. &
                & ASSOCIATED(this%DirichletRandbedingung) ) THEN
                CALL WARN('RandbedingungClass::pruefeAssoziation','Achtung RB wird überschrieben')
            ENDIF
       END SUBROUTINE pruefeAssoziation
   
       SUBROUTINE initialisiereRandbedingungDirichlet(this,DirichletRB)  
           TYPE(Randbedingung_t), INTENT(inout) :: this
           TYPE(DirichletRandbedingung_t),INTENT(in), TARGET :: DirichletRB
   
           CALL pruefeAssoziation(this)
           this%DirichletRandbedingung => DirichletRB
           NULLIFY(this%CauchyRandbedingung)
           NULLIFY(this%VolumenRandbedingung)
           NULLIFY(this%NeumannRandbedingung)
           CALL INFO('RandbedingungClass', 'DirichletRandbedingung',level=8)
       END SUBROUTINE initialisiereRandbedingungDirichlet
 
       SUBROUTINE initialisiereRandbedingungNeumann(this,NeumannRB)  
           TYPE(Randbedingung_t), INTENT(inout) :: this
           TYPE(NeumannRandbedingung_t),INTENT(in), TARGET :: NeumannRB
   
           CALL pruefeAssoziation(this)
           this%NeumannRandbedingung => NeumannRB
           NULLIFY(this%CauchyRandbedingung)
           NULLIFY(this%VolumenRandbedingung)
           NULLIFY(this%DirichletRandbedingung)
           CALL INFO('RandbedingungClass','NeumannRandbedingung',level=8)
       END SUBROUTINE initialisiereRandbedingungNeumann
   
       SUBROUTINE initialisiereRandbedingungCauchy(this, CauchyRB)
           TYPE(Randbedingung_t), INTENT(inout) :: this
           TYPE(CauchyRandbedingung_t), INTENT(IN), TARGET :: CauchyRB
           CALL pruefeAssoziation(this)
           this%CauchyRandbedingung => CauchyRB
           NULLIFY(this%NeumannRandbedingung)
           NULLIFY(this%VolumenRandbedingung)
           NULLIFY(this%DirichletRandbedingung)
           CALL INFO('RandbedingungClass', 'CauchyRandbedingung initialisiert', level=8)
       END SUBROUTINE initialisiereRandbedingungCauchy
   
       SUBROUTINE initialisiereRandbedingungVolumen(this, VolumenRB)
           TYPE(Randbedingung_t), INTENT(inout) :: this
           TYPE(VolumenRandbedingung_t), INTENT(in), TARGET :: VolumenRB
           CALL pruefeAssoziation(this)
           this%VolumenRandbedingung => VolumenRB
           NULLIFY(this%CauchyRandbedingung)
           NULLIFY(this%NeumannRandbedingung)
           NULLIFY(this%DirichletRandbedingung)
           CALL INFO('RandbedingungClass', 'VolumenRandbedingung', level=8)
       END SUBROUTINE initialisiereRandbedingungVolumen

       LOGICAL FUNCTION istVolumen(this)
            TYPE(Randbedingung_t), INTENT(in) :: this
            IF(ASSOCIATED(this%VolumenRandbedingung)) THEN
                istVolumen = .true.
            ELSE IF(.not.ASSOCIATED(this%VolumenRandbedingung)) THEN
                istVolumen = .false.
            ELSE 
                CALL FATAL('RandbedingungClass', 'istVolumen kann Pointer Status &
                    & nicht bestimmen.')
            ENDIF
            Return 
       END FUNCTION istVolumen

       LOGICAL FUNCTION istNeumann(this)
            TYPE(Randbedingung_t), INTENT(in) :: this
            IF(ASSOCIATED(this%NeumannRandbedingung)) THEN
                istNeumann = .true.
            ELSE IF(.not.ASSOCIATED(this%NeumannRandbedingung)) THEN
                istNeumann = .false.
            ELSE 
                CALL FATAL('RandbedingungClass', 'istNeumann kann Pointer Status &
                    & nicht bestimmen.')
            ENDIF
            Return 
       END FUNCTION istNeumann

       LOGICAL FUNCTION istCauchy(this)
           TYPE(Randbedingung_t), INTENT(in) :: this
           IF(ASSOCIATED(this%CauchyRandbedingung)) THEN
               istCauchy = .true.
           ELSE IF(.not.ASSOCIATED(this%CauchyRandbedingung)) THEN
               istCauchy = .false.
           ELSE 
               CALL FATAL('RandbedingungClass', 'istCauchy kann Pointer Status &
                   & nicht bestimmen.')
           ENDIF
           Return 
       END FUNCTION istCauchy

       LOGICAL FUNCTION istDirichlet(this)
           TYPE(Randbedingung_t), INTENT(in) :: this
           IF(ASSOCIATED(this%DirichletRandbedingung)) THEN
               istDirichlet = .true.
           ELSE IF(.not.ASSOCIATED(this%DirichletRandbedingung)) THEN
               istDirichlet = .false.
           ELSE 
               CALL FATAL('RandbedingungClass', 'istDirichlet kann Pointer Status &
                   & nicht bestimmen.')
           ENDIF
           Return 
       END FUNCTION istDirichlet

       FUNCTION gibRandbedingungDirichlet(this)
            Type(Randbedingung_t), INTENT(in) :: this
            REAL(KIND=dp), DIMENSION(:), POINTER :: gibRandbedingungDirichlet
            IF(istDirichlet(this))THEN
                gibRandbedingungDirichlet => gibWert(this%DirichletRandbedingung)
            ELSE
                CALL WARN('RandbedingungClass:gibRandbedingungDirichlet',&
                    & 'Keine Dirichlet Randbedingung für Element &
                    & definiert.')
            ENDIF
            RETURN
       END FUNCTION gibRandbedingungDirichlet

       FUNCTION gibRandbedingungNeumann(this) RESULT(WertNeumann)
            Type(Randbedingung_t), INTENT(in) :: this
            REAL(KIND=dp), DIMENSION(:), POINTER :: WertNeumann

            IF(istNeumann(this)) THEN
                WertNeumann => gibWert(this%NeumannRandbedingung)
            ELSE
                WertNeumann => NULL()
            ENDIF
            RETURN

       END FUNCTION gibRandbedingungNeumann
       
       FUNCTION gibRandbedingungCauchy(this,i) RESULT(WertCauchy)
           Type(Randbedingung_t), INTENT(in) :: this
           INTEGER, INTENT(in) :: i
           REAL(KIND=dp), DIMENSION(:), POINTER ::WertCauchy
           IF(istCauchy(this)) THEN
               WertCauchy => gibWert(this%CauchyRandbedingung,i)
           ELSE
               WertCauchy => NULL()
           ENDIF
           RETURN
       END FUNCTION gibRandbedingungCauchy

       FUNCTION gibRandbedingungVolumen(this)
            Type(Randbedingung_t), INTENT(in) :: this
            REAL(Kind=dp), DIMENSION(:),POINTER :: gibRandbedingungVolumen

            IF (ASSOCIATED(this%VolumenRandbedingung)) THEN
                gibRandbedingungVolumen => gibWert(this%VolumenRandbedingung)
            ELSE 
                gibRandbedingungVolumen => NULL()
            ENDIF
            RETURN
       END FUNCTION gibRandbedingungVolumen
 
           !> setzt Dirichlet-, Neumann- oder VolumenRandbedingung 
       SUBROUTINE setzeRandbedingungSkalar(this, wert)
           Type(Randbedingung_t), INTENT(inout) :: this 
           REAL(kind=dp), INTENT(in), DIMENSION(:), target :: Wert
           WRITE(MESSAGE,*)'setzeRandbedingungSkalar ', Wert
           CALL INFO('RandbedingungClass', MESSAGE, level=9)
           IF(istNeumann(this)) THEN
               CALL setze(this%NeumannRandbedingung,Wert)
           ELSE IF(istDirichlet(this)) THEN
               CALL setze(this%DirichletRandbedingung,Wert)
           ELSE IF(istVolumen(this)) THEN
               CALL setze(this%VolumenRandbedingung, Wert)
               !this%VolumenRandbedingung%Wert=>Wert
           ELSE IF(istCauchy(this)) THEN
               CALL FATAL('RandbedingungClass:setzeRandbedingungSkalar', 'Skalar&
                   & reicht nicht aus für Cauchy Randbedingungen. Nutze Methode &
                   & setzeRandbedingungCauchy.')
           ELSE 
               CALL FATAL('RandbedingungClass:setzeRandbedingungSkalar', &
                   & 'Keine Randbedingungen für Element initialisiert. Kann &
                   & Wert nicht zuordnen.')
           ENDIF
       END SUBROUTINE setzeRandbedingungSkalar
           !> Setzt die nötigen Werte für Cauchy Randbedingungen
           !!
           !> Die nötigen Werte sind ein Skalierungsfaktor, hier 
           !! \f$S\f$, und ein Referenzwert \f$X_{ref}\f$. der
           !! Skalierungsfaktor korreliert die Wertdifferenz
           !! \f$X_{ref}-X\f$ mit einem Fluss \f$q\f$ durch den
           !! Gebietsrand.
           !! \f$ q=S\left(X_{ref}-X\right)\f$
           !! \param this Randbedingung 
           !! \param WertRef \f$X_{ref}\f$
           !! \param S Skalierungsfaktor
       SUBROUTINE setzeRandbedingungCauchy(this, WertRef, S)
           Type(Randbedingung_t), INTENT(inout) :: this 
           REAL(kind=dp), INTENT(in), DIMENSION(:), target :: WertRef, S
           IF(ASSOCIATED(this%CauchyRandbedingung)) THEN
               CALL setze(this%CauchyRandbedingung, WertRef, S)
           ELSE IF(.not.ASSOCIATED(this%CauchyRandbedingung)) THEN
               CALL WARN('RandbedingungClass:setzeRandbedingungCauchy', &
                   & 'Element hat keine Cauchy Randbedingung' )
           ENDIF
       END SUBROUTINE setzeRandbedingungCauchy
 END MODULE RandbedingungClass
