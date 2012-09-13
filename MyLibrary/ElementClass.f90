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
    !> \defgroup ElementClass ElementClass
    !> \ingroup Elementeigenschaften 
    !!\{
    !> \class ElementClass ElementClass.f90 "ElementClass.f90"
    !> \brief sammelt und verarbeitet elementbasierte
    !! Informationen
    !! 
    !> Elmers Simulationsablauf ist elementbasiert. Das heißt für jedes
    !! Element werden die Systemmatrizen berechnet und durch Elmers Routinen
    !! in die Gesammtmatrizen eingepflegt. Diese Klasse sammelt und verarbeitet
    !! die für die Komposition der Lokalen Matrizen notwendigen
    !! Informationen, wie 
    !!  - Randbedingung
    !!  - Elementmatrizen
    !!  - Algorithmus zum berechnen der Elementmatrizen 
    !> \todo
    !!  - Geometrie und Materialwerte nur beim ersten Aufruf pro 
    !!    Zeitschritt setzen 
 MODULE ElementClass
      USE DefUtils
      USE LokaleGekoppelteMatrizenClass 
      USE MaterialSammlung_Class
      USE geometrischeEigenschaften_Class
      USE RandbedingungClass
      USE DirichletRandbedingungClass
      USE CauchyRandbedingungClass
      USE NeumannRandbedingungClass
      USE VolumenRandbedingungClass
      USE Eriksson2006DGL_Class
      USE MyStdModules
      implicit none

          !> Globale Klassenvariable
          !!
          !> Diese Variablen sind für alle Objekte vom Typ 
          !! BauteilElement_t les- und schreibbar
          !> \param elementezaehler zählt die Anzahl der erzeugten 
          !! Objekte 
          !> \param durchlaufzaehler zählt die Anzahl der Durchläufe 
          !! durch die linked list um festzustellen wann alle Elemente 
          !! berechnet wurden wird elementezaehler mit durchlaufzaehler
          !! verglichen
      INTEGER, private :: elementezaehler=0, durchlaufzaehler=1
      LOGICAL, private :: found

          !! \struct ElementClass::BauteilElement_t
          !! Type Struktur als Kollektion der Elementeigenschaften
          !!
          !> sammelt Elementeigenschaften
          !> \param elementNumber vielleicht überflüssig
          !! \param istGebietselement seperiert in 
          !!        Gebietselemente und Randgebietselemente
          !> \todo 
          !!  - Aufbau des structs überdenken
          !!  - DGL
      TYPE BauteilElement_t
          INTEGER:: elementNummer
          LOGICAL:: istGebietselement
          INTEGER :: elementDoF
          TYPE(LokaleGekoppelteMatrizen_t):: elementMatrizen 
          Type(BauteilElement_t), POINTER:: naechstesElement
          TYPE(geometrischeEigenschaften_t):: geometrie
          TYPE(materialSammlung_t) :: materialWerte
          TYPE(Randbedingung_t):: temperaturRandbedingung
          TYPE(Randbedingung_t):: feuchtigkeitRandbedingung
          TYPE(Eriksson2006DGL_t) :: ErikssonDGL
          TYPE(Element_t), POINTER :: element !< Elmers Elementinfo
          TYPE(ValueList_t), Pointer :: material=>NULL() !< Elmers Materialinfo
          TYPE(Solver_t) :: Solver
          TYPE(Model_t) :: Model
          REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: vorherigeLsg

          ! Probleme mit weglaufenden Pointern zwingen mich alle 
          ! für Randbedingungen nötigen Variablen hier zu definieren
          REAL(kind=dp), DIMENSION(:), POINTER, PRIVATE :: moisup, &
              & heatsup, moiflu, heatflu, mcsurf, tsurf, smc, st
          TYPE( NeumannRandbedingung_t ), PRIVATE :: NeumannRBTemp
          TYPE( CauchyRandbedingung_t ) , PRIVATE :: CauchyRBTemp
          TYPE( VolumenRandbedingung_t ), PRIVATE :: VolumenRBTemp
          TYPE( NeumannRandbedingung_t ), PRIVATE :: NeumannRBFeucht
          TYPE( CauchyRandbedingung_t ) , PRIVATE :: CauchyRBFeucht
          TYPE( VolumenRandbedingung_t ), PRIVATE :: VolumenRBFeucht

      END TYPE

      PRIVATE:: initialisiereBauteilElement, berechneElement, &
          & initialisiereMaterial, initialisiereRandbedingungen, &
          & holeMaterialVektor

      interface neu
          MODULE PROCEDURE initialisiereBauteilElement
      END interface
      interface berechne 
          MODULE PROCEDURE berechneElement
      END interface
      interface holeMaterial 
          MODULE PROCEDURE holeMaterialVektor
          MODULE PROCEDURE holeMaterialTensor
      END interface
  CONTAINS
      SUBROUTINE initialisiereBauteilElement(this, ElemNr, Solver, Model, &
          & istGebietselement)
          TYPE(BauteilElement_t):: this
          INTEGER :: ElemNr
          TYPE(Model_t) :: Model
          TYPE(Solver_t) :: Solver
          LOGICAL :: istGebietselement 

          this%elementNummer=ElemNr 
          this%solver = Solver
          this%Model = Model
          this%istGebietselement = istGebietselement

          IF(this%istGebietselement) THEN
              this%element => GetActiveElement(this%elementNummer,&
                  & usolver= this%solver)
              this%material => GetMaterial(this%element)
          ELSE IF(.not.this%istGebietselement) THEN
              this%element => GetBoundaryElement(&
                  & this%ElementNummer, usolver=this%solver)
          ELSE
              CALL FATAL('ElementClass: ', 'Elementtype ist unbestimmt.')
          ENDIF
          
          CALL initialisiereGeometrie(this)

          ! GetElementNOFDOFs gibt Größe der Submatrizen
          this%elementDoF = 2*GetElementNOFDOFs(this%Element)
          ALLOCATE(this%vorherigeLsg(2,GetElementNOFNodes(this%Element)))
          CALL GetVectorLocalSolution(this%vorherigeLsg, UElement=this%Element)


          !Definiere neue Matrizen für dieses Element
          CALL neu(this%elementMatrizen, this%elementDoF)
          
              ! Pointer zeigt möglicherweise ins Nirvana.
              ! Daher setze ich ihn erstmal auf NULL
          NULLIFY(this%naechstesElement)

          !print *,'Elementezähler: ', elementezaehler
          ! initialisier Randbedingung entsprechend istGebietselement
          CAll initialisiereRandbedingungen(this)
          ! Randbedingung bekommen ihre Werte über initialisiereRandbedingungen
          IF (this%istGebietselement) THEN
              ! initialisere Materialkennwerte für Gebietselemente
              CALL initialisiereMaterial(this)
          ENDIF

          ! am Ende wird noch die DGL intitialisiert
          CALL initialisiereDGL(this)
      END SUBROUTINE initialisiereBauteilElement

      SUBROUTINE setzeElementzaehler(anzahlElemente)
          TYPE(BauteilElement_t):: this
          INTEGER :: anzahlElemente
          elementezaehler = AnzahlElemente
      END SUBROUTINE setzeElementzaehler

          ! DOCUMENTATION
          !> initialisiert die Materialparameter für 
          !! die Berechnung der Eriksson2006DGL_Class
          !!
          !> Dabei sind folgende Materialeigenschaften für die
          !! Berechnung der Matrizen notwendig
          !! materialliste:
          !! - density \f$\rho\f$
          !! - dry density \f$\rho_s\f$
          !! - Activation energy for bound water \f$E_b\f$
          !! - Feuchtigkeit
          !!   - Mass transfer coeff
          !!   - diffusion coeff \f$D_{\omega}(X,T)\f$
          !! - Temperatur 
          !!   - Konduktivität \f$D_T(X,T)\f$
          !!   - specific heat \f$c\f$
      SUBROUTINE initialisiereMaterial(this)
          TYPE( BauteilElement_t ), INTENT(inout) :: this

          CALL holeMaterial(this, this%materialWerte%dichte, 'density')
          CALL holeMaterial(this, this%materialWerte%spezifischeWaerme,'specific heat capacity')

          CALL holeMaterial(this, this%materialWerte%diffusionskoeffizientTensor, 'Diffusivity')
          CALL holeMaterial(this, this%materialWerte%konduktivitaetTensor, 'Conductivity')

          CALL holeMaterial(this, this%materialWerte%E_b, 'Eb')
          CALL holeMaterial(this, this%materialWerte%R, 'R')
          CALL holeMaterial(this, this%materialWerte%relativeHumidity, 'RH')
          CALL holeMaterial(this, this%materialWerte%dwdH, 'dwdH')

          Write(MESSAGE, *) 'density',this%materialWerte%dichte
          CALL INFO('ElementClass:initialisiereMaterial', MESSAGE, level=9)
          Write(MESSAGE, *) &
              & 'spezifischeWaerme',this%materialWerte%spezifischeWaerme
          CALL INFO('ElementClass:initialisiereMaterial', MESSAGE, level=9)
          Write(MESSAGE, *) 'Diffusivity',&
              & this%materialWerte%diffusionskoeffizientTensor(:,:,1)
          CALL INFO('ElementClass:initialisiereMaterial', MESSAGE, level=9)
          Write(MESSAGE, *) 'Conductivity',&
              & this%materialWerte%konduktivitaetTensor(:,:,1)
          CALL INFO('ElementClass:initialisiereMaterial', MESSAGE, level=9)
          Write(MESSAGE, *) 'Eb',this%materialWerte%E_b
          CALL INFO('ElementClass:initialisiereMaterial', MESSAGE, level=9)
          Write(MESSAGE, *) 'R',this%materialWerte%R
          CALL INFO('ElementClass:initialisiereMaterial', MESSAGE, level=9)
          Write(MESSAGE, *) 'RH',this%materialWerte%relativeHumidity
          CALL INFO('ElementClass:initialisiereMaterial', MESSAGE, level=9)
          Write(MESSAGE, *) 'dwdH',this%materialWerte%dwdH
          CALL INFO('ElementClass:initialisiereMaterial', MESSAGE, level=9)
      END SUBROUTINE initialisiereMaterial

          ! DOCUMENTATION
          !> initialisiert die Werte für die Geometrie für 
          !! die Berechnung der Eriksson2006DGL_Class
          !!
          !> Dabei sind folgende Werte für die
          !! Berechnung der Matrizen notwendig:
          !! - detJ
          !! - ds
          !! - Number of Integrationpoints
          !! - Formfunktionen \f$N\f$
          !! - abgeleitete Formfunktionen \f$B\f$
          !! - anzahlKnoten
      SUBROUTINE initialisiereGeometrie(this)
          TYPE( BauteilElement_t ), INTENT(inout) :: this

          TYPE(ValueList_t), POINTER :: SolverParams


          ! setze Element auch in Geometrieeigenschaften_t für Eriksson2006DGL_Class
          this%geometrie%Element => this%element
          this%geometrie%anzahlKnoten=GetElementNOFNodes(this%Element)

          SolverParams => GetSolverParams(this%solver)
          
          this%geometrie%Knotenfreiwerte = &
              & ListGetInteger(SolverParams,'Variable DOFs', found, 1,3)
          IF(.not.found) CALL FATAL('ElementClass: ', 'Knotenfreiwerte &
              & unbekannt. Kann Variable DOFS in Solver Bereich der &
              & .sif Datei nicht finden.')
          this%geometrie%DIM = CoordinateSystemDimension()
          CALL GetElementNodes( this%geometrie%Nodes,usolver=this%Solver,&
              & UElement=this%element)
          this%geometrie%IP = GaussPoints( this%Element )
          ALLOCATE(&
             & this%geometrie%N(1,this%geometrie%anzahlKnoten), &
             & this%geometrie%B(this%geometrie%anzahlKnoten, &
             &                  this%geometrie%DIM), &
             & this%geometrie%dBdx(this%geometrie%anzahlKnoten, &
             &                      this%geometrie%DIM, &
             &                      this%geometrie%DIM ) )
      END SUBROUTINE initialisiereGeometrie

      SUBROUTINE initialisiereRandbedingungen(this)
          TYPE( BauteilElement_t ), INTENT(inout) :: this

          !TYPE( NeumannRandbedingung_t ), TARGET, SAVE :: NeumannRBTemp
          !TYPE( CauchyRandbedingung_t ), TARGET , SAVE :: CauchyRBTemp
          !TYPE( VolumenRandbedingung_t ), TARGET, SAVE :: VolumenRBTemp
          !TYPE( NeumannRandbedingung_t ), TARGET, SAVE :: NeumannRBFeucht
          !TYPE( CauchyRandbedingung_t ), TARGET , SAVE :: CauchyRBFeucht
          !TYPE( VolumenRandbedingung_t ), TARGET, SAVE :: VolumenRBFeucht

          TYPE(ValueList_t), POINTER :: Liste
          CHARACTER(LEN=MAX_NAME_LEN) :: BoundaryType


          WRITE(MESSAGE, *)'initialisiereRandbedingungen:: Element ', this%elementNummer
          CALL INFO('ElementClass', MESSAGE, level=5)

          ! GEBIETSELEMENT
          IF (this%istGebietselement) THEN
              ! Schaue ob Volumenlast vorhanden
              Liste => GetBodyForce(this%Element, found)
              IF (ASSOCIATED(Liste).and.found) THEN
                  ! schaue ob für das Element Feuchtigkeitslast vorhanden
                  ALLOCATE(this%moisup(this%geometrie%anzahlKnoten))
                  this%moisup = & 
                      & GetReal(Liste, 'MoistureSupply', &
                          & found,UElement=this%element)
                  IF(found) THEN
                      CALL neu(this%feuchtigkeitRandbedingung, this%VolumenRBFeucht)
                      CALL setze(this%feuchtigkeitRandbedingung, this%moisup)
                      CALL INFO('ElementClass','setze BodyForce Moisture', level=6)
                  ENDIF
                  ! oder ob Temperaturvolumenlast vorhanden ist
                  ALLOCATE(this%heatsup(this%geometrie%anzahlKnoten))
                  this%heatsup = & 
                      & GetReal(Liste, 'HeatSource', found, UElement=this%element)
                  IF(found) THEN
                      CALL neu(this%temperaturRandbedingung, this%VolumenRBTemp)
                      CALL setze(this%temperaturRandbedingung, this%heatsup)
                      CALL INFO('ElementClass','setze BodyForce Temperatur', level=6)
                  ENDIF
              ENDIF

          ! RANDELEMENT
          ELSE IF (.not.this%istGebietselement) THEN
              ! randbedingungen für element prüfen und setzen
              Liste=>GetBC(UElement=this%element)
              !> \todo Werte für randbedingungen einlesen
              !Feuchtigkeit Neumann oder Cauchy
              BoundaryType = GetString(Liste, 'Moisture Boundary Type', found)
              IF(found) THEN
                  IF(BoundaryType.eq.'neumann') THEN
                      ALLOCATE(this%moiflu(this%geometrie%anzahlKnoten))
                      this%moiflu = GetREAL(Liste, 'moistureflux',found)
                      IF(found) THEN
                          CALL neu(this%feuchtigkeitRandbedingung, this%NeumannRBFeucht)
                          CALL Setze(this%feuchtigkeitRandbedingung,this%moiflu)
                          CALL INFO('ElementClass','setze Neumann Moisture', level=6)
                      ELSE
                          WRITE(MESSAGE,*)'In Boundary Condition ', GetBCID(), &
                               & ' wurde der nötige Wert "MoistureFlux" nicht &
                               & gefunden. Dieser ist nötig, wenn "neumann" als&
                               & "Moisture Boundary Type" gesetzt ist.'
                          CALL WARN('ElementClass:initialisiereRandbedingungen',&
                              &MESSAGE)
                      ENDIF
                  ELSE IF(BoundaryType.eq.'cauchy') THEN
                      ALLOCATE(this%mcsurf(this%geometrie%anzahlKnoten))
                      this%mcsurf = GetReal(Liste, 'MC_SURF', found)
                      IF(found) THEN
                          ALLOCATE(this%smc(this%geometrie%anzahlKnoten))
                          this%smc = GetReal(Liste, 'S_MC', found)
                          IF(found) THEN
                              CALL neu(this%feuchtigkeitRandbedingung, this%CauchyRBFeucht)
                              CALL setze(this%feuchtigkeitRandbedingung, &
                                  & this%mcsurf,this%smc) 
                              CALL INFO('ElementClass','setze Cauchy Moisture', level=6)
                          ELSE
                              WRITE(MESSAGE,*) 'S_MC fehlt bei &
                                  & Cauchy-Randbedingung in Boundary Condition &
                                  &', GetBCID()
                                  CALL WARN('ElementClass:initialisiereRandbedingungen',&
                                      & MESSAGE)
                          ENDIF
                      ELSE
                          WRITE(MESSAGE, *) 'MC_Surf fehlt bei Cauchy-Randbedingung&
                              & in Boundary Condition ', GetBCID()
                          CALL WARN('ElementClass:initialisiereRandbedingungen',&
                              & MESSAGE)
                      ENDIF
                  ELSE 
                      WRITE(MESSAGE, *) 'Unbekannter Randbedingungstyp ', &
                          & BoundaryType, ' für Feuchtigkeits&
                          &randbedingung "Moisture Boundary Type" in Boundary &
                          & Condition ', GetBCID(), ' in der sif Datei &
                          & definiert. Bisher ist nur Neumann und Cauchy für &
                          & diesen Solver unterstützt.'
                      CALL WARN('ElementClass:initialisiereRandbedingungen',&
                          & MESSAGE)
                  ENDIF
              ENDIF
              !Temperatur Neumann oder Cauchy
              BoundaryType = GetString(Liste, 'Temperature Boundary Type', found)
              IF(found) THEN
                  IF(BoundaryType.eq.'neumann') THEN
                      ALLOCATE(this%heatflu(this%geometrie%anzahlKnoten))
                      this%heatflu = GetREAL(Liste, 'heatflux',found)
                      IF(found) THEN
                          CALL neu(this%temperaturRandbedingung, this%NeumannRBTemp)
                          CALL Setze(this%temperaturRandbedingung, this%heatflu)
                          CALL INFO('ElementClass','setze Neumann Heat', level=6)
                      ELSE
                          WRITE(MESSAGE,*)'In Boundary Condition ', GetBCID(), &
                               & ' wurde der nötige Wert "HeatFlux" nicht &
                               & gefunden. Dieser ist nötig, wenn "neumann" als&
                               & "Moisture Boundary Type" gesetzt ist.'
                          CALL WARN('ElementClass:initialisiereRandbedingungen',&
                              &MESSAGE)
                      ENDIF
                  ELSE IF(BoundaryType.eq.'cauchy') THEN
                      ALLOCATE(this%tsurf(this%geometrie%anzahlKnoten))
                      this%tsurf = GetReal(Liste, 'T_SURF', found, &
                          & UElement=this%element)
                      IF(found) THEN
                          ALLOCATE(this%st(this%geometrie%anzahlKnoten))
                          this%st = GetReal(Liste, 'S_T', found,&
                              & UElement=this%element)
                          IF(found) THEN
                              CALL neu(this%temperaturRandbedingung, this%CauchyRBTemp)
                              CALL setze(this%temperaturRandbedingung, &
                                  & this%tsurf,this%st) 
                              CALL INFO('ElementClass','setze Cauchy HEat', level=6)
                          ELSE
                              WRITE(MESSAGE,*) 'S_T fehlt bei &
                                  & Cauchy-Randbedingung in Boundary Condition &
                                  &', GetBCID()
                                  CALL WARN('ElementClass:initialisiereRandbedingungen',&
                                      & MESSAGE)
                          ENDIF
                      ELSE
                          WRITE(MESSAGE, *) 'T_Surf fehlt bei Cauchy-Randbedingung&
                              & in Boundary Condition ', GetBCID()
                          CALL WARN('ElementClass:initialisiereRandbedingungen',&
                              & MESSAGE)
                      ENDIF
                  ELSE 
                      WRITE(MESSAGE, *) 'Unbekannter Randbedingungstyp ', &
                          & BoundaryType, ' für Temperatur&
                          &randbedingung "Temperature Boundary Type" in Boundary &
                          & Condition ', GetBCID(), ' in der sif Datei &
                          & definiert. Bisher ist nur Neumann und Cauchy für &
                          & diesen Solver unterstützt.'
                      CALL WARN('ElementClass:initialisiereRandbedingungen',&
                          & MESSAGE)
                  ENDIF
              ENDIF
              !Check for Dirichlet RB --> Werden nach der Assemblierung 
              ! wegen dieser Prüfung ist dieses Package nicht 
              ! in der Lage dynamisch Elemente an- und abzuschalten
              !> \todo diese Abfrage optimieren
              ! von BauteilClass über eine Elmermethode gesetzt
          ENDIF
      END SUBROUTINE initialisiereRandbedingungen

          !> Übergibt die für dieses Element spezifischen Material- und
          !! Geometrieeigenschaften an die DGL, sodass diese ihre Gleichungen
          !! für dieses Element optimieren kann
      SUBROUTINE initialisiereDGL(this)
          TYPE(BauteilElement_t), INTENT(inout) :: this
          !> \todo Methodenaufruf überdenken
          IF(this%istGebietselement) THEN
              ! Nutzt überladen Funktion um entsprechend
              ! für ein Gebietselement zu initilisieren
              CALL neu(this%ErikssonDGL,&
                      &this%elementMatrizen,&
                      &this%materialWerte,&
                      &this%feuchtigkeitRandbedingung,&
                      &this%temperaturRandbedingung, &
                      &this%Geometrie)
          ELSE
              ! nutzt überladene Funktion um entsprechend
              ! für ein Randgebietselemente zu initilisieren
              ! beachte die andere Argumentliste im vergleich
              ! zum oberen Aufruf
              CALL neu(this%ErikssonDGL,&
                      &this%elementMatrizen,&
                      &this%feuchtigkeitRandbedingung,&
                      &this%temperaturRandbedingung, &
                      &this%Geometrie)
          ENDIF
      END SUBROUTINE initialisiereDGL

          ! DOCUMENTATION
          !> \brief Gibt mithilfe von Elmers Routine den Wert elmername zurück
          !! und macht zusätzliche eine Fehlerprüfung
      SUBROUTINE holeMaterialVektor(this, materialWert, elmername)
          TYPE(BauteilElement_t), INTENT(inout) :: this
          REAL(kind=dp), DIMENSION(:), POINTER :: materialWert
          CHARACTER(*), INTENT(in) :: elmername 

          ALLOCATE(materialWert(this%geometrie%anzahlKnoten))
          materialWert=&
              & GetReal(this%material, elmername, Found, UElement=this%element)
          IF (.not.Found) THEN
              print *, 'Fehlende Materialeigenschaft für Body ', &
                  & this%element%BodyID
              WRITE(MESSAGE, *)elmername, ' nicht in sif Datei gefunden!'
              CALL FATAL('ElementClass',MESSAGE)
          ENDIF
      END SUBROUTINE holeMaterialVektor

      SUBROUTINE holeMaterialTensor(this, materialWert, elmername)
          TYPE(BauteilElement_t), INTENT(inout) :: this
          REAL(kind=dp), DIMENSION(:,:,:), POINTER, INTENT(inout) :: materialWert
          CHARACTER(*), INTENT(in) :: elmername

          REAL(kind=dp), DIMENSION(:,:,:), POINTER:: dummy

          CALL ListGetRealArray(this%material,elmername,dummy, &
              & 8,this%element%NodeIndexes &
              & , found)
          ALLOCATE(materialWert(SIZE(dummy,1),SIZE(dummy,2), SIZE(dummy,3)))
          CALL ListGetRealArray(this%material,elmername,materialWert, &
              & 8,this%element%NodeIndexes &
              & , found)
          IF (.not.Found) THEN
              print *, 'Fehlende Materialeigenschaft für Body ', &
                  & this%element%BodyID
              WRITE(MESSAGE, *)elmername, ' nicht in sif Datei gefunden!'
              CALL FATAL('ElementClass',MESSAGE)
          ENDIF
      END SUBROUTINE holeMaterialTensor

          !> Setze den Pointer auf das Nachbarelement
          !!
          !> Normalerweise ist der Pointer auf das 
          !! nächste Element auf NULL gesetzt und 
          !! muss erst gesetzt werden, um die Möglichkeit
          !! zu geben durch alle Elementen als eine 
          !! Linked List zu iterieren.
      SUBROUTINE setzePointerAufNaechstesElement(this, PointerAufNaechstesElement)
          TYPE(BauteilElement_t) :: this
          TYPE(BauteilElement_t), POINTER :: PointerAufNaechstesElement

          this%naechstesElement => PointerAufNaechstesElement
          !print *, 'Element ', this%elementNummer, 'zeigt zu Element ', &
              !& this%naechstesElement%elementNummer
      END SUBROUTINE setzePointerAufNaechstesElement
      
          !> Löse die DGL für dieses Element und
          !! folge dem POINTER auf das naechste Element 
          !!
          !> Da jedes Element die Möglichkeit hat sein 
          !! Nachbarelement zu kennen kann eine Linked
          !! List erzeugt werden, die es ermöglicht über 
          !! alle Elemente zu iterieren. Der Pointer wird
          !! durch die SUBROUTINE setzePointerAufNaechstesElement
          !! gesetzt. Wird diese SUBROUTINE aufgerufen werden
          !! zuerst die Systemmatrizen für dieses Element 
          !! berechnet und dann wird zum nächsten Element 
          !! gegangen. Dort werden die Matrizen berechnet, usw.,
          !! bis das Letzte Element keinen Pointer mehr auf 
          !! ein weiteres Element mehr hat.
          !> \sa setzePointerAufNaechstesElement
      RECURSIVE SUBROUTINE iteriereUeberAlleElemente(this)
          TYPE(BauteilElement_t) :: this

          WRITE(MESSAGE, *) 'Element ', this%elementNummer, 'is happy to welcome &
              & you at RECURSIVE SUBROUTINE iteriereUeberAlleElemente mit &
              & dem durchlaufzaehler von : ', durchlaufzaehler, '/',elementezaehler
          CALL INFO('ElementClass', MESSAGE, level=10)
          CALL berechne(this)
          durchlaufzaehler = durchlaufzaehler + 1
          !print *, 'durchlaufzaehler', durchlaufzaehler, 'elementezaehler',&
              !& elementezaehler
          IF(.not.(durchlaufzaehler.GT.elementezaehler)) THEN
              IF(ASSOCIATED(this%naechstesElement)) THEN
                  CALL iteriereUeberAlleElemente(this%naechstesElement)
              ELSE
                  CALL FATAL('ElementClass', 'Die Linked List ist nicht &
                      & geschlossen')
              ENDIF
          ELSE
              ! setze durchlaufzaehler zurück auf 0 um für die nächste
              ! Berechnungsdurchläufe wieder alle Elemente einzubeziehen
              durchlaufzaehler = 1
          ENDIF
          RETURN
      END SUBROUTINE iteriereUeberAlleElemente
          !> Löse die DGL für dieses Element und
          !! folge dem POINTER auf das naechste Element 
          !!
          !> Da jedes Element die Möglichkeit hat sein 
          !! Nachbarelement zu kennen kann eine Linked
          !! List erzeugt werden, die es ermöglicht über 
          !! alle Elemente zu iterieren. Der Pointer wird
          !! durch die SUBROUTINE setzePointerAufNaechstesElement
          !! gesetzt. Wird diese SUBROUTINE aufgerufen werden
          !! zuerst die Systemmatrizen für dieses Element 
          !! berechnet und dann wird zum nächsten Element 
          !! gegangen. Dort werden die Matrizen berechnet, usw.,
          !! bis das Letzte Element keinen Pointer mehr auf 
          !! ein weiteres Element mehr hat.
          !> \sa setzePointerAufNaechstesElement
      RECURSIVE SUBROUTINE iteriereUeberAlleElementeStaggered1(this)
          TYPE(BauteilElement_t) :: this

          !print *, 'Element ', this%elementNummer, 'is happy to welcome &
          !    & you at RECURSIVE SUBROUTINE iteriereUeberAlleElemente mit &
          !    & dem durchlaufzaehler von : ', durchlaufzaehler
          CALL berechne(this)
          durchlaufzaehler = durchlaufzaehler + 1
          !print *, 'durchlaufzaehler', durchlaufzaehler, 'elementezaehler',&
              !& elementezaehler
          IF(.not.(durchlaufzaehler.GE.elementezaehler)) THEN
              IF(ASSOCIATED(this%naechstesElement)) THEN
                  CALL iteriereUeberAlleElemente(this%naechstesElement)
              ELSE
                  CALL FATAL('ElementClass', 'Die Linked List ist nicht &
                      & geschlossen')
              ENDIF
          ELSE
              ! setze durchlaufzaehler zurück auf 0 um für die nächste
              ! Berechnungsdurchläufe wieder alle Elemente einzubeziehen
              durchlaufzaehler = 1
          ENDIF
          RETURN
      END SUBROUTINE iteriereUeberAlleElementeStaggered1
          !> Löse die DGL für dieses Element und
          !! folge dem POINTER auf das naechste Element 
          !!
          !> Da jedes Element die Möglichkeit hat sein 
          !! Nachbarelement zu kennen kann eine Linked
          !! List erzeugt werden, die es ermöglicht über 
          !! alle Elemente zu iterieren. Der Pointer wird
          !! durch die SUBROUTINE setzePointerAufNaechstesElement
          !! gesetzt. Wird diese SUBROUTINE aufgerufen werden
          !! zuerst die Systemmatrizen für dieses Element 
          !! berechnet und dann wird zum nächsten Element 
          !! gegangen. Dort werden die Matrizen berechnet, usw.,
          !! bis das Letzte Element keinen Pointer mehr auf 
          !! ein weiteres Element mehr hat.
          !> \sa setzePointerAufNaechstesElement
      RECURSIVE SUBROUTINE iteriereUeberAlleElementeStaggered2(this)
          TYPE(BauteilElement_t) :: this

          !print *, 'Element ', this%elementNummer, 'is happy to welcome &
          !    & you at RECURSIVE SUBROUTINE iteriereUeberAlleElemente mit &
          !    & dem durchlaufzaehler von : ', durchlaufzaehler
          CALL berechne(this)
          durchlaufzaehler = durchlaufzaehler + 1
          !print *, 'durchlaufzaehler', durchlaufzaehler, 'elementezaehler',&
              !& elementezaehler
          IF(.not.(durchlaufzaehler.GE.elementezaehler)) THEN
              IF(ASSOCIATED(this%naechstesElement)) THEN
                  CALL iteriereUeberAlleElemente(this%naechstesElement)
              ELSE
                  CALL FATAL('ElementClass', 'Die Linked List ist nicht &
                      & geschlossen')
              ENDIF
          ELSE
              ! setze durchlaufzaehler zurück auf 0 um für die nächste
              ! Berechnungsdurchläufe wieder alle Elemente einzubeziehen
              durchlaufzaehler = 1
          ENDIF
          RETURN
      END SUBROUTINE iteriereUeberAlleElementeStaggered2

          !> Löse die definierte DGL für dieses Element
      SUBROUTINE berechneElement(this)
           TYPE(BauteilElement_t), INTENT(inout) :: this

           WRITE(MESSAGE, *) 'Berechne Elementmatrizen von Element: ', this%elementNummer
           CALL INFO('ElementClass', MESSAGE, level=5)

           ! nullifiziere Elementmatrizen vor der Berechnung
           CALL nullifiziere(this%elementMatrizen)
           IF(.not.this%istGebietselement)THEN
               !print *, 'ELEMENT ', this%elementNummer, &
               !    & ActiveBoundaryElement(UElement=this%element), &
               !    & GetElementFamily(this%element)
               IF(.not.ActiveBoundaryElement(&
                   & UElement=this%element )&
                   & .or. GetElementFamily(this%element) == 1)THEN
                   CALL INFO('ElementClass', 'nichts zu tun für Randgebietselemente',&
                      & level = 10)
               ELSE
               !IF(.true.)THEN
                   CALL berechne(this%ErikssonDGL,&
                       & this%vorherigeLsg,&
                       & this%istGebietselement)

                     !IF(.not.(MAXVAL(ABS(MASS(this%elementMatrizen)))<=0.1d-300))THEN 
                     !    CALL printArray(&
                     !        & MASS(this%elementMatrizen), 'MASS')
                     !ENDIF
                     !IF(.not.(MAXVAL(ABS(STIFF(this%elementMatrizen)))<=0.1d-300))THEN 
                     !    CALL printArray(&
                     !        & STIFF(this%elementMatrizen), 'STIFF')
                     !ENDIF
                     !IF(.not.(MAXVAL(ABS(FORCE(this%elementMatrizen)))<=0.1d-300))THEN 
                     !    CALL printArray(&
                     !        & FORCE(this%elementMatrizen), 'FORCE')
                     !ENDIF
                   
                    CALL schuettleMatrizen(this%elementMatrizen)
                 
                   !> \todo allow static Simulationsablauf
                    IF ( .TRUE. ) THEN
                      CALL Default1stOrderTime( MASS(this%elementMatrizen), &
                          & STIFF(this%elementMatrizen),&
                          & FORCE(this%elementMatrizen), &
                          & UElement=this%element)
                      !CALL Default1stOrderTime( MASS(this%elementMatrizen), &
                      !    & STIFF(this%elementMatrizen),&
                      !    & FORCE(this%elementMatrizen), &
                      !    & UElement=this%element)
                    END IF
                      CALL DefaultUpdateEquations( STIFF(this%elementMatrizen), &
                          & FORCE(this%ElementMatrizen), UElement=this%element )
                ENDIF
            ELSE IF(this%istGebietselement) THEN 
               CALL berechne(this%ErikssonDGL,&
                   & this%vorherigeLsg,&
                   & this%istGebietselement)
                !IF(.not.(MAXVAL(ABS(MASS(this%elementMatrizen)))<=0.1d-300))THEN 
                !    CALL printArray(&
                !    & MASS(this%elementMatrizen), 'MASS')
                !ENDIF
                !IF(.not.(MAXVAL(ABS(STIFF(this%elementMatrizen)))<=0.1d-300))THEN 
                !CALL printArray(&
                !    & STIFF(this%elementMatrizen), 'STIFF')
                !ENDIF
                !IF(.not.(MAXVAL(ABS(FORCE(this%elementMatrizen)))<=0.1d-300))THEN 
                !CALL printArray(&
                !    & FORCE(this%elementMatrizen), 'FORCE')
                !ENDIF
                CALL schuettleMatrizen(this%elementMatrizen)
               !> \todo allow static Simulationsablauf
               IF ( .TRUE. ) THEN
                 !CALL DefaultUpdateMass(MASS(this%elementMatrizen))
                 CALL Default1stOrderTime( MASS(this%elementMatrizen), &
                     & STIFF(this%elementMatrizen),&
                     & FORCE(this%elementMatrizen), &
                     & UElement=this%element)
               END IF
                 CALL DefaultUpdateEquations( STIFF(this%elementMatrizen), &
                     & FORCE(this%ElementMatrizen), UElement=this%element )

            ENDIF
      END SUBROUTINE berechneElement
          !> Löse die definierte DGL für dieses Element
      SUBROUTINE berechneElementStaggered1(this)
           TYPE(BauteilElement_t), INTENT(inout) :: this

           WRITE(MESSAGE, *) 'Berechne Elementmatrizen von Element: ', this%elementNummer
           CALL INFO('ElementClass', MESSAGE, level=5)

           ! nullifiziere Elementmatrizen vor der Berechnung
           CALL nullifiziere(this%elementMatrizen)
           IF(.not.this%istGebietselement)THEN
               !IF(ActiveBoundaryElement(this%element))THEN
               IF(.true.)THEN
                   CALL berechne(this%ErikssonDGL,&
                       & this%vorherigeLsg,&
                       & this%istGebietselement)
                !IF(.not.(MAXVAL(MASS(this%elementMatrizen))<=0.1d-300))THEN 
                !    CALL printArray(&
                !        & MASS(this%elementMatrizen), 'MASS')
                !ENDIF
                !IF(.not.(MAXVAL(STIFF(this%elementMatrizen))<=0.1d-300))THEN 
                !    CALL printArray(&
                !        & STIFF(this%elementMatrizen), 'STIFF')
                !ENDIF
                !IF(.not.(MAXVAL(FORCE(this%elementMatrizen))<=0.1d-300))THEN 
                !    CALL printArray(&
                !        & FORCE(this%elementMatrizen), 'FORCE')
                !ENDIF
                   
                   !> \todo allow static Simulationsablauf
                    IF ( .TRUE. ) THEN
                      CALL Default1stOrderTime( MASS(this%elementMatrizen), &
                          & STIFF(this%elementMatrizen),&
                          & FORCE(this%elementMatrizen), &
                          & UElement=this%element)
                    END IF
                      CALL DefaultUpdateEquations( STIFF(this%elementMatrizen), &
                          & FORCE(this%ElementMatrizen), UElement=this%element )
                ENDIF
            ELSE IF(this%istGebietselement) THEN 
               CALL berechne(this%ErikssonDGL,&
                   & this%vorherigeLsg,&
                   & this%istGebietselement)
                !CALL printArray(this%Geometrie%B, 'B')
                !IF(.not.(MAXVAL(MASS(this%elementMatrizen))<=0.1d-300))THEN 
                !    CALL printArray(&
                !    & MASS(this%elementMatrizen), 'MASS')
                !ENDIF
                !IF(.not.(MAXVAL(STIFF(this%elementMatrizen))<=0.1d-300))THEN 
                !CALL printArray(&
                !    & STIFF(this%elementMatrizen), 'STIFF')
                !ENDIF
                !IF(.not.(MAXVAL(FORCE(this%elementMatrizen))<=0.1d-300))THEN 
                !CALL printArray(&
                !    & FORCE(this%elementMatrizen), 'FORCE')
                !ENDIF
               !> \todo allow static Simulationsablauf
               IF ( .TRUE. ) THEN
                 CALL Default1stOrderTime( MASS(this%elementMatrizen), &
                     & STIFF(this%elementMatrizen),&
                     & FORCE(this%elementMatrizen), &
                     & UElement=this%element)
               END IF
                 CALL DefaultUpdateEquations( STIFF(this%elementMatrizen), &
                     & FORCE(this%ElementMatrizen), UElement=this%element )

            ENDIF
      END SUBROUTINE berechneElementStaggered1
          !> Löse die definierte DGL für dieses Element
      SUBROUTINE berechneElementStaggered2(this)
           TYPE(BauteilElement_t), INTENT(inout) :: this

           WRITE(MESSAGE, *) 'Berechne Elementmatrizen von Element: ', this%elementNummer
           CALL INFO('ElementClass', MESSAGE, level=5)

           ! nullifiziere Elementmatrizen vor der Berechnung
           CALL nullifiziere(this%elementMatrizen)
           IF(.not.this%istGebietselement)THEN
               !IF(ActiveBoundaryElement(this%element))THEN
               IF(.true.)THEN
                   CALL berechne(this%ErikssonDGL,&
                       & this%vorherigeLsg,&
                       & this%istGebietselement)
                !IF(.not.(MAXVAL(MASS(this%elementMatrizen))<=0.1d-300))THEN 
                !    CALL printArray(&
                !        & MASS(this%elementMatrizen), 'MASS')
                !ENDIF
                !IF(.not.(MAXVAL(STIFF(this%elementMatrizen))<=0.1d-300))THEN 
                !    CALL printArray(&
                !        & STIFF(this%elementMatrizen), 'STIFF')
                !ENDIF
                !IF(.not.(MAXVAL(FORCE(this%elementMatrizen))<=0.1d-300))THEN 
                !    CALL printArray(&
                !        & FORCE(this%elementMatrizen), 'FORCE')
                !ENDIF
                   
                   !> \todo allow static Simulationsablauf
                    IF ( .TRUE. ) THEN
                      CALL Default1stOrderTime( MASS(this%elementMatrizen), &
                          & STIFF(this%elementMatrizen),&
                          & FORCE(this%elementMatrizen), &
                          & UElement=this%element)
                    END IF
                      CALL DefaultUpdateEquations( STIFF(this%elementMatrizen), &
                          & FORCE(this%ElementMatrizen), UElement=this%element )
                ENDIF
            ELSE IF(this%istGebietselement) THEN 
               CALL berechne(this%ErikssonDGL,&
                   & this%vorherigeLsg,&
                   & this%istGebietselement)
                !CALL printArray(this%Geometrie%B, 'B')
                !IF(.not.(MAXVAL(MASS(this%elementMatrizen))<=0.1d-300))THEN 
                !    CALL printArray(&
                !    & MASS(this%elementMatrizen), 'MASS')
                !ENDIF
                !IF(.not.(MAXVAL(STIFF(this%elementMatrizen))<=0.1d-300))THEN 
                !CALL printArray(&
                !    & STIFF(this%elementMatrizen), 'STIFF')
                !ENDIF
                !IF(.not.(MAXVAL(FORCE(this%elementMatrizen))<=0.1d-300))THEN 
                !CALL printArray(&
                !    & FORCE(this%elementMatrizen), 'FORCE')
                !ENDIF
               !> \todo allow static Simulationsablauf
               IF ( .TRUE. ) THEN
                 CALL Default1stOrderTime( MASS(this%elementMatrizen), &
                     & STIFF(this%elementMatrizen),&
                     & FORCE(this%elementMatrizen), &
                     & UElement=this%element)
               END IF
                 CALL DefaultUpdateEquations( STIFF(this%elementMatrizen), &
                     & FORCE(this%ElementMatrizen), UElement=this%element )

            ENDIF
      END SUBROUTINE berechneElementStaggered2
          ! DOCUMENTATION
          !> This routine is taken from Elmers Stress.f90 and modified
          !!
          !> 
    !------------------------------------------------------------------------------
      SUBROUTINE InputTensor( Tensor, IsScalar, Name, Material, n, NodeIndexes )
    !------------------------------------------------------------------------------
          REAL(KIND=dp) :: Tensor(:,:,:)
          INTEGER :: n, NodeIndexes(:)
          LOGICAL :: IsScalar
          CHARACTER(LEN=*) :: Name
          TYPE(ValueList_t), POINTER :: Material
    !------------------------------------------------------------------------------
          LOGICAL :: FirstTime = .TRUE., stat
          REAL(KIND=dp), POINTER :: Hwrk(:,:,:)

          INTEGER :: i,j

          SAVE FirstTime, Hwrk
    !------------------------------------------------------------------------------
          IF ( FirstTime ) THEN
             NULLIFY( Hwrk )
             FirstTime = .FALSE.
          END IF

          Tensor = 0.0d0
          IsScalar = .TRUE.

          CALL ListGetRealArray( Material, Name, Hwrk, n, NodeIndexes, stat )
          IF ( .NOT. stat ) RETURN

          IsScalar = SIZE(HWrk,1) == 1 .AND. SIZE(HWrk,2) == 1
          IF ( IsScalar ) THEN
            DO i=1,SIZE(Tensor,1)
              Tensor(i,i,1:n) = Hwrk(1,1,1:n)
            END DO
          ELSE
            IF ( SIZE(Hwrk,1) == 1 ) THEN
               DO i=1,MIN(6,SIZE(HWrk,2) )
                  Tensor( i,i,1:n ) = Hwrk( 1,1,1:n )
               END DO
            ELSE IF ( SIZE(Hwrk,2) == 1 ) THEN
               DO i=1,MIN(6,SIZE(Hwrk,1))
                  Tensor( i,i,1:n ) = Hwrk( i,1,1:n )
               END DO
            ELSE
              DO i=1,MIN(6,SIZE(Hwrk,1))
                 DO j=1,MIN(6,SIZE(Hwrk,2))
                    Tensor( i,j,1:n ) = Hwrk( i,j,1:n )
                 END DO
              END DO
            END IF
          END IF

    !------------------------------------------------------------------------------
       END SUBROUTINE InputTensor
    !------------------------------------------------------------------------------
 End MODULE 
 !!\}
