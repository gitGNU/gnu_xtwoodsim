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
    !> \addtogroup BauteilClass
    !> \ingroup BauteilClass 
    !! \{
    !> \brief BauteilClass speichert die vorhergehende Lösung und verwaltet den
    !! Simulationsablauf
    !!
    !! Da die Simulation nach einem sogenannten Staggerd Iteration Scheme 
    !! ablaufen soll muss jeweils der vorhergehende Lösungsschritt 
    !! gespeichert werden. Weiterhin verwaltet BauteilClass seine Elemente
    !! vom Typ ElementeClass
MODULE BauteilClass
      USE DefUtils 
      USE LokaleGekoppelteMatrizenClass
      USE ElementClass
      implicit none


      LOGICAL, PRIVATE :: found
          ! DOCUMENTATION
          !> Array bestehend aus Pointern auf BauteilElemente::BauteilElement_t
          !!
          !> Diese Form der Linked List erlaubt es Pointer auf derived TYPE
          !! Definitionen in einem Array abzuspeichern, welches relativ
          !! klein ist, da es nur die Pointer und nicht die BauteilElement 
          !! selbst enthält
          !> \sa Ellies1994 Seite 573-575
          !> \param PointerAufElement verweist auf ein BauteilElement_t
          !! \param PointerAufNaechstesElement verweist auf das nächste 
          !! BauteilElement_t um schnell über alle Bauteile iterieren zu 
          !! können
      TYPE, PRIVATE :: ElementePointer
          TYPE(BauteilElement_t), POINTER:: PointerAufElement
          TYPE(BauteilElement_t), POINTER:: PointerAufNaechstesElement
      END TYPE
          ! DOCUMENTATION
          !> \struct BauteilClass::Bauteil_t
          !> bietet die Voraussetzung für das Staggerd Iteration Scheme
          !!
          !! Type Struktur zur Sammlung aller wichtiger Bauteileigeschaften
          !> indem der Lösungsvector für die Skalarfelder
          !! zwischengespeichert wird
      TYPE, Public :: Bauteil_t 
          PRIVATE
          REAL(KIND=dp), POINTER:: vorhergehendeLoesung
          REAL(KIND=dp) :: NonlinearTol=0.0d0
          INTEGER:: anzahlGebietselemente=0, anzahlRandgebietselemente=0,&
              & problemDim=1, NonlinearIter=0
          Type(ElementePointer), Allocatable, DIMENSION(:) :: ElementePointerListe
          TYPE(Solver_t) :: Solver
      END TYPE

      PRIVATE:: initialisiereBauteil, &
          & testLokaleMatrizenClass, testLokaleGekoppelteMatrizenClass, &
          & fuegeElementZuPointerlisteHinzu, berechneBauteil
          
      INTERFACE neu
          MODULE PROCEDURE initialisiereBauteil
      END INTERFACE
      INTERFACE berechne 
          MODULE PROCEDURE berechneBauteil
      END INTERFACE
  CONTAINS
          ! DOCUMENTATION
          !> Füllen der Werte des Bauteil_t
          !! 
          !> Initialisiert ein Bauteil für die Simulation.
          !! Dazu werden die für die global wichtigen Parameter aus 
          !! dem Solver_t und Model_t gezogen, welche Elmer bereitstellt.
          !! Weiterhin werden die Elemente (Gebiets- und Randgebietselemente)
          !! initialisiert
      SUBROUTINE initialisiereBauteil(this, Solver, Model)
          TYPE(Bauteil_t), INTENT(inout) :: this
          Type(Solver_t), INTENT(in) :: Solver
          TYPE(Model_t), INTENT(in) :: Model

          INTEGER:: AnzahlElemente=0, i, j
          TYPE(BauteilElement_t), POINTER :: PointerAufNaechstesElement
          TYPE(ValueList_t), POINTER :: SolverParams

          NULLIFY(PointerAufNaechstesElement)

          this%ProblemDim = CoordinateSystemDimension()
          this%anzahlGebietselemente = GetNOFActive() !Solver % NumberOfActiveElements
          this%anzahlRandgebietselemente = GetNOFBoundaryElements() !Solver % Mesh % NumberOfBoundaryElements
          this%Solver = Solver
          SolverParams => GetSolverParams(this%Solver)
          this%NonlinearIter = GetInteger(SolverParams, 'Nonlinear System Max&
              & Iterations', found)
          IF(.not.found) THEN
              CALL WARN('BauteilClass:initialisiereBauteil', 'Nonlinear System &
                  & Max Iterations nicht im sif file gefunden. Setze Wert auf &
                  & 1')
                  this%NonlinearIter = 1
          ENDIF
          this%NonlinearTol = GetConstReal(SolverParams, 'Nonlinear System&
              & Convergence Tolerance', found)
          IF(.not.found) THEN
              CALL WARN('BauteilClass:initialisiereBauteil', 'Nonlinear &
                  & System Convergence Tolerance not found in sif. Setze &
                  & Wert auf 1.0d-3')
              this%NonlinearTol = 1.0d-03
          ENDIF
          
          AnzahlElemente = this%anzahlGebietselemente + &
              & this%anzahlRandgebietselemente
              ! ElementClass braucht die Anzahl der Elemente 
              ! als globale Klassenvariable
          CALL setzeElementzaehler(AnzahlElemente)

          IF(ALLOCATED(this%ElementePointerListe))THEN
              IF(.not.(SIZE(this%ElementePointerListe).eq.AnzahlElemente))THEN
                  DEALLOCATE(this%ElementePointerListe)
                  ALLOCATE(this%ElementePointerListe(AnzahlElemente))
              ENDIF
          ELSE
              ALLOCATE(this%ElementePointerListe(AnzahlElemente))
          ENDIF

              ! NULLIFY sicherheitshalber die Pointer der Liste,
              ! damit sie nicht auf Nonsense verweisen.
          DO i=1, AnzahlElemente
              NULLIFY(this%ElementePointerListe(i)%PointerAufElement)
              NULLIFY(this%ElementePointerListe(i)%PointerAufNaechstesElement)
          ENDDO
              ! da zwei verschiedene elementtypen Initialisiert
              ! werden sollen und ich es nicht für notwendig und
              ! vorteilhaft halte in jedem SChleifendurchlauf 
              ! einen weiteren Boolean Vergleich zu machen werden die 
              ! Elemente in zwei Schleifen intialisiert. 
              ! Bringt das einen Geschwindigkeitsvorteil?
          DO i=1, this%anzahlGebietselemente
              !print *, 'intialisiert Element: ', i
              CALL fuegeElementZuPointerlisteHinzu(this, i,i, .true., Solver, Model)
          ENDDO
          DO i=1, this%anzahlRandgebietselemente
              j = this%anzahlGebietselemente+i
              !print *, 'Initialisiert Randelement: ', j
              CALL fuegeElementZuPointerlisteHinzu(this,j, i,.false., Solver, Model)
          ENDDO
              ! Das Erste Element muss noch den Pointer auf das Letzte bekommen,
              ! um den Kreis zu schließen und alle Elemente mit Pointern auf
              ! ein Nachbarelement zu versehen
          CALL setzePointerAufNaechstesElement( &
              & this%ElementePointerListe(1)%PointerAufElement, &
              &  this%ElementePointerListe(AnzahlElemente)%PointerAufElement)
      END SUBROUTINE initialisiereBauteil

          !> Reserviere Speicherplatz und erzeuge ein neues BauteilElemente-
          !! Objekt
          !!
          !> Speicherplatz für ein BauteilElement_t wird reserviert, 
          !! initialisiert ein neues Objekt vom Typ BauteilElement_t, fügt
          !! dieses in das Array ElementePointerListe hinzu und
          !! setzt bei diesem den Pointer zum virtuellen Nachbarelement. 
      SUBROUTINE fuegeElementZuPointerlisteHinzu(this, listenpos, ElemNr, &
          & istGebietselement, Solver, Model)
          TYPE(Bauteil_t), INTENT(inout) :: this
          INTEGER, INTENT(in) :: ElemNr, listenpos
          LOGICAL :: istGebietselement
          TYPE(Solver_t), INTENT(in) :: Solver
          TYPE(Model_t), INTENT(in) :: Model

          LOGICAL, SAVE :: ersterAufruf=.true.
          
          TYPE(BauteilElement_t), POINTER, SAVE :: PointerAufNaechstesElement
          IF(ersterAufruf) THEN 
              NULLIFY(PointerAufNaechstesElement)
              PointerAufNaechstesElement => NULL()
          ENDIF
          ersterAufruf = .false.

          ALLOCATE(this%ElementePointerListe(listenpos)%PointerAufElement)
          ! intialisiere dieses neue Element
          CALL neu(this%ElementePointerListe(listenpos)%PointerAufElement, &
              &ElemNr, Solver, Model, istGebietselement)
          ! und setze den Pointer auf das nächste Element
          CALL setzePointerAufNaechstesElement( &
              & this%ElementePointerListe(listenpos)%PointerAufElement, &
              & PointerAufNaechstesElement )
          ! speichere den aktuellen Pointer für das nächste
          ! Element 
          PointerAufNaechstesElement => &
              & this%ElementePointerListe(listenpos)%PointerAufElement
      END SUBROUTINE fuegeElementZuPointerlisteHinzu

          ! DOCUMENTATION
          !> intialisiert und steuert die Eriksson2006Solve berechnung
          !! für dieses Mesh 
          !! 
          !> \todo ausbauen
      SUBROUTINE berechneBauteil(this) 
          TYPE(Bauteil_t), INTENT(inout):: this
          INTEGER :: iter
          LOGICAL :: converged
          REAL(KIND=dp), SAVE :: Norm, PrevNorm=0, RelativChange
          DO iter=1, this%NonlinearIter
              converged = .false.
              WRITE(MESSAGE,'(A,I5,A,I5)') 'Nonlinear iteration no.', &
                  & iter, 'of max. ', this%NonlinearIter
              CALL INFO('Eriksson2006', MESSAGE, level=2)
              CALL DefaultInitialize(this%Solver)

              CALL iteriereUeberAlleElemente( &
                  & this%ElementePointerListe(1)%PointerAufElement)

              CALL DefaultFinishAssembly()
              CALL DefaultDirichletBCs()

              Norm=DefaultSolve()
              !IF(PrevNorm + Norm /=0.0d0) THEN
              !    RelativChange = 2.0d0 + ABS(PrevNorm-Norm)/(PrevNorm+Norm)
              !ELSE
              !    RelativChange = 0.0d0
              !ENDIF
              RelativChange = this%Solver%Variable%NonlinChange
              WRITE(MESSAGE, *) 'Result Norm :', Norm
              CALL INFO('Eriksson2006', MESSAGE, level=2)
              WRITE(MESSAGE, *) 'Relativ Change :', RelativChange
              CALL INFO('Eriksson2006', MESSAGE, level=2)
              IF(RelativChange < this%NonlinearTol) THEN
                  converged =.TRUE.
                  EXIT
              ELSE
                  PrevNorm = Norm
              ENDIF
          ENDDO
      END SUBROUTINE berechneBauteil
          ! DOCUMENTATION
          !> intialisiert und steuert die Eriksson2006Solve berechnung
          !! für dieses Mesh 
          !! 
          !> \todo ausbauen
      SUBROUTINE berechneBauteilStaggered(this) 
          TYPE(Bauteil_t), INTENT(inout):: this
          INTEGER :: iter
          LOGICAL :: converged
          REAL(KIND=dp), SAVE :: Norm, PrevNorm=0, RelativChange
          DO iter=1, this%NonlinearIter
              converged = .false.
              WRITE(MESSAGE,'(A,I5,A,I5)') 'Nonlinear iteration no.', &
                  & iter, 'of max. ', this%NonlinearIter
              CALL INFO('Eriksson2006', MESSAGE, level=2)
              CALL DefaultInitialize(this%Solver)

              CALL iteriereUeberAlleElementeStaggered1( &
                  & this%ElementePointerListe(1)%PointerAufElement)

              CALL DefaultFinishAssembly()
              CALL DefaultDirichletBCs()

              Norm=DefaultSolve()
              !IF(PrevNorm + Norm /=0.0d0) THEN
              !    RelativChange = 2.0d0 + ABS(PrevNorm-Norm)/(PrevNorm+Norm)
              !ELSE
              !    RelativChange = 0.0d0
              !ENDIF
              RelativChange = this%Solver%Variable%NonlinChange
              WRITE(MESSAGE, *) 'Result Norm :', Norm
              CALL INFO('Eriksson2006', MESSAGE, level=2)
              WRITE(MESSAGE, *) 'Relativ Change :', RelativChange
              CALL INFO('Eriksson2006', MESSAGE, level=2)
              IF(RelativChange < this%NonlinearTol) THEN
                  converged =.TRUE.
                  EXIT
              ELSE
                  PrevNorm = Norm
              ENDIF
          ENDDO
          DO iter=1, this%NonlinearIter
              converged = .false.
              WRITE(MESSAGE,'(A,I5,A,I5)') 'Nonlinear iteration no.', &
                  & iter, 'of max. ', this%NonlinearIter
              CALL INFO('Eriksson2006', MESSAGE, level=2)
              CALL DefaultInitialize(this%Solver)

              CALL iteriereUeberAlleElementeStaggered2( &
                  & this%ElementePointerListe(1)%PointerAufElement)

              CALL DefaultFinishAssembly()
              CALL DefaultDirichletBCs()

              Norm=DefaultSolve()
              !IF(PrevNorm + Norm /=0.0d0) THEN
              !    RelativChange = 2.0d0 + ABS(PrevNorm-Norm)/(PrevNorm+Norm)
              !ELSE
              !    RelativChange = 0.0d0
              !ENDIF
              RelativChange = this%Solver%Variable%NonlinChange
              WRITE(MESSAGE, *) 'Result Norm :', Norm
              CALL INFO('Eriksson2006', MESSAGE, level=2)
              WRITE(MESSAGE, *) 'Relativ Change :', RelativChange
              CALL INFO('Eriksson2006', MESSAGE, level=2)
              IF(RelativChange < this%NonlinearTol) THEN
                  converged =.TRUE.
                  EXIT
              ELSE
                  PrevNorm = Norm
              ENDIF
          ENDDO
          DO iter=1, this%NonlinearIter
              converged = .false.
              WRITE(MESSAGE,'(A,I5,A,I5)') 'Nonlinear iteration no.', &
                  & iter, 'of max. ', this%NonlinearIter
              CALL INFO('Eriksson2006', MESSAGE, level=2)
              CALL DefaultInitialize(this%Solver)

              CALL iteriereUeberAlleElementeStaggered1( &
                  & this%ElementePointerListe(1)%PointerAufElement)

              CALL DefaultFinishAssembly()
              CALL DefaultDirichletBCs()

              Norm=DefaultSolve()
              !IF(PrevNorm + Norm /=0.0d0) THEN
              !    RelativChange = 2.0d0 + ABS(PrevNorm-Norm)/(PrevNorm+Norm)
              !ELSE
              !    RelativChange = 0.0d0
              !ENDIF
              RelativChange = this%Solver%Variable%NonlinChange
              WRITE(MESSAGE, *) 'Result Norm :', Norm
              CALL INFO('Eriksson2006', MESSAGE, level=2)
              WRITE(MESSAGE, *) 'Relativ Change :', RelativChange
              CALL INFO('Eriksson2006', MESSAGE, level=2)
              IF(RelativChange < this%NonlinearTol) THEN
                  converged =.TRUE.
                  EXIT
              ELSE
                  PrevNorm = Norm
              ENDIF
          ENDDO
      END SUBROUTINE berechneBauteilStaggered
      
          ! ---- SUBROUTINE berechneRelativeFeuchte 
          !> berechne relative Feuchtigkeit aus dem Feuchtegehalt
          !! und der Temperatur
          !! 
          !> Die DGL gib den Feuchtegehalt als Lösung
          !! dieser muss für Quellung noch in die
          !! relative Feuchtigkeit umgerechnet werden.
          !! Das geschieht mithilfe von sogenannten
          !! Sorptionsisothermen. 
          !> \todo 
          !! - Formeln für relative Feuchtigkeit einbetten und 
          !!   Lösungen im Postprozessor als Variable verfügbar 
          !!   machen
      SUBROUTINE berechneRelativeFeuchte(this)
          TYPE(Bauteil_t), INTENT(inout) :: this
          !this%vorhergehendeLoesung = 0
      END SUBROUTINE berechneRelativeFeuchte

      ! ############# TESTROUTINEN #################
      SUBROUTINE testLokaleMatrizenClass(this)
          TYPE( Bauteil_t ), INTENT(inout) :: this

          !> deklariere ein Objekt LokaleMatrizen der Klasse
          !! LokaleMatrizen_t
          TYPE( LokaleMatrizen_t ) :: LokaleMatrizen, newMa

          !> deklarieren der STIFF, MASS und FORCE Matrizen
          !! für Testzwecke
          REAL(KIND=dp), DIMENSION(4,4) :: aSTIFF, aMASS,a
          REAL(KIND=dp), DIMENSION(4)   :: aFORCE,b

          CALL neu(LokaleMatrizen, aFORCE, aSTIFF, aMASS)
          b = 3.5d0
          CALL setze(LokaleMatrizen,wertFORCE=b)
          
          a=8.5d0
          CALL setze(LokaleMatrizen,wertMASS=a)

          CALL neu(newMa, aFORCE,aSTIFF,aMASS)
          newMa = LokaleMatrizen + LokaleMatrizen !< teste + operator auf LokaleMatrizen_t angewendet
      END SUBROUTINE testLokaleMatrizenClass
          ! DOCUMENTATION
          !> Testroutine, um LokaleGekoppelteMatrizenClass zu testen
      SUBROUTINE testLokaleGekoppelteMatrizenClass(this)
          TYPE( Bauteil_t ), INTENT(inout) :: this
          
          !> deklariere ein Objekt LokaleGekoppelteMatrizen der Klasse
          !! LokaleGekoppelteMatrizenClass
          TYPE(LokaleGekoppelteMatrizen_t) :: M1,M2
          !> deklarieren der STIFF, MASS und FORCE Matrizen
          !! für Testzwecke
          REAL(KIND=dp), DIMENSION(4,4) :: STIFFa, MASSa, a
          REAL(KIND=dp), DIMENSION(4)   :: FORCEa, b

          b = 3.5d0
          CALL neu(M1, FORCEa, STIFFa, MASSa)
          CALL setze(M1,wertFORCE=b)
          a = 8.5d0
          CALL setze(M1,wertMASS=a)

          CALL neu(M2, FORCEa,STIFFa,MASSa)
          M2 = M1 + M1 !< teste + operator auf LokaleMatrizen_t angewendet
          CALL setze(M1, wertFORCE=b)
      END SUBROUTINE testLokaleGekoppelteMatrizenClass
END MODULE
!> \}

