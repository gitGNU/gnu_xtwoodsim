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
    ! ---- Eriksson2006Solve_Init
    !> \defgroup Eriksson2006Solve
    !! \{
    !> \ingroup Eriksson2006Solve 
    !! Eriksson2006Solve.f90 beinhaltet die SUBROUTINE Eriksson2006Solve
    !! welche durch den Elmer Kern in jeder Steady state Iteration 
    !! aufgerufen wird. 

    !> @brief Initialisiere einige Standardparameter für die Simulation

    !> Wird automatisch beim Aufruf der SUBROUTINE Eriksson2006Solve
    !!aufgerufen.
    !> @param Model bereitgestellt von Elmer
    !> @param Solver bereitgestellt von Elmer
    !> @param dt Schrittweite, bereitgestellt von Elmer
    !> @param Transient bereitgestellt von Elmer
SUBROUTINE Eriksson2006Solve_Init( Model, Solver, dt, Transient )
      USE DefUtils

      TYPE(Model_t) :: Model
      TYPE(Solver_t) :: Solver
      REAL(KIND=dp) :: dt
      LOGICAL :: Transient

      INTEGER :: dim
      TYPE(ValueList_t), POINTER :: SolverParams


      SolverParams => GetSolverParams()

     IF ( .NOT. ListCheckPresent( SolverParams,'Variable') ) THEN
      dim = CoordinateSystemDimension()
      !CALL ListAddInteger( SolverParams, 'Variable DOFs', dim )
      !CALL ListAddString( SolverParams, 'Variable', 'Displacement' )
     END IF
     !CALL ListAddInteger( SolverParams, 'Time derivative order', 2 )
End SUBROUTINE Eriksson2006Solve_Init

    ! ---- Eriksson2006Solve
    !> @brief Eriksson2006Solve organisiert den Wärme-Feuchtigkeitstransportsolver

    !> Aufruf erfolgt durch Elmer 
    !! Der Solver übernimmt die Assemblierung der Elementmatrizen für 
    !! Rand- und Gebietselemente, berücksichtigt dabei die Randbedingungen
    !! für Dirichlet-, Neumann- und Cauchyrandbedingungen und prüft die 
    !! Fehlerschranken für die nichtlineare Iterationsschleife, nachdem 
    !! das Gleichungssystem durch, von Elmer bereitgestellte, Routinen
    !! gelöst wurde
    !! \sa http://www.elmerfem.org/doxygen
    !> @param Model bereitgestellt von Elmer
    !> @param Solver bereitgestellt von Elmer
    !> @param dt Schrittweite, bereitgestellt von Elmer
    !> @param Transient bereitgestellt von Elmer
    !> \todo Deactivierung oder spätere Aktivierung
    !! funktioniert bisher nicht! Alle Elemente müssen
    !! während der gesammten Simulationsdauer aktiv sein.
    !! Ursache sind Beschränkungen in ElementClass, in der 
    !! SUBROUTINE initialisiereBauteilElement, wo geprüft wird
    !! ob aktuell der erste Zeitschritt durchlaufen wird.
SUBROUTINE Eriksson2006Solve( Model, Solver, dt, Transient )
      USE DefUtils 
      USE BauteilClass

      TYPE(Model_t) :: Model
      TYPE(Solver_t) :: Solver
      REAL(KIND=dp) :: dt
      LOGICAL :: Transient

      !> deklariere ein neues Bauteil
      TYPE( Bauteil_t ) :: Bauteil
      CALL neu(Bauteil, Solver, Model)
      CALL berechne(Bauteil)
      !CALL berechneRelativeFeuchte(Bauteil)
      !CALL testLokaleMatrizenClass(Bauteil)
      !CALL testLokaleGekoppelteMatrizenClass(Bauteil)
End SUBROUTINE Eriksson2006Solve
!> \}
