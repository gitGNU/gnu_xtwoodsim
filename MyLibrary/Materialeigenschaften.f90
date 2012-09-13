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
!> struct enhält die für diese Simulation nötigen materialWerte
!!
!> Alle hier gelisteten Materialwerte werden momentan in ElementClass
!! gesetzt und durch Eriksson2006DGL_Class in der Berechnung verwendet
!!
!> \todo mache MaterialSammlung_Class eigenständig durch eigene
!! Funktionen und als linked List unabhängig von ElementClass, 
!! Eriksson2006DGL_Class, etc
MODULE MaterialSammlung_Class
 USE DefUtils
 IMPLICIT NONE
  TYPE :: materialSammlung_t
      REAL(Kind=dp), DIMENSION(:), POINTER :: dichte=>NULL()
      REAL(Kind=dp), DIMENSION(:), POINTER :: trockendichte=>NULL()
      REAL(Kind=dp), DIMENSION(:), POINTER :: diffusionskoeffizient=>NULL()
      REAL(Kind=dp), DIMENSION(:), POINTER :: Konduktivitaet=>NULL()
      REAL(Kind=dp), DIMENSION(:), POINTER :: spezifischeWaerme=>NULL()

      REAL(Kind=dp), DIMENSION(:), POINTER :: E_b=>NULL() !< Activation energy
          !! of bound water
      REAL(Kind=dp), DIMENSION(:), POINTER :: R=>NULL() !< gas konstant
      REAL(Kind=dp), DIMENSION(:), POINTER :: relativeHumidity=>NULL()
      REAL(Kind=dp), DIMENSION(:), POINTER :: dwdH=>NULL()

      REAL(Kind=dp), DIMENSION(:,:,:), POINTER :: dichteTensor=>NULL()
      REAL(Kind=dp), DIMENSION(:,:,:), POINTER :: trockendichteTensor=>NULL()
      REAL(Kind=dp), DIMENSION(:,:,:), POINTER :: diffusionskoeffizientTensor=>NULL()
      REAL(Kind=dp), DIMENSION(:,:,:), POINTER :: KonduktivitaetTensor=>NULL()
      REAL(Kind=dp), DIMENSION(:,:,:), POINTER :: spezifischeWaermeTensor=>NULL()
  END TYPE materialSammlung_t
END MODULE MaterialSammlung_Class

    !> struct dieses Moduls enthält die Werte für die geometrie 
    !! eines Elements
    !!
    !> Alle hier gelisteten Werte werden momentan in ElementClass
    !! gesetzt und durch Eriksson2006DGL_Class in der Berechnung verwendet
    !!
    !> \todo mache geometrischeEigenschaften_Class eigenständig durch eigene
    !! Funktionen und als linked List unabhängig von ElementClass, 
    !! Eriksson2006DGL_Class, etc
MODULE geometrischeEigenschaften_Class
USE DefUtils
IMPLICIT NONE
    !> geometrischeEigenschaften_t 
    !! Sammlung der für die Simulation nötigen geometrischen 
    !! Elementeigenschaften
  TYPE :: geometrischeEigenschaften_t
      REAL(Kind=dp) :: detJ
      !> Achtung beim setzen der DIMENSION für N und B
      REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: N
      REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: B
      REAL(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: dBdx
      TYPE(GaussIntegrationPoints_t) :: IP
      INTEGER:: anzahlKnoten
      INTEGER:: Knotenfreiwerte
      TYPE(Nodes_t) :: Nodes
      INTEGER :: DIM
      TYPE(Element_t), POINTER :: Element
  END TYPE
END MODULE geometrischeEigenschaften_Class
