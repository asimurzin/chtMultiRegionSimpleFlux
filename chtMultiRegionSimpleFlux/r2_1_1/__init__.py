#!/usr/bin/env python

#--------------------------------------------------------------------------------------
## pythonFlu - Python wrapping for OpenFOAM C++ API
## Copyright (C) 2010- Alexey Petrov
## Copyright (C) 2009-2010 Pebble Bed Modular Reactor (Pty) Limited (PBMR)
## 
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
## 
## See http://sourceforge.net/projects/pythonflu
##
## Author : Alexey PETROV, Andrey SIMURZIN
##


#------------------------------------------------------------------------------------
from Foam import ref, man


#------------------------------------------------------------------------------------
def main_standalone( argc, argv ):

    args = ref.setRootCase( argc, argv )

    runTime = man.createTime( args )
    
    rp = ref.compressible.regionProperties( runTime )
    
    from fluid import createFluidMeshes
    fluidRegions = createFluidMeshes( rp, runTime )
    
    from solid import createSolidMeshes,createSolidField
    solidRegions = createSolidMeshes( rp,runTime )
    
    from fluid import createFluidFields
    thermoFluid, rhoFluid, kappaFluid, UFluid, phiFluid, gFluid, turbulence, initialMassFluid, ghFluid, ghfFluid, \
                 p_rghFluid, radiation, pRefCellFluid, pRefValueFluid, rhoMax, rhoMin = createFluidFields( fluidRegions, runTime )

    from solid import createSolidField
    thermos = createSolidField( solidRegions, runTime )
    
    from fluid import initContinuityErrs
    cumulativeContErr = initContinuityErrs( fluidRegions )
    
    
    while runTime.loop() :
        ref.ext_Info()<< "Time = " << runTime.timeName() << ref.nl << ref.nl
                
        for i in range( fluidRegions.__len__() ):
            ref.ext_Info() << "\nSolving for fluid region " << fluidRegions[ i ].name() << ref.nl

            from fluid import setRegionFluidFields
            mesh, thermo, rho, kappa, U, phi, turb, p, psi, h, initialMass, p_rgh, gh, ghf, rad, pRefCell, pRefValue = \
                   setRegionFluidFields( i, fluidRegions, thermoFluid, rhoFluid, kappaFluid, UFluid, \
                                        phiFluid, turbulence, initialMassFluid, ghFluid, ghfFluid, p_rghFluid, radiation, pRefCellFluid, pRefValueFluid )
                
            from fluid import readFluidMultiRegionSIMPLEControls
            simple, nNonOrthCorr, momentumPredictor, transonic = readFluidMultiRegionSIMPLEControls( mesh ) 
                
            from fluid import solveFluid
            cumulativeContErr[ i ] = solveFluid( runTime, i, mesh, thermo, rad, thermoFluid, rho, kappa, U, phi, h, turb, p, psi, initialMass, p_rgh, gh, ghf, \
                                           simple, transonic, nNonOrthCorr, momentumPredictor, cumulativeContErr[ i ], rhoMax, rhoMin, pRefCell, pRefValue )
                
            pass
                
        for i in range( solidRegions.__len__() ):
            ref.ext_Info() << "\nSolving for solid region " << solidRegions[ i ].name() << ref.nl
               
            from solid import setRegionSolidFields
            mesh, thermo, rho, cp, kappa, T = setRegionSolidFields( i, solidRegions, thermos )

            from solid import readSolidMultiRegionSIMPLEControls
            simple, nNonOrthCorr = readSolidMultiRegionSIMPLEControls( mesh )
               
            from solid import solveSolid
            solveSolid( mesh, kappa, T, nNonOrthCorr )
            pass                
        runTime.write()

        ref.ext_Info()<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s" \
                      << "  ClockTime = " << runTime.elapsedClockTime() << " s" \
                      << ref.nl << ref.nl    

    ref.ext_Info() << "End\n"
    
    import os
    return os.EX_OK

    
#--------------------------------------------------------------------------------------
argv = None
import sys, os
from Foam import FOAM_VERSION
if FOAM_VERSION( ">=", "020101" ):
    if __name__ == "__main__" :
        argv = sys.argv
        os._exit( main_standalone( len( argv ), argv ) )
        pass
else:
    ref.ext_Info() << "\n\n To use this solver, it is necessary to SWIG OpenFOAM-2.1.1 or higher \n"    
    pass


#--------------------------------------------------------------------------------------

