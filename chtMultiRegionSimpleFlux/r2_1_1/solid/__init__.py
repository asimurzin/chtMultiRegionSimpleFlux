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


#-----------------------------------------------------------------
from Foam import ref, man


#-----------------------------------------------------------------
def createSolidMeshes( rp, runTime ):
    
    solidRegions = list()
    for index in range( rp.solidRegionNames().size() ):
        ref.ext_Info()<< "Create solid mesh for region " << rp.solidRegionNames()[ index ] \
            << " for time = " << runTime.timeName() << ref.nl << ref.nl
        
        solidRegions.append( man.fvMesh( man.IOobject ( rp.solidRegionNames()[ index ],
                                                        ref.fileName( runTime.timeName() ),
                                                        runTime,
                                                        ref.IOobject.MUST_READ ) ) )
        pass

    return solidRegions


#---------------------------------------------------------------------
def createSolidField( solidRegions, runTime ):

    thermos = list()
    for index in range( solidRegions.__len__() ):
        ref.ext_Info() << "*** Reading solid mesh thermophysical properties for region " \
                       << solidRegions[ index ].name() << ref.nl << ref.nl

        ref.ext_Info() << "    Adding to thermos\n" << ref.nl
        thermos.append( man.basicSolidThermo.New( solidRegions[ index ] ) )
        pass
   
    return thermos


#-----------------------------------------------------------------------------------------------------
def readSolidMultiRegionSIMPLEControls( mesh ):
    simple = mesh.solutionDict().subDict( ref.word( "SIMPLE" ) )
    nNonOrthCorr = simple.lookupOrDefault( ref.word( "nNonOrthogonalCorrectors" ), 0 )
    
    return simple, nNonOrthCorr


#-------------------------------------------------------------------------------------------------------
def setRegionSolidFields( i, solidRegions, thermos ):
    mesh = solidRegions[ i ]
    thermo = thermos[ i ]

    rho = man.volScalarField( thermo.rho(), man.Deps( thermo ) )

    cp = man.volScalarField( thermo.Cp(), man.Deps( thermo ) )
    
    kappa = man.volScalarField( thermo.ext_K(), man.Deps( thermo ) )
    # tmp<volSymmTensorField> tK = thermo.directionalK();
    
    # const volSymmTensorField& K = tK();
    T = man.volScalarField( thermo.T(), man.Deps( thermo ) )
    
    return mesh, thermo, rho, cp, kappa, T


#-------------------------------------------------------------------------------------------------------
def solveSolid( mesh, kappa, T, nNonOrthCorr ):
    for index in range( nNonOrthCorr + 1 ):
       TEqn = - ref.fvm.laplacian( kappa, T )
       TEqn.relax()
       TEqn.solve()
       pass
  
    ref.ext_Info()<< "Min/max T:" << T.ext_min().value() << ' ' << T.ext_max().value() << ref.nl
    
    pass
