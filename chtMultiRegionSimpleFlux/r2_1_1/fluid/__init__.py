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
def createFluidMeshes( rp, runTime ) :
    
    fluidRegions = list()
    for index in range( rp.fluidRegionNames().size() ) :
        ref.ext_Info()<< "Create fluid mesh for region " << rp.fluidRegionNames()[ index ] \
                  << " for time = " << runTime.timeName() << ref.nl << ref.nl
        mesh = man.fvMesh( man.IOobject( rp.fluidRegionNames()[ index ],
                                         ref.fileName( runTime.timeName() ),
                                         runTime,
                                         ref.IOobject.MUST_READ ) )
        fluidRegions.append( mesh )
        pass
    
    return fluidRegions
    
    
#-------------------------------------------------------------------
def createFluidFields( fluidRegions, runTime ) :
    
    # Initialise fluid field pointer lists
    thermoFluid = list() 
    rhoFluid = list()
    kappaFluid = list()
    UFluid = list()
    phiFluid = list()
    gFluid = list()
    turbulence =list()
    p_rghFluid = list()
    ghFluid = list()
    ghfFluid = list()
    radiation =list()
    initialMassFluid = list()
    pRefCellFluid = list()
    pRefValueFluid = list()

    rhoMax = list()
    rhoMin = list()
    
    #Populate fluid field pointer lists

    for index in range( fluidRegions.__len__() ) :
        ref.ext_Info() << "*** Reading fluid mesh thermophysical properties for region " \
            << fluidRegions[ index ].name() << ref.nl << ref.nl

        ref.ext_Info()<< "    Adding to thermoFluid\n" << ref.nl
        
        thermo = man.basicRhoThermo.New( fluidRegions[ index ] )
        thermoFluid.append( thermo )
        
        ref.ext_Info()<< "    Adding to rhoFluid\n" << ref.nl
        rhoFluid.append( man.volScalarField( man.IOobject( ref.word( "rho" ), 
                                                           ref.fileName( runTime.timeName() ), 
                                                           fluidRegions[ index ], 
                                                           ref.IOobject.NO_READ, 
                                                           ref.IOobject.AUTO_WRITE ),
                                              man.volScalarField( thermoFluid[ index ].rho(), man.Deps( thermoFluid[ index ] ) ) ) )
        
        ref.ext_Info()<< "    Adding to kappaFluid\n" << ref.nl
        kappaFluid.append( man.volScalarField( man.IOobject( ref.word( "kappa" ),
                                                             ref.fileName( runTime.timeName() ),
                                                             fluidRegions[ index ],
                                                             ref.IOobject.NO_READ,
                                                             ref.IOobject.NO_WRITE ),
                                               man.volScalarField( thermoFluid[ index ].Cp() * thermoFluid[ index ].alpha(), 
                                                                   man.Deps( thermoFluid[ index ] ) ) ) )
                                                       
        ref.ext_Info()<< "    Adding to UFluid\n" << ref.nl
        UFluid.append( man.volVectorField( man.IOobject( ref.word( "U" ),
                                                         ref.fileName( runTime.timeName() ),
                                                         fluidRegions[ index ],
                                                         ref.IOobject.MUST_READ,
                                                         ref.IOobject.AUTO_WRITE ),
                                           fluidRegions[ index ] ) )
        
        ref.ext_Info()<< "    Adding to phiFluid\n" << ref.nl
        phiFluid.append( man.surfaceScalarField( man.IOobject( ref.word( "phi" ),
                                                               ref.fileName( runTime.timeName() ),
                                                               fluidRegions[ index ],
                                                               ref.IOobject.READ_IF_PRESENT,
                                                               ref.IOobject.AUTO_WRITE),
                                                  man.linearInterpolate( rhoFluid[ index ] * UFluid[ index ] ) & 
                                                  man.surfaceVectorField( fluidRegions[ index ].Sf(), man.Deps( fluidRegions[ index ] ) ) ) )
        
        ref.ext_Info()<< "    Adding to gFluid\n" << ref.nl
        gFluid.append( man.uniformDimensionedVectorField( man.IOobject( ref.word( "g" ),
                                                                        ref.fileName( runTime.constant() ),
                                                                        fluidRegions[ index ],
                                                                        ref.IOobject.MUST_READ,
                                                                        ref.IOobject.NO_WRITE ) ) )        
        
        ref.ext_Info()<< "    Adding to turbulence\n" << ref.nl
        turbulence.append( man.compressible.turbulenceModel.New( rhoFluid[ index ],
                                                                 UFluid[ index ],
                                                                 phiFluid[ index ],
                                                                 thermoFluid[ index ] ) )
        ref.ext_Info() << "    Adding to ghFluid\n" << ref.nl
        ghFluid.append( man.volScalarField( ref.word( "gh" ) , 
                                            gFluid[ index ] & man.volVectorField( fluidRegions[ index ].C(), man.Deps( fluidRegions[ index ] ) ) ) )

        ref.ext_Info() << "    Adding to ghfFluid\n" << ref.nl
        ghfFluid.append( man.surfaceScalarField( ref.word( "ghf" ), 
                                                 gFluid[ index ] & man.surfaceVectorField( fluidRegions[ index ].Cf(), man.Deps( fluidRegions[ index ] ) ) ) )

        p_rghFluid.append( man.volScalarField( man.IOobject( ref.word( "p_rgh" ),
                                                             ref.fileName( runTime.timeName() ),
                                                             fluidRegions[ index ],
                                                             ref.IOobject.MUST_READ,
                                                             ref.IOobject.AUTO_WRITE ),
                                               fluidRegions[ index ] ) )
        # Force p_rgh to be consistent with p
        p_rghFluid[ index ] << thermoFluid[ index ].p()() - rhoFluid[ index ] * ghFluid[ index ]
        
        radiation.append( man.radiation.radiationModel.New( man.volScalarField( thermoFluid[ index ].T(), man.Deps( thermoFluid[ index ] ) ) ) )
        
        initialMassFluid.append( ref.fvc.domainIntegrate( rhoFluid[ index ] ).value()  )
        
        pRefCellFluid.append( 0 )
        pRefValueFluid.append( 0.0 )
        ref.setRefCell( thermoFluid[ index ].p(), p_rghFluid[ index ], 
                        fluidRegions[ index ].solutionDict().subDict( ref.word( "SIMPLE" ) ), 
                        pRefCellFluid[ index ], pRefValueFluid[ index ] )
        
        rhoMax.append( ref.dimensionedScalar( fluidRegions[ index ].solutionDict().subDict( ref.word( "SIMPLE" ) ).lookup( ref.word( "rhoMax" ) ) ) )
        rhoMin.append( ref.dimensionedScalar( fluidRegions[ index ].solutionDict().subDict( ref.word( "SIMPLE" ) ).lookup( ref.word( "rhoMin" ) ) ) )
        pass

    
    return thermoFluid, rhoFluid, kappaFluid, UFluid, phiFluid, gFluid, turbulence, \
           initialMassFluid, ghFluid, ghfFluid, p_rghFluid, radiation, pRefCellFluid, pRefValueFluid, rhoMax, rhoMin
        

#-----------------------------------------------------------------------------------------------------------------------
def readFluidMultiRegionSIMPLEControls( mesh ) :
    
    simple = mesh.solutionDict().subDict( ref.word( "SIMPLE" ) )
    nNonOrthCorr = simple.lookupOrDefault( ref.word( "nNonOrthogonalCorrectors" ), 0 )
    momentumPredictor = simple.lookupOrDefault( ref.word( "momentumPredictor" ), ref.Switch( True ) )
    transonic = simple.lookupOrDefault( ref.word( "transonic" ), ref.Switch( False ) )
    
    return simple, nNonOrthCorr, momentumPredictor, transonic


#--------------------------------------------------------------------------------------------------------------------------
def setRegionFluidFields( i, fluidRegions, thermoFluid, rhoFluid, kappaFluid, UFluid, phiFluid, turbulence, \
                          initialMassFluid, ghFluid, ghfFluid, p_rghFluid, radiation, pRefCellFluid, pRefValueFluid ):
    mesh = fluidRegions[ i ]

    thermo = thermoFluid[ i ]
    rho = rhoFluid[ i ]
    kappa = kappaFluid[ i ]
    U = UFluid[ i ]
    phi = phiFluid[ i ]
    
    turb = turbulence[ i ]

    p = thermo.p()
    psi = thermo.psi()
    h = thermo.h()
    
    p_rgh = p_rghFluid[ i ]
    gh = ghFluid[ i ]
    ghf = ghfFluid[ i ]
    
    rad = radiation[ i ]
    
    initialMass = ref.dimensionedScalar( ref.word( "initialMass" ), ref.dimMass, initialMassFluid[ i ] )
    
    pRefCell = pRefCellFluid[ i ]
    pRefValue = pRefValueFluid[ i ]

    return mesh, thermo, rho, kappa, U, phi, turb, p, psi, h, initialMass, p_rgh, gh, ghf, rad, pRefCell, pRefValue


#--------------------------------------------------------------------------------------------------------------------------
def initContinuityErrs( fluidRegions ):
    cumulativeContErr = list()
    for index in range( fluidRegions.__len__() ) :
        cumulativeContErr.append( 0.0 )
        pass
    
    return cumulativeContErr


#--------------------------------------------------------------------------------------------------------------------------
def fun_UEqn( rho, U, phi, ghf, p_rgh, turb, mesh ):
    # Solve the Momentum equation
    
    UEqn = man.fvm.div( phi, U ) + man.fvVectorMatrix( turb.divDevRhoReff( U ), man.Deps( turb, U ) )

    UEqn.relax()
    
    ref.solve( UEqn() == ref.fvc.reconstruct( ( -ghf * ref.fvc.snGrad( rho ) - ref.fvc.snGrad( p_rgh ) ) * mesh.magSf() ) )
    
    return UEqn


#--------------------------------------------------------------------------------------------------------------------------
def fun_hEqn( U, rho, h, phi, turb, thermo, rad, mesh ):
    
    hEqn = ( ( ref.fvm.div( phi, h ) -  ref.fvm.Sp( ref.fvc.div( phi ), h ) - ref.fvm.laplacian( turb.alphaEff(), h ) ) \
               == - ref.fvc.div( phi, 0.5 * U.magSqr(), ref.word( "div(phi,K)" ) ) + rad.Sh( thermo() )  )  # mixed calculation
    
    hEqn.relax()
    hEqn.solve() 
   
    thermo.correct()
    rad.correct()
    
    ref.ext_Info()<< "Min/max T:" << thermo.T().ext_min().value() << ' ' \
        << thermo.T().ext_max().value() << ref.nl
        
    pass


#---------------------------------------------------------------------------------------------------------    
def fun_pEqn( runTime, i, mesh, p, rho, turb, thermo, thermoFluid, kappa, UEqn, U, phi, psi, \
              initialMass, p_rgh, gh, ghf, nNonOrthCorr, cumulativeContErr, rhoMax, rhoMin, pRefCell, pRefValue ):
    
    rho << thermo.rho()
    rho << rho.ext_max( rhoMin[ i ] )
    rho << rho.ext_min( rhoMax[ i ] )
    rho.relax()

    rAU = 1.0 / UEqn.A()
    rhorAUf = ref.surfaceScalarField( ref.word( "(rho*(1|A(U)))" ), ref.fvc.interpolate( rho * rAU ) )

    U << rAU * UEqn.H()
    #UEqn.clear()

    phi << ref.fvc.interpolate( rho ) * ( ref.fvc.interpolate( U ) & mesh.Sf() )
    closedVolume = ref.adjustPhi( phi, U, p_rgh )
    compressibility = ref.fvc.domainIntegrate( psi )
    compressible = ( compressibility.value() > ref.SMALL)

    buoyancyPhi = rhorAUf * ghf * ref.fvc.snGrad( rho ) * mesh.magSf()
    phi -= buoyancyPhi

    # Solve pressure
    for nonOrth in range( nNonOrthCorr + 1 ):
        p_rghEqn = ref.fvm.laplacian( rhorAUf, p_rgh ) == ref.fvc.div( phi ) 

        if compressible:
            tmp = ref.getRefCellValue(p_rgh, pRefCell)
            pass
        else:
            tmp = pRefValue
            pass
        p_rghEqn.setReference( pRefCell, tmp )

        p_rghEqn.solve()

        if nonOrth == nNonOrthCorr:
            # Calculate the conservative fluxes
            phi -= p_rghEqn.flux()

            # Explicitly relax pressure for momentum corrector
            p_rgh.relax()

            # Correct the momentum source with the pressure gradient flux
            # calculated from the relaxed pressure
            U -= rAU * ref.fvc.reconstruct( ( buoyancyPhi + p_rghEqn.flux() ) / rhorAUf )
            U.correctBoundaryConditions()
            pass
        pass

    p << p_rgh + rho * gh
    
    cumulativeContErr = ref.ContinuityErrs( phi, runTime, mesh, cumulativeContErr )

    # For closed-volume cases adjust the pressure level
    # to obey overall mass continuity
    if closedVolume and compressible:
        p += ( initialMass - ref.fvc.domainIntegrate( thermo.rho() ) ) / compressibility
        p_rgh << p() - rho * gh
        pass

    rho << thermo.rho();
    rho << rho.ext_max( rhoMin[ i ] )
    rho << rho.ext_min( rhoMax[ i ] )
    rho.relax()

    ref.ext_Info() << "Min/max rho:" << rho.ext_min().value() << ' ' << rho.ext_max().value() << ref.nl

    # Update thermal conductivity
    kappa << thermo.Cp() * turb.alphaEff()
    
    return cumulativeContErr


#--------------------------------------------------------------------------------------------------------------------------
def solveFluid( runTime, i, mesh, thermo, rad, thermoFluid, rho, kappa, U, phi, h, turb, p, psi, initialMass, p_rgh, gh, ghf, \
                pimple, transonic, nNonOrthCorr, momentumPredictor, cumulativeContErr, rhoMax, rhoMin, pRefCell, pRefValue ) :
    p_rgh.storePrevIter()
    rho.storePrevIter()
    
    UEqn = fun_UEqn( rho, U, phi, ghf, p_rgh, turb, mesh )
    fun_hEqn( U, rho, h, phi, turb, thermo, rad, mesh )
    
    # --- PISO loop
    cumulativeContErr =  fun_pEqn( runTime, i, mesh, p, rho, turb, thermo, thermoFluid, kappa, UEqn, U, phi, psi, \
                                   initialMass, p_rgh, gh, ghf, nNonOrthCorr, cumulativeContErr, rhoMax, rhoMin, pRefCell, pRefValue )
    
    turb.correct()
    
    return cumulativeContErr


#-----------------------------------------------------------------------------------------------------------------
