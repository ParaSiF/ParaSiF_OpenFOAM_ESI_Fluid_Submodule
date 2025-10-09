/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    interMUIFoam

Group
    grpMultiphaseSolvers

Description
    Solver for two incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing with two-way data exchange with Python script undet Numpy data
    format on cell center values.

\*---------------------------------------------------------------------------*/

/*****************************************************************************
* Parallel Partitioned Multi-Physics Simulation Framework (ParaSiF)          *
*                                                                            *
* Copyright (C) 2025 The ParaSiF Development Team                            *
* All rights reserved                                                        *
*                                                                            *
* This software is licensed under the GNU General Public License version 3   *
*                                                                            *
* ** GNU General Public License, version 3 **                                *
*                                                                            *
* This program is free software: you can redistribute it and/or modify       *
* it under the terms of the GNU General Public License as published by       *
* the Free Software Foundation, either version 3 of the License, or          *
* (at your option) any later version.                                        *
*                                                                            *
* This program is distributed in the hope that it will be useful,            *
* but WITHOUT ANY WARRANTY; without even the implied warranty of             *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
* GNU General Public License for more details.                               *
*                                                                            *
* You should have received a copy of the GNU General Public License          *
* along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
*****************************************************************************/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
// * * * * * * * * * * * MUI Include * * * * * * * * * * * * * * * * * //
    #include "mui.h"
    #include "muiFiles/mui_config.h"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible, isothermal immiscible fluids"
        " using VOF phase-fraction based interface capturing.\n"
        "With optional mesh motion and mesh topology changes including"
        " adaptive re-meshing with two-way data exchange with Python"
        " script undet Numpy data format on cell center values."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createAlphaFluxes.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * MUI Initialise * * * * * * * * * * * * * * * * * //

    Info<< "\nCreate MUI unifaces\n" << endl;
    // Define the name of MUI interfaceName
    List<word> interfaceNames;
    // Define the name of MUI domain
    word domainName;
    // Define the name of push/fetch values
    word nameFetch = "prghDouble";
    word namePush = "prgh";

    interfaceNames.append(static_cast<word>("MLDataInterface"));
    domainName = static_cast<word>("OpenFOAM");

    // Define the name of MUI interfaces
    std::vector<std::string> interfaces;

    forAll(interfaceNames, iN)
    {
        interfaces.emplace_back(interfaceNames[iN]);
    }

    // Declare MUI objects using MUI configure file
    auto ifs = mui::create_uniface<mui::mui_config>( domainName, interfaces );

    const vectorField& centres = mesh.C();

    // Initialise min and max vectors
    vector minCoord(VGREAT, VGREAT, VGREAT);
    vector maxCoord(-VGREAT, -VGREAT, -VGREAT);

    forAll(centres, cellI)
    {
        const vector& c = centres[cellI];

        // component-wise comparison
        minCoord.x() = Foam::min(minCoord.x(), c.x());
        minCoord.y() = Foam::min(minCoord.y(), c.y());
        minCoord.z() = Foam::min(minCoord.z(), c.z());

        maxCoord.x() = Foam::max(maxCoord.x(), c.x());
        maxCoord.y() = Foam::max(maxCoord.y(), c.y());
        maxCoord.z() = Foam::max(maxCoord.z(), c.z());
    }

    label nProcs = Pstream::nProcs();

    ifs[0]->push( "nProcs", nProcs );

    word namePushSize = "cellSize_" + name(Pstream::myProcNo());

    ifs[0]->push( namePushSize, centres.size() );

    Pout << "{OpenFOAM} rank " << Pstream::myProcNo() << " has " << centres.size() << " cells." << nl;
/* 
    // Get the communicator OpenFOAM is using
    MPI_Comm comm = Pstream::mpiCommunicator();

    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    Info << "Rank " << rank << " of " << size 
         << " using OpenFOAM communicator." << nl; */

   // annouce send span
    mui::geometry::box<mui::mui_config> send_region( {minCoord.x(), minCoord.y(), minCoord.z()}, {maxCoord.x(), maxCoord.y(), maxCoord.z()} );
    mui::geometry::box<mui::mui_config> recv_region( {minCoord.x(), minCoord.y(), minCoord.z()}, {maxCoord.x(), maxCoord.y(), maxCoord.z()} );
    printf( "{OpenFOAM} send region for: %lf %lf %lf - %lf %lf %lf\n", minCoord.x(), minCoord.y(), minCoord.z(), maxCoord.x(), maxCoord.y(), maxCoord.z() );
    ifs[0]->announce_send_span( 0, label((runTime.endTime().value() - runTime.startTime().value())/runTime.deltaT().value() + 0.5), send_region );
    ifs[0]->announce_recv_span( 0, label((runTime.endTime().value() - runTime.startTime().value())/runTime.deltaT().value() + 0.5), recv_region );

    // define spatial and temporal samplers
    scalar r    = 1.0;  // search radius
    mui::sampler_pseudo_n2_linear<mui::mui_config> s1(r);
    mui::temporal_sampler_exact<mui::mui_config> s2;

    // commit ZERO step
    ifs[0]->commit(0);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                mesh.update();

                if (mesh.changing())
                {
                    // Do not apply previous time-step mesh compression flux
                    // if the mesh topology changed
                    if (mesh.topoChanging())
                    {
                        talphaPhi1Corr0.clear();
                    }

                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);

                        mixture.correct();
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            mixture.correct();

            if (pimple.frozenFlow())
            {
                continue;
            }

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
