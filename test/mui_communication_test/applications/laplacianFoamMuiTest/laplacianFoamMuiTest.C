/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    laplacianFoamMuiTest

Description
    MUI coupling test based on laplace equation solver for a scalar quantity.

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

#include "mui.h"
#include "mui_config.h"
#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "MUI coupling test based on laplace equation solver for a scalar quantity."
    );

    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"

    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCreate MUI unifaces\n" << endl;
    // Define the name of MUI interfaceName
    List<word> interfaceNames;
    // Define the name of MUI domain
    word domainName;
    
    interfaceNames.append(static_cast<word>("testInterface"));
    domainName = static_cast<word>("OF0");

    // Define the name of MUI interfaces
    std::vector<std::string> interfaces;

    forAll(interfaceNames, iN)
    {
        interfaces.emplace_back(interfaceNames[iN]);
    }

    // Declare MUI objects using MUI configure file
    auto ifs = mui::create_uniface<mui::mui_config>( domainName, interfaces );

    // setup parameters
    constexpr static int    Nx        = 41; // number of grid points in x axis
    constexpr static int    Ny        = 5; // number of grid points in y axis
    constexpr static int    Nz        = 5; // number of grid points in z axis
    const char* name_fetchX = "dispX";
    const char* name_fetchY = "dispY";
    const char* name_fetchZ = "dispZ";
    const char* name_pushX0 = "forceX0";
    const char* name_pushY0 = "forceY0";
    const char* name_pushZ0 = "forceZ0";
    const char* name_pushX1 = "forceX1";
    const char* name_pushY1 = "forceY1";
    const char* name_pushZ1 = "forceZ1";
    double r    = 1.0;                      // search radius
    int Nt = Nx * Ny * Nz; // total time steps
    int steps = 1000; // total time steps
    int nSubIter = 1; // total time steps
    double timeStepSize    = 0.1;
    double local_x0 = 0.; // local origin
    double local_y0 = 0.;
    double local_z0 = 0.;
    double local_x1 = 20.;
    double local_y1 = 2.;
    double local_z1 = 2.;
    double local_x2 = 0.; // local origin
    double local_y2 = 0.;
    double local_z2 = 0.;
    double local_x3 = 20.;
    double local_y3 = 2.;
    double local_z3 = 2.;
    double monitorX = 20.;
    double monitorY = 0.;
    double monitorZ = 0.;
    double pp[Nx][Ny][Nz][3], pf[Nx][Ny][Nz][3];
    double force_pushX[Nx][Ny][Nz],force_pushY[Nx][Ny][Nz],force_pushZ[Nx][Ny][Nz];
    double displacement_fetchX[Nx][Ny][Nz], displacement_fetchY[Nx][Ny][Nz], displacement_fetchZ[Nx][Ny][Nz];
    double displacement_fetchX_Store[Nx][Ny][Nz], displacement_fetchY_Store[Nx][Ny][Nz], displacement_fetchZ_Store[Nx][Ny][Nz];

    // Push points generation and evaluation
    for ( int i = 0; i < Nx; ++i ) {
        for ( int j = 0; j < Ny; ++j ) {
            for ( int k = 0; k < Nz; ++k ) {
                double x = local_x0+(i*(local_x1-local_x0)/(Nx-1));
                double y = local_y0+(j*(local_y1-local_y0)/(Ny-1));
                double z = local_z0+(k*(local_z1-local_z0)/(Nz-1));
                pp[i][j][k][0] = x;
                pp[i][j][k][1] = y;
                pp[i][j][k][2] = z;
                force_pushX[i][j][k] = 0.;
                force_pushY[i][j][k] = 0.;
                force_pushZ[i][j][k] = 0.;
            }
        }
    }

    // Fetch points generation and evaluation
    for ( int i = 0; i < Nx; ++i ) {
        for ( int j = 0; j < Ny; ++j ) {
            for ( int k = 0; k < Nz; ++k ) {
                double x = local_x2+(i*(local_x3-local_x2)/(Nx-1));
                double y = local_y2+(j*(local_y3-local_y2)/(Ny-1));
                double z = local_z2+(k*(local_z3-local_z2)/(Nz-1));
                pf[i][j][k][0] = x;
                pf[i][j][k][1] = y;
                pf[i][j][k][2] = z;
                displacement_fetchX[i][j][k] = 0.0;
                displacement_fetchY[i][j][k] = 0.0;
                displacement_fetchZ[i][j][k] = 0.0;
                displacement_fetchX_Store[i][j][k] = 0.0;
                displacement_fetchY_Store[i][j][k] = 0.0;
                displacement_fetchZ_Store[i][j][k] = 0.0;
            }
        }
    }

   // annouce send span
    mui::geometry::box<mui::mui_config> send_region( {local_x0, local_y0, local_z0}, {local_x1, local_y1, local_z1} );
    mui::geometry::box<mui::mui_config> recv_region( {local_x2, local_y2, local_z2}, {local_x3, local_y3, local_z3} );
    printf( "{OF0} send region for: %lf %lf %lf - %lf %lf %lf\n", local_x0, local_y0, local_z0, local_x1, local_y1, local_z1 );
    ifs[0]->announce_send_span( 0, steps*10, send_region );
    ifs[0]->announce_recv_span( 0, steps*10, recv_region );

    // define spatial and temporal samplers
    mui::sampler_pseudo_n2_linear<mui::mui_config> s1(r);
    mui::temporal_sampler_exact<mui::mui_config> s2;

    if (Pstream::myProcNo() == 0) {
        mui::point3d locp0( 13, 12, 11 );
        ifs[0]->push( "OF1", locp0, 888.78 );
    } else {
        mui::point3d locp0( 12, 11, 10 );
        ifs[0]->push( "OF1", locp0, 888.78 );
    }

    // commit ZERO step
    ifs[0]->commit(0);
    mui::point3d locp1( 23, 22, 21 );
    double receive = ifs[0]->fetch( "PUSHER1", locp1, 0, s1, s2 );
    std::cout << "{OF0} 0 receive: " <<  receive << std::endl;
   // Commit ZERO step
   ifs[0]->commit(0);
   Info << "MUI interface " << interfaceNames[0] << " commit zero step" << nl << endl;

   // Begin time loops
   for ( int n = 1; n <= steps; ++n ) {

       printf("\n");
       printf("{OF0} %d Step \n", n);

       // Begin iteration loops
       for ( int iter = 1; iter <= nSubIter; ++iter ) {

           printf("{OF0} %d iteration \n", iter);

           int totalIter = ( (n - 1) * nSubIter ) + iter;
           double total_force_Y=0.0;
           // push data to the other solver
           for ( int i = 0; i < Nx; ++i ) {
               for ( int j = 0; j < Ny; ++j ) {
                   for ( int k = 0; k < Nz; ++k ) {
                       force_pushX[i][j][k] = 0.;
                       force_pushY[i][j][k] = -((n*timeStepSize)*(10)/7.0);
                       force_pushZ[i][j][k] = 0.;

                       if (std::abs(pp[i][j][k][0] - 20.0) <= 0.00001 ){
                           mui::point3d locp( pp[i][j][k][0], pp[i][j][k][1], pp[i][j][k][2] );
                            if (Pstream::myProcNo() == 0) {
                                ifs[0]->push( name_pushX0, locp, force_pushX[i][j][k] );
                                ifs[0]->push( name_pushY0, locp, force_pushY[i][j][k] );
                                ifs[0]->push( name_pushZ0, locp, force_pushZ[i][j][k] );
                            } else {
                                ifs[0]->push( name_pushX1, locp, force_pushX[i][j][k] );
                                ifs[0]->push( name_pushY1, locp, force_pushY[i][j][k] );
                                ifs[0]->push( name_pushZ1, locp, force_pushZ[i][j][k] );
                            }
                            std::cout << "{OF0} push point: " <<  locp[0] << ", " <<  locp[1] << ", "<<  locp[2] << " with value " << force_pushY[i][j][k] << std::endl;
                            total_force_Y += force_pushY[i][j][k];
                       }
                   }
               }
           }
           printf( "{OF0} total_force_Y: %lf at time: %f [s]\n", total_force_Y, (n*timeStepSize));
           int sent = ifs[0]->commit( totalIter );
           if ((totalIter-1)>=1){
               // push data to the other solver
               for ( int i = 0; i < Nx; ++i ) {
                   for ( int j = 0; j < Ny; ++j ) {
                       for ( int k = 0; k < Nz; ++k ) {
                           mui::point3d locf( pf[i][j][k][0], pf[i][j][k][1], pf[i][j][k][2] );
                           displacement_fetchX[i][j][k] = ifs[0]->fetch( name_fetchX, locf,
                               (totalIter-1),
                               s1,
                               s2 );
                           displacement_fetchY[i][j][k] = ifs[0]->fetch( name_fetchY, locf,
                               (totalIter-1),
                               s1,
                               s2 );
                           displacement_fetchZ[i][j][k] = ifs[0]->fetch( name_fetchZ, locf,
                               (totalIter-1),
                               s1,
                               s2 );
                           if ((std::abs(pf[i][j][k][0] - 20.0) < 0.0001) && (std::abs(pf[i][j][k][1] - 0.0) < 0.0001) && (std::abs(pf[i][j][k][2] - 0.0) < 0.0001)) {
                               printf( "{OF0} fetch disp : %lf, %lf, %lf at time: %f [s], i: %d; j: %d; k: %d; pf[i][j][k][0]: %lf\n", displacement_fetchX[i][j][k], displacement_fetchY[i][j][k], displacement_fetchZ[i][j][k], (n*timeStepSize), i, j, k, pf[i][j][k][0]);
                           }
                       }
                   }
               }
           }
       }

   }
    
    Info<< "\nCalculating temperature distribution\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                fvm::ddt(T) - fvm::laplacian(DT, T)
             ==
                fvOptions(T)
            );

            fvOptions.constrain(TEqn);
            TEqn.solve();
            fvOptions.correct(T);
        }

        #include "write.H"

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;
    mui::mpi_finalize_after_split();
    return 0;
}


// ************************************************************************* //
