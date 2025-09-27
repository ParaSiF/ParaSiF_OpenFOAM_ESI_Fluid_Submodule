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

/*
 * 3D_CPP_PUSHER_FETCHER_1.cpp
 *
 *  Created on: 10 Jan 2019
 *      Author: Wendi Liu
 */

#include "mui.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include "pusher_fetcher_config.h"

int main(int argc, char ** argv) {
    using namespace mui;

    MPI_Comm  world = mui::mpi_split_by_app();

    std::ofstream displacementOutFile("dispCpp.txt");

    // Define the name of MUI interfaces
    std::vector<std::string> interfaces;
    std::string domainName="PUSHER_FETCHER_1";
    std::string appName="threeDInterface0";

    interfaces.emplace_back(appName);

    // Declare MUI objects using MUI configure file
    auto ifs = mui::create_uniface<mui::pusher_fetcher_config>( domainName, interfaces );

    int rank, size;
    MPI_Comm_rank( world, &rank );
    MPI_Comm_size( world, &size );

    // setup parameters
    constexpr static int    Nx        = 41; // number of grid points in x axis
    constexpr static int    Ny        = 5; // number of grid points in y axis
    constexpr static int    Nz        = 5; // number of grid points in z axis
    const char* name_fetchX = "forceX";
    const char* name_fetchY = "forceY";
    const char* name_fetchZ = "forceZ";
    const char* name_pushX = "dispX";
    const char* name_pushY = "dispY";
    const char* name_pushZ = "dispZ";
    double r    = 1.0;                      // search radius
    int Nt = Nx * Ny * Nz; // total time steps
    int steps = 100; // total time steps
    int nSubIter = 5; // total time steps
    double timeStepSize    = 0.025;
    double local_x0 = 0.45; // local origin
    double local_y0 = 0.;
    double local_z0 = 0.;
    double local_x1 = 0.55;
    double local_y1 = 0.2;
    double local_z1 = -0.18;
    double local_x2 = 0.45; // local origin
    double local_y2 = 0.;
    double local_z2 = 0.;
    double local_x3 = 0.55;
    double local_y3 = 0.2;
    double local_z3 = -0.18;
    double monitorX = 0.5;
    double monitorY = 0.2;
    double monitorZ = 0.;
    double pp[Nx][Ny][Nz][3], pf[Nx][Ny][Nz][3];
    double displacement_pushX[Nx][Ny][Nz],displacement_pushY[Nx][Ny][Nz],displacement_pushZ[Nx][Ny][Nz];
    double force_fetchX[Nx][Ny][Nz], force_fetchY[Nx][Ny][Nz], force_fetchZ[Nx][Ny][Nz];

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
                displacement_pushX[i][j][k] = 0.;
                displacement_pushY[i][j][k] = 0.;
                displacement_pushZ[i][j][k] = 0.;
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
                force_fetchX[i][j][k] = 0.0;
                force_fetchY[i][j][k] = 0.0;
                force_fetchZ[i][j][k] = 0.0;
            }
        }
    }

   // annouce send span
    mui::geometry::box<mui::pusher_fetcher_config> send_region( {local_x0, local_y0, local_z0}, {local_x1, local_y1, local_z1} );
    mui::geometry::box<mui::pusher_fetcher_config> recv_region( {local_x2, local_y2, local_z2}, {local_x3, local_y3, local_z3} );
    ifs[0]->announce_send_span( 0, steps*10, send_region );
    ifs[0]->announce_recv_span( 0, steps*10, recv_region );

    // define spatial and temporal samplers
    mui::sampler_pseudo_n2_linear<mui::pusher_fetcher_config> s1(r);
    mui::temporal_sampler_exact<mui::pusher_fetcher_config> s2;

    // commit ZERO step
    ifs[0]->commit(0);
    int totalIter = 0; 
    // Begin time loops
    for ( int n = 1; n <= steps; ++n ) {

        printf("\n");
        printf("{PUSHER_FETCHER_1} Time Step %d \n", n);

        // Begin iteration loops
        for ( int iter = 1; iter <= nSubIter; ++iter ) {

            printf("{PUSHER_FETCHER_1} %d iteration \n", iter);

            totalIter = ( (n - 1) * nSubIter ) + iter;
            // push data to the other solver
            for ( int i = 0; i < Nx; ++i ) {
                for ( int j = 0; j < Ny; ++j ) {
                    for ( int k = 0; k < Nz; ++k ) {
                        displacement_pushX[i][j][k] = (0.625 * (pp[i][j][k][1] * pp[i][j][k][1])) * sin(2 * M_PI * 0.5 * (timeStepSize * n));
                        displacement_pushY[i][j][k] = 0.;
                        displacement_pushZ[i][j][k] = 0.;

                        point3d locp( pp[i][j][k][0], pp[i][j][k][1], pp[i][j][k][2] );
                        ifs[0]->push( name_pushX, locp, displacement_pushX[i][j][k] );
                        ifs[0]->push( name_pushY, locp, displacement_pushY[i][j][k] );
                        ifs[0]->push( name_pushZ, locp, displacement_pushZ[i][j][k] );
                        std::cout << "{PUSHER_FETCHER_1} push point: " <<  locp[0] << ", " <<  locp[1] << ", "<<  locp[2] << " with disp: " << displacement_pushX[i][j][k] << ", " << displacement_pushY[i][j][k] << ", " << displacement_pushZ[i][j][k] << std::endl;
                    }
                }
            }
            int sent = ifs[0]->commit( totalIter );
            printf( "{PUSHER_FETCHER_1} commit at: %d \n", (totalIter));
            if ((totalIter-0)>=1){
                // push data to the other solver
                for ( int i = 0; i < Nx; ++i ) {
                    for ( int j = 0; j < Ny; ++j ) {
                        for ( int k = 0; k < Nz; ++k ) {
                            point3d locf( pf[i][j][k][0], pf[i][j][k][1], pf[i][j][k][2] );
                            force_fetchX[i][j][k] = ifs[0]->fetch( name_fetchX, locf,
                                (totalIter-0),
                                s1,
                                s2 );
                            force_fetchY[i][j][k] = ifs[0]->fetch( name_fetchY, locf,
                                (totalIter-0),
                                s1,
                                s2 );
                            force_fetchZ[i][j][k] = ifs[0]->fetch( name_fetchZ, locf,
                                (totalIter-0),
                                s1,
                                s2 );
                            if ((std::abs(pf[i][j][k][0] - monitorX) < 0.0001) && (std::abs(pf[i][j][k][1] - monitorY) < 0.0001) && (std::abs(pf[i][j][k][2] - monitorZ) < 0.0001)) {
                                std::cout << "{PUSHER_FETCHER_1} fetch force: " <<  force_fetchX[i][j][k] << ", " << force_fetchY[i][j][k] << ", " << force_fetchZ[i][j][k] << " at iteration " << (totalIter-0) << std::endl;
                            }
                        }
                    }
                }
            }

            for ( int i = 0; i < Nx; ++i ) {
                for ( int j = 0; j < Ny; ++j ) {
                    for ( int k = 0; k < Nz; ++k ) {
                        if ((std::abs(pf[i][j][k][0] - monitorX) < 0.0001) && (std::abs(pf[i][j][k][1] - monitorY) < 0.0001) && (std::abs(pf[i][j][k][2] - monitorZ) < 0.0001) && (iter == nSubIter)) {
                            displacementOutFile.open("dispCpp.txt", std::ios_base::app);
                            displacementOutFile << n*timeStepSize << " " << force_fetchY[i][j][k] << std::endl;
                            displacementOutFile.close();
                        }
                    }
                }
            }
        }
    }

    MPI_Finalize();

    return 0;
}