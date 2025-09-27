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
#include "pusher_fetcher_config.h"

int main(int argc, char ** argv) {
    using namespace mui;

    MPI_Comm  world = mui::mpi_split_by_app();

    std::ofstream displacementOutFile("dispCpp.txt");

    // Define the name of MUI interfaces
    std::vector<std::string> interfaces;
    std::string domainName="PUSHER_FETCHER_1";
    std::string appName="testInterface";

    interfaces.emplace_back(appName);

    // Declare MUI objects using MUI configure file
    auto ifs = mui::create_uniface<mui::config_3d>( domainName, interfaces );

    int rank, size;
    MPI_Comm_rank( world, &rank );
    MPI_Comm_size( world, &size );

    // setup parameters
    constexpr static int    Nx        = 41; // number of grid points in x axis
    constexpr static int    Ny        = 5; // number of grid points in y axis
    constexpr static int    Nz        = 5; // number of grid points in z axis
    const char* name_fetchX0 = "dispX0";
    const char* name_fetchY0 = "dispY0";
    const char* name_fetchZ0 = "dispZ0";
    const char* name_fetchX1 = "dispX1";
    const char* name_fetchY1 = "dispY1";
    const char* name_fetchZ1 = "dispZ1";
    const char* name_pushX = "forceX";
    const char* name_pushY = "forceY";
    const char* name_pushZ = "forceZ";
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
    mui::geometry::box<mui::config_3d> send_region( {local_x0, local_y0, local_z0}, {local_x1, local_y1, local_z1} );
    mui::geometry::box<mui::config_3d> recv_region( {local_x2, local_y2, local_z2}, {local_x3, local_y3, local_z3} );
    ifs[0]->announce_send_span( 0, steps*10, send_region );
    ifs[0]->announce_recv_span( 0, steps*10, recv_region );

    // define spatial and temporal samplers
    mui::sampler_pseudo_n2_linear<mui::config_3d> s1(r);
    mui::temporal_sampler_exact<mui::config_3d> s2;

    point3d locp0( 23, 22, 21 );
    ifs[0]->push( "PUSHER1", locp0, 888.78 );
    // commit ZERO step
    ifs[0]->commit(0);
    point3d locr0( 13, 12, 11 );
    point3d locr1( 12, 11, 10 );
    double receive0 = ifs[0]->fetch( "OF1", locr0, 0, s1, s2 );
    double receive1 = ifs[0]->fetch( "OF1", locr1, 0, s1, s2 );
    std::cout << "{PUSHER_FETCHER_1} 0 receive0: " <<  receive0 << std::endl;
    std::cout << "{PUSHER_FETCHER_1} 0 receive1: " <<  receive1 << std::endl;
    // Begin time loops
    for ( int n = 1; n <= steps; ++n ) {

        printf("\n");
        printf("{PUSHER_FETCHER_1} %d Step \n", n);

        // Begin iteration loops
        for ( int iter = 1; iter <= nSubIter; ++iter ) {

            printf("{PUSHER_FETCHER_1} %d iteration \n", iter);

            int totalIter = ( (n - 1) * nSubIter ) + iter;
            double total_force_Y=0.0;
            // push data to the other solver
            for ( int i = 0; i < Nx; ++i ) {
                for ( int j = 0; j < Ny; ++j ) {
                    for ( int k = 0; k < Nz; ++k ) {
                        if (((n*timeStepSize)<=7.)&& (i==(Nx-1))){
                            force_pushX[i][j][k] = 0.;
                            force_pushY[i][j][k] = -((n*timeStepSize)*(20)/7.0);
                            force_pushZ[i][j][k] = 0.;
                        }else{
                            force_pushX[i][j][k] = 0.;
                            force_pushY[i][j][k] = 0.;
                            force_pushZ[i][j][k] = 0.;
                        }

                        if (std::abs(pp[i][j][k][0] - 20.0) <= 0.00001 ){
                            point3d locp( pp[i][j][k][0], pp[i][j][k][1], pp[i][j][k][2] );
                            ifs[0]->push( name_pushX, locp, force_pushX[i][j][k] );
                            ifs[0]->push( name_pushY, locp, force_pushY[i][j][k] );
                            ifs[0]->push( name_pushZ, locp, force_pushZ[i][j][k] );
                            std::cout << "{PUSHER_FETCHER_1} push point: " <<  locp[0] << ", " <<  locp[1] << ", "<<  locp[2] << std::endl;
                            total_force_Y += force_pushY[i][j][k];
                        }
                    }
                }
            }
            printf( "{PUSHER_FETCHER_1} total_force_Y: %lf at time: %f [s]\n", total_force_Y, (n*timeStepSize));
            int sent = ifs[0]->commit( totalIter );
            if ((totalIter-1)>=1){
                // push data to the other solver
                for ( int i = 0; i < Nx; ++i ) {
                    for ( int j = 0; j < Ny; ++j ) {
                        for ( int k = 0; k < Nz; ++k ) {
                            point3d locf( pf[i][j][k][0], pf[i][j][k][1], pf[i][j][k][2] );
                            displacement_fetchX[i][j][k] = ifs[0]->fetch( name_fetchX0, locf,
                                (totalIter-1),
                                s1,
                                s2 );
                            displacement_fetchY[i][j][k] = ifs[0]->fetch( name_fetchY0, locf,
                                (totalIter-1),
                                s1,
                                s2 );
                            displacement_fetchZ[i][j][k] = ifs[0]->fetch( name_fetchZ0, locf,
                                (totalIter-1),
                                s1,
                                s2 );
                            if ((std::abs(pf[i][j][k][0] - 20.0) < 0.0001) && (std::abs(pf[i][j][k][1] - 0.0) < 0.0001) && (std::abs(pf[i][j][k][2] - 0.0) < 0.0001)) {
                                printf( "{PUSHER_FETCHER_1} fetch 0 disp : %lf, %lf, %lf at time: %f [s], i: %d; j: %d; k: %d; pf[i][j][k][0]: %lf\n", displacement_fetchX[i][j][k], displacement_fetchY[i][j][k], displacement_fetchZ[i][j][k], (n*timeStepSize), i, j, k, pf[i][j][k][0]);
                            }

                            displacement_fetchX[i][j][k] = ifs[0]->fetch( name_fetchX1, locf,
                                (totalIter-1),
                                s1,
                                s2 );
                            displacement_fetchY[i][j][k] = ifs[0]->fetch( name_fetchY1, locf,
                                (totalIter-1),
                                s1,
                                s2 );
                            displacement_fetchZ[i][j][k] = ifs[0]->fetch( name_fetchZ1, locf,
                                (totalIter-1),
                                s1,
                                s2 );
                            if ((std::abs(pf[i][j][k][0] - 20.0) < 0.0001) && (std::abs(pf[i][j][k][1] - 0.0) < 0.0001) && (std::abs(pf[i][j][k][2] - 0.0) < 0.0001)) {
                                printf( "{PUSHER_FETCHER_1} fetch 1 disp : %lf, %lf, %lf at time: %f [s], i: %d; j: %d; k: %d; pf[i][j][k][0]: %lf\n", displacement_fetchX[i][j][k], displacement_fetchY[i][j][k], displacement_fetchZ[i][j][k], (n*timeStepSize), i, j, k, pf[i][j][k][0]);
                            }
                        }
                    }
                }
            }

            for ( int i = 0; i < Nx; ++i ) {
                for ( int j = 0; j < Ny; ++j ) {
                    for ( int k = 0; k < Nz; ++k ) {
                        if ((pf[i][j][k][0] == monitorX) && (pf[i][j][k][1] == monitorY) && (pf[i][j][k][2] == monitorZ)) {
                            displacementOutFile.open("dispCpp.txt", std::ios_base::app);
                            displacementOutFile << n*timeStepSize << " " << displacement_fetchY[i][j][k] << std::endl;
                            displacementOutFile.close();
                        }
                    }
                }
            }

        }

    }

    return 0;
}