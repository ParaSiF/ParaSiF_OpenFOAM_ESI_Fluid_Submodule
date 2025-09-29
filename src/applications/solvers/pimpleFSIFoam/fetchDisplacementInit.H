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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // Debug switch for fetchDisplacementInit.H
    bool debugFetchDisplacementInit = false;
    
    // Declare MUI FSI residual output.
    fileName muiResidualFile("muiFSIResidual.txt");
    OFstream muiResidual(runTime.system()/muiResidualFile);

    pointVectorField* pointDisplacement = 
        const_cast<pointVectorField*>
        (
            &mesh.lookupObject<pointVectorField>
            (
                "pointDisplacement"
            )
        );

    List<List<vectorField>> dispVectorValsList;
    List<List<vectorField>> dispVectorOldValsList;
    List<List<vectorField>> oriCoordFetchPointsVectorValsList;

    bool rankHasPatchPointsForfetch = false;

    forAll(interfaceNames, iN)
    {

        List<vectorField> dispVectorValsListLocal;
        List<vectorField> dispVectorOldValsListLocal;
        List<vectorField> oriCoordFetchPointsVectorValsListLocal;

        forAll(patchIDs[iN], pI)
        {

            if(pointDisplacement->boundaryFieldRef()[patchIDs[iN][pI]].size() != 0)
            {
                rankHasPatchPointsForfetch = true;
            }

            // Declare a vectorField to store the fetched displacements.
            vectorField dispVectorValsLocal(pointDisplacement->boundaryFieldRef()[patchIDs[iN][pI]].size());
            vectorField dispVectorOldValsLocal(pointDisplacement->boundaryFieldRef()[patchIDs[iN][pI]].size());
            // Declare a vectorField to store the original point coordinates.
            vectorField oriCoordFetchPointsVectorValsLocal(pointDisplacement->boundaryFieldRef()[patchIDs[iN][pI]].size());

            dispVectorValsListLocal.append(dispVectorValsLocal);
            dispVectorOldValsListLocal.append(dispVectorOldValsLocal);
            oriCoordFetchPointsVectorValsListLocal.append(oriCoordFetchPointsVectorValsLocal);
        }

        dispVectorValsList.append(dispVectorValsListLocal);
        dispVectorOldValsList.append(dispVectorOldValsListLocal);
        oriCoordFetchPointsVectorValsList.append(oriCoordFetchPointsVectorValsListLocal);
    }

        scalar recvXMin = VGREAT;
        scalar recvYMin = VGREAT;
        scalar recvZMin = VGREAT;

        scalar recvXMax = -VGREAT;
        scalar recvYMax = -VGREAT;
        scalar recvZMax = -VGREAT;

    forAll(interfaceNames, iN)
    {

        forAll(patchIDs[iN], pI)
        {
            // Loop over points to Store the original point coordinates
            forAll(dispVectorValsList[iN][pI], pointI)
            {
                oriCoordFetchPointsVectorValsList[iN][pI][pointI].component(vector::X) =
                    mesh.boundaryMesh()[patchIDs[iN][pI]].localPoints()[pointI].component(vector::X);
                oriCoordFetchPointsVectorValsList[iN][pI][pointI].component(vector::Y) =
                    mesh.boundaryMesh()[patchIDs[iN][pI]].localPoints()[pointI].component(vector::Y);
                oriCoordFetchPointsVectorValsList[iN][pI][pointI].component(vector::Z) =
                    mesh.boundaryMesh()[patchIDs[iN][pI]].localPoints()[pointI].component(vector::Z);

                dispVectorValsList[iN][pI][pointI].component(vector::X) = 0.0;
                dispVectorValsList[iN][pI][pointI].component(vector::Y) = 0.0;
                dispVectorValsList[iN][pI][pointI].component(vector::Z) = 0.0;

                dispVectorOldValsList[iN][pI][pointI].component(vector::X) = 0.0;
                dispVectorOldValsList[iN][pI][pointI].component(vector::Y) = 0.0;
                dispVectorOldValsList[iN][pI][pointI].component(vector::Z) = 0.0;
            }
        }
    }

    const pointField& meshPoints = mesh.points();


    if (rankHasPatchPointsForfetch)
    {
        forAll(interfaceNames, iN)
        {

            forAll(patchIDs[iN], pI)
            {
                forAll(oriCoordFetchPointsVectorValsList[iN][pI], pointI)
                {
                    recvXMin =
                        ((oriCoordFetchPointsVectorValsList[iN][pI][pointI].component(vector::X)) < recvXMin) ? (oriCoordFetchPointsVectorValsList[iN][pI][pointI].component(vector::X)) : recvXMin;
                    recvYMin =
                        ((oriCoordFetchPointsVectorValsList[iN][pI][pointI].component(vector::Y)) < recvYMin) ? (oriCoordFetchPointsVectorValsList[iN][pI][pointI].component(vector::Y)) : recvYMin;
                    recvZMin =
                        ((oriCoordFetchPointsVectorValsList[iN][pI][pointI].component(vector::Z)) < recvZMin) ? (oriCoordFetchPointsVectorValsList[iN][pI][pointI].component(vector::Z)) : recvZMin;

                    recvXMax =
                        ((oriCoordFetchPointsVectorValsList[iN][pI][pointI].component(vector::X)) > recvXMax) ? (oriCoordFetchPointsVectorValsList[iN][pI][pointI].component(vector::X)) : recvXMax;
                    recvYMax =
                        ((oriCoordFetchPointsVectorValsList[iN][pI][pointI].component(vector::Y)) > recvYMax) ? (oriCoordFetchPointsVectorValsList[iN][pI][pointI].component(vector::Y)) : recvYMax;
                    recvZMax =
                        ((oriCoordFetchPointsVectorValsList[iN][pI][pointI].component(vector::Z)) > recvZMax) ? (oriCoordFetchPointsVectorValsList[iN][pI][pointI].component(vector::Z)) : recvZMax;
                }
            }
        }
    } else
    {
        forAll(meshPoints, pts)
        {
            recvXMin =
                (meshPoints[pts][0] < recvXMin) ? meshPoints[pts][0] : recvXMin;
            recvYMin =
                (meshPoints[pts][1] < recvYMin) ? meshPoints[pts][1] : recvYMin;
            recvZMin =
                (meshPoints[pts][2] < recvZMin) ? meshPoints[pts][2] : recvZMin;

            recvXMax =
                (meshPoints[pts][0] > recvXMax) ? meshPoints[pts][0] : recvXMax;
            recvYMax =
                (meshPoints[pts][1] > recvYMax) ? meshPoints[pts][1] : recvYMax;
            recvZMax =
                (meshPoints[pts][2] > recvZMax) ? meshPoints[pts][2] : recvZMax;
        }
    }

    // create on the stack
    mui::algo_aitken<mui::mui_config> aitken_x(static_cast<double>(initUndRelxCpl), 1.0);
    mui::algo_aitken<mui::mui_config> aitken_y(static_cast<double>(initUndRelxCpl), 1.0);
    mui::algo_aitken<mui::mui_config> aitken_z(static_cast<double>(initUndRelxCpl), 1.0);

    mui::algo_fixed_relaxation<mui::mui_config> fr_x(static_cast<double>(initUndRelxCpl));
    mui::algo_fixed_relaxation<mui::mui_config> fr_y(static_cast<double>(initUndRelxCpl));
    mui::algo_fixed_relaxation<mui::mui_config> fr_z(static_cast<double>(initUndRelxCpl));

    // define spatial and temporal samplers
    mui::sampler_pseudo_n2_linear<mui::mui_config> spatial_sampler(rSampler);
    mui::temporal_sampler_exact<mui::mui_config> chrono_sampler;

    std::vector<std::pair<mui::point3d, double>> ptsVluInit;
    mui::point3d rcv_point;
    int count=0;

    forAll(interfaceNames, iN)
    {
        forAll(patchIDs[iN], pI)
        {
            // Loop over points to receive displacement values and assigned to dispVectorVals
            //- Blocking fetch for value passed through interface from coupled solver

            forAll(dispVectorValsList[iN][pI], pointI)
            {
                rcv_point[0] = oriCoordFetchPointsVectorValsList[iN][pI][pointI].component(vector::X);
                rcv_point[1] = oriCoordFetchPointsVectorValsList[iN][pI][pointI].component(vector::Y);
                rcv_point[2] = oriCoordFetchPointsVectorValsList[iN][pI][pointI].component(vector::Z);

                if (rankHasPatchPointsForfetch)
                {
                    ptsVluInit.push_back(std::make_pair(rcv_point,0.0));
                    ++count;
                }
            }
        }
    }

    int totPointSize = 0;

    forAll(interfaceNames, iN)
    {

        forAll(patchIDs[iN], pI)
        {

            totPointSize += pointDisplacement->boundaryFieldRef()[patchIDs[iN][pI]].size();
        }

        if (cplMethod == "Aitken")
        {
            Info << "{OpenFOAM} interface " << interfaceNames[iN] << " Aitken coupling method selected!" << endl;

            muiResidual << "MUI CouplingLib Utility" << "\t" << "Aitken method" << endl;
            muiResidual << "interfaceNames"<< "\t" << "Time" << "\t" << "Time-Step" << "\t" << "Sub-Iteration" << "\t"
                << "UnderRelx-Factor" << "\t" << "pointSize" << "\t"
                << "squareSumResidual" << "\t" << "maxResidualL-2Norm" << "\t"
                << "residualL2Norm" << "\t" << "forceX"
                        << "\t" << "forceY" << "\t" << "forceZ"<< endl;

        } else if (cplMethod == "FixedRelaxation")
        {
            Info << "{OpenFOAM} interface " << interfaceNames[iN] << " Fixed Relaxation coupling method selected!" << endl;

            muiResidual << "MUI CouplingLib Utility" << "\t" << "FixedRelaxation method" << endl;
            muiResidual << "interfaceNames"<< "\t" << "Time" << "\t" << "Time-Step" << "\t" << "Sub-Iteration" << "\t"
                << "UnderRelx-Factor" << "\t" << "pointSize" << "\t"
                << "squareSumResidual" << "\t" << "maxResidualL-2Norm" << "\t"
                << "residualL2Norm" << "\t" << "forceX"
                        << "\t" << "forceY" << "\t" << "forceZ"<< endl;

        } else if (cplMethod == "Loose")
        {
            Info << "{OpenFOAM} Loose coupling method selected!" << endl;
        } else
        {
            FatalError << "cplFSIMethod: " << cplMethod << " in fsiDict is not recognized!" << exit(FatalError);
        }
    }

    if (debugFetchDisplacementInit == true)
    {
        Pout << "recvXMin = " << recvXMin << endl;
        Pout << "recvYMin = " << recvYMin << endl;
        Pout << "recvZMin = " << recvZMin << endl;

        Pout << "recvXMax = " << recvXMax << endl;
        Pout << "recvYMax = " << recvYMax << endl;
        Pout << "recvZMax = " << recvZMax << endl;
    }

    // Annouce send span
    mui::geometry::box<mui::mui_config> recv_region( {recvXMin, recvYMin, recvZMin}, {recvXMax, recvYMax, recvZMax} );

    forAll(interfaceNames, iN)
    {
        ifs[iN]->announce_recv_span( runTime.startTime().value(), runTime.endTime().value(), recv_region );

        // Commit ZERO step
        ifs[iN]->commit(0);
        Info << "{OpenFOMA} : MUI interface " << interfaceNames[iN] << " commit zero step" << nl << endl;
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// *************************  FILE END  ************************************ //
