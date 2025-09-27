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

    // Debug switch for fetchDisplacement.H
    bool debugFetchDisplacement = false;
    double res_x = 0.0,res_y = 0.0,res_z = 0.0;

    mui::point3d rcv_point;
    Info << "{OpenFOAM} Fetch iteration " << (totalCurrentIter - 0) << nl << endl;
    if((totalCurrentIter-0) >= 0)
    {
        forAll(interfaceNames, iN)
        {
            int pointCounter = 0;

            bool patchIDToFetch = false;

            forAll(patch6DoFIDs[iN], pI)
            {
                forAll(patchIDs[iN], ptI)
                {
                    if(patch6DoFIDs[iN][pI] == patchIDs[iN][ptI])
                    {
                        patchIDToFetch = true;
                        break;
                    }
                }
                // Loop over points to receive displacement values and assigned to dispVectorVals
                //- Blocking fetch for value passed through interface from coupled solver
                if(patchIDToFetch == true)
                {
                    forAll(dispVectorValsList[iN][pI], pointI)
                    {
                        ++pointCounter;

                        rcv_point[0] = oriCoordFetchPointsVectorValsList[iN][pI][pointI].component(vector::X);
                        rcv_point[1] = oriCoordFetchPointsVectorValsList[iN][pI][pointI].component(vector::Y);
                        rcv_point[2] = oriCoordFetchPointsVectorValsList[iN][pI][pointI].component(vector::Z);

                        if (rankHasPatchPointsForfetch)
                        {
                            double fetchTempX =0.0, fetchTempY = 0.0, fetchTempZ = 0.0;
                            if (cplMethod == "Aitken")
                            {
                                    fetchTempX =
                                            ifs[iN]->fetch(fetchNameX, rcv_point, static_cast<double>(totalCurrentIter-0), spatial_sampler, chrono_sampler,aitken_x);
                                    fetchTempY =
                                            ifs[iN]->fetch(fetchNameY, rcv_point, static_cast<double>(totalCurrentIter-0), spatial_sampler, chrono_sampler,aitken_y);
                                    fetchTempZ =
                                            ifs[iN]->fetch(fetchNameZ, rcv_point, static_cast<double>(totalCurrentIter-0), spatial_sampler, chrono_sampler,aitken_z);
                                    res_x = aitken_x.get_residual_L2_Norm(static_cast<double>(totalCurrentIter-0));
                                    res_y = aitken_y.get_residual_L2_Norm(static_cast<double>(totalCurrentIter-0));
                                    res_z = aitken_z.get_residual_L2_Norm(static_cast<double>(totalCurrentIter-0));

                            } else if (cplMethod == "FixedRelaxation")
                            {
                                    fetchTempX =
                                            ifs[iN]->fetch(fetchNameX, rcv_point, static_cast<double>(totalCurrentIter-0), spatial_sampler, chrono_sampler,fr_x);
                                    fetchTempY =
                                            ifs[iN]->fetch(fetchNameY, rcv_point, static_cast<double>(totalCurrentIter-0), spatial_sampler, chrono_sampler,fr_y);
                                    fetchTempZ =
                                            ifs[iN]->fetch(fetchNameZ, rcv_point, static_cast<double>(totalCurrentIter-0), spatial_sampler, chrono_sampler,fr_z);

                                    res_x = fr_x.get_residual_L2_Norm(static_cast<double>(totalCurrentIter-0));
                                    res_y = fr_y.get_residual_L2_Norm(static_cast<double>(totalCurrentIter-0));
                                    res_z = fr_z.get_residual_L2_Norm(static_cast<double>(totalCurrentIter-0));

                            }  else if (cplMethod == "Loose")
                            {
                                    fetchTempX =
                                            ifs[iN]->fetch(fetchNameX, rcv_point, static_cast<double>(totalCurrentIter-0), spatial_sampler, chrono_sampler);
                                    fetchTempY =
                                            ifs[iN]->fetch(fetchNameY, rcv_point, static_cast<double>(totalCurrentIter-0), spatial_sampler, chrono_sampler);
                                    fetchTempZ =
                                            ifs[iN]->fetch(fetchNameZ, rcv_point, static_cast<double>(totalCurrentIter-0), spatial_sampler, chrono_sampler);
                            } else
                            {
                                    FatalError << "cplFSIMethod: " << cplMethod << " in fsiDict is not recognized!" << exit(FatalError);
                            }




                            double angle2DRotationTemp = rotation2DAmplitude_*Foam::sin(2*3.1416*rotation2DFrq_*runTime.value());

                            dispVectorValsList[iN][pI][pointI].component(vector::X) = fetchTempX * Foam::cos(angle2DRotationTemp) - fetchTempY * Foam::sin(angle2DRotationTemp);
                            dispVectorValsList[iN][pI][pointI].component(vector::Y) = fetchTempX * Foam::sin(angle2DRotationTemp) + fetchTempY * Foam::cos(angle2DRotationTemp);
                            dispVectorValsList[iN][pI][pointI].component(vector::Z) = fetchTempZ;

                        }

                        if((dispVectorValsList[iN][pI][pointI].component(vector::X) != 0.0) && (debugFetchDisplacement == true))
                        {
                            Pout << "dispVectorValsList[iN][pI][pointI].component(vector::X): " << dispVectorValsList[iN][pI][pointI].component(vector::X) << nl << endl;
                        }
                        if((dispVectorValsList[iN][pI][pointI].component(vector::Y) != 0.0) && (debugFetchDisplacement == true))
                        {
                            Pout << "dispVectorValsList[iN][pI][pointI].component(vector::Y): " << dispVectorValsList[iN][pI][pointI].component(vector::Y) << nl << endl;
                        }
                        if((dispVectorValsList[iN][pI][pointI].component(vector::Z) != 0.0) && (debugFetchDisplacement == true))
                        {
                            Pout << "dispVectorValsList[iN][pI][pointI].component(vector::Z): " << dispVectorValsList[iN][pI][pointI].component(vector::Z) << nl << endl;
                        }
                    }
                }
            }
        }

        Info << "{OpenFOAM} : Coupling Relaxation Residuals at time  "<<runTime.value() << " and sub iter " << (totalCurrentIter-0) << " is " 
                << res_x << "  "<< res_y << "  "<< res_z  << endl;

    }
    else
    {
        forAll(interfaceNames, iN)
        {
            forAll(patch6DoFIDs[iN], pI)
            {
                // Loop over points to receive displacement values and assigned to dispVectorVals
                //- Blocking fetch for value passed through interface from coupled solver
                forAll(dispVectorValsList[iN][pI], pointI)
                {
                    rcv_point[0] = oriCoordFetchPointsVectorValsList[iN][pI][pointI].component(vector::X);
                    rcv_point[1] = oriCoordFetchPointsVectorValsList[iN][pI][pointI].component(vector::Y);
                    rcv_point[2] = oriCoordFetchPointsVectorValsList[iN][pI][pointI].component(vector::Z);

                    dispVectorValsList[iN][pI][pointI].component(vector::X) = 0.0;
                    dispVectorValsList[iN][pI][pointI].component(vector::Y) = 0.0;
                    dispVectorValsList[iN][pI][pointI].component(vector::Z) = 0.0;

                    if((dispVectorValsList[iN][pI][pointI].component(vector::X) != 0.0) && (debugFetchDisplacement == true))
                    {
                        Pout << "dispVectorValsList[iN][pI][pointI].component(vector::X): " << dispVectorValsList[iN][pI][pointI].component(vector::X) << nl << endl;
                    }
                    if((dispVectorValsList[iN][pI][pointI].component(vector::Y) != 0.0) && (debugFetchDisplacement == true))
                    {
                        Pout << "dispVectorValsList[iN][pI][pointI].component(vector::Y): " << dispVectorValsList[iN][pI][pointI].component(vector::Y) << nl << endl;
                    }
                    if((dispVectorValsList[iN][pI][pointI].component(vector::Z) != 0.0) && (debugFetchDisplacement == true))
                    {
                        Pout << "dispVectorValsList[iN][pI][pointI].component(vector::Z): " << dispVectorValsList[iN][pI][pointI].component(vector::Z) << nl << endl;
                    }
                }
            }
        }
    }

    forAll(interfaceNames, iN)
    {
        int pointCounter = 0;

        bool patchIDToFetch = false;

        forAll(patch6DoFIDs[iN], pI)
        {

            fixedValuePointPatchVectorField& PointPointerLocal =
                refCast<fixedValuePointPatchVectorField>
                (
                    pointDisplacement->boundaryFieldRef()[patch6DoFIDs[iN][pI]]
                );

            // Loop over points to assign dispVectorVals to cellDisplacement boundaryFieldRef, i.e. move the nodes
            forAll(dispVectorValsList[iN][pI], pointI)
            {
                ++pointCounter;
                mui::point3d rot_pointTemp;

                rot_pointTemp[0] = oriCoordFetchPointsVectorValsList[iN][pI][pointI].component(vector::X);
                rot_pointTemp[1] = oriCoordFetchPointsVectorValsList[iN][pI][pointI].component(vector::Y);
                rot_pointTemp[2] = oriCoordFetchPointsVectorValsList[iN][pI][pointI].component(vector::Z);

                double angle2DRotationTemp = rotation2DAmplitude_*Foam::sin(2*3.1416*rotation2DFrq_*runTime.value());

                PointPointerLocal[pointI].component(vector::X) = Foam::cos(angle2DRotationTemp) * (rot_pointTemp[0] - rotation2DOriginX_) - Foam::sin(angle2DRotationTemp) * (rot_pointTemp[1] - rotation2DOriginY_) + rotation2DOriginX_ - rot_pointTemp[0];
                PointPointerLocal[pointI].component(vector::Y) = Foam::sin(angle2DRotationTemp) * (rot_pointTemp[0] - rotation2DOriginX_) + Foam::cos(angle2DRotationTemp) * (rot_pointTemp[1] - rotation2DOriginY_) + rotation2DOriginY_ - rot_pointTemp[1];


                PointPointerLocal[pointI].component(vector::X) += rigidMotionAmplitudeX_*Foam::sin(2*3.1416*rigidMotionFrqX_*runTime.value());
                PointPointerLocal[pointI].component(vector::Y) += rigidMotionAmplitudeY_*Foam::sin(2*3.1416*rigidMotionFrqY_*runTime.value());
                PointPointerLocal[pointI].component(vector::Z) = rigidMotionAmplitudeZ_*Foam::sin(2*3.1416*rigidMotionFrqZ_*runTime.value());

                forAll(patchIDs[iN], ptI)
                {
                    if(patch6DoFIDs[iN][pI] == patchIDs[iN][ptI])
                    {
                        patchIDToFetch = true;
                        break;
                    }
                }

                if(patchIDToFetch == true)
                {
                    ++pointCounter;

                    if (applyFetchX)
                    {
                        PointPointerLocal[pointI].component(vector::X) += dispVectorValsList[iN][pI][pointI].component(vector::X);
                    }
                    if (applyFetchY)
                    {
                        PointPointerLocal[pointI].component(vector::Y) += dispVectorValsList[iN][pI][pointI].component(vector::Y);
                    }
                    if (applyFetchZ)
                    {
                        PointPointerLocal[pointI].component(vector::Z) += dispVectorValsList[iN][pI][pointI].component(vector::Z);
                    }

                    if((PointPointerLocal[pointI].component(vector::X) != 0.0) && (debugFetchDisplacement == true))
                    {
                        Pout << "PointPointerLocal[pointI].component(vector::X): " << PointPointerLocal[pointI].component(vector::X) << nl << endl;
                    }
                    if((PointPointerLocal[pointI].component(vector::Y) != 0.0) && (debugFetchDisplacement == true))
                    {
                        Pout << "PointPointerLocal[pointI].component(vector::Y): " << PointPointerLocal[pointI].component(vector::Y) << nl << endl;
                    }
                    if((PointPointerLocal[pointI].component(vector::Z) != 0.0) && (debugFetchDisplacement == true))
                    {
                        Pout << "PointPointerLocal[pointI].component(vector::Z): " << PointPointerLocal[pointI].component(vector::Z) << nl << endl;
                    }

                    if(Pstream::master())
                    {
                        if (applyFetchX || applyFetchY || applyFetchZ)
                        {
                            if (cplMethod == "Loose")
                            {
                                // Do nothing for Loose coupling method
                            }
                            else
                            {
                                muiResidual << interfaceNames[iN] << "\t" << runTime.timeName() << "\t" << timeSteps << "\t" << subIter << "\t"
                                            << "\t" << forceTotal[0].component(vector::Y) << "\t" << forceTotal[0].component(vector::Z) << endl;
                            }
                        }
                    }
                }
            }
            patchIDToFetch = false;
        }
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// *************************  FILE END  ************************************ //
