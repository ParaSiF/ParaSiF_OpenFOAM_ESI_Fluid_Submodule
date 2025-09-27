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

    mui::point3d rcv_point;

    scalar t = mesh.time().value();
 
    Info << "{OpenFOAM} Fetch iteration " << (totalCurrentIter - 0) << nl << endl;
    if((totalCurrentIter-0) >= 0)
    {                    
        forAll(interfaceNames, iN)
        {
            int pointCounter = 0;

            forAll(patchIDs[iN], pI)
            {
                // Loop over points to receive displacement values and assigned to dispVectorVals
                //- Blocking fetch for value passed through interface from coupled solver  

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
                        } else if (cplMethod == "FixedRelaxation")
                        {
                            fetchTempX =
                                    ifs[iN]->fetch(fetchNameX, rcv_point, static_cast<double>(totalCurrentIter-0), spatial_sampler, chrono_sampler,fr_x);
                            fetchTempY =
                                    ifs[iN]->fetch(fetchNameY, rcv_point, static_cast<double>(totalCurrentIter-0), spatial_sampler, chrono_sampler,fr_y);
                            fetchTempZ =
                                    ifs[iN]->fetch(fetchNameZ, rcv_point, static_cast<double>(totalCurrentIter-0), spatial_sampler, chrono_sampler,fr_z);
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

                            dispVectorValsList[iN][pI][pointI].component(vector::X) = fetchTempX;
                            dispVectorValsList[iN][pI][pointI].component(vector::Y) = fetchTempY;
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
    else
    {
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

        forAll(patchIDs[iN], pI)
        {

            fixedValuePointPatchVectorField& PointPointerLocal = 
                refCast<fixedValuePointPatchVectorField>
                (
                    pointDisplacement->boundaryFieldRef()[patchIDs[iN][pI]]
                );

            // Loop over points to assign dispVectorVals to cellDisplacement boundaryFieldRef, i.e. move the nodes
            forAll(dispVectorValsList[iN][pI], pointI)
            {
                ++pointCounter;

                PointPointerLocal[pointI].component(vector::X) = rigidMotionAmplitudeX_*Foam::sin(2*3.1416*rigidMotionFrqX_*t);
                PointPointerLocal[pointI].component(vector::Y) = rigidMotionAmplitudeY_*Foam::sin(2*3.1416*rigidMotionFrqY_*t);
                PointPointerLocal[pointI].component(vector::Z) = rigidMotionAmplitudeZ_*Foam::sin(2*3.1416*rigidMotionFrqZ_*t);

                PointPointerLocal[pointI].component(vector::X) += dispVectorValsList[iN][pI][pointI].component(vector::X);
                PointPointerLocal[pointI].component(vector::Y) += dispVectorValsList[iN][pI][pointI].component(vector::Y);
                PointPointerLocal[pointI].component(vector::Z) += dispVectorValsList[iN][pI][pointI].component(vector::Z);


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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// *************************  FILE END  ************************************ //
