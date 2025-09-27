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

    // Debug switch for pushForce.H
    bool debugPushForce = false;

    // Obtain the normal vectors on the boundary faces
    const surfaceVectorField::Boundary& Sfb = mesh.Sf().boundaryField();

    // Obtain the mag normal vectors on the boundary faces
    const surfaceScalarField::Boundary& magSfb = mesh.magSf().boundaryField();

    // Obtain pressure fields, velocity fields, density and viscosity
    const volScalarField& pForce = mesh.lookupObject<volScalarField>("p");
    const volVectorField& UForce = mesh.lookupObject<volVectorField>("U");
    const volScalarField& rhoForce = mesh.lookupObject<volScalarField>("rho");
    const volScalarField& alphaWater = mesh.lookupObject<volScalarField>("alpha.water");
    const dimensionedScalar rhoForceWater(transportProperties.subDict("water").lookup("rho"));
    const dimensionedScalar rhoForceAir(transportProperties.subDict("air").lookup("rho"));
    const dimensionedScalar nuForceWater(transportProperties.subDict("water").lookup("nu"));
    const dimensionedScalar nuForceAir(transportProperties.subDict("air").lookup("nu"));

    Info << "{OpenFOAM} : rho Water: " << rhoForceWater << " nu Water: " << nuForceWater << nl << endl;
    Info << "{OpenFOAM} : rho Air: " << rhoForceAir << " nu Air: " << nuForceAir << nl << endl;

    //to avoid negative values
    const volScalarField limitedAlphaWater
    (
        "limitedAlphaWater",
        min(max(alphaWater, scalar(0)), scalar(1))
    );

    tmp<volSymmTensorField> tdevRhoReffWater = (-rhoForce*nuForceWater*dev(twoSymm(fvc::grad(UForce))));
    const volSymmTensorField::Boundary& devRhoReffbWater
        = tdevRhoReffWater().boundaryField();

    tmp<volSymmTensorField> tdevRhoReffAir = (-rhoForce*nuForceWater*dev(twoSymm(fvc::grad(UForce))));
    const volSymmTensorField::Boundary& devRhoReffbAir
        = tdevRhoReffAir().boundaryField();

    // Calculate the total forces on each face: Total_force = Pressure_force + Share_stress.
    // Note: incompressible solvers only.
    forAll(interfaceNames, iN)
    {
        forAll(patchIDs[iN], pI)
        {
            faceForces->boundaryFieldRef()[patchIDs[iN][pI]] =
                (
                    Sfb[patchIDs[iN][pI]] 
                    * 
                    pForce.boundaryField()[patchIDs[iN][pI]]

                )
                +
                (
                    (Sfb[patchIDs[iN][pI]]
                    )
                    & 
                    devRhoReffbWater[patchIDs[iN][pI]]
                );

            // Calculate the total traction on each face: Total_traction = Total_force / Face_area_magnitude.
            if (twoDimCorrection)
            {
                faceForces2D->boundaryFieldRef()[patchIDs[iN][pI]] =
                    faceForces->boundaryFieldRef()[patchIDs[iN][pI]]
                    /
                    twoDimCellHeight;

                faceTraction->boundaryFieldRef()[patchIDs[iN][pI]] = 
                    faceForces2D->boundaryFieldRef()[patchIDs[iN][pI]]
                    /
                    magSfb[patchIDs[iN][pI]];
            } else
            {
                faceTraction->boundaryFieldRef()[patchIDs[iN][pI]] = 
                    faceForces->boundaryFieldRef()[patchIDs[iN][pI]]
                    /
                    magSfb[patchIDs[iN][pI]];
            }
        }
    }

    List<vector> forceTotal;
    List<vector> force2DTotal;
    List<vector> tractionTotal;
    
    // Loop over all the interfaces
    forAll(interfaceNames, iN)
    {    
        vector forceTotalLocal = vector::zero;
        vector force2DTotalLocal = vector::zero;
        vector tractionTotalLocal = vector::zero;
        forceTotal.append(forceTotalLocal);
        force2DTotal.append(force2DTotalLocal);
        tractionTotal.append(tractionTotalLocal);
    }

    mui::point3d send_point;

    // Loop over all the interfaces
    forAll(interfaceNames, iN)
    {
        forAll(patchIDs[iN], pI)
        {
            forAll(oriCoordPushPointsVectorValsList[iN][pI], faceI)
            {
                vector facePush;

                if(pushForce)
                {
                    facePush = faceForces->boundaryFieldRef()[patchIDs[iN][pI]][faceI];
                } 
                else
                {
                    facePush = faceTraction->boundaryFieldRef()[patchIDs[iN][pI]][faceI];
                }

                //- Create 3D sending point at current face center
                // ToDo: change it to map  solution
                send_point[0] = static_cast<double>(static_cast<float>(oriCoordPushPointsVectorValsList[iN][pI][faceI].component(vector::X)));
                send_point[1] = static_cast<double>(static_cast<float>(oriCoordPushPointsVectorValsList[iN][pI][faceI].component(vector::Y)));
                send_point[2] = static_cast<double>(static_cast<float>(oriCoordPushPointsVectorValsList[iN][pI][faceI].component(vector::Z)));

                if (rankHasPatchPointsForPush)
                {
                    if (pushDelayIterations >= totalCurrentIter)
                    {
                        double tempZero = 0.0;

                        //- Push value of three dimensional tractions with label "tractionX/tractionY/tractionZ" at each face center
                        ifs[iN]->push(sendNameX, send_point, tempZero);
                        ifs[iN]->push(sendNameY, send_point, tempZero);
                        ifs[iN]->push(sendNameZ, send_point, tempZero);
                    } else
                    {

                        double angle2DRotationTemp = rotation2DAmplitude_*Foam::sin(2*3.1416*rotation2DFrq_*runTime.value());

                        double tempPushX = facePush.component(vector::X) * Foam::cos(angle2DRotationTemp) +  facePush.component(vector::Y) * Foam::sin(angle2DRotationTemp);
                        double tempPushY = facePush.component(vector::Y) * Foam::cos(angle2DRotationTemp) -  facePush.component(vector::X) * Foam::sin(angle2DRotationTemp);
                        double tempPushZ = facePush.component(vector::Z); 

                        //- Push value of three dimensional tractions with label "tractionX/tractionY/tractionZ" at each face center
                        ifs[iN]->push(sendNameX, send_point, tempPushX);
                        ifs[iN]->push(sendNameY, send_point, tempPushY);
                        ifs[iN]->push(sendNameZ, send_point, tempPushZ);
                    }
                }

                if (((facePush.component(vector::X) != 0.0) ||
                    (facePush.component(vector::Y) != 0.0) || 
                    (facePush.component(vector::Z) != 0.0)) &&
                    (debugPushForce == true))
                {
                    Pout << "send_point[0]= " << send_point[0] << nl <<endl;
                    Pout << "send_point[1]= " << send_point[1] << nl <<endl;
                    Pout << "send_point[2]= " << send_point[2] << nl <<endl;
                    if(pushForce)
                    {
                        Pout << "forceX= " << facePush.component(vector::X) << nl << endl;
                        Pout << "forceY= " << facePush.component(vector::Y) << nl << endl;
                        Pout << "forceZ= " << facePush.component(vector::Z) << nl << endl;
                    }
                    else
                    {
                        Pout << "tractionX= " << facePush.component(vector::X) << nl << endl;
                        Pout << "tractionY= " << facePush.component(vector::Y) << nl << endl;
                        Pout << "tractionZ= " << facePush.component(vector::Z) << nl << endl;
                    }
                }

                if (twoDimCorrection)
                {
                    force2DTotal[iN].component(vector::X) += faceForces2D->boundaryFieldRef()[patchIDs[iN][pI]][faceI].component(vector::X);
                    force2DTotal[iN].component(vector::Y) += faceForces2D->boundaryFieldRef()[patchIDs[iN][pI]][faceI].component(vector::Y);
                    force2DTotal[iN].component(vector::Z) += faceForces2D->boundaryFieldRef()[patchIDs[iN][pI]][faceI].component(vector::Z);
                }

                if(pushForce)
                {
                    tractionTotal[iN].component(vector::X) += faceTraction->boundaryFieldRef()[patchIDs[iN][pI]][faceI].component(vector::X);
                    tractionTotal[iN].component(vector::Y) += faceTraction->boundaryFieldRef()[patchIDs[iN][pI]][faceI].component(vector::Y);
                    tractionTotal[iN].component(vector::Z) += faceTraction->boundaryFieldRef()[patchIDs[iN][pI]][faceI].component(vector::Z);

                    forceTotal[iN].component(vector::X) += facePush.component(vector::X);
                    forceTotal[iN].component(vector::Y) += facePush.component(vector::Y);
                    forceTotal[iN].component(vector::Z) += facePush.component(vector::Z);
                }
                else
                {
                    forceTotal[iN].component(vector::X) += faceForces->boundaryFieldRef()[patchIDs[iN][pI]][faceI].component(vector::X);
                    forceTotal[iN].component(vector::Y) += faceForces->boundaryFieldRef()[patchIDs[iN][pI]][faceI].component(vector::Y);
                    forceTotal[iN].component(vector::Z) += faceForces->boundaryFieldRef()[patchIDs[iN][pI]][faceI].component(vector::Z);

                    tractionTotal[iN].component(vector::X) += facePush.component(vector::X);
                    tractionTotal[iN].component(vector::Y) += facePush.component(vector::Y);
                    tractionTotal[iN].component(vector::Z) += facePush.component(vector::Z);
                }
            }
        }

        //- Commit the pushed values to the interface at current time step
        ifs[iN]->commit(static_cast<scalar>(totalCurrentIter));

        if ((totalCurrentIter-2) > 0)
        {
            ifs[iN]->forget(static_cast<scalar>(0), static_cast<scalar>(totalCurrentIter-2));
            Info << "MUI interface " << iN << " forget from 0 to " << (totalCurrentIter-2) << " iterations" << nl << endl;
        } 

        Info << "{OpenFOMA}: MUI interface " << interfaceNames[iN] << " Finshed commiting forces at " << totalCurrentIter << " iterations" << nl << endl;
    }

    forAll(interfaceNames, iN)
    {
        // Integrate the forces of the patch
        reduce(forceTotal[iN], sumOp<vector>()); 
        reduce(tractionTotal[iN], sumOp<vector>());

        tractionTotal[iN] /= totalFaceNList[iN];

        Info << "{OpenFOAM} : " << interfaceNames[iN] << " forceX = " << forceTotal[iN].component(vector::X) << endl;
        Info << "{OpenFOAM} : " << interfaceNames[iN] << " forceY = " << forceTotal[iN].component(vector::Y) << endl;
        Info << "{OpenFOAM} : " << interfaceNames[iN] << " forceZ = " << forceTotal[iN].component(vector::Z) << endl;

        if (twoDimCorrection)
        {
            reduce(force2DTotal[iN], sumOp<vector>()); 
            Info << "{OpenFOAM} : " << interfaceNames[iN] << " force2DX = " << force2DTotal[iN].component(vector::X) << endl;
            Info << "{OpenFOAM} : " << interfaceNames[iN] << " force2DY = " << force2DTotal[iN].component(vector::Y) << endl;
            Info << "{OpenFOAM} : " << interfaceNames[iN] << " force2DZ = " << force2DTotal[iN].component(vector::Z) << endl;
        }

        Info << "{OpenFOAM} : " << interfaceNames[iN] << " tractionX = " << tractionTotal[iN].component(vector::X) << endl;
        Info << "{OpenFOAM} : " << interfaceNames[iN] << " tractionY = " << tractionTotal[iN].component(vector::Y) << endl;
        Info << "{OpenFOAM} : " << interfaceNames[iN] << " tractionZ = " << tractionTotal[iN].component(vector::Z) << endl;
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// *************************  FILE END  ************************************ //
