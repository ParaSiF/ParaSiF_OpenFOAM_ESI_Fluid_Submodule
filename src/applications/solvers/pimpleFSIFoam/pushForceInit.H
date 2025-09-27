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

    // Time Steps counter
    int timeSteps = 0;

    // Total iteration counter
    int totalCurrentIter = 0;

    // Define the name of MUI domain
    word domainName;

    // Define the name of MUI interfaceName
    List<word> interfaceNames;
    List<List<word>> patchNames;

    int subIterationNum;
    int forgetTimeSteps;

    word sendNameX;
    word sendNameY;
    word sendNameZ;

    word fetchNameX;
    word fetchNameY;
    word fetchNameZ;

    double rSampler;
 
    //- Tight/Loose coupling method
    word cplMethod;

    //- Initial under relaxation factor of the coupling
    double initUndRelxCpl;

    //- Whether to push force or traction
    bool pushForce;

    //- Whether apply the fetched X value or not(i.e. one-way coupling)
    bool applyFetchX;
    //- Whether apply the fetched Y value or not(i.e. one-way coupling)
    bool applyFetchY;
    //- Whether apply the fetched Z value or not(i.e. one-way coupling)
    bool applyFetchZ;

    //- Maximum value of the under relaxation number (valid only for IQN-ILS)
    double maxUndRelxCpl;

    //- Number of Aitken iterations to at the beginning of each IQN-ILS time-steps (valid only for IQN-ILS)
    int numAitkenIterations;

    //- Number of time steps that will not push forces/tractions through MUI.
    int pushDelay;

    //- Whether it is a 2D case, i.e. whether to have a 2D force correction before push? 
    //- A 2D case is the case that the front and back boundaries are empty boundary condition. 
    //- There usually has one cell between the front and back boundaries.
    bool twoDimCorrection;

    //- The cell height of the only cell between the two empty boundaries [m].
    //- Valid only when twoDimCorrection is true.
    double twoDimCellHeight;

    scalar rigidMotionFrqX_;
    scalar rigidMotionFrqY_;
    scalar rigidMotionFrqZ_;

    scalar rigidMotionAmplitudeX_;
    scalar rigidMotionAmplitudeY_;
    scalar rigidMotionAmplitudeZ_;

    //- Whether to change the sub-iteration number
    bool changeSubIter;

    //- When to change the sub-iteration number [s].
    double changeSubIterTime;

    //- New sub-iteration number.
    int subIterationNumNew;  

    string fsiFile("fsiDict");

    IOdictionary fsiDict
    (
        IOobject
        (
            fsiFile,
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    //Read fsi dictionary file if it exists
    if (fsiDict.headerOk())
    {

        word domainFSIName;
        if (!(fsiDict.readIfPresent("domainFSIName", domainFSIName)))
        {
            FatalIOErrorIn("", fsiDict)
                           << "Missing domainFSIName entry" << exit(FatalIOError);
        }
        domainName = domainFSIName;

        const PtrList<entry> interfaceConfigurations
        (
            fsiDict.lookup("interfaceLists")
        );

        //Iterate through interface lists
        forAll(interfaceConfigurations, cC)
        {
            const entry& interfaceI = interfaceConfigurations[cC];
            const dictionary& interfaceIDict = interfaceI.dict();

            interfaceNames.append(static_cast<word>(interfaceIDict.dictName()));
            
            List<word> patchNamesLocal(interfaceIDict.lookup("PatchNames"));

            patchNames.append(patchNamesLocal);

        }

        word cplFSIMethod;
        if (!(fsiDict.readIfPresent("cplFSIMethod", cplFSIMethod)))
        {
            FatalIOErrorIn("", fsiDict)
                           << "Missing cplFSIMethod entry" << exit(FatalIOError);
        }
        cplMethod = cplFSIMethod;

        int subIterationNumber;
        if (!(fsiDict.readIfPresent("subIterationNumber", subIterationNumber)))
        {
            FatalIOErrorIn("", fsiDict)
                           << "Missing subIterationNumber entry" << exit(FatalIOError);
        }
        subIterationNum = subIterationNumber;

        int forgetFSITimeSteps;
        if (!(fsiDict.readIfPresent("forgetFSITimeSteps", forgetFSITimeSteps)))
        {
            FatalIOErrorIn("", fsiDict)
                           << "Missing forgetFSITimeSteps entry" << exit(FatalIOError);
        }
        forgetTimeSteps = forgetFSITimeSteps;

        word sendFSINameX;
        word sendFSINameY;
        word sendFSINameZ;
        word fetchFSINameX;
        word fetchFSINameY;
        word fetchFSINameZ;

        if (!(fsiDict.readIfPresent("sendFSINameX", sendFSINameX)))
        {
            FatalIOErrorIn("", fsiDict)
                           << "Missing sendFSINameX entry" << exit(FatalIOError);
        }
        if (!(fsiDict.readIfPresent("sendFSINameY", sendFSINameY)))
        {
            FatalIOErrorIn("", fsiDict)
                           << "Missing sendFSINameY entry" << exit(FatalIOError);
        }
        if (!(fsiDict.readIfPresent("sendFSINameZ", sendFSINameZ)))
        {
            FatalIOErrorIn("", fsiDict)
                           << "Missing sendFSINameZ entry" << exit(FatalIOError);
        }
        if (!(fsiDict.readIfPresent("fetchFSINameX", fetchFSINameX)))
        {
            FatalIOErrorIn("", fsiDict)
                           << "Missing fetchFSINameX entry" << exit(FatalIOError);
        }
        if (!(fsiDict.readIfPresent("fetchFSINameY", fetchFSINameY)))
        {
            FatalIOErrorIn("", fsiDict)
                           << "Missing fetchFSINameY entry" << exit(FatalIOError);
        }
        if (!(fsiDict.readIfPresent("fetchFSINameZ", fetchFSINameZ)))
        {
            FatalIOErrorIn("", fsiDict)
                           << "Missing fetchFSINameZ entry" << exit(FatalIOError);
        }
        sendNameX = sendFSINameX;
        sendNameY = sendFSINameY;
        sendNameZ = sendFSINameZ;
        fetchNameX = fetchFSINameX;
        fetchNameY = fetchFSINameY;
        fetchNameZ = fetchFSINameZ;

        bool pushFSIForce;
        if (!(fsiDict.readIfPresent("pushFSIForce", pushFSIForce)))
        {
            FatalIOErrorIn("", fsiDict)
                           << "Missing pushFSIForce entry" << exit(FatalIOError);
        }
        pushForce = pushFSIForce;

        bool applyFSIFetchX;
        if (!(fsiDict.readIfPresent("applyFSIFetchX", applyFSIFetchX)))
        {
            FatalIOErrorIn("", fsiDict)
                           << "Missing applyFSIFetchX entry" << exit(FatalIOError);
        }
        applyFetchX = applyFSIFetchX;

        bool applyFSIFetchY;
        if (!(fsiDict.readIfPresent("applyFSIFetchY", applyFSIFetchY)))
        {
            FatalIOErrorIn("", fsiDict)
                           << "Missing applyFSIFetchY entry" << exit(FatalIOError);
        }
        applyFetchY = applyFSIFetchY;

        bool applyFSIFetchZ;
        if (!(fsiDict.readIfPresent("applyFSIFetchZ", applyFSIFetchZ)))
        {
            FatalIOErrorIn("", fsiDict)
                           << "Missing applyFSIFetchZ entry" << exit(FatalIOError);
        }
        applyFetchZ = applyFSIFetchZ;

        double searchAreaSampler;
        if (!(fsiDict.readIfPresent("searchAreaSampler", searchAreaSampler)))
        {
            FatalIOErrorIn("", fsiDict)
                           << "Missing searchAreaSampler entry" << exit(FatalIOError);
        }
        rSampler = searchAreaSampler;

        if (cplMethod != "Loose")
        {
            double initFSIUndRelxCpl;
            if (!(fsiDict.readIfPresent("initFSIUndRelxCpl", initFSIUndRelxCpl)))
            {
                FatalIOErrorIn("", fsiDict)
                               << "Missing initFSIUndRelxCpl entry" << exit(FatalIOError);
            }
            initUndRelxCpl = initFSIUndRelxCpl;
        }

        if (cplMethod == "IQN_ILS")
        {
            double maxFSIUndRelxCpl;
            if (!(fsiDict.readIfPresent("maxFSIUndRelxCpl", maxFSIUndRelxCpl)))
            {
                FatalIOErrorIn("", fsiDict)
                               << "Missing maxFSIUndRelxCpl entry" << exit(FatalIOError);
            }
            maxUndRelxCpl = maxFSIUndRelxCpl;

            int numFSIAitkenIterations;
            if (!(fsiDict.readIfPresent("numFSIAitkenIterations", numFSIAitkenIterations)))
            {
                FatalIOErrorIn("", fsiDict)
                               << "Missing numFSIAitkenIterations entry" << exit(FatalIOError);
            }
            numAitkenIterations = numFSIAitkenIterations;
        }

        int pushFSIDelay;
        if (!(fsiDict.readIfPresent("pushFSIDelay", pushFSIDelay)))
        {
            FatalIOErrorIn("", fsiDict)
                           << "Missing pushFSIDelay entry" << exit(FatalIOError);
        }
        pushDelay = pushFSIDelay;

        bool twoDimFSICorrection;
        if (!(fsiDict.readIfPresent("twoDimFSICorrection", twoDimFSICorrection)))
        {
            FatalIOErrorIn("", fsiDict)
                           << "Missing twoDimFSICorrection entry" << exit(FatalIOError);
        }
        twoDimCorrection = twoDimFSICorrection;

        if (twoDimCorrection)
        {
            double twoDimFSICellHeight;
            if (!(fsiDict.readIfPresent("twoDimFSICellHeight", twoDimFSICellHeight)))
            {
                FatalIOErrorIn("", fsiDict)
                               << "Missing twoDimFSICellHeight entry" << exit(FatalIOError);
            }
            twoDimCellHeight = twoDimFSICellHeight;
        }
        
        double rigidMotionFrqX;
        if (!(fsiDict.readIfPresent("rigidMotionFrqX", rigidMotionFrqX)))
        {
            FatalIOErrorIn("", fsiDict)
                           << "Missing rigidMotionFrqX entry" << exit(FatalIOError);
        }
        rigidMotionFrqX_ = rigidMotionFrqX;
        
        double rigidMotionFrqY;
        if (!(fsiDict.readIfPresent("rigidMotionFrqY", rigidMotionFrqY)))
        {
            FatalIOErrorIn("", fsiDict)
                           << "Missing rigidMotionFrqY entry" << exit(FatalIOError);
        }
        rigidMotionFrqY_ = rigidMotionFrqY;
        
        double rigidMotionFrqZ;
        if (!(fsiDict.readIfPresent("rigidMotionFrqZ", rigidMotionFrqZ)))
        {
            FatalIOErrorIn("", fsiDict)
                           << "Missing rigidMotionFrqZ entry" << exit(FatalIOError);
        }
        rigidMotionFrqZ_ = rigidMotionFrqZ;
        
        double rigidMotionAmplitudeX;
        if (!(fsiDict.readIfPresent("rigidMotionAmplitudeX", rigidMotionAmplitudeX)))
        {
            FatalIOErrorIn("", fsiDict)
                           << "Missing rigidMotionAmplitudeX entry" << exit(FatalIOError);
        }
        rigidMotionAmplitudeX_ = rigidMotionAmplitudeX;
        
        double rigidMotionAmplitudeY;
        if (!(fsiDict.readIfPresent("rigidMotionAmplitudeY", rigidMotionAmplitudeY)))
        {
            FatalIOErrorIn("", fsiDict)
                           << "Missing rigidMotionAmplitudeY entry" << exit(FatalIOError);
        }
        rigidMotionAmplitudeY_ = rigidMotionAmplitudeY;
        
        double rigidMotionAmplitudeZ;
        if (!(fsiDict.readIfPresent("rigidMotionAmplitudeZ", rigidMotionAmplitudeZ)))
        {
            FatalIOErrorIn("", fsiDict)
                           << "Missing rigidMotionAmplitudeZ entry" << exit(FatalIOError);
        }
        rigidMotionAmplitudeZ_ = rigidMotionAmplitudeZ;

        bool changeFSISubIter;
        if (!(fsiDict.readIfPresent("changeFSISubIter", changeFSISubIter)))
        {
            FatalIOErrorIn("", fsiDict)
                           << "Missing changeFSISubIter entry" << exit(FatalIOError);
        }
        changeSubIter = changeFSISubIter;
        
        double changeFSISubIterTime;
        if (!(fsiDict.readIfPresent("changeFSISubIterTime", changeFSISubIterTime)))
        {
            FatalIOErrorIn("", fsiDict)
                           << "Missing changeFSISubIterTime entry" << exit(FatalIOError);
        }
        changeSubIterTime = changeFSISubIterTime;

        int newFSISubIterationNum;
        if (!(fsiDict.readIfPresent("newFSISubIterationNum", newFSISubIterationNum)))
        {
            FatalIOErrorIn("", fsiDict)
                           << "Missing newFSISubIterationNum entry" << exit(FatalIOError);
        }
        subIterationNumNew = newFSISubIterationNum;

    }
    else
    {
            FatalIOErrorIn("", fsiDict)
                           << "Missing fsiDict" << exit(FatalIOError);
    }

    List<List<int>> patchIDs;

    forAll(interfaceNames, iN)
    {

        List<int> localPatchIDs;

        forAll(patchNames[iN], pN)
        {
            int localPatchID = mesh.boundaryMesh().findPatchID(patchNames[iN][pN]);

        Info << "{OpenFOMA} : PatchName = " << patchNames[iN][pN] << " has ID = " << localPatchID << endl;

            // Check if patch exists
            if(localPatchID == -1)
            {
            FatalError << "{OpenFOMA} : PatchName " << patchNames[iN][pN] << " not found!" << exit(FatalError);
            }

            localPatchIDs.append(localPatchID);

        }

        patchIDs.append(localPatchIDs);

    }
    
    volVectorField * faceForces = new volVectorField
    (
        IOobject
        (
            "faceForces",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "faceForces",
            dimensionSet(1,1,-2,0,0,0,0),
            vector::zero          
        )
    );

    volVectorField * faceForces2D = new volVectorField
    (
        IOobject
        (
            "faceForces2D",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "faceForces2D",
            dimensionSet(1,1,-2,0,0,0,0),
            vector::zero          
        )
    );

    volVectorField * faceTraction = new volVectorField
    (
        IOobject
        (
            "faceTraction",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "faceTraction",
            dimensionSet(1,-1,-2,0,0,0,0),
            vector::zero          
        )
    );

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    bool rankHasPatchPointsForPush = false;

    List<scalar> totalFaceNList;
    // Declare a vectorField to store the original point coordinates.
    List<List<vectorField>> oriCoordPushPointsVectorValsList;

    forAll(interfaceNames, iN)
    {

        scalar totalFaceNLocal;
        List<vectorField> oriCoordPushPointsVectorValsListLocal;

        forAll(patchIDs[iN], pI)
        {
            if(mesh.boundaryMesh()[patchIDs[iN][pI]].size() != 0)
            {
                rankHasPatchPointsForPush = true;
            }

            totalFaceNLocal += mesh.boundaryMesh()[patchIDs[iN][pI]].size();
            reduce(totalFaceNLocal, sumOp<scalar>());
            vectorField oriCoordPushPointsVectorValsLocal(mesh.Cf().boundaryField()[patchIDs[iN][pI]].size());
            oriCoordPushPointsVectorValsListLocal.append(oriCoordPushPointsVectorValsLocal);
        }
        totalFaceNList.append(totalFaceNLocal);
        oriCoordPushPointsVectorValsList.append(oriCoordPushPointsVectorValsListLocal);
    }

    // Define the name of MUI interfaces
    std::vector<std::string> interfaces;

    forAll(interfaceNames, iN)
    {
        interfaces.emplace_back(interfaceNames[iN]);
    }

    // Declare MUI objects using MUI configure file
    auto ifs = mui::create_uniface<mui::mui_config>( domainName, interfaces );

    //- If this is not a parallel run then need to first call MPI_Init for MUI (otherwise this is called during PStream creation)
    if (!args.parRunControl().parRun())
    {
        MPI_Comm local_world = mui::mpi_split_by_app(argc, argv);
    }

    scalar sendXMin = VGREAT;
    scalar sendYMin = VGREAT;
    scalar sendZMin = VGREAT;
    
    scalar sendXMax = -VGREAT;
    scalar sendYMax = -VGREAT;
    scalar sendZMax = -VGREAT;

    int const pushDelayIterations = pushDelay * subIterationNum;

    // Loop over faces to Store the original face center coordinates
    forAll(interfaceNames, iN)
    {
        forAll(patchIDs[iN], pI)
        {
            // ToDo: use map function to link face ID with coordinates
            forAll(oriCoordPushPointsVectorValsList[iN][pI], faceI)
            {
               oriCoordPushPointsVectorValsList[iN][pI][faceI].component(vector::X) =
                    mesh.Cf().boundaryField()[patchIDs[iN][pI]][faceI].component(vector::X);
                oriCoordPushPointsVectorValsList[iN][pI][faceI].component(vector::Y) =
                    mesh.Cf().boundaryField()[patchIDs[iN][pI]][faceI].component(vector::Y);
                oriCoordPushPointsVectorValsList[iN][pI][faceI].component(vector::Z) =
                    mesh.Cf().boundaryField()[patchIDs[iN][pI]][faceI].component(vector::Z);
            }
        }
    }

    if (rankHasPatchPointsForPush)
    {

        forAll(interfaceNames, iN)
        {

            forAll(patchIDs[iN], pI)
            {
                forAll(oriCoordPushPointsVectorValsList[iN][pI], faceI)
                {
                    sendXMin =
                        ( (oriCoordPushPointsVectorValsList[iN][pI][faceI].component(vector::X)) < sendXMin) ? (oriCoordPushPointsVectorValsList[iN][pI][faceI].component(vector::X)) : sendXMin;
                    sendYMin =
                        ((oriCoordPushPointsVectorValsList[iN][pI][faceI].component(vector::Y)) < sendYMin) ? (oriCoordPushPointsVectorValsList[iN][pI][faceI].component(vector::Y)) : sendYMin;
                    sendZMin =
                        ((oriCoordPushPointsVectorValsList[iN][pI][faceI].component(vector::Z)) < sendZMin) ? (oriCoordPushPointsVectorValsList[iN][pI][faceI].component(vector::Z)) : sendZMin;

                    sendXMax =
                        ((oriCoordPushPointsVectorValsList[iN][pI][faceI].component(vector::X)) > sendXMax) ? (oriCoordPushPointsVectorValsList[iN][pI][faceI].component(vector::X)) : sendXMax;
                    sendYMax =
                        ((oriCoordPushPointsVectorValsList[iN][pI][faceI].component(vector::Y)) > sendYMax) ? (oriCoordPushPointsVectorValsList[iN][pI][faceI].component(vector::Y)) : sendYMax;
                    sendZMax =
                        ((oriCoordPushPointsVectorValsList[iN][pI][faceI].component(vector::Z)) > sendZMax) ? (oriCoordPushPointsVectorValsList[iN][pI][faceI].component(vector::Z)) : sendZMax;
                }
            }

        }
    } else
    {

        const pointField& meshPoints = mesh.points();

        forAll(meshPoints, pts)
        {
            sendXMin =
                (meshPoints[pts][0] < sendXMin) ? meshPoints[pts][0] : sendXMin;
            sendYMin =
                (meshPoints[pts][1] < sendYMin) ? meshPoints[pts][1] : sendYMin;
            sendZMin =
                (meshPoints[pts][2] < sendZMin) ? meshPoints[pts][2] : sendZMin;

            sendXMax =
                (meshPoints[pts][0] > sendXMax) ? meshPoints[pts][0] : sendXMax;
            sendYMax =
                (meshPoints[pts][1] > sendYMax) ? meshPoints[pts][1] : sendYMax;
            sendZMax =
                (meshPoints[pts][2] > sendZMax) ? meshPoints[pts][2] : sendZMax;
        }
    }


    Pout << "sendXMin = " << sendXMin << endl;
    Pout << "sendYMin = " << sendYMin << endl;
    Pout << "sendZMin = " << sendZMin << endl;

    Pout << "sendXMax = " << sendXMax << endl;
    Pout << "sendYMax = " << sendYMax << endl;
    Pout << "sendZMax = " << sendZMax << endl;

    // Annouce send span
    mui::geometry::box<mui::mui_config> send_region( {sendXMin, sendYMin, sendZMin}, {sendXMax, sendYMax, sendZMax} );
    // Set the maximum time step size
    int totalTimeSteps = ceil((runTime.endTime().value() - runTime.startTime().value())/runTime.deltaT().value());

    forAll(interfaceNames, iN)
    {
        ifs[iN]->announce_send_span( 0, totalTimeSteps, send_region );
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// *************************  FILE END  ************************************ //