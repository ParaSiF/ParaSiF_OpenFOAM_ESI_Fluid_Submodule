{
    if (correctPhi)
    {
        rAU.ref() = 1.0/UEqn.A();
    }
    else
    {
        rAU = 1.0/UEqn.A();
    }

    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));
    volVectorField HbyA(constrainHbyA(rAU()*UEqn.H(), U, p_rgh));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + MRF.zeroFilter(fvc::interpolate(rho*rAU())*fvc::ddtCorr(U, phi, Uf))
    );
    MRF.makeRelative(phiHbyA);

    if (p_rgh.needReference())
    {
        fvc::makeRelative(phiHbyA, U);
        adjustPhi(phiHbyA, U, p_rgh);
        fvc::makeAbsolute(phiHbyA, U);
    }

    surfaceScalarField phig
    (
        (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

    while (pimple.correctNonOrthogonal())
    {

    // * * * * * * * * * * * MUI Initialise * * * * * * * * * * * * * * * * * //

        Info << "{OpenFOAM} non-orthogonal iteration: " << pimple.corrNonOrtho() << " with nNonOrthCorr " << pimple.nNonOrthCorr() << nl;

        #include "clockTime.H"

        // Start timer
        Foam::clockTime loopTimer;

        scalar pressureDoubleFetch = 0.0;
        scalar pressureSend = 0.0;

        label totalIter = (runTime.timeIndex() - 1) * pimple.nCorrPIMPLE() * pimple.nCorrPISO() * (pimple.nNonOrthCorr()+1)
            + (pimple.corr() - 1) * pimple.nCorrPISO() * (pimple.nNonOrthCorr()+1)
            + (pimple.corrPISO() - 1) * (pimple.nNonOrthCorr()+1)
            + pimple.corrNonOrtho();

        const Foam::volScalarField& pField = p_rgh;

        forAll(pField, cellI)
        {
            const scalar pVal = pField[cellI];
            const vector& coord = centres[cellI];

            mui::point3d locp( coord.x(), coord.y(), coord.z() );

            ifs[0]->push( "prgh", locp, pVal );

            pressureSend = pVal;
        }

        Pout << "{OpenFOAM} Sending - p_rgh: " << pressureSend << nl;

        // commit step
        ifs[0]->commit(totalIter);

        Info << "{OpenFOAM} commit at " << totalIter << nl;

        forAll(pField, cellI)
        {
            const vector& coord = centres[cellI];
            mui::point3d locf( coord.x(), coord.y(), coord.z() );

            pressureDoubleFetch = ifs[0]->fetch( nameFetch, locf, totalIter, s1, s2 );

        }

        Pout << "{OpenFOAM} Receiving - p_rgh Double: " << pressureDoubleFetch << nl;

        // Print elapsed time
        Info << "{OpenFOAM} Two-way Data passing in: " 
            << loopTimer.elapsedTime() << " s" << nl;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(p_rgh.select(pimple.finalInnerIter()));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p_rghEqn.flux();

            p_rgh.relax();

            U = HbyA + rAU()*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    #include "continuityErrs.H"

    // Correct Uf if the mesh is moving
    fvc::correctUf(Uf, U, phi);

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }

    if (!correctPhi)
    {
        rAU.clear();
    }
}
