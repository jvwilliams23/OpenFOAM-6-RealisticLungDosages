/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "dynamicLagrangianDev.H"
#include "fvOptions.H"
#include "fvcGrad.H"
#include "fvm.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void dynamicLagrangianDev<BasicTurbulenceModel>::correctNut
(
    const tmp<volTensorField>& gradU
)
{
    this->nut_ = (flm_/fmm_)*sqr(this->delta())*mag(dev(symm(gradU)));
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
void dynamicLagrangianDev<BasicTurbulenceModel>::correctNut()
{
    correctNut(fvc::grad(this->U_));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
dynamicLagrangianDev<BasicTurbulenceModel>::dynamicLagrangianDev
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    LESeddyViscosity<BasicTurbulenceModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    flm_
    (
        IOobject
        (
            IOobject::groupName("flm", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("flm", pow(dimLength, 4.0)/pow(dimTime, 4.0), VSMALL),
        "zeroGradient"
    ),
    fmm_
    (
        IOobject
        (
            IOobject::groupName("fmm", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("fmm", pow(dimLength, 4.0)/pow(dimTime, 4.0), VSMALL),
        "zeroGradient"
    ),
    theta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "theta",
            this->coeffDict_,
            1.5
        )
    ),
    k_
    (
        IOobject
        (
            IOobject::groupName("k", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("epsilon", sqr(dimLength)/sqr(dimTime), VSMALL)
    ),
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", this->U_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("epsilon", sqr(dimLength)/sqr(dimTime)/dimTime, VSMALL)
    ),
    Ce_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ce",
            this->coeffDict_,
            1.0
        )
    ),
    Cs_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cs",
            this->coeffDict_,
            0.065
        )
    ),
    simpleFilter_(U.mesh()),
    filterPtr_(LESfilter::New(U.mesh(), this->coeffDict())),
    filter_(filterPtr_()),

    flm0_("flm0", flm_.dimensions(), 0.0),
    fmm0_("fmm0", fmm_.dimensions(), vSmall)
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool dynamicLagrangianDev<BasicTurbulenceModel>::read()
{
    if (LESeddyViscosity<BasicTurbulenceModel>::read())
    {
        filter_.read(this->coeffDict());
        theta_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void dynamicLagrangianDev<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    LESeddyViscosity<BasicTurbulenceModel>::correct();

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();

    volSymmTensorField S(dev(symm(gradU)));
    volScalarField magS(mag(S));

    volVectorField Uf(filter_(U));
    volSymmTensorField Sf(dev(symm(fvc::grad(Uf))));
    volScalarField magSf(mag(Sf));

    volSymmTensorField L(dev(filter_(sqr(U)) - (sqr(filter_(U)))));
    volSymmTensorField M
    (
        2.0*sqr(this->delta())*(filter_(magS*S) - 4.0*magSf*Sf)
    );

    volScalarField invT
    (
        alpha*rho*(1.0/(theta_.value()*this->delta()))*pow(flm_*fmm_, 1.0/8.0)
    );

    volScalarField LM(L && M);

    fvScalarMatrix flmEqn
    (
        fvm::ddt(alpha, rho, flm_)
      + fvm::div(alphaRhoPhi, flm_)
     ==
        invT*LM
      - fvm::Sp(invT, flm_)
      + fvOptions(alpha, rho, flm_)
    );

    flmEqn.relax();
    fvOptions.constrain(flmEqn);
    flmEqn.solve();
    fvOptions.correct(flm_);
    bound(flm_, flm0_);

    volScalarField MM(M && M);

    fvScalarMatrix fmmEqn
    (
        fvm::ddt(alpha, rho, fmm_)
      + fvm::div(alphaRhoPhi, fmm_)
     ==
        invT*MM
      - fvm::Sp(invT, fmm_)
      + fvOptions(alpha, rho, fmm_)
    );

    fmmEqn.relax();
    fvOptions.constrain(fmmEqn);
    fmmEqn.solve();
    fvOptions.correct(fmm_);
    bound(fmm_, fmm0_);

    // JW 05/04/21
    //epsilon_ = sqr(Cs_*this->delta())*pow(magS,3.0);
    //k_ = k();
    k_ = pow(2.0*flm_/fmm_, 2.0/3.0)
              * pow(this->Ce_, -2.0/3.0)
              * sqr(this->delta())*magSqr(dev(symm(gradU)));

    epsilon_ = Ce_*pow(k_, 1.5)/this->delta();
    //k_ = Ce_*pow(this->delta()*epsilon_, 2.0/3.0);

    correctNut(gradU);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
