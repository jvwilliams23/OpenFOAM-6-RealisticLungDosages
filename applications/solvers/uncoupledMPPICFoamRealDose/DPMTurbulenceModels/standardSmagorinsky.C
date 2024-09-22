/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "standardSmagorinsky.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void standardSmagorinsky<BasicTurbulenceModel>::correctNut
(
    const tmp<volTensorField>& gradU
)
{
    this->nut_ = max(Cs_*sqr(this->delta())*mag(dev(symm(gradU))),-1.0*this->nu());
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
void standardSmagorinsky<BasicTurbulenceModel>::correctNut()
{
    correctNut(fvc::grad(this->U_));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
standardSmagorinsky<BasicTurbulenceModel>::standardSmagorinsky
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

    k_
    (
        IOobject
        (
            IOobject::groupName("k", this->U_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
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

    Cs_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cs",
            this->coeffDict_,
            0.065
        )
    ),

    Ce_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ce",
            this->coeffDict_,
            1.0
        )
    )

{
//    bound(k_, this->kMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool standardSmagorinsky<BasicTurbulenceModel>::read()
{
    if (LESeddyViscosity<BasicTurbulenceModel>::read())
    {

        return true;
    }
    else
    {
        return false;
    }
}

/*template<class BasicTurbulenceModel>
tmp<volScalarField> standardSmagorinsky<BasicTurbulenceModel>::epsilon() const
{
     volScalarField k = this->k_; //(this->k(fvc::grad(this->U_)));
 
     return tmp<volScalarField>
     (
         new volScalarField
         (
             IOobject
             (
                 IOobject::groupName("epsilon", this->alphaRhoPhi_.group()),
                 this->runTime_.timeName(),
                 this->mesh_,
                 IOobject::NO_READ,
                 IOobject::AUTO_WRITE //NO_WRITE
             ),
             this->Ce_*pow(mag(k),1.5)/this->delta()
         )
     );
}*/

template<class BasicTurbulenceModel>
void standardSmagorinsky<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    //const surfaceScalarField& phi = this->phi_;
    const volVectorField& U = this->U_;
    //fv::options& fvOptions(fv::options::New(this->mesh_));

    LESeddyViscosity<BasicTurbulenceModel>::correct();

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();

    //const volTensorField tau = this->mu() * (gradU + gradU.T());
    //k_ = 0.5 * tr(tau);

    volSymmTensorField S(dev(symm(gradU)));
    volScalarField magS(mag(S));
    //volTensorField S( 0.5*(gradU+gradU) );
    //volScalarField magS(

    epsilon_ = sqr(Cs_*this->delta())*pow(magS,3.0);
    k_ = Ce_*pow(this->delta()*epsilon_, 2.0/3.0);

    k_.correctBoundaryConditions();
    //epsilon_.correctBoundaryConditions();

    correctNut(gradU);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
