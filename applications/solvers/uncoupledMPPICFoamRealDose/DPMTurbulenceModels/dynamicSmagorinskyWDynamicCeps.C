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

#include "dynamicSmagorinskyWDynamicCeps.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
volScalarField dynamicSmagorinskyWDynamicCeps<BasicTurbulenceModel>::Ce
(
    const volSymmTensorField& D,
    const volScalarField& KK
) const
{
    const volScalarField Ce
    (
        simpleFilter_(this->nuEff()*(filter_(magSqr(D)) - magSqr(filter_(D))))
       /simpleFilter_(pow(KK, 1.5)/(2.0*this->delta()))
    );

    tmp<volScalarField> tfld = 0.5*(mag(Ce) + Ce);
    return tfld();
}


template<class BasicTurbulenceModel>
volScalarField dynamicSmagorinskyWDynamicCeps<BasicTurbulenceModel>::CeDynamic() const
{
    const volSymmTensorField D(dev(symm(fvc::grad(this->U_))));

    volScalarField KK
    (
        0.5*(filter_(magSqr(this->U_)) - magSqr(filter_(this->U_)))
    );
    KK.max(dimensionedScalar("small", KK.dimensions(), small));

    return Ce(D, KK);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void dynamicSmagorinskyWDynamicCeps<BasicTurbulenceModel>::correctNut
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
void dynamicSmagorinskyWDynamicCeps<BasicTurbulenceModel>::correctNut()
{
    correctNut(fvc::grad(this->U_));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
dynamicSmagorinskyWDynamicCeps<BasicTurbulenceModel>::dynamicSmagorinskyWDynamicCeps
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
        IOobject
        (
            IOobject::groupName("Cs", this->U_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar ("Cs", dimless,SMALL)
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

    simpleFilter_(U.mesh()),
    filterPtr_(LESfilter::New(U.mesh(), this->coeffDict())),
    filter_(filterPtr_())
{
//    bound(k_, this->kMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool dynamicSmagorinskyWDynamicCeps<BasicTurbulenceModel>::read()
{
    if (LESeddyViscosity<BasicTurbulenceModel>::read())
    {
        filter_.read(this->coeffDict());        

        return true;
    }
    else
    {
        return false;
    }
}

template<class BasicTurbulenceModel>
tmp<volScalarField> dynamicSmagorinskyWDynamicCeps<BasicTurbulenceModel>::epsilon() const
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
            //CeDynamic()*mag(k)*sqrt(mag(k))/this->delta()

            this->Ce_*mag(k)*sqrt(mag(k))/this->delta()
         )
     );
}

template<class BasicTurbulenceModel>
void dynamicSmagorinskyWDynamicCeps<BasicTurbulenceModel>::correct()
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

    volSymmTensorField S(dev(symm(gradU)));
    volScalarField magS(mag(S));

    volVectorField Uf(filter_(U));

    volSymmTensorField Sf(filter_(S));  
//    volSymmTensorField Sf(dev(symm(fvc::grad(Uf))));
    
    volScalarField magSf(mag(Sf));
          
    
    volSymmTensorField LL =
    (dev(filter_(sqr(U)) - (sqr(filter_(U)))));

    volSymmTensorField MM
    (
        sqr(this->delta())*(filter_(magS*S) - 4.0*magSf*Sf)
    );
    
    volScalarField MMMM = fvc::average(magSqr(MM));
    MMMM.max(VSMALL);

    Cs_= 0.5* fvc::average(LL && MM)/MMMM;

    volScalarField KK =
    0.5*(filter_(magSqr(U)) - magSqr(filter_(U)));
    
    volScalarField mm
    (
        sqr(this->delta())*(4.0*sqr(mag(Sf)) - filter_(sqr(magS)))
       
    );

    volScalarField mmmm = fvc::average(magSqr(mm));
    mmmm.max(VSMALL);

    //k_ = fvc::average(KK*mm)/mmmm * sqr(this->delta())*magSqr(S);
    //epsilon_ = epsilon();

    // JW 26/04/21
    epsilon_ = mag(Cs_) * sqr(this->delta())*pow(magS,3.0);
    //k_ = CeDynamic()*pow(this->delta()*mag(epsilon_), 2.0/3.0);
    k_ = Ce_*pow(this->delta()*mag(epsilon_), 2.0/3.0);

    k_.correctBoundaryConditions();
    //epsilon_.correctBoundaryConditions();
    correctNut(gradU);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
