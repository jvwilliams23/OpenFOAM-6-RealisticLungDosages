/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::LESModels::dynamicSmagorinskyWDynamicCeps

Group
    grpLESTurbulence

Description
    Dynamic SGS model with local averaging, based on Germano and Lilly works.
    
    An older version of this code had been developed for version 2.3.1 and older by Alberto Passalacqua.
    
    The current version has been developed for version 4.1, Ehsan Asgari, e.asgari@aut.ac.ir

    The makeTurbulentFluidThermoModels.C file has been updated for OpenFOAM version 6, Eric Daymo, eadaymo@tonkomo.com 

SourceFiles
    dynamicSmagorinskyWDynamicCeps.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicSmagorinskyWDynamicCeps_H
#define dynamicSmagorinskyWDynamicCeps_H

#include "LESeddyViscosity.H"
#include "simpleFilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                       Class dynamicSmagorinskyWDynamicCeps Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class dynamicSmagorinskyWDynamicCeps
:
    public LESeddyViscosity<BasicTurbulenceModel>
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        dynamicSmagorinskyWDynamicCeps(const dynamicSmagorinskyWDynamicCeps&);
        void operator=(const dynamicSmagorinskyWDynamicCeps&);


protected:

    // Protected data

        volScalarField k_;
        volScalarField epsilon_;
        volScalarField Cs_;
        const dimensionedScalar Ce_;

//       volScalarField beta_;
   
        simpleFilter simpleFilter_;
        autoPtr<LESfilter> filterPtr_;
        LESfilter& filter_;

      
    


    // Protected Member Functions
        //- Calculate Ce by filtering the velocity field U
        volScalarField Ce
        (
            const volSymmTensorField& D,
            const volScalarField& KK
        ) const;

        volScalarField CeDynamic() const;
        

        //- Update sub-grid eddy-viscosity
        void correctNut(const tmp<volTensorField>& gradU);

        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;

    //- Runtime type information
    TypeName("dynamicSmagorinskyWDynamicCeps");


    // Constructors

        //- Construct from components
        dynamicSmagorinskyWDynamicCeps
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~dynamicSmagorinskyWDynamicCeps()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", this->nut_ + this->nu())
            );
        }

         //- Return sub-grid disipation rate
         virtual tmp<volScalarField> epsilon() const;

        //- Correct Eddy-Viscosity and related properties
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "dynamicSmagorinskyWDynamicCeps.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
