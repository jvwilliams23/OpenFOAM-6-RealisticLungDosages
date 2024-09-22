/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::CollidingRealDosageCloud

Description
    Adds coolisions to kinematic clouds

SourceFiles
    CollidingRealDosageCloudI.H
    CollidingRealDosageCloud.C

\*---------------------------------------------------------------------------*/

#ifndef CollidingRealDosageCloud_H
#define CollidingRealDosageCloud_H

#include "particle.H"
#include "Cloud.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class CollisionModel;

/*---------------------------------------------------------------------------*\
                       Class CollidingRealDosageCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CollidingRealDosageCloud
:
    public CloudType
{
public:

    // Public typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of parcel the cloud was instantiated for
        typedef typename CloudType::particleType parcelType;

        //- Convenience typedef for this cloud type
        typedef CollidingRealDosageCloud<CloudType> collidingCloudType;


private:

    // Private data

        //- Cloud copy pointer
        autoPtr<CollidingRealDosageCloud<CloudType>> cloudCopyPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        CollidingRealDosageCloud(const CollidingRealDosageCloud&);

        //- Disallow default bitwise assignment
        void operator=(const CollidingRealDosageCloud&);


protected:

    // Protected data

        //- Thermo parcel constant properties
        typename parcelType::constantProperties constProps_;


        // References to the cloud sub-models

            //- Collision model
            autoPtr<CollisionModel<CollidingRealDosageCloud<CloudType>>>
                collisionModel_;


        // Initialisation

            //- Set cloud sub-models
            void setModels();


        // Cloud evolution functions

            //- Move-collide particles
            template<class TrackCloudType>
            void moveCollide
            (
                TrackCloudType& cloud,
                typename parcelType::trackingData& td,
                const scalar deltaT
            );

            //- Reset state of cloud
            void cloudReset(CollidingRealDosageCloud<CloudType>& c);


public:

    // Constructors

        //- Construct given carrier gas fields
        CollidingRealDosageCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu,
            const dimensionedVector& g,
            bool readFields = true
        );

        //- Copy constructor with new name
        CollidingRealDosageCloud
        (
            CollidingRealDosageCloud<CloudType>& c,
            const word& name
        );

        //- Copy constructor with new name - creates bare cloud
        CollidingRealDosageCloud
        (
            const fvMesh& mesh,
            const word& name,
            const CollidingRealDosageCloud<CloudType>& c
        );

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> clone(const word& name)
        {
            return autoPtr<Cloud<parcelType>>
            (
                new CollidingRealDosageCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> cloneBare(const word& name) const
        {
            return autoPtr<Cloud<parcelType>>
            (
                new CollidingRealDosageCloud(this->mesh(), name, *this)
            );
        }


    //- Destructor
    virtual ~CollidingRealDosageCloud();


    // Member Functions

        // Access

            //- Return a reference to the cloud copy
            inline const CollidingRealDosageCloud& cloudCopy() const;

            //- Return the constant properties
            inline const typename parcelType::constantProperties&
                constProps() const;


            // Sub-models

                //- Return const access to the collision model
                inline const CollisionModel<CollidingRealDosageCloud<CloudType>>&
                    collision() const;

                //- Return reference to the collision model
                inline CollisionModel<CollidingRealDosageCloud<CloudType>>&
                    collision();

        // Check

            //- Total rotational kinetic energy in the system
            inline scalar rotationalKineticEnergyOfSystem() const;


        // Cloud evolution functions

            //- Store the current cloud state
            void storeState();

            //- Reset the current cloud to the previously stored state
            void restoreState();

            //- Evolve the cloud
            void evolve();

            //- Particle motion
            template<class TrackCloudType>
            void motion
            (
                TrackCloudType& cloud,
                typename parcelType::trackingData& td
            );


        // I-O

            //- Print cloud information
            void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CollidingRealDosageCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CollidingRealDosageCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
