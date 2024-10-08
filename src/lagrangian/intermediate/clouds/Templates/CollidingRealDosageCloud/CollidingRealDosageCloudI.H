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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::CollidingRealDosageCloud<CloudType>&
Foam::CollidingRealDosageCloud<CloudType>::cloudCopy() const
{
    return cloudCopyPtr_();
}


template<class CloudType>
inline const typename CloudType::particleType::constantProperties&
Foam::CollidingRealDosageCloud<CloudType>::constProps() const
{
    return constProps_;
}


template<class CloudType>
inline const Foam::CollisionModel<Foam::CollidingRealDosageCloud<CloudType>>&
Foam::CollidingRealDosageCloud<CloudType>::collision() const
{
    return collisionModel_();
}


template<class CloudType>
inline Foam::CollisionModel<Foam::CollidingRealDosageCloud<CloudType>>&
Foam::CollidingRealDosageCloud<CloudType>::collision()
{
    return collisionModel_();
}


template<class CloudType>
inline Foam::scalar
Foam::CollidingRealDosageCloud<CloudType>::rotationalKineticEnergyOfSystem() const
{
    scalar rotationalKineticEnergy = 0.0;

    forAllConstIter(typename CollidingRealDosageCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();

        rotationalKineticEnergy +=
            p.nParticle()*0.5*p.momentOfInertia()*(p.omega() & p.omega());
    }

    return rotationalKineticEnergy;
}


// ************************************************************************* //
