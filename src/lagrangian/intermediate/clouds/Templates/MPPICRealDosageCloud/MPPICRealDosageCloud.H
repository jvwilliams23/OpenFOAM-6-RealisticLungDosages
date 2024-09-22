/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::MPPICRealDosageCloud

Description
    Adds MPPIC modelling to kinematic clouds

SourceFiles
    MPPICRealDosageCloudI.H
    MPPICRealDosageCloud.C

\*---------------------------------------------------------------------------*/

#ifndef MPPICRealDosageCloud_H
#define MPPICRealDosageCloud_H

#include "particle.H"
#include "Cloud.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "fvMesh.H"
#include "volFields.H"
#include "AveragingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class PackingModel;

template<class CloudType>
class DampingModel;

template<class CloudType>
class IsotropyModel;

/*---------------------------------------------------------------------------*\
                       Class MPPICRealDosageCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class MPPICRealDosageCloud
:
    public CloudType
{
public:

    // Public typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of parcel the cloud was instantiated for
        typedef typename CloudType::parcelType parcelType;

        //- Convenience typedef for this cloud type
        typedef MPPICRealDosageCloud<CloudType> MPPICCloudType;


private:

    // Private data

        //- Cloud copy pointer
        autoPtr<MPPICRealDosageCloud<CloudType> > cloudCopyPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        MPPICRealDosageCloud(const MPPICRealDosageCloud&);

        //- Disallow default bitwise assignment
        void operator=(const MPPICRealDosageCloud&);
	
	/*
	template<class TrackData>
	void updateParcelCharge(TrackData& td);

	template<class TrackData>
	void setChargeCellValues(TrackData& td);	
	*/
	
protected:

    // Protected data

        // References to the cloud sub-models

            //- Packing model
            autoPtr<PackingModel<MPPICRealDosageCloud<CloudType> > > packingModel_;

            //- Damping model
            autoPtr<DampingModel<MPPICRealDosageCloud<CloudType> > >
                dampingModel_;

            //- Exchange model
            autoPtr<IsotropyModel<MPPICRealDosageCloud<CloudType> > >
                isotropyModel_;


        // Initialisation

            //- Set cloud sub-models
            void setModels();
	
	/*
	//charge density
	volScalarField& rhoe_;
	
	//mean charge
	volScalarField charge_;
	
	//charge deficit
	volScalarField chargeDef_;
	
	//mean work function
	volScalarField workFunction_;
	
	//volume fraction
	volScalarField volumeFraction_;
	
	//particle velocity
	volVectorField us_;
	
	volScalarField granularTemperature_;
			
	const volVectorField& electricField_;
	*/
	

public:

    // Constructors
	
        //- Construct given carrier gas fields
        MPPICRealDosageCloud
        (
	    const word& cloudName,
	    const volScalarField& rho,
	    const volVectorField& U,
	    const volScalarField& mu,
        const dimensionedVector& g,
//		const scalar epsilon0,
//	    const volVectorField& electricField,
//		volScalarField& rhoe,
//		volScalarField& surfCharge,
            bool readFields = true
        );

        //- Copy constructor with new name
        MPPICRealDosageCloud
        (
            MPPICRealDosageCloud<CloudType>& c,
            const word& name
        );

        //- Copy constructor with new name - creates bare cloud
        MPPICRealDosageCloud
        (
            const fvMesh& mesh,
            const word& name,
            const MPPICRealDosageCloud<CloudType>& c
        );

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType> > clone(const word& name)
        {
            return autoPtr<Cloud<parcelType> >
            (
                new MPPICRealDosageCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType> > cloneBare(const word& name) const
        {
            return autoPtr<Cloud<parcelType> >
            (
                new MPPICRealDosageCloud(this->mesh(), name, *this)
            );
        }

    //- Destructor
    virtual ~MPPICRealDosageCloud();

    // Member Functions

        // Access
	    
	    /*
	    inline volScalarField& surfaceCharge();    
	    volScalarField& charge();
	    volScalarField& chargeDef();
	    volScalarField& workFunction();
	    volScalarField& granularTemperature();
	    volScalarField& volumeFraction();
	    volVectorField& us();
	    	    
	    inline const volVectorField& electricField() const;
	    */
	    
            //- Return a reference to the cloud copy
            inline const MPPICRealDosageCloud& cloudCopy() const;

            //- Return const access to the packing model
            inline const PackingModel<MPPICRealDosageCloud<CloudType> >&
                packingModel() const;

            //- Return a reference to the packing model
            inline PackingModel<MPPICRealDosageCloud<CloudType> >& packingModel();

            //- Return condt access to the damping model
            inline const DampingModel<MPPICRealDosageCloud<CloudType> >&
                dampingModel() const;

            //- Return a reference to the damping model
            inline DampingModel<MPPICRealDosageCloud<CloudType> >& dampingModel();

            //- Return condt access to the isotropy model
            inline const IsotropyModel<MPPICRealDosageCloud<CloudType> >&
                isotropyModel() const;

            //- Return a reference to the isotropy model
            inline IsotropyModel<MPPICRealDosageCloud<CloudType> >& isotropyModel();


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

	    //- Return the charge density
	    inline const tmp<volScalarField> rhoe() const; 
	     
		inline void updateRhoe();

        //- I-O

            //- Print cloud information
            void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MPPICRealDosageCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "MPPICRealDosageCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
