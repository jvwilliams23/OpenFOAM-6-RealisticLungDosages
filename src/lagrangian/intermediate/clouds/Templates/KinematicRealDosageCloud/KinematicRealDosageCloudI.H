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

#include "fvmSup.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::KinematicRealDosageCloud<CloudType>&
Foam::KinematicRealDosageCloud<CloudType>::cloudCopy() const
{
    return cloudCopyPtr_();
}


template<class CloudType>
inline const Foam::fvMesh& Foam::KinematicRealDosageCloud<CloudType>::mesh() const
{
    return mesh_;
}


template<class CloudType>
inline const Foam::IOdictionary&
Foam::KinematicRealDosageCloud<CloudType>::particleProperties() const
{
    return particleProperties_;
}


template<class CloudType>
inline const Foam::IOdictionary&
Foam::KinematicRealDosageCloud<CloudType>::outputProperties() const
{
    return outputProperties_;
}


template<class CloudType>
inline Foam::IOdictionary& Foam::KinematicRealDosageCloud<CloudType>::outputProperties()
{
    return outputProperties_;
}


template<class CloudType>
inline const Foam::cloudSolution&
Foam::KinematicRealDosageCloud<CloudType>::solution() const
{
    return solution_;
}


template<class CloudType>
inline Foam::cloudSolution& Foam::KinematicRealDosageCloud<CloudType>::solution()
{
    return solution_;
}


template<class CloudType>
inline const typename CloudType::particleType::constantProperties&
Foam::KinematicRealDosageCloud<CloudType>::constProps() const
{
    return constProps_;
}


template<class CloudType>
inline typename CloudType::particleType::constantProperties&
Foam::KinematicRealDosageCloud<CloudType>::constProps()
{
    return constProps_;
}


template<class CloudType>
inline const Foam::dictionary&
Foam::KinematicRealDosageCloud<CloudType>::subModelProperties() const
{
    return subModelProperties_;
}


template<class CloudType>
inline const Foam::volScalarField& Foam::KinematicRealDosageCloud<CloudType>::rho() const
{
    return rho_;
}


template<class CloudType>
inline const Foam::volVectorField& Foam::KinematicRealDosageCloud<CloudType>::U() const
{
    return U_;
}


template<class CloudType>
inline const Foam::volScalarField& Foam::KinematicRealDosageCloud<CloudType>::mu() const
{
    return mu_;
}


template<class CloudType>
inline const Foam::dimensionedVector& Foam::KinematicRealDosageCloud<CloudType>::g() const
{
    return g_;
}


template<class CloudType>
inline Foam::scalar Foam::KinematicRealDosageCloud<CloudType>::pAmbient() const
{
    return pAmbient_;
}


template<class CloudType>
inline Foam::scalar& Foam::KinematicRealDosageCloud<CloudType>::pAmbient()
{
    return pAmbient_;
}


template<class CloudType>
//inline const typename CloudType::parcelType::forceType&
inline const typename Foam::KinematicRealDosageCloud<CloudType>::forceType&
Foam::KinematicRealDosageCloud<CloudType>::forces() const
{
    return forces_;
}


template<class CloudType>
inline typename Foam::KinematicRealDosageCloud<CloudType>::forceType&
Foam::KinematicRealDosageCloud<CloudType>::forces()
{
    return forces_;
}


template<class CloudType>
inline typename Foam::KinematicRealDosageCloud<CloudType>::functionType&
Foam::KinematicRealDosageCloud<CloudType>::functions()
{
    return functions_;
}


template<class CloudType>
inline const Foam::InjectionModelList<Foam::KinematicRealDosageCloud<CloudType>>&
Foam::KinematicRealDosageCloud<CloudType>::injectors() const
{
    return injectors_;
}


template<class CloudType>
inline Foam::InjectionModelList<Foam::KinematicRealDosageCloud<CloudType>>&
Foam::KinematicRealDosageCloud<CloudType>::injectors()
{
    return injectors_;
}


template<class CloudType>
inline const Foam::DispersionModel<Foam::KinematicRealDosageCloud<CloudType>>&
Foam::KinematicRealDosageCloud<CloudType>::dispersion() const
{
    return dispersionModel_;
}


template<class CloudType>
inline Foam::DispersionModel<Foam::KinematicRealDosageCloud<CloudType>>&
Foam::KinematicRealDosageCloud<CloudType>::dispersion()
{
    return dispersionModel_();
}


template<class CloudType>
inline const Foam::PatchInteractionModel<Foam::KinematicRealDosageCloud<CloudType>>&
Foam::KinematicRealDosageCloud<CloudType>::patchInteraction() const
{
    return patchInteractionModel_;
}


template<class CloudType>
inline Foam::PatchInteractionModel<Foam::KinematicRealDosageCloud<CloudType>>&
Foam::KinematicRealDosageCloud<CloudType>::patchInteraction()
{
    return patchInteractionModel_();
}


template<class CloudType>
inline const Foam::StochasticCollisionModel<Foam::KinematicRealDosageCloud<CloudType>>&
Foam::KinematicRealDosageCloud<CloudType>::stochasticCollision() const
{
    return stochasticCollisionModel_();
}


template<class CloudType>
inline Foam::StochasticCollisionModel<Foam::KinematicRealDosageCloud<CloudType>>&
Foam::KinematicRealDosageCloud<CloudType>::stochasticCollision()
{
    return stochasticCollisionModel_();
}


template<class CloudType>
inline const Foam::SurfaceFilmModel<Foam::KinematicRealDosageCloud<CloudType>>&
Foam::KinematicRealDosageCloud<CloudType>::surfaceFilm() const
{
    return surfaceFilmModel_();
}


template<class CloudType>
inline Foam::SurfaceFilmModel<Foam::KinematicRealDosageCloud<CloudType>>&
Foam::KinematicRealDosageCloud<CloudType>::surfaceFilm()
{
    return surfaceFilmModel_();
}


template<class CloudType>
inline const Foam::integrationScheme&
Foam::KinematicRealDosageCloud<CloudType>::UIntegrator() const
{
    return UIntegrator_;
}


template<class CloudType>
inline Foam::label Foam::KinematicRealDosageCloud<CloudType>::nParcels() const
{
    return this->size();
}


template<class CloudType>
inline Foam::scalar Foam::KinematicRealDosageCloud<CloudType>::massInSystem() const
{
    scalar sysMass = 0.0;
    forAllConstIter(typename KinematicRealDosageCloud<CloudType>, *this, iter)
    {
         const parcelType& p = iter();
         sysMass += p.nParticle()*p.mass();
    }

    return sysMass;
}


template<class CloudType>
inline Foam::vector
Foam::KinematicRealDosageCloud<CloudType>::linearMomentumOfSystem() const
{
    vector linearMomentum(Zero);

    forAllConstIter(typename KinematicRealDosageCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();

        linearMomentum += p.nParticle()*p.mass()*p.U();
    }

    return linearMomentum;
}


template<class CloudType>
inline Foam::scalar
Foam::KinematicRealDosageCloud<CloudType>::linearKineticEnergyOfSystem() const
{
    scalar linearKineticEnergy = 0.0;

    forAllConstIter(typename KinematicRealDosageCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();

        linearKineticEnergy += p.nParticle()*0.5*p.mass()*(p.U() & p.U());
    }

    return linearKineticEnergy;
}


template<class CloudType>
inline Foam::scalar Foam::KinematicRealDosageCloud<CloudType>::Dij
(
    const label i,
    const label j
) const
{
    scalar si = 0.0;
    scalar sj = 0.0;
    forAllConstIter(typename KinematicRealDosageCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        si += p.nParticle()*pow(p.d(), i);
        sj += p.nParticle()*pow(p.d(), j);
    }

    reduce(si, sumOp<scalar>());
    reduce(sj, sumOp<scalar>());
    sj = max(sj, vSmall);

    return si/sj;
}


template<class CloudType>
inline Foam::scalar Foam::KinematicRealDosageCloud<CloudType>::Dmax() const
{
    scalar d = -great;
    forAllConstIter(typename KinematicRealDosageCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        d = max(d, p.d());
    }

    reduce(d, maxOp<scalar>());

    return max(0.0, d);
}


template<class CloudType>
inline Foam::Random& Foam::KinematicRealDosageCloud<CloudType>::rndGen() const
{
    return rndGen_;
}


template<class CloudType>
inline Foam::List<Foam::DynamicList<typename CloudType::particleType*>>&
Foam::KinematicRealDosageCloud<CloudType>::cellOccupancy()
{
    if (cellOccupancyPtr_.empty())
    {
        buildCellOccupancy();
    }

    return cellOccupancyPtr_();
}


template<class CloudType>
inline const Foam::scalarField&
Foam::KinematicRealDosageCloud<CloudType>::cellLengthScale() const
{
    return cellLengthScale_;
}


template<class CloudType>
inline Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::KinematicRealDosageCloud<CloudType>::UTrans()
{
    return UTrans_();
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::KinematicRealDosageCloud<CloudType>::UTrans() const
{
    return UTrans_();
}


template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::KinematicRealDosageCloud<CloudType>::UCoeff()
{
    return UCoeff_();
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::KinematicRealDosageCloud<CloudType>::UCoeff() const
{
    return UCoeff_();
}


template<class CloudType>
inline Foam::tmp<Foam::fvVectorMatrix>
Foam::KinematicRealDosageCloud<CloudType>::SU(volVectorField& U) const
{
    if (debug)
    {
        Info<< "UTrans min/max = " << min(UTrans()).value() << ", "
            << max(UTrans()).value() << nl
            << "UCoeff min/max = " << min(UCoeff()).value() << ", "
            << max(UCoeff()).value() << endl;
    }

    if (solution_.coupled())
    {
        if (solution_.semiImplicit("U"))
        {
            const volScalarField::Internal
                Vdt(mesh_.V()*this->db().time().deltaT());

            return UTrans()/Vdt - fvm::Sp(UCoeff()/Vdt, U) + UCoeff()/Vdt*U;
        }
        else
        {
            tmp<fvVectorMatrix> tfvm(new fvVectorMatrix(U, dimForce));
            fvVectorMatrix& fvm = tfvm.ref();

            fvm.source() = -UTrans()/(this->db().time().deltaT());

            return tfvm;
        }
    }

    return tmp<fvVectorMatrix>(new fvVectorMatrix(U, dimForce));
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::KinematicRealDosageCloud<CloudType>::vDotSweep() const
{
    tmp<volScalarField> tvDotSweep
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":vDotSweep",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimTime, 0.0),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );

    volScalarField& vDotSweep = tvDotSweep.ref();
    forAllConstIter(typename KinematicRealDosageCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label celli = p.cell();

        if (p.active() == 1)
        {
            vDotSweep[celli] += p.nParticle()*p.areaP()*mag(p.U() - U_[celli]);
        }
    }

    vDotSweep.primitiveFieldRef() /= mesh_.V();
    vDotSweep.correctBoundaryConditions();

    return tvDotSweep;
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::KinematicRealDosageCloud<CloudType>::theta() const
{
    tmp<volScalarField> ttheta
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":theta",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );

    volScalarField& theta = ttheta.ref();
    forAllConstIter(typename KinematicRealDosageCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label celli = p.cell();

        if (p.active() == 1)
        {
            theta[celli] += p.nParticle()*p.volume();
        }
    }

    theta.primitiveFieldRef() /= mesh_.V();
    theta.correctBoundaryConditions();

    return ttheta;
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::KinematicRealDosageCloud<CloudType>::alpha() const
{
    tmp<volScalarField> talpha
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":alpha",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    scalarField& alpha = talpha.ref().primitiveFieldRef();
    forAllConstIter(typename KinematicRealDosageCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label celli = p.cell();

        if (p.active() == 1)
        {
            alpha[celli] += p.nParticle()*p.mass();
        }
    }

    alpha /= (mesh_.V()*rho_);

    return talpha;
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::KinematicRealDosageCloud<CloudType>::rhoEff() const
{
    tmp<volScalarField> trhoEff
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":rhoEff",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimDensity, 0.0)
        )
    );

    scalarField& rhoEff = trhoEff.ref().primitiveFieldRef();
    forAllConstIter(typename KinematicRealDosageCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label celli = p.cell();

        if (p.active() == 1)
        {
            rhoEff[celli] += p.nParticle()*p.mass();
        }
    }

    rhoEff /= mesh_.V();

    return trhoEff;
}


// ************************************************************************* //
