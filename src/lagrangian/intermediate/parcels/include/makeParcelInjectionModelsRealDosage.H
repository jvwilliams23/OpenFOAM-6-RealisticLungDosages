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

#ifndef makeParcelInjectionModels_H
#define makeParcelInjectionModels_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CellZoneInjection.H"
#include "ConeInjection.H"
#include "ConeNozzleInjection.H"
#include "FieldActivatedInjection.H"
#include "InflationInjection.H"
#include "KinematicLookupTableInjection.H"
#include "ManualInjection.H"
#include "NoInjection.H"
#include "PatchInjection.H"
#include "PatchFlowRateInjection.H"
//#include "ManualInjectionRealDosage.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeParcelInjectionModels(CloudType)                                   \
                                                                               \
    makeInjectionModel(CloudType);                                             \
                                                                               \
    makeInjectionModelType(CellZoneInjection, CloudType);                      \
    makeInjectionModelType(ConeInjection, CloudType);                          \
    makeInjectionModelType(ConeNozzleInjection, CloudType);                    \
    makeInjectionModelType(FieldActivatedInjection, CloudType);                \
    makeInjectionModelType(InflationInjection, CloudType);                     \
    makeInjectionModelType(KinematicLookupTableInjection, CloudType);          \
    makeInjectionModelType(ManualInjection, CloudType);                        \
    makeInjectionModelType(NoInjection, CloudType);                            \
    makeInjectionModelType(PatchInjection, CloudType);                         \
    makeInjectionModelType(PatchFlowRateInjection, CloudType);                 \
    //makeInjectionModelType(ManualInjectionRealDosage, CloudType);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
