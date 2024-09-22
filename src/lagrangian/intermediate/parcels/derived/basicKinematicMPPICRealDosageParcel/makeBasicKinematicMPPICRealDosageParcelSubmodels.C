/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2015 OpenFOAM Foundation
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

#include "basicKinematicMPPICRealDosageCloud.H"

#include "makeParcelCloudFunctionObjects.H"

// Kinematic sub-models
//#include "makeParcelForcesRealDosage.H"
#include "makeParcelForces.H"
#include "makeParcelDispersionModels.H"
//#include "makeParcelInjectionModelsRealDosage.H"
#include "makeParcelInjectionModels.H"
//#include "makeParcelPatchInteractionModelsRealDosage.H"
#include "makeParcelPatchInteractionModels.H"
#include "makeParcelStochasticCollisionModels.H"
#include "makeParcelSurfaceFilmModels.H"

// MPPIC sub-models
#include "makeMPPICParcelDampingModels.H"
#include "makeMPPICParcelIsotropyModels.H"
#include "makeMPPICParcelPackingModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeParcelCloudFunctionObjects(basicKinematicMPPICRealDosageCloud);

// Kinematic sub-models
makeParcelForces(basicKinematicMPPICRealDosageCloud);

makeParcelDispersionModels(basicKinematicMPPICRealDosageCloud);
makeParcelInjectionModels(basicKinematicMPPICRealDosageCloud);
makeParcelPatchInteractionModels(basicKinematicMPPICRealDosageCloud);
makeParcelStochasticCollisionModels(basicKinematicMPPICRealDosageCloud);
makeParcelSurfaceFilmModels(basicKinematicMPPICRealDosageCloud);

// MPPIC sub-models
makeMPPICParcelDampingModels(basicKinematicMPPICRealDosageCloud);
makeMPPICParcelIsotropyModels(basicKinematicMPPICRealDosageCloud);
makeMPPICParcelPackingModels(basicKinematicMPPICRealDosageCloud);

// ************************************************************************* //
