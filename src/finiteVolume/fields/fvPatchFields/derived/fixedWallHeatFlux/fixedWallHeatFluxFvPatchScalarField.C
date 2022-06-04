/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "fixedWallHeatFluxFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedWallHeatFluxFvPatchScalarField::fixedWallHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    q_(p.size(), 0)
{}


Foam::fixedWallHeatFluxFvPatchScalarField::fixedWallHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict, false),
    q_("q", dict, p.size())
{
    if (dict.found("value") && dict.found("gradient"))
    {
        fvPatchField<scalar>::operator=(scalarField("value", dict, p.size()));
        gradient() = scalarField("gradient", dict, p.size());
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = Zero;
    }
}


Foam::fixedWallHeatFluxFvPatchScalarField::fixedWallHeatFluxFvPatchScalarField
(
    const fixedWallHeatFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    q_(p.size(), 0)
{}


Foam::fixedWallHeatFluxFvPatchScalarField::fixedWallHeatFluxFvPatchScalarField
(
    const fixedWallHeatFluxFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF),
    q_(wbppsf.q_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//void Foam::fixedWallHeatFluxFvPatchScalarField::updateCoeffs
//(
//    const scalarField& snGradp
//)
//{
//    if (updated())
//    {
//        return;
//    }
//
//    curTimeIndex_ = this->db().time().timeIndex();
//
//    gradient() = snGradp;
//    fixedGradientFvPatchScalarField::updateCoeffs();
//}


void Foam::fixedWallHeatFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

	const scalarField& alphaEff = 
		patch().lookupPatchField<volScalarField, scalar>("alphaEff");
	gradient() = q_/alphaEff;

    //if (curTimeIndex_ != this->db().time().timeIndex())
    //{
    //    FatalErrorInFunction
    //        << "updateCoeffs(const scalarField& snGradp) MUST be called before"
    //           " updateCoeffs() or evaluate() to set the boundary gradient."
    //        << exit(FatalError);
    //}
}


void Foam::fixedWallHeatFluxFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedWallHeatFluxFvPatchScalarField
    );
}


// ************************************************************************* //
