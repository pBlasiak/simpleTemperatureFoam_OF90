/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::fixedWallHeatFluxFvPatchScalarField

Description
    This boundary condition sets the pressure gradient to the provided value
    such that the flux on the boundary is that specified by the velocity
    boundary condition.

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fixedWallHeatFlux;
    }
    \endverbatim

See also
    Foam::fixedGradientFvPatchField

SourceFiles
    fixedWallHeatFluxFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedWallHeatFluxFvPatchScalarField_H
#define fixedWallHeatFluxFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class fixedWallHeatFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fixedWallHeatFluxFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private Data

        //- Wall heat flux, W/m2
        scalarField q_;


public:

    //- Runtime type information
    TypeName("fixedWallHeatFlux");


    // Constructors

        //- Construct from patch and internal field
        fixedWallHeatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedWallHeatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedWallHeatFluxFvPatchScalarField onto
        //  a new patch
        fixedWallHeatFluxFvPatchScalarField
        (
            const fixedWallHeatFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        fixedWallHeatFluxFvPatchScalarField
        (
            const fixedWallHeatFluxFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        fixedWallHeatFluxFvPatchScalarField
        (
            const fixedWallHeatFluxFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedWallHeatFluxFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        //- Update the patch pressure gradient field from the given snGradp
        virtual void updateCoeffs(const scalarField& snGradp);

        //- Update the patch pressure gradient field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "volFields.H"
//
//namespace Foam
//{
//    template<class GradBC>
//    inline void setSnGrad
//    (
//        volScalarField::Boundary& bf,
//        const FieldField<fvsPatchField, scalar>& snGrad
//    )
//    {
//        forAll(bf, patchi)
//        {
//            if (isA<GradBC>(bf[patchi]))
//            {
//                refCast<GradBC>(bf[patchi]).updateCoeffs(snGrad[patchi]);
//            }
//        }
//    }
//
//    template<class GradBC>
//    inline void setSnGrad
//    (
//        volScalarField::Boundary& bf,
//        const tmp<FieldField<fvsPatchField, scalar>>& tsnGrad
//    )
//    {
//        setSnGrad<GradBC>(bf, tsnGrad());
//    }
//}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
