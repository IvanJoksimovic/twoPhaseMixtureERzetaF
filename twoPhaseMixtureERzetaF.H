/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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
    Foam::RASModels::ERzetaF

Description
    Lien and Kalitzin's v2-f turbulence model for incompressible and
    compressible flows, with a limit imposed on the turbulent viscosity given
    by Davidson et al.

    The model solves for turbulence kinetic energy k and turbulence dissipation
    rate epsilon, with additional equations for the turbulence stress normal to
    streamlines, v2, and elliptic damping function, f.

    The variant implemented employs N=6, such that f=0 on walls.

    Wall boundary conditions are:

        k       = kLowReWallFunction
        epsilon = epsilonWallFunction
        v2      = v2WallFunction
        f       = fWallFunction

    These are applicable to both low- and high-Reynolds number flows.

    Inlet values can be approximated by:

        v2      = 2/3 k
        f       = zero-gradient

    References:
    \verbatim
        Lien, F. S., & Kalitzin, G. (2001).
        Computations of transonic flow with the ERzetaF turbulence model.
        International Journal of Heat and Fluid Flow, 22(1), 53-61.

        Davidson, L., Nielsen, P., & Sveningsson, A. (2003).
        Modifications of the v2-f model for computing the flow in a
        3D wall jet.
        Turbulence, Heat and Mass Transfer, 4, 577-584
    \endverbatim

    The default model coefficients are
    \verbatim
        ERzetaFCoeffs
        {
            Cmu         0.22;
            CmuKEps     0.09;
            C1          1.4;
            C2          0.3;
            CL          0.23;
            Ceta        70;
            Ceps2       1.9;
            Ceps3       -0.33;
            sigmaEps    1.3;
            sigmaK      1;
        }
    \endverbatim

Note
    If the kLowReWallFunction is employed, a velocity variant of the turbulent
    viscosity wall function should be used, e.g. nutUWallFunction.  Turbulence
    k variants (nutk...) for this case will not behave correctly.

SourceFiles
    ERzetaF.C
\*---------------------------------------------------------------------------*/
#ifndef twoPhaseMixtureERzetaF_H
#define twoPhaseMixtureERzetaF_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class twoPhaseMixtureERzetaF Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class twoPhaseMixtureERzetaF
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
    // Private data

        mutable twoPhaseMixtureERzetaF<BasicTurbulenceModel> *liquidTurbulencePtr_;


    // Private Member Functions

        //- No copy construct
        twoPhaseMixtureERzetaF(const twoPhaseMixtureERzetaF&) = delete;

        //- No copy assignment
        void operator=(const twoPhaseMixtureERzetaF&) = delete;

        //- Return the turbulence model for the other phase
        twoPhaseMixtureERzetaF<BasicTurbulenceModel>& liquidTurbulence() const;


protected:

    // Protected data

        // Model coefficients

        // Model coefficients
            dimensionedScalar Cmu_;
            dimensionedScalar COmega2_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar sigmaK_;
            dimensionedScalar sigmaOmega_;
            dimensionedScalar sigmaCDv_;
            dimensionedScalar sigmaCDt_;
            dimensionedScalar sigmaZeta_;
            dimensionedScalar CTau_;
            dimensionedScalar CL_;
            dimensionedScalar CEta_;
            dimensionedScalar Csas_;
            dimensionedScalar CT2_;
            dimensionedScalar Clim_;

        // Fields

            volScalarField k_;
            volScalarField omega_;
            volScalarField epsilon_;
            volScalarField zeta_;
            volScalarField delta_;
            volScalarField f_;
            const volScalarField& yr_;

            //- Lower limit of zeta
            dimensionedScalar mTSmall_;
            dimensionedScalar zetaMin_;
            //- Lower limit of f
            dimensionedScalar fMin_;

        // Mixture fields

            autoPtr<volScalarField> Ct2_;
            autoPtr<volScalarField> rhom_;
            autoPtr<volScalarField> km_;
            //autoPtr<volScalarField> epsilonm_;
            autoPtr<volScalarField> omegam_;
            autoPtr<volScalarField> zetam_;
            autoPtr<volScalarField> deltam_;
            autoPtr<volScalarField> fm_;


    // Protected Member Functions

        wordList epsilonBoundaryTypes(const volScalarField& epsilon) const;

        void correctInletOutlet
        (
            volScalarField& vsf,
            const volScalarField& refVsf
        ) const;

        void initMixtureFields();

        virtual void correctNut();

        tmp<volScalarField> Ct2() const;

        tmp<volScalarField> rholEff() const;
        tmp<volScalarField> rhogEff() const;
        tmp<volScalarField> rhom() const;

        tmp<volScalarField> mix
        (
            const volScalarField& fc,
            const volScalarField& fd
        ) const;

        tmp<volScalarField> mixU
        (
            const volScalarField& fc,
            const volScalarField& fd
        ) const;

        tmp<surfaceScalarField> mixFlux
        (
            const surfaceScalarField& fc,
            const surfaceScalarField& fd
        ) const;

        tmp<volScalarField> bubbleG() const;
        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> epsilonSource() const;
/*
        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff(const volScalarField& nutm) const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkEff",
                    nutm/sigmak_
                )
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff(const volScalarField& nutm) const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DepsilonEff",
                    nutm/sigmaEps_
                )
            );
        }
*/
public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("twoPhaseMixtureERzetaF");


    // Constructors

        //- Construct from components
        twoPhaseMixtureERzetaF
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
    virtual ~twoPhaseMixtureERzetaF() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkEff",
                    this->nut_/sigmaK_ + this->nu()
                )
            );
        }

            //- Return the effective diffusivity for omega
            tmp<volScalarField> DomegaEff() const
            {
                return tmp<volScalarField>
                (
                    new volScalarField
                    (
                        "DomegaEff",
                        this->nut_/sigmaOmega_ + this->nu()
                    )
                );
            }

            //
            tmp<volScalarField> DzetaEff() const
            {
                return tmp<volScalarField>
                (
                    new volScalarField
                    (
                        "DzetaEff",
                        this->nut_/sigmaZeta_ + this->nu()
                    )
                );
            }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Return turbulence stress normal to streamlines
        virtual tmp<volScalarField> zeta() const
        {
            return zeta_;
        }

        //- Return the damping function
        virtual tmp<volScalarField> f() const
        {
            return f_;
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "twoPhaseMixtureERzetaF.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
