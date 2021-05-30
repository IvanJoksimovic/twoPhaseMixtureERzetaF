/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "twoPhaseMixtureERzetaF.H"
#include "fvOptions.H"
#include "bound.H"
#include "twoPhaseSystem.H"
#include "dragModel.H"
#include "virtualMassModel.H"
#include "fixedValueFvPatchFields.H"
#include "inletOutletFvPatchFields.H"
#include "fvmSup.H"
#include "wallDist.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
twoPhaseMixtureERzetaF<BasicTurbulenceModel>::twoPhaseMixtureERzetaF
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    liquidTurbulencePtr_(nullptr),

    Cmu_		(dimensioned<scalar>::lookupOrAddToDict("Cmu", 			this->coeffDict_, 0.22)),
    COmega2_		(dimensioned<scalar>::lookupOrAddToDict("COmega2", 		this->coeffDict_, 0.9)),
    C1_			(dimensioned<scalar>::lookupOrAddToDict("C1", 			this->coeffDict_, 0.4)),
    C2_			(dimensioned<scalar>::lookupOrAddToDict("C2", 			this->coeffDict_, 0.65)),
    sigmaK_		(dimensioned<scalar>::lookupOrAddToDict("sigmaK", 		this->coeffDict_, 1.1)),
    sigmaOmega_		(dimensioned<scalar>::lookupOrAddToDict("sigmaOmega", 		this->coeffDict_, 1.1)),
    sigmaCDv_		(dimensioned<scalar>::lookupOrAddToDict("sigmaCDv", 		this->coeffDict_, 1.2)),
    sigmaCDt_		(dimensioned<scalar>::lookupOrAddToDict("sigmaCDt", 		this->coeffDict_, 1.6)),
    sigmaZeta_		(dimensioned<scalar>::lookupOrAddToDict("sigmaZeta", 		this->coeffDict_, 1.2)),
    CTau_		(dimensioned<scalar>::lookupOrAddToDict("CTau", 		this->coeffDict_, 6.0)),
    CL_			(dimensioned<scalar>::lookupOrAddToDict("CL", 			this->coeffDict_, 0.36)),
    CEta_		(dimensioned<scalar>::lookupOrAddToDict("CEta", 		this->coeffDict_, 85)),
    Csas_		(dimensioned<scalar>::lookupOrAddToDict("Csas", 		this->coeffDict_, 4)),
    CT2_		(dimensioned<scalar>::lookupOrAddToDict("CT2", 			this->coeffDict_, 1.0)),
    Clim_		(dimensioned<scalar>::lookupOrAddToDict("Clim", 		this->coeffDict_, 0.0)),

    k_ 		(IOobject (IOobject::groupName("k", 		alphaRhoPhi.group()), this->runTime_.timeName(),this->mesh_, IOobject::MUST_READ, IOobject::AUTO_WRITE ), this->mesh_),
    omega_ 	(IOobject (IOobject::groupName("omega", 	alphaRhoPhi.group()), this->runTime_.timeName(),this->mesh_, IOobject::MUST_READ, IOobject::AUTO_WRITE ), this->mesh_),
    epsilon_ 	(IOobject (IOobject::groupName("epsilon", 	alphaRhoPhi.group()), this->runTime_.timeName(),this->mesh_, IOobject::NO_READ, IOobject::NO_WRITE ), omega_*k_),
    zeta_ 	(IOobject (IOobject::groupName("zeta", 		alphaRhoPhi.group()), this->runTime_.timeName(),this->mesh_, IOobject::MUST_READ, IOobject::AUTO_WRITE ), this->mesh_),

    delta_
    (
        IOobject
        (
            IOobject::groupName("delta", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("delta", dimLength, 1e-10)
    ),

    f_ 		(IOobject (IOobject::groupName("f", 		alphaRhoPhi.group()), this->runTime_.timeName(),this->mesh_, IOobject::MUST_READ, IOobject::AUTO_WRITE ), this->mesh_),
    yr_(wallDist::New(this->mesh_).y()),
    mTSmall_("mTSmall", dimensionSet(0, 0, -1, 0, 0, 0, 0),1e-10),
    zetaMin_("zetaMin", dimless, 1e-10),
    fMin_("fMin", dimless/dimTime, 1e-10)

  {
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);
    bound(zeta_, zetaMin_);
    bound(f_, fMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
  }


/*
template<class BasicTurbulenceModel>
wordList twoPhaseMixtureERzetaF<BasicTurbulenceModel>::epsilonBoundaryTypes
(
    const volScalarField& epsilon
) const
{
    const volScalarField::Boundary& ebf = epsilon.boundaryField();

    wordList ebt = ebf.types();

    forAll(ebf, patchi)
    {
        if (isA<fixedValueFvPatchScalarField>(ebf[patchi]))
        {
            ebt[patchi] = fixedValueFvPatchScalarField::typeName;
        }
    }

    return ebt;
}
*/

template<class BasicTurbulenceModel>
void twoPhaseMixtureERzetaF<BasicTurbulenceModel>::correctInletOutlet
(
    volScalarField& vsf,
    const volScalarField& refVsf
) const
{
    volScalarField::Boundary& bf = vsf.boundaryFieldRef();
    const volScalarField::Boundary& refBf =
        refVsf.boundaryField();

    forAll(bf, patchi)
    {
        if
        (
            isA<inletOutletFvPatchScalarField>(bf[patchi])
         && isA<inletOutletFvPatchScalarField>(refBf[patchi])
        )
        {
            refCast<inletOutletFvPatchScalarField>
            (bf[patchi]).refValue() =
            refCast<const inletOutletFvPatchScalarField>
            (refBf[patchi]).refValue();
        }
    }
}


template<class BasicTurbulenceModel>
void twoPhaseMixtureERzetaF<BasicTurbulenceModel>::initMixtureFields()
{
    if (rhom_) return;

    // Local references to gas-phase properties
    const volScalarField& kg = this->k_;
    const volScalarField& zetag = this->zeta_;
    const volScalarField& omegag = this->omega_;
    const volScalarField& fg = this->f_;

    // Local references to liquid-phase properties
    twoPhaseMixtureERzetaF<BasicTurbulenceModel>& turbc = this->liquidTurbulence();
    const volScalarField& kl = turbc.k_;
    const volScalarField& zetal = turbc.zeta_;
    const volScalarField& omegal = turbc.omega_;
    const volScalarField& fl = turbc.f_;

    word startTimeName
    (
        this->runTime_.timeName(this->runTime_.startTime().value())
    );

    Ct2_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "Ct2",
                startTimeName,
                this->mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            Ct2()
        )
    );

    rhom_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "rhom",
                startTimeName,
                this->mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            rhom()
        )
    );

    km_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "km",
                startTimeName,
                this->mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mix(kl, kg),
            kl.boundaryField().types()
        )
    );
    correctInletOutlet(km_(), kl);

    zetam_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "zetam",
                startTimeName,
                this->mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mix(zetal, zetag),
            zetal.boundaryField().types()
        )
    );
    correctInletOutlet(zetam_(),zetal);

    omegam_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "omegam",
                startTimeName,
                this->mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mix(omegal, omegag),
            omegal.boundaryField().types()
        )
    );
    correctInletOutlet(omegam_(), omegal);

    fm_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "fm",
                startTimeName,
                this->mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mix(fl,fg),
            fl.boundaryField().types()
        )
    );
    correctInletOutlet(fm_(), fl);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool twoPhaseMixtureERzetaF<BasicTurbulenceModel>::read()
{
  if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
  {
      Cmu_.readIfPresent(this->coeffDict());
      COmega2_.readIfPresent(this->coeffDict());
      C1_.readIfPresent(this->coeffDict());
      C2_.readIfPresent(this->coeffDict());
      sigmaK_.readIfPresent(this->coeffDict());
      sigmaOmega_.readIfPresent(this->coeffDict());
      sigmaCDv_.readIfPresent(this->coeffDict());
      sigmaCDt_.readIfPresent(this->coeffDict());
      sigmaZeta_.readIfPresent(this->coeffDict());
      CTau_.readIfPresent(this->coeffDict());
      CL_.readIfPresent(this->coeffDict());
      CEta_.readIfPresent(this->coeffDict());
      Csas_.readIfPresent(this->coeffDict());
      CT2_.readIfPresent(this->coeffDict());
      Clim_.readIfPresent(this->coeffDict());


      return true;
  }
  else
  {
      return false;
  }
}

/*
template<class BasicTurbulenceModel>
void twoPhaseMixtureERzetaF<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = Cmu_*sqr(k_)/epsilon_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}
*/

template<class BasicTurbulenceModel>
twoPhaseMixtureERzetaF<BasicTurbulenceModel>&
twoPhaseMixtureERzetaF<BasicTurbulenceModel>::liquidTurbulence() const
{
    if (!liquidTurbulencePtr_)
    {
        const volVectorField& U = this->U_;

        const transportModel& gas = this->transport();
        const twoPhaseSystem& fluid =
            refCast<const twoPhaseSystem>(gas.fluid());
        const transportModel& liquid = fluid.otherPhase(gas);

        liquidTurbulencePtr_ =
           &const_cast<twoPhaseMixtureERzetaF<BasicTurbulenceModel>&>
            (
                U.db().lookupObject<twoPhaseMixtureERzetaF<BasicTurbulenceModel>>
                (
                    IOobject::groupName
                    (
                        turbulenceModel::propertiesName,
                        liquid.name()
                    )
                )
            );
    }

    return *liquidTurbulencePtr_;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> twoPhaseMixtureERzetaF<BasicTurbulenceModel>::Ct2() const
{
    const twoPhaseMixtureERzetaF<BasicTurbulenceModel>& liquidTurbulence =
        this->liquidTurbulence();

    const transportModel& gas = this->transport();
    const twoPhaseSystem& fluid = refCast<const twoPhaseSystem>(gas.fluid());
    const transportModel& liquid = fluid.otherPhase(gas);

    const volScalarField& alphag = this->alpha_;

    volScalarField magUr(mag(liquidTurbulence.U() - this->U()));

    volScalarField beta
    (
        (6*this->Cmu_/(4*sqrt(3.0/2.0)))
       *fluid.Kd()/liquid.rho()
       *(liquidTurbulence.k_/liquidTurbulence.epsilon_)
    );
    volScalarField Ct0((3 + beta)/(1 + beta + 2*gas.rho()/liquid.rho()));
    volScalarField fAlphad((180 + (-4.71e3 + 4.26e4*alphag)*alphag)*alphag);

    return sqr(1 + (Ct0 - 1)*exp(-fAlphad));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> twoPhaseMixtureERzetaF<BasicTurbulenceModel>::rholEff() const
{
    const transportModel& gas = this->transport();
    const twoPhaseSystem& fluid = refCast<const twoPhaseSystem>(gas.fluid());
    return fluid.otherPhase(gas).rho();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> twoPhaseMixtureERzetaF<BasicTurbulenceModel>::rhogEff() const
{
    const transportModel& gas = this->transport();
    const twoPhaseSystem& fluid = refCast<const twoPhaseSystem>(gas.fluid());
    const virtualMassModel& virtualMass =
        fluid.lookupSubModel<virtualMassModel>(gas, fluid.otherPhase(gas));
    return
        gas.rho()
      + virtualMass.Cvm()*fluid.otherPhase(gas).rho();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> twoPhaseMixtureERzetaF<BasicTurbulenceModel>::rhom() const
{
    const volScalarField& alphag = this->alpha_;
    const volScalarField& alphal = this->liquidTurbulence().alpha_;

    return alphal*rholEff() + alphag*rhogEff();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> twoPhaseMixtureERzetaF<BasicTurbulenceModel>::mix
(
    const volScalarField& fc,
    const volScalarField& fd
) const
{
    const volScalarField& alphag = this->alpha_;
    const volScalarField& alphal = this->liquidTurbulence().alpha_;

    return (alphal*rholEff()*fc + alphag*rhogEff()*fd)/rhom_();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> twoPhaseMixtureERzetaF<BasicTurbulenceModel>::mixU
(
    const volScalarField& fc,
    const volScalarField& fd
) const
{
    const volScalarField& alphag = this->alpha_;
    const volScalarField& alphal = this->liquidTurbulence().alpha_;

    return
        (alphal*rholEff()*fc + alphag*rhogEff()*Ct2_()*fd)
       /(alphal*rholEff() + alphag*rhogEff()*Ct2_());
}


template<class BasicTurbulenceModel>
tmp<surfaceScalarField> twoPhaseMixtureERzetaF<BasicTurbulenceModel>::mixFlux
(
    const surfaceScalarField& fc,
    const surfaceScalarField& fd
) const
{
    const volScalarField& alphag = this->alpha_;
    const volScalarField& alphal = this->liquidTurbulence().alpha_;

    surfaceScalarField alphalf(fvc::interpolate(alphal));
    surfaceScalarField alphagf(fvc::interpolate(alphag));

    surfaceScalarField rholEfff(fvc::interpolate(rholEff()));
    surfaceScalarField rhogEfff(fvc::interpolate(rhogEff()));

    return
       (alphalf*rholEfff*fc + alphagf*rhogEfff*fvc::interpolate(Ct2_())*fd)
      /(alphalf*rholEfff + alphagf*rhogEfff*fvc::interpolate(Ct2_()));
}

/*
template<class BasicTurbulenceModel>
tmp<volScalarField> twoPhaseMixtureERzetaF<BasicTurbulenceModel>::bubbleG() const
{
    const twoPhaseMixtureERzetaF<BasicTurbulenceModel>& liquidTurbulence =
        this->liquidTurbulence();

    const transportModel& gas = this->transport();
    const twoPhaseSystem& fluid = refCast<const twoPhaseSystem>(gas.fluid());
    const transportModel& liquid = fluid.otherPhase(gas);

    const dragModel& drag = fluid.lookupSubModel<dragModel>(gas, liquid);

    volScalarField magUr(mag(liquidTurbulence.U() - this->U()));

    // Lahey model
    tmp<volScalarField> bubbleG
    (
        Cp_
       *liquid*liquid.rho()
       *(
            pow3(magUr)
          + pow(drag.CdRe()*liquid.nu()/gas.d(), 4.0/3.0)
           *pow(magUr, 5.0/3.0)
        )
       *gas
       /gas.d()
    );

    // Simple model
    // tmp<volScalarField> bubbleG
    // (
    //     Cp_*liquid*drag.K()*sqr(magUr)
    // );

    return bubbleG;
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> twoPhaseMixtureERzetaF<BasicTurbulenceModel>::kSource() const
{
    return fvm::Su(bubbleG()/rhom_(), km_());
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> twoPhaseMixtureERzetaF<BasicTurbulenceModel>::epsilonSource() const
{
    return fvm::Su(C3_*epsilonm_()*bubbleG()/(rhom_()*km_()), epsilonm_());
}

*/

template<class BasicTurbulenceModel>
void twoPhaseMixtureERzetaF<BasicTurbulenceModel>::correct()
{
    const transportModel& gas = this->transport();
    const twoPhaseSystem& fluid = refCast<const twoPhaseSystem>(gas.fluid());

    // Only solve the mixture turbulence for the gas-phase
    if (&gas != &fluid.phase1())
    {
        // This is the liquid phase but check the model for the gas-phase
        // is consistent
        this->liquidTurbulence();

        return;
    }

    if (!this->turbulence_)
    {
        return;
    }

    // Initialise the mixture fields if they have not yet been constructed
    initMixtureFields();

    // Local references to gas-phase properties
    tmp<surfaceScalarField> phig = this->phi();
    const volVectorField& Ug = this->U_;
    const volScalarField& alphag = this->alpha_;
    volScalarField& kg = this->k_;
    volScalarField& omegag = this->omega_;
    volScalarField& zetag = this->zeta_;
    volScalarField& fg = this->f_;
    //volScalarField& epsilong = this->epsilon_;
    volScalarField& nutg = this->nut_;

    // Local references to liquid-phase properties
    twoPhaseMixtureERzetaF<BasicTurbulenceModel>& liquidTurbulence = this->liquidTurbulence();
    tmp<surfaceScalarField> phil = liquidTurbulence.phi();
    const volVectorField& Ul = liquidTurbulence.U_;
    const volScalarField& alphal = liquidTurbulence.alpha_;
    volScalarField& kl = liquidTurbulence.k_;
    volScalarField& omegal = liquidTurbulence.omega_;
    volScalarField& zetal = liquidTurbulence.zeta_;
    volScalarField& fl = liquidTurbulence.f_;
    //volScalarField& epsilonl = liquidTurbulence.epsilon_;
    volScalarField& nutl = liquidTurbulence.nut_;

    // Local references to mixture properties
    volScalarField& rhom = rhom_();
    volScalarField& km = km_();
    //volScalarField& epsilonm = epsilonm_();

    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    // Update the effective mixture density
    rhom = this->rhom();

    // Mixture flux
    surfaceScalarField phim("phim", mixFlux(phil, phig));

    // Mixture velocity divergence
    volScalarField divUm
    (
        mixU
        (
            fvc::div(fvc::absolute(phil, Ul)),
            fvc::div(fvc::absolute(phig, Ug))
        )
    );

    const volTensorField gradUl(fvc::grad(Ul));
    const volScalarField S2l(2*magSqr(dev(symm(gradUl))));
    const volScalarField Gl(this->GName(), nutl*S2l);
    /*
    const volScalarField Tl(this->Tau());
    const volScalarField Ll(this->L());



    /*

    tmp<volScalarField> Gc;
    {
        tmp<volTensorField> tgradUl = fvc::grad(Ul);
        Gc = tmp<volScalarField>
        (
            new volScalarField
            (
                this->GName(),
                nutl*(tgradUl() && dev(twoSymm(tgradUl())))
            )
        );
        tgradUl.clear();

        // Update k, epsilon and G at the wall
        kl.boundaryFieldRef().updateCoeffs();
        epsilonl.boundaryFieldRef().updateCoeffs();

        Gc.ref().checkOut();
    }

    tmp<volScalarField> Gd;
    {
        tmp<volTensorField> tgradUg = fvc::grad(Ug);
        Gd = tmp<volScalarField>
        (
            new volScalarField
            (
                this->GName(),
                nutg*(tgradUg() && dev(twoSymm(tgradUg())))
            )
        );
        tgradUg.clear();

        // Update k, epsilon and G at the wall
        kg.boundaryFieldRef().updateCoeffs();
        epsilong.boundaryFieldRef().updateCoeffs();

        Gd.ref().checkOut();
    }

    // Mixture turbulence generation
    volScalarField Gm(mix(Gc, Gd));

    // Mixture turbulence viscosity
    volScalarField nutm(mixU(nutl, nutg));

    // Update the mixture k and epsilon boundary conditions
    km == mix(kl, kg);
    bound(km, this->kMin_);
    epsilonm == mix(epsilonl, epsilong);
    bound(epsilonm, this->epsilonMin_);

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilonm)
      + fvm::div(phim, epsilonm)
      - fvm::Sp(fvc::div(phim), epsilonm)
      - fvm::laplacian(DepsilonEff(nutm), epsilonm)
     ==
        C1_*Gm*epsilonm/km
      - fvm::SuSp(((2.0/3.0)*C1_)*divUm, epsilonm)
      - fvm::Sp(C2_*epsilonm/km, epsilonm)
      + epsilonSource()
      + fvOptions(epsilonm)
    );

    epsEqn.ref().relax();
    fvOptions.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilonm.boundaryFieldRef());
    solve(epsEqn);
    fvOptions.correct(epsilonm);
    bound(epsilonm, this->epsilonMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kmEqn
    (
        fvm::ddt(km)
      + fvm::div(phim, km)
      - fvm::Sp(fvc::div(phim), km)
      - fvm::laplacian(DkEff(nutm), km)
     ==
        Gm
      - fvm::SuSp((2.0/3.0)*divUm, km)
      - fvm::Sp(epsilonm/km, km)
      + kSource()
      + fvOptions(km)
    );

    kmEqn.ref().relax();
    fvOptions.constrain(kmEqn.ref());
    solve(kmEqn);
    fvOptions.correct(km);
    bound(km, this->kMin_);
    km.correctBoundaryConditions();

    volScalarField Cc2(rhom/(alphal*rholEff() + alphag*rhogEff()*Ct2_()));
    kl = Cc2*km;
    kl.correctBoundaryConditions();
    epsilonl = Cc2*epsilonm;
    epsilonl.correctBoundaryConditions();
    liquidTurbulence.correctNut();

    Ct2_() = Ct2();
    kg = Ct2_()*kl;
    kg.correctBoundaryConditions();
    epsilong = Ct2_()*epsilonl;
    epsilong.correctBoundaryConditions();
    nutg = Ct2_()*(liquidTurbulence.nu()/this->nu())*nutl;
    */
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
