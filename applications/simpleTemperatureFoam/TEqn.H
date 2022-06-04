
alphat = turbulence->nut()/Prt;
alphat.correctBoundaryConditions();

volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + alphat);

fvScalarMatrix TEqn
(
	fvm::div(phi, T) == fvm::laplacian(alphaEff, T)
);

TEqn.relax();
TEqn.solve();
