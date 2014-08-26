from proteus import *
from proteus.default_n import *
from Classic_Transport_p import *
from Classic_Transport import *

DT = T/100.
nDTout = 100
#timeIntegration = BackwardEuler_cfl
#DT = 0.1
#runCFL = 0.3
#timeIntegration = BackwardEuler
#timeIntegration = ForwardEuler
timeIntegration = FLCBDF
stepController = FLCBDF_controller
systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep

rtol_u[0] = 1.0e-6
atol_u[0] = 1.0e-6
rtol_u[1] = 1.0e-6
atol_u[1] = 1.0e-6

subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,lag=True)
shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=1.0,lag=True)

#numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,1:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,4)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

nn = 101
nLevels = 1

massLumping = True

multilevelNonlinearSolver  = Newton
levelNonlinearSolver = Newton
#maxNonlinearIts = 25
#maxLineSearches = 2
fullNewtonFlag = True
#nl_atol_res = 1.0e-6
#l_atol_res = 1.0e-9
#tolFac = 1.0e-6