from proteus import *
from proteus.default_n import *
from Dilute_Transport_p import *
from Dilute_Transport import *

DT = T/500.
nDTout = 500

#timeIntegration = BackwardEuler
#timeIntegration = ForwardEuler
#timeIntegration = FLCBDF
#stepController = FLCBDF_controller
#systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep

#rtol_u[0] = 1.0e-3
#atol_u[0] = 1.0e-3
#rtol_u[1] = 1.0e-3
#atol_u[1] = 1.0e-3


femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,1:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,4)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

nn = 1001
nLevels = 1

massLumping = True

multilevelNonlinearSolver  = Newton
levelNonlinearSolver = Newton
maxNonlinearIts = 25
maxLineSearches = 2
fullNewtonFlag = True
nl_atol_res = 1.0e-6
l_atol_res = 1.0e-9
tolFac = 1.0e-6