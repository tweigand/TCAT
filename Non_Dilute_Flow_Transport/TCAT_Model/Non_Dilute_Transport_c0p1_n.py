from proteus import *
from proteus.default_n import *
from Non_Dilute_Transport_p import *
from Non_Dilute_Transport import *

DT = T/100.
nDTout = 100

timeIntegration = BackwardEuler
#timeIntegration = FLCBDF
#stepController = FLCBDF_controller
#systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep

#rtol_u[0] = 1.0e-6
#atol_u[0] = 1.0e-6
#rtol_u[1] = 1.0e-6
#atol_u[1] = 1.0e-6
#rtol_u[2] = 1.0e-6
#atol_u[2] = 1.0e-6

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,1:C0_AffineLinearOnSimplexWithNodalBasis,2:C0_AffineLinearOnSimplexWithNodalBasis}
#femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis,1:DG_AffineLinearOnSimplexWithNodalBasis,2:DG_AffineLinearOnSimplexWithNodalBasis}
elementQuadrature = SimplexGaussQuadrature(nd,4)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

#numericalFluxType = Advection_DiagonalUpwind

nn=201
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