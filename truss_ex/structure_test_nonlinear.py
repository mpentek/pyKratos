from __future__ import print_function, absolute_import, division 
import sys
sys.path.append("..")
#print(sys.path)

from numpy import *
from pyKratos import *

# add variables to be allocated from the list in variables.py
solution_step_variables = [
    DISPLACEMENT_X,     #Verschiebung
    DISPLACEMENT_Y,
    ROTATION,
    IS_LAGRANGIAN,
    EXTERNAL_FORCE_X,
    EXTERNAL_FORCE_Y, 
    EXTERNAL_MOMENT,  
]

#       1 - 3
#     / |   |
#   0   2   4

property_list = { 
    0: {YOUNG_MODULUS: 210e9,      #in N/m^2
        SECTION_TYPE:  0.001885,    # in m^2 -> Rundstab d = 8mm
        }
}

node_list = {
    0: array([0.0, 0.0]),
    1: array([5, -5]),
    2: array([5, 0]),
    3: array([10, -5]),
    4: array([10, 0]),
}

element_connectivities = {
    1: [0, [0, 1]],     ## 0 is property pointer
    2: [0, [1, 2]],     ## über prop müssste qx übergeben werden
    3: [0, [1, 3]],
    4: [0, [3, 4]],
}

condition_connectivities = {
    1: [0, [1]],
    2: [0, [3]],
}


nodal_values = {        #Randbedingungen
DISPLACEMENT_X: [
        [0, True, 0.0], # first column is Id of node, second col if fixity, third is imposed vlaue ggf. vorverformung
        [2, True, 0.0],
        [4, True, 0.0],
    ],
DISPLACEMENT_Y: [
    [0, True, 0.0],
    [2, True, 0.0],
    [4, True, 0.0],
    ],
EXTERNAL_FORCE_Y: [
    [1, True, 40000],
    [3, True, 40000],
],
EXTERNAL_FORCE_X: [
    [1, True, 1000],
]
}

### ERSTELLT Modell
buffer_size = 1  # store current step and 2 in the past         
model_part = model_part.ModelPart(buffer_size, solution_step_variables)
model_part.AddNodes(node_list)
model_part.AddProperties(property_list)
model_part.AddElements("truss_element_nonlinear", element_connectivities)
model_part.AddConditions("point_condition_2d", condition_connectivities)
model_part.AddNodalValues(nodal_values)

print(model_part)

## Static Scheme
scheme = static_scheme.StaticScheme(model_part)

# Builder And Solver
builder_and_solver = builder_and_solver.BuilderAndSolver(
    model_part, scheme)

## Solving Strategy
strategy = solving_strategy.SolvingStrategy(
    model_part, scheme, builder_and_solver)
strategy.Initialize()

for i in range(0,2):    
    strategy.Solve()

# Output
# zero_based_indices_for_nodes = True
# GiDIO = gid_io.GidIO("gid_out",zero_based_indices_for_nodes)
# GiDIO.WriteMesh(model_part,"outmesh")

# time = 0  
# GiDIO.WriteNodalResults(DISPLACEMENT,model_part.NodeIterators(), time)

# import plot_contour
# plot_contour.PlotContour(model_part.NodeIterators(), DISPLACEMENT_X, "disp_DISPLACEMENT_X.png" )
# plot_contour.PlotContour(model_part.NodeIterators(), DISPLACEMENT_Y, "disp_DISPLACEMENT_Y.png" )

import plot_system
scale = 100
plot_system.PlotSystem(model_part.NodeIterators(), model_part.ElementIterators(), DISPLACEMENT, "plot_DISPLACEMENT.png", scale)

# import plot_internal_forces
# plot_internal_forces.PlotInternalForces(model_part.ElementIterators())




print('finished')