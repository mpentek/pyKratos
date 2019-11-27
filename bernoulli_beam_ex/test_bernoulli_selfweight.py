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
    BODY_FORCE_X,
    BODY_FORCE_Y,
    EXTERNAL_FORCE_X,
    EXTERNAL_FORCE_Y,
    EXTERNAL_MOMENT,   
    LINE_LOAD_X,
    LINE_LOAD_Y,
]


#      1
#     / \
#   0    2       


node_list = {
    0: array([0.0, 0.0]), # Koordinaten in m
    1: array([3, -2]),
    2: array([6, 0]),
}

property_list = { 
    0: {YOUNG_MODULUS: 210e9,           # in N/m^2 (steel 210000 N/mm²) 
        SECTION_TYPE:  84.5e-4,         # in m^2   (IPE 400)     
        MOMENT_INERTIA_AREA: 23130e-8,  # in m^4 (Iy from IPE 400)
        DENSITY: 7850,                  # in kg/m^3 (steel)
        BODY_FORCE_Y: 9.81,             # in m/s^2  
        }
}

element_connectivities = {
    0: [0, [0, 1]],     ## 0 is property pointer
    1: [0, [1, 2]],     
}

condition_connectivities = {
    0: [0, [1] ],
}

# element_conditions = {
#     # 0: [0, [0]],
# }


nodal_values = {        #Randbedingungen
    DISPLACEMENT_X: [
        [0, True, 0.0], # first column is Id of nodeB, second col if fixity, third is imposed value ggf. vorverformung in [m]
        [2, True, 0.0],
    ],
    DISPLACEMENT_Y: [
        [0, True, 0.0],
        [2, True, 0.0],
    ],
    ROTATION: [
        [0, True, 0.0],
        [2, True, 0.0],    
    ],
    EXTERNAL_FORCE_X: [
        [1, True, 200000],    # in N
    ],
    EXTERNAL_FORCE_Y: [
        [1, True, 100000],    # in N
    ],
}

# ## Constant Line Loads
# element_values = {
#     LINE_LOAD_X: [
         
#     ],
#     LINE_LOAD_Y: [
#          #[0, True, 0],    # in N
#     ],
# }

### ERSTELLT Modell
buffer_size = 1  # store current step and 2 in the past         
model_part = model_part.ModelPart(buffer_size, solution_step_variables)
model_part.AddNodes(node_list)
model_part.AddProperties(property_list)
model_part.AddElements("bernoulli_beam", element_connectivities)
model_part.AddConditions("point_condition_2d", condition_connectivities)
# model_part.AddElementConditions("line_load_condition_2d", element_conditions)
model_part.AddNodalValues(nodal_values)
# model_part.AddElementValues(element_values)

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

strategy.Solve()

# zero_based_indices_for_nodes = True
# GiDIO = gid_io.GidIO("gid_out",zero_based_indices_for_nodes)
# GiDIO.WriteMesh(model_part,"outmesh")

# time = 0  
# GiDIO.WriteNodalResults(DISPLACEMENT,model_part.NodeIterators(), time)

import plot_system
scale = 1000
plot_system.PlotSystem(model_part.NodeIterators(), model_part.ElementIterators(), DISPLACEMENT, "plot_DISPLACEMENT.png", scale)

# import plot_internal_forces
# plot_internal_forces.PlotInternalForces(model_part.ElementIterators())

print('finished')