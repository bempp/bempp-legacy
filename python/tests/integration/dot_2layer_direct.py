import sys
sys.path.append("../..")


from bempp import lib as blib
#from bempp import visualization as vis
import numpy as np
import tempfile
import os
import subprocess
import math

def evalBoundaryData(point):
    return 1

def evalNullData(point):
    return 0

def Keijzer(n):
    th = math.asin(1.0/n)
    costh = math.fabs(math.cos(th))
    R0 = ((n-1.0)*(n-1.0)) / ((n+1.0)*(n+1.0))
    return (2.0/(1.0-R0) - 1.0 + costh*costh*costh) / (1.0 - costh*costh)

# Define physical parameters

c = .3
freq = 100e6
omega = 2*np.pi*freq*1E-12
refind = 1.4
alpha = Keijzer(refind)

# Outer region

mua1 = .01
mus1 = 1.
kappa1 = 1./(3.*(mua1+mus1))
w1 = np.sqrt(mua1/kappa1+1j*omega/(c*kappa1))

# Inner region
mua2 = .02
mus2 = .5
kappa2 = 1./(3.*(mua2+mus2))
w2 = np.sqrt(mua2/kappa2+1j*omega/(c*kappa2))


# We consider two spheres. One has radius r1, the other radius r2<r1. We want to look at
# the low-rank interaction between the spheres. Gmsh is used to create the meshes

r1 = 25.
r2 = 7.5
element_size1 = 1.
element_size2 = element_size1*r2/r1
gmsh_command = "gmsh"
sphere_definition = "../../../examples/meshes/sphere.txt"

sphere_def = open(sphere_definition,'r').read()

# Construct two Gmsh meshes with the required parameters

s1_geo, s1_geo_name = tempfile.mkstemp(suffix='.geo',dir=os.getcwd(),text=True)
s2_geo, s2_geo_name = tempfile.mkstemp(suffix='.geo',dir=os.getcwd(),text=True)
s1_msh_name = os.path.splitext(s1_geo_name)[0]+".msh"
s2_msh_name = os.path.splitext(s2_geo_name)[0]+".msh"

s1_geo_f = os.fdopen(s1_geo,"w")
s2_geo_f = os.fdopen(s2_geo,"w")

s1_geo_f.write("rad = "+str(r1)+";\nlc = "+str(element_size1)+";\n"+sphere_def)
s2_geo_f.write("rad = "+str(r2)+";\nlc = "+str(element_size2)+";\n"+sphere_def)
s1_geo_f.close()
s2_geo_f.close()

# Use Gmsh to create meshes

subprocess.check_call(gmsh_command+" -2 "+s1_geo_name,shell=True)
subprocess.check_call(gmsh_command+" -2 "+s2_geo_name,shell=True)

# Read the meshes into BEM++ Objects

sphere1 = blib.createGridFactory().importGmshGrid("triangular",s1_msh_name)
sphere2 = blib.createGridFactory().importGmshGrid("triangular",s2_msh_name)

# Clean up the temporary files

os.remove(s1_geo_name)
os.remove(s2_geo_name)
os.remove(s1_msh_name)
os.remove(s2_msh_name)

# Create Context

accuracy_options = blib.createAccuracyOptions()
# 1 orders higher than default accuracy for regular integrals
accuracy_options.doubleRegular.setRelativeQuadratureOrder(1)
# 0 orders higher than default accuracy for regular integrals
accuracy_options.doubleSingular.setRelativeQuadratureOrder(0)
strategy = blib.createNumericalQuadratureStrategy("float64", "complex128", accuracy_options)
options = blib.createAssemblyOptions()
aca_options = blib.createAcaOptions()
aca_options.eps=1E-6
options.switchToAca(aca_options)
context = blib.createContext(strategy, options)


# Create the spaces

sphere1_plc = blib.createPiecewiseLinearContinuousScalarSpace(context,sphere1)
sphere2_plc = blib.createPiecewiseLinearContinuousScalarSpace(context,sphere2)


# Now create the operators
slp11 = blib.createModifiedHelmholtz3dSingleLayerBoundaryOperator(context,sphere1_plc,sphere1_plc,sphere1_plc,w1)
dlp11 = blib.createModifiedHelmholtz3dDoubleLayerBoundaryOperator(context,sphere1_plc,sphere1_plc,sphere1_plc,w1)
id11  = blib.createIdentityOperator(context,sphere1_plc,sphere1_plc,sphere1_plc)

slp22_w1 = blib.createModifiedHelmholtz3dSingleLayerBoundaryOperator(context,sphere2_plc,sphere2_plc,sphere2_plc,w1)
dlp22_w1 = blib.createModifiedHelmholtz3dDoubleLayerBoundaryOperator(context,sphere2_plc,sphere2_plc,sphere2_plc,w1)
id22  = blib.createIdentityOperator(context,sphere2_plc,sphere2_plc,sphere2_plc)

slp22_w2 = blib.createModifiedHelmholtz3dSingleLayerBoundaryOperator(context,sphere2_plc,sphere2_plc,sphere2_plc,w2)
dlp22_w2 = blib.createModifiedHelmholtz3dDoubleLayerBoundaryOperator(context,sphere2_plc,sphere2_plc,sphere2_plc,w2)

slp12 = blib.createModifiedHelmholtz3dSingleLayerBoundaryOperator(context,sphere2_plc,sphere1_plc,sphere1_plc,w1)
dlp12 = blib.createModifiedHelmholtz3dDoubleLayerBoundaryOperator(context,sphere2_plc,sphere1_plc,sphere1_plc,w1)

slp21 = blib.createModifiedHelmholtz3dSingleLayerBoundaryOperator(context,sphere1_plc,sphere2_plc,sphere2_plc,w1)
dlp21 = blib.createModifiedHelmholtz3dDoubleLayerBoundaryOperator(context,sphere1_plc,sphere2_plc,sphere2_plc,w1)


scale = 1.0/(2.0*alpha*kappa1)
lhs_k11 = 0.5*id11 + dlp11 + scale*slp11
lhs_k12 = -1.0*dlp12
lhs_k13 = -(1.0/kappa1)*slp12
lhs_k21 = dlp21 + scale*slp21
lhs_k22 = 0.5*id22 - dlp22_w1
lhs_k23 = -(1.0/kappa1)*slp22_w1
# lhs_k31 -- empty
lhs_k32 = 0.5*id22 + dlp22_w2
lhs_k33 = (1.0/kappa2) * slp22_w2

structure = blib.createBlockedOperatorStructure(context)
structure.setBlock(0, 0, lhs_k11)
structure.setBlock(0, 1, lhs_k12)
structure.setBlock(0, 2, lhs_k13)
structure.setBlock(1, 0, lhs_k21)
structure.setBlock(1, 1, lhs_k22)
structure.setBlock(1, 2, lhs_k23)
# structure.setBlock(2, 0, ...); -- empty
structure.setBlock(2, 1, lhs_k32)
structure.setBlock(2, 2, lhs_k33)

lhsOp = blib.createBlockedBoundaryOperator(context,structure)

rhs1 = scale*slp11
rhs2 = scale*slp21

boundaryData1 = rhs1 * blib.createGridFunction(
    context, sphere1_plc, sphere1_plc, evalBoundaryData)
boundaryData2 = rhs2 * blib.createGridFunction(
    context, sphere1_plc, sphere1_plc, evalBoundaryData)
boundaryData3 = blib.createGridFunction(
    context, sphere2_plc, sphere2_plc, evalNullData)

rhs = [boundaryData1, boundaryData2, boundaryData3]

solver = blib.createDefaultIterativeSolver(lhsOp)
params = blib.defaultGmresParameterList(1e-8)
solver.initializeSolver(params)
solution = solver.solve(rhs)
u0 = solution.gridFunction(0)
u1 = solution.gridFunction(1)
v1 = solution.gridFunction(2)

# write out VTK files
u0.exportToVtk("vertex_data", "u0", "u0")
u1.exportToVtk("vertex_data", "u1", "u1")
v1.exportToVtk("vertex_data", "v1", "v1")
