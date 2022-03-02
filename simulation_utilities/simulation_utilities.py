# The idea here is to optimize the density of the sphere such that given:
# - the touching area (defined as flat)
# - the mechanical properties and terrestrial acceleration
# we can find a solution such that
# - the stress at the edge of the touching area has an average z-component of zero
#
# The last condition is the natural condition for a flat continuation of the surface, that is, the object at equilibrium
# on its touching area


from ngsolve import *
from ngsolve.bvp import BVP
import netgen.meshing as nm
import math

def epsilon(u):
    """Strain from deformation field, symmetric"""
    return  0.5*(Grad(u)+Grad(u).trans)

def sigma(strain,E,nu):
    mu = E / 2 / (1 + nu)
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    return lam*Trace(strain)*Id(round(math.sqrt(strain.dim))) + 2*mu*strain




def getElasticBilinearForm(fes,E,nu,compression):
    a = BilinearForm(fes, symmetric=True, condense=True)
    v = fes.TestFunction()
    u = fes.TrialFunction()
    e=epsilon(v)
    s=sigma(epsilon(u),E,nu)
    a+=(s[0]*e[0]+s[1]*e[1]+s[2]*e[2]+s[3]*e[3]+s[4]*e[4]+s[5]*e[5]+s[6]*e[6]+s[7]*e[7]+s[8]*e[8])*dx
    return a

# This is like the additional thermal expansion element described in this paper
    # IEEE TRANSACTIONS ON MEDICAL IMAGING, VOL. 25, NO. 11, NOVEMBER 2006, https://pubmed.ncbi.nlm.nih.gov/17117771/
def GetIsotropicCompressionLinearContribution(fes, E,nu,compression):
    v = fes.TestFunction()
    e=-Id(v.dim) * compression
    s=sigma(epsilon(v),E,nu)
    return (s[0]*e[0]+s[1]*e[1]+s[2]*e[2]+s[3]*e[3]+s[4]*e[4]+s[5]*e[5]+s[6]*e[6]+s[7]*e[7]+s[8]*e[8])*dx







