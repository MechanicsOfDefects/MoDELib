# Example Python mobility for MoDELib.
#
# Full mobility (DislocationMobilityPy), selected in the material file with
#     mobility_<slipSystem>=python <path/to/this/file.py>;
# MoDELib constructs MobilitySolver() and calls velocityPy at each quadrature point.
#
# Inputs:  stress [GPa] (3x3), burgers, tangent, normal (3-vectors), temp [K], dL, dt (code units)
# Returns: magnitude of the dislocation velocity in m/s

class MobilitySolver:
    def velocityPy(self, stress, burgers, tangent, normal, temp, dL, dt):
        return 1.0 # m/s
