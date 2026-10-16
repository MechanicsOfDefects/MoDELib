# Python mobility component (screw) for W-xRe: Gaussian-process regression of log10(velocity)
# on temperature, resolved shear stress and Re content.
#
# Material file:  mobility_a/2<111>{110}s=python ../../../python/MobilityLaws/WRe_GPR/DDMobilityPy_WRe_GPR.py <Re at.%>;
# Interface:      MobilitySolver().velocityPy(stress[GPa], burgers, normal, temp[K], dL, dt, soluteContent[at.%])
# Returns:        magnitude of the dislocation velocity in m/s

import os
import numpy as np
import pandas as pd
import joblib

_here = os.path.dirname(os.path.abspath(__file__))
loaded_model = joblib.load(os.path.join(_here, 'WRe_GPR_3features_log10Velocity_Matern.joblib'))

class MobilitySolver():
        def velocityPy(self, stress, burgers, normal, temp, dL, dt, soluteContent):
            stress = np.asarray(stress, dtype=float)
            burgers = np.asarray(burgers, dtype=float).reshape(-1)
            normal = np.asarray(normal, dtype=float).reshape(-1)

            slip_direction = burgers / np.linalg.norm(burgers)  # Unit slip direction
            plane_normal = normal / np.linalg.norm(normal)      # Unit plane normal

            tau = float(slip_direction @ stress @ plane_normal) # Resolved shear stress: tau = b_hat . stress . n
            tau_model = abs(tau)                                # Positive shear-stress magnitudes

            features = pd.DataFrame(
                {
                    "temp": [float(temp)],
                    "tau_yz": [tau_model],
                    "solute_percent": [float(soluteContent)],
                }
            )

            log10_velocity = float(loaded_model.predict(features)[0])
            velocity_mps = 10.0**log10_velocity
            return velocity_mps
