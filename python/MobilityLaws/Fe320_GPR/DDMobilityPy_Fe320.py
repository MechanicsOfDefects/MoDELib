# Full Python mobility for BCC Fe at 320 K: Gaussian-process regression trained on MD
# (8 features: sorted stress eigenvalues, |projection| of their eigenvectors on the line
# tangent, temperature, character angle; the model predicts sqrt(velocity)).
#
# Material file:  mobility_a/2<111>{110}=python ../../../python/MobilityLaws/Fe320_GPR/DDMobilityPy_Fe320.py;
# Interface:      MobilitySolver().velocityPy(stress[GPa], burgers, tangent, normal, temp[K], dL, dt)
# Returns:        magnitude of the dislocation velocity in m/s

import os
import numpy as np
import pandas as pd
import joblib

_here = os.path.dirname(os.path.abspath(__file__))
loaded_model = joblib.load(os.path.join(_here, 'Fe320_GPR_8features_sqrtVelocity_RBF_1001.joblib'))

def sort_eigenvalues_and_vectors(eigenvalues, eigenvectors):
    idx = eigenvalues.argsort()[::-1]  # Sort in descending order
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    return eigenvalues, eigenvectors

def angle_atan2(b_unit, n_unit, t_unit, f_unit):
    # Compute the angle using atan2 for better numerical stability
    cos_theta = np.sum(b_unit*t_unit, axis=-1)
    sin_theta = np.sum(np.cross(b_unit, t_unit)*n_unit, axis=-1)
    th_deg = np.degrees(np.arctan2(sin_theta, cos_theta))
    return (th_deg % 180 + 180) % 180  # Ensure angle is in [0, 180)

class MobilitySolver():
        def velocityPy(self,stress,burgers,tangent,normal,temp,dL,dt):
                # Character angle in degrees
                f_pk = np.cross(stress@burgers, tangent)
                f_g = f_pk - np.dot(f_pk,normal)*normal
                f_g_norm = f_g / np.linalg.norm(f_g)
                miss_angle = angle_atan2(burgers, normal, tangent, f_g_norm)

                # Eigenvalues and eigenvectors, project eigenvectors along tangent, absolute values
                eigValue, eigVector = np.linalg.eigh(0.5*(stress+stress.T)) # stress is symmetric: real eigenvalues
                eigValue, eigVector = sort_eigenvalues_and_vectors(eigValue, eigVector)
                projection = eigVector.T @ tangent
                anev = np.abs(projection)

                # Feature set
                features = pd.DataFrame({
                    'e1': eigValue[0],
                    'e2': eigValue[1],
                    'e3': eigValue[2],
                    'anev1': anev[0],
                    'anev2': anev[1],
                    'anev3': anev[2],
                    'temp': temp,
                    'miss': miss_angle
                }, index=[0])

                velocity_pred = loaded_model.predict(features)
                velocity_mps = float(velocity_pred[0])**2

                screw_char = abs(np.dot(burgers, tangent))  # 1 = screw, 0 = edge
                screw_cutoff = np.cos(np.radians(10.0))     # only within 10 deg of screw
                if screw_char >= screw_cutoff:
                    velocity_mps *= 0.5

                return velocity_mps
