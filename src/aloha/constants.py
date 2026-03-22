import numpy as np
from scipy.constants import c, e, epsilon_0, m_e, mu_0, pi

# Vacuum impedance and admittance
Z0 = np.sqrt(mu_0 / epsilon_0)
Y0 = 1 / Z0
