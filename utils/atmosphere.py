import numpy as np
from utils.helpers import logger, memoize, filter, isreal, dump, load, plot_attempts
import os.path
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import multiprocessing
num_cores = multiprocessing.cpu_count()


# find index of the last item in z_arr that is bigger than is 'le' to z
arr_len = None
def theFind(arr, z):
  global arr_len
  if (not arr_len): arr_len = len(arr)
  for idx in range(arr_len-1, -1, -1):
    if arr[idx] <= z:
      return idx
  assert False, 'it should never reach this line'

class Atmosphere:
  def __init__(self):
    self.T0_arr = [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65]
    self.p0_arr = [101325, 22632.1, 5474.89, 868.02, 110.91, 66.94, 3.96]
    self.B_arr = [-0.0065, 0, 0.001, 0.0028, 0.0, -0.0028, -0.002]
    self.z_arr = [0, 11000, 20000, 32000, 47000, 51000, 71000]
    self.g0 = 9.81
    self.R0 = 6370
    self.pMin = 0.000000001
    self.Tmin = 3
    self.R = 287

  def _Dynamics(self, z):
    # idx = max(find(z >= self.z_arr)) # find index of the last item in z_arr that is bigger than is 'le' to z
    idx = theFind(self.z_arr, z)
    T0 = self.T0_arr[idx]
    p0 = self.p0_arr[idx]
    B = self.B_arr[idx]
    z0 = self.z_arr[idx]
    return [T0, p0, B, z0]

  def ParametersAtAltitude(self, z):
    # Get the dynamics variables
    [T0, p0, B, z0] = self._Dynamics(z)

    # Now calculate gravity
    g = self.g0 * ((self.R0 / (self.R0 + z/1000)) ** 2)

    # Above critical alt?
    if z > 180000:
      p = self.pMin
      T = self.Tmin
    else:
      # First calculate temperature
      T = T0 + B * (z - z0)
      if T < self.Tmin: #or not isreal(T): # WARN
        T = self.Tmin
      # end

      # Now calculate pressure
      if B == 0:
        p = p0 * np.exp(-g * (z - z0) / (self.R * T0))
      else:
        p = p0 * ((T0 / T) ** (g / (self.R * B)))
      # end

      # Keep to a min
      if p < self.pMin: #or not isreal(p):
        p = self.pMin
      # end
    # end

    # Then calculate density
    rho = p / (self.R * T)
  # end
    return [p, rho, T, g]
# end

if __name__ == "__main__":
  # plot atmosphere vs alt
  # do we need to compute it?
  dt = 1e-3
  target_alt = 2e+5
  fname = f'out/atmosphere-dt{dt}'
  print('preparing the data..')
  if not os.path.isfile(fname):
    atm = Atmosphere()
    history = {
      'pa': [],
      'rho': [],
      't': [],
      'g': [],
    }
    count = target_alt / dt
    zs = np.linspace(1, 200000, count)
    values = Parallel(n_jobs=num_cores)(delayed(atm.ParametersAtAltitude)(z) for z in zs)
    for [Pa, rho, T, g] in values:
      history['pa'].append(Pa)
      history['rho'].append(rho)
      history['t'].append(T)
      history['g'].append(g)
    dump(fname, history)
  else:
    history = load(fname)

  print('plotting..')
  plot_attempts(history['pa'])
  plt.show()
