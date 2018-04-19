import numpy as np

def filter(pred, lst):
  # Filters lst and returns the matching indices
  matches = []
  for index in range(len(lst)-1, -1, -1):
    if pred(lst[index]):
      matches.append(index)
  return matches

# mimic matlab's isreal()
def isreal(arr):
  # if input is not a list
  return np.isreal(arr)
  # if inp is a list:
  # return all(it == True for it in np.isreal(2))

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
    # idx = max(find(z >= self.z_arr)) #TODO find index of the last item in z_arr that is bigger than is 'le' to z
    idx = max(filter(lambda x: z >= x, self.z_arr))
    T0 = self.T0_arr[idx]
    p0 = self.p0_arr[idx]
    B = self.B_arr[idx]
    z0 = self.z_arr[idx]
    return [T0, p0, B, z0]

  def ParametersAtAltitude(self, z):
    # Get the dynamics variables
    [T0, p0, B, z0] = self._Dynamics(z)

    # Now calculate gravity
    g = np.power(self.g0 * (self.R0 / (self.R0 + z/1000)), 2)

    # Above critical alt?
    if z > 180000:
      p = self.pMin
      T = self.Tmin
    else:
      # First calculate temperature
      T = T0 + B * (z - z0)
      if T < self.Tmin or  not isreal(T): #TODO
        T = self.Tmin
      # end

      # Now calculate pressure
      if B == 0:
        p = p0 * exp(-g * (z - z0) / (self.R * T0))
      else:
        p = p0 * np.power((T0 / T), (g / (self.R * B)))
      # end

      # Keep to a min
      if p < self.pMin or  not isreal(p):
        p = self.pMin
      # end
    # end

    # Then calculate density
    rho = p / (self.R * T)
  # end
    return [p, rho, T, g]
# end

