class Atmosphere:
  def __init__(self):
    T0_arr = [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65]
    p0_arr = [101325, 22632.1, 5474.89, 868.02, 110.91, 66.94, 3.96]
    B_arr = [-0.0065, 0, 0.001, 0.0028, 0.0, -0.0028, -0.002]
    z_arr = [0, 11000, 20000, 32000, 47000, 51000, 71000]
    g0 = 9.81
    R0 = 6370
    pMin = 0.000000001
    Tmin = 3
    R = 287

  def _Dynamics(self, z):
    idx = max(find(z >= self.z_arr)) #TODO
    T0 = self.T0_arr[idx]
    p0 = self.p0_arr[idx]
    B = self.B_arr[idx]
    z0 = self.z_arr[idx]
    return [T0, p0, B, z0]

  def ParametersAtAltitude(self, z):
    # Get the dynamics variables
    [T0, p0, B, z0] = self._Dynamics(z)

    # Now calculate gravity
    g = math.pow(self.g0 * (self.R0 / (self.R0 + z/1000)), 2)

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
        p = p0 * math.pow((T0 / T), (g / (self.R * B)))
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

