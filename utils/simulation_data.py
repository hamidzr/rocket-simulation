import numpy as np

# WARN TODO arrays start at 1 in matlab, array indices are refrenced w/ ()
class SimulationData:
  def __init__(self, dt, tMax):
    self.dt = 0.01
    self.t = []
    self.Fd = []
    self.Ft = []
    self.m = []
    self.theta = []
    self.z = []
    self.vz = []
    self.az = []
    self.x = []
    self.vx = []
    self.ax = []
    self.I = 2

    self.dt = dt
    self.t = range(0, tMax, dt) # 0:dt:tMax
    self.Fd = [0] * len(self.t)
    self.Ft = [0] * len(self.t)
    self.m = [0] * len(self.t)
    self.theta = [0] * len(self.t)
    self.z = [0] * len(self.t)
    self.vz = [0] * len(self.t)
    self.az = [0] * len(self.t)
    self.x = [0] * len(self.t)
    self.vx = [0] * len(self.t)
    self.ax = [0] * len(self.t)

    def SetInitialConditions(self, mFull):
      self.z[1] = 0
      self.vz[1] = 0
      self.az[1] = -9.81
      self.x[1] = 0
      self.vx[1] = 0
      self.ax[1] = 0
      self.m[1] = mFull; # Start full
      self.Fd[1] = 0
      self.Ft[1] = 0

    def ResetFrame(self):
      self.I = 2

    def Thrust(self, mdot, vjet, Pe, Pa, Ae, N):
      self.m[self.I] = self.m[self.I-1] - mdot * self.dt * N
      self.Ft[self.I] = (mdot * vjet + (Pe - Pa) * Ae) * N


    def NoThrust(self):
      self.m[self.I] = self.m[self.I-1]
      self.Ft[self.I] = 0


    def RelativeWind(self):
      zpart = abs(self.vz[self.I-1])
      xpart = abs(self.vx[self.I-1])
      relWind = abs(cos(self.theta[self.I-1]) * zpart) + abs(sin(self.theta[self.I-1]) * xpart)
      return relWind

    def UpdateVZ(self, dv_axial, g):
      self.vz[self.I] = self.vz[self.I-1] + dv_axial * cos(self.theta[self.I-1]) - g*self.dt


    def UpdateAngle(self, dv_axial, thetaWindow, thetaAngle):
      # Migrate to the angle over first few seconds
      if self.t[self.I] < thetaWindow:
        self.theta[self.I] = thetaAngle * (self.t[self.I] / thetaWindow)
        self.vx[self.I] = tan(self.theta[self.I-1]) * self.vz[self.I-1]
      else:
        self.vx[self.I] = self.vx[self.I-1] + dv_axial * sin(self.theta[self.I-1])
        self.theta[self.I] = atan(self.vx[self.I]/self.vz[self.I])

    def CalculateNewXandZ(self):
      # Net acceleration
      dvx = self.vx[self.I] - self.vx[self.I-1]
      self.ax[self.I] = dvx / self.dt
      self.x[self.I] = self.x[self.I-1] + 0.5 * dvx * self.dt + self.vx[self.I-1] * self.dt
      dvz = self.vz[self.I] - self.vz[self.I-1]
      self.az[self.I] = dvz / self.dt
      self.z[self.I] = self.z[self.I-1] + 0.5 * dvz * self.dt + self.vz[self.I-1] * self.dt

    def NextFrame(self):
      self.I = self.I + 1

    def RemainingMass(self, m0): #TODO
      mass = self.m[self.I] - m0
      return mass

    def velocity(self, offset):
      # v = (self.vz[self.I + offset].^2 + self.vx[self.I + offset].^2).^0.5
      t1 = math.pow(self.vz[self.I + offset], 2) # WARN changes
      t2 = math.pow(self.vx[self.I + offset], 2)
      v = t1 + math.sqrt(t2)
      return v

