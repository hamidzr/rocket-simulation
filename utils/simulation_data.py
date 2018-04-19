import numpy as np
from utils.helpers import logger

# WARN TODO arrays start at 1 in matlab, array indices are refrenced w/ ()
class SimulationData:
  def __init__(self, dt, tMax):
    # self.dt = 0.01
    # self.t = []
    # self.Fd = []
    # self.Ft = []
    # self.m = []
    # self.theta = []
    # self.z = []
    # self.vz = []
    # self.az = []
    # self.x = []
    # self.vx = []
    # self.ax = []

    self.I = 1 # WARN changed 2 to 1

    self.dt = dt
    self.t = np.arange(0, tMax, dt) # 0:dt:tMax
    self.Fd = np.zeros(len(self.t))
    self.Ft = np.zeros(len(self.t))
    self.m = np.zeros(len(self.t))
    self.theta = np.zeros(len(self.t))
    self.z = np.zeros(len(self.t))
    self.vz = np.zeros(len(self.t))
    self.az = np.zeros(len(self.t))
    self.x = np.zeros(len(self.t))
    self.vx = np.zeros(len(self.t))
    self.ax = np.zeros(len(self.t))

  def SetInitialConditions(self, mFull):
    self.z[0] = 0
    self.vz[0] = 0
    self.az[0] = -9.81
    self.x[0] = 0
    self.vx[0] = 0
    self.ax[0] = 0
    self.m[0] = mFull; # Start full
    self.Fd[0] = 0
    self.Ft[0] = 0

  # TODO what does this do?
  def ResetFrame(self):
    self.I = 1 # WARN changed 2 to 1

  def Thrust(self, mdot, vjet, Pe, Pa, Ae, N):
    self.m[self.I] = self.m[self.I-1] - mdot * self.dt * N
    self.Ft[self.I] = (mdot * vjet + (Pe - Pa) * Ae) * N


  def NoThrust(self):
    self.m[self.I] = self.m[self.I-1]
    self.Ft[self.I] = 0


  def RelativeWind(self):
    zpart = abs(self.vz[self.I-1])
    xpart = abs(self.vx[self.I-1])
    relWind = abs(np.cos(self.theta[self.I-1]) * zpart) + abs(np.sin(self.theta[self.I-1]) * xpart)
    return relWind

  def UpdateVZ(self, dv_axial, g):
    vz =  self.vz[self.I-1] + dv_axial * np.cos(self.theta[self.I-1]) - g*self.dt
    self.vz[self.I] = vz


  def UpdateAngle(self, dv_axial, thetaWindow, thetaAngle):
    # Migrate to the angle over first few seconds
    if self.t[self.I] < thetaWindow:
      self.theta[self.I] = thetaAngle * (self.t[self.I] / thetaWindow)
      self.vx[self.I] = np.tan(self.theta[self.I-1]) * self.vz[self.I-1]
    else:
      self.vx[self.I] = self.vx[self.I-1] + dv_axial * np.sin(self.theta[self.I-1])
      self.theta[self.I] = np.arctan(self.vx[self.I]/self.vz[self.I]) # WARN atan

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
    t1 = np.power(self.vz[self.I + offset], 2) # WARN changes
    t2 = np.power(self.vx[self.I + offset], 2)
    v = t1 + np.sqrt(t2)
    return v

