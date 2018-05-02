from utils.atmosphere import Atmosphere
from utils.simulation_data import SimulationData
from utils.helpers import logger, plot_attempts, dump, load, plot_batch_attempts, sort_two_lists
import numpy as np
import matplotlib.pyplot as plt
import argparse
import time

parser = argparse.ArgumentParser()
parser.add_argument("--dt", type=float, default=0.001, help="simulation tick")
parser.add_argument("--max_attempts", type=int, default=10, help="max number of attempts")
parser.add_argument("--save", type=str, default='out/default.out', help="save results to file")
parser.add_argument("--load", type=str, help="load from file insted of computing")
args = parser.parse_args()


## Rocket Model
# Grant Fennessy
# Hamid Zare


## Global parameters
# global atmosphere LandingVelocityAllowance TargetAltitude mSecondStage
# global m0_propellant mEmpty m0_full Cd_up Cd_down Pe Ae mdot vjet A
atmosphere = Atmosphere()
mSecondStage = 111500 # kg
m0_propellant = 433100 #kg
mEmpty = 22200 # kg
m0_full = m0_propellant + mSecondStage + mEmpty # kg
Cd_up = 0.3
Cd_down = 0.82
N = 9 # Number of engines
A = (3.7 ** 2) / 4 # Diameter of 3.7m
Thrust_vacuum = 8227000 # N
Thrust_ground = 7607000 # N
Isp_vacuum = 311 # s
Isp_ground = 282 # s
vjet = 3000
TargetAltitude = 200000 # m
errorAllowance = 0.0001 # target altitude error allowance
LandingVelocityAllowance = 6 # m/s


# Parameters being solved
MAX_LAUNCH_RATIO = 0.99999999 # should never happen, but prevents >= 1.
launch_ratio = 0.94
launch_ratio_BEST = launch_ratio
MaxDiscoveredVx = 0
takeoff_angle = 0.00034 # Start attempt angle
takeoff_angle_BEST = None
landingBurnAltitude = 25000 # Pick an abitrary start
landingBurnAltitude_BEST = None
I_apex_BEST = 0 # WARN was 1
I_land_BEST = 0

# Create simulation data
dt = args.dt
data = SimulationData(dt, 2000)

# keep track of experiments
attempts_history = [] # [(angles, alts), ...] len of maxAttempts
ratios = []

def SolveThrust(vjet, T_sea, T_vac, Isp_sea, Isp_vac, N):
  # Vacuum first
  mdot_vac = T_vac / (Isp_vac * 9.81)
  pe_Ae = T_vac - mdot_vac * vjet


  # Sea level
  mdot_sea = T_sea / (Isp_sea * 9.81)
  Ae = (T_sea - mdot_sea * vjet - pe_Ae) / -101325
  pe = pe_Ae / Ae

  # Per Engine
  Ae = Ae / N
  mdot = (mdot_vac + mdot_sea) / (2 * N)
  return [pe, Ae, mdot]

[Pe, Ae, mdot] = SolveThrust(vjet, Thrust_ground, Thrust_vacuum, Isp_ground, Isp_vacuum, N)

## Solve Takeoff Angle
# Responsible for solving the takeoff angle given a launch ratio.
# Will start at 0 to ensure the desired altitude can even be reached.
def SolveTakeoffAngle(launch_ratio, data):
  # Bring in globals
  global atmosphere, TargetAltitude, mSecondStage, A, takeoff_angle
  global m0_propellant, mEmpty, m0_full, Cd_up, Pe, Ae, mdot, vjet

  # Determine when to stop
  N = 9 # Use all engines
  stopMass = mEmpty + mSecondStage + m0_propellant * (1 - launch_ratio)
  data.SetInitialConditions(m0_full)

  # Iterate the times until there
  # use the best takeoff angle so far
  takeoff_angle = takeoff_angle_BEST if takeoff_angle_BEST else takeoff_angle
  adjustment =  0.000005 # Start adjustment
  adjustmentDir = 1
  attempt = 1

  # keep track of what values you tried
  angles = []
  errors = []

  while True: # Keep looping until terminated
    # Initial conditions
    data.theta[0] = 0

    print(f'testing angle {takeoff_angle}')
    # Loop until at peak altitude
    data.ResetFrame()
    while data.vz[data.I-1] >= 0:
      # Update atmosphere
      [Pa, rho, T, g] = atmosphere.ParametersAtAltitude(data.z[data.I-1])

      # Update mass
      if data.m[data.I-1] > stopMass:
        data.Thrust(mdot, vjet, Pe, Pa, Ae, N)
      else:
        data.NoThrust()
      # end

      # Angle of attack + Drag
      relWind = data.RelativeWind()
      # data.Fd[data.I] = 0.5 * Cd_up * A * rho * relWind.^2
      data.Fd[data.I] = 0.5 * Cd_up * A * rho * (relWind ** 2)

      # Axial -> X/Y
      a_axial = (data.Ft[data.I] - data.Fd[data.I]) / data.m[data.I]
      dv_axial = a_axial * data.dt
      data.UpdateVZ(dv_axial, g)

      # Update the angle
      data.UpdateAngle(dv_axial, 3, takeoff_angle)

      # Net acceleration
      data.CalculateNewXandZ()

      # Next frame
      data.NextFrame()
    # end

    # Rewind a frame
    data.I = data.I - 1

    # Evaluate
    attempt = attempt + 1
    error = abs(data.z[data.I] - TargetAltitude)

    angles.append(takeoff_angle)
    errors.append(error)

    if error <= TargetAltitude * errorAllowance:
      success = 1
      print("--- Error = {}m after {} attempts.\n".format(error, attempt))
      return [takeoff_angle, success, (angles, errors)]
    # end

    # logger.debug(data.z[data.I])
    # If didn't reach altitude...
    # logger.debug('z', data.z[0:10])
    # logger.debug('x', data.x[0:10])
    if data.z[data.I] < TargetAltitude:
      # If takeoff angle is 0, not enough fuel!
      # QUESTOIN why not set this at the begining and check this?
      if takeoff_angle == 0:
        success = 0
        print(f'wasted {attempt} to figureout there is not enough fuel')
        return [takeoff_angle, success, (angles, errors)]
      # end

      # Need to decrease the angle
      if adjustmentDir == 1:
        adjustment = adjustment / 2
      # end
      takeoff_angle = takeoff_angle - adjustment
      adjustmentDir = -1

    else:
      # Need to increase angle
      if adjustmentDir == -1:
        # If was moving other way, reduce adjustment amount
        adjustment = adjustment / 2
      # end
      takeoff_angle = takeoff_angle + adjustment
      adjustmentDir = 1
    # end
  # end
  return [takeoff_angle, success, (angles, errors)]
# end




## Solve Landing Burn
# Responsible for solving the takeoff angle given a launch ratio.
# Will start at 0 to ensure the desired altitude can even be reached.
def SolveLandingBurn(data, I_apex):
  # Bring in globals
  global atmosphere, LandingVelocityAllowance, A, landingBurnAltitude
  global mEmpty, Cd_down, Pe, Ae, mdot, vjet

  # Determine when to stop
  N = 3 # Just use three engines
  stopMass = mEmpty # Stop when no more propellant

  # Iterate the times until there
  # use the best takeoff angle so far
  landingBurnAltitude = landingBurnAltitude_BEST if landingBurnAltitude_BEST else landingBurnAltitude
  aAdjDir = 0
  aAdjVal = 2000
  adjRatio = 0.7
  SHUTDOWN_ALT_RATIO = 1.8 # how much should the altitude be changed based on the engine shutoff alt empirical

  # dir +1 or -1; constant is False or a value
  def adjust_landing_alt(direction, constant=False, reason=None):
    nonlocal aAdjDir, aAdjVal, adjRatio
    global landingBurnAltitude
    assert direction == 1 or direction == -1, 'illegal input'

    if aAdjDir == -1 * direction: # if changing dir
      aAdjVal = aAdjVal * adjRatio # binary search effect

    if constant:
      effective_change = constant * direction
    else:
      effective_change = aAdjVal * direction


    print(f'{round(effective_change, 1)} alt because {reason}')
    landingBurnAltitude = landingBurnAltitude + effective_change
    aAdjDir = direction # keep track of cur direction



  # keep track of attempts
  altitudes = []

  # for attempt = 1:40
  for attempt in range(0, 40): #WARN
    # Reset I to the apex
    data.I = I_apex
    burnInitiated_I = 0 # Set to index of burn initialization
    burnStopped_I = 0 # Set when burn stops
    manualBurnShutdown = 0 # When burn is manually stopped due to sufficient slowdown
    outOfFuel = 0 # When burn expires due to lack of fuel
    altitudes.append(landingBurnAltitude)

    # Loop until landed
    while data.z[data.I-1] > 0:
      # Update atmosphere
      [Pa, rho, T, g] = atmosphere.ParametersAtAltitude(data.z[data.I-1])

      # Initiate burn?
      if burnInitiated_I == 0 and data.z[data.I-1] <= landingBurnAltitude:
        burnInitiated_I = data.I
      # end

      # Done with burn?
      if burnInitiated_I > 0 and burnStopped_I == 0:
        # Out of fuel?
        if data.m[data.I-1] <= stopMass:
          burnStopped_I = data.I
          outOfFuel = 1
        elif data.velocity(-1) < LandingVelocityAllowance/4:
          burnStopped_I = data.I
          manualBurnShutdown = 1
        # end
      # end

      # Update mass - keep burning while engines on
      if burnInitiated_I > 0 and burnStopped_I == 0:
        data.Thrust(mdot, vjet, Pe, Pa, Ae, N)
      else:
        data.NoThrust()
      # end

      # Angle of attack + Drag
      relWind = data.RelativeWind()
      # data.Fd[data.I] = 0.5 * Cd_down * A * rho * relWind.^2
      data.Fd[data.I] = 0.5 * Cd_down * A * rho * (relWind ** 2)

      # Axial -> X/Y
      a_axial = (data.Ft[data.I] + data.Fd[data.I]) / data.m[data.I]
      dv_axial = a_axial * data.dt
      data.UpdateVZ(dv_axial, g)

      # Update the angle
      data.UpdateAngle(dv_axial, 0, 0)

      # Net acceleration
      data.CalculateNewXandZ()

      # Next frame
      data.NextFrame()
    # end landed at this point


    # Rewind a frame
    data.I = data.I - 1

    # Determine if success or not
    landingVelocity = data.velocity(0)
    landedSafely = landingVelocity < LandingVelocityAllowance
    remainingFuel = data.m[data.I] - stopMass

    # If landed safely, then we found appropriate altitude to burn
    # engine. Stop the simulation here.
    if landedSafely:
      print("--- Landed safely after #{} attempts.\n".format(attempt))
      print("--- Landing burn altitude = {}m\n".format(landingBurnAltitude))
      return landingBurnAltitude, altitudes
    # end

    # Ok well did we stop burning early?
    if manualBurnShutdown > 0:
      # Let's just reduce the burn altitude by that much
      # landingBurnAltitude = landingBurnAltitude - data.z[burnStopped_I]
      shutdown_alt = data.z[burnStopped_I]
      adjust_landing_alt(-1, constant=shutdown_alt*SHUTDOWN_ALT_RATIO, reason=f'manual shutdown at {round(shutdown_alt)}m alt')

    # Did we run out of fuel?
    elif outOfFuel:
      LOW_ALT = 5
      # If we ran out of fuel at ~ LOW_ALT meters or less, and didn't reach nearly
      # safe velocities, then chances are we just didn't have enough fuel
      # regardless of when we burn. Fail
      if data.z[burnStopped_I] <= LOW_ALT:
        print("--- Ran out of fuel at low altitude.\n")
        return landingBurnAltitude, altitudes
      # end

      # Otherwise we should reduce the altitude
      shutdown_alt = data.z[burnStopped_I]
      adjust_landing_alt(-1, constant=shutdown_alt*SHUTDOWN_ALT_RATIO, reason="out of fuel at high alt")

    # Crashed while still burning!
    else:
      # Increase altitude based on amount of fuel remaining
      assert remainingFuel > 0
      adjust_landing_alt(+1, constant=remainingFuel/3, reason=f'crashed w/ {round(remainingFuel)}kg remaining fuel')
    # end
  # end
  print("--- Number of attempts has expired.\n")
  return landingBurnAltitude, altitudes
# end

## The Model
def simulate():
  global MAX_LAUNCH_RATIO, launch_ratio, launch_ratio_BEST, MaxDiscoveredVx
  global takeoff_angle, takeoff_angle_BEST, landingBurnAltitude
  global landingBurnAltitude_BEST, I_apex_BEST, I_land_BEST
  print("Starting search......\n")
  ratioAdjustment = 0.01
  ratioAdjustmentDir = -1
  MaxAttempts = args.max_attempts
  for attempt in range(0, MaxAttempts): # WARN 
    ratios.append(launch_ratio)
    # Print attempt
    print("##########\nStarting attempt #{}/{}\n".format(attempt, MaxAttempts))
    if launch_ratio > MAX_LAUNCH_RATIO:
      launch_ratio = MAX_LAUNCH_RATIO
    # end
    print("Launch fuel ratio = {}\n".format(launch_ratio))

    # Solve the takeoff angle
    print("1) Solving takeoff angle......\n")
    [takeoff_angle, success, (angles, errors)] = SolveTakeoffAngle(launch_ratio, data)
    I_apex = data.I
    if success == 0:
      print("--- Takeoff angle search failed.  Increasing launch ratio.\n")
      if ratioAdjustmentDir == -1:
        ratioAdjustment = ratioAdjustment/2
      # end
      launch_ratio = launch_ratio + ratioAdjustment
      ratioAdjustmentDir = 1
      continue
    # end
    print("--- Found takeoff angle = {}\n".format(takeoff_angle))

    # Detach second stage
    print("2) Detaching second stage.\n")
    data.m[I_apex] = data.m[I_apex] - mSecondStage

    # At this point we're at apex.  Time to solve best landing burn
    print("3) Solving landing burn......\n")
    landingBurnAltitude, altitudes = SolveLandingBurn(data, I_apex+1)
    I_land = data.I

    # Print some stats
    remainingMass = data.RemainingMass(mEmpty)
    landingSpeed = data.velocity(0)
    print("--- Landed with {}kg of fuel at {}m/s.\n".format(remainingMass, landingSpeed))

    # Is this optimal?
    if data.vx[I_apex] >= MaxDiscoveredVx and landingSpeed < LandingVelocityAllowance:
      print("*** This is a new optimal spec!\n")
      MaxDiscoveredVx = data.vx[I_apex]
      launch_ratio_BEST = launch_ratio
      takeoff_angle_BEST = takeoff_angle
      landingBurnAltitude_BEST = landingBurnAltitude
      I_apex_BEST = I_apex
      I_land_BEST = I_land
    # end

    # Increase or decrease fuel burn ratio
    if remainingMass > 0 and abs(landingSpeed) < LandingVelocityAllowance:
      print("Increasing launch fuel ratio.\n")
      if ratioAdjustmentDir == -1:
        ratioAdjustment = ratioAdjustment/2
      # end
      launch_ratio = launch_ratio + ratioAdjustment
      ratioAdjustmentDir = 1


    else:
      print("Reducing launch fuel ratio due to land failure.\n")
      if ratioAdjustmentDir == 1:
        ratioAdjustment = ratioAdjustment/2
      # end
      launch_ratio = launch_ratio - ratioAdjustment
      ratioAdjustmentDir = -1
    # end

    attempts_history.append((angles, altitudes))
# end

if not args.load:
  print('Simulation started at:', time.ctime())
  simulate()
  print('Simulation finished at:', time.ctime())
  # First need to trim the results
  t = data.t[0:I_land_BEST]
  z = data.z[0:I_land_BEST]
  vz = data.vz[0:I_land_BEST]
  az = data.az[0:I_land_BEST]
  x = data.x[0:I_land_BEST]
  vx = data.vx[0:I_land_BEST]
  Fd = data.Fd[0:I_land_BEST]
  print("##########\nAlgorithm complete!\n")
  print("Best launch ratio = {}\n".format(launch_ratio_BEST))
  print("Required takeoff angle = {}\n".format(takeoff_angle_BEST))
  print("Achieved horizontal velocity = {} m/s\n".format(MaxDiscoveredVx))
  print("Landing Burn Altitude = {}\n".format(landingBurnAltitude_BEST))
else:
  # load the data
  print('loadig simulation data..')
  (results, how, ratios, attempts_history) = load(args.load)

  dt = results['dt']

  t = how['t']
  z = how['z']
  vz = how['vz']
  az = how['az']
  x = how['x']
  vx = how['vx']
  Fd = how['Fd']

  print('loaded saved simulation data')
  print(results)

if (not args.load and args.save):
  print('saving simulation data..')
  how = {
    't': data.t[0:I_land_BEST],
    'z': data.z[0:I_land_BEST],
    'vz': data.vz[0:I_land_BEST],
    'az': data.az[0:I_land_BEST],
    'x': data.x[0:I_land_BEST],
    'vx': data.vx[0:I_land_BEST],
    'Fd': data.Fd[0:I_land_BEST]
  }
  results = {
    'dt': dt,
    'launch_ratio': launch_ratio_BEST,
    'i_land': I_land_BEST,
    'takeoff_angle': takeoff_angle_BEST,
    'max_discovered_vx': MaxDiscoveredVx,
    'landing_burn_altitude': landingBurnAltitude_BEST,
  }
  data = (results, how, ratios, attempts_history)
  dump(f'{args.save}', data)


########## Plots ###############

plt.plot(t, z)
plt.xlabel('Time')
plt.ylabel("Altitude")
plt.title("Altitude vs Time")
plt.savefig('figs/alt-time.jpg')
plt.close()

plt.plot(t, vx)
plt.xlabel('Time')
plt.ylabel("Horizontal Position")
plt.title("X vs Time")
plt.savefig('figs/x-time.jpg')
plt.close()

plt.plot(t, Fd)
plt.xlabel('Time')
plt.ylabel("Drag")
plt.title("Fd vs Time")
plt.savefig('figs/fd-time.jpg')
plt.close()

plt.plot(t, az)
plt.xlabel('Time')
plt.ylabel("Acceleration")
plt.title("Az vs Time")
plt.savefig('figs/vertical-acceleration-time.jpg')
plt.close()

plt.plot(t, z, label="z")
plt.plot(t, x, label="x")
plt.xlabel('Time')
plt.ylabel("Position")
plt.title("Position vs Time")
plt.legend()
plt.savefig('figs/position-time.jpg')
plt.close()

plt.plot(x, z)
plt.xlabel('Horizontal Pos')
plt.ylabel("Vertical Pos")
plt.title("X vs Z")
plt.legend()
plt.savefig('figs/x-z.jpg')
plt.close()

plt.plot(t, vz, label="vz")
plt.plot(t, vx, label="vx")
plt.xlabel('Time')
plt.ylabel("Speed")
plt.title("Speed vs Time")
plt.legend()
plt.savefig('figs/speed-time.jpg')
plt.close()

# launch ratios vs takeoff angles
best_angles = []
best_altitudes = []
for attempt, (angles, altitudes) in enumerate(attempts_history):
  best_angle = angles[-1] *10000
  best_alt = altitudes[-1]
  best_angles.append(best_angle)
  best_altitudes.append(best_alt)

r1 = ratios[:]
r1, best_altitude = sort_two_lists(r1, best_altitudes)
plt.plot(r1, best_altitudes)
plt.scatter(r1, best_altitudes)
plt.xlabel('Launch ratio')
plt.ylabel("Burn altitude (m)")
plt.title("Launch Ratio vs Landing Burn ltitudes")
plt.savefig('figs/lr-altitude.jpg')
plt.close()

r2 = ratios[:]
r2, best_angles = sort_two_lists(r2, best_angles)
plt.plot(r2, best_angles)
plt.scatter(r2, best_angles)
plt.xlabel('Launch ratio')
plt.ylabel("Takeoff angle (*1e+4)")
plt.title("Launch Ratio vs Takeoff Angle ")
plt.savefig('figs/lr-angle.jpg')
plt.close()

plot_attempts(ratios, f'dt{dt}-ratios.jpg',
             ylabel='landing burn altitude', title='Explored ratios')
plt.close()


# plot individual attempts
for attempt, (angles, altitudes) in enumerate(attempts_history):
  # plot attempts
  angles = np.array(angles) * 10000
  plot_attempts(angles, fname=f'dt{dt}-a{attempt}-angles.jpg',
               ylabel='takeoff angle (*1e+4)',
               title=f'Explored angles at LR {ratios[attempt]} - Attempt #{attempt+1}',
               line_label='Angles')
  altitudes = np.array(altitudes) / 1000
  plot_attempts(altitudes, fname=f'dt{dt}-a{attempt}-altitudes.jpg',
               ylabel='landing burn altitude (1km)',
               title=f'Explored altitudes at LR {ratios[attempt]} - Attempt #{attempt+1}',
               line_label="Altitutdes")


angles_arr = []
altitudes_arr = []
for attempt, (angles, altitudes) in enumerate(attempts_history):
  angles_arr.append(np.array(angles) * 1e+4)
  altitudes_arr.append(np.array(altitudes) / 1000)

plot_batch_attempts(altitudes_arr, ratios, fname=f'dt{dt}-altitudes.jpg',
               ylabel='Landing burn altitude (km)',
               title=f'Landing Burn Altitude State Space Exploration')

plot_batch_attempts(angles_arr, ratios, fname=f'dt{dt}-angles.jpg',
               ylabel='Takeoff angle (*1e+4)',
               title=f'Takeoff Angle State Space Exploration')
