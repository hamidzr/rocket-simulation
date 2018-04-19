from utils.atmosphere import Atmosphere
from utils.simulation_data import SimulationData
import math


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
A = math.pow(3.7, 2) / 4 # Diameter of 3.7m
Thrust_vacuum = 8227000 # N
Thrust_ground = 7607000 # N
Isp_vacuum = 311 # s
Isp_ground = 282 # s
vjet = 3000
TargetAltitude = 200000 # m
LandingVelocityAllowance = 6 # m/s
[Pe, Ae, mdot] = SolveThrust(vjet, Thrust_ground, Thrust_vacuum, Isp_ground, Isp_vacuum, N)


# Parameters being solved
MAX_LAUNCH_RATIO = 0.99999999 # should never happen, but prevents >= 1.
launch_ratio = 0.94
launch_ratio_BEST = launch_ratio
MaxDiscoveredVx = 0
takeoff_angle = 0.00001
takeoff_angle_BEST = takeoff_angle
landingBurnAltitude = 6000 # m
landingBurnAltitude_BEST = landingBurnAltitude
I_apex_BEST = 1
I_land_BEST = 1

# Create simulation data
data = SimulationData(0.001, 2000)


## The Model
print("Starting search......\n")
ratioAdjustment = 0.01
ratioAdjustmentDir = -1
MaxAttempts = 10
for attempt in range(0, MaxAttempts): # WARN 
  # Print attempt
  print("##########\nStarting attempt #%d/%d\n", attempt, MaxAttempts)
  if launch_ratio > MAX_LAUNCH_RATIO:
    launch_ratio = MAX_LAUNCH_RATIO
  # end
  print("Launch fuel ratio = %f\n", launch_ratio)

  # Solve the takeoff angle
  print("1) Solving takeoff angle......\n")
  [takeoff_angle, success] = SolveTakeoffAngle(launch_ratio, data)
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
  print("--- Found takeoff angle = %f\n", takeoff_angle)

  # Detach second stage
  print("2) Detaching second stage.\n")
  data.m(I_apex) = data.m(I_apex) - mSecondStage

  # At this point we're at apex.  Time to solve best landing burn
  print("3) Solving landing burn......\n")
  landingBurnAltitude = SolveLandingBurn(data, I_apex+1)
  I_land = data.I

  # Print some stats
  remainingMass = data.RemainingMass(mEmpty)
  landingSpeed = data.velocity(0)
  print("--- Landed with %fkg of fuel at %fm/s.\n", remainingMass, landingSpeed)

  # Is this optimal?
  if data.vx(I_apex) >= MaxDiscoveredVx and landingSpeed < LandingVelocityAllowance:
    print("*** This is a new optimal spec!\n")
    MaxDiscoveredVx = data.vx(I_apex)
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
# end

print("##########\nAlgorithm complete!\n")
print("Best launch ratio = %f\n", launch_ratio_BEST)
print("Required takeoff angle = %f\n", takeoff_angle_BEST)
print("Achieved horizontal velocity = %f m/s\n", MaxDiscoveredVx)
print("Landing Burn Altitude = %f\n", landingBurnAltitude_BEST)


## Plots
# First need to trim the results
t = data.t(0:I_land_BEST)
z = data.z(0:I_land_BEST)
vz = data.vz(0:I_land_BEST)
az = data.az(0:I_land_BEST)
x = data.x(0:I_land_BEST)
vx = data.vx(0:I_land_BEST)
Fd = data.Fd(0:I_land_BEST)

# TODO plots
# figure(0)
# plot(t, vz)
# xlabel('Time')
# ylabel("Altitude")
# title("Altitude vs Time")

# figure(2)
# plot(t, vx)
# xlabel('Time')
# ylabel("Horizontal Position")
# title("X vs Time")

# figure(3)
# plot(t, Fd)
# xlabel('Time')
# ylabel("Drag")
# title("Fd vs Time")




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


## Solve Takeoff Angle
# Responsible for solving the takeoff angle given a launch ratio.
# Will start at 0 to ensure the desired altitude can even be reached.
def SolveTakeoffAngle(launch_ratio, data):
  # Bring in globals
  global atmosphere TargetAltitude mSecondStage A
  global m0_propellant mEmpty m0_full Cd_up Pe Ae mdot vjet

  # Determine when to stop
  N = 9 # Use all engines
  stopMass = mEmpty + mSecondStage + m0_propellant * (1 - launch_ratio)
  data.SetInitialConditions(m0_full)

  # Iterate the times until there
  takeoff_angle = 0.000333 # Start attempt angle
  adjustment =  0.000001 # Start adjustment
  adjustmentDir = 1
  attempt = 1
  errorAllowance = 0.001
  while True: # Keep looping until terminated
    # Initial conditions
    data.theta(0) = 0

    # Loop until at peak altitude
    data.ResetFrame()
    while data.vz(data.I-1) >= 0:
      # Update atmosphere
      [Pa, rho, T, g] = atmosphere.ParametersAtAltitude(data.z(data.I-1))

      # Update mass
      if data.m(data.I-1) > stopMass:
        data.Thrust(mdot, vjet, Pe, Pa, Ae, N)
      else:
        data.NoThrust()
      # end

      # Angle of attack + Drag
      relWind = data.RelativeWind()
      data.Fd(data.I) = 0.5 * Cd_up * A * rho * relWind.^2

      # Axial -> X/Y
      a_axial = (data.Ft(data.I) - data.Fd(data.I)) / data.m(data.I)
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
    error = abs(data.z(data.I) - TargetAltitude)
    if error <= TargetAltitude * errorAllowance:
      success = 1
      print("--- Error = %fm after %d attempts.\n", error, attempt)
      return
    # end

    # If didn't reach altitude...
    if data.z(data.I) < TargetAltitude:
      # If takeoff angle is 0, not enough fuel!
      if takeoff_angle == 0:
        success = 0
        return
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
  return [takeoff_angle, success]
# end




## Solve Landing Burn
# Responsible for solving the takeoff angle given a launch ratio.
# Will start at 0 to ensure the desired altitude can even be reached.
function landingBurnAltitude = SolveLandingBurn(data, I_apex)
  # Bring in globals
  global atmosphere LandingVelocityAllowance A
  global mEmpty Cd_down Pe Ae mdot vjet

  # Determine when to stop
  N = 3 # Just use three engines
  stopMass = mEmpty # Stop when no more propellant

  # Iterate the times until there
  landingBurnAltitude = 27000 # Pick an abitrary start
  for attempt = 1:40
    # Reset I to the apex
    data.I = I_apex
    burnInitiated_I = 0 # Set to index of burn initialization
    burnStopped_I = 0 # Set when burn stops
    manualBurnShutdown = 0 # When burn is manually stopped due to sufficient slowdown
    outOfFuel = 0 # When burn expires due to lack of fuel

    # Loop until landed
    while data.z(data.I-1) > 0:
      # Update atmosphere
      [Pa, rho, T, g] = atmosphere.ParametersAtAltitude(data.z(data.I-1))

      # Initiate burn?
      if burnInitiated_I == 0 and data.z(data.I-1) <= landingBurnAltitude:
        burnInitiated_I = data.I
      # end

      # Done with burn?
      if burnInitiated_I > 0 and burnStopped_I == 0:
        # Out of fuel?
        if data.m(data.I-1) <= stopMass:
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
      data.Fd(data.I) = 0.5 * Cd_down * A * rho * relWind.^2

      # Axial -> X/Y
      a_axial = (data.Ft(data.I) + data.Fd(data.I)) / data.m(data.I)
      dv_axial = a_axial * data.dt
      data.UpdateVZ(dv_axial, g)

      # Update the angle
      data.UpdateAngle(dv_axial, 0, 0)

      # Net acceleration
      data.CalculateNewXandZ()

      # Next frame
      data.NextFrame()
    # end

    # Rewind a frame
    data.I = data.I - 1

    # Determine if success or not
    landingVelocity = data.velocity(0)
    landedSafely = landingVelocity < LandingVelocityAllowance
    remainingFuel = data.m(data.I) - stopMass

    # If landed safely, then we found appropriate altitude to burn
    # engine. Stop the simulation here.
    if landedSafely:
      print("--- Landed safely after #%d attempts.\n", attempt)
      print("--- Landing burn altitude = %dm\n", landingBurnAltitude)
      return
    # end

    # Ok well did we stop burning early?
    if manualBurnShutdown > 0:
      # Let's just reduce the burn altitude by that much
      landingBurnAltitude = landingBurnAltitude - data.z(burnStopped_I)

    # Did we run out of fuel?
    elif outOfFuel:
      # If we ran out of fuel at ~500m or less, and didn't reach nearly
      # safe velocities, then chances are we just didn't have enough fuel
      # regardless of when we burn. Fail
      if data.z(burnStopped_I) <= 500:
        print("--- Ran out of fuel at low altitude.\n")
        return
      # end

      # Otherwise we should reduce the altitude by bulk.
      landingBurnAltitude = landingBurnAltitude - 500

    # Crashed while still burning!
    else:
      # Increase altitude based on amount of fuel remaining
      landingBurnAltitude = landingBurnAltitude + remainingFuel/3
    # end
  # end
  print("--- Number of attempts has expired.\n")
# end
