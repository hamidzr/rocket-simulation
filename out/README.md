# Simulation data output
Since running the simulation is costly here we provide the data gathered by the simulaiton so that it can be easily loaded and inspected.

The output data inclues: 

```
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
# ratios: fuel ratios tested.
data = (results, how, ratios, attempts_history)

```
