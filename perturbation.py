import orekit
vm = orekit.initVM()

from orekit.pyhelpers import setup_orekit_curdir
setup_orekit_curdir()

from org.orekit.propagation.analytical.tle import TLE, TLEPropagator
from org.orekit.propagation.analytical import KeplerianPropagator
from org.orekit.time import TimeScalesFactory
from org.orekit.frames import FramesFactory
import numpy as np
import matplotlib.pyplot as plt
import time

#International Space Station's TLE
ISS_tle1 = "1 25544U 98067A   26165.80282117  .00008449  00000+0  16005-3 0  9998"
ISS_tle2 = "2 25544  51.6336 311.8940 0004856 186.3522 173.7405 15.49231655571384"

#Hubble Telescope's TLE
hub_tle1 = "1 20580U 90037B   26165.14417512  .00006344  00000+0  19999-3 0  9995"
hub_tle2 = "2 20580  28.4715 101.3706 0001912 110.8714 249.2086 15.30716351788104"

#geo satellite's TLE
goes19_tle1 = "1 60133U 24119A   26165.81642218 -.00000243  00000+0  00000+0 0  9990"
goes19_tle2 = "2 60133   0.0120  97.6670 0000378  89.6211 294.5445  1.00272177  6947"

ISS_tle = TLE(ISS_tle1,ISS_tle2)
hub_tle = TLE(hub_tle1,hub_tle2)
goes19_tle = TLE(goes19_tle1,goes19_tle2)

#TLE_Propagation
ISS_propagator = TLEPropagator.selectExtrapolator(ISS_tle)
hub_propagator = TLEPropagator.selectExtrapolator(hub_tle)
goes19_propagator = TLEPropagator.selectExtrapolator(goes19_tle)

#Keplerain Orbit Propagation
ISS_kep_orbit = ISS_propagator.propagate(ISS_tle.getDate()).getOrbit()
ISS_kep_prop = KeplerianPropagator(ISS_kep_orbit)

hub_kep_orbit = hub_propagator.propagate(hub_tle.getDate()).getOrbit()
hub_kep_prop = KeplerianPropagator(hub_kep_orbit)

goes19_kep_orbit = goes19_propagator.propagate(goes19_tle.getDate()).getOrbit()
goes19_kep_prop = KeplerianPropagator(goes19_kep_orbit)

#Calculation Part
def orbitcalc (tle_prop , kep_prop , startdate , duration_days = 7.0 , step_secs = 600.0): #Step size is 10 mins
    
    times = []
    dist = []
    enddate = startdate.shiftedBy(duration_days * 24.0 * 60.0 * 60.0)
    currentdate = startdate
    while (currentdate.compareTo(enddate) <= 0.0):
        #prop TLE
        state_tle = tle_prop.propagate(currentdate)
        #get pos
        pos_tle = state_tle.getPVCoordinates().getPosition()
        
        #prop kep
        state_kep = kep_prop.propagate(currentdate)
        #get pos
        pos_kep = state_kep.getPVCoordinates().getPosition()

        #calc and store the distance between TLE and KEP
        distance = pos_kep.distance(pos_tle) / 1000.0
        dist.append(distance)

        #Elasped days

        elasped_days = currentdate.durationFrom(startdate) / 86400.0
        times.append(elasped_days) 

        currentdate = currentdate.shiftedBy(step_secs)

    return times,dist

#calling the func
print("Simulating ISS...")
time.sleep(1)
iss_times, iss_dist = orbitcalc(ISS_propagator, ISS_kep_prop,ISS_tle.getDate())

print("Simulating Hubble...")
time.sleep(1)
hub_times, hub_dist = orbitcalc(hub_propagator, hub_kep_prop,hub_tle.getDate())

print("Simulating GOES-19...")
time.sleep(1)
goes19_times, goes19_dist = orbitcalc(goes19_propagator, goes19_kep_prop,goes19_tle.getDate())

#Plotting the results
plt.figure(figsize=(12,7))
plt.plot(iss_times,iss_dist, label = "International Space Station (LEO)", color='red', linewidth =2)
plt.plot(hub_times,hub_dist, label = "Hubble (LEO)", color='blue', linewidth =2)
plt.plot(goes19_times,goes19_dist, label = "GOES-19 (GEO)", color='green', linewidth =2)

plt.yscale
plt.xlabel("Elasped Time in Days", fontsize = 10)
plt.ylabel("Distance b/w the Satellites", fontsize = 10)
plt.title("Perturbation Comparison of Different Orbits at Different Altitudes", fontsize = 10)
plt.grid(True,which="both",linestyle= "--",alpha = 0.5)
plt.legend(fontsize=11)
plt.tight_layout()
plt.show()