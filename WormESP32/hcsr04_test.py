import time
import hcsr04

#create an instance of the HCSR04 distance sensor.
# trigger=32, echo=33
distance_sensor = hcsr04.HCSR04(32, 33)
  
while True: 
    print(distance_sensor.distance_cm())
    time.sleep(3)
