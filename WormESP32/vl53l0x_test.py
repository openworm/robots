import time, display, network
from machine import PWM, I2C, Pin
from math import pi, radians
from microWebSrv import MicroWebSrv
from BOTS import *  #import all custom code for the board
import vl53l0x

#create an instance of the I2C bus
I2C_sensor_bus = I2C(1, speed=400000, sda=32, scl=33)

#create an instance of the VL53L0X distance sensor.
distance_sensor = vl53l0x.VL53L0X(I2C_sensor_bus)
  
while True: 
    print(distance_sensor.read())
    time.sleep(3)



