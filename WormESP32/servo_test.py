from machine import PWM, I2C, Pin
from math import pi, radians
from BOTS import *  #import all custom code for the board

#create an instance of the I2C bus
I2C_bus = I2C(0, sda=27, scl=26)

#create instance of the servo driver
robot = Servo(I2C_bus)

#set servos to test positions.
robot.set_servo (radians(45), 0)
robot.set_servo (radians(-45), 1)
robot.set_servo (radians(45), 2)
robot.set_servo (radians(-45), 3)
robot.set_servo (radians(45), 4)
robot.set_servo (radians(-45), 5)
robot.set_servo (radians(45), 6)
robot.set_servo (radians(-45), 7)

I2C_bus.deinit()
