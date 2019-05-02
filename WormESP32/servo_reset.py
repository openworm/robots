from machine import PWM, I2C, Pin
from math import pi, radians
from BOTS import *  #import all custom code for the board

#create an instance of the I2C bus
I2C_servo_bus = I2C(0, speed=100000, sda=27, scl=26)

#create instance of the servo driver
robot = Servo(I2C_servo_bus)

#set the servos to position 0
robot.set_servo (0, 0)
robot.set_servo (0, 1)
robot.set_servo (0, 2)
robot.set_servo (0, 3)
robot.set_servo (0, 4)
robot.set_servo (0, 5)
robot.set_servo (0, 6)
robot.set_servo (0, 7)

I2C_servo_bus.deinit()
