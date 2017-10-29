from Servo import Servo #import servo class

robot_servo = Servo()  #create instance of Servo

#set first 8 servos to position 60 degree and set next 8 servos to -60 degrees
robot_servo.set_servo (0, 0)
robot_servo.set_servo (1, 0)
robot_servo.set_servo (2, 0)
robot_servo.set_servo (3, 0)
robot_servo.set_servo (4, 0)
robot_servo.set_servo (5, 0)
robot_servo.set_servo (6, 0)
robot_servo.set_servo (7, 0)
robot_servo.set_servo (8, 0)
robot_servo.set_servo (9, 0)
robot_servo.set_servo (10, 0)
robot_servo.set_servo (11, 0)
robot_servo.set_servo (12, 0)
robot_servo.set_servo (13, 0)
robot_servo.set_servo (14, 0)
robot_servo.set_servo (15, 0)

# do beaware I haven't created a shutdown toutine so even after the program termintes 
# the PCA9685 will continue to hold these positions until power down.
