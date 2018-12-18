# Run with segment angles file.
# Options: [<input file name (defaults to segment_angles.txt)>]
#from machine import PWM, I2C, Pin
from math import pi, radians
#from BOTS import *  #import all custom code for the board
import sys
from time import sleep
#import sonar

# Parameters
# ping_angle should be set such that left and right pings should provide
# food object sensing coverage extensively yet exclusively on their respective sides.
#angle_amplifier = 1.25
angle_amplifier = 1.0
servo_delay = .05
ping_angle = 40
ping_stride = 2
ping_before_pause = 1
ping_after_pause = 1
angle_bias = 1.75
left_turn_angle_delta = 2.0
right_turn_angle_delta = -2.0
arrival_distance = 10.0
head_swing_step_restart = 5

# Read input file.
angles_array = []
filename = "segment_angles.txt"
if len(sys.argv) == 2:
    filename = sys.argv[2]
with open(filename, "r") as f:
    for line in f:
        vals = (line[:-1]).split(" ")
        angles = []
        for angle in vals:
            angles.append(int(angle))
        angles_array.append(angles);

#create an instance of the I2C bus
#I2C_bus = I2C(0, sda=27, scl=26)

#create instance of the servo driver
#robot = Servo(I2C_bus)

# Run.
left_food_distance = -1.0
right_food_distance = -1.0
left_check = False
right_check = False
head_swing_count = 0
right_ping_count = 0
left_ping_count = 0
right_ping_toggle = True
left_ping_toggle = False
turn_angle_delta = 0.0
angles_array_len = len(angles_array)
step = 0
while step < angles_array_len:
    #print("step=", step)
    for segment in range(8):
        position = (int)(((float)(angles_array[step][segment]) * angle_amplifier) + angle_bias + turn_angle_delta)
        #robot.set_servo (radians(position), 7 - segment)
        #print("servo=", segment, "position=", position)
    sleep(servo_delay)

    if step > 0 and step < (angles_array_len - 1):
        prev_angle = (int)(((float)(angles_array[step - 1][0]) * angle_amplifier) + angle_bias + turn_angle_delta)
        curr_angle = (int)(((float)(angles_array[step][0]) * angle_amplifier) + angle_bias + turn_angle_delta)
        next_angle = (int)(((float)(angles_array[step + 1][0]) * angle_amplifier) + angle_bias + turn_angle_delta)
        
        # Restart steps?
        if curr_angle > 0 and next_angle <= 0:
            head_swing_count += 1
            if head_swing_step_restart != -1 and head_swing_count >= head_swing_step_restart:
                step = -1
                head_swing_count = 0
                print("step restart")

        # Check food distance.
        change = False
        if right_check == False and prev_angle <= 0 and curr_angle < 0 and next_angle < 0:
            curr_angle_delta = abs(-curr_angle - ping_angle)
            if abs(-prev_angle - ping_angle) >= curr_angle_delta and abs(-next_angle - ping_angle) >= curr_angle_delta:
               right_check = True
               left_check = False
               if right_ping_count > 0 and (right_ping_count % ping_stride) == 0:
                   if right_ping_toggle == True:
                       right_ping_toggle = False
                       sleep(ping_before_pause)
                       right_food_distance = float(input("enter right distance..."))
                       #right_food_distance = float(sonar.ping())
                       print("right_food_distance=", right_food_distance)
                       sleep(ping_after_pause)
                       change = True
                       #turn_angle_delta = int(input("right, turn_angle_delta="))
                   else:
                       right_ping_toggle = True
               right_ping_count = right_ping_count + 1
        if left_check == False and prev_angle >= 0 and curr_angle > 0 and next_angle > 0:
            curr_angle_delta = abs(curr_angle - ping_angle)
            if abs(prev_angle - ping_angle) >= curr_angle_delta and abs(next_angle - ping_angle) >= curr_angle_delta:
               left_check = True
               right_check = False
               if left_ping_count > 0 and (left_ping_count % ping_stride) == 0:
                   if left_ping_toggle == True:
                       left_ping_toggle = False
                       sleep(ping_before_pause)
                       left_food_distance = float(input("enter left distance..."))
                       #left_food_distance = float(sonar.ping())
                       print("left_food_distance=", left_food_distance)
                       sleep(ping_after_pause)
                       change = True
                       #turn_angle_delta = int(input("left, turn_angle_delta="))
                   else:
                       left_ping_toggle = True
               left_ping_count = left_ping_count + 1
               
        # Determine direction.
        if change == True and left_food_distance != -1.0 and right_food_distance != -1.0:
            if left_food_distance <= arrival_distance or right_food_distance <= arrival_distance:
                print("arrived")
                break
            elif left_food_distance < right_food_distance:
                print("food to left")
                turn_angle_delta = left_turn_angle_delta
            elif right_food_distance < left_food_distance:
                print("food to right")
                turn_angle_delta = right_turn_angle_delta
            else:
                turn_angle_delta = 0.0
    step += 1
    
#I2C_bus.deinit()
