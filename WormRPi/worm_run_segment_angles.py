# Run with segment angles file.
# Options: [<input file name (defaults to segment_angles.txt)>]

import sys
#from Servo import Servo
from time import sleep
#import sonar

# Parameters
# ping_angle should be set such that left and right pings should provide
# food object sensing coverage extensively yet exclusively on their respective sides.
angle_amplifier = 1.25
servo_delay = .01
ping_angle = 40
ping_stride = 2
ping_before_pause = 1
ping_after_pause = 1
angle_bias = -2.0
left_turn_angle_delta = -2.0
right_turn_angle_delta = 4.0
arrival_distance = 10.0
head_swing_step_reset = -1

# Read input file.
angles_array = []
filename = "segment_angles.txt"
if len(sys.argv) == 2:
    filename = sys.argv[2]
input_file = open(filename, 'r')
lines = input_file.readlines()
input_file.close()
for line in lines:
    vals = (line[:-1]).split(" ")
    angles = []
    for angle in vals:
        angles.append(int(angle))
    angles_array.append(angles);

# Set up servos.
#robot = Servo()

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
    print("step=", step)
    for segment in range(8):
        position = (int)(((float)(angles_array[step][segment]) * angle_amplifier) + angle_bias + turn_angle_delta)
        #robot.set_servo ((11 - segment), position)
        #print("servo=", (11 - segment), "position=", position)
    sleep(servo_delay)

    # Check food distance.
    if step > 0 and step < (angles_array_len - 1):
        prev_angle = (int)(((float)(angles_array[step - 1][0]) * angle_amplifier) + angle_bias + turn_angle_delta)
        curr_angle = (int)(((float)(angles_array[step][0]) * angle_amplifier) + angle_bias + turn_angle_delta)
        next_angle = (int)(((float)(angles_array[step + 1][0]) * angle_amplifier) + angle_bias + turn_angle_delta)
        if curr_angle > 0 and next_angle <= 0:
            head_swing_count += 1
            if head_swing_step_reset != -1 and head_swing_count >= head_swing_step_reset:
                step = -1
                head_swing_count = 0
                robot = Servo()
                print("step reset")
        change = False
        if right_check == False and prev_angle >= 0 and curr_angle > 0 and next_angle > 0:
            curr_angle_delta = abs(curr_angle - ping_angle)
            if abs(prev_angle - ping_angle) >= curr_angle_delta and abs(next_angle - ping_angle) >= curr_angle_delta:
               right_check = True
               left_check = False
               if right_ping_count > 0 and (right_ping_count % ping_stride) == 0:
                   if right_ping_toggle == True:
                       right_ping_toggle = False
                       sleep(ping_before_pause)
                       #dummy = int(input("enter 0 to read right distance..."))
                       #right_food_distance = sonar.ping()
                       print("right_food_distance=", right_food_distance)
                       sleep(ping_after_pause)
                       #change = True
                       #turn_angle_delta = int(input("right, turn_angle_delta="))
                   else:
                       right_ping_toggle = True
               right_ping_count = right_ping_count + 1
        if left_check == False and prev_angle <= 0 and curr_angle < 0 and next_angle < 0:
            curr_angle_delta = abs(-curr_angle - ping_angle)
            if abs(-prev_angle - ping_angle) >= curr_angle_delta and abs(-next_angle - ping_angle) >= curr_angle_delta:
               left_check = True
               right_check = False
               if left_ping_count > 0 and (left_ping_count % ping_stride) == 0:
                   if left_ping_toggle == True:
                       left_ping_toggle = False
                       sleep(ping_before_pause)
                       #dummy = int(input("enter 0 to read left distance..."))
                       #left_food_distance = sonar.ping()
                       print("left_food_distance=", left_food_distance)
                       sleep(ping_after_pause)
                       #change = True
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
