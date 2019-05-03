# Run with segment angles file "segment_angles.txt".
# An ultrasonic distance sensor simulates food foraging by locating and moving toward the nearest object.
# Maximum steps can optionally be set by creating file
# "max_steps.txt" containing value.

from machine import PWM, I2C, Pin
from math import pi, radians
from BOTS import *  #import all custom code for the board
import sys
from time import sleep
import hcsr04

# Parameters.
angle_amplifier = 1.0  # Increases wave amplitude
servo_delay = .05
right_ping_angles = [ -45.0 ]  # Angles to ping objects
left_ping_angles = [ 45.0 ]
ping_stride = 1
ping_before_pause = .5
ping_after_pause = 0
max_ping_tries = 5
angle_bias = 1.75
left_turn_angle_delta = 3.327  # Increase to make sharper turn
right_turn_angle_delta = -1.5
forward_angle_delta = .327
head_swing_step_index_restart = 10
goal_distance = 25.0
max_valid_food_distance = 200

# Load segment angles file.
segment_angles_file = "segment_angles.txt"
print("loading segment angles...", end='')
angles_array = []
with open(segment_angles_file, "r") as f:
    for line in f:
        vals = (line[:-1]).split(" ")
        angles = []
        for angle in vals:
            angles.append(int(angle))
        angles_array.append(angles);
print("done")

# Load optional maximum steps from file.
max_steps = -1
max_steps_file = "max_steps.txt"
try:
    with open(max_steps_file, "r") as f:
        max_steps = int(next(f))
        print("max_steps=", max_steps)
except:
    pass

# Straighten robot.
exec(open("servo_reset.py").read())
sleep(2)

# Initialize servos.
I2C_servo_bus = I2C(0, speed=100000, sda=27, scl=26)
robot = Servo(I2C_servo_bus)

# Initialize sonic  distance sensor.
# trigger=32, echo=33
distance_sensor = hcsr04.HCSR04(32, 33)

# Run.
right_ping_angles_len = len(right_ping_angles)
right_ping_distances = [ max_valid_food_distance + 1 ] * right_ping_angles_len
left_ping_angles_len = len(left_ping_angles)
left_ping_distances = [ max_valid_food_distance + 1 ] * left_ping_angles_len
head_swing_side = 'neutral'
head_swing_count = 0
right_swing_count = 0
left_swing_count = 0
turn_angle_delta = forward_angle_delta
angles_array_len = len(angles_array)
step = 0
step_index = 0
while step_index < angles_array_len and (max_steps == -1 or step < max_steps):
    #print("step=", step)
    for segment in range(8):
        position = (int)(((float)(angles_array[step_index][segment]) * angle_amplifier) + angle_bias + turn_angle_delta)
        try:
            robot.set_servo (radians(position), 7 - segment)              
            #print("servo=", segment, "position=", position)
        except:
            print("error setting servo angle")
    sleep(servo_delay)

    if step_index > 0 and step_index < (angles_array_len - 1):
        prev_angle = (int)(((float)(angles_array[step_index - 1][0]) * angle_amplifier) + angle_bias + turn_angle_delta)
        curr_angle = (int)(((float)(angles_array[step_index][0]) * angle_amplifier) + angle_bias + turn_angle_delta)
        next_angle = (int)(((float)(angles_array[step_index + 1][0]) * angle_amplifier) + angle_bias + turn_angle_delta)
        
        # Restart step index?
        if curr_angle > 0 and next_angle <= 0:
            head_swing_count += 1
            if head_swing_step_index_restart != -1 and head_swing_count >= head_swing_step_index_restart:
                step_index = -1
                head_swing_count = 0
                head_swing_side = 'neutral'
                right_swing_count = left_swing_count = 0
                #print("step restart")

        # Check right food distance.
        if prev_angle <= 0 and curr_angle < prev_angle:
            if head_swing_side == 'neutral' or head_swing_side == 'left':
                if head_swing_side == 'left':
                    left_swing_count = left_swing_count + 1
                    #print("left food distance=", min(left_ping_distances))
                head_swing_side = 'right'
            if (right_swing_count % ping_stride) == 0:
                for i in range(right_ping_angles_len):
                    ping_angle = right_ping_angles[i]
                    curr_angle_delta = abs(curr_angle - ping_angle)
                    if abs(prev_angle - ping_angle) >= curr_angle_delta and abs(next_angle - ping_angle) >= curr_angle_delta:
                        prev_distance = right_ping_distances[i]
                        for ping in range(max_ping_tries):
                            sleep(ping_before_pause)
                            try:
                                right_ping_distances[i] = distance_sensor.distance_cm()
                            except:
                                print("error sensing distance")
                            sleep(ping_after_pause)
                            if right_ping_distances[i] < 0 or right_ping_distances[i] > max_valid_food_distance:
                                #print("out of bounds right=", right_ping_distances[i], " ping=", ping)
                                distance_sensor = hcsr04.HCSR04(32, 33)
                                if right_ping_distances[i] < 0:
                                    right_ping_distances[i] = prev_distance
                                else:
                                    right_ping_distances[i] = max_valid_food_distance + 1
                            else:
                                break
                        print("right ping: angle=", curr_angle, "distance=", right_ping_distances[i])

        # Check left food distance.
        if prev_angle >= 0 and curr_angle > prev_angle:
            if head_swing_side == 'neutral' or head_swing_side == 'right':
                if head_swing_side == 'right':
                    right_swing_count = right_swing_count + 1
                    #print("right food distance=", min(right_ping_distances))
                head_swing_side = 'left'
            if (left_swing_count % ping_stride) == 0:
                for i in range(left_ping_angles_len):
                    ping_angle = left_ping_angles[i]
                    curr_angle_delta = abs(curr_angle - ping_angle)
                    if abs(prev_angle - ping_angle) >= curr_angle_delta and abs(next_angle - ping_angle) >= curr_angle_delta:
                        prev_distance = left_ping_distances[i]
                        for ping in range(max_ping_tries):
                            sleep(ping_before_pause)
                            try:
                                left_ping_distances[i] = distance_sensor.distance_cm()
                            except:
                                print("error sensing distance")
                            sleep(ping_after_pause)
                            if left_ping_distances[i] < 0 or left_ping_distances[i] > max_valid_food_distance:
                                #print("out of bounds left=", left_ping_distances[i], " ping=", ping)
                                distance_sensor = hcsr04.HCSR04(32, 33)
                                if left_ping_distances[i] < 0:
                                    left_ping_distances[i] = prev_distance
                                else:
                                    left_ping_distances[i] = max_valid_food_distance + 1
                            else:
                                break
                        print("left ping: angle=", curr_angle, "distance=", left_ping_distances[i])
               
        # Determine direction.
        min_right_ping_distance = min(right_ping_distances)
        #print("right_ping_distance=", min_right_ping_distance) ###
        min_left_ping_distance = min(left_ping_distances)
        #print("left_ping_distance=", min_left_ping_distance) ###
        prev_turn_angle_delta = turn_angle_delta
        if min_right_ping_distance <= goal_distance or min_left_ping_distance <= goal_distance:
            print("food found")
            break
        elif min_right_ping_distance < min_left_ping_distance:
            turn_angle_delta = right_turn_angle_delta
            if turn_angle_delta != prev_turn_angle_delta:
                print("food right")
        elif min_left_ping_distance < min_right_ping_distance:
            turn_angle_delta = left_turn_angle_delta
            if turn_angle_delta != prev_turn_angle_delta:
                print("food left")           
        else:
            turn_angle_delta = forward_angle_delta
            if turn_angle_delta != prev_turn_angle_delta:
                print("food forward")
    step_index += 1
    step += 1

I2C_servo_bus.deinit()
print("exiting")
sys.exit(0)




