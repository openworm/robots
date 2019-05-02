# Run with segment angles file.
# Options: [<input file name (defaults to segment_angles.txt)>]
from machine import PWM, I2C, Pin
from math import pi, radians
from BOTS import *  #import all custom code for the board
import sys
from time import sleep
import vl53l0x   # Laser distance sensor

# Parameters.
#angle_amplifier = 1.25
angle_amplifier = 1.0
servo_delay = .05
right_ping_angles = [ -35.0 ]
left_ping_angles = [ 35.0 ]
ping_stride = 2
ping_before_pause = 5
ping_after_pause = 3
angle_bias = 1.75
left_turn_angle_delta = 2.0
right_turn_angle_delta = -2.0
head_swing_step_restart = 5
goal_distance = 10.0
max_valid_food_distance = 6000

# Load segment angles file.
print("Loading segment angles...", end='')
angles_array = []
filename = "segment_angles.txt"
if len(sys.argv) == 2:
    filename = sys.argv[1]
with open(filename, "r") as f:
    for line in f:
        vals = (line[:-1]).split(" ")
        angles = []
        for angle in vals:
            angles.append(int(angle))
        angles_array.append(angles);
print("done")

# Initialize servos.
exec(open("servo_reset.py").read())
sleep(2)

# Run.
right_ping_angles_len = len(right_ping_angles)
right_ping_distances = [ max_valid_food_distance + 1 ] * right_ping_angles_len
left_ping_angles_len = len(left_ping_angles)
left_ping_distances = [ max_valid_food_distance + 1 ] * left_ping_angles_len
head_swing_side = 'neutral'
head_swing_count = 0
right_swing_count = 0
left_swing_count = 0
turn_angle_delta = 0.0
angles_array_len = len(angles_array)
step = 0
while step < angles_array_len:
    #print("step=", step)
    I2C_servo_bus = I2C(0, speed=100000, sda=27, scl=26)
    robot = Servo(I2C_servo_bus)
    for segment in range(8):
        position = (int)(((float)(angles_array[step][segment]) * angle_amplifier) + angle_bias + turn_angle_delta)
        try:
            robot.set_servo (radians(position), 7 - segment)
            #print("servo=", segment, "position=", position)
        except:
            print("error setting servo angle")
    I2C_servo_bus.deinit()
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
                right_ping_distances = [ max_valid_food_distance + 1 ] * right_ping_angles_len
                left_ping_distances = [ max_valid_food_distance + 1 ] * left_ping_angles_len
                head_swing_side = 'neutral'
                right_swing_count = left_swing_count = 0
                turn_angle_delta = 0.0
                print("step restart")

        # Check right food distance.
        if prev_angle <= 0 and curr_angle < prev_angle:
            if head_swing_side == 'neutral' or head_swing_side == 'left':
                if head_swing_side == 'left':
                    left_swing_count = left_swing_count + 1
                    print("left food distance=", min(left_ping_distances))
                head_swing_side = 'right'
            if (right_swing_count % ping_stride) == 0:
                for i in range(right_ping_angles_len):
                    ping_angle = right_ping_angles[i]
                    curr_angle_delta = abs(curr_angle - ping_angle)
                    if abs(prev_angle - ping_angle) >= curr_angle_delta and abs(next_angle - ping_angle) >= curr_angle_delta:
                        sleep(ping_before_pause)  ###
                        I2C_sensor_bus = I2C(0, speed=100000, sda=32, scl=33)
                        distance_sensor = vl53l0x.VL53L0X(I2C_sensor_bus)
                        sleep(ping_before_pause)
                        try:
                            right_ping_distances[i] = distance_sensor.read()
                        except:
                            print("error sensing distance")
                        I2C_sensor_bus.deinit()
                        print("ping: angle=", curr_angle, "distance=", right_ping_distances[i])  ###
                        if right_ping_distances[i] > max_valid_food_distance:
                            right_ping_distances[i] = max_valid_food_distance + 1
                        sleep(ping_after_pause)

        # Check left food distance.
        if prev_angle >= 0 and curr_angle > prev_angle:
            if head_swing_side == 'neutral' or head_swing_side == 'right':
                if head_swing_side == 'right':
                    right_swing_count = right_swing_count + 1
                    print("right food distance=", min(right_ping_distances))
                head_swing_side = 'left'
            if (left_swing_count % ping_stride) == 0:
                for i in range(left_ping_angles_len):
                    ping_angle = left_ping_angles[i]
                    curr_angle_delta = abs(curr_angle - ping_angle)
                    if abs(prev_angle - ping_angle) >= curr_angle_delta and abs(next_angle - ping_angle) >= curr_angle_delta:
                        sleep(ping_before_pause) ###
                        I2C_sensor_bus = I2C(0, speed=100000, sda=32, scl=33)
                        distance_sensor = vl53l0x.VL53L0X(I2C_sensor_bus)
                        sleep(ping_before_pause)
                        try:
                            left_ping_distances[i] = distance_sensor.read()
                        except:
                            print("error sensing distance")
                        I2C_sensor_bus.deinit()
                        print("ping: angle=", curr_angle, "distance=", left_ping_distances[i])  ###
                        if left_ping_distances[i] > max_valid_food_distance:
                            left_ping_distances[i] = max_valid_food_distance + 1
                        sleep(ping_after_pause)
               
        # Determine direction.
        min_right_ping_distance = min(right_ping_distances)
        min_left_ping_distance = min(left_ping_distances)
        min_right_ping_distance = min_left_ping_distance ###       
        if min_right_ping_distance <= goal_distance or min_left_ping_distance <= goal_distance:
            print("food found")
            break
        elif min_right_ping_distance < min_left_ping_distance:
            print("food to right")
            turn_angle_delta = right_turn_angle_delta
        elif min_left_ping_distance < min_right_ping_distance:
            print("food to left")
            turn_angle_delta = left_turn_angle_delta
        else:
            turn_angle_delta = 0.0
    step += 1

print("exiting")



