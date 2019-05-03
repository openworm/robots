# Run with segment angles file "segment_angles.txt".
# Maximum steps can optionally be set by creating file
# "max_steps.txt" containing value.

from machine import PWM, I2C, Pin
from math import pi, radians
from BOTS import *  #import all custom code for the board
import sys
from time import sleep

# Parameters.
angle_amplifier = 1.0  # Increases wave amplitude
servo_delay = .05
angle_bias = 1.75
forward_angle_delta = .327
head_swing_step_index_restart = 10

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

# Run.
head_swing_count = 0
angles_array_len = len(angles_array)
step = 0
step_index = 0
while step_index < angles_array_len and (max_steps == -1 or step < max_steps):
    #print("step=", step)
    for segment in range(8):
        position = (int)(((float)(angles_array[step_index][segment]) * angle_amplifier) + angle_bias + forward_angle_delta)
        try:
            robot.set_servo (radians(position), 7 - segment)              
            #print("servo=", segment, "position=", position)
        except:
            print("error setting servo angle")
    sleep(servo_delay)

    if step_index > 0 and step_index < (angles_array_len - 1):
        curr_angle = (int)(((float)(angles_array[step_index][0]) * angle_amplifier) + angle_bias + forward_angle_delta)
        next_angle = (int)(((float)(angles_array[step_index + 1][0]) * angle_amplifier) + angle_bias + forward_angle_delta)
        
        # Restart step index?
        if curr_angle > 0 and next_angle <= 0:
            head_swing_count += 1
            if head_swing_step_index_restart != -1 and head_swing_count >= head_swing_step_index_restart:
                step_index = -1
                head_swing_count = 0
                #print("step restart")

    step_index += 1
    step += 1

I2C_servo_bus.deinit()
print("exiting")
sys.exit(0)




