# Generate segment angles file.
# Options: <number of steps (formats: n or m:n where m leading steps are excluded)> [<output file name (defaults to segment_angles.txt)>]

import sys
sys.path.insert(1, 'WormSim/Model/x64/Release')
import wormpy

# Parameters
amplifier = 1
stride = 1

if len(sys.argv) == 1:
    print("Options: <number of steps (formats: n or m:n where m leading steps are excluded)> [<output file name> (defaults to segment_angles.txt)]")
    sys.exit(1)
steps = 0
exclude_steps = 0
vals = sys.argv[1].split(":")
if len(vals) == 1:
    steps = (int)(vals[0])
else:
    exclude_steps = (int)(vals[0])
    steps = (int)(vals[1])

# Open output file.
filename = "segment_angles.txt"
if len(sys.argv) == 3:
    filename = sys.argv[2]
output_file = open(filename, "w")

motor_file = open("motors.txt", "w")
muscle_file = open("muscles.txt", "w")

wormpy.init()
step_count = 0
for step in range(steps):
    wormpy.step(0.0)
    for s in range(8):
        motor_dorsal = wormpy.get_dorsal_motor_activation(s)
        if motor_dorsal < 0.0:
            motor_dorsal_str = "{0:.2f}".format(motor_dorsal)
        else:
            motor_dorsal_str = "+{0:.2f}".format(motor_dorsal)
        motor_ventral = wormpy.get_ventral_motor_activation(s)
        if motor_ventral < 0.0:
            motor_ventral_str = "{0:.2f}".format(motor_ventral)
        else:
            motor_ventral_str = "+{0:.2f}".format(motor_ventral)
        motor_file.write(motor_dorsal_str)
        motor_file.write("/")
        motor_file.write(motor_ventral_str)
        motor_file.write(" ")
        muscle_dorsal = wormpy.get_dorsal_muscle_activation(s)
        if muscle_dorsal < 0.0:
            muscle_dorsal_str = "{0:.2f}".format(muscle_dorsal)
        else:
            muscle_dorsal_str = "+{0:.2f}".format(muscle_dorsal)
        muscle_ventral = wormpy.get_ventral_muscle_activation(s)
        if muscle_ventral < 0.0:
            muscle_ventral_str = "{0:.2f}".format(muscle_ventral)
        else:
            muscle_ventral_str = "+{0:.2f}".format(muscle_ventral)
        muscle_file.write(muscle_dorsal_str)
        muscle_file.write("/")
        muscle_file.write(muscle_ventral_str)
        muscle_file.write(" ")
    motor_file.write("\n")
    muscle_file.write("\n")
    if (step % stride) == 0:
        if step_count >= exclude_steps:
            for segment in range(8) :
                angle = (int)(wormpy.get_segment_angle(segment) * amplifier)
                output_file.write(str(angle))
                if segment < 7:
                    output_file.write(" ")
            output_file.write("\n")
        step_count = step_count + 1

wormpy.term()
output_file.close()
muscle_file.close()
motor_file.close()
sys.exit(0)
