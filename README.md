# robots
C. elegans robots.

The goal of this project is twofold:
1. To build a robot that simulates the neuromuscular function of a C. elegans nematode worm.
2. To specify a kit of parts and instructions that will allow a student to also build the robot.

The robot is controlled by a Raspberry Pi/ESP32 processor that contains a recorded simulation of the worm's
neuromuscular system (see references). The robot's body is a sequence of segments that mutually exert simulated muscle
contractions impemented by servos.

Folders:
1. WormHost: PC code to communicate with the onboard Raspberry Pi.
2. WormRPi: Raspberry Pi onboard code.
3. WormESP32: ESP32 onboard code.
4. WormSim: C. elegans neuromuscular simulator.
5. assembly: parts list, 3D printing .stl shape files, and assembly instructions.
6. docs: robot description.

References:
Boyle, Berri and Cohen, "Gait modulation in C. elegans: an integrated neuromechanical model", Front. Comput. Neurosci., 2012.
Eduardo J. Izquierdo and Randall D. Beer, "An Integrated Neuromechanical Model of Steering in C. elegans", ECAL15