# Worm simulation.
import wormpy
wormpy.init()
wormpy.step()
angles = []
activations = []
for segment in range(12) :
    angles.append(wormpy.get_segment_angle(segment))
    activations.append(wormpy.get_ventral_muscle_activation(segment))
print angles
print activations
wormpy.term()
