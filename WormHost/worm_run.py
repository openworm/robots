# Worm simulation.
import wormpy
wormpy.init()
wormpy.step()
activations = []
for segment in range(12) :
    activations.append(wormpy.get_ventral_muscle_activation(segment));
print activations
wormpy.term()
