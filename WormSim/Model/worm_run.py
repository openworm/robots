# Worm simulation.
import sys
sys.path.insert(1, 'Release')
import wormpy
wormpy.init()
wormpy.step()
activations = []
for segment in range(12) :
    activations.append(wormpy.get_ventral_muscle_activation(segment));
print activations
wormpy.term()
