# Worm simulation test.
import sys
sys.path.insert(1, 'Release')
import wormpy
wormpy.init_test()
go = 1
while go == 1 :
   go = wormpy.step()
   activations = []
   for segment in range(12) :
        activations.append(wormpy.get_ventral_muscle_activation(segment));
wormpy.term()

