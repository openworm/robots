# Worm simulation test.
import sys
sys.path.insert(1, 'Release')
import wormpy
wormpy.init_test()
wormpy.amplify_forward(1)
#wormpy.amplify_turn(1.25)
go = 1
while go == 1 :
   go = wormpy.step()
   angles = []
   activations = []
   for segment in range(12) :
        angles.append(wormpy.get_segment_angle(segment))
        activations.append(wormpy.get_ventral_muscle_activation(segment))
   print(angles)
   #print activations
wormpy.term()

