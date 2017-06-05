// WormSim API.
int init();
int init_test();
int step();
float get_dorsal_muscle_activation(int segment);
float get_ventral_muscle_activation(int segment);
void term();

// Motor neuron forward amplification.
// 1=no amplification.
void amplify_forward(float amplifier);

// Motor neuron turn amplification.
// 1=no amplification.
void amplify_turn(float amplifier);
