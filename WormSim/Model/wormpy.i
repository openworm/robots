/* File: wormpy.i */
%module wormpy

%{
#define SWIG_FILE_WITH_INIT
#include "wormsim.h"
%}

// wormsim API.
int init();
int init_test();
int step();
float get_segment_angle(int segment);
float get_dorsal_muscle_activation(int segment);
float get_ventral_muscle_activation(int segment);
void term();
void amplify_forward(float amplifier);
void amplify_turn(float amplifier);


