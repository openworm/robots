# Listen for food audio.

import os
from scipy.io.wavfile import read
import numpy

def listen():

    # Record audio.
    result = os.system("arecord -d 1 -D dmic_sv -c2 -r 48000 -f S32_LE -t wav -V mono food.wav >/dev/null 2>/dev/null")

    # Extract amplitudes.
    fs, data = read('food.wav')
    data = numpy.array(data).astype(float)
    data_size = len(data)

    # Extract mean of max amplitudes.
    amplitude_sample_size = 50
    amplitudes = []
    idx = 0
    while idx < data_size:
        if data[idx][0] < 0:
            amplitudes.append((float)(abs(data[idx][0])))
        idx += 1
    amplitudes.sort(key=float)
    amplitude_sum = 0.0
    amplitudes_size = len(amplitudes)
    idx = amplitudes_size - amplitude_sample_size
    while idx < amplitudes_size:
        amplitude_sum += amplitudes[idx]
        idx += 1
    return ((int)((amplitude_sum / (float)(amplitude_sample_size)) / 1000000))

if __name__ == '__main__':
    print listen()