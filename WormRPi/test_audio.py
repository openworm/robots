#!/usr/bin/env python3
import pyaudio
import wave

CHUNK = 4096
FORMAT = pyaudio.paInt32
CHANNELS = 1
RATE = 16000 
RECORD_SECONDS = 10
WAVE_OUTPUT_FILENAME = "test.wav"

p = pyaudio.PyAudio()

stream = p.open(format=FORMAT,
            channels=CHANNELS,
            rate=RATE,
            input=True,
            frames_per_buffer=CHUNK) 

frames = []

for i in range(0, int(RATE / CHUNK * RECORD_SECONDS)):
     data = stream.read(CHUNK)
     #print(data)
     #frames.append(data)

stream.stop_stream()
stream.close()
p.terminate()

