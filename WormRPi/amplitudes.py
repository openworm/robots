"""Usage:
    amplitudes WAV_FILE

    Returns the (linear) amplitude of the signal inside the wav file, sample by sample.
"""
from __future__ import division
import docopt

import scipy.io.wavfile

MAX_WAV16_AMP = 32767  # = 2**15-1  # because wav16 is signed (wav8 isn't)


def main():
    args = docopt.docopt(__doc__)

    one_pos_to_neg_crossing_path = r"W:\materials\audio\one_pos_to_neg_crossing.wav"
    if False: args['WAV_FILE'] = one_pos_to_neg_crossing_path

    rate, amp_arr = scipy.io.wavfile.read(args['WAV_FILE'])
    for amp in (amp_arr / MAX_WAV16_AMP):
        print(amp)

if __name__ == '__main__':
    main()
