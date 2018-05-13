C elegans worm robot Raspberry Pi Zero W code.

To enable http communications with host:
1. Create $HOME/robot directory on RPi.
2. Copy python and json files to $HOME/robot.
3. Start flask web server: python flask_app.py

Audio microphone installation instructions:
https://learn.adafruit.com/adafruit-i2s-mems-microphone-breakout/
https://pinout.xyz/pinout/i2c#

To set audio volume sensitivity:
Copy asoundrc.txt to $HOME/.asoundrc on the robot.
Run arecord first to "prime the pump"
Then either:
  alsamixer
    F6 to select snd_rpi_simple_card
    F4 to set volume
or
  amixer -M set Capture 50%
Finally to make sensitivity persistent:
  sudo alsactl store 1
