C elegans worm robot Raspberry Pi Zero W code.

Download Raspbian to SD boot card.

To enable wireless, create wpa_supplicant.conf in boot directory with this content:
country=<country code, e.g. US>
ctrl_interface=DIR=/var/run/wpa_supplicant GROUP=netdev
update_config=1

network={
    ssid="<network name>"
    psk="<network password>"
}

To enable ssh, create file named "ssh" in boot directory.

To find RPi IP: nmap -sn 192.168.1.0/24

Enable I2C:
https://www.abelectronics.co.uk/kb/article/1/i2c--smbus-and-raspbian-linux

Login with ssh:
ssh <RPi IP>
login=pi, password=raspberry

Setup:
1. Create $HOME/robot directory on RPi:
mkdir $HOME/robot
2. Copy python and segment_angles.txt files to $HOME/robot.

Run scripts: servo_reset.py (straighten robot), worm_run.py, worm_run_angles.py, worm_run_segment_angles.py

Ultrasonic distance sensor:

Parallax Ping Ultrasonic distance sensor:
https://www.parallax.com/sites/default/files/downloads/28015-PING-Sensor-Product-Guide-v2.0.pdf
See: sonar.py

Deprecated audio sensor:
 
Audio microphone installation instructions:
https://learn.adafruit.com/adafruit-i2s-mems-microphone-breakout/
https://pinout.xyz/pinout/i2c#
See: food_listener.py

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

To enable http communications with host:
1. Copy json files to $HOME/robot.
2. Start flask web server: python flask_app.py
