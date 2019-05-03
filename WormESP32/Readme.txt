C elegans worm robot with ESP32.
A distance sensor simulates food foraging by locating and moving toward the nearest object.

Power up the ESP32 and connect to its network: Robot_Worm
ESP32 IP address should be 192.169.4.1

Transfer python and segment_angles.txt to ESP32 using ftp.

FTP:
For linux:
In ftp connection, issue "passive" command to avoid attempt to open a server-to-client data connection.
This connection will be refused.
Cygwin:
Get ftp by installing inetutils package (also gets telnet).
FileZilla instructions:
https://github.com/loboris/MicroPython_ESP32_psRAM_LoBo/wiki/ftpserver

Login with telnet: 
telnet 192.168.4.1
login=micro, password=python

Run scripts:
1. To straighten robot:
exec(open("servo_reset.py").read())
2. To test servos:
exec(open("servo_test.py").read())
3. Run no food sensor:
exec(open("worm_run_segment_angles.py").read())
4. Run with laser distance food sensor:
exec(open("worm_run_segment_angles_laser.py").read())
5. Run with ultrasonic distance food sensor:
exec(open("worm_run_segment_angles_sonic.py").read())

The maximim number of steps can be limited by creating a file named "max_steps.txt" containing the value.

Note: power cycle board if you encounter python errors after repeated runs.

