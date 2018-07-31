# Sonar distance to nearby object.
# Using Parallax PING))) ultrasonic distance sensor.
# https://www.parallax.com/sites/default/files/downloads/28015-PING-Sensor-Product-Guide-v2.0.pdf

import time
import RPi.GPIO as GPIO

GPIO.setmode(GPIO.BOARD)

timeout = 0.020

def ping():

        GPIO.setup(11, GPIO.OUT)
        #cleanup output
        GPIO.output(11, 0)

        time.sleep(0.000002)

        #send signal
        GPIO.output(11, 1)

        time.sleep(0.000005)

        GPIO.output(11, 0)

        GPIO.setup(11, GPIO.IN)
        
        goodread=True
        watchtime=time.time()
        while GPIO.input(11)==0 and goodread:
                starttime=time.time()
                if (starttime-watchtime > timeout):
                        goodread=False

        if goodread:
                watchtime=time.time()
                while GPIO.input(11)==1 and goodread:
                        endtime=time.time()
                        if (endtime-watchtime > timeout):
                                goodread=False
        
        if goodread:
                duration=endtime-starttime
                distance=duration*34000/2
                return distance
        else:
                return -1.0

if __name__ == '__main__':
    print ping()