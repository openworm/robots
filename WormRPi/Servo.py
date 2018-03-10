# Copyright (c) 2017 Out of the BOTS
# http://outofthebots.com.au/
# Author: Shane Gingell
#
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, 
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, 
#    this list of conditions and the following disclaimer in the documentation 
#    and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
# OR TORT INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



# Servo class provides very basic functionality for controlling 
# sevoing using the PCA9685 chipset with a RPi. 


import smbus #import the System Mangement package for I2C coms
import time

class Servo :

  global bus
  global address
  bus = smbus.SMBus(1) #RPi3 uses bus 1 but if using bus 0 then change here
  address = 0x40 #factory standard I2C slave address change if needed
  
  @classmethod
  def __init__(self): #initizing will turn on all servos and set to poiition 0
    bus.write_byte_data(address, 0x01, 0x04)
    bus.write_byte_data(address, 0x00, 0x01)
    time.sleep(0.005)

    oldreg0 = bus.read_byte_data(address, 0x00)
    newreg0 = (oldreg0 & 0x7F) | 0x10
    bus.write_byte_data(address, 0, newreg0)

    prescale = int((25000000/(4096*50))-1)
    bus.write_byte_data(address, 0xFE, prescale)

    bus.write_byte_data(address, 0, oldreg0)
    time.sleep(0.005)
    bus.write_byte_data(address, 0, oldreg0 | 0x80)

    bus.write_byte_data(address, 0xFA, 0)
    bus.write_byte_data(address, 0xFB, 0)
    bus.write_byte_data(address, 0xFC, 0x33)
    bus.write_byte_data(address, 0xFD, 0x01)  
  
  def set_servo (self,servo_number, position): #routine for setting servo positions. no error checking
	step = 202/90
	if position < -90 : position = -90
	if position > 90 : position = 90
	if position == 0 : offtime = 307	
	else  : offtime = 307 + (step * position)    
	startreg = servo_number * 4 + 8    
	bus.write_byte_data(address, startreg, offtime & 0xff)    
	bus.write_byte_data(address, startreg+1, offtime >> 8)

#still much funtionality missing especally a shutdown down function.
  




