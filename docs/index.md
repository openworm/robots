---
title: OpenWorm Robot build instructions
layout: default
---

# OpenWorm Crawling Robot v 0.1

By Shane Gingell (Shane@OutoftheBOTS.com.au), Tom Portegys
(portegys@openworm.org) and the OpenWorm Foundation
(info@openworm.org) Version 1.0 August 25th, 2017. 

We have an explainer video [here on Youtube](https://youtu.be/6gbrbjGFhD4).

![Finished product image](https://guides.github.com/activities/hello-world/branching.png)

## Expected Costs

This can be made for around US$150 dollars.

## Purchased Parts

* __Wheels__ : This part isn’t critical and easily exchanged for what is
  easily available locally. We used 25mm RC aircraft wheels and then
  drilled out the axle holes to 3mm purchased [from here](https://www.aliexpress.com/item/100pcs-Small-Light-Foam-Tail-Wheel-Diam-25mm-Thickness-11mm-Shaft-hole-2mm-For-RC-Remote/32705923540.html?spm=2114.01010208.3.10.b4OseI&ws_ab_test=searchweb0_0,searchweb201602_4_10152_10065_10151_10130_10068_5010017_10136_10157_10137_10060_10138_10131_10155_10062_10132_10156_10133_437_10154_10056_10055_10054_10059_303_100031_10099_10103_10102_10096_10147_10052_10053_10050_10107_10142_10051_5020019_10084_10083_10080_10082_10081_10110_519_10175_10111_10112_10113_10114_10037_10182_10185_10032_10078_10079_10077_10073_10123_142,searchweb201603_13,ppcSwitch_4&btsid=c1dbe6fc-0b51-4f03-8d57-6b99288429d4&algo_expid=840bed72-95b2-42cf-bef2-a14471322b64-1&algo_pvid=840bed72-95b2-42cf-bef2-a14471322b64).
* __Nuts and bolts__ : We used 3mm x 20mm bolts with nyloc nuts. It is critical that nyloc nuts are used.
* __PMW board__ : This drives all the servos. This is one [from
  Adafruit](https://www.adafruit.com/product/815), but we used a model
  with the same chip so the same software could be used to drive
  it. [This
  model](http://www.ebay.com.au/itm/PCA9685-16-Channel-12-bit-PWM-Servo-motor-Driver-I2C-Module-For-Arduino-Robot-P6/182391317074?_trkparms=aid%3D222007%26algo%3DSIM.MBE%26ao%3D2%26asc%3D20140106155344%26meid%3D55e89eed8cb54264bcdedb85b836058c%26pid%3D100005%26rk%3D6%26rkt%3D6%26sd%3D262714356063&_trksid=p2047675.c100005.m1851)
  is better built and cheaper. An alternative can be [found here](http://www.ebay.com/itm/like/191931347026?chn=ps&dispItem=1).
* __Battery__ : We used a RC 7.4v 2S battery (you will need a
  discharge rate of min 3c but most RC batteries are like 25c so
  generally lots more than needed) the capacity doesn’t matter so much
  but we used 2000mah as it seemed to be a nice balance between size
  and lasting long enough before going flat. [This
  model](http://www.ebay.com.au/itm/1x-2x-2000mAh-7-4V-25C-LiPo-Battery-2S-Balance-Charger-for-Syma-X8HG-Drone-/332136206449?var=&hash=item4d54dc5071:m:mdi5sJNL8wo9I-nlZyge7tg)
  was chosen because it also comes with a charger. An alternative can
  be [found
  here](http://www.ebay.com/itm/like/351738012331?chn=ps&dispItem=1).
* __Voltage step down converters__ : We had to use __two different
  models__ for these. __The first__ is [an adjustable micro
  360](http://www.ebay.com.au/itm/5pcs-Mini-360-DC-DC-Buck-Converter-Step-Down-Module-4-75V-23V-to-1V-17V-/122098217064?hash=item1c6d9ef068:g:4vAAAOSwV0RXvWXl)
  for the RPi and adjusted it to 5v. An alternative source for this
  part can be [found
  here](http://www.ebay.com/itm/like/351738012331?chn=ps&dispItem=1). __The
  second__ is for the servos. We used a larger amp adjustable converter
  and adjusted it to 6v. This model can be [found here](http://www.ebay.com.au/itm/XL4015-5A-DC-DC-Step-Down-Buck-Converter-Module-Power-Supply-LED-Lithium-Charger/192234308138?var=492161169440&_trkparms=aid%3D222007%26algo%3DSIM.MBE%26ao%3D2%26asc%3D20140106155344%26meid%3D538f67598ae24dbba1f6594d01d6020b%26pid%3D100005%26rk%3D6%26rkt%3D6%26sd%3D332197058401&_trksid=p2047675.c100005.m1851). An alternative source can be [found here](http://www.ebay.com/itm/like/112237134540?chn=ps&dispItem=1).
* __Servos__ : We used the Corona brand model CS-939MG brought
  directly [from the factory
  here](https://wxrgdz.en.alibaba.com/product/60608230195-804446832/Corona_CS939MG_Analog_Metal_Gear_Servo_2_5kg_0_14sec_12_5g_for_RC_toys.html). An
  alternative source for this part can be [found from Hobby King
  here](https://hobbyking.com/en_us/corona-939mg-metal-gear-servo-2-5kg-0-14sec-12-5g.html?countrycode=US&gclid=Cj0KCQjwoZTNBRCWARIsAOMZHmGEeijs9xE-hcR2DAJ2OFpjulCqenN-gff_-upFRulK0nbGEl4q_wQaAmjZEALw_wcB&gclsrc=aw.ds).
* __Low battery alarm__ : From [Ebay](http://www.ebay.com.au/itm/LED-RC-Lipo-Li-ion-Battery-Low-Voltage-Meters-Alarms-Test-Buzzer-Monitor-KB/232225019917?_trkparms=aid%3D555017%26algo%3DPL.CASSINI%26ao%3D1%26asc%3D20160706105120%26meid%3D94998925138d4c3380cf4f04fcb59951%26pid%3D100508%26rk%3D1%26rkt%3D1%26&_trksid=p2045573.c100508.m3226). Alternative source from [Amazon](https://www.amazon.com/dp/B00SCJOITA/ref=asc_df_B00SCJOITA5145422/?tag=hyprod-20&creative=395033&creativeASIN=B00SCJOITA&linkCode=df0&hvadid=198096709148&hvpos=1o1&hvnetw=g&hvrand=9385563450414889614&hvpone=&hvptwo=&hvqmt=&hvdev=c&hvdvcmdl=&hvlocint=&hvlocphy=9033288&hvtargid=pla-349312979893).
* __Raspberry Pi Zero W__ : [From AdaFruit](https://www.adafruit.com/product/3410).
* __Wires and Miscellany__ : You will also need some __dupont jumper wires__ and some
  __3 amp power wires__ and __soldering iron__ and __multimeter__ for
  adjusting the voltage on buck converters.


## Step 1 -- Print the 3D printer Parts

**NOTE: Important stuff to keep in mind**

If you want to print these yourself, we recommend buying the [Tevo
Slimbot here](https://www.3dprintersbay.com/tevo-slimbot-dual-head-large-bed-auto-leveling-3d-printer-diy-kit)

![Branching](https://guides.github.com/activities/hello-world/branching.png)

1. __Worm body__ : [STL file from Github here](https://github.com/openworm/robots/blob/master/3D%20Printing%20Shapefiles/Worm%20Body%2011-8-17.stl).
1. __Worm head__ : [STL file from Github here](https://github.com/openworm/robots/blob/master/3D%20Printing%20Shapefiles/Worm_Head.stl).
1. __Worm PWM stand__ : [STL file from Github here](https://github.com/openworm/robots/blob/master/3D%20Printing%20Shapefiles/Worm_PWM_Stand.stl).
1. __Worm RPi Stand__ : [STL file from Github here](https://github.com/openworm/robots/blob/master/3D%20Printing%20Shapefiles/Worm_RPi_Stand.stl).
1. __Worm battery stand__ : [STL file from Github here](https://github.com/openworm/robots/blob/master/3D%20Printing%20Shapefiles/Worm_battery_stand.stl).

## Step 2 -- Build the body

![Branching](https://guides.github.com/activities/hello-world/branching.png)

1. Subtask 1
1. Subtask 2
1. Subtask 3
1. Subtask 4
1. Subtask 5

## Step 3 -- Add the electronics

![Branching](https://guides.github.com/activities/hello-world/branching.png)

1. Subtask 1
1. Subtask 2
1. Subtask 3
1. Subtask 4
1. Subtask 5

## Step 4 -- Upload the code 

![Branching](https://guides.github.com/activities/hello-world/branching.png)

1. Subtask 1
1. Subtask 2
1. Subtask 3
1. Subtask 4
1. Subtask 5

## Step 5 -- Integrate & Run

![Branching](https://guides.github.com/activities/hello-world/branching.png)

1. Subtask 1
1. Subtask 2
1. Subtask 3
1. Subtask 4
1. Subtask 5
