import flask
from flask import Flask, json, jsonify, request, abort
import json
from Servo import Servo
from time import sleep

robot = Servo()
servo_list = [0, 1, 2, 3, 4, 5, 6, 7]
servo_delay = 1

app = Flask(__name__)

@app.route('/')
def index():
    return 'Routes: robot/sensors, and robot/muscles'

@app.route('/robot/sensors', methods=['GET'])
def get_sensors():
    with open("sensors.json") as sensors_file:
        sensors_data = json.load(sensors_file)
    return jsonify({'sensors': sensors_data})

@app.route('/robot/muscles', methods=['POST'])
def put_muscles():
    if not request.json or not 'muscles' in request.json:
        abort(400)
    muscles = request.json.get('muscles', "")
    for i in servo_list:
        robot.set_servo ((11 - i), muscles[i])
    #sleep(servo_delay)
    with open('muscles.json', "a") as muscles_file:
        muscles_file.write('{ "muscles" : %s }\n' % muscles)
    return "OK", 201

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0')
