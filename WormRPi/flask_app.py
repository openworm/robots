import flask
from flask import Flask, json, jsonify, request
from Servo import Servo

app = Flask(__name__)

robot_servo = Servo()
servo_list = [0, 1, 2, 3, 4, 5, 6, 7]

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
    muscles = {
        'muscles': request.json.get('muscles', "")
    }
    with open('muscles.json', 'w') as muscles_file:
        json.dump(muscles, muscles_file)
        for i in servo_list:
            robot_servo.set_servo ((11 - i), muscles[i])
        return jsonify({'muscles': muscles}), 201

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0')

