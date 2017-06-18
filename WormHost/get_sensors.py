import requests, json

url = "http://192.168.0.25:5000/robot/sensors"
r = requests.get(url)
print r.json()

