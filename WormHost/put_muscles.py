import requests, json

url = "http://192.168.0.25:5000/robot/muscles"
payload = json.load(open("muscles.json"))
headers = {'content-type': 'application/json', 'Accept-Charset': 'UTF-8'}
r = requests.post(url, data=json.dumps(payload), headers=headers)
print r.json()

