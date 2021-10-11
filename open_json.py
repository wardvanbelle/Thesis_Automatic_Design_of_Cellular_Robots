import json

with open('building_blocks_database.json') as f:
	data = json.load(f)

for state in data:
	print(state)