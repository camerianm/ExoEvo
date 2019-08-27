from mineralDB import minerals
import getall as get

def planets_from_summary():
	file = 'cumulative.csv'
	f = open(file,'r')
	lines = f.readlines()
	f.close()

	files = {}

	for line in lines:
		temp = line.split(',')
		if temp[0] not in files:
			files[temp[0]] = {}
		files[temp[0]][temp[1]] = float(temp[2])

	for file in files:
		files[file]['composition'] = {}
		for mineral in minerals.keys():
			if (mineral in files[file]) and (files[file][mineral]>0):
					files[file]['composition'][mineral] = files[file][mineral]
			try:
				del files[file][mineral]
			except KeyError:
				pass
		files[file]['composition'] = get.adds_up(files[file]['composition'])
		# print("\n" + file + "\n" + "\n".join("{}: {}".format(k, v) for k, v in files[file]['composition'].items()) + "}")
	return(files)
