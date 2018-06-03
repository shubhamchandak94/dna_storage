import itertools
import math

for n in [2,3,4,5,6,7,8,9,10]:
	chars = 'ACGT'
	num_valid = 0
	for item in itertools.product(chars, repeat=n):
		s = ''.join(item)
		if s[0] == s[1]:
			continue
		if s[-2] == s[-1]:
			continue
		if s.find('AAA') != -1:
			continue
		if s.find('CCC') != -1:
			continue
		if s.find('GGG') != -1:
			continue
		if s.find('TTT') != -1:
			continue
		num_valid += 1
	print n
	print num_valid
	print 1.0*num_valid/(4**n)
	print math.floor(math.log(num_valid,2))/n		

for n in [2,3,4,5,6,7,8,9,10,11]:
	chars = 'ACGT'
	num_valid = 0
	for item in itertools.product(chars, repeat=n):
		s = ''.join(item)
		if s.find('AAA') != -1:
			continue
		if s.find('CCC') != -1:
			continue
		if s.find('GGG') != -1:
			continue
		if s.find('TTT') != -1:
			continue
		num_valid += 1
	print n
	print 0.75*num_valid
	print 0.75*num_valid/(4**n)
	print math.floor(math.log(0.75*num_valid,2))/n
