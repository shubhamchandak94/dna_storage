import math

counts = [3,12];
for i in range(120):
	counts.append(3*counts[-1]+3*counts[-2])
	print 'Number of bases: ' + str(i+3)
	print 'Number of valid strings: ' + str(counts[-1])
	print 'Max number of bits: ' + str(math.floor(math.log(counts[-1],2)))
	print 'Rate in bits/base: ' + str(math.floor(math.log(counts[-1],2))*1.0/(i+3))
	print
