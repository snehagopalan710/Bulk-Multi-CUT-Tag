import linecache
import argparse
import os.path

# Argument parsing:
def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-odir', type=str, required=True,
						help='output folder path')
	parser.add_argument('-f', type=str, required=True,
						help='filename')
	parser.add_argument('-endbp', type=int, required=True,
						help='end range that count as duplicates')
	parser.add_argument('-perc_cutoff', type=float, required=True,
						help='keep ab pairs with more than this percentage')
	parser.add_argument('-aba', nargs='+', required=True,
						help='antibody name + barcodes for antibody a')
	parser.add_argument('-abb', nargs='+', required=True,
						help='antibody name + barcodes for antibody b')
	args = parser.parse_args()
	return args

# Add new dup to the dictionary for storing duplicates info
def adddupdict(dupdict, umi):
	info = umi.split(':')
	dupdict['start'].add(info[1])
	dupdict['end'].add(int(info[2]))
	dupdict['cb'].add(info[3])

# is this line duplicates?
def isdup(dupdict,umi,endbp):
	decision = False
	info = umi.split(":")
	start = info[1]
	end = int(info[2])
	cb = info[3]
	if cb in dupdict['cb'] and start in dupdict['start']:
		if end in dupdict['end']:
			decision = True
		else:
			endextent = []
			for i in range(1,endbp+1):
				endextent += [x+i for x in dupdict['end']]
				endextent += [x-i for x in dupdict['end']]
			if end in endextent:
				decision = True
	return decision

def selectFrag(dups, perc_cutoff):
	loccb = {}
	ab = {}
	for umi, count in dups.items():
		info = umi.split(':')
		umiL = '%s:%s:%s:%s' % (info[0],info[1],info[2],info[3])
		umiR = '%s:%s' % (info[4],info[5])
		if loccb.has_key(umiL):
			loccb[umiL] += count
		else:
			loccb[umiL] = count
		if ab.has_key(umiR):
			ab[umiR] += count
		else:
			ab[umiR] = count
	loccb_most = max(loccb, key=loccb.get)
	ab_most = max(ab, key=ab.get)
	if ab[ab_most] > perc_cutoff*sum(ab.values()):
		fragment = loccb_most.split(':') + ab_most.split(':')
	else:
		fragment = None
	return (fragment, ab.values())

def abwrite(fragment, abaout, abbout, abalist, abblist):
	ab1 = fragment[4]
	ab2 = fragment[5]
	cut1 = [fragment[i] for i in [0,1,3]]
	cut2 = [fragment[i] for i in [0,2,3]]
	cuts = {ab1:cut1, ab2:cut2}
	for ab, cut in cuts.items():
		if ab in abalist:
			abaout.write('%s\n' % '\t'.join(cut))
		elif ab in abblist:
			abbout.write('%s\n' % '\t'.join(cut))

def run(args):
	odir = args.odir
	inFile = args.f
	inf = open(inFile, 'r')
	endbp = args.endbp
	perc_cutoff = args.perc_cutoff
	abaname = args.aba[0]
	abalist = args.aba[1:]
	abbname = args.abb[0]
	abblist = args.abb[1:]

	# fragment
	outFile = os.path.join(odir, inFile.replace('.txt','_dedup.txt'))
	out = open(outFile,'wb')
	# antibody barcodes frequency within each dups
	freqFile = os.path.join(odir, inFile.replace('.txt','_freq.txt'))
	freq = open(freqFile, 'wb')
	# cutsite tracks for antibody, barcodes as file names.
	abaFile = os.path.join(odir, inFile.replace('.txt','_'+abaname+'.txt'))
	abaout = open(abaFile, 'wb')
	abbFile = os.path.join(odir, inFile.replace('.txt','_'+abbname+'.txt'))
	abbout = open(abbFile, 'wb')

	# dedup
	dupdict = {'start':set(),'end':set(),'cb':set()}
	dups = {}
	# first line
	line = inf.readline().split()
	count = int(line[0])
	umi = line[1]
	adddupdict(dupdict,umi)
	dups[umi] = count

	while 1:
		line = inf.readline().split()
		if not line:
			break
		count = int(line[0])
		umi = line[1]
		if isdup(dupdict, umi, endbp):
			adddupdict(dupdict,umi)
			dups[umi] = count
		else:
			# write the previous dup to file.
			(fragment, frequency) = selectFrag(dups, perc_cutoff)
			freq.write('%s\n' % '\t'.join(map(str,frequency)))
			if fragment != None:
				out.write('%s\n' % '\t'.join(fragment))
				# write citsite tracks
				abwrite(fragment, abaout, abbout, abalist, abblist)

			# start a new dups
			dupdict = {'start':set(),'end':set(),'cb':set()}
			dups = {}
			adddupdict(dupdict,umi)
			dups[umi] = count

	# write the last dups
	(fragment, frequency) = selectFrag(dups, perc_cutoff)
	freq.write('%s\n' % '\t'.join(map(str,frequency)))
	if fragment != None:
		out.write('%s\n' % '\t'.join(fragment))

	out.close()
	freq.close()
	abaout.close()
	abbout.close()
	return

if __name__ == "__main__":
	args = get_args()
	run(args)

