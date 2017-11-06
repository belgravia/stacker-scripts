import sys, csv, collections

try:
	psl_in = open(sys.argv[1])
	psl_out = sys.argv[2]
except:
	sys.stderr.write('usage: script.py psl_in psl_out\n')
	sys.exit(1)

Edge = collections.namedtuple('Edge', ['node', 'weight'])

class Node:
	""" A node object that can be used for building directed graphs. Node
	objects have weighted pointers to other nodes that can be added using the
	addNext method. """
	def __init__(self, label, score=None):
		""" Sets the label of the node to an attribute and initializes the list
		of nodes that can reach the current node. """
		self.label = label
		self.nextNodes = {}  # keys are nodes following out edges, values are weight

	def addPrevious(self, prevNode, weight):
		""" Appends another node with an outgoing edge to this node as well as
		the weight of the edge as a tuple to self.previous. """
		self.previous += [Edge(prevNode, weight=weight)]

	def addNext(self, nextNode, weight=1, new=False):
		if nextNode not in self.nextNodes:
			self.nextNodes[nextNode] = 0
			new = True
		self.nextNodes[nextNode] += weight
		return new

class DAG:
	""" A directed acyclic graph. The current methods support finding the
	longest path through the graph. """
	def __init__(self, start, end):
		"""  """
		self.nodes = {}  # {nodelabel: nodeobject}
		self.numPaths = str(1)
		self.start, self.end = start, end
		self.source = Node(str(start))
		self.nodes[str(start)] = self.source

	def addEdges(self, edges):
		""" Edges is a zip of starts and edges. """
		prevEnd = 0
		novelPath = False
		for start, end in edges:
			if prevEnd:
				if not novelPath:
					print(prevEnd, start)
					found, prevEnd2, start2 =  self.labelsFound(prevEnd, start, True)
				if (novelPath or not found) or (found and self.nodes[start2] not in self.nodes[prevEnd2].nextNodes):
					novelPath = True
					start += '-' + self.numPaths
				else:
					prevEnd = prevEnd2
					start = start2
				self.addEdge(prevEnd, start)
				if not novelPath:
					found, start2, end2 = self.labelsFound(start, end, True)
				if (novelPath or not found) or (found and self.nodes[end2] not in self.nodes[start2].nextNodes):
					end += '-' + self.numPaths
					novelPath = True
				else:
					end = end2
				self.addEdge(start, end)
			else:
				self.addEdge(self.source.label, end)
			prevEnd = end
		if novelPath:
			self.numPaths = str(1 +int(self.numPaths))

	def addEdge(self, labelA, labelB, weight=1):
		""" Adds a directed edge from the node with labelA to the node with
		labelB. Creates new nodes if they aren't in the graph already. Returns 
		True if a novel Edge was created. """
		# if labelA in self.nodes:
		A = self.nodes[labelA]
		# else:
		# 	A = Node(labelA)
		# 	self.nodes[labelA] = A
		# 	sys.stderr.write('this part of the code should never be executed\n')
		if labelB in self.nodes:
			B = self.nodes[labelB]
		else:
			B = Node(labelB)
			self.nodes[labelB] = B
			print('new node ' + labelB)
		return A.addNext(B)

	def labelsFound(self, labelA, labelB, outputOne=False):
		""" Set outputOne to true to return   """
		if labelA in self.nodes:
			labelAfound = True
		else:
			try:
				i = [label[:-2] for label in self.nodes].index(labelA)
				labelA = list(self.nodes.keys())[i]
				labelAfound = True
			except:
				labelAfound = False
		if labelB in self.nodes:
			labelBfound = True
		else:
			try:
				i = [label[:-2] for label in self.nodes].index(labelB)
				labelB = list(self.nodes.keys())[i]
				labelBfound = True
			except:
				labelBfound = False
		if outputOne:
			print('both {} {} nodes are found == {}'.format(labelA, labelB, labelAfound and labelBfound))
			return labelAfound and labelBfound, labelA, labelB
		else:
			return labelAfound, labelBfound, labelA, labelB

	def refine(self, wiggle=6):
		""" I think this is not entirely necessary right now. """
		print('not refining')

	def updateWeight(self, nodeA, nodeB):
		""" Finds both nod"""
		hA, hB = nodeA.label.find('-'), nodeB.label.find('-')
		labelA = nodeA.label[:-hA] if hA >= 0 else nodeA.label
		labelB = nodeB.label[:-hB] if hB >= 0 else nodeB.label

		labelAfound, labelBfound, labelA, labelB = self.labelsFound(labelA, labelB)

		if labelAfound:
			A = self.nodes[labelA]
			if labelBfound:
				B = self.nodes[labelB]
			else:
				B = Node(labelB)
				self.nodes[labelB] = B
		else:
			sys.stderr.write('how often does this happen\n')
			A = Node(labelA)
			B = Node(labelB)
			self.nodes[labelA] = A
			self.nodes[labelB] = B
		A.addNext(B)

	def mergeDAG(self, node, novelPath=False, wiggle=6):
		""" Provide the source node of a DAG to merge with this DAG. """
		for nextNode in node.nextNodes:
			print('merging', node.label, nextNode.label)
			self.updateWeight(node, nextNode)
			self.mergeDAG(nextNode)

	def findPaths(self, node, path='', allpaths=[]):
		for nextNode in node.nextNodes:
			allpaths = self.findPaths(nextNode, path + '>' + node.label, allpaths)
		if not node.nextNodes:
			allpaths += [path + '>' + node.label]
			return allpaths
		return allpaths if allpaths else []

	def pathListToPSL(self, pathlist, chrom='chr1'):
		entries = []
		isoform_num = 0
		for path in pathlist:
			path = path.split('>')[1:]
			blockstarts = []
			blocksizes = []
			add = True
			for label in path:
				if '-' in label:
					label = label[:label.find('-')]
				if add:
					blockstarts += [label]
					prevStart = label
					add = False
				else:
					blocksizes += [int(label) - int(prevStart)]
					add = True
			rng = 'Isoform' + str(isoform_num)
			blocksizes = [str(b) for b in blocksizes]
			entries += [ [0]*9 + [rng] + [0]*3 + [chrom] + [0] + [path[0], path[-1], len(blocksizes),
							','.join(blocksizes)+',', 'ignore', ','.join(blockstarts)+','] ]
			isoform_num += 1
		return entries


isoforms = {}
for line in psl_in:  # collect all possible paths
	line = line.rstrip().split('\t')
	chrom, start, end, name = line[13], int(line[15]), int(line[16]), line[9]
	blockstarts = line[20].split(',')[:-1]
	blocksizes = [int(size) for size in line[18].split(',')[:-1]]
	blockends = [str(int(_start) + size) for _start,size in zip(blockstarts, blocksizes)]
	p_left, p_right = name.split('_')[-2], name.split('_')[-1]
	if chrom not in isoforms:
		isoforms[chrom] = {}
	if start not in isoforms[chrom]:
		isoforms[chrom][start] = {}
		isoforms[chrom][start]['start_score'] = 0
	isoforms[chrom][start]['start_score'] += 2 if p_left else 1
	if end not in isoforms[chrom][start]:
		isoforms[chrom][start][end] = {}
		isoforms[chrom][start][end]['end_score'] = 0
		isoforms[chrom][start][end]['dag'] = DAG(start, end)
	isoforms[chrom][start][end]['end_score'] += 2 if p_right else 1
	print(line[9])
	isoforms[chrom][start][end]['dag'].addEdges(zip(blockstarts, blockends))

print(isoforms)
for chrom in isoforms:  # refine paths 
	for start in isoforms[chrom]:
		for end in isoforms[chrom][start]:
			if end == 'start_score':
				continue
			dag = isoforms[chrom][start][end]['dag']
			dag.refine()
			dag.findPaths(dag.source, allpaths=[])

for chrom in isoforms:  # only keep most frequent TSS/TES
	starts = sorted(isoforms[chrom])
	TSS_group = {}
	TES_group = {}
	for start in starts:
		if not TSS_group:
			TSS_group[start] = isoforms[chrom][start]['start_score']
			if start != starts[-1]:
				continue
		elif abs(start - sum(list(TSS_group.keys()))/len(TSS_group)) < 20:
			TSS_group[start] = isoforms[chrom][start]['start_score']
			continue
		bestTSS = sorted(TSS_group.items(), key=lambda tup: tup[1])[-1][0]
		bestTES, bestEndScore = 0, 0
		for ES in isoforms[chrom][bestTSS]:
			if ES != 'start_score' and isoforms[chrom][bestTSS][ES]['end_score'] > bestEndScore:
				bestTES, bestEndScore = ES, isoforms[chrom][bestTSS][ES]['end_score']
		for end in isoforms[chrom][bestTSS].copy():
			if end != bestTES and end != 'start_score':
				allEnds = isoforms[chrom][bestTSS]
				allEnds[end]['dag'].source.label = str(bestTSS)
				allEnds[bestTES]['dag'].mergeDAG(allEnds[end]['dag'].source)
				allEnds.pop(end)
		for tss in TSS_group:
			if tss != bestTSS:
				for end in isoforms[chrom][tss].copy():
					if end == 'start_score':
						continue
					bestDAG = isoforms[chrom][bestTSS][bestTES]['dag']
					otherDAG = isoforms[chrom][tss][end]['dag']
					otherDAG.source.label = str(bestTSS)
					bestDAG.mergeDAG(otherDAG.source)
					isoforms[chrom].pop(tss)
		TSS_group = {}

with open(psl_out, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for chrom in isoforms:
		for start in isoforms[chrom]:
			for end in isoforms[chrom][start]:
				if end == 'start_score':
					continue
				dag = isoforms[chrom][start][end]['dag'] 
				paths = dag.findPaths(dag.source, allpaths=[])
				print(paths)
				entries = dag.pathListToPSL(paths)
				for e in entries:
					writer.writerow(e)