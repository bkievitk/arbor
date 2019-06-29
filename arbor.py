from Bio import SeqIO
from Bio.Seq import Seq

class Features:
	def __init__(self, length, minDistance, symbolBucket):
		self.length = length
		self.minDistance = minDistance
		self.symbolBucket = symbolBucket
	def toString(self):
		return "length: {}\nminDistance: {}\nsymbolBucket: {}".format(self.length, self.minDistance, self.symbolBucket)

def process():		
	proteins = {}
	files = ["GCF_002014815.1_ASM201481v1_genomic.gbff", "GCF_001239625.1_7068_7_24_genomic.gbff"]
	for file in files:
		for seqIO in SeqIO.parse(file, "genbank"):
			# Extract repeat regions.
			repeatRegions = []
			for feature in seqIO.features:
				if feature.type == "repeat_region":
					repeatRegions.append(feature.location)
						
			# For every protein
			for feature in seqIO.features:
				if feature.type == "CDS":
					# Find the minimum distance to a repeate region.
					minDistance = float('inf')
					for repeat in repeatRegions:
						if feature.location.start > repeat.start:
							minDistance = min(minDistance, feature.location.start - repeat.end)
						else:
							minDistance = min(minDistance, repeat.start - feature.location.end)

					# Get id and translation and add to feature map.
					id = feature.qualifiers.get("protein_id", "");
					translation = feature.qualifiers.get("translation", "")
					if len(id) > 0 and len(translation) > 0:
						length = len(translation[0])
						
						# Calculate percent of sequence that is composed of each given character.
						cs = {}
						for c in translation[0]:
							if c in cs:
								cs[c] += 1
							else:
								cs[c] = 1
						for c in cs:
							cs[c] /= length
						proteins[id[0]] = Features(length, minDistance, cs)						
	
	# Get symbol buckets for target proteins.
	sb1 = proteins["WP_010922251.1"].symbolBucket
	sb2 = proteins["WP_053019794.1"].symbolBucket
	
	c = 0
	# Check each protein.
	for proteinId in proteins:
		protein = proteins[proteinId]

		# Get distance with both target proteins.
		sum1 = distance(sb1, protein.symbolBucket)
		sum2 = distance(sb2, protein.symbolBucket)
		
		# Protein must be long enough to count and we can't cheat and use the distance with itself.
		if protein.length > 1000 and ((sum1 > 0 and sum1 < .2) or (sum2 > 0 and sum2 < .2)): 
			c += 1
		
		# Target proteins are close to repeat regins and are long.
		if protein.minDistance < 5000 and protein.length > 1000:
			print("Matching protein:")
			print(proteins[proteinId].toString())
			print()
	
	# Show the % of proteins that match within the bucket limit of .2 and have > 1000 elements.
	print("Symbol Bucket Matches:")
	print("{}%\n".format(c / len(proteins) * 100))

	# Print out information on our target proteins by name to compare.
	print("Target Proteins:")
	print(proteins["WP_053019794.1"].toString())
	print(proteins["WP_010922251.1"].toString())

def distance(a, b):
	sum = 0
	for s in a:
		sum += abs(a.get(s,0) - b.get(s,0))
	return sum
			
process()				
