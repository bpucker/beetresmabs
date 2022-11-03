### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python filter_VCF_by_goldstandard.py
					--in <INPUT_VCF>
					--gold <GOLDSTANDARD_VCF>
					--out <OUTPUT_VCF>
					"""

import sys, os

# --- end of imports --- #

def load_vcf( vcf ):
	"""! @brief load variants from VCF file """
	
	variants = {}
	with open( vcf, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if not ',' in parts[4]:
					if parts[6] == "PASS":
						variants.update( { parts[0] + "_%_" + parts[1].zfill(8): None } )
			line = f.readline()
	return variants


def main( arguments ):
	"""! @brief run everything """
	
	in_file = arguments[ arguments.index('--in')+1 ]
	gold_file = arguments[ arguments.index('--gold')+1 ]
	out_file = arguments[ arguments.index('--out')+1 ]
	
	gold = load_vcf( gold_file )
	
	pass_counter = 0
	fail_counter = 0
	
	with open( out_file, "w" ) as out:
		with open( in_file, "r" ) as f:
			line = f.readline()
			while line:
				if line[0] != '#':
					parts = line.strip().split('\t')
					try:
						gold[ parts[0] + "_%_" + parts[1].zfill(8) ]
						out.write( line )
						pass_counter += 1
					except KeyError:
						fail_counter += 1
				else:
					out.write( line )
				line = f.readline()
	
	print "PASS variants: " + str( pass_counter )
	print "FAIL variants: " + str( fail_counter )


if '--in' in sys.argv and '--out' in sys.argv and '--gold' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
