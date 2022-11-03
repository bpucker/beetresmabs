### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.25 ###

__usage__ = """
					python filter_parent_variants.py
					--p1 <FULL_PATH_TO_P1_VCF>
					--p2 <FULL_PATH_TO_P2_VCF>
					--f1 <FULL_PATH_TO_F1_VCF>
					--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import os, sys

# --- end of imports --- #


def load_heterozygous_variants( vcf_file, min_allele_freq, max_allele_freq ):
	"""! @brief load only heterozygous variants from default GATK output VCF """
	
	heterozygous_variants = {}
	coverages = []
	
	with open( vcf_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				try:
					reads = map( float, parts[-1].split(':')[1].split(',') )
					if len( reads ) == 2:
						if parts[6] == "PASS":
							try:
								allel_ratio = reads[0] / ( reads[0] + reads[1] )
								coverages.append( reads[0] + reads[1] )
								if allel_ratio > 0.3 and allel_ratio < 0.7:
									heterozygous_variants.update( { parts[0].replace('_pilon', "") + '_%_' + parts[1].zfill( 9 ):  { 'ref': parts[3], 'alt': parts[4], 'cov': reads[0] + reads[1], 'line': line } } )
							except ZeroDivisionError:
								pass
				except IndexError:
					pass
			line = f.readline()
	average = sum( coverages ) / len( coverages )
	return heterozygous_variants, average


def load_hom_variants( vcf_file ):
	"""! @brief load all variants from VCF file (GATK default expected) """
	
	variants = {}
	
	with open( vcf_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[-1][:3] == "1/1":
					if parts[6] == "PASS":
						variants.update( { parts[0].replace('_pilon', "") + '_%_' + parts[1].zfill( 9 ): { 'ref': parts[3], 'alt': parts[4], 'line': line } } )
			line = f.readline()
	return variants


def load_all_variants( vcf_file ):
	"""! @brief load all variants from VCF file (GATK default expected) """
	
	variants = {}
	
	with open( vcf_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				variants.update( { parts[0].replace('_pilon', "") + '_%_' + parts[1].zfill( 9 ): { 'ref': parts[3], 'alt': parts[4] } } )
			line = f.readline()
	return variants


def write_variants_to_file( filename, variants ):
	"""! @brief write all variants of given dictionary into output file """
	
	with open( filename, "w" ) as out:
		out.write( "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE\n" )
		for key in sorted( variants.keys() ):
			out.write( variants[ key ][ 'line' ] )


def main( arguments ):
	"""! @brief runs everything """
	
	prefix = arguments[ arguments.index( '--out' )+1 ]
	
	vcf_p1 = arguments[ arguments.index( '--p1' )+1 ]
	vcf_p2 = arguments[ arguments.index( '--p2' )+1 ]
	vcf_f1 = arguments[ arguments.index( '--f1' )+1 ]
	
	if prefix[-1] != "/":
		prefix += "/"
	
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	
	output_vcf_p1 = prefix + "P1_hom.vcf"
	output_vcf_p2 = prefix + "P2_hom.vcf"
	output_vcf_f1 = prefix + "F1_het.vcf"
	output_gold = prefix + "gold_standard.vcf"
	
	if '--min_het_freq' in arguments:
		min_allele_freq = float( arguments[ arguments.index( '--min_het_freq' )+1 ] )
	else:
		min_allele_freq = 0.4
	
	if '--max_het_freq' in arguments:
		max_allele_freq = float( arguments[ arguments.index( '--max_het_freq' )+1 ] )
	else:
		max_allele_freq = 0.6
	
	cov_deviation = 0.2
	
	# --- get heterozygous F1 variants --- #
	f1_hetero_variants, average = load_heterozygous_variants( vcf_f1, min_allele_freq, max_allele_freq )
	print "number of heterozygous variants in F1: " + str( len( f1_hetero_variants.keys() ) )
	f1_keys = sorted( f1_hetero_variants.keys() )
	write_variants_to_file( output_vcf_f1, f1_hetero_variants )
		
	# --- compare P1 and P2 variants --- #
	print "loading P1 variants ..."
	p1_variants = load_hom_variants( vcf_p1 )
	all_p1_variants = load_all_variants( vcf_p1 )
	write_variants_to_file( output_vcf_p1, p1_variants )
	
	print "loading P2 variants .... "
	p2_variants = load_hom_variants( vcf_p2 )
	all_p2_variants = load_all_variants( vcf_p2 )
	write_variants_to_file( output_vcf_p2, p2_variants )
	
	print "homozygous P1 variants: " + str( len( p1_variants.keys() ) )
	print "homozygous P2 variants: " + str( len( p2_variants.keys() ) )
		
	print "constructing total variant list ... "
	total_keys = list( set( p1_variants.keys() + p2_variants.keys() ) )
	
	print "checking for surviving variants ... "
	print "number of variants to check: " + str( len( total_keys ) )
	
	remaining_keys = []
	min_cov = average - ( average * cov_deviation )
	max_cov = average + ( average * cov_deviation )
	
	print "coverage borders: " + str( min_cov ) + " - " + str( max_cov )
	for key in total_keys:
		try:
			if min_cov <= f1_hetero_variants[ key ]['cov'] <= max_cov:
				status = True
				#in P1 clean, but not in P2 at all
				try:
					p1_variants[ key ]
					try:
						all_p2_variants[ key ]
						status = False
					except KeyError:
						pass
				
				#in P2 clean, but not in P1 at all
				except KeyError:
					p2_variants[ key ]
					try:
						all_p1_variants[ key ]
						status = False
					except KeyError:
						pass
				if status:
					remaining_keys.append( key )
		except KeyError:
			pass
	
	print "number of gold variants: " + str( len( remaining_keys ) )
	
	# --- writing remaining P1 and P2 variants into output files --- #
	clean_variants = {}
	for key in remaining_keys:
		clean_variants.update( { key: f1_hetero_variants[ key ] } )
	write_variants_to_file( output_gold, clean_variants )


if '--p1' in sys.argv and '--p2' in sys.argv and '--f1' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )

