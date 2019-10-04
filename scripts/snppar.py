#!/usr/bin/env python
#
# 
# ****Function Description***
# SNPPar: Parallel/homoplasic SNP finder
#
# Author(s) 
#			D. J. Edwards (David.Edwards@monash.edu) 
#			K. E. Holt
#			B. J. Pope
#			S. Duchene
#
# Example command:
'''
snppar -s snps.csv -g genbank.gb -t tree.tre
'''
#
# Last modified - 4/10/2019
# Recent Changes:	changed default reporting to homoplasic, not parallel
#					change of some input commands as a result
#					added user command to log output
# To add (v0.2dev):	missingness report - highest SNP, isolate, overall missingness
#

import os,sys,subprocess,string,re,random,collections,operator,argparse
from operator import itemgetter
from argparse import ArgumentParser
from Bio import SeqIO, AlignIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, generic_dna
from Bio.SeqRecord import SeqRecord
from Bio.Seq import _dna_complement_table as dna_complement_table
from Bio.Data.CodonTable import TranslationError
from Bio.Align import MultipleSeqAlignment
from ete3 import Tree
from datetime import datetime

# Constants declaration
version = 'V0.1dev'
genefeatures = 'CDS'
excludefeatures = 'gene,misc_feature,repeat_region,mobile_element'
nt = ['A','C','G','T']

def parseArguments():
	parser = ArgumentParser(description='\nSNPPar: Parallel/homoplasic SNP Finder '+ version)
	# Inputs
	parser.add_argument('-s', '--snptable', type=str, help='SNP table (i.e. RedDog output)')
	parser.add_argument('-m', '--mfasta', type=str, help='SNPs in MFASTA format')
	parser.add_argument('-l', '--snp_position_list', type=str, help='SNP position list (required for MFASTA input)')
	parser.add_argument('-t', '--tree', type=str, required=True, help='Phylogenetic tree (required)')
	parser.add_argument('-g', '--genbank', type=str, required=True, help='Genbank reference (required)')
	# Optional Inputs
	parser.add_argument('-d', '--directory', type=str, default='', help='Output directory')
	parser.add_argument('-p', '--prefix', type=str, default='', help='Prefix to add to output files')
	# Optional Flags
	parser.add_argument('-P', '--parallel', default=False, action="store_true", help='Flag for reporting of parallel calls')
	parser.add_argument('-S', '--strict', default=False, action="store_true", help='Flag to output strict parallel calls (for testing, sets \'-P\' to True")')
	parser.add_argument('-C', '--convergent', default=False, action="store_true", help='Flag for reporting of convergent calls')
	parser.add_argument('-R', '--revertant', default=False, action="store_true", help='Flag for reporting of revertant calls')
	parser.add_argument('-a', '--no_all_calls', default=False, action="store_true", help='Flag to turn off reporting of all events at each call position (homoplasic reporting)')
	parser.add_argument('-h', '--no_homoplasic', default=False, action="store_true", help='Flag to turn off homoplasic calls output')
	parser.add_argument('-e', '--no_all_events', default=False, action="store_true", help='Flag to turn off reporting of all mutation events')
	parser.add_argument('-c', '--counting', default=False, action="store_true", help='Flag to display counts during SNP testing - warning: slow with large data sets')
	parser.add_argument('-u', '--no_clean_up', default=False, action="store_true", help='Flag to turn off deletion of intermediate files on completion of run')
	parser.add_argument('-f', '--fastml', default=False, action="store_true", help='Flag to use fastML for ASR (default ASR: TreeTime)')
	# Further optional inputs
	parser.add_argument('-x', '--fastml_execute', type=str, default="fastml", help='Command to execute fastML (default command: "fastml" i.e. on PATH)')
	return parser.parse_args()

def executeCommand(command, log):
	# adapted, with permission, from code by Stephen Watts
	result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	# Manually decoding as subprocess.run decoding replaces \r with \n
	result.stdout = result.stdout.decode()
	result.stderr = result.stderr.decode()
	if result.returncode != 0:
		message = 'Failed to run command: ' + command
		logPrint(log, message, "CRITICAL")
		message = 'stdout: ' + result.stdout
		log(log, message, "CRITICAL")
		message = 'stderr: ' + result.stderr
		log(log, message, "CRITICAL")
		sys.exit(1)
	else:
		message = 'stdout: ' + result.stdout
		logPrint(log, message, "INFO")
	return

def getOutputDirectory(output_directory, fastml):
	if output_directory != '':
		if output_directory[-1] != '/':
			output_directory += '/'
		if not os.path.isdir(output_directory):
			os.mkdir(output_directory)
	if fastml and not os.path.isdir(output_directory+"fastml_out/"):
		os.mkdir(output_directory+"fastml_out/")
	return output_directory

def setLog(directory):
	now = datetime.now()
	timestamp = datetime.timestamp(now)
	dt_object = datetime.fromtimestamp(timestamp)
	time = dt_object.strftime("%Y-%m-%d_%H-%M-%S.%f")
	log = directory + time + "_log.txt"
	message = "\nLog started: ("+ log +")"
	print(message)
	message = time + ": INFO : Log started: ("+ log +")\n"
	outputToFile(log, message)
	return log

def logPrint(log, message, message_type):
	print(message)
	while message.startswith('\n'):
		message = message.lstrip('\n')
	now = datetime.now()
	timestamp = datetime.timestamp(now)
	dt_object = datetime.fromtimestamp(timestamp)
	time = dt_object.strftime("%Y-%m-%d_%H:%M:%S.%f")	
	message = time + ": " + message_type + " : " + message
	if not message.endswith('\n'):
		message += '\n'
	appendToFile(log, message)
	return

def log(log, message, message_type):
	while message.startswith('\n'):
		message = message.lstrip('\n')
	now = datetime.now()
	timestamp = datetime.timestamp(now)
	dt_object = datetime.fromtimestamp(timestamp)
	time = dt_object.strftime("%Y-%m-%d_%H:%M:%S.%f")	
	message = time + ": " + message_type + " : " + message
	if not message.endswith('\n'):
		message += '\n'
	appendToFile(log, message)
	return

def echoUserCommand(command,log):
	message = "\nUser command: " + " ".join(command) + '\n'
	logPrint(log, message, 'INFO')
	return

def readInput(input_file_name):
	input_file_handle = open(input_file_name, 'r')
	input_file = input_file_handle.readlines()
	input_file_handle.close()
	return input_file

def outputToFile(output_file_name, output):
	output_file_handle = open(output_file_name, "w")
	output_file_handle.write(output)
	output_file_handle.close()
	return

def appendToFile(output_file_name, output):
	output_file_handle = open(output_file_name, "a")
	output_file_handle.write(output)
	output_file_handle.close()
	return

def isVariable(snp_calls):
	return len(set(snp_calls.upper()).intersection(nt)) > 1

def readSNPTable(infile, log):
	message = "\nReading SNP table from " + infile
	logPrint(log, message, 'INFO')
	snptable = []
	strains = [] # strains from header
	strainlist = []
	ignored =[]
	fields = []
	count = 0
	keep = []
	keep_ingroup = []
	snp_list = []
	lines = readInput(infile)
	for i in range(len(lines)):
		if fields == []:
			fields = lines[i].rstrip().split(',')
		if len(strains)==0:
			strains = fields
			if len(strainlist) == 0:
				for j in range(1,len(strains)):
					strainlist.append(strains[j]) 
					keep.append(j)
					keep_ingroup.append(j)
		else:
			j=0
			snp = ''
			while lines[i][j] != ',':
				snp += lines[i][j]
				j+=1
			# create list of in-group snp calls
			snp_calls_ingroup = ''
			for k in keep_ingroup:
				snp_calls_ingroup += lines[i][(j+2*(k-1)+1)].upper() 
			if isVariable(snp_calls_ingroup):
				# create list of all snp calls
				snp_calls = ''
				if len(keep) == len(keep_ingroup):
					snp_calls = snp_calls_ingroup
				else:
					for k in keep:
						snp_calls += lines[i][(j+2*(k-1)+1)].upper() 
				snptable.append([snp, snp_calls])
				snp_list.append(int(snp))
				count +=1
			else:
				ignored.append(snp)
				logPrint(log, "invariant position removed: " + snp, 'INFO')
	strains.pop(0) # remove SNP column header
	strains_used = []
	for strain in strains:
		if strain in strainlist:
			strains_used.append(strain)
	message = "\nFinished reading " + str(len(snptable) + len(ignored)) + " SNPs in total\n"
	message += "...keeping " + str(len(snptable)) + " variable SNPs and ignoring " + str(len(ignored)) +" SNPs\n"
	message += "that are non-variable among the " +  str(len(strains_used)) + " isolates"			
	logPrint(log, message, 'INFO')
	return(snptable, strains_used, len(snptable), len(strains_used), snp_list)

def readMFASTA(infile_mfasta, infile_snp_positions, log):
	message = "\nReading SNPs from " + infile_mfasta + " and SNP positions from " + infile_snp_positions
	logPrint(log, message, 'INFO')
	snptable = []
	strains = [] # strains from headers in MFASTA
	strainlist = []
	ignored =[]
	keep = []
	keep_ingroup = []
	snp_list = []
	all_calls = []
	snp_positions = []
	snp_position_lines = readInput(infile_snp_positions)
	for line in snp_position_lines:
		snp_positions.append(line.rstrip())
		all_calls.append('')
	mfasta_lines = readInput(infile_mfasta)
	for line in mfasta_lines:
		if line.startswith('>'):
			name = line.lstrip('>')
			name = name.rstrip()
			strains.append(name)
		else:
			for i in range(len(line.rstrip())):
				all_calls[i]+=line[i].upper()

	for i in range(0,len(strains)):
		strainlist.append(strains[i]) 
		keep.append(i)
		keep_ingroup.append(i)
	strains_used = []
	for strain in strains:
		if strain in strainlist:
			strains_used.append(strain)
	for i in range(len(all_calls)):
		# create list of in-group snp calls
		snp_calls_ingroup = ''
		for j in keep_ingroup:
			snp_calls_ingroup += all_calls[i][j]
		if isVariable(snp_calls_ingroup):
			# create list of all snp calls (currently ingroup and all calls are the same)
			snp_calls = ''
			if len(keep) == len(keep_ingroup):
				snp_calls = snp_calls_ingroup
			else:
				for j in keep:
					snp_calls += all_calls[i][j] 
			snptable.append([snp_positions[i], snp_calls])
			snp_list.append(int(snp_positions[i]))
		else:
			ignored.append(snp_positions[i])
	message = "\nFinished reading " + str(len(snptable) + len(ignored)) + " SNPs in total\n"
	message += "...keeping " + str(len(snptable)) + " variable SNPs and ignoring " + str(len(ignored)) +" SNPs\n"
	message += "that are non-variable among the " +  str(len(strains_used)) + " isolates"			
	logPrint(log, message, 'INFO')
	return(snptable, strains_used, len(snptable), len(strains_used), snp_list)

def readTree(new_tree, log):
	message = "\nReading tree from " + new_tree + "..."
	logPrint(log, message, 'INFO')
	isolate_names = []
	node_names = []
	tree = Tree(new_tree,format=1) # tree for parsing mapping results
	edge = 1
	for node in tree.traverse("preorder"):
		if not node.is_leaf():
			node.name = "N%d" %edge
			node_names.append("N%d" %edge)
			edge += 1
		else: 
			isolate_names.append(node.name)
	return tree, isolate_names, node_names

def sameStrains(strains1,strains2):
	if strains1 == [] and strains2 == []:
		return True
	if len(strains1)!=len(strains2):
		return False
	test_score = 0
	for item in strains1:
		if item in strains2:
			test_score += 1
	if test_score!=len(strains1):
		return False
	test_score = 0
	for item in strains2:
		if item in strains1:
			test_score += 1
	if test_score!=len(strains2):
		return False
	return True

def countCalls(snp,total_isolate_count,snptable):
	#call_counts: A, C, G, T, NA
	call_counts = [0,0,0,0,0]
	for j in range(total_isolate_count):
		call = snptable[snp][1][j].lower()
		if (call == "a"):
			call_counts[0] += 1
		elif (call=="c"):
			call_counts[1] += 1
		elif (call=="g"):
			call_counts[2] += 1
		elif (call=="t"):
			call_counts[3] += 1
		else:
			call_counts[4] += 1
	return call_counts

def getStrainsWithCalls(snp, total_isolate_count, snptable, call1, call2, strains):
	snp_set = []
	alt_set = []
	na_set = []
	for j in range(total_isolate_count):
		call = snptable[snp][1][j].lower()
		if call == call1:
			snp_set.append(strains[j])
		elif call == call2:
			alt_set.append(strains[j])
		else:
			na_set.append(strains[j])
	return snp_set, alt_set, na_set

def setSequence(snp, snptable, node_sequences, tree, snp_set, snp_call, ancestor_call, node_names, na_nodes):
	node_to_test = tree.get_common_ancestor(snp_set)
	descendant_nodes = node_to_test.get_descendants()
	descendant_node_names = []
	descendant_node_names.append(node_to_test.name)
	for node in descendant_nodes:
		descendant_node_names.append(node.name)
	node_sequences.append([snptable[snp][0],""])
	for i in range(len(node_names)):
		if node_names[i] in descendant_node_names:
			node_sequences[-1][1] += snp_call.upper()
		else:
			node_sequences[-1][1] += ancestor_call.upper()
	return node_sequences

def testSNPPattern(snp_set,alt_set,na_set,tree,node_names):
	test_tree, removed_nodes = getNANodes(tree, na_set, node_names)
	if not(test_tree.check_monophyly(values=snp_set, target_attr="name",unrooted=False)[0]):
		if not(test_tree.check_monophyly(values=alt_set, target_attr="name",unrooted=False)[0]):
			return True, False, False, removed_nodes
		else:
			return False, False, True, removed_nodes
	else:
		return False, True, False, removed_nodes

def makeSNPPatterns(snp,snp_set,alt_set,na_set,tree,snps_to_map,monophyletic_snps,snp_call,alt_call,node_names, monophyletic_node_sequences, biallelic_homoplasic, snptable):
	is_paraphyletic, snp_set_monophylic, alt_set_monophylic, na_nodes  = testSNPPattern(snp_set,alt_set,na_set,tree,node_names)
	snp_pattern = [[snp_set, alt_set, na_set, is_paraphyletic, snp_set_monophylic, alt_set_monophylic, na_nodes]]
	if is_paraphyletic:
		snps_to_map.append(snp)
		biallelic_homoplasic += 1
	elif snp_set_monophylic:
		monophyletic_snps.append([snp, alt_call.upper(), snp_call.upper(), snp_set])
		monophyletic_node_sequences = setSequence(snp,snptable,monophyletic_node_sequences, tree, snp_set, snp_call, alt_call, node_names, na_nodes)
	elif alt_set_monophylic:
		monophyletic_snps.append([snp, snp_call.upper(), alt_call.upper(), alt_set])
		monophyletic_node_sequences = setSequence(snp,snptable,monophyletic_node_sequences, tree, alt_set, alt_call, snp_call, node_names, na_nodes)
	return snp_pattern, snps_to_map, monophyletic_snps, monophyletic_node_sequences, biallelic_homoplasic

def addToSNPPatterns(snp,snp_pattern,snp_set,alt_set,na_set,tree,snps_to_map,monophyletic_snps,snp_call,alt_call,node_names,monophyletic_node_sequences, biallelic_homoplasic, snptable):
	match = -1
	k = 0
	while k < len(snp_pattern) and not(match==k):
		if (sameStrains(snp_pattern[k][0],snp_set)):
			if (sameStrains(snp_pattern[k][2],na_set)):
				match = k
		k += 1
	if match == -1:
		is_paraphyletic, snp_set_monophylic, alt_set_monophylic, na_nodes = testSNPPattern(snp_set,alt_set,na_set,tree,node_names)
		new_pattern = [snp_set,alt_set,na_set,is_paraphyletic, snp_set_monophylic, alt_set_monophylic, na_nodes]
		snp_pattern.append(new_pattern)
		if is_paraphyletic:
			snps_to_map.append(snp)
			biallelic_homoplasic += 1
		elif snp_set_monophylic:
			monophyletic_snps.append([snp, alt_call.upper(), snp_call.upper(), snp_set])
			monophyletic_node_sequences = setSequence(snp,snptable,monophyletic_node_sequences, tree, snp_set, snp_call, alt_call, node_names,na_nodes)
		elif alt_set_monophylic:
			monophyletic_snps.append([snp, snp_call.upper(), alt_call.upper(), alt_set])
			monophyletic_node_sequences = setSequence(snp,snptable,monophyletic_node_sequences, tree, alt_set, alt_call, snp_call, node_names, na_nodes)
	else:
		if snp_pattern[match][3]:
			snps_to_map.append(snp)
			biallelic_homoplasic += 1
		elif snp_pattern[match][4]:
			monophyletic_snps.append([snp, alt_call.upper(), snp_call.upper(), snp_set])
			monophyletic_node_sequences = setSequence(snp,snptable,monophyletic_node_sequences, tree, snp_set, snp_call, alt_call, node_names, snp_pattern[match][6])
		elif snp_pattern[match][5]:
			monophyletic_snps.append([snp, snp_call.upper(), alt_call.upper(), alt_set])
			monophyletic_node_sequences = setSequence(snp,snptable,monophyletic_node_sequences, tree, alt_set, alt_call, snp_call, node_names, snp_pattern[match][6])
	return snp_pattern, snps_to_map, monophyletic_snps, monophyletic_node_sequences, biallelic_homoplasic

def getNANodes(tree, na_set, node_names):
	removed_nodes = []
	test_tree = tree.copy()
	if na_set:
		for isolate in na_set:
			if test_tree.search_nodes(name=isolate):
				test_tree.search_nodes(name=isolate)[0].detach()
			else:
				print(isolate + " not found")
# ugly reiterrative bit - needs to be changed to check parent of each removed node to see if too should be removed 
# (i.e. when both child nodes are removed then the internal node also needs to be removed)
		removal_finished = False
		while not removal_finished:
			remove = False
			for leaf in test_tree:
				if leaf.name in node_names:
					removed_nodes.append(leaf.name)
					leaf.name = 'remove'
					remove = True
			if remove:
				test_tree.search_nodes(name='remove')[0].detach()
			else:
				removal_finished = True
	return test_tree, removed_nodes

def getSNPsToTest(total_snp_count,total_isolate_count,snptable,strains,tree,node_names,counting,log):
	message = "\nParsing SNPs to find bi-, tri- and quadallelic SNPs...\n"
	message += "Also testing if biallelic SNPs are homoplasic"
	logPrint(log, message, 'INFO')
	snps_to_map = []
	monophyletic_snps = []
	monophyletic_node_sequences = []
	other_snps = []
	snp_pattern = []
	biallelic = 0
	biallelic_homoplasic = 0
	triallelic = 0
	quadallelic = 0
	if counting:
		message = "\nSNP no.\\patterns\\homoplasic"
		print(message)
	for i in range(total_snp_count):
		call_counts = countCalls(i,total_isolate_count,snptable)
		if ((call_counts[0] + call_counts[1] + call_counts[4]) == total_isolate_count):
			if (call_counts[0] > 1 and call_counts[1] > 1):
				biallelic += 1
				if (call_counts[0] >= call_counts[1]):
					snp_set, alt_set, na_set = getStrainsWithCalls(i, total_isolate_count, snptable, "c", "a", strains)
					snp_call = "c"
					alt_call = "a"
				else:
					snp_set, alt_set, na_set = getStrainsWithCalls(i, total_isolate_count, snptable, "a", "c", strains)
					snp_call = "a"
					alt_call = "c"
				if not snp_pattern:
					snp_pattern, snps_to_map, monophyletic_snps, monophyletic_node_sequences, biallelic_homoplasic = makeSNPPatterns(i,snp_set,alt_set,na_set,tree,snps_to_map,monophyletic_snps,snp_call,alt_call,node_names, monophyletic_node_sequences, biallelic_homoplasic, snptable)
				else:
					snp_pattern, snps_to_map, monophyletic_snps, monophyletic_node_sequences, biallelic_homoplasic = addToSNPPatterns(i,snp_pattern,snp_set,alt_set,na_set,tree,snps_to_map,monophyletic_snps,snp_call,alt_call,node_names, monophyletic_node_sequences, biallelic_homoplasic, snptable)
			else:
				if call_counts[0] == 1:
					snp_set, alt_set, na_set = getStrainsWithCalls(i, total_isolate_count, snptable, "a", "c", strains)
					if len(na_set) > 1:
						test_tree, na_nodes = getNANodes(tree,na_set,node_names)
					else:
						na_nodes = []
					monophyletic_snps.append([i, "C", "A", snp_set])
					monophyletic_node_sequences.append([snptable[i][0],""])
					for j in range(len(node_names)):
						monophyletic_node_sequences[-1][1] += "C"
				elif call_counts[1] == 1:
					snp_set, alt_set, na_set = getStrainsWithCalls(i, total_isolate_count, snptable, "c", "a", strains)
					if len(na_set) > 1:
						test_tree, na_nodes = getNANodes(tree,na_set,node_names)
					else:
						na_nodes = []
					monophyletic_snps.append([i, "A", "C", snp_set])
					monophyletic_node_sequences.append([snptable[i][0],""])
					for j in range(len(node_names)):
						monophyletic_node_sequences[-1][1] += "A"
		elif ((call_counts[0] + call_counts[2] + call_counts[4]) == total_isolate_count):
			if (call_counts[0] > 1 and call_counts[2] > 1):
				biallelic += 1
				if (call_counts[0] >= call_counts[2]):
					snp_set, alt_set, na_set = getStrainsWithCalls(i, total_isolate_count, snptable, "g", "a", strains)
					snp_call = "g"
					alt_call = "a"
				else:
					snp_set, alt_set, na_set = getStrainsWithCalls(i, total_isolate_count, snptable, "a", "g", strains)
					snp_call = "a"
					alt_call = "g"
				if not snp_pattern:
					snp_pattern, snps_to_map, monophyletic_snps, monophyletic_node_sequences, biallelic_homoplasic = makeSNPPatterns(i,snp_set,alt_set,na_set,tree,snps_to_map,monophyletic_snps,snp_call,alt_call,node_names, monophyletic_node_sequences, biallelic_homoplasic, snptable)
				else:
					snp_pattern, snps_to_map, monophyletic_snps, monophyletic_node_sequences, biallelic_homoplasic = addToSNPPatterns(i,snp_pattern,snp_set,alt_set,na_set,tree,snps_to_map,monophyletic_snps,snp_call,alt_call,node_names, monophyletic_node_sequences, biallelic_homoplasic, snptable)
			else:
				if call_counts[0] == 1:
					snp_set, alt_set, na_set = getStrainsWithCalls(i, total_isolate_count, snptable, "a", "g", strains)
					if len(na_set) > 1:
						test_tree, na_nodes = getNANodes(tree,na_set,node_names)
					else:
						na_nodes = []
					monophyletic_snps.append([i, "G", "A", snp_set])
					monophyletic_node_sequences.append([snptable[i][0],""])
					for j in range(len(node_names)):
						monophyletic_node_sequences[-1][1] += "G"
				elif call_counts[2] == 1:
					snp_set, alt_set, na_set = getStrainsWithCalls(i, total_isolate_count, snptable, "g", "a", strains)
					if len(na_set) > 1:
						test_tree, na_nodes = getNANodes(tree,na_set,node_names)
					else:
						na_nodes = []	
					monophyletic_snps.append([i, "A", "G", snp_set])
					monophyletic_node_sequences.append([snptable[i][0],""])
					for j in range(len(node_names)):
						monophyletic_node_sequences[-1][1] += "A"
		elif ((call_counts[0] + call_counts[3] + call_counts[4]) == total_isolate_count):
			if (call_counts[0] > 1 and call_counts[3] > 1):
				biallelic += 1
				if (call_counts[0] >= call_counts[3]):
					snp_set, alt_set, na_set = getStrainsWithCalls(i, total_isolate_count, snptable, "t", "a", strains)
					snp_call = "t"
					alt_call = "a"
				else:
					snp_set, alt_set, na_set = getStrainsWithCalls(i, total_isolate_count, snptable, "a", "t", strains)
					snp_call = "a"
					alt_call = "t"
				if not snp_pattern:
					snp_pattern, snps_to_map, monophyletic_snps, monophyletic_node_sequences, biallelic_homoplasic = makeSNPPatterns(i,snp_set,alt_set,na_set,tree,snps_to_map,monophyletic_snps,snp_call,alt_call,node_names, monophyletic_node_sequences, biallelic_homoplasic, snptable)
				else:
					snp_pattern, snps_to_map, monophyletic_snps, monophyletic_node_sequences, biallelic_homoplasic = addToSNPPatterns(i,snp_pattern,snp_set,alt_set,na_set,tree,snps_to_map,monophyletic_snps,snp_call,alt_call,node_names, monophyletic_node_sequences, biallelic_homoplasic, snptable)
			else:
				if call_counts[0] == 1:
					snp_set, alt_set, na_set = getStrainsWithCalls(i, total_isolate_count, snptable, "a", "t", strains)
					if len(na_set) > 1:
						test_tree, na_nodes = getNANodes(tree,na_set,node_names)
					else:
						na_nodes = []
					monophyletic_snps.append([i, "T", "A", snp_set])
					monophyletic_node_sequences.append([snptable[i][0],""])
					for j in range(len(node_names)):
						monophyletic_node_sequences[-1][1] += "T"
				elif call_counts[3] == 1:
					snp_set, alt_set, na_set = getStrainsWithCalls(i, total_isolate_count, snptable, "t", "a", strains)
					if len(na_set) > 1:
						test_tree, na_nodes = getNANodes(tree,na_set,node_names)
					else:
						na_nodes = []
					monophyletic_snps.append([i, "A", "T", snp_set])
					monophyletic_node_sequences.append([snptable[i][0],""])
					for j in range(len(node_names)):
						monophyletic_node_sequences[-1][1] += "A"
		elif ((call_counts[1] + call_counts[2] + call_counts[4]) == total_isolate_count):
			if (call_counts[1] > 1 and call_counts[2] > 1):
				biallelic += 1
				if (call_counts[1] >= call_counts[2]):
					snp_set, alt_set, na_set = getStrainsWithCalls(i, total_isolate_count, snptable, "g", "c", strains)
					snp_call = "g"
					alt_call = "c"
				else:
					snp_set, alt_set, na_set = getStrainsWithCalls(i, total_isolate_count, snptable, "c", "g", strains)
					snp_call = "c"
					alt_call = "g"
				if not snp_pattern:
					snp_pattern, snps_to_map, monophyletic_snps, monophyletic_node_sequences, biallelic_homoplasic = makeSNPPatterns(i,snp_set,alt_set,na_set,tree,snps_to_map,monophyletic_snps,snp_call,alt_call,node_names, monophyletic_node_sequences, biallelic_homoplasic, snptable)
				else:
					snp_pattern, snps_to_map, monophyletic_snps, monophyletic_node_sequences, biallelic_homoplasic = addToSNPPatterns(i,snp_pattern,snp_set,alt_set,na_set,tree,snps_to_map,monophyletic_snps,snp_call,alt_call,node_names, monophyletic_node_sequences, biallelic_homoplasic, snptable)
			else:
				if call_counts[1] == 1:
					snp_set, alt_set, na_set = getStrainsWithCalls(i, total_isolate_count, snptable, "c", "g", strains)
					if len(na_set) > 1:
						test_tree, na_nodes = getNANodes(tree,na_set,node_names)
					else:
						na_nodes = []
					monophyletic_snps.append([i, "G", "C", snp_set])
					monophyletic_node_sequences.append([snptable[i][0],""])
					for j in range(len(node_names)):
						monophyletic_node_sequences[-1][1] += "G"
				elif call_counts[2] == 1:
					snp_set, alt_set, na_set = getStrainsWithCalls(i, total_isolate_count, snptable, "g", "c", strains)
					if len(na_set) > 1:
						test_tree, na_nodes = getNANodes(tree,na_set,node_names)
					else:
						na_nodes = []
					monophyletic_snps.append([i, "C", "G", snp_set])
					monophyletic_node_sequences.append([snptable[i][0],""])
					for j in range(len(node_names)):
						monophyletic_node_sequences[-1][1] += "C"
		elif ((call_counts[1] + call_counts[3] + call_counts[4]) == total_isolate_count):
			if (call_counts[1] > 1 and call_counts[3] > 1):
				biallelic += 1
				if (call_counts[1] >= call_counts[3]):
					snp_set, alt_set, na_set = getStrainsWithCalls(i, total_isolate_count, snptable, "t", "c", strains)
					snp_call = "t"
					alt_call = "c"
				else:
					snp_set, alt_set, na_set = getStrainsWithCalls(i, total_isolate_count, snptable, "c", "t", strains)
					snp_call = "c"
					alt_call = "t"
				if not snp_pattern:
					snp_pattern, snps_to_map, monophyletic_snps, monophyletic_node_sequences, biallelic_homoplasic = makeSNPPatterns(i,snp_set,alt_set,na_set,tree,snps_to_map,monophyletic_snps,snp_call,alt_call,node_names, monophyletic_node_sequences, biallelic_homoplasic, snptable)
				else:
					snp_pattern, snps_to_map, monophyletic_snps, monophyletic_node_sequences, biallelic_homoplasic = addToSNPPatterns(i,snp_pattern,snp_set,alt_set,na_set,tree,snps_to_map,monophyletic_snps,snp_call,alt_call,node_names, monophyletic_node_sequences, biallelic_homoplasic, snptable)
			else:
				if call_counts[1] == 1:
					snp_set, alt_set, na_set = getStrainsWithCalls(i, total_isolate_count, snptable, "c", "t", strains)
					if len(na_set) > 1:
						test_tree, na_nodes = getNANodes(tree,na_set,node_names)
					else:
						na_nodes = []
					monophyletic_snps.append([i, "T", "C", snp_set])
					monophyletic_node_sequences.append([snptable[i][0],""])
					for j in range(len(node_names)):
						monophyletic_node_sequences[-1][1] += "T"
				elif call_counts[3] == 1:
					snp_set, alt_set, na_set = getStrainsWithCalls(i, total_isolate_count, snptable, "t", "c", strains)
					if len(na_set) > 1:
						test_tree, na_nodes = getNANodes(tree,na_set,node_names)
					else:
						na_nodes = []
					monophyletic_snps.append([i, "C", "T", snp_set])
					monophyletic_node_sequences.append([snptable[i][0],""])
					for j in range(len(node_names)):
						monophyletic_node_sequences[-1][1] += "C"
		elif ((call_counts[2] + call_counts[3] + call_counts[4]) == total_isolate_count):
			if (call_counts[2] > 1 and call_counts[3] > 1):
				biallelic += 1
				if (call_counts[2] >= call_counts[3]):
					snp_set, alt_set, na_set = getStrainsWithCalls(i, total_isolate_count, snptable, "t", "g", strains)
					snp_call = "t"
					alt_call = "g"
				else:
					snp_set, alt_set, na_set = getStrainsWithCalls(i, total_isolate_count, snptable, "g", "t", strains)
					snp_call = "g"
					alt_call = "t"
				if not snp_pattern:
					snp_pattern, snps_to_map, monophyletic_snps, monophyletic_node_sequences, biallelic_homoplasic = makeSNPPatterns(i,snp_set,alt_set,na_set,tree,snps_to_map,monophyletic_snps,snp_call,alt_call,node_names, monophyletic_node_sequences, biallelic_homoplasic, snptable)
				else:
					snp_pattern, snps_to_map, monophyletic_snps, monophyletic_node_sequences, biallelic_homoplasic = addToSNPPatterns(i,snp_pattern,snp_set,alt_set,na_set,tree,snps_to_map,monophyletic_snps,snp_call,alt_call,node_names, monophyletic_node_sequences, biallelic_homoplasic, snptable)
			else:
				if call_counts[2] == 1:
					snp_set, alt_set, na_set = getStrainsWithCalls(i, total_isolate_count, snptable, "g", "t", strains)
					if len(na_set) > 1:
						test_tree, na_nodes = getNANodes(tree,na_set,node_names)
					else:
						na_nodes = []
					monophyletic_snps.append([i, "T", "G", snp_set])
					monophyletic_node_sequences.append([snptable[i][0],""])
					for j in range(len(node_names)):
						monophyletic_node_sequences[-1][1] += "T"
				elif call_counts[3] == 1:
					snp_set, alt_set, na_set = getStrainsWithCalls(i, total_isolate_count, snptable, "t", "g", strains)
					if len(na_set) > 1:
						test_tree, na_nodes = getNANodes(tree,na_set,node_names)
					else:
						na_nodes = []
					monophyletic_snps.append([i, "G", "T", snp_set])
					monophyletic_node_sequences.append([snptable[i][0],""])
					for j in range(len(node_names)):
						monophyletic_node_sequences[-1][1] += "G"
		elif ((call_counts[0] + call_counts[1] + call_counts[2] + call_counts[4]) == total_isolate_count):
			if (call_counts[0] > 0 and call_counts[1] > 0 and call_counts[2] > 0):
				other_snps.append(i)
				triallelic += 1
		elif ((call_counts[0] + call_counts[1] + call_counts[3] + call_counts[4]) == total_isolate_count):
			if (call_counts[0] > 0 and call_counts[1] > 0 and call_counts[3] > 0):
				other_snps.append(i)
				triallelic += 1
		elif ((call_counts[0] + call_counts[2] + call_counts[3] + call_counts[4]) == total_isolate_count):
			if (call_counts[0] > 0 and call_counts[2] > 0 and call_counts[3] > 0):
				other_snps.append(i)
				triallelic += 1
		elif ((call_counts[1] + call_counts[2] + call_counts[3] + call_counts[4]) == total_isolate_count):
			if (call_counts[1] > 0 and call_counts[2] > 0 and call_counts[3] > 0):
				other_snps.append(i)
				triallelic += 1
		elif ((call_counts[0] + call_counts[1] + call_counts[2] + call_counts[3] + call_counts[4]) == total_isolate_count):
			if (call_counts[0] > 0 and call_counts[1] > 0 and call_counts[2] > 0 and call_counts[3] > 0):
				other_snps.append(i)
				quadallelic += 1
		if counting:
			print((str(i+1)+'\\'+str(len(snp_pattern))+'\\'+str(biallelic_homoplasic)), end='\r')
	if counting:
		print()
	message = "\nBiallelic SNPs (>1 one isolate): " + str(biallelic)
	message += "\nBiallelic SNP patterns tested: " + str(len(snp_pattern))
	message += "\nParaphyletic SNPs found: " + str(biallelic_homoplasic)
	message += "\nTriallelic SNPs found: " + str(triallelic)
	message += "\nQuadallelic SNPs found: " + str(quadallelic)
	snps_to_map = snps_to_map + other_snps
	snps_to_map.sort()
	message += "\n\nTotal SNPs for mapping: " + str(len(snps_to_map))
	logPrint(log, message, 'INFO')
	return snps_to_map, monophyletic_snps, monophyletic_node_sequences

def readSeqs(output_seqs,strains,snps_to_map,snptable):
	sequences = readInput(output_seqs)
	nodes = []
	node_snptable = []
	for i in range(len(snps_to_map)):
		node_snptable.append([snptable[snps_to_map[i]][0],''])
	i = 0
	while i < len(sequences):
		test_name = sequences[i][1:-1]
		if not(test_name in strains):
			if not(test_name in nodes):
				nodes.append(test_name)
			i += 1
			calls = sequences[i].rstrip()
			for j in range(len(calls)):
				node_snptable[j][1] = ''.join([node_snptable[j][1],calls[j]])
			i += 1
		else:
			i += 2
	return node_snptable, nodes

def readSeqsTT(output_seqs,strains,snps_to_map,snptable):
	sequences = readInput(output_seqs)
	nodes = []
	node_snptable = []
	for i in range(len(snps_to_map)):
		node_snptable.append([snptable[snps_to_map[i]][0],''])
	i = 0
	while i < len(sequences):
		test_name = sequences[i][1:-1]
		if not(test_name in strains):
			test_name = 'N' + str(int(test_name[5:])+1)
			if not(test_name in nodes):
				nodes.append(test_name)
			i += 1
			calls = ''
			while i < len(sequences) and not(sequences[i].startswith('>')):
				calls += sequences[i][:-1]
				i += 1
			for j in range(len(calls)):
				node_snptable[j][1] = ''.join([node_snptable[j][1],calls[j]])
		else:
			i+=1
			while i < len(sequences) and not(sequences[i].startswith('>')):
				i += 1
	return node_snptable, nodes

def getSnpsOnBranch(parent, child, newtable, strains, node_snptable, nodes, snps_mapped,tree, snps_to_map, snptable):
	parent_seq = ''
	parent_index = nodes.index(parent)
	for i in range(len(node_snptable)):
		parent_seq += node_snptable[i][1][parent_index]
	child_seq = ''
	if child in nodes:
		child_index = nodes.index(child)
		for i in range(len(node_snptable)):
			child_seq += node_snptable[i][1][child_index]
	else:
		child_index = strains.index(child)
		for i in range(len(newtable)):
			child_seq += newtable[i][1][child_index]
	for i in range(len(parent_seq)):
		if parent_seq[i].upper() in nt and child_seq[i].upper() in nt:
			if parent_seq[i].upper() != child_seq[i].upper():
				snps_mapped.append([int(snptable[snps_to_map[i]][0]), parent,child, parent_seq[i].upper(), child_seq[i].upper()])
	return snps_mapped

def mapSNPs(snps_to_map, snptable, strains, tree_name, prefix,log):
	aln_file_name = prefix + "snps_to_map.mfasta"
	newtable = []
	for i in snps_to_map:
		newtable.append(snptable[i])
	writeMFASTA(aln_file_name,newtable,strains)
	output_seqs = prefix+"fastml_out/fastml_seq.joint.fasta"
	fastml_tree_name = prefix+"fastml_out/fastml_tree.newick.txt"
	fastml_command = " ".join([fastml_execute,"-t",tree_name,"-s",aln_file_name,
						"-x",fastml_tree_name,
						"-y",prefix+"fastml_out/fastml_tree.ancestor.txt",
						"-j",output_seqs,
						"-k",prefix+"fastml_out/fastml_seq.marginal.fasta",
						"-d",prefix+"fastml_out/fastml_prob.joint.txt",
						"-y",prefix+"fastml_out/fastml_prob.marginal.txt",
						"-mh","-qf","-b"])
	message = "\nRunning fastml: " + fastml_command+"\n"
	logPrint(log, message, 'INFO')
	executeCommand(fastml_command, log)
	fastml_tree = Tree(fastml_tree_name,format=1)
	message = "\nExtracting internal node sequences from ASR results..."
	logPrint(log, message, 'INFO')
	node_snptable, nodes = readSeqs(output_seqs,strains,snps_to_map,snptable)
	snps_mapped = []
	message = "\nExtracting mutation events from ASR results..."
	logPrint(log, message, 'INFO')
	for node in fastml_tree.traverse("preorder"):
		if not node.is_leaf():
			for child in node.children:
				snps_mapped = getSnpsOnBranch(node.name, child.name, newtable, strains, node_snptable, nodes, snps_mapped, fastml_tree, snps_to_map, snptable)
	snps_mapped = sorted(snps_mapped, key=itemgetter(0))
	return snps_mapped, fastml_tree, node_snptable, nodes

def mapSNPsTT(snps_to_map, snptable, strains, tree_name, directory, tree, prefix,log):
	aln_file_name = prefix + "snps_to_map.mfasta"
	newtable = []
	for i in snps_to_map:
		newtable.append(snptable[i])
	writeMFASTAaddNs(aln_file_name,newtable,strains)
	output_dir = directory+"treetime_out/"
	treetime_command = " ".join(["treetime ancestral --aln",aln_file_name,"--tree",tree_name,"--outdir",output_dir,"--report-ambiguous --verbose 2"])
	message = "\nRunning TreeTime: " + treetime_command+"\n"
	logPrint(log, message, 'INFO')
	executeCommand(treetime_command, log)
	message = "\nExtracting mutation events from ASR results..."
	logPrint(log, message, 'INFO')
	snps_mapped = readMappedSNPs(output_dir+'annotated_tree.nexus',tree,snps_to_map,snptable)
	output_seqs = directory+"treetime_out/ancestral_sequences.fasta"
	message = "\nExtracting internal node sequences from ASR results..."
	logPrint(log, message, 'INFO')
	node_snptable, nodes = readSeqsTT(output_seqs,strains,snps_to_map,snptable)
	return snps_mapped, node_snptable, nodes

def readGenbank(genbank,log):
	message = "\nReading Genbank file from " + genbank
	logPrint(log, message, 'INFO')
	handle = open(genbank,"r")
	record = SeqIO.read(handle, "genbank")
	sequence = record.seq
	geneannot = record.features
	return record, sequence, geneannot

def makeGeneIndex(geneannot,sequenceLength,log):
	feature_list = []
	feature_count = 0
	for feature in geneannot:
		if (feature.type != "source" 
			and feature.type not in excludefeatures 
			and feature.type in genefeatures):
			if not str(feature.location).startswith('join'):
				strand = feature.location.strand
				if strand:
					start = feature.location.nofuzzy_start
					stop = feature.location.nofuzzy_end + 1
				else:
					start = feature.location.nofuzzy_start + 1
					stop = feature.location.nofuzzy_end				
				tag = feature
				feature_list.append([start,stop,feature_count,strand])
			else:
				message = "Split gene found in Genbank reference\nThis gene will not be included in the results: " + str(feature.location)
				logPrint(log, message, 'WARNING')
		feature_count += 1
	feature_slice = []
	if len(feature_list) > 0:
		slice_size = sequenceLength//len(feature_list)+1
		for slice in range((sequenceLength//slice_size)+2):
			feature_slice.append([])
	else:
		slice_size = len(record) +1
		feature_slice.append([])
		feature_slice.append([])
	feature_count=0
	for feature in feature_list:
		slice1 = feature_list[feature_count][0]//slice_size
		slice2 = feature_list[feature_count][1]//slice_size
		feature_slice[slice1].append([feature_list[feature_count][0],feature_list[feature_count][1],feature_list[feature_count][2],feature_list[feature_count][3]])
		while slice1 < slice2:
			slice1 += 1
			feature_slice[slice1].append([feature_list[feature_count][0],feature_list[feature_count][1],feature_list[feature_count][2],feature_list[feature_count][3]])
		feature_count += 1
	message = "Found "+str(len(feature_list))+" genes"
	logPrint(log, message, 'INFO')
	return slice_size, feature_slice, feature_list

def findNearestGenes(x, feature_slice, slice_id, geneannot, sequenceLength):
	gene_up = []
	gene_down = []
	if len(feature_slice[slice_id]) == 0:
		offset = 1
	else:
		offset = 0
	while not(gene_up and gene_down):
		if not(offset):
			if len(feature_slice[slice_id]) == 1:
				strand = feature_slice[slice_id][0][3]
				try:
					locus_tag = geneannot[feature_slice[slice_id][0][2]].qualifiers['locus_tag'][0]
				except:
					locus_tag = "Tag_" + str(feature_slice[slice_id][0][1])

				if strand == 1:
					if x < feature_slice[slice_id][0][0]:
						gene_down = [locus_tag, strand, feature_slice[slice_id][0][0] + 1 - x]
					else:
						gene_up = [locus_tag, strand, x + 1 - feature_slice[slice_id][0][1]]
				else:
					if x <= feature_slice[slice_id][0][0]:
						gene_down = [locus_tag, strand, feature_slice[slice_id][0][0] + 1 - x]
					else:
						gene_up = [locus_tag, strand, x + 1 - feature_slice[slice_id][0][1]]
			elif len(feature_slice[slice_id]) > 1:
				start_distances = []
				stop_distances = []
				for i in range(len(feature_slice[slice_id])):

					strand = feature_slice[slice_id][i][3]
					if strand == 1:
						start_distances.append(feature_slice[slice_id][i][0] + 1 - x)
						stop_distances.append(x + 1 - feature_slice[slice_id][i][1])
					else:
						start_distances.append(feature_slice[slice_id][i][0] + 1 - x)
						stop_distances.append(x + 1 - feature_slice[slice_id][i][1])
				start_index = -1
				stop_index = -1
				start_min = -1
				stop_min = -1
				for i in range(len(feature_slice[slice_id])):
					if start_distances[i] > 0:
						if start_min == -1:
							start_min = start_distances[i]
							start_index = i
						else:
							if start_distances[i] < start_min:
								start_min = start_distances[i]
								start_index = i
					if stop_distances[i] > 0:
						if stop_min == -1:
							stop_min = stop_distances[i]
							stop_index = i
						else:
							if stop_distances[i] < stop_min:
								stop_min = stop_distances[i]
								stop_index = i
				if start_min > 0:
					strand = feature_slice[slice_id][start_index][3]
					try:
						locus_tag = geneannot[feature_slice[slice_id][start_index][2]].qualifiers['locus_tag'][0]
					except:
						locus_tag = "Tag_" + str(feature_slice[slice_id][start_index][1])
					gene_down = [locus_tag, strand, start_min]
				if stop_min > 0:
					strand = feature_slice[slice_id][stop_index][3]
					try:
						locus_tag = geneannot[feature_slice[slice_id][stop_index][2]].qualifiers['locus_tag'][0]
					except:
						locus_tag = "Tag_" + str(feature_slice[slice_id][stop_index][1])
					gene_up = [locus_tag, strand, stop_min]
			offset = 1
		else:
			if not(gene_up):
				cross_break = False
				index_offset = slice_id - offset
				if index_offset < 0:
					index_offset = len(feature_slice) + index_offset
					cross_break = True
				if len(feature_slice[index_offset]) == 1:
					strand = feature_slice[index_offset][0][3]
					try:
						locus_tag = geneannot[feature_slice[index_offset][0][2]].qualifiers['locus_tag'][0]
					except:
						locus_tag = "Tag_" + str(feature_slice[index_offset][0][2])
					if strand == 1:
						if not(cross_break):
							gene_up = [locus_tag, strand, x + 1 - feature_slice[index_offset][0][1]]
						else:
							gene_up = [locus_tag, strand, sequenceLength - feature_slice[index_offset][0][1] + x + 1]
					else:
						if not(cross_break):
							gene_up = [locus_tag, strand, x + 1 - feature_slice[index_offset][0][1]]
						else:
							gene_up = [locus_tag, strand, sequenceLength - feature_slice[index_offset][0][1] + x + 1]
				elif len(feature_slice[index_offset]) > 1:
					stop_distances = []
					for i in range(len(feature_slice[index_offset])):
						strand = feature_slice[index_offset][i][3]
						if strand == 1:
							if not(cross_break):
								stop_distances.append(x + 1 - feature_slice[index_offset][i][1])
							else:
								stop_distances.append(sequenceLength - feature_slice[index_offset][0][1] + x + 1)
						else:
							if not(cross_break):
								stop_distances.append(x + 1 - feature_slice[index_offset][i][1])
							else:
								stop_distances.append(sequenceLength - feature_slice[index_offset][0][1] + x + 1)
					stop_index = -1
					stop_min = -1
					for i in range(len(feature_slice[index_offset])):
						if stop_distances[i] > 0:
							if stop_min == -1:
								stop_min = stop_distances[i]
								stop_index = i
							else:
								if stop_distances[i] < stop_min:
									stop_min = stop_distances[i]
									stop_index = i
					if stop_min > 0 :
						strand = feature_slice[index_offset][stop_index][3]
						try:
							locus_tag = geneannot[feature_slice[index_offset][stop_index][2]].qualifiers['locus_tag'][0]
						except:
							locus_tag = "Tag_" + str(feature_slice[index_offset][stop_index][2])
						gene_up = [locus_tag, strand, stop_min]
			if not(gene_down):
				cross_break = False
				index_offset = slice_id + offset
				if index_offset >= len(feature_slice):
					index_offset = index_offset - len(feature_slice)
					cross_break = True
				if len(feature_slice[index_offset]) == 1:
					strand = feature_slice[index_offset][0][3]
					try:
						locus_tag = geneannot[feature_slice[index_offset][0][2]].qualifiers['locus_tag'][0]
					except:
						locus_tag = "Tag_" + str(feature_slice[index_offset][0][2])
					if strand == 1:
						if not(cross_break):
							gene_down = [locus_tag, strand, feature_slice[index_offset][0][0] + 1 - x]
						else:
							gene_down = [locus_tag, strand, (sequenceLength-x+feature_slice[index_offset][0][0] + 1) ]
					else:
						if not(cross_break):
							gene_down = [locus_tag, strand, feature_slice[index_offset][0][0] + 1 - x]
						else:
							gene_down = [locus_tag, strand, (sequenceLength-x+feature_slice[index_offset][0][0]+1)]
				elif len(feature_slice[index_offset]) > 1:
					start_distances = []
					for i in range(len(feature_slice[index_offset])):
						strand = feature_slice[index_offset][i][3]
						if strand == 1:
							if not(cross_break):
								start_distances.append(feature_slice[index_offset][i][0] + 1 - x)
							else:
								start_distances.append(sequenceLength-x+feature_slice[index_offset][i][0] + 1)
						else:
							if not(cross_break):
								start_distances.append(feature_slice[index_offset][i][0] - x + 1)
							else:
								start_distances.append(sequenceLength-x+feature_slice[index_offset][i][0]+1)
					start_index = -1
					start_min = -1
					for i in range(len(feature_slice[index_offset])):
						if start_distances[i] > 0:
							if start_min == -1:
								start_min = start_distances[i]
								start_index = i
							else:
								if start_distances[i] < start_min:
									start_min = start_distances[i]
									start_index = i
					if start_min > 0:
						strand = feature_slice[index_offset][start_index][3]
						try:
							locus_tag = geneannot[feature_slice[index_offset][start_index][2]].qualifiers['locus_tag'][0]
						except:
							locus_tag = "Tag_" + str(feature_slice[index_offset][start_index][2])
						gene_down = [locus_tag, strand, start_min]
			offset +=1
	intergenic_genes = [x, "intergenic"] + gene_up + gene_down
	return intergenic_genes

def assignSNPs(snps_to_test, snptable, slice_size, feature_slice, geneannot, sequenceLength, output, snp_type,log):
	message = "\nAssigning " + snp_type + " to genes"
	logPrint(log, message, 'INFO')
	intergenic = 0
	intragenic = 0
	for i in range(len(snps_to_test)):
		x = int(snptable[snps_to_test[i]][0])
		slice_id = x//slice_size
		if feature_slice[slice_id] == []:
			intergenic_genes = findNearestGenes(x, feature_slice, slice_id, geneannot, sequenceLength)
			output.append(intergenic_genes)
			intergenic += 1
		else:
			intragenic_flag = False
			for i in range(len(feature_slice[slice_id])):
				add = False
				if feature_slice[slice_id][i][3] == 1:
					if x > feature_slice[slice_id][i][0] and x < feature_slice[slice_id][i][1]:
						intragenic_flag = True
						add = True
						y = x - feature_slice[slice_id][i][0] - 1
				else:
					if x > feature_slice[slice_id][i][0] and x < feature_slice[slice_id][i][1]:
						intragenic_flag = True
						add = True
						y = feature_slice[slice_id][i][1] - x -1
				if add:
					try:
						locus_tag = geneannot[feature_slice[slice_id][i][2]].qualifiers['locus_tag'][0]
					except:
						locus_tag = "Tag_" + str(feature_slice[slice_id][i][2])
					strand = feature_slice[slice_id][i][3]
					codon = y//3 + 1
					codon_position = y - (y//3*3) + 1
					output.append([x,"intragenic",locus_tag,strand,int(codon),int(codon_position)]) 
					intragenic += 1
			if not(intragenic_flag):
				intergenic_genes = findNearestGenes(x, feature_slice, slice_id, geneannot, sequenceLength)
				output.append(intergenic_genes)
				intergenic += 1
	message = "\nOf the " + str(len(snps_to_test)) +" "+ snp_type +" for testing, " + str(intergenic) + " are intergenic, and " + str(intragenic) + " are intragenic"
	logPrint(log, message, 'INFO')
	if (intergenic + intragenic) > len(snps_to_test):
		message = str(intergenic + intragenic - len(snps_to_test)) + " SNPs occur in overlapping genes"
		logPrint(log, message, 'INFO')
	return output

def getMonophyleticSNPs(snp_list, snptable, slice_size, feature_slice, geneannot, sequenceLength, input_data, snp_type,log):
	snps = []
	output = []
	for item in snp_list:
		snps.append(item[0])
	mono_output = assignSNPs(snps, snptable, slice_size, feature_slice, geneannot, sequenceLength, input_data, snp_type,log)
	i = 0
	j = 0
	while i < len(snp_list) and j < len(mono_output):
		if int(snptable[snp_list[i][0]][0]) == mono_output[j][0]:
			output.append(mono_output[j] + snp_list[i][1:])
			j += 1
		elif int(snptable[snp_list[i][0]][0]) > mono_output[j][0]:
			j += 1
		elif mono_output[j][0] > int(snptable[snp_list[i][0]][0]):
			i += 1
	return output

def writeMFASTA(output_file,snptable,strains):
	output = ""
	for strain in range(len(strains)): # cycle over strains
		output = output + ">" + strains[strain] + "\n"
		seq = ''
		for snp in range(len(snptable)): # cycle over SNPs
			seq += snptable[snp][1][strain]
		output = output + seq + "\n"
	outputToFile(output_file, output)	
	return

def writeMFASTAaddNs(output_file,snptable,strains):
	output = ""
	for strain in range(len(strains)): # cycle over strains
		output = output + ">" + strains[strain] + "\n"
		seq = ''
		for snp in range(len(snptable)): # cycle over SNPs
			call = snptable[snp][1][strain]
			if call in nt:
				seq +=  call
			else:
				seq += 'N'
		output = output + seq + "\n"
	outputToFile(output_file, output)	
	return

def makeMFASTA(geneannot,gene_list,gene_snps,intergenic,snptable,strains,prefix):
	gene_tags = []
	if not os.path.isdir(prefix+ "mfasta/"):
			os.mkdir(prefix+ "mfasta/")
	for i in range(len(gene_list)):
		if gene_snps[i]:
			try:
				locus_tag = geneannot[gene_list[i]].qualifiers['locus_tag'][0]
			except:
				locus_tag = "Tag_" + str(gene_list[i])
			gene_tags.append(locus_tag)
			output_file_name = prefix+ "mfasta/" + locus_tag + ".mfasta"
			newtable = []
			for j in gene_snps[i]:
				newtable.append(snptable[j])
			writeMFASTA(output_file_name,newtable,strains)
	output_file_name = prefix+ "mfasta/intergenic.mfasta"
	newtable = []
	for j in intergenic:
		newtable.append(snptable[j])
	writeMFASTA(output_file_name,newtable,strains)
	return gene_tags

def combineParallelResults(parallel_output, snps_mapped,log):
	message = "\nProcessing mapped mutation events..."
	logPrint(log, message, 'INFO')
	parallel_mapped = []
	for item in parallel_output:
		SNP = item[0]
		for result in snps_mapped:
			if result[0] == SNP:
				parallel_mapped.append(item + result[1:])
	return parallel_mapped

def convertMonophyleticResults(monophyletic_output,tree,log):
	message = "\nProcessing monophyletic mutation events..."
	logPrint(log, message, 'INFO')
	monophyletic_mapped = []
	for item in monophyletic_output:
		if len(item[-1])==1:
			child = item[-1][0]
		else:
			isolates = item[-1]
			child_node = tree.get_common_ancestor(isolates)
			child = child_node.name
		try:
			parent_node = tree.search_nodes(name = child)[0].up
			parent = parent_node.name
			monophyletic_mapped.append(item[:-3] + [parent,child,item[-3],item[-2]])
		except:
			print(item)
	return monophyletic_mapped

def combineTables(unordered_table, unordered_nodes, ordered_table, ordered_nodes,log):
	message = "\nCombining monophyletic and mapped mutation events"
	message += "\nAlso combining monophyletic and mapped node sequences..."
	logPrint(log, message, 'INFO')
	order = []
	for item in unordered_nodes:
		order.append(ordered_nodes.index(item))
	for i in range(len(unordered_table)):
		new_calls = []
		for j in range(len(order)):
			new_calls.append("")
		for j in range(len(order)):
			new_calls[order[j]] = unordered_table[i][1][j]
		new_call = ""
		for item in new_calls:
			new_call += item
		unordered_table[i][1] = new_call
	combined_snptable = ordered_table + unordered_table
	for i in range(len(combined_snptable)):
		combined_snptable[i][0] = int(combined_snptable[i][0])
	combined_snptable = sorted(combined_snptable, key=itemgetter(0))
	for i in range(len(combined_snptable)):
		combined_snptable[i][0] = str(combined_snptable[i][0])
	return combined_snptable

def isItemInList(item,the_list):
	result = []
	try:
		index = the_list.index(item)
	except:
		index = -1
	if not index == -1:
		result.append(True)
	else:
		result.append(False)
	result.append(index)
	return result

def checkAncestralBaseForward(codon_base,codonseq,codon,snp_list,node_snptable,ancestral_group):
	testSNP = isItemInList(codon[codon_base],snp_list)
	if testSNP[0]:
		codonseq[codon_base] = node_snptable[testSNP[1]][1][int(ancestral_group[1:])-1]
	return codonseq

def checkAncestralBaseReverse(codon_base,codonseq,codon,snp_list,node_snptable,ancestral_group):
	testSNP = isItemInList(codon[codon_base],snp_list)
	if testSNP[0]:
		codonseq[codon_base] = node_snptable[testSNP[1]][1][int(ancestral_group[1:])-1].translate(dna_complement_table)
	return codonseq

def getAncestralBasesForward(codonseq,positionInCodon,codon,snp_list,node_snptable,ancestral_group):
	if positionInCodon == 3:
		codonseq = checkAncestralBaseForward(0,codonseq,codon,snp_list,node_snptable,ancestral_group)
		codonseq = checkAncestralBaseForward(1,codonseq,codon,snp_list,node_snptable,ancestral_group)
	elif positionInCodon == 1:
		codonseq = checkAncestralBaseForward(1,codonseq,codon,snp_list,node_snptable,ancestral_group)
		codonseq = checkAncestralBaseForward(2,codonseq,codon,snp_list,node_snptable,ancestral_group)
	elif positionInCodon == 2:
		codonseq = checkAncestralBaseForward(0,codonseq,codon,snp_list,node_snptable,ancestral_group)
		codonseq = checkAncestralBaseForward(2,codonseq,codon,snp_list,node_snptable,ancestral_group)
	return codonseq

def getAncestralBasesReverse(codonseq,positionInCodon,codon,snp_list,node_snptable,ancestral_group):
	if positionInCodon == 3:
		codonseq = checkAncestralBaseReverse(0,codonseq,codon,snp_list,node_snptable,ancestral_group)
		codonseq = checkAncestralBaseReverse(1,codonseq,codon,snp_list,node_snptable,ancestral_group)
	elif positionInCodon == 1:
		codonseq = checkAncestralBaseReverse(1,codonseq,codon,snp_list,node_snptable,ancestral_group)
		codonseq = checkAncestralBaseReverse(2,codonseq,codon,snp_list,node_snptable,ancestral_group)
	elif positionInCodon == 2:
		codonseq = checkAncestralBaseReverse(0,codonseq,codon,snp_list,node_snptable,ancestral_group)
		codonseq = checkAncestralBaseReverse(2,codonseq,codon,snp_list,node_snptable,ancestral_group)
	return codonseq

def checkDerivedBaseForward(codon_base,codonseq,codon,snp_list,node_snptable,derived_group,strains,snptable):
	testSNP = isItemInList(codon[codon_base],snp_list)
	if testSNP[0]:
		testGroup = isItemInList(derived_group,strains)
		if testGroup[0]:
			codonseq[codon_base] = snptable[testSNP[1]][1][testGroup[1]]
		else:
			codonseq[codon_base] = node_snptable[testSNP[1]][1][int(derived_group[1:])-1]
	return codonseq

def checkDerivedBaseReverse(codon_base,codonseq,codon,snp_list,node_snptable,derived_group,strains,snptable):
	testSNP = isItemInList(codon[codon_base],snp_list)
	if testSNP[0]:
		testGroup = isItemInList(derived_group,strains)
		if testGroup[0]:
			codonseq[codon_base] = snptable[testSNP[1]][1][testGroup[1]].translate(dna_complement_table)
		else:
			codonseq[codon_base] = node_snptable[testSNP[1]][1][int(derived_group[1:])-1].translate(dna_complement_table)
	return codonseq

def getDerivedBasesForward(codonseq,positionInCodon,codon,snp_list,node_snptable,derived_group,strains,snptable):
	if positionInCodon == 3:
		codonseq = checkDerivedBaseForward(0,codonseq,codon,snp_list,node_snptable,derived_group,strains,snptable)
		codonseq = checkDerivedBaseForward(1,codonseq,codon,snp_list,node_snptable,derived_group,strains,snptable)
	elif positionInCodon == 1:
		codonseq = checkDerivedBaseForward(1,codonseq,codon,snp_list,node_snptable,derived_group,strains,snptable)
		codonseq = checkDerivedBaseForward(2,codonseq,codon,snp_list,node_snptable,derived_group,strains,snptable)
	elif positionInCodon == 2:
		codonseq = checkDerivedBaseForward(0,codonseq,codon,snp_list,node_snptable,derived_group,strains,snptable)
		codonseq = checkDerivedBaseForward(2,codonseq,codon,snp_list,node_snptable,derived_group,strains,snptable)
	return codonseq

def getDerivedBasesReverse(codonseq,positionInCodon,codon,snp_list,node_snptable,derived_group,strains,snptable):
	if positionInCodon == 3:
		codonseq = checkDerivedBaseReverse(0,codonseq,codon,snp_list,node_snptable,derived_group,strains,snptable)
		codonseq = checkDerivedBaseReverse(1,codonseq,codon,snp_list,node_snptable,derived_group,strains,snptable)
	elif positionInCodon == 1:
		codonseq = checkDerivedBaseReverse(1,codonseq,codon,snp_list,node_snptable,derived_group,strains,snptable)
		codonseq = checkDerivedBaseReverse(2,codonseq,codon,snp_list,node_snptable,derived_group,strains,snptable)
	elif positionInCodon == 2:
		codonseq = checkDerivedBaseReverse(0,codonseq,codon,snp_list,node_snptable,derived_group,strains,snptable)
		codonseq = checkDerivedBaseReverse(2,codonseq,codon,snp_list,node_snptable,derived_group,strains,snptable)
	return codonseq

def getCodons(positionInCodon,genestrand,snpPosition,derived,ancestral,derived_group,ancestral_group,sequence,snptable,strains,node_snptable,tree_nodes,snp_list,log):
	codon = ()
	# determine coordinates of codon within genome
	if genestrand == 1:
	# note genestart is in -1 offset space, snp is not
		if positionInCodon == 3:
			codon = (snpPosition-2,snpPosition-1,snpPosition)
		elif positionInCodon == 1:
			codon = (snpPosition,snpPosition+1,snpPosition+2)
		elif positionInCodon == 2:
			codon = (snpPosition-1,snpPosition,snpPosition+1)
		else:
			message = "Unrecognised position in codon: " + positionInCodon
			logPrint(log, message, "CRITICAL")
			sys.exit(1)
	elif genestrand == -1:
	# note genestop is not in -1 offset space
		if positionInCodon == 3:
			codon = (snpPosition+2,snpPosition+1,snpPosition)
		elif positionInCodon == 1:
			codon = (snpPosition,snpPosition-1,snpPosition-2)
		elif positionInCodon == 2:
			codon = (snpPosition+1,snpPosition,snpPosition-1)
		else:
			message = "Unrecognised position in codon: " + positionInCodon
			logPrint(log, message, "CRITICAL")
			sys.exit(1)
	else:
		message = "Unrecognised gene strand:" + genestrand
		logPrint(log, message, "CRITICAL")
		sys.exit(1)
	# extract codon sequence from reference genome
	codonseq = [str(sequence[codon[0]-1]), str(sequence[codon[1]-1]), str(sequence[codon[2]-1])]
	if genestrand == -1:
		# complement the reverse strand
		codonseq = [s.translate(dna_complement_table) for s in codonseq]
	# insert ancestral bases
	# check if other 2 bases in codon are in snp_list, and change where appropriate
	if genestrand == 1:
		codonseq[positionInCodon-1] = ancestral # replace snp within codon
		codonseq = getAncestralBasesForward(codonseq,positionInCodon,codon,snp_list,node_snptable,ancestral_group)
	elif genestrand == -1:
		codonseq[positionInCodon-1] = ancestral.translate(dna_complement_table)
		codonseq = getAncestralBasesReverse(codonseq,positionInCodon,codon,snp_list,node_snptable,ancestral_group)
	ancestral_codon = Seq(''.join(codonseq),IUPAC.unambiguous_dna)
	# mutate with current SNP
	# again checking other bases in codon
	if genestrand == 1:
		codonseq[positionInCodon-1] = derived # replace snp within codon
		codonseq = getDerivedBasesForward(codonseq,positionInCodon,codon,snp_list,node_snptable,derived_group,strains,snptable)
	elif genestrand == -1:
		codonseq[positionInCodon-1] = derived.translate(dna_complement_table)
		codonseq = getDerivedBasesReverse(codonseq,positionInCodon,codon,snp_list,node_snptable,derived_group,strains,snptable)
	derived_codon = Seq(''.join(codonseq),IUPAC.unambiguous_dna)

	# Translate codons; codons containing ambigous bases cannot always be translated, 
	# in these cases set amino acid (AA) to None
	try:
		ancestralAA = ancestral_codon.translate()
	except TranslationError:
		ancestralAA = None
	try:
		derivedAA = derived_codon.translate()
	except TranslationError:
		derivedAA = None
	result = [ancestral_codon,derived_codon,ancestralAA,derivedAA]
	return result

def getConsequences(mapped, snptable, strains, node_snptable, nodes, sequence, snp_list,log):
	message = "\nAssigning coding consequenses to the mutation events..."
	logPrint(log, message, 'INFO')
	syn_count = 0
	ns_count = 0
	ambiguous_count = 0
	intergenic_count = 0
	for i in range(len(mapped)):
		if mapped[i][1] == 'intragenic':
			result = getCodons(mapped[i][5],mapped[i][3],mapped[i][0],mapped[i][9],mapped[i][8],mapped[i][7],mapped[i][6],sequence,snptable,strains,node_snptable,nodes,snp_list,log)
			for item in result:
				mapped[i].append(str(item))
			if result[2] and result[2] == result[3]:
 				mapped[i].append("S")
 				syn_count += 1
			elif result[2] and result[3]:
				mapped[i].append("NS")
				ns_count += 1
			else:
				mapped[i].append("ambiguous")
				ambiguous_count += 1
		else:
			intergenic_count += 1
	message =  str(ns_count)+" non-synonymous changes, "+str(syn_count)+" synonymous changes"
	if ambiguous_count:
		message += (", "+str(ambiguous_count)+" ambiguous changes, and "+str(intergenic_count)+" SNPs in non-coding regions")
	else:
		message += (", and "+str(intergenic_count)+" SNP events in non-coding regions")
	logPrint(log, message, 'INFO')
	return mapped

def findEvents(GetStrict, noAllCalls, getParallel, getConvergent, getRevertant, noHomoplasic, noAllEvents, mapped, prefix, tree, isolate_names,log):
	message = "\nParsing all mutation events..."
	logPrint(log, message, 'INFO')
	last_item = -1
	convergent = []
	revertant = []
	parallel = []
	strict_parallel = []
	strict_count = 0
	i = 0
	while i < len(mapped):
		if not(last_item == -1) and mapped[i][0] == mapped[last_item][0]:		# same SNP position
			# find other mutation events with same position (sorted list!)
			same = 1
			while (i+same) < len(mapped) and mapped[i][0] == mapped[i+same][0]:
				same += 1
			# test them pairwise
			for j in range(last_item,last_item+same):
				for k in range(i,last_item+same+1):
					if not(j==k):
						if mapped[j][2] == mapped[k][2]:		# same gene
							if mapped[j][1] == "intergenic":		
								if mapped[j][10] == mapped[k][10] and mapped[j][11] == mapped[k][11]:		#same ancestor base and same base change
									if j not in parallel:
										parallel.append(j)
										if GetStrict:
											strict = [mapped[j][0], mapped[j][8], mapped[j][9], mapped[j][10], mapped[j][11]]
											if strict not in strict_parallel:
												strict_parallel.append(strict)
									if k not in parallel:
										parallel.append(k)
										if GetStrict:
											strict = [mapped[k][0], mapped[k][8], mapped[k][9], mapped[k][10], mapped[k][11]]
											if strict not in strict_parallel:
												strict_parallel.append(strict)
												strict_count += 1
								elif mapped[j][11] == mapped[k][11]:	#same base change and different ancestor base
									if j not in convergent:
										convergent.append(j)
									if k not in convergent:
										convergent.append(k)
								elif (getRevertant or getHomoplastic):
									if mapped[j][10] == mapped[k][11]: #base change and ancestor match in pair
										test_nodes = []
										node = tree.search_nodes(name=mapped[j][9])[0]
										for item in node.get_descendants():
											test_nodes.append(item.name)
										for name in test_nodes:
											if name == mapped[k][8]:
												if j not in revertant:
													revertant.append(j)
												if k not in revertant:
													revertant.append(k)
									elif mapped[j][11] == mapped[k][10]: #opposite base change and ancestor match in pair
										test_nodes = []
										node = tree.search_nodes(name=mapped[k][9])[0]
										for item in node.get_descendants():
											test_nodes.append(item.name)
										for name in test_nodes:
											if name == mapped[j][8]:
												if j not in revertant:
													revertant.append(j)
												if k not in revertant:
													revertant.append(k)
							else:		# i.e. intragenic
								if mapped[j][8] == mapped[k][8] and mapped[j][9] == mapped[k][9]:		#same ancestor base and same base change
									if j not in parallel:
										parallel.append(j)
										if GetStrict:
											strict = [mapped[j][0], mapped[j][6], mapped[j][7], mapped[j][8], mapped[j][9]]
											if strict not in strict_parallel:
												strict_parallel.append(strict)
									if k not in parallel:
										parallel.append(k)
										if GetStrict:
											strict = [mapped[k][0], mapped[k][6], mapped[k][7], mapped[k][8], mapped[k][9]]
											if strict not in strict_parallel:
												strict_parallel.append(strict)
												strict_count += 1
								elif (getConvergent or getHomoplastic) and mapped[j][9] == mapped[k][9]:	#same base change and different ancestor base
									if j not in convergent:
										convergent.append(j)
									if k not in convergent:
										convergent.append(k)
								elif (getRevertant or getHomoplastic):
									if mapped[j][8] == mapped[k][9]: #base change and ancestor match in pair
										test_nodes = []
										node = tree.search_nodes(name=mapped[j][7])[0]
										for item in node.get_descendants():
											test_nodes.append(item.name)
										for name in test_nodes:
											if name == mapped[k][6]:
												if j not in revertant:
													revertant.append(j)
												if k not in revertant:
													revertant.append(k)
									elif mapped[j][9] == mapped[k][8]: #base change and ancestor match in pair
										test_nodes = []
										node = tree.search_nodes(name=mapped[k][7])[0]
										for item in node.get_descendants():
											test_nodes.append(item.name)
										for name in test_nodes:
											if name == mapped[j][6]:
												if j not in revertant:
													revertant.append(j)
												if k not in revertant:
													revertant.append(k)
			# reset i
			i = last_item + same
			last_item = i
			i += 1
		else: # reset
			last_item = i
			i += 1
	if not noHomoplasic:
		homoplasic = []		
		for item in parallel:
			homoplasic.append(item)
		for item in convergent:
			if item not in homoplasic:
				homoplasic.append(item)
		for item in revertant:
			if item not in homoplasic:
				homoplasic.append(item)
		homoplasic = sorted(homoplasic)
		message = "\nFound " + str(len(homoplasic)) + " mutation events that are homoplasic"
		logPrint(log, message, 'INFO')
		if not noAllCalls:
			message = "\nWriting all calls at homoplasic event positions to " + prefix + "homoplasic_events_all_calls.tsv"
			logPrint(log, message, 'INFO')
			outputAllCallsTSV(prefix+"homoplasic_events_all_calls.tsv",mapped,homoplasic,log)
		else:
			message = "\nWriting homoplasic events to " + prefix + "homoplasic_events.tsv"
			logPrint(log, message, 'INFO')
			outputTSV(prefix+"homoplasic_events.tsv",mapped,homoplasic,log)
	if getParallel:
		message = "\nFound " + str(len(parallel)) + " mutation events that are parallel"
		logPrint(log, message, 'INFO')
		if not noAllCalls:
			message = "\nWriting all calls at parallel event positions to " + prefix + "parallel_events_all_calls.tsv"
			logPrint(log, message, 'INFO')
			outputAllCallsTSV(prefix+"parallel_events_all_calls.tsv",mapped,parallel,log)
		else:
			message = "\nWriting parallel events to " + prefix + "parallel_events.tsv"
			logPrint(log, message, 'INFO')
			outputTSV(prefix+"parallel_events.tsv",mapped,parallel,log)
	if getConvergent:
		message = "\nFound " + str(len(convergent)) + " mutation events that are convergent"
		logPrint(log, message, 'INFO')
		if not noAllCalls:
			message = "\nWriting all calls at convergent event positions to " + prefix + "convergent_events_all_calls.tsv"
			logPrint(log, message, 'INFO')
			outputAllCallsTSV(prefix+"convergent_events_all_calls.tsv",mapped,convergent,log)
		else:
			message = "\nWriting convergent events to " + prefix + "convergent_events.tsv"
			logPrint(log, message, 'INFO')
			outputTSV(prefix+"convergent_events.tsv",mapped,convergent,log)
	if getRevertant:
		message = "\nFound " + str(len(revertant)) + " mutation events that are revertant"
		logPrint(log, message, 'INFO')
		if not noAllCalls:
			message = "\nWriting all calls at revertant event positions to " + prefix + "revertant_events_all_calls.tsv"
			logPrint(log, message, 'INFO')
			outputAllCallsTSV(prefix+"revertant_events_all_calls.tsv",mapped,revertant,log)
		else:
			message = "\nWriting revertant events to " + prefix + "revertant_events.tsv"
			logPrint(log, message, 'INFO')
			outputTSV(prefix+"revertant_events.tsv",mapped,revertant,log)
	if GetStrict:
		message = "\nOf the parallel mutation events, " + str(len(strict_parallel)) + " mutation events are parallel when ignoring genes (i.e. 'strict')"
		message += "\n\nWriting strict parallel events to " + prefix + "strict_parallel_events.tsv"
		logPrint(log, message, 'INFO')
		outputList(prefix+"strict_parallel_events.tsv",strict_parallel)
	if not noAllEvents:
		message = "\nWriting all mutation events to " + prefix + "all_mutation_events.tsv"
		logPrint(log, message, 'INFO')
		outputAllEventsTSV(prefix+"all_mutation_events.tsv",mapped,log)
	total_node_count = adjustRootCounts(countEvents(mapped,getUniqueIndex(mapped)),tree)
	if getParallel:
		parallel_node_count = dropRootCounts(countEvents(mapped,parallel),tree)
	else:
		parallel_node_count = []
	if getConvergent:
		convergent_node_count = dropRootCounts(countEvents(mapped,convergent),tree)
	else:
		convergent_node_count = []
	if getRevertant:
		revertant_node_count = dropRootCounts(countEvents(mapped,revertant),tree)
	else:
		revertant_node_count = []
	if not noHomoplasic:
		homoplasic_node_count = dropRootCounts(countEvents(mapped,homoplasic),tree)
	else:
		homoplasic_node_count = []
	message = "\nMapping mutation events to NEXUS tree: " + prefix + "node_labelled_nexus.tre"
	logPrint(log, message, 'INFO')
	tree_out = mapEventsToTree(tree,total_node_count,parallel_node_count,convergent_node_count,revertant_node_count,homoplasic_node_count,isolate_names)
	outputToFile(prefix+'node_labelled_nexus.tre',tree_out)
	message = "\nAlso mapping mutation events to Newick NHX tree: " + prefix + "node_labelled_newick.tre"
	logPrint(log, message, 'INFO')
	tree_out = mapEventsToNHXTree(tree,total_node_count,parallel_node_count,convergent_node_count,revertant_node_count,homoplasic_node_count,isolate_names)
	outputToFile(prefix+'node_labelled_newick.tre',tree_out)
	return 

def mapEventsToTree(tree,total_node_count,parallel_node_count,convergent_node_count,revertant_node_count,homoplasic_node_count,isolate_names):
	full_split = []
	labels = []
	new_tree = tree.write(format=1)
	new_tree = new_tree.split(',')
	for item in new_tree:
		new_entry = item.split(')')
		for item in new_entry:
			full_split.append(item)
	for i in range(len(full_split)):
		name = full_split[i]
		while name.startswith('('):
			name = name.lstrip('(')
		label = name.split(':')[0]
		if not(label == ';'):
			labels.append(label)
		else:
			labels.append('N1')
	taxa_index_counter = 0
	tree_out = '#NEXUS\nBegin taxa;\n\tDimensions ntax=' + str(len(isolate_names)) + ';\n\tTaxlabels\n'
	tree_out_2 = 'Begin trees;\n\tTranslate\n'  
	for i in range(len(full_split)):
		if labels[i] in isolate_names:
			tree_out += '\t\t' + labels[i] + '\n'
			taxa_index_counter += 1
		if full_split[i] == ';':
			full_split[i] = labels[i] + '[&total=0'
			if homoplasic_node_count:
				full_split[i] += ',homoplasic=0'
			if parallel_node_count:
				full_split[i] += ',parallel=0'
			if convergent_node_count:
				full_split[i] += ',convergent=0'
			if revertant_node_count:
				full_split[i] += ',revertant=0'
			full_split[i] += '];'
		else:
			if labels[i] in isolate_names:
				if taxa_index_counter < len(isolate_names):
					tree_out_2 +=  ('\t\t' + str(taxa_index_counter) + ' ' + labels[i] + ',\n')
				else:
					tree_out_2 +=  ('\t\t' + str(taxa_index_counter) + ' ' + labels[i] + '\n')
			if labels[i] in total_node_count[1]:
				new_index = total_node_count[1].index(labels[i])
				entry = full_split[i].split(':')
				if labels[i] in isolate_names:
					j=0
					while entry[0][j] == '(':
						j+=1
					full_split[i] = ''
					if j:
						for k in range(j):
							full_split[i] += '('								
					full_split[i] += str(taxa_index_counter) + '[&total=' + str(total_node_count[2][new_index])
				else:
					full_split[i] = labels[i] + '[&total=' + str(total_node_count[2][new_index]) 
				if homoplasic_node_count:
					if labels[i] in homoplasic_node_count[1]:
						new_index = homoplasic_node_count[1].index(labels[i])
						full_split[i] += (',homoplasic=' + str(homoplasic_node_count[2][new_index]))
					else:
						full_split[i] += ',homoplasic=0'
				if parallel_node_count:
					if labels[i] in parallel_node_count[1]:
						new_index = parallel_node_count[1].index(labels[i])
						full_split[i] += (',parallel=' + str(parallel_node_count[2][new_index]))
					else:
						full_split[i] += ',parallel=0'
				if convergent_node_count:
					if labels[i] in convergent_node_count[1]:
						new_index = convergent_node_count[1].index(labels[i])
						full_split[i] += (',convergent=' + str(convergent_node_count[2][new_index]))
					else:
						full_split[i] += ',convergent=0'
				if revertant_node_count:
					if labels[i] in revertant_node_count[1]:
						new_index = revertant_node_count[1].index(labels[i])
						full_split[i] += (',revertant=' + str(revertant_node_count[2][new_index]))
					else:
						full_split[i] += ',revertant=0'
				full_split[i] += ']:' + entry[1]
			else: 
				entry = full_split[i].split(':')
				j=0
				while entry[0][j] == '(':
					j+=1
				full_split[i] = ''
				if j:
					for k in range(j):
						full_split[i] += '('
				if labels[i] in isolate_names:
					full_split[i] += str(taxa_index_counter) + '[&total=0'
				else:
					full_split[i] += labels[i] + '[&total=0'
				if homoplasic_node_count:
					full_split[i] += ',homoplasic=0'
				if parallel_node_count:
					full_split[i] += ',parallel=0'
				if convergent_node_count:
					full_split[i] += ',convergent=0'
				if revertant_node_count:
					full_split[i] += ',revertant=0'
				full_split[i] += ']:' + entry[1]
	tree_out += ('\t\t;\nEnd;\n\n' + tree_out_2 + '\t\t;\ntree TREE1 = [&R] ' + full_split[0]) 
	for i in range(1,len(full_split)):
		if labels[i] in isolate_names:
			tree_out += ',' + full_split[i]
		else:
			tree_out += ')' + full_split[i]
	tree_out += '\nEnd;\n'
	return tree_out

def mapEventsToNHXTree(tree,total_node_count,parallel_node_count,convergent_node_count,revertant_node_count,homoplasic_node_count,isolate_names):
	full_split = []
	labels = []
	new_tree = tree.write(format=1)
	new_tree = new_tree.split(',')
	for item in new_tree:
		new_entry = item.split(')')
		for item in new_entry:
			full_split.append(item)
	for i in range(len(full_split)):
		name = full_split[i]
		while name.startswith('('):
			name = name.lstrip('(')
		label = name.split(':')[0]
		if not(label == ';'):
			labels.append(name.split(':')[0])
		else:
			labels.append('N1')
	new_full_split = []
	for i in range(len(full_split)):
		if full_split[i] == ';':
			full_split[i] = 'N1[&&NHX:total=0' 
			if parallel_node_count:
				full_split[i] += ':parallel=0'
			if convergent_node_count:
				full_split[i] += ':convergent=0'
			if revertant_node_count:
				full_split[i] += ':revertant=0'
			if homoplasic_node_count:
				full_split[i] += ':homoplasic=0'
			full_split[i]+= '];'
		else:
			if labels[i] in total_node_count[1]:
				new_index = total_node_count[1].index(labels[i])
				entry = full_split[i].split(':')
				full_split[i] = entry[0]+ '[&&NHX:total=' + str(total_node_count[2][new_index])
				if parallel_node_count:
					if labels[i] in parallel_node_count[1]:
						new_index = parallel_node_count[1].index(labels[i])
						full_split[i] += (':parallel=' + str(parallel_node_count[2][new_index]))
					else:
						full_split[i] += ':parallel=0'
				if convergent_node_count:
					if labels[i] in convergent_node_count[1]:
						new_index = convergent_node_count[1].index(labels[i])
						full_split[i] += (':convergent=' + str(convergent_node_count[2][new_index]))
					else:
						full_split[i] += ':convergent=0'
				if revertant_node_count:
					if labels[i] in revertant_node_count[1]:
						new_index = revertant_node_count[1].index(labels[i])
						full_split[i] += (':revertant=' + str(revertant_node_count[2][new_index]))
					else:
						full_split[i] += ':revertant=0'
				if homoplasic_node_count:
					if labels[i] in homoplasic_node_count[1]:
						new_index = homoplasic_node_count[1].index(labels[i])
						full_split[i] += (':homoplasic=' + str(homoplasic_node_count[2][new_index]))
					else:
						full_split[i] += ':homoplasic=0'
				full_split[i] += ']:' + entry[1]
			else:
				entry = full_split[i].split(':')
				full_split[i] = entry[0] + '[&&NHX:total=0'
				if parallel_node_count:
					full_split[i] += ':parallel=0'
				if convergent_node_count:
					full_split[i] += ':convergent=0'
				if revertant_node_count:
					full_split[i] += ':revertant=0'
				if homoplasic_node_count:
					full_split[i] += ':homoplasic=0'
				full_split[i] += ']:' + entry[1]
	tree_out = full_split[0]
	for i in range(1,len(full_split)):
		if labels[i] in isolate_names:
			tree_out += ',' + full_split[i]
		else:
			tree_out += ')' + full_split[i]
	return tree_out

def adjustRootCounts(node_counts,tree):
	root_indexes = []
	root_distances = []
	root_counts = []
	children = []
	for i in range(len(node_counts[0])):
		if node_counts[0][i] == 'N1':
			children.append(node_counts[1][i])
			root_indexes.append(i)
			root_node = tree.get_tree_root()
			root_distances.append(root_node.get_distance(tree.search_nodes(name=node_counts[1][i])[0]))
			root_counts.append(node_counts[2][i])
	if len(children) < 2:
		root_node = tree.get_tree_root()
		for child in root_node.children:
			if child.name not in children:
				children.append(child.name)
				root_distances.append(root_node.get_distance(child))
	total_distance = 0
	for i in root_distances:
		total_distance += i
	total_counts = 0
	for i in root_counts:
		total_counts += i
	if root_indexes:
		new_counts = []
		new_counts.append(round(total_counts*root_distances[0]/total_distance))
		new_counts.append(total_counts - new_counts[0])
		if len(root_indexes) == 2:
			for i in range(len(root_indexes)):
				node_counts[2][root_indexes[i]] = new_counts[i]
		else:
			node_counts[2][root_indexes[0]] = new_counts[0]
			if new_counts[1] > 0:
				node_counts[0].append('N1')
				node_counts[1].append(children[1])
				node_counts[2].append(new_counts[1])
	return node_counts

def dropRootCounts(node_counts,tree):
	new_node_counts = [[],[],[]]
	for i in range(len(node_counts[0])):
		if not(node_counts[0][i] == 'N1'):
			new_node_counts[0].append(node_counts[0][i])
			new_node_counts[1].append(node_counts[1][i])
			new_node_counts[2].append(node_counts[2][i])
	return new_node_counts

def getUniqueIndex(mapped):
	list_index = []
	for i in range(len(mapped)):
		list_index.append(i)
	last_item = -1
	i = 0
	last_item_tested = False
	while i < len(mapped):
		if not(last_item == -1) and mapped[i][0] == mapped[last_item][0]:		# same SNP position
			# find other mutation events with same position (sorted list!)
			same = 1
			while (i+same) < len(mapped) and mapped[i][0] == mapped[i+same][0]:
				same += 1
			# test them pairwise
			for j in range(last_item,last_item+same):
				for k in range(j+1,last_item+same+1):
					if mapped[j][2] == mapped[k][2]:
						if j not in list_index:
							list_index.append(j)
						if k not in list_index:
							list_index.append(j)
			# reset i
			i = last_item + same
			last_item = i
			i += 1
			last_item_tested = True
		else:
			if not last_item_tested and last_item not in list_index:
				list_index.append(last_item)
			last_item = i
			i += 1
			last_item_tested = False
	if not last_item_tested and last_item not in list_index:
		list_index.append(last_item)
	return list_index

def countEvents(mapped, list_index):
	node_counts = [[],[],[]]
	for i in list_index:
		done = False
		j = 0
		while j < len(node_counts[0]) and not done:
			if mapped[i][1] == 'intergenic':
				if mapped[i][8] == node_counts[0][j]:
					if mapped[i][9] == node_counts[1][j]:
						node_counts[2][j] += 1
						done = True
			else:
				if mapped[i][6] == node_counts[0][j]:
					if mapped[i][7] == node_counts[1][j]:
						node_counts[2][j] += 1
						done = True
			j+=1
		if not done:
			if mapped[i][1] == 'intergenic':
				node_counts[0].append(mapped[i][8])
				node_counts[1].append(mapped[i][9])
				node_counts[2].append(1)
			else:
				node_counts[0].append(mapped[i][6])
				node_counts[1].append(mapped[i][7])
				node_counts[2].append(1)
	return node_counts

def outputTSV(fileName,mapping,list_index,log):
	if list_index:
		output = 'Position\tType\tAncestor_Node\tDerived_Node\tAncestor_Call\tDerived_Call\tGene\tStrand\tCodon\tCodon_Position\tAncestor_Codon\tDerived_Codon\tAncestor_A.A.\tDerived_A.A.\tChange\tUp_Gene\tUp_Gene_Strand\tUp_Gene_Distance\tDown_Gene\tDown_Gene_Strand\tDown_Gene_Distance\n'
		for item in list_index:
			if mapping[item][1] == "intergenic":
				output += str(mapping[item][0]) + '\t' + str(mapping[item][1]) + '\t' + str(mapping[item][8]) + '\t' + str(mapping[item][9]) + '\t' + str(mapping[item][10]) + '\t' + str(mapping[item][11]) + '\t'
				output += '-\t-\t-\t-\t-\t-\t-\t-\t-\t' 
				output += str(mapping[item][2]) + '\t' + str(mapping[item][3]) + '\t' + str(mapping[item][4]) + '\t' + str(mapping[item][5]) + '\t' + str(mapping[item][6]) + '\t' + str(mapping[item][7]) + '\n'
			else:
				output += str(mapping[item][0]) + '\t' + str(mapping[item][1]) + '\t' + str(mapping[item][6]) + '\t' + str(mapping[item][7]) + '\t' + str(mapping[item][8]) + '\t' + str(mapping[item][9]) + '\t'
				output += str(mapping[item][2]) + '\t' + str(mapping[item][3]) + '\t' + str(mapping[item][4]) + '\t' + str(mapping[item][5]) + '\t' + str(mapping[item][10]) + '\t' + str(mapping[item][11]) + '\t' + str(mapping[item][12]) + '\t' + str(mapping[item][13]) + '\t' + str(mapping[item][14]) + '\t'
				output += '-\t-\t-\t-\t-\t-\n'
		outputToFile(fileName,output)
	else:
		message = "No output available... writing of file aborted!"
		logPrint(log, message, 'INFO')
		print(message)
	return

def outputAllCallsTSV(fileName,mapping,list_index,log):
	new_list_index = []
	for item in list_index:
		for i in range(len(mapping)):
			if (mapping[i][0] == mapping[item][0]) and not (mapping[item] == mapping[i]):
				if i not in new_list_index:
					new_list_index.append(i)
	new_list_index = sorted(new_list_index)
	outputTSV(fileName,mapping,new_list_index,log)
	return

def outputAllEventsTSV(fileName,mapping,log):
	list_index = []
	for i in range(len(mapping)):
		list_index.append(i)
	outputTSV(fileName,mapping,list_index,log)
	return

def outputList(output_file_name,a_list):
	output_file_handle = open(output_file_name,"w")
	for item in a_list:
		output = ''
		for i in range(len(item)):
			output += str(item[i])
			if i < len(item)-1:
				output += '\t'
		output += '\n' 
		output_file_handle.write(output)
	output_file_handle.close()
	return

def readMappedSNPs(tt_nexus_tree,the_tree,snps_to_map,snptable):
	output = []
	tree_file = readInput(tt_nexus_tree)
	items = tree_file[6][12:-1].split('],')
	further_items = []
	for item in items:
		more_items = item.split(')')
		for item in more_items:
			further_items.append(item)
	even_more_items = []
	for item in further_items:
		even_more = item.split(',(')
		for item in even_more:
			even_more_items.append(item)
	for item in even_more_items:
		if not item.find("mutations") == -1:
			parts = item.split('[&mutations="')
			if parts[-1][-1] == "]":
				parts[-1] = parts[-1].rstrip(']')
			if parts[-1][-1] == '"':
				parts[-1] = parts[-1].rstrip('"')
			mutations = parts[-1].split(',')
			if not parts[0].find(',') == -1:
				parts[0] = parts[0].split(',')[1]
			parts[0] = parts[0].split(':')[0]
			while parts[0].startswith('('):
				parts[0] = parts[0].lstrip('(')
			if parts[0].startswith('NODE_'):
				parts[0] = 'N' + str(int(parts[0][5:])+1)
			derived_node = the_tree.search_nodes(name=parts[0])[0]
			derived = derived_node.name
			ancestor_node = derived_node.up
			ancestor = ancestor_node.name
			for mutation in mutations:
				snp_index = int(mutation[1:-1])-1
				ancestor_call = mutation[0]
				derived_call = mutation[-1]
				output.append([int(snptable[snps_to_map[snp_index]][0]),ancestor,derived,ancestor_call,derived_call])
	output = sorted(output, key=itemgetter(0))
	return output

def cleanUp(directory, prefix, fastml,log):
	#clean up temp files from fastml step incl. snps_to_map.mfasta, prob.marginal.txt and log.txt
	message = "\nCleaning up intermediate files..."
	logPrint(log, message, 'INFO')
	os.system('rm ' + prefix + 'snps_to_map.mfasta')
	if fastml:
		os.system('rm log.txt')
		os.system('rm prob.marginal.txt')
		os.system('rm -rf ' + directory + 'fastml_out')
	else:
		os.system('rm -rf ' + directory + 'treetime_out')
	return

def printEnd(log):
	message_int = random.randint(1,6)
	message = "\n...Finished\n"
	if message_int == 3:
		message = "\n...Terminado\n"
	elif message_int == 4:
		message = "\n...Fertig\n"
	elif message_int == 5:
		message = "\n...Selesai\n"
	elif message_int == 6:
		message = "\n...Terminé\n"
	logPrint(log, message, 'INFO')
	return

def main():
	# Parse arguments
	arguments = parseArguments()
	directory = getOutputDirectory(arguments.directory, arguments.fastml)
	log = setLog(directory)
	if arguments.strict:
		arguments.parallel = True
	message = '\nSNPPar:\tParallel SNP Finder '+ version
	if arguments.fastml:
		ASR = 'FastML'
	else:
		ASR = 'TreeTime'
	message += '\n\t(utilising ' + ASR + " for ASR)"
	logPrint(log, message, 'INFO')
	echoUserCommand(sys.argv,log)
	if arguments.prefix:
		prefix = directory + arguments.prefix
	else:
		prefix = directory
	# read in SNP table
	if arguments.snptable:
		snptable, strains, total_snp_count, total_isolate_count,snpPositionList = readSNPTable(arguments.snptable, log)
	elif arguments.mfasta and arguments.snp_position_list:
		snptable, strains, total_snp_count, total_isolate_count,snpPositionList = readMFASTA(arguments.mfasta,arguments.snp_position_list, log)
	elif arguments.mfasta:
		message = "No SNP positions list included... terminating run."
		logPrint(log, message, 'CRITICAL')
		sys.exit(1)
	elif arguments.snp_position_list:
		message = "No MFASTA included... terminating run."
		logPrint(log, message, 'CRITICAL')
		sys.exit(1)
	else:
		message = "No SNPs provided... terminating run."
		logPrint(log, message, 'CRITICAL')
		sys.exit(1)
	# read in genbank file
	record, sequence, geneannot = readGenbank(arguments.genbank, log)
	# index gene features
	slice_size, feature_slice, feature_list = makeGeneIndex(geneannot,len(sequence),log)
	# read in tree
	tree, tree_strains, tree_nodes = readTree(arguments.tree,log)
	#check same isolates are found tree and alignment
	if not sameStrains(strains,tree_strains):
		message = "\nEach isolate in the tree should also be found in the SNP table\nand visa-versa - your data does not match this requirement!"
		logPrint(log, message, 'CRITICAL')
		sys.exit(1)
	else:
		message = "Tree and SNP table have same isolates"
		logPrint(log, message, 'INFO')
	# find biallelic SNPs that occur in more than one isolate
	# and test whether their SNP pattern is paraphyletic
	snps_to_map, monophyletic_snps, monophyletic_node_sequences = getSNPsToTest(total_snp_count,total_isolate_count,snptable,strains,tree,tree_nodes,arguments.counting,log)
	# check if snps to be mapped are intra- or intergenic
	monophyletic_output = getMonophyleticSNPs(monophyletic_snps, snptable, slice_size, feature_slice, geneannot, len(sequence), [], "monophyletic SNPs",log)
	parallel_output = assignSNPs(snps_to_map, snptable, slice_size, feature_slice, geneannot, len(sequence), [], "SNPs to map",log)
	# map paraphyletic biallelic [and other (tri- and quadallelic) snps] to tree to get ancestral genotypes
	if arguments.fastml:
		snps_mapped, tree_with_nodes, mapped_node_sequences, node_names_mapped = mapSNPs(snps_to_map, snptable, strains, arguments.tree, directory,log)
	else:
		snps_mapped, mapped_node_sequences, node_names_mapped = mapSNPsTT(snps_to_map,snptable,strains,arguments.tree,directory,tree,prefix,log)
	# process mapped SNPs to get mutation events
	parallel_mapped = combineParallelResults(parallel_output,snps_mapped,log)
	monophyletic_mapped = convertMonophyleticResults(monophyletic_output,tree,log)
	mapped = sorted(parallel_mapped + monophyletic_mapped, key=itemgetter(0))
	node_snptable = combineTables(mapped_node_sequences, node_names_mapped, monophyletic_node_sequences, tree_nodes,log)
	message = "\nWriting node sequences to " + prefix + "node_sequences.fasta"
	logPrint(log, message, 'INFO')
	writeMFASTA(prefix+"node_sequences.fasta",node_snptable,tree_nodes)
	mapped = getConsequences(mapped, snptable, strains, node_snptable, tree_nodes, sequence, snpPositionList,log)
	findEvents(arguments.strict, arguments.no_all_calls, arguments.parallel, arguments.convergent, arguments.revertant, arguments.no_homoplasic, arguments.no_all_events, mapped, prefix, tree, tree_strains,log)
	if not arguments.no_clean_up:
		cleanUp(directory,prefix,arguments.fastml,log)
	printEnd(log)
	return

if __name__ == '__main__':
	main()
