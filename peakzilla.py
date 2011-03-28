#!/usr/bin/python

# Copyright (c) Jonas Steinmann, 2011
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License version 3 as 
# published by the Free Software Foundation.

try:
	from numpy import median, convolve, ones, std
except ImportError:
	sys.stderr.write("Failed to import from numpy, please install numpy!\n")
	sys.exit(1)

import sys
import csv
from operator import add, mul
from time import strftime, localtime
from collections import deque
from array import array
from optparse import OptionParser
from math import fabs

def main():
	# option parser
	usage = 'python peakzilla.py [OPTIONS] chip.bed control.bed > results.tsv'
	parser = OptionParser(usage=usage)
	
	parser.add_option("-s", "--fragment_size",\
	type = "int", dest="fragment_size", default="200",\
	help = "upper limit of fragment size in bp: default=200")
	
	parser.add_option("-m", "--model_threshold",\
	type = "float", dest="model_threshold", default="120",\
	help = "fold enrichment threshold for modeling: default=120")
	
	parser.add_option("-t", "--peak_threshold",\
	type = "float", dest="peak_threshold", default="40",\
	help = "fold enrichment threshold for peak finding: default=40")
	
	parser.add_option("-f", "--fdr",\
	type = "float", dest="fdr", default='1',\
	help = "cutoff for the estimated FDR value: default = 1")
	
	parser.add_option("-q", "--quiet",\
	action = "store_false", dest="verbose", default=True,\
	help = "don't print status messages")
	
	# read arguments and options
	(options, args) = parser.parse_args()
	if len(args) != 2:
		# return help message if argment number is incorrect
		parser.print_help()
		sys.exit(0)
	ip = args[0]
	control = args[1]
	
	# load tags
	print_status('Loading tags ...', options.verbose)
	ip_tags = TagContainer()
	control_tags = TagContainer()
	ip_tags(ip)
	control_tags(control)

	# first attempt of modeling peak size
	print_status('Modeling peak size and shift ...', options.verbose)
	peak_model = PeakShiftModel(ip_tags, options.fragment_size, options.model_threshold)
	
	# change model threshold until it yields a reasonable number of peaks
	while peak_model.peaks_incorporated < 200 or peak_model.peaks_incorporated > 500:
		if peak_model.peaks_incorporated < 500:
			options.model_threshold = options.model_threshold / 2
			print_status('Model threshold was set too high, trying: %.1f'  % options.model_threshold, options.verbose)
			peak_model = PeakShiftModel(ip_tags, options.fragment_size, options.model_threshold)
		else:
			options.model_threshold = options.model_threshold * 1.25
			print_status('Model threshold was set too low, trying: %.1f' % options.model_threshold, options.verbose)
			peak_model = PeakShiftModel(ip_tags, options.fragment_size, options.model_threshold)
	print_status('Used best %d peaks for modeling ...' % peak_model.peaks_incorporated, options.verbose)
	print_status('Peak size is %d bp' % peak_model.peak_size, options.verbose)
	print_status('Peak size SD is %d bp' % peak_model.peak_size_std, options.verbose)
			
	# first attempt to find candidate peaks in control sample
	print_status('Finding potential false positives ...', options.verbose)
	control_peaks = PeakContainer(control_tags, ip_tags, peak_model.peak_size, options.peak_threshold, peak_model.plus_model, peak_model.minus_model)
	
	# change peak threshold until it yields a reasonable number of peaks
	while control_peaks.peak_count < 2000 or control_peaks.peak_count > 10000:
		if control_peaks.peak_count < 2000:
			options.peak_threshold = options.peak_threshold / 2
			print_status('Peak threshold was set too high, trying: %.2f'  % options.peak_threshold, options.verbose)
			control_peaks = PeakContainer(control_tags, ip_tags, peak_model.peak_size, options.peak_threshold, peak_model.plus_model, peak_model.minus_model)
		else:
			options.peak_threshold = options.peak_threshold * 1.25
			print_status('Peak threshold was set too low, trying: %.2f' % options.peak_threshold, options.verbose)
			control_peaks = PeakContainer(control_tags, ip_tags, peak_model.peak_size, options.peak_threshold, peak_model.plus_model, peak_model.minus_model)
	print_status('%d potential false positives found' % control_peaks.peak_count, options.verbose)
	
	# find candidate peaks in IP sample
	print_status('Finding peak candidates ...', options.verbose)
	ip_peaks = PeakContainer(ip_tags, control_tags, peak_model.peak_size, options.peak_threshold, peak_model.plus_model, peak_model.minus_model)
	print_status('%d candidate peaks found' % ip_peaks.peak_count, options.verbose)
	
	# model tag distribution 
	print_status('Modeling tag distribution ...', options.verbose)
	plus_model = ip_peaks.model_tag_distribution()[0]
	minus_model = ip_peaks.model_tag_distribution()[1]
	

	# output model in an R friendly fashion
	print 'plus = c',
	print tuple(plus_model)
	print 'minus = c',
	print tuple(minus_model)

	# find peaks using emirical model
	print_status('Finding peaks with empirical peak model ...', options.verbose)
	control_peaks = PeakContainer(control_tags, ip_tags, peak_model.peak_size, options.peak_threshold, plus_model, minus_model)
	print_status('%d potential false positives found' % control_peaks.peak_count, options.verbose)
	ip_peaks = PeakContainer(ip_tags, control_tags, peak_model.peak_size, options.peak_threshold, plus_model, minus_model)
	print_status('%d candidate peaks found' % ip_peaks.peak_count, options.verbose)
	

	
	print_status('Calculating tag distribution scores ...', options.verbose)
	ip_peaks.determine_distribution_scores(plus_model, minus_model)
	control_peaks.determine_distribution_scores(plus_model, minus_model)
	 
	# calculate FDR
	print_status('Calculating FDR ...', options.verbose)
	ip_peaks.calculate_fdr(control_peaks.peaks)
	
	# write output as bed files
	print_status('Writing results to file ...', options.verbose)
	ip_peaks.write_to_stdout(options)
	print_status('Done!', options.verbose)

class TagContainer:
	# class for loading, storing and manipulating sequence tags
	def __init__(self):
		# intitialize an empty object
		self.tags = {}
		self.tag_number = 0
		self.sorted = False

	def __call__(self, bed_file):
		# when called like a function load bed file and return self
		self.load_bed(bed_file)
		self.sort_tags()
		return self
	
	def add_tag(self, chrom, strand, fiveprime):
		# add tag to dictionary
		if not chrom in self.tags:
			self.tags[chrom] = {}
			# store tags in an array of unsigned integers (4 bytes)
			self.tags[chrom]['+'] = array('i',[])
			self.tags[chrom]['-'] = array('i',[])
			self.tags[chrom][strand].append(fiveprime)
		else:
			self.tags[chrom][strand].append(fiveprime)
		# keep track of total number of tags added
		self.tag_number += 1
		
	def load_bed(self, bed_file):
		# parse a bed file and add contents to self
		for i in csv.reader(open(bed_file), delimiter='\t'):
			chrom = i[0]
			start = int(i[1])
			end = int(i[2])
			strand = i[5]
			# determine five prime end
			if strand == '+':
				fiveprime = start
			else:
				fiveprime = end
			# add tag to container
			self.add_tag(chrom, strand, fiveprime)
			
	def sort_tags(self):
		# sort all tags while preserving the array
		for chrom in self.tags.keys():
			# as sorted returns conversion back to array is required
			self.tags[chrom]['+'] = array('i', sorted(self.tags[chrom]['+']))
			self.tags[chrom]['-'] = array('i', sorted(self.tags[chrom]['-']))
		# change sorted flag to true
		self.sorted = True
	
	def get_chrom_size(self, chrom):
		# chromosome size to consider for scanning of both strands
		if not self.sorted:
			self.sort_tags()
		if self.tags[chrom]['+'] and self.tags[chrom]['-']:
			chrom_size = self.tags[chrom]['-'][-1]
			return chrom_size
		else:
			return 0

	def genome_size(self):
		# genome size to consider for scanning of both strands
		genome_size = 0
		for chrom in self.tags.keys():
			genome_size += self.get_chrom_size(chrom)
		return genome_size
			
	def get_tags(self, chrom, strand):
		# return the whole array of tags
		if chrom in self.tags:
			return self.tags[chrom][strand]
		else:
			return []
			
	def get_tag_number(self, chrom, strand):
		# find out how many tags are mapped to a particular comsomome and strand
		return len(self.tags[chrom][strand])
		
	def get_chrom_names(self):
		# retreive a sorted list of all chromosome names
		return self.tags.keys()

class PeakShiftModel:
	# class for modeling peak size and strand shift
	def __init__(self, tags, fragment_size, fold_threshold):
		self.tags = tags
		self.window_size = fragment_size / 2
		self.fold_threshold = fold_threshold
		self.tag_threshold = tags.tag_number / float(tags.genome_size()) * fragment_size / 2 * fold_threshold
		self.peak_shifts = []
		self.peak_shift = None
		self.peak_size = None
		self.peak_size_std = None
		self.plus_model = None
		self.minus_model = None
		self.peaks_incorporated = 0
		self.build()

	def build(self):
		# for all chromosomes look for shifted peaks
		for chrom in self.tags.get_chrom_names():
			plus_peaks = self.find_simple_peaks(chrom, '+')
			minus_peaks = self.find_simple_peaks(chrom, '-')
			self.determine_shifts(plus_peaks, minus_peaks)
		# calculate the meidan peak_shift
		if self.peak_shifts:
			self.peak_shift = int(median(self.peak_shifts))
			# peak size is 2 * shift size + 1
			self.peak_size = self.peak_shift * 2 + 1
			self.plus_model = [1] * self.peak_shift  + [0] * (self.peak_shift + 1)
			self.minus_model = [0] * (self.peak_shift + 1) + [1] * self.peak_shift
			self.peak_size_std = std([shift * 2 + 1 for shift in self.peak_shifts])

	def find_simple_peaks(self, chrom, strand):
		# return maxima of tag counts in regions with more tags than threshold
		tags = self.tags.get_tags(chrom, strand)
		window = deque([])
		peak_position_list = []
		peak_region = []
		for tag in tags:
			# add a new tag to the window and reposition it
			window.append(tag)
			window_start = tag - self.window_size
			# get rid of all the tags not fitting in the window
			while window[0] < window_start:
				window.popleft()
			# identify maxima of enriched regions
			tag_count = len(window)
			if tag_count > self.tag_threshold:
				position = tag - self.window_size / 2
				peak_region.append((tag_count, position))
			elif peak_region:
				peak_position_list.append(max(peak_region)[1])
				peak_region = []
		return peak_position_list
	
	def determine_shifts(self, plus_peaks, minus_peaks):
		# looks for minus peaks upstream of plus peaks within fragment size
		minus_peaks = deque(minus_peaks)
		for plus_peak in plus_peaks:
			while minus_peaks:
				minus_peak = minus_peaks[0]
				if minus_peak > plus_peak:
					peak_shift = minus_peak - plus_peak
					if peak_shift < self.window_size * 2:
						# only append if in agreement with max fragment size
						self.peak_shifts.append(peak_shift)
						self.peaks_incorporated += 1
					break
				minus_peaks.popleft()


class Peak:
	# class for peak related infromation and fuctions
	def __init__(self):
		self.size = 0
		self.shift = 0
		self.position = None
		self.tags = ([],[])
		self.score = 0
		self.background = 0
		self.fold_enrichment = 0
		self.plus_freq_dist = None
		self.minus_freq_dist = None
		self.fdr = None
		self.dist_score = None
		self.survivals = 0

	def __len__(self):
		# for truth testing and number of tags
		return int(self.score)
	
	def calc_fold_enrichment(self, tag_threshold, avg_tag_density):
		# calculate fold enrichment minus local background over avg background
		if self.background > tag_threshold * 2: # can do this smarter?
			# if background is really hight dont consider this region
			self.fold_enrichment = 0
		else:
			self.fold_enrichment = (self.score - self.background) / avg_tag_density
		
	def get_score(self):
		# return tag score
		return self.score * self.dist_score
		
	def determine_tag_distribution(self, filter_width):
		# return smoothed frequency distribution position of tags
		# normalize tags for position
		plus_tags = [tags - self.position + self.shift for tags in self.tags[0]]
		minus_tags = [tags - self.position + self.shift for tags in self.tags[1]]
		plus_dist = [0] * (self.size)
		minus_dist = [0] * (self.size)
		# project tags to list
		for i in plus_tags:
			plus_dist[i] += 1
		for i in minus_tags:
			minus_dist[i] += 1
		# use a flat moving window to improve S/N ratio
		window=ones(filter_width,'d')
		# smooth by convolution of the singal with the window
		self.plus_freq_dist = list(convolve(window/window.sum(), plus_dist, mode='same'))
		self.minus_freq_dist = list(convolve(window/window.sum(), minus_dist, mode='same'))
		# normalize distribution height
		norm_factor = (sum(self.plus_freq_dist) + sum(self.minus_freq_dist)) / self.size
		for i in range(self.size):
			self.plus_freq_dist[i] = self.plus_freq_dist[i]/norm_factor
			self.minus_freq_dist[i] = self.minus_freq_dist[i]/norm_factor
		
	def calc_distribution_score(self, plus_model, minus_model):
		# calculate the percentage of similarity between tag dist and model
		overlaps = []
		for i in range(len(plus_model)):
			peak_freq = self.plus_freq_dist[i]
			model_freq = plus_model[i]
			overlap = 1 - (fabs(model_freq - peak_freq)/(model_freq + peak_freq))
			overlaps.append(overlap)
		self.plus_dist_score = sum(overlaps) / len(overlaps)
		
		overlaps = []
		for i in range(len(minus_model)):
			peak_freq = self.minus_freq_dist[i]
			model_freq = minus_model[i]
			overlap = 1 - (fabs(model_freq - peak_freq)/(model_freq + peak_freq))
			overlaps.append(overlap)
		self.minus_dist_score = sum(overlaps) / len(overlaps)
		
		self.dist_score = (self.plus_dist_score + self.minus_dist_score) / 2

class PeakContainer:
	# a class to identify and classify potential peaks
	def __init__(self, ip_tags, control_tags, peak_size, fold_threshold, plus_model, minus_model):
		self.ip_tags = ip_tags
		self.control_tags = control_tags
		self.peak_size = peak_size
		self.peak_shift = (peak_size - 1) / 2
		self.tag_threshold = ip_tags.tag_number / float(ip_tags.genome_size()) * peak_size * fold_threshold
		self.avg_tag_density = ip_tags.tag_number / float(ip_tags.genome_size()) * peak_size
		self.cov_coefficient = ip_tags.tag_number / float(control_tags.tag_number)
		self.plus_model = plus_model
		self.minus_model = minus_model
		self.peaks = {}
		self.peak_count = 0
		self.plus_window = deque([])
		self.minus_window = deque([])
		self.position = 0
		self.build()

	def build(self):
		# perform main peak finding tasks
		for chrom in self.ip_tags.get_chrom_names():
			self.find_peaks(chrom)
			self.measure_background(chrom)
			self.determine_fold_enrichment(chrom)

	def calculate_score(self):
		score = 0
		for tag in self.plus_window:
			score += self.plus_model[tag - self.position + self.peak_shift]
		for tag in self.minus_window:
			score += self.minus_model[tag - self.position + self.peak_shift]
		return score
	
	def find_peaks(self, chrom):
		# identify peak candidates on chromosome
		self.peaks[chrom] = []
		# convert tag arrays to deque for fast appending and popping
		plus_tags = deque(self.ip_tags.get_tags(chrom, '+'))
		minus_tags = deque(self.ip_tags.get_tags(chrom, '-'))
		# initalize windows and stuff
		score_buffer = deque([])
		peak_candidate = Peak()
		# reset scanning windows and position on chromosome
		self.plus_window = deque([])
		self.minus_window = deque([])
		self.position = 0
		while plus_tags and minus_tags:
			# fill windows
			while plus_tags and plus_tags[0] <= (self.position + self.peak_shift):
				self.plus_window.append(plus_tags.popleft())
			while minus_tags and minus_tags[0] <= (self.position + self.peak_shift):
				self.minus_window.append(minus_tags.popleft())
			# get rid of old tags not fitting in the window any more
			while self.plus_window and self.plus_window[0] < (self.position - self.peak_shift):
				self.plus_window.popleft()
			while self.minus_window and self.minus_window[0] < (self.position - self.peak_shift):
				self.minus_window.popleft()
			# add position to region if over threshold
			score = self.calculate_score()
			if score > self.tag_threshold:
				# save all scores in buffer
				score_buffer.append(score)
				# get rid of old scores that are outside of the filter
				if len(score_buffer) > self.peak_size:
					score_buffer.popleft()
				# if current score is as big or bigger, consider it instead
				if score >= peak_candidate.score:
					peak_candidate.size = self.peak_size
					peak_candidate.shift = self.peak_shift
					peak_candidate.score = score
					peak_candidate.tags = (list(self.plus_window), list(self.minus_window))
					peak_candidate.survivals = 0
					if self.position >= 0:
						peak_candidate.position = self.position
					else:
						peak_candidate.position = 0
				# candidate survives if current score is smaller
				else:
					peak_candidate.survivals += 1
				# if candidate survives long enough do the expensive lookup
				if peak_candidate.survivals == self.peak_shift:
					# check score buffer to see whether candidate is a maximum
					# candidate is in the middle of the buffer now
					if peak_candidate.score == max(score_buffer):
						self.add_peak(peak_candidate, chrom)
					# consider current score next, reset survivals
					peak_candidate = Peak()
				# while in enriched region move windows in 1 bp steps
				self.position += 1
			else:
				# if we still have a candidate check whether its a max and add
				if peak_candidate:
					if peak_candidate.score == max(score_buffer):
						self.add_peak(peak_candidate, chrom)
					peak_candidate = Peak()
					score_buffer = deque([])
				# determine the next informative position in the genome and move there
				if plus_tags and minus_tags:
					distance_to_next = plus_tags[0] - self.position + 1
				self.position += distance_to_next

	def add_peak(self, peak, chrom):
		# calculate tag distribution frequency and add peak to container
		peak.determine_tag_distribution(11) # 11 is emp optimal window width
		self.peaks[chrom].append(peak)
		self.peak_count += 1
	
	def measure_background(self, chrom):
		# for every peak check background level
		plus_tags = deque(self.control_tags.get_tags(chrom, '+'))
		minus_tags = deque(self.control_tags.get_tags(chrom, '-'))
		# convert to deque for super fast and efficient popleft
		self.plus_window = deque([])
		self.minus_window = deque([])
		for peak in self.peaks[chrom]:
			# fill windows
			while plus_tags and plus_tags[0] <= (peak.position + self.peak_shift):
				self.plus_window.append(plus_tags.popleft())
			while minus_tags and minus_tags[0] <= (peak.position + self.peak_shift):
				self.minus_window.append(minus_tags.popleft())
			# get rid of old tags not fitting in the window any more
			while self.plus_window and self.plus_window[0] < (peak.position - self.peak_shift):
				self.plus_window.popleft()
			while self.minus_window and self.minus_window[0] < (peak.position - self.peak_shift):
				self.minus_window.popleft()
			# calculate normalized background level
			# add position to region if over threshold
			self.position = peak.position
			background = self.calculate_score()
			nrom_background = background * self.cov_coefficient
			peak.background = nrom_background
	
	def determine_fold_enrichment(self, chrom):
		# for evey peak calculate fold enrichment
		for chrom in self.peaks.keys():
			for peak in self.peaks[chrom]:
				peak.calc_fold_enrichment(self.tag_threshold, self.avg_tag_density)
	
	def model_tag_distribution(self):
		# pick top 200 peak positions and build model
		ranked_peak_tags = []
		for chrom in self.peaks.keys():
			for peak in self.peaks[chrom]:
				ranked_peak_tags.append((peak.score, (peak.plus_freq_dist, peak.minus_freq_dist)))
		# find the tag count of the 200th largest peak
		tag_threshold = sorted(ranked_peak_tags)[-200][0]
		# add tags from highest peaks to the model
		top_tags = [i[1] for i in ranked_peak_tags if i[0] > tag_threshold]
		plus_model = [0] * self.peak_size
		minus_model = [0] * self.peak_size
		for tags in top_tags:
			plus_tags = tags[0]
			minus_tags = tags[1]
			plus_model = map(add, plus_tags, plus_model)
			minus_model = map(add, minus_tags, minus_model)
		# nromalize model for number of total peaks
		norm_factor = (sum(plus_model) + sum(minus_model)) / self.peak_size
		for i in range(self.peak_size):
			plus_model[i] = plus_model[i]/norm_factor
			minus_model[i] = minus_model[i]/norm_factor
		return (plus_model, minus_model)

	def determine_distribution_scores(self, plus_model, minus_model):
		# calculate distribution similarity of every peak to model
		for chrom in self.peaks.keys():
			for peak in self.peaks[chrom]:
				peak.calc_distribution_score(plus_model, minus_model)

	def calculate_fdr(self, control_peaks):
		# create a dictionary to correlate scores with FDR values
		score2fdr = {}
		ip_scores = []
		control_scores = []
		for chrom in self.peaks.keys():
			for peak in self.peaks[chrom]:
				ip_scores.append(peak.get_score())
		for chrom in control_peaks.keys():	
			for peak in control_peaks[chrom]:
				control_scores.append(peak.get_score())
		ip_scores = deque(sorted(ip_scores, reverse=True))
		control_scores = deque(sorted(control_scores, reverse=True))
		# calculate FDR at all relevant cutoffs
		ip_count = float(0)
		control_count = float(0)
		while ip_scores:
			ip_score = ip_scores.popleft()
			ip_count += 1
			while control_scores and control_scores[0] >= ip_score:
				control_scores.popleft()
				control_count +=1
			ip_fdr = control_count / ip_count * 100
			score2fdr[str(ip_score)] = ip_fdr
		# add fdr to each peak object
		for chrom in self.peaks.keys():
			for peak in self.peaks[chrom]:
				peak.fdr = score2fdr[str(peak.get_score())]

	def write_to_stdout(self, options):
		# write results to stdout
		sys.stdout.write('Chromosome\tStart\tEnd\tName\tSummit\tScore\tFoldEnrichmen\tDistributionScore\tFDR\n')
		peak_count = 0
		for chrom in sorted(self.peaks.keys()):
			for peak in self.peaks[chrom]:
				if peak.fdr <= options.fdr:
					peak_count += 1
					summit = peak.position
					start = summit - self.peak_shift
					end = summit + self.peak_shift
					name = chrom + '_Peak_' + str(peak_count)
					score = peak.get_score()
					enrichment = peak.fold_enrichment
					dist_score = peak.dist_score
					fdr = peak.fdr
					output = (chrom, start, end, name, summit, score, enrichment, dist_score, fdr)
					sys.stdout.write('%s\t%d\t%d\t%s\t%d\t%.2f\t%.2f\t%.2f\t%.2f\n' % output)
		print_status('%d peaks detected at FDR %.1f%% \n' % (peak_count, options.fdr), options.verbose)

def print_status(string, boolean):
	# switchable printing to stderror
	if boolean:
		sys.stderr.write('%s %s\n' % (strftime("%H:%M:%S", localtime()), string))
	
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("Program canceled by user!\n")
        sys.exit(0)
