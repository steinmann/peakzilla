#!/usr/bin/env python

# Copyright (c) Jonas Steinmann, 2010-2011
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License version 3 as 
# published by the Free Software Foundation.

import sys
import csv
import os
from cPickle import dump
from operator import add
from time import strftime, localtime
from collections import deque
from array import array
from optparse import OptionParser
from copy import copy
from math import exp, sqrt, pi

def main():
	# option parser
	usage = 'python peakzilla.py [OPTIONS] chip.bed control.bed > results.tsv'
	parser = OptionParser(usage=usage)
	
	parser.add_option("-s", "--fragment_size",\
	type = "int", dest="fragment_size", default="200",\
	help = "upper limit of fragment size in bp: default=200")
	
	parser.add_option("-f", "--fdr",\
	type = "float", dest="fdr", default='1',\
	help = "cutoff for the estimated FDR value: default = 1")
	
	parser.add_option("-c", "--score_cutoff",\
	type = "float", dest="score_cutoff", default='2',\
	help = "minimum cutoff for score (~ fold enrichment): default = 2")
	
	parser.add_option("-q", "--quiet",\
	action = "store_false", dest="verbose", default=True,\
	help = "don't print status messages")
	
	# read arguments and options
	(options, args) = parser.parse_args()
	if len(args) != 2:
		# return help message if argment number is incorrect
		parser.print_help()
		sys.exit(0)
	ip_file = args[0]
	control_file = args[1]
	
	# load tags
	print_status('Loading tags ...', options.verbose)
	ip_tags = TagContainer()
	control_tags = TagContainer()
	ip_tags(ip_file)
	control_tags(control_file)
	
	# report tag number
	print_status('Tags in IP: %d' % ip_tags.tag_number, options.verbose)
	print_status('Tags in control: %d' % control_tags.tag_number, options.verbose)

	# first attempt of modeling peak size
	print_status('Modeling peak size and shift ...', options.verbose)
	model_threshold = 120
	peak_model = PeakShiftModel(ip_tags, options.fragment_size, model_threshold)
	
	# change model threshold until it yields a reasonable number of peaks
	while peak_model.peaks_incorporated < 200 or peak_model.peaks_incorporated > 500:
		if peak_model.peaks_incorporated < 500:
			model_threshold = model_threshold / 2
			peak_model = PeakShiftModel(ip_tags, options.fragment_size, model_threshold)
		else:
			model_threshold = model_threshold * 1.25
			peak_model = PeakShiftModel(ip_tags, options.fragment_size, model_threshold)
	print_status('Peak size is %d bp' % peak_model.peak_size, options.verbose)
	
	# find peaks for modeling
	print_status('Finding peaks for modeling ...', options.verbose)
	ip_peaks = PeakContainer(ip_tags, control_tags, peak_model.peak_size, peak_model.plus_model, peak_model.minus_model)
	
	# model tag distribution 
	print_status('Modeling tag distribution ...', options.verbose)
	plus_model = ip_peaks.model_tag_distribution()[0]
	minus_model = ip_peaks.model_tag_distribution()[1]
	
	# plot model using R
	print_status('Plotting the model ...', options.verbose)
	create_model_plot(plus_model, minus_model, ip_file)
	
	# serialize model
	print_status('Serializing the model for later use ...', options.verbose)
	filename = ip_file[:-4] + '_model.pkl'
	dump((plus_model, minus_model), open(filename, 'w'))

	# find peaks using emirical model
	print_status('Finding peaks with empirical peak model ...', options.verbose)
	ip_peaks = PeakContainer(ip_tags, control_tags, peak_model.peak_size, plus_model, minus_model)
	
	# find peaks using emirical model in control sample
	print_status('Finding peaks in control sample ...', options.verbose)
	control_peaks = PeakContainer(control_tags, ip_tags, peak_model.peak_size, plus_model, minus_model)
	
	# calculate distribution scores
	print_status('Calculating tag distribution scores ...', options.verbose)
	ip_peaks.determine_distribution_scores(plus_model, minus_model)
	control_peaks.determine_distribution_scores(plus_model, minus_model)
	 
	# calculate FDR
	print_status('Calculating FDR ...', options.verbose)
	ip_peaks.calculate_fdr(control_peaks.peaks)
	
	# write output as bed files
	ip_peaks.write_to_stdout(options)
	
	# write peaks in input to file
	print_status('Writing input peaks to %s' % control_file[:-4] + '_peaks.tsv', options.verbose)
	control_peaks.write_artifact_peaks(control_file)
	print_status('Done!', options.verbose)


def create_model_plot(plus_model, minus_model, ip_file_name):
	# output model in an R friendly fashion
		filename = ip_file_name[:-4] + '_model.R'
		pdf_name = ip_file_name[:-4] + '_model.pdf'
		try:
			f = open(filename, 'w')
			f.close()
			f = open(filename, 'a')
			f.write('plus = c',)
			f.write(str(tuple(plus_model))+'\n')
			f.write('minus = c',)
			f.write(str(tuple(minus_model))+'\n')
			f.write('pdf(file=\'%s\')\n' % pdf_name)
			f.write('shift = (length(plus) - 1) / 2\n')
			f.write('plot(seq(-shift,shift), plus,type=\'l\', col=\'red\')\n')
			f.write('lines(seq(-shift,shift), minus, col=\'blue\')\n')
			f.write('dev.off()')
			f.close()
		except:
			print_status('WARNING: You dont have write permissions in this folder! No modeldata saved!', True)
		try:
			os.system('R --slave < %s > /dev/null' % filename)
		except:
			print_status('No plot created, need to install R for that!', True)


def print_status(string, boolean):
	# switchable printing to stderror
	if boolean:
		sys.stderr.write('%s %s\n' % (strftime("%H:%M:%S", localtime()), string))


def chisquare(f_obs, f_exp):
	# calculates a one-way chi-square for observed versus exprected frequencies
	chisq = 0
	df = len(f_obs)-1
	for i in range(len(f_obs)):
		chisq = chisq + (f_obs[i]-f_exp[i]) ** 2 / float(f_exp[i])
	# return chi-sqare value and associated p-value for f_obs == f_exp
	return chisq, chisqprob(chisq, df)


def chisqprob(chisq, df):
	# returns chisquare probability (works only for high degrees of freedom)
	if df < 30:
		raise ValueError('Function does not work for df < 30!')
	if chisq < 15:
		return 1.0
	a = 0.5 * chisq
	y = exp(-a)
	chisq = 0.5 * (df - 1.0)
	if df % 2 == 0:
		e = 1.0
		z = 1.0
	else:
		e = 1.0 / sqrt(pi) / sqrt(a)
		z = 0.5
	c = 0.0
	while (z <= chisq):
		e = e * (a/float(z))
		c = c + e
		z = z + 1.0
	return (c*y)


def median(numlist):
	# calculate median
    s = sorted(numlist)
    l = len(numlist)
    if l == 0:
		return float('nan')
    if l%2 == 0:
        return (s[l/2] + s[l/2-1]) / 2.0
    else:
        return float(s[l/2])


def convolve(signal, filter_width):
	# smooth signal with a flat scanning window of filter_width
	filter_width = float(filter_width)
	overhang = int((filter_width-1) / 2)
	window = deque([])
	result = []
	for i in signal + overhang * [0]:
		window.append(i)
		while len(window) > filter_width:
			window.popleft()
		result.append(sum(window)/ filter_width)
	return result[overhang:]


class TagContainer:

	# class for loading, storing and manipulating sequence tags
	def __init__(self):
		# intitialize an empty object
		self.tags = {}
		self.tag_number = 0

	def __call__(self, bed_file):
		# when called like a function load bed file and return self
		self.load_bed(bed_file)
		self.sort_tags()
		return self
	
	def add_tag(self, chrom, strand, fiveprime):
		# add tag to dictionary
		if not chrom in self.tags:
			self.tags[chrom] = {}
			# store tags as an array of unsigned integers (4 bytes)
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
			try:
				chrom = i[0]
				start = int(i[1])
				end = int(i[2])
				strand = i[5]
				# determine five prime end
				if strand == '+':
					fiveprime = start
				elif strand == '-':
					fiveprime = end
				# add tag to container
				self.add_tag(chrom, strand, fiveprime)
			except:
				sys.stderr.write("Input file is not in BED format!\n")
				sys.exit(1)
				
			
	def sort_tags(self):
		# sort all tags while preserving the array
		for chrom in self.tags.keys():
			# as sorted returns conversion back to array is required
			self.tags[chrom]['+'] = array('i', sorted(self.tags[chrom]['+']))
			self.tags[chrom]['-'] = array('i', sorted(self.tags[chrom]['-']))
	
	def get_chrom_size(self, chrom):
		# chromosome size to consider for scanning of both strands
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
					if peak_shift < self.window_size * 2 and peak_shift > self.window_size / 2 :
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
		self.median_score_ip = None
		self.median_score_control = None
		self.fold_enrichment = 0
		self.plus_freq_dist = None
		self.minus_freq_dist = None
		self.fdr = None
		self.dist_score = None
		self.survivals = 0
		self.plus_reg_tags_ip = None
		self.plus_reg_tags_control = None

	def __len__(self):
		# for truth testing and number of tags
		return int(self.score)
	
	def calc_fold_enrichment(self):
		# normalize by median score count in 4 kbp window around peak summit
		cov_norm_factor = self.median_score_ip / self.median_score_control
		# if background is abnormaly low use local score median for fold enrichment calculation
		if self.background < self.median_score_ip:
			self.fold_enrichment = (self.score - self.background * cov_norm_factor) / self.median_score_ip
		# else use median normalized background at exactly the same position in control sample
		else:
			self.fold_enrichment = self.score / (self.background * cov_norm_factor)
	
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
		# smooth by convolution of the singal with the window
		self.plus_freq_dist = convolve(plus_dist, filter_width)
		self.minus_freq_dist = convolve(minus_dist, filter_width)
		# normalize distribution height
		norm_factor = (sum(self.plus_freq_dist) + sum(self.minus_freq_dist)) / self.size
		for i in range(self.size):
			self.plus_freq_dist[i] = self.plus_freq_dist[i]/norm_factor
			self.minus_freq_dist[i] = self.minus_freq_dist[i]/norm_factor

	def calc_distribution_score(self, plus_model, minus_model):
		# concatenate plus and minus distributions and models for testing
		model = plus_model[:self.shift] + minus_model[-self.shift:]
		freq_dist = self.plus_freq_dist[:self.shift] + self.minus_freq_dist[-self.shift:]
		# dist score is the p-value returned by the chi-square test
		self.dist_score = chisquare(freq_dist, model)[1]
		
	def get_score(self):
		# final score is fold enrichment times goodness of fit to model
		return self.fold_enrichment * self.dist_score


class PeakContainer:
	# a class to identify and classify potential peaks
	def __init__(self, ip_tags, control_tags, peak_size, plus_model, minus_model):
		self.ip_tags = ip_tags
		self.control_tags = control_tags
		self.peak_size = peak_size
		self.peak_shift = (peak_size - 1) / 2
		self.score_threshold = 10
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
			self.measure_local_median(chrom)
			self.measure_background(chrom)
			self.determine_fold_enrichment(chrom)

	def calculate_score(self):
		# calculate score
		score = 0
		tag_shift = self.peak_shift - self.position
		plus = self.plus_model
		minus = self.minus_model
		for tag in self.plus_window:
			score += plus[tag + tag_shift]
		for tag in self.minus_window:
			score += minus[tag + tag_shift]
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
			# if number of candidates found is high readjust threshold
			if self.peak_count > 35000:
				self.adjust_threshold()
			# add position to region if over threshold
			score = self.calculate_score()
			if score > self.score_threshold:
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

	def adjust_threshold(self):
		# allows for dynamic adjustment of peak calling threshold
		# restricts the number of candidate peaks to investigate to 30000
		peak_scores = []
		for chrom in self.peaks.keys():
			for peak in self.peaks[chrom]:
				peak_scores.append(peak.score)
		# set score to 30000th
		self.score_threshold = sorted(peak_scores)[-30000]
		# remove peaks below threshold
		for chrom in self.peaks.keys():
			self.peaks[chrom] = [peak for peak in self.peaks[chrom] if peak.score >= self.score_threshold]
		# recount peaks
		self.peak_count = 0
		for chrom in self.peaks.keys():
			self.peak_count += len(self.peaks[chrom])

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
			peak.background = self.calculate_score()
		
	def measure_local_median(self, chrom):
		# measure median score arround a peak +/- 2000 bp for ip tags
		plus_tags = deque(self.ip_tags.get_tags(chrom, '+'))
		minus_tags = deque(self.ip_tags.get_tags(chrom, '-'))
		local_plus_tags = deque([])
		local_minus_tags = deque([])
		for peak in self.peaks[chrom]:
			# fill windows
			while plus_tags and plus_tags[0] <= (peak.position + 2000):
				local_plus_tags.append(plus_tags.popleft())
			while minus_tags and minus_tags[0] <= (peak.position + 2000):
				local_minus_tags.append(minus_tags.popleft())
			# get rid of old tags not fitting in the window any more
			while local_plus_tags  and local_plus_tags [0] < (peak.position - 2000):
				local_plus_tags .popleft()
			while local_minus_tags and local_minus_tags[0] < (peak.position - 2000):
				local_minus_tags.popleft()
			# calculate normalized background level
			# add position to region if over threshold
			peak.plus_reg_tags_ip = len(local_plus_tags) + len(local_minus_tags)
			peak.median_score_ip = self.get_median(copy(local_plus_tags), copy(local_minus_tags), peak.position)	
		# measure median score arround a peak +/- 2000 bp for control tags
		plus_tags = deque(self.control_tags.get_tags(chrom, '+'))
		minus_tags = deque(self.control_tags.get_tags(chrom, '-'))
		local_plus_tags = deque([])
		local_minus_tags = deque([])
		for peak in self.peaks[chrom]:
			# fill windows
			while plus_tags and plus_tags[0] <= (peak.position + 2000):
				local_plus_tags.append(plus_tags.popleft())
			while minus_tags and minus_tags[0] <= (peak.position + 2000):
				local_minus_tags.append(minus_tags.popleft())
			# get rid of old tags not fitting in the window any more
			while local_plus_tags  and local_plus_tags [0] < (peak.position - 2000):
				local_plus_tags .popleft()
			while local_minus_tags and local_minus_tags[0] < (peak.position - 2000):
				local_minus_tags.popleft()
			# calculate normalized background level
			# add position to region if over threshold
			peak.plus_reg_tags_control = len(local_plus_tags) + len(local_minus_tags)
			peak.median_score_control = self.get_median(copy(local_plus_tags), copy(local_minus_tags), peak.position)
			
	def get_median(self, plus_tags, minus_tags, peak_position):
		# calculate score every 10 bp and return median score
		score_list = []
		self.plus_window = deque([])
		self.minus_window = deque([])
		self.position = peak_position - 2000 + self.peak_shift
		region_end = peak_position + 2000 - self.peak_shift
		while self.position <= region_end:
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
			if score > 0:
				# dont include region that is not mappable
				score_list.append(score)
			self.position += 10
		score_median = median(score_list)
		# if no tags are in the region at all just return 1 as median score
		if score_median > 0:
			return score_median
		else:
			return 1

	def determine_fold_enrichment(self, chrom):
		# for evey peak calculate fold enrichment
		for chrom in self.peaks.keys():
			for peak in self.peaks[chrom]:
				peak.calc_fold_enrichment()
	
	def model_tag_distribution(self):
		# use tags from top 200 peaks to build the distribution model
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
		sys.stdout.write('Chromosome\tStart\tEnd\tName\tSummit\tScore\tRawScore\tLocalScoreMedian\tLocalScoreMedianBg\tBackground\tFoldEnrichment\tDistributionScore\tFDR\n')
		peak_count = 0
		for chrom in sorted(self.peaks.keys()):
			for peak in self.peaks[chrom]:
				score = peak.get_score()
				if peak.fdr <= options.fdr and score >= options.score_cutoff:
					peak_count += 1
					summit = peak.position
					start = summit - self.peak_shift
					if start < 0:
						start = 0
					end = summit + self.peak_shift
					name = chrom + '_Peak_' + str(peak_count)
					raw_score = peak.score
					loc_med_score = peak.median_score_ip
					loc_med_bg_score = peak.median_score_control
					background = peak.background
					enrichment = peak.fold_enrichment
					dist_score = peak.dist_score
					fdr = peak.fdr
					output = (chrom, start, end, name, summit, score, raw_score, loc_med_score, loc_med_bg_score, background, enrichment, dist_score, fdr)
					sys.stdout.write('%s\t%d\t%d\t%s\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % output)
		print_status('%d peaks detected at FDR %.1f%% and %.1f fold enrichment cutoff' % (peak_count, options.fdr, options.score_cutoff), options.verbose)
	
	def write_artifact_peaks(self, control_file_name):
		# write peaks found in input to file
		filename = control_file_name[:-4] + '_peaks.tsv'
		f = open(filename, 'w')
		f.write('Chromosome\tStart\tEnd\tName\tSummit\tScore\tRawScore\tLocalScoreMedian\tLocalScoreMedianBg\tBackground\tFoldEnrichment\tDistributionScore\n')
		f.close()
		peak_count = 0
		for chrom in sorted(self.peaks.keys()):
			for peak in self.peaks[chrom]:
				if peak.get_score() > 1:
					peak_count += 1
					summit = peak.position
					start = summit - self.peak_shift
					if start < 0:
						start = 0
					end = summit + self.peak_shift
					name = chrom + '_Peak_' + str(peak_count)
					score = peak.get_score()
					raw_score = peak.score
					loc_med_score = peak.median_score_ip
					loc_med_bg_score = peak.median_score_control
					background = peak.background
					enrichment = peak.fold_enrichment
					dist_score = peak.dist_score
					output = (chrom, start, end, name, summit, score, raw_score, loc_med_score, loc_med_bg_score, background, enrichment, dist_score)
					f = open(filename, 'a')
					f.write('%s\t%d\t%d\t%s\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % output)
					f.close()


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("Program canceled by user!\n")
        sys.exit(0)
