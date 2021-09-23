#!/usr/bin/python3

import sys
from typing import *
from argparse import ArgumentParser

import mne
import numpy as np
import matplotlib.pyplot as plt
from mne.preprocessing import ICA
from mne.io import Raw, read_raw

from .common import *

ICA_METHODS = [
	'fastica', 'infomax', 'picard',
]

DATA_TYPES = [
	'eeg', 'meg', 'ieeg',
]

MONTAGES = {
	'standard_1020': [
		 'Fp1', 'Fp2',   'F9',   'F7',  'F3',  'Fz',
		  'F4',  'F8',  'F10',  'FT9', 'FT7', 'FT8',
		'FT10',  'T9',   'T7',   'C3',  'Cz',  'C4',
		  'T8', 'T10',  'TP9',  'TP7', 'CPz', 'TP8',
		'TP10',  'P9',   'P7',   'P3',  'Pz',  'P4',
		  'P8', 'P10',  'PO7',  'PO8',  'O1',  'O2',
		  'A1',  'A2', 'ECG1', 'ECG2',
	],
}

INTERP_MODES = [
	'accurate', 'fast'
]

INTERP_METHODS = [
	'spline', 'MNE'
]

def make_argument_parser() -> ArgumentParser:
	MONTAGE_NAMES = [k for k in MONTAGES.keys()]
	parse_notch_filter = parse_optional(parse_with_mapping(float, {'europe': 50, 'usa': 60}))
	
	parser = ArgumentParser(description='M/EEG signal data preprocessing tool')
	parser.add_argument('input', metavar='INPUT', help='input file to analyze')
	parser.add_argument('-t', '--type', help='type of analized data', **use_first_as_default(DATA_TYPES))
	parser.add_argument('-i', '--interactive', action='store_true', help='show windows for interactive editing')
	parser.add_argument('--ica', action='store_true', help='perform ICA decomposition')
	parser.add_argument('--resample', metavar='FREQ', type=int, help='resample data at given frequency')
	add_output_options(parser)
	add_logging_options(parser)
	
	# Montage options
	montage_opts = parser.add_argument_group('montage options')
	montage_opts.add_argument('-m', '--montage', help='apply selected montage', **use_first_as_default(MONTAGE_NAMES))
	montage_opts.add_argument('-c', '--channels', metavar='CHANNEL', nargs='+', help='include selected channels')
	montage_opts.add_argument('-e', '--exclude', metavar='CHANNEL', nargs='+', help='exclude selected channels')
	montage_opts.add_argument('-r', '--reference', metavar='CHANNEL', nargs='+', default='average', help='reference channels to use')
	montage_opts.add_argument('--ecg', metavar='CHANNEL', nargs='+', help='mark channels as ECG channels')
	montage_opts.add_argument('--eog', metavar='CHANNEL', nargs='+', help='mark channels as EOG channels')
	montage_opts.add_argument('--stim', metavar='CHANNEL', nargs='+', help='mark channels as STIM channels')
	montage_opts.add_argument('--bads', metavar='CHANNEL', nargs='+', help='mark channels as bad channels')
	
	# Filter options
	filter_opts = parser.add_argument_group('filter options')
	filter_opts.add_argument('-fl', '--filter-low', metavar='FREQ', type=parse_optional(float), default=.3, help='low frequency filter in Hz')
	filter_opts.add_argument('-fh', '--filter-high', metavar='FREQ', type=parse_optional(float), default=70, help='high frequency filter in Hz')
	filter_opts.add_argument('-fn', '--filter-notch', metavar='FREQ', type=parse_notch_filter, default=50, help='notch frequency filter in Hz')
	
	# Interpolation options
	interp_opts = parser.add_argument_group('channel interpolation options')
	interp_opts.add_argument('--interpolate-bads', action='store_true', default=True, help='interpolate channels marked as bad')
	interp_opts.add_argument('--interpolate-mode', help='channel interpolation mode', **use_first_as_default(INTERP_MODES))
	interp_opts.add_argument('--interpolate-method', choices=INTERP_METHODS, help = 'channel interpolation method to use')
	
	# ICA options
	ica_opts = parser.add_argument_group('ICA options')
	pca_opts = ica_opts.add_mutually_exclusive_group()
	pca_opts.add_argument('--pca-count', metavar='COUNT', type=int, help='number of PCA components to keep')
	pca_opts.add_argument('--pca-variance', metavar='PERCENTAGE', type=float, help='keep as many PCA components as to explain given variance')
	ica_opts.add_argument('--ica-seed', metavar='SEED', type=int, help='seed for deterministic ICA results')
	ica_opts.add_argument('--ica-filter-low', metavar='FREQ', type=float, default=1, help='low frequency filter before ICA analysis')
	ica_opts.add_argument('--ica-method', help='ICA method to use', **use_first_as_default(ICA_METHODS))
	ica_opts.add_argument('--ica-exclude', metavar='INDEX', type=int, nargs='*', help='ICA components to exclude from result')
	ica_opts.add_argument('--ica-exclude-ecg', action='store_true', help='exclude ICA sources detected as ECG')
	ica_opts.add_argument('--ica-exclude-eog', action='store_true', help='exclude ICA sources detected as EOG')
	ica_opts.add_argument('--ica-show-pca-variance', action='store_true', help='show cumulative PCA variance')
	ica_opts.add_argument('--ica-show-sources', action='store_true', help='show ICA components')
	return parser

def normalize_channel_names(raw: Raw) -> Raw:
	mapping = {ch: ch.replace('FP', 'Fp') for ch in raw.ch_names
	                                      if ch.startswith('FP')}
	if 'OI' in raw.ch_names: mapping['OI'] = 'A1'
	if 'OD' in raw.ch_names: mapping['OD'] = 'A2'
	if mapping: raw.rename_channels(mapping)
	return raw

def load_file(fname: str) -> Optional[Raw]:
	raw = read_raw(fname, preload=True)
	normalize_channel_names(raw)
	return raw

def apply_signal_filters(raw: Raw, low_filter: float = None, high_filter: float = None, notch_filter: float = None) -> None:
	filtered = raw.filter(l_freq=low_filter, h_freq=high_filter, picks='data')
	if notch_filter is not None:
		filtered.notch_filter(notch_filter, picks='data')

def perform_ica_analysis(raw: Raw, method: str, pca_components: float = None, seed: float = None) -> ICA:
	ica = ICA(max_iter='auto', n_components=pca_components, random_state=seed)
	return ica.fit(raw)

def is_ica_analysis_required(args: Mapping[str, any]) -> bool:
	return args.ica             or \
	       args.ica_exclude     or \
	       args.ica_exclude_ecg or \
	       args.ica_exclude_eog

def plot_pca_cumulative_variance(ica: ICA, show: bool = True, block: bool = False) -> plt.Figure:
	pca_total_variance = sum(ica.pca_explained_variance_)
	x = np.arange(1, len(ica.pca_explained_variance_) + 1)
	y = np.cumsum(ica.pca_explained_variance_) / pca_total_variance
	plt.xlabel('number of components')
	plt.ylabel('cumulative explained variance')
	plt.plot(x, y)
	
	if show:
		plt.show(block=block)

def main() -> Optional[int]:
	parser = make_argument_parser()
	if len(sys.argv) == 1:
		parser.print_usage(file=sys.stderr)
		return -1
	
	args = parser.parse_args()
	process_logging_options(args)
	
	raw = load_file(args.input)
	channels = MONTAGES.get(args.montage)
	raw.pick_channels(channels, ordered=True)
	
	if args.ecg is None:
		raw.set_channel_types({ch: 'ecg' for ch in raw.ch_names
		                                 if ch.startswith('ECG')})
	else:
		raw.set_channel_types({ch: 'ecg' for ch in args.ecg})
	
	if args.eog is not None:
		raw.set_channel_types({ch: 'eog' for ch in args.eog})
	if args.stim is not None:
		raw.set_channel_types({ch: 'stim' for ch in args.stim})
	if args.bads is not None:
		raw.info['bads'].extend(args.bads)
	
	if args.channels:
		raw.pick(args.channels)
	if args.exclude:
		raw.drop_channels(args.exclude)
	
	average_proj = args.reference == 'average'
	raw.set_montage(args.montage, on_missing='warn')
	raw.set_eeg_reference(args.reference, average_proj)
	if average_proj:
		raw.apply_proj()

	if args.filter_notch:
		freq = args.filter_notch
		fnotch_freqs = np.arange(freq, 5 * freq + 1, freq)
	else:
		fnotch_freqs = None
	
	apply_signal_filters(raw, low_filter=args.filter_low,
	                          high_filter=args.filter_high,
	                          notch_filter=fnotch_freqs)

	if args.resample:
		raw.resample(sfreq=args.resample)
	
	if args.interactive:
		raw.plot(block=True)
	
	if args.interpolate_bads:
		raw.interpolate_bads(mode=args.interpolate_mode, method=args.interpolate_method)
	
	if is_ica_analysis_required(args):
		filtered = raw.copy().pick('data')
		filtered.filter(args.ica_filter_low, None)
		pca_components = args.pca_count or args.pca_variance
		ica = perform_ica_analysis(filtered, method=args.ica_method,
		                           pca_components=pca_components,
		                           seed=args.ica_seed)

		if args.ica_show_pca_variance:
			plot_pca_cumulative_variance(ica, block=True)
		
		if args.ica_exclude:
			ica.exclude.extend(args.ica_exclude)
		if args.ica_exclude_ecg:
			ecg_ids, _ = ica.find_bads_ecg(raw)
			ica.exclude.extend(ecg_ids)
		if args.ica_exclude_eog:
			eog_ids, _ = ica.find_bads_eog(raw)
			ica.exclude.extend(eog_ids)
		
		if args.ica_show_sources:
			ica.plot_sources(inst=filtered, block=True)
			plt.show()
			
		if args.interactive:
			ica.plot_components(inst=filtered, show=False)
		
		ica.apply(raw, n_pca_components=pca_components)
	
	if args.output:
		out_file = make_output_filename(args.output, args.type, compress=args.gzip)
		raw.save(out_file, overwrite=True)


if __name__ == '__main__':
	try:
		res = main()
		sys.exit(res)
	
	except KeyboardInterrupt:
		print('Interrupted.', file=sys.stderr)
		sys.exit(-1)
	
	except Exception as e:
		print(f'Error: {e}.', file=sys.stderr)
		sys.exit(-1)
