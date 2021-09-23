import os
from typing import *
from argparse import ArgumentParser

import mne

def add_logging_options(parser: ArgumentParser) -> ArgumentParser:
	loglevel_opts = parser.add_mutually_exclusive_group()
	loglevel_opts.add_argument('-q', '--quiet', action='store_true', help='do not show messages')
	loglevel_opts.add_argument('-v', '--verbose', action='count', help='show additional info')

def add_freesurfer_options(parser: ArgumentParser) -> ArgumentParser:
	parser.add_argument('-s', '--subject', metavar='ID', required=True, help='FreeSurfer subject ID')
	parser.add_argument('-d', '--subjects-dir', metavar='DIR', help='FreeSurfer subjects folder path')
	return parser

def add_output_options(parser: ArgumentParser) -> ArgumentParser:
	parser.add_argument('-z', '--gzip', action='store_true', help='use gzip to compress output data')
	parser.add_argument('-o', '--output', metavar='NAME', required=True, help='used as template for output files')
	return parser

def make_output_filename(name: str, type: str = None, ext: str='fif', compress: bool = False) -> str:
	fname = name
	if type:
		fname += f'_{type}'
	if ext:
		fname += f'.{ext}'
	if compress:
		fname += '.gz'
	return fname

def parse_optional(t: Callable[[str], any]) -> Callable[[str], any]:
	return lambda x: None if x == 'none' else t(x)

def parse_with_mapping(t: Callable[[str], any], mapping: Mapping[str, any]) -> Callable[[str], any]:
	return lambda x: mapping.get(x) or t(x)

def parse_either_as(*ts: Iterable[Callable[[str], any]]) -> Callable[[str], any]:
	
	def __helper(x: str) -> any:
		for t in ts:
			try:
				return t(x)
			except:
				pass
		return None
	
	return __helper

def process_freesurfer_options(args: dict) -> None:
	subjects_dir = os.getenv('SUBJECTS_DIR')
	if args.subjects_dir:
		mne.set_config('SUBJECTS_DIR', args.subjects_dir)
		subjects_dir = args.subjects_dir

def process_logging_options(args: dict) -> None:
	if args.quiet:
		mne.set_log_level(False)
	elif args.verbose:
		if args.verbose == 1:
			mne.set_log_level('INFO')
		elif args.verbose >= 2:
			mne.set_log_level('DEBUG')

def use_first_as_default(choices: List[str]) -> dict:
	return {'choices': choices, 'default': choices[0]}
