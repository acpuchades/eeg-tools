import os
from typing import *
from argparse import ArgumentParser

import mne

def add_global_options(parser: ArgumentParser) -> ArgumentParser:
	parser.add_argument('-q', '--quiet', action='store_true', help='do not show messages')

def add_freesurfer_options(parser: ArgumentParser) -> ArgumentParser:
	parser.add_argument('-s', '--subject', metavar='ID', required=True, help='FreeSurfer subject ID')
	parser.add_argument('-d', '--subjects-dir', metavar='DIR', help='FreeSurfer subjects folder path')
	return parser

def add_output_options(parser: ArgumentParser) -> ArgumentParser:
	parser.add_argument('-z', '--gzip', action='store_true', help='use gzip to compress output data')
	parser.add_argument('-o', '--output', metavar='NAME', required=True, help='used as template for output files')
	return parser

def make_output_filename(name: str, type: str = None, compress: bool = False) -> str:
	fname = name
	if type:
		fname += f'_{type}'
	fname += '.fif'
	if compress:
		fname += '.gz'
	return fname

def parse_optional(t: Callable[[str], any]):
	return lambda x: None if x == 'none' else t(x)

def parse_with_mapping(t: Callable[[str], any], mapping: Mapping[str, any]):
	return lambda x: mapping.get(x) or t(x)

def process_freesurfer_options(args: dict) -> None:
	subjects_dir = os.getenv('SUBJECTS_DIR')
	if args.subjects_dir:
		mne.set_config('SUBJECTS_DIR', args.subjects_dir)
		subjects_dir = args.subjects_dir

def process_global_options(args: dict) -> None:
	if args.quiet:
		mne.set_log_level('warning')

def use_first_as_default(choices: List[str]) -> dict:
	return {'choices': choices, 'default': choices[0]}
