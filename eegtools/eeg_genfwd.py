#!/usr/bin/python3

import os
import sys
from typing import *
from argparse import ArgumentParser

import mne
from mne.source_space import SourceSpaces

from .common import *

BEM_METHODS = [
	'watershed', 'flash'
]

SOURCE_TYPES = [
	'volume', 'surface'
]

def make_argument_parser() -> ArgumentParser:
	parser = ArgumentParser(description='Generate forward solution to use with DSL')
	parser.add_argument('input', metavar='FIF_FILE', help='input file to compute forward solution for')
	
	add_freesurfer_options(parser)
	parser.add_argument('-b', '--bem-method', help='method to use for generating BEM surfaces', **use_first_as_default(BEM_METHODS))
	parser.add_argument('-j', '--source-spacing', type=parse_either_as(float, str), help='spacing for source spaces')
	parser.add_argument('-p', '--source-type', help='type of source space to generate', **use_first_as_default(SOURCE_TYPES))
	parser.add_argument('-t', '--transform', required=True, metavar='TRANS_FILE', type=parse_optional(str), help='path to transformation FIF file')
	parser.add_argument('-r', '--regenerate', action='store_true', help='recreate intermediate files')
	add_output_options(parser)
	add_logging_options(parser)
	
	return parser

def generate_bem_solution(subject: str, **kwargs) -> mne.bem.ConductorModel:
	surfaces = mne.make_bem_model(subject=subject, **kwargs)
	bem = mne.make_bem_solution(surfaces)
	return bem

def generate_bem_surfaces(subject: str, method: str, **kwargs) -> None:
	if method == 'watershed':
		mne.bem.make_watershed_bem(subject, **kwargs)
	elif method == 'flash':
		mne.bem.make_flash_bem(subject, **kwargs)

def generate_source_space(subject: str, type: str, **kwargs) -> mne.SourceSpaces:
	if type == 'surface':
		return mne.setup_source_space(subject, **kwargs)
	elif type == 'volume':
		return mne.setup_volume_source_space(subject, **kwargs)

def main() -> Optional[int]:
	parser = make_argument_parser()
	if len(sys.argv) == 1:
		parser.print_usage(file=sys.stderr)
		return -1
	
	args = parser.parse_args()
	process_logging_options(args)
	process_freesurfer_options(args)
	
	raw = mne.io.read_raw(args.input)
	
	bem = None
	bem_type = f'{args.bem_method}_bem'
	bem_file = make_output_filename(args.output, bem_type, compress=args.gzip)
	if not args.regenerate:
		try:
			bem = mne.read_bem_solution(bem_file)
		except Exception:
			pass
	
	if not bem:
		subjects_dir = args.subjects_dir or os.getenv('SUBJECTS_DIR')
		ws_dir = os.path.join(subjects_dir, args.subject, 'bem', args.bem_method)
		if not os.path.exists(ws_dir) or args.regenerate:
			generate_bem_surfaces(args.subject, args.bem_method, overwrite=True)
		bem = generate_bem_solution(args.subject)
		mne.write_bem_solution(bem_file, bem, overwrite=True)
	
	src_space = None
	src_type = f'{args.source_spacing}_src'
	src_file = make_output_filename(args.output, src_type, compress=args.gzip)
	if not args.regenerate:
		try:
			src_space = mne.read_source_spaces(src_file)
		except:
			pass
	
	if not src_space:
		if args.source_type == 'surface':
			src_space = generate_source_space(args.subject, 'surface',
			                                  spacing=args.source_spacing)
		if args.source_type == 'volume':
			src_space = generate_source_space(args.subject, 'volume')
		mne.write_source_spaces(src_file, src_space, overwrite=True)
	
	fwd_type = f'{args.source_spacing}_fwd'
	fwd = mne.make_forward_solution(raw.info, args.transform, src_space, bem)
	output_file = make_output_filename(args.output, fwd_type, compress=args.gzip)
	mne.write_forward_solution(output_file, fwd, overwrite=True)


if __name__ == '__main__':
	try:
		res = main()
		sys.exit(res)
	
	except InterruptedError:
		print('Interrupted.', file=sys.stderr)
		sys.exit(-1)
		
	except Exception as e:
		print(f'Error: {e}.', file=sys.stderr)
		sys.exit(-1)
