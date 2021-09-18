from typing import *


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
