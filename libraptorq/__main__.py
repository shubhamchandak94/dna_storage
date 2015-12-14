#!/usr/bin/env python2
from __future__ import print_function

import itertools as it, operator as op, functools as ft
from os.path import dirname, basename, exists, isdir, join, abspath
import os, sys, types, math, json, base64


try: import libraptorq
except ImportError:
	# Make sure tool works from a checkout
	if __name__ != '__main__': raise
	pkg_root = abspath(dirname(__file__))
	for pkg_root in pkg_root, dirname(pkg_root):
		if isdir(join(pkg_root, 'libraptorq'))\
				and exists(join(pkg_root, 'setup.py')):
			sys.path.insert(0, dirname(__file__))
			try: import libraptorq
			except ImportError: pass
			else: break
	else: raise ImportError('Failed to find/import "libraptorq" module')

from libraptorq import RQEncoder, RQDecoder, RQError


sys.stdout, sys.stderr = (
	os.fdopen(s.fileno(), 'wb', 0) for s in [sys.stdout, sys.stderr] )
p = lambda fmt,*a,**k:\
	print(*( [fmt.format(*a,**k)]\
		if isinstance(fmt, types.StringTypes) and (a or k)
		else [[fmt] + list(a), k] ), file=sys.stderr)

b64_encode = base64.urlsafe_b64encode
b64_decode = lambda s:\
	base64.urlsafe_b64decode(bytes(s))\
		if '-' in s or '_' in s else bytes(s).decode('base64')


def main(args=None, error_func=None):
	import argparse
	parser = argparse.ArgumentParser(
		description='Encode/decode data using RaptorQ rateless'
			' erasure encoding ("fountain code") algorithm, using libRaptorQ through CFFI.')
	parser.add_argument('--debug', action='store_true', help='Verbose operation mode.')
	cmds = parser.add_subparsers( dest='cmd',
		title='Supported operations (have their own suboptions as well)' )


	cmd = cmds.add_parser('encode',
		help='Encode file into chunks and dump these along with OTI parameters as a JSON structure.')
	cmd.add_argument('path_src', nargs='?',
		help='Path to a file which contents should be encoded. Stdin will be used, if not specified.')
	cmd.add_argument('path_dst', nargs='?',
		help='Path to write resulting JSON to. Will be dumped to stdout, if not specified.')

	cmd.add_argument('-j', '--threads',
		type=int, metavar='n',
		help='Number of encoder threads to use. 0 to scale to all cpus (default).')
	# cmd.add_argument('-t', '--rq-type',
	# 	default='32', metavar='{ NONE | 8 | 16 | 32 | 64 }',
	# 	help='No idea what it means, see RFC6330. Default: %(default)s')
	cmd.add_argument('-k', '--min-subsymbol-size',
		type=int, default=8, metavar='bytes',
		help='No idea what it means, see RFC6330. Default: %(default)s')
	cmd.add_argument('-s', '--symbol-size',
		type=int, default=16, metavar='bytes',
		help='No idea what it means, see RFC6330. Default: %(default)s')
	cmd.add_argument('-m', '--max-memory',
		type=int, default=200, metavar='megabytes?',
		help='Uh... max memory in megs? Not sure. Default: %(default)s')

	cmd.add_argument('-n', '--repair-symbols-rate',
		default=0, type=float, metavar='float',
		help='Number of extra symbols to generate above what is required'
				' to reassemble to file as a fraction of that "required" count.'
			' For example, if 100 symbols are required, 0.5 will generate 150 symbols.'
			' Default is to only generate required amount (i.e. "-n 0").')

	cmd.add_argument('-d', '--drop-rate',
		default=0, type=float, metavar='float',
		help='Drop specified randomly-picked fraction'
			' of resulting chunks (incl. ones for repair),'
			' i.e. just discard these right before output. Mainly useful for testing.')


	cmd = cmds.add_parser('decode', help='Decode lines of base64 into a file.')
	cmd.add_argument('path_src', nargs='?',
		help='Path to a file with JSON structure, such as produced by "encode" operation.'
			' Stdin will be used, if not specified.')
	cmd.add_argument('path_dst', nargs='?',
		help='Path to write assembled file to. Will be dumped to stdout, if not specified.')


	opts = parser.parse_args(sys.argv[1:] if args is None else args)

	global log
	import logging
	logging.basicConfig(
		format='%(asctime)s :: %(levelname)s :: %(message)s',
		datefmt='%Y-%m-%d %H:%M:%S',
		level=logging.DEBUG if opts.debug else logging.WARNING )
	log = logging.getLogger()

	src = sys.stdin if not opts.path_src else open(opts.path_src, 'rb')
	try: data = src.read()
	finally: src.close()


	if opts.cmd == 'encode':
		with RQEncoder( data,
				opts.min_subsymbol_size, opts.symbol_size, opts.max_memory ) as enc:
			enc.precompute(opts.threads)

			symbols, n_block_reps_total = dict(), 0

			n_blocks = enc.blocks
			for sbn in xrange(n_blocks):

				n_block_syms = enc.symbols(sbn)
				for esi in xrange(n_block_syms):
					sym_id = enc.sym_id(esi, sbn)
					symbols[sym_id] = enc.encode(sym_id)

				n_block_reps = int(math.ceil(n_block_syms * opts.repair_symbols_rate))
				n_block_reps_max = enc.max_repair(sbn)
				if n_block_reps_max < n_block_reps:
					log.warn( 'Number of repair symbols is'
						' above max allowed: %s > %s', n_block_reps, n_block_reps_max )
				n_block_reps = min(n_block_reps, n_block_reps_max)
				n_block_reps_total += n_block_reps
				for esi in xrange(n_block_reps):
					sym_id = enc.sym_id(n_block_syms + esi, sbn)
					symbols[sym_id] = enc.encode(sym_id)
			oti_scheme, oti_common = enc.oti_scheme, enc.oti_common

		n_dropped = 0
		if opts.drop_rate > 0:
			import random
			n_dropped = int(round(len(symbols) * opts.drop_rate, 0))
			for n in xrange(n_dropped):
				k = random.choice(symbols.keys())
				del symbols[k]

		data = json.dumps(
			dict( oti_scheme=oti_scheme, oti_common=oti_common,
				symbols=dict((k, b64_encode(v)) for k,v in symbols.viewitems()) ),
			sort_keys=True, indent=2, separators=(',', ': ') )
		log.debug(
			'Encoded %s block(s), %s symbol(s) total (%s for repair). Dropped %s symbol(s).',
			n_blocks, len(symbols) + n_dropped, n_block_reps_total, n_dropped )


	elif opts.cmd == 'decode':
		data = json.loads(data)
		n_symbols, n_discarded = len(data['symbols']), 0
		with RQDecoder(data['oti_common'], data['oti_scheme']) as dec:
			for sym_id, sym in data['symbols'].viewitems():
				sym_id, sym = int(sym_id), b64_decode(sym)
				try: dec.add_symbol(sym, sym_id)
				except Exception as err:
					# log.debug('Failed to add symbol - %s', err) # XXX: not sure if this should happen
					n_discarded += 1
			data = dec.decode()
		log.debug( 'Decoded %sB of data from %s symbols'
			' (total, discarded: %s)', len(data), n_symbols, n_discarded )


	else: raise NotImplementedError(opts.cmd)


	dst = sys.stdout if not opts.path_dst else open(opts.path_dst, 'wb')
	try: dst.write(data)
	finally: dst.close()


if __name__ == '__main__': sys.exit(main())