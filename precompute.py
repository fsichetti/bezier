### Bezier subdivision pre-computations ###
from polyfem_indices import index_set_polyfem
import os

## Setup ##
COMBINATIONS = [
	(1,1,1), (1,1,2), (1,1,3), (1,1,4), (1,1,5),
	(2,1,1), (2,1,2),
	(2,2,1), (2,2,2), (2,2,3),
	(3,1,1),
	(3,2,1),
	(3,3,1), (3,3,2),
	(2,2,4),
	(3,3,3),
	(3,3,4),
]
COMBINATIONS = [
	(3,3,4),
]

EXPRESSION_CHUNK_SIZE = {(3,3,4) : 50}

DRY_RUN = False
POLYFEM_ORDER = True
DYNAMIC = True
STATIC = True

WRITE_MATRICES = False
WRITE_LAGVEC = False
WRITE_CORNERS = False
WRITE_LAGPOLY = True

INFO_ORDER = False
INFO_JAC_ORDER = False
INFO_LAGBASIS = False
INFO_LAGVECTOR = False

TO_ORIGIN = False

## Imports

from itertools import product
from sympy import \
	symbols, prod, Matrix, det, Poly, expand, collect, \
	pprint, numbered_symbols, cse, QQ, ZZ
from sympy.printing.cxx import CXX11CodePrinter
from sympy.polys.matrices import DomainMatrix
from concurrent.futures import ProcessPoolExecutor

## Custom printer for rationals ##
class MyCodePrinter(CXX11CodePrinter):
	def _print_Rational(self, expr):
		# return f'R_{expr.p}_{expr.q}'
		return f'R({expr.p}, {expr.q})'
	def _print_Pow(self, expr):
		base, exp = expr.as_base_exp()
		return f"powi({self._print(base)}, {exp})"
		


## Subscripting ##
def subscripts(b, i):
	return symbols([f'{b}[{j}]' for j in i])

## Index set ##
def index_set(n, s, p):
	assert(0 < s <= n)
	if POLYFEM_ORDER: return index_set_polyfem(n,s,p)
	else: return [tup for tup in product(range(p+1), repeat=n)
		if sum(tup[:s]) <= p and all(x <= p for x in tup[s:])]

## Corner indices set ##
def corner_indices_set(n, s, p, dynamic):
	simplex_ord = n * p - s
	tensor_ord = n * p - 1
	assert(simplex_ord >= 0 and tensor_ord >= 0)
	time_ord = n
	en = enumerate(index_set_J(n, s, p, dynamic))
	if dynamic:
		return [ind for ind,tup in en
			if all(x == simplex_ord or x == 0 for x in tup[:s]) \
			and all(x == tensor_ord or x == 0 for x in tup[s:n]) \
			and tup[-1] == time_ord]
	else:
		return [ind for ind,tup in en
			if all(x == simplex_ord or x == 0 for x in tup[:s]) \
			and all(x == tensor_ord or x == 0 for x in tup[s:])]

## Univariate Lagrange polynomial ##
def lagrange_uni(m, p, i, x):
	k_values = [k for k in range(0, p+1) if k != i]
	return prod((m * x - k) / (i - k) for k in k_values)

## Multivariate Lagrange polynomial ##
def lagrange(n, s, p, i, x_name):
	assert(len(i) == n)
	x = subscripts(x_name, range(n))
	i_slack = p - sum(i[:s])
	assert(i_slack >= 0)
	x_slack = 1 - sum(x[:s])
	return expand(
		lagrange_uni(p, i_slack, i_slack, x_slack) * \
		prod([lagrange_uni(p, i[k], i[k], x[k]) for k in range(s)]) * \
		prod([lagrange_uni(p, p, i[k], x[k]) for k in range(s,n)])
	)

## Multivariate Lagrange polynomial ##
def lagrange_J(n, s, p, i, x_name):
	assert(len(i) == n)
	Ps = n*p - s
	Pt = n*p - 1
	x = subscripts(x_name, range(n))
	i_slack = Ps - sum(i[:s])
	assert(i_slack >= 0)
	x_slack = 1 - sum(x[:s])
	return expand(
		lagrange_uni(Ps, i_slack, i_slack, x_slack) * \
		prod([lagrange_uni(Ps, i[k], i[k], x[k]) for k in range(s)]) * \
		prod([lagrange_uni(Pt, Pt, i[k], x[k]) for k in range(s,n)])
	)

## Geometric map component ##
def geo_map_component(n, s, p, t, x_name, c_name, ni):
	indices = index_set(n, s, p)
	if t < 0:
		if (__debug__): print(f"Static geometric map ({n} {s} {p}), component {ni}")
		cp = subscripts(c_name, range(ni, len(indices)*n, n))
	else:
		if (__debug__): print(f"Geometric map ({n} {s} {p}) at time {t}, component {ni}")
		cp = subscripts(c_name, range(t + 2*ni, len(indices)*2*n, 2*n))
	if TO_ORIGIN: cp[0] = 0
	lag = [lagrange(n, s, p, k, x_name) for k in indices]
	return expand(sum(cp[k] * lag[k] for k in range(len(indices))))
	
## Geometric map ##
def aux_map_component(args):
	n, s, p, t, x_name, c_name, i = args
	return geo_map_component(n, s, p, t, x_name, c_name, i)

def aux_combine_results(args):
	i, T, results_0, results_1 = args
	return (1 - T) * results_0[i] + T * results_1[i]

def geo_map(n, s, p, x_name, c_name, t_name=""):
	# Dynamic case
	if len(t_name) > 0:
		T = symbols(t_name)

		with ProcessPoolExecutor() as executor:
			args_0 = [(n, s, p, 0, x_name, c_name, i) for i in range(n)]
			args_1 = [(n, s, p, 1, x_name, c_name, i) for i in range(n)]

			# Submit tasks for t=0 and t=1 components
			futures_0 = [executor.submit(aux_map_component, arg) for arg in args_0]
			futures_1 = [executor.submit(aux_map_component, arg) for arg in args_1]

			results_0 = [future.result() for future in futures_0]
			results_1 = [future.result() for future in futures_1]

		with ProcessPoolExecutor() as executor:
			# Prepare arguments for combining results
			args_combined = [(i, T, results_0, results_1) for i in range(n)]

			futures_combined = [executor.submit(aux_combine_results, arg)
				for arg in args_combined]
			results_combined = [future.result() for future in futures_combined]

		return results_combined + [T]

	# Static case
	else:
		with ProcessPoolExecutor() as executor:
			args = [(n, s, p, -1, x_name, c_name, i) for i in range(n)]
			futures = [executor.submit(aux_map_component, arg) for arg in args]
			results = [future.result() for future in futures]
		return results


## Domain points for the Jacobian determinant ##
def index_set_J(n, s, p, dynamic):
	simplex_ord = n * p - s
	tensor_ord = n * p - 1
	assert(simplex_ord >= 0 and tensor_ord >= 0)
	if dynamic:
		time_ord = n
		return [tup for tup in product(range(n*p+1), repeat=n+1)
			if sum(tup[:s]) <= simplex_ord \
			and all(x <= tensor_ord for x in tup[s:n]) \
			and tup[-1] <= time_ord]
	else:
		return [tup for tup in product(range(n*p+1), repeat=n)
			if sum(tup[:s]) <= simplex_ord \
			and all(x <= tensor_ord for x in tup[s:])]

def domain_pts_J(n, s, p, dynamic):
	simplex_den = max(1, n * p - s)
	tensor_den = max(1, n * p - 1)
	if dynamic:
		time_den = n
		return [
			tuple(
				QQ(elem, simplex_den) if i < s else
				QQ(elem, time_den) if i == n else
				QQ(elem, tensor_den)
				for i, elem in enumerate(tup)
			)
			for tup in index_set_J(n, s, p, True)]
	else:
		return [
			tuple(
				QQ(elem, simplex_den) if i < s else
				QQ(elem, tensor_den)
				for i, elem in enumerate(tup)
			)
			for tup in index_set_J(n, s, p, False)]

## Determinant ##
def determinant(J, dynamic):
	off = 1 if dynamic else 0
	if J.shape == (1 + off,1 + off):
		return J[0,0]
	elif J.shape == (2 + off,2 + off):
		return (J[0,0] * J[1,1]) - (J[0,1] * J[1,0])
	elif J.shape == (3 + off,3 + off):
		return \
			(J[0,0] * ((J[1,1] * J[2,2]) - (J[1,2] * J[2,1]))) +\
			(J[0,1] * ((J[1,2] * J[2,0]) - (J[1,0] * J[2,2]))) +\
			(J[0,2] * ((J[1,0] * J[2,1]) - (J[1,1] * J[2,0])))
	else: raise ValueError

## Lagrange vector ##
def aux_evaluate_det(args):
	pt, J, xt, n = args
	if (__debug__): print(f"\t Det eval {pt} start")
	# if J.shape == (n+1,n+1):
	# 	assert(all(elem == 0 for elem in J[-1, :-1]) and J[-1, -1] == 1)
	# 	J = J[:-1, :-1]
	# assert(J.shape == (n,n))
	# d = J.subs({xt[k]: pt[k] for k in range(len(xt))}).det(method='berkowitz')
	dyn = J.shape == (n+1,n+1)
	d = determinant(J.subs({xt[k]: pt[k] for k in range(len(xt))}), dynamic=dyn)
	if (__debug__): print(f"\t Det eval {pt} end")
	return d

def lagrange_vector(n, s, p, c_name, dynamic):
	x_name = 'u'
	t_name = ''
	xt = subscripts(x_name, range(n))
	if dynamic:
		t_name = 't'
		xt += [symbols(t_name)]

	if (__debug__): print("Jac start")
	gmap = geo_map(n, s, p, x_name, c_name, t_name)
	J =	Matrix([[expand(poly).diff(v) for v in xt] for poly in gmap])
	if (__debug__): print("Jac end")

	pts = domain_pts_J(n, s, p, dynamic)

	with ProcessPoolExecutor() as executor:
		futures = [
			executor.submit(aux_evaluate_det, (pt, J, xt, n)) for pt in pts
		]
		lag_vec = [future.result() for future in futures]

	return lag_vec



## Time subdivision maps ##
def time_subdiv_map(pt, shift):
	h = 1 if shift else 0
	return(pt[:-1] + ((pt[-1] + h) * QQ(1,2),))

## Space subdivision maps ##
def space_subdiv_map(pt, n, s, q):
	assert(1 <= s <= n)
	assert(n <= len(pt) <= n+1)
	assert(0 <= q < 2**len(pt))
	res = list(pt)
	qs = q % (2**s)
	# Simplex part
	if s == 1:
		res[0] = pt[0] + qs
	elif s == 2:
		if qs == 0:
			res[0] = pt[0]
			res[1] = pt[1]
		elif qs == 1:
			res[0] = pt[0] + 1
			res[1] = pt[1]
		elif qs == 2:
			res[0] = pt[0]
			res[1] = pt[1] + 1
		elif qs == 3:
			res[0] = 1 - pt[0]
			res[1] = 1 - pt[1]
	elif s == 3:
		if qs == 0:
			res[0] = pt[0]
			res[1] = pt[1]
			res[2] = pt[2]
		elif qs == 1:
			res[0] = pt[0] + 1
			res[1] = pt[1]
			res[2] = pt[2]
		elif qs == 2:
			res[0] = pt[0]
			res[1] = pt[1] + 1
			res[2] = pt[2]
		elif qs == 3:
			res[0] = pt[0]
			res[1] = pt[1]
			res[2] = pt[2] + 1
		elif qs == 4:
			res[0] = 1 - pt[1] - pt[2]
			res[1] = pt[1]
			res[2] = pt[0] + pt[1] + pt[2]
		elif qs == 5:
			res[0] = 1 - pt[1]
			res[1] = pt[0] + pt[1]
			res[2] = pt[1] + pt[2]
		elif qs == 6:
			res[0] = pt[0] + pt[1]
			res[1] = 1 - pt[0]
			res[2] = pt[2]
		elif qs == 7:
			res[0] = pt[0]
			res[1] = pt[1] + pt[2]
			res[2] = 1 - pt[0] - pt[1]
	else:
		raise Exception('Unsupported simplex dimension')

	# Tensor part
	for k in range(s, len(pt)):
		qt = (q >> k) & 1
		res[k] = pt[k] + qt

	# Divide and return
	return tuple(r * QQ(1,2) for r in res)

def time_subdiv_maps(pts, shift):
	return [time_subdiv_map(pt, shift) for pt in pts]
def space_subdiv_maps(pts, n, s, q):
	return [space_subdiv_map(pt, n, s, q) for pt in pts]

## Compress matrix ##
def mat_compress(m):
	l = m.rows
	assert(m.cols == l)
	s = str(l)
	for i in range(l):
		for j in range(l):
			if m[i, j] != 0:
				s += f', {i}, {j}, {m[i,j].p}, {m[i,j].q}'
	return s

## Binomial and multinomial coefficients ##
def multinomial(ii):
	res = 1
	k = sum(ii)
	i0 = ii.index(max(ii))
	for a in ii[:i0] + ii[i0+1:]:
		for j in range(1,a+1):
			res *= k
			res //= j
			k -= 1
	return res

## Subdivision matrices ##
def eval_basis(args):
	xt, basis, pts = args
	return DomainMatrix.from_Matrix(Matrix(
	[[b.eval({xt[k]: pt[k] for k in range(len(xt))})
	for b in basis] for pt in pts])).to_field()
def aux_matmul(args): return args[0] * args[1]

def subdiv_matrices(n, s, p, dynamic):
	simplex_ord = n * p - s
	tensor_ord = n * p - 1
	time_ord = n
	xt = subscripts('u', range(n))
	if dynamic: xt += [symbols('t')]
	xt_full = xt + [1 - sum(xt[:s])] + [1 - k for k in xt[s:]]
	exponents = index_set_J(n, s, p, dynamic)
	if dynamic: 
		exponents_full = [list(e) + [simplex_ord - sum(e[:s])] + \
			[tensor_ord - e[k] for k in range(s, n)] + [time_ord - e[n]]
			for e in exponents]
	else: 
		exponents_full = [list(e) + [simplex_ord - sum(e[:s])] + \
			[tensor_ord - e[k] for k in range(s, n)]
			for e in exponents]

	r_full = range(len(xt_full))
	if dynamic: 
		basis = [Poly(
				multinomial(list(e[:s]) + [simplex_ord - sum(e[:s])]) *
				prod(multinomial((ee, tensor_ord - ee)) for ee in e[s:n]) *
				multinomial((e[n], time_ord - e[n])) *
				prod([xt_full[k]**exponents_full[i][k] for k in r_full]),
			xt, domain='QQ') for i,e in enumerate(exponents)]
	else: 
		basis = [Poly(
				multinomial(list(e[:s]) + [simplex_ord - sum(e[:s])]) *
				prod(multinomial((ee, tensor_ord - ee)) for ee in e[s:n]) *
				prod([xt_full[k]**exponents_full[i][k] for k in r_full]),
			xt, domain='QQ') for i,e in enumerate(exponents)]
	assert(sum(basis) == 1)
	if (__debug__): print('Basis computed')
	pts = domain_pts_J(n, s, p, dynamic)

	with ProcessPoolExecutor() as executor:
		if dynamic: 
			args_combined = [(xt, basis, pts)] + \
				[(xt, basis, time_subdiv_maps(pts, t)) for t in (False, True)] + \
				[(xt, basis, space_subdiv_maps(pts, n, s, q)) for q in range(2**(n+1))]
		else: 
			args_combined = [(xt, basis, pts)] + \
				[(xt, basis, space_subdiv_maps(pts, n, s, q)) for q in range(2**n)]

		futures_combined = [executor.submit(eval_basis, args)
			for args in args_combined]
		results_combined = [future.result() for future in futures_combined]
	if (__debug__): print('Prepared subdivision matrices')

	l2b = results_combined[0].inv()
	if dynamic: 
		tsd = [l2b * m for m in results_combined[1:3]]
		ssd = [l2b * m for m in results_combined[3:]]
		if (__debug__): print('Computed subdivision matrices')
		return [
			l2b.to_Matrix(),
			[m.to_Matrix() for m in tsd],
			[m.to_Matrix() for m in ssd]
		]
	else:
		ssd = [l2b * m for m in results_combined[1:]]
		if (__debug__): print('Computed subdivision matrices')
		return [
			l2b.to_Matrix(),
			[m.to_Matrix() for m in ssd]
		]

## Format matrices ##
def matrices_formatted(n, s, p, dynamic):
	lines = []
	if dynamic:
		lines += ['template<>\n' +
			f'void initMatricesT<{n}, {s}, {p}> ' +
			'(Matrix<Interval> &l2b, ' +
			'span<Matrix<Interval>> tsd, ' +
			'span<Matrix<Interval>> ssd) ' +
			'{']
	else:
		lines += ['template<>\n' +
			f'void initMatrices<{n}, {s}, {p}> ' +
			'(Matrix<Interval> &l2b, ' +
			f'span<Matrix<Interval>> ssd) ' +
			'{']
	
	mm = subdiv_matrices(n, s, p, dynamic)
	lines.append('l2b.fill({' + mat_compress(mm[0]) + '});')
	if dynamic:
		lines.append('tsd[0].fill({' + mat_compress(mm[1][0]) + '});')
		lines.append('tsd[1].fill({' + mat_compress(mm[1][1]) + '});')
	for q,a in enumerate(mm[-1]):
		lines.append(f'ssd[{q}].fill(' + '{' + mat_compress(a) + '});')
	return '\n\t'.join(lines) + '\n}'

def unique_rational_coefficients(polynomials):
	unique_coeffs = set()
	
	for poly in polynomials:
		coeffs = Poly(poly, domain='QQ').coeffs()
		for coeff in coeffs:
			if coeff.is_rational and coeff.denominator != 1:
				unique_coeffs.add(abs(coeff))
	
	return list(unique_coeffs)

## Subdivide a list in chunks
def chunkize(data, target_size):
	assert(target_size > 0)
	L = len(data)
	num_chunks = (L-1) // target_size + 1
	chunk_size = L // num_chunks
	leftovers = L - (num_chunks * chunk_size)
	assert leftovers < num_chunks
	chunks = []
	for i in range(leftovers):
		chunks.append(data[i*(chunk_size+1) : (i+1)*(chunk_size+1)])
	offset = leftovers * (chunk_size+1)
	for i in range(num_chunks-leftovers):
		chunks.append(data[offset + i*chunk_size : offset + (i+1)*chunk_size])
	return chunks

## Format lag vector ##
def lag_vec_formatted(file, n, s, p, dynamic):
	my_ccode = MyCodePrinter().doprint
	func_name = 'lagrangeVector' + ('T' if dynamic else '')
	cfunc_name = 'chunk' + ('T' if dynamic else '')

	NCP = len(index_set(n,s,p)) * n
	if dynamic: NCP *= 2

	expr = lagrange_vector(n, s, p, 'cp', dynamic)
	chunk_size = len(expr)
	if (n,s,p) in EXPRESSION_CHUNK_SIZE:
		chunk_size = EXPRESSION_CHUNK_SIZE[(n,s,p)]
	chunks = chunkize(expr, chunk_size)

	f.write('#include "../lagrangeVector.hpp"\n\n')
	f.write('#define R(p, q) (Interval(p) / q)\n')
	f.write('#define I const Interval \n\n')

	if (len(chunks) > 1):
		if (__debug__): print(f'Splitting vector in {len(chunks)} chunks')
		if (__debug__): print(f'Writing chunk function imports')
		file.write(f'#include "{func_name}_{n}_{s}_{p}_chunks.hpp"\n')
		file.write('\n')

	f.write('namespace element_validity {\n')


	file.write(
		'template<>\n' +
		f'void {func_name}<{n}, {s}, {p}>' +
		'(const span<const fp_t> cpFP, const span<Interval> out) {'
	)
	file.write(f'\n\tassert(cpFP.size() == {NCP});')
	file.write(f'\n\tassert(out.size() == {len(expr)});')
	file.write(f'\n\tstd::array<Interval, {NCP}> cp;')
	file.write(f'\n\tfor (int i = 0; i < {NCP}; ++i) cp[i] = cpFP[i];')
	if TO_ORIGIN:
		file.write(f'\n\tfor (int i = S-1; i >= 0; --i) cp[i] -= cp[i%{n}];')

	if (len(chunks) == 1):
		if (__debug__): print('CSE start')
		temp_expr = cse(expr, numbered_symbols('tmp_'), order='canonical')
		if (__debug__): print('CSE end')
		assert len(expr) == len(temp_expr[1])

		# if (__debug__): print('Unique coefficients listing start')
		# temp_coef = unique_rational_coefficients(expr)
		# if (__debug__): print('Unique coefficients listing end')

		if (__debug__): print('Writing temporaries')
		# for qq in temp_coef:
		# 	num, den = qq.as_numer_denom()
		# 	file.write(f'\n\tI R_{num}_{den} = R({num},{den});')
		for helper in temp_expr[0]:
			file.write('\n\tI ' + my_ccode(helper[1], helper[0]))
		if (__debug__): print('Writing expressions')
		for c,result in enumerate(temp_expr[1]):
			file.write('\n\t' + my_ccode(result, f'out[{c}]'))
		if (__debug__): print('Finished writing')
	else:
		for c in range(len(chunks)):
			file.write(f'\n\t{cfunc_name}_{n}_{s}_{p}_{c}(cp, out);')

		# Create new files to write chunks
		chunk_file_directory = os.path.dirname(f.name)
		chunk_file_base_name, chunk_file_extension = \
			os.path.splitext(os.path.basename(f.name))
		
		chunks_head_name = os.path.join(chunk_file_directory,
			f"{chunk_file_base_name}_chunks.hpp")

		with open(chunks_head_name, 'w') as cfile:
			cfile.write(
				'#pragma once\n' +
				'#include "numbers/Interval.hpp"\n' +
				'#include <array>\n' +
				'namespace element_validity {\n'
			)
			for c,chunk in enumerate(chunks):
				cfile.write(f'void {cfunc_name}_{n}_{s}_{p}_{c}' +
					f'(const std::array<Interval, {NCP}> &cp, ' + 
					'const span<Interval> &out);\n')
			cfile.write('}\n')
			

		offset = 0
		for c,chunk in enumerate(chunks):
			if (__debug__): print(f'Processing chunk {c}')
			chunk_file_name = os.path.join(chunk_file_directory,
				f"{chunk_file_base_name}_chunk{c}{chunk_file_extension}")
			with open(chunk_file_name, 'w') as cfile:
				cfile.write(
					f'#include "{chunk_file_base_name}_chunks.hpp"\n' +
					'#define R(p, q) (Interval(p) / q)\n' +
					'#define I const Interval \n\n' +
					'namespace element_validity {\n' +
					f'void {cfunc_name}_{n}_{s}_{p}_{c}' +
					f'(const std::array<Interval, {NCP}> &cp, ' + 
					'const span<Interval> &out) {'
				)

				if (__debug__): print('CSE start')
				temp_expr = cse(chunk, numbered_symbols('tmp_'), order='canonical')
				if (__debug__): print('CSE end')
				len(chunk) == len(temp_expr[1])

				# if (__debug__): print('Unique coefficients listing start')
				# temp_coef = unique_rational_coefficients(expr)
				# if (__debug__): print('Unique coefficients listing end')

				if (__debug__): print('Writing temporaries')
				# for qq in temp_coef:
				# 	num, den = qq.as_numer_denom()
				# 	cfile.write(f'\n\tI R_{num}_{den} = R({num},{den});')
				for helper in temp_expr[0]:
					cfile.write('\n\tI ' + my_ccode(helper[1], helper[0]))
				if (__debug__): print('Writing expressions')
				for i,result in enumerate(temp_expr[1]):
					cfile.write('\n\t' + my_ccode(result, f'out[{offset+i}]'))
				cfile.write('\n}\n}\n')
				cfile.write('#undef R\n')
				cfile.write('#undef I')
			if (__debug__): print(f'Finished writing chunk {c} to {chunk_file_name}')
			offset += len(chunk)
	f.write('\n}\n}\n')
	f.write('#undef R\n')
	f.write('#undef I')
			


## Format lag basis ##
def lag_eval_formatted(file, n, s, p, dynamic):
	func_name = 'lagrangeEvaluate' + ('T' if dynamic else '')
	indices = index_set_J(n, s, p, dynamic)
	lag = [lagrange_J(n, s, p, k, 'x') for k in indices]
	if dynamic:
		tlag = [lagrange_uni(n, n, i, symbols(f'x[{n}]')) for i in range(n+1)]
		expr = [pa * pb for pa in lag for pb in tlag]
	else: expr = lag
	if (__debug__): print('CSE start')
	temp_expr = cse(expr, numbered_symbols('tmp_'), order='canonical')
	if (__debug__): print('CSE end')

	# if (__debug__): print('Unique coefficients listing start')
	# temp_coef = unique_rational_coefficients(expr)
	# if (__debug__): print('Unique coefficients listing end')

	file.write(
		'template<>\n' +
		f'Interval {func_name}<{n}, {s}, {p}>(\n' +
		'\tconst span<const fp_t> xFP,\n' +
		'\tconst span<const Interval> lagVec\n' + 
		') {'
	)
	nx = n
	if dynamic: nx += 1
	file.write(f'\n\tassert(xFP.size() == {nx});')
	file.write(f'\n\tassert(lagVec.size() == {len(expr)});')
	assert len(expr) == len(temp_expr[1])
	file.write(f'\n\tstd::array<Interval, {nx}> x;')
	file.write(f'\n\tfor (int i = 0; i < {nx}; ++i) x[i] = xFP[i];')
	file.write(f'\n\tInterval acc = 0.;')
	my_ccode = MyCodePrinter().doprint
	if (__debug__): print('Writing temporaries')
	# for qq in temp_coef:
	# 	num, den = qq.as_numer_denom()
	# 	file.write(f'\n\tI R_{num}_{den} = R({num},{den});')
	for helper in temp_expr[0]:
		file.write('\n\tI ' + my_ccode(helper[1], helper[0]))
	if (__debug__): print('Writing expressions')
	for i,result in enumerate(temp_expr[1]):
		file.write(f'\n\tacc += lagVec[{i}] * ({my_ccode(result)});')
	file.write(f'\n\treturn acc;')
	if (__debug__): print('Finished writing')
	file.write('\n}')

## Format corner indices
def corners_formatted(n, s, p, dynamic):
	func_name = 'cornerIndices' + ('T' if dynamic else '')
	ind = map(str,corner_indices_set(n, s, p, dynamic))
	return f'template<>\nvoid {func_name}<{n}, {s}, {p}>(std::vector<int> &v) ' + \
		'{ v = {' + ','.join(ind) + '}; }'

## Write ##
# if WRITE_CMAKE and not DRY_RUN:
# 	comb = [
# 		('lagvec', 'lagrangeVector'),
# 		('transmat', 'transMatrices'),
# 		('lagvecT', 'lagrangeVectorT'),
# 		('transmatT', 'transMatricesT'),
# 		('lageval', 'lagrangeEvaluate'),
# 		# ('lagevalT', 'lagrangeEvaluateT'),
# 	]
# 	for cc in comb:
# 		with open(f'src/validity/{cc[0]}/CMakeLists.txt', 'w') as f:
# 			f.write('set(SOURCES\n')
# 			for n,s,p in COMBINATIONS:
# 				f.write(f'\t{cc[1]}_{n}_{s}_{p}.cpp\n')
# 			f.write(')\n\n')

# 			f.write('source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})\n')
# 			f.write('target_sources(bezier PRIVATE ${SOURCES})')

if WRITE_MATRICES:
	print('Writing matrices...')
	path = lambda n,s,p,d: f'src/validity/transmat{"T" if d else ""}/transMatrices{"T" if d else ""}_{n}_{s}_{p}.cpp'
	if DRY_RUN: path = lambda n,s,p,d: '/dev/null'
	for n,s,p in COMBINATIONS:
		if DYNAMIC:
			with open(path(n,s,p,True), 'w') as f:
				f.write('#include "../transMatrices.hpp"\n\n')
				f.write('namespace element_validity {\n')
				if (__debug__): print(f'({n} {s} {p})')
				f.write(matrices_formatted(n, s, p, True))
				f.write('\n}\n')
		if STATIC:
			with open(path(n,s,p,False), 'w') as f:
				f.write('#include "../transMatrices.hpp"\n\n')
				f.write('namespace element_validity {\n')
				if (__debug__): print(f'({n} {s} {p})')
				f.write(matrices_formatted(n, s, p, False))
				f.write('\n}\n')
	print()
	print('Done writing matrices.')
else:
	print('Not writing matrices.')

if WRITE_LAGVEC:
	print('Writing Lagrange vectors...')
	path = lambda n,s,p,d: f'src/validity/lagvec{"T" if d else ""}/lagrangeVector{"T" if d else ""}_{n}_{s}_{p}.cpp'
	if DRY_RUN: path = lambda n,s,p,d: '/dev/null'
	for n,s,p in COMBINATIONS:
		if DYNAMIC:
			with open(path(n,s,p,True), 'w') as f:
				lag_vec_formatted(f,n,s,p,True)
		if STATIC:
			with open(path(n,s,p,False), 'w') as f:
				lag_vec_formatted(f,n,s,p,False)
	print('Done writing Lagrange vectors.')
else:
	print('Not writing Lagrange vectors.')

if WRITE_LAGPOLY:
	print('Writing Lagrange basis polynomials...')
	path = lambda n,s,p,d: f'src/validity/lageval{"T" if d else ""}/lagrangeEvaluate{"T" if d else ""}_{n}_{s}_{p}.cpp'
	if DRY_RUN: path = lambda n,s,p,d: '/dev/null'
	for n,s,p in COMBINATIONS:
		# if DYNAMIC:
		# 	with open(path(n,s,p,True), 'w') as f:
		# 		f.write('#include "../lagrangeEvaluate.hpp"\n\n')
		# 		f.write('#define R(p, q) (Interval(p) / q)\n')
		# 		f.write('#define I const Interval \n\n')
		# 		f.write('namespace element_validity {\n')
		# 		lag_eval_formatted(f,n,s,p,True)
		# 		f.write('}\n')
		# 		f.write('#undef R\n')
		# 		f.write('#undef I')
		if STATIC:
			with open(path(n,s,p,False), 'w') as f:
				f.write('#include "../lagrangeEvaluate.hpp"\n\n')
				f.write('#define R(p, q) (Interval(p) / q)\n')
				f.write('#define I const Interval \n\n')
				f.write('namespace element_validity {\n')
				lag_eval_formatted(f,n,s,p,False)
				f.write('}\n')
				f.write('#undef R\n')
				f.write('#undef I')
	print('Done writing Lagrange basis.')
else:
	print('Not writing Lagrange basis.')

if WRITE_CORNERS:
	print('Writing corner indices...')
	path = 'src/validity/cornerIndices.cpp'
	if DRY_RUN: path = '/dev/null'
	with open(path, 'w') as f:
		f.write('#include "cornerIndices.hpp"\n\n')
		f.write('namespace element_validity {\n')

		for n,s,p in COMBINATIONS:
			f.write(corners_formatted(n,s,p,False))
			f.write('\n\n')
			f.write(corners_formatted(n,s,p,True))
			f.write('\n\n')

		f.write('\n}\n')
	print('Done writing corner indices.')
else:
	print('Not writing corner indices.')

if INFO_ORDER:
	for n,s,p in COMBINATIONS:
		print(f'This is the order of geometric map control points for elements of type ({n},{s},{p}):')
		ind = index_set(n, s, p)
		for i,t in enumerate(ind):
			print(f'\t{i}: {t}')
		print()

if INFO_JAC_ORDER:
	for n,s,p in COMBINATIONS:
		print(f'This is the order of jacobian control points for elements of type ({n},{s},{p}):')
		ind = domain_pts_J(n, s, p, False)
		for i,t in enumerate(ind):
			print(f'\t{i}: {t}')
		print(f'This is the order of dynamic jacobian control points for elements of type ({n},{s},{p}):')
		indT = domain_pts_J(n, s, p, True)
		for i,t in enumerate(indT):
			print(f'\t{i}: {t}')
		print()

if INFO_LAGBASIS:
	for n,s,p in COMBINATIONS:
		print(f'This is the Lagrange basis for elements of type ({n},{s},{p}):')
		indices = index_set(n, s, p)
		lag = [lagrange(n, s, p, k, 'x') for k in indices]
		for i,poly in enumerate(lag):
			print(f'\t{i}: {poly}')
		print()

if INFO_LAGVECTOR:
	for n,s,p in COMBINATIONS:
		print(f'This is the Lagrange vector for elements of type ({n},{s},{p}):')
		lag = lagrange_vector(n, s, p, 'c', False)
		for i,poly in enumerate(lag):
			print(f'\t{i}: {poly}')
		print()
