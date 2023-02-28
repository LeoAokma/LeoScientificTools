import numpy as np
import time
mass = {
	'H' : 1.007553, 'Li' : 7.016283, 'C' : 12, 'N' : 14.00308, 'O' : 15.99491, 'Na' : 22.99004,
	'Cl-' : 34.96912, 'Cl+' : 36.96672,	'K': 38.96397217,
}
atoms = ['C', 'H', 'N', 'O', 'Cl', 'Na', 'K']
# print('C', 'H', 'N', 'O', 'Cl35', 'Na', 'Error', sep='\t')

def analyse(ms_weight):
	result = []
	start_time = time.time()
	# H, C, N, O, Na, Cl35, Cl37
	for c_number in range(0, round(ms_weight / mass['C']) + 1):
		res_c = ms_weight - c_number * mass['C']
		for n_number in range(0, round(res_c / mass['N']) + 1):
			res_n = res_c - n_number * mass['N']
			for o_number in range(0, round(res_n / mass['O']) + 1):
				res_o = res_n - o_number * mass['O']
				for cl35 in range(0, round(res_o / mass['Cl-']) + 1):
					res_cl = res_o - cl35 * mass['Cl-']
					for na_number in range(0, 2):
						res_k = res_cl - na_number * mass['Na']
						for k_number in range(0, 2):
							res = res_k - k_number * mass['K']
							if res > 0 :
								h_number = round(res / mass['H'])
								res = res - h_number * mass['H']
							error = res
							if abs(error) < 0.2:
								omega = c_number + 1 - (h_number + cl35 - n_number)/2
								# saturation constrain
								if -2 < omega < 8:
									if k_number + na_number < 2:
										result.append([c_number, h_number, n_number, o_number, cl35, na_number, k_number, error])
									#print("{}\t{}\t{}\t{}\t{}\t{}\t{:.3%}".format(c_number, h_number, n_number, o_number, cl35, na_number, LgE))
	end_time = time.time()
	duration = end_time - start_time
	print('Calculation Complete, found {} match, in {:.2f} s'.format(len(result), duration))
	return result, duration

def constrain(req, result):
	cfm = req.strip().split(' ')
	cfm_res = result[:]
	for _ in range(5):
		if cfm[_].upper() != 'N':
			if cfm[_][0] == '=':
				cfm_res = [x for x in cfm_res if x[_] == int(cfm[_][1:]) ]
			elif cfm[_][0] == '<':
				cfm_res = [x for x in cfm_res if x[_] < int(cfm[_][1:]) ]
			elif cfm[_][0] == '>':
				cfm_res = [x for x in cfm_res if x[_] > int(cfm[_][1:]) ]
	
	sorted_res = sorted(cfm_res, key = lambda x: abs(x[-1]), reverse=False)
	print('......Constrain Index Complete......')
	print('Found {} result(s).'.format(len(cfm_res)))
	print('匹配的化学式\t残余化学式\t绝对误差')
	for _ in sorted_res:
		# print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3f}".format(*_))
		formula = ''
		for ind in range(len(_) - 1):
			if _[ind] >0:
				formula += atoms[ind]
				if _[ind] >1:
					formula += str(_[ind])
		print(formula, end='\t')
		# 计算其残余分子式
		cal_for = np.array(_[:-1])
		augment = cal_for - np.array(tar_for)
		plus = '+'
		minus = '-'
		for num in range(len(augment)):
			if augment[num] >0:
				plus += atoms[num]
				if augment[num] >1:
					plus += str(augment[num])
			elif augment[num] <0:
				minus += atoms[num]
				minus += str(abs(augment[num]))
		if len(minus) == 1:
			minus = ""
		if len(plus) == 1:
			plus = ""	
		print("{} {}\t{:.3f}".format(plus, minus, _[-1]))
	return sorted_res


ms_signals = input('Please enter the ms signal(m/z)\nDivided by space.\n').strip().split(' ')
signals = [float(x) for x in ms_signals]

cfm_input = input('Input the general constrain condition(C H N O Cl) split by space, N for no constrain, use >, < and = before numbers.\n')

tar_formula = input('请输入目标化学式，用空格分割(C, H, N, O, Cl)\n').strip().split(' ')
tar_for = [int(x) for x in tar_formula]
tar_for.append(0)
tar_for.append(0)

for signal in signals:
	print('-------正在分析信号{}-------'.format(signal))
	res, dur = analyse(signal)
	out = constrain(cfm_input, res)
	
print('全部任务已完成')
