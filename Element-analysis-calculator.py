element = {'C':12.011, 'H':1.00797, 'N': 14.007, 'O':15.999, 'S':32.06,
					'Eu':151.964, 'I':126.90447, 'Br':79.904, 'Cl':35.45, 'F':18.998403163,
					'P':30.973761998, 'La':138.90547, 'Ce':140.116
						}

def formula_check(form):
	'''
	input string
	'''
	figs = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0']
	atoms = []
	atom_nums = []
	num = ''
	atom = ''
	atom_cal = 0
	for _ in form:
		if _ in figs:
			num += _
		elif _ == ' ':
			if num == '':
				atom_nums.append(1)
			else:
				atom_nums.append(int(num))
			if atom != '':
				atoms.append(atom)
		elif _.isupper():
			if atom_cal == 0:
				atom += _
				atom_cal =+ 1
			elif atom_cal != 0:
				if atom != '':
					atoms.append(atom.rstrip())
				if num != '':
					atom_nums.append(int(num))
					num = ''
				else:
					atom_nums.append(1)
				atom = _
				if atom_cal == 2:
					atom_cal = 0
		else:
			atom_cal += 1
			atom += _
			atoms.append(atom.rstrip())
			atom=''
			'''
			if num != '':
				atom_nums.append(int(num))
				num = ''
			else:
				atom_nums.append(1)
			'''
	return atoms, atom_nums

def cal_percentage(formula_atom, formula_num, report=True, atom_list=['C', 'N', 'H', 'S']):
	'''
	:param report: Boolean, choose whether to report the results or not
	'''
	mw = 0
	for _ in range(len(formula_atom)):
		mw += element[formula_atom[_]] * formula_num[_]
	
	results = {}
	for atom in atom_list:
		if atom in formula_atom:
			percentage = element[atom] * formula_num[formula_atom.index(atom)] / mw
			results[atom] = percentage
		else:
			results[atom] = 0
	if report:
		print('Type\t', end='')
		for _ in results:
			print('{}%'.format(_).ljust(6), end='\t')
		print('\n', end='')
		print('Value\t', end='')
		for _ in results:
			print('{:.2%}'.format(results[_]), end='\t')
		print('\n', end='')
	return results

formula_input = input('Please enter the chemical formula(e.g. C2H6O)\nNo seperations between atoms and numbers are needed, and mind the capital alphabets!\n')
formula = formula_input.strip() + ' '
print('Successfully read formula as:', end='\n')
form_read = formula_check(formula)
print(form_read)
cal_percentage(*form_read)
