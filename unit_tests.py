from model import *

# NOTE: run with pytest unit_tests.py

def test_ode_model():
	"""
	Test if function ode_model is working properly by comparing it with a known result.
	"""
	ode_test1 = ode_model(10,1,-1,1,1,0)
	ode_test2 = ode_model(10,2.7,-1,1,1,0)
	ode_test3 = ode_model(10,-1000,-1,1,1,0)
	ode_test4 = ode_model(0,0,0,0,0,0)
	ode_test5 = ode_model(-1,-1,-1,-1,-1,-1)
	ode_test6 = ode_model(1000,1000,1000,1000,1000,1000)

	ode_sln1 = 0
	ode_sln2 = -1.7
	ode_sln3 = 1001
	ode_sln4 = 0
	ode_sln5 = -1
	ode_sln6 = -1000000

	assert ode_test1 == ode_sln1
	assert round(ode_test2,1) == ode_sln2
	assert ode_test3 == ode_sln3
	assert ode_test4 == ode_sln4
	assert ode_test5 == ode_sln5
	assert ode_test6 == ode_sln6

def test_solve_ode():
	"""
	Test if function solve_ode is working properly by comparing it with a known result.
	"""
	t_test,solve_ode_test = solve_ode(ode_model,0,10,0.1,0,[-1,-1,-1,0])

	t_sln = np.arange(0,10.1,0.1)
	ode_sln = 1-np.e**(-1*t_sln)

	for i in range(10):
		assert (solve_ode_test[i] - ode_sln[i]) < 1.e-8
