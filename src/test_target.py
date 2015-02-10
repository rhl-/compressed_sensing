import unittest
import common
import numpy as np
import cmath
import gen_target as target

class target_test(unittest.TestCase):
	def assertAlmostEqualDigits(self,a,b):
		self.assertAlmostEqual(a,b,self.digits)

	def get_sample(self,t=common.TARGET_TYPES.HPSD):
		return target.getTargetSample(self.n, self.r, t)

	def setUp(self):
		self.digits=13
		self.n = 100
		self.r = 25

	def test_RPSD(self):
		A = self.get_sample(common.TARGET_TYPES.RPSD)
		self.assertAlmostEqualDigits(np.linalg.norm(A - A.T), 0)
		evals = np.linalg.eigvalsh(A)
		ell = np.amin(evals)
		# test that eigenvalues are bigger than zero
		self.assertAlmostEqualDigits(ell, 0)
		# test that nonzero eigenvalues are close to 1
		self.assertAlmostEqualDigits(np.amax(evals-1),0)

	def test_RSYM(self):
		A = self.get_sample(common.TARGET_TYPES.RSYM)
		self.assertAlmostEqualDigits(np.linalg.norm(A - A.T), 0)
		evals = np.linalg.eigvalsh(A)
		evals = abs(evals)
		max_diff_from_1= np.amax(abs(evals[np.where(evals) > 1e-3] - 1))
		self.assertAlmostEqualDigits(max_diff_from_1, 0)

	def test_HPSD(self):
		A = self.get_sample(common.TARGET_TYPES.HPSD)
		self.assertAlmostEqualDigits(np.linalg.norm(A - A.conj().transpose()), 0)
		evals = np.linalg.eigvalsh(A)
		ell = np.amin(evals)
		# test that eigenvalues are bigger than zero
		self.assertAlmostEqualDigits(ell, 0)
		# test that nonzero eigenvalues are close to 1
		self.assertAlmostEqualDigits(np.amax(evals-1),0)

	def test_HERM(self):
		A = self.get_sample(common.TARGET_TYPES.HERM)
		self.assertAlmostEqualDigits(np.linalg.norm(A - A.conj().transpose()), 0)
		evals = np.linalg.eigvalsh(A)
		evals = abs(evals)
		max_diff_from_1= np.amax(abs(evals[np.where(evals) > 1e-3] - 1))
		self.assertAlmostEqualDigits(max_diff_from_1, 0)






if __name__ == '__main__':
    unittest.main()
