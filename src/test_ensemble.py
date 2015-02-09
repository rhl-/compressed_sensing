import unittest
import common
import numpy as np
import cmath
import gen_ensemble as ensemble

class ensemble_test(unittest.TestCase):
	def assertAlmostEqualDigits(self,a,b):
		self.assertAlmostEqual(a,b,self.digits)

	def is_perm_format(self,A):
		nnz = np.count_nonzero( A)
 		self.assertEquals( nnz, self.m*(self.n))
		for matrix in A:
			self.assertEquals( np.count_nonzero( matrix), self.n)
			for row in np.reshape(matrix, (self.n,self.n)):
				self.assertEquals( np.count_nonzero( row), 1)
			for column in np.reshape(matrix, (self.n,self.n)).T:
				self.assertEquals( np.count_nonzero( column), 1)

	def get_sample(self,t1, t2=common.TARGET_TYPES.HPSD):
		return ensemble.getEnsembleSample(self.m, self.n, t1, t2)

	def setUp(self):
		self.digits=4
		self.m = 128
		self.n = 128

	def test_entry(self):
		A = self.get_sample(common.ENSEMBLE_TYPES.ENTRY)
		nnz = np.count_nonzero(A)
 		self.assertEquals(nnz, self.m)
		for matrix in A:
			self.assertEquals(np.count_nonzero(matrix), 1)
		
	def test_perm(self):
		A = self.get_sample(common.ENSEMBLE_TYPES.PERM)
		self.is_perm_format( A)

	def test_rsperm(self):
		A = self.get_sample(common.ENSEMBLE_TYPES.RSPERM)
		self.is_perm_format( A)
		B = self.n*np.multiply(A,A);	
		self.assertAlmostEqualDigits( sum(sum(B)), np.count_nonzero( B))

	def test_csperm(self):
		A = self.get_sample(common.ENSEMBLE_TYPES.CSPERM)
		self.is_perm_format( A)
		B = np.multiply(A,A)
		B = np.multiply(B,B)
		B = (self.n**2)*B
		self.assertAlmostEqualDigits( sum(sum(B)), np.count_nonzero( B))

	def test_rgperm(self):
		A = self.get_sample(common.ENSEMBLE_TYPES.RGPERM)
		self.is_perm_format( A)
		B = A[ np.nonzero( A)].flatten()
		self.assertAlmostEqual( np.mean(B), 0.0, 2)
		self.assertAlmostEqual( np.mean( np.multiply(B,B).flatten()), 1.0/self.n, 2)

	def test_cgperm(self):
		A = self.get_sample(common.ENSEMBLE_TYPES.CGPERM)
		self.is_perm_format( A)
		B = np.multiply(A,A)
		B = np.multiply(B,B)
		B = (self.n**2)*B
		B = A[ np.nonzero( A)].flatten()
		self.assertAlmostEqual( np.mean( B), 0.0, 2)
		self.assertAlmostEqual( np.mean( np.multiply( abs(B), abs(B))), 1.0/(self.n**2), 1)

	def test_rdirac(self):
		A = self.get_sample(common.ENSEMBLE_TYPES.RDIRAC)
		for row in A:
			self.assertAlmostEqualDigits(np.linalg.norm( row), 1)

	def test_cdirac(self):
		A = self.get_sample(common.ENSEMBLE_TYPES.CDIRAC)
		for row in A:
			self.assertAlmostEqualDigits(np.linalg.norm( row),1)
	
	def test_rgauss(self):
		A = self.get_sample(common.ENSEMBLE_TYPES.RGAUSS)
		B = A.flatten();
		self.assertAlmostEqualDigits( np.mean( B), 0)
		self.assertAlmostEqualDigits( np.mean( np.multiply( abs(B), abs(B))), 1.0/(self.n**2))
	
	def test_cgauss(self):
		A = self.get_sample(common.ENSEMBLE_TYPES.CGAUSS)
		B = A.flatten();
		self.assertAlmostEqualDigits( np.mean( B), 0)
		self.assertAlmostEqualDigits( np.mean( np.multiply( abs(B), abs(B))), 1.0/(self.n**2))

if __name__ == '__main__':
    unittest.main()
