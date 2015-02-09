import unittest
import common
import numpy as np
import cmath
import gen_ensemble as ensemble


class ensemble_test(unittest.TestCase):
	def is_perm_format(self,A):
		nnz = np.count_nonzero( A)
 		self.assertEquals( nnz, self.m*self.n*self.n)
		for matrix in A:
			self.assertTrue( np.count_nonzero( matrix) == self.n**2)
			for row in np.reshape(matrix, (self.n,self.n)):
				self.assertEquals( np.count_nonzero( row), 1)
			for column in np.reshape(matrix, (self.n,self.n)).T:
				self.assertEquals( np.count_nonzero( column), 1)

	def get_sample(self,t1, t2=common.TARGET_TYPES.HPSD):
		return ensemble.getEnsembleSample(self.m, self.n, t1, t2)

	def setUp(self):
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
		B = n*np.multiply(A,A);	
		self.assertAlmostEqual( sum(sum(B)), np.count_nonzero( B))

	def test_csperm(self):
		A = self.get_sample(common.ENSEMBLE_TYPES.CSPERM)
		self.is_perm_format( A)
		B = np.multiply(A,A)
		B = np.multiply(B,B)
		B = (n**2)*B
		self.assertAlmostEqual( sum(sum(B)), np.count_nonzero( B))

	def test_rgperm(self):
		A = self.get_sample(common.ENSEMBLE_TYPES.RGPERM)
		self.is_perm_format( A)
		B = n*np.multiply(A,A)
		self.assertAlmostEqual( sum(sum(B)), np.count_nonzero( B))

	def test_cgperm(self):
		A = self.get_sample(common.ENSEMBLE_TYPES.CGPERM)
		self.is_perm_format( A)
		B = np.multiply(A,A)
		B = np.multiply(B,B)
		B = (n**2)*B
		self.assertAlmostEqual( sum(sum(B)), np.count_nonzero( B))
		B = np.nonzeroflatten( A);
		self.assertAlmostEqual( np.mean( B), 0)
		self.assertAlmostEqual( np.mean( np.multiply( abs(B), abs(B))), 1.0/(self.n**2))

	def test_rdirac(self):
		A = self.get_sample(common.ENSEMBLE_TYPES.RDIRAC)
		for row in A:
			self.assertEqual(np.linalg.norm( row), 1)

	def test_cdirac(self):
		A = self.get_sample(common.ENSEMBLE_TYPES.CDIRAC)
		for row in A:
			self.assertEqual(np.linalg.norm( row),1)
	
	def test_rgauss(self):
		A = self.get_sample(common.ENSEMBLE_TYPES.RGAUSS)
		B = A.flatten();
		self.assertAlmostEqual( np.mean( B), 0)
		self.assertAlmostEqual( np.mean( np.multiply( abs(B), abs(B))), 1.0/(self.n**2))
	
	def test_cgauss(self):
		A = self.get_sample(common.ENSEMBLE_TYPES.CGAUSS)
		B = A.flatten();
		self.assertAlmostEqual( np.mean( B), 0)
		self.assertAlmostEqual( np.mean( np.multiply( abs(B), abs(B))), 1.0/(self.n**2))

if __name__ == '__main__':
    unittest.main()
