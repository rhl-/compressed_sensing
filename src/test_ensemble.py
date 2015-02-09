import unittest
import common
import numpy as np
import cmath
import gen_ensemble as ensemble


class ensemble_test(unittest.TestCase):

	def get_sample(self,t1, t2=common.TARGET_TYPES.HPSD):
		return ensemble.getEnsembleSample(self.m, self.n, t1, t2)
	def setUp(self):
		self.m = 128
		self.n = 128

	def test_entry(self):
		A = self.get_sample(common.ENSEMBLE_TYPES.ENTRY)
		nnz = np.count_nonzero(A)
 		self.assertTrue(nnz == self.m)
		for matrix in A:
			self.assertTrue(np.count_nonzero(matrix == 1))
		
	def test_perm(self):
		A = self.get_sample(common.ENSEMBLE_TYPES.PERM)
		nnz = np.count_nonzero( A)
 		self.assertTrue( nnz == self.m*self.n*self.n)
		for matrix in A:
			self.assertTrue( np.count_nonzero( matrix) == self.n**2)
			for row in np.reshape(matrix, (n,n)):
				self.assertTrue( np.count_nonzero( row) == 1)
			for column in np.reshape(matrix, (n,n)).T:
				self.assertTrue( np.count_nonzero( column) == 1)

	def test_rsperm(self):
		A = self.get_sample(common.ENSEMBLE_TYPES.RSPERM)
		nnz = np.count_nonzero( A)
 		self.assertTrue( nnz == self.m*self.n*self.n)
		for matrix in A:
			self.assertTrue( np.count_nonzero( matrix) == self.n**2)
			for row in np.reshape(matrix, (n,n)):
				self.assertTrue( np.count_nonzero( row) == 1)
			for column in np.reshape(matrix, (n,n)).T:
				self.assertTrue( np.count_nonzero( column) == 1)
		B = n*np.multiply(A,A);	
		self.assertTrue( abs(sum(sum(B)) - np.count_nonzero( B)) < 1e-4)

	def test_csperm(self):
		A = self.get_sample(common.ENSEMBLE_TYPES.CSPERM)
		nnz = np.count_nonzero( A)
 		self.assertTrue( nnz == self.m*self.n*self.n)
		for matrix in A:
			self.assertTrue( np.count_nonzero( matrix) == self.n**2)
			for row in np.reshape(matrix, (n,n)):
				self.assertTrue( np.count_nonzero( row) == 1)
			for column in np.reshape(matrix, (n,n)).T:
				self.assertTrue( np.count_nonzero( column) == 1)
		B = np.multiply(A,A)
		B = np.multiply(B,B)
		B = (n**2)*B
		self.assertTrue( abs(sum(sum(B)) - np.count_nonzero( B)) < 1e-4)

	def test_rgperm(self):
		A = self.get_sample(common.ENSEMBLE_TYPES.RGPERM)
		nnz = np.count_nonzero( A)
 		self.assertTrue( nnz == self.m*self.n*self.n)
		for matrix in A:
			self.assertTrue( np.count_nonzero( matrix) == self.n**2)
			for row in np.reshape(matrix, (n,n)):
				self.assertTrue( np.count_nonzero( row) == 1)
			for column in np.reshape(matrix, (n,n)).T:
				self.assertTrue( np.count_nonzero( column) == 1)
		B = n*np.multiply(A,A)
		self.assertTrue( abs(sum(sum(B)) - np.count_nonzero( B)) < 1e-4)

	def test_cgperm(self):
		A = self.get_sample(common.ENSEMBLE_TYPES.CGPERM)
		nnz = np.count_nonzero( A)
 		self.assertTrue( nnz == self.m*self.n*self.n)
		for matrix in A:
			self.assertTrue( np.count_nonzero( matrix) == self.n**2)
			for row in np.reshape(matrix, (self.n,self.n)):
				self.assertTrue( np.count_nonzero( row) == 1)
			for column in np.reshape(matrix, (self.n,self.n)).T:
				self.assertTrue( np.count_nonzero( column) == 1)
		B = np.multiply(A,A)
		B = np.multiply(B,B)
		B = (n**2)*B
		self.assertTrue( abs(sum(sum(B)) - np.count_nonzero( B)) < 1e-4)
		B = np.nonzeroflatten( A);
		self.assertTrue( np.mean( B) < 1e-3)
		self.assertTrue( abs(np.mean( numpy.multiplies(B,B)) - 1.0/n) < 1e-3)

	def test_rdirac(self):
		A = self.get_sample(common.ENSEMBLE_TYPES.RDIRAC)
		for row in A:
			self.assertTrue(np.linalg.norm( row) == 1)

	def test_cdirac(self):
		A = self.get_sample(common.ENSEMBLE_TYPES.CDIRAC)
		for row in A:
			self.assertTrue(np.linalg.norm( row) == 1)
	
	def test_rgauss(self):
		A = self.get_sample(common.ENSEMBLE_TYPES.RGAUSS)
		B = A.flatten();
		self.assertTrue( abs(np.mean( B)) < 1e-3)
		self.assertTrue( abs(np.mean( np.multiply( abs(B), abs(B))) - 1.0/(self.n**2)) < 1e-3)
	
	def test_cgauss(self):
		A = self.get_sample(common.ENSEMBLE_TYPES.CGAUSS)
		B = A.flatten();
		self.assertTrue( abs(np.mean( B)) < 1e-3)
		self.assertTrue( abs(np.mean( np.multiply( abs(B), abs(B))) - 1.0/(self.n**2)) < 1e-3)

if __name__ == '__main__':
    unittest.main()
