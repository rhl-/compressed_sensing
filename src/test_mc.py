import unittest
import common
import mc_driver as driver
import numpy as np
import cmath

class target_test(unittest.TestCase):
	def setUp(self):
		self.digits=5

	def test_1(self):
		myDriver = driver.MCDriver()
		myDriver.run_trials()

if __name__ == '__main__':
	unittest.main()
