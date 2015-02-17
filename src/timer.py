import time
class Timer:
	def tic(self):
	     self.t = time.time()
	def toc(self):
	     return time.time() - self.t
