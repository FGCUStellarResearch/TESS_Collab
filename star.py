import numpy as np

class Star(object):
	def __init__(self, time, flux):
		self.time = time
		self.flux = flux
		self.init_stats()

	def init_stats(self):
		# Mean flux
		self.fmean = np.mean(self.flux)
		
		# Standard deviation of twice-differenced (whitened) time series
		self.fstd = np.std(np.diff(np.diff(self.flux)))
		
		# Relative 5-95 percentile flux range
		frange = np.percentile(self.flux,95) - np.percentile(self.flux,5)
		frange = frange / np.mean(self.flux)
		self.frange = abs(frange)

		# Relative differenced (whitened) standard deviation
		srange = np.std(np.diff(self.flux)) / np.mean(self.flux)
		self.srange = abs(srange)
