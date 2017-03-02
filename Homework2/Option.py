import numpy as np

class Option(object):

	def __init__(self, S0, K, sigma, T, r):
		self._S0 = S0
		self._K = K
		self._sigma = sigma
		self._T = T
		self._r = r
		self._N = 50
		self._dt = self._T / self._N
		self._sqrt_dt = np.sqrt(self._dt)
		self._u = np.exp(self._sigma * self._sqrt_dt)
		self._d = 1 / self._u
		self._p = (np.exp(self._r * self._dt) - self._d) / (self._u - self._d)
		self._discount_factor = np.exp(-self._r * self._T)

	def getU(self):
		return self._u

	def getD(self):
		return self._d

	def getS0(self):
		return self._S0

	def setS0(self, S0):
		self._S0 = S0

	def getK(self):
		return self._K

	def setK(self, K):
		self._K = K

	def getVol(self):
		return self._sigma

	def setVol(self, sigma):
		self._sigma = sigma

	def getMaturity(self):
		return self._T

	def setMaturity(self, T):
		self._T = T

	def price(self, t, mu, W):
		return self._S0 * np.exp(mu * t + 
			self._sigma * self._sqrt_dt * W)

	def payoff(self, S):
		return max(S - self._K, 0)



