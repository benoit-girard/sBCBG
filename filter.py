import numpy as np

def lowpass(series, fc, b):

	N = int(np.ceil((4 / b)))
	if not N % 2: N += 1  # Make sure that N is odd.
	n = np.arange(N)
	 
	# Compute sinc filter.
	h = np.sinc(2 * fc * (n - (N - 1) / 2.))
	 
	# Compute Blackman window.
	w = 0.42 - 0.5 * np.cos(2 * np.pi * n / (N - 1)) + 0.08 * np.cos(4 * np.pi * n / (N - 1))
	 
	# Multiply sinc filter with window.
	h = h * w
	 
	# Normalize to get unity gain.
	h = h / np.sum(h)

	return np.convolve(series, h)

def main():
	pass

if __name__ == '__main__':
  main()
