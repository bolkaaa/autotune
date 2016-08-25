#ifndef MAYER_H
#define MAYER_H

#define REAL float

extern "C" {
	void mayer_realfft(int n, REAL *real);
	void mayer_realifft(int n, REAL *real);
}

#endif
