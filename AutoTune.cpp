/* AutoTune.cpp

   Ported to SuperCollider UGen by Bolka 2015
   
   A pitch-correcting LADSPA plugin.
   Free software by Thomas A. Baran.
   http://web.mit.edu/tbaran/www/autotalent.html
   VERSION 0.2
   March 20, 2010

   This program is free software; you can redistribute it and/or modify        
   it under the terms of the GNU General Public License as published by        
   the Free Software Foundation; either version 2 of the License, or           
   (at your option) any later version.                                         
	                                                                            
   This program is distributed in the hope that it will be useful,             
   but WITHOUT ANY WARRANTY; without even the implied warranty of              
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               
   GNU General Public License for more details.                                
	                                                                            
   You should have received a copy of the GNU General Public License           
   along with this program; if not, write to the Free Software                 
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  
   
*/
	/*****************************************************************************/

	#include <stdlib.h>
	#include <string.h>
	#include <math.h>
	#include <stdio.h>
	
	#include "SC_PlugIn.h"
	#include "mayer_fft.h"
	
	#define PI (float) 3.14159265358979323846
	#define L2SC (float) 3.32192809488736218171
	
	enum {
		AT_INPUT1,
		AT_TUNE,
		AT_FIXED,
		AT_PULL,
		AT_A,
		AT_Bb,
		AT_B,
		AT_C,
		AT_Db,
		AT_D,
		AT_Eb,
		AT_E,
		AT_F,
		AT_Gb,
		AT_G,
		AT_Ab,
		AT_AMOUNT,
		AT_SMOOTH,
		AT_SHIFT,
		AT_SCWARP,
		AT_LFOAMP,
		AT_LFORATE,
		AT_LFOSHAPE,
		AT_LFOSYMM,
		AT_LFOQUANT,
		AT_FCORR,
		AT_FWARP,
		AT_MIX
	};

	static InterfaceTable *ft;

	typedef struct {
	  int nfft;        // size of FFT
	  int numfreqs;    // number of frequencies represented (nfft/2 + 1)
	  float* fft_data; // array for writing/reading to/from FFT function
	} fft_vars;

	// oh god

	struct AutoTune : public Unit {
	  
	  // Params
	  float m_pfTune;
	  float m_pfFixed;
	  float m_pfPull;
	  int m_pfA;
	  int m_pfBb;
	  int m_pfB;
	  int m_pfC;
	  int m_pfDb;
	  int m_pfD;
	  int m_pfEb;
	  int m_pfE;
	  int m_pfF;
	  int m_pfGb;
	  int m_pfG;
	  int m_pfAb;
	  float m_pfAmount;
	  float m_pfSmooth;
	  float m_pfShift;
	  int m_pfScwarp;
	  float m_pfLfoamp;
	  float m_pfLforate;
	  float m_pfLfoshape;
	  float m_pfLfosymm;
	  int m_pfLfoquant;
	  int m_pfFcorr;
	  float m_pfFwarp;
	  float m_pfMix;
	  // End params
	  
	  fft_vars* fmembvars; // member variables for fft routine

	  unsigned long fs; // Sample rate

	  unsigned long cbsize; // size of circular buffer
	  unsigned long corrsize; // cbsize/2 + 1
	  unsigned long cbiwr;
	  unsigned long cbord;
	  float* cbi; // circular input buffer
	  float* cbf; // circular formant correction buffer
	  float* cbo; // circular output buffer

	  float* cbwindow; // hann of length N/2, zeros for the rest
	  float* acwinv; // inverse of autocorrelation of window
	  float* hannwindow; // length-N hann
	  int noverlap;

	  float* ffttime;
	  float* fftfreqre;
	  float* fftfreqim;

	  // VARIABLES FOR LOW-RATE SECTION
	  float aref; // A tuning reference (Hz)
	  float inpitch; // Input pitch (semitones)
	  float conf; // Confidence of pitch period estimate (between 0 and 1)
	  float outpitch; // Output pitch (semitones)
	  float vthresh; // Voiced speech threshold
	  
	  float pmax; // Maximum allowable pitch period (seconds)
	  float pmin; // Minimum allowable pitch period (seconds)
	  unsigned long nmax; // Maximum period index for pitch prd est
	  unsigned long nmin; // Minimum period index for pitch prd est

	  float lrshift; // Shift prescribed by low-rate section
	  int ptarget; // Pitch target, between 0 and 11
	  float sptarget; // Smoothed pitch target

	  float lfophase;

	  // VARIABLES FOR PITCH SHIFTER
	  float phprdd; // default (unvoiced) phase period
	  double inphinc; // input phase increment
	  double outphinc; // input phase increment
	  double phincfact; // factor determining output phase increment
	  double phasein;
	  double phaseout;
	  float* frag; // windowed fragment of speech
	  unsigned long fragsize; // size of fragment in samples

	  // VARIABLES FOR FORMANT CORRECTOR
	  int ford;
	  float falph;
	  float flamb;
	  float* fk;
	  float* fb;
	  float* fc;
	  float* frb;
	  float* frc;
	  float* fsig;
	  float* fsmooth;
	  float fhp;
	  float flp;
	  float flpa;
	  float** fbuff;
	  float* ftvec;
	  float fmute;
	  float fmutealph;
	};

	fft_vars* fft_con(AutoTune *unit, int nfft) 
	{
	  fft_vars* membvars = (fft_vars*) RTAlloc(unit->mWorld, sizeof(fft_vars));

	  membvars->nfft = nfft;
	  membvars->numfreqs = nfft / 2 + 1;

	  membvars->fft_data = (float*) RTAlloc(unit->mWorld, nfft * sizeof(float));
	  memset(membvars->fft_data, 0, nfft * sizeof(float));
	  
	  return membvars;
	}

	void fft_des(AutoTune *unit, fft_vars *membvars) 
	{
	  RTFree(unit->mWorld, membvars->fft_data);
	  RTFree(unit->mWorld, membvars);
	}

	void fft_forward(fft_vars *membvars, float *input, float *output_re, float *output_im)
	{
	  int ti;
	  int nfft;
	  int hnfft;
	  int numfreqs;

	  nfft = membvars->nfft;
	  hnfft = nfft/2;
	  numfreqs = membvars->numfreqs;

	  for (ti=0; ti<nfft; ti++) {
		membvars->fft_data[ti] = input[ti];
	  }

	  mayer_realfft(nfft, membvars->fft_data);

	  output_im[0] = 0;
	  for (ti=0; ti<hnfft; ti++) {
		output_re[ti] = membvars->fft_data[ti];
		output_im[ti+1] = membvars->fft_data[nfft-1-ti];
	  }
	  output_re[hnfft] = membvars->fft_data[hnfft];
	  output_im[hnfft] = 0;
	}

	void fft_inverse(fft_vars *membvars, float *input_re, float *input_im, float *output)
	{
	  int ti;
	  int nfft;
	  int hnfft;
	  int numfreqs;

	  nfft = membvars->nfft;
	  hnfft = nfft/2;
	  numfreqs = membvars->numfreqs;

	  for (ti = 0; ti<hnfft; ti++) {
		membvars->fft_data[ti] = input_re[ti];
		membvars->fft_data[nfft-1-ti] = input_im[ti+1];
	  }
	  
	  membvars->fft_data[hnfft] = input_re[hnfft];

	  mayer_realifft(nfft, membvars->fft_data);

	  for (ti = 0; ti<nfft; ti++) {
		output[ti] = membvars->fft_data[ti];
	  }
	}

	extern "C" {

		static void AutoTune_next(AutoTune *unit, int inNumSamples);
		static void AutoTune_Ctor(AutoTune *unit);
		static void AutoTune_Dtor(AutoTune *unit); 
	
	}
	
   /********************
    *  THE CONSTRUCTOR *
    ********************/

	void AutoTune_Ctor(AutoTune *unit) {
	
	  SETCALC(AutoTune_next);
	  unsigned long ti;
	  
	  unit->m_pfAmount = IN0(AT_AMOUNT);
	  unit->m_pfSmooth = IN0(AT_SMOOTH) * 0.8;
	  unit->m_pfTune = IN0(AT_TUNE);
	  unit->m_pfA = IN0(AT_A);
	  unit->m_pfBb = IN0(AT_Bb);
	  unit->m_pfB = IN0(AT_B);
	  unit->m_pfC = IN0(AT_C);
	  unit->m_pfDb = IN0(AT_Db);
	  unit->m_pfD = IN0(AT_D);
	  unit->m_pfEb = IN0(AT_Eb);
	  unit->m_pfE = IN0(AT_E);
	  unit->m_pfF = IN0(AT_F);
	  unit->m_pfGb = IN0(AT_Gb);
	  unit->m_pfG = IN0(AT_G);
	  unit->m_pfAb = IN0(AT_Ab);
	  unit->m_pfFixed = IN0(AT_FIXED);
	  unit->m_pfPull = IN0(AT_PULL);
	  unit->m_pfShift = IN0(AT_SHIFT);
	  unit->m_pfScwarp = IN0(AT_SCWARP);
	  unit->m_pfLfoamp = IN0(AT_LFOAMP);
	  unit->m_pfLforate = IN0(AT_LFORATE);
	  unit->m_pfLfoshape = IN0(AT_LFOSHAPE);
	  unit->m_pfLfosymm = IN0(AT_LFOSYMM);
	  unit->m_pfLfoquant = IN0(AT_LFOQUANT);
	  unit->m_pfFcorr = IN0(AT_FCORR);
	  unit->m_pfFwarp = IN0(AT_FWARP);
	  unit->m_pfMix = IN0(AT_MIX);

	  unit->aref = 440;
	  unit->fs = SAMPLERATE;

	  unit->inpitch = 0.0; // !!!

	  if (SAMPLERATE >= 88200) {
		unit->cbsize = 4096;
	  }
	  else {
		unit->cbsize = 2048;
	  }
	  
	  unit->corrsize = unit->cbsize / 2 + 1;

	  unit->pmax = 1/(float)70;  // max and min periods (ms)
	  unit->pmin = 1/(float)700; // eventually may want to bring these out as sliders

	  unit->nmax = (unsigned long)(SAMPLERATE * unit->pmax);
	  
	  if (unit->nmax > unit->corrsize) {
		unit->nmax = unit->corrsize;
	  }
	  
	  unit->nmin = (unsigned long)(SAMPLERATE * unit->pmin);

	  // allocate ...
	  
	  unit->cbi = (float *)RTAlloc(unit->mWorld, unit->cbsize * sizeof(float));
	  memset(unit->cbi, 0, unit->cbsize * sizeof(float));
	  
	  unit->cbf = (float *)RTAlloc(unit->mWorld, unit->cbsize * sizeof(float));
	  memset(unit->cbf, 0, unit->cbsize * sizeof(float));
	  
	  unit->cbo = (float *)RTAlloc(unit->mWorld, unit->cbsize * sizeof(float));
	  memset(unit->cbo, 0, unit->cbsize * sizeof(float));
	  
	  unit->cbiwr = 0;
	  unit->cbord = 0;
	  unit->lfophase = 0;

	  // Initialize formant corrector
	  unit->ford = 7; // should be sufficient to capture formants
	  unit->falph = pow(0.001, (float) 80 / (SAMPLERATE));
	  unit->flamb = -(0.8517*sqrt(atan(0.06583*SAMPLERATE))-0.1916); // or about -0.88 @ 44.1kHz
	  
	  // allocate ...
	  
	  unit->fk = (float *)RTAlloc(unit->mWorld, unit->ford *sizeof(float));
	  memset(unit->fk, 0, unit->ford *sizeof(float));
	  
	  unit->fb = (float *)RTAlloc(unit->mWorld, unit->ford * sizeof(float));
	  memset(unit->fb, 0, unit->ford * sizeof(float));
	  
	  unit->fc = (float *)RTAlloc(unit->mWorld, unit->ford * sizeof(float));
	  memset(unit->fc, 0, unit->ford * sizeof(float));
	  
	  unit->frb = (float *)RTAlloc(unit->mWorld, unit->ford *sizeof(float));
	  memset(unit->frb, 0, unit->ford *sizeof(float));
	  
	  unit->frc = (float *)RTAlloc(unit->mWorld, unit->ford * sizeof(float));
	  memset(unit->frc, 0, unit->ford * sizeof(float));
	  
	  unit->fsig = (float *)RTAlloc(unit->mWorld, unit->ford * sizeof(float));
	  memset(unit->fsig, 0, unit->ford * sizeof(float));
	  
	  unit->fsmooth = (float *)RTAlloc(unit->mWorld, unit->ford * sizeof(float));
	  memset(unit->fsmooth, 0, unit->ford * sizeof(float));
	  
	  // 
	  
	  unit->fhp = 0;
	  unit->flp = 0;
	  unit->flpa = pow(0.001, (float) 10 / (SAMPLERATE));
	  unit->fbuff = (float**) RTAlloc(unit->mWorld, (unit->ford) * sizeof(float*));
	  
	  for (ti=0; ti<unit->ford; ti++) {
		unit->fbuff[ti] = (float *)RTAlloc(unit->mWorld, unit->cbsize * sizeof(float));
		memset(unit->fbuff[ti], 0, unit->cbsize * sizeof(float));
	  }
	  
	  unit->ftvec = (float *)RTAlloc(unit->mWorld, unit->ford * sizeof(float));
	  memset(unit->ftvec, 0, unit->ford * sizeof(float));
	  
	  unit->fmute = 1;
	  unit->fmutealph = pow(0.001, (float)1 / (SAMPLERATE));

	  // Standard raised cosine window, max height at N/2
	  unit->hannwindow = (float *)RTAlloc(unit->mWorld, unit->cbsize * sizeof(float));
	  memset(unit->hannwindow, 0, unit->cbsize * sizeof(float));
	  
	  for (ti=0; ti<unit->cbsize; ti++) {
		unit->hannwindow[ti] = -0.5*cos(2*PI*ti/unit->cbsize) + 0.5;
	  }

	  // Generate a window with a single raised cosine from N/4 to 3N/4
	  unit->cbwindow = (float *)RTAlloc(unit->mWorld, unit->cbsize * sizeof(float));
	  memset(unit->cbwindow, 0, unit->cbsize * sizeof(float));
	  
	  for (ti=0; ti<(unit->cbsize / 2); ti++) {
		unit->cbwindow[ti+unit->cbsize/4] = -0.5*cos(4*PI*ti/(unit->cbsize - 1)) + 0.5;
	  }

	  unit->noverlap = 4;

	  unit->fmembvars = fft_con(unit, unit->cbsize);

	  unit->ffttime = (float *)RTAlloc(unit->mWorld, unit->cbsize * sizeof(float));
	  memset(unit->ffttime, 0, unit->cbsize * sizeof(float));
	  
	  unit->fftfreqre = (float *)RTAlloc(unit->mWorld, unit->corrsize * sizeof(float));
	  memset(unit->fftfreqre, 0, unit->corrsize * sizeof(float));
	  
	  unit->fftfreqim = (float *)RTAlloc(unit->mWorld, unit->corrsize * sizeof(float));
 	  memset(unit->fftfreqim, 0, unit->corrsize * sizeof(float));
	  
	  // ---- Calculate autocorrelation of window ----
	  
	  unit->acwinv = (float *)RTAlloc(unit->mWorld, unit->cbsize * sizeof(float));
	  memset(unit->acwinv, 0, unit->cbsize * sizeof(float));
	  
	  for (ti=0; ti<unit->cbsize; ti++) {
		unit->ffttime[ti] = unit->cbwindow[ti];
	  }
	  
	  fft_forward(unit->fmembvars, unit->cbwindow, unit->fftfreqre, unit->fftfreqim);
	  
	  for (ti=0; ti<unit->corrsize; ti++) {
		unit->fftfreqre[ti] = (unit->fftfreqre[ti])*(unit->fftfreqre[ti]) + (unit->fftfreqim[ti])*(unit->fftfreqim[ti]);
		unit->fftfreqim[ti] = 0;
	  }
	  
	  fft_inverse(unit->fmembvars, unit->fftfreqre, unit->fftfreqim, unit->ffttime);
	  
	  for (ti=1; ti<unit->cbsize; ti++) {
		unit->acwinv[ti] = unit->ffttime[ti]/unit->ffttime[0];
		if (unit->acwinv[ti] > 0.000001) {
		  unit->acwinv[ti] = (float)1/unit->acwinv[ti];
		}
		else {
		  unit->acwinv[ti] = 0;
		}
	  }
	  
	  unit->acwinv[0] = 1;
	  
	  // ---- END Calculate autocorrelation of window ----
	  
	  unit->lrshift = 0;
	  unit->ptarget = 0;
	  unit->sptarget = 0;

	  unit->vthresh = 0.7;  //  The voiced confidence (unbiased peak) threshold level

	  // Pitch shifter initialization
	  
	  unit->phprdd = 0.01; // Default period
	  unit->inphinc = (float)1/(unit->phprdd * SAMPLERATE);
	  unit->phincfact = 1;
	  unit->phasein = 0;
	  unit->phaseout = 0;
	
	  unit->frag = (float *)RTAlloc(unit->mWorld, unit->cbsize * sizeof(float));
	  memset(unit->frag, 0, unit->cbsize * sizeof(float));
	
	  unit->fragsize = 0;
	  	  
	  AutoTune_next(unit, 1);
	  
	  //printf("AutoTune memory allocation finished\n");
	
	}

	// Called every time we get a new chunk of audio
	void AutoTune_next(AutoTune *unit, int inNumSamples) 
	{
	  float *pfInput = IN(0);
	  
	  float *pfOutput = OUT(0);
	  float *pfOutpitch = OUT(1);
	  float *pfOutconf = OUT(2);
	  
	  float fAmount;
	  float fSmooth;
	  int iNotes[12];
	  int iPitch2Note[12];
	  int iNote2Pitch[12];
	  int numNotes;
	  float fTune;
	  float fFixed;
	  float fPull;
	  float fShift;
	  int iScwarp;
	  float fLfoamp;
	  float fLforate;
	  float fLfoshape;
	  float fLfosymm;
	  int iLfoquant;
	  int iFcorr;
	  float fFwarp;
	  float fMix;
	  
	  unsigned long lSampleIndex;

	  long int N;
	  long int Nf;
	  long int fs;
	  float pmin;
	  float pmax;
	  unsigned long nmin;
	  unsigned long nmax;

	  long int ti;
	  long int ti2;
	  long int ti3;
	  long int ti4;
	  float tf;
	  float tf2;

	  // Variables for cubic spline interpolator
	  float indd;
	  int ind0;
	  int ind1;
	  int ind2;
	  int ind3;
	  float vald;
	  float val0;
	  float val1;
	  float val2;
	  float val3;

	  int lowersnap;
	  int uppersnap;
	  float lfoval;

	  float pperiod;
	  float inpitch;
	  float conf;
	  float outpitch;
	  float aref;
	  float fa;
	  float fb;
	  float fc;
	  float fk;
	  float flamb;
	  float frlamb;
	  float falph;
	  float foma;
	  float f1resp;
	  float f0resp;
	  float flpa;
	  int ford;
	  
	  fAmount = IN0(AT_AMOUNT);
	  fSmooth = IN0(AT_SMOOTH) * 0.8;
	  fTune = IN0(AT_TUNE);
	  iNotes[0] = IN0(AT_A);
	  iNotes[1] = IN0(AT_Bb);
	  iNotes[2] = IN0(AT_B);
	  iNotes[3] = IN0(AT_C);
	  iNotes[4] = IN0(AT_Db);
	  iNotes[5] = IN0(AT_D);
	  iNotes[6] = IN0(AT_Eb);
	  iNotes[7] = IN0(AT_E);
	  iNotes[8] = IN0(AT_F);
	  iNotes[9] = IN0(AT_Gb);
	  iNotes[10] = IN0(AT_G);
	  iNotes[11] = IN0(AT_Ab);
	  fFixed = IN0(AT_FIXED);
	  fPull = IN0(AT_PULL);
	  fShift = IN0(AT_SHIFT);
	  iScwarp = IN0(AT_SCWARP);
	  fLfoamp = IN0(AT_LFOAMP);
	  fLforate = IN0(AT_LFORATE);
	  fLfoshape = IN0(AT_LFOSHAPE);
	  fLfosymm = IN0(AT_LFOSYMM);
	  iLfoquant = IN0(AT_LFOQUANT);
	  iFcorr = IN0(AT_FCORR);
	  fFwarp = IN0(AT_FWARP);
	  fMix = IN0(AT_MIX);

	  // Some logic for the semitone->scale and scale->semitone conversion
	  // If no notes are selected as being in the scale, instead snap to all notes
	  
	  ti2 = 0;
	  for (ti=0; ti<12; ti++) {
		if (iNotes[ti]>=0) {
		  iPitch2Note[ti] = ti2;
		  iNote2Pitch[ti2] = ti;
		  ti2 = ti2 + 1;
		}
		else {
		  iPitch2Note[ti] = -1;
		}
	  }
	  numNotes = ti2;
	  while (ti2<12) {
		iNote2Pitch[ti2] = -1;
		ti2 = ti2 + 1;
	  }
	  if (numNotes==0) {
		for (ti=0; ti<12; ti++) {
		  iNotes[ti] = 1;
		  iPitch2Note[ti] = ti;
		  iNote2Pitch[ti] = ti;
		}
		numNotes = 12;
	  }
	  
	  iScwarp = (iScwarp + numNotes*5)%numNotes;

	  ford = unit->ford;
	  falph = unit->falph;
	  foma = (float)1 - falph;
	  flpa = unit->flpa;
	  flamb = unit->flamb;
	  tf = pow((float)2,fFwarp/2)*(1+flamb)/(1-flamb);
	  frlamb = (tf - 1)/(tf + 1);

	  unit->aref = (float)fTune;

	  N = unit->cbsize;
	  Nf = unit->corrsize;
	  fs = unit->fs;

	  pmax = unit->pmax;
	  pmin = unit->pmin;
	  nmax = unit->nmax;
	  nmin = unit->nmin;

	  aref = unit->aref;
	  pperiod = unit->pmax;
	  inpitch = unit->inpitch;
	  conf = unit->conf;
	  outpitch = unit->outpitch;

	  /*******************
	   *  MAIN DSP LOOP  *
	   *******************/
	   
	   for (lSampleIndex = 0; lSampleIndex < inNumSamples; lSampleIndex++)  {
    
    // load data into circular buffer
    tf = (float) *(pfInput++);
    ti4 = unit->cbiwr;
    unit->cbi[ti4] = tf;

    if (iFcorr>=1) {
      // Somewhat experimental formant corrector
      //  formants are removed using an adaptive pre-filter and
      //  re-introduced after pitch manipulation using post-filter
      // tf is signal input
      fa = tf - unit->fhp; // highpass pre-emphasis filter
      unit->fhp = tf;
      fb = fa;
      for (ti=0; ti<ford; ti++) {
	unit->fsig[ti] = fa*fa*foma + unit->fsig[ti]*falph;
	fc = (fb-unit->fc[ti])*flamb + unit->fb[ti];
	unit->fc[ti] = fc;
	unit->fb[ti] = fb;
	fk = fa*fc*foma + unit->fk[ti]*falph;
	unit->fk[ti] = fk;
	tf = fk/(unit->fsig[ti] + 0.000001);
	tf = tf*foma + unit->fsmooth[ti]*falph;
	unit->fsmooth[ti] = tf;
	unit->fbuff[ti][ti4] = tf;
	fb = fc - tf*fa;
	fa = fa - tf*fc;
      }
      unit->cbf[ti4] = fa;
      // Now hopefully the formants are reduced
      // More formant correction code at the end of the DSP loop
    }
    else {
      unit->cbf[ti4] = tf;
    }


    // Input write pointer logic
    unit->cbiwr++;
    if (unit->cbiwr >= N) {
      unit->cbiwr = 0;
    }


    // ********************
    // * Low-rate section *
    // ********************

    // Every N/noverlap samples, run pitch estimation / manipulation code
    if ((unit->cbiwr)%(N/unit->noverlap) == 0) {

      // ---- Obtain autocovariance ----

      // Window and fill FFT buffer
      ti2 = unit->cbiwr;
      for (ti=0; ti<N; ti++) {
	unit->ffttime[ti] = (float)(unit->cbi[(ti2-ti+N)%N]*unit->cbwindow[ti]);
      }

      // Calculate FFT
      fft_forward(unit->fmembvars, unit->ffttime, unit->fftfreqre, unit->fftfreqim);

      // Remove DC
      unit->fftfreqre[0] = 0;
      unit->fftfreqim[0] = 0;

      // Take magnitude squared
      for (ti=1; ti<Nf; ti++) {
	unit->fftfreqre[ti] = (unit->fftfreqre[ti])*(unit->fftfreqre[ti]) + (unit->fftfreqim[ti])*(unit->fftfreqim[ti]);
	unit->fftfreqim[ti] = 0;
      }

      // Calculate IFFT
      fft_inverse(unit->fmembvars, unit->fftfreqre, unit->fftfreqim, unit->ffttime);

      // Normalize
      tf = (float)1/unit->ffttime[0];
      for (ti=1; ti<N; ti++) {
	unit->ffttime[ti] = unit->ffttime[ti] * tf;
      }
      unit->ffttime[0] = 1;

      //  ---- END Obtain autocovariance ----


      //  ---- Calculate pitch and confidence ----

      // Calculate pitch period
      //   Pitch period is determined by the location of the max (biased)
      //     peak within a given range
      //   Confidence is determined by the corresponding unbiased height
      tf2 = 0;
      pperiod = pmin;
      for (ti=nmin; ti<nmax; ti++) {
	ti2 = ti-1;
	ti3 = ti+1;
	if (ti2<0) {
	  ti2 = 0;
	}
	if (ti3>Nf) {
	  ti3 = Nf;
	}
	tf = unit->ffttime[ti];

	if (tf>unit->ffttime[ti2] && tf>=unit->ffttime[ti3] && tf>tf2) {
	  tf2 = tf;
	  ti4 = ti;
	}
      }
      if (tf2>0) {
	conf = tf2*unit->acwinv[ti4];
	if (ti4>0 && ti4<Nf) {
	  // Find the center of mass in the vicinity of the detected peak
	  tf = unit->ffttime[ti4-1]*(ti4-1);
	  tf = tf + unit->ffttime[ti4]*(ti4);
	  tf = tf + unit->ffttime[ti4+1]*(ti4+1);
	  tf = tf/(unit->ffttime[ti4-1] + unit->ffttime[ti4] + unit->ffttime[ti4+1]);
	  pperiod = tf/fs;
	}
	else {
	  pperiod = (float)ti4/fs;
	}
      }

      // Convert to semitones
      tf = (float) -12*log10((float)aref*pperiod)*L2SC;
      if (conf>=unit->vthresh) {
	inpitch = tf;
	unit->inpitch = tf; // update pitch only if voiced
      }
      unit->conf = conf;

      //  ---- END Calculate pitch and confidence ----


      //  ---- Modify pitch in all kinds of ways! ----

      outpitch = inpitch;

      // Pull to fixed pitch
      outpitch = (1-fPull)*outpitch + fPull*fFixed;

      // -- Convert from semitones to scale notes --
      ti = (int)(outpitch/12 + 32) - 32; // octave
      tf = outpitch - ti*12; // semitone in octave
      ti2 = (int)tf;
      ti3 = ti2 + 1;
      // a little bit of pitch correction logic, since it's a convenient place for it
      if (iNotes[ti2%12]<0 || iNotes[ti3%12]<0) { // if between 2 notes that are more than a semitone apart
	lowersnap = 1;
	uppersnap = 1;
      }
      else {
	lowersnap = 0;
	uppersnap = 0;
	if (iNotes[ti2%12]==1) { // if specified by user
	  lowersnap = 1;
	}
	if (iNotes[ti3%12]==1) { // if specified by user
	  uppersnap = 1;
	}
      }
      // (back to the semitone->scale conversion)
      // finding next lower pitch in scale
      while (iNotes[(ti2+12)%12]<0) {
      	ti2 = ti2 - 1;
      }
      // finding next higher pitch in scale
      while (iNotes[ti3%12]<0) {
      	ti3 = ti3 + 1;
      }
      tf = (tf-ti2)/(ti3-ti2) + iPitch2Note[(ti2+12)%12];
      if (ti2<0) {
      	tf = tf - numNotes;
      }
      outpitch = tf + numNotes*ti;
      // -- Done converting to scale notes --

      // The actual pitch correction
      ti = (int)(outpitch+128) - 128;
      tf = outpitch - ti - 0.5;
      ti2 = ti3-ti2;
      if (ti2>2) { // if more than 2 semitones apart, put a 2-semitone-like transition halfway between
	tf2 = (float)ti2/2;
      }
      else {
	tf2 = (float)1;
      }
      if (fSmooth<0.001) {
	tf2 = tf*tf2/0.001;
      }
      else {
	tf2 = tf*tf2/fSmooth;
      }
      if (tf2<-0.5) tf2 = -0.5;
      if (tf2>0.5) tf2 = 0.5;
      tf2 = 0.5*sin(PI*tf2) + 0.5; // jumping between notes using horizontally-scaled sine segment
      tf2 = tf2 + ti;
      if ( (tf<0.5 && lowersnap) || (tf>=0.5 && uppersnap) ) {
	outpitch = fAmount*tf2 + ((float)1-fAmount)*outpitch;
      }

      // Add in pitch shift
      outpitch = outpitch + fShift;

      // LFO logic
      tf = fLforate*N/(unit->noverlap*fs);
      if (tf>1) tf=1;
      unit->lfophase = unit->lfophase + tf;
      if (unit->lfophase>1) unit->lfophase = unit->lfophase-1;
      lfoval = unit->lfophase;
      tf = (fLfosymm + 1)/2;
      if (tf<=0 || tf>=1) {
	if (tf<=0) lfoval = 1-lfoval;
      }
      else {
	if (lfoval<=tf) {
	  lfoval = lfoval/tf;
	}
	else {
	  lfoval = 1 - (lfoval-tf)/(1-tf);
	}
      }
      if (fLfoshape>=0) {
	// linear combination of cos and line
	lfoval = (0.5 - 0.5*cos(lfoval*PI))*fLfoshape + lfoval*(1-fLfoshape);
	lfoval = fLfoamp*(lfoval*2 - 1);
      }
      else {
	// smoosh the sine horizontally until it's squarish
	tf = 1 + fLfoshape;
	if (tf<0.001) {
	  lfoval = (lfoval - 0.5)*2/0.001;
	}
	else {
	  lfoval = (lfoval - 0.5)*2/tf;
	}
	if (lfoval>1) lfoval = 1;
	if (lfoval<-1) lfoval = -1;
	lfoval = fLfoamp*sin(lfoval*PI*0.5);
      }
      // add in quantized LFO
      if (iLfoquant>=1) {
	outpitch = outpitch + (int)(numNotes*lfoval + numNotes + 0.5) - numNotes;
      }


      // Convert back from scale notes to semitones
      outpitch = outpitch + iScwarp; // output scale rotate implemented here
      ti = (int)(outpitch/numNotes + 32) - 32;
      tf = outpitch - ti*numNotes;
      ti2 = (int)tf;
      ti3 = ti2 + 1;
      outpitch = iNote2Pitch[ti3%numNotes] - iNote2Pitch[ti2];
      if (ti3>=numNotes) {
	outpitch = outpitch + 12;
      }
      outpitch = outpitch*(tf - ti2) + iNote2Pitch[ti2];
      outpitch = outpitch + 12*ti;
      outpitch = outpitch - (iNote2Pitch[iScwarp] - iNote2Pitch[0]); //more scale rotation here

      // add in unquantized LFO
      if (iLfoquant<=0) {
	outpitch = outpitch + lfoval*2;
      }


      if (outpitch<-36) outpitch = -48;
      if (outpitch>24) outpitch = 24;

      unit->outpitch = outpitch;

      //  ---- END Modify pitch in all kinds of ways! ----

      // Compute variables for pitch shifter that depend on pitch
      unit->inphinc = aref*pow(2,inpitch/12)/fs;
      unit->outphinc = aref*pow(2,outpitch/12)/fs;
      unit->phincfact = unit->outphinc/unit->inphinc;
    }
    // ************************
    // * END Low-Rate Section *
    // ************************



    // *****************
    // * Pitch Shifter *
    // *****************

    // Pitch shifter (kind of like a pitch-synchronous version of Fairbanks' technique)
    //   Note: pitch estimate is naturally N/2 samples old
    unit->phasein = unit->phasein + unit->inphinc;
    unit->phaseout = unit->phaseout + unit->outphinc;

    //   When input phase resets, take a snippet from N/2 samples in the past
    if (unit->phasein >= 1) {
      unit->phasein = unit->phasein - 1;
      ti2 = unit->cbiwr - N/2;
      for (ti=-N/2; ti<N/2; ti++) {
	unit->frag[(ti+N)%N] = unit->cbf[(ti + ti2 + N)%N];
      }
    }

    //   When output phase resets, put a snippet N/2 samples in the future
    if (unit->phaseout >= 1) {
      unit->fragsize = unit->fragsize*2;
      if (unit->fragsize > N) {
	unit->fragsize = N;
      }
      unit->phaseout = unit->phaseout - 1;
      ti2 = unit->cbord + N/2;
      ti3 = (long int)(((float)unit->fragsize) / unit->phincfact);
      if (ti3>=N/2) {
	ti3 = N/2 - 1;
      }
      for (ti=-ti3/2; ti<(ti3/2); ti++) {
	tf = unit->hannwindow[(long int)N/2 + ti*(long int)N/ti3];
	// 3rd degree polynomial interpolator - based on eqns from Hal Chamberlin's book
	indd = unit->phincfact*ti;
	ind1 = (int)indd;
	ind2 = ind1+1;
	ind3 = ind1+2;
	ind0 = ind1-1;
	val0 = unit->frag[(ind0+N)%N];
	val1 = unit->frag[(ind1+N)%N];
	val2 = unit->frag[(ind2+N)%N];
	val3 = unit->frag[(ind3+N)%N];
	vald = 0;
	vald = vald - (float)0.166666666667 * val0 * (indd - ind1) * (indd - ind2) * (indd - ind3);
	vald = vald + (float)0.5 * val1 * (indd - ind0) * (indd - ind2) * (indd - ind3);
	vald = vald - (float)0.5 * val2 * (indd - ind0) * (indd - ind1) * (indd - ind3);
	vald = vald + (float)0.166666666667 * val3 * (indd - ind0) * (indd - ind1) * (indd - ind2);
	unit->cbo[(ti + ti2 + N)%N] = unit->cbo[(ti + ti2 + N)%N] + vald*tf;
      }
      unit->fragsize = 0;
    }
    unit->fragsize++;

    //   Get output signal from buffer
    tf = unit->cbo[unit->cbord]; // read buffer

    unit->cbo[unit->cbord] = 0; // erase for next cycle
    unit->cbord++; // increment read pointer
    if (unit->cbord >= N) {
      unit->cbord = 0;
    }

    // *********************
    // * END Pitch Shifter *
    // *********************

    ti4 = (unit->cbiwr + 2)%N;
    if (iFcorr>=1) {
      // The second part of the formant corrector
      // This is a post-filter that re-applies the formants, designed
      //   to result in the exact original signal when no pitch
      //   manipulation is performed.
      // tf is signal input
      // gotta run it 3 times because of a pesky delay free loop
      //  first time: compute 0-response
      tf2 = tf;
      fa = 0;
      fb = fa;
      for (ti=0; ti<ford; ti++) {
	fc = (fb-unit->frc[ti])*frlamb + unit->frb[ti];
	tf = unit->fbuff[ti][ti4];
	fb = fc - tf*fa;
	unit->ftvec[ti] = tf*fc;
	fa = fa - unit->ftvec[ti];
      }
      tf = -fa;
      for (ti=ford-1; ti>=0; ti--) {
	tf = tf + unit->ftvec[ti];
      }
      f0resp = tf;
      //  second time: compute 1-response
      fa = 1;
      fb = fa;
      for (ti=0; ti<ford; ti++) {
	fc = (fb-unit->frc[ti])*frlamb + unit->frb[ti];
	tf = unit->fbuff[ti][ti4];
	fb = fc - tf*fa;
	unit->ftvec[ti] = tf*fc;
	fa = fa - unit->ftvec[ti];
      }
      tf = -fa;
      for (ti=ford-1; ti>=0; ti--) {
	tf = tf + unit->ftvec[ti];
      }
      f1resp = tf;
      //  now solve equations for output, based on 0-response and 1-response
      tf = (float)2*tf2;
      tf2 = tf;
      tf = ((float)1 - f1resp + f0resp);
      if (tf!=0) {
	tf2 = (tf2 + f0resp) / tf;
      }
      else {
	tf2 = 0;
      }
      //  third time: update delay registers
      fa = tf2;
      fb = fa;
      for (ti=0; ti<ford; ti++) {
	fc = (fb-unit->frc[ti])*frlamb + unit->frb[ti];
	unit->frc[ti] = fc;
	unit->frb[ti] = fb;
	tf = unit->fbuff[ti][ti4];
	fb = fc - tf*fa;
	fa = fa - tf*fc;
      }
      tf = tf2;
      tf = tf + flpa*unit->flp;  // lowpass post-emphasis filter
      unit->flp = tf;
      // Bring up the gain slowly when formant correction goes from disabled
      // to enabled, while things stabilize.
      if (unit->fmute>0.5) {
	tf = tf*(unit->fmute - 0.5)*2;
      }
      else {
	tf = 0;
      }
      tf2 = unit->fmutealph;
      unit->fmute = (1-tf2) + tf2*unit->fmute;
      // now tf is signal output
      // ...and we're done messing with formants
    }
    else {
      unit->fmute = 0;
    }

    // Write audio to output of plugin
    // Mix (blend between original (delayed) =0 and processed =1)
    *(pfOutput++) = (float) fMix*tf + (1.0-fMix)*unit->cbi[ti4];
    *(pfOutpitch++) = (float) outpitch; 
    *(pfOutconf++) = (float) conf;

  }

	}

	/********************
	 *  THE DESTRUCTOR! *
	 ********************/
	void AutoTune_Dtor(AutoTune *unit) 
	{
	  int ti;
	  fft_des(unit, unit->fmembvars);
	  
	  RTFree(unit->mWorld, unit->cbi);
	  RTFree(unit->mWorld, unit->cbf);
	  RTFree(unit->mWorld, unit->cbo);
	  RTFree(unit->mWorld, unit->cbwindow);
	  RTFree(unit->mWorld, unit->hannwindow);
	  RTFree(unit->mWorld, unit->acwinv);
	  RTFree(unit->mWorld, unit->frag);
	  RTFree(unit->mWorld, unit->ffttime);
	  RTFree(unit->mWorld, unit->fftfreqre);
	  RTFree(unit->mWorld, unit->fftfreqim);
	  RTFree(unit->mWorld, unit->fk);
	  RTFree(unit->mWorld, unit->fb);
	  RTFree(unit->mWorld, unit->fc);
	  RTFree(unit->mWorld, unit->frb);
	  RTFree(unit->mWorld, unit->frc);
	  RTFree(unit->mWorld, unit->fsmooth);
	  RTFree(unit->mWorld, unit->fsig);
	  
	  for (ti=0; ti<unit->ford; ti++) {
		RTFree(unit->mWorld, unit->fbuff[ti]);
	  }
	  
	  RTFree(unit->mWorld, unit->fbuff);
	  RTFree(unit->mWorld, unit->ftvec);
	}

	PluginLoad(AutoTune)
	{ 
		ft = inTable; 
    	DefineDtorUnit(AutoTune);
	};

	// All done
