
	AutoTune : MultiOutUGen {

		*ar { arg in, tune = 440.0, fixed = 0, pull = 0,
				  scale = #[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
				  amount = 0, smooth = 0, shift = 0, scwarp = 0,
				  lfoamp = 0, lforate = 0.5, lfosymm = 0, lfoshape = 0,
				  lfoquant = 0, fcorr = 0, fwarp = 0, fmix = 1,
				  mul = 1.0, add = 0.0;

			^this.multiNew('audio', in, tune, fixed, pull,
				  scale[0], scale[1], scale[2], scale[3],
				  scale[4], scale[5], scale[6], scale[7],
				  scale[8], scale[9], scale[10], scale[11],
				  amount, smooth, shift, scwarp,
				  lfoamp, lforate, lfoshape, lfosymm, lfoquant,
				  fcorr, fwarp, fmix).madd(mul, add)
		}

		init {arg ... theInputs;
			inputs = theInputs;
			channels = [
				OutputProxy(rate, this, 0),
				OutputProxy(rate, this, 1),
				OutputProxy(rate, this, 2)
			];
			^channels
		}

	}

    // USAGE
    
    /*
        ~tuneback = { arg tune1 = 0, tune2 = 0;
            TuneBack.ar(In.ar(0, 1), In.ar(1, 1), tune1, tune2);
        }.play;
        
        try different settings:
        
        ~tuneback.set(\tune1, 1);
        ~tuneback.set(\tune1, 0.5);
        ~tuneback.set(\tune2, 0.75);
        ~tuneback.set(\tune2, 1);
        ~tuneback.set(\tune2, 0);
        
    */

	TuneBack : UGen {

		*ar { arg in1, in2, tune1 = 0.0, tune2 = 0.0;

			var sig1, sig2, pitch1 = 0, pitch2 = 0, conf1 = 0, conf2 = 0;
			var fb;

			fb = LocalIn.ar(4, 0);

			#sig1, pitch1, conf1 = AutoTune.ar(in1, fixed: fb[2], pull: tune1, fcorr:1);
			#sig2, pitch2, conf2 = AutoTune.ar(in2, fixed: fb[0], pull: tune2, fcorr:1);

			LocalOut.ar([pitch1, conf1, pitch2, conf2]);

			^[sig1, sig2];
		}

	}
