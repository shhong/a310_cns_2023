TITLE Voltage-dependent spike time-dependent plasticity
COMMENT

Written by Sungo Hong, CNS Unit, OIST, Japan
March 2017

ENDCOMMENT

NEURON {
	POINT_PROCESS CCSTDPSyn
	NONSPECIFIC_CURRENT i
	RANGE u1, u2, u3, vca
	RANGE taum, taup, taux, tauca
	RANGE Am, Ap, thetam, thetap, wmin, wmax, tau1, tau2
}

UNITS {
      (nA) = (nanoamp)
      (mV) = (millivolt)
      (uS) = (microsiemens)
}

PARAMETER {
	e = 0	(mV)
	tau1 = 0.1 (ms) <1e-9,1e9>
	tau2 = 1.5 (ms) <1e-9,1e9>

	taum = 6 (ms)      : time constant for filtering membrane potential
	taup = 2 (ms)       : time constant for post filter for potentiation
	taux = 28 (ms)      : time constant for low-pass r
	tauca = 5 (ms)
	Am = 0.29e-6 (uS/mV)     : amplitude for depression
	Ap = 0.28e-8 (uS/mV/mV)  : amplitude for potentiation
	thetam = -70.6 (mV) : threshold for depression, theta_minus
	thetap = -45.3 (mV)   : threshold for potentiation, theta_plus

	wmin = 0 (uS)
	wmax = 0.5 (uS)

}

ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	vca (mV)

	factor
}

STATE {
	A (uS)
	B (uS)

	u1 (mV)
	u2 (mV)
	u3
}

INITIAL {
	LOCAL tp
	if (tau1/tau2 > 0.9999) {
		tau1 = 0.9999*tau2
	}
	if (tau1/tau2 < 1e-9) {
		tau1 = tau2*1e-9
	}

	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor


	A = 0
	B = 0
	g = 0

	vca = v
	u1 = v
	u2 = v
	u3 = 0

	: Start monitoring a postsynaptic spike
	net_send(0, 1)
}

BREAKPOINT {
	SOLVE state METHOD cnexp

	g = B - A
	i = g*(v - e)
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2

	vca = (v - vca)/tauca   : delayed v
	u1' = (vca - u1)/taum   : u_- low-pass filtered vca
	u2' = (vca - u2)/taup   : u_+ low-pass filtered vca
	u3' = Ap*rectified(v-thetap)*rectified(u2-thetam) : potentiation eligibility
}

NET_RECEIVE(w (uS), x, tpre (ms)) {
	INITIAL { x=0 tpre=-1e9}

	if (flag==0) { : Presynaptic spike
		: printf("entry flag=%g t=%g w=%g x=%g tpre=%g\n", flag, t, w, x, tpre)

		: Activate the synapse
		A = A + w*factor : TODO can add w*factor
		B = B + w*factor : TODO can add w*factor
		:printf("entry flag=%g t=%g A=%g B=%g\n", flag, t, A, B)

		: Decay of the spike count
		x = exp(-(t-tpre)/taux)*x

:		if (spiking==0) { : is a cell spiking
		: otherwise presynaptic spike induces depression
		: printf("DEP t=%g w=%g dw=%g x=%g tpre=%g\n", t, w, -Am*rectified(u1-thetam), x, tpre)
		w = w- Am*rectified(u1-thetam)
		: But not below wmin
		if (w<wmin) { w = wmin }
:		}

		: Accumulate presynaptic spike count
		x = x + 1
		tpre = t

	} else if (flag==2) { : Postsynaptic spike begin
		: printf("entry flag=%g t=%g\n", flag, t)
		: Wait until the spike ends
		WATCH (v<thetap) 3

	} else if (flag==3) { : Postsynaptic spike end

		: printf("entry flag=%g t=%g\n", flag, t)

		: Potentiate all netcons now
		FOR_NETCONS(wi, xi, tprei) {
			: printf("POT t=%g w=%g dw=%g, x=%g tpre=%g\n", t, wi, u3*xi*exp(-(t-tprei)/taux), xi, tprei)
			wi = wi + u3*xi*exp(-(t-tprei)/taux)

			: unless they are not too stong
			if (wi>wmax) { wi = wmax }
		}

		: Reset acculumated potentiation
		u3 = 0

		: Everything back to normal
		net_send(0, 1)

	} else {
		: watch for a postsynaptic spike
		WATCH (v>thetap)  2
	}
}

UNITSOFF
FUNCTION rectified (v (mV)) {
	if (v>0) {
		rectified = v
	} else {
		rectified = 0
	}
}

UNITSON
