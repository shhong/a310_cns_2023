
from neuron import h
from cell_template import Cell

class Mainen(Cell):
    def __init__(self, i, rho=150, c=10):
        # We initialize the cell as in the `Cell` class.
        super().__init__()
        h.v_init = -70
        h.celsius = 37

        self.id = i
        self.set_rho(rho)
        self.set_c(c)

    def create_sections(self):
        self.soma = h.Section(name="soma", cell=self)
        self.dend = h.Section(name="dend", cell=self)

    def build_topology(self):
        self.dend.connect(self.soma(0.5), 0)

    def build_subsets(self):
        pass

    def define_geometry(self):
        self.soma.diam = 10/h.PI
        self.soma.L = 10
        self.dend.diam = 10/h.PI
        self.dend.L = 200

    def define_biophysics(self):
        self.soma.insert("na")
        self.soma.insert("kv")

        self.dend.insert("pas")
        self.dend.insert("na")
        self.dend.insert("ca")
        self.dend.insert("km")
        self.dend.insert("kca")
        self.dend.insert("cad")

        self.soma.ek = -90
        self.soma.ena = 60
        self.soma.gbar_na = 30000
        self.soma.gbar_kv = 1500

        self.dend.ek = -90
        self.dend.ena = 60
        self.dend.eca = 140
        h.ion_style("ca_ion", 0, 1, 0, 0, 0, sec=self.dend)

        self.dend.g_pas = 1/30000
        self.dend.gbar_na = 15
        self.dend.gbar_ca = 0.3
        self.dend.gbar_kca = 3
        self.dend.gbar_km = 0.1

        self.soma.cm = 0.75
        self.dend.cm = 0.75

    def set_rho(self, rho):
        self.rho = rho
        self.dend.L = self.rho*self.soma.L

    def set_c(self, c):
        self.c = c
        self.dend.Ra = self.dend.Ra*self.c/h.ri(.5, sec=self.dend)

    def create_synapses(self):
        self.synlist.append(h.Exp2Syn(self.dend(0.5))) # AMPA
        self.synlist[-1].tau1 = 0.1
        self.synlist[-1].tau2 = 2.0
        self.synlist[-1].e = 0

        self.synlist.append(h.Exp2Syn(self.dend(0.5))) # NMDA-like
        self.synlist[-1].tau1 = 0.5
        self.synlist[-1].tau2 = 20.0
        self.synlist[-1].e = 0

        self.synlist.append(h.Exp2Syn(self.dend(0.5))) # GABA
        self.synlist[-1].e = -75
        self.synlist[-1].tau1 = 0.5
        self.synlist[-1].tau2 = 4.0

        self.synlist.append(h.Exp2Syn(self.soma(0.5))) # GABA
        self.synlist[-1].e = -75
        self.synlist[-1].tau1 = 0.5
        self.synlist[-1].tau2 = 4.0

class Mease(Cell):
    def __init__(self, i):
        super().__init__()
        h.v_init = -70
        h.celsius = 37

        self.id = i
        
    def create_sections(self):
        self.soma = h.Section(name="soma", cell=self)

    def build_topology(self):
        pass

    def build_subsets(self):
        pass

    def define_geometry(self):
        self.soma.diam = 15/h.PI
        self.soma.L = 15

    def define_biophysics(self):
        self.soma.insert("mease")
        # self.soma.insert("hh")
        self.soma.insert("pas")
        self.set_gpas(0*1e-5)
        self.set_epas(10)
        
    def set_gpas(self, gpas):
        self.soma.g_pas = gpas
        
    def set_epas(self, epas):
        self.soma.e_pas = epas

    def create_synapses(self):
        self.synlist.append(h.Exp2Syn(self.soma(0.5))) # AMPA
        self.synlist[-1].tau1 = 0.1
        self.synlist[-1].tau2 = 20.0
        self.synlist[-1].e = 0

        # self.synlist.append(h.Exp2Syn(self.soma(0.5))) # NMDA-like
        # self.synlist[-1].tau1 = 1.0
        # self.synlist[-1].tau2 = 5.0
        # self.synlist[-1].e = 0

        # self.synlist.append(h.Exp2Syn(self.soma(0.5))) # GABA
        # self.synlist[-1].e = -75
        # self.synlist[-1].tau1 = 0.5
        # self.synlist[-1].tau2 = 5.0
