""" master mix components"""
##components and dilution factors

##Basic buffermixes
wilham = {"5XPBS":0.2, "2M NaCl":0.085, "0.1M EDTA":0.01, "1mM ThT":0.01}
salvadores = {"1M Tris-HCl pH7.4": 0.1, "1mM ThT": 0.005}

##Surfactants
surfactants = ["F-127", "F-68", "PE-P84", "F108"]
surfact_concs = [0.001, 0.0011, 0.005]

##substrate protein batches and dilution factors
sb_dil_fact = {"HaFLPrP 'M'":0.222,
               "abeta#1":0.1}

##physical_parameters of various assays
rtquic =    {"plate_shaker":"BMG-floustar/omega","rpm":900, "config":"double orbital", "shake on":87, "shake off": 33,   "read every":15, "plate_type":96, "reaction_vol":100, "max_time": 100, "temp":42, "excite":435, "emit": 485}
abetaPMCA = {"plate_shaker":"Thermomixer",       "rpm":500, "config":"double orbital", "shake on":60, "shake off": 1740, "read every":60, "plate_type":96, "reaction_vol":100, "max_time": 200, "temp":22, "excite":435, "emit": 485}

##molecular_weights
mw = {"HamFLPrP": 22803.23,
      "abeta_42":4514.10,
      "abeta_40:":4329.86}

class Master_mix(object):
    def __init__(self, mm_components, numReps, react_vol, seed_vol, sub_batch):
        self.components = mm_components
        self.numReps = numReps
        self.react_vol = react_vol
        self.seed_vol = seed_vol
        self.sub_batch = sub_batchd
        self.mm_vol = self.water_vol = ((self.react_vol-self.seed_vol)\
                                        * numReps * 1.1)
        self.recipe = {}
        self.recipe["Seed vol"] = self.seed_vol
        for name, dil_fact in self.components.items():
            self.recipe[name] = dil_fact * self.mm_vol
            self.water_vol -= self.recipe[name]
        self.recipe[sub_batch] = sb_dil_fact[sub_batch]*self.mm_vol           
        self.water_vol -= self.recipe[sub_batch]
        assert self.water_vol > 0
        self.recipe["Water"] = self.water_vol
    def __str__(self):
        readout ="MM vol..... "+str(round(self.mm_vol,2))+"\n"
        readout += "No. reps..... "+str(self.numReps)+"\n"
        for name, vol in self.recipe.items():
            readout += (name+"..... "+str(round(vol, 2))+ "\n")
        return readout

##rtquic1 = Master_mix(wilham_mm, 12, 100, 2, "HaFLPrP 'M'")
##abeta = Master_mix(salvadores_mm, 12, 200, 40, "abeta#1")
##print(abetaPMCA)
            
class Phys_params(object):
     pass

class RTQuICplate(object):
    row = ["A","B","C","D","E","F","G","H"]
    def __init__(self, seeds, buffers, substrates):
        self.seeds = seeds 
        self.buffers = buffers
        self.substrates = substrates
        self.mastermixes = []
        ##mastermixes of each buffer comp will be split into submastermixes
        ##for each additive dilution
        self.subMastermixes = []
        self.plate = {}
        self.setPlate()
        self.setMasterMixes()
        
    def setPlate(self):
        i = 0
        if i < 96:
            for sub in self.substrates:           
                for buff in self.buffers:
                    for dil in buff.getDil():
                        for seed in self.seeds:
                            for r in range(seed.get_numReps()):
                                self.plate[RTQuICplate.row[i//12] + \
str((i%12)+1)+ ": "] = sub.get_name() + buff.get_name() + " " + \
str(dil) + seed.get_name()+ " " + str(seed.get_seedVol()) +"ul"
                                i += 1
        else:
            raise ValueError('full')
        
    def setMasterMixes(self, seed_vol = 2, additive_vol = 10, well_vol = 100):
        numReps = 0 
        for seed in self.seeds:
            numReps += len(self.seeds) * len(self.buffers) * len(self.substrates)
        for sub in self.substrates:           
            mm_vols = {}
            mm_vols["Seed vol (ul per well)"] = seed_vol
            mm_vols["MM Vol before additive"] = \
(well_vol - seed_vol - additive_vol)* numReps * 1.2 
            water_vol = mm_vols["MM Vol before additive"] 
            mm_vols[sub.get_name()+" vol"] = sub.get_dilFac() * numReps *  1.2 * well_vol
            water_vol -= mm_vols[sub.get_name()+" vol"]
            print (buff.get_bufferBase())
            for buff_part, fact in buff.get_bufferBase().items():
                    mm_vols[buff_part] = fact * numReps * 1.2 * well_vol
                    water_vol -= mm_vols[buff_part] 
            mm_vols["Water"] = water_vol 
            self.mastermixes.append(mm_vols)
            for buff in self.buffers:
                for dil in buff.getDil():
                    subMM = {}
                    subMM["sMM Vol before additive"]= (well_vol - seed_vol - additive_vol) * (numReps/len(self.buffers)) *1.1  
                subMM["Vol Additive stock"] = additive_vol * (numReps/len(self.buffers)) *1.1 * (dil/buff.get_stock())
                subMM["Top up water vol"] = additive_vol - subMM["Vol Additive stock"]
                self.subMastermixes.append(subMM)
                    
    def get_mastermixes(self):
        print("Mastermixes and submastermixes")
        print("Mastermixes")
        for mix in self.mastermixes:
            for key, value in mix.items():
                print(key + ": " + str(round(value,1)))
            print ("\n=========\n")
        print("Submastermixes")
        for mix in self.subMastermixes:
            for key, value in mix.items():
                print(key + ": " + str(round(value,1)))
            print ("\n=========\n")
        return self.mastermixes, self.subMastermixes
    def getPlate(self):
        return self.plate
    def __str__(self):
        string = "Plate setup\n"
        for well, contents in self.plate.items():
            string += (well + contents + "\n")
        return string

class Component(object):
    def __init__(self, name):
        self.name = name
    def get_name(self):
        return self.name
        

class Seed (Component):
    def __init__(self, name, numReps, seedVol):
        Component.__init__(self, name)
        self.numReps = numReps
        self.seedVol = seedVol
    def get_numReps(self):
        return self.numReps
    def get_seedVol(self):
        return self.seedVol                                  


class Substrate(Component):
    def __init__(self, name, dilFac):
        Component.__init__(self, name)
        self.dilFac = dilFac
    def get_dilFac(self):
        return self.dilFac


class Additive(Component):
    def __init__(self, name, stock_conc, concs):
        Component.__init__(self, name)
        self.stock_conc = stock_conc
        self.concs = concs
    def get_stock_conc(self):
        return self.stock_conc
    def get_concs(self):
        return self.concs



##for surfactant in surfactants:
##    for conc in surfactant_concs:
##        component[ 

class Buffer(Component):
    def __init__(self, name, bufferBase, additive, additive_dilutions, additive_stock):
        Component.__init__(self, name)
        self.bufferBase  = bufferBase
        self.additive  = additive
        self.additive_dilutions  = additive_dilutions
        self.additive_stock = additive_stock
    def get_bufferBase(self):
        return self.bufferBase
    def getDil(self):
        return self.additive_dilutions
    def get_stock(self):
        return self.additive_stock

buffers = []

for surf in surfactants[:2]:
    buff = Buffer(surf, wilham, surf, surfact_concs, 0.005)
    buffers.append(buff)


seeds_M = [Seed("Water",9,2), Seed("MM1 10-3",3,2)]
substrates_M = [Substrate("HSFLPRP",0.1)]
ex5 = RTQuICplate(seeds_M, buffers, substrates_M)
print(ex5)
print(ex5.get_mastermixes())
    
    
