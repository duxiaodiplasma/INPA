import pickle

with open('22mm.pickle') as f:
     ini22,geo22,res22 = pickle.load(f)

with open('55mm.pickle') as f:
     ini55,geo55,res55 = pickle.load(f)

with open('12mm.pickle') as f:
     ini12,geo12,res12 = pickle.load(f)

with open('12mm_halffoil.pickle') as f:
     ini12h,geo12h,res12h = pickle.load(f)

with open('55mm_halffoil.pickle') as f:
     ini55h,geo55h,res55h = pickle.load(f)
