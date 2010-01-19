# This file was created automatically by SWIG.
import MEADc
class Coord:
    def __init__(self,*args):
        self.this = apply(MEADc.new_Coord,args)
        self.thisown = 1

    def __del__(self,MEADc=MEADc):
        if self.thisown == 1 :
            MEADc.delete_Coord(self)
    def write(*args):
        val = apply(MEADc.Coord_write,args)
        return val
    def cross(*args):
        val = apply(MEADc.Coord_cross,args)
        if val: val = CoordPtr(val) ; val.thisown = 1
        return val
    def dot(*args):
        val = apply(MEADc.Coord_dot,args)
        return val
    def __iadd__(*args):
        val = apply(MEADc.Coord___iadd__,args)
        if val: val = CoordPtr(val) ; val.thisown = 1
        return val
    def __isub__(*args):
        val = apply(MEADc.Coord___isub__,args)
        if val: val = CoordPtr(val) ; val.thisown = 1
        return val
    def __imul__(*args):
        val = apply(MEADc.Coord___imul__,args)
        if val: val = CoordPtr(val) ; val.thisown = 1
        return val
    def __idiv__(*args):
        val = apply(MEADc.Coord___idiv__,args)
        if val: val = CoordPtr(val) ; val.thisown = 1
        return val
    def __neg__(*args):
        val = apply(MEADc.Coord___neg__,args)
        if val: val = CoordPtr(val) ; val.thisown = 1
        return val
    def __gt__(*args):
        val = apply(MEADc.Coord___gt__,args)
        return val
    def __lt__(*args):
        val = apply(MEADc.Coord___lt__,args)
        return val
    def __add__(*args):
        val = apply(MEADc.Coord___add__,args)
        if val: val = CoordPtr(val) ; val.thisown = 1
        return val
    def __sub__(*args):
        val = apply(MEADc.Coord___sub__,args)
        if val: val = CoordPtr(val) ; val.thisown = 1
        return val
    def __mul__(*args):
        val = apply(MEADc.Coord___mul__,args)
        if val: val = CoordPtr(val) ; val.thisown = 1
        return val
    def __rmul__(*args):
        val = apply(MEADc.Coord___rmul__,args)
        if val: val = CoordPtr(val) ; val.thisown = 1
        return val
    def __mul____L86(*args):
        val = apply(MEADc.Coord___mul____L86,args)
        return val
    def __div__(*args):
        val = apply(MEADc.Coord___div__,args)
        if val: val = CoordPtr(val) ; val.thisown = 1
        return val
    def __eq__(*args):
        val = apply(MEADc.Coord___eq__,args)
        return val
    def __ne__(*args):
        val = apply(MEADc.Coord___ne__,args)
        return val
    def __pos__(*args):
        val = apply(MEADc.Coord___pos__,args)
        if val: val = CoordPtr(val) ; val.thisown = 1
        return val
    def __cmp__(*args):
        val = apply(MEADc.Coord___cmp__,args)
        return val
    __setmethods__ = {
        "x" : MEADc.Coord_x_set,
        "y" : MEADc.Coord_y_set,
        "z" : MEADc.Coord_z_set,
    }
    def __setattr__(self,name,value):
        if (name == "this") or (name == "thisown"): self.__dict__[name] = value; return
        method = Coord.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value
    __getmethods__ = {
        "x" : MEADc.Coord_x_get,
        "y" : MEADc.Coord_y_get,
        "z" : MEADc.Coord_z_get,
    }
    def __getattr__(self,name):
        method = Coord.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C Coord instance at %s>" % (self.this,)
class CoordPtr(Coord):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = Coord






class ChargeDist_lett:
    def __init__(self,this):
        self.this = this

    def __repr__(self):
        return "<C ChargeDist_lett instance at %s>" % (self.this,)
class ChargeDist_lettPtr(ChargeDist_lett):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = ChargeDist_lett



class DielectricEnvironment_lett:
    def __init__(self,this):
        self.this = this

    def get_cuberep(*args):
        val = apply(MEADc.DielectricEnvironment_lett_get_cuberep,args)
        return val
    def __repr__(self):
        return "<C DielectricEnvironment_lett instance at %s>" % (self.this,)
class DielectricEnvironment_lettPtr(DielectricEnvironment_lett):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = DielectricEnvironment_lett



class ElectrolyteEnvironment_lett:
    def __init__(self,this):
        self.this = this

    def get_cuberep(*args):
        val = apply(MEADc.ElectrolyteEnvironment_lett_get_cuberep,args)
        return val
    def __repr__(self):
        return "<C ElectrolyteEnvironment_lett instance at %s>" % (self.this,)
class ElectrolyteEnvironment_lettPtr(ElectrolyteEnvironment_lett):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = ElectrolyteEnvironment_lett



class ElstatPot_lett:
    def __init__(self,this):
        self.this = this

    def get_cuberep(*args):
        val = apply(MEADc.ElstatPot_lett_get_cuberep,args)
        return val
    def __add__(*args):
        val = apply(MEADc.ElstatPot_lett___add__,args)
        if val: val = ElstatPotCombinationPtr(val) ; val.thisown = 1
        return val
    def __add____L80(*args):
        val = apply(MEADc.ElstatPot_lett___add____L80,args)
        if val: val = ElstatPotCombinationPtr(val) ; val.thisown = 1
        return val
    def __repr__(self):
        return "<C ElstatPot_lett instance at %s>" % (self.this,)
class ElstatPot_lettPtr(ElstatPot_lett):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = ElstatPot_lett



class AnalyticEP(ElstatPot_lett):
    def __init__(self,this):
        self.this = this

    def __repr__(self):
        return "<C AnalyticEP instance at %s>" % (self.this,)
class AnalyticEPPtr(AnalyticEP):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = AnalyticEP



class Atom:
    def __init__(self,*args):
        self.this = apply(MEADc.new_Atom,args)
        self.thisown = 1

    def __del__(self,MEADc=MEADc):
        if self.thisown == 1 :
            MEADc.delete_Atom(self)
    def write(*args):
        val = apply(MEADc.Atom_write,args)
        return val
    def __eq__(*args):
        val = apply(MEADc.Atom___eq__,args)
        return val
    def __ne__(*args):
        val = apply(MEADc.Atom___ne__,args)
        return val
    def __lt__(*args):
        val = apply(MEADc.Atom___lt__,args)
        return val
    def __gt__(*args):
        val = apply(MEADc.Atom___gt__,args)
        return val
    def __cmp__(*args):
        val = apply(MEADc.Atom___cmp__,args)
        return val
    __setmethods__ = {
        "coord" : MEADc.Atom_coord_set,
        "rad" : MEADc.Atom_rad_set,
        "charge" : MEADc.Atom_charge_set,
        "atname" : MEADc.Atom_atname_set,
        "resname" : MEADc.Atom_resname_set,
        "chainid" : MEADc.Atom_chainid_set,
        "resnum" : MEADc.Atom_resnum_set,
    }
    def __setattr__(self,name,value):
        if (name == "this") or (name == "thisown"): self.__dict__[name] = value; return
        method = Atom.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value
    __getmethods__ = {
        "coord" : lambda x : CoordPtr(MEADc.Atom_coord_get(x)),
        "rad" : MEADc.Atom_rad_get,
        "charge" : MEADc.Atom_charge_get,
        "atname" : MEADc.Atom_atname_get,
        "resname" : MEADc.Atom_resname_get,
        "chainid" : MEADc.Atom_chainid_get,
        "resnum" : MEADc.Atom_resnum_get,
    }
    def __getattr__(self,name):
        method = Atom.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C Atom instance at %s>" % (self.this,)
class AtomPtr(Atom):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = Atom





class AtomID:
    def __init__(self,*args):
        self.this = apply(MEADc.new_AtomID,args)
        self.thisown = 1

    def __del__(self,MEADc=MEADc):
        if self.thisown == 1 :
            MEADc.delete_AtomID(self)
    def __eq__(*args):
        val = apply(MEADc.AtomID___eq__,args)
        return val
    def __ne__(*args):
        val = apply(MEADc.AtomID___ne__,args)
        return val
    def __lt__(*args):
        val = apply(MEADc.AtomID___lt__,args)
        return val
    def __gt__(*args):
        val = apply(MEADc.AtomID___gt__,args)
        return val
    def __cmp__(*args):
        val = apply(MEADc.AtomID___cmp__,args)
        return val
    __setmethods__ = {
        "resnum" : MEADc.AtomID_resnum_set,
        "atname" : MEADc.AtomID_atname_set,
        "chainid" : MEADc.AtomID_chainid_set,
    }
    def __setattr__(self,name,value):
        if (name == "this") or (name == "thisown"): self.__dict__[name] = value; return
        method = AtomID.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value
    __getmethods__ = {
        "resnum" : MEADc.AtomID_resnum_get,
        "atname" : MEADc.AtomID_atname_get,
        "chainid" : MEADc.AtomID_chainid_get,
    }
    def __getattr__(self,name):
        method = AtomID.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C AtomID instance at %s>" % (self.this,)
class AtomIDPtr(AtomID):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = AtomID








class AtomSet:
    def __init__(self,*args):
        self.this = apply(MEADc.new_AtomSet,args)
        self.thisown = 1

    def __del__(self,MEADc=MEADc):
        if self.thisown == 1 :
            MEADc.delete_AtomSet(self)
    def read(*args):
        val = apply(MEADc.AtomSet_read,args)
        return val
    def read__L21(*args):
        val = apply(MEADc.AtomSet_read__L21,args)
        return val
    def geom_cent(*args):
        val = apply(MEADc.AtomSet_geom_cent,args)
        if val: val = CoordPtr(val) ; val.thisown = 1
        return val
    def set_coords_to(*args):
        val = apply(MEADc.AtomSet_set_coords_to,args)
        return val
    def build_from_vectors(*args):
        val = apply(MEADc.AtomSet_build_from_vectors,args)
        return val
    def insert(*args):
        val = apply(MEADc.AtomSet_insert,args)
        return val
    def clear(*args):
        val = apply(MEADc.AtomSet_clear,args)
        return val
    def copy(*args):
        val = apply(MEADc.AtomSet_copy,args)
        if val: val = AtomSetPtr(val) 
        return val
    def has_key(*args):
        val = apply(MEADc.AtomSet_has_key,args)
        return val
    def __len__(*args):
        val = apply(MEADc.AtomSet___len__,args)
        return val
    def __delitem__(*args):
        val = apply(MEADc.AtomSet___delitem__,args)
        return val
    def update(*args):
        val = apply(MEADc.AtomSet_update,args)
        return val
    def keys(*args):
        val = apply(MEADc.AtomSet_keys,args)
        return val
    def values(*args):
        val = apply(MEADc.AtomSet_values,args)
        return val
    def __getitem__(*args):
        val = apply(MEADc.AtomSet___getitem__,args)
        if val: val = AtomPtr(val) 
        return val
    def __setitem__(*args):
        val = apply(MEADc.AtomSet___setitem__,args)
        return val
    __setmethods__ = {
        "name" : MEADc.AtomSet_name_set,
    }
    def __setattr__(self,name,value):
        if (name == "this") or (name == "thisown"): self.__dict__[name] = value; return
        method = AtomSet.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value
    __getmethods__ = {
        "name" : MEADc.AtomSet_name_get,
    }
    def __getattr__(self,name):
        method = AtomSet.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C AtomSet instance at %s>" % (self.this,)
class AtomSetPtr(AtomSet):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = AtomSet







class PointCharge:
    def __init__(self,*args):
        self.this = apply(MEADc.new_PointCharge,args)
        self.thisown = 1

    def __del__(self,MEADc=MEADc):
        if self.thisown == 1 :
            MEADc.delete_PointCharge(self)
    def __eq__(*args):
        val = apply(MEADc.PointCharge___eq__,args)
        return val
    __setmethods__ = {
        "charge" : MEADc.PointCharge_charge_set,
        "coord" : MEADc.PointCharge_coord_set,
    }
    def __setattr__(self,name,value):
        if (name == "this") or (name == "thisown"): self.__dict__[name] = value; return
        method = PointCharge.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value
    __getmethods__ = {
        "charge" : MEADc.PointCharge_charge_get,
        "coord" : lambda x : CoordPtr(MEADc.PointCharge_coord_get(x)),
    }
    def __getattr__(self,name):
        method = PointCharge.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C PointCharge instance at %s>" % (self.this,)
class PointChargePtr(PointCharge):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = PointCharge





class AtomChargeSet(ChargeDist_lett,AtomSet):
    def __init__(self,*args):
        self.this = apply(MEADc.new_AtomChargeSet,args)
        self.thisown = 1

    def __del__(self,MEADc=MEADc):
        if self.thisown == 1 :
            MEADc.delete_AtomChargeSet(self)
    def total_charge(*args):
        val = apply(MEADc.AtomChargeSet_total_charge,args)
        return val
    def has_charges(*args):
        val = apply(MEADc.AtomChargeSet_has_charges,args)
        return val
    def vacuum_coulomb(*args):
        val = apply(MEADc.AtomChargeSet_vacuum_coulomb,args)
        return val
    def number_points(*args):
        val = apply(MEADc.AtomChargeSet_number_points,args)
        return val
    def __mul__(*args):
        val = apply(MEADc.AtomChargeSet___mul__,args)
        if val: val = AtomChargeSetPtr(val) 
        return val
    def __rmul__(*args):
        val = apply(MEADc.AtomChargeSet___rmul__,args)
        if val: val = AtomChargeSetPtr(val) 
        return val
    def __div__(*args):
        val = apply(MEADc.AtomChargeSet___div__,args)
        if val: val = AtomChargeSetPtr(val) 
        return val
    def __add__(*args):
        val = apply(MEADc.AtomChargeSet___add__,args)
        if val: val = AtomChargeSetPtr(val) 
        return val
    def __add____L68(*args):
        val = apply(MEADc.AtomChargeSet___add____L68,args)
        if val: val = ManyPointChargePtr(val) 
        return val
    def __add____L69(*args):
        val = apply(MEADc.AtomChargeSet___add____L69,args)
        if val: val = ManyPointChargePtr(val) 
        return val
    def pointcharges(*args):
        val = apply(MEADc.AtomChargeSet_pointcharges,args)
        return val
    def __mul____L113(*args):
        val = apply(MEADc.AtomChargeSet___mul____L113,args)
        return val
    def __mul____L114(*args):
        val = apply(MEADc.AtomChargeSet___mul____L114,args)
        return val
    __setmethods__ = {
        "name" : MEADc.AtomChargeSet_name_set,
    }
    def __setattr__(self,name,value):
        if (name == "this") or (name == "thisown"): self.__dict__[name] = value; return
        method = AtomChargeSet.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value
    __getmethods__ = {
        "name" : MEADc.AtomChargeSet_name_get,
    }
    def __getattr__(self,name):
        method = AtomChargeSet.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C AtomChargeSet instance at %s>" % (self.this,)
class AtomChargeSetPtr(AtomChargeSet):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = AtomChargeSet






class ManyPointCharge(ChargeDist_lett):
    def __init__(self,*args):
        self.this = apply(MEADc.new_ManyPointCharge,args)
        self.thisown = 1

    def __del__(self,MEADc=MEADc):
        if self.thisown == 1 :
            MEADc.delete_ManyPointCharge(self)
    def total_charge(*args):
        val = apply(MEADc.ManyPointCharge_total_charge,args)
        return val
    def has_charges(*args):
        val = apply(MEADc.ManyPointCharge_has_charges,args)
        return val
    def vacuum_coulomb(*args):
        val = apply(MEADc.ManyPointCharge_vacuum_coulomb,args)
        return val
    def number_points(*args):
        val = apply(MEADc.ManyPointCharge_number_points,args)
        return val
    def __mul__(*args):
        val = apply(MEADc.ManyPointCharge___mul__,args)
        if val: val = ManyPointChargePtr(val) 
        return val
    def __rmul__(*args):
        val = apply(MEADc.ManyPointCharge___rmul__,args)
        if val: val = ManyPointChargePtr(val) 
        return val
    def __div__(*args):
        val = apply(MEADc.ManyPointCharge___div__,args)
        if val: val = ManyPointChargePtr(val) 
        return val
    def __add__(*args):
        val = apply(MEADc.ManyPointCharge___add__,args)
        if val: val = ManyPointChargePtr(val) 
        return val
    def __add____L55(*args):
        val = apply(MEADc.ManyPointCharge___add____L55,args)
        if val: val = ManyPointChargePtr(val) 
        return val
    def __add____L56(*args):
        val = apply(MEADc.ManyPointCharge___add____L56,args)
        if val: val = ManyPointChargePtr(val) 
        return val
    def pointcharges(*args):
        val = apply(MEADc.ManyPointCharge_pointcharges,args)
        return val
    def append(*args):
        val = apply(MEADc.ManyPointCharge_append,args)
        return val
    def __len__(*args):
        val = apply(MEADc.ManyPointCharge___len__,args)
        return val
    def extend(*args):
        val = apply(MEADc.ManyPointCharge_extend,args)
        return val
    def count(*args):
        val = apply(MEADc.ManyPointCharge_count,args)
        return val
    def index(*args):
        val = apply(MEADc.ManyPointCharge_index,args)
        return val
    def insert(*args):
        val = apply(MEADc.ManyPointCharge_insert,args)
        return val
    def remove(*args):
        val = apply(MEADc.ManyPointCharge_remove,args)
        return val
    def __getitem__(*args):
        val = apply(MEADc.ManyPointCharge___getitem__,args)
        if val: val = PointChargePtr(val) ; val.thisown = 1
        return val
    def __setitem__(*args):
        val = apply(MEADc.ManyPointCharge___setitem__,args)
        return val
    def __delitem__(*args):
        val = apply(MEADc.ManyPointCharge___delitem__,args)
        return val
    def __mul____L119(*args):
        val = apply(MEADc.ManyPointCharge___mul____L119,args)
        return val
    def __mul____L120(*args):
        val = apply(MEADc.ManyPointCharge___mul____L120,args)
        return val
    def __repr__(self):
        return "<C ManyPointCharge instance at %s>" % (self.this,)
class ManyPointChargePtr(ManyPointCharge):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = ManyPointCharge





class OnePointCharge(ChargeDist_lett):
    def __init__(self,*args):
        self.this = apply(MEADc.new_OnePointCharge,args)
        self.thisown = 1

    def __del__(self,MEADc=MEADc):
        if self.thisown == 1 :
            MEADc.delete_OnePointCharge(self)
    def vacuum_coulomb(*args):
        val = apply(MEADc.OnePointCharge_vacuum_coulomb,args)
        return val
    def total_charge(*args):
        val = apply(MEADc.OnePointCharge_total_charge,args)
        return val
    def has_charges(*args):
        val = apply(MEADc.OnePointCharge_has_charges,args)
        return val
    def number_points(*args):
        val = apply(MEADc.OnePointCharge_number_points,args)
        return val
    def get_charge(*args):
        val = apply(MEADc.OnePointCharge_get_charge,args)
        return val
    def get_coord(*args):
        val = apply(MEADc.OnePointCharge_get_coord,args)
        if val: val = CoordPtr(val) ; val.thisown = 1
        return val
    def pointcharge(*args):
        val = apply(MEADc.OnePointCharge_pointcharge,args)
        if val: val = PointChargePtr(val) ; val.thisown = 1
        return val
    def __mul__(*args):
        val = apply(MEADc.OnePointCharge___mul__,args)
        if val: val = OnePointChargePtr(val) 
        return val
    def __rmul__(*args):
        val = apply(MEADc.OnePointCharge___rmul__,args)
        if val: val = OnePointChargePtr(val) 
        return val
    def __div__(*args):
        val = apply(MEADc.OnePointCharge___div__,args)
        if val: val = OnePointChargePtr(val) 
        return val
    def __add__(*args):
        val = apply(MEADc.OnePointCharge___add__,args)
        if val: val = ManyPointChargePtr(val) 
        return val
    def __add____L66(*args):
        val = apply(MEADc.OnePointCharge___add____L66,args)
        if val: val = ManyPointChargePtr(val) 
        return val
    def __add____L69(*args):
        val = apply(MEADc.OnePointCharge___add____L69,args)
        if val: val = ManyPointChargePtr(val) 
        return val
    def __mul____L123(*args):
        val = apply(MEADc.OnePointCharge___mul____L123,args)
        return val
    def __mul____L124(*args):
        val = apply(MEADc.OnePointCharge___mul____L124,args)
        return val
    def __repr__(self):
        return "<C OnePointCharge instance at %s>" % (self.this,)
class OnePointChargePtr(OnePointCharge):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = OnePointCharge




class UniformDielectric(DielectricEnvironment_lett):
    def __init__(self,*args):
        self.this = apply(MEADc.new_UniformDielectric,args)
        self.thisown = 1

    def __del__(self,MEADc=MEADc):
        if self.thisown == 1 :
            MEADc.delete_UniformDielectric(self)
    def epsext_value(*args):
        val = apply(MEADc.UniformDielectric_epsext_value,args)
        return val
    def __repr__(self):
        return "<C UniformDielectric instance at %s>" % (self.this,)
class UniformDielectricPtr(UniformDielectric):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = UniformDielectric



class DielectricSphere(DielectricEnvironment_lett):
    def __init__(self,*args):
        self.this = apply(MEADc.new_DielectricSphere,args)
        self.thisown = 1

    def __del__(self,MEADc=MEADc):
        if self.thisown == 1 :
            MEADc.delete_DielectricSphere(self)
    def epsext_value(*args):
        val = apply(MEADc.DielectricSphere_epsext_value,args)
        return val
    def radius_value(*args):
        val = apply(MEADc.DielectricSphere_radius_value,args)
        return val
    def epsin_value(*args):
        val = apply(MEADc.DielectricSphere_epsin_value,args)
        return val
    def get_center(*args):
        val = apply(MEADc.DielectricSphere_get_center,args)
        if val: val = CoordPtr(val) ; val.thisown = 1
        return val
    def __repr__(self):
        return "<C DielectricSphere instance at %s>" % (self.this,)
class DielectricSpherePtr(DielectricSphere):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = DielectricSphere



class DielectricSlab(DielectricEnvironment_lett):
    def __init__(self,*args):
        self.this = apply(MEADc.new_DielectricSlab,args)
        self.thisown = 1

    def __del__(self,MEADc=MEADc):
        if self.thisown == 1 :
            MEADc.delete_DielectricSlab(self)
    def epsslab_value(*args):
        val = apply(MEADc.DielectricSlab_epsslab_value,args)
        return val
    def epsext_value(*args):
        val = apply(MEADc.DielectricSlab_epsext_value,args)
        return val
    def zupper_value(*args):
        val = apply(MEADc.DielectricSlab_zupper_value,args)
        return val
    def zlower_value(*args):
        val = apply(MEADc.DielectricSlab_zlower_value,args)
        return val
    def __repr__(self):
        return "<C DielectricSlab instance at %s>" % (self.this,)
class DielectricSlabPtr(DielectricSlab):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = DielectricSlab



class TwoValueDielectricByAtoms(DielectricEnvironment_lett):
    def __init__(self,*args):
        self.this = apply(MEADc.new_TwoValueDielectricByAtoms,args)
        self.thisown = 1

    def __del__(self,MEADc=MEADc):
        if self.thisown == 1 :
            MEADc.delete_TwoValueDielectricByAtoms(self)
    def epsext_value(*args):
        val = apply(MEADc.TwoValueDielectricByAtoms_epsext_value,args)
        return val
    def epsin_value(*args):
        val = apply(MEADc.TwoValueDielectricByAtoms_epsin_value,args)
        return val
    def __repr__(self):
        return "<C TwoValueDielectricByAtoms instance at %s>" % (self.this,)
class TwoValueDielectricByAtomsPtr(TwoValueDielectricByAtoms):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = TwoValueDielectricByAtoms






class ThreeValueDielectricByAtoms(DielectricEnvironment_lett):
    def __init__(self,*args):
        self.this = apply(MEADc.new_ThreeValueDielectricByAtoms,args)
        self.thisown = 1

    def __del__(self,MEADc=MEADc):
        if self.thisown == 1 :
            MEADc.delete_ThreeValueDielectricByAtoms(self)
    def epsext_value(*args):
        val = apply(MEADc.ThreeValueDielectricByAtoms_epsext_value,args)
        return val
    def __repr__(self):
        return "<C ThreeValueDielectricByAtoms instance at %s>" % (self.this,)
class ThreeValueDielectricByAtomsPtr(ThreeValueDielectricByAtoms):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = ThreeValueDielectricByAtoms






class TwoValueDielMembAtoms(DielectricEnvironment_lett):
    def __init__(self,*args):
        self.this = apply(MEADc.new_TwoValueDielMembAtoms,args)
        self.thisown = 1

    def __del__(self,MEADc=MEADc):
        if self.thisown == 1 :
            MEADc.delete_TwoValueDielMembAtoms(self)
    def epsext_value(*args):
        val = apply(MEADc.TwoValueDielMembAtoms_epsext_value,args)
        return val
    def __repr__(self):
        return "<C TwoValueDielMembAtoms instance at %s>" % (self.this,)
class TwoValueDielMembAtomsPtr(TwoValueDielMembAtoms):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = TwoValueDielMembAtoms




class ThreeValueDielMembAtomsAtoms(DielectricEnvironment_lett):
    def __init__(self,*args):
        self.this = apply(MEADc.new_ThreeValueDielMembAtomsAtoms,args)
        self.thisown = 1

    def __del__(self,MEADc=MEADc):
        if self.thisown == 1 :
            MEADc.delete_ThreeValueDielMembAtomsAtoms(self)
    def epsext_value(*args):
        val = apply(MEADc.ThreeValueDielMembAtomsAtoms_epsext_value,args)
        return val
    def __repr__(self):
        return "<C ThreeValueDielMembAtomsAtoms instance at %s>" % (self.this,)
class ThreeValueDielMembAtomsAtomsPtr(ThreeValueDielMembAtomsAtoms):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = ThreeValueDielMembAtomsAtoms



class UniformElectrolyte(ElectrolyteEnvironment_lett):
    def __init__(self,*args):
        self.this = apply(MEADc.new_UniformElectrolyte,args)
        self.thisown = 1

    def __del__(self,MEADc=MEADc):
        if self.thisown == 1 :
            MEADc.delete_UniformElectrolyte(self)
    def ionic_strength(*args):
        val = apply(MEADc.UniformElectrolyte_ionic_strength,args)
        return val
    def __repr__(self):
        return "<C UniformElectrolyte instance at %s>" % (self.this,)
class UniformElectrolytePtr(UniformElectrolyte):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = UniformElectrolyte



class ElySphere(ElectrolyteEnvironment_lett):
    def __init__(self,*args):
        self.this = apply(MEADc.new_ElySphere,args)
        self.thisown = 1

    def __del__(self,MEADc=MEADc):
        if self.thisown == 1 :
            MEADc.delete_ElySphere(self)
    def ionic_strength(*args):
        val = apply(MEADc.ElySphere_ionic_strength,args)
        return val
    def get_radius(*args):
        val = apply(MEADc.ElySphere_get_radius,args)
        return val
    def get_center(*args):
        val = apply(MEADc.ElySphere_get_center,args)
        if val: val = CoordPtr(val) ; val.thisown = 1
        return val
    def __repr__(self):
        return "<C ElySphere instance at %s>" % (self.this,)
class ElySpherePtr(ElySphere):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = ElySphere



class ElectrolyteByAtoms(ElectrolyteEnvironment_lett):
    def __init__(self,*args):
        self.this = apply(MEADc.new_ElectrolyteByAtoms,args)
        self.thisown = 1

    def __del__(self,MEADc=MEADc):
        if self.thisown == 1 :
            MEADc.delete_ElectrolyteByAtoms(self)
    def ionic_strength(*args):
        val = apply(MEADc.ElectrolyteByAtoms_ionic_strength,args)
        return val
    def __repr__(self):
        return "<C ElectrolyteByAtoms instance at %s>" % (self.this,)
class ElectrolyteByAtomsPtr(ElectrolyteByAtoms):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = ElectrolyteByAtoms






class SolvAccVol:
    def __init__(self,*args):
        self.this = apply(MEADc.new_SolvAccVol,args)
        self.thisown = 1

    def __del__(self,MEADc=MEADc):
        if self.thisown == 1 :
            MEADc.delete_SolvAccVol(self)
    def anal_calc(*args):
        val = apply(MEADc.SolvAccVol_anal_calc,args)
        return val
    def accessible(*args):
        val = apply(MEADc.SolvAccVol_accessible,args)
        return val
    def check_is_calculated(*args):
        val = apply(MEADc.SolvAccVol_check_is_calculated,args)
        return val
    def write_top_in_binary(*args):
        val = apply(MEADc.SolvAccVol_write_top_in_binary,args)
        return val
    def write_top_in_ascii(*args):
        val = apply(MEADc.SolvAccVol_write_top_in_ascii,args)
        return val
    def read_top_in_binary(*args):
        val = apply(MEADc.SolvAccVol_read_top_in_binary,args)
        return val
    def read_top_in_ascii(*args):
        val = apply(MEADc.SolvAccVol_read_top_in_ascii,args)
        return val
    def write_top_in_binary__L30(*args):
        val = apply(MEADc.SolvAccVol_write_top_in_binary__L30,args)
        return val
    def write_top_in_ascii__L31(*args):
        val = apply(MEADc.SolvAccVol_write_top_in_ascii__L31,args)
        return val
    def read_top_in_binary__L32(*args):
        val = apply(MEADc.SolvAccVol_read_top_in_binary__L32,args)
        return val
    def read_top_in_ascii__L33(*args):
        val = apply(MEADc.SolvAccVol_read_top_in_ascii__L33,args)
        return val
    def get_cuberep(*args):
        val = apply(MEADc.SolvAccVol_get_cuberep,args)
        return val
    def tag_points(*args):
        val = apply(MEADc.SolvAccVol_tag_points,args)
        return val
    def __repr__(self):
        return "<C SolvAccVol instance at %s>" % (self.this,)
class SolvAccVolPtr(SolvAccVol):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = SolvAccVol






class AccTag_enum:
    interior = MEADc.AccTag_enum_interior
    exterior = MEADc.AccTag_enum_exterior
    undecided = MEADc.AccTag_enum_undecided
    in_tube = MEADc.AccTag_enum_in_tube
    def __init__(self,this):
        self.this = this

    def __repr__(self):
        return "<C AccTag_enum instance at %s>" % (self.this,)
class AccTag_enumPtr(AccTag_enum):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = AccTag_enum



class Poly:
    def __init__(self,*args):
        self.this = apply(MEADc.new_Poly,args)
        self.thisown = 1

    def derivative(*args):
        val = apply(MEADc.Poly_derivative,args)
        if val: val = PolyPtr(val) ; val.thisown = 1
        return val
    def degree(*args):
        val = apply(MEADc.Poly_degree,args)
        return val
    def size(*args):
        val = apply(MEADc.Poly_size,args)
        return val
    def coefficients(*args):
        val = apply(MEADc.Poly_coefficients,args)
        return val
    def set_varstring(*args):
        val = apply(MEADc.Poly_set_varstring,args)
        return val
    def varstring(*args):
        val = apply(MEADc.Poly_varstring,args)
        return val
    def output(*args):
        val = apply(MEADc.Poly_output,args)
        return val
    def write(*args):
        val = apply(MEADc.Poly_write,args)
        return val
    def __call__(*args):
        val = apply(MEADc.Poly___call__,args)
        return val
    def __neg__(*args):
        val = apply(MEADc.Poly___neg__,args)
        if val: val = PolyPtr(val) ; val.thisown = 1
        return val
    def __add__(*args):
        val = apply(MEADc.Poly___add__,args)
        if val: val = PolyPtr(val) ; val.thisown = 1
        return val
    def __sub__(*args):
        val = apply(MEADc.Poly___sub__,args)
        if val: val = PolyPtr(val) ; val.thisown = 1
        return val
    def __mul__(*args):
        val = apply(MEADc.Poly___mul__,args)
        if val: val = PolyPtr(val) ; val.thisown = 1
        return val
    def __mul____L82(*args):
        val = apply(MEADc.Poly___mul____L82,args)
        if val: val = PolyPtr(val) ; val.thisown = 1
        return val
    def __div__(*args):
        val = apply(MEADc.Poly___div__,args)
        if val: val = PolyPtr(val) ; val.thisown = 1
        return val
    def __eq__(*args):
        val = apply(MEADc.Poly___eq__,args)
        return val
    def __rmul__(*args):
        val = apply(MEADc.Poly___rmul__,args)
        if val: val = PolyPtr(val) ; val.thisown = 1
        return val
    def __pos__(*args):
        val = apply(MEADc.Poly___pos__,args)
        if val: val = PolyPtr(val) ; val.thisown = 1
        return val
    def __repr__(self):
        return "<C Poly instance at %s>" % (self.this,)
class PolyPtr(Poly):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = Poly








class Legendre(Poly):
    def __init__(self,*args):
        self.this = apply(MEADc.new_Legendre,args)
        self.thisown = 1

    def ell(*args):
        val = apply(MEADc.Legendre_ell,args)
        return val
    def __repr__(self):
        return "<C Legendre instance at %s>" % (self.this,)
class LegendrePtr(Legendre):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = Legendre






class CubeLatSpec:
    def __init__(self,*args):
        self.this = apply(MEADc.new_CubeLatSpec,args)
        self.thisown = 1

    def __del__(self,MEADc=MEADc):
        if self.thisown == 1 :
            MEADc.delete_CubeLatSpec(self)
    def resolve(*args):
        val = apply(MEADc.CubeLatSpec_resolve,args)
        return val
    def get_grid_dim(*args):
        val = apply(MEADc.CubeLatSpec_get_grid_dim,args)
        return val
    def get_spacing(*args):
        val = apply(MEADc.CubeLatSpec_get_spacing,args)
        return val
    def is_resolved(*args):
        val = apply(MEADc.CubeLatSpec_is_resolved,args)
        return val
    def get_centering_style(*args):
        val = apply(MEADc.CubeLatSpec_get_centering_style,args)
        return val
    def write(*args):
        val = apply(MEADc.CubeLatSpec_write,args)
        return val
    def get_center(*args):
        val = apply(MEADc.CubeLatSpec_get_center,args)
        if val: val = CoordPtr(val) ; val.thisown = 1
        return val
    def __eq__(*args):
        val = apply(MEADc.CubeLatSpec___eq__,args)
        return val
    def __lt__(*args):
        val = apply(MEADc.CubeLatSpec___lt__,args)
        return val
    def __gt__(*args):
        val = apply(MEADc.CubeLatSpec___gt__,args)
        return val
    def __cmp__(*args):
        val = apply(MEADc.CubeLatSpec___cmp__,args)
        return val
    def __repr__(self):
        return "<C CubeLatSpec instance at %s>" % (self.this,)
class CubeLatSpecPtr(CubeLatSpec):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = CubeLatSpec






class CenteringStyle_enum:
    ON_ORIGIN = MEADc.CenteringStyle_enum_ON_ORIGIN
    ON_CENT_OF_INTR = MEADc.CenteringStyle_enum_ON_CENT_OF_INTR
    ON_GEOM_CENT = MEADc.CenteringStyle_enum_ON_GEOM_CENT
    SPECIFIC = MEADc.CenteringStyle_enum_SPECIFIC
    def __init__(self,this):
        self.this = this

    def __repr__(self):
        return "<C CenteringStyle_enum instance at %s>" % (self.this,)
class CenteringStyle_enumPtr(CenteringStyle_enum):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = CenteringStyle_enum



class PhysCond:
    def __init__(self,*args):
        self.this = apply(MEADc.new_PhysCond,args)
        self.thisown = 1

    def write(*args):
        val = apply(MEADc.PhysCond_write,args)
        return val
    __setmethods__ = {
        "epsext" : MEADc.PhysCond_epsext_set,
        "solrad" : MEADc.PhysCond_solrad_set,
        "sterln" : MEADc.PhysCond_sterln_set,
        "ionicstr" : MEADc.PhysCond_ionicstr_set,
        "T" : MEADc.PhysCond_T_set,
        "kBolt" : MEADc.PhysCond_kBolt_set,
        "conconv" : MEADc.PhysCond_conconv_set,
        "econv" : MEADc.PhysCond_econv_set,
        "bohr_radius" : MEADc.PhysCond_bohr_radius_set,
        "proton_charge" : MEADc.PhysCond_proton_charge_set,
        "hueck" : MEADc.PhysCond_hueck_set,
        "kappasq" : MEADc.PhysCond_kappasq_set,
        "ln10kT" : MEADc.PhysCond_ln10kT_set,
    }
    def __setattr__(self,name,value):
        if (name == "this") or (name == "thisown"): self.__dict__[name] = value; return
        method = PhysCond.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value
    __getmethods__ = {
        "epsext" : MEADc.PhysCond_epsext_get,
        "solrad" : MEADc.PhysCond_solrad_get,
        "sterln" : MEADc.PhysCond_sterln_get,
        "ionicstr" : MEADc.PhysCond_ionicstr_get,
        "T" : MEADc.PhysCond_T_get,
        "kBolt" : MEADc.PhysCond_kBolt_get,
        "conconv" : MEADc.PhysCond_conconv_get,
        "econv" : MEADc.PhysCond_econv_get,
        "bohr_radius" : MEADc.PhysCond_bohr_radius_get,
        "proton_charge" : MEADc.PhysCond_proton_charge_get,
        "hueck" : MEADc.PhysCond_hueck_get,
        "kappasq" : MEADc.PhysCond_kappasq_get,
        "ln10kT" : MEADc.PhysCond_ln10kT_get,
    }
    def __getattr__(self,name):
        method = PhysCond.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C PhysCond instance at %s>" % (self.this,)
class PhysCondPtr(PhysCond):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = PhysCond



class Blab:
    def __init__(self,*args):
        self.this = apply(MEADc.new_Blab,args)
        self.thisown = 1

    def get_level(*args):
        val = apply(MEADc.Blab_get_level,args)
        return val
    def set_level(*args):
        val = apply(MEADc.Blab_set_level,args)
        return val
    __setmethods__ = {
        "level" : MEADc.Blab_level_set,
    }
    def __setattr__(self,name,value):
        if (name == "this") or (name == "thisown"): self.__dict__[name] = value; return
        method = Blab.__setmethods__.get(name,None)
        if method: return method(self,value)
        self.__dict__[name] = value
    __getmethods__ = {
        "level" : MEADc.Blab_level_get,
    }
    def __getattr__(self,name):
        method = Blab.__getmethods__.get(name,None)
        if method: return method(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C Blab instance at %s>" % (self.this,)
class BlabPtr(Blab):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = Blab



class FinDiffMethod:
    def __init__(self,*args):
        self.this = apply(MEADc.new_FinDiffMethod,args)
        self.thisown = 1

    def __del__(self,MEADc=MEADc):
        if self.thisown == 1 :
            MEADc.delete_FinDiffMethod(self)
    def read(*args):
        val = apply(MEADc.FinDiffMethod_read,args)
        return val
    def add_level(*args):
        val = apply(MEADc.FinDiffMethod_add_level,args)
        return val
    def resolve(*args):
        val = apply(MEADc.FinDiffMethod_resolve,args)
        return val
    def is_resolved(*args):
        val = apply(MEADc.FinDiffMethod_is_resolved,args)
        return val
    def add_level__L29(*args):
        val = apply(MEADc.FinDiffMethod_add_level__L29,args)
        return val
    def write(*args):
        val = apply(MEADc.FinDiffMethod_write,args)
        return val
    def levels(*args):
        val = apply(MEADc.FinDiffMethod_levels,args)
        return val
    def __repr__(self):
        return "<C FinDiffMethod instance at %s>" % (self.this,)
class FinDiffMethodPtr(FinDiffMethod):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = FinDiffMethod





class AnalySphere(AnalyticEP):
    def __init__(self,*args):
        self.this = apply(MEADc.new_AnalySphere,args)
        self.thisown = 1

    def __del__(self,MEADc=MEADc):
        if self.thisown == 1 :
            MEADc.delete_AnalySphere(self)
    def solve(*args):
        val = apply(MEADc.AnalySphere_solve,args)
        return val
    def value(*args):
        val = apply(MEADc.AnalySphere_value,args)
        return val
    def field(*args):
        val = apply(MEADc.AnalySphere_field,args)
        if val: val = CoordPtr(val) ; val.thisown = 1
        return val
    def displacement(*args):
        val = apply(MEADc.AnalySphere_displacement,args)
        if val: val = CoordPtr(val) ; val.thisown = 1
        return val
    def get_rad_diel(*args):
        val = apply(MEADc.AnalySphere_get_rad_diel,args)
        return val
    def get_rad_ely(*args):
        val = apply(MEADc.AnalySphere_get_rad_ely,args)
        return val
    def get_maxterm(*args):
        val = apply(MEADc.AnalySphere_get_maxterm,args)
        return val
    def __mul__(*args):
        val = apply(MEADc.AnalySphere___mul__,args)
        return val
    def __mul____L69(*args):
        val = apply(MEADc.AnalySphere___mul____L69,args)
        if val: val = ElstatPotCombinationPtr(val) ; val.thisown = 1
        return val
    def __rmul__(*args):
        val = apply(MEADc.AnalySphere___rmul__,args)
        if val: val = ElstatPotCombinationPtr(val) ; val.thisown = 1
        return val
    def __div__(*args):
        val = apply(MEADc.AnalySphere___div__,args)
        if val: val = ElstatPotCombinationPtr(val) ; val.thisown = 1
        return val
    def __repr__(self):
        return "<C AnalySphere instance at %s>" % (self.this,)
class AnalySpherePtr(AnalySphere):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = AnalySphere




class AnalySlab(AnalyticEP):
    def __init__(self,*args):
        self.this = apply(MEADc.new_AnalySlab,args)
        self.thisown = 1

    def __del__(self,MEADc=MEADc):
        if self.thisown == 1 :
            MEADc.delete_AnalySlab(self)
    def solve(*args):
        val = apply(MEADc.AnalySlab_solve,args)
        return val
    def value(*args):
        val = apply(MEADc.AnalySlab_value,args)
        return val
    def field(*args):
        val = apply(MEADc.AnalySlab_field,args)
        if val: val = CoordPtr(val) ; val.thisown = 1
        return val
    def displacement(*args):
        val = apply(MEADc.AnalySlab_displacement,args)
        if val: val = CoordPtr(val) ; val.thisown = 1
        return val
    def get_zlower(*args):
        val = apply(MEADc.AnalySlab_get_zlower,args)
        return val
    def get_zupper(*args):
        val = apply(MEADc.AnalySlab_get_zupper,args)
        return val
    def get_maxterm(*args):
        val = apply(MEADc.AnalySlab_get_maxterm,args)
        return val
    def __mul__(*args):
        val = apply(MEADc.AnalySlab___mul__,args)
        return val
    def __mul____L71(*args):
        val = apply(MEADc.AnalySlab___mul____L71,args)
        if val: val = ElstatPotCombinationPtr(val) ; val.thisown = 1
        return val
    def __rmul__(*args):
        val = apply(MEADc.AnalySlab___rmul__,args)
        if val: val = ElstatPotCombinationPtr(val) ; val.thisown = 1
        return val
    def __div__(*args):
        val = apply(MEADc.AnalySlab___div__,args)
        if val: val = ElstatPotCombinationPtr(val) ; val.thisown = 1
        return val
    def __repr__(self):
        return "<C AnalySlab instance at %s>" % (self.this,)
class AnalySlabPtr(AnalySlab):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = AnalySlab




class Debye(AnalyticEP):
    def __init__(self,*args):
        self.this = apply(MEADc.new_Debye,args)
        self.thisown = 1

    def __del__(self,MEADc=MEADc):
        if self.thisown == 1 :
            MEADc.delete_Debye(self)
    def solve(*args):
        val = apply(MEADc.Debye_solve,args)
        return val
    def value(*args):
        val = apply(MEADc.Debye_value,args)
        return val
    def field(*args):
        val = apply(MEADc.Debye_field,args)
        if val: val = CoordPtr(val) ; val.thisown = 1
        return val
    def displacement(*args):
        val = apply(MEADc.Debye_displacement,args)
        if val: val = CoordPtr(val) ; val.thisown = 1
        return val
    def get_kappa(*args):
        val = apply(MEADc.Debye_get_kappa,args)
        return val
    def __mul__(*args):
        val = apply(MEADc.Debye___mul__,args)
        return val
    def __mul____L47(*args):
        val = apply(MEADc.Debye___mul____L47,args)
        if val: val = ElstatPotCombinationPtr(val) ; val.thisown = 1
        return val
    def __rmul__(*args):
        val = apply(MEADc.Debye___rmul__,args)
        if val: val = ElstatPotCombinationPtr(val) ; val.thisown = 1
        return val
    def __div__(*args):
        val = apply(MEADc.Debye___div__,args)
        if val: val = ElstatPotCombinationPtr(val) ; val.thisown = 1
        return val
    def __repr__(self):
        return "<C Debye instance at %s>" % (self.this,)
class DebyePtr(Debye):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = Debye



class FinDiffElstatPot(ElstatPot_lett):
    def __init__(self,*args):
        self.this = apply(MEADc.new_FinDiffElstatPot,args)
        self.thisown = 1

    def __del__(self,MEADc=MEADc):
        if self.thisown == 1 :
            MEADc.delete_FinDiffElstatPot(self)
    def solve(*args):
        val = apply(MEADc.FinDiffElstatPot_solve,args)
        return val
    def solve_using_coarse_init(*args):
        val = apply(MEADc.FinDiffElstatPot_solve_using_coarse_init,args)
        return val
    def write_coarse_field(*args):
        val = apply(MEADc.FinDiffElstatPot_write_coarse_field,args)
        return val
    def coarse_lattice_spec(*args):
        val = apply(MEADc.FinDiffElstatPot_coarse_lattice_spec,args)
        if val: val = CubeLatSpecPtr(val) ; val.thisown = 1
        return val
    def value(*args):
        val = apply(MEADc.FinDiffElstatPot_value,args)
        return val
    def field(*args):
        val = apply(MEADc.FinDiffElstatPot_field,args)
        if val: val = CoordPtr(val) ; val.thisown = 1
        return val
    def displacement(*args):
        val = apply(MEADc.FinDiffElstatPot_displacement,args)
        if val: val = CoordPtr(val) ; val.thisown = 1
        return val
    def get_cuberep(*args):
        val = apply(MEADc.FinDiffElstatPot_get_cuberep,args)
        return val
    def __mul__(*args):
        val = apply(MEADc.FinDiffElstatPot___mul__,args)
        return val
    def __mul____L97(*args):
        val = apply(MEADc.FinDiffElstatPot___mul____L97,args)
        if val: val = ElstatPotCombinationPtr(val) ; val.thisown = 1
        return val
    def __rmul__(*args):
        val = apply(MEADc.FinDiffElstatPot___rmul__,args)
        if val: val = ElstatPotCombinationPtr(val) ; val.thisown = 1
        return val
    def __div__(*args):
        val = apply(MEADc.FinDiffElstatPot___div__,args)
        if val: val = ElstatPotCombinationPtr(val) ; val.thisown = 1
        return val
    def __repr__(self):
        return "<C FinDiffElstatPot instance at %s>" % (self.this,)
class FinDiffElstatPotPtr(FinDiffElstatPot):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = FinDiffElstatPot




class ElstatPotCombination:
    def __init__(self,*args):
        self.this = apply(MEADc.new_ElstatPotCombination,args)
        self.thisown = 1

    def solve(*args):
        val = apply(MEADc.ElstatPotCombination_solve,args)
        return val
    def value(*args):
        val = apply(MEADc.ElstatPotCombination_value,args)
        return val
    def field(*args):
        val = apply(MEADc.ElstatPotCombination_field,args)
        if val: val = CoordPtr(val) ; val.thisown = 1
        return val
    def displacement(*args):
        val = apply(MEADc.ElstatPotCombination_displacement,args)
        if val: val = CoordPtr(val) ; val.thisown = 1
        return val
    def __del__(self,MEADc=MEADc):
        if self.thisown == 1 :
            MEADc.delete_ElstatPotCombination(self)
    def __iadd__(*args):
        val = apply(MEADc.ElstatPotCombination___iadd__,args)
        if val: val = ElstatPotCombinationPtr(val) ; val.thisown = 1
        return val
    def __iadd____L90(*args):
        val = apply(MEADc.ElstatPotCombination___iadd____L90,args)
        if val: val = ElstatPotCombinationPtr(val) ; val.thisown = 1
        return val
    def __imul__(*args):
        val = apply(MEADc.ElstatPotCombination___imul__,args)
        if val: val = ElstatPotCombinationPtr(val) ; val.thisown = 1
        return val
    def __idiv__(*args):
        val = apply(MEADc.ElstatPotCombination___idiv__,args)
        if val: val = ElstatPotCombinationPtr(val) ; val.thisown = 1
        return val
    def __add__(*args):
        val = apply(MEADc.ElstatPotCombination___add__,args)
        if val: val = ElstatPotCombinationPtr(val) ; val.thisown = 1
        return val
    def __add____L98(*args):
        val = apply(MEADc.ElstatPotCombination___add____L98,args)
        if val: val = ElstatPotCombinationPtr(val) ; val.thisown = 1
        return val
    def __mul__(*args):
        val = apply(MEADc.ElstatPotCombination___mul__,args)
        if val: val = ElstatPotCombinationPtr(val) ; val.thisown = 1
        return val
    def __rmul__(*args):
        val = apply(MEADc.ElstatPotCombination___rmul__,args)
        if val: val = ElstatPotCombinationPtr(val) ; val.thisown = 1
        return val
    def __div__(*args):
        val = apply(MEADc.ElstatPotCombination___div__,args)
        if val: val = ElstatPotCombinationPtr(val) ; val.thisown = 1
        return val
    def __mul____L106(*args):
        val = apply(MEADc.ElstatPotCombination___mul____L106,args)
        return val
    def __repr__(self):
        return "<C ElstatPotCombination instance at %s>" % (self.this,)
class ElstatPotCombinationPtr(ElstatPotCombination):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = ElstatPotCombination





class Moments:
    def __init__(self,*args):
        self.this = apply(MEADc.new_Moments,args)
        self.thisown = 1

    def ellmax(*args):
        val = apply(MEADc.Moments_ellmax,args)
        return val
    def get_imag(*args):
        val = apply(MEADc.Moments_get_imag,args)
        return val
    def set(*args):
        val = apply(MEADc.Moments_set,args)
        return val
    def __neg__(*args):
        val = apply(MEADc.Moments___neg__,args)
        if val: val = MomentsPtr(val) ; val.thisown = 1
        return val
    def __add__(*args):
        val = apply(MEADc.Moments___add__,args)
        if val: val = MomentsPtr(val) ; val.thisown = 1
        return val
    def __sub__(*args):
        val = apply(MEADc.Moments___sub__,args)
        if val: val = MomentsPtr(val) ; val.thisown = 1
        return val
    def __mul__(*args):
        val = apply(MEADc.Moments___mul__,args)
        if val: val = MomentsPtr(val) ; val.thisown = 1
        return val
    def __rmul__(*args):
        val = apply(MEADc.Moments___rmul__,args)
        if val: val = MomentsPtr(val) ; val.thisown = 1
        return val
    def __pos__(*args):
        val = apply(MEADc.Moments___pos__,args)
        if val: val = MomentsPtr(val) ; val.thisown = 1
        return val
    def __repr__(self):
        return "<C Moments instance at %s>" % (self.this,)
class MomentsPtr(Moments):
    def __init__(self,this):
        self.this = this
        self.thisown = 0
        self.__class__ = Moments





#-------------- FUNCTION WRAPPERS ------------------

def order_one_Legendre(*args, **kwargs):
    val = apply(MEADc.order_one_Legendre,args,kwargs)
    if val: val = LegendrePtr(val); val.thisown = 1
    return val

def next_Legendre_using_prev2(*args, **kwargs):
    val = apply(MEADc.next_Legendre_using_prev2,args,kwargs)
    if val: val = LegendrePtr(val); val.thisown = 1
    return val

create_legendre_series = MEADc.create_legendre_series

def momentsOfElstatPot(*args, **kwargs):
    val = apply(MEADc.momentsOfElstatPot,args,kwargs)
    if val: val = MomentsPtr(val); val.thisown = 1
    return val

print_moments = MEADc.print_moments

compare_moments = MEADc.compare_moments

def momentsOfChargeDist(*args, **kwargs):
    val = apply(MEADc.momentsOfChargeDist,args,kwargs)
    if val: val = MomentsPtr(val); val.thisown = 1
    return val



#-------------- VARIABLE WRAPPERS ------------------



#-------------- USER INCLUDE -----------------------

import exceptions

class Error(exceptions.Exception):
   def __init__(self, args=None):
      self.args = args
import types, re
class __Dummy:
   def __init__(self,*args):
      trydefault = 0
      try:
         if len(args)==1:
            self.this = apply(MEADc.new_Coord__L17,args)
         elif len(args)==3:
            self.this = apply(MEADc.new_Coord__L16,args)
         else:
            trydefault = 1
            self.this = apply(MEADc.new_Coord,args)
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               self.this = apply(MEADc.new_Coord,args)
            except StandardError, e:
               error = 'No Coord methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               self.thisown = 0
               raise Error, error
         else:
            self.thisown = 0
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No Coord methods take ' + nargs
               raise Error, error
            else:
               raise e
      self.thisown = 1
      return
   def __mul__(*args):
      trydefault = 0
      try:
         if len(args)==2:
            val = apply(MEADc.Coord___mul____L86,args)
            if val and isinstance(val, types.StringType) and re.search('_p_float', val):
               val = floatPtr(val) ; val.thisown = 1
         else:
            trydefault = 1
            val = apply(MEADc.Coord___mul__,args)
            if val and isinstance(val, types.StringType):
               mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
               if mo:
                  rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                  val = eval(rettype) ; val.thisown = 1
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               val = apply(MEADc.Coord___mul__,args)
               if val and isinstance(val, types.StringType):
                  mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
                  if mo:
                     rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                     val = eval(rettype) ; val.thisown = 1
            except StandardError, e:
               error = 'No Coord.__mul__ methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               raise Error, error
         else:
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No Coord.__mul__ methods take ' + nargs
               raise Error, error
            else:
               raise e
      return val
   def __deepcopy__(self, memo = None):
      val = MEADc.new_Coord__deepcopy(self)
      if val: val = CoordPtr(val) ; val.thisown = 1
      return val

Coord.__init__ = __Dummy.__dict__['__init__']
Coord.__mul__ = __Dummy.__dict__['__mul__']
del Coord.__dict__['_Coord__mul____L86']
Coord.__deepcopy__ = __Dummy.__dict__['__deepcopy__']

import types, re
import types, re
import types, re
import types, re
class __Dummy:
   def __add__(*args):
      trydefault = 0
      try:
         if len(args)==2:
            val = apply(MEADc.ElstatPot_lett___add____L80,args)
            if val and isinstance(val, types.StringType) and re.search('_p_ElstatPotCombination', val):
               val = ElstatPotCombinationPtr(val) ; val.thisown = 1
         else:
            trydefault = 1
            val = apply(MEADc.ElstatPot_lett___add__,args)
            if val and isinstance(val, types.StringType):
               mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
               if mo:
                  rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                  val = eval(rettype) ; val.thisown = 1
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               val = apply(MEADc.ElstatPot_lett___add__,args)
               if val and isinstance(val, types.StringType):
                  mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
                  if mo:
                     rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                     val = eval(rettype) ; val.thisown = 1
            except StandardError, e:
               error = 'No ElstatPot_lett.__add__ methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               raise Error, error
         else:
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No ElstatPot_lett.__add__ methods take ' + nargs
               raise Error, error
            else:
               raise e
      return val

ElstatPot_lett.__add__ = __Dummy.__dict__['__add__']
del ElstatPot_lett.__dict__['_ElstatPot_lett__add____L80']

import types, re
class __Dummy:
   def __init__(self,*args):
      trydefault = 0
      try:
         if len(args)==1:
            self.this = apply(MEADc.new_Atom__L23,args)
         else:
            trydefault = 1
            self.this = apply(MEADc.new_Atom,args)
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               self.this = apply(MEADc.new_Atom,args)
            except StandardError, e:
               error = 'No Atom methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               self.thisown = 0
               raise Error, error
         else:
            self.thisown = 0
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No Atom methods take ' + nargs
               raise Error, error
            else:
               raise e
      self.thisown = 1
      return
   def __deepcopy__(self, memo = None):
      val = MEADc.new_Atom__deepcopy(self)
      if val: val = AtomPtr(val) ; val.thisown = 1
      return val

Atom.__init__ = __Dummy.__dict__['__init__']
Atom.__deepcopy__ = __Dummy.__dict__['__deepcopy__']

import types, re
class __Dummy:
   def __init__(self,*args):
      trydefault = 0
      try:
         if len(args)==1:
            try:
               self.this = apply(MEADc.new_AtomID__L17,args)
            except:
               self.this = apply(MEADc.new_AtomID__L18,args)
         elif len(args)==2:
            self.this = apply(MEADc.new_AtomID__L41,args)
         elif len(args)==3:
            self.this = apply(MEADc.new_AtomID__L16,args)
         else:
            trydefault = 1
            self.this = apply(MEADc.new_AtomID,args)
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               self.this = apply(MEADc.new_AtomID,args)
            except StandardError, e:
               error = 'No AtomID methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               self.thisown = 0
               raise Error, error
         else:
            self.thisown = 0
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No AtomID methods take ' + nargs
               raise Error, error
            else:
               raise e
      self.thisown = 1
      return
   def __deepcopy__(self, memo = None):
      val = MEADc.new_AtomID__deepcopy(self)
      if val: val = AtomIDPtr(val) ; val.thisown = 1
      return val

AtomID.__init__ = __Dummy.__dict__['__init__']
AtomID.__deepcopy__ = __Dummy.__dict__['__deepcopy__']

import types, re
class __Dummy:
   def read(*args):
      trydefault = 0
      try:
         if len(args)==2:
            val = apply(MEADc.AtomSet_read__L21,args)
            if val and isinstance(val, types.StringType) and re.search('_p_void', val):
               val = voidPtr(val) ; val.thisown = 1
         else:
            trydefault = 1
            val = apply(MEADc.AtomSet_read,args)
            if val and isinstance(val, types.StringType):
               mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
               if mo:
                  rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                  val = eval(rettype) ; val.thisown = 1
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               val = apply(MEADc.AtomSet_read,args)
               if val and isinstance(val, types.StringType):
                  mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
                  if mo:
                     rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                     val = eval(rettype) ; val.thisown = 1
            except StandardError, e:
               error = 'No AtomSet.read methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               raise Error, error
         else:
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No AtomSet.read methods take ' + nargs
               raise Error, error
            else:
               raise e
      return val
   def __init__(self,*args):
      trydefault = 0
      try:
         if len(args)==1:
            try:
               self.this = apply(MEADc.new_AtomSet__L16,args)
            except:
               try:
                  self.this = apply(MEADc.new_AtomSet__L17,args)
               except:
                  self.this = apply(MEADc.new_AtomSet__L18,args)
         else:
            trydefault = 1
            self.this = apply(MEADc.new_AtomSet,args)
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               self.this = apply(MEADc.new_AtomSet,args)
            except StandardError, e:
               error = 'No AtomSet methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               self.thisown = 0
               raise Error, error
         else:
            self.thisown = 0
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No AtomSet methods take ' + nargs
               raise Error, error
            else:
               raise e
      self.thisown = 1
      return
   def __deepcopy__(self, memo = None):
      val = MEADc.new_AtomSet__deepcopy(self)
      if val: val = AtomSetPtr(val) ; val.thisown = 1
      return val

AtomSet.read = __Dummy.__dict__['read']
del AtomSet.__dict__['read__L21']
AtomSet.__init__ = __Dummy.__dict__['__init__']
AtomSet.__deepcopy__ = __Dummy.__dict__['__deepcopy__']

import types, re
class __Dummy:
   def __init__(self,*args):
      trydefault = 0
      try:
         if len(args)==2:
            try:
               self.this = apply(MEADc.new_PointCharge__L16,args)
            except:
               self.this = apply(MEADc.new_PointCharge__L17,args)
         else:
            trydefault = 1
            self.this = apply(MEADc.new_PointCharge,args)
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               self.this = apply(MEADc.new_PointCharge,args)
            except StandardError, e:
               error = 'No PointCharge methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               self.thisown = 0
               raise Error, error
         else:
            self.thisown = 0
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No PointCharge methods take ' + nargs
               raise Error, error
            else:
               raise e
      self.thisown = 1
      return

PointCharge.__init__ = __Dummy.__dict__['__init__']

import types, re
class __Dummy:
   def __init__(self,*args):
      trydefault = 0
      try:
         if len(args)==1:
            try:
               self.this = apply(MEADc.new_AtomChargeSet__L16,args)
            except:
               self.this = apply(MEADc.new_AtomChargeSet__L17,args)
         else:
            trydefault = 1
            self.this = apply(MEADc.new_AtomChargeSet,args)
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               self.this = apply(MEADc.new_AtomChargeSet,args)
            except StandardError, e:
               error = 'No AtomChargeSet methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               self.thisown = 0
               raise Error, error
         else:
            self.thisown = 0
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No AtomChargeSet methods take ' + nargs
               raise Error, error
            else:
               raise e
      self.thisown = 1
      return
   def __mul__(*args):
      trydefault = 0
      try:
         if len(args)==2:
            try:
               val = apply(MEADc.AtomChargeSet___mul____L113,args)
               if val and isinstance(val, types.StringType) and re.search('_p_float', val):
                  val = floatPtr(val) ; val.thisown = 1
            except:
               val = apply(MEADc.AtomChargeSet___mul____L114,args)
               if val and isinstance(val, types.StringType) and re.search('_p_float', val):
                  val = floatPtr(val) ; val.thisown = 1
         else:
            trydefault = 1
            val = apply(MEADc.AtomChargeSet___mul__,args)
            if val and isinstance(val, types.StringType):
               mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
               if mo:
                  rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                  val = eval(rettype) ; val.thisown = 1
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               val = apply(MEADc.AtomChargeSet___mul__,args)
               if val and isinstance(val, types.StringType):
                  mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
                  if mo:
                     rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                     val = eval(rettype) ; val.thisown = 1
            except StandardError, e:
               error = 'No AtomChargeSet.__mul__ methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               raise Error, error
         else:
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No AtomChargeSet.__mul__ methods take ' + nargs
               raise Error, error
            else:
               raise e
      return val
   def __add__(*args):
      trydefault = 0
      try:
         if len(args)==2:
            try:
               val = apply(MEADc.AtomChargeSet___add____L68,args)
               if val and isinstance(val, types.StringType) and re.search('_p_ManyPointCharge', val):
                  val = ManyPointChargePtr(val) ; val.thisown = 1
            except:
               val = apply(MEADc.AtomChargeSet___add____L69,args)
               if val and isinstance(val, types.StringType) and re.search('_p_ManyPointCharge', val):
                  val = ManyPointChargePtr(val) ; val.thisown = 1
         else:
            trydefault = 1
            val = apply(MEADc.AtomChargeSet___add__,args)
            if val and isinstance(val, types.StringType):
               mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
               if mo:
                  rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                  val = eval(rettype) ; val.thisown = 1
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               val = apply(MEADc.AtomChargeSet___add__,args)
               if val and isinstance(val, types.StringType):
                  mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
                  if mo:
                     rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                     val = eval(rettype) ; val.thisown = 1
            except StandardError, e:
               error = 'No AtomChargeSet.__add__ methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               raise Error, error
         else:
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No AtomChargeSet.__add__ methods take ' + nargs
               raise Error, error
            else:
               raise e
      return val
   def __deepcopy__(self, memo = None):
      val = MEADc.new_AtomChargeSet__deepcopy(self)
      if val: val = AtomChargeSetPtr(val) ; val.thisown = 1
      return val

AtomChargeSet.__init__ = __Dummy.__dict__['__init__']
AtomChargeSet.__mul__ = __Dummy.__dict__['__mul__']
del AtomChargeSet.__dict__['_AtomChargeSet__mul____L113']
del AtomChargeSet.__dict__['_AtomChargeSet__mul____L114']
AtomChargeSet.__add__ = __Dummy.__dict__['__add__']
del AtomChargeSet.__dict__['_AtomChargeSet__add____L68']
del AtomChargeSet.__dict__['_AtomChargeSet__add____L69']
AtomChargeSet.__deepcopy__ = __Dummy.__dict__['__deepcopy__']

import types, re
class __Dummy:
   def __mul__(*args):
      trydefault = 0
      try:
         if len(args)==2:
            try:
               val = apply(MEADc.ManyPointCharge___mul____L119,args)
               if val and isinstance(val, types.StringType) and re.search('_p_float', val):
                  val = floatPtr(val) ; val.thisown = 1
            except:
               val = apply(MEADc.ManyPointCharge___mul____L120,args)
               if val and isinstance(val, types.StringType) and re.search('_p_float', val):
                  val = floatPtr(val) ; val.thisown = 1
         else:
            trydefault = 1
            val = apply(MEADc.ManyPointCharge___mul__,args)
            if val and isinstance(val, types.StringType):
               mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
               if mo:
                  rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                  val = eval(rettype) ; val.thisown = 1
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               val = apply(MEADc.ManyPointCharge___mul__,args)
               if val and isinstance(val, types.StringType):
                  mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
                  if mo:
                     rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                     val = eval(rettype) ; val.thisown = 1
            except StandardError, e:
               error = 'No ManyPointCharge.__mul__ methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               raise Error, error
         else:
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No ManyPointCharge.__mul__ methods take ' + nargs
               raise Error, error
            else:
               raise e
      return val
   def __init__(self,*args):
      trydefault = 0
      try:
         if len(args)==1:
            self.this = apply(MEADc.new_ManyPointCharge__L16,args)
         elif len(args)==3:
            self.this = apply(MEADc.new_ManyPointCharge__L17,args)
         else:
            trydefault = 1
            self.this = apply(MEADc.new_ManyPointCharge,args)
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               self.this = apply(MEADc.new_ManyPointCharge,args)
            except StandardError, e:
               error = 'No ManyPointCharge methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               self.thisown = 0
               raise Error, error
         else:
            self.thisown = 0
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No ManyPointCharge methods take ' + nargs
               raise Error, error
            else:
               raise e
      self.thisown = 1
      return
   def __add__(*args):
      trydefault = 0
      try:
         if len(args)==2:
            try:
               val = apply(MEADc.ManyPointCharge___add____L55,args)
               if val and isinstance(val, types.StringType) and re.search('_p_ManyPointCharge', val):
                  val = ManyPointChargePtr(val) ; val.thisown = 1
            except:
               val = apply(MEADc.ManyPointCharge___add____L56,args)
               if val and isinstance(val, types.StringType) and re.search('_p_ManyPointCharge', val):
                  val = ManyPointChargePtr(val) ; val.thisown = 1
         else:
            trydefault = 1
            val = apply(MEADc.ManyPointCharge___add__,args)
            if val and isinstance(val, types.StringType):
               mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
               if mo:
                  rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                  val = eval(rettype) ; val.thisown = 1
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               val = apply(MEADc.ManyPointCharge___add__,args)
               if val and isinstance(val, types.StringType):
                  mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
                  if mo:
                     rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                     val = eval(rettype) ; val.thisown = 1
            except StandardError, e:
               error = 'No ManyPointCharge.__add__ methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               raise Error, error
         else:
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No ManyPointCharge.__add__ methods take ' + nargs
               raise Error, error
            else:
               raise e
      return val

ManyPointCharge.__mul__ = __Dummy.__dict__['__mul__']
del ManyPointCharge.__dict__['_ManyPointCharge__mul____L119']
del ManyPointCharge.__dict__['_ManyPointCharge__mul____L120']
ManyPointCharge.__init__ = __Dummy.__dict__['__init__']
ManyPointCharge.__add__ = __Dummy.__dict__['__add__']
del ManyPointCharge.__dict__['_ManyPointCharge__add____L55']
del ManyPointCharge.__dict__['_ManyPointCharge__add____L56']

import types, re
class __Dummy:
   def __mul__(*args):
      trydefault = 0
      try:
         if len(args)==2:
            try:
               val = apply(MEADc.OnePointCharge___mul____L123,args)
               if val and isinstance(val, types.StringType) and re.search('_p_float', val):
                  val = floatPtr(val) ; val.thisown = 1
            except:
               val = apply(MEADc.OnePointCharge___mul____L124,args)
               if val and isinstance(val, types.StringType) and re.search('_p_float', val):
                  val = floatPtr(val) ; val.thisown = 1
         else:
            trydefault = 1
            val = apply(MEADc.OnePointCharge___mul__,args)
            if val and isinstance(val, types.StringType):
               mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
               if mo:
                  rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                  val = eval(rettype) ; val.thisown = 1
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               val = apply(MEADc.OnePointCharge___mul__,args)
               if val and isinstance(val, types.StringType):
                  mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
                  if mo:
                     rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                     val = eval(rettype) ; val.thisown = 1
            except StandardError, e:
               error = 'No OnePointCharge.__mul__ methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               raise Error, error
         else:
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No OnePointCharge.__mul__ methods take ' + nargs
               raise Error, error
            else:
               raise e
      return val
   def __add__(*args):
      trydefault = 0
      try:
         if len(args)==2:
            try:
               val = apply(MEADc.OnePointCharge___add____L66,args)
               if val and isinstance(val, types.StringType) and re.search('_p_ManyPointCharge', val):
                  val = ManyPointChargePtr(val) ; val.thisown = 1
            except:
               val = apply(MEADc.OnePointCharge___add____L69,args)
               if val and isinstance(val, types.StringType) and re.search('_p_ManyPointCharge', val):
                  val = ManyPointChargePtr(val) ; val.thisown = 1
         else:
            trydefault = 1
            val = apply(MEADc.OnePointCharge___add__,args)
            if val and isinstance(val, types.StringType):
               mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
               if mo:
                  rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                  val = eval(rettype) ; val.thisown = 1
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               val = apply(MEADc.OnePointCharge___add__,args)
               if val and isinstance(val, types.StringType):
                  mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
                  if mo:
                     rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                     val = eval(rettype) ; val.thisown = 1
            except StandardError, e:
               error = 'No OnePointCharge.__add__ methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               raise Error, error
         else:
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No OnePointCharge.__add__ methods take ' + nargs
               raise Error, error
            else:
               raise e
      return val
   def __init__(self,*args):
      trydefault = 0
      try:
         if len(args)==2:
            self.this = apply(MEADc.new_OnePointCharge__L16,args)
         else:
            trydefault = 1
            self.this = apply(MEADc.new_OnePointCharge,args)
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               self.this = apply(MEADc.new_OnePointCharge,args)
            except StandardError, e:
               error = 'No OnePointCharge methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               self.thisown = 0
               raise Error, error
         else:
            self.thisown = 0
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No OnePointCharge methods take ' + nargs
               raise Error, error
            else:
               raise e
      self.thisown = 1
      return

OnePointCharge.__mul__ = __Dummy.__dict__['__mul__']
del OnePointCharge.__dict__['_OnePointCharge__mul____L123']
del OnePointCharge.__dict__['_OnePointCharge__mul____L124']
OnePointCharge.__add__ = __Dummy.__dict__['__add__']
del OnePointCharge.__dict__['_OnePointCharge__add____L66']
del OnePointCharge.__dict__['_OnePointCharge__add____L69']
OnePointCharge.__init__ = __Dummy.__dict__['__init__']

import types, re
import types, re
import types, re
import types, re
class __Dummy:
   def __init__(self,*args):
      trydefault = 0
      try:
         if len(args)==4:
            self.this = apply(MEADc.new_ThreeValueDielectricByAtoms__L30,args)
         elif len(args)==5:
            self.this = apply(MEADc.new_ThreeValueDielectricByAtoms__L31,args)
         elif len(args)==7:
            self.this = apply(MEADc.new_ThreeValueDielectricByAtoms__L29,args)
         else:
            trydefault = 1
            self.this = apply(MEADc.new_ThreeValueDielectricByAtoms,args)
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               self.this = apply(MEADc.new_ThreeValueDielectricByAtoms,args)
            except StandardError, e:
               error = 'No ThreeValueDielectricByAtoms methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               self.thisown = 0
               raise Error, error
         else:
            self.thisown = 0
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No ThreeValueDielectricByAtoms methods take ' + nargs
               raise Error, error
            else:
               raise e
      self.thisown = 1
      return

ThreeValueDielectricByAtoms.__init__ = __Dummy.__dict__['__init__']

class __Dummy:
   def __init__(self,*args):
      trydefault = 0
      try:
         if len(args)==2:
            self.this = apply(MEADc.new_TwoValueDielectricByAtoms__L18,args)
         elif len(args)==3:
            self.this = apply(MEADc.new_TwoValueDielectricByAtoms__L19,args)
         elif len(args)==4:
            self.this = apply(MEADc.new_TwoValueDielectricByAtoms__L17,args)
         else:
            trydefault = 1
            self.this = apply(MEADc.new_TwoValueDielectricByAtoms,args)
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               self.this = apply(MEADc.new_TwoValueDielectricByAtoms,args)
            except StandardError, e:
               error = 'No TwoValueDielectricByAtoms methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               self.thisown = 0
               raise Error, error
         else:
            self.thisown = 0
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No TwoValueDielectricByAtoms methods take ' + nargs
               raise Error, error
            else:
               raise e
      self.thisown = 1
      return

TwoValueDielectricByAtoms.__init__ = __Dummy.__dict__['__init__']

import types, re
class __Dummy:
   def __init__(self,*args):
      trydefault = 0
      try:
         if len(args)==6:
            self.this = apply(MEADc.new_TwoValueDielMembAtoms__L16,args)
         else:
            trydefault = 1
            self.this = apply(MEADc.new_TwoValueDielMembAtoms,args)
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               self.this = apply(MEADc.new_TwoValueDielMembAtoms,args)
            except StandardError, e:
               error = 'No TwoValueDielMembAtoms methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               self.thisown = 0
               raise Error, error
         else:
            self.thisown = 0
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No TwoValueDielMembAtoms methods take ' + nargs
               raise Error, error
            else:
               raise e
      self.thisown = 1
      return

TwoValueDielMembAtoms.__init__ = __Dummy.__dict__['__init__']

import types, re
import types, re
import types, re
class __Dummy:
   def __init__(self,*args):
      trydefault = 0
      try:
         if len(args)==1:
            self.this = apply(MEADc.new_ElectrolyteByAtoms__L31,args)
         elif len(args)==2:
            self.this = apply(MEADc.new_ElectrolyteByAtoms__L32,args)
         elif len(args)==3:
            self.this = apply(MEADc.new_ElectrolyteByAtoms__L18,args)
         else:
            trydefault = 1
            self.this = apply(MEADc.new_ElectrolyteByAtoms,args)
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               self.this = apply(MEADc.new_ElectrolyteByAtoms,args)
            except StandardError, e:
               error = 'No ElectrolyteByAtoms methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               self.thisown = 0
               raise Error, error
         else:
            self.thisown = 0
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No ElectrolyteByAtoms methods take ' + nargs
               raise Error, error
            else:
               raise e
      self.thisown = 1
      return

ElectrolyteByAtoms.__init__ = __Dummy.__dict__['__init__']

# Class to implement the AccTag enum
# Define the enum values in the dict from a SWIG generated AccTag_enum class
# Don't allow these enum attributes to be set
# (But allow other attributes to be added)

class AccTag_enum(AccTag_enum):
   def __init__(self):
      self.this = "None"
      self.__dict__.update({'interior' : AccTag_enum.interior, 'exterior' : AccTag_enum.exterior, 'undecided' : AccTag_enum.undecided, 'in_tube' : AccTag_enum.in_tube})
   def __setattr__(self, name, value):
      if (name == "interior") or (name == "exterior") or (name == "undecided") or (name == "in_tube"):
         error = "Can't set constant enum attribute " + name
         raise AttributeError, error
      else:
         self.__dict__[name] = value

# Create an instance
AccTag = AccTag_enum()
import types, re
class __Dummy:
   def read_top_in_binary(*args):
      trydefault = 0
      try:
         if len(args)==2:
            val = apply(MEADc.SolvAccVol_read_top_in_binary__L32,args)
            if val and isinstance(val, types.StringType) and re.search('_p_int', val):
               val = intPtr(val) ; val.thisown = 1
         else:
            trydefault = 1
            val = apply(MEADc.SolvAccVol_read_top_in_binary,args)
            if val and isinstance(val, types.StringType):
               mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
               if mo:
                  rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                  val = eval(rettype) ; val.thisown = 1
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               val = apply(MEADc.SolvAccVol_read_top_in_binary,args)
               if val and isinstance(val, types.StringType):
                  mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
                  if mo:
                     rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                     val = eval(rettype) ; val.thisown = 1
            except StandardError, e:
               error = 'No SolvAccVol.read_top_in_binary methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               raise Error, error
         else:
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No SolvAccVol.read_top_in_binary methods take ' + nargs
               raise Error, error
            else:
               raise e
      return val
   def write_top_in_ascii(*args):
      trydefault = 0
      try:
         if len(args)==2:
            val = apply(MEADc.SolvAccVol_write_top_in_ascii__L31,args)
            if val and isinstance(val, types.StringType) and re.search('_p_void', val):
               val = voidPtr(val) ; val.thisown = 1
         else:
            trydefault = 1
            val = apply(MEADc.SolvAccVol_write_top_in_ascii,args)
            if val and isinstance(val, types.StringType):
               mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
               if mo:
                  rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                  val = eval(rettype) ; val.thisown = 1
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               val = apply(MEADc.SolvAccVol_write_top_in_ascii,args)
               if val and isinstance(val, types.StringType):
                  mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
                  if mo:
                     rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                     val = eval(rettype) ; val.thisown = 1
            except StandardError, e:
               error = 'No SolvAccVol.write_top_in_ascii methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               raise Error, error
         else:
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No SolvAccVol.write_top_in_ascii methods take ' + nargs
               raise Error, error
            else:
               raise e
      return val
   def read_top_in_ascii(*args):
      trydefault = 0
      try:
         if len(args)==2:
            val = apply(MEADc.SolvAccVol_read_top_in_ascii__L33,args)
            if val and isinstance(val, types.StringType) and re.search('_p_int', val):
               val = intPtr(val) ; val.thisown = 1
         else:
            trydefault = 1
            val = apply(MEADc.SolvAccVol_read_top_in_ascii,args)
            if val and isinstance(val, types.StringType):
               mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
               if mo:
                  rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                  val = eval(rettype) ; val.thisown = 1
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               val = apply(MEADc.SolvAccVol_read_top_in_ascii,args)
               if val and isinstance(val, types.StringType):
                  mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
                  if mo:
                     rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                     val = eval(rettype) ; val.thisown = 1
            except StandardError, e:
               error = 'No SolvAccVol.read_top_in_ascii methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               raise Error, error
         else:
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No SolvAccVol.read_top_in_ascii methods take ' + nargs
               raise Error, error
            else:
               raise e
      return val
   def __init__(self,*args):
      trydefault = 0
      try:
         if len(args)==1:
            self.this = apply(MEADc.new_SolvAccVol__L17,args)
         elif len(args)==2:
            self.this = apply(MEADc.new_SolvAccVol__L16,args)
         else:
            trydefault = 1
            self.this = apply(MEADc.new_SolvAccVol,args)
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               self.this = apply(MEADc.new_SolvAccVol,args)
            except StandardError, e:
               error = 'No SolvAccVol methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               self.thisown = 0
               raise Error, error
         else:
            self.thisown = 0
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No SolvAccVol methods take ' + nargs
               raise Error, error
            else:
               raise e
      self.thisown = 1
      return
   def write_top_in_binary(*args):
      trydefault = 0
      try:
         if len(args)==2:
            val = apply(MEADc.SolvAccVol_write_top_in_binary__L30,args)
            if val and isinstance(val, types.StringType) and re.search('_p_void', val):
               val = voidPtr(val) ; val.thisown = 1
         else:
            trydefault = 1
            val = apply(MEADc.SolvAccVol_write_top_in_binary,args)
            if val and isinstance(val, types.StringType):
               mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
               if mo:
                  rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                  val = eval(rettype) ; val.thisown = 1
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               val = apply(MEADc.SolvAccVol_write_top_in_binary,args)
               if val and isinstance(val, types.StringType):
                  mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
                  if mo:
                     rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                     val = eval(rettype) ; val.thisown = 1
            except StandardError, e:
               error = 'No SolvAccVol.write_top_in_binary methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               raise Error, error
         else:
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No SolvAccVol.write_top_in_binary methods take ' + nargs
               raise Error, error
            else:
               raise e
      return val
   def __deepcopy__(self, memo = None):
      val = MEADc.new_SolvAccVol__deepcopy(self)
      if val: val = SolvAccVolPtr(val) ; val.thisown = 1
      return val

SolvAccVol.read_top_in_binary = __Dummy.__dict__['read_top_in_binary']
del SolvAccVol.__dict__['read_top_in_binary__L32']
SolvAccVol.write_top_in_ascii = __Dummy.__dict__['write_top_in_ascii']
del SolvAccVol.__dict__['write_top_in_ascii__L31']
SolvAccVol.read_top_in_ascii = __Dummy.__dict__['read_top_in_ascii']
del SolvAccVol.__dict__['read_top_in_ascii__L33']
SolvAccVol.__init__ = __Dummy.__dict__['__init__']
SolvAccVol.write_top_in_binary = __Dummy.__dict__['write_top_in_binary']
del SolvAccVol.__dict__['write_top_in_binary__L30']
SolvAccVol.__deepcopy__ = __Dummy.__dict__['__deepcopy__']

import types, re
class __Dummy:
   def __mul__(*args):
      trydefault = 0
      try:
         if len(args)==2:
            val = apply(MEADc.Poly___mul____L82,args)
            if val: val = PolyPtr(val) ; val.thisown = 1
         else:
            trydefault = 1
            val = apply(MEADc.Poly___mul__,args)
            if val and isinstance(val, types.StringType):
               mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
               if mo:
                  rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                  val = eval(rettype) ; val.thisown = 1
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               val = apply(MEADc.Poly___mul__,args)
               if val and isinstance(val, types.StringType):
                  mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
                  if mo:
                     rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                     val = eval(rettype) ; val.thisown = 1
            except StandardError, e:
               error = 'No Poly.__mul__ methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               raise Error, error
         else:
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No Poly.__mul__ methods take ' + nargs
               raise Error, error
            else:
               raise e
      return val
   def __init__(self,*args):
      trydefault = 0
      try:
         if len(args)==1:
            try:
               self.this = apply(MEADc.new_Poly__L19,args)
            except:
               try:
                  self.this = apply(MEADc.new_Poly__L18,args)
               except:
                  self.this = apply(MEADc.new_Poly__L16,args)
         elif len(args)==2:
            self.this = apply(MEADc.new_Poly__L17,args)
         else:
            trydefault = 1
            self.this = apply(MEADc.new_Poly,args)
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               self.this = apply(MEADc.new_Poly,args)
            except StandardError, e:
               error = 'No Poly methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               self.thisown = 0
               raise Error, error
         else:
            self.thisown = 0
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No Poly methods take ' + nargs
               raise Error, error
            else:
               raise e
      self.thisown = 1
      return
   def __deepcopy__(self, memo = None):
      val = MEADc.new_Poly__deepcopy(self)
      if val: val = PolyPtr(val) ; val.thisown = 1
      return val

Poly.__mul__ = __Dummy.__dict__['__mul__']
del Poly.__dict__['_Poly__mul____L82']
Poly.__init__ = __Dummy.__dict__['__init__']
Poly.__deepcopy__ = __Dummy.__dict__['__deepcopy__']

class __Dummy:
   def __init__(self,*args):
      trydefault = 0
      try:
         if len(args)==1:
            try:
               self.this = apply(MEADc.new_Legendre__L41,args)
            except:
               self.this = apply(MEADc.new_Legendre__L42,args)
         else:
            trydefault = 1
            self.this = apply(MEADc.new_Legendre,args)
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               self.this = apply(MEADc.new_Legendre,args)
            except StandardError, e:
               error = 'No Legendre methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               self.thisown = 0
               raise Error, error
         else:
            self.thisown = 0
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No Legendre methods take ' + nargs
               raise Error, error
            else:
               raise e
      self.thisown = 1
      return
   def __deepcopy__(self, memo = None):
      val = MEADc.new_Legendre__deepcopy(self)
      if val: val = LegendrePtr(val) ; val.thisown = 1
      return val

Legendre.__init__ = __Dummy.__dict__['__init__']
Legendre.__deepcopy__ = __Dummy.__dict__['__deepcopy__']

# Class to implement the CenteringStyle enum
# Define the enum values in the dict from a SWIG generated CenteringStyle_enum class.
# Don't allow these enum attributes to be set
# (But allow other attributes to be added)

class CenteringStyle_enum(CenteringStyle_enum):
   def __init__(self):
      self.this = "None"
      self.__dict__.update({'ON_ORIGIN' : CenteringStyle_enum.ON_ORIGIN, 'ON_CENT_OF_INTR' : CenteringStyle_enum.ON_CENT_OF_INTR, 'ON_GEOM_CENT' : CenteringStyle_enum.ON_GEOM_CENT, 'SPECIFIC' : CenteringStyle_enum.SPECIFIC})
   def __setattr__(self, name, value):
      if (name == "ON_ORIGIN") or (name == "ON_CENT_OF_INTR") or (name == "ON_GEOM_CENT") or (name == "SPECIFIC"):
         error = "Can't set constant enum attribute " + name
         raise AttributeError, error
      else:
         self.__dict__[name] = value

# Create an instance
CenteringStyle = CenteringStyle_enum()
import types, re
class __Dummy:
   def __init__(self,*args):
      trydefault = 0
      try:
         if len(args)==1:
            self.this = apply(MEADc.new_CubeLatSpec__L17,args)
         elif len(args)==3:
            self.this = apply(MEADc.new_CubeLatSpec__L33,args)
         else:
            trydefault = 1
            self.this = apply(MEADc.new_CubeLatSpec,args)
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               self.this = apply(MEADc.new_CubeLatSpec,args)
            except StandardError, e:
               error = 'No CubeLatSpec methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               self.thisown = 0
               raise Error, error
         else:
            self.thisown = 0
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No CubeLatSpec methods take ' + nargs
               raise Error, error
            else:
               raise e
      self.thisown = 1
      return
   def __deepcopy__(self, memo = None):
      val = MEADc.new_CubeLatSpec__deepcopy(self)
      if val: val = CubeLatSpecPtr(val) ; val.thisown = 1
      return val

CubeLatSpec.__init__ = __Dummy.__dict__['__init__']
CubeLatSpec.__deepcopy__ = __Dummy.__dict__['__deepcopy__']

# Create an instance of the PhysCond static class
# Create an instance of the Blab class.

PhysCond = PhysCond()
Blab = Blab()
import types, re
import types, re
class __Dummy:
   def add_level(*args):
      trydefault = 0
      try:
         if len(args)==4:
            val = apply(MEADc.FinDiffMethod_add_level__L29,args)
            if val and isinstance(val, types.StringType) and re.search('_p_void', val):
               val = voidPtr(val) ; val.thisown = 1
         else:
            trydefault = 1
            val = apply(MEADc.FinDiffMethod_add_level,args)
            if val and isinstance(val, types.StringType):
               mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
               if mo:
                  rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                  val = eval(rettype) ; val.thisown = 1
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               val = apply(MEADc.FinDiffMethod_add_level,args)
               if val and isinstance(val, types.StringType):
                  mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
                  if mo:
                     rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                     val = eval(rettype) ; val.thisown = 1
            except StandardError, e:
               error = 'No FinDiffMethod.add_level methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               raise Error, error
         else:
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No FinDiffMethod.add_level methods take ' + nargs
               raise Error, error
            else:
               raise e
      return val
   def __init__(self,*args):
      trydefault = 0
      try:
         if len(args)==1:
            self.this = apply(MEADc.new_FinDiffMethod__L16,args)
         else:
            trydefault = 1
            self.this = apply(MEADc.new_FinDiffMethod,args)
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               self.this = apply(MEADc.new_FinDiffMethod,args)
            except StandardError, e:
               error = 'No FinDiffMethod methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               self.thisown = 0
               raise Error, error
         else:
            self.thisown = 0
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No FinDiffMethod methods take ' + nargs
               raise Error, error
            else:
               raise e
      self.thisown = 1
      return
   def __deepcopy__(self, memo = None):
      val = MEADc.new_FinDiffMethod__deepcopy(self)
      if val: val = FinDiffMethodPtr(val) ; val.thisown = 1
      return val

FinDiffMethod.add_level = __Dummy.__dict__['add_level']
del FinDiffMethod.__dict__['add_level__L29']
FinDiffMethod.__init__ = __Dummy.__dict__['__init__']
FinDiffMethod.__deepcopy__ = __Dummy.__dict__['__deepcopy__']

import types, re
class __Dummy:
   def __mul__(*args):
      trydefault = 0
      try:
         if len(args)==2:
            val = apply(MEADc.AnalySphere___mul____L69,args)
            if val and isinstance(val, types.StringType) and re.search('_p_ElstatPotCombination', val):
               val = ElstatPotCombinationPtr(val) ; val.thisown = 1
         else:
            trydefault = 1
            val = apply(MEADc.AnalySphere___mul__,args)
            if val and isinstance(val, types.StringType):
               mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
               if mo:
                  rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                  val = eval(rettype) ; val.thisown = 1
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               val = apply(MEADc.AnalySphere___mul__,args)
               if val and isinstance(val, types.StringType):
                  mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
                  if mo:
                     rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                     val = eval(rettype) ; val.thisown = 1
            except StandardError, e:
               error = 'No AnalySphere.__mul__ methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               raise Error, error
         else:
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No AnalySphere.__mul__ methods take ' + nargs
               raise Error, error
            else:
               raise e
      return val
   def __init__(self,*args):
      trydefault = 0
      try:
         if len(args)==3:
            self.this = apply(MEADc.new_AnalySphere__L33,args)
         else:
            trydefault = 1
            self.this = apply(MEADc.new_AnalySphere,args)
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               self.this = apply(MEADc.new_AnalySphere,args)
            except StandardError, e:
               error = 'No AnalySphere methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               self.thisown = 0
               raise Error, error
         else:
            self.thisown = 0
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No AnalySphere methods take ' + nargs
               raise Error, error
            else:
               raise e
      self.thisown = 1
      return

AnalySphere.__mul__ = __Dummy.__dict__['__mul__']
del AnalySphere.__dict__['_AnalySphere__mul____L69']
AnalySphere.__init__ = __Dummy.__dict__['__init__']

import types, re
class __Dummy:
   def __init__(self,*args):
      trydefault = 0
      try:
         if len(args)==3:
            self.this = apply(MEADc.new_AnalySlab__L35,args)
         else:
            trydefault = 1
            self.this = apply(MEADc.new_AnalySlab,args)
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               self.this = apply(MEADc.new_AnalySlab,args)
            except StandardError, e:
               error = 'No AnalySlab methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               self.thisown = 0
               raise Error, error
         else:
            self.thisown = 0
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No AnalySlab methods take ' + nargs
               raise Error, error
            else:
               raise e
      self.thisown = 1
      return
   def __mul__(*args):
      trydefault = 0
      try:
         if len(args)==2:
            val = apply(MEADc.AnalySlab___mul____L71,args)
            if val and isinstance(val, types.StringType) and re.search('_p_ElstatPotCombination', val):
               val = ElstatPotCombinationPtr(val) ; val.thisown = 1
         else:
            trydefault = 1
            val = apply(MEADc.AnalySlab___mul__,args)
            if val and isinstance(val, types.StringType):
               mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
               if mo:
                  rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                  val = eval(rettype) ; val.thisown = 1
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               val = apply(MEADc.AnalySlab___mul__,args)
               if val and isinstance(val, types.StringType):
                  mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
                  if mo:
                     rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                     val = eval(rettype) ; val.thisown = 1
            except StandardError, e:
               error = 'No AnalySlab.__mul__ methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               raise Error, error
         else:
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No AnalySlab.__mul__ methods take ' + nargs
               raise Error, error
            else:
               raise e
      return val

AnalySlab.__init__ = __Dummy.__dict__['__init__']
AnalySlab.__mul__ = __Dummy.__dict__['__mul__']
del AnalySlab.__dict__['_AnalySlab__mul____L71']

import types, re
class __Dummy:
   def __mul__(*args):
      trydefault = 0
      try:
         if len(args)==2:
            val = apply(MEADc.Debye___mul____L47,args)
            if val and isinstance(val, types.StringType) and re.search('_p_ElstatPotCombination', val):
               val = ElstatPotCombinationPtr(val) ; val.thisown = 1
         else:
            trydefault = 1
            val = apply(MEADc.Debye___mul__,args)
            if val and isinstance(val, types.StringType):
               mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
               if mo:
                  rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                  val = eval(rettype) ; val.thisown = 1
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               val = apply(MEADc.Debye___mul__,args)
               if val and isinstance(val, types.StringType):
                  mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
                  if mo:
                     rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                     val = eval(rettype) ; val.thisown = 1
            except StandardError, e:
               error = 'No Debye.__mul__ methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               raise Error, error
         else:
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No Debye.__mul__ methods take ' + nargs
               raise Error, error
            else:
               raise e
      return val

Debye.__mul__ = __Dummy.__dict__['__mul__']
del Debye.__dict__['_Debye__mul____L47']

import types, re
class __Dummy:
   def __mul__(*args):
      trydefault = 0
      try:
         if len(args)==2:
            val = apply(MEADc.FinDiffElstatPot___mul____L97,args)
            if val and isinstance(val, types.StringType) and re.search('_p_ElstatPotCombination', val):
               val = ElstatPotCombinationPtr(val) ; val.thisown = 1
         else:
            trydefault = 1
            val = apply(MEADc.FinDiffElstatPot___mul__,args)
            if val and isinstance(val, types.StringType):
               mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
               if mo:
                  rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                  val = eval(rettype) ; val.thisown = 1
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               val = apply(MEADc.FinDiffElstatPot___mul__,args)
               if val and isinstance(val, types.StringType):
                  mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
                  if mo:
                     rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                     val = eval(rettype) ; val.thisown = 1
            except StandardError, e:
               error = 'No FinDiffElstatPot.__mul__ methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               raise Error, error
         else:
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No FinDiffElstatPot.__mul__ methods take ' + nargs
               raise Error, error
            else:
               raise e
      return val
   def __init__(self,*args):
      trydefault = 0
      try:
         if len(args)==3:
            self.this = apply(MEADc.new_FinDiffElstatPot__L17,args)
         else:
            trydefault = 1
            self.this = apply(MEADc.new_FinDiffElstatPot,args)
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               self.this = apply(MEADc.new_FinDiffElstatPot,args)
            except StandardError, e:
               error = 'No FinDiffElstatPot methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               self.thisown = 0
               raise Error, error
         else:
            self.thisown = 0
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No FinDiffElstatPot methods take ' + nargs
               raise Error, error
            else:
               raise e
      self.thisown = 1
      return

FinDiffElstatPot.__mul__ = __Dummy.__dict__['__mul__']
del FinDiffElstatPot.__dict__['_FinDiffElstatPot__mul____L97']
FinDiffElstatPot.__init__ = __Dummy.__dict__['__init__']

import types, re
class __Dummy:
   def __mul__(*args):
      trydefault = 0
      try:
         if len(args)==2:
            val = apply(MEADc.ElstatPotCombination___mul____L106,args)
            if val and isinstance(val, types.StringType) and re.search('_p_float', val):
               val = floatPtr(val) ; val.thisown = 1
         else:
            trydefault = 1
            val = apply(MEADc.ElstatPotCombination___mul__,args)
            if val and isinstance(val, types.StringType):
               mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
               if mo:
                  rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                  val = eval(rettype) ; val.thisown = 1
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               val = apply(MEADc.ElstatPotCombination___mul__,args)
               if val and isinstance(val, types.StringType):
                  mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
                  if mo:
                     rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                     val = eval(rettype) ; val.thisown = 1
            except StandardError, e:
               error = 'No ElstatPotCombination.__mul__ methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               raise Error, error
         else:
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No ElstatPotCombination.__mul__ methods take ' + nargs
               raise Error, error
            else:
               raise e
      return val
   def __init__(self,*args):
      trydefault = 0
      try:
         if len(args)==1:
            self.this = apply(MEADc.new_ElstatPotCombination__L16,args)
         else:
            trydefault = 1
            self.this = apply(MEADc.new_ElstatPotCombination,args)
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               self.this = apply(MEADc.new_ElstatPotCombination,args)
            except StandardError, e:
               error = 'No ElstatPotCombination methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               self.thisown = 0
               raise Error, error
         else:
            self.thisown = 0
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No ElstatPotCombination methods take ' + nargs
               raise Error, error
            else:
               raise e
      self.thisown = 1
      return
   def __add__(*args):
      trydefault = 0
      try:
         if len(args)==2:
            val = apply(MEADc.ElstatPotCombination___add____L98,args)
            if val: val = ElstatPotCombinationPtr(val) ; val.thisown = 1
         else:
            trydefault = 1
            val = apply(MEADc.ElstatPotCombination___add__,args)
            if val and isinstance(val, types.StringType):
               mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
               if mo:
                  rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                  val = eval(rettype) ; val.thisown = 1
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               val = apply(MEADc.ElstatPotCombination___add__,args)
               if val and isinstance(val, types.StringType):
                  mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
                  if mo:
                     rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                     val = eval(rettype) ; val.thisown = 1
            except StandardError, e:
               error = 'No ElstatPotCombination.__add__ methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               raise Error, error
         else:
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No ElstatPotCombination.__add__ methods take ' + nargs
               raise Error, error
            else:
               raise e
      return val
   def __iadd__(*args):
      trydefault = 0
      try:
         if len(args)==2:
            val = apply(MEADc.ElstatPotCombination___iadd____L90,args)
            if val: val = ElstatPotCombinationPtr(val) ; val.thisown = 1
         else:
            trydefault = 1
            val = apply(MEADc.ElstatPotCombination___iadd__,args)
            if val and isinstance(val, types.StringType):
               mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
               if mo:
                  rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                  val = eval(rettype) ; val.thisown = 1
      except StandardError, e:
         if not trydefault:
            laste = e
            try:
               val = apply(MEADc.ElstatPotCombination___iadd__,args)
               if val and isinstance(val, types.StringType):
                  mo = re.search(r'^_[\da-fA-F]*_p_(.*)', val)
                  if mo:
                     rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'
                     val = eval(rettype) ; val.thisown = 1
            except StandardError, e:
               error = 'No ElstatPotCombination.__iadd__ methods apply using the given arguments.\n'
               if re.search('arguments', e.args[0]):
                  error = error + 'Possible error is ' + laste.args[0]
               else:
                  error = 'One possible error is ' + e.args[0] + '\n'
                  error = error + 'Another possible error is ' + laste.args[0]
               raise Error, error
         else:
            if re.search('arguments', e.args[0]):
               if len(args)==1:
                  nargs = '1 argument'
               else:
                  nargs = str(len(args)) + ' arguments'
               error = 'No ElstatPotCombination.__iadd__ methods take ' + nargs
               raise Error, error
            else:
               raise e
      return val
   def __deepcopy__(self, memo = None):
      val = MEADc.new_ElstatPotCombination__deepcopy(self)
      if val: val = ElstatPotCombinationPtr(val) ; val.thisown = 1
      return val

ElstatPotCombination.__mul__ = __Dummy.__dict__['__mul__']
del ElstatPotCombination.__dict__['_ElstatPotCombination__mul____L106']
ElstatPotCombination.__init__ = __Dummy.__dict__['__init__']
ElstatPotCombination.__add__ = __Dummy.__dict__['__add__']
del ElstatPotCombination.__dict__['_ElstatPotCombination__add____L98']
ElstatPotCombination.__iadd__ = __Dummy.__dict__['__iadd__']
del ElstatPotCombination.__dict__['_ElstatPotCombination__iadd____L90']
ElstatPotCombination.__deepcopy__ = __Dummy.__dict__['__deepcopy__']

import types, re
