'''
---------------------------
 Licensing and Distribution
---------------------------

Program name: TorsiFlex
Version     : 2021.2
License     : MIT/x11

Copyright (c) 2021, David Ferro Costas (david.ferro@usc.es) and
Antonio Fernandez Ramos (qf.ramos@usc.es)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the Software
is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.
---------------------------

*----------------------------------*
| Module     :  modtorsiflex       |
| Sub-module :  tpespoint          |
| Last Update:  2020/12/21 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains the TorPESpoint
class for TorsiFlex
'''

#===================================================#
import common.fncs as fncs
#===================================================#

#===================================================#
#       Torsional vector in different formats       #
#===================================================#
def ivec2fvec(ivec): return [float(angle)      for angle in ivec]
def fvec2ivec(fvec): return [int(round(angle)) for angle in fvec]
def svec2fvec(svec): return [float(angle)      for angle in svec.split("_")]
def svec2ivec(svec): return [int(angle)        for angle in svec.split("_")]
def ivec2svec(ivec): return "_".join([ "%003i"%(round(angle)%360)  for angle in ivec])
def fvec2svec(fvec): return "_".join(["%003.f"%(round(angle)%360)  for angle in fvec])
#===================================================#

class TorPESpoint():
      ''' angles in degrees '''

      def __init__(self,angles,cycle=None):
          '''
          _svec: string version of the vector
          _ivec: vector of integers
          _fvec: vector of floats
          '''

          if type(angles) == str: self._fvec = svec2fvec(angles)
          else                  : self._fvec = ivec2fvec(angles)
          # dimension of the vector
          self._dim  = len(self._fvec)
          # limits for each angle
          if cycle is None: self._cycle = [360 for idx in range(self._dim)]
          else            : self._cycle = cycle
          self._fvec = [vi%li for vi,li in zip(self._fvec,self._cycle)]
          # Now, get the other versions of the vector
          self._ivec = fvec2ivec(self._fvec)
          self._svec = ivec2svec(self._ivec)

      def __str__(self): return self._svec

      def fvec_in_180(self):
         fvec = [angle if angle <= 180 else angle-360.0 for angle in self._fvec]
         return fvec

      def abs_diff(self,point):
          diff = [None for idx in range(self._dim)]
          for idx in range(self._dim):
              args   = (self._fvec[idx],point._fvec[idx])
              kwargs = {"u": "deg", "limit": self._cycle[idx]}
              diff[idx] = fncs.angular_dist(*args,**kwargs)
          return diff

      def distance_if_smaller(self,point,ref=float("inf")):
          refe2 = ref**2
          dist2 = 0.0
          for idx in range(self._dim):
              args   = (self._fvec[idx],point._fvec[idx])
              kwargs = {"u": "deg", "limit": self._cycle[idx]}
              dist2 += fncs.angular_dist(*args,**kwargs)**2
              if dist2 > refe2: return float("inf")
          return dist2**0.5

      def distance(self,point):
          return self.distance_if_smaller(point,ref=float("inf"))

      def max1Ddistance(self,point):
          return max(self.abs_diff(point))

      def is_same(self,point,eps):
          for idx in range(self._dim):
              args   = (self._fvec[idx],point._fvec[idx])
              kwargs = {"u": "deg", "limit": self._cycle[idx]}
              diff1D = fncs.angular_dist(*args,**kwargs)
              if diff1D > eps: return False
          return True

      def is_inside(self,point,eps):
          return self.is_same(point,eps)

      def is_in_list(self,points,eps):
          for point in points:
              if self.is_same(point,eps): return True, point
          return False, None

      def is_in_domain(self,tdomain):
          for angle,domain in zip(self._fvec,tdomain):
              if not fncs.float_in_domain(angle,domain): return False
          return True

      def closest(self,points):
          if len(points) == 0: return None
          # initialize
          mindist = float("inf")
          closest = None
          # compare with all points (brute force)
          for point in points:
              dist = self.distance_if_smaller(point,mindist)
              if dist < mindist: mindist,closest = dist,point
          return closest

      def nomenclature(self):
          name = ""
          for angle in self._ivec:
              angle = angle%360
              # return name for this angle
              if          angle ==   0: name += "C"
              elif   0  < angle <=  30: name += "C+"
              elif  30  < angle <=  75: name += "G+"
              elif  75  < angle <=  90: name += "g+"
              elif  90  < angle <= 105: name += "a+"
              elif 105  < angle <= 150: name += "A+"
              elif 150  < angle <  180: name += "T+"
              elif        angle == 180: name += "T"
              elif 180  < angle <  210: name += "T-"
              elif 210 <= angle <  255: name += "A-"
              elif 255 <= angle <  270: name += "a-"
              elif 270 <= angle <  285: name += "g-"
              elif 285 <= angle <  330: name += "G-"
              elif 330 <= angle <  360: name += "C-"
          return name


def tests():
    cycle = (360,360,120)
    pp1 = TorPESpoint( (90,180,+50) ,cycle)
    print("pp1:", pp1)
    print("nomenclature:",pp1.nomenclature())
    print("")
    pp2 = TorPESpoint( (85,190,-40),cycle )
    print("pp2:", pp2)
    print("")
    print("pp1 - pp2:", pp1.abs_diff(pp2))
    print("pp2 - pp1:", pp2.abs_diff(pp1))
    print("max diff:",pp1.max1Ddistance(pp2))
    print("distance:",pp1.distance(pp2))
    print("distance:",pp1.distance_if_smaller(pp2,1))
    print("")
    print("")
    pp3 = TorPESpoint( (270,180,-50),cycle )
    print("diff:", pp1.abs_diff(pp3))
    suma = [a+b for a,b in zip(pp1._fvec,pp3._fvec)]
    pp4 = TorPESpoint(suma,cycle)
    print("sum: ",pp4)



if __name__ == '__main__': tests()

