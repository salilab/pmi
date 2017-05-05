class Segments(object):
    def __init__(self,index):
        if type(index) is int:
            self.segs=[[index]]
        elif type(index) is list:
            self.segs=[[index[0]]]
            for i in index[1:]:
                self.add(i)

    def add(self,index):
        mergeleft=None
        mergeright=None
        for n,s in enumerate(self.segs):
            if index in s:
                pass
            else:
                if s[0]-index==1:
                    mergeleft=n
                if index-s[-1]==1:
                    mergeright=n
        if mergeright is None and mergeleft is None:
            self.segs.append([index])
        if not mergeright is None and mergeleft is None:
            self.segs[mergeright].append(index)
        if not mergeleft is None and  mergeright is None:
            self.segs[mergeleft]=[index]+self.segs[mergeleft]
        if not mergeleft is None and not mergeright is None:
            self.segs[mergeright]=self.segs[mergeright]+[index]+self.segs[mergeleft]
            del self.segs[mergeleft]

        for n in range(len(self.segs)):
            self.segs[n].sort()

        self.segs.sort(key=lambda tup: tup[0])

    def remove(self,index):
        for n,s in enumerate(self.segs):
            if index in s:
                if s[0]==index:
                    self.segs[n]=s[1:]
                elif s[-1]==index:
                    self.segs[n]=s[:-1]
                else:
                    i=self.segs[n].index(index)
                    self.segs[n]=s[:i]
                    self.segs.append(s[i+1:])
        for n in range(len(self.segs)):
            self.segs[n].sort()
        self.segs.sort(key=lambda tup: tup[0])


s=Segments(1)
print s.segs
s.add(2)
print s.segs
s.add(0)
print s.segs
s.add(4)
print s.segs
s.add(-3)
print s.segs
s.add(3)
print s.segs
s.add(-1)
print s.segs
s.add(-2)
print s.segs
s.remove(0)
print s.segs
s.remove(1)
print s.segs
s.remove(-1)
print s.segs
s.remove(3)
print s.segs
s.remove(-3)
print s.segs

s=Segments([1,2,-1,4,9,-2])
print s.segs
