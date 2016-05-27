"""@namespace IMP.pmi.mmcif
   @brief Support for the mmCIF file format.
"""

from __future__ import print_function

class _LineWriter(object):
    def __init__(self, writer, line_len=80, multi_line_len=70):
        self.writer = writer
        self.line_len = line_len
        self.multi_line_len = multi_line_len
        self.column = 0
    def write(self, val):
        if isinstance(val, str) and len(val) > self.multi_line_len:
            self.writer.fh.write("\n;")
            for i in range(0, len(val), self.multi_line_len):
                self.writer.fh.write(val[i:i+self.multi_line_len])
                self.writer.fh.write("\n")
            self.writer.fh.write(";\n")
            self.column = 0
            return
        val = self.writer._repr(val)
        if self.column > 0:
            if self.column + len(val) + 1 > self.line_len:
                self.writer.fh.write("\n")
                self.column = 0
            else:
                self.writer.fh.write(" ")
                self.column += 1
        self.writer.fh.write(val)
        self.column += len(val)


class CifCategoryWriter(object):
    def __init__(self, writer, category):
        self.writer = writer
        self.category = category
    def write(self, **kwargs):
        self.writer._write(self.category, kwargs)
    def __enter__(self):
        return self
    def __exit__(self, exc_type, exc_value, traceback):
        pass


class CifLoopWriter(object):
    def __init__(self, writer, category, keys):
        self.writer = writer
        self.category = category
        self.keys = keys
    def write(self, **kwargs):
        l = _LineWriter(self.writer)
        for k in self.keys:
            l.write(kwargs.get(k, self.writer.omitted))
        self.writer.fh.write("\n")
    def __enter__(self):
        f = self.writer.fh
        f.write("#\nloop_\n")
        for k in self.keys:
            f.write("%s.%s\n" % (self.category, k))
        return self
    def __exit__(self, exc_type, exc_value, traceback):
        self.writer.fh.write("#\n")


class CifWriter(object):
    omitted = '.'
    unknown = '?'

    def __init__(self, fh):
        self.fh = fh
    def category(self, category):
        return CifCategoryWriter(self, category)
    def loop(self, category, keys):
        return CifLoopWriter(self, category, keys)
    def _write(self, category, kwargs):
        for key in kwargs:
            self.fh.write("%s.%s %s\n" % (category, key,
                                          self._repr(kwargs[key])))
    def _repr(self, obj):
        if isinstance(obj, str) and '"' not in obj \
           and "'" not in obj and " " not in obj:
            return obj
        elif isinstance(obj, float):
            return "%.3f" % obj
        else:
            return repr(obj)
