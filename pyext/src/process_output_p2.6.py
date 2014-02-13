#! /usr/bin/env python

from optparse import OptionParser
from optparse import Option, OptionValueError

import difflib
import cPickle

class MultipleOption(Option):
    ACTIONS = Option.ACTIONS + ("extend",)
    STORE_ACTIONS = Option.STORE_ACTIONS + ("extend",)
    TYPED_ACTIONS = Option.TYPED_ACTIONS + ("extend",)
    ALWAYS_TYPED_ACTIONS = Option.ALWAYS_TYPED_ACTIONS + ("extend",)

    def take_action(self, action, dest, opt, value, values, parser):
        if action == "extend":
            values.ensure_value(dest, []).append(value)
        else:
            Option.take_action(self, action, dest, opt, value, values, parser)

parser = OptionParser(option_class=MultipleOption,usage='Process output data file saved as dictionaries. It has two modality: print selected fields for all lines or print a particular line where a filed has a given value. Example of usage: process_output.py --soft -s To -s E -s S -f log.3.native-2-no-red. process_output.py --soft --search_field EV0 --search_value 5.67750116023 -f log.3.native-2-no-red')
parser.add_option('-f', action="store", dest="filename", help="file name to process" )
parser.add_option('-s', dest="fields", action="extend", type="string", help="Specify all fields to be printed. Multiple flags will append a list of fields to be printed" )
parser.add_option('-t', dest="single_column_field", help="Specify a single column field to be printed. It will be printed as a column. If the field name is not complete, it will print all fields whose name contain the queried string." )
parser.add_option('-p', action="store_true", dest="print_fields", default=False, help="print the fields contained in the file")
parser.add_option('-n', action="store", dest="print_raw_number", help="print the selected raw" )
parser.add_option('--soft', action="store_true", dest="soft_match", default=False, help="Soft match. Closest matching field will be printed, e.g. S will give Step_Number, En will give energy, etc. ")
parser.add_option('--search_field', dest="search_field", help="Search a line from the file. Specify the field to be searched for. ")
parser.add_option('--search_value', dest="search_value", help="Search a line from the file. Specify the value to be searched for. ")

(result,args)=parser.parse_args()

#open the file
if result.filename!=None:
   f=open(result.filename,"r")
else:
   print "Error: No file name provided. Use -h for help"
   exit()


#get the keys from the first line
for line in f.readlines():
    d=eval(line)
    klist=d.keys()
    klist.sort()
    break
f.close()    
    
#print the keys    
if result.print_fields:    
    for key in klist:
        print key

#the field string matching is by default strict, i.e., the input string must be the same as the one in the file
match_strictness=1.0
if result.soft_match: match_strictness=0.1

#print the queried fields
if result.fields!=None:
   field_list=[]
   #check whether the fields exist and convert them to best maching existing field names
   for field in result.fields:
      found_entries=difflib.get_close_matches(field,klist,1,match_strictness)
      if len(found_entries)==0:
         print "Error: field "+field+" non found"
         exit()
      else:
         field_list.append(found_entries[0])
   
   #print comment line   
   s0=' '.join(["%20s" % (field) for field in field_list])
   print "# "+s0
   
   #print fields values
   f=open(result.filename,"r")
   line_number=0
   for line in f.readlines():
      line_number+=1
      try:
         d=eval(line)
      except:
         print "# Warning: skipped line number " + str(line_number) + " not a valid line"
         continue
      s0=' '.join(["%20s" % (str(d[field])) for field in field_list])
      print "> "+s0
   f.close()

if result.single_column_field!=None:
   field_list=[]
   for k in klist:
       if result.single_column_field in k:
           field_list.append(k)
              
   f=open(result.filename,"r")
   line_number=0
   for line in f.readlines():
      try:
         d=eval(line)
      except:
         print "# Warning: skipped line number " + str(line_number) + " not a valid line"
         continue   
      for key in field_list:
         print key, d[key]
      print " "               
   f.close()      

if (result.search_field!=None and result.search_value!=None):
   #check whether the fields exist and convert them to best maching existing field names
   found_entries=difflib.get_close_matches(result.search_field,klist,1,match_strictness)
   if len(found_entries)==0:
       print "Error: field "+results.search_field+" non found"
       exit()
   else:
       corrected_field=found_entries[0]
   #print fields values
   f=open(result.filename,"r")
   line_number=0
   for line in f.readlines():
      line_number+=1
      try:
         d=eval(line)
      except:
         print "# Warning: skipped line number " + str(line_number) + " not a valid line"
         continue
      if (str(d[corrected_field])==result.search_value):
         for key in klist:
             print key, d[key]
   f.close()

if (result.print_raw_number!=None):
   #check whether the fields exist and convert them to best maching existing field names
   f=open(result.filename,"r")
   line_number=0
   for line in f.readlines():
      line_number+=1
      if (line_number==int(result.print_raw_number)):
        try:
          d=eval(line)
        except:
          print "# Warning: skipped line number " + str(line_number) + " not a valid line"
          break
        for key in klist:
           print key, d[key]
   f.close()
       


