import ROOT
import os

def loadTree(mytree, directory, category, year ):

   if (year == '_12022'):
      mytree.Add(directory+'snapshot_mc-11'+year+category+'.root')
      mytree.Add(directory+'snapshot_mc-12'+year+category+'.root')
      mytree.Add(directory+'snapshot_mc-13'+year+category+'.root')
      mytree.Add(directory+'snapshot_mc-14'+year+category+'.root')                  

   if (year == '_22022'):
      mytree.Add(directory+'snapshot_mc-15'+year+category+'.root')
      mytree.Add(directory+'snapshot_mc-16'+year+category+'.root')
      mytree.Add(directory+'snapshot_mc-17'+year+category+'.root')

   if (year == '_12023'):
      mytree.Add(directory+'snapshot_mc-21'+year+category+'.root')
      mytree.Add(directory+'snapshot_mc-22'+year+category+'.root')
      mytree.Add(directory+'snapshot_mc-23'+year+category+'.root')
      mytree.Add(directory+'snapshot_mc-24'+year+category+'.root')

   if (year == '_22023'):
      mytree.Add(directory+'snapshot_mc-31'+year+category+'.root')
      mytree.Add(directory+'snapshot_mc-32'+year+category+'.root')      
      
   mytree.Add(directory+'snapshot_mc10'+year+category+'.root')
   mytree.Add(directory+'snapshot_mc11'+year+category+'.root')
   mytree.Add(directory+'snapshot_mc12'+year+category+'.root')
   mytree.Add(directory+'snapshot_mc13'+year+category+'.root')
   mytree.Add(directory+'snapshot_mc14'+year+category+'.root')
   mytree.Add(directory+'snapshot_mc15'+year+category+'.root')      

   mytree.Add(directory+'snapshot_mc100'+year+category+'.root') # DY
   mytree.Add(directory+'snapshot_mc101'+year+category+'.root') # EWK 
   mytree.Add(directory+'snapshot_mc102'+year+category+'.root') # TT12L
   
   mytree.Add(directory+'snapshot_mc201'+year+category+'.root') # VV
   mytree.Add(directory+'snapshot_mc202'+year+category+'.root') # VV
   mytree.Add(directory+'snapshot_mc203'+year+category+'.root') # VV
   mytree.Add(directory+'snapshot_mc204'+year+category+'.root') # VV   
   
   return mytree
