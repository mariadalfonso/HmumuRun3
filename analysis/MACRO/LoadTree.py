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
      mytree.Add(directory+'snapshot_mc-23'+year+category+'.root')
      mytree.Add(directory+'snapshot_mc-24'+year+category+'.root')

   if (year == '_22023'):
      mytree.Add(directory+'snapshot_mc-31'+year+category+'.root')
      mytree.Add(directory+'snapshot_mc-32'+year+category+'.root')

   if (year == '_2024'):
      mytree.Add(directory+'snapshot_mc-41'+year+category+'.root')
      mytree.Add(directory+'snapshot_mc-42'+year+category+'.root')
      mytree.Add(directory+'snapshot_mc-43'+year+category+'.root')
      mytree.Add(directory+'snapshot_mc-44'+year+category+'.root')
      mytree.Add(directory+'snapshot_mc-45'+year+category+'.root')
      mytree.Add(directory+'snapshot_mc-46'+year+category+'.root')
      mytree.Add(directory+'snapshot_mc-47'+year+category+'.root')
      mytree.Add(directory+'snapshot_mc-48'+year+category+'.root')
      mytree.Add(directory+'snapshot_mc-49'+year+category+'.root')
      mytree.Add(directory+'snapshot_mc-50'+year+category+'.root')
      mytree.Add(directory+'snapshot_mc-51'+year+category+'.root')
      mytree.Add(directory+'snapshot_mc-52'+year+category+'.root')
      mytree.Add(directory+'snapshot_mc-53'+year+category+'.root')
      mytree.Add(directory+'snapshot_mc-54'+year+category+'.root')

   #Hmumu
   mytree.Add(directory+'snapshot_mc10'+year+category+'.root')
   mytree.Add(directory+'snapshot_mc11'+year+category+'.root')
   mytree.Add(directory+'snapshot_mc12'+year+category+'.root')
   mytree.Add(directory+'snapshot_mc13'+year+category+'.root')
   mytree.Add(directory+'snapshot_mc14'+year+category+'.root')
   mytree.Add(directory+'snapshot_mc15'+year+category+'.root')

   #H Zgamma
   mytree.Add(directory+'snapshot_mc20'+year+category+'.root')
   mytree.Add(directory+'snapshot_mc21'+year+category+'.root')
   mytree.Add(directory+'snapshot_mc22'+year+category+'.root')
   mytree.Add(directory+'snapshot_mc23'+year+category+'.root')
   mytree.Add(directory+'snapshot_mc24'+year+category+'.root')

   if (year == '_12022' or year == '_22022' or year == '_12023' or year == '_22023'):
      mytree.Add(directory+'snapshot_mc25'+year+category+'.root')
   elif year == '_2024':
      mytree.Add(directory+'snapshot_mc26'+year+category+'.root')

#      mytree.Add(directory+'snapshot_mc101'+year+category+'.root') # EWK

   if (year == '_12022' or year == '_22022' or year == '_12023' or year == '_22023'):
      mytree.Add(directory+'snapshot_mc100'+year+category+'.root') # DY
   if year == '_2024':
      mytree.Add(directory+'snapshot_mc103'+year+category+'.root') # DY
      mytree.Add(directory+'snapshot_mc104'+year+category+'.root') # DY

   mytree.Add(directory+'snapshot_mc102'+year+category+'.root') # TT12L

   mytree.Add(directory+'snapshot_mc201'+year+category+'.root') # VV
   mytree.Add(directory+'snapshot_mc202'+year+category+'.root') # VV
   mytree.Add(directory+'snapshot_mc203'+year+category+'.root') # VV
   mytree.Add(directory+'snapshot_mc204'+year+category+'.root') # VV
   mytree.Add(directory+'snapshot_mc205'+year+category+'.root') # VV

   mytree.Add(directory+'snapshot_mc211'+year+category+'.root') # VV
   mytree.Add(directory+'snapshot_mc212'+year+category+'.root') # VV
   mytree.Add(directory+'snapshot_mc213'+year+category+'.root') # VV
   mytree.Add(directory+'snapshot_mc214'+year+category+'.root') # VV

   return mytree
