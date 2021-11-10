import ihm.dictionary
import urllib.request

d = ihm.dictionary.Dictionary()
for dn in 'mmcif_pdbx_v50.dic', 'mmcif_ihm.dic':
     url = 'http://mmcif.wwpdb.org/dictionaries/ascii/' + dn
     with urllib.request.urlopen(url) as fh:
         d += ihm.dictionary.read(fh)

with open('threading_NPC_inner_ring.cif') as fh:
     d.validate(fh)
