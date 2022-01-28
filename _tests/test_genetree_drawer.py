#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pytest
import os.path as op
import tempfile
import shutil
from collections import namedtuple


from genetree_drawer import *

logger.setLevel(logging.DEBUG)

FileForTesting = namedtuple('TestFile', ['name', 'content'])
TEMPDIRNAME = tempfile.mkdtemp(prefix='test_genetree_drawer')

def newtestfile(basename, content):
    return FileForTesting(op.join(TEMPDIRNAME, basename), content)


def test_can_use_gui():
    import matplotlib as mpl
    try:
        mpl.use('TkAgg')
    except ImportError:
        mpl.use('Qt5Agg')

### Dependency modules to test first:
# - dendro.reconciled

### Unit features
class Test_read_params_fromfile:
    """When the --fromfile option is used."""
    pass

### Simple data

SIMPLE_PHYLTREE = newtestfile("simplePhylTree.nwk",
"""
(
  Ggo:9.0631,
  (
	Hsa:6.6509,
	(
	  Ppa:2.82006,
	  Ptr:2.82006
	)Pan:3.83085
  )HomoPan:2.41219
)Homininae:6.69907;
""")

SIMPLE_TREEBEST = newtestfile(
        "treebest_onedup.nwk",
        "((Ppa_G0:1[&&NHX:S=Ppa],Ppa_G1:1[&&NHX:S=Ppa]):1[&&NHX:S=Ppa:D=2],Ptr_G0:1[&&NHX:S=Ptr]):1[&&NHX:S=Pan:D=0];")
# Will fail because Ptr is not an expected child of Ppa.
SIMPLE_TREEBEST_T = newtestfile(
        "treebest_transfer.nwk",
        "((Ppa_G0[&&NHX:S=Ppa],Ptr_G1[&&NHX:S=Ptr]):1[&&NHX:S=Ppa:D=11],Ptr_G0:1[&&NHX:S=Ptr]):1[&&NHX:S=Pan:D=0];")
SIMPLE_TREEBEST_T2 = newtestfile(
        "treebest_transfer2.nwk",
        "((Ppa_G0:1[&&NHX:S=Ppa],Ptr_G1:1[&&NHX:S=Ptr:T=-1]):1[&&NHX:S=Ppa:T=1],Ptr_G0:1[&&NHX:S=Ptr])[&&NHX:S=Pan:D=0];")

### Complex data to test the entire program result

# Example data from:
# - ~/ws7/lateralite/
# - ~/ws7/DUPLI_data93/alignments_analysis/duprates/
# - ~/ws7/DUPLI_data93/alignments/ENSGT*

# Example command (NODAL Amniota tree based on Ensembl 93 with edition cutoff=0.35):
# cd ~/ws7/DUPLI_data93/alignments/ENSGT00910000143982/subtreesGoodQualO2/
# genetree_drawer.py -a3 -t -e 93 -p ../../PhylTree.TimeTree201901.Ensembl93-like.goodQual.nwk -l AmniotaENSGT00910000143982.A.b -i Amniota,Theria,Euarchontoglires,Neopterygii,Carnivora,Rodentia,Primates,Neognathae -- ~/ws7/lateralite/data93/NODAL_035GoodQual.pdf:16x13 subtrees035GoodQual/EuteleostomiENSGT00910000143982.A.b.nwk
#
# genetree_drawer.py -c 'Artiofabula,Sauria,Myotis lucifugus,Zooamata' -i 'Neopterygii,Neognathae,Amniota,Theria,Euarchontoglires,Primates,Carnivora' -a3 -e93 -p ~/ws7/DUPLI_data93/PhylTree.TimeTree201901.Ensembl93-like.goodQual.nwk -l AmniotaENSGT00390000008796.b -- data93/LMLN2_035GoodQual.pdf:16x16 ~/ws7/DUPLI_data93/alignments/ENSGT00390000008796/subtrees035GoodQual/EuteleostomiENSGT00390000008796.b.nwk,LMLN2

#### DATA

# PhylTree.TimeTree201901.Ensembl93-like.goodQual.nwk
PHYLTREE_NEWICK = newtestfile(
    'PhylTree.nwk',
    """(
  Saccharomyces cerevisiae|Yeast|saccharomyces_cerevisiae:1105.06,
  (
    (
      Drosophila melanogaster|Fruitfly|drosophila_melanogaster:743,
      Caenorhabditis elegans|caenorhabditis_elegans|C_elegans:743
    )Ecdysozoa:53.5556,
    (
      (
        (
          Lepisosteus oculatus|Spotted gar|lepisosteus_oculatus:314.684,
          (
            (
              Danio rerio|danio_rerio|Zebrafish:152.6,
              Astyanax mexicanus|astyanax_mexicanus|Cave fish:152.6
            )Otophysi:77.303,
            (
              (
                Oreochromis niloticus|Tilapia|oreochromis_niloticus:119,
                (
                  Oryzias latipes|Medaka|oryzias_latipes:93.1723,
                  (
                    Poecilia formosa|poecilia_formosa|Amazon molly:57,
                    Xiphophorus maculatus|Platyfish|xiphophorus_macula:57
                  )Poeciliinae:36.1723
                )Atherinomorphae:25.8277
              )Ovalentaria:9,
              (
                Gasterosteus aculeatus|gasterosteus_aculeatus|Stickleback:110,
                (
                  Tetraodon nigroviridis|tetraodon_nigroviridis|Tetraodon:51.6,
                  Takifugu rubripes|Fugu|takifugu_rubripes:51.6
                )Tetraodontidae:58.4
              )Eupercaria:18
            )Percomorphaceae:101.903
          )Clupeocephala:84.7817
        )Neopterygii:120.64,
        (
          (
            Monodelphis domestica|monodelphis_domestica|Opossum:158.598,
            (
              Loxodonta africana|Elephant|loxodonta_africana:105.46,
              (
                (
                  Myotis lucifugus|Microbat|myotis_lucifugus:78.5287,
                  (
                    (
                      Sus scrofa|sus_scrofa|Pig:61.966,
                      (
                        Bos taurus|Cow|bos_taurus:24.6,
                        (
                          Capra hircus|capra_hircus|Goat:9.75,
                          Ovis aries|Sheep|ovis_aries:9.75
                        )Caprinae:14.85
                      )Bovidae:37.366
                    )Artiofabula:15.7889,
                    (
                      (
                        (
                          Felis catus|felis_catus|Cat:15.175,
                          (
                            Panthera tigris altaica|Tiger|panthera_tigris_altaica:7.41234,
                            Panthera pardus|panthera_pardus|leopard:7.41234
                          )Panthera:7.76265
                        )Felidae:39.1464,
                        (
                          Canis lupus familiaris|canis_familiaris|Dog:45.5283,
                          Mustela putorius furo|Ferret|Mustela_putorius_furo:45.5283
                        )Caniformia:8.79309
                      )Carnivora:22.9052,
                      Equus caballus|equus_caballus|Horse:77.2266
                    )Zooamata:0.528362
                  )Fereuungulata:0.773786
                )Scrotifera:17.9337,
                (
                  (
                    Oryctolagus cuniculus|oryctolagus_cuniculus|Rabbit:82.1408,
                    (
                      Ictidomys tridecemlineatus|Squirrel|ictidomys_tridecemlineatus:72.8767,
                      (
                        Fukomys damarensis|fukomys_damarensis|Damara mole rat:43.3585,
                        (
                          (
                            Cavia aperea|Brazilian-Guinea pig|cavia_aperea:5.65679,
                            Cavia porcellus|Guinea pig|cavia_porcellus:5.65679
                          )Cavia:30.1222,
                          (
                            Octodon degus|Degu|octodon_degus:32.8554,
                            Chinchilla lanigera|Long-tailed chinchilla|chinchilla_lanigera:32.8554
                          )ChinchillaOctodon:2.92363
                        )Caviomorpha:7.57947
                      )Hystricomorpha:29.5182,
                      (
                        Dipodomys ordii|dipodomys_ordii|Kangaroo rat:69.8984,
                        (
                          Jaculus jaculus|jaculus_jaculus|Lesser Egyptian jerboa:54.8024,
                          (
                            Nannospalax galili|Upper Galilee mountains blind mole rat|nannospalax_galili:45.2692,
                            (
                              (
                                Peromyscus maniculatus bairdii|peromyscus_maniculatus_bairdii|Northern American deer mouse:28.8,
                                Microtus ochrogaster|microtus_ochrogaster|Prairie vole:28.8,
                                Mesocricetus auratus|mesocricetus_auratus|Golden Hamster:28.8
                              )Cricetidae:3.86363,
                              (
                                Rattus norvegicus|Rat|rattus_norvegicus:20.8874,
                                (
                                  Mus pahari|mus_pahari|Shrew mouse:8.29083,
                                  (
                                    Mus caroli|mus_caroli|Ryukyu mouse:7.4126,
                                    (
                                      Mus musculus|mus_musculus|Mouse:3.06548,
                                      Mus spretus|spretus|mus_spretus:3.06548
                                    )Mus_B:4.34712
                                  )Mus_A:0.878228
                                )Mus:12.5966
                              )Murinae:11.7762
                            )Eumuroida:12.6056
                          )Muroidea:9.53321
                        )Myomorpha:15.096
                      )MyodontaCastorimorpha:2.97833
                    )Rodentia:9.26409
                  )Glires:7.68239,
                  (
                    (
                      Otolemur garnettii|Bushbaby|otolemur_garnettii:59.3244,
                      (
                        Propithecus coquereli|propithecus_coquereli|Coquerels sifaka:37.7826,
                        Microcebus murinus|Mouse lemur|microcebus_murinus:37.7826
                      )Lemuriformes:21.5417
                    )Strepsirrhini:14.5125,
                    (
                      (
                        Aotus nancymaae|aotus_nancymaae|Mas night monkey:19.6811,
                        Callithrix jacchus|callithrix_jacchus|Marmoset:19.6811,
                        (
                          Saimiri boliviensis boliviensis|saimiri_boliviensis|Bolivian squirrel monkey:16.0705,
                          Cebus capucinus imitator|cebus_capucinus|Capuchin:16.0705
                        )Cebidae:3.61059
                      )Platyrrhini:23.4702,
                      (
                        (
                          Nomascus leucogenys|Gibbon|nomascus_leucogenys:20.1892,
                          (
                            Pongo abelii|pongo_abelii|Orangutan:15.7622,
                            (
                              Gorilla gorilla gorilla|gorilla_gorilla|Gorilla:9.0631,
                              (
                                Homo sapiens|Human|homo_sapiens:6.6509,
                                (
                                  Pan paniscus|Bonobo|pygmy chimpanzee|pan_paniscus:2.82006,
                                  Pan troglodytes|Chimpanzee|pan_troglodytes:2.82006
                                )Pan:3.83085
                              )HomoPan:2.41219
                            )Homininae:6.69907
                          )Hominidae:4.42705
                        )Hominoidea:9.25233,
                        (
                          Colobus angolensis palliatus|Angola colobus|colobus_palliatus:19.4223,
                          (
                            (
                              Cercocebus atys|Sooty mangabey|cercocebus_atys:12.4,
                              Mandrillus leucophaeus|Drill|mandrillus_leucophaeus:12.4,
                              Papio anubis|Olive baboon|papio_anubis:12.4,
                              (
                                Macaca fascicularis|macaca_fascicularis|Crab-eating macaque:5.27817,
                                Macaca mulatta|macaca_mulatta|Rhesus:5.27817,
                                Macaca nemestrina|Pig-tailed macaque|macaca_nemestrina:5.27816
                              )Macaca:7.12184
                            )Papionini:1.34957,
                            Chlorocebus sabaeus|Vervet Monkey|chlorocebus_sabaeus:13.7496
                          )Cercopithecinae:5.67273
                        )Cercopithecidae:10.0192
                      )Catarrhini:13.7097
                    )Simiiformes:30.6856
                  )Primates:15.9863
                )Euarchontoglires:6.6392
              )Boreoeutheria:8.99739
            )Eutheria:53.1378
          )Theria:153.306,
          (
            Anolis carolinensis|anolis_carolinensis|Lizard:279.657,
            (
              (
                Ficedula albicollis|ficedula_albicollis|Flycatcher:43.7,
                Taeniopygia guttata|taeniopygia_guttata|Zebra finch:43.7
              )Passeriformes:54.3429,
              (
                Meleagris gallopavo|Turkey|meleagris_gallopavo:37.2,
                Gallus gallus|Chicken|gallus_gallus:37.2
              )Phasianidae:60.8429
            )Neognathae:181.614
          )Sauria:32.2469
        )Amniota:123.42
      )Euteleostomi:241.124,
      (
        Ciona savignyi|Ciona_sav|ciona_savignyi:100,
        Ciona intestinalis|Ciona_int|ciona_intestinalis:100
      )Ciona:576.448
    )Chordata:120.108
  )Bilateria:308.504
)Opisthokonta;
"""
)

# ENSGT00910000143982/subtrees035GoodQual/EuteleostomiENSGT00910000143982.A.b.nwk
NODAL_NEWICK = newtestfile(
    'NODAL.nwk',
    (
"((((((((ENSTRUG00000012437:0.088844[&&NHX:ID=8879:S=Takifugu.rubripes:D=N],"
"ENSTNIG00000015847:0.126717[&&NHX:ID=8528:S=Tetraodon.nigroviridis:D=N])"
"TetraodontidaeENSGT00910000143982.A.b.c.a:0.121732[&&NHX:ID=8092999:B=99:S=Tetraodontidae:D=N],"
"ENSGACG00000008499:0.10048[&&NHX:ID=8875:S=Gasterosteus.aculeatus:D=N])"
"EupercariaENSGT00910000143982.A.b.c.a:0.000528[&&NHX:ID=8092994:B=0:S=Eupercaria:D=N],"
"(((ENSPFOG00000008931:0.012784[&&NHX:ID=9181:S=Poecilia.formosa:D=N],"
"ENSXMAG00000017893:0.040258[&&NHX:ID=8217:S=Xiphophorus.maculatus:D=N])"
"PoeciliinaeENSGT00910000143982.A.b.c.a:0.076586[&&NHX:ID=8092989:B=100:S=Poeciliinae:D=N],"
"ENSORLG00000006553:0.136817[&&NHX:ID=9152:S=Oryzias.latipes:D=N])"
"AtherinomorphaeENSGT00910000143982.A.b.c.a:0.033357[&&NHX:ID=8092984:B=27:S=Atherinomorphae:D=N],"
"ENSONIG00000013352:0.066249[&&NHX:ID=8268:S=Oreochromis.niloticus:D=N])"
"OvalentariaENSGT00910000143982.A.b.c.a:0.042389[&&NHX:ID=8092980:B=8:S=Ovalentaria:D=N])"
"PercomorphaceaeENSGT00910000143982.A.b.c.a:0.20458[&&NHX:ID=8092977:B=0:S=Percomorphaceae:D=N],"
"(ENSAMXG00000021167:0.145127[&&NHX:ID=8326:S=Astyanax.mexicanus:D=N],"
"ENSDARG00000014309:0.276025[&&NHX:ID=9055:S=Danio.rerio:D=N])"
"OtophysiENSGT00910000143982.A.b.c.a:0.105416[&&NHX:ID=8092964:B=78:S=Otophysi:D=N])"
"ClupeocephalaENSGT00910000143982.A.b.c.a:0.034524[&&NHX:ID=8092961:B=17:S=Clupeocephala:D=N],"
"((ENSAMXG00000008977:0.112421[&&NHX:ID=8743:S=Astyanax.mexicanus:D=N],"
"ENSDARG00000057096:0.270078[&&NHX:ID=8971:S=Danio.rerio:D=N])"
"OtophysiENSGT00910000143982.A.b.c.b:0.023617[&&NHX:ID=8092919:B=61:S=Otophysi:D=N],"
"(((ENSTRUG00000010779:0.055247[&&NHX:ID=8316:S=Takifugu.rubripes:D=N],"
"ENSTNIG00000013237:0.057252[&&NHX:ID=8258:S=Tetraodon.nigroviridis:D=N])"
"TetraodontidaeENSGT00910000143982.A.b.c.b:0.072501[&&NHX:ID=8092956:B=100:S=Tetraodontidae:D=N],"
"ENSGACG00000017712:0.15561[&&NHX:ID=8288:S=Gasterosteus.aculeatus:D=N])"
"EupercariaENSGT00910000143982.A.b.c.b:0.016135[&&NHX:ID=8092952:B=85:S=Eupercaria:D=N],"
"(((ENSXMAG00000002355:0.037475[&&NHX:ID=8128:S=Xiphophorus.maculatus:D=N],"
"ENSPFOG00000006036:0.044127[&&NHX:ID=8842:S=Poecilia.formosa:D=N])"
"PoeciliinaeENSGT00910000143982.A.b.c.b:0.140594[&&NHX:ID=8092944:B=100:S=Poeciliinae:D=N],"
"ENSORLG00000011275:0.232772[&&NHX:ID=9067:S=Oryzias.latipes:D=N])"
"AtherinomorphaeENSGT00910000143982.A.b.c.b:0.011004[&&NHX:ID=8092937:B=6:S=Atherinomorphae:D=N],"
"ENSONIG00000004012:0.144392[&&NHX:ID=8546:S=Oreochromis.niloticus:D=N])"
"OvalentariaENSGT00910000143982.A.b.c.b:0.056251[&&NHX:ID=8092932:B=4:S=Ovalentaria:D=N])"
"PercomorphaceaeENSGT00910000143982.A.b.c.b:0.131936[&&NHX:ID=8092928:B=3:S=Percomorphaceae:D=N])"
"ClupeocephalaENSGT00910000143982.A.b.c.b:0.104302[&&NHX:ID=8092916:B=5:S=Clupeocephala:D=N])"
"ClupeocephalaENSGT00910000143982.A.b.c:0.076526[&&NHX:ID=8092910:B=7:S=Clupeocephala:SIS=1.0:D=Y],"
"(((((ENSXMAG00000015109:0.015042[&&NHX:ID=8236:S=Xiphophorus.maculatus:D=N],"
"ENSPFOG00000005318:0.023904[&&NHX:ID=8983:S=Poecilia.formosa:D=N])"
"PoeciliinaeENSGT00910000143982.A.b.d:0.132076[&&NHX:ID=8092510:B=100:S=Poeciliinae:D=N],"
"ENSORLG00000009098:0.168817[&&NHX:ID=9068:S=Oryzias.latipes:D=N])"
"AtherinomorphaeENSGT00910000143982.A.b.d:0.038942[&&NHX:ID=8092504:B=45:S=Atherinomorphae:D=N],"
"ENSONIG00000016114:0.157221[&&NHX:ID=8201:S=Oreochromis.niloticus:D=N])"
"OvalentariaENSGT00910000143982.A.b.d:0.025644[&&NHX:ID=8092499:B=22:S=Ovalentaria:D=N],"
"((ENSTRUG00000012942:0.051852[&&NHX:ID=8935:S=Takifugu.rubripes:D=N],"
"ENSTNIG00000005578:0.061915[&&NHX:ID=8017:S=Tetraodon.nigroviridis:D=N])"
"TetraodontidaeENSGT00910000143982.A.b.d:0.072113[&&NHX:ID=8092493:B=100:S=Tetraodontidae:D=N],"
"ENSGACG00000002333:0.113672[&&NHX:ID=9080:S=Gasterosteus.aculeatus:D=N])"
"EupercariaENSGT00910000143982.A.b.d:0.025686[&&NHX:ID=8092490:B=94:S=Eupercaria:D=N])"
"PercomorphaceaeENSGT00910000143982.A.b.d:0.285711[&&NHX:ID=8092487:B=35:S=Percomorphaceae:D=N],"
"(ENSAMXG00000015719:0.168314[&&NHX:ID=9091:S=Astyanax.mexicanus:D=N],"
"ENSDARG00000101279:0.341961[&&NHX:ID=9163:S=Danio.rerio:D=N])"
"OtophysiENSGT00910000143982.A.b.d:0.154284[&&NHX:ID=8092516:B=94:S=Otophysi:D=N])"
"ClupeocephalaENSGT00910000143982.A.b.d:0.154284[&&NHX:ID=100003014:S=Clupeocephala:_r=1:D=N])"
"ClupeocephalaENSGT00910000143982.A.b:0.076526[&&NHX:ID=100003013:B=92:_f=_8092479_:S=Clupeocephala:_r=1:D=Y],"
"(((ENSLOCG00000016067:0.005352[&&NHX:ID=8958:S=Lepisosteus.oculatus:D=N],"
"ENSLOCG00000016069:0.008034[&&NHX:ID=1609847:S=Lepisosteus.oculatus:D=N])"
"Lepisosteus.oculatusENSGT00910000143982.A.b.e.a:0.010815[&&NHX:ID=8092901:B=3:S=Lepisosteus.oculatus:SIS=1.0:D=Y],"
"ENSLOCG00000016070:0.004197[&&NHX:ID=8053:S=Lepisosteus.oculatus:D=N])"
"Lepisosteus.oculatusENSGT00910000143982.A.b.e:0.160233[&&NHX:ID=8092898:B=3:S=Lepisosteus.oculatus:SIS=1.0:D=Y],"
"ENSLOCG00000005102:0.386068[&&NHX:ID=8103:S=Lepisosteus.oculatus:D=N])"
"Lepisosteus.oculatusENSGT00910000143982.A.b:0.160233[&&NHX:ID=100003012:S=Lepisosteus.oculatus:_r=1:D=Y])"
"NeopterygiiENSGT00910000143982.A.b:0.131748[&&NHX:ID=100003011:B=64:_f=_8092895_.8092476_:S=Neopterygii:_r=1:D=N],"
"((ENSACAG00000008399:0.577071[&&NHX:ID=8620:S=Anolis.carolinensis:D=N],"
"((ENSTGUG00000004739:0.008434[&&NHX:ID=8479:S=Taeniopygia.guttata:D=N],"
"ENSFALG00000006750:0.037295[&&NHX:ID=8102:S=Ficedula.albicollis:D=N])"
"PasseriformesENSGT00910000143982.A.b:0.150789[&&NHX:ID=8093061:B=99:S=Passeriformes:D=N],"
"(ENSGALG00000003209:0.092701[&&NHX:ID=8055:S=Gallus.gallus:D=N],"
"(ENSMGAG00000016191:0.248647[&&NHX:ID=1619661:S=Meleagris.gallopavo:D=N],"
"ENSMGAG00000002207:0.04777[&&NHX:ID=9034:S=Meleagris.gallopavo:D=N])"
"Meleagris.gallopavoENSGT00910000143982.A.b:0.04777[&&NHX:ID=100003025:S=Meleagris.gallopavo:_r=1:D=Y])"
"PhasianidaeENSGT00910000143982.A.b:0.118956[&&NHX:ID=100003024:B=100:_f=_8093056_:S=Phasianidae:_r=1:D=N])"
"NeognathaeENSGT00910000143982.A.b:0.127569[&&NHX:ID=100003022:B=73:_f=_8093043_:S=Neognathae:_r=1:D=N])"
"SauriaENSGT00910000143982.A.b:0.097792[&&NHX:ID=100003019:S=Sauria:_r=1:D=N],"
"(((((((ENSPCOG00000016462:0.019737[&&NHX:ID=8899:S=Propithecus.coquereli:D=N],"
"ENSMICG00000015080:0.04425[&&NHX:ID=9119:S=Microcebus.murinus:D=N])"
"LemuriformesENSGT00910000143982.A.b:0.025094[&&NHX:ID=8092613:B=91:S=Lemuriformes:D=N],"
"ENSOGAG00000005716:0.065971[&&NHX:ID=8926:S=Otolemur.garnettii:D=N])"
"StrepsirrhiniENSGT00910000143982.A.b:0.02211[&&NHX:ID=8092608:B=22:S=Strepsirrhini:D=N],"
"(((((ENSGGOG00000002581:0.00424[&&NHX:ID=8749:S=Gorilla.gorilla.gorilla:D=N],"
"((ENSPPAG00000032411:0[&&NHX:ID=8365:S=Pan.paniscus:D=N],"
"ENSPTRG00000002592:0.004237[&&NHX:ID=1871691:S=Pan.troglodytes:D=N])"
"PanENSGT00910000143982.A.b:0.001058[&&NHX:ID=8092706:B=86:S=Pan:D=N],"
"ENSG00000156574:0.007445[&&NHX:ID=9159:S=Homo.sapiens:D=N])"
"HomoPanENSGT00910000143982.A.b:0.001058[&&NHX:ID=100003034:S=HomoPan:_r=1:D=N])"
"HomininaeENSGT00910000143982.A.b:0.01213[&&NHX:ID=8092698:B=80:_f=_8092702_:S=Homininae:D=N],"
"ENSPPYG00000002370:0.007086[&&NHX:ID=8892:S=Pongo.abelii:D=N])"
"HominidaeENSGT00910000143982.A.b:0[&&NHX:ID=8092690:B=2:S=Hominidae:D=N],"
"ENSNLEG00000016126:0.004949[&&NHX:ID=8897:S=Nomascus.leucogenys:D=N])"
"HominoideaENSGT00910000143982.A.b:0.006869[&&NHX:ID=8092687:B=2:S=Hominoidea:D=N],"
"((ENSCSAG00000008299:0.010628[&&NHX:ID=8041:S=Chlorocebus.sabaeus:D=N],"
"((ENSMMUG00000023171:0.002133[&&NHX:ID=1716192:S=Macaca.mulatta:D=N],"
"ENSMFAG00000041504:0.006382[&&NHX:ID=8341:S=Macaca.fascicularis:D=N],"
"ENSMNEG00000036384:0.004252[&&NHX:ID=1716193:S=Macaca.nemestrina:D=N])"
"MacacaENSGT00910000143982.A.b:0.001056[&&NHX:ID=8092676:B=64:_f=_8092681_:S=Macaca:D=N],"
"ENSPANG00000022260:0.00324[&&NHX:ID=8812:S=Papio.anubis:D=N],"
"ENSCATG00000037917:0.001061[&&NHX:ID=8742:S=Cercocebus.atys:D=N],"
"ENSMLEG00000031653:0.046494[&&NHX:ID=8412:S=Mandrillus.leucophaeus:D=N])"
"PapioniniENSGT00910000143982.A.b:0.001056[&&NHX:ID=100003035:S=Papionini:_r=1:D=N])"
"CercopithecinaeENSGT00910000143982.A.b:0[&&NHX:ID=8092657:B=0:_f=_8092668_.8092665_.8092660_:S=Cercopithecinae:D=N],"
"ENSCANG00000040789:0.008949[&&NHX:ID=8886:S=Colobus.angolensis.palliatus:D=N])"
"CercopithecidaeENSGT00910000143982.A.b:0.008293[&&NHX:ID=8092644:B=6:S=Cercopithecidae:D=N])"
"CatarrhiniENSGT00910000143982.A.b:0.00764[&&NHX:ID=8092638:B=53:S=Catarrhini:D=N],"
"(ENSCJAG00000016288:0.016444[&&NHX:ID=8930:S=Callithrix.jacchus:D=N],"
"ENSANAG00000023960:0.009755[&&NHX:ID=8639:S=Aotus.nancymaae:D=N],"
"(ENSSBOG00000035806:0.015751[&&NHX:ID=9171:S=Saimiri.boliviensis.boliviensis:D=N],"
"ENSCCAG00000032672:0.007007[&&NHX:ID=8312:S=Cebus.capucinus.imitator:D=N])"
"CebidaeENSGT00910000143982.A.b:0.007007[&&NHX:ID=100003036:S=Cebidae:_r=1:D=N])"
"PlatyrrhiniENSGT00910000143982.A.b:0.014124[&&NHX:ID=8092622:B=67:_f=_8092631_.8092627_:S=Platyrrhini:D=N])"
"SimiiformesENSGT00910000143982.A.b:0.040431[&&NHX:ID=8092619:B=78:S=Simiiformes:D=N])"
"PrimatesENSGT00910000143982.A.b:0.012485[&&NHX:ID=100003032:B=0:_f=_8092606_:S=Primates:_r=1:D=N],"
"(((ENSFDAG00000014759:0.065764[&&NHX:ID=8686:S=Fukomys.damarensis:D=N],"
"((ENSCLAG00000016370:0.013697[&&NHX:ID=8016:S=Chinchilla.lanigera:D=N],"
"ENSODEG00000013687:0.038807[&&NHX:ID=8509:S=Octodon.degus:D=N])"
"ChinchillaOctodonENSGT00910000143982.A.b:0.00758[&&NHX:ID=8092745:B=55:S=ChinchillaOctodon:D=N],"
"(ENSCPOG00000025772:0[&&NHX:ID=8159:S=Cavia.porcellus:D=N],"
"ENSCAPG00000009829:0.053007[&&NHX:ID=8824:S=Cavia.aperea:D=N])"
"CaviaENSGT00910000143982.A.b:0.048738[&&NHX:ID=8092752:B=100:S=Cavia:D=N])"
"CaviomorphaENSGT00910000143982.A.b:0.00758[&&NHX:ID=100003038:S=Caviomorpha:_r=1:D=N])"
"HystricomorphaENSGT00910000143982.A.b:0.039553[&&NHX:ID=8092727:B=0:_f=_8092730_:S=Hystricomorpha:D=N],"
"ENSSTOG00000012946:0.064069[&&NHX:ID=8855:S=Ictidomys.tridecemlineatus:D=N],"
"(((((ENSMAUG00000012281:0.061062[&&NHX:ID=8080:S=Mesocricetus.auratus:D=N],"
"ENSPEMG00000020412:0.050469[&&NHX:ID=9044:S=Peromyscus.maniculatus.bairdii:D=N],"
"ENSMOCG00000002962:0.047725[&&NHX:ID=8551:S=Microtus.ochrogaster:D=N])"
"CricetidaeENSGT00910000143982.A.b:0.008259[&&NHX:ID=8092773:B=81:_f=_8092776_:S=Cricetidae:D=N],"
"((((ENSMUSG00000037171:0.002118[&&NHX:ID=7998:S=Mus.musculus:D=N],"
"MGP_SPRETEiJ_G0016176:0.002121[&&NHX:ID=1895151:S=Mus.spretus:D=N])"
"Mus_BENSGT00910000143982.A.b:0.004289[&&NHX:ID=8092808:B=98:S=Mus_B:D=N],"
"MGP_CAROLIEiJ_G0015366:0.010784[&&NHX:ID=8997:S=Mus.caroli:D=N])"
"Mus_AENSGT00910000143982.A.b:0.004174[&&NHX:ID=8092803:B=97:S=Mus_A:D=N],"
"MGP_PahariEiJ_G0030782:0.019662[&&NHX:ID=8755:S=Mus.pahari:D=N])"
"MusENSGT00910000143982.A.b:0.006865[&&NHX:ID=8092796:B=86:S=Mus:D=N],"
"ENSRNOG00000000556:0.022028[&&NHX:ID=8266:S=Rattus.norvegicus:D=N])"
"MurinaeENSGT00910000143982.A.b:0.043659[&&NHX:ID=8092793:B=87:S=Murinae:D=N])"
"EumuroidaENSGT00910000143982.A.b:0.045774[&&NHX:ID=8092769:B=91:S=Eumuroida:D=N],"
"ENSNGAG00000007985:0.055041[&&NHX:ID=8184:S=Nannospalax.galili:D=N])"
"MuroideaENSGT00910000143982.A.b:0.006881[&&NHX:ID=8092763:B=37:S=Muroidea:D=N],"
"ENSJJAG00000019905:0.128339[&&NHX:ID=8227:S=Jaculus.jaculus:D=N])"
"MyomorphaENSGT00910000143982.A.b:0.016972[&&NHX:ID=8092757:B=22:S=Myomorpha:D=N],"
"ENSDORG00000011307:0.067933[&&NHX:ID=9023:S=Dipodomys.ordii:D=N])"
"MyodontaCastorimorphaENSGT00910000143982.A.b:0.016972[&&NHX:ID=100003037:S=MyodontaCastorimorpha:_r=1:D=N])"
"RodentiaENSGT00910000143982.A.b:0.020971[&&NHX:ID=8092716:B=1:_f=_8092724_.8092719_:S=Rodentia:D=N],"
"ENSOCUG00000008685:0.053052[&&NHX:ID=8460:S=Oryctolagus.cuniculus:D=N])"
"GliresENSGT00910000143982.A.b:0.007468[&&NHX:ID=100003031:S=Glires:_r=1:D=N])"
"EuarchontogliresENSGT00910000143982.A.b:0.014937[&&NHX:ID=100003030:B=0:_f=_8092604_:S=Euarchontoglires:_r=1:D=N],"
"(ENSMLUG00000015297:0.056412[&&NHX:ID=8568:S=Myotis.lucifugus:D=N],"
"((ENSSSCG00000039399:0.040227[&&NHX:ID=8583:S=Sus.scrofa:D=N],"
"((ENSOARG00000006065:0.002254[&&NHX:ID=8767:S=Ovis.aries:D=N],"
"ENSCHIG00000000383:0.008416[&&NHX:ID=8180:S=Capra.hircus:D=N])"
"CaprinaeENSGT00910000143982.A.b:0.013956[&&NHX:ID=8092865:B=100:S=Caprinae:D=N],"
"ENSBTAG00000013090:0.021119[&&NHX:ID=9188:S=Bos.taurus:D=N])"
"BovidaeENSGT00910000143982.A.b:0.054003[&&NHX:ID=8092862:B=99:S=Bovidae:D=N])"
"ArtiofabulaENSGT00910000143982.A.b:0.022534[&&NHX:ID=100003044:S=Artiofabula:_r=1:D=N],"
"(ENSECAG00000017055:0.038551[&&NHX:ID=8353:S=Equus.caballus:D=N],"
"(((ENSPPRG00000014162:0[&&NHX:ID=8354:S=Panthera.pardus:D=N],"
"ENSPTIG00000014186:0[&&NHX:ID=8164:S=Panthera.tigris.altaica:D=N])"
"PantheraENSGT00910000143982.A.b:0.004166[&&NHX:ID=8092598:B=100:S=Panthera:D=N],"
"ENSFCAG00000001230:0.010901[&&NHX:ID=9108:S=Felis.catus:D=N])"
"FelidaeENSGT00910000143982.A.b:0.056277[&&NHX:ID=8092594:B=100:S=Felidae:D=N],"
"(ENSCAFG00000014052:0.044554[&&NHX:ID=8410:S=Canis.lupus.familiaris:D=N],"
"ENSMPUG00000004293:0.051655[&&NHX:ID=8378:S=Mustela.putorius.furo:D=N])"
"CaniformiaENSGT00910000143982.A.b:0.006687[&&NHX:ID=100003046:B=49:_f=_8092590_:S=Caniformia:_r=1:D=N])"
"CarnivoraENSGT00910000143982.A.b:0.004428[&&NHX:ID=100003045:B=22:_f=_8092587_:S=Carnivora:_r=1:D=N])"
"ZooamataENSGT00910000143982.A.b:0.003155[&&NHX:ID=100003042:B=0:_f=_8092873_:S=Zooamata:_r=1:D=N])"
"FereuungulataENSGT00910000143982.A.b:0.006311[&&NHX:ID=100003041:S=Fereuungulata:_r=1:D=N])"
"ScrotiferaENSGT00910000143982.A.b:0.025676[&&NHX:ID=100003040:B=0:_f=_8092843_.8092840_:S=Scrotifera:_r=1:D=N])"
"BoreoeutheriaENSGT00910000143982.A.b:0.024543[&&NHX:ID=100003027:B=0:_f=_8092570_.8092580_:S=Boreoeutheria:_r=1:D=N],"
"ENSLAFG00000021867:0.135493[&&NHX:ID=8905:S=Loxodonta.africana:D=N])"
"EutheriaENSGT00910000143982.A.b:0.067982[&&NHX:ID=8092552:B=0:_f=_8092818_.8092577_.8092558_.8092556_:S=Eutheria:SIS=0.0:D=N],"
"ENSMODG00000012158:0.281295[&&NHX:ID=8700:S=Monodelphis.domestica:D=N])"
"TheriaENSGT00910000143982.A.b:0.516509[&&NHX:ID=8092536:B=95:S=Theria:D=N])"
"AmniotaENSGT00910000143982.A.b:0.161705[&&NHX:ID=100003017:B=24:_f=_8093071_.8092530_:S=Amniota:_r=1:D=N])"
"EuteleostomiENSGT00910000143982.A.b:0.121325[&&NHX:ID=100003009:B=2:_f=_8092890_.8092473_.8092471_:S=Euteleostomi:SIS=0.1398:_r=1:D=N];"
)
)


# ~/ws7/DUPLI_data93/alignments/ENSGT00390000008796/subtrees035GoodQual/EuteleostomiENSGT00390000008796.b.nwk
LMLN2_NEWICK = newtestfile(
    'LMLN2.nwk',
    (
"(((ENSMODG00000004607:0.262975[&&NHX:D=N:ID=1048241:S=Monodelphis.domestica],"
"(((((((ENSNLEG00000031249:0.014263[&&NHX:D=N:ID=1048240:S=Nomascus.leucogenys],"
"((ENSGGOG00000043009:0.003211[&&NHX:D=N:ID=1048324:S=Gorilla.gorilla.gorilla],"
"((ENSPTRG00000044020:0.001068[&&NHX:D=N:ID=1865044:S=Pan.troglodytes],"
"ENSPPAG00000041297:0.001068[&&NHX:D=N:ID=1865045:S=Pan.paniscus])"
"PanENSGT00390000008796.b:0[&&NHX:D=N:B=98:ID=21037925:S=Pan],"
"ENSG00000283654:0[&&NHX:D=N:ID=1048253:S=Homo.sapiens])"
"HomoPanENSGT00390000008796.b:0[&&NHX:_r=1:D=N:ID=100371848:S=HomoPan])"
"HomininaeENSGT00390000008796.b:0.00598853[&&NHX:_f=_21037921_:D=N:B=83:ID=21037917:S=Homininae])"
"HominidaeENSGT00390000008796.b:0.00395747[&&NHX:reinserted=True:S=Hominidae])"
"HominoideaENSGT00390000008796.b:0.001308[&&NHX:D=N:B=99:ID=21037914:S=Hominoidea],"
"(ENSCANG00000032248:0.040587[&&NHX:D=N:ID=1048286:S=Colobus.angolensis.palliatus],"
"(((ENSMNEG00000039815:0.004285[&&NHX:D=N:ID=1048271:S=Macaca.nemestrina],"
"ENSMFAG00000039515:0.0161[&&NHX:D=N:ID=1048300:S=Macaca.fascicularis])"
"MacacaENSGT00390000008796.b:0[&&NHX:D=N:B=11:ID=21037952:S=Macaca],"
"ENSPANG00000033367:0.004408[&&NHX:D=N:ID=1048284:S=Papio.anubis],"
"ENSMLEG00000035762:0.006452[&&NHX:D=N:ID=1048226:S=Mandrillus.leucophaeus],"
"ENSCATG00000044967:0.005353[&&NHX:D=N:ID=1048238:S=Cercocebus.atys])"
"PapioniniENSGT00390000008796.b:0.00101628[&&NHX:_f=_21037947_.21037944_:D=N:B=13:ID=21037942:S=Papionini])"
"CercopithecinaeENSGT00390000008796.b:0.00427172[&&NHX:reinserted=True:S=Cercopithecinae])"
"CercopithecidaeENSGT00390000008796.b:0.009623[&&NHX:D=N:B=10:ID=21037930:S=Cercopithecidae])"
"CatarrhiniENSGT00390000008796.b:0.006894[&&NHX:D=N:B=22:ID=21037911:S=Catarrhini],"
"(ENSCJAG00000044082:0.016559[&&NHX:D=N:ID=1048310:S=Callithrix.jacchus],"
"ENSANAG00000017674:0.045132[&&NHX:D=N:ID=1048297:S=Aotus.nancymaae],"
"(ENSSBOG00000023270:0.007511[&&NHX:D=N:ID=1048316:S=Saimiri.boliviensis.boliviensis],"
"ENSCCAG00000030416:0.028002[&&NHX:D=N:ID=1048273:S=Cebus.capucinus.imitator])"
"CebidaeENSGT00390000008796.b:0.007511[&&NHX:_r=1:D=N:ID=100371849:S=Cebidae])"
"PlatyrrhiniENSGT00390000008796.b:0.029506[&&NHX:_f=_21037904_.21037902_:D=N:B=85:ID=21037899:S=Platyrrhini])"
"SimiiformesENSGT00390000008796.b:0.041215[&&NHX:D=N:B=87:ID=21037897:S=Simiiformes],"
"((ENSPCOG00000024150:0.024068[&&NHX:D=N:ID=1048225:S=Propithecus.coquereli])"
"LemuriformesENSGT00390000008796.b:0.0137224[&&NHX:reinserted=True:S=Lemuriformes])"
"StrepsirrhiniENSGT00390000008796.b:0.00924464[&&NHX:reinserted=True:S=Strepsirrhini])"
"PrimatesENSGT00390000008796.b:0.00897[&&NHX:D=N:B=1:ID=21037891:S=Primates],"
"(((((((ENSRNOG00000057827:0.040595[&&NHX:D=N:ID=1048237:S=Rattus.norvegicus],"
"(((ENSMUSG00000114865:0.0110201[&&NHX:D=N:ID=1048263:S=Mus.musculus])"
"Mus_BENSGT00390000008796.b:0.0156274[&&NHX:reinserted=True:S=Mus_B])"
"Mus_AENSGT00390000008796.b:0.00315714[&&NHX:reinserted=True:S=Mus_A])"
"MusENSGT00390000008796.b:0.0452834[&&NHX:reinserted=True:S=Mus])"
"MurinaeENSGT00390000008796.b:0.0210307[&&NHX:D=N:B=100:ID=21037958:S=Murinae])"
"EumuroidaENSGT00390000008796.b:0.0225119[&&NHX:reinserted=True:S=Eumuroida])"
"MuroideaENSGT00390000008796.b:0.017025[&&NHX:reinserted=True:S=Muroidea])"
"MyomorphaENSGT00390000008796.b:0.0269594[&&NHX:reinserted=True:S=Myomorpha])"
"MyodontaCastorimorphaENSGT00390000008796.b:0.00531889[&&NHX:reinserted=True:S=MyodontaCastorimorpha])"
"RodentiaENSGT00390000008796.b:0.0165444[&&NHX:reinserted=True:S=Rodentia])"
"GliresENSGT00390000008796.b:0.0137197[&&NHX:reinserted=True:S=Glires])"
"EuarchontogliresENSGT00390000008796.b:0.00897[&&NHX:_r=1:D=N:ID=100371847:S=Euarchontoglires],"
"(((ENSECAG00000012205:0.075502[&&NHX:D=N:ID=1048215:S=Equus.caballus],"
"(((ENSPTIG00000013858:0[&&NHX:D=N:ID=1890960:S=Panthera.tigris.altaica],"
"ENSPPRG00000007607:0.001074[&&NHX:D=N:ID=1048254:S=Panthera.pardus])"
"PantheraENSGT00390000008796.b:0.007366[&&NHX:D=N:B=100:ID=21037889:S=Panthera],"
"ENSFCAG00000045976:0.002402[&&NHX:D=N:ID=1048235:S=Felis.catus])"
"FelidaeENSGT00390000008796.b:0.0291095[&&NHX:D=N:B=100:ID=21037888:S=Felidae])"
"CarnivoraENSGT00390000008796.b:0.0170325[&&NHX:reinserted=True:S=Carnivora])"
"ZooamataENSGT00390000008796.b:0.00081727[&&NHX:D=N:B=45:ID=21037886:S=Zooamata])"
"FereuungulataENSGT00390000008796.b:0.00119689[&&NHX:reinserted=True:S=Fereuungulata])"
"ScrotiferaENSGT00390000008796.b:0.0277398[&&NHX:reinserted=True:S=Scrotifera])"
"BoreoeutheriaENSGT00390000008796.b:0.0202141[&&NHX:_f=_21037885_:D=N:SIS=0.0:B=6:ID=21037884:S=Boreoeutheria])"
"EutheriaENSGT00390000008796.b:0.119381[&&NHX:reinserted=True:S=Eutheria])"
"TheriaENSGT00390000008796.b:0.414697[&&NHX:D=N:B=81:ID=21037882:S=Theria])"
"AmniotaENSGT00390000008796.b:0.333856[&&NHX:reinserted=True:S=Amniota],"
"(((ENSAMXG00000003887:0.262097[&&NHX:D=N:ID=1048244:S=Astyanax.mexicanus],"
"ENSDARG00000098915:0.362318[&&NHX:D=N:ID=1048272:S=Danio.rerio])"
"OtophysiENSGT00390000008796.b:0.130069[&&NHX:D=N:B=100:ID=21037878:S=Otophysi],"
"(((ENSORLG00000007295:0.281804[&&NHX:D=N:ID=1048216:S=Oryzias.latipes],"
"(ENSPFOG00000009690:0.13399[&&NHX:D=N:ID=1048276:S=Poecilia.formosa])"
"PoeciliinaeENSGT00390000008796.b:0.0850306[&&NHX:reinserted=True:S=Poeciliinae])"
"AtherinomorphaeENSGT00390000008796.b:0.0361841[&&NHX:D=N:B=100:ID=21037879:S=Atherinomorphae])"
"OvalentariaENSGT00390000008796.b:0.0126088[&&NHX:reinserted=True:S=Ovalentaria])"
"PercomorphaceaeENSGT00390000008796.b:0.142764[&&NHX:reinserted=True:S=Percomorphaceae])"
"ClupeocephalaENSGT00390000008796.b:0.205049[&&NHX:D=N:B=100:ID=21037877:S=Clupeocephala])"
"NeopterygiiENSGT00390000008796.b:0.291773[&&NHX:reinserted=True:S=Neopterygii])"
"EuteleostomiENSGT00390000008796.b:0.762584[&&NHX:D=N:B=8:ID=21037876:S=Euteleostomi];"
)
)


# ALE output data
#
# Directory: ENSGT00390000001638/subtreesGoodQualO2/realign/
# or:        ENSGT00910000143978/subtreesGoodQualO2/realign/
#
#PhylTree-Simii.TimeTree201901.Ensembl93-like.goodQual.shortnames.nwk_SimiiformesENSGT00390000001638_pb.treelist.shortnames.ale.uml_rec
#PhylTree-Simii.TimeTree201901.Ensembl93-like.goodQual.shortnames.nwk_SimiiformesENSGT00390000001638_pb.treelist.shortnames.ale.uTs
#PhylTree.TimeTree201901.Ensembl93-like.goodQual.binary.basic.shortnames.nwk_SimiiformesENSGT00390000001638_pb.treelist.shortnames.ale.cons_tree
#PhylTree.TimeTree201901.Ensembl93-like.goodQual.binary.basic.shortnames.nwk_SimiiformesENSGT00390000001638_pb.treelist.shortnames.ale.ml_rec
#PhylTree.TimeTree201901.Ensembl93-like.goodQual.binary.basic.shortnames.nwk_SimiiformesENSGT00390000001638_pb.treelist.shortnames.ale.Ts

# Reconciled tree from ALE (SimiiformesENSGT00390000001638):
ALE_PHYLTREE = newtestfile(
    'ALE_PhylTree.nwk',
"(Scere:1,((Dmela:0.672362,Celeg:0.672362)69:0.0484642,(((Loc:0.284767,((Dar:0.138092,Amx:0.138092)61:0.0699537,((Oni:0.107686,(Orl:0.0843143,(Pfo:0.0515809,Xma:0.0515809)41:0.0327333)53:0.0233722)59:0.00814435,(Gac:0.0995421,(Tni:0.0466943,Tru:0.0466943)38:0.0528478)58:0.0162887)60:0.0922149)63:0.076721)66:0.10917,((Mod:0.14352,(Laf:0.0954336,((Mlu:0.0710629,((Ssc:0.0560748,(Bta:0.0222612,(Chi:0.00882305,Oar:0.00882305)11:0.0134382)24:0.0338136)43:0.0142878,(((Fca:0.0137323,(Pti:0.00670764,Ppr:0.00670764)7:0.00702465)16:0.0354247,(Caf:0.0411998,Mpu:0.0411998)37:0.00795712)39:0.0207276,Eca:0.0698845)48:0.000478103)49:0.000700214)50:0.0162287,((Ocu:0.0743315,((Sto:0.0658577,(Fda:0.0392363,((Cap:0.00511899,Cpo:0.00511899)5:0.0272584,(Ode:0.0297318,Cla:0.0297318)29:0.00264566)30:0.00685889)34:0.0266214)45:9.05111e-05,(Dor:0.063253,(Jja:0.0495923,(Nga:0.0409654,((Pem:0.0260619,(Moc:0.0259714,Mau:0.0259714)25:9.04928e-05)26:0.0034963,(Rno:0.0189016,(_pa:0.00750261,(_ca:0.00670787,(Mus:0.00277404,_sp:0.00277404)2:0.00393383)8:0.000794734)9:0.011399)23:0.0106566)28:0.0114072)36:0.00862686)40:0.0136608)44:0.00269516)46:0.00838333)51:0.006952,((Oga:0.0536843,(Pco:0.0341905,Mic:0.0341905)32:0.0194937)42:0.0131328,((Ana:0.01781,(Cja:0.0177195,(Sbo:0.0145427,Cca:0.0145427)18:0.00317684)20:9.04951e-05)21:0.0212388,((Nle:0.0182698,(Ppy:0.0142636,(Ggo:0.00820145,(Hsa:0.00601859,(Ppa:0.00255195,Ptr:0.00255195)1:0.00346664)6:0.00218286)10:0.00606219)17:0.00400615)22:0.00837268,(Can:0.0175758,((Cat:0.0112211,((Mle:0.0110401,Pan:0.0110401)12:9.04951e-05,((Mfa:0.00468587,Mmu:0.00468587)3:9.04883e-05,Mne:0.00477636)4:0.00635426)13:9.04917e-05)14:0.00122128,Csa:0.0124424)15:0.00513341)19:0.00906668)27:0.0124063)33:0.0277683)47:0.0144665)52:0.00600802)54:0.00814209)57:0.048086)62:0.138731,(Aca:0.253069,((Fal:0.0395454,Tgu:0.0395454)35:0.0491764,(Mga:0.0336633,Gal:0.0336633)31:0.0550585)55:0.164348)64:0.0291811)65:0.111686)67:0.2182,(Csav:0.0904928,Cin:0.0904928)56:0.521644)68:0.108689)70:0.279174)71;"
)

# With transfers marked as `T@`
ALE_NEWICK0 = newtestfile(
    'ALE_reconciled_genetree.nwk',
    "((Mic_ENSMICG00000008869:0.00291241,Pco_ENSPCOG00000010109:0.00238338).37.27:0.007607,(Can_ENSCANG00000006755@13|Can:0.000698673,((Ana_ENSANAG00000021473:0.00861211,(Cja_ENSCJAG00000000450:0.00355522,(Cca_ENSCCAG00000023954@9|Cca:0.000775261,Sbo_ENSSBOG00000025199@0|Sbo:0.00149833)T@12|-1:0.00162328)T@15|Cja:0.000966641).17:0.00166647,((Csa_ENSCSAG00000002856:0.00226518,(Mle_ENSMLEG00000036246:0.000714659,(Cat_ENSCATG00000039328:0.00112995,(Pan_ENSPANG00000019447@2|Pan:0.00382531,(Mne_ENSMNEG00000044003@4|Mne:0.0023071,(Mfa_ENSMFAG00000031618:0.000551396,Mmu_ENSMMUG00000005983:0.00232172)@4|3.3:0.00180815)T@7|-1:0.00256803)T@8|-1:0.000509848)T@9|Cat:0.00110414).11:0.000494).16.12:0.00100459,(Nle_ENSNLEG00000029729:0.00636055,(Ppy_ENSPPYG00000010198:0.00733296,(Ggo_ENSGGOG00000013433:0.00386273,(Hsa_ENSG00000105464:0.000460522,(Ppa_ENSPPAG00000034600:0.000663474,Ptr_ENSPTRG00000048202:0.0182578).1:0.000854063).5:0.00011525).9:0.000992671).14:0.00109895).18:0.00189119).22:0.00102454).28:0.00270523)T@38|28:0.007607).41:0;")

ALE_TREEBEST = newtestfile(
    'ALE_treebest.nwk',
"((ENSMICG00000008869:0.00291241,ENSPCOG00000010109:0.00238338).37.27:0.007607,(ENSCANG00000006755@13|Can:0.000698673,((ENSANAG00000021473:0.00861211,(ENSCJAG00000000450:0.00355522,(ENSCCAG00000023954@9|Cca:0.000775261,ENSSBOG00000025199@0|Sbo:0.00149833)T@12|-1:0.00162328)T@15|Cja:0.000966641).17:0.00166647,((ENSCSAG00000002856:0.00226518,(ENSMLEG00000036246:0.000714659,(ENSCATG00000039328:0.00112995,(ENSPANG00000019447@2|Pan:0.00382531,(ENSMNEG00000044003@4|Mne:0.0023071,(ENSMFAG00000031618:0.000551396,ENSMMUG00000005983:0.00232172)@4|3.3:0.00180815)T@7|-1:0.00256803)T@8|-1:0.000509848)T@9|Cat:0.00110414).11:0.000494).16.12:0.00100459,(ENSNLEG00000029729:0.00636055,(ENSPPYG00000010198:0.00733296,(ENSGGOG00000013433:0.00386273,(ENSG00000105464:0.000460522,(ENSPPAG00000034600:0.000663474,ENSPTRG00000048202:0.0182578).1:0.000854063).5:0.00011525).9:0.000992671).14:0.00109895).18:0.00189119).22:0.00102454).28:0.00270523)T@38|28:0.007607).41:0;"
)

# Another ALE experiment

ALE2_PHYLTREE = newtestfile("ALE2_Stree.nwk",
"((((((ANASP:0.171429,ANAVT:0.171429)6:0.342857,NOSP7:0.514286)18:0.371429,TRIEI:0.885714)31:0.0285714,((((CYAA5:0.428571,CYAP8:0.428571)15:0.285714,(CYAP7:0.685714,MICAN:0.685714)24:0.0285714)25:0.114286,SYNY3:0.828571)29:0.0285714,SYNP2:0.857143)30:0.0571429)32:0.0285714,((THEEB:0.742857,CYAP4:0.742857)26:0.0285714,ACAM1:0.771429)27:0.171429)33:0.0571429,(((SYNR3:0.657143,(((PROMM:0.114286,PROM3:0.114286)4:0.485714,((PRMAR1:0.4,PROM4:0.4)14:0.0857143,((PROM1:0.0571429,PROMT:0.0571429)2:0.4,((PROM9:0.257143,(PROM2:0.228571,(PROM0:0.2,PROMS:0.2)7:0.0285714)8:0.0285714)9:0.0285714,(PROM5:0.142857,PROMP:0.142857)5:0.142857)10:0.171429)16:0.0285714)17:0.114286)21:0.0285714,((SYNPW:0.314286,SYNS3:0.314286)11:0.228571,(SYNPX:0.371429,(SYNS9:0.342857,SYNSC:0.342857)12:0.0285714)13:0.171429)19:0.0857143)22:0.0285714)23:0.142857,(SYNE7:0.0285714,SYNP6:0.0285714)1:0.771429)28:0.171429,((SYNJA:0.0857143,SYNJB:0.0857143)3:0.485714,GLVIO1:0.571429)20:0.4)34:0.0285714)35;")

ALE2_REC_TREEBEST = newtestfile("ALE2_rec_treebest.nwk", 
        "((((PE950:0.0567907[&&NHX:S=SYNJB:D=0],PE2030:0.0454848[&&NHX:S=SYNJA:D=0])3:0.157165[&&NHX:S=3:D=0]):0.157165[&&NHX:S=20:D=0]):0.157165[&&NHX:S=34:D=0],(((PE420:0.266217[&&NHX:S=CYAP4:D=0],(PE1494:0.334417[&&NHX:S=THEEB:D=0],(PE1169:0.362213[&&NHX:S=GLVIO1:D=0]):0.362213[&&NHX:S=GLVIO1:D=0:T=-1])THEEB:0.131324[&&NHX:S=THEEB:D=11:T=1])26:0.099468[&&NHX:S=26:D=0]):0.099468[&&NHX:S=27:D=0],((PE668:0.263141[&&NHX:S=TRIEI:D=0],(PE4607:0.0631294[&&NHX:S=NOSP7:D=0],(PE931:0.0147708[&&NHX:S=ANAVT:D=0],PE2969:0.0205304[&&NHX:S=ANASP:D=0])6:0.0482404[&&NHX:S=6:D=0])18:0.147972[&&NHX:S=18:D=0])31:0.0665971[&&NHX:S=31:D=0],(((PE135:0.0836723[&&NHX:S=ACAM1:D=0]):0.0836723[&&NHX:S=ACAM1:D=0:T=-1]):0.0836723[&&NHX:S=SYNP2:D=11:T=1],(((PE594:0.221589[&&NHX:S=CYAA5:D=0],PE2632:0.176689[&&NHX:S=CYAP8:D=0])15:0.054886[&&NHX:S=15:D=0],(PE2280:0.17352[&&NHX:S=CYAP7:D=0],PE2716:0.244643[&&NHX:S=MICAN:D=0])24:0.0380066[&&NHX:S=24:D=0])25:0.0613626[&&NHX:S=25:D=0],(PE2642:0.22923[&&NHX:S=SYNY3:D=0],(((PE221:0.00291872[&&NHX:S=SYNE7:D=0],PE1291:0.00296591[&&NHX:S=SYNP6:D=0])1:0.140013[&&NHX:S=1:D=0]):0.140013[&&NHX:S=1:D=0:T=-1],(PE944:0.670514[&&NHX:S=SYNR3:D=0],(((PE774:0.182687[&&NHX:S=SYNPW:D=0],PE1833:0.162686[&&NHX:S=SYNS3:D=0])11:0.056033[&&NHX:S=11:D=0]):0.056033[&&NHX:S=11:D=0:T=-1],(((PE1481:0.117437[&&NHX:S=SYNPX:D=0],(PE1008:0.161326[&&NHX:S=SYNSC:D=0],PE925:0.146404[&&NHX:S=SYNS9:D=0])12:0.0190145[&&NHX:S=12:D=0])13:0.0669585[&&NHX:S=13:D=0]):0.0669585[&&NHX:S=19:D=0],((PE1845:0.00426106[&&NHX:S=PROM3:D=0],PE495:0.0101808[&&NHX:S=PROMM:D=0])4:0.144436[&&NHX:S=4:D=0],(((PE651:0.00569711[&&NHX:S=PROM1:D=0],PE635:0.0158877[&&NHX:S=PROMT:D=0])2:0.18974[&&NHX:S=2:D=0]):0.18974[&&NHX:S=16:D=0],(PE1063:0.316481[&&NHX:S=PRMAR1:D=0],(PE1056:0.246757[&&NHX:S=PROM4:D=0],((PE662:0.079947[&&NHX:S=PROM5:D=0]):0.079947[&&NHX:S=PROM5:D=0:T=-1],((PE651:0.0644765[&&NHX:S=PROMP:D=0]):0.0644765[&&NHX:S=5:D=0],(PE623:0.0631102[&&NHX:S=PROM9:D=0],((PE651:0.00945917[&&NHX:S=PROMS:D=0],(PE623:0.0104173[&&NHX:S=PROM0:D=0],(PE676:0.0140147[&&NHX:S=PROM2:D=0]):0.0140147[&&NHX:S=PROM2:D=0:T=-1])PROM0:0.00690004[&&NHX:S=PROM0:D=11:T=1])7:0.0111632[&&NHX:S=7:D=0]):0.0111632[&&NHX:S=8:D=0])9:0.0287127[&&NHX:S=9:D=0])10:0.0655906[&&NHX:S=10:D=0])10:0.547024[&&NHX:S=10:D=12:T=0])PROM4:0.116101[&&NHX:S=PROM4:D=11:T=1])14:0.0755022[&&NHX:S=14:D=0])17:0.191706[&&NHX:S=17:D=0])21:0.0627178[&&NHX:S=21:D=0])22:0.053263[&&NHX:S=22:D=0])22:0.0960767[&&NHX:S=22:D=11:T=1])23:0.810502[&&NHX:S=23:D=0])23:0.125576[&&NHX:S=23:D=12:T=0])SYNY3:0.118545[&&NHX:S=SYNY3:D=11:T=1])29:0.0848267[&&NHX:S=29:D=0])30:0.0585885[&&NHX:S=30:D=0])32:0.0795774[&&NHX:S=32:D=0])33:0.471496[&&NHX:S=33:D=0]):1[&&NHX:S=35:D=0];")

ALE2_REC_TREEBEST_TEST1 = newtestfile("ALE2_rec_treebest_test1.nwk",
        "((((PE950:0.0567907[&&NHX:S=SYNJB:D=0],PE2030:0.0454848[&&NHX:S=SYNJA:D=0])3:0.157165[&&NHX:S=3:D=0]):0.157165[&&NHX:S=20:D=0]):0.157165[&&NHX:S=34:D=0],(((PE420:0.266217[&&NHX:S=CYAP4:D=0],(PE1494:0.334417[&&NHX:S=THEEB:D=0],(PE1169:0.362213[&&NHX:S=GLVIO1:D=0]):0.362213[&&NHX:S=GLVIO1:D=0:T=-1])THEEB:0.131324[&&NHX:S=THEEB:D=11:T=1])26:0.099468[&&NHX:S=26:D=0]):0.099468[&&NHX:S=27:D=0],((PE668:0.263141[&&NHX:S=TRIEI:D=0],(PE4607:0.0631294[&&NHX:S=NOSP7:D=0],(PE931:0.0147708[&&NHX:S=ANAVT:D=0],PE2969:0.0205304[&&NHX:S=ANASP:D=0])6:0.0482404[&&NHX:S=6:D=0])18:0.147972[&&NHX:S=18:D=0])31:0.0665971[&&NHX:S=31:D=0],((PE135:0.125508[&&NHX:S=ACAM1:D=0]):0.125508[&&NHX:S=SYNP2@26|ACAM1:D=11:T=1],(((PE594:0.221589[&&NHX:S=CYAA5:D=0],PE2632:0.176689[&&NHX:S=CYAP8:D=0])15:0.054886[&&NHX:S=15:D=0],(PE2280:0.17352[&&NHX:S=CYAP7:D=0],PE2716:0.244643[&&NHX:S=MICAN:D=0])24:0.0380066[&&NHX:S=24:D=0])25:0.0613626[&&NHX:S=25:D=0],(PE2642:0.22923[&&NHX:S=SYNY3:D=0],(((PE221:0.00291872[&&NHX:S=SYNE7:D=0],PE1291:0.00296591[&&NHX:S=SYNP6:D=0])1:0.140013[&&NHX:S=1:D=0]):0.140013[&&NHX:S=1:D=0:T=-1],(PE944:0.670514[&&NHX:S=SYNR3:D=0],(((PE774:0.182687[&&NHX:S=SYNPW:D=0],PE1833:0.162686[&&NHX:S=SYNS3:D=0])11:0.056033[&&NHX:S=11:D=0]):0.056033[&&NHX:S=11:D=0:T=-1],(((PE1481:0.117437[&&NHX:S=SYNPX:D=0],(PE1008:0.161326[&&NHX:S=SYNSC:D=0],PE925:0.146404[&&NHX:S=SYNS9:D=0])12:0.0190145[&&NHX:S=12:D=0])13:0.0669585[&&NHX:S=13:D=0]):0.0669585[&&NHX:S=19:D=0],((PE1845:0.00426106[&&NHX:S=PROM3:D=0],PE495:0.0101808[&&NHX:S=PROMM:D=0])4:0.144436[&&NHX:S=4:D=0],(((PE651:0.00569711[&&NHX:S=PROM1:D=0],PE635:0.0158877[&&NHX:S=PROMT:D=0])2:0.18974[&&NHX:S=2:D=0]):0.18974[&&NHX:S=16:D=0],(PE1063:0.316481[&&NHX:S=PRMAR1:D=0],(PE1056:0.246757[&&NHX:S=PROM4:D=0],((PE662:0.079947[&&NHX:S=PROM5:D=0]):0.079947[&&NHX:S=PROM5:D=0:T=-1],((PE651:0.0644765[&&NHX:S=PROMP:D=0]):0.0644765[&&NHX:S=5:D=0],(PE623:0.0631102[&&NHX:S=PROM9:D=0],((PE651:0.00945917[&&NHX:S=PROMS:D=0],(PE623:0.0104173[&&NHX:S=PROM0:D=0],(PE676:0.0140147[&&NHX:S=PROM2:D=0]):0.0140147[&&NHX:S=PROM2:D=0:T=-1])PROM0:0.00690004[&&NHX:S=PROM0:D=11:T=1])7:0.0111632[&&NHX:S=7:D=0]):0.0111632[&&NHX:S=8:D=0])9:0.0287127[&&NHX:S=9:D=0])10:0.0655906[&&NHX:S=10:D=0])10:0.547024[&&NHX:S=10:D=12:T=0])PROM4:0.116101[&&NHX:S=PROM4:D=11:T=1])14:0.0755022[&&NHX:S=14:D=0])17:0.191706[&&NHX:S=17:D=0])21:0.0627178[&&NHX:S=21:D=0])22:0.053263[&&NHX:S=22:D=0])22:0.0960767[&&NHX:S=22:D=11:T=1])23:0.810502[&&NHX:S=23:D=0])23:0.125576[&&NHX:S=23:D=12:T=0])SYNY3:0.118545[&&NHX:S=SYNY3:D=11:T=1])29:0.0848267[&&NHX:S=29:D=0])30:0.0585885[&&NHX:S=30:D=0])32:0.0795774[&&NHX:S=32:D=0])33:0.471496[&&NHX:S=33:D=0]):1[&&NHX:S=35:D=0];")

# SIMULATED DATA (simul-scrollsaw)

EUKA_2SUPERGROUPS = newtestfile("euka_2supergroups.nwk", "(Amorphea:2000,Diphoda:2000)Eukaryota:100;")

EUKA_10SPECIES = newtestfile("euka_10species.nwk",
        "(((Acanthamoeba:900,Acytostelium:900)Amoebozoa:600,((Basidiobolus:1100,Nematostella:1100)Opisthokonta:200,(Thecamonas:900)Apusozoa:400)Obazoa:200)Amorphea:500,(((Naeglaria:1000,Euglena:1000)Excavata:500)Discoba:400,((Saprolegnia:1100,Plasmodiophora:1100)SAR:400,(Klebsormidium:1200,Physcomistrella:1200)Archaeplastida:300,(Emiliana:700,Chrysochromulina:700)Haptophyta:800)Diaphoretickes:400)Diphoda:100)Eukaryota;")

SIM5 = newtestfile("sim5.nwk",
        "(Amorphea_G001:-5839[&&NHX:S=Amorphea],(Diphoda_lost:0[&&NHX:S=Diphoda:D=1],(Diphoda_G005:304.952[&&NHX:S=Diphoda],Diphoda_G006:304.952[&&NHX:S=Diphoda]):0[&&NHX:S=Diphoda:D=1]):-5839[&&NHX:S=Diphoda:D=1])r:1[&&NHX:S=Eukaryota];"
"((Amorphea_G005:0[&&NHX:S=Amorphea],Amorphea_G006:0[&&NHX:S=Amorphea]):-199.303[&&NHX:S=Amorphea:D=1],((Diphoda_G007:708.515[&&NHX:S=Diphoda],(Diphoda_G010:818.834[&&NHX:S=Diphoda],Diphoda_G011:818.834[&&NHX:S=Diphoda]):708.515[&&NHX:S=Diphoda:D=1]):0[&&NHX:S=Diphoda:D=1],(((Diphoda_G015:616.67[&&NHX:S=Diphoda],Diphoda_G016:616.67[&&NHX:S=Diphoda]):457.423[&&NHX:S=Diphoda:D=1],Diphoda_G012:457.423[&&NHX:S=Diphoda]):565.501[&&NHX:S=Diphoda:D=1],Diphoda_G008:565.501[&&NHX:S=Diphoda]):0[&&NHX:S=Diphoda:D=1]):-199.303[&&NHX:S=Diphoda:D=1])r:1[&&NHX:S=Eukaryota];"
"(Amorphea_G002:-636.601[&&NHX:S=Amorphea],Diphoda_G001:-636.601[&&NHX:S=Diphoda])r:1[&&NHX:S=Eukaryota];"
"(Amorphea_G003:-2259.48[&&NHX:S=Amorphea],(Diphoda_G002:0[&&NHX:S=Diphoda],Diphoda_G003:0[&&NHX:S=Diphoda]):-2259.48[&&NHX:S=Diphoda:D=1])r:1[&&NHX:S=Eukaryota];"
"(Amorphea_G004:-712.099[&&NHX:S=Amorphea],((Diphoda_G009:1427.3[&&NHX:S=Diphoda],(Diphoda_G013:129.762[&&NHX:S=Diphoda],Diphoda_G014:129.762[&&NHX:S=Diphoda]):1427.3[&&NHX:S=Diphoda:D=1]):0[&&NHX:S=Diphoda:D=1],Diphoda_G004:0[&&NHX:S=Diphoda]):-712.099[&&NHX:S=Diphoda:D=1])r:1[&&NHX:S=Eukaryota];")
SIM5_TRANSFERS = newtestfile("sim5_t.nwk",
        "((((Amorphea_G014:319.798[&&NHX:S=Amorphea],Amorphea_G015:319.798[&&NHX:S=Amorphea]):196.525[&&NHX:S=Amorphea:D=1],Amorphea_G009:196.525[&&NHX:S=Amorphea]):0[&&NHX:S=Amorphea:D=1],(Amorphea_G010:120.726[&&NHX:S=Amorphea],Amorphea_G011:120.726[&&NHX:S=Amorphea]):0[&&NHX:S=Amorphea:D=1]):-1723.21[&&NHX:S=Amorphea:D=1],((Diphoda_G004:1234.17[&&NHX:S=Diphoda],Amorphea_G005:1234.17[&&NHX:S=Amorphea:D=0:T=-1]):0[&&NHX:S=Diphoda:D=11:T=1],(Amorphea_G006:776.043[&&NHX:S=Amorphea],Diphoda_G005:776.043[&&NHX:S=Diphoda:D=0:T=-1]):0[&&NHX:S=Amorphea:D=11:T=1]):-1723.21[&&NHX:S=Diphoda:D=11:T=1])r:1[&&NHX:S=Eukaryota];"
"(Amorphea_G001:-1259.05[&&NHX:S=Amorphea],Diphoda_lost:-1259.05[&&NHX:S=Diphoda:D=1])r:1[&&NHX:S=Eukaryota];"
"((Amorphea_G003:0[&&NHX:S=Amorphea],(Amorphea_lost:976.962[&&NHX:S=Amorphea:D=1],(Amorphea_lost:106.352[&&NHX:S=Amorphea:D=1],Amorphea_G016:106.352[&&NHX:S=Amorphea]):976.962[&&NHX:S=Amorphea:D=1]):0[&&NHX:S=Amorphea:D=1]):-4621.3[&&NHX:S=Amorphea:D=1],Diphoda_G001:-4621.3[&&NHX:S=Diphoda])r:1[&&NHX:S=Eukaryota];"
"(Amorphea_G002:-207.341[&&NHX:S=Amorphea],(Diphoda_G003:0[&&NHX:S=Diphoda],(Amorphea_G007:1533.5[&&NHX:S=Amorphea],Amorphea_G008:1533.5[&&NHX:S=Amorphea]):0[&&NHX:S=Amorphea:D=1:T=-1]):-207.341[&&NHX:S=Diphoda:D=11:T=1])r:1[&&NHX:S=Eukaryota];"
"((Amorphea_G004:0[&&NHX:S=Amorphea],(Amorphea_G012:142.255[&&NHX:S=Amorphea],Amorphea_G013:142.255[&&NHX:S=Amorphea]):0[&&NHX:S=Amorphea:D=1]):-1722.5[&&NHX:S=Amorphea:D=1],Diphoda_G002:-1722.5[&&NHX:S=Diphoda])r:1[&&NHX:S=Eukaryota];")
SIM5_10SP = newtestfile("sim5_10sp.nwk",
        "((((Acanthamoeba_G005:53.4677[&&NHX:S=Acanthamoeba],Acytostelium_G004:53.4677[&&NHX:S=Acytostelium]):0[&&NHX:S=Amoebozoa],(Acanthamoeba_G006:53.4677[&&NHX:S=Acanthamoeba],(Acytostelium_G008:0[&&NHX:S=Acytostelium],Acytostelium_G009:0[&&NHX:S=Acytostelium]):53.4677[&&NHX:S=Acytostelium:D=1]):0[&&NHX:S=Amoebozoa]):-167.071[&&NHX:S=Amoebozoa:D=1],((Basidiobolus_G001:-2816.57[&&NHX:S=Basidiobolus],Nematostella_G001:-2816.57[&&NHX:S=Nematostella]):-3.00033[&&NHX:S=Opisthokonta],(Thecamonas_G001:-456.71[&&NHX:S=Thecamonas]):-3.00033[&&NHX:S=Apusozoa]):-167.071[&&NHX:S=Obazoa]):-128.887[&&NHX:S=Amorphea],(((Naeglaria_G001:-3322.52[&&NHX:S=Naeglaria],Euglena_G001:-3322.52[&&NHX:S=Euglena]):-568.899[&&NHX:S=Excavata]):-5026.72[&&NHX:S=Discoba],((Saprolegnia_G001:-65.754[&&NHX:S=Saprolegnia],Plasmodiophora_G001:-65.754[&&NHX:S=Plasmodiophora]):-3615.69[&&NHX:S=SAR],((Klebsormidium_G005:0[&&NHX:S=Klebsormidium],Klebsormidium_lost:0[&&NHX:S=Klebsormidium:D=1]):-1176.41[&&NHX:S=Klebsormidium:D=1],(Physcomistrella_G006:0[&&NHX:S=Physcomistrella],(Physcomistrella_G011:0[&&NHX:S=Physcomistrella],Physcomistrella_G012:0[&&NHX:S=Physcomistrella]):0[&&NHX:S=Physcomistrella:D=1]):-1176.41[&&NHX:S=Physcomistrella:D=1]):-3615.69[&&NHX:S=Archaeplastida],((Emiliana_G005:-3321.2[&&NHX:S=Emiliana],Chrysochromulina_G005:-3321.2[&&NHX:S=Chrysochromulina]):0[&&NHX:S=Haptophyta],((Emiliana_G006:0[&&NHX:S=Emiliana],Emiliana_G007:0[&&NHX:S=Emiliana]):-1309.7[&&NHX:S=Emiliana:D=1],Chrysochromulina_G006:-1309.7[&&NHX:S=Chrysochromulina]):0[&&NHX:S=Haptophyta]):-3615.69[&&NHX:S=Haptophyta:D=1]):-5026.72[&&NHX:S=Diaphoretickes]):-128.887[&&NHX:S=Diphoda])r:1[&&NHX:S=Eukaryota];"
"(((Acanthamoeba_G001:-877.041[&&NHX:S=Acanthamoeba],(Acytostelium_G005:0[&&NHX:S=Acytostelium],Acytostelium_G006:0[&&NHX:S=Acytostelium]):-877.041[&&NHX:S=Acytostelium:D=1]):-6626.17[&&NHX:S=Amoebozoa],Obazoa_lost:-6626.17[&&NHX:S=Obazoa:D=1]):-1561.82[&&NHX:S=Amorphea],(((Naeglaria_G002:-519.752[&&NHX:S=Naeglaria],Euglena_G002:-519.752[&&NHX:S=Euglena]):-20.2037[&&NHX:S=Excavata]):-3588.03[&&NHX:S=Discoba],((Saprolegnia_G002:-1167.64[&&NHX:S=Saprolegnia],Plasmodiophora_G002:-1167.64[&&NHX:S=Plasmodiophora]):-173.232[&&NHX:S=SAR],((Klebsormidium_G006:0[&&NHX:S=Klebsormidium],Klebsormidium_G007:0[&&NHX:S=Klebsormidium]):-843.24[&&NHX:S=Klebsormidium:D=1],Physcomistrella_G001:-843.24[&&NHX:S=Physcomistrella]):-173.232[&&NHX:S=Archaeplastida],((Emiliana_G003:-531.166[&&NHX:S=Emiliana],Chrysochromulina_G003:-531.166[&&NHX:S=Chrysochromulina]):0[&&NHX:S=Haptophyta],(Emiliana_G004:-7330.87[&&NHX:S=Emiliana],Chrysochromulina_G004:-7330.87[&&NHX:S=Chrysochromulina]):0[&&NHX:S=Haptophyta]):-173.232[&&NHX:S=Haptophyta:D=1]):-3588.03[&&NHX:S=Diaphoretickes]):-1561.82[&&NHX:S=Diphoda])r:1[&&NHX:S=Eukaryota];"
"((((Acanthamoeba_G004:-1227.29[&&NHX:S=Acanthamoeba],Acytostelium_G003:-1227.29[&&NHX:S=Acytostelium]):56.2937[&&NHX:S=Amoebozoa],((Basidiobolus_G003:-1304.22[&&NHX:S=Basidiobolus],Nematostella_G003:-1304.22[&&NHX:S=Nematostella]):-4618.02[&&NHX:S=Opisthokonta],Apusozoa_lost:-4618.02[&&NHX:S=Apusozoa:D=1]):56.2937[&&NHX:S=Obazoa]):0[&&NHX:S=Amorphea],(((Acanthamoeba_lost:0[&&NHX:S=Acanthamoeba:D=1],Acanthamoeba_G007:0[&&NHX:S=Acanthamoeba]):-1765.95[&&NHX:S=Acanthamoeba:D=1],((Acytostelium_G010:167.125[&&NHX:S=Acytostelium],Acytostelium_G011:167.125[&&NHX:S=Acytostelium]):0[&&NHX:S=Acytostelium:D=1],Acytostelium_G007:0[&&NHX:S=Acytostelium]):-1765.95[&&NHX:S=Acytostelium:D=1]):56.2937[&&NHX:S=Amoebozoa],((Basidiobolus_G004:-2420.87[&&NHX:S=Basidiobolus],Nematostella_G004:-2420.87[&&NHX:S=Nematostella]):-2549.11[&&NHX:S=Opisthokonta],(Thecamonas_G003:-2749.77[&&NHX:S=Thecamonas]):-2549.11[&&NHX:S=Apusozoa]):56.2937[&&NHX:S=Obazoa]):0[&&NHX:S=Amorphea]):-1808.89[&&NHX:S=Amorphea:D=1],((((Naeglaria_G005:-3630.03[&&NHX:S=Naeglaria],Euglena_G005:-3630.03[&&NHX:S=Euglena]):0[&&NHX:S=Excavata],(Naeglaria_G006:-3460.47[&&NHX:S=Naeglaria],Euglena_G006:-3460.47[&&NHX:S=Euglena]):0[&&NHX:S=Excavata]):-1483.25[&&NHX:S=Excavata:D=1]):-118.138[&&NHX:S=Discoba],(((Saprolegnia_G005:0[&&NHX:S=Saprolegnia],Saprolegnia_G006:0[&&NHX:S=Saprolegnia]):-3523.73[&&NHX:S=Saprolegnia:D=1],(Plasmodiophora_G005:0[&&NHX:S=Plasmodiophora],Plasmodiophora_G006:0[&&NHX:S=Plasmodiophora]):-3523.73[&&NHX:S=Plasmodiophora:D=1]):-9563.07[&&NHX:S=SAR],(((Klebsormidium_G003:-188.753[&&NHX:S=Klebsormidium],(Physcomistrella_G007:0[&&NHX:S=Physcomistrella],Physcomistrella_G008:0[&&NHX:S=Physcomistrella]):-188.753[&&NHX:S=Physcomistrella:D=1]):13.6299[&&NHX:S=Archaeplastida],(Klebsormidium_lost:-541.196[&&NHX:S=Klebsormidium:D=1],Physcomistrella_G004:-541.196[&&NHX:S=Physcomistrella]):13.6299[&&NHX:S=Archaeplastida]):0[&&NHX:S=Archaeplastida:D=1],(Klebsormidium_G002:-119.089[&&NHX:S=Klebsormidium],Physcomistrella_G003:-119.089[&&NHX:S=Physcomistrella]):0[&&NHX:S=Archaeplastida]):-9563.07[&&NHX:S=Archaeplastida:D=1],(Emiliana_G001:-3202.5[&&NHX:S=Emiliana],Chrysochromulina_G001:-3202.5[&&NHX:S=Chrysochromulina]):-9563.07[&&NHX:S=Haptophyta]):-118.138[&&NHX:S=Diaphoretickes]):-1808.89[&&NHX:S=Diphoda])r:1[&&NHX:S=Eukaryota];"
"(((Acanthamoeba_G002:-504.426[&&NHX:S=Acanthamoeba],Acytostelium_G001:-504.426[&&NHX:S=Acytostelium]):-2251.59[&&NHX:S=Amoebozoa],(((Basidiobolus_G005:85.2748[&&NHX:S=Basidiobolus],(Nematostella_G006:0[&&NHX:S=Nematostella],Nematostella_G007:0[&&NHX:S=Nematostella]):85.2748[&&NHX:S=Nematostella:D=1]):0[&&NHX:S=Opisthokonta],(Basidiobolus_G006:85.2748[&&NHX:S=Basidiobolus],Nematostella_G005:85.2748[&&NHX:S=Nematostella]):0[&&NHX:S=Opisthokonta]):-1944.65[&&NHX:S=Opisthokonta:D=1],(Thecamonas_G002:-1790.22[&&NHX:S=Thecamonas]):-1944.65[&&NHX:S=Apusozoa]):-2251.59[&&NHX:S=Obazoa]):-1103.33[&&NHX:S=Amorphea],(((Naeglaria_G003:-3474.03[&&NHX:S=Naeglaria],Euglena_G003:-3474.03[&&NHX:S=Euglena]):-903.248[&&NHX:S=Excavata]):-1828.31[&&NHX:S=Discoba],((Saprolegnia_lost:-11825.3[&&NHX:S=Saprolegnia:D=1],Plasmodiophora_G003:-11825.3[&&NHX:S=Plasmodiophora]):-699.355[&&NHX:S=SAR],(Klebsormidium_G001:-1658.37[&&NHX:S=Klebsormidium],((Physcomistrella_G013:0[&&NHX:S=Physcomistrella],Physcomistrella_G014:0[&&NHX:S=Physcomistrella]):0[&&NHX:S=Physcomistrella:D=1],Physcomistrella_lost:0[&&NHX:S=Physcomistrella:D=1]):-1658.37[&&NHX:S=Physcomistrella:D=1]):-699.355[&&NHX:S=Archaeplastida],Haptophyta_lost:-699.355[&&NHX:S=Haptophyta:D=1]):-1828.31[&&NHX:S=Diaphoretickes]):-1103.33[&&NHX:S=Diphoda])r:1[&&NHX:S=Eukaryota];"
"(((Acanthamoeba_G003:-4689.58[&&NHX:S=Acanthamoeba],Acytostelium_G002:-4689.58[&&NHX:S=Acytostelium]):-891.348[&&NHX:S=Amoebozoa],((Basidiobolus_G002:-1486.97[&&NHX:S=Basidiobolus],Nematostella_G002:-1486.97[&&NHX:S=Nematostella]):-3512.4[&&NHX:S=Opisthokonta],((Thecamonas_G004:0[&&NHX:S=Thecamonas],Thecamonas_G005:0[&&NHX:S=Thecamonas]):-5718.96[&&NHX:S=Thecamonas:D=1]):-3512.4[&&NHX:S=Apusozoa]):-891.348[&&NHX:S=Obazoa]):-2617.51[&&NHX:S=Amorphea],(((Naeglaria_G004:-707.106[&&NHX:S=Naeglaria],Euglena_G004:-707.106[&&NHX:S=Euglena]):-6269.53[&&NHX:S=Excavata]):-157.751[&&NHX:S=Discoba],(((Saprolegnia_G003:-1014.73[&&NHX:S=Saprolegnia],Plasmodiophora_G004:-1014.73[&&NHX:S=Plasmodiophora]):326.538[&&NHX:S=SAR],(Klebsormidium_lost:-408.7[&&NHX:S=Klebsormidium:D=1],Physcomistrella_G002:-408.7[&&NHX:S=Physcomistrella]):326.538[&&NHX:S=Archaeplastida],(Emiliana_G002:-6137.69[&&NHX:S=Emiliana],Chrysochromulina_lost:-6137.69[&&NHX:S=Chrysochromulina:D=1]):326.538[&&NHX:S=Haptophyta]):0[&&NHX:S=Diaphoretickes],(((Saprolegnia_G004:0[&&NHX:S=Saprolegnia],(Saprolegnia_G007:0[&&NHX:S=Saprolegnia],Saprolegnia_G008:0[&&NHX:S=Saprolegnia]):0[&&NHX:S=Saprolegnia:D=1]):-1642.72[&&NHX:S=Saprolegnia:D=1],(Plasmodiophora_G007:0[&&NHX:S=Plasmodiophora],Plasmodiophora_G008:0[&&NHX:S=Plasmodiophora]):-1642.72[&&NHX:S=Plasmodiophora:D=1]):326.538[&&NHX:S=SAR],(((Klebsormidium_G008:0[&&NHX:S=Klebsormidium],Klebsormidium_G009:0[&&NHX:S=Klebsormidium]):55.6654[&&NHX:S=Klebsormidium:D=1],(Physcomistrella_G009:0[&&NHX:S=Physcomistrella],Physcomistrella_G010:0[&&NHX:S=Physcomistrella]):55.6654[&&NHX:S=Physcomistrella:D=1]):0[&&NHX:S=Archaeplastida],(Klebsormidium_G004:55.6654[&&NHX:S=Klebsormidium],Physcomistrella_G005:55.6654[&&NHX:S=Physcomistrella]):0[&&NHX:S=Archaeplastida]):326.538[&&NHX:S=Archaeplastida:D=1],(((Emiliana_G008:394.334[&&NHX:S=Emiliana],Emiliana_G009:394.334[&&NHX:S=Emiliana]):0[&&NHX:S=Emiliana:D=1],(Emiliana_G010:247.002[&&NHX:S=Emiliana],Emiliana_G011:247.002[&&NHX:S=Emiliana]):0[&&NHX:S=Emiliana:D=1]):-1666.66[&&NHX:S=Emiliana:D=1],Chrysochromulina_G002:-1666.66[&&NHX:S=Chrysochromulina]):326.538[&&NHX:S=Haptophyta]):0[&&NHX:S=Diaphoretickes]):-157.751[&&NHX:S=Diaphoretickes:D=1]):-2617.51[&&NHX:S=Diphoda])r:1[&&NHX:S=Eukaryota];")
SIM5_10SP_TRANSFERS = newtestfile("sim5_10sp_t.nwk",
        """(((((Acanthamoeba_G004:-30.3998[&&NHX:S=Acanthamoeba],Acytostelium_G002:-30.3998[&&NHX:S=Acytostelium]):194.331[&&NHX:S=Amoebozoa],(((Basidiobolus_G006:0[&&NHX:S=Basidiobolus],(Basidiobolus_G009:85.3993[&&NHX:S=Basidiobolus],Basidiobolus_G010:85.3993[&&NHX:S=Basidiobolus]):0[&&NHX:S=Basidiobolus:D=1]):-1578.2[&&NHX:S=Basidiobolus:D=1],Nematostella_G004:-1578.2[&&NHX:S=Nematostella]):-1937.82[&&NHX:S=Opisthokonta],(Thecamonas_G002:-3122.9[&&NHX:S=Thecamonas]):-1937.82[&&NHX:S=Apusozoa]):194.331[&&NHX:S=Obazoa]):5.35171[&&NHX:S=Amorphea],((Saprolegnia_G002:-694.679[&&NHX:S=Saprolegnia],Plasmodiophora_G004:-694.679[&&NHX:S=Plasmodiophora]):194.331[&&NHX:S=SAR],(Klebsormidium_G004:-4230.09[&&NHX:S=Klebsormidium],Physcomistrella_G004:-4230.09[&&NHX:S=Physcomistrella]):194.331[&&NHX:S=Archaeplastida],Haptophyta_lost:194.331[&&NHX:S=Haptophyta:D=1]):5.35171[&&NHX:S=Diaphoretickes:D=0:T=-1]):0[&&NHX:S=Amorphea:D=11:T=1],(((Acanthamoeba_G005:-2409.34[&&NHX:S=Acanthamoeba],Acytostelium_G003:-2409.34[&&NHX:S=Acytostelium]):67.8463[&&NHX:S=Amoebozoa],(((Basidiobolus_G007:0[&&NHX:S=Basidiobolus],(Basidiobolus_G011:129.168[&&NHX:S=Basidiobolus],Basidiobolus_G012:129.168[&&NHX:S=Basidiobolus]):0[&&NHX:S=Basidiobolus:D=1]):-16.6164[&&NHX:S=Basidiobolus:D=1],(Nematostella_lost:0[&&NHX:S=Nematostella:D=1],Acanthamoeba_G009:0[&&NHX:S=Acanthamoeba:D=0:T=-1]):-16.6164[&&NHX:S=Nematostella:D=11:T=1]):-2143.37[&&NHX:S=Opisthokonta],((Thecamonas_G007:0[&&NHX:S=Thecamonas],Thecamonas_G008:0[&&NHX:S=Thecamonas]):-1300.64[&&NHX:S=Thecamonas:D=1]):-2143.37[&&NHX:S=Apusozoa]):67.8463[&&NHX:S=Obazoa]):131.836[&&NHX:S=Amorphea],((Saprolegnia_G003:-6361.88[&&NHX:S=Saprolegnia],Plasmodiophora_G005:-6361.88[&&NHX:S=Plasmodiophora]):67.8463[&&NHX:S=SAR],(Klebsormidium_G005:-1968.66[&&NHX:S=Klebsormidium],(Physcomistrella_G008:0[&&NHX:S=Physcomistrella],Physcomistrella_G009:0[&&NHX:S=Physcomistrella]):-1968.66[&&NHX:S=Physcomistrella:D=1]):67.8463[&&NHX:S=Archaeplastida],(Emiliana_G003:-2321.6[&&NHX:S=Emiliana],Chrysochromulina_G002:-2321.6[&&NHX:S=Chrysochromulina]):67.8463[&&NHX:S=Haptophyta]):131.836[&&NHX:S=Diaphoretickes:D=0:T=-1]):0[&&NHX:S=Amorphea:D=11:T=1]):-2050.28[&&NHX:S=Amorphea:D=1],(((((Naeglaria_G003:0[&&NHX:S=Naeglaria],(Physcomistrella_G016:0[&&NHX:S=Physcomistrella],Physcomistrella_G017:0[&&NHX:S=Physcomistrella]):0[&&NHX:S=Physcomistrella:D=1:T=-1]):-424.594[&&NHX:S=Naeglaria:D=11:T=1],Euglena_G002:-424.594[&&NHX:S=Euglena]):368.974[&&NHX:S=Excavata]):0[&&NHX:S=Discoba],(((((Naeglaria_G010:0[&&NHX:S=Naeglaria],Acanthamoeba_G008:0[&&NHX:S=Acanthamoeba:D=0:T=-1]):0[&&NHX:S=Naeglaria:D=11:T=1],Naeglaria_G004:0[&&NHX:S=Naeglaria]):-9778.15[&&NHX:S=Naeglaria:D=1],(Euglena_G004:0[&&NHX:S=Euglena],Euglena_G005:0[&&NHX:S=Euglena]):-9778.15[&&NHX:S=Euglena:D=1]):0[&&NHX:S=Excavata],Physcomistrella_G007:0[&&NHX:S=Physcomistrella:D=0:T=-1]):368.974[&&NHX:S=Excavata:D=11:T=1]):0[&&NHX:S=Discoba]):-2783.1[&&NHX:S=Discoba:D=1],((Saprolegnia_G001:-2739.15[&&NHX:S=Saprolegnia],Plasmodiophora_G001:-2739.15[&&NHX:S=Plasmodiophora]):-1235.74[&&NHX:S=SAR],(Klebsormidium_G001:-3965.52[&&NHX:S=Klebsormidium],Physcomistrella_G001:-3965.52[&&NHX:S=Physcomistrella]):-1235.74[&&NHX:S=Archaeplastida],(Emiliana_G001:-736.96[&&NHX:S=Emiliana],((Chrysochromulina_G009:25.8763[&&NHX:S=Chrysochromulina],Chrysochromulina_G010:25.8763[&&NHX:S=Chrysochromulina]):0[&&NHX:S=Chrysochromulina:D=1],Chrysochromulina_G008:0[&&NHX:S=Chrysochromulina]):-736.96[&&NHX:S=Chrysochromulina:D=1]):-1235.74[&&NHX:S=Haptophyta]):-2783.1[&&NHX:S=Diaphoretickes]):-2050.28[&&NHX:S=Diphoda])r:1[&&NHX:S=Eukaryota];
(((Acanthamoeba_G001:-2310.38[&&NHX:S=Acanthamoeba],Acytostelium_lost:-2310.38[&&NHX:S=Acytostelium:D=1]):-636[&&NHX:S=Amoebozoa],(((Basidiobolus_G003:-304.821[&&NHX:S=Basidiobolus],(Nematostella_G007:0[&&NHX:S=Nematostella],Nematostella_G008:0[&&NHX:S=Nematostella]):-304.821[&&NHX:S=Nematostella:D=1]):56.3398[&&NHX:S=Opisthokonta],(Thecamonas_G004:-2468.43[&&NHX:S=Thecamonas]):56.3398[&&NHX:S=Apusozoa]):0[&&NHX:S=Obazoa],((Basidiobolus_lost:-2637.18[&&NHX:S=Basidiobolus:D=1],Nematostella_G005:-2637.18[&&NHX:S=Nematostella]):56.3398[&&NHX:S=Opisthokonta],(Thecamonas_G005:-4231.96[&&NHX:S=Thecamonas]):56.3398[&&NHX:S=Apusozoa]):0[&&NHX:S=Obazoa]):-636[&&NHX:S=Obazoa:D=1]):-248.737[&&NHX:S=Amorphea],((((((Naeglaria_G011:389.262[&&NHX:S=Naeglaria],Naeglaria_G012:389.262[&&NHX:S=Naeglaria]):0[&&NHX:S=Naeglaria:D=1],Naeglaria_G005:0[&&NHX:S=Naeglaria]):-354.609[&&NHX:S=Naeglaria:D=1],(Euglena_G008:0[&&NHX:S=Euglena],Acanthamoeba_G011:0[&&NHX:S=Acanthamoeba:D=0:T=-1]):-354.609[&&NHX:S=Euglena:D=11:T=1]):0[&&NHX:S=Excavata],((Naeglaria_G006:0[&&NHX:S=Naeglaria],Naeglaria_G007:0[&&NHX:S=Naeglaria]):-1243.7[&&NHX:S=Naeglaria:D=1],Euglena_G003:-1243.7[&&NHX:S=Euglena]):0[&&NHX:S=Excavata]):-1250.8[&&NHX:S=Excavata:D=1]):-819.691[&&NHX:S=Discoba],((((Saprolegnia_G012:0[&&NHX:S=Saprolegnia],(Saprolegnia_G015:406.902[&&NHX:S=Saprolegnia],Saprolegnia_G016:406.902[&&NHX:S=Saprolegnia]):0[&&NHX:S=Saprolegnia:D=1]):-1066.48[&&NHX:S=Saprolegnia:D=1],((Plasmodiophora_G017:205.906[&&NHX:S=Plasmodiophora],Plasmodiophora_G018:205.906[&&NHX:S=Plasmodiophora]):0[&&NHX:S=Plasmodiophora:D=1],Plasmodiophora_G011:0[&&NHX:S=Plasmodiophora]):-1066.48[&&NHX:S=Plasmodiophora:D=1]):0[&&NHX:S=SAR],((((Saprolegnia_G007:-1014.71[&&NHX:S=Saprolegnia],Plasmodiophora_G007:-1014.71[&&NHX:S=Plasmodiophora]):34.8277[&&NHX:S=SAR],(Saprolegnia_G008:-559.823[&&NHX:S=Saprolegnia],(Plasmodiophora_G012:0[&&NHX:S=Plasmodiophora],Plasmodiophora_G013:0[&&NHX:S=Plasmodiophora]):-559.823[&&NHX:S=Plasmodiophora:D=1]):34.8277[&&NHX:S=SAR]):1.89771[&&NHX:S=SAR:D=1],(Saprolegnia_G006:-437.256[&&NHX:S=Saprolegnia],Plasmodiophora_lost:-437.256[&&NHX:S=Plasmodiophora:D=1]):1.89771[&&NHX:S=SAR]):120.717[&&NHX:S=SAR:D=1],(Saprolegnia_G005:-2751.07[&&NHX:S=Saprolegnia],Plasmodiophora_lost:-2751.07[&&NHX:S=Plasmodiophora:D=1]):120.717[&&NHX:S=SAR]):0[&&NHX:S=SAR:D=1]):-2074.79[&&NHX:S=SAR:D=1],(((Klebsormidium_G010:0[&&NHX:S=Klebsormidium],Klebsormidium_G011:0[&&NHX:S=Klebsormidium]):-6794.25[&&NHX:S=Klebsormidium:D=1],Physcomistrella_G006:-6794.25[&&NHX:S=Physcomistrella]):0[&&NHX:S=Archaeplastida],(Klebsormidium_G006:-1321.5[&&NHX:S=Klebsormidium],((Physcomistrella_G013:0[&&NHX:S=Physcomistrella],Physcomistrella_G014:0[&&NHX:S=Physcomistrella]):0[&&NHX:S=Physcomistrella:D=1],(Physcomistrella_G015:0[&&NHX:S=Physcomistrella],Nematostella_G009:0[&&NHX:S=Nematostella:D=0:T=-1]):0[&&NHX:S=Physcomistrella:D=11:T=1]):-1321.5[&&NHX:S=Physcomistrella:D=1]):0[&&NHX:S=Archaeplastida]):-2074.79[&&NHX:S=Archaeplastida:D=1],(((Emiliana_G008:0[&&NHX:S=Emiliana],Plasmodiophora_G014:0[&&NHX:S=Plasmodiophora:D=0:T=-1]):-884.366[&&NHX:S=Emiliana:D=11:T=1],Chrysochromulina_G003:-884.366[&&NHX:S=Chrysochromulina]):0[&&NHX:S=Haptophyta],(Saprolegnia_G004:-2612.74[&&NHX:S=Saprolegnia],((Plasmodiophora_G015:290.138[&&NHX:S=Plasmodiophora],Plasmodiophora_G016:290.138[&&NHX:S=Plasmodiophora]):0[&&NHX:S=Plasmodiophora:D=1],Plasmodiophora_G010:0[&&NHX:S=Plasmodiophora]):-2612.74[&&NHX:S=Plasmodiophora:D=1]):0[&&NHX:S=SAR:D=0:T=-1]):-2074.79[&&NHX:S=Haptophyta:D=11:T=1]):-819.691[&&NHX:S=Diaphoretickes]):-248.737[&&NHX:S=Diphoda])r:1[&&NHX:S=Eukaryota];
((((Acanthamoeba_G003:-2751.5[&&NHX:S=Acanthamoeba],Acytostelium_G001:-2751.5[&&NHX:S=Acytostelium]):383.962[&&NHX:S=Amoebozoa],((Basidiobolus_G001:-1168.45[&&NHX:S=Basidiobolus],Nematostella_G003:-1168.45[&&NHX:S=Nematostella]):-1546.62[&&NHX:S=Opisthokonta],Apusozoa_lost:-1546.62[&&NHX:S=Apusozoa:D=1]):383.962[&&NHX:S=Obazoa]):0[&&NHX:S=Amorphea],((((Acanthamoeba_G012:0[&&NHX:S=Acanthamoeba],Acanthamoeba_G013:0[&&NHX:S=Acanthamoeba]):-9401.43[&&NHX:S=Acanthamoeba:D=1],Acytostelium_G004:-9401.43[&&NHX:S=Acytostelium]):321.812[&&NHX:S=Amoebozoa],((Basidiobolus_G002:-392.95[&&NHX:S=Basidiobolus],Nematostella_lost:-392.95[&&NHX:S=Nematostella:D=1]):-1726.54[&&NHX:S=Opisthokonta],(Thecamonas_G003:-2924.29[&&NHX:S=Thecamonas]):-1726.54[&&NHX:S=Apusozoa]):321.812[&&NHX:S=Obazoa]):62.15[&&NHX:S=Amorphea],((((Saprolegnia_G013:0[&&NHX:S=Saprolegnia],Acanthamoeba_G010:0[&&NHX:S=Acanthamoeba:D=0:T=-1]):0[&&NHX:S=Saprolegnia:D=11:T=1],(Saprolegnia_G014:0[&&NHX:S=Saprolegnia],Saprolegnia_lost:0[&&NHX:S=Saprolegnia:D=1]):0[&&NHX:S=Saprolegnia:D=1]):-122.056[&&NHX:S=Saprolegnia:D=1],Plasmodiophora_G006:-122.056[&&NHX:S=Plasmodiophora]):321.812[&&NHX:S=SAR],(((Klebsormidium_G012:0[&&NHX:S=Klebsormidium],Klebsormidium_G013:0[&&NHX:S=Klebsormidium]):0[&&NHX:S=Klebsormidium:D=1],(Klebsormidium_G014:0[&&NHX:S=Klebsormidium],Thecamonas_G011:0[&&NHX:S=Thecamonas:D=0:T=-1]):0[&&NHX:S=Klebsormidium:D=11:T=1]):-1229.44[&&NHX:S=Klebsormidium:D=1],Physcomistrella_G005:-1229.44[&&NHX:S=Physcomistrella]):321.812[&&NHX:S=Archaeplastida],((Emiliana_G004:-536.633[&&NHX:S=Emiliana],Chrysochromulina_G004:-536.633[&&NHX:S=Chrysochromulina]):0[&&NHX:S=Haptophyta],((Emiliana_G005:-727.526[&&NHX:S=Emiliana],Chrysochromulina_G006:-727.526[&&NHX:S=Chrysochromulina]):0[&&NHX:S=Haptophyta],((Emiliana_G006:187.99[&&NHX:S=Emiliana],Chrysochromulina_G007:187.99[&&NHX:S=Chrysochromulina]):0[&&NHX:S=Haptophyta],Physcomistrella_lost:0[&&NHX:S=Physcomistrella:D=1:T=-1]):0[&&NHX:S=Haptophyta:D=11:T=1]):0[&&NHX:S=Haptophyta:D=1]):321.812[&&NHX:S=Haptophyta:D=1]):62.15[&&NHX:S=Diaphoretickes:D=0:T=-1]):0[&&NHX:S=Amorphea:D=11:T=1]):-1708.2[&&NHX:S=Amorphea:D=1],(((Naeglaria_G001:-594.644[&&NHX:S=Naeglaria],(Euglena_G006:0[&&NHX:S=Euglena],Euglena_G007:0[&&NHX:S=Euglena]):-594.644[&&NHX:S=Euglena:D=1]):-1414.13[&&NHX:S=Excavata]):-4035.09[&&NHX:S=Discoba],((Saprolegnia_lost:-4053.54[&&NHX:S=Saprolegnia:D=1],(Plasmodiophora_G008:0[&&NHX:S=Plasmodiophora],Plasmodiophora_G009:0[&&NHX:S=Plasmodiophora]):-4053.54[&&NHX:S=Plasmodiophora:D=1]):-2445.94[&&NHX:S=SAR],(((Klebsormidium_G008:0[&&NHX:S=Klebsormidium],Klebsormidium_G009:0[&&NHX:S=Klebsormidium]):-5999.98[&&NHX:S=Klebsormidium:D=1],(Physcomistrella_G010:0[&&NHX:S=Physcomistrella],Physcomistrella_G011:0[&&NHX:S=Physcomistrella]):-5999.98[&&NHX:S=Physcomistrella:D=1]):0[&&NHX:S=Archaeplastida],(Klebsormidium_G007:-2378.19[&&NHX:S=Klebsormidium],(Physcomistrella_G012:0[&&NHX:S=Physcomistrella],Euglena_G009:0[&&NHX:S=Euglena:D=0:T=-1]):-2378.19[&&NHX:S=Physcomistrella:D=11:T=1]):0[&&NHX:S=Archaeplastida]):-2445.94[&&NHX:S=Archaeplastida:D=1],((Emiliana_lost:0[&&NHX:S=Emiliana:D=1],Emiliana_G007:0[&&NHX:S=Emiliana]):-562.308[&&NHX:S=Emiliana:D=1],Chrysochromulina_lost:-562.308[&&NHX:S=Chrysochromulina:D=1]):-2445.94[&&NHX:S=Haptophyta]):-4035.09[&&NHX:S=Diaphoretickes]):-1708.2[&&NHX:S=Diphoda])r:1[&&NHX:S=Eukaryota];
((((Acanthamoeba_G006:-2538.14[&&NHX:S=Acanthamoeba],Acytostelium_G005:-2538.14[&&NHX:S=Acytostelium]):0[&&NHX:S=Amoebozoa],(Acanthamoeba_G007:-3194.05[&&NHX:S=Acanthamoeba],((Acytostelium_G008:170.975[&&NHX:S=Acytostelium],(Acytostelium_G009:214.885[&&NHX:S=Acytostelium],Acytostelium_G010:214.885[&&NHX:S=Acytostelium]):170.975[&&NHX:S=Acytostelium:D=1]):0[&&NHX:S=Acytostelium:D=1],(Acanthamoeba_G014:102.148[&&NHX:S=Acanthamoeba],Acanthamoeba_G015:102.148[&&NHX:S=Acanthamoeba]):0[&&NHX:S=Acanthamoeba:D=1:T=-1]):-3194.05[&&NHX:S=Acytostelium:D=11:T=1]):0[&&NHX:S=Amoebozoa]):-4037.72[&&NHX:S=Amoebozoa:D=1],(((((Basidiobolus_G013:353.823[&&NHX:S=Basidiobolus],Basidiobolus_G014:353.823[&&NHX:S=Basidiobolus]):0[&&NHX:S=Basidiobolus:D=1],Basidiobolus_G008:0[&&NHX:S=Basidiobolus]):0[&&NHX:S=Basidiobolus:D=1],Acytostelium_G007:0[&&NHX:S=Acytostelium:D=0:T=-1]):-2118.35[&&NHX:S=Basidiobolus:D=11:T=1],Nematostella_G001:-2118.35[&&NHX:S=Nematostella]):-1971.16[&&NHX:S=Opisthokonta],(Thecamonas_G001:-2435.61[&&NHX:S=Thecamonas]):-1971.16[&&NHX:S=Apusozoa]):-4037.72[&&NHX:S=Obazoa]):-2434.85[&&NHX:S=Amorphea],(((((Naeglaria_G008:0[&&NHX:S=Naeglaria],Naeglaria_G009:0[&&NHX:S=Naeglaria]):0[&&NHX:S=Naeglaria:D=1],Naeglaria_G002:0[&&NHX:S=Naeglaria]):-894.161[&&NHX:S=Naeglaria:D=1],Euglena_G001:-894.161[&&NHX:S=Euglena]):-2019.14[&&NHX:S=Excavata]):-3684.33[&&NHX:S=Discoba],(((Saprolegnia_G010:0[&&NHX:S=Saprolegnia],Saprolegnia_G011:0[&&NHX:S=Saprolegnia]):-439.354[&&NHX:S=Saprolegnia:D=1],Plasmodiophora_G002:-439.354[&&NHX:S=Plasmodiophora]):-5760.87[&&NHX:S=SAR],(Klebsormidium_G002:-744.237[&&NHX:S=Klebsormidium],Physcomistrella_G002:-744.237[&&NHX:S=Physcomistrella]):-5760.87[&&NHX:S=Archaeplastida],(Emiliana_G002:-2037.81[&&NHX:S=Emiliana],Chrysochromulina_G001:-2037.81[&&NHX:S=Chrysochromulina]):-5760.87[&&NHX:S=Haptophyta]):-3684.33[&&NHX:S=Diaphoretickes]):-2434.85[&&NHX:S=Diphoda])r:1[&&NHX:S=Eukaryota];
(((Acanthamoeba_G002:-1425.3[&&NHX:S=Acanthamoeba],Acytostelium_lost:-1425.3[&&NHX:S=Acytostelium:D=1]):-2522.19[&&NHX:S=Amoebozoa],(((Basidiobolus_G004:0[&&NHX:S=Basidiobolus],Basidiobolus_G005:0[&&NHX:S=Basidiobolus]):-1029.4[&&NHX:S=Basidiobolus:D=1],Nematostella_G002:-1029.4[&&NHX:S=Nematostella]):-496.731[&&NHX:S=Opisthokonta],((Thecamonas_G006:-1234.42[&&NHX:S=Thecamonas]):0[&&NHX:S=Apusozoa],(Thecamonas_lost:-3270.23[&&NHX:S=Thecamonas:D=1]):0[&&NHX:S=Apusozoa]):-496.731[&&NHX:S=Apusozoa:D=1]):-2522.19[&&NHX:S=Obazoa]):-3672.96[&&NHX:S=Amorphea],(Discoba_lost:-2187.81[&&NHX:S=Discoba:D=1],(((Saprolegnia_G009:0[&&NHX:S=Saprolegnia],((Thecamonas_G009:0[&&NHX:S=Thecamonas],Thecamonas_G010:0[&&NHX:S=Thecamonas]):27.3765[&&NHX:S=Thecamonas:D=1]):0[&&NHX:S=Apusozoa:D=0:T=-1]):-1452.5[&&NHX:S=Saprolegnia:D=11:T=1],Plasmodiophora_G003:-1452.5[&&NHX:S=Plasmodiophora]):-852.386[&&NHX:S=SAR],(Klebsormidium_G003:-761.738[&&NHX:S=Klebsormidium],Physcomistrella_G003:-761.738[&&NHX:S=Physcomistrella]):-852.386[&&NHX:S=Archaeplastida],(((Emiliana_G009:0[&&NHX:S=Emiliana],Emiliana_G010:0[&&NHX:S=Emiliana]):-1700.99[&&NHX:S=Emiliana:D=1],Chrysochromulina_G005:-1700.99[&&NHX:S=Chrysochromulina]):0[&&NHX:S=Haptophyta],(Nematostella_G006:0[&&NHX:S=Nematostella],Acytostelium_G006:0[&&NHX:S=Acytostelium:D=0:T=-1]):0[&&NHX:S=Nematostella:D=11:T=1]):-852.386[&&NHX:S=Haptophyta:D=11:T=1]):-2187.81[&&NHX:S=Diaphoretickes]):-3672.96[&&NHX:S=Diphoda])r:1[&&NHX:S=Eukaryota];""")

SIM50_TRANSFERS = newtestfile("sim50_t.nwk",
        """(Amorphea_G001:2000[&&NHX:S=Amorphea],Diphoda_G001:2000[&&NHX:S=Diphoda])r:100[&&NHX:S=Eukaryota];
(((Amorphea_G051:393.063[&&NHX:S=Amorphea],(Amorphea_G068:393.063[&&NHX:S=Amorphea],Amorphea_G069:393.063[&&NHX:S=Amorphea]):0[&&NHX:S=Amorphea:D=1]):1606.94[&&NHX:S=Amorphea:D=1],Diphoda_G018:2000[&&NHX:S=Diphoda]):33.7677[&&NHX:S=Eukaryota],((Amorphea_G052:1690.98[&&NHX:S=Amorphea],((Amorphea_G080:1237.25[&&NHX:S=Amorphea],Amorphea_G081:1237.25[&&NHX:S=Amorphea]):453.738[&&NHX:S=Amorphea:D=1],(Diphoda_G080:1237.25[&&NHX:S=Diphoda],Diphoda_G081:1237.25[&&NHX:S=Diphoda]):453.738[&&NHX:S=Diphoda:D=1:T=-1]):0[&&NHX:S=Amorphea:D=11:T=1]):309.017[&&NHX:S=Amorphea:D=1],Diphoda_G019:2000[&&NHX:S=Diphoda]):33.7677[&&NHX:S=Eukaryota])r:66.2323[&&NHX:S=Eukaryota:D=1];
((Amorphea_G023:1556.94[&&NHX:S=Amorphea],Amorphea_G024:1556.94[&&NHX:S=Amorphea]):443.065[&&NHX:S=Amorphea:D=1],Diphoda_lost:1509.28[&&NHX:S=Diphoda:D=1])r:100[&&NHX:S=Eukaryota];
(Amorphea_G002:2000[&&NHX:S=Amorphea],Diphoda_G002:2000[&&NHX:S=Diphoda])r:100[&&NHX:S=Eukaryota];
((Amorphea_G025:988.389[&&NHX:S=Amorphea],Amorphea_G026:988.389[&&NHX:S=Amorphea]):1011.61[&&NHX:S=Amorphea:D=1],(Diphoda_G021:911.604[&&NHX:S=Diphoda],Diphoda_G022:911.604[&&NHX:S=Diphoda]):1088.4[&&NHX:S=Diphoda:D=1])r:100[&&NHX:S=Eukaryota];
(Amorphea_G003:2000[&&NHX:S=Amorphea],((Diphoda_G046:1339.68[&&NHX:S=Diphoda],Amorphea_G055:1339.68[&&NHX:S=Amorphea:D=0:T=-1]):0[&&NHX:S=Diphoda:D=11:T=1],Diphoda_G023:1339.68[&&NHX:S=Diphoda]):660.32[&&NHX:S=Diphoda:D=1])r:100[&&NHX:S=Eukaryota];
(Amorphea_G004:2000[&&NHX:S=Amorphea],Diphoda_G003:2000[&&NHX:S=Diphoda])r:100[&&NHX:S=Eukaryota];
(Amorphea_lost:420.564[&&NHX:S=Amorphea:D=1],((Diphoda_G047:768.248[&&NHX:S=Diphoda],Amorphea_G056:768.248[&&NHX:S=Amorphea:D=0:T=-1]):0[&&NHX:S=Diphoda:D=11:T=1],(Diphoda_G048:768.248[&&NHX:S=Diphoda],Amorphea_G057:768.248[&&NHX:S=Amorphea:D=0:T=-1]):0[&&NHX:S=Diphoda:D=11:T=1]):1231.75[&&NHX:S=Diphoda:D=1])r:100[&&NHX:S=Eukaryota];
(Amorphea_G005:2000[&&NHX:S=Amorphea],(Diphoda_G024:1474.08[&&NHX:S=Diphoda],Amorphea_G020:1474.08[&&NHX:S=Amorphea:D=0:T=-1]):525.917[&&NHX:S=Diphoda:D=11:T=1])r:100[&&NHX:S=Eukaryota];
((Amorphea_G027:1010.24[&&NHX:S=Amorphea],Amorphea_G028:1010.24[&&NHX:S=Amorphea]):989.755[&&NHX:S=Amorphea:D=1],Diphoda_G004:2000[&&NHX:S=Diphoda])r:100[&&NHX:S=Eukaryota];
(Amorphea_lost:1120[&&NHX:S=Amorphea:D=1],Diphoda_G005:2000[&&NHX:S=Diphoda])r:100[&&NHX:S=Eukaryota];
((Amorphea_G029:91.1246[&&NHX:S=Amorphea],Amorphea_G030:91.1246[&&NHX:S=Amorphea]):1908.88[&&NHX:S=Amorphea:D=1],Diphoda_G006:2000[&&NHX:S=Diphoda])r:100[&&NHX:S=Eukaryota];
((Amorphea_G031:829.245[&&NHX:S=Amorphea],(Diphoda_G070:829.245[&&NHX:S=Diphoda],((Diphoda_G084:425.331[&&NHX:S=Diphoda],(Diphoda_G090:153.764[&&NHX:S=Diphoda],Diphoda_G091:153.764[&&NHX:S=Diphoda]):271.567[&&NHX:S=Diphoda:D=1]):5.60178[&&NHX:S=Diphoda:D=1],(Diphoda_G085:425.331[&&NHX:S=Diphoda],Diphoda_G086:425.331[&&NHX:S=Diphoda]):5.60178[&&NHX:S=Diphoda:D=1]):398.313[&&NHX:S=Diphoda:D=1]):0[&&NHX:S=Diphoda:D=1:T=-1]):1170.75[&&NHX:S=Amorphea:D=11:T=1],(Diphoda_G025:618.483[&&NHX:S=Diphoda],Amorphea_G021:618.483[&&NHX:S=Amorphea:D=0:T=-1]):1381.52[&&NHX:S=Diphoda:D=11:T=1])r:100[&&NHX:S=Eukaryota];
((Amorphea_G032:846.585[&&NHX:S=Amorphea],Amorphea_G033:846.585[&&NHX:S=Amorphea]):1153.41[&&NHX:S=Amorphea:D=1],Diphoda_G007:2000[&&NHX:S=Diphoda])r:100[&&NHX:S=Eukaryota];
((Amorphea_G034:259.46[&&NHX:S=Amorphea],Amorphea_G035:259.46[&&NHX:S=Amorphea]):1740.54[&&NHX:S=Amorphea:D=1],(Diphoda_G026:1366.97[&&NHX:S=Diphoda],Diphoda_G027:1366.97[&&NHX:S=Diphoda]):633.032[&&NHX:S=Diphoda:D=1])r:100[&&NHX:S=Eukaryota];
(Amorphea_G006:2000[&&NHX:S=Amorphea],(Diphoda_G028:1648.88[&&NHX:S=Diphoda],(Diphoda_G049:1648.88[&&NHX:S=Diphoda],Diphoda_G050:1648.88[&&NHX:S=Diphoda]):0[&&NHX:S=Diphoda:D=1]):351.123[&&NHX:S=Diphoda:D=1])r:100[&&NHX:S=Eukaryota];
(Amorphea_G007:2000[&&NHX:S=Amorphea],(Diphoda_G029:1100.82[&&NHX:S=Diphoda],((Amorphea_G072:992.398[&&NHX:S=Amorphea],Diphoda_G074:992.398[&&NHX:S=Diphoda:D=0:T=-1]):108.425[&&NHX:S=Amorphea:D=11:T=1],(Diphoda_G075:992.398[&&NHX:S=Diphoda],(Amorphea_G082:482.39[&&NHX:S=Amorphea],Amorphea_G083:482.39[&&NHX:S=Amorphea]):510.008[&&NHX:S=Amorphea:D=1:T=-1]):108.425[&&NHX:S=Diphoda:D=11:T=1]):0[&&NHX:S=Amorphea:D=11:T=1]):899.177[&&NHX:S=Diphoda:D=11:T=1])r:100[&&NHX:S=Eukaryota];
(((Amorphea_G060:1264[&&NHX:S=Amorphea],Amorphea_G061:1264[&&NHX:S=Amorphea]):0[&&NHX:S=Amorphea:D=1],Diphoda_lost:0[&&NHX:S=Diphoda:D=1:T=-1]):736.003[&&NHX:S=Amorphea:D=11:T=1],((Diphoda_G051:1830.1[&&NHX:S=Diphoda],(Diphoda_G076:720.717[&&NHX:S=Diphoda],Diphoda_G077:720.717[&&NHX:S=Diphoda]):1109.38[&&NHX:S=Diphoda:D=1]):0[&&NHX:S=Diphoda:D=1],Diphoda_G030:1830.1[&&NHX:S=Diphoda]):169.899[&&NHX:S=Diphoda:D=1])r:100[&&NHX:S=Eukaryota];
((Amorphea_G036:756.826[&&NHX:S=Amorphea],Amorphea_G037:756.826[&&NHX:S=Amorphea]):1243.17[&&NHX:S=Amorphea:D=1],Diphoda_lost:495.971[&&NHX:S=Diphoda:D=1])r:100[&&NHX:S=Eukaryota];
(Amorphea_G008:2000[&&NHX:S=Amorphea],Diphoda_lost:481.268[&&NHX:S=Diphoda:D=1])r:100[&&NHX:S=Eukaryota];
r:21.4943[&&NHX:S=Eukaryota:D=1];
((Amorphea_G038:1959.16[&&NHX:S=Amorphea],(Amorphea_G062:1959.16[&&NHX:S=Amorphea],(Amorphea_G077:1306.68[&&NHX:S=Amorphea],(Amorphea_G086:966.673[&&NHX:S=Amorphea],Amorphea_G087:966.673[&&NHX:S=Amorphea]):340.003[&&NHX:S=Amorphea:D=1]):652.487[&&NHX:S=Amorphea:D=1]):0[&&NHX:S=Amorphea:D=1]):40.8372[&&NHX:S=Amorphea:D=1],Diphoda_G008:2000[&&NHX:S=Diphoda])r:100[&&NHX:S=Eukaryota];
(Amorphea_G009:2000[&&NHX:S=Amorphea],((((Diphoda_G082:1399.96[&&NHX:S=Diphoda],Diphoda_G083:1399.96[&&NHX:S=Diphoda]):76.0948[&&NHX:S=Diphoda:D=1],((Diphoda_G087:597.798[&&NHX:S=Diphoda],Diphoda_G088:597.798[&&NHX:S=Diphoda]):802.162[&&NHX:S=Diphoda:D=1],((Diphoda_G092:200.929[&&NHX:S=Diphoda],Amorphea_G088:200.929[&&NHX:S=Amorphea:D=0:T=-1]):396.869[&&NHX:S=Diphoda:D=11:T=1],Diphoda_G089:597.798[&&NHX:S=Diphoda]):802.162[&&NHX:S=Diphoda:D=1]):76.0948[&&NHX:S=Diphoda:D=1]):314.2[&&NHX:S=Diphoda:D=1],((Amorphea_G084:523.711[&&NHX:S=Amorphea],Amorphea_G085:523.711[&&NHX:S=Amorphea]):952.344[&&NHX:S=Amorphea:D=1],Amorphea_G073:1476.06[&&NHX:S=Amorphea]):314.2[&&NHX:S=Amorphea:D=1:T=-1]):0[&&NHX:S=Diphoda:D=11:T=1],(Diphoda_G052:1790.26[&&NHX:S=Diphoda],Amorphea_G058:1790.26[&&NHX:S=Amorphea:D=0:T=-1]):0[&&NHX:S=Diphoda:D=11:T=1]):209.745[&&NHX:S=Diphoda:D=1])r:100[&&NHX:S=Eukaryota];
(Amorphea_G010:2000[&&NHX:S=Amorphea],(Diphoda_lost:0[&&NHX:S=Diphoda:D=1],(Diphoda_G053:1914.99[&&NHX:S=Diphoda],Diphoda_G054:1914.99[&&NHX:S=Diphoda]):0[&&NHX:S=Diphoda:D=1]):85.0131[&&NHX:S=Diphoda:D=1])r:100[&&NHX:S=Eukaryota];
(Amorphea_lost:804.071[&&NHX:S=Amorphea:D=1],(Diphoda_G031:1140.44[&&NHX:S=Diphoda],Diphoda_G032:1140.44[&&NHX:S=Diphoda]):859.559[&&NHX:S=Diphoda:D=1])r:100[&&NHX:S=Eukaryota];
(Amorphea_G011:2000[&&NHX:S=Amorphea],Diphoda_G009:2000[&&NHX:S=Diphoda])r:100[&&NHX:S=Eukaryota];
(Amorphea_G012:2000[&&NHX:S=Amorphea],((Diphoda_lost:192.883[&&NHX:S=Diphoda:D=1],Diphoda_G055:1442.36[&&NHX:S=Diphoda]):0[&&NHX:S=Diphoda:D=1],(Diphoda_G056:1442.36[&&NHX:S=Diphoda],Diphoda_G057:1442.36[&&NHX:S=Diphoda]):0[&&NHX:S=Diphoda:D=1]):557.639[&&NHX:S=Diphoda:D=1])r:100[&&NHX:S=Eukaryota];
(Amorphea_G013:2000[&&NHX:S=Amorphea],Diphoda_G010:2000[&&NHX:S=Diphoda])r:100[&&NHX:S=Eukaryota];
(Amorphea_lost:467.152[&&NHX:S=Amorphea:D=1],(Diphoda_G033:11.0533[&&NHX:S=Diphoda],Diphoda_G034:11.0533[&&NHX:S=Diphoda]):1988.95[&&NHX:S=Diphoda:D=1])r:100[&&NHX:S=Eukaryota];
((Amorphea_G039:1381.5[&&NHX:S=Amorphea],Amorphea_G040:1381.5[&&NHX:S=Amorphea]):618.495[&&NHX:S=Amorphea:D=1],Diphoda_G011:2000[&&NHX:S=Diphoda])r:100[&&NHX:S=Eukaryota];
((Amorphea_G041:1750.84[&&NHX:S=Amorphea],(Amorphea_G063:1750.84[&&NHX:S=Amorphea],(Amorphea_G078:1650.96[&&NHX:S=Amorphea],Diphoda_G079:1650.96[&&NHX:S=Diphoda:D=0:T=-1]):99.8736[&&NHX:S=Amorphea:D=11:T=1]):0[&&NHX:S=Amorphea:D=1]):249.164[&&NHX:S=Amorphea:D=1],(Diphoda_G035:1658.57[&&NHX:S=Diphoda],(Diphoda_G058:1658.57[&&NHX:S=Diphoda],Diphoda_G059:1658.57[&&NHX:S=Diphoda]):0[&&NHX:S=Diphoda:D=1]):341.427[&&NHX:S=Diphoda:D=1])r:100[&&NHX:S=Eukaryota];
((Amorphea_G042:582.038[&&NHX:S=Amorphea],(Amorphea_G064:582.038[&&NHX:S=Amorphea],Diphoda_G071:582.038[&&NHX:S=Diphoda:D=0:T=-1]):0[&&NHX:S=Amorphea:D=11:T=1]):1417.96[&&NHX:S=Amorphea:D=1],((Diphoda_lost:953.339[&&NHX:S=Diphoda:D=1],Diphoda_G060:1311.37[&&NHX:S=Diphoda]):0[&&NHX:S=Diphoda:D=1],Diphoda_G036:1311.37[&&NHX:S=Diphoda]):688.631[&&NHX:S=Diphoda:D=1])r:100[&&NHX:S=Eukaryota];
((Amorphea_G043:1421.5[&&NHX:S=Amorphea],Amorphea_G044:1421.5[&&NHX:S=Amorphea]):578.501[&&NHX:S=Amorphea:D=1],Diphoda_lost:1999.87[&&NHX:S=Diphoda:D=1])r:100[&&NHX:S=Eukaryota];
((((Amorphea_G070:1341.29[&&NHX:S=Amorphea],Amorphea_G071:1341.29[&&NHX:S=Amorphea]):0[&&NHX:S=Amorphea:D=1],Amorphea_G053:1341.29[&&NHX:S=Amorphea]):658.714[&&NHX:S=Amorphea:D=1],(Diphoda_G042:1250.45[&&NHX:S=Diphoda],Diphoda_G043:1250.45[&&NHX:S=Diphoda]):749.551[&&NHX:S=Diphoda:D=1]):2.04415[&&NHX:S=Eukaryota],((Amorphea_G054:1336.48[&&NHX:S=Amorphea],Diphoda_G045:1336.48[&&NHX:S=Diphoda:D=0:T=-1]):663.52[&&NHX:S=Amorphea:D=11:T=1],Diphoda_G020:2000[&&NHX:S=Diphoda]):2.04415[&&NHX:S=Eukaryota])r:97.9559[&&NHX:S=Eukaryota:D=1];
(Amorphea_lost:372.356[&&NHX:S=Amorphea:D=1],Diphoda_G012:2000[&&NHX:S=Diphoda])r:100[&&NHX:S=Eukaryota];
(Amorphea_G014:2000[&&NHX:S=Amorphea],(Diphoda_G037:441.792[&&NHX:S=Diphoda],(Diphoda_G061:441.792[&&NHX:S=Diphoda],(Amorphea_G074:425.551[&&NHX:S=Amorphea],Amorphea_G075:425.551[&&NHX:S=Amorphea]):16.2415[&&NHX:S=Amorphea:D=1:T=-1]):0[&&NHX:S=Diphoda:D=11:T=1]):1558.21[&&NHX:S=Diphoda:D=1])r:100[&&NHX:S=Eukaryota];
(((Amorphea_lost:37.759[&&NHX:S=Amorphea:D=1],Amorphea_G065:1018.02[&&NHX:S=Amorphea]):0[&&NHX:S=Amorphea:D=1],Diphoda_lost:0[&&NHX:S=Diphoda:D=1:T=-1]):981.981[&&NHX:S=Amorphea:D=11:T=1],Diphoda_G013:2000[&&NHX:S=Diphoda])r:100[&&NHX:S=Eukaryota];
((Amorphea_G045:1650.7[&&NHX:S=Amorphea],Diphoda_G044:1650.7[&&NHX:S=Diphoda:D=0:T=-1]):349.298[&&NHX:S=Amorphea:D=11:T=1],Diphoda_G014:2000[&&NHX:S=Diphoda])r:100[&&NHX:S=Eukaryota];
(Amorphea_lost:429.545[&&NHX:S=Amorphea:D=1],((Diphoda_G062:784.338[&&NHX:S=Diphoda],Diphoda_G063:784.338[&&NHX:S=Diphoda]):0[&&NHX:S=Diphoda:D=1],Amorphea_G022:784.338[&&NHX:S=Amorphea:D=0:T=-1]):1215.66[&&NHX:S=Diphoda:D=11:T=1])r:100[&&NHX:S=Eukaryota];
(Amorphea_G015:2000[&&NHX:S=Amorphea],(Diphoda_G038:1647.62[&&NHX:S=Diphoda],(Diphoda_G064:1647.62[&&NHX:S=Diphoda],Diphoda_G065:1647.62[&&NHX:S=Diphoda]):0[&&NHX:S=Diphoda:D=1]):352.379[&&NHX:S=Diphoda:D=1])r:100[&&NHX:S=Eukaryota];
(Amorphea_G016:2000[&&NHX:S=Amorphea],Diphoda_lost:12.8015[&&NHX:S=Diphoda:D=1])r:100[&&NHX:S=Eukaryota];
(Amorphea_G017:2000[&&NHX:S=Amorphea],(Diphoda_G039:1921.39[&&NHX:S=Diphoda],(Diphoda_G066:1921.39[&&NHX:S=Diphoda],Amorphea_G059:1921.39[&&NHX:S=Amorphea:D=0:T=-1]):0[&&NHX:S=Diphoda:D=11:T=1]):78.6076[&&NHX:S=Diphoda:D=1])r:100[&&NHX:S=Eukaryota];
(Amorphea_G018:2000[&&NHX:S=Amorphea],Diphoda_G015:2000[&&NHX:S=Diphoda])r:100[&&NHX:S=Eukaryota];
((Amorphea_lost:0[&&NHX:S=Amorphea:D=1],Amorphea_G046:837.128[&&NHX:S=Amorphea]):1162.87[&&NHX:S=Amorphea:D=1],Diphoda_G016:2000[&&NHX:S=Diphoda])r:100[&&NHX:S=Eukaryota];
((Amorphea_lost:0[&&NHX:S=Amorphea:D=1],((Amorphea_lost:57.0726[&&NHX:S=Amorphea:D=1],Amorphea_G079:1345.87[&&NHX:S=Amorphea]):571.05[&&NHX:S=Amorphea:D=1],Amorphea_lost:571.05[&&NHX:S=Amorphea:D=1]):0[&&NHX:S=Amorphea:D=1]):83.0768[&&NHX:S=Amorphea:D=1],Diphoda_lost:804.268[&&NHX:S=Diphoda:D=1])r:100[&&NHX:S=Eukaryota];
((Amorphea_G047:1387.92[&&NHX:S=Amorphea],(Amorphea_G066:1387.92[&&NHX:S=Amorphea],Amorphea_G067:1387.92[&&NHX:S=Amorphea]):0[&&NHX:S=Amorphea:D=1]):612.079[&&NHX:S=Amorphea:D=1],Diphoda_lost:377.153[&&NHX:S=Diphoda:D=1])r:100[&&NHX:S=Eukaryota];
(Amorphea_lost:1697.43[&&NHX:S=Amorphea:D=1],(Diphoda_G040:1180.09[&&NHX:S=Diphoda],Diphoda_G041:1180.09[&&NHX:S=Diphoda]):819.912[&&NHX:S=Diphoda:D=1])r:100[&&NHX:S=Eukaryota];
(Amorphea_G019:2000[&&NHX:S=Amorphea],(((Diphoda_G078:108.464[&&NHX:S=Diphoda],Amorphea_G076:108.464[&&NHX:S=Amorphea:D=0:T=-1]):1857.35[&&NHX:S=Diphoda:D=11:T=1],Diphoda_G067:1965.82[&&NHX:S=Diphoda]):0[&&NHX:S=Diphoda:D=1],(Diphoda_G068:1965.82[&&NHX:S=Diphoda],Diphoda_G069:1965.82[&&NHX:S=Diphoda]):0[&&NHX:S=Diphoda:D=1]):34.1817[&&NHX:S=Diphoda:D=1])r:100[&&NHX:S=Eukaryota];
((Amorphea_G048:1881.67[&&NHX:S=Amorphea],(Diphoda_G072:1881.67[&&NHX:S=Diphoda],Diphoda_G073:1881.67[&&NHX:S=Diphoda]):0[&&NHX:S=Diphoda:D=1:T=-1]):118.329[&&NHX:S=Amorphea:D=11:T=1],(Diphoda_lost:0[&&NHX:S=Diphoda:D=1],Diphoda_lost:0[&&NHX:S=Diphoda:D=1]):483.293[&&NHX:S=Diphoda:D=1])r:100[&&NHX:S=Eukaryota];
((Amorphea_G049:396.868[&&NHX:S=Amorphea],Amorphea_G050:396.868[&&NHX:S=Amorphea]):1603.13[&&NHX:S=Amorphea:D=1],Diphoda_G017:2000[&&NHX:S=Diphoda])r:100[&&NHX:S=Eukaryota];
""")
DESIGNED_FALSE_POSITIVE_EUKA2 = newtestfile("designed_false_positive.nwk",
        "((Amorphea_G001:1000[&&NHX:S=Amorphea],(Diphoda_G002:990[&&NHX:S=Diphoda]):10[&&NHX:S=Diphoda:T=-1:D=11]):1000[&&NHX:S=Amorphea:T=1:D=11],Diphoda_G001:2000[&&NHX:S=Diphoda])r:100[&&NHX:S=Eukaryota];")

# #genetree_drawer.py -p _tests/euka_2supergroups.nwk -t - sim5.nwk
# #genetree_drawer.py -p _tests/euka_2supergroups.nwk -t - <(head -1 sim5.nwk)
# genetree_drawer.py -g leaf -p _tests/euka_2supergroups.nwk -t -  sim5.nwk
# genetree_drawer.py -A -g leaf -p _tests/euka_2supergroups.nwk -t - sim5.nwk
# genetree_drawer.py -a3 -A -g leaf -p _tests/euka_2supergroups.nwk -t - sim5.nwk
# genetree_drawer.py -a2 -A -g leaf -p _tests/euka_2supergroups.nwk -t - sim5.nwk
# genetree_drawer.py -a2 -A -g leaf -p _tests/euka_2supergroups.nwk -t - sim5_t.nwk
# genetree_drawer.py -A -g leaf -p _tests/euka_2supergroups.nwk -t - sim5_10sp.nwk # Should FAIL
# genetree_drawer.py -A -g leaf -p _tests/euka_10species.nwk -t - sim5_10sp.nwk
# genetree_drawer.py -a1 -A -g leaf -p _tests/euka_10species.nwk -t - sim5_10sp.nwk
# genetree_drawer.py -a2 -A -g leaf -p _tests/euka_10species.nwk -t - sim5_10sp.nwk
# genetree_drawer.py -a3 -A -g leaf -p _tests/euka_10species.nwk -t - sim5_10sp.nwk
# genetree_drawer.py -a4 -A -g leaf -p _tests/euka_10species.nwk -t - sim5_10sp.nwk
# genetree_drawer.py -a5 -A -g leaf -p _tests/euka_10species.nwk -t - sim5_10sp.nwk
# genetree_drawer.py -a2 -A -g leaf -p _tests/euka_10species.nwk -t - sim5_10sp_t.nwk
# genetree_drawer.py -k -a2 -A -g leaf -p _tests/euka_10species.nwk -t - sim5_10sp_t.nwk
# genetree_drawer.py -a2 -A -p _tests/euka_2supergroups.nwk -t -g leaf -- - <(fgrep Diphoda_G037 sim50_t.nwk)
# genetree_drawer.py -a2 -A -p _tests/euka_2supergroups.nwk -t -g leaf -T "T" -- - <(fgrep Diphoda_G037 sim50_t.nwk)
# genetree_drawer.py -a2 -A -p _tests/euka_2supergroups.nwk -t -g leaf,dup -T "T" -- - <(fgrep Diphoda_G037 sim50_t.nwk)
# genetree_drawer.py -a2 -k -A -p _tests/euka_2supergroups.nwk -t -g leaf,dup -T "T" -- - <(fgrep Diphoda_G037 sim50_t.nwk)
# genetree_drawer.py -a2 -A -p _tests/euka_2supergroups.nwk -t -g leaf,dup -T "T" -- - <(fgrep Amorphea_G052 sim50_t.nwk)
# genetree_drawer.py -a2 -A -p _tests/euka_2supergroups.nwk -t -g leaf,dup -T "T" -- - <(fgrep Amorphea_G048 sim50_t.nwk)
# genetree_drawer.py -a2 -A -p _tests/euka_2supergroups.nwk -t -g leaf,dup -T "T" -- - _tests/designed_false_positive.nwk
# genetree_drawer.py -a3 -A -p _tests/euka_2supergroups.nwk -t -g leaf,dup -T "T" -- - _tests/designed_false_positive.nwk


def setup_module():
    for testfile in (PHYLTREE_NEWICK, NODAL_NEWICK, LMLN2_NEWICK,
            SIMPLE_PHYLTREE, SIMPLE_TREEBEST, SIMPLE_TREEBEST_T, SIMPLE_TREEBEST_T2,
            ALE2_PHYLTREE, ALE2_REC_TREEBEST, ALE2_REC_TREEBEST_TEST1,
            EUKA_2SUPERGROUPS, EUKA_10SPECIES, SIM5, SIM5_TRANSFERS,
            SIM5_10SP, SIM5_10SP_TRANSFERS, SIM50_TRANSFERS, DESIGNED_FALSE_POSITIVE_EUKA2):
        with open(testfile.name, 'w') as out:
            out.write(testfile.content)


#def teardown_module():
#    shutil.rmtree(TEMPDIRNAME)

class Test_Simple_TreeBest:
    kwargs = dict(outfile='-', angle_style=3, phyltreefile=SIMPLE_PHYLTREE.name, treebest=True)

    def test_onedup(self):
        gd = run([SIMPLE_TREEBEST.name], [], **self.kwargs)
    def test_onetransfer_fails_with_KeyError(self):
        with pytest.raises(KeyError):
            gd = run([SIMPLE_TREEBEST_T.name], [], **self.kwargs)
    def test_onetransfer2(self):
        gd = run([SIMPLE_TREEBEST_T2.name], [], **self.kwargs)
    def test_colorize_clade(self):
        gd = run([SIMPLE_TREEBEST.name], [], colorize_clades=['Pan'], **self.kwargs)



def test_NODAL():
    displayed_internal = ['Amniota','Theria','Euarchontoglires','Neopterygii','Carnivora','Rodentia','Primates','Neognathae']
    colorized_descent = ['AmniotaENSGT00910000143982.A.b']

    gd = run([NODAL_NEWICK.name], [], outfile='-', ensembl_version=93, angle_style=3,
             phyltreefile=PHYLTREE_NEWICK.name,
             colorize_descent=colorized_descent, internal=displayed_internal,
             )


#def test_ALE_newick():
#    gd = run([ALE_FILENAME], [], outfile='-', angle_style=3,
#             phyltreefile=ALE_PHYLTREE_FILENAME,
#            )

class Test_ALE:
    kwargs = dict(outfile='-', angle_style=3, phyltreefile=ALE2_PHYLTREE.name, treebest=True, genenames='dup,transfer', tags='T,id')
    def test_rec_treebest(self):
        gd = run([ALE2_REC_TREEBEST.name], [], colorize_clades=['19'], **self.kwargs)
        # Ok but the line in the species branch leading to SYNPW clade (node 19) should be black if speciation, else labelled as transfer received.
        # Also one transfer source has wrong y position on the branch.
    def test_rec_treebest_test1_fails_KeyError(self):
        with pytest.raises(KeyError):
            gd = run([ALE2_REC_TREEBEST_TEST1.name], [], **self.kwargs)
        # KeyError: "SYNP2@26|ACAM1"

class Test_Euka_Sim:
    kwargs = dict(outfile='-', angle_style=3, treebest=True, genenames='leaf,dup,transfer', tags='T,id')
    kwargs_2sg = dict(phyltreefile=EUKA_2SUPERGROUPS.name, **kwargs)
    kwargs_10sp = dict(phyltreefile=EUKA_10SPECIES.name, **kwargs)

    def test_sim5(self):
        gd = run([SIM5.name], [], **self.kwargs_2sg)

    def test_sim5_transfers(self):
        gd = run([SIM5_TRANSFERS.name], [], **self.kwargs_2sg)
        # The drawn lines are messy (children_ys inverted "ERROR: Children's relative Y not sorted"). FIXME.
        #assert False
        #gd.branchings

    def test_sim5_10sp_wrong_speciestree_fails_AssertionError(self):
        with pytest.raises(AssertionError):
            gd = run([SIM5_10SP.name], [], **self.kwargs_2sg)
        #FIXME: hard to diagnostic error (line: assert expected_children_taxa, children_taxa
        # AssertionError: {'Acanthamoeba', 'Acytostelium'}

    def test_sim5_10sp(self):
        gd = run([SIM5_10SP.name], [], **self.kwargs_10sp)

    def test_sim5_10sp_transfers(self):
        gd = run([SIM5_10SP_TRANSFERS.name], [], **self.kwargs_10sp)

    def test_false_positive(self):
        gd = run([DESIGNED_FALSE_POSITIVE_EUKA2.name], [], **self.kwargs_2sg)
