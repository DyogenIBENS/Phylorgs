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

def setup_module():
    for testfile in (PHYLTREE_NEWICK, NODAL_NEWICK, LMLN2_NEWICK,
            SIMPLE_PHYLTREE, SIMPLE_TREEBEST, SIMPLE_TREEBEST_T, SIMPLE_TREEBEST_T2):
        with open(testfile.name, 'w') as out:
            out.write(testfile.content)


#def teardown_module():
#    shutil.rmtree(TEMPDIRNAME)

class Test_Simple_TreeBest:
    kwargs = dict(outfile='-', angle_style=3, phyltreefile=SIMPLE_PHYLTREE.name, treebest=True)

    def test_onedup(self):
        gd = run([SIMPLE_TREEBEST.name], [], **self.kwargs)
    def test_onetransfer(self):
        with pytest.raises(KeyError):
            gd = run([SIMPLE_TREEBEST_T.name], [], **self.kwargs)
    def test_onetransfer2(self):
        gd = run([SIMPLE_TREEBEST_T2.name], [], **self.kwargs)


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
