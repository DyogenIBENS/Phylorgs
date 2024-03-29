#!/usr/bin/env python3


"""Module to help converting identifiers for proteins/genes/species from Genomicus/Ensembl.
"""

from __future__ import print_function

from sys import stderr, stdin
import os.path as op
import re
from functools import partial
from copy import deepcopy
import logging
logger = logging.getLogger(__name__)
#ch = logging.StreamHandler()
#ch.setFormatter(logging.Formatter("%(levelname)s:%(module)s l.%(lineno)d:%(message)s"))
#logger.addHandler(ch)


ENSEMBL_VERSION = 85

CELEGANS_REG = re.compile(r'^([234]R|A[CH]|B[0E]|C[0-5C-E]|cT|[DFHKMRTWY][0-9CHYT]|E[0_EG]|JC|LL|SS|Z[CK]|V[BCFHMWYZ]|P[AD]).*')

PROT2SP = {85:
           { #'Y': 'Saccharomyces cerevisiae',  # there are C.elegans prot with 'Y' too
            'Q0': 'Saccharomyces cerevisiae',
            'FB': 'Drosophila melanogaster',
            #'WBGene0': 'Caenorhabditis elegans',  # No consensus
            '': 'Caenorhabditis elegans',  # just so that the values reflect all species.
            'ENSCINP': 'Ciona intestinalis',
            'ENSCSAV': 'Ciona savignyi',
            'ENSCSAP': 'Chlorocebus sabaeus',
            'ENSPMAP': 'Petromyzon marinus',
            'ENSXETP': 'Xenopus tropicalis',
            'ENSPSIP': 'Pelodiscus sinensis',
            'ENSGALP': 'Gallus gallus',
            'ENSMGAP': 'Meleagris gallopavo',
            'ENSTGUP': 'Taeniopygia guttata',
            'ENSFALP': 'Ficedula albicollis',
            'ENSAPLP': 'Anas platyrhynchos',
            'ENSACAP': 'Anolis carolinensis',
            'ENSOANP': 'Ornithorhynchus anatinus',
            'ENSMEUP': 'Macropus eugenii',
            'ENSSHAP': 'Sarcophilus harrisii',
            'ENSMODP': 'Monodelphis domestica',
            'ENSLAFP': 'Loxodonta africana',
            'ENSETEP': 'Echinops telfairi',
            'ENSPCAP': 'Procavia capensis',
            'ENSDNOP': 'Dasypus novemcinctus',
            'ENSCHOP': 'Choloepus hoffmanni',
            'ENSSARP': 'Sorex araneus',
            'ENSEEUP': 'Erinaceus europaeus',
            'ENSMLUP': 'Myotis lucifugus',
            'ENSPVAP': 'Pteropus vampyrus',
            'ENSTTRP': 'Tursiops truncatus',
            'ENSBTAP': 'Bos taurus',
            'ENSOARP': 'Ovis aries',
            'ENSVPAP': 'Vicugna pacos',
            'ENSSSCP': 'Sus scrofa',
            'ENSECAP': 'Equus caballus',
            'ENSMPUP': 'Mustela putorius furo',
            'ENSAMEP': 'Ailuropoda melanoleuca',
            'ENSCAFP': 'Canis lupus familiaris',
            'ENSFCAP': 'Felis catus',
            'ENSTBEP': 'Tupaia belangeri',
            'ENSPANP': 'Papio anubis',
            'ENSMMUP': 'Macaca mulatta',
            'ENSPPYP': 'Pongo abelii',
            'ENSGGOP': 'Gorilla gorilla gorilla',
            'ENSPTRP': 'Pan troglodytes',
            'ENSP000':    'Homo sapiens',       # ENSG
            'ENSNLEP': 'Nomascus leucogenys',
            'ENSCJAP': 'Callithrix jacchus',
            'ENSTSYP': 'Tarsius syrichta',
            'ENSOGAP': 'Otolemur garnettii',
            'ENSMICP': 'Microcebus murinus',
            'ENSOPRP': 'Ochotona princeps',
            'ENSOCUP': 'Oryctolagus cuniculus',
            'ENSCPOP': 'Cavia porcellus',
            'ENSRNOP': 'Rattus norvegicus',
            'ENSMUSP': 'Mus musculus',
            'ENSSTOP': 'Ictidomys tridecemlineatus',
            'ENSDORP': 'Dipodomys ordii',
            'ENSLACP': 'Latimeria chalumnae',
            'ENSLOCP': 'Lepisosteus oculatus',
            'ENSGACP': 'Gasterosteus aculeatus',
            'ENSTNIP': 'Tetraodon nigroviridis',
            'ENSTRUP': 'Takifugu rubripes',
            'ENSONIP': 'Oreochromis niloticus',
            'ENSORLP': 'Oryzias latipes',
            'ENSPFOP': 'Poecilia formosa',
            'ENSXMAP': 'Xiphophorus maculatus',
            'ENSGMOP': 'Gadus morhua',
            'ENSAMXP': 'Astyanax mexicanus',
            'ENSDARP': 'Danio rerio'}
           }

PROT2SP[86] = deepcopy(PROT2SP[85])
PROT2SP[86].update({'MGP_SPR': 'Mus spretus'})  # TODO: check the protein id

PROT2SP[87] = deepcopy(PROT2SP[86])
PROT2SP[87].update({'ENSTSYP': 'Carlito syrichta'})

PROT2SP[88] = PROT2SP[87]

PROT2SP[90] = deepcopy(PROT2SP[88])
PROT2SP[90].update({'ENSMEUP': 'Notamacropus eugenii',  # Update of Macropus e.
                    'ENSCAPP': 'Cavia aperea',
                    'ENSCLAP': 'Chinchilla lanigera',
                    #'ENSCGRP00001': 'Cricetulus griseus CHOK1GS',
                    #'ENSCGRP00000': 'Cricetulus griseus Crigri',  # !!!!
                    'ENSCGRP': 'Cricetulus griseus',  # !!!!
                    'ENSFDAP':      'Fukomys damarensis',
                    #'ENSHGLP00000': 'Heterocephalus glaber female',
                    #'ENSHGLP00100': 'Heterocephalus glaber male',
                    'ENSHGLP': 'Heterocephalus glaber',
                    'ENSJJAP': 'Jaculus jaculus',
                    'ENSMAUP': 'Mesocricetus auratus',
                    'ENSMOCP': 'Microtus ochrogaster',
                    'MGP_CAR': 'Mus caroli',
                    'MGP_Pah': 'Mus pahari',
                    'ENSNGAP': 'Nannospalax galili',
                    'ENSODEP': 'Octodon degus',
                    'ENSPEMP': 'Peromyscus maniculatus bairdii'})

PROT2SP[91] = deepcopy(PROT2SP[90])
PROT2SP[91].update({'ENSANAP': 'Aotus nancymaae',
                    'ENSCCAP': 'Cebus capucinus imitator',  # 'Cebus capucinus'
                    'ENSCATP': 'Cercocebus atys',
                    'ENSCANP': 'Colobus angolensis palliatus',
                    'ENSMFAP': 'Macaca fascicularis',
                    'ENSMNEP': 'Macaca nemestrina',
                    'ENSMLEP': 'Mandrillus leucophaeus',
                    'ENSPPAP': 'Pan paniscus',
                    'ENSPCOP': 'Propithecus coquereli',
                    'ENSRBIP': 'Rhinopithecus bieti',
                    'ENSRROP': 'Rhinopithecus roxellana',
                    'ENSSBOP': 'Saimiri boliviensis boliviensis'})

PROT2SP[92] = deepcopy(PROT2SP[91])
PROT2SP[92].update({'ENSCHIP':'Capra hircus'})

PROT2SP[93] = deepcopy(PROT2SP[92])
PROT2SP[93].update({'ENSEBUP': 'Eptatretus burgeri',
                    'ENSPPRP': 'Panthera pardus',
                    'ENSPTIP': 'Panthera tigris altaica'})

GENE2SP = {}

GENE2SP[85] = {'Y': 'Saccharomyces cerevisiae',
               'Q': 'Saccharomyces cerevisiae',
               'FBgn': 'Drosophila melanogaster',
               'WBGene0': 'Caenorhabditis elegans',
               'ENSCING': 'Ciona intestinalis',
               'ENSCSAV': 'Ciona savignyi',
               'ENSCSAG': 'Chlorocebus sabaeus',
               'ENSPMAG': 'Petromyzon marinus',
               'ENSXETG': 'Xenopus tropicalis',
               'ENSPSIG': 'Pelodiscus sinensis',
               'ENSGALG': 'Gallus gallus',
               'ENSMGAG': 'Meleagris gallopavo',
               'ENSTGUG': 'Taeniopygia guttata',
               'ENSFALG': 'Ficedula albicollis',
               'ENSAPLG': 'Anas platyrhynchos',
               'ENSACAG': 'Anolis carolinensis',
               'ENSOANG': 'Ornithorhynchus anatinus',
               'ENSMEUG': 'Macropus eugenii',
               'ENSSHAG': 'Sarcophilus harrisii',
               'ENSMODG': 'Monodelphis domestica',
               'ENSLAFG': 'Loxodonta africana',
               'ENSETEG': 'Echinops telfairi',
               'ENSPCAG': 'Procavia capensis',
               'ENSDNOG': 'Dasypus novemcinctus',
               'ENSCHOG': 'Choloepus hoffmanni',
               'ENSSARG': 'Sorex araneus',
               'ENSEEUG': 'Erinaceus europaeus',
               'ENSMLUG': 'Myotis lucifugus',
               'ENSPVAG': 'Pteropus vampyrus',
               'ENSTTRG': 'Tursiops truncatus',
               'ENSBTAG': 'Bos taurus',
               'ENSOARG': 'Ovis aries',
               'ENSVPAG': 'Vicugna pacos',
               'ENSSSCG': 'Sus scrofa',
               'ENSECAG': 'Equus caballus',
               'ENSMPUG': 'Mustela putorius furo',
               'ENSAMEG': 'Ailuropoda melanoleuca',
               'ENSCAFG': 'Canis lupus familiaris',
               'ENSFCAG': 'Felis catus',
               'ENSTBEG': 'Tupaia belangeri',
               'ENSPANG': 'Papio anubis',
               'ENSMMUG': 'Macaca mulatta',
               'ENSPPYG': 'Pongo abelii',
               'ENSGGOG': 'Gorilla gorilla gorilla',
               'ENSPTRG': 'Pan troglodytes',
               'ENSG000':    'Homo sapiens',       # ENSG
               'ENSNLEG': 'Nomascus leucogenys',
               'ENSCJAG': 'Callithrix jacchus',
               'ENSTSYG': 'Tarsius syrichta',
               'ENSOGAG': 'Otolemur garnettii',
               'ENSMICG': 'Microcebus murinus',
               'ENSOPRG': 'Ochotona princeps',
               'ENSOCUG': 'Oryctolagus cuniculus',
               'ENSCPOG': 'Cavia porcellus',
               'ENSRNOG': 'Rattus norvegicus',
               'ENSMUSG': 'Mus musculus',
               'ENSSTOG': 'Ictidomys tridecemlineatus',
               'ENSDORG': 'Dipodomys ordii',
               'ENSLACG': 'Latimeria chalumnae',
               'ENSLOCG': 'Lepisosteus oculatus',
               'ENSGACG': 'Gasterosteus aculeatus',
               'ENSTNIG': 'Tetraodon nigroviridis',
               'ENSTRUG': 'Takifugu rubripes',
               'ENSONIG': 'Oreochromis niloticus',
               'ENSORLG': 'Oryzias latipes',
               'ENSPFOG': 'Poecilia formosa',
               'ENSXMAG': 'Xiphophorus maculatus',
               'ENSGMOG': 'Gadus morhua',
               'ENSAMXG': 'Astyanax mexicanus',
               'ENSDARG': 'Danio rerio'}

GENE2SP[86] = deepcopy(GENE2SP[85])
GENE2SP[86].update(**{'MGP_SPR': 'Mus spretus'})

GENE2SP[87] = deepcopy(GENE2SP[86])
GENE2SP[87].update({'ENSTSYG': 'Carlito syrichta'})
GENE2SP[88] = GENE2SP[87]

GENE2SP[90] = deepcopy(GENE2SP[88])
GENE2SP[90].update({'ENSMEUG': 'Notamacropus eugenii',  # Update of Macropus e.
                    'ENSCAPG': 'Cavia aperea',
                    'ENSCLAG': 'Chinchilla lanigera',
                    #'ENSCGRG00001': 'Cricetulus griseus CHOK1GS',
                    #'ENSCGRG00000': 'Cricetulus griseus Crigri',  # !!!!
                    'ENSCGRG': 'Cricetulus griseus',
                    'ENSFDAG':      'Fukomys damarensis',
                    #'ENSHGLG00000': 'Heterocephalus glaber female',
                    #'ENSHGLG00100': 'Heterocephalus glaber male',
                    'ENSHGLG': 'Heterocephalus glaber',
                    'ENSJJAG': 'Jaculus jaculus',
                    'ENSMAUG': 'Mesocricetus auratus',
                    'ENSMOCG': 'Microtus ochrogaster',
                    'MGP_CAR': 'Mus caroli',
                    'MGP_Pah': 'Mus pahari',
                    'ENSNGAG': 'Nannospalax galili',
                    'ENSODEG': 'Octodon degus',
                    'ENSPEMG': 'Peromyscus maniculatus bairdii'})

GENE2SP[91] = deepcopy(GENE2SP[90])
GENE2SP[91].update({'ENSANAG': 'Aotus nancymaae',
                    'ENSCCAG': 'Cebus capucinus imitator',  # 'Cebus capucinus'
                    'ENSCATG': 'Cercocebus atys',
                    'ENSCANG': 'Colobus angolensis palliatus',
                    'ENSMFAG': 'Macaca fascicularis',
                    'ENSMNEG': 'Macaca nemestrina',
                    'ENSMLEG': 'Mandrillus leucophaeus',
                    'ENSPPAG': 'Pan paniscus',
                    'ENSPCOG': 'Propithecus coquereli',
                    'ENSRBIG': 'Rhinopithecus bieti',
                    'ENSRROG': 'Rhinopithecus roxellana',
                    'ENSSBOG': 'Saimiri boliviensis boliviensis'})

GENE2SP[92] = deepcopy(GENE2SP[91])
GENE2SP[92].update({'ENSCHIG':'Capra hircus'})

GENE2SP[93] = deepcopy(GENE2SP[92])
GENE2SP[93].update({'ENSEBUG': 'Eptatretus burgeri',
                    'ENSPPRG': 'Panthera pardus',
                    'ENSPTIG': 'Panthera tigris altaica'})


SP2GENEID = {}
for version, conversion in GENE2SP.items():
    SP2GENEID[version] = {}
    for geneid, species in conversion.items():
        if species in SP2GENEID[version]:
            SP2GENEID[version][species] += '|' + geneid
        else:
            SP2GENEID[version][species] = geneid
#from tabletools import inverse

#@myTools.memoize
def get_available_ensembl_version(version, available):
    # Search sorted and get the element after.
    for i, av in enumerate(sorted(available)):
        if version <= av:
            if version != av:
                logger.debug('Iter %d', i)
                logger.warning('No precomputed conversion dict for '
                               'version %s, fall back to %s.',
                               version, av)
            return av

    logger.warning('All precomputed versions are older than %s. '
                   'Fall back to the most recent one: %s.',
                   version, av)
    return av


class fallbackContent(object):

    def __init__(self, srcdict, fallbackfunc):  # warn=True
        self.available = srcdict
        self._fallbackfunc = fallbackfunc
        self.fallbacks = {}

    def __getitem__(self, key):
        try:
            return self.available[key]
        except KeyError:
            return self.available[self.set_fallback(key)]
    
    def get(self, key):
        try:
            return self.__getitem__(key)
        except KeyError:
            return None

    def set_fallback(self, key):
        try:
            return self.fallbacks[key]
        except KeyError:
            fallback_key = self._fallbackfunc(key)
            self.fallbacks[key] = fallback_key
            return fallback_key


class fallback_gene2species(fallbackContent):
    def __init__(self):
        fallbackContent.__init__(self,
            GENE2SP,
            partial(get_available_ensembl_version, available=GENE2SP.keys())
            )

class fallback_prot2species(fallbackContent):
    def __init__(self):
        fallbackContent.__init__(self,
            PROT2SP,
            partial(get_available_ensembl_version, available=PROT2SP.keys())
            )


GENE2SP_F = fallback_gene2species()
PROT2SP_F = fallback_prot2species()


def convert_prot2species(modernID, ensembl_version=ENSEMBL_VERSION, default=None):
    ensembl_version = PROT2SP_F.set_fallback(ensembl_version)
    if ensembl_version >= 90:
        try:
            return PROT2SP[ensembl_version][modernID[:12]]
        except KeyError:
            pass
    try:
        return PROT2SP[ensembl_version][modernID[:7]]
    except KeyError:
        try:
            # Saccharomyces cerevisiae (Q0) or Drosophila melanogaster
            return PROT2SP[ensembl_version][modernID[:2]]
        except KeyError as err:
            if re.match('Y[A-Z]', modernID):
                return 'Saccharomyces cerevisiae'
            elif CELEGANS_REG.match(modernID):
                return 'Caenorhabditis elegans'
            elif default is not None:
                return default
            else:
                raise KeyError(err.args[0] + ' (protein: %s, Ens.%d)' % \
                                (modernID, ensembl_version))



def convert_gene2species(modernID, ensembl_version=ENSEMBL_VERSION):
    """Give the species of an **extant** sequence id (from ENSEMBL)."""
    ensembl_version = GENE2SP_F.set_fallback(ensembl_version)
    gene2sp = GENE2SP[ensembl_version]
    if ensembl_version >= 90:
        try:
            return gene2sp[modernID[:12]]
        except KeyError:
            pass
    try:
        return gene2sp[modernID[:7]]
    except KeyError:
        try:
            return gene2sp[modernID[:4]]
        except KeyError:
            try:
                return gene2sp[modernID[0]]
            except KeyError:
                raise KeyError("%s can't be assigned to a *species* (from "
                               "Ensembl %d)" % (modernID, ensembl_version))


def convert_prot2gene(protID, gene_info, cprot=2, cgene=0, shorten_species=False,
                      ensembl_version=ENSEMBL_VERSION):
    warning.warn('Bad (inefficient) function (too much IO). You should load all conversions at once in a dict.')
    sp = convert_prot2species(protID, ensembl_version)
    if shorten_species:
        spsplit = sp.split()
        sp2 = spsplit[0][0].lower() + spsplit[-1]
    else:
        sp2 = sp.replace(' ', '.')
    return grep_prot(gene_info % sp2, protID, cprot=cprot, cgene=cgene)


def load_assembly2species(filename="~/ws2/UCSC_genome_releases_full.tsv",
                          fromcol=2, tocol=0):
    conversion = {}
    with open(op.expanduser(filename)) as stream:
        header = next(stream)
        for line in stream:
            fields = line.rstrip().split('\t')
            conversion[fields[fromcol]] = fields[tocol]
    return conversion

def try_load_assembly2species(ucsc_conv_filename='~/ws7/UCSC_genome_releases_full.tsv'):
    try:
        return load_assembly2species(ucsc_conv_filename)
    except FileNotFoundError:
        logger.warning("Conversion file not found: %r", ucsc_conv_filename)
        return {}


class OnDemandData(object):
    """A dictionary-like object that is filled only when someone queries it.
    Read-only.

    Allows to save loading time if unused (e.g. Disk I/O...).
    """
    #data = None  # Default *class* attribute needed during unpickling.
    def __init__(self, loader):
        """loader: function returning a dictionary"""
        self.data = None
        self.loader = loader
        self.loaded = False

    def load(self):
        """Load or reload"""
        self.data = self.loader()

    def _load_on_error(method):
        """Tries first, or load data"""
        def newmethod(self, *args, **kwargs):
            try:
                # For example, try: dict.__getitem__(None, 'a')  -> TypeError
                return method(self, *args, **kwargs)
            except AttributeError:  # Because data is None
                self.data = self.loader()
                self.loaded = True
                return method(self, *args, **kwargs)
        return newmethod

    @_load_on_error
    def __getitem__(self, key):
        return self.data.__getitem__(key)

    @_load_on_error
    def __iter__(self):
        return self.data.__iter__()

    @_load_on_error
    def __len__(self):
        return self.data.__len__()  # Why two definitions of length...

    @_load_on_error
    def get(self, key, default):
        return self.data.get(key, default)

    @_load_on_error
    def items(self):
        return self.data.items()

    @_load_on_error
    def keys(self):
        return self.data.keys()

    @_load_on_error
    def __contains__(self, key):
        return self.data.__contains__(key)

    def __repr__(self):
        return '%s(loader=%s)' % ((self.__class__.__name__, self.loader))


assembly2species = OnDemandData(try_load_assembly2species)


def ultimate_seq2sp(seqname, ensembl_version=ENSEMBL_VERSION):
    """From an extant sequence name, find the corresponding species.
    Recognizes Ensembl gene IDs, Ensembl protein IDs, and also UCSC assembly
    names such as 'loxAfr3'"""
    try:
        sp = convert_gene2species(seqname, ensembl_version)
    except KeyError:
        try:
            sp = convert_prot2species(seqname, ensembl_version)
        except KeyError:
            assembly = re.match('[A-Za-z0-9]+', seqname).group()
            sp = assembly2species[assembly]
    return sp



def main():
    import argparse as ap
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('infile', nargs='?', default=stdin, type=ap.FileType('r'),
                        help='One identifier per line.')
    parser.add_argument('-e', '--ensembl-version', type=int, default=93)
    # Return the corresponding species
    args = parser.parse_args()

    for line in args.infile:
        print(ultimate_seq2sp(line.rstrip(), args.ensembl_version))


if __name__ == '__main__':
    main()
