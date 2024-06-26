from setuptools import setup


def readme():
    with open('README.md') as f:
        return f.read()

setup(name='Phylorgs',
      version='0.1.0-b',
      description='Tools for phylogenomics with focus on molecular dating and Ensembl database',
      long_description=readme(),
      url='https://github.com/DyogenIBENS/Phylorgs',
      author='Guillaume Louvel',
      author_email='guillaume.louvel__bioinfo'+('@')+'normalesup.org',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'License :: OSI Approved :: CEA CNRS Inria Logiciel Libre License, version 2.1 (CeCILL-2.1)',
          'Natural Language :: English',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'Environment :: Console',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Utilities'
      ],
      keywords='genomics phylogeny phylogenomics bioinformatics alignment clock',
      packages=['datasci', 'dendro', 'ensembltools', 'genchron',
                'genomicustools', 'pamliped', 'seqtools',
                'tabletools', 'taxtools', 'UItools'],
      py_modules=['archiparse', 'argparse_custom', 'genetree_drawer', 'IOtools',
                  'prep_genetree_drawing', 'objectools', 'run_environment', 'siphon'],
      #packages_data={'mypkg': 'data/*.dat'},
      #data_files=[('destdirectory', ['file1', 'file2'])]  # files to be installed relative to sys.prefix
      install_requires=[ #find_packages()
          'numpy',
          'pandas',
          'scipy',
          'matplotlib',
          'seaborn',
          'biopython',
          'ete3',
          'scikit-learn',
          'statsmodels',
          'requests',
          'markdown'
          #'LibsDyogen_py3',
          #'ToolsDyogen_py3'
          #'evosite3D'
          #'fluidcondor'
      #],
      #dependency_links=[
      #    'LibsDyogen @ git+https://github.com/DyogenIBENS/LibsDyogen_py3',
      #    'ToolsDyogen @ git+https://github.com/DyogenIBENS/ToolsDyogen_py3',
      #    'evosite3D @ git+https://github.com/romainstuder/evosite3D',
      #    'fluidcondor @ git+https://github.com/GullumLuvl/fluidcondor',
      ],
      extras_require={}, # pandas dataframe styling: Jinja2
      entry_points = {
        'console_scripts': [
            'vizdensity=datasci.vizdensity:main',
            'ProtTree_SortChildren=dendro.ProtTree_SortChildren:main',
            'ProtTree_SortLeaves=dendro.ProtTree_SortLeaves:main',
            'ProtTree_cleaner=dendro.ProtTree_cleaner:main',
            'ProtTree_fuseTips=dendro.ProtTree_fuseTips:main',
            'ProtTree_prune_species=dendro.ProtTree_prune_species:main',
            'add_branch=dendro.add_branch:main',
            'ale2treebest=dendro.ale2treebest:main',
            'annotate_phyltree_synonyms_qualities=dendro.annotate_phyltree_synonyms_qualities:main',
            'cladecollapsum=dendro.cladecollapsum:main',
            'cladematch=dendro.cladematch:main',
            'cladeplot=dendro.cladeplot:main',
            'dendroformat=dendro.formats:main',
            'dendroprint=dendro.printer:main',
            'flatten_polyphylies=dendro.flatten_polyphylies:main',
            'genetree_prune_species=dendro.genetree_prune_species:main',
            'genomicus2newick=dendro.genomicus2newick:main',
            'genomicus2treebest=dendro.genomicus2treebest:main',
            'indent_nwk=dendro.indent_nwk:main',
            'interactree=dendro.interactree:main',
            'meandist2age=dendro.meandist2age:main',
            'newick_SortLeaves=dendro.newick_SortLeaves:main',
            'numbernodes=dendro.numbernodes:main',
            'phyltree_prune_badqualspecies=dendro.phyltree_prune_badqualspecies:main',
            'reconciledtree2paralogytree=dendro.reconciledtree2paralogytree:main',
            'rename_nodes=dendro.rename_nodes:main',
            'simforest=dendro.simforest:main',
            'tabulate_nhx=dendro.tabulate_nhx:main',
            'time_fromspeciestree=dendro.time_fromspeciestree:main',
            'set_branchlengths=dendro.set_branchlengths:main',
            'transform_branchlengths=dendro.transform_branchlengths:main',
            'tree_len2node=dendro.tree_len2node:main',
            'treebest2genomicus=dendro.treebest2genomicus:main',
            'treed=dendro.treed:main',
            'unroot_binarise=dendro.unroot_binarise:main',
            'basalize=dendro.basalize:main',
            'robinsonfoulds=dendro.robinsonfoulds:main',
            'validate_monophylies=dendro.validate_monophylies:main'
            'prot2gene=ensembltools.prot2gene:main',
            'request_biomart=ensembltools.request_biomart:main',
            'allocate_codeml_memory=genchron.allocate_codeml_memory:main',
            'dSvisualizor=genchron.analyse.dSvisualizor:main',
            'generate_dNdStable=genchron.analyse.generate_dNdStable:main',
            'del_edited_branches=genchron.del_edited_branches:main',
            'del_splitgenes_in_subtrees=genchron.del_splitgenes_in_subtrees:main',
            'find_non_overlapping_codeml_results=genchron.find_non_overlapping_codeml_results:main',
            'find_omitted_sequences=genchron.find_omitted_sequences:main',
            'mark_all_branches=genchron.mark_all_branches:main',
            'MPL=genchron.MPL:main',
            'prune2family=genchron.prune2family:main',
            'reshape_tree=genchron.reshape_tree:main',
            'subtrees_stats=genchron.subtrees_stats:main',
            'genetree_drawer=genetree_drawer:main',
            'count_robusts=genomicustools.count_robusts:main',
            'forest=genomicustools.forest:main',
            'identify=genomicustools.identify:main',
            'iterate_genefamilies=genomicustools.iterate_genefamilies:main',
            'codeml_nuc2fasta=pamliped.codeml_nuc2fasta:main',
            'codeml_trees2nw=pamliped.codeml_trees2nw:main',
            'prep_genetree_drawing=prep_genetree_drawing.py',
            'al_concat=seqtools.al_concat:main',
            'al_discard_stops=seqtools.al_discard_stops:main',
            'al_len=seqtools.al_len:main',
            'al_merge_splits=seqtools.al_merge_splits:main',
            'al_order=seqtools.al_order:main',
            'backtransX=seqtools.backtransX:main',
            'bed_join_consecutive=seqtools.bed_join_consecutive:main',
            'bed_reverseminus=seqtools.bed_reverseminus:main',
            'compo_freq=seqtools.compo_freq:main',
            'compo_hetero=seqtools.compo_heterogeneity:main',
            'fasta_join_revcomp=seqtools.fasta_join_revcomp:main',
            'seq_translate=seqtools.seq_translate:main',
            'fillpositions=seqtools.fillpositions:main',
            'nei_gojobori_dNdS=phylorg.nei_gojobori_dNdS:main',
            'plotal=seqtools.plotal:main',
            'printal=seqtools.printal:main',
            'seq_from_axt=seqtools.seq_from_axt:main',
            'seq_from_maf=seqtools.seq_from_maf:main',
            'seqconv=seqtools.seqconv:main',
            'seqname_grep=seqtools.seqname_grep:main',
            'specify=seqtools.specify:main',
            'ungap=seqtools.ungap:main',
            'inverse_relation=tabletools.inverse_relation:main',
            'reprecise=tabletools.reprecise:main',
            'FindClade2=taxtools.FindClade2:main',
            'gettaxonsynonyms=taxtools.gettaxonsynonyms:main',
            'gettaxontree=taxtools.gettaxontree:main',
            'ott_speciestree=taxtools.ott_speciestree:main',
            'timetree_taxonify=taxtools.timetree_taxonify:main',
            'nbargs=UItools.nbargs:main',
            'analyse_duprates=duprates.analyse_duprates:main',
            ]
        #TODO: snake files?
          },
      scripts=[
            'dendro/check_newick.R',
            'dendro/phyltree_nwk_standardize.sh',
            'treebestpipe/exonerate2cds.pl',
            'duprates/convert-and-specify.sh',
            'duprates/generate_generax_familyfile.sh',
            'duprates/generax_family1by1.sh',
            'duprates/drop_outgroup.py',
            'duprates/ale-master.sh',
            'duprates/cafe-master.sh',
            'ensembltools/query_align_cds.pl',
            'ensembltools/query_gene_info.pl',
            'genchron/analyse/date_dup.R',
            'genchron/gather_job_characteristics.sh',
            'genchron/hmmcleaner_codons.sh',
            'genchron/hmmcleaner_codonzilla.sh',
            'genchron/hmmcleaner_fsa.sh',
            'genchron/reshape_tree.sh',
            'genchron/sim_relaxed.R',
            'genchron/beastwrap.sh',
            'pamliped/run_codeml_separatedir.sh',
            'seqtools/fa2phy.sed',
            'seqtools/index_msa.sh',
            'seqtools/query_ucsc_aligned.sh',
              ],
      tests_require=['pytest'],
      #test_suite='pytest',
      include_package_data=True,
      zip_safe=False)

# For data files, use
# import pkg_resources
# my_data = pkg_resources.resource_string(__name__, "foo.dat")
