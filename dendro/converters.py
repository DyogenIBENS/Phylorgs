
"""Module to convert Python tree objects between different classes.

(Ete3, LibsDyogen, Bio.Phylo, Bio.Nexus, scipy.hclust.linkage...)

Also see dendro.formats to switch between *file* formats.
"""


from dendro.bates import rev_dfw_descendants, dfw_descendants_generalized
import logging


logger = logging.getLogger(__name__)


def get_data(tree, nodedist):
    return tree.data.get(nodedist[0], [])


def ete3_to_linkage(tree):
    ### PB: distances are identical on each side of any fork.
    Z = []
    node_is = {}
    leaf_i = 0
    for node_dist, children_dists in \
                                rev_dfw_descendants(tree,
                                                get_data,
                                                include_leaves=True,
                                                queue=[(tree.root, tree.dist)]):
        if not children_dists:
            node_is[node_dist[0]] = leaf_i
            leaf_i += 1
            continue

        children_is = [node_is.pop(ch) for ch,_ in children_dists]
        Z.append([])


def ProtTree_to_ete3(prottree):
    import ete3
    tree = ete3.Tree(name=prottree.root, dist=0)
    current_nodes = [tree]
    while current_nodes:
        current = current_nodes.pop()
        current.add_features(**prottree.info[current.name])
        for child, dist in prottree.data.get(current.name, []):
            current_nodes.append(current.add_child(name=child, dist=dist))
    return tree


def PhylTree_to_ete3(phyltree, nosinglechild=False):
    import ete3
    tree = ete3.Tree(name=phyltree.root, dist=getattr(phyltree, 'rootlength', 0))
    tree.add_features(treename=getattr(phyltree, 'name', ''))
    current_nodes = [tree]

    phyltree_nodedata = []
    for ft in ('ages', 'commonNames', 'fileName', 'indBranches', 'indNames'):
        try:
            phyltree_nodedata.append((ft, getattr(phyltree, ft)))
        except AttributeError:
            logger.warning('phyltree has no attribute %r.', ft)

    while current_nodes:
        current = current_nodes.pop()
        extra_features = {}
        current.add_features(**{ft: data.get(current.name) for ft, data in phyltree_nodedata})
        childlist = phyltree.items.get(current.name, [])
        if not childlist:
            current.add_features(Esp2X=(current.name in phyltree.lstEsp2X),
                                 Esp6X=(current.name in phyltree.lstEsp6X),
                                 EspFull=(current.name in phyltree.lstEspFull))
            # dicGenes
            # dicGenomes
        while nosinglechild and (len(childlist) == 1):
            current.name = childlist[0][0]
            current.dist += childlist[0][1]
            childlist = phyltree.items.get(current.name, [])

        for child, dist in childlist:
            current_nodes.append(current.add_child(name=child, dist=dist))
    return tree


def PhylTree_draw(phyltree, images=None):
    """Draw tree in a gui window using Ete3"""
    import ete3
    ptree = PhylTree_to_ete3(phyltree)
    ts = ete3.TreeStyle()
    ts.scale = 2
    ts.show_branch_length = True
    nameface = ete3.faces.AttrFace("name", fsize=12)
    def mylayout(node):
        if node.is_leaf():
            # add species image
            pass
        else:
            ete3.faces.add_face_to_node(nameface, node, column=0,
                                        position='branch-top')
    ptree.show(layout=mylayout, tree_style=ts, name='Species tree')


def ete3_to_ete3(tree):
    return tree

def BioPhylo_to_BioNexusTree(tree):
    from Bio.Nexus.Trees import Tree, Nodes, NodeData
    from dendro.any import BioPhylo

    ntree = Tree(name=tree.name, rooted=tree.rooted, weight=tree.weight)
    root = ntree.node(0)
    root.set_data(NodeData(taxon=tree.root.name,
                           branchlength=tree.root.branch_length,
                           support=tree.root.confidence,
                           comment=tree.root.comment))

    node2ids = {tree.root: 0}
    for parent, children in dfw_descendants_generalized(tree, BioPhylo.get_children):
        parent_id = node2ids[parent]
        for child in children:
            new = Nodes.Node(NodeData(taxon=child.name,
                                      branchlength=child.branch_length,
                                      support=child.confidence,
                                      comment=child.comment))

            node2ids[child] = ntree.add(new, parent_id)
    return ntree


def BioNexusTrees_to_BioPhylo(ntrees, translate=None):
    """If a translate dict is given, encode node labels by integer ids"""
    if translate:
        untranslate = {taxon: str(tax_id) for tax_id, taxon in translate.items()}
    from Bio.Phylo.BaseTree import Clade, Tree
    trees = []
    for idx, ntree in enumerate(ntrees):
        nroot = ntree.node(ntree.root)
        if nroot.data.taxon and translate:
            nodename = untranslate[nroot.data.taxon]
        else:
            nodename = nroot.data.taxon
        root = Clade(branch_length=nroot.data.branchlength,
                     name=nodename,
                     confidence=nroot.data.support)
        tree = Tree(root, rooted=ntree.rooted, id=idx, name=ntree.name)
        tree.weight = ntree.weight
        matching_clades = {nroot: root}  # nexus node -> Phylo.BaseTree.Clade
        queue = [nroot]
        while queue:
            nnode = queue.pop(0)
            node = matching_clades.pop(nnode)
            nchildren = [ntree.node(ch_id) for ch_id in nnode.succ]
            for nchild in nchildren:
                if nchild.data.taxon and translate:
                    nodename = untranslate[nchild.data.taxon]
                else:
                    nodename = nchild.data.taxon
                child = Clade(branch_length=nchild.data.branchlength,
                             name=nodename,
                             confidence=nchild.data.support)
                child.comment = nchild.data.comment
                matching_clades[nchild] = child
                node.clades.append(child)
                queue.append(nchild)
        trees.append(tree)
    return trees


def One_BioNexusTree_to_BioPhylo(ntree):
    return BioNexusTrees_to_BioPhylo([ntree])[0]


converterchoice = {'PhylTree': {'Ete3': PhylTree_to_ete3},
                   'ProtTree': {'Ete3': ProtTree_to_ete3},
                   'ProteinTree': {'Ete3': ProtTree_to_ete3},
                   'Ete3': {'Ete3': ete3_to_ete3},
                   'BioPhylo': {'BioNexus': BioPhylo_to_BioNexusTree},
                   'BioNexus': {'BioPhylo': One_BioNexusTree_to_BioPhylo}}
for from_, todict in list(converterchoice.items()):
    for to_, func in list(todict.items()):
        todict[to_.lower()] = func
    converterchoice[from_.lower()] = todict
converterchoice['phyltree']['ete3_f1'] = converterchoice['phyltree']['ete3']
converterchoice['prottree']['ete3_f1'] = converterchoice['prottree']['ete3']
converterchoice['ete3_f1'] = converterchoice['ete3']


VALID_CONVERSIONS = """
Valid conversions:

 * PhylTree|ProtTree -> Ete3  (the reverse is not useful because PhylTree/ProtTree have newick parsers)
 * BioPhylo <-> BioNexus      (although BioNexus cannot be used to output newick)
"""


def _eval_scalar(val):
    if val.lower() in ('none', 'false', 'true'):
        return eval(val.capitalize())
    try:
        return int(val)
    except ValueError:
        try:
            return float(val)
        except ValueError:
            # val will remain a string.
            return val


def main():
    from sys import stdin
    import argparse as ap
    from dendro.parsers import chooseparser
    from dendro.any import methodchoice

    parser = ap.ArgumentParser(description=__doc__, epilog=VALID_CONVERSIONS,
                               formatter_class=ap.RawDescriptionHelpFormatter)
    parser.add_argument('inputtree', default=stdin)
    parser.add_argument('-p', '--parser', default='phyltree', help='[%(default)s]')
    parser.add_argument('-o', '--outformat', default='ete3', help='[%(default)s]')
    parser.add_argument('-w', '--out-kwargs', nargs='+', default=[],
                        help=('list of key=value arguments to give to '
                              '`print_newick`. Ex: "features=[]" for the ete3 '
                              'output format.'))

    args = parser.parse_args()

    parse_tree = chooseparser(args.parser)
    treetype = args.parser.split(':', 1)[0].strip().lower()
    outformat = args.outformat.lower()
    convert = converterchoice[treetype][outformat]
    print_newick = methodchoice[outformat].print_newick

    kwargs = {}
    for keyval in args.out_kwargs:
        key, val = keyval.split('=', 1)
        if val[0] == '[' and val[-1] == ']':
            newval = []
            for elem in val[1:-1].split(','):
                newval.append(_eval_scalar(elem.strip()))
        else:
            newval = _eval_scalar(val)
        kwargs[key] = newval

    for tree in parse_tree(args.inputtree):
        print_newick(convert(tree), **kwargs)


if __name__ == '__main__':
    main()
