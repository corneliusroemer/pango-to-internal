#! /usr/local/Caskroom/mambaforge/base/envs/nextstrain/bin/python

# Get all insertions and code as binary vector: Use metadata for this
# Load tree
# Reconstruct using tree time
#%%
import copy
import json
from collections import defaultdict

import augur.ancestral as ancestral
import augur.utils
import Bio.Align
import click
import pandas as pd
from BCBio import GFF
from mergedeep import merge
import pango_designation
import itertools
from treetime import TreeAnc

from Bio import Phylo


# @click.command()
# @click.option("--designations", default="test_data/designations.csv", type=str)
# @click.option("--tree", default="test_data/tree.nwk", type=str)
# @click.option("--alias", default="test_data/alias.json", type=str)
# @click.option("--output", type=click.File("w"))
# def main(designations, tree, alias, output, genemap):
#%%
class Aliasor:
    def __init__(self, alias_file):
        import pandas as pd

        aliases = pd.read_json(alias_file)

        self.alias_dict = {}
        for column in aliases.columns:
            if column.startswith("X"):
                self.alias_dict[column] = column
            else:
                self.alias_dict[column] = aliases[column][0]

        self.alias_dict["A"] = "A"
        self.alias_dict["B"] = "B"

        self.realias_dict = {v: k for k, v in self.alias_dict.items()}

    def compress(self, name):
        name_split = name.split(".")
        if len(name_split) < 5:
            return name
        letter = self.realias_dict[".".join(name_split[0:4])]
        if len(name_split) == 5:
            return letter + "." + name_split[4]
        else:
            return letter + "." + ".".join(name_split[4:])

    def uncompress(self, name):
        name_split = name.split(".")
        letter = name_split[0]
        unaliased = self.alias_dict[letter]
        if len(name_split) == 1:
            return name
        if len(name_split) == 2:
            return unaliased + "." + name_split[1]
        else:
            return unaliased + "." + ".".join(name_split[1:])


#%%
# designations = "test_data/designations.csv"
designations = "test_data/real_designations.csv"
tree = "test_data/tree.nwk"
alias = "test_data/alias.json"
#%%
aliasor = Aliasor(alias)
#%%
meta = pd.read_csv(designations, index_col=0)

#%%
# Filter meta down to strains in tree
tree = Phylo.read("test_data/nextstrain__tree.nwk", format="newick")
# print(tree)
unpruned_tree = copy.deepcopy(tree)
#%%
tips = list(map(lambda x: x.name, tree.get_terminals()))

#%%
meta = meta[meta.index.isin(tips)]
#%%
meta['unaliased'] = meta['lineage'].apply(aliasor.uncompress)

#%%
# Find all lineages so we can map lineage to position
# Unalias
# Make sure to add intermediary lineages at all levels of the hierarchy
# Then get dicts to map lineage <-> position

#%%

lineages_unaliased = set(meta.unaliased.unique())
lineages_unaliased
#%%
def get_lineage_hierarchy(lineage):
    """
    Takes lineage and returns list including all parental lineages
    >>> get_lineage_hierarchy("A.B.C.D.E")
    ['A', 'A.B', 'A.B.C', 'A.B.C.D', 'A.B.C.D.E']
    """
    lineage_split = lineage.split(".")
    hierarchy = []
    if len(lineage_split) == 1:
        return lineage
    for i in range(len(lineage_split)):
        hierarchy.append(".".join(lineage_split[0:i+1]))
    return hierarchy
    
#%%
lineage_hierachy = set(itertools.chain.from_iterable(map(get_lineage_hierarchy, lineages_unaliased)))
lineage_hierachy
#%%

# Construct of insertions
characters = sorted(lineage_hierachy)

# Dictionary that maps property to position in pseudosequence
mapping = {}
for pos, character in enumerate(characters):
    mapping[character] = pos

inverse_mapping = {v: k for k, v in mapping.items()}
#%%

# Take lineage and turn into binary vector
# Convention: A = 0, G = 1
def lineage_to_vector(lineage):
    """Turn lineage into binary vector"""
    character_list = get_lineage_hierarchy(lineage)
    vector = ["A"] * len(characters)
    for character in character_list:
        vector[mapping[character]] = "G"
    return "".join(vector)

meta = meta.assign(trait_vector=meta['unaliased'].apply(lineage_to_vector))

meta
#%%

# Print mapping
# print(mapping)
# Print all sequences with insertions
# print(meta.insertions_vector[meta.insertions_vector.apply(lambda x: "G" in x) > 0])

# Transform binary vector into BioSeq alignment
alignment = Bio.Align.MultipleSeqAlignment(
    [
        Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(vector), id=strain, name=strain)
        for strain, vector in meta.trait_vector.iteritems()
    ]
)
print(alignment)

#%%
# Prune tree to only include designated tips



#%%
# Prune tree
# for node in tree.get_terminals():
#     if node.name not in meta.index:
#         tree.prune(node)
# print(tree)

#%%
tt_tree = copy.deepcopy(tree)
tt_tree
#%%
tt = TreeAnc(tree=tree, aln=alignment)
# tt.infer_ancestral_sequences(method='fitch',reconstruct_tip_states=True)
tt.infer_ancestral_sequences(method='ml',reconstruct_tip_states=True)
nodes_with_mutations = ancestral.collect_mutations_and_sequences(tt)

#%%
# Get reconstructed alignment
aln = tt.get_tree_dict()
print(aln)
#%%
# Convert reconstructed alignment to lineage
def seq_to_lineage_list(seq):
    """
    Takes a BioSeq object and returns a list of lineages
    >>> seq_to_lineage_list("GGAAG")
    ['B', 'B.1', 'B.1.3']
    """
    lineage_list = []
    for index, character in enumerate(seq):
        if character == "G":
            lineage_list.append(inverse_mapping[index])
    return lineage_list

def lineage_list_to_lineage(lineage_list):
    return lineage_list[-1]

# seq_to_lineage_list("GGAAAAAGAAA")
#%%
meta['reconstructed'] = meta['unaliased']
meta
#%%
# TODO: Infer reconstructed tips too
internal_names = [node.name for node in tree.get_nonterminals()]
for rec in aln:
    # if rec.id in internal_names:
    meta.loc[rec.id,'reconstructed_full'] = str(seq_to_lineage_list(rec.seq))
    meta.loc[rec.id,'reconstructed'] = lineage_list_to_lineage(seq_to_lineage_list(rec.seq))

meta
#%%
meta.reconstructed_full.value_counts().to_csv("test_data/reconstructed_full.csv")

#%%
# Walk tree from root to tips and inferring lineage based on parents
# root = 'M'
# for clade in unpruned_tree.find_clades():
#     if clade.name not in meta.index:
#         path = unpruned_tree.get_path(target=clade) 
#         for i, node in enumerate(path):
#             if node.name not in meta.index:
#                 meta.loc[node.name,'reconstructed'] = meta.loc[path[i-1].name,'reconstructed']

# meta
#%%
meta['realiased'] = meta['reconstructed'].apply(aliasor.compress)
meta
#%%
# meta.realiased.value_counts().to_csv("test_data/realiased.csv")
#%%
output = output = {"nodes": {}}
for row in meta.itertuples():
    output["nodes"][row.Index] = {"pango_lineage": row.realiased}

json.dump(output, open("pango_reconstructed.json","w"), indent=4)

#%%
meta['realiased'].to_csv("test_data/realiased.csv")
#%%
augur.utils.write_json({"nodes":meta['realiased'].rename('inferred_lineage').to_frame().to_dict(orient='index')}, "test_data/lineage.json")

#%%
# if trait == "insertions":
#     gm = list(GFF.parse(genemap))[0]
#     def position_to_gene(position):
#         """Convert position to gene using genemap"""
#         gene_name = codon_number = reading_frame = None
#         for feature in gm.features:
#             if position in feature.location:
#                 gene_name = feature.qualifiers["gene_name"][0]
#                 codon_number = (position - feature.location.start) // 3 + 1
#                 reading_frame = (position - feature.location.start - 1) % 3
#         return {"gene_name": gene_name, "codon_number": codon_number, "reading_frame": reading_frame}

# def mut_to_str(mut: str) -> dict:
#     """Convert mutation to string"""
#     if trait == "insertions":
#         indel_type = "ins" if mut[0] == "A" else "rev ins"
#         insertion = inverse_mapping[int(mut[1:-1]) - 1].split(":")
#         insertion_position = insertion[0]
#         inserted_nucleotides = insertion[1]
#         gene_position = position_to_gene(int(insertion_position))
#         if gene_position["gene_name"] is None:
#             position_string = ""
#         else:
#             frame_string = f"/{gene_position['reading_frame']}" if gene_position["reading_frame"] != 0 else ""
#             position_string = f" ({gene_position['gene_name']}:{gene_position['codon_number']}{frame_string})"
#         return (indel_type, insertion_position + inserted_nucleotides + position_string)

#     else:
#         frameshift_type = "frameshift" if mut[0] == "A" else "rev frameshift"
#         frameshift = inverse_mapping[int(mut[1:-1]) - 1]
#         return (frameshift_type, frameshift)

# for node in nodes_with_mutations["nodes"].values():
#     list_of_pairs = list(map(mut_to_str, node["muts"]))
#     node["muts"] = defaultdict(list)
#     for a,b in list_of_pairs:
#         node["muts"][a].append(b)
#     node["muts"] = dict(node["muts"])

# for key, value in nodes_with_mutations["nodes"].items():
#     if value["muts"] != {}:
#         print(f"{key}: {value['muts']}")

# merge(result, nodes_with_mutations)

# #     json.dump(result , output)


# # if __name__ == "__main__":
# #     main()
# %%

# %%
