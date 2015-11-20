#FILTERING COMPLEXES

Whilst there will undoubtedly be changes, at present the pipeline is as follows:

1. Load tablex.out. This includes all the complexes with chains over length x residues.

2. Get all the mouse chains for these complexes. Only accept a complex if +80% of its chains map to a mouse protein. This is important not just because it increases the amount of available data per complex, but also because if only a few mouse proteins map to a complex, there's a good chance that the complex doesn't exist, or is in a different form in the mouse.

3. For each complex, get all the interface data. For every single interface in the complex, map the pdb chains to individual mouse proteins. This gives a situation where a pair of mouse proteins will often have multiple interfaces associated with it. Of these choose the biggest.

4. Assign PFAM domains to all proteins. If no PFAM domain can't be assigned, give it a unique random number. We don't want to get rid of these proteins because they are more likely to be unique proteins than they are to be paralagous, and data is precious.

5. Cluster all your proteins based on similarity of protein chains. Currently I'm doing this based on DBSCAN algorithm. If 2 complexes share >= 2/3rds (>66%) the same uniprot IDS (based on the larger of the two complexes), then add them to the same cluster.

7. For each cluster, sort in order of the following criteria:
    1. Number of mapped mouse proteins. A possible alternative here is to favour not the total number, but diversity in terms of N/E ratios. U we don't really care about so much.
    2. Number of different interfaces available.
    3. Average sequence identity of the complex.
Then pick the top hit from each cluster.

8. Do a final cluster and pass to remove proteins in which all chains have identical PFAM compositions to the chains of another complex.

9. Remove complexes with almost fully paralagous protein domains (70% similarity) across the length of the whole complex.

10. Remove ribosomes (though perhaps manually keep 40s and 60s subunits?)

11. Remove Immunoglobin proteins from V/C1/C2-set families.

