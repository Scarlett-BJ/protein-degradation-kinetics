#!/usr/bin/env python3

import filtercomp as fc
from ixntools import structure

pool = fc.load_table(1, 2)
pool = fc.get_chains(pool, 80)
# pool = fc.get_interfaces(pool)
fc.get_interfaces2(pool)

# cmap = structure.chain_map()
# # print(cmap['1fs2'])
# # # print(cmap['1go4'])
# for s in pool:
#     if s.name != '5bsa':
#         continue
#     for i in s.interfaces:
#         print(s, i, s.interfaces[i])

# cmap = structure.similar_chains()
# print(cmap)
# for s in cmap:
#     struc = cmap[s]
#     if len(set(struc.keys())) != len(set(struc.values())):
#         print(s, struc)
#         print()