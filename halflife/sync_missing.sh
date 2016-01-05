#!/bin/bash

for f in $(cat "missing.tmp");
do
    $(rsync marshlab:/media/jonwells/Mmus_coexpressdb/Mmu.v13-01.G20959-S31479.rma.mrgeo.d/$f ../../expression/expression/data/coexpressdb/)
done
