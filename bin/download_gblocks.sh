#! /bin/bash

if [ -f "bin/Gblocks" ]
then
    echo "Gblocks already installed."
else
    wget http://molevol.cmima.csic.es/castresana/Gblocks/Gblocks_Linux64_0.91b.tar.Z
    tar -zxf Gblocks_Linux64_0.91b.tar.Z
    mv Gblocks_0.91b/Gblocks bin/
    rm Gblocks_Linux64_0.91b.tar.Z
    rm -r Gblocks_0.91b
fi
