#! /bin/bash

if [ -f "bin/revtrans.py" ]
then
    echo "RevTrans already installed."
else
    wget http://www.cbs.dtu.dk/services/RevTrans/releases/revtrans-1.4.tgz
    tar -zxf revtrans-1.4.tgz 
    mv RevTrans-1.4/* bin/
    rm revtrans-1.4.tgz
    rmdir RevTrans-1.4
fi
