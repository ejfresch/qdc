#!/bin/bash

printf "Installing base packages...\n"
sudo apt-get install python-matplotlib python-qt4 python2.7-dev libxml2 libxml2-dev unzip wget g++ make swig
if [ $? -ne 0 ]
then
exit
fi
printf "[Done]\n"

printf "Downloading libsbml..."
rm -rf libsbml*
wget http://sourceforge.net/projects/sbml/files/libsbml/4.0.1/libsbml-4.0.1-src.zip
if [ $? -ne 0 ]
then
exit
fi
printf "[Done]\n"

printf "Decompressing libsbml..."
unzip libsbml-4.0.1-src.zip
if [ $? -ne 0 ]
then
exit
fi
printf "[Done]\n"

printf "Compiling libsbml..."
mkdir libsbml
if [ $? -ne 0 ]
then
exit
fi

DIR=`pwd`

cd libsbml-4.0.1
if [ $? -ne 0 ]
then
exit
fi



./configure --with-python --prefix="$DIR/libsbml" --with-swig
if [ $? -ne 0 ]
then
exit
fi

make 
if [ $? -ne 0 ]
then
exit
fi

make install
if [ $? -ne 0 ]
then
exit
fi
printf "[Done]\n"

printf "Configuring qdc..."

mv ../libsbml/lib/python* ../libsbml/lib/python

cd ../../

make parser

mkdir src

printf "[Done]\n"
