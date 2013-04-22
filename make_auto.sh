#!/bin/bash

. ~/.bashrc

make all_carbon 
create_carbon_network

run_carbon data/my_net.xml data/zone.xml out

make_figure_data.sh

cd idl
idl make_idl_figures
mv *.eps ~
cd ..


