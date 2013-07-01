#!/bin/bash

FILE=out

print_properties ${FILE} time t9 rho > idl/prop
print_mass_fractions_in_zones ${FILE} n h1c h2 he2c li3c be4c b5c c6c n7c o8c o8r nn > idl/mass_fractions

#make_flow_data.sh

compute_carbon_flows ${FILE} '' > flows
sed_flow_names.sh
