#!/bin/bash

FILE=out

print_properties ${FILE} time t9 rho > idl/prop
print_mass_fractions_in_zones ${FILE} n h1 h2 he2 li3 be4 b5 c6C n7C o8C o8R nn > idl/mass_fractions

make_flow_data.sh

compute_carbon_flows ${FILE} '' > flows
sed_flow_names.sh
