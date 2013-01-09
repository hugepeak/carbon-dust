#!/bin/bash

FILE=out

print_properties ${FILE} time t9 rho > idl/prop
print_mass_fractions_in_zones ${FILE} n h1 h2 he2 li3 be4 b5 c6 n7 o8C o8R > idl/mass_fractions
print_carbon_flows ${FILE} '' "[(product = 'he2' and product = 'gamma') or (product = 'he2' and product = 'n')]" > idl/c2_flows
print_carbon_flows ${FILE} '' "[(product = 'h2' and product = 'gamma') or (reactant = 'h2' and reactant = 'gamma') or (reactant = 'h2' and reactant = 'electron')]" > idl/co_flows

compute_carbon_flows ${FILE} '' > flows
sed_flow_names.sh
