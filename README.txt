Put this small package in nucnet-tools-code/my_projects.

>cd ...my_projects/carbon-dust

>./make all_carbon

To make carbon network data:

>./make carbon_data

To run the carbon oxygen netowrk:

>./run_carbon data/my_net.xml data/zone.xml out

By default, the network is up to n=8 for carbon chain and carbon ring.

To make data for idl plotting, do

>./make_figure_data.sh
