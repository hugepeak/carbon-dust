Put this small package in nucnet_projects directory.

>cd carbon-dust

To make carbon executables:

>make all_carbon

To make carbon network data:

>make carbon_data

To run the carbon oxygen netowrk:

>./run_carbon data/my_net.xml data/zone.xml out

By default, the network is up to n=8 for carbon chain and carbon ring.

To make the plot:

>cd python
>python plot_figures.py

To run the generic netowrk:

>./run_generic data/generic_net.xml data/generic_zone.xml out

By default, the network is up to n=8 for carbon chain and carbon ring.

