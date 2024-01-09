# Alter OpenPNM functions to run the codes
Two function scripts need to be altered to allow for the Hydraulic conductance simulations: 
	a.	Transport.
		Location:  OpenPNM – algorithms – transport
	b.	Reactive transport. 
		Location: OpenPNM- algorithms – reactive transport

The easiest way to correctly alter these is by copy-pasting the function scripts defined in this folder into these inherent OpenPNM scripts. 
By changing to another (random) OpenPNM version via Gitkraken and subsequently returning to OpenPNM v3.0.0, the original transport and reactive transport function scripts are retrieved. 
