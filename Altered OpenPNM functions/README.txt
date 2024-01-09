Two function scripts need to be altered to allow for the Hydraulic conductance of Larachi: 
	a.	Transport.
		i.	 Location:  OpenPNM – algorithms – transport
	b.	Reactive transport. 
		i.	Location: OpenPNM- algorithms – reactive transport


The easiest way to correctly alter these is by copy-pasting the function scripts defined in this folder into these inherent OpenPNM scripts:
	a. OpenPNM_Transport_altered for the transport function script of OpenPNM.
	b. OpenPNM_Reactive_transport_altered for the reactive transport function script of OpenPNM.

Just select CTRL+A, CTRL+C in the altered function scripts and CTRL+A, CTRL+V in the inherent OpenPNM function scripts. 

By changing to another (random) OpenPNM version via Gitkraken and subsequently returning to OpenPNM version 3.0, the original transport and reactive transport function scripts are retrieved. 
Hence, you dont need to save a copy of them.