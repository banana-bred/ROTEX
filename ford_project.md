---
project: ROTEX
summary: Calculates electron-impact rotational excitation cross sections for linear and nonlinear molecules in the rigid rotor approximation.
extensions: f
fixed_extensions:
max_frontpage_items: 6
display: public
license: MIT
---

The program (main driver) can be found in main/[[rotex.f(file)]].
It calls a small number of driver routines from [[rotex__drivers]], which do (or call other routines that do) most of the work.

{!README.html!}
