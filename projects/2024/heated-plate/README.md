# Heated Plate

UMass Lowell Fall 2024 <br>
Dept. of Chemical Engineering, Nuclear Program <br>
Engy-4390: Nuclear Systems and Design Analysis 

View the project on `NBViewer`: [![NBViewer](https://raw.githubusercontent.com/jupyter/design/master/logos/Badges/nbviewer_badge.svg)](https://nbviewer.jupyter.org/github/dpploy/engy-4390/blob/main/projects/2024/heated-plate)

Run the project on `Binder`: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/dpploy/engy-4390/HEAD?filepath=projects%2F2024%2Fheated-plate)

  >**Student:** [Noah Johnson](https://github.com/Noah-R_Johnson) <br>
  >**Mentor/Collaborator:** [Prof. Valmor F. de Almeida](https://github.com/dealmeidavf) <br>
  >[Dept. of Chemical Engineering (Nuclear Energy Program)](https://www.uml.edu/Engineering/Chemical/faculty/de-Almeida-Valmor.aspx) <br>
  >University of Massachusetts Lowell, USA <br>

This is a 1-D model of the fuel plate and cooling channel of the UMass Lowell Research Reactor. This problem help students to think about heat transfer in
multi-domains. The interfacial boundary conditions used are standard. The temperature profile and heat flux are calculated using two different implementations of the 
Finite Element method.

|  |
|:---:|
| <img width="900" src="pics/domain.png" title="pyvista"> |
| <p style="text-align:center;"><b>Heat transfer between a heated plate, a fluid channel, and an unheated plate.</b></p> |

|  |
|:---:|
| <img width="900" src="pics/results.png" title="My result"> |
| <p style="text-align:center;"><b>No-flux Neumann boundary conditions on the left with a heat source in the left plate and an equivalent heat sink in the fluid channel with a Robin boundary condition on the far right.</b></p> |

This Semester I have spent my time improving the projects technical aspects clarifying the data presented and creating a comparison to gold standard provided by Professor Walmor de Almeida. The past groups didn't have this data and with this data future groups will be able to model potentially more tests with different boundary conditions and understanding how they can check their work after they're done. 

References:

 + [Eng-5330: Computational Transport Phenomena: course notes](https://github.com/dpploy/engy-5330)
